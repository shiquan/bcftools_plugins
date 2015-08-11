/*   select.c  --
 *   This is a plugin of bcftools for selecting fields in defined format from VCF/BCF file
 *
 *   The extracting rules are most similar with `bcftools-query -f`, but still
 *   there are two different rules in this plugin. 
 *   1) each sample can be exported per line by specify %SAMPLE in the INFO region
 *   2) split entry by specified tag  '--split NONE,ALT,GT,SAMPLE,HGVS,ALL'
 *
 *   Copyright (C) 2015 BGI Research
 *
 *   Author : SHI Quan (shiquan@genomics.cn)
 *
 *   License : MIT
 *
 *   This program is inspired by pd3's vcfquery.c .
 *
 * 
 * Demo:
 *
 * bcftools +select -f '%CHROM\t%POS\t%REF[\t%TGT\t%DP]' demo.vcf.gz
 *
 * bcftools +select -f '%BED\t%SAMPLE\t%REF\t%ALT\t[%DP]' demo.vcf.gz
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <getopt.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <ctype.h>
#include "bcftools.h"

#define SPLIT_NONE    1          //  bcftools query mode, one bcf1_t per line
#define SPLIT_TAG     0
#define SPLIT_SAMPLE  2          //  flag for spliting by sample
#define SPLIT_ALT     4          //  flag for split by allele, each alt per line
#define SPLIT_GT      2          //  flag for split by genotype (0/0, 0/1, 1/1), it's same with split by samples ^ SPLIT_ALT
#define SPLIT_TRANS    (1<<4)    //  flag for split by TRANS name, if use this flag any tags related with transcripts should be split 
#define SPLIT_DEFAULT 0
#define SPLIT_ALL     ( SPLIT_SAMPLE | SPLIT_ALT | SPLIT_TRANS )

#define IS_SEP    1
#define IS_FIX    1
#define IS_SAM    2
#define IS_ALT    4
#define IS_GT     (1<<3)
#define IS_TRANS   (1<<4)

#define S_CHROM     1
#define S_POS       2
#define S_ID        3
#define S_REF       4
#define S_ALT       5
#define S_BED       6 // %BED is only used with type & SPLIT_ALT , bed calculated from ref and alt position
#define S_QUAL      7
#define S_FILTER    8
#define S_INFO      9
#define S_FORMAT    10
#define S_SAMPLE    11
#define S_SEP       12
#define S_TYPE      13
#define S_GT        14 // GT is useless for us, export TGT directly
#define S_TGT       14
#define S_IUPAC_GT  15
#define S_FIRST_ALT 16
#define S_TRANS      17

#define MEMPOOL     2621440

/* split flag */
static int split_flag = SPLIT_DEFAULT;

bcf_hdr_t * header;

/* the fmt_t and convert_t are adapted from pd3's convert.c 
 * if we want export the sample info and allele info per line
 * we must rewrite the process_* functions in convert.c
 */
typedef struct _convert convert_t;
typedef struct _tags     tags_t;

typedef struct _fmt
{
    int id, type, is_gtf;
    char *key;
    bcf_fmt_t *fmt;
    void (*handler)(bcf1_t *, struct _fmt *, int iala, int isample, tags_t *);
}
fmt_t;

struct _convert
{
    int mfmt, nfmt;
    int max_unpack;
    fmt_t *fmt;
    char *format_str;
};

typedef struct
{
    int type;
    int m, n;    // the array size of a[], usually $m==1
    char **a; 
}
mval_t;

struct _tags
{
    int n, m; // n : fields , m : genotypesn*samples
    int k, l; // the sizes of matrix
    mval_t **trans;
};

static int skip_ref = 0;
static int print_header = 0;
void clear_tag(tags_t *t)
{
    int i, j, k;
    for (i=0; i<t->n; ++i)
    {
	for (j=0; j<t->m; ++j)
	{
	    for (k=0; k < t->trans[i][j].m; ++k) free(t->trans[i][j].a[k]);
	}
	free(t->trans[i]);
    }
    free(t->trans);
}
int tags2str(tags_t *t, kstring_t *s)
{
    int i, j, k, d;
    int l_ori = s->l;
    for (i=0; i<t->n; ++i) // allele
    {
	d=1; k=0;
	for (; k<d;k++) // trans
	{
	    for (j=0; j<t->m; ++j) // rows
	    {
		mval_t f= t->trans[i][j];
		
		if ((f.type & IS_TRANS) && f.m > 1)
		{
		    if (d != 1 && d != f.m) error("ERROR: TRANS tags have different transcripts!\n");
		    d=f.m;
		}
		if (f.type & IS_TRANS)
		{
		    kputs(f.a[k], s);
		}
		else
		{
		    fprintf(stderr, "%s\n", f.a[0]);
		    kputs(f.a[0], s);
		}
	    }
	    kputc('\n', s);
	}
    }
    return s->l -l_ori;
}

static convert_t *convert = NULL;

int convert_header(kstring_t *str)
{
    int l_ori = str->l;
    int i;
    str->l=0;
    for ( i=0; i<convert->nfmt; ++i )
    {
	if ( convert->fmt[i].is_gtf && !(split_flag & SPLIT_SAMPLE ))
	{
	    int j = i, js, k;
	    while ( convert->fmt[j].is_gtf ) j++;
	    for ( js = 0; js < bcf_hdr_nsamples(header); ++js )
	    {
		for ( k=i; k<j; ++k )
		{
		    if ( convert->fmt[k].type & IS_SEP )
		    {
			if ( convert->fmt[k].key ) kputs( convert->fmt[k].key, str);
		    }
		    else
			ksprintf(str, "%s:%s", header->samples[js], convert->fmt[k].key);
		}
	    }
	    i = j -1; // skip sample fields
	    continue;
	}
	if ( convert->fmt[i].type & IS_SAM ) {
	    if ( !(split_flag & SPLIT_SAMPLE) ) error("split tag has inconsistent format\n");
	    //if ( nsam > 1 ) error("more than 1 sample name field\n");
	    //else nsam++;
	    kputs("SAMPLE", str);
	    continue;
	}
	// Fixed fields
	if ( convert->fmt[i].key ) kputs(convert->fmt[i].key, str);
    }    
    return str->l - l_ori;
}
void init_tags(tags_t *t, bcf1_t *line)
{
    int i, j;
    t->l = split_flag & SPLIT_ALT ? line->n_allele : 1;
    if ( split_flag & SPLIT_SAMPLE ) t->l *= header->n[BCF_DT_SAMPLE];
    t->k = convert->nfmt;
    t->trans = (mval_t**)calloc(t->l, sizeof(mval_t*));
    for (i=0; i<t->l; ++i)
    {
	t->trans[i] = (mval_t*)calloc(t->k, sizeof(mval_t));
	for (j=0; j<t->k; ++j)
	{
	    t->trans[i][j].m = 0;
	    t->trans[i][j].n = 1;
	    t->trans[i][j].type = 0;
	    t->trans[i][j].a = (char**)calloc(t->trans[i][j].n, sizeof(char*));
	}
    }
    for (i=0; i<convert->nfmt; ++i)
    {
	if (convert->fmt[i].id == -1) continue;
	for (j=0; j<(int)line->n_fmt; j++)
	{
	    if ( line->d.fmt[j].id==convert->fmt[i].id ) { convert->fmt[i].fmt = &line->d.fmt[j]; break; }
	}
    }
}
int convert_line(bcf1_t *line, kstring_t *str)
{
    int l_ori = str->l;
    bcf_unpack(line, convert->max_unpack);
    int i, k=0;
    tags_t tag = { 0, 0, 0};
    init_tags(&tag, line);
    int isample=-1, nsample=0;
    int row = 0;
    if ( split_flag & SPLIT_SAMPLE )
    {
	isample = 0;
	nsample = header->n[BCF_DT_SAMPLE];
    }
    for ( ; isample<nsample; ++isample )
    {
	int iala=-1;
	if ( split_flag & SPLIT_ALT )
	{
	    bcf_fmt_t *fgt = bcf_get_fmt(header, line, "GT");
#define BRANCH(type_t, vector_end) do {					\
		type_t *ptr = (type_t*) (fgt->p + isample*fgt->size);	\
		iala = ptr[k]&1;					\
	    } while(0)
	    fprintf(stderr, "geno: %d\n", fgt->n);
	    for ( k=0; k<fgt->n; ++k)
	    {
		switch (fgt->type) {
		    case BCF_BT_INT8:  BRANCH(int8_t, bcf_int8_vector_end); break;
		    case BCF_BT_INT16: BRANCH(int16_t, bcf_int16_vector_end); break;
		    case BCF_BT_INT32: BRANCH(int32_t, bcf_int32_vector_end); break;
		    default: fprintf(stderr, "FIXME: type %d in bcf_format_gt?\n", fgt->type); abort(); break;
		}
		if ( skip_ref && !iala ) continue;
		tag.m = row++;
		for ( i=0; i<convert->nfmt; i++ )
		{
		    tag.n = i;
		    fprintf(stderr,"m: %d\tn:%d\tl:%d\tk:%d\n",tag.m, tag.n, tag.l, tag.k);
		    if ( convert->fmt[i].handler ) convert->fmt[i].handler(line, &convert->fmt[i], iala, isample, &tag);
		    //fprintf(stderr, "key: %s id:%d\n", tag.trans[k][i].a[0], convert->fmt[i].id);
		}
	    }
#undef BRANCH
	}
	else
	{
	    tag.m = row++;	    
	    for ( i=0; i<convert->nfmt; i++ )
	    {
		tag.n = i;
		//fprintf(stderr,"m: %d\tn:%d\tl:%d\tk:%d\n",tag.m, tag.n, tag.l, tag.k);
		if ( convert->fmt[i].handler ) convert->fmt[i].handler(line, &convert->fmt[i], iala, isample, &tag);
	    }
	}
    }
    tags2str(&tag, str);
    //fprintf(stderr, "l: %d\tm: %d\t%s\n", str->l, str->m, str->s);
    clear_tag(&tag);
    return str->l - l_ori;
}

void free_list(char **s, int n)
{
    int i;
    for (i = 0; i < n; ++i) free(s[i]);
    free(s);
}
int init_split_flag(char *s)
{
    int n, i;
    int flag = 0;
    char **list = hts_readlist(s, 0, &n);

    for (i = 0; i < n; ++i)
    {
	if (strcmp(list[i], "NONE") == 0) { free_list(list, n); return SPLIT_NONE; }
	else if (strcmp(list[i], "ALT") == 0) flag |= SPLIT_ALT;
	else if (strcmp(list[i], "SAMPLE") == 0) flag |= SPLIT_SAMPLE;
	else if (strcmp(list[i], "TRANS") == 0) flag |= (SPLIT_TRANS|SPLIT_ALT);
	else if (strcmp(list[i], "ALL") == 0) flag |= SPLIT_ALL;
	else {
	    fprintf(stderr, "cannot recongize this tag : %s\n", list[i]);
	    exit(1);
	}
    }
    free_list(list, n);
    return flag;
}
static void process_first_alt(bcf1_t *line, fmt_t *fmt, int iala, int isample, tags_t *tag)
{
    kstring_t str = { 0, 0, 0 };
    mval_t *pv= &tag->trans[tag->m][tag->n];
    pv->type = fmt->type;
    if ( line->n_allele == 1 ) kputc('.', &str);
    else kputs(line->d.allele[1], &str);
    if (pv->n == pv->m)
    {
	pv->n+= 2;
	pv->a = (char**)realloc(pv->a, pv->n*sizeof(char*));
    }
    pv->a[pv->m++] = str.s;

}
static void process_chrom(bcf1_t *line, fmt_t *fmt, int iala, int isample, tags_t *tag)
{
    kstring_t str = { 0, 0, 0 };
    kputs(header->id[BCF_DT_CTG][line->rid].key, &str);
    mval_t *pv= &tag->trans[tag->m][tag->n];
    pv->type = fmt->type;
    if (pv->n == pv->m)
    {
	pv->n+= 2;
	pv->a = (char**)realloc(pv->a, pv->n*sizeof(char*));
    }
    pv->a[pv->m++] = str.s;

}
static void process_pos(bcf1_t *line, fmt_t *fmt, int iala, int isample, tags_t *tag)
{
    kstring_t str = { 0, 0, 0 };
    kputw(line->pos+1, &str);
    mval_t *pv = &tag->trans[tag->m][tag->n];
    pv->type = fmt->type;
    if (pv->n == pv->m)
    {
	pv->n+= 2;
	pv->a = (char**)realloc(pv->a, pv->n*sizeof(char*));
    }
    pv->a[pv->m++] = str.s;
}
static void process_bed(bcf1_t *line, fmt_t *fmt, int iala, int isample, tags_t *tag)
{
    
}
static void process_ref(bcf1_t *line, fmt_t *fmt, int iala, int isample, tags_t *tag)
{
    kstring_t str = { 0, 0, 0 };
    kputs(line->d.allele[0], &str);
    mval_t *pv= &tag->trans[tag->m][tag->n];
    pv->type = fmt->type;
    if (pv->n == pv->m)
    {
	pv->n+= 2;
	pv->a = (char**)realloc(pv->a, pv->n*sizeof(char*));
    }
    pv->a[pv->m++] = str.s;

}
static void process_id(bcf1_t *line, fmt_t *fmt, int iala, int isample, tags_t *tag)
{
}
static void process_alt(bcf1_t *line, fmt_t *fmt, int iala, int isample, tags_t *tag)
{
    kstring_t str = { 0, 0, 0 };
    mval_t *pv= &tag->trans[tag->m][tag->n];
    pv->type = fmt->type;
    if ( line->n_allele == 1 ) kputc('.', &str);
    else if ( iala == -1 )
    {
	int i;
	for (i=1; i<line->n_allele; ++i)
	{
	    if ( i>1 ) kputc(',', &str);
	    kputs(line->d.allele[i], &str);
	}
    }
    else kputs(line->d.allele[iala], &str);
    if (pv->n == pv->m)
    {
	pv->n+= 2;
	pv->a = (char**)realloc(pv->a, pv->n*sizeof(char*));
    }
    pv->a[pv->m++] = str.s;

}
static void process_qual(bcf1_t *line, fmt_t *fmt, int iala, int isample, tags_t *tag)
{
    kstring_t str = { 0, 0, 0 };
    if ( bcf_float_is_missing(line->qual) ) kputc('.', &str);
    else ksprintf(&str, "%g", line->qual);
    mval_t *pv= &tag->trans[tag->m][tag->n];
    pv->type = fmt->type;
        if (pv->n == pv->m)
    {
	pv->n+= 2;
	pv->a = (char**)realloc(pv->a, pv->n*sizeof(char*));
    }
    pv->a[pv->m++] = str.s;

}
static void process_filter(bcf1_t *line, fmt_t *fmt, int iala, int isample, tags_t *tag)
{
    kstring_t str = { 0, 0, 0 };
    int i;
    if ( line->d.n_flt )
    {
	for (i=0; i<line->d.n_flt; i++)
	{
	    if (i) kputc(';', &str);
	    kputs(header->id[BCF_DT_ID][line->d.flt[i]].key, &str);
	}
    }
    else kputc('.', &str);
    
    mval_t *pv= &tag->trans[tag->m][tag->n];
    pv->type = fmt->type;
        if (pv->n == pv->m)
    {
	pv->n+= 2;
	pv->a = (char**)realloc(pv->a, pv->n*sizeof(char*));
    }
    pv->a[pv->m++] = str.s;

    
}
static void process_tgt(bcf1_t *line, fmt_t *fmt, int iala, int isample, tags_t *tag)
{
    if (iala != -1 ) error("%TGT only used with SPLIT_GT \n");
    kstring_t str = { 0, 0, 0 };
    if ( fmt->fmt==NULL ) kputc('.', &str);
    else
    {
	assert(fmt->fmt->type==BCF_BT_INT8);
	int l;
	int8_t *x = (int8_t*)(fmt->fmt->p + isample*fmt->fmt->size);
	for (l=0; l<fmt->fmt->n && x[l]!=bcf_int8_vector_end; ++l)
	{
	    if (l) kputc("/|"[x[l]>>1], &str);
	    if (x[l]>>1) kputs(line->d.allele[(x[l]>>1)-1], &str);
	    else kputc('.', &str);
	}
	if (l==0) kputc('.', &str);
    }
    mval_t *pv= &tag->trans[tag->m][tag->n];
    pv->type = fmt->type;
    if (pv->n == pv->m)
    {
	pv->n+= 2;
	pv->a = (char**)realloc(pv->a, pv->n*sizeof(char*));
    }
    pv->a[pv->m++] = str.s;

}
static void process_trans(bcf1_t *line, fmt_t *fmt, int iala, int isample, tags_t *tag)
{
}
static void process_info(bcf1_t *line, fmt_t *fmt, int iala, int isample, tags_t *tag)
{
    kstring_t str = { 0, 0, 0 };
    mval_t *pv= &tag->trans[tag->m][tag->n];
    pv->type = fmt->type;
    int i;
    for (i=0; i<line->n_info; i++)
	if ( line->d.info[i].key == fmt->id ) break;
    if ( i==line->n_info )
    {
	kputc('.', &str);
	if (pv->n == pv->m)
	{
	    pv->n+= 2;
	    pv->a = (char**)realloc(pv->a, pv->n*sizeof(char*));
	}
	pv->a[pv->m++] = str.s;

	return;
    }

    bcf_info_t *info = &line->d.info[i];
    if ( info->len < 0)
    {
	kputc('1', &str);
	if (pv->n == pv->m)
	{
	    pv->n+= 2;
	    pv->a = (char**)realloc(pv->a, pv->n*sizeof(char*));
	}
	pv->a[pv->m++] = str.s;

	return;
    }

    if ( info->len == 1 )
    {
	switch (info->type)
	{
	case BCF_BT_INT8: if ( info->v1.i==bcf_int8_missing ) kputc('.', &str); else kputw(info->v1.i, &str); break;
	case BCF_BT_INT16: if ( info->v1.i==bcf_int16_missing ) kputc('.', &str); else kputw(info->v1.i, &str); break;
	case BCF_BT_INT32: if ( info->v1.i==bcf_int32_missing ) kputc('.', &str); else kputw(info->v1.i, &str); break;
	case BCF_BT_FLOAT: if ( bcf_float_is_missing(info->v1.f) ) kputc('.', &str); else ksprintf(&str, "%g", info->v1.f); break;
	case BCF_BT_CHAR: kputc(info->v1.i, &str); break;
	default: fprintf(stderr,"todo: type %d\n", info->type); exit(1); break;
	}
    }
    else
	bcf_fmt_array(&str, info->len, info->type, info->vptr);
    if (pv->n == pv->m)
    {
	pv->n+= 2;
	pv->a = (char**)realloc(pv->a, pv->n*sizeof(char*));
    }
    pv->a[pv->m++] = str.s;

}
static void process_format(bcf1_t *line, fmt_t *fmt, int iala, int isample, tags_t *tag)
{
    kstring_t str = { 0, 0, 0};
    mval_t *pv= &tag->trans[tag->m][tag->n];
    pv->type = fmt->type;

    if ( fmt->fmt == NULL ) kputc('.', &str);
    else
	bcf_fmt_array(&str, fmt->fmt->n, fmt->fmt->type, fmt->fmt->p + isample*fmt->fmt->size);
    if (pv->n == pv->m)
    {
	pv->n+= 2;
	pv->a = (char**)realloc(pv->a, pv->n*sizeof(char*));
    }
    pv->a[pv->m++] = str.s;


}
static void process_sample(bcf1_t *line, fmt_t *fmt, int iala, int isample, tags_t *tag)
{
    assert(isample > -1);
    kstring_t str = { 0, 0, 0};
    mval_t *pv= &tag->trans[tag->m][tag->n];
    pv->type = fmt->type;
    kputs(header->samples[isample], &str);
        if (pv->n == pv->m)
    {
	pv->n+= 2;
	pv->a = (char**)realloc(pv->a, pv->n*sizeof(char*));
    }
    pv->a[pv->m++] = str.s;

}
static void process_sep(bcf1_t *line, fmt_t *fmt, int iala, int isample, tags_t *tag)
{
    if (!fmt->key) error("FIXME: This is an empty fmt\n");
    kstring_t str = { 0, 0, 0 };
    mval_t *pv= &tag->trans[tag->m][tag->n];
    pv->type = fmt->type;
    kputs(fmt->key,&str);
    if (pv->n == pv->m)
    {
	pv->n+= 2;
	pv->a = (char**)realloc(pv->a, pv->n*sizeof(char*));
    }
    pv->a[pv->m++] = str.s;

}

/* The register_tag and parse_tag functions are adapted from pd3's convert.c 
 * use S-tags instead of T-tags 
 */
fmt_t *register_tag(int type, char *key, int is_gtf)
{
    convert->nfmt++;
    if (convert->nfmt == convert->mfmt)
    {
	convert->mfmt += 10;
	convert->fmt = (fmt_t*)realloc(convert->fmt, convert->mfmt*sizeof(fmt_t));
    }
    fmt_t *fmt = &convert->fmt[ convert->nfmt-1 ];
    fmt->type = type;
    fmt->is_gtf = is_gtf;
    if ( key == NULL ) error("BUGS: empty key\n");
    fmt->key = strdup(key);
    if (fmt->type == S_SEP)
    {
	fmt->id = -1;
    }
    else
    {
	fmt->id = bcf_hdr_id2int(header, BCF_DT_ID, fmt->key);
	if ( fmt->id == -1 && !fmt->is_gtf)
	{
	    if ( !strcmp("CHROM", key) ) fmt->type = S_CHROM;
	    else if ( !strcmp("POS", key) ) fmt->type = S_POS;
	    else if ( !strcmp("BED", key) ) fmt->type = S_BED;
	    else if ( !strcmp("ID", key) ) fmt->type = S_ID;
	    else if ( !strcmp("REF", key) ) fmt->type = S_REF;
	    else if ( !strcmp("ALT", key) ) fmt->type = S_ALT;
	    else if ( !strcmp("FIRST_ALT", key) ) fmt->type = S_FIRST_ALT;
	    else if ( !strcmp("QUAL", key) ) fmt->type = S_QUAL;
	    else if ( !strcmp("FILTER", key) ) fmt->type = S_FILTER;
	    else if ( !strcmp("SAMPLE", key) ) fmt->type = S_SAMPLE;
	    else  error("No such tag in the header %d\n", key);
    	}
	else
	{
	    if ( !strcmp("HGVS", key) ) fmt->type = S_TRANS;
	    else if ( !strcmp("SAMPLE", key) ) fmt->type = S_SAMPLE;
	    else if ( fmt->id >= 0 && bcf_hdr_idinfo_exists(header, BCF_HL_INFO, fmt->id) )
	    {
    		fmt->type = S_INFO;
    		fprintf(stderr, "Warning: Assuming INFO %s\n", key);
    	    }
	    else
	    {
		fmt->type = S_FORMAT;
	    }
	}
    }
    switch ( fmt->type )
    {
    case S_FIRST_ALT: fmt->handler = &process_first_alt; break;
    case S_CHROM: fmt->handler = &process_chrom; break;
    case S_POS: fmt->handler = &process_pos; break;
    case S_BED: fmt->handler = &process_bed; break;
    case S_ID: fmt->handler = &process_id; break;
    case S_REF: fmt->handler = &process_ref; break;
    case S_ALT: fmt->handler = &process_alt; break;
    case S_QUAL: fmt->handler = &process_qual; break;
    case S_FILTER: fmt->handler = &process_filter; convert->max_unpack |= BCF_UN_FLT; break;
    case S_INFO: fmt->handler = &process_info; convert->max_unpack |= BCF_UN_INFO; break;	
    case S_TRANS: fmt->handler = &process_trans; convert->max_unpack |= BCF_UN_INFO; break;
    case S_FORMAT: fmt->handler = &process_format; convert->max_unpack |= BCF_UN_FMT; break;
    case S_SAMPLE: fmt->handler = &process_sample; break;
    case S_SEP: fmt->handler = &process_sep; break;
    case S_TGT: fmt->handler = &process_tgt; convert->max_unpack |= BCF_UN_FMT; break;
    default: error("TODO: handler for type %d\n", fmt->type);
    }
    /* if ( key ) */
    /* { */
    /* 	if ( fmt->type==S_INFO ) */
    /* 	{ */
    /* 	    fmt->id = bcf_hdr_id2int(header, BCF_DT_ID, key); */
    /* 	    if ( fmt->id==-1 ) error("Error: no such tag defined in the VCF header: INFO/%s\n", key); */
    /* 	} */
    /* } */
    return fmt;
}
static char *parse_tag(char *p, int is_gtf)
{
    char *q = ++p;
    while ( *q && (isalnum(*q) || *q=='_' || *q=='.')) q++;
    kstring_t str = { 0, 0, 0};
    if ( q-p==0 ) error("Could not parse format string: %s\n", convert->format_str);
    kputsn(p, q-p, &str);
    if ( is_gtf )
    {
	if ( !strcmp(str.s, "SAMPLE") )
	{
	    if ( !(split_flag & SPLIT_SAMPLE) ) register_tag(S_SAMPLE, "SAMPLE", is_gtf);
	    else error("SAMPLE tag should be defined in FMT region, if split by samples");
	}
	else if ( !strcmp(str.s, "GT") || !strcmp(str.s, "TGT")) register_tag(S_TGT, "GT", is_gtf);
	else if ( !strcmp(str.s, "IUPACGT") ) register_tag(S_IUPAC_GT, "GT", is_gtf);
	else
	{
	    register_tag(S_FORMAT, str.s, is_gtf);
	}
    }
    else
    {
	if ( !strcmp(str.s, "CHROM") ) register_tag(S_CHROM, "CHROM", is_gtf);
	else if ( !strcmp(str.s, "POS") ) register_tag(S_POS, "POS", is_gtf);
	else if ( !strcmp(str.s, "BED") ) register_tag(S_BED, "BED", is_gtf);
	else if ( !strcmp(str.s, "POS") ) register_tag(S_CHROM, "POS", is_gtf);
	else if ( !strcmp(str.s, "REF") ) register_tag(S_REF, "REF", is_gtf);
	else if ( !strcmp(str.s, "ALT") ) register_tag(S_ALT, "ALT", is_gtf);
	else if ( !strcmp(str.s, "HGVS") ) register_tag(S_TRANS, "HGVS", is_gtf);
	else if ( !strcmp(str.s, "FUNC") ) register_tag(S_TRANS, "FUNC", is_gtf);
	else if ( !strcmp(str.s, "SAMPLE") ) { split_flag |= SPLIT_SAMPLE; register_tag(S_SAMPLE, "SAMPLE", is_gtf); }
	else if ( !strcmp(str.s, "FIRST_ALT") ) register_tag(S_FIRST_ALT, "FIRST_ALT", is_gtf);
	else if ( !strcmp(str.s, "QUAL") ) register_tag(S_QUAL, "QUAL", is_gtf);
	else if ( !strcmp(str.s, "INFO") )
	{
	    if ( *q=='/' ) error("Could not parse format string: %s\n", convert->format_str);
	    p = ++q;
	    str.l = 0;
	    while ( *q && (isalnum(*q) || *q=='_' || *q=='.') ) q++;
	    if ( q-p==0 ) error("Could not parse format string: %s\n", convert->format_str);
	    kputsn(p, q-p, &str);
	    register_tag(S_INFO, str.s, is_gtf);
	}
	else
	{
	    register_tag(S_INFO, str.s, is_gtf);
	}
    }
    free(str.s);
    return q;
}
static char *parse_sep(char *p, int is_gtf)
{
    char *q = p;
    kstring_t str = { 0, 0, 0};
    while ( *q && *q!='[' && *q!=']' && *q!='%' )
    {
	if ( *q=='\\' )
	{
	    q++;
	    if ( *q=='n' ) kputc('\n', &str);
	    else if ( *q == 't') kputc('\t', &str);
	    else kputc(*q, &str);
	}
	else kputc(*q, &str);
	q++;
    }
    if ( !str.l ) error("Could not parse format string: %s\n", convert->format_str);
    register_tag(S_SEP, str.s, is_gtf);
    free(str.s);
    return q;
}
void convert_init(char *s)
{
    convert = (convert_t*)malloc(sizeof(convert_t));
    convert->format_str = strdup(s);
    convert->nfmt = 0;
    convert->mfmt = 2;
    convert->fmt = (fmt_t*)calloc(convert->mfmt, sizeof(fmt_t)); 
    convert->max_unpack = 0;
    int is_gtf = 0;
    char *p = convert->format_str;
    while ( *p )
    {
	switch (*p)
	{
	    case '[': is_gtf = 1; p++; break;
	    case ']': is_gtf = 0; p++; break;
	    case '%': p = parse_tag(p, is_gtf); break;
	    default:  p = parse_sep(p, is_gtf); break;
	}
    }
}
const char *about(void)
{
    return "Select tags in pre-defined format.\n";
}

const char *usage(void)
{
    return
	"\n"
	"About : Select tags from VCF/BCF file.\n"
	"Usage:\n"
	"Standalone mode:\n"
	"\tbcfselect [Options] in.vcf.gz\n"
	"Options:\n"
	"\t-f, --format   see man page for deatils.\n"
	"BCFtools plugin mode:\n"
	"\tbcftools +select [General Options] -- [Plugin Options]\n"
	"General Options:\n"
	"\trun \"bcftools plugin\" for a list of common options.\n"
	"Plugin Options: same with standalone mode.\n"
	"\n"
	"Website:\n";
}

kstring_t *mempool = NULL;

int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    static struct option const long_opts[] =
	{
	    {"format", required_argument, NULL, 'f'},
	    {"skip-ref", no_argument, NULL, 'r'},
	    {"split", required_argument, NULL, 's'},
	    {0, 0, 0, 0}
	};
    char c;
    char * format = NULL;
    while ((c = getopt_long(argc, argv, "f:sh?", long_opts, NULL)) >= 0)
    {
	switch (c)
	{
	case 'f': format = strdup(optarg); break;
	case 'r': skip_ref = 1; break;
	case 's': split_flag = init_split_flag(optarg); break;
	case 'h':
	case '?':
	default: fprintf(stderr, "%s", usage()); break;
	}
    }
    header = in;
    convert_init(format);
    free(format);
    mempool = (kstring_t*)malloc(sizeof(kstring_t));
    mempool->m = mempool->l = 0;
    print_header = 1;
    return 0;
}


bcf1_t *process(bcf1_t *rec)
{
    convert_line(rec, mempool);
    //fprintf(stderr, "%s\n", mempool->s);
    mempool->l = 0;
    return 0;
}

void destroy(void)
{
}
#ifdef _SELECT_MAIN

int main(int argc, char **argv)
{
    return run(argc, argv);
}
#endif
