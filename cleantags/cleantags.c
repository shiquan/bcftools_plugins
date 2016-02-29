#include <stdio.h>
#include <stdlib.h>
#include <htslib/vcf.h>

bcf_hdr_t *header = NULL;

const char * about(void)
{
    return
	"Simple clean all the INFO tags and FMT tags except GT.\n";
}

static void remove_hdr_lines(bcf_hdr_t *hdr, int type)
{
    int i = 0, nrm = 0;
    while ( i<hdr->nhrec )
    {
        if ( hdr->hrec[i]->type!=type ) { i++; continue; }
        bcf_hrec_t *hrec = hdr->hrec[i];
        if ( type==BCF_HL_FMT )
        {
            // everything except FORMAT/GT
            int id = bcf_hrec_find_key(hrec, "ID");
            if ( id>=0 && !strcmp(hrec->vals[id],"GT") ) { i++; continue; }
        }
        nrm++;
        hdr->nhrec--;
        if ( i < hdr->nhrec )
            memmove(&hdr->hrec[i],&hdr->hrec[i+1],(hdr->nhrec-i)*sizeof(bcf_hrec_t*));
        bcf_hrec_destroy(hrec);
    }
    if ( nrm ) bcf_hdr_sync(hdr);
}

int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    header = in;
    bcf_hdr_t *out_hdr = out;
    remove_hdr_lines(out_hdr, BCF_HL_FLT);
    remove_hdr_lines(out_hdr, BCF_HL_INFO);
    remove_hdr_lines(out_hdr, BCF_HL_FMT);
    return 0;
}


bcf1_t *process(bcf1_t *rec)
{
    int i;
    bcf_unpack(rec, BCF_UN_ALL);
    // remove filter
    bcf_update_filter(header, rec, NULL, 0);
    // remove quality
    bcf_float_set_missing(rec->qual);

    // remove info
    for (i=0; i<rec->n_info; i++)
    {
        bcf_info_t *inf = &rec->d.info[i];
        if ( inf->vptr_free )
        {
            free(inf->vptr - inf->vptr_off);
            inf->vptr_free = 0;
        }
        rec->d.shared_dirty |= BCF1_DIRTY_INF;
        inf->vptr = NULL;
    }

    // remove format
    for (i=0; i<rec->n_fmt; i++)
    {
        bcf_fmt_t *fmt = &rec->d.fmt[i];
        const char *key = bcf_hdr_int2id(header,BCF_DT_ID,fmt->id);
        if ( key[0]=='G' && key[1]=='T' && !key[2] ) continue;

        if ( fmt->p_free )
        {
            free(fmt->p - fmt->p_off);
            fmt->p_free = 0;
        }
        rec->d.indiv_dirty = 1;
        fmt->p = NULL;
    }
    return rec;
}

void destroy(void)
{
    return;
}
