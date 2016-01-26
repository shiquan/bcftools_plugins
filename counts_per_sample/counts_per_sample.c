#include <stdio.h>
#include <stdlib.h>
#include <htslib/vcf.h>

int nsamples = 0;

/*
    This short description is used to generate the output of `bcftools plugin -l`.
*/
const char *about(void)
{
    return
    "A minimal plugin which counts number of non-ref positions for each position in FORMAT level,\n"
    "and average depth and average genotype quality score for these positions";
}

/*
    Called once at startup, allows to initialize local variables.
    Return 1 to suppress VCF/BCF header from printing, 0 otherwise.
*/
unsigned int *ncounts = NULL;
unsigned int *depths = NULL;
unsigned int *gqualitys = NULL;
int gt_id = -1;
int dp_id = -1;
int gq_id = -1;
bcf_hdr_t *header = NULL;

int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    nsamples = bcf_hdr_nsamples(in);
    ncounts = (unsigned int*)calloc(nsamples, sizeof(unsigned int));
    depths = (unsigned int*)calloc(nsamples, sizeof(unsigned int));
    gqualitys = (unsigned int*)calloc(nsamples, sizeof(unsigned int));
    memset(ncounts, 0, nsamples);
    memset(depths, 0, nsamples);
    memset(gqualitys, 0, nsamples);
    gt_id = bcf_hdr_id2int(in, BCF_DT_ID, "GT");
    dp_id = bcf_hdr_id2int(in, BCF_DT_ID, "DP");
    gq_id = bcf_hdr_id2int(in, BCF_DT_ID, "GQ");

    if (!bcf_hdr_idinfo_exists(in, BCF_HL_FMT, gt_id) ) {
	fprintf(stderr, "No GT in the header!\n");
	exit(-1);
    }
    if (!bcf_hdr_idinfo_exists(in, BCF_HL_FMT, dp_id) ) {
	fprintf(stderr, "No FMT/DP in the header!\n");
	exit(-1);
    }
    if (!bcf_hdr_idinfo_exists(in, BCF_HL_FMT, gq_id) ) {
	fprintf(stderr, "No FMT/GQ in the header!\n");
	exit(-1);
    }
	
    header = in;
    return 1;
}


/*
    Called for each VCF record. Return rec to output the line or NULL
    to suppress output.
*/
bcf1_t *process(bcf1_t *rec)
{
    int type = bcf_get_variant_types(rec);
    //assert(fmt1->type == BCF_BT_INT8);
    if ( type ^ VCF_REF ) {
	    bcf_fmt_t *fmt = bcf_get_fmt_id(rec, gt_id);
	    bcf_fmt_t *fmt1 = bcf_get_fmt_id(rec, dp_id);
	    bcf_fmt_t *fmt2 = bcf_get_fmt_id(rec, gq_id);
	    if (fmt2 == NULL) return NULL;
	    assert(fmt->type == BCF_BT_INT8);
	    assert(fmt2->type == BCF_BT_INT8);
	    int i;
	    for (i=0; i<nsamples; ++i) {
		int8_t *x = (int8_t*)(fmt->p + i*(fmt->size));
		int8_t *x2 = (int8_t*)(fmt2->p + i*(fmt2->size));	    
	    // stat snps
	    if (fmt1->type == BCF_BT_INT8) {
		int8_t *x1 = (int8_t*)(fmt1->p + i*(fmt1->size));

		int j;
		for (j=1; j<fmt->n && x[j] != bcf_int8_vector_end; ++j) {
		    if (x[j]>>1 > 1) {
			ncounts[i]++;
			depths[i] += (unsigned int)x1[0];
			gqualitys[i] += (unsigned int)x2[0];
			break;
		    }
		}
	    } else {
		assert(fmt1->type == BCF_BT_INT16);
		int16_t *x1 = (int16_t*)(fmt1->p + i*(fmt1->size));
		int j;
		for (j=1; j<fmt->n && x1[j] != bcf_int16_vector_end; ++j) {
		    if (x[j]>>1 > 1) {
			ncounts[i]++;
			depths[i] += (unsigned int)x1[0];
			gqualitys[i] += (unsigned int)x2[0];
			break;
		    }
		}
	    }
	}
    }
    return NULL;
}


/*
    Clean up.
*/
void destroy(void)
{
    int i;
    printf("#SAMPLE\tCOUNTS\tAVG_DEPTH\tAVG_GQ\n");
    for (i=0; i<nsamples; ++i) {
	printf("%s\t%u\t%.2f\t%u\n", header->samples[i], ncounts[i], (float)depths[i]/ncounts[i], gqualitys[i]/ncounts[i]);
    }
    free(ncounts);
    free(depths);
}


