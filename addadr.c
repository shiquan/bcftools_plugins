/* plugins alleledepthratio.c -- calculate the allele depth ratio
   
 */

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <math.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>

static int32_t *ala = NULL;
static int mala = 0;
static bcf_hdr_t *in_hdr, *out_hdr;

const char *about(void)
{
    return "Add ADR tag, allele depth ratio.\n";
}

const char *usage(void)
{
    return
	"\n"
	"About: Add ADR tag in the VCF file.\n"
	"Usage: bcftools +addadr [General Options] -- [Plugin Options]\n"
	"Options:\n"
	"   run \"bcftools plugin\" for a list of common optionsn"
	"\n"
	"Example:\n"
	"   bcftools +addadr in.vcf\n"
	"\n";
}

int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    in_hdr = in;
    out_hdr = out;
    bcf_hdr_append(out_hdr, "##FORMAT=<ID=ADR,Number=G,Type=Float,Description=\"Allele Depth Ratio\">");
    return 0;
}

bcf1_t *process(bcf1_t *rec)
{
    int i, n, m;
    int depth;
    m = bcf_get_format_int32(in_hdr, rec, "DP", &ala, &mala);
    if (m != 1)  depth = 0;
    else depth = ala[0];
    n = bcf_get_format_int32(in_hdr, rec, "AD", &ala, &mala);
    if (n < 1) return rec;
    float *ratios = (float*)calloc(n, sizeof(float));
    if (depth == 0) {
	i = 0;
	while (i < n)
	{
	    depth += ala[i];
	}
    }
    for (i = 0; i < n; i++)
    {
	ratios[i] = (float)ala[i]/depth;
    }
    bcf_update_format_float(out_hdr, rec, "ADR", ratios, n);
    free(ratios);
    return rec;
}

void destroy(void)
{
    free(ala);
}
