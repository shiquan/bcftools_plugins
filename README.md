Plugins:

addadr -- add ADR (Allele Depth Ratio) tag in your VCF files

ADR = (float)AD/DP


To compile these plugins:

** download the lastest bcftools
** tar zxvf bcftools-xx.tar.gz
** cd bcftools-xx/plugins
** git clone https://github.com/shiquan/bcftools_plugins.git
** cd ..; make
** mv plugins/*.so $BCFTOOLS_PLUGINS/
