# To whom it may concern

This repo is a collection of bcftools plugins. By the time  I write all these plugins, I am working at BGI research as a bioinformatician. Thanks to HTSlib, bcftools has become much more flexible and robust than ever before. I am happy to write plugins to extent the features of bcftools and meanwhile I can manage my daily jobs by bcftools and these plugins. I would love to share all my codes with people who may need them. Please let me know if you have any problems when you try to use these plugins.



# Introduce to these plugins

Here is the brief introduction, try to enter each folder for details. Some folders may be empty, I will add the detailed manuals for them soon.

### addadr

This is a simple demo plugin to add the ADR tag in the vcf FORMAT field.

We denote `ADR = (float)AD/DP`

This plugin will be used as an example in my unofficial BCFtools guide (not release yet).



### select

This plugin is a variant of bcftools-query program, vcfquery can extract each field from VCF/BCF file line by line. But sometime we may want to split alleles and samples into different lines, this plugin can extract fields from VCF/BCF file more flexible.

Here is a demo of different format string of bcftools query and select:

``` 
bcftools query  -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' file.vcf.gz
bcftools +select -f '%BED\t%SAMPLE\t%REF\t%ALT\n' file.vcf.gz
```



# How to Compile:

There are three different ways to compiles these plugins.

### First way , compile all these plugins with latest bcftools. ***The basic way.***



1)  try to  download the lastest bcftools

2)  tar zxvf bcftools-xx.tar.gz; cd bcftools-xx/

3)  git clone https://github.com/shiquan/bcftools_plugins.git

4) mv bcftools_plugins/\*/\*.c  plugins/

5) make

6) mv plugins/*.so $BCFTOOLS_PLUGINS/



Note:  the environment variable BCFTOOLS_PLUGINS must be set properly , if you don’t know what I am talking about, please try to read the bcftools manual (https://samtools.github.io/bcftools/bcftools.html) or contact your system administer.



### Second way, compile into shared library directly, ***for high level user only***.

`gcc -g -Wall -Wc++-compat -I htslib_include_dir -fPIC -shared -o pluginX.so pluginX.c -L htslib_dir -lhts`

Note: htslib should be installed in your system and the version of bcftools and htslib should be consistent.



### Third way, compile the plugin into a standalone program, ***perhaps the simplest way***.

you should try to compile it with  `make` in each folder, or compiler it from scratch, here is an example.

`gcc -g -O2 -Ihtslib_include -D_SELECT_MAIN -o bcfselect plugins/select.c -Lhtslib_dir -lhts`