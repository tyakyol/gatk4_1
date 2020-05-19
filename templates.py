def bwa_index(fa, amb, ann, bwt, pac, sa):
    '''
    Template for genome indexing.
    '''
    inputs = [fa]
    outputs = [amb, ann, bwt, pac, sa]
    options = {}
    spec = '''
bwa index {fa}
    '''.format(fa=fa)
    return inputs, outputs, options, spec


def picard_dict(fa, dict_):
    '''
    Template for creating sequence dictionary. HaplotypeCaller from GATK needs
    this file.
    '''
    inputs = [fa]
    outputs = [dict_]
    options = {}
    spec = '''
java -jar ./libexec/picard/picard.jar CreateSequenceDictionary \
R={fa} \
O={dict_}   
    '''.format(fa=fa, dict_=dict_)
    return inputs, outputs, options, spec


def samtools_faidx(fa, fai):
    '''
    HaplotypeCaller from GATK needs this file.
    '''
    inputs = [fa]
    outputs = [fai]
    options = {}
    spec = '''
samtools faidx {fa} 
    '''.format(fa=fa)
    return inputs, outputs, options, spec


def bwa_map(fa, fq1, fq2, output):
    '''
    Template for mapping short reads to the genome.
    '''
    inputs = [fa, fq1, fq2,
              '{}.amb'.format(fa),
              '{}.ann'.format(fa),
              '{}.bwt'.format(fa),
              '{}.pac'.format(fa),
              '{}.sa'.format(fa)]
    outputs = [output]
    options = {}
    spec = '''
bwa mem {fa} {fq1} {fq2} | samtools view -Sb > {output}
'''.format(fa=fa, fq1=fq1, fq2=fq2, output=output)
    return inputs, outputs, options, spec


def picard_rg(mapped, ref, rgroup, name):
    '''
    Template for adding read groups.
    '''
    inputs = [mapped, ref, name]
    outputs = [rgroup]
    options = {}
    spec = '''
java -jar ./libexec/picard/picard.jar AddOrReplaceReadGroups \
I={mapped} \
O={rgroup} \
RGLB=lib1 \
RGPL=illumina \
RGPU=unit1 \
RGSM={name} \
R={ref} \
CREATE_INDEX=true \
MAX_RECORDS_IN_RAM=4000000 \
USE_JDK_INFLATER=true \
USE_JDK_DEFLATER=true
    '''.format(mapped=mapped, ref=ref, rgroup=rgroup, name=name)
    return inputs, outputs, options, spec


def gatk_md(rgroup, dups, bai, sbi):
    '''
    Template for marking duplicates. This code not only marks duplicates but
    also removes them, I am not sure if this is correct.
    '''
    inputs = [rgroup]
    outputs = [dups, bai, sbi]
    options = {}
    spec = '''
gatk MarkDuplicatesSpark \
-I {rgroup} \
-O {dups} \
--remove-sequencing-duplicates    
    '''.format(rgroup=rgroup, dups=dups)
    return inputs, outputs, options, spec


def gatk_haplotypecaller(fa, fai, dict_, bam, gvcf, idx):
    '''
    Template for generating gvcf files.
    '''
    inputs = [fa, fai, dict_, bam]
    outputs = [gvcf, idx]
    options = {}
    spec = '''
gatk --java-options '-Xmx4g' HaplotypeCaller \
-R {fa} \
-I {bam} \
-O {gvcf} \
-ERC GVCF
'''.format(fa=fa, bam=bam, gvcf=gvcf)
    return inputs, outputs, options, spec


def gvcf_list(gvcf, tsv):
    '''
    Template for creating map file of g.vcf files.
    '''
    inputs = []
    for g in gvcf:
        inputs.append(g)
    outputs = [tsv]
    options = {}
    spec = '''
python3.7 sample_map.py {tsv}   
    '''.format(tsv=tsv)
    return inputs, outputs, options, spec


def gatk_genomicsdbimport(samples, db):
    '''
    Template for combining g.vcf files.
    '''
    inputs = [samples]
    outputs = [db]
    options = {}
    spec = '''
gatk --java-options '-Xmx4g -Xms4g' GenomicsDBImport \
--sample-name-map {samples} \
--genomicsdb-workspace-path {db} \
-L LjG1.1_chr1 -L LjG1.1_chr2 
    '''.format(samples=samples, db=db)
    return inputs, outputs, options, spec


def gatk_genotypegvcfs(fa, db, vcf):
    '''
    Template for joint genotyping.
    '''
    inputs = [fa, db]
    outputs = [vcf]
    options = {}
    spec = '''
gatk --java-options '-Xmx4g' GenotypeGVCFs \
-R {fa} \
-V gendb://{db} \
-O {vcf}
    '''.format(fa=fa, db=db, vcf=vcf)
    return inputs, outputs, options, spec
