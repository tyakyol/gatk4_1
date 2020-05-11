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
java -jar ./share/picard-2.22.4-0/picard.jar CreateSequenceDictionary \
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


def picard_rg(mapped, rgroup, name):
    '''
    Template for adding read groups.
    '''
    inputs = [mapped, name]
    outputs = [rgroup]
    options = {}
    spec = '''
java -jar ./share/picard-2.22.4-0/picard.jar AddOrReplaceReadGroups \
I={mapped} \
O={rgroup} \
RGLB=lib1 \
RGPL=illumina \
RGPU=unit1 \
RGSM={name}
    '''.format(mapped=mapped, rgroup=rgroup, name=name)
    return inputs, outputs, options, spec


def picard_md(rgroup, dups, bai, sbi):
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


####### Base (Quality Score) Recalibration #######


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


def gatk_combinegvcfs(fa, gvcf0, gvcf1, gvcf2, cohort):
    '''
    Template for combining multiple gvcf files.
    GenomicsDBImport might be a better solution.
    '''
    inputs = [fa, gvcf0, gvcf1, gvcf2]
    outputs = [cohort]
    options = {}
    spec = '''
gatk CombineGVCFs \
-R {fa} \
--variant {gvcf0} \
--variant {gvcf1} \
--variant {gvcf2} \
-O {cohort}
    '''.format(fa=fa, gvcf0=gvcf0, gvcf1=gvcf1, gvcf2=gvcf2, cohort=cohort)
    return inputs, outputs, options, spec


def gatk_genotypegvcfs(fa, gvcf, vcf):
    '''
    Template for joint genotyping.
    '''
    inputs = [fa, gvcf]
    outputs = [vcf]
    options = {}
    spec = '''
gatk --java-options '-Xmx4g' GenotypeGVCFs \
-R {fa} \
-V {gvcf} \
-O {vcf}  
    '''.format(fa=fa, gvcf=gvcf, vcf=vcf)
    return inputs, outputs, options, spec
