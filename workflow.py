from gwf import Workflow
from templates import *
import glob

files1 = sorted(glob.glob('data/*_R1.fastq'))
files2 = sorted(glob.glob('data/*_R2.fastq'))
files = list(zip(files1, files2))
rg = 'data/refGen.fa'

gwf = Workflow()

gwf.target_from_template('bwaIndex',
                         bwa_index(fa=rg,
                                   amb=rg+'.amb',
                                   ann=rg+'.ann',
                                   bwt=rg+'.bwt',
                                   pac=rg+'.pac',
                                   sa=rg+'.sa'
                                   ))

gwf.target_from_template('picardDict',
                         picard_dict(fa=rg,
                                     dict_=rg[0:11]+'.dict'
                                     ))

gwf.target_from_template('samtoolsFaidx',
                         samtools_faidx(fa=rg,
                                        fai=rg+'.fai'))

for i in range(len(files)):
    gwf.target_from_template('bwaMapping_{}'.format(files[i][0][11:16]),
                             bwa_map(fa=rg,
                                     fq1=files[i][0],
                                     fq2=files[i][1],
                                     output='results/mapped_{}.bam'.format(
                                         files[i][0][11:16])
                                     ))

    gwf.target_from_template('picardRG_{}'.format(files[i][0][11:16]),
                             picard_rg(mapped='results/mapped_{}.bam'.format(
                                 files[i][0][11:16]),
                                       rgroup='results/rg_{}.bam'.format(
                                           files[i][0][11:16]),
                                       name=files[i][0][11:16]
                             ))

    gwf.target_from_template('picardMD_{}'.format(files[i][0][11:16]),
                             picard_md(rgroup='results/rg_{}.bam'.format(
                                 files[i][0][11:16]),
                                       dups='results/md_{}.bam'.format(
                                           files[i][0][11:16]),
                                       bai='results/md_{}.bam.bai'.format(
                                           files[i][0][11:16]),
                                       sbi='results/md_{}.bam.sbi'.format(
                                           files[i][0][11:16])
                                       ))

    gwf.target_from_template('gatkHaplotypeCaller_{}'.format(files[i][0][11:16]),
                             gatk_haplotypecaller(fa=rg,
                                                  dict_=rg[0:11]+'.dict',
                                                  fai=rg+'.fai',
                                                  bam='results/md_{}.bam'.format(
                                                      files[i][0][11:16]),
                                                  gvcf='results/{}.g.vcf'.format(
                                                      files[i][0][11:16]),
                                                  idx='results/{}.g.vcf.idx'.format(
                                                      files[i][0][11:16])
                                                  ))

gwf.target_from_template('gatkCombinegvcfs',
                         gatk_combinegvcfs(fa=rg,
                                           gvcf0='results/Gifu_.g.vcf',
                                           gvcf1='results/MG001.g.vcf',
                                           gvcf2='results/MG002.g.vcf',
                                           cohort='results/All_accessions.g.vcf'
                                           ))

gwf.target_from_template('gatkGenotypegvcfs',
                         gatk_genotypegvcfs(fa=rg,
                                            gvcf='results/All_accessions.g.vcf',
                                            vcf='results/All_accessions.vcf'
                                            ))
