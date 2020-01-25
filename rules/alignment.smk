################################################################################ 
# STAR alignment
################################################################################ 
if config['options']['pdx'] == False:


    # input fastq files
    def align_inputs(wildcards):
        fastq_path = path.join(raw_dir, "{sample}{{pair}}_trimmed.fastq.gz".format(
            sample=wildcards.sample))
        pairs = ["_R1", "_R2"] if config['options']['paired'] else [""]
        return expand(fastq_path, pair = pairs)
    
    
    # Using shared memory with genomeLoad
    rule star_preload:
        params:
            index = config['star_align']['index']
        output:
            temp(touch(path.join(bam_dir, 'star_preload.done')))
        conda:
            '../envs/star.yaml'
        log:
            path.join(log_dir, 'star_align', 'genome_preload.log')
        shell:
            'STAR --genomeLoad LoadAndExit --genomeDir {params.index} 2> {log}'
    
    
    ##############################################
    # 2-PASS alignemt for novel junction discovery
    ##############################################
    
    if config['options']['star_2pass']:
    
    # to run 2-pass alignment from the own data, refer to :
    # https://groups.google.com/forum/#!topic/rna-star/Cpsf-_rLK9I
    
    # Or you can also use CTAT junction information and run usual 1-pass alignment
    # if the data is driven from cancer
    # https://github.com/NCIP/Trinity_CTAT/wiki
    
    # Because this is TCGA cancer data, junction information from TCGA tumors
    # is included in CTAT library
    
    # normal alignment process to generate SJ.out.tab
        rule star_align_1st:
            input:
                fq = align_inputs,
                tmp = path.join(bam_dir, 'star_preload.done')
            output:
                path.join(bam_dir, 'star_1', '{sample}', 'SJ.out.tab')
            params:
                index = config['star_align']['index'],
                options = format_options(config['star_align']['options']),
                out_prefix = path.join(bam_dir, 'star_1', '{sample}/')  # '/' should be at the end of folder name: '{sample}/'
            conda:
                '../envs/star.yaml'
            threads:
                config['star_align']['threads']
            log:
                path.join(log_dir, 'star_align', '{sample}_1.log')
            shell:
                'STAR {params.options} --genomeDir {params.index} '
                '--outFileNamePrefix {params.out_prefix} --runThreadN {threads} '
                '--readFilesCommand gunzip -c --outSAMtype BAM SortedByCoordinate '
                '--genomeLoad LoadAndKeep --limitBAMsortRAM 30000000000 '
                '--readFilesIn {input.fq} 2> {log}'
    
    # Genome unloading
        rule star_unload_1st:
            input:
                expand(path.join(bam_dir, 'star_1', '{sample}','Aligned.sortedByCoord.out.bam'), sample = all_samples)
            output:
                temp(touch(path.join(bam_dir, 'star_unload_1.done')))
            params:
                index=config['star_align']['index']
            conda:
                '../envs/star.yaml'
            log:
                path.join(log_dir, 'star_align', 'genome_remove.log')
            shell:
                'STAR --genomeLoad Remove --genomeDir {params.index} > {log}'
        
    # Filtering out the junctions with low reliability to reduce the size of SJ list
    # remove junctions with # of uniquely mapped reads crossing the junction <=2
    # also remove annotated junctions because they will be included later stage
        rule filter_SJ:
            input:
                sj = expand(path.join(bam_dir, 'star_1', '{sample}', 'SJ.out.tab'), sample = all_samples),
                tmp = path.join(bam_dir, 'star_unload_1.done')
            output:
                path.join(bam_dir, 'SJ_db', 'SJ.out.comb.tab')
            params:
                path.join(bam_dir, 'star_1')
            shell:
                "cat {input.sj} | awk '($7>2 && $6==0)' | cut -f1-6 | sort | uniq > {output} "
                "&& rm -rf {params}"
        
    # indexing genome with annotations and combined SJ list
        rule star_SJ_indexing:
            input:
                path.join(bam_dir, 'SJ_db', 'SJ.out.comb.tab')
            output:
                path.join(bam_dir, 'SJ_db', 'sjdbList.out.tab')
            params:
                sjdb = path.join(bam_dir, 'SJ_db'),
                fa = config['star_align']['fasta'],
                gtf = config['star_align']['gtf'],
                options = format_options(config['star_align']['options'])
            conda:
                '../envs/star.yaml'
            threads:
                config['star_align']['threads_indexing']
            log:
                path.join(log_dir, 'star_align', '{sample}_SJ_indexing.log')
            shell:
                'STAR --runMode genomeGenerate {params.options} --genomeDir {params.sjdb} '
                '--genomeFastaFiles {params.fa} --sjdbGTFfile {params.gtf} '
                '--runThreadN {threads} --sjdbFileChrStartEnd {input} 2>{log}'
    
    # Using shared memory with genomeLoad
        rule star_preload_2nd:
            input:
                path.join(bam_dir, 'SJ_db', 'sjdbList.out.tab')
            output:
                temp(touch(path.join(bam_dir, 'star_preload_2.done')))
            params:
                index = path.join(bam_dir, 'SJ_db')
            conda:
                '../envs/star.yaml'
            log:
                path.join(log_dir, 'star_align', 'genome_preload_2.log')
            shell:
                'STAR --genomeLoad LoadAndExit --genomeDir {params.index} 2> {log}'
    
    # second alignment    
        rule star_align_2nd:
            input:
                fq = align_inputs,
                tmp = path.join(bam_dir, 'star_preload_2.done'),
                sj = path.join(bam_dir, 'SJ_db', 'sjdbList.out.tab')
            output: 
                protected(path.join(bam_dir, '{sample}','Aligned.sortedByCoord.out.bam'))
            params:
                index = path_join(bam_dir, 'SJ_db'),
                options = format_options(config['star_align']['options']),
                out_prefix = path.join(bam_dir, '{sample}/'),  # '/' should be at the end of folder name: '{sample}/'
                qc_prefix = path.join(qc_dir, 'star_align','{sample}'),
                qc_dir = path.join(qc_dir, 'star_align')
            conda:
                '../envs/star.yaml'
            threads:
                config['star_align']['threads']
            log:
                path.join(log_dir, 'star_align', '{sample}_2.log')
            shell:
                'STAR {params.options} --genomeDir {params.index} '
                '--outFileNamePrefix {params.out_prefix} --runThreadN {threads} '
                '--readFilesCommand gunzip -c --outSAMtype BAM SortedByCoordinate '
                '--genomeLoad LoadAndKeep --limitBAMsortRAM 30000000000 '
                '--readFilesIn {input.fq} 2> {log} '
                '&& mkdir -p {params.qc_dir} '
                '&& mv {params.out_prefix}Log.final.out {params.qc_prefix}.Log.final.out'
    
    # Genome unloading
        rule star_unload:
            input:
                expand(path.join(bam_dir, '{sample}','Aligned.sortedByCoord.out.bam'),
                    sample=all_samples)
            output:
                touch(path.join(bam_dir, 'star_unload.done'))
            params:
                index = path.join(bam_dir, 'SJ_db') 
            conda:
                '../envs/star.yaml'
            log:
                path.join(log_dir, 'star_align', 'genome_remove.log')
            shell:
                'STAR --genomeLoad Remove --genomeDir {params.index} > {log}'
    
    
    ##############################################
    # 1-PASS alignment 
    ##############################################
    
    else:
    
        rule star_align:
            input:
                fq = align_inputs,
                tmp = path.join(bam_dir, 'star_preload.done')
            output:
                protected(path.join(bam_dir, '{sample}','Aligned.sortedByCoord.out.bam'))
            params:
                index = config['star_align']['index'],
                options = format_options(config['star_align']['options']),
                out_prefix = path.join(bam_dir, '{sample}/'),  # '/' should be at the end of folder name: '{sample}/'
                qc_prefix = path.join(qc_dir,'star_align','{sample}'),
                qc_dir = path.join(qc_dir,'star_align')
            conda:
                '../envs/star.yaml'
            threads:
                config['star_align']['threads']
            log:
                path.join(log_dir, 'star_align', '{sample}.log')
            shell:
                'STAR {params.options} --genomeDir {params.index} '
                '--outFileNamePrefix {params.out_prefix} --runThreadN {threads} '
                '--readFilesCommand gunzip -c --outSAMtype BAM SortedByCoordinate '
                '--genomeLoad LoadAndKeep --limitBAMsortRAM 30000000000 '
                '--readFilesIn {input.fq} 2> {log} '
                '&& mkdir -p {params.qc_dir} '
                '&& mv {params.out_prefix}Log.final.out {params.qc_prefix}.Log.final.out'
        
        # Genome unloading
        rule star_unload:
            input:
                expand(path.join(bam_dir, '{sample}','Aligned.sortedByCoord.out.bam'),
                    sample=all_samples)
            output:
                touch(path.join(bam_dir, 'star_unload.done'))
            params:
                index=config['star_align']['index']
            conda:
                '../envs/star.yaml'
            log:
                path.join(log_dir, 'star_align', 'genome_remove.log')
            shell:
                'STAR --genomeLoad Remove --genomeDir {params.index} > {log}'
        
    
    
    ################################################################################ 
    # indexing the bam files
    ################################################################################ 
    
    rule sambamba_index:
        input:
            path.join(bam_dir, '{sample}','Aligned.sortedByCoord.out.bam')
        output:
            path.join(bam_dir, '{sample}','Aligned.sortedByCoord.out.bam.bai')
        conda:
            '../envs/sambamba.yaml'
        threads:
            config['sambamba_index']['threads']
        shell:
            'sambamba index -t {threads} {input} {output}'
