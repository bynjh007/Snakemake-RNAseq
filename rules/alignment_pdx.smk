# PDX data

if config['options']['pdx']:

    # input fastq files
    def align_inputs(wildcards):
        fastq_path = path.join(raw_dir, "{sample}{{pair}}_trimmed.fastq.gz".format(
            sample=wildcards.sample))
        pairs = ["_R1", "_R2"] if config['options']['paired'] else [""]
        return expand(fastq_path, pair = pairs)

    # Using shared memory with genomeLoad
    rule star_preload:
        params:
            index_host = config['star_align']['index_host'],
            index_graft = config['star_align']['index_graft']
        output:
            temp(touch(path.join(bam_dir, 'star_preload.done')))
        conda:
            '../envs/star.yaml'
        log:
            path.join(log_dir, 'star_align', 'genome_preload.log')
        shell:
            'STAR --genomeLoad LoadAndExit --genomeDir {params.index_host}'
            '&& STAR --genomeLoad LoadAndExit --genomeDir {params.index_graft} 2> {log}'
    
    # Alignment to host genome
    rule star_align_host:
        input:
            fq = align_inputs,
            tmp = path.join(bam_dir, 'star_preload.done')
        output:
            temp(path.join(bam_dir, '{sample}_host','Aligned.out.bam'))
        params:
            index = config['star_align']['index_host'],
            options = format_options(config['star_align']['options']),
            out_prefix = path.join(bam_dir, '{sample}_host/'),  
            # '/' should be at the end of folder name: '{sample}/'
            qc_prefix = path.join(qc_dir,'star_align','{sample}_host'),
            qc_dir = path.join(qc_dir,'star_align_host')
        conda:
            '../envs/star.yaml'
        threads:
            config['star_align']['threads']
        log:
            path.join(log_dir, 'star_align', '{sample}_host.log')
        shell:
            'STAR {params.options} --genomeDir {params.index} '
            '--outFileNamePrefix {params.out_prefix} --runThreadN {threads} '
            '--readFilesCommand gunzip -c --outSAMtype BAM Unsorted '
            '--genomeLoad LoadAndKeep --limitBAMsortRAM 30000000000 '
            '--readFilesIn {input.fq} 2> {log}'
        
    # Alignment to graft genome
    rule star_align_graft:
        input:
            fq = align_inputs,
            tmp = path.join(bam_dir, 'star_preload.done')
        output:
            temp(path.join(bam_dir, '{sample}_graft','Aligned.out.bam'))
        params:
            index = config['star_align']['index_graft'],
            options = format_options(config['star_align']['options']),
            out_prefix = path.join(bam_dir, '{sample}_graft/'),  
            # '/' should be at the end of folder name: '{sample}/'
            qc_prefix = path.join(qc_dir,'star_align_graft','{sample}_graft'),
            qc_dir = path.join(qc_dir,'star_align_graft')
        conda:
            '../envs/star.yaml'
        threads:
            config['star_align']['threads']
        log:
            path.join(log_dir, 'star_align', '{sample}_graft.log')
        shell:
            'STAR {params.options} --genomeDir {params.index} '
            '--outFileNamePrefix {params.out_prefix} --runThreadN {threads} '
            '--readFilesCommand gunzip -c --outSAMtype BAM Unsorted '
            '--genomeLoad LoadAndKeep --limitBAMsortRAM 30000000000 '
            '--readFilesIn {input.fq} 2> {log} '
            '&& mkdir -p {params.qc_dir} '
            '&& mv {params.out_prefix}Log.final.out {params.qc_prefix}.Log.final.out'

    rule star_unload:
        input:
            expand(path.join(bam_dir, '{sample}_{organism}','Aligned.out.bam'),
                sample = all_samples, organism = ['host', 'graft'])
        output:
            touch(path.join(bam_dir, 'star_unload.done'))
        params:
            index_host = config['star_align']['index_host'],
            index_graft = config['star_align']['index_graft']
        conda:
            '../envs/star.yaml'
        log:
            path.join(log_dir, 'star_align', 'genome_remove.log')
        shell:
            'STAR --genomeLoad Remove --genomeDir {params.index_host} '
            '&& STAR --genomeLoad Remove --genomeDir {params.index_graft} '
            '2> {log}'

    rule sambamba_sort_qname:
        input:
            path.join(bam_dir, '{sample}_{organism}', 'Aligned.out.bam')
        output:
            temp(path.join(bam_dir, '{sample}_{organism}', 'Aligned.sorted.bam'))
        conda:
            '../envs/sambamba.yaml'
        threads:
            config['sambamba_sort']['threads']
        shell:
            'sambamba sort --natural-sort -t {threads} '
            '-o {output} {input}'

    rule disambiguate:
        input:
            host = path.join(bam_dir, '{sample}_host', 'Aligned.sorted.bam'),
            graft = path.join(bam_dir, '{sample}_graft', 'Aligned.sorted.bam')
        output:
            host_amb = temp(path.join(bam_dir, '{sample}_host', 'Aligned.amb.bam')),
            graft_amb = temp(path.join(bam_dir, '{sample}_graft', 'Aligned.amb.bam')),
            host_disamb = temp(path.join(bam_dir, '{sample}_host', 'Aligned.disamb.bam')),
            graft_disamb = temp(path.join(bam_dir, '{sample}_graft', 'Aligned.disamb.bam'))
        params:
            algorithm = config['disambiguate']['algorithm'],
            prefix = '',
            out_dir = path.join(bam_dir, '{sample}'),
            options = config['disambiguate']['options']
        conda:
            '../envs/disambiguate.yaml'
        shell:
            'ngs_disambiguate {params.options} -o {params.out_dir} -s '
            '{params.prefix} -a {params.algorithm} {input.host} {input.graft} '
            '&& mv {params.out_dir}/ambiguousSpeciesA.bam {output.host_amb} '
            '&& mv {params.out_dir}/ambiguousSpeciesB.bam {output.graft_amb} '
            '&& mv {params.out_dir}/disambiguatedSpeciesA.bam {output.host_disamb} '
            '&& mv {params.out_dir}/disambiguatedSpeciesB.bam {output.graft_disamb}'

    rule sambamba_sort:
        input:
            path.join(bam_dir, '{sample}_graft', 'Aligned.disamb.bam')
        output:
            path.join(bam_dir, '{sample}', 'Aligned.sortedByCoord.out.bam')
        conda: 
            '../envs/sambamba.yaml'
        params:
            options = config['sambamba_sort']['options']
        shell:
            'sambamba sort -t {threads} {params.options} -o {output} {input}'

    rule sambamba_index:
        input:
            path.join(bam_dir, '{sample}', 'Aligned.sortedByCoord.out.bam')
        output:
            path.join(bam_dir, '{sample}', 'Aligned.sortedByCoord.out.bam.bai')
        conda: 
            '../envs/sambamba.yaml'
        params:
            rm_dir = path.join(bam_dir, '{sample}_host')
        shell:
            'sambamba index -t {threads} {input} {output} '
            '&& rm -rf {params.rm_dir}_host'


