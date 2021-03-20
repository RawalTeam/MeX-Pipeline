import os

# noinspection PyUnresolvedReferences
rule rule_run_fastp:
    input:
        fq1 = config['args']['fastq1']
    params:
        fq2 = config['args']['fastq2'],
        paired = config['params']['paired'],
        output1_file = config['outputs']['fastp']['fqp1'],
        output2_file = config['outputs']['fastp']['fqp2'],
        log_file = config['logs']['fastp'],
        mex_path = config['params']['mex_path']
    threads: 2
    output:
        config['outputs']['fastp']['fqp1']
    conda: config['envs']['preprocessing']
    script: config['scripts']['fastp']

# noinspection PyUnresolvedReferences
rule rule_run_fastqc:
    input:
        fq1 = config['outputs']['fastp']['fqp1']
    params:
        fq2 = config['outputs']['fastp']['fqp2'],
        paired = config['params']['paired'],
        output_dir = config['outputs']['fastqc'],
        log_file = config['logs']['fastqc'],
        mex_path = config['params']['mex_path']
    threads: 2
    output:
        config['targets']['fastqc']
    conda: config['envs']['preprocessing']
    script: config['scripts']['fastqc']

# noinspection PyUnresolvedReferences
rule rule_run_ngs_te_mapper:
    input:
        fq1 = config['outputs']['fastp']['fqp1'],
        genome = config['args']['genome'],
        te_fasta = config['args']['te'],
    params:
        fq2 = config['outputs']['fastp']['fqp2'],
        paired = config['params']['paired'],
        tool_config = config['params']['tools_config'],
        output_dir = config['outputs']['ngs_te_mapper'],
        log_file = config['logs']['ngs_te_mapper'],
        mex_path = config['params']['mex_path']
    threads: config['args']['threads']
    output:
        config['targets']['ngs_te_mapper'][0],
        config['targets']['ngs_te_mapper'][1]
    conda: config['envs']['ngs_te_mapper']
    script: config['scripts']['ngs_te_mapper']

rule rule_run_ngs_te_mapper2:
    input:
        fq1 = config['outputs']['fastp']['fqp1'],
        genome = config['args']['genome'],
        te_fasta = config['args']['te'],
    params:
        fq2 = config['outputs']['fastp']['fqp2'],
        paired = config['params']['paired'],
        tool_config = config['params']['tools_config'],
        output_dir = config['outputs']['ngs_te_mapper2'],
        log_file = config['logs']['ngs_te_mapper2'],
        mex_path = config['params']['mex_path']
    threads: config['args']['threads']
    output:
        config['targets']['ngs_te_mapper2'][0]
    conda: config['envs']['ngs_te_mapper2']
    script: config['scripts']['ngs_te_mapper2']

# noinspection PyUnresolvedReferences
rule rule_run_vep:
    input:
        vcf = config['targets']['ngs_te_mapper'],
    params:
        log_file = config['logs']['vep'],
        mex_path = config['params']['mex_path'],
        assembly = config['args']['assembly']
    threads: config['args']['threads']
    output:
        config['targets']['vep']
    conda: config['envs']['vep']
    script: config['scripts']['vep']