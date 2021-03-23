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
    threads: config['args']['processes']
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
    threads: config['args']['processes']
    output:
        config['targets']['fastqc']
    conda: config['envs']['preprocessing']
    script: config['scripts']['fastqc']


# noinspection PyUnresolvedReferences
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
    threads: config['args']['processes']
    output:
        config['targets']['ngs_te_mapper2'][0],
        config['targets']['ngs_te_mapper2'][1]
    conda: config['envs']['ngs_te_mapper2']
    script: config['scripts']['ngs_te_mapper2']

# noinspection PyUnresolvedReferences
rule rule_run_vep:
    input:
        ngs_te_mapper2_vcfs = config['targets']['ngs_te_mapper2'],
    params:
        log_file = config['logs']['vep'],
        mex_path = config['params']['mex_path'],
        assembly = config['args']['assembly']
    threads: config['args']['processes']
    output:
        config['targets']['vep'][0], config['targets']['vep'][1],
    conda: config['envs']['vep']
    script: config['scripts']['vep']