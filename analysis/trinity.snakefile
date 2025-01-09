"""
Trinity assembly
"""

if RUN["trinity"]:
    '''
        Run trinity to assemble Illumina reads
    '''
    # Add final output files from this module to 'all_outputs' from
    # the main Snakefile scope. SAMPLES is also from the main Snakefile scope.
    trinity = expand(OUTDIR+"assembly/trinity/unprocessed/{sample}_trinity.fasta",
            sample=SAMPLES)

    all_outputs.extend(trinity)

    trinity_config = config["trinity"]

    rule trinity:
        input:
            forward_pair = short_read_input().forward_pair,
            reverse_pair = short_read_input().reverse_pair
        output:
            assembly=OUTDIR+"assembly/trinity/unprocessed/{sample}/Trinity.fasta"
        conda:
            ENVDIR+"trinity.yaml"
        params:
            folder=OUTDIR+"assembly/trinity/unprocessed/{sample}",
            tmpfolder=OUTDIR+"assembly/trinity/unprocessed/{sample}-trinity", # Trinity needs to include trinity in output folder
            sample="{sample}",
            memory=trinity_config["memory"],
            min_contig_length=trinity_config["min_contig_length"],
            options=trinity_config["options"],
        message: "Assembly using trinity:\nR1:\n{input.forward_pair}\nR2:\n{input.reverse_pair}\n"
        log: LOGDIR+"assembly/trinity/{sample}.log"
        benchmark:
            BENCHMARKDIR+"trinity/{sample}.bmk"
        threads: THREADS
        shell:
            """
                Trinity \
                --seqType fq \
                --left {input.forward_pair} \
                --right {input.reverse_pair} \
                --CPU {threads} \
                --min_contig_length  {params.min_contig_length} \
                --output {params.tmpfolder} \
                --max_memory {params.memory} \
                {params.options} \
                2>&1 | tee {log}

                mv -T {params.tmpfolder} {params.folder}
            """

    rule trinity_cp:
        input: assembly = rules.trinity.output.assembly
        output: assembly = OUTDIR+"assembly/trinity/unprocessed/{sample}_trinity.fasta"
        shell:
            """
                cp {input.assembly} {output.assembly}
            """
