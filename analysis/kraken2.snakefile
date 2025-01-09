"""
Annotate contigs using megablast
"""
import itertools

if RUN["megablast"]:
    '''
        Annotate contigs using megablast
    '''

    # Use only active databases as specified in the config
    active_databases = config['megablast']['active_databases']
    sample_assembler_combinations = expand(OUTDIR + "annotation/megablast/{sample}_{assembler}_{db}_megablast.txt",
                                        sample=SAMPLES,
                                        db=active_databases,
                                        assembler=get_assembly_list())
    all_outputs.extend(sample_assembler_combinations)


    megablast_config = config["megablast"]

    rule megablast:
        input:
            fasta=get_final_fasta_input() #get_fasta_input()
        output:
            megablast=OUTDIR+ "annotation/megablast/{sample}_{assembler}_{db}_megablast.txt" #output file
        params:
            prefix="{sample}_{assembler}_{db}",
            db_path=lambda wildcards: config['megablast']['databases'][wildcards.db]['path'],
            outfmt=megablast_config['outfmt'],
            culling_limit=megablast_config['culling_limit'],
            evalue=megablast_config['evalue'],
            options=megablast_config['options'],

        message: "Analysing using megablast on {wildcards.db} database:\nSample id: {params.prefix}"
        threads: 6 #THREADS
        conda: ENVDIR+"blast.yaml"
        shell:
            """
            blastn -task megablast \
            -query {input.fasta}  \
            -db {params.db_path}   \
            -outfmt {params.outfmt}  \
            -culling_limit {params.culling_limit}  \
            -num_threads {threads}  \
            -evalue {params.evalue}  \
            -out {output.megablast}
            """
