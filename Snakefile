#
# Copyright (c) 2014 Hyeshik Chang
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
# - Hyeshik Chang <hyeshik@snu.ac.kr>
#

shell.prefix('set -e; set -o pipefail; ')
shell.executable(os.popen('which bash').read().strip()) # pipefail is supported by bash only.

SAMPLES = ['Control', 'TUT247KD']

TRIM_FIRST_CYCLES = 15 # to suppress mismatches from mismatched primer
MINIMUM_LENGTH = 15 + TRIM_FIRST_CYCLES
ADAPTER_SEQUENCE = 'TGGAATTCTCGGGTGCCAAGG'

SRR_ACCESSIONS = {
    'Control': 'SRR1734356',
    'TUT247KD': 'SRR1734357',
}

subworkflow reference_preparation:
    workdir: 'reference'
    snakefile: 'reference/Snakefile'

rule all:
    input:
        expand('sequences/{sample}.fa.gz', sample=SAMPLES),
        expand('alignments/{sample}.psl.gz', sample=SAMPLES),
        expand('stats/{sample}.mods.txt.gz', sample=SAMPLES),
        expand('geosubmission/{sample}.txt.gz', sample=SAMPLES)

rule download_fastq_from_sra:
    output: 'sequences/{sample}.fq.gz'
    params: srrno=lambda wc: SRR_ACCESSIONS[wc.sample]
    shell: 'cd sequences && fastq-dump --gzip {params.srrno} && \
            mv {params.srrno}.fastq.gz {wildcards.sample}.fq.gz'

rule trim_adapters:
    input: 'sequences/{sample}.fq.gz'
    output: 'sequences/{sample}.fa.gz'
    shell: 'zcat {input} | \
            cutadapt --trimmed-only --minimum-length={MINIMUM_LENGTH} \
                     --cut={TRIM_FIRST_CYCLES} -a {ADAPTER_SEQUENCE} - | \
            fastx_collapser -Q33 - | gzip -c - > {output}'

rule align_reads_by_blat:
    input: seq='sequences/{sample}.fa.gz', ref=reference_preparation('hairpins.fa')
    output: 'alignments/{sample}.psl.gz'
    shell: 'blat -noTrimA -tileSize=8 -stepSize=4 -noHead -out=pslx \
                -minIdentity=70 {input.ref} {input.seq} /dev/stdout | \
            grep "^[0-9]" | gzip -c - > {output}'

rule align_reads_by_blast:
    input: seq='sequences/{sample}.fa.gz', ref=reference_preparation('hairpins.fa')
    output: 'alignments/{sample}.blast.gz'
    params: refdb='reference/hairpins'
    threads: 32
    shell: 'zcat {input.seq} | \
            blastn -db {params.refdb} -outfmt 6 -num_threads {threads} -word_size 8 \
                -strand plus -query /dev/stdin | \
            gzip -c - > {output}'

rule find_untemplated_modifications:
    input: seq='sequences/{sample}.fa.gz', ref=reference_preparation('hairpins.fa'), \
           aln='alignments/{sample}.psl.gz'
    output: 'stats/{sample}.mods.txt.gz'
    shell: 'tools/psl2untemplatemods.py --read={input.seq} --untmpl-t-rescue \
                --reference={input.ref} --alignments={input.aln} --output - | \
            gzip -c - > {output}'

rule generate_GEO_sequence_tabular_list:
    input: 'stats/{sample}.mods.txt.gz'
    output: 'geosubmission/{sample}.txt.gz'
    shell: "zcat {input} | awk '{{printf \"%s\\t%d\\n\", $4, $3;}}' | \
            sed -e '1s/^.*$/SEQUENCE\tCOUNT/' | gzip -c - > {output}"

# ex: syntax=snakemake
