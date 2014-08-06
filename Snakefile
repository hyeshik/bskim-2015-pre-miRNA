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

SAMPLES = """Control TUT247KD""".split()

PEAR_SUFFIXES = ['assembled', 'discarded', 'unassembled.forward', 'unassembled.reverse']

subworkflow reference_preparation:
    workdir: 'reference'
    snakefile: 'reference/Snakefile'

rule all:
    input:  expand('sequences/{sample}.assembled.fq.gz', sample=SAMPLES), \
            reference_preparation('hairpins/hairpins.version')

rule uncompress_fastq:
    input: 'sequences/{filename}.fq.gz'
    output: temp('scratch/{filename}.fq')
    shell: 'zcat {input} > {output}'

rule join_paired_end_reads:
    input: expand('scratch/{{sample}}_{read}.fq', read=['R1', 'R2'])
    output: expand(temp('scratch/{{sample}}.{suffix}.fastq'), suffix=PEAR_SUFFIXES)
    params: outputprefix='scratch/{sample}'
    threads: 32
    shell: 'tools/pear -j {threads} -y 8G -f {input[0]} -r {input[1]} \
                -o {params.outputprefix}'

rule compress_pear_outputs:
    input: 'scratch/{sample,[^.]+}.{suffix}.fastq'
    output: 'sequences/{sample}.{suffix}.fq.gz'
    threads: 32
    shell: 'pigz -p {threads} -c {input} > {output}'

