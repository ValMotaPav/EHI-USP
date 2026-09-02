# Login + site do Earth Hologenome Iniciative
ssh -J valentina.pavelecini@marfim.lad.pucrs.br valentina.pavelecini@pantanal.lad.pucrs.br  
cd /labgenomaarea2/valentina.pavelecini  
conda activate EHI  

https://www.earthhologenome.org/bioinformatics/index.html


# 31/08/2026

## Realizando o primeiro passo da pipeline do EHI: filtro de qualidade

```
2.1 Quality-filtering
Raw sequencing data require an initial preprocessing to get rid off low-quality nucleotides and reads, as well as any remains of sequencing adaptors that can mess around in the downstream analyses. An efficient way to do so is to use the software fastp, which can perform all above-mentioned operations in a single go and directly on compressed files.
fastp \
    --in1 {input.r1i} --in2 {input.r2i} \
    --out1 {output.r1o} --out2 {output.r2o} \
    --trim_poly_g \
    --trim_poly_x \
    --low_complexity_filter \
    --n_base_limit 5 \
    --qualified_quality_phred 20 \
    --length_required 60 \
    --thread {threads} \
    --html {output.fastp_html} \
    --json {output.fastp_json} \
    --adapter_sequence {params.adapter1} \
    --adapter_sequence_r2 {params.adapter2}
```

Primeiro faremos um teste com apenas um arquivo de Raw_Data.  
Para os adaptadores, clonamos o workflow da EHI utilizando ```git clone https://github.com/earthhologenome/EHI_bioinformatics.git``` e procuramos pelos adaptadores padrão Illumina utilizados na pipeline com
```grep -R "adapter1\|adapter2" /labgenomaarea2/valentina.pavelecini/EHI/EHI_bioinformatics```, que nos revela que os adaptadores são ```AGATCGGAAGAGCACACGTCTGAACTCCAGTCA``` para R1 e ```AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT``` para R2:

```
/labgenomaarea2/valentina.pavelecini/EHI/EHI_bioinformatics/0_Code/configs/1_Preprocess_QC_config.yaml:adapter1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
/labgenomaarea2/valentina.pavelecini/EHI/EHI_bioinformatics/0_Code/configs/1_Preprocess_QC_config.yaml:adapter2: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
/labgenomaarea2/valentina.pavelecini/EHI/EHI_bioinformatics/0_Code/1_Preprocess_QC.snakefile:        adapter1 = expand("{adapter1}", adapter1=config['adapter1']),
/labgenomaarea2/valentina.pavelecini/EHI/EHI_bioinformatics/0_Code/1_Preprocess_QC.snakefile:        adapter2 = expand("{adapter2}", adapter2=config['adapter2'])
/labgenomaarea2/valentina.pavelecini/EHI/EHI_bioinformatics/0_Code/1_Preprocess_QC.snakefile:            --adapter_sequence {params.adapter1} \
```

E então rodamos o comando:
```
fastp \
    --in1 /labgenomaarea2/valentina.pavelecini/raw_data/TF-2587-PM-1-A_S5_L001_R1_001.fastq.gz \
    --in2 /labgenomaarea2/valentina.pavelecini/raw_data/TF-2587-PM-1-A_S5_L001_R2_001.fastq.gz \
    --out1 /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-PM-1-A_R1.fastq.gz \
    --out2 /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-PM-1-A_R2.fastq.gz \
    --trim_poly_g \
    --trim_poly_x \
    --low_complexity_filter \
    --n_base_limit 5 \
    --qualified_quality_phred 20 \
    --length_required 60 \
    --thread 8 \
    --html /labgenomaarea2/valentina.pavelecini/EHI/TF-2587-PM-1-A_fastp.html \
    --json /labgenomaarea2/valentina.pavelecini/EHI/TF-2587-PM-1-A_fastp.json \
    --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
```
Em que:
- ```fastp \``` executa o programa fastp, que irá remover adaptadores, filtrar reads de baixa qualidade, remover reads muito curtos, tratar sequências de baixa complexidade, fazer trimming Poly-G/Poly-X, e gerar relatórios.
- ```--in1 /labgenomaarea2/valentina.pavelecini/raw_data/TF-2587-PM-1-A_S5_L001_R1_001.fastq.gz \``` indica o arquivo de entrada do Read 1
- ```--in2 /labgenomaarea2/valentina.pavelecini/raw_data/TF-2587-PM-1-A_S5_L001_R2_001.fastq.gz \``` indica o arquivo de entrada do Read 2
- ```--out1 /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-PM-1-A_R1.fastq.gz \``` indica o arquivo de saída do Read 1, ou seja, o Read filtrado
- ```--out2 /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-PM-1-A_R2.fastq.gz \``` indica o arquivo de saída do Read 2, ou seja, o Read filtrado
- ```--trim_poly_g \``` realiza o trimming de Poly-G, que são sequências com muitos G consecutivos
- ```--trim_poly_x \``` realiza trimming de Poly-X, que são sequências com muitas bases iguais consecutivas
- ```--low_complexity_filter \``` é similar ao trimming de Poly-X, mas com sequências longas com pouca diversidade de bases, e portanto, pouca informação para algumas análises
- ```--n_base_limit 5 \``` estabelece o limite para filtragem bases desconhecidas, significa que reads com mais de 5 bases desconhecidas podem ser filtrados
- ```--qualified_quality_phred 20 \``` estabelece que 20 é o limite de qualidade da escala Phred das bases, ou seja, que a probabilidade de erro da base é menor (~1%) se o limite fosse 10 (10%)
- ```--length_required 60 \``` indica o comprimento mínimo de pares de bases permitido para um read após o trimming, caso contrário, ele é descartado
- ```--thread 8 \``` indica que 8 threads/CPUs poderão trabalhar paralelamente
- ```--html /labgenomaarea2/valentina.pavelecini/EHI/TF-2587-PM-1-A_fastp.html \``` indica o diretório em que será criado o relatório em formato HTML e seu nome, que poderá ser aberto em um navegador e fará análises como: qualidade antes/depois, quantidade de reads, bases, adaptadores, trimming, duplicação, distribuição de comprimento, etc.
- ```--json /labgenomaarea2/valentina.pavelecini/EHI/TF-2587-PM-1-A_fastp.json \``` indica o diretório em que será criado o relatório em formato JSON e seu nome, que poderá ser utilizado com o MultiQC para gerar um relatório geral para todas as sequências, em seguida.
- ```--adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \``` indica qual sequência adaptadora estamos utilizando para R1
- - ```--adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT``` indica qual sequência adaptadora estamos utilizando para R2


Êxito:


```
Reduce worker threads to 2 due to CPU cores limit
Read1 before filtering:
total reads: 6589707
total bases: 995045757
Q20 bases: 974772737(97.9626%)
Q30 bases: 945056606(94.9762%)
Q40 bases: 0(0%)

Read2 before filtering:
total reads: 6589707
total bases: 995045757
Q20 bases: 966744610(97.1558%)
Q30 bases: 924231698(92.8833%)
Q40 bases: 0(0%)

Read1 after filtering:
total reads: 6442594
total bases: 914808225
Q20 bases: 904229986(98.8437%)
Q30 bases: 879472457(96.1374%)
Q40 bases: 0(0%)

Read2 after filtering:
total reads: 6442594
total bases: 914758959
Q20 bases: 899433314(98.3246%)
Q30 bases: 864863106(94.5455%)
Q40 bases: 0(0%)

Filtering result:
reads passed filter: 12885188
reads failed due to low quality: 130486
reads failed due to too many N: 1874
reads failed due to too short: 153640
reads failed due to low complexity: 4338
reads failed due to adapter dimer: 3888
reads with adapter trimmed: 3653530
bases trimmed due to adapters: 115675898
reads with polyX in 3' end: 43496
bases trimmed in polyX tail: 783414

Duplication rate: 18.4939%

Insert size peak (evaluated by paired-end reads): 151

JSON report: /labgenomaarea2/valentina.pavelecini/EHI/TF-2587-PM-1-A_fastp.json
HTML report: /labgenomaarea2/valentina.pavelecini/EHI/TF-2587-PM-1-A_fastp.html

fastp --in1 /labgenomaarea2/valentina.pavelecini/raw_data/TF-2587-PM-1-A_S5_L001_R1_001.fastq.gz --in2 /labgenomaarea2/valentina.pavelecini/raw_data/TF-2587-PM-1-A_S5_L001_R2_001.fastq.gz --out1 /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-PM-1-A_R1.fastq.gz --out2 /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-PM-1-A_R2.fastq.gz --trim_poly_g --trim_poly_x --low_complexity_filter --n_base_limit 5 --qualified_quality_phred 20 --length_required 60 --thread 8 --html /labgenomaarea2/valentina.pavelecini/EHI/TF-2587-PM-1-A_fastp.html --json /labgenomaarea2/valentina.pavelecini/EHI/TF-2587-PM-1-A_fastp.json --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
fastp v1.3.6, time used: 203 seconds
```

Observou-se que as threads foram diminuídas para 2 devido ao limite da CPU, e por isso esse valor será utilizado daqui pra frente.

______________________________________________________

Dado que o teste foi um sucesso, prosseguiremos com as outras 5 amostras.  
Vamos iniciar criando um script .sh com ```nano EHI_filtrados``` para automatizar esta etapa:  

```
#!/bin/bash

set -e

RAW="/labgenomaarea2/valentina.pavelecini/raw_data"
OUT="/labgenomaarea2/valentina.pavelecini/EHI/filtrados"
REPORT="/labgenomaarea2/valentina.pavelecini/EHI/relat_fastp"

mkdir -p "$OUT"
mkdir -p "$REPORT"

ADAPTER1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPTER2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

THREADS=2

samples=(
    "TF-2587-PM-5-A_S7_L001"
    "TF-2587-RF-1-B_S1_L001"
    "TF-2587-RF-2-B_S2_L001"
    "TF-2587-RF-4-B_S3_L001"
    "TF-2587-RF-5-B_S4_L001"
)

for sample in "${samples[@]}"; do

    echo "=========================================="
    echo "Processando: $sample"
    echo "Início: $(date)"
    echo "=========================================="

    fastp \
        --in1 "$RAW/${sample}_R1_001.fastq.gz" \
        --in2 "$RAW/${sample}_R2_001.fastq.gz" \
        --out1 "$OUT/${sample}_R1.fastq.gz" \
        --out2 "$OUT/${sample}_R2.fastq.gz" \
        --trim_poly_g \
        --trim_poly_x \
        --low_complexity_filter \
        --n_base_limit 5 \
        --qualified_quality_phred 20 \
        --length_required 60 \
        --thread "$THREADS" \
        --html "$REPORT/${sample}_fastp.html" \
        --json "$REPORT/${sample}_fastp.json" \
        --adapter_sequence "$ADAPTER1" \
        --adapter_sequence_r2 "$ADAPTER2"

    echo "Finalizado: $sample"
    echo "Fim: $(date)"
    echo ""

done

echo "=========================================="
echo "TODAS AS 5 AMOSTRAS FORAM PROCESSADAS"
echo "Fim: $(date)"
echo "=========================================="
```

Em que:
- ```#!/bin/bash``` indica que os comandos no arquivos devem ser executados utilizando o interpretador Bash
- ```set -e``` indica que se algum comando apresentar erro, o script será interrompido
- ```RAW="/labgenomaarea2/valentina.pavelecini/raw_data"``` cria uma variável chamada RAW, que se refere ao caminho do diretório "/labgenomaarea2/valentina.pavelecini/raw_data"
- ```OUT="/labgenomaarea2/valentina.pavelecini/EHI/filtrados"``` cria uma variável chamada OUT, que se refere ao caminho do diretório "/labgenomaarea2/valentina.pavelecini/EHI/filtrados", o qual será usado para colocar os arquivos pós-filtragem
- ```REPORT="/labgenomaarea2/valentina.pavelecini/EHI/relat_fastp"``` cria uma variável chamada REPORT, que se refere ao caminho do diretório "/labgenomaarea2/valentina.pavelecini/EHI/filtrados", o qual será usado para colocar os relatórios
- ```mkdir -p "$OUT"``` cria um diretório do correspondente á variável OUT. -p indica que ele só será criado caso já não existir
- ```mkdir -p "$REPORT"``` cria um diretório do correspondente á variável REPORT. -p indica que ele só será criado caso já não existir
- ```ADAPTER1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"``` cria a variável ADAPTER1 que indica o adaptador de R1
- ```ADAPTER2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"``` cria a variável ADAPTER2 que indica o adaptador de R2
- ```THREADS=2``` limita as threads/CPU para 2, como feito automaticamente antes
-
```
samples=(
    "TF-2587-PM-5-A_S7_L001"
    "TF-2587-RF-1-B_S1_L001"
    "TF-2587-RF-2-B_S2_L001"
    "TF-2587-RF-4-B_S3_L001"
    "TF-2587-RF-5-B_S4_L001"
)
```
indica as 5 amostras restantes para analisar  

- ```
  for sample in "${samples[@]}"; do

    echo "=========================================="
    echo "Processando: $sample"
    echo "Início: $(date)"
    echo "=========================================="

    fastp \
        --in1 "$RAW/${sample}_R1_001.fastq.gz" \
        --in2 "$RAW/${sample}_R2_001.fastq.gz" \
        --out1 "$OUT/${sample}_R1.fastq.gz" \
        --out2 "$OUT/${sample}_R2.fastq.gz" \
        --trim_poly_g \
        --trim_poly_x \
        --low_complexity_filter \
        --n_base_limit 5 \
        --qualified_quality_phred 20 \
        --length_required 60 \
        --thread "$THREADS" \
        --html "$REPORT/${sample}_fastp.html" \
        --json "$REPORT/${sample}_fastp.json" \
        --adapter_sequence "$ADAPTER1" \
        --adapter_sequence_r2 "$ADAPTER2"

    echo "Finalizado: $sample"
    echo "Fim: $(date)"
    echo ""
```
done

echo "=========================================="
echo "TODAS AS 5 AMOSTRAS FORAM PROCESSADAS"
echo "Fim: $(date)"
echo "=========================================="
```
indica que para cada amostra em samples, fastp será executado seguindo o comando que realizamos no teste, substituindo cada parâmetro pela variável correspondente. Além disso, serão ecoados no terminal indicadores de início/fim do(s) processo(s), acompanhado da data e hora de início e fim.


#02/09/2026

Rodando o script automatizado com ```bash EHI_filtragem.sh```, que retornou:
```
==========================================
Processando: TF-2587-PM-5-A_S7_L001
Início: Wed Sep  2 08:50:08 -03 2026
==========================================
```
