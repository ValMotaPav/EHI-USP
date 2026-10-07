# Login + site do Earth Hologenome Iniciative + Fluxograma
```
ssh -J valentina.pavelecini@marfim.lad.pucrs.br valentina.pavelecini@pantanal.lad.pucrs.br
```
```
cd /labgenomaarea2/valentina.pavelecini  
conda activate EHI  
```

Pipeline: https://www.earthhologenome.org/bioinformatics/index.html  

Fluxograma: https://www.canva.com/design/DAHUDcTIaEc/Z9jiEwyVO9jZrfU_UPkavA/edit?ui=e30  

# Objetivo
Reconstruir um genoma procariótico a partir dos reads metagenômicos de amostra de água, produzindo um MAG (Metagenome-Assembled Genome)

# Realizando o primeiro passo da pipeline do EHI: filtro de qualidade

## 31/08/2026

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


## 02/09/2026

Rodando o script automatizado com ```bash EHI_filtragem.sh```, que retornou:
```
==========================================
Processando: TF-2587-PM-5-A_S7_L001
Início: Wed Sep  2 08:50:08 -03 2026
==========================================
Read1 before filtering:
total reads: 9985759
total bases: 1507849609
Q20 bases: 1480418421(98.1808%)
Q30 bases: 1436979655(95.2999%)
Q40 bases: 0(0%)

Read2 before filtering:
total reads: 9985759
total bases: 1507849609
Q20 bases: 1464950014(97.1549%)
Q30 bases: 1400006136(92.8479%)
Q40 bases: 0(0%)

Read1 after filtering:
total reads: 9798845
total bases: 1413505653
Q20 bases: 1397961173(98.9003%)
Q30 bases: 1360623963(96.2588%)
Q40 bases: 0(0%)

Read2 after filtering:
total reads: 9798845
total bases: 1413424118
Q20 bases: 1388147223(98.2117%)
Q30 bases: 1331828682(94.2271%)
Q40 bases: 0(0%)

Filtering result:
reads passed filter: 19597690
reads failed due to low quality: 183926
reads failed due to too many N: 2848
reads failed due to too short: 180998
reads failed due to low complexity: 4658
reads failed due to adapter dimer: 1398
reads with adapter trimmed: 4164630
bases trimmed due to adapters: 141649715
reads with polyX in 3' end: 60749
bases trimmed in polyX tail: 1034931

Duplication rate: 42.0141%

Insert size peak (evaluated by paired-end reads): 152

JSON report: /labgenomaarea2/valentina.pavelecini/EHI/fastp_reports/TF-2587-PM-5-A_S7_L001_fastp.json
HTML report: /labgenomaarea2/valentina.pavelecini/EHI/fastp_reports/TF-2587-PM-5-A_S7_L001_fastp.html

fastp --in1 /labgenomaarea2/valentina.pavelecini/raw_data/TF-2587-PM-5-A_S7_L001_R1_001.fastq.gz --in2 /labgenomaarea2/valentina.pavelecini/raw_data/TF-2587-PM-5-A_S7_L001_R2_001.fastq.gz --out1 /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-PM-5-A_S7_L001_R1.fastq.gz --out2 /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-PM-5-A_S7_L001_R2.fastq.gz --trim_poly_g --trim_poly_x --low_complexity_filter --n_base_limit 5 --qualified_quality_phred 20 --length_required 60 --thread 2 --html /labgenomaarea2/valentina.pavelecini/EHI/fastp_reports/TF-2587-PM-5-A_S7_L001_fastp.html --json /labgenomaarea2/valentina.pavelecini/EHI/fastp_reports/TF-2587-PM-5-A_S7_L001_fastp.json --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
fastp v1.3.6, time used: 317 seconds
Finalizado: TF-2587-PM-5-A_S7_L001
Fim: Wed Sep  2 08:55:26 -03 2026

==========================================
Processando: TF-2587-RF-1-B_S1_L001
Início: Wed Sep  2 08:55:26 -03 2026
==========================================
Read1 before filtering:
total reads: 27153134
total bases: 4100123234
Q20 bases: 4032254733(98.3447%)
Q30 bases: 3914129614(95.4637%)
Q40 bases: 0(0%)

Read2 before filtering:
total reads: 27153134
total bases: 4100123234
Q20 bases: 4000837905(97.5785%)
Q30 bases: 3838615384(93.622%)
Q40 bases: 0(0%)

Read1 after filtering:
total reads: 26644982
total bases: 3820107019
Q20 bases: 3779236655(98.9301%)
Q30 bases: 3680215874(96.338%)
Q40 bases: 0(0%)

Read2 after filtering:
total reads: 26644982
total bases: 3820003784
Q20 bases: 3760197317(98.4344%)
Q30 bases: 3621180897(94.7952%)
Q40 bases: 0(0%)

Filtering result:
reads passed filter: 53289964
reads failed due to low quality: 502976
reads failed due to too many N: 7526
reads failed due to too short: 464238
reads failed due to low complexity: 25710
reads failed due to adapter dimer: 15854
reads with adapter trimmed: 12754373
bases trimmed due to adapters: 433310480
reads with polyX in 3' end: 172548
bases trimmed in polyX tail: 3464229

Duplication rate: 39.1692%

Insert size peak (evaluated by paired-end reads): 151

JSON report: /labgenomaarea2/valentina.pavelecini/EHI/fastp_reports/TF-2587-RF-1-B_S1_L001_fastp.json
HTML report: /labgenomaarea2/valentina.pavelecini/EHI/fastp_reports/TF-2587-RF-1-B_S1_L001_fastp.html

fastp --in1 /labgenomaarea2/valentina.pavelecini/raw_data/TF-2587-RF-1-B_S1_L001_R1_001.fastq.gz --in2 /labgenomaarea2/valentina.pavelecini/raw_data/TF-2587-RF-1-B_S1_L001_R2_001.fastq.gz --out1 /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-RF-1-B_S1_L001_R1.fastq.gz --out2 /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-RF-1-B_S1_L001_R2.fastq.gz --trim_poly_g --trim_poly_x --low_complexity_filter --n_base_limit 5 --qualified_quality_phred 20 --length_required 60 --thread 2 --html /labgenomaarea2/valentina.pavelecini/EHI/fastp_reports/TF-2587-RF-1-B_S1_L001_fastp.html --json /labgenomaarea2/valentina.pavelecini/EHI/fastp_reports/TF-2587-RF-1-B_S1_L001_fastp.json --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
fastp v1.3.6, time used: 844 seconds
Finalizado: TF-2587-RF-1-B_S1_L001
Fim: Wed Sep  2 09:09:30 -03 2026

==========================================
Processando: TF-2587-RF-2-B_S2_L001
Início: Wed Sep  2 09:09:30 -03 2026
==========================================
Read1 before filtering:
total reads: 42723059
total bases: 6451181909
Q20 bases: 6335105148(98.2007%)
Q30 bases: 6137850356(95.143%)
Q40 bases: 0(0%)

Read2 before filtering:
total reads: 42723059
total bases: 6451181909
Q20 bases: 6258058716(97.0064%)
Q30 bases: 5978075705(92.6664%)
Q40 bases: 0(0%)

Read1 after filtering:
total reads: 41834992
total bases: 5937359191
Q20 bases: 5870570308(98.8751%)
Q30 bases: 5709384285(96.1603%)
Q40 bases: 0(0%)

Read2 after filtering:
total reads: 41834992
total bases: 5937233856
Q20 bases: 5834018888(98.2616%)
Q30 bases: 5600183995(94.3231%)
Q40 bases: 0(0%)

Filtering result:
reads passed filter: 83669984
reads failed due to low quality: 846650
reads failed due to too many N: 11762
reads failed due to too short: 865214
reads failed due to low complexity: 41700
reads failed due to adapter dimer: 10808
reads with adapter trimmed: 22628308
bases trimmed due to adapters: 801743789
reads with polyX in 3' end: 333342
bases trimmed in polyX tail: 6086987

Duplication rate: 48.5783%

Insert size peak (evaluated by paired-end reads): 151

JSON report: /labgenomaarea2/valentina.pavelecini/EHI/fastp_reports/TF-2587-RF-2-B_S2_L001_fastp.json
HTML report: /labgenomaarea2/valentina.pavelecini/EHI/fastp_reports/TF-2587-RF-2-B_S2_L001_fastp.html

fastp --in1 /labgenomaarea2/valentina.pavelecini/raw_data/TF-2587-RF-2-B_S2_L001_R1_001.fastq.gz --in2 /labgenomaarea2/valentina.pavelecini/raw_data/TF-2587-RF-2-B_S2_L001_R2_001.fastq.gz --out1 /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-RF-2-B_S2_L001_R1.fastq.gz --out2 /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-RF-2-B_S2_L001_R2.fastq.gz --trim_poly_g --trim_poly_x --low_complexity_filter --n_base_limit 5 --qualified_quality_phred 20 --length_required 60 --thread 2 --html /labgenomaarea2/valentina.pavelecini/EHI/fastp_reports/TF-2587-RF-2-B_S2_L001_fastp.html --json /labgenomaarea2/valentina.pavelecini/EHI/fastp_reports/TF-2587-RF-2-B_S2_L001_fastp.json --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
fastp v1.3.6, time used: 1328 seconds
Finalizado: TF-2587-RF-2-B_S2_L001
Fim: Wed Sep  2 09:31:38 -03 2026

==========================================
Processando: TF-2587-RF-4-B_S3_L001
Início: Wed Sep  2 09:31:38 -03 2026
==========================================
Read1 before filtering:
total reads: 10318813
total bases: 1558140763
Q20 bases: 1524124785(97.8169%)
Q30 bases: 1475318019(94.6845%)
Q40 bases: 0(0%)

Read2 before filtering:
total reads: 10318813
total bases: 1558140763
Q20 bases: 1506594561(96.6918%)
Q30 bases: 1436789946(92.2118%)
Q40 bases: 0(0%)

Read1 after filtering:
total reads: 9897059
total bases: 1370410315
Q20 bases: 1355701841(98.9267%)
Q30 bases: 1320333450(96.3458%)
Q40 bases: 0(0%)

Read2 after filtering:
total reads: 9897059
total bases: 1370391929
Q20 bases: 1348468879(98.4002%)
Q30 bases: 1297633798(94.6907%)
Q40 bases: 0(0%)

Filtering result:
reads passed filter: 19794118
reads failed due to low quality: 200108
reads failed due to too many N: 2738
reads failed due to too short: 625134
reads failed due to low complexity: 9930
reads failed due to adapter dimer: 5598
reads with adapter trimmed: 6910285
bases trimmed due to adapters: 258719439
reads with polyX in 3' end: 72472
bases trimmed in polyX tail: 1395069

Duplication rate: 46.7192%

Insert size peak (evaluated by paired-end reads): 40

JSON report: /labgenomaarea2/valentina.pavelecini/EHI/fastp_reports/TF-2587-RF-4-B_S3_L001_fastp.json
HTML report: /labgenomaarea2/valentina.pavelecini/EHI/fastp_reports/TF-2587-RF-4-B_S3_L001_fastp.html

fastp --in1 /labgenomaarea2/valentina.pavelecini/raw_data/TF-2587-RF-4-B_S3_L001_R1_001.fastq.gz --in2 /labgenomaarea2/valentina.pavelecini/raw_data/TF-2587-RF-4-B_S3_L001_R2_001.fastq.gz --out1 /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-RF-4-B_S3_L001_R1.fastq.gz --out2 /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-RF-4-B_S3_L001_R2.fastq.gz --trim_poly_g --trim_poly_x --low_complexity_filter --n_base_limit 5 --qualified_quality_phred 20 --length_required 60 --thread 2 --html /labgenomaarea2/valentina.pavelecini/EHI/fastp_reports/TF-2587-RF-4-B_S3_L001_fastp.html --json /labgenomaarea2/valentina.pavelecini/EHI/fastp_reports/TF-2587-RF-4-B_S3_L001_fastp.json --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
fastp v1.3.6, time used: 305 seconds
Finalizado: TF-2587-RF-4-B_S3_L001
Fim: Wed Sep  2 09:36:43 -03 2026

==========================================
Processando: TF-2587-RF-5-B_S4_L001
Início: Wed Sep  2 09:36:43 -03 2026
==========================================
Read1 before filtering:
total reads: 25842990
total bases: 3902291490
Q20 bases: 3823251041(97.9745%)
Q30 bases: 3702119866(94.8704%)
Q40 bases: 0(0%)

Read2 before filtering:
total reads: 25842990
total bases: 3902291490
Q20 bases: 3774302904(96.7202%)
Q30 bases: 3609168137(92.4884%)
Q40 bases: 0(0%)

Read1 after filtering:
total reads: 25206565
total bases: 3599243529
Q20 bases: 3555169366(98.7755%)
Q30 bases: 3456327177(96.0293%)
Q40 bases: 0(0%)

Read2 after filtering:
total reads: 25206565
total bases: 3599027607
Q20 bases: 3535995165(98.2486%)
Q30 bases: 3398878324(94.4388%)
Q40 bases: 0(0%)

Filtering result:
reads passed filter: 50413130
reads failed due to low quality: 736940
reads failed due to too many N: 6974
reads failed due to too short: 499956
reads failed due to low complexity: 23826
reads failed due to adapter dimer: 5154
reads with adapter trimmed: 12814364
bases trimmed due to adapters: 428453030
reads with polyX in 3' end: 225747
bases trimmed in polyX tail: 4285570

Duplication rate: 26.1427%

Insert size peak (evaluated by paired-end reads): 151

JSON report: /labgenomaarea2/valentina.pavelecini/EHI/fastp_reports/TF-2587-RF-5-B_S4_L001_fastp.json
HTML report: /labgenomaarea2/valentina.pavelecini/EHI/fastp_reports/TF-2587-RF-5-B_S4_L001_fastp.html

fastp --in1 /labgenomaarea2/valentina.pavelecini/raw_data/TF-2587-RF-5-B_S4_L001_R1_001.fastq.gz --in2 /labgenomaarea2/valentina.pavelecini/raw_data/TF-2587-RF-5-B_S4_L001_R2_001.fastq.gz --out1 /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-RF-5-B_S4_L001_R1.fastq.gz --out2 /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-RF-5-B_S4_L001_R2.fastq.gz --trim_poly_g --trim_poly_x --low_complexity_filter --n_base_limit 5 --qualified_quality_phred 20 --length_required 60 --thread 2 --html /labgenomaarea2/valentina.pavelecini/EHI/fastp_reports/TF-2587-RF-5-B_S4_L001_fastp.html --json /labgenomaarea2/valentina.pavelecini/EHI/fastp_reports/TF-2587-RF-5-B_S4_L001_fastp.json --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
fastp v1.3.6, time used: 807 seconds
Finalizado: TF-2587-RF-5-B_S4_L001
Fim: Wed Sep  2 09:50:10 -03 2026

==========================================
TODAS AS 5 AMOSTRAS FORAM PROCESSADAS
Fim: Wed Sep  2 09:50:10 -03 2026
==========================================
```
Êxito.

# Próxima etapa: Avaliação da complexidade taxonômica com Nonpareil

```
nonpareil \
    -s {input.non_host_r1} \
    -f fastq \
    -T kmer \
    -t {threads} \
    -b {wildcards.sample}

#Script to extract nonpareil values of interest
Rscript {config[codedir]}/scripts/nonpareil_table.R {output.npo} {output.npstats}
```

OBSERVAÇÃO: Pularemos completamente a etapa da pipeline do EHI de separar dados do hospedeiro, afinal estamos trabalhando com amostras de água, sem um hospedeiro conhecido específico.  

Começaremos instalando o nonparail com ```conda install -c bioconda nonpareil```. Depois, criamos um diretório chamado ```nonpareil``` com ```mkdir nonpareil```, onde irão os arquivos de saída.  

Iremos fazer o primeiro teste com o arquivo ```TF-2587-PM-1-A_S5_L001_R1.fastq.gz```, o mesmo que utilizamos para teste na primeira etapa, porém filtrado. Depois, automatizaremos.

```
nonpareil \
    -s /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-PM-1-A_S5_L001_R1.fastq.gz \
    -f fastq \
    -T kmer \
    -t 2 \
    -b TF-2587-PM-1-A_S5_L001
```
Obs: este código está sendo executado em ```/labgenomaarea2/valentina.pavelecini/EHI/nonpareil```, pois é onde quero que fiquem os arquivo de saída.  
Retornou erro ```Fatal error: Segmentation fault (core dumped)``` porque coloquei o nome errado, é ```TF-2587-PM-1-A_R1.fastq.gz```

```
nonpareil \
    -s /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-PM-1-A_R1.fastq.gz \
    -f fastq \
    -T kmer \
    -t 2 \
    -b TF-2587-PM-1-A_S5_L001
```
Êxito:
```
Nonpareil v3.5.5
 [      0.4]   The file /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-PM-1-A_R1.fastq.gz.enve-tmp.500452 was created
 [      0.4]  Reading /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-PM-1-A_R1.fastq.gz.enve-tmp.500452
 [      0.4]   Picking 10000 random sequences
 [      0.4]   Counting kmers
 [      2.2]  Read file with 6442594 sequences
 [      2.2]  Average read length is 141.993772bp
 [      2.2]  Sub-sampling library
 [      2.4]  Evaluating consistency
 [      2.4]  Everything seems correct
```

_______________________________________________________

Agora vamos automatizar para as outras 5 amostras. Começamos criando o arquivo .sh com ```nano nonpareil.sh``` e criando o comando:
```
#!/bin/bash
set -e

FILTRADOS="/labgenomaarea2/valentina.pavelecini/EHI/filtrados"
OUT="/labgenomaarea2/valentina.pavelecini/EHI/nonpareil"

mkdir -p "$OUT"

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

    nonpareil \
        -s "$FILTRADOS/${sample}_R1.fastq.gz" \
        -f fastq \
        -T kmer \
        -t "$THREADS" \
        -b "$OUT/$sample"

    echo "Finalizado: $sample"
    echo "Fim: $(date)"
    echo ""

done

echo "=========================================="
echo "TODAS AS 5 AMOSTRAS FORAM PROCESSADAS"
echo "Fim: $(date)"
echo "=========================================="
```

E iniciando com ```bash nonpareil.sh```:

```
==========================================
Processando: TF-2587-PM-5-A_S7_L001
Início: Wed Sep  2 10:42:02 -03 2026
==========================================
Nonpareil v3.5.5
 [      0.6]   The file /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-PM-5-A_S7_L001_R1.fastq.gz.enve-tmp.500778 was created
 [      0.6]  Reading /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-PM-5-A_S7_L001_R1.fastq.gz.enve-tmp.500778
 [      0.6]   Picking 10000 random sequences
 [      0.6]   Counting kmers
 [      3.4]  Read file with 9798845 sequences
 [      3.4]  Average read length is 144.252272bp
 [      3.4]  Sub-sampling library
 [      3.7]  Evaluating consistency
 [      3.7]  Everything seems correct
Finalizado: TF-2587-PM-5-A_S7_L001
Fim: Wed Sep  2 10:46:05 -03 2026

==========================================
Processando: TF-2587-RF-1-B_S1_L001
Início: Wed Sep  2 10:46:05 -03 2026
==========================================
Nonpareil v3.5.5
 [      1.7]   The file /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-RF-1-B_S1_L001_R1.fastq.gz.enve-tmp.501007 was created
 [      1.7]  Reading /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-RF-1-B_S1_L001_R1.fastq.gz.enve-tmp.501007
 [      1.7]   Picking 10000 random sequences
 [      1.7]   Counting kmers
 [      9.1]  Read file with 26644982 sequences
 [      9.1]  Average read length is 143.370599bp
 [      9.1]  Sub-sampling library
 [      9.4]  Evaluating consistency
 [      9.4]  Everything seems correct
Finalizado: TF-2587-RF-1-B_S1_L001
Fim: Wed Sep  2 10:56:29 -03 2026

==========================================
Processando: TF-2587-RF-2-B_S2_L001
Início: Wed Sep  2 10:56:29 -03 2026
==========================================
Nonpareil v3.5.5
 [      2.5]   The file /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-RF-2-B_S2_L001_R1.fastq.gz.enve-tmp.501362 was created
 [      2.5]  Reading /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-RF-2-B_S2_L001_R1.fastq.gz.enve-tmp.501362
 [      2.5]   Picking 10000 random sequences
 [      2.5]   Counting kmers
 [     14.5]  Read file with 41834992 sequences
 [     14.5]  Average read length is 141.923278bp
 [     14.5]  Sub-sampling library
 [     14.7]  Evaluating consistency
 [     14.7]  Everything seems correct
Finalizado: TF-2587-RF-2-B_S2_L001
Fim: Wed Sep  2 11:12:55 -03 2026

==========================================
Processando: TF-2587-RF-4-B_S3_L001
Início: Wed Sep  2 11:12:55 -03 2026
==========================================
Nonpareil v3.5.5
 [      0.6]   The file /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-RF-4-B_S3_L001_R1.fastq.gz.enve-tmp.502151 was created
 [      0.6]  Reading /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-RF-4-B_S3_L001_R1.fastq.gz.enve-tmp.502151
 [      0.6]   Picking 10000 random sequences
 [      0.6]   Counting kmers
 [      3.3]  Read file with 9897059 sequences
 [      3.3]  Average read length is 138.466419bp
 [      3.3]  Sub-sampling library
 [      3.5]  Evaluating consistency
 [      3.5]  Everything seems correct
Finalizado: TF-2587-RF-4-B_S3_L001
Fim: Wed Sep  2 11:16:48 -03 2026

==========================================
Processando: TF-2587-RF-5-B_S4_L001
Início: Wed Sep  2 11:16:48 -03 2026
==========================================
Nonpareil v3.5.5
 [      1.5]   The file /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-RF-5-B_S4_L001_R1.fastq.gz.enve-tmp.502341 was created
 [      1.5]  Reading /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-RF-5-B_S4_L001_R1.fastq.gz.enve-tmp.502341
 [      1.5]   Picking 10000 random sequences
 [      1.5]   Counting kmers
 [      8.6]  Read file with 25206565 sequences
 [      8.6]  Average read length is 142.789925bp
 [      8.6]  Sub-sampling library
 [      8.8]  Evaluating consistency
 [      8.8]  Everything seems correct
Finalizado: TF-2587-RF-5-B_S4_L001
Fim: Wed Sep  2 11:26:38 -03 2026

==========================================
TODAS AS 5 AMOSTRAS FORAM PROCESSADAS
Fim: Wed Sep  2 11:26:38 -03 2026
==========================================
```

## 11/09/2026

A documentação disponível no site da Earth Hologenome Initiative apresenta uma etapa adicional após a execução do Nonpareil, na qual o script nonpareil_table.R é utilizado para extrair valores de interesse dos arquivos .npo e gerar arquivos .npstats:
```
#Script to extract nonpareil values of interest
Rscript {config[codedir]}/scripts/nonpareil_table.R {output.npo} {output.npstats}
```
Entretanto, como indicado após clonarmos o repositório mais recente do EHI, essa etapa não está presente na versão atual do 1_Preprocess_QC.snakefile que estamos utilizando.   
Na versão atual do workflow, os arquivos .npo gerados pelo Nonpareil são utilizados diretamente pela regra de geração do relatório, que cria o arquivo nonpareil_metadata.tsv a partir desses resultados. Portanto, não executaremos a etapa nonpareil_table.R descrita na página do site, pois ela não faz parte da versão do workflow que estamos reproduzindo.

Desta forma, a etapa de Avaliação da complexidade taxonômica com Nonpareil está completa.


# Próxima etapa: Avaliação da fração procariótica

Esta etapa avaliará a proporção, diversidade, e atividade de bactérias e archeas através de diferentes amostras ambientais. Essa etapa incluir duas análises diferentes:
```
#Run singlem pipe
singlem pipe \
    -1 {input.non_host_r1} \
    -2 {input.non_host_r2} \
    --otu-table {params.pipe_uncompressed} \
    --taxonomic-profile {output.condense} \
    --threads {threads}
```
Identifica/classifica marcadores taxonômicos e produz um perfil taxonômico. Ou seja, a abundância relativa dos microrganismos presentes na amostra ambiental ("Quem está lá?" e "Em qual quantidade?")

```
#Run singlem read_fraction
singlem read_fraction \
    -1 {input.non_host_r1} \
    -2 {input.non_host_r2} \
    --input-profile {output.condense} \
    --output-tsv {output.read_fraction} \
    --output-per-taxon-read-fractions {params.read_fraction_taxa}
```
Estima a fração de reads atribuída a cada táxon.   

___________________

Vamos começar criando um diretório específico com ```mkdir SingleM```, e então instalando o SingleM com ```conda install -c conda-forge -c bioconda --strict-channel-priority singlem```, em que:
- ```-c conda-forge```: indica o canal conda-forge em que está o pacote singlem
- ```-c bioconda```: indica o canal bioconda em que está o pacote singlem
- ```--strict-channel-priority```: evita que o Conda misture versões de pacotes dos canais diferentes de maneira desnecessária
- ```singlem```: pacote desejado


Com o SingleM instalado, vamos primeiro produzir o perfil taxonômico:
```
singlem pipe \
    -1 /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-PM-1-A_R1.fastq.gz \
    -2 /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-PM-1-A_R2.fastq.gz \
    --otu-table /labgenomaarea2/valentina.pavelecini/EHI/SingleM/PM-1-A_otu_table.tsv \
    --taxonomic-profile /labgenomaarea2/valentina.pavelecini/EHI/SingleM/PM-1-A_taxonomic_profile.tsv \
    --threads 2
```

Erro:
```
(/labgenomaarea2/valentina.pavelecini/conda_envs/EHI) valentina.pavelecini@pantanal:/labgenomaarea2/valentina.pavelecini/EHI/SingleM$ singlem pipe \
>     -1 /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-PM-1-A_R1.fastq.gz \
>     -2 /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-PM-1-A_R2.fastq.gz \
>     --otu-table /labgenomaarea2/valentina.pavelecini/EHI/SingleM/PM-1-A_otu_table.tsv \
>     --taxonomic-profile /labgenomaarea2/valentina.pavelecini/EHI/SingleM/PM-1-A_taxonomic_profile.tsv \
>     --threads 2
2026/09/11 10:00:47 AM INFO: SingleM v0.21.4
Traceback (most recent call last):
  File "/labgenomaarea2/valentina.pavelecini/conda_envs/EHI/bin/singlem", line 8, in <module>
    sys.exit(main())
             ^^^^^^
  File "/labgenomaarea2/valentina.pavelecini/conda_envs/EHI/lib/python3.12/site-packages/singlem/main.py", line 780, in main
    singlem.pipe.SearchPipe().run(
  File "/labgenomaarea2/valentina.pavelecini/conda_envs/EHI/lib/python3.12/site-packages/singlem/pipe.py", line 67, in run
    metapackage = self._parse_packages_or_metapackage(**kwargs)
                  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/labgenomaarea2/valentina.pavelecini/conda_envs/EHI/lib/python3.12/site-packages/singlem/pipe.py", line 117, in _parse_packages_or_metapackage
    return Metapackage.acquire_default()
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/labgenomaarea2/valentina.pavelecini/conda_envs/EHI/lib/python3.12/site-packages/singlem/metapackage.py", line 156, in acquire_default
    backpack = Metapackage.acquire_default_backpack()
               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/labgenomaarea2/valentina.pavelecini/conda_envs/EHI/lib/python3.12/site-packages/singlem/metapackage.py", line 113, in acquire_default_backpack
    raise Exception("The {} environment variable, which points to the default data directory, is not set. To download the default SingleM metapackage, use 'singlem data'. The metapackage can also be downloaded manually from https://doi.org/{}".format(DATA_ENVIRONMENT_VARIABLE, DATA_DOI))
Exception: The SINGLEM_METAPACKAGE_PATH environment variable, which points to the default data directory, is not set. To download the default SingleM metapackage, use 'singlem data'. The metapackage can also be downloaded manually from https://doi.org/10.5281/zenodo.5739611
```

Não temos o metapackage, ou seja, o banco de referência para o SingleM analisar os reads. Começaremos criando uma pasta com ```mkdir metapackage```. O SingleM 0.21.4 possui uma função própria para baixar o metapackage padrão, então usaremos ela:
```
singlem data \
    --output-directory /labgenomaarea2/valentina.pavelecini/EHI/SingleM/metapackage
```
Output:
```
2026/09/11 10:33:51 AM INFO: Extracting files from archive...
2026/09/11 10:39:01 AM INFO: Verifying version and checksums...
2026/09/11 10:40:54 AM INFO: Verification success.
2026/09/11 10:40:54 AM INFO: Finished downloading data
2026/09/11 10:40:54 AM INFO: The environment variable SINGLEM_METAPACKAGE_PATH can now be set to /labgenomaarea2/valentina.pavelecini/EHI/SingleM/metapackage
2026/09/11 10:40:54 AM INFO: For instance, the following can be included in your .bashrc (requires logout and login after inclusion):
2026/09/11 10:40:54 AM INFO: export SINGLEM_METAPACKAGE_PATH='/labgenomaarea2/valentina.pavelecini/EHI/SingleM/metapackage/S6.5.0.GTDB_r232.metapackage_20260319.smpkg.zb'
```

Feito isso, definiremos a variável de ambiente "SINGLEM_METAPACKAGE_PATH" como indicado pelo próprio SingleM na mensagem após a instalação
```
export SINGLEM_METAPACKAGE_PATH='/labgenomaarea2/valentina.pavelecini/EHI/SingleM/metapackage/S6.5.0.GTDB_r232.metapackage_20260319.smpkg.zb'
```
_NOTA: a variável precisa ser definida em toda sessão. Ou seja, toda vez que deslogar do LAD e logar de novo, esse comando precisa ser rodado novamente._

E então tentaremos produzir o perfil taxonômico de novo:
```
singlem pipe \
    -1 /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-PM-1-A_R1.fastq.gz \
    -2 /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-PM-1-A_R2.fastq.gz \
    --otu-table /labgenomaarea2/valentina.pavelecini/EHI/SingleM/PM-1-A_otu_table.tsv \
    --taxonomic-profile /labgenomaarea2/valentina.pavelecini/EHI/SingleM/PM-1-A_taxonomic_profile.tsv \
    --threads 2
```

_____________________________

## 18/09/2026

Output do SingleM Pipe:
```
2026/09/18 09:03:18 AM INFO: SingleM v0.21.4
2026/09/18 09:03:18 AM INFO: Retrieval successful. Location of backpack is: /labgenomaarea2/valentina.pavelecini/EHI/SingleM/metapackage/S6.5.0.GTDB_r232.metapackage_20260319.smpkg.zb
2026/09/18 09:03:19 AM INFO: Loaded 59 SingleM packages
2026/09/18 09:03:25 AM INFO: Using as input 1 different pairs of sequence files e.g. /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-PM-1-A_R1.fastq.gz & /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-PM-1-A_R2.fastq.gz
2026/09/18 09:03:25 AM INFO: DIAMOND version: diamond version 2.2.6
2026/09/18 09:03:25 AM INFO: Filtering sequence files through DIAMOND blastx
2026/09/18 09:03:25 AM INFO: Filtering TF-2587-PM-1-A_R1.fastq.gz
2026/09/18 09:11:55 AM INFO: Found 27634 hits for TF-2587-PM-1-A_R1.fastq.gz
2026/09/18 09:11:55 AM INFO: Filtering TF-2587-PM-1-A_R2.fastq.gz
2026/09/18 09:20:20 AM INFO: Found 27544 hits for TF-2587-PM-1-A_R2.fastq.gz
2026/09/18 09:20:20 AM INFO: Finished DIAMOND prefilter phase
2026/09/18 09:20:20 AM INFO: Assigning sequences to SingleM packages with DIAMOND ..
2026/09/18 09:20:21 AM INFO: Extracting reads from 1 sample(s) across 59 package(s) using 2 thread(s)
Extracting reads: 100%|█████████████████████████████████████████████████████████| 59/59 [00:59<00:00,  1.01s/sample×pkg]
2026/09/18 09:21:21 AM INFO: After read extraction, 26645 sequence(s) remain
2026/09/18 09:21:21 AM INFO: Running taxonomic assignment ..
2026/09/18 09:21:21 AM INFO: Assigning taxonomy by singlem query ..
Querying taxonomy: 100%|█████████████████████████████████████████████████████████| 5419/5419 [01:49<00:00, 49.45 seqs/s]
2026/09/18 09:23:11 AM INFO: Finished running singlem query-based taxonomic assignment, now running diamond using 2 thread(s) ..
2026/09/18 09:23:11 AM INFO: Assigning taxonomy with DIAMOND blastx to 7750 OTUs (The 81.6% that were not assigned by smafa) ..
DIAMOND taxonomy: 100%|████████████████████████████████████████████████████████████| 118/118 [35:45<00:00, 18.18s/chunk]
2026/09/18 09:58:57 AM INFO: Finished running taxonomic assignment
2026/09/18 09:59:58 AM INFO: Finished
2026/09/18 09:59:58 AM INFO: Writing /labgenomaarea2/valentina.pavelecini/EHI/SingleM/PM-1-A_otu_table.tsv
2026/09/18 09:59:58 AM INFO: Writing taxonomic profile to /labgenomaarea2/valentina.pavelecini/EHI/SingleM/PM-1-A_taxonomic_profile.tsv
2026/09/18 09:59:58 AM INFO: Using minimum taxon coverage of 0.35
2026/09/18 09:59:58 AM INFO: Removing off-target OTUs from TF-2587-PM-1-A_R1
2026/09/18 09:59:58 AM INFO: Found 12555.54 assigned and 0.00 unassigned OTU coverage units
2026/09/18 09:59:58 AM INFO: After removing off-target OTUs, found 10763.80 assigned and 0.00 unassigned OTU coverage units
2026/09/18 09:59:58 AM INFO: Total OTU coverage by query: 2100.4268488180187
2026/09/18 09:59:58 AM INFO: Total OTU coverage by diamond: 8663.371088033391
2026/09/18 09:59:58 AM INFO: Applying species-wise expectation maximization algorithm to OTU table
2026/09/18 09:59:58 AM INFO: Found 2 species uniquely hitting >= 10 marker genes
2026/09/18 10:00:00 AM INFO: Species-wise EM converged in 55 steps
2026/09/18 10:00:00 AM INFO: Gathering equivalence classes
2026/09/18 10:00:00 AM INFO: Demultiplexing OTU table
2026/09/18 10:00:01 AM INFO: Finished expectation maximization
2026/09/18 10:00:01 AM INFO: Converting DIAMOND IDs to taxons
2026/09/18 10:00:36 AM INFO: Converted 3718 Diamond-assigned OTU taxon_ids to taxon strings
2026/09/18 10:00:36 AM INFO: Applying genus-wise expectation maximization algorithm to OTU table
2026/09/18 10:00:39 AM INFO: Genus-wise EM converged in 50 steps
2026/09/18 10:00:39 AM INFO: Gathering equivalence classes
2026/09/18 10:00:39 AM INFO: Demultiplexing OTU table
2026/09/18 10:00:40 AM INFO: Finished genus expectation maximization
2026/09/18 10:00:40 AM INFO: Total profile coverage after condense domain to species: 298.77297293837285
2026/09/18 10:00:40 AM INFO: Total profile coverage after push down: 298.7729729383729
2026/09/18 10:00:40 AM INFO: Taxonomic level coverage:
2026/09/18 10:00:40 AM INFO: kingdom:   6.09%   2 taxons
2026/09/18 10:00:40 AM INFO: phylum:    5.86%   20 taxons
2026/09/18 10:00:40 AM INFO: class:     15.77%  37 taxons
2026/09/18 10:00:40 AM INFO: order:     16.47%  57 taxons
2026/09/18 10:00:40 AM INFO: family:    21.72%  66 taxons
2026/09/18 10:00:40 AM INFO: genus:     31.74%  51 taxons
2026/09/18 10:00:40 AM INFO: species:   2.35%   7 taxons
2026/09/18 10:00:40 AM INFO: Finished condense
```

Como a análise do SingleM é demorada (~1h), não iremos criar um bash para rodar todas as amostras de uma vez só, e sim rodar cada uma manualmente:

Já fizemos PM-1-A

### PM-5-A
```
singlem pipe \
    -1 /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-PM-5-A_S7_L001_R1.fastq.gz \
    -2 /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-PM-5-A_S7_L001_R2.fastq.gz \
    --otu-table /labgenomaarea2/valentina.pavelecini/EHI/SingleM/PM-5-A_otu_table.tsv \
    --taxonomic-profile /labgenomaarea2/valentina.pavelecini/EHI/SingleM/PM-5-A_taxonomic_profile.tsv \
    --threads 2
```

Output:
```
2026/09/18 10:23:11 AM INFO: SingleM v0.21.4
2026/09/18 10:23:11 AM INFO: Retrieval successful. Location of backpack is: /labgenomaarea2/valentina.pavelecini/EHI/SingleM/metapackage/S6.5.0.GTDB_r232.metapackage_20260319.smpkg.zb
2026/09/18 10:23:12 AM INFO: Loaded 59 SingleM packages
2026/09/18 10:23:18 AM INFO: Using as input 1 different pairs of sequence files e.g. /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-PM-5-A_S7_L001_R1.fastq.gz & /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-PM-5-A_S7_L001_R2.fastq.gz
2026/09/18 10:23:18 AM INFO: DIAMOND version: diamond version 2.2.6
2026/09/18 10:23:18 AM INFO: Filtering sequence files through DIAMOND blastx
2026/09/18 10:23:18 AM INFO: Filtering TF-2587-PM-5-A_S7_L001_R1.fastq.gz
2026/09/18 10:36:23 AM INFO: Found 38396 hits for TF-2587-PM-5-A_S7_L001_R1.fastq.gzHI/SingleM/metapackage
2026/09/18 10:36:24 AM INFO: Filtering TF-2587-PM-5-A_S7_L001_R2.fastq.gz
2026/09/18 10:49:22 AM INFO: Found 37997 hits for TF-2587-PM-5-A_S7_L001_R2.fastq.gz
2026/09/18 10:49:22 AM INFO: Finished DIAMOND prefilter phase
2026/09/18 10:49:22 AM INFO: Assigning sequences to SingleM packages with DIAMOND ..
2026/09/18 10:49:23 AM INFO: Extracting reads from 1 sample(s) across 59 package(s) using 2 thread(s)
Extracting reads: 100%|█████████████████████████████████████████████████████████████████████████████████████| 59/59 [01:13<00:00,  1.24s/sample×pkg]
2026/09/18 10:50:37 AM INFO: After read extraction, 35425 sequence(s) remain
2026/09/18 10:50:37 AM INFO: Running taxonomic assignment ..
2026/09/18 10:50:37 AM INFO: Assigning taxonomy by singlem query ..
Querying taxonomy: 100%|█████████████████████████████████████████████████████████████████████████████████████| 4450/4450 [01:37<00:00, 45.78 seqs/s]
2026/09/18 10:52:15 AM INFO: Finished running singlem query-based taxonomic assignment, now running diamond using 2 thread(s) ..
2026/09/18 10:52:15 AM INFO: Assigning taxonomy with DIAMOND blastx to 9745 OTUs (The 76.1% that were not assigned by smafa) ..
DIAMOND taxonomy: 100%|████████████████████████████████████████████████████████████████████████████████████████| 118/118 [45:06<00:00, 22.94s/chunk]
2026/09/18 11:37:22 AM INFO: Finished running taxonomic assignment
2026/09/18 11:38:27 AM INFO: Finished
2026/09/18 11:38:27 AM INFO: Writing /labgenomaarea2/valentina.pavelecini/EHI/SingleM/PM-5-A_otu_table.tsv
2026/09/18 11:38:27 AM INFO: Writing taxonomic profile to /labgenomaarea2/valentina.pavelecini/EHI/SingleM/PM-5-A_taxonomic_profile.tsv
2026/09/18 11:38:27 AM INFO: Using minimum taxon coverage of 0.35
2026/09/18 11:38:27 AM INFO: Removing off-target OTUs from TF-2587-PM-5-A_S7_L001_R1
2026/09/18 11:38:27 AM INFO: Found 17498.27 assigned and 0.00 unassigned OTU coverage units
2026/09/18 11:38:27 AM INFO: After removing off-target OTUs, found 14662.68 assigned and 0.00 unassigned OTU coverage units
2026/09/18 11:38:27 AM INFO: Total OTU coverage by query: 3853.890936463486
2026/09/18 11:38:27 AM INFO: Total OTU coverage by diamond: 10808.791954260652
2026/09/18 11:38:27 AM INFO: Applying species-wise expectation maximization algorithm to OTU table
2026/09/18 11:38:27 AM INFO: Found 3 species uniquely hitting >= 10 marker genes
2026/09/18 11:38:28 AM INFO: Species-wise EM converged in 63 steps
2026/09/18 11:38:28 AM INFO: Gathering equivalence classes
2026/09/18 11:38:28 AM INFO: Demultiplexing OTU table
2026/09/18 11:38:29 AM INFO: Finished expectation maximization
2026/09/18 11:38:29 AM INFO: Converting DIAMOND IDs to taxons
2026/09/18 11:39:01 AM INFO: Converted 2911 Diamond-assigned OTU taxon_ids to taxon strings
2026/09/18 11:39:01 AM INFO: Applying genus-wise expectation maximization algorithm to OTU table
2026/09/18 11:39:04 AM INFO: Genus-wise EM converged in 59 steps
2026/09/18 11:39:04 AM INFO: Gathering equivalence classes
2026/09/18 11:39:04 AM INFO: Demultiplexing OTU table
2026/09/18 11:39:04 AM INFO: Finished genus expectation maximization
2026/09/18 11:39:04 AM INFO: Total profile coverage after condense domain to species: 408.681532696004
2026/09/18 11:39:04 AM INFO: Total profile coverage after push down: 408.68153269600407
2026/09/18 11:39:04 AM INFO: Taxonomic level coverage:
2026/09/18 11:39:04 AM INFO: kingdom:   6.78%   2 taxons
2026/09/18 11:39:04 AM INFO: phylum:    5.67%   19 taxons
2026/09/18 11:39:04 AM INFO: class:     12.14%  32 taxons
2026/09/18 11:39:04 AM INFO: order:     11.61%  41 taxons
2026/09/18 11:39:04 AM INFO: family:    13.33%  51 taxons
2026/09/18 11:39:04 AM INFO: genus:     44.03%  51 taxons
2026/09/18 11:39:04 AM INFO: species:   6.44%   12 taxons
2026/09/18 11:39:04 AM INFO: Finished condense
```


## 23/09/2026

Mudei um pouco a estratégia hoje: vou usar uma janela do tmux pra rodar as próximas amostras. O tmux é uma janela em que os processos nela continuam rodando mesmo que o computador desligue (tipo um job só que mais interativo). Comecei criando uma janela para essas análises do SingleM com ```tmux new -s SingleM-EHI```.  Defini os parâmetros dessa sessão com ```srun -N 2 -t 01:30:00 -n 8 --pty bash -i```, em que ```-N 2``` separa 2 máquinas pro processo, ```-t 01:30:00``` aloca essa sessão por 1 hora e meia, ```-n``` define o limite pra 8 threads (ao invés de 2), e ```--pty bash -i``` indica que esses parâmetros vão valer para qualquer comando realizado nessa janela na próxima 1h30. O objetivo é diminuir um pouco o tempo de análise. Além disso, pra iniciar cada processo, sempre lembremos de rodar ```export SINGLEM_METAPACKAGE_PATH='/labgenomaarea2/valentina.pavelecini/EHI/SingleM/metapackage/S6.5.0.GTDB_r232.metapackage_20260319.smpkg.zb'``` para cada nova sessão.


Além disso, vamos adicionar mais uma linha aos comandos seguintes: ```&> /labgenomaarea2/valentina.pavelecini/conda_envs/EHI/SingleM/outputs/output-XX-X.log```, que salva os outputs dos processos em um arquivo .log.

Por isso, os próximos comandos serão um pouco diferentes dos anteriores:


### RF-1-B
```
singlem pipe \
    -1 /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-RF-1-B_S1_L001_R1.fastq.gz \
    -2 /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-RF-1-B_S1_L001_R2.fastq.gz \
    --otu-table /labgenomaarea2/valentina.pavelecini/EHI/SingleM/RF-1-B_otu_table.tsv \
    --taxonomic-profile /labgenomaarea2/valentina.pavelecini/EHI/SingleM/RF-1-B_taxonomic_profile.tsv \
    --threads 2 \
    &> /labgenomaarea2/valentina.pavelecini/EHI/SingleM/outputs/output-RF-1.log
```

Output:
```
2026/09/25 08:35:23 AM INFO: SingleM v0.21.4
2026/09/25 08:35:23 AM INFO: Retrieval successful. Location of backpack is: /labgenomaarea2/valentina.pavelecini/EHI/SingleM/metapackage/S6.5.0.GTDB_r232.metapackage_20260319.smpkg.zb
2026/09/25 08:35:24 AM INFO: Loaded 59 SingleM packages
2026/09/25 08:35:29 AM INFO: Using as input 1 different pairs of sequence files e.g. /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-RF-1-B_S1_L001_R1.fastq.gz & /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-RF-1-B_S1_L001_R2.fastq.gz
2026/09/25 08:35:29 AM INFO: DIAMOND version: diamond version 2.2.6
2026/09/25 08:35:29 AM INFO: Filtering sequence files through DIAMOND blastx
2026/09/25 08:35:29 AM INFO: Filtering TF-2587-RF-1-B_S1_L001_R1.fastq.gz
2026/09/25 09:10:54 AM INFO: Found 101749 hits for TF-2587-RF-1-B_S1_L001_R1.fastq.gz
2026/09/25 09:10:55 AM INFO: Filtering TF-2587-RF-1-B_S1_L001_R2.fastq.gz
2026/09/25 09:46:07 AM INFO: Found 101149 hits for TF-2587-RF-1-B_S1_L001_R2.fastq.gz
2026/09/25 09:46:08 AM INFO: Finished DIAMOND prefilter phase
2026/09/25 09:46:08 AM INFO: Assigning sequences to SingleM packages with DIAMOND ..
2026/09/25 09:46:10 AM INFO: Extracting reads from 1 sample(s) across 59 package(s) using 2 thread(s)
Extracting reads: 100%|█████████████████████████████████████████████████████████| 59/59 [02:44<00:00,  2.78s/sample×pkg]
2026/09/25 09:48:56 AM INFO: After read extraction, 91500 sequence(s) remain
2026/09/25 09:48:56 AM INFO: Running taxonomic assignment ..
2026/09/25 09:48:56 AM INFO: Assigning taxonomy by singlem query ..
Querying taxonomy: 100%|███████████████████████████████████████████████████████| 10348/10348 [02:04<00:00, 82.94 seqs/s]
2026/09/25 09:51:02 AM INFO: Finished running singlem query-based taxonomic assignment, now running diamond using 2 thread(s) ..
2026/09/25 09:51:02 AM INFO: Assigning taxonomy with DIAMOND blastx to 26615 OTUs (The 79.4% that were not assigned by smafa) ..
DIAMOND taxonomy: 100%|██████████████████████████████████████████████████████████| 118/118 [2:09:12<00:00, 65.70s/chunk]
2026/09/25 12:00:16 PM INFO: Finished running taxonomic assignment
2026/09/25 12:01:57 PM INFO: Finished
2026/09/25 12:01:58 PM INFO: Writing /labgenomaarea2/valentina.pavelecini/EHI/SingleM/RF-1-B_otu_table.tsv
2026/09/25 12:01:58 PM INFO: Writing taxonomic profile to /labgenomaarea2/valentina.pavelecini/EHI/SingleM/RF-1-B_taxonomic_profile.tsv
2026/09/25 12:01:58 PM INFO: Using minimum taxon coverage of 0.35
2026/09/25 12:01:58 PM INFO: Removing off-target OTUs from TF-2587-RF-1-B_S1_L001_R1
2026/09/25 12:01:58 PM INFO: Found 45275.84 assigned and 0.00 unassigned OTU coverage units
2026/09/25 12:01:58 PM INFO: After removing off-target OTUs, found 38591.91 assigned and 0.00 unassigned OTU coverage units
2026/09/25 12:01:59 PM INFO: Total OTU coverage by query: 8280.254504643764
2026/09/25 12:01:59 PM INFO: Total OTU coverage by diamond: 30311.65820502048
2026/09/25 12:01:59 PM INFO: Applying species-wise expectation maximization algorithm to OTU table
2026/09/25 12:01:59 PM INFO: Found 1 species uniquely hitting >= 10 marker genes
2026/09/25 12:02:03 PM INFO: Species-wise EM converged in 86 steps
2026/09/25 12:02:03 PM INFO: Gathering equivalence classes
2026/09/25 12:02:03 PM INFO: Demultiplexing OTU table
2026/09/25 12:02:03 PM INFO: Finished expectation maximization
2026/09/25 12:02:03 PM INFO: Converting DIAMOND IDs to taxons
2026/09/25 12:02:50 PM INFO: Converted 7134 Diamond-assigned OTU taxon_ids to taxon strings
2026/09/25 12:02:50 PM INFO: Applying genus-wise expectation maximization algorithm to OTU table
2026/09/25 12:02:56 PM INFO: Genus-wise EM converged in 59 steps
2026/09/25 12:02:56 PM INFO: Gathering equivalence classes
2026/09/25 12:02:56 PM INFO: Demultiplexing OTU table
2026/09/25 12:02:57 PM INFO: Finished genus expectation maximization
2026/09/25 12:02:57 PM INFO: Total profile coverage after condense domain to species: 1076.3839064635652
2026/09/25 12:02:57 PM INFO: Total profile coverage after push down: 1076.3839064635658
2026/09/25 12:02:57 PM INFO: Taxonomic level coverage:
2026/09/25 12:02:57 PM INFO: kingdom:   3.07%   2 taxons
2026/09/25 12:02:57 PM INFO: phylum:    3.58%   22 taxons
2026/09/25 12:02:57 PM INFO: class:     8.07%   46 taxons
2026/09/25 12:02:57 PM INFO: order:     9.89%   76 taxons
2026/09/25 12:02:57 PM INFO: family:    17.28%  122 taxons
2026/09/25 12:02:57 PM INFO: genus:     50.00%  100 taxons
2026/09/25 12:02:57 PM INFO: species:   8.10%   16 taxons
2026/09/25 12:02:57 PM INFO: Finished condense
```

## 30/09/2026

Usar tmux não deu certo então azar. Voltando a fazer normal:


### RF-4-B
```
singlem pipe \
    -1 /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-RF-4-B_S3_L001_R1.fastq.gz \
    -2 /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-RF-4-B_S3_L001_R2.fastq.gz \
    --otu-table /labgenomaarea2/valentina.pavelecini/EHI/SingleM/RF-4-B_otu_table.tsv \
    --taxonomic-profile /labgenomaarea2/valentina.pavelecini/EHI/SingleM/RF-4-B_taxonomic_profile.tsv \
    --threads 2
```
Output:
```
2026/09/30 08:46:40 AM INFO: SingleM v0.21.4
2026/09/30 08:46:40 AM INFO: Retrieval successful. Location of backpack is: /labgenomaarea2/valentina.pavelecini/EHI/SingleM/metapackage/S6.5.0.GTDB_r232.metapackage_20260319.smpkg.zb
2026/09/30 08:46:41 AM INFO: Loaded 59 SingleM packages
2026/09/30 08:46:46 AM INFO: Using as input 1 different pairs of sequence files e.g. /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-RF-4-B_S3_L001_R1.fastq.gz & /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-RF-4-B_S3_L001_R2.fastq.gz
2026/09/30 08:46:47 AM INFO: DIAMOND version: diamond version 2.2.6
2026/09/30 08:46:47 AM INFO: Filtering sequence files through DIAMOND blastx
2026/09/30 08:46:47 AM INFO: Filtering TF-2587-RF-4-B_S3_L001_R1.fastq.gz
2026/09/30 08:58:58 AM INFO: Found 35535 hits for TF-2587-RF-4-B_S3_L001_R1.fastq.gz
2026/09/30 08:58:58 AM INFO: Filtering TF-2587-RF-4-B_S3_L001_R2.fastq.gz
2026/09/30 09:11:32 AM INFO: Found 35744 hits for TF-2587-RF-4-B_S3_L001_R2.fastq.gz
2026/09/30 09:11:32 AM INFO: Finished DIAMOND prefilter phase
2026/09/30 09:11:32 AM INFO: Assigning sequences to SingleM packages with DIAMOND ..
2026/09/30 09:11:33 AM INFO: Extracting reads from 1 sample(s) across 59 package(s) using 2 thread(s)
Extracting reads: 100%|█████████████████████████████████████████████████████████| 59/59 [01:06<00:00,  1.13s/sample×pkg]
2026/09/30 09:12:40 AM INFO: After read extraction, 33278 sequence(s) remain
2026/09/30 09:12:40 AM INFO: Running taxonomic assignment ..
2026/09/30 09:12:40 AM INFO: Assigning taxonomy by singlem query ..
Querying taxonomy: 100%|█████████████████████████████████████████████████████████| 2891/2891 [01:30<00:00, 31.80 seqs/s]
2026/09/30 09:14:12 AM INFO: Finished running singlem query-based taxonomic assignment, now running diamond using 2 thread(s) ..
2026/09/30 09:14:12 AM INFO: Assigning taxonomy with DIAMOND blastx to 9289 OTUs (The 80.6% that were not assigned by smafa) ..
DIAMOND taxonomy: 100%|██████████████████████████████████████████████████████████| 118/118 [1:05:40<00:00, 33.39s/chunk]
2026/09/30 10:19:53 AM INFO: Finished running taxonomic assignment
2026/09/30 10:20:55 AM INFO: Finished
2026/09/30 10:20:56 AM INFO: Writing /labgenomaarea2/valentina.pavelecini/EHI/SingleM/RF-4-B_otu_table.tsv
2026/09/30 10:20:56 AM INFO: Writing taxonomic profile to /labgenomaarea2/valentina.pavelecini/EHI/SingleM/RF-4-B_taxonomic_profile.tsv
2026/09/30 10:20:56 AM INFO: Using minimum taxon coverage of 0.35
2026/09/30 10:20:56 AM INFO: Removing off-target OTUs from TF-2587-RF-4-B_S3_L001_R1
2026/09/30 10:20:56 AM INFO: Found 15457.02 assigned and 0.00 unassigned OTU coverage units
2026/09/30 10:20:56 AM INFO: After removing off-target OTUs, found 13080.66 assigned and 0.00 unassigned OTU coverage units
2026/09/30 10:20:56 AM INFO: Total OTU coverage by query: 2687.6175338801227
2026/09/30 10:20:56 AM INFO: Total OTU coverage by diamond: 10393.037546219548
2026/09/30 10:20:56 AM INFO: Applying species-wise expectation maximization algorithm to OTU table
2026/09/30 10:20:56 AM INFO: Found 0 species uniquely hitting >= 10 marker genes
2026/09/30 10:20:57 AM INFO: Species-wise EM converged in 56 steps
2026/09/30 10:20:57 AM INFO: Gathering equivalence classes
2026/09/30 10:20:57 AM INFO: Demultiplexing OTU table
2026/09/30 10:20:58 AM INFO: Finished expectation maximization
2026/09/30 10:20:58 AM INFO: Converting DIAMOND IDs to taxons
2026/09/30 10:21:19 AM INFO: Converted 1900 Diamond-assigned OTU taxon_ids to taxon strings
2026/09/30 10:21:19 AM INFO: Applying genus-wise expectation maximization algorithm to OTU table
2026/09/30 10:21:23 AM INFO: Genus-wise EM converged in 113 steps
2026/09/30 10:21:23 AM INFO: Gathering equivalence classes
2026/09/30 10:21:23 AM INFO: Demultiplexing OTU table
2026/09/30 10:21:23 AM INFO: Finished genus expectation maximization
2026/09/30 10:21:23 AM INFO: Total profile coverage after condense domain to species: 354.8656978489279
2026/09/30 10:21:23 AM INFO: Total profile coverage after push down: 354.865697848928
2026/09/30 10:21:23 AM INFO: Taxonomic level coverage:
2026/09/30 10:21:23 AM INFO: kingdom:   12.19%  2 taxons
2026/09/30 10:21:23 AM INFO: phylum:    8.27%   14 taxons
2026/09/30 10:21:23 AM INFO: class:     14.42%  18 taxons
2026/09/30 10:21:23 AM INFO: order:     10.41%  24 taxons
2026/09/30 10:21:23 AM INFO: family:    12.20%  27 taxons
2026/09/30 10:21:23 AM INFO: genus:     36.89%  30 taxons
2026/09/30 10:21:23 AM INFO: species:   5.61%   8 taxons
2026/09/30 10:21:23 AM INFO: Finished condense
```


### RF-5-B
```
singlem pipe \
    -1 /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-RF-5-B_S4_L001_R1.fastq.gz \
    -2 /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-RF-5-B_S4_L001_R2.fastq.gz \
    --otu-table /labgenomaarea2/valentina.pavelecini/EHI/SingleM/RF-5-B_otu_table.tsv \
    --taxonomic-profile /labgenomaarea2/valentina.pavelecini/EHI/SingleM/RF-5-B_taxonomic_profile.tsv \
    --threads 2
```

Output:
```
2026/09/28 02:42:10 PM INFO: SingleM v0.21.4
2026/09/28 02:42:10 PM INFO: Retrieval successful. Location of backpack is: /labgenomaarea2/valentina.pavelecini/EHI/SingleM/metapackage/S6.5.0.GTDB_r232.metapackage_20260319.smpkg.zb
2026/09/28 02:42:10 PM INFO: Loaded 59 SingleM packages
2026/09/28 02:42:15 PM INFO: Using as input 1 different pairs of sequence files e.g. /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-RF-5-B_S4_L001_R1.fastq.gz & /labgenoma>
2026/09/28 02:42:15 PM INFO: DIAMOND version: diamond version 2.2.6
2026/09/28 02:42:15 PM INFO: Filtering sequence files through DIAMOND blastx
2026/09/28 02:42:15 PM INFO: Filtering TF-2587-RF-5-B_S4_L001_R1.fastq.gz
^M                                                                                ^M2026/09/28 03:05:38 PM INFO: Found 81757 hits for TF-2587-RF-5-B_S4_L001_R1.fastq.gz
2026/09/28 03:05:38 PM INFO: Filtering TF-2587-RF-5-B_S4_L001_R2.fastq.gz
^M                                                                                ^M2026/09/28 03:29:00 PM INFO: Found 81504 hits for TF-2587-RF-5-B_S4_L001_R2.fastq.gz
2026/09/28 03:29:00 PM INFO: Finished DIAMOND prefilter phase
2026/09/28 03:29:00 PM INFO: Assigning sequences to SingleM packages with DIAMOND ..
2026/09/28 03:29:00 PM INFO: Extracting reads from 1 sample(s) across 59 package(s) using 2 thread(s)
2026/09/28 03:30:10 PM INFO: Finished extracting reads for chunk 1 of 59
2026/09/28 03:30:10 PM INFO: Finished extracting reads for chunk 2 of 59
2026/09/28 03:30:13 PM INFO: Finished extracting reads for chunk 3 of 59
2026/09/28 03:30:14 PM INFO: Finished extracting reads for chunk 4 of 59
2026/09/28 03:30:16 PM INFO: Finished extracting reads for chunk 5 of 59
2026/09/28 03:30:17 PM INFO: Finished extracting reads for chunk 6 of 59
2026/09/28 03:30:19 PM INFO: Finished extracting reads for chunk 7 of 59
2026/09/28 03:30:20 PM INFO: Finished extracting reads for chunk 8 of 59
2026/09/28 03:30:22 PM INFO: Finished extracting reads for chunk 9 of 59
2026/09/28 03:30:23 PM INFO: Finished extracting reads for chunk 10 of 59
2026/09/28 03:30:25 PM INFO: Finished extracting reads for chunk 11 of 59
2026/09/28 03:30:26 PM INFO: Finished extracting reads for chunk 12 of 59
2026/09/28 03:30:28 PM INFO: Finished extracting reads for chunk 13 of 59
2026/09/28 03:30:29 PM INFO: Finished extracting reads for chunk 14 of 59
2026/09/28 03:30:30 PM INFO: Finished extracting reads for chunk 15 of 59
2026/09/28 03:30:32 PM INFO: Finished extracting reads for chunk 16 of 59
2026/09/28 03:30:32 PM INFO: Finished extracting reads for chunk 17 of 59
2026/09/28 03:30:33 PM INFO: Finished extracting reads for chunk 18 of 59
2026/09/28 03:30:33 PM INFO: Finished extracting reads for chunk 19 of 59
2026/09/28 03:30:33 PM INFO: Finished extracting reads for chunk 20 of 59
2026/09/28 03:30:34 PM INFO: Finished extracting reads for chunk 21 of 59
2026/09/28 03:30:35 PM INFO: Finished extracting reads for chunk 22 of 59
2026/09/28 03:30:35 PM INFO: Finished extracting reads for chunk 23 of 59
2026/09/28 03:30:36 PM INFO: Finished extracting reads for chunk 24 of 59
2026/09/28 03:30:36 PM INFO: Finished extracting reads for chunk 25 of 59
2026/09/28 03:30:36 PM INFO: Finished extracting reads for chunk 26 of 59
2026/09/28 03:30:37 PM INFO: Finished extracting reads for chunk 27 of 59
2026/09/28 03:30:36 PM INFO: Finished extracting reads for chunk 26 of 59
2026/09/28 03:30:37 PM INFO: Finished extracting reads for chunk 27 of 59
2026/09/28 03:30:37 PM INFO: Finished extracting reads for chunk 28 of 59
2026/09/28 03:30:38 PM INFO: Finished extracting reads for chunk 29 of 59
2026/09/28 03:30:38 PM INFO: Finished extracting reads for chunk 30 of 59
2026/09/28 03:30:39 PM INFO: Finished extracting reads for chunk 31 of 59
2026/09/28 03:30:39 PM INFO: Finished extracting reads for chunk 32 of 59
2026/09/28 03:30:40 PM INFO: Finished extracting reads for chunk 33 of 59
2026/09/28 03:30:40 PM INFO: Finished extracting reads for chunk 34 of 59
2026/09/28 03:30:41 PM INFO: Finished extracting reads for chunk 35 of 59
2026/09/28 03:30:41 PM INFO: Finished extracting reads for chunk 36 of 59
2026/09/28 03:30:42 PM INFO: Finished extracting reads for chunk 37 of 59
2026/09/28 03:30:44 PM INFO: Finished extracting reads for chunk 38 of 59
2026/09/28 03:30:46 PM INFO: Finished extracting reads for chunk 39 of 59
2026/09/28 03:30:47 PM INFO: Finished extracting reads for chunk 40 of 59
2026/09/28 03:30:49 PM INFO: Finished extracting reads for chunk 41 of 59
2026/09/28 03:30:50 PM INFO: Finished extracting reads for chunk 42 of 59
2026/09/28 03:30:52 PM INFO: Finished extracting reads for chunk 43 of 59
2026/09/28 03:30:53 PM INFO: Finished extracting reads for chunk 44 of 59
2026/09/28 03:30:55 PM INFO: Finished extracting reads for chunk 45 of 59
2026/09/28 03:30:55 PM INFO: Finished extracting reads for chunk 46 of 59
2026/09/28 03:30:57 PM INFO: Finished extracting reads for chunk 47 of 59
2026/09/28 03:30:58 PM INFO: Finished extracting reads for chunk 48 of 59
2026/09/28 03:31:00 PM INFO: Finished extracting reads for chunk 49 of 59
2026/09/28 03:31:01 PM INFO: Finished extracting reads for chunk 50 of 59
2026/09/28 03:31:02 PM INFO: Finished extracting reads for chunk 51 of 59
2026/09/28 03:31:04 PM INFO: Finished extracting reads for chunk 52 of 59
2026/09/28 03:31:05 PM INFO: Finished extracting reads for chunk 53 of 59
2026/09/28 03:31:07 PM INFO: Finished extracting reads for chunk 54 of 59
2026/09/28 03:31:08 PM INFO: Finished extracting reads for chunk 55 of 59
2026/09/28 03:31:10 PM INFO: Finished extracting reads for chunk 56 of 59
2026/09/28 03:31:11 PM INFO: Finished extracting reads for chunk 57 of 59
2026/09/28 03:31:13 PM INFO: Finished extracting reads for chunk 58 of 59
2026/09/28 03:31:13 PM INFO: Finished extracting reads for chunk 59 of 59
2026/09/28 03:31:13 PM INFO: After read extraction, 73093 sequence(s) remain
2026/09/28 03:31:13 PM INFO: Running taxonomic assignment ..
2026/09/28 03:31:13 PM INFO: Assigning taxonomy by singlem query ..
2026/09/28 03:31:30 PM INFO: Progress: 189/7893 (2%)
2026/09/28 03:31:36 PM INFO: Progress: 1843/7893 (23%)
2026/09/28 03:31:43 PM INFO: Progress: 3538/7893 (44%)
2026/09/28 03:31:48 PM INFO: Progress: 5022/7893 (63%)
2026/09/28 03:31:53 PM INFO: Progress: 6579/7893 (83%)
2026/09/28 03:31:58 PM INFO: Finished running singlem query-based taxonomic assignment, now running diamond using 2 thread(s) ..
2026/09/28 03:31:58 PM INFO: Assigning taxonomy with DIAMOND blastx to 16790 OTUs (The 63.6% that were not assigned by smafa) ..
2026/09/28 03:34:44 PM INFO: Finished DIAMOND blastx chunk 1 of 118
2026/09/28 03:35:18 PM INFO: Finished DIAMOND blastx chunk 2 of 118
2026/09/28 03:36:11 PM INFO: Finished DIAMOND blastx chunk 3 of 118
2026/09/28 03:37:00 PM INFO: Finished DIAMOND blastx chunk 4 of 118
2026/09/28 03:39:05 PM INFO: Finished DIAMOND blastx chunk 5 of 118
2026/09/28 03:40:27 PM INFO: Finished DIAMOND blastx chunk 6 of 118
2026/09/28 03:41:46 PM INFO: Finished DIAMOND blastx chunk 7 of 118
2026/09/28 03:43:39 PM INFO: Finished DIAMOND blastx chunk 8 of 118
2026/09/28 03:44:29 PM INFO: Finished DIAMOND blastx chunk 9 of 118
2026/09/28 03:46:08 PM INFO: Finished DIAMOND blastx chunk 10 of 118
2026/09/28 03:47:48 PM INFO: Finished DIAMOND blastx chunk 11 of 118
2026/09/28 03:49:23 PM INFO: Finished DIAMOND blastx chunk 12 of 118
2026/09/28 03:50:24 PM INFO: Finished DIAMOND blastx chunk 13 of 118
2026/09/28 03:51:31 PM INFO: Finished DIAMOND blastx chunk 14 of 118
2026/09/28 03:53:40 PM INFO: Finished DIAMOND blastx chunk 15 of 118
2026/09/28 03:54:39 PM INFO: Finished DIAMOND blastx chunk 16 of 118
2026/09/28 03:55:20 PM INFO: Finished DIAMOND blastx chunk 17 of 118
2026/09/28 03:56:04 PM INFO: Finished DIAMOND blastx chunk 18 of 118
2026/09/28 03:58:02 PM INFO: Finished DIAMOND blastx chunk 19 of 118
2026/09/28 03:58:37 PM INFO: Finished DIAMOND blastx chunk 20 of 118
2026/09/28 04:00:42 PM INFO: Finished DIAMOND blastx chunk 21 of 118
2026/09/28 04:00:50 PM INFO: Finished DIAMOND blastx chunk 22 of 118
2026/09/28 04:02:37 PM INFO: Finished DIAMOND blastx chunk 23 of 118
2026/09/28 04:02:40 PM INFO: Finished DIAMOND blastx chunk 24 of 118
2026/09/28 04:04:21 PM INFO: Finished DIAMOND blastx chunk 25 of 118
2026/09/28 04:04:47 PM INFO: Finished DIAMOND blastx chunk 26 of 118
2026/09/28 04:06:00 PM INFO: Finished DIAMOND blastx chunk 27 of 118
2026/09/28 04:06:22 PM INFO: Finished DIAMOND blastx chunk 28 of 118
2026/09/28 04:07:51 PM INFO: Finished DIAMOND blastx chunk 29 of 118
2026/09/28 04:08:22 PM INFO: Finished DIAMOND blastx chunk 30 of 118
2026/09/28 04:10:57 PM INFO: Finished DIAMOND blastx chunk 31 of 118
2026/09/28 04:10:57 PM INFO: Finished DIAMOND blastx chunk 31 of 118
2026/09/28 04:10:57 PM INFO: Finished DIAMOND blastx chunk 32 of 118
2026/09/28 04:10:58 PM INFO: Finished DIAMOND blastx chunk 33 of 118
2026/09/28 04:10:59 PM INFO: Finished DIAMOND blastx chunk 34 of 118
2026/09/28 04:11:00 PM INFO: Finished DIAMOND blastx chunk 35 of 118
2026/09/28 04:11:00 PM INFO: Finished DIAMOND blastx chunk 36 of 118
2026/09/28 04:11:01 PM INFO: Finished DIAMOND blastx chunk 37 of 118
2026/09/28 04:11:02 PM INFO: Finished DIAMOND blastx chunk 38 of 118
2026/09/28 04:11:02 PM INFO: Finished DIAMOND blastx chunk 39 of 118
2026/09/28 04:11:02 PM INFO: Finished DIAMOND blastx chunk 40 of 118
2026/09/28 04:11:03 PM INFO: Finished DIAMOND blastx chunk 41 of 118
2026/09/28 04:11:04 PM INFO: Finished DIAMOND blastx chunk 42 of 118
2026/09/28 04:11:05 PM INFO: Finished DIAMOND blastx chunk 43 of 118
2026/09/28 04:11:05 PM INFO: Finished DIAMOND blastx chunk 44 of 118
2026/09/28 04:11:07 PM INFO: Finished DIAMOND blastx chunk 45 of 118
2026/09/28 04:11:07 PM INFO: Finished DIAMOND blastx chunk 46 of 118
2026/09/28 04:11:07 PM INFO: Finished DIAMOND blastx chunk 47 of 118
2026/09/28 04:11:07 PM INFO: Finished DIAMOND blastx chunk 48 of 118
2026/09/28 04:11:07 PM INFO: Finished DIAMOND blastx chunk 49 of 118
2026/09/28 04:11:08 PM INFO: Finished DIAMOND blastx chunk 50 of 118
2026/09/28 04:11:08 PM INFO: Finished DIAMOND blastx chunk 51 of 118
2026/09/28 04:11:08 PM INFO: Finished DIAMOND blastx chunk 52 of 118
2026/09/28 04:11:08 PM INFO: Finished DIAMOND blastx chunk 53 of 118
2026/09/28 04:11:08 PM INFO: Finished DIAMOND blastx chunk 54 of 118
2026/09/28 04:11:08 PM INFO: Finished DIAMOND blastx chunk 55 of 118
2026/09/28 04:11:09 PM INFO: Finished DIAMOND blastx chunk 56 of 118
2026/09/28 04:11:09 PM INFO: Finished DIAMOND blastx chunk 57 of 118
2026/09/28 04:11:09 PM INFO: Finished DIAMOND blastx chunk 58 of 118
2026/09/28 04:11:09 PM INFO: Finished DIAMOND blastx chunk 59 of 118
2026/09/28 04:11:09 PM INFO: Finished DIAMOND blastx chunk 60 of 118
2026/09/28 04:11:09 PM INFO: Finished DIAMOND blastx chunk 61 of 118
2026/09/28 04:11:10 PM INFO: Finished DIAMOND blastx chunk 62 of 118
2026/09/28 04:11:10 PM INFO: Finished DIAMOND blastx chunk 63 of 118
2026/09/28 04:11:10 PM INFO: Finished DIAMOND blastx chunk 64 of 118
2026/09/28 04:11:10 PM INFO: Finished DIAMOND blastx chunk 65 of 118
2026/09/28 04:11:10 PM INFO: Finished DIAMOND blastx chunk 66 of 118
2026/09/28 04:11:11 PM INFO: Finished DIAMOND blastx chunk 67 of 118
2026/09/28 04:11:11 PM INFO: Finished DIAMOND blastx chunk 68 of 118
2026/09/28 04:11:11 PM INFO: Finished DIAMOND blastx chunk 69 of 118
2026/09/28 04:11:11 PM INFO: Finished DIAMOND blastx chunk 70 of 118
2026/09/28 04:11:11 PM INFO: Finished DIAMOND blastx chunk 71 of 118
2026/09/28 04:11:12 PM INFO: Finished DIAMOND blastx chunk 72 of 118
2026/09/28 04:11:13 PM INFO: Finished DIAMOND blastx chunk 73 of 118
2026/09/28 04:11:14 PM INFO: Finished DIAMOND blastx chunk 74 of 118
2026/09/28 04:14:44 PM INFO: Finished DIAMOND blastx chunk 75 of 118
2026/09/28 04:14:50 PM INFO: Finished DIAMOND blastx chunk 76 of 118
2026/09/28 04:16:34 PM INFO: Finished DIAMOND blastx chunk 77 of 118
2026/09/28 04:16:47 PM INFO: Finished DIAMOND blastx chunk 78 of 118
2026/09/28 04:18:49 PM INFO: Finished DIAMOND blastx chunk 79 of 118
2026/09/28 04:18:51 PM INFO: Finished DIAMOND blastx chunk 80 of 118
2026/09/28 04:21:25 PM INFO: Finished DIAMOND blastx chunk 81 of 118
2026/09/28 04:21:35 PM INFO: Finished DIAMOND blastx chunk 82 of 118
2026/09/28 04:22:25 PM INFO: Finished DIAMOND blastx chunk 83 of 118
2026/09/28 04:22:28 PM INFO: Finished DIAMOND blastx chunk 84 of 118
2026/09/28 04:24:39 PM INFO: Finished DIAMOND blastx chunk 85 of 118
2026/09/28 04:25:28 PM INFO: Finished DIAMOND blastx chunk 86 of 118
2026/09/28 04:25:37 PM INFO: Finished DIAMOND blastx chunk 87 of 118
2026/09/28 04:26:29 PM INFO: Finished DIAMOND blastx chunk 88 of 118
2026/09/28 04:28:16 PM INFO: Finished DIAMOND blastx chunk 89 of 118
2026/09/28 04:28:35 PM INFO: Finished DIAMOND blastx chunk 90 of 118
2026/09/28 04:29:03 PM INFO: Finished DIAMOND blastx chunk 91 of 118
2026/09/28 04:29:25 PM INFO: Finished DIAMOND blastx chunk 92 of 118
2026/09/28 04:30:58 PM INFO: Finished DIAMOND blastx chunk 93 of 118
2026/09/28 04:31:07 PM INFO: Finished DIAMOND blastx chunk 94 of 118
2026/09/28 04:32:27 PM INFO: Finished DIAMOND blastx chunk 95 of 118
2026/09/28 04:32:48 PM INFO: Finished DIAMOND blastx chunk 96 of 118
2026/09/28 04:33:30 PM INFO: Finished DIAMOND blastx chunk 97 of 118
2026/09/28 04:33:51 PM INFO: Finished DIAMOND blastx chunk 98 of 118
2026/09/28 04:35:13 PM INFO: Finished DIAMOND blastx chunk 99 of 118
2026/09/28 04:35:39 PM INFO: Finished DIAMOND blastx chunk 100 of 118
2026/09/28 04:35:45 PM INFO: Finished DIAMOND blastx chunk 101 of 118
2026/09/28 04:36:12 PM INFO: Finished DIAMOND blastx chunk 102 of 118
2026/09/28 04:38:02 PM INFO: Finished DIAMOND blastx chunk 103 of 118
2026/09/28 04:38:41 PM INFO: Finished DIAMOND blastx chunk 104 of 118
2026/09/28 04:39:05 PM INFO: Finished DIAMOND blastx chunk 105 of 118
2026/09/28 04:39:40 PM INFO: Finished DIAMOND blastx chunk 106 of 118
2026/09/28 04:43:38 PM INFO: Finished DIAMOND blastx chunk 107 of 118
2026/09/28 04:44:26 PM INFO: Finished DIAMOND blastx chunk 108 of 118
2026/09/28 04:46:37 PM INFO: Finished DIAMOND blastx chunk 109 of 118
2026/09/28 04:47:10 PM INFO: Finished DIAMOND blastx chunk 110 of 118
2026/09/28 04:51:06 PM INFO: Finished DIAMOND blastx chunk 111 of 118
2026/09/28 04:52:18 PM INFO: Finished DIAMOND blastx chunk 112 of 118
2026/09/28 04:53:02 PM INFO: Finished DIAMOND blastx chunk 113 of 118
2026/09/28 04:53:30 PM INFO: Finished DIAMOND blastx chunk 114 of 118
2026/09/28 04:54:26 PM INFO: Finished DIAMOND blastx chunk 115 of 118
2026/09/28 04:55:07 PM INFO: Finished DIAMOND blastx chunk 116 of 118
2026/09/28 04:56:49 PM INFO: Finished DIAMOND blastx chunk 117 of 118
2026/09/28 04:57:28 PM INFO: Finished DIAMOND blastx chunk 118 of 118
2026/09/28 04:57:28 PM INFO: Finished running taxonomic assignment
2026/09/28 04:58:02 PM INFO: Finished
2026/09/28 04:58:02 PM INFO: Writing /labgenomaarea2/valentina.pavelecini/EHI/SingleM/RF-5-B_otu_table.tsv
2026/09/28 04:58:02 PM INFO: Writing taxonomic profile to /labgenomaarea2/valentina.pavelecini/EHI/SingleM/RF-5-B_taxonomic_profile.tsv
2026/09/28 04:58:02 PM INFO: Using minimum taxon coverage of 0.35
2026/09/28 04:58:02 PM INFO: Removing off-target OTUs from TF-2587-RF-5-B_S4_L001_R1
2026/09/28 04:58:02 PM INFO: Found 35791.95 assigned and 0.00 unassigned OTU coverage units
2026/09/28 04:58:02 PM INFO: After removing off-target OTUs, found 30091.35 assigned and 0.00 unassigned OTU coverage units
2026/09/28 04:58:02 PM INFO: Total OTU coverage by query: 11334.89855749553
2026/09/28 04:58:02 PM INFO: Total OTU coverage by diamond: 18756.45272374535
2026/09/28 04:58:02 PM INFO: Applying species-wise expectation maximization algorithm to OTU table
2026/09/28 04:58:02 PM INFO: Found 6 species uniquely hitting >= 10 marker genes
2026/09/28 04:58:06 PM INFO: Species-wise EM converged in 96 steps
2026/09/28 04:58:06 PM INFO: Gathering equivalence classes
2026/09/28 04:58:06 PM INFO: Demultiplexing OTU table
2026/09/28 04:58:06 PM INFO: Finished expectation maximization
2026/09/28 04:58:06 PM INFO: Converting DIAMOND IDs to taxons
2026/09/28 04:58:34 PM INFO: Converted 4603 Diamond-assigned OTU taxon_ids to taxon strings
2026/09/28 04:58:34 PM INFO: Applying genus-wise expectation maximization algorithm to OTU table
2026/09/28 04:58:38 PM INFO: Genus-wise EM converged in 79 steps
2026/09/28 04:58:38 PM INFO: Gathering equivalence classes
2026/09/28 04:58:38 PM INFO: Demultiplexing OTU table
2026/09/28 04:58:38 PM INFO: Finished genus expectation maximization
2026/09/28 04:58:39 PM INFO: Total profile coverage after condense domain to species: 844.656681119038
2026/09/28 04:58:39 PM INFO: Total profile coverage after push down: 844.6566811190384
2026/09/28 04:58:39 PM INFO: Taxonomic level coverage:
2026/09/28 04:58:39 PM INFO: kingdom:   3.68%   2 taxons
2026/09/28 04:58:39 PM INFO: phylum:    2.37%   20 taxons
2026/09/28 04:58:39 PM INFO: class:     5.40%   31 taxons
2026/09/28 04:58:39 PM INFO: order:     7.14%   47 taxons
2026/09/28 04:58:39 PM INFO: family:    11.16%  60 taxons
2026/09/28 04:58:39 PM INFO: genus:     46.45%  66 taxons
2026/09/28 04:58:39 PM INFO: species:   23.78%  27 taxons
2026/09/28 04:58:39 PM INFO: Finished condense
```

### RF-2-B
```
singlem pipe \
    -1 /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-RF-2-B_S2_L001_R1.fastq.gz \
    -2 /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-RF-2-B_S2_L001_R2.fastq.gz \
    --otu-table /labgenomaarea2/valentina.pavelecini/EHI/SingleM/RF-2-B_otu_table.tsv \
    --taxonomic-profile /labgenomaarea2/valentina.pavelecini/EHI/SingleM/RF-2-B_taxonomic_profile.tsv \
    --threads 2
```

Output:
```
2026/09/28 08:41:31 AM INFO: SingleM v0.21.4
2026/09/28 08:41:31 AM INFO: Retrieval successful. Location of backpack is: /labgenomaarea2/valentina.pavelecini/EHI/SingleM/metapackage/S6.5.0.GTDB_r232.metapackage_20260319.smpkg.zb
2026/09/28 08:41:31 AM INFO: Loaded 59 SingleM packages
2026/09/28 08:41:37 AM INFO: Using as input 1 different pairs of sequence files e.g. /labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-RF-2-B_S2_L001_R1.fastq.gz & /labgenoma>
2026/09/28 08:41:37 AM INFO: DIAMOND version: diamond version 2.2.6
2026/09/28 08:41:37 AM INFO: Filtering sequence files through DIAMOND blastx
2026/09/28 08:41:37 AM INFO: Filtering TF-2587-RF-2-B_S2_L001_R1.fastq.gz
^M                                                                                ^M2026/09/28 09:19:59 AM INFO: Found 150227 hits for TF-2587-RF-2-B_S2_L001_R1.fastq.gz
2026/09/28 09:19:59 AM INFO: Filtering TF-2587-RF-2-B_S2_L001_R2.fastq.gz
^M                                                                                ^M2026/09/28 09:58:28 AM INFO: Found 150981 hits for TF-2587-RF-2-B_S2_L001_R2.fastq.gz
2026/09/28 09:58:28 AM INFO: Finished DIAMOND prefilter phase
2026/09/28 09:58:28 AM INFO: Assigning sequences to SingleM packages with DIAMOND ..
2026/09/28 09:58:28 AM INFO: Extracting reads from 1 sample(s) across 59 package(s) using 2 thread(s)
2026/09/28 10:00:33 AM INFO: Finished extracting reads for chunk 1 of 59
2026/09/28 10:00:34 AM INFO: Finished extracting reads for chunk 2 of 59
2026/09/28 10:00:38 AM INFO: Finished extracting reads for chunk 3 of 59
2026/09/28 10:00:40 AM INFO: Finished extracting reads for chunk 4 of 59
2026/09/28 10:00:44 AM INFO: Finished extracting reads for chunk 5 of 59
2026/09/28 10:00:45 AM INFO: Finished extracting reads for chunk 6 of 59
2026/09/28 10:00:49 AM INFO: Finished extracting reads for chunk 7 of 59
2026/09/28 10:00:51 AM INFO: Finished extracting reads for chunk 8 of 59
2026/09/28 10:00:54 AM INFO: Finished extracting reads for chunk 9 of 59
2026/09/28 10:00:57 AM INFO: Finished extracting reads for chunk 10 of 59
2026/09/28 10:01:00 AM INFO: Finished extracting reads for chunk 11 of 59
2026/09/28 10:01:03 AM INFO: Finished extracting reads for chunk 12 of 59
2026/09/28 10:01:05 AM INFO: Finished extracting reads for chunk 13 of 59
2026/09/28 10:01:08 AM INFO: Finished extracting reads for chunk 14 of 59
2026/09/28 10:01:10 AM INFO: Finished extracting reads for chunk 15 of 59
2026/09/28 10:01:12 AM INFO: Finished extracting reads for chunk 16 of 59
2026/09/28 10:01:12 AM INFO: Finished extracting reads for chunk 17 of 59
2026/09/28 10:01:15 AM INFO: Finished extracting reads for chunk 18 of 59
2026/09/28 10:01:15 AM INFO: Finished extracting reads for chunk 19 of 59
2026/09/28 10:01:16 AM INFO: Finished extracting reads for chunk 20 of 59
2026/09/28 10:01:17 AM INFO: Finished extracting reads for chunk 21 of 59
2026/09/28 10:01:18 AM INFO: Finished extracting reads for chunk 22 of 59
2026/09/28 10:01:19 AM INFO: Finished extracting reads for chunk 23 of 59
2026/09/28 10:01:20 AM INFO: Finished extracting reads for chunk 24 of 59
2026/09/28 10:01:21 AM INFO: Finished extracting reads for chunk 25 of 59
2026/09/28 10:01:21 AM INFO: Finished extracting reads for chunk 26 of 59
2026/09/28 10:01:22 AM INFO: Finished extracting reads for chunk 27 of 59
2026/09/28 10:01:23 AM INFO: Finished extracting reads for chunk 28 of 59
2026/09/28 10:01:24 AM INFO: Finished extracting reads for chunk 29 of 59
2026/09/28 10:01:25 AM INFO: Finished extracting reads for chunk 30 of 59
2026/09/28 10:01:26 AM INFO: Finished extracting reads for chunk 31 of 59
2026/09/28 10:01:27 AM INFO: Finished extracting reads for chunk 32 of 59
2026/09/28 10:01:28 AM INFO: Finished extracting reads for chunk 33 of 59
2026/09/28 10:01:28 AM INFO: Finished extracting reads for chunk 34 of 59
2026/09/28 10:01:30 AM INFO: Finished extracting reads for chunk 35 of 59
2026/09/28 10:01:30 AM INFO: Finished extracting reads for chunk 36 of 59
2026/09/28 10:01:33 AM INFO: Finished extracting reads for chunk 37 of 59
2026/09/28 10:01:37 AM INFO: Finished extracting reads for chunk 38 of 59
2026/09/28 10:01:40 AM INFO: Finished extracting reads for chunk 39 of 59
2026/09/28 10:01:41 AM INFO: Finished extracting reads for chunk 40 of 59
2026/09/28 10:01:46 AM INFO: Finished extracting reads for chunk 41 of 59
2026/09/28 10:01:46 AM INFO: Finished extracting reads for chunk 42 of 59
2026/09/28 10:01:50 AM INFO: Finished extracting reads for chunk 43 of 59
2026/09/28 10:01:51 AM INFO: Finished extracting reads for chunk 44 of 59
2026/09/28 10:01:56 AM INFO: Finished extracting reads for chunk 45 of 59
2026/09/28 10:01:56 AM INFO: Finished extracting reads for chunk 46 of 59
2026/09/28 10:02:00 AM INFO: Finished extracting reads for chunk 47 of 59
2026/09/28 10:02:02 AM INFO: Finished extracting reads for chunk 48 of 59
2026/09/28 10:02:05 AM INFO: Finished extracting reads for chunk 49 of 59
2026/09/28 10:02:07 AM INFO: Finished extracting reads for chunk 50 of 59
2026/09/28 10:02:09 AM INFO: Finished extracting reads for chunk 51 of 59
2026/09/28 10:02:12 AM INFO: Finished extracting reads for chunk 52 of 59
2026/09/28 10:02:14 AM INFO: Finished extracting reads for chunk 53 of 59
2026/09/28 10:02:17 AM INFO: Finished extracting reads for chunk 54 of 59
2026/09/28 10:02:20 AM INFO: Finished extracting reads for chunk 55 of 59
2026/09/28 10:02:23 AM INFO: Finished extracting reads for chunk 56 of 59
2026/09/28 10:02:24 AM INFO: Finished extracting reads for chunk 57 of 59
2026/09/28 10:02:28 AM INFO: Finished extracting reads for chunk 58 of 59
2026/09/28 10:02:28 AM INFO: Finished extracting reads for chunk 59 of 59
2026/09/28 10:02:28 AM INFO: After read extraction, 139297 sequence(s) remain
2026/09/28 10:02:28 AM INFO: Running taxonomic assignment ..
2026/09/28 10:02:28 AM INFO: Assigning taxonomy by singlem query ..
2026/09/28 10:02:55 AM INFO: Progress: 330/12696 (2%)
2026/09/28 10:03:04 AM INFO: Progress: 640/12696 (5%)
2026/09/28 10:03:13 AM INFO: Progress: 923/12696 (7%)
2026/09/28 10:03:22 AM INFO: Progress: 1267/12696 (9%)
2026/09/28 10:03:31 AM INFO: Progress: 1581/12696 (12%)
2026/09/28 10:03:40 AM INFO: Progress: 1885/12696 (14%)
2026/09/28 10:03:50 AM INFO: Progress: 2192/12696 (17%)
2026/09/28 10:03:59 AM INFO: Progress: 2536/12696 (19%)
2026/09/28 10:04:08 AM INFO: Progress: 2815/12696 (22%)
2026/09/28 10:04:17 AM INFO: Progress: 3133/12696 (24%)
2026/09/28 10:04:26 AM INFO: Progress: 3463/12696 (27%)
2026/09/28 10:04:35 AM INFO: Progress: 3709/12696 (29%)
2026/09/28 10:04:44 AM INFO: Progress: 4015/12696 (31%)
2026/09/28 10:04:53 AM INFO: Progress: 4331/12696 (34%)
2026/09/28 10:05:03 AM INFO: Progress: 4571/12696 (36%)
2026/09/28 10:05:13 AM INFO: Progress: 4867/12696 (38%)
2026/09/28 10:05:25 AM INFO: Progress: 5451/12696 (42%)
2026/09/28 10:05:41 AM INFO: Progress: 5577/12696 (43%)
2026/09/28 10:05:51 AM INFO: Progress: 5656/12696 (44%)
2026/09/28 10:05:58 AM INFO: Progress: 5723/12696 (45%)
2026/09/28 10:06:05 AM INFO: Progress: 5804/12696 (45%)
2026/09/28 10:06:13 AM INFO: Progress: 5838/12696 (45%)
2026/09/28 10:06:19 AM INFO: Progress: 5900/12696 (46%)
2026/09/28 10:06:28 AM INFO: Progress: 6247/12696 (49%)
2026/09/28 10:06:37 AM INFO: Progress: 6558/12696 (51%)
2026/09/28 10:06:46 AM INFO: Progress: 6881/12696 (54%)
2026/09/28 10:06:56 AM INFO: Progress: 7204/12696 (56%)
2026/09/28 10:07:05 AM INFO: Progress: 7536/12696 (59%)
2026/09/28 10:07:14 AM INFO: Progress: 7820/12696 (61%)
2026/09/28 10:07:23 AM INFO: Progress: 8018/12696 (63%)
2026/09/28 10:07:33 AM INFO: Progress: 8330/12696 (65%)
2026/09/28 10:07:42 AM INFO: Progress: 8649/12696 (68%)
2026/09/28 10:07:51 AM INFO: Progress: 8890/12696 (70%)
2026/09/28 10:08:01 AM INFO: Progress: 9153/12696 (72%)
2026/09/28 10:08:10 AM INFO: Progress: 9466/12696 (74%)
2026/09/28 10:08:19 AM INFO: Progress: 9749/12696 (76%)
2026/09/28 10:08:28 AM INFO: Progress: 10044/12696 (79%)
2026/09/28 10:08:38 AM INFO: Progress: 10303/12696 (81%)
2026/09/28 10:08:47 AM INFO: Progress: 10598/12696 (83%)
2026/09/28 10:08:57 AM INFO: Progress: 10917/12696 (85%)
2026/09/28 10:09:06 AM INFO: Progress: 11251/12696 (88%)
2026/09/28 10:09:15 AM INFO: Progress: 11579/12696 (91%)
2026/09/28 10:09:24 AM INFO: Progress: 11846/12696 (93%)
2026/09/28 10:09:33 AM INFO: Progress: 12096/12696 (95%)
2026/09/28 10:09:42 AM INFO: Progress: 12387/12696 (97%)
2026/09/28 10:09:52 AM INFO: Progress: 12696/12696 (100%)
2026/09/28 10:09:52 AM INFO: Finished running singlem query-based taxonomic assignment, now running diamond using 2 thread(s) ..
2026/09/28 10:09:52 AM INFO: Assigning taxonomy with DIAMOND blastx to 42014 OTUs (The 81.6% that were not assigned by smafa) ..
2026/09/28 10:16:56 AM INFO: Finished DIAMOND blastx chunk 1 of 118
2026/09/28 10:17:06 AM INFO: Finished DIAMOND blastx chunk 2 of 118
2026/09/28 10:20:12 AM INFO: Finished DIAMOND blastx chunk 3 of 118
2026/09/28 10:20:18 AM INFO: Finished DIAMOND blastx chunk 4 of 118
2026/09/28 10:27:23 AM INFO: Finished DIAMOND blastx chunk 5 of 118
2026/09/28 10:27:37 AM INFO: Finished DIAMOND blastx chunk 6 of 118
2026/09/28 10:34:05 AM INFO: Finished DIAMOND blastx chunk 7 of 118
2026/09/28 10:34:35 AM INFO: Finished DIAMOND blastx chunk 8 of 118
2026/09/28 10:39:15 AM INFO: Finished DIAMOND blastx chunk 9 of 118
2026/09/28 10:40:43 AM INFO: Finished DIAMOND blastx chunk 10 of 118
2026/09/28 10:45:57 AM INFO: Finished DIAMOND blastx chunk 11 of 118
2026/09/28 10:47:45 AM INFO: Finished DIAMOND blastx chunk 12 of 118
2026/09/28 10:52:14 AM INFO: Finished DIAMOND blastx chunk 13 of 118
2026/09/28 10:53:59 AM INFO: Finished DIAMOND blastx chunk 14 of 118
2026/09/28 11:03:09 AM INFO: Finished DIAMOND blastx chunk 15 of 118
2026/09/28 11:07:23 AM INFO: Finished DIAMOND blastx chunk 16 of 118
2026/09/28 11:07:56 AM INFO: Finished DIAMOND blastx chunk 17 of 118
2026/09/28 11:13:07 AM INFO: Finished DIAMOND blastx chunk 18 of 118
2026/09/28 11:15:38 AM INFO: Finished DIAMOND blastx chunk 19 of 118
2026/09/28 11:19:14 AM INFO: Finished DIAMOND blastx chunk 20 of 118
2026/09/28 11:21:26 AM INFO: Finished DIAMOND blastx chunk 21 of 118
2026/09/28 11:24:25 AM INFO: Finished DIAMOND blastx chunk 22 of 118
2026/09/28 11:24:33 AM INFO: Finished DIAMOND blastx chunk 23 of 118
2026/09/28 11:28:44 AM INFO: Finished DIAMOND blastx chunk 24 of 118
2026/09/28 11:28:55 AM INFO: Finished DIAMOND blastx chunk 25 of 118
2026/09/28 11:32:27 AM INFO: Finished DIAMOND blastx chunk 26 of 118
2026/09/28 11:33:39 AM INFO: Finished DIAMOND blastx chunk 27 of 118
2026/09/28 11:35:55 AM INFO: Finished DIAMOND blastx chunk 28 of 118
2026/09/28 11:39:56 AM INFO: Finished DIAMOND blastx chunk 29 of 118
2026/09/28 11:40:49 AM INFO: Finished DIAMOND blastx chunk 30 of 118
2026/09/28 11:45:30 AM INFO: Finished DIAMOND blastx chunk 31 of 118
2026/09/28 11:45:31 AM INFO: Finished DIAMOND blastx chunk 32 of 118
2026/09/28 11:45:31 AM INFO: Finished DIAMOND blastx chunk 33 of 118
2026/09/28 11:45:31 AM INFO: Finished DIAMOND blastx chunk 34 of 118
2026/09/28 11:45:32 AM INFO: Finished DIAMOND blastx chunk 35 of 118
2026/09/28 11:45:32 AM INFO: Finished DIAMOND blastx chunk 36 of 118
2026/09/28 11:45:33 AM INFO: Finished DIAMOND blastx chunk 37 of 118
2026/09/28 11:45:33 AM INFO: Finished DIAMOND blastx chunk 38 of 118
2026/09/28 11:45:34 AM INFO: Finished DIAMOND blastx chunk 39 of 118
2026/09/28 11:45:34 AM INFO: Finished DIAMOND blastx chunk 40 of 118
2026/09/28 11:45:34 AM INFO: Finished DIAMOND blastx chunk 41 of 118
2026/09/28 11:45:34 AM INFO: Finished DIAMOND blastx chunk 42 of 118
2026/09/28 11:45:36 AM INFO: Finished DIAMOND blastx chunk 43 of 118
2026/09/28 11:45:36 AM INFO: Finished DIAMOND blastx chunk 44 of 118
2026/09/28 11:45:37 AM INFO: Finished DIAMOND blastx chunk 45 of 118
2026/09/28 11:45:37 AM INFO: Finished DIAMOND blastx chunk 46 of 118
2026/09/28 11:45:38 AM INFO: Finished DIAMOND blastx chunk 47 of 118
2026/09/28 11:45:38 AM INFO: Finished DIAMOND blastx chunk 48 of 118
2026/09/28 11:45:39 AM INFO: Finished DIAMOND blastx chunk 49 of 118
2026/09/28 11:45:39 AM INFO: Finished DIAMOND blastx chunk 50 of 118
2026/09/28 11:45:39 AM INFO: Finished DIAMOND blastx chunk 51 of 118
2026/09/28 11:45:39 AM INFO: Finished DIAMOND blastx chunk 52 of 118
2026/09/28 11:45:40 AM INFO: Finished DIAMOND blastx chunk 53 of 118
2026/09/28 11:45:40 AM INFO: Finished DIAMOND blastx chunk 54 of 118
2026/09/28 11:45:40 AM INFO: Finished DIAMOND blastx chunk 55 of 118
2026/09/28 11:45:40 AM INFO: Finished DIAMOND blastx chunk 56 of 118
2026/09/28 11:45:40 AM INFO: Finished DIAMOND blastx chunk 57 of 118
2026/09/28 11:45:40 AM INFO: Finished DIAMOND blastx chunk 58 of 118
2026/09/28 11:45:41 AM INFO: Finished DIAMOND blastx chunk 59 of 118
2026/09/28 11:45:41 AM INFO: Finished DIAMOND blastx chunk 60 of 118
2026/09/28 11:45:41 AM INFO: Finished DIAMOND blastx chunk 61 of 118
2026/09/28 11:45:41 AM INFO: Finished DIAMOND blastx chunk 62 of 118
2026/09/28 11:45:41 AM INFO: Finished DIAMOND blastx chunk 63 of 118
2026/09/28 11:45:41 AM INFO: Finished DIAMOND blastx chunk 64 of 118
2026/09/28 11:45:42 AM INFO: Finished DIAMOND blastx chunk 65 of 118
2026/09/28 11:45:42 AM INFO: Finished DIAMOND blastx chunk 66 of 118
2026/09/28 11:45:43 AM INFO: Finished DIAMOND blastx chunk 67 of 118
2026/09/28 11:45:43 AM INFO: Finished DIAMOND blastx chunk 68 of 118
2026/09/28 11:45:43 AM INFO: Finished DIAMOND blastx chunk 69 of 118
2026/09/28 11:45:43 AM INFO: Finished DIAMOND blastx chunk 70 of 118
2026/09/28 11:45:44 AM INFO: Finished DIAMOND blastx chunk 71 of 118
2026/09/28 11:45:44 AM INFO: Finished DIAMOND blastx chunk 72 of 118
2026/09/28 11:45:48 AM INFO: Finished DIAMOND blastx chunk 73 of 118
2026/09/28 11:45:48 AM INFO: Finished DIAMOND blastx chunk 74 of 118
2026/09/28 11:54:05 AM INFO: Finished DIAMOND blastx chunk 75 of 118
2026/09/28 11:55:24 AM INFO: Finished DIAMOND blastx chunk 76 of 118
2026/09/28 12:00:05 PM INFO: Finished DIAMOND blastx chunk 77 of 118
2026/09/28 12:00:36 PM INFO: Finished DIAMOND blastx chunk 78 of 118
2026/09/28 12:06:36 PM INFO: Finished DIAMOND blastx chunk 79 of 118
2026/09/28 12:07:04 PM INFO: Finished DIAMOND blastx chunk 80 of 118
2026/09/28 12:12:45 PM INFO: Finished DIAMOND blastx chunk 81 of 118
2026/09/28 12:13:47 PM INFO: Finished DIAMOND blastx chunk 82 of 118
2026/09/28 12:16:41 PM INFO: Finished DIAMOND blastx chunk 83 of 118
2026/09/28 12:17:12 PM INFO: Finished DIAMOND blastx chunk 84 of 118
2026/09/28 12:20:28 PM INFO: Finished DIAMOND blastx chunk 85 of 118
2026/09/28 12:21:10 PM INFO: Finished DIAMOND blastx chunk 86 of 118
2026/09/28 12:23:34 PM INFO: Finished DIAMOND blastx chunk 87 of 118
2026/09/28 12:24:01 PM INFO: Finished DIAMOND blastx chunk 88 of 118
2026/09/28 12:29:48 PM INFO: Finished DIAMOND blastx chunk 89 of 118
2026/09/28 12:30:23 PM INFO: Finished DIAMOND blastx chunk 90 of 118
2026/09/28 12:31:38 PM INFO: Finished DIAMOND blastx chunk 91 of 118
2026/09/28 12:32:03 PM INFO: Finished DIAMOND blastx chunk 92 of 118
2026/09/28 12:34:50 PM INFO: Finished DIAMOND blastx chunk 93 of 118
2026/09/28 12:35:34 PM INFO: Finished DIAMOND blastx chunk 94 of 118
2026/09/28 12:39:30 PM INFO: Finished DIAMOND blastx chunk 95 of 118
2026/09/28 12:40:24 PM INFO: Finished DIAMOND blastx chunk 96 of 118
2026/09/28 12:42:26 PM INFO: Finished DIAMOND blastx chunk 97 of 118
2026/09/28 12:43:01 PM INFO: Finished DIAMOND blastx chunk 98 of 118
2026/09/28 12:46:15 PM INFO: Finished DIAMOND blastx chunk 99 of 118
2026/09/28 12:47:20 PM INFO: Finished DIAMOND blastx chunk 100 of 118
2026/09/28 12:47:31 PM INFO: Finished DIAMOND blastx chunk 101 of 118
2026/09/28 12:48:40 PM INFO: Finished DIAMOND blastx chunk 102 of 118
2026/09/28 12:51:58 PM INFO: Finished DIAMOND blastx chunk 103 of 118
2026/09/28 12:52:52 PM INFO: Finished DIAMOND blastx chunk 104 of 118
2026/09/28 12:54:30 PM INFO: Finished DIAMOND blastx chunk 105 of 118
2026/09/28 12:55:04 PM INFO: Finished DIAMOND blastx chunk 106 of 118
2026/09/28 01:03:15 PM INFO: Finished DIAMOND blastx chunk 107 of 118
2026/09/28 01:04:17 PM INFO: Finished DIAMOND blastx chunk 108 of 118
2026/09/28 01:09:50 PM INFO: Finished DIAMOND blastx chunk 109 of 118
2026/09/28 01:10:04 PM INFO: Finished DIAMOND blastx chunk 110 of 118
2026/09/28 01:20:12 PM INFO: Finished DIAMOND blastx chunk 111 of 118
2026/09/28 01:20:56 PM INFO: Finished DIAMOND blastx chunk 112 of 118
2026/09/28 01:22:34 PM INFO: Finished DIAMOND blastx chunk 113 of 118
2026/09/28 01:23:14 PM INFO: Finished DIAMOND blastx chunk 114 of 118
2026/09/28 01:25:45 PM INFO: Finished DIAMOND blastx chunk 115 of 118
2026/09/28 01:25:49 PM INFO: Finished DIAMOND blastx chunk 116 of 118
2026/09/28 01:32:35 PM INFO: Finished DIAMOND blastx chunk 117 of 118
2026/09/28 01:33:23 PM INFO: Finished DIAMOND blastx chunk 118 of 118
2026/09/28 01:33:23 PM INFO: Finished running taxonomic assignment
2026/09/28 01:34:20 PM INFO: Finished
2026/09/28 01:34:20 PM INFO: Writing /labgenomaarea2/valentina.pavelecini/EHI/SingleM/RF-2-B_otu_table.tsv
2026/09/28 01:34:20 PM INFO: Writing taxonomic profile to /labgenomaarea2/valentina.pavelecini/EHI/SingleM/RF-2-B_taxonomic_profile.tsv
2026/09/28 01:34:20 PM INFO: Using minimum taxon coverage of 0.35
2026/09/28 01:34:20 PM INFO: Removing off-target OTUs from TF-2587-RF-2-B_S2_L001_R1
2026/09/28 01:34:20 PM INFO: Found 69521.07 assigned and 0.00 unassigned OTU coverage units
2026/09/28 01:34:20 PM INFO: After removing off-target OTUs, found 59802.04 assigned and 0.00 unassigned OTU coverage units
2026/09/28 01:34:20 PM INFO: Total OTU coverage by query: 11717.857394122675
2026/09/28 01:34:20 PM INFO: Total OTU coverage by diamond: 48084.18592057608
2026/09/28 01:34:20 PM INFO: Applying species-wise expectation maximization algorithm to OTU table
2026/09/28 01:34:20 PM INFO: Found 4 species uniquely hitting >= 10 marker genes
2026/09/28 01:34:23 PM INFO: Species-wise EM converged in 57 steps
2026/09/28 01:34:23 PM INFO: Gathering equivalence classes
2026/09/28 01:34:23 PM INFO: Demultiplexing OTU table
2026/09/28 01:34:23 PM INFO: Finished expectation maximization
2026/09/28 01:34:23 PM INFO: Converting DIAMOND IDs to taxons
2026/09/28 01:35:10 PM INFO: Converted 9156 Diamond-assigned OTU taxon_ids to taxon strings
2026/09/28 01:35:10 PM INFO: Applying genus-wise expectation maximization algorithm to OTU table
2026/09/28 01:35:24 PM INFO: Genus-wise EM converged in 217 steps
2026/09/28 01:35:24 PM INFO: Gathering equivalence classes
2026/09/28 01:35:24 PM INFO: Demultiplexing OTU table
2026/09/28 01:35:24 PM INFO: Finished genus expectation maximization
2026/09/28 01:35:24 PM INFO: Total profile coverage after condense domain to species: 1663.716928486934
2026/09/28 01:35:24 PM INFO: Total profile coverage after push down: 1663.7169284869337
2026/09/28 01:35:24 PM INFO: Taxonomic level coverage:
2026/09/28 01:35:24 PM INFO: kingdom:   2.33%   2 taxons
2026/09/28 01:35:24 PM INFO: phylum:    2.86%   22 taxons
2026/09/28 01:35:24 PM INFO: class:     7.56%   46 taxons
2026/09/28 01:35:24 PM INFO: order:     8.74%   79 taxons
2026/09/28 01:35:24 PM INFO: family:    15.66%  133 taxons
2026/09/28 01:35:24 PM INFO: genus:     54.84%  126 taxons
2026/09/28 01:35:24 PM INFO: species:   8.01%   15 taxons
2026/09/28 01:35:24 PM INFO: Finished condense
```

Agora, terminamos o pré-processamento dos dados e iremos prosseguir para a terceira parte: Assembly and Binning.

# Agrupamento e categorização: Agrupamento metagenômico

```
megahit \
    -t {threads} \
    --verbose \
    --min-contig-len 1500 \
    -1 {input.r1} -2 {input.r2} \
    -f \
    -o {config[workdir]}/{wildcards.PRB}_{wildcards.EHI}_assembly
```

Neste passo, iremos reconstruir o genoma dos microorganismos presentes nas amostras. O material genético de todos os organismos presentes nessa amostra serão alvo do sequenciamento, para que a comunidade microbiana seja estudada por completo simultâneamente.

Estaremos usando o Megahit, devido a ser o agrupador padrão da pipeline do EHI, sua simplicidade e aos requisitos de memória altos demais do Metaspades. Instalamos ele com ```conda install -c bioconda megahit```

Temos duas opções: rodar todas as amostras no Megahit, uma de cada vez individualmente (individual assembly), ou agrupá-las e então rodar no Megahit (coassembly). Como não fizemos a etapa de remover genoma hospedeiro, que é necessário para o coassembly, iremos em compensação fazer com o reads filtrados diretamente. Se não der certo, vamos fazer individual assembly.

Primeiro, abrimos uma janela do tmux com ```tmux new -s megahit``` e definimos os parâmetros da sessão para ```srun -N 1 -n 1 -c 24 -t 24:00:00 --pty bash -i```, em que:
- ```-N 1``` = 1 máquina
- ```-n 1``` = 1 tarefa
- ```-c 24``` 24 CPUs/cores para essa tarefa
- ```-t 24:00:00``` = 24 horas
```--pty bash -i``` = abre um shell interativo

Depois, definimos as variáveis R1 e R2 como, respectivamente, todos os reads filtrados R1 e todos os reads filtrados R2:
```
R1="/labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-RF-1-B_S1_L001_R1.fastq.gz,/labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-RF-2-B_S2_L001_R1.fastq.gz,/labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-RF-4-B_S3_L001_R1.fastq.gz,/labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-RF-5-B_S4_L001_R1.fastq.gz,/labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-PM-1-A_R1.fastq.gz,/labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-PM-5-A_S7_L001_R1.fastq.gz"

R2="/labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-RF-1-B_S1_L001_R2.fastq.gz,/labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-RF-2-B_S2_L001_R2.fastq.gz,/labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-RF-4-B_S3_L001_R2.fastq.gz,/labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-RF-5-B_S4_L001_R2.fastq.gz,/labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-PM-1-A_R2.fastq.gz,/labgenomaarea2/valentina.pavelecini/EHI/filtrados/TF-2587-PM-5-A_S7_L001_R2.fastq.gz"
```

Depois, criamos uma parta dentro de Megahit chamada ```coassembly```, em que será salvo o output do comando, que vai ser mais ou menos assim:

```
megahit \
    -t 24 \
    --verbose \
    --min-contig-len 1500 \
    -1 "$R1" \
    -2 "$R2" \
    -f \
    -o /labgenomaarea2/valentina.pavelecini/EHI/Megahit/coassembly \
    &> /labgenomaarea2/valentina.pavelecini/EHI/Megahit/coassembly/megahit.log
```
Em que:
- ```-t 24``` = usa os 24 CPUs que reservamos.
- ```-1 "$R1"``` = usa as seis amostras R1 que colocamos na variável.
- ```-2 "$R2"``` = usa as seis amostras R2 que colocamos na variável.
- ```--min-contig-len 1500``` = mantém contigs (sequências de DNA reconstruídas a partir de vários fragmentos menores/reads) maior que 1500 bp.
- ```-f``` = força a execução caso o diretório de saída já exista.
- ```-o /labgenomaarea2/valentina.pavelecini/EHI/Megahit/coassembly \``` = coloca os resultados aqui.
- ```&> /labgenomaarea2/valentina.pavelecini/EHI/Megahit/coassembly/megahit.log``` = salva output no arquivo .log.

Output:
Arquivo foi copiado para o computador com o comando abaixo e está disponível neste repositório em ```megahit.log```
```
scp -J valentina.pavelecini@marfim.lad.pucrs.br `
valentina.pavelecini@pantanal.lad.pucrs.br:/labgenomaarea2/valentina.pavelecini/EHI/Megahit/coassembly/megahit.log `
C:\Users\labgenoma\Downloads\
```


Arquivos criados:
- ```final.contigs.fa``` = resultado final do coassembly.
- ```done``` = marcador criado pelo MEGAHIT indicando que a execução terminou.
- ```megahit.log``` = log da execução para verificar se houve erros.
- ```log``` = informações da execução.
- ```options.json``` = parâmetros usados pelo MEGAHIT.
- ```checkpoints.txt``` = informações sobre o progresso/checkpoints.
- ```intermediate_contigs``` = pasta com resultados intermediários produzidos durante as diferentes etapas de assembly.



## 05/10/2026

Vamos fazer também a etapa do QUAST, para obtermos algumas estatísticas gerais do assembly:
```
quast \
    -o {output.report} \
    --threads {threads} \
    {input.assembly}

```

Vamos fazer isso na janela do tmux que criamos antes, entrando nela com ```tmux at```. Vamos pedir ```srun -N 1 -n 1 -c 24 -t 6:00:00 --pty bash -i``` (1 máquina, 1 tarefa, 24 CPUs/cores para essa tarefa, por 6 horas). Menos tempo dessa vez porque é um comando mais leve. Vamos criar também um diretório específico para o QUAST com ```mkdir QUAST```. Então, vamos baixar o software com ```conda install -c bioconda quast```.

E no nosso caso, para rodar o comando, utilizamos o ```final.contigs.fa```. E claro, os outputs serão salvos em ```quast.log```.

```
quast \
    -o /labgenomaarea2/valentina.pavelecini/EHI/Megahit/QUAST \
    --threads 24 \
    /labgenomaarea2/valentina.pavelecini/EHI/Megahit/coassembly/final.contigs.fa \
    &> /labgenomaarea2/valentina.pavelecini/EHI/Megahit/QUAST/quast.log
```

## 07/10/2026

Interpretando os dados de report.tsv:

Número total de contigs    63.170                             A montagem é bastante fragmentada
Comprimento total          256.725.534 bp (~256,7 Mb)         Quantidade total de sequência montada
Maior contig	           542.717 bp	                      O maior fragmento montado tem ~543 kb
N50                        5.267 bp                           Metade da montagem está em contigs ≥ 5,3 kb
N90                        1.736 bp                           90% da montagem está em contigs ≥ 1,7 kb
Contigs ≥ 5 kb             8.203                              8.203 contigs têm pelo menos 5 kb
Contigs ≥ 10 kb            3.333	                          3.333 contigs têm pelo menos 10 kb
Contigs ≥ 25 kb            1.044                              1.044 contigs têm pelo menos 25 kb
Contigs ≥ 50 kb            423	                              423 contigs têm pelo menos ≥50 kb
GC                         47,24%	                          Conteúdo Guanina/Citosina médio da montagem é 47,24%
N's/100 kbp                0                                  Não há bases ambíguas N na montagem


# Mapeamento do agrupamento:

Assim, a próxima etapa é mapear gerar informações de cobertura dos contigs do agrupamento (coassembly), que serão usadas pelos programas de binning para agrupar contigs em possíveis MAGs.

_______________________________________________________________________________________

Na primeira etapa, organizaremos os contigs do agrupamento em um índice com o bowtie2:

```
bowtie2-build \
    --large-index \
    --threads {threads} \
    {input.contigs} {output.contigs}
```

Para isso, primeiro vamos abrir a janela do tmux com ```tmux new -s bowtie2```, e iniciar os parâmetros da sessão com ```srun -N 1 -n 1 -c 24 -t 6:00:00 --pty bash -i```, os mesmos de antes. Também precisamos baixar o bowtie2 e fazemos isso com ```conda install -c bioconda bowtie2```. Depois disso, já podemos rodar o comando:

```
bowtie2-build \
    --large-index \
    --threads 24 \
    /labgenomaarea2/valentina.pavelecini/EHI/Megahit/coassembly/final.contigs.fa \
    /labgenomaarea2/valentina.pavelecini/EHI/bowtie2/3.2.bowtie2.contigs \
    &> /labgenomaarea2/valentina.pavelecini/EHI/bowtie2/bowtie2contigs.log
```

Em que ```3.2``` se refere á etapa da pipeline do EHI, ```bowtie2``` ao programa usado e ```contigs``` ao conteúdo
