# Login + site do Earth Hologenome Iniciative + Fluxograma
```
ssh -J valentina.pavelecini@marfim.lad.pucrs.br valentina.pavelecini@pantanal.lad.pucrs.br  
cd /labgenomaarea2/valentina.pavelecini  
conda activate EHI  
```

https://www.earthhologenome.org/bioinformatics/index.html  
https://www.canva.com/design/DAHUDcTIaEc/Z9jiEwyVO9jZrfU_UPkavA/edit?ui=e30  

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

Começaremos instalando o nonparail com ```conda install bioconda::nonpareil```. Depois, criamos um diretório chamado ```nonpareil``` com ```mkdir nonpareil```, onde irão os arquivos de saída.  

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
