# Previsão da mortalidade com o modelo de Lee–Carter: uma introdução

Este repositório reúne o código em **R** utilizado na aplicação empírica do artigo *Previsão da mortalidade com o modelo de Lee–Carter: uma introdução*, desenvolvido para os **Cadernos de Estatística da UFF**. O objetivo é ilustrar a estimação bayesiana do modelo de Lee–Carter, a avaliação do ajuste e a previsão das taxas de mortalidade dos Estados Unidos, incluindo uma aplicação à esperança de vida temporária.

## Dados e período de estudo

As taxas centrais de mortalidade são provenientes da [Human Mortality Database (HMD)](https://www.mortality.org/), para a população **total (ambos os sexos)** dos Estados Unidos. Utilizam-se **21 faixas etárias** — `0`, `1-4`, `5-9`, ..., `95-99` — e observações anuais entre **1970 e 2010**.

O período **1970–2000** é utilizado para estimar o modelo (31 anos), enquanto **2001–2010** é reservado à avaliação das previsões fora da amostra. Os dados posteriores a 2000 não entram na estimação.

O arquivo de entrada esperado é `Mx_USA.txt`, contendo, entre outras, as colunas `Year`, `Age` e `Total` no formato da HMD. Consulte as condições de acesso e utilização da base no site da HMD.

## Arquivos principais

```text
.
├── README.md
├── Lee_Carter_Analise_Completa_USA.R
├── Mx_USA.txt                         # arquivo de entrada (obtido na HMD)
├── Figura_Motivacao_USA_1970_2010_estiloR_branco.pdf
├── Figura_2_USA_fundo_branco.pdf
├── Resumo_Variacao_Anual_USA_1971_2010_R.csv
└── resultados_LC_USA/               # criado automaticamente pelo script
```

O arquivo `.R` é o script principal. As Figuras 1 e 2 e o resumo exploratório são gravados no diretório de trabalho; os resultados da estimação, projeção e diagnósticos são gravados em `resultados_LC_USA/`.

## Requisitos

- **R**, com funções básicas para leitura de dados, cálculos e gráficos.
- **BayesMortalityPlus**, para o ajuste bayesiano e as previsões do modelo.
- **coda**, para os diagnósticos ESS e Geweke (opcional, mas recomendado).
- **fields**, para a barra de cores do mapa de resíduos apresentada no script.

Instale os pacotes necessários uma única vez:

```r
install.packages(c("BayesMortalityPlus", "coda", "fields"))
```

## Como reproduzir a análise

1. Coloque `Mx_USA.txt` e `Lee_Carter_Analise_Completa_USA.R` no mesmo diretório.
2. Abra o R nesse diretório e execute:

   ```r
   source("Lee_Carter_Analise_Completa_USA.R")
   ```

3. Consulte os arquivos PDF, CSV e RDS produzidos. O ajuste bayesiano pode levar algum tempo, dependendo do computador.

**Nota sobre o nome do arquivo:** se o script baixado estiver com o sufixo `(1)` — `Lee_Carter_Analise_Completa_USA(1).R` —, renomeie-o para `Lee_Carter_Analise_Completa_USA.R` ou utilize o nome efetivo na chamada a `source()`.

## Modelo e configuração computacional

O ajuste utiliza a implementação bayesiana do modelo de Lee–Carter disponibilizada pela função `BayesMortalityPlus::blc()`, aplicada ao logaritmo das taxas centrais de mortalidade. O modelo considera um índice temporal latente com evolução estocástica e produz amostras posteriores para os parâmetros e estados.

No script fornecido, são utilizadas as seguintes configurações:

| Configuração | Valor |
|---|---:|
| Semente aleatória | `20261105` |
| Número total de iterações (`M`) | `20000` |
| Aquecimento / burn-in (`bn`) | `10000` |
| Intervalo de armazenamento (`thin`) | `4` |
| Número de cadeias | `1` |

A configuração corresponde a **2.500 amostras posteriores armazenadas após o burn-in**, conforme a convenção de armazenamento adotada no script. Para avaliar a amostragem, o código calcula o tamanho efetivo da amostra (ESS) e a estatística de Geweke, além de gerar exemplos de gráficos de trajetória. Esses diagnósticos são complementares e, isoladamente, não garantem convergência.

## Figuras produzidas

| Figura | Arquivo | Conteúdo |
|---|---|---|
| 1 | `Figura_Motivacao_USA_1970_2010_estiloR_branco.pdf` | Evolução das log-taxas de mortalidade em seis faixas etárias selecionadas, 1970–2010 |
| 2 | `Figura_2_USA_fundo_branco.pdf` | Mediana anual das variações de log-mortalidade e dispersão entre faixas etárias (MAD) |
| 3 | `resultados_LC_USA/Figura_3_Posteriores_LeeCarter_USA.pdf` | Medianas e intervalos de credibilidade de 95% de αₓ, βₓ e κₜ |
| 4 | `resultados_LC_USA/Figura_4_Previsoes_LC_USA.pdf` | Taxas observadas, medianas preditivas e intervalos preditivos de 95% para 2001–2010 |
| 5 | `resultados_LC_USA/Figura_5_Residuos_LC_USA.pdf` | Mapa dos resíduos padronizados e evolução da variância residual por ano |
| 6 | `resultados_LC_USA/Figura_6_Esperanca_Vida_Temporaria_USA.pdf` | Esperança de vida temporária até os 100 anos ao nascer e aos 60 anos: observado e previsto |

As seis faixas etárias destacadas nas Figuras 1 e 4 são `0`, `20-24`, `40-44`, `60-64`, `80-84` e `90-94`. O ajuste utiliza todas as 21 faixas.

## Resultados adicionais

A execução também produz arquivos auxiliares, como:

- `Resumo_Variacao_Anual_USA_1971_2010_R.csv`: medidas exploratórias anuais.
- `resultados_LC_USA/outputBLC.rds`: objeto do ajuste e configurações da simulação.
- `resultados_LC_USA/diagnosticos_mcmc.csv`: ESS e estatísticas de Geweke, quando `coda` está instalado.
- `resultados_LC_USA/Diagnosticos_tracos_LC_USA.pdf`: exemplos de trajetórias MCMC.
- `resultados_LC_USA/previsoes_LC_2001_2010.csv`: taxas observadas, medianas, intervalos preditivos e indicação de cobertura.
- `resultados_LC_USA/avaliacao_previsoes.csv`: MAE, RMSE, cobertura empírica e largura média dos intervalos.
- `resultados_LC_USA/residuos_padronizados_mediana.csv` e `variancia_residual_anual.csv`: resumos dos resíduos.
- `resultados_LC_USA/esperanca_vida_temporaria_previsoes.csv`: valores observados e resumos preditivos da esperança de vida temporária.

Os arquivos `.rds` também preservam objetos e amostras posteriores para análises complementares.

## Esperança de vida temporária

A esperança de vida é calculada a partir das taxas centrais observadas e preditivas, assumindo **força de mortalidade constante dentro de cada intervalo etário**, com larguras de 1, 4 e 5 anos. A medida é **temporária até os 100 anos**: ao nascer, corresponde aos anos esperados vividos até essa idade; aos 60 anos, corresponde aos anos adicionais esperados entre 60 e 100 anos.

Essa medida **não é a esperança de vida completa**, pois não inclui sobrevivência acima dos 100 anos, e corresponde a uma medida **de período**, baseada no perfil das taxas do ano considerado. O cálculo é aplicado individualmente às amostras preditivas para produzir medianas e intervalos de 95%.

## Observações sobre a reprodutibilidade

O código realiza uma única cadeia de amostragem; por isso, não apresenta o diagnóstico entre cadeias R-hat. As bandas das Figuras 3 e 4 possuem interpretações distintas: na Figura 3, representam **intervalos de credibilidade dos parâmetros**; na Figura 4, **intervalos preditivos das taxas de mortalidade**.

Para executar o script integralmente em modo não interativo, verifique a chamada `plot(ajuste1, paramter = "fitted")` no código enviado: o nome esperado do argumento é `parameter`, com **e**. Essa correção não altera o modelo nem os cálculos.

## Referências

- **Lee, R. D.; Carter, L. R. (1992).** Modeling and forecasting U.S. mortality. *Journal of the American Statistical Association*, 87(419), 659–671.
- **Pedroza, C. (2006).** A Bayesian forecasting model: predicting U.S. male mortality. *Biostatistics*, 7(4), 530–550.
- **Human Mortality Database.** https://www.mortality.org/
- **BayesMortalityPlus:** pacote R utilizado para a implementação bayesiana do Lee–Carter.

## Citação do artigo

Este repositório acompanha o artigo *Previsão da mortalidade com o modelo de Lee–Carter: uma introdução*, preparado para os **Cadernos de Estatística da UFF**. Os dados bibliográficos definitivos (autoria, ano, volume e DOI, se houver) devem ser acrescentados após a publicação.
