# ARTIGO LEE-CARTER — TODAS AS FIGURAS (1 a 6)
# Arquivo de entrada: Mx_USA.txt no diretorio de trabalho.
# Figuras 1 e 2: R base, sem pacote extra.
# Figuras 3 a 6: ajuste bayesiano, BayesMortalityPlus.
# Execucao: source("Lee_Carter_Analise_Completa_USA.R")

# Cadernos de Estatística da UFF
# Lee-Carter — EUA, 1970–2010: figuras 1 e 2 e resumo exploratório
# Reproduzível apenas com funções do R base.
# Coloque Mx_USA.txt na pasta de trabalho e execute source("Lee_Carter_Analise_Completa_USA.R").

arquivo <- "Mx_USA.txt"
if (!file.exists(arquivo)) stop("Arquivo Mx_USA.txt não encontrado no diretório de trabalho.")

Mx <- read.table(arquivo, header = TRUE, na.strings = ".", stringsAsFactors = FALSE)
faixas <- c("0", "1-4", paste0(seq(5, 95, 5), "-", seq(9, 99, 5)))
anos <- 1970:2010

dados <- subset(Mx, Year %in% anos & Age %in% faixas,
                select = c("Year", "Age", "Total"))
if (anyNA(dados$Total) || any(dados$Total <= 0)) {
  stop("Há taxas ausentes ou não positivas no período e faixas selecionados.")
}
dados <- dados[order(match(dados$Age, faixas), dados$Year), ]
stopifnot(nrow(dados) == length(faixas) * length(anos))

# Matriz com idades nas linhas e anos nas colunas.
mx <- matrix(dados$Total, nrow = length(anos), ncol = length(faixas))
mx <- t(mx)
dimnames(mx) <- list(faixas, as.character(anos))
logmx <- log(mx)

# Figura 1: seis faixas selecionadas, gráficos em duas linhas e três colunas.
selecionadas <- c("0", "20-24", "40-44", "60-64", "80-84", "90-94")
pdf("Figura_Motivacao_USA_1970_2010_estiloR_branco.pdf", width = 10, height = 5.7,
    bg = "white")
par(mfrow = c(2, 3), mar = c(3.5, 5.5, 2.0, 0.6) + 0.3, mgp = c(2.8, 0.7, 0),las = 1, bty = "o", bg = "white")
for (i in seq_along(selecionadas)) {
  fx <- selecionadas[i]
  plot(anos, logmx[fx, ],type = "l", col = "black", lwd = 1.1, main = paste("idade:", fx),cex.main = 1.2,
       cex.lab = 1.4,cex.axis = 1.4, xlab = if (i > 3) "ano" else "",
       ylab = if (i %% 3 == 1) expression(log(m[xt])) else "",xlim = c(1970, 2010),xaxt = "n")
  axis(1, at = c(1970, 1980, 1990, 2000, 2010), labels = i > 3, cex.axis = 1.4)
  abline(v = 2000, lty = 2, col = "gray50")
}
dev.off()

# Diferenças anuais das log-taxas: 21 idades x 40 anos.
delta <- logmx[, -1, drop = FALSE] - logmx[, -ncol(logmx), drop = FALSE]
M <- apply(delta, 2, median)
D <- apply(delta, 2, function(z) median(abs(z - median(z))))
anos_delta <- 1971:2010
resumo <- data.frame(ano = anos_delta, mediana = M, dispersao_MAD = D)
write.csv(resumo, "Resumo_Variacao_Anual_USA_1971_2010_R.csv", row.names = FALSE)

# Figura 2: dois painéis horizontais; FUNDO BRANCO, SEM GRID, SEM TÍTULOS.
pdf("Figura_2_USA_fundo_branco.pdf", width = 9.2, height = 3.2,
    bg = "white")
par(mfrow = c(1, 2), mar = c(4.5, 6.5, 1.5, 1.0),mgp = c(4.2, 0.8, 0), las = 1, bty = "o", bg = "white")
plot(anos_delta, M, type = "l", lwd = 1.1, col = "black",
     xlab = "ano", ylab = "mediana da variação anual",
     xlim = c(1970, 2010), xaxt = "n", cex.lab=1.2, cex=1.2, cex.axis=1.2)
axis(1, at = c(1970, 1980, 1990, 2000, 2010), cex.axis=1.2)
plot(anos_delta, D, type = "l", lwd = 1.1, col = "black",
     xlab = "ano", ylab = "dispersão entre idades",
     xlim = c(1970, 2010), xaxt = "n", cex.lab=1.2, cex=1.2, cex.axis=1.2)
axis(1, at = c(1970, 1980, 1990, 2000, 2010), cex.axis=1.2)
dev.off()

# ========================================================
# ESTIMACAO, PREVISAO, RESIDUOS E ESPERANCA DE VIDA
# ========================================================
# Cadernos de Estatistica da UFF
# Secao 4.2 — Estimacao bayesiana do Lee-Carter com BayesMortalityPlus
# Arquivo de entrada: Mx_USA.txt (HMD)
# Apenas 1970–2000 e usado no ajuste; 2001–2010 fica fora da amostra.
# Saidas: figura, tabelas CSV, objetos RDS e diagnosticos MCMC.

if (!requireNamespace("BayesMortalityPlus", quietly = TRUE)) {
  stop("Instale o pacote primeiro: install.packages('BayesMortalityPlus')")
}

set.seed(20261105)
arquivo <- "Mx_USA.txt"
if (!file.exists(arquivo)) stop("Nao encontrei Mx_USA.txt no diretorio de trabalho.")
dir.create("resultados_LC_USA", showWarnings = FALSE)
saida <- "resultados_LC_USA"

mx_us <- read.table(arquivo, header = TRUE, na.strings = ".", stringsAsFactors = FALSE)
idades <- c("0", "1-4", paste0(seq(5, 95, 5), "-", seq(9, 99, 5)))
anos_ajuste <- 1970:2000
anos_previsao <- 2001:2010

# Y: matriz 21 faixas etarias x 31 anos — log da taxa central Total.
# Usar 'match' evita reordenacao incorreta de faixas como '10-14' e '5-9'.
Y <- t(vapply(idades, function(faixa) {
  z <- mx_us[mx_us$Age == faixa & mx_us$Year %in% anos_ajuste, ]
  z <- z[match(anos_ajuste, z$Year), ]
  if (nrow(z) != length(anos_ajuste) || anyNA(z$Total) || any(z$Total <= 0)) {
    stop("Dados faltantes, nao positivos ou anos incompletos para a faixa ", faixa)
  }
  log(z$Total)
}, numeric(length(anos_ajuste))))
dimnames(Y) <- list(idades, as.character(anos_ajuste))
stopifnot(identical(dim(Y), c(21L, 31L)), all(is.finite(Y)))

# M e o total de iteracoes antes do descarte; bn e o burn-in.
M <- 20000
bn <- 10000
thin <- 4

ajuste1 <- BayesMortalityPlus::blc(Y = Y, M = M, bn = bn, thin = thin)

saveRDS(list(cadeia1 = ajuste1,
             anos_ajuste = anos_ajuste, idades = idades,
             M = M, bn = bn, thin = thin),
        file.path(saida, "outputBLC.rds"))

plot(ajuste1, paramter= "fitted")
plot(ajuste1) ## alpha , beta , kappa


# Mediana posterior e intervalos de credibilidade de 95%
posterior <- function(mat, nomes) {
  stopifnot(is.matrix(mat), nrow(mat) == length(nomes))
  data.frame(indice = nomes, mediana = apply(mat, 1, median),
             li95 = apply(mat, 1, quantile, probs = 0.025),
             ls95 = apply(mat, 1, quantile, probs = 0.975),
             row.names = NULL)
}

alpha <- posterior(as.matrix(ajuste1$alpha), idades)
beta  <- posterior(as.matrix(ajuste1$beta), idades)
kappa <- posterior(as.matrix(ajuste1$kappa), as.character(anos_ajuste))

# Parametros adicionais: drift theta; precises phiw e phiv.
# Para obter variancias, transformar CADA sorteio de precisao antes de resumir.
# Parâmetros adicionais: drift theta e variância de evolução
resumir <- function(x) c(mediana = median(x), li95 = unname(quantile(x, .025)), ls95 = unname(quantile(x, .975)))

resumo_estaticos <- rbind(
  theta = resumir(as.numeric(ajuste1$theta)),
  sigma2_w = resumir(1 / as.numeric(ajuste1$phiw))
)

print(resumo_estaticos)

# Figura dos parametros (mediana e intervalos de credibilidade de 95%).
# Base R, fundo branco, sem grid, estilo igual as figuras exploratorias.
# Figura 3: parâmetros do Lee-Carter (mediana e IC de 95%)
plot_intervalo <- function(z, eixo_x, eixo_y, eixo_x_nome, rotulos = NULL) {
  ylim <- range(c(z$li95, z$ls95))
  folga <- diff(range(eixo_x)) * 0.06
  plot(eixo_x, z$mediana, type = "n", ylim = ylim, xlim = range(eixo_x) + c(-folga, folga), xlab = eixo_x_nome, ylab = eixo_y, xaxt = "n", cex.lab = 2, cex.axis = 2, cex = 2)
  if (!is.null(rotulos)) {
    pos <- seq(1, length(eixo_x), by = 4)
    axis(1, at = eixo_x[pos], labels = rotulos[pos], las = 1, cex.axis = 2)
  } else {
    axis(1, at = seq(1970, 2000, by = 10), cex.axis = 2)
  }
  polygon(c(eixo_x, rev(eixo_x)), c(z$li95, rev(z$ls95)), col = "gray85", border = NA)
  lines(eixo_x, z$mediana, col = "black", lwd = 1.3)
}

pdf(file.path(saida, "Figura_3_Posteriores_LeeCarter_USA.pdf"), width = 14.3, height = 3.5, bg = "white")

par(mfrow = c(1, 3), mar = c(4.5, 6.8, 1.5, 1.5) +0.1, mgp = c(3.2, 0.8, 0), las = 1, bty = "o", bg = "white")

plot_intervalo(alpha, seq_along(idades), expression(alpha[x]), "faixa etária", rotulos = idades)
plot_intervalo(beta, seq_along(idades), expression(beta[x]), "faixa etária", rotulos = idades)
abline(h = 0, lty = 2, col = "gray60")
plot_intervalo(kappa, anos_ajuste, expression(kappa[t]), "ano")

dev.off()

# Diagnostico visual das cadeias (exemplos: theta e kappa em 1985).
pdf(file.path(saida, "Diagnosticos_tracos_LC_USA.pdf"),
    width = 9.0, height = 5.2, bg = "white")
par(mfrow = c(2,1), mar = c(3.4,4.2,1.0,0.8), bg = "white")
plot(ajuste1$theta, type = "l", col = "black", xlab = "iteracao salva",
     ylab = expression(theta))
pos <- match("1985", colnames(Y))
plot(ajuste1$kappa[pos, ], type = "l", col = "black", xlab = "iteracao salva",
     ylab = expression(kappa[1985]))
dev.off()

# R-hat e tamanho efetivo por componente, se o pacote coda estiver instalado.
# Diagnósticos MCMC para uma única cadeia: ESS e Geweke
if (requireNamespace("coda", quietly = TRUE)) {
  
  pnames <- c(paste0("alpha_", idades), paste0("beta_", idades),
              paste0("kappa_", anos_ajuste), "theta", "sigma2_w")
  
  x <- cbind(t(ajuste1$alpha), t(ajuste1$beta), t(ajuste1$kappa),
             theta = as.numeric(ajuste1$theta),
             sigma2_w = 1 / as.numeric(ajuste1$phiw))
  
  colnames(x) <- pnames
  cadeia <- coda::mcmc(x)
  
  ess <- coda::effectiveSize(cadeia)
  geweke <- coda::geweke.diag(cadeia)$z
  
  diag <- data.frame(parametro = names(ess), ESS = as.numeric(ess),
                     Geweke = as.numeric(geweke[names(ess)]))
  
  write.csv(diag, file.path(saida, "diagnosticos_mcmc.csv"), row.names = FALSE)
  
  cat("ESS mínimo:", min(diag$ESS), "\n")
  cat("Percentual de |Geweke| > 1.96:", mean(abs(diag$Geweke) > 1.96, na.rm = TRUE) * 100, "%\n")
  
} else {
  message("Para ESS e Geweke, instale: install.packages('coda')")
}
############################################################
## 4.3 — PROJECAO POSTERIOR DA MORTALIDADE (2001-2010)
## 4.4 — DIAGNOSTICOS DE RESIDUOS DO LEE-CARTER
## Utiliza os ajustes salvos pela etapa anterior.
## NENHUM ano posterior a 2000 participa da estimacao.
############################################################

# A versao do pacote pode ser consultada por packageVersion("BayesMortalityPlus").
# Esta secao simula as distribuicoes preditivas a partir das cadeias,
# para reproduzir graficos no estilo do manuscrito mortality jump.
# Distintamente dos graficos com 3 modelos, aqui mostramos APENAS Lee-Carter.

# Extrair draws orientados no formato [idade ou ano, amostra].
# A dimensao das amostras do objeto BLC segue sua documentacao.
# Extrair amostras posteriores do modelo Lee-Carter
extrair_matriz <- function(ajuste, nome, nlinhas) {
  v <- ajuste[[nome]]
  if (is.null(dim(v)) || length(dim(v)) != 2L)
    stop("Esperava matriz em ", nome, ". Inspecione str(ajuste$", nome, ").")
  if (nrow(v) == nlinhas) return(v)
  if (ncol(v) == nlinhas) return(t(v))
  stop("Dimensoes inesperadas de ", nome, ": ", paste(dim(v), collapse = " x "))
}

if (!exists("ajuste1"))
  stop("Execute a estimacao do modelo antes da previsao.")

A <- extrair_matriz(ajuste1, "alpha", length(idades))
B <- extrair_matriz(ajuste1, "beta", length(idades))
K <- extrair_matriz(ajuste1, "kappa", length(anos_ajuste))
PHI_V <- extrair_matriz(ajuste1, "phiv", length(idades))
THETA <- as.numeric(ajuste1$theta)
PHI_W <- as.numeric(ajuste1$phiw)

S <- ncol(A)
stopifnot(ncol(B) == S, ncol(K) == S, ncol(PHI_V) == S,
          length(THETA) == S, length(PHI_W) == S,
          all(PHI_V > 0), all(PHI_W > 0))

# Dados observados apenas para comparar as previsoes
extrair_taxas <- function(anos) {
  z <- t(vapply(idades, function(faixa) {
    d <- mx_us[mx_us$Age == faixa & mx_us$Year %in% anos,]
    d <- d[match(anos, d$Year),]
    if (nrow(d) != length(anos) || anyNA(d$Total) || any(d$Total <= 0))
      stop("Taxas observadas ausentes para a faixa ", faixa)
    d$Total
  }, numeric(length(anos))))
  dimnames(z) <- list(idades, as.character(anos))
  z
}

MX_TESTE <- extrair_taxas(anos_previsao)
MX_TREINO <- exp(Y)

# Previsao nativa do BayesMortalityPlus
pred_nativa_1 <- predict(ajuste1, h = length(anos_previsao))
saveRDS(pred_nativa_1, file.path(saida, "previsoesBayesMortalityPlus.rds"))

# Simulacao posterior preditiva coerente com:
# kappa[t+h] = kappa[t+h-1] + theta + N(0, 1/phiw)
# log(mx)[x,t+h] = alpha[x] + beta[x]*kappa[t+h] + N(0, 1/phiv[x])
# Portanto, considera incerteza de parametros, do estado futuro e da observacao.

H <- length(anos_previsao)
MX_PRED <- array(NA_real_, c(length(idades), H, S),
                 dimnames=list(idades, as.character(anos_previsao), NULL))
kappa_futuro <- K[nrow(K), ]
for (h in seq_len(H)) {
  kappa_futuro <- kappa_futuro + THETA + rnorm(S, 0, sqrt(1/PHI_W))
  log_pred <- sweep(B, 2, kappa_futuro, `*`) + A
  log_pred <- log_pred + matrix(rnorm(length(idades)*S),length(idades),S) *
                         sqrt(1/PHI_V)
  MX_PRED[,h,] <- exp(log_pred)
}
saveRDS(MX_PRED,file.path(saida,"amostras_preditivas_mx_2001_2010.rds"))

# Resumir as amostras preditivas por idade x ano.
pred_med <- apply(MX_PRED, c(1,2), median)
pred_li  <- apply(MX_PRED, c(1,2), quantile, probs=.025)
pred_ls  <- apply(MX_PRED, c(1,2), quantile, probs=.975)

previsao_csv <- expand.grid(idade=idades, ano=anos_previsao,
                            KEEP.OUT.ATTRS=FALSE, stringsAsFactors=FALSE)
previsao_csv$observado <- as.vector(MX_TESTE)
previsao_csv$mediana <- as.vector(pred_med)
previsao_csv$li95 <- as.vector(pred_li)
previsao_csv$ls95 <- as.vector(pred_ls)
previsao_csv$coberto <- with(previsao_csv, observado>=li95 & observado<=ls95)
write.csv(previsao_csv,file.path(saida,"previsoes_LC_2001_2010.csv"),row.names=FALSE)

# Grafico no estilo dos exemplos: observado, mediana e faixa de previsao.
# Idades selecionadas apenas para visualizacao; o modelo usa todas as 21.
#idades_grafico <- c("20-24", "40-44", "55-59", "60-64", "80-84", "90-94")
idades_grafico <- c("0", "20-24", "40-44", "60-64", "80-84", "90-94")

cores <- list(linha="black", mediana="#205b92", faixa="gray84")
pdf(file.path(saida,"Figura_4_Previsoes_LC_USA.pdf"),
    width = 10, height = 5.7,bg="white")
# par(mfrow=c(2,3),mar=c(3.2,3.4,1.7,.5),mgp=c(2.0,.55,0),
#     las=1,bty="o",bg="white",cex.axis=.77,cex.lab=.85)
par(mfrow=c(2,3),mar=c(3.2,3.4,1.7,.5),mgp=c(2.0,.55,0),
    las=1,bty="o",bg="white",cex.axis=.77,cex.lab=.85)

for (faixa in idades_grafico) {
  i <- match(faixa, idades)
  mx_obs <- c(MX_TREINO[i,], MX_TESTE[i,])
  yr_all <- c(anos_ajuste, anos_previsao)
  lim <- range(mx_obs, pred_li[i,], pred_ls[i,])
  
  plot(yr_all, mx_obs, type = "n", xlab = "ano", ylab = expression(m[x,t]), ylim = lim, main = paste("idade", faixa), xlim = range(yr_all), cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.4)
  par(xpd = FALSE)
  
  polygon(c(anos_previsao, rev(anos_previsao)), c(pred_li[i,], rev(pred_ls[i,])), border = NA, col = cores$faixa)
  lines(yr_all, mx_obs, col = cores$linha, lwd = 1, lty = 2)
  lines(anos_previsao, pred_med[i,], col = cores$mediana, lwd = 1.6)
  
  segments(2000.5, par("usr")[3], 2000.5, par("usr")[4], lty = 2, lwd = 0.8, col = "black")
  
  if (faixa == idades_grafico[1]) {
    legend("bottomleft", inset = c(0.02, 0.02), bty = "n", cex = 0.85,
           legend = c("observado", "mediana preditiva", "IP 95%"),
           lty = c(2, 1, NA), lwd = c(1, 1.6, NA),
           col = c("black", cores$mediana, NA),
           fill = c(NA, NA, cores$faixa), border = NA)
  }
}
dev.off()

############################################################
## RESIDUOS PADRONIZADOS (AJUSTE 1970-2000)
############################################################
# Em cada amostra posterior s:
# r[x,t,s] = {log(mx)[x,t] - alpha[x,s] - beta[x,s]*kappa[t,s]}
#            / sigma[x,s], com sigma^2[x,s] = 1/phiv[x,s].
# E[x,t] eh mediana posterior do residuo padronizado.
# V[t,s] eh a variancia amostral entre as 21 idades dos residuos padronizados.
# Esta eh uma estatistica de diagnostico, NAO e sigma_x^2/lambda_t do DynFull.
############################################################
# FIGURA 5: ANALISE DOS RESIDUOS DO LEE-CARTER
############################################################

X <- length(idades)
Tn <- length(anos_ajuste)

# Residuos padronizados e variancia residual por ano
R_all <- array(NA_real_, dim = c(X, Tn, S))
V <- matrix(NA_real_, nrow = Tn, ncol = S)

for (s in seq_len(S)) {
  fitted_s <- sweep(outer(B[,s], K[,s]), 1, A[,s], "+")
  resid_s <- sweep(Y - fitted_s, 1, sqrt(PHI_V[,s]), "*")
  R_all[,,s] <- resid_s
  V[,s] <- apply(resid_s, 2, var)
}

# Resumos posteriores
R_med <- apply(R_all, c(1, 2), median)
V_med <- apply(V, 1, median)
V_li <- apply(V, 1, quantile, probs = 0.025)
V_ls <- apply(V, 1, quantile, probs = 0.975)

# Mapa de calor: azul negativo, branco zero, vermelho positivo
mapa_lim <- 3
pal <- grDevices::colorRampPalette(c("#5e94be", "white", "#d58e78"))(101)
Z <- t(pmax(pmin(R_med, mapa_lim), -mapa_lim))

stopifnot(identical(dim(Z), c(Tn, X)))

# Figura 5: mapa de calor e variância residual lado a lado
# Figura 5: análise dos resíduos do Lee-Carter

pdf(file.path(saida, "Figura_5_Residuos_LC_USA.pdf"), width = 12, height = 4.5, bg = "white")
old_par <- par(no.readonly = TRUE)

par(mfrow = c(1, 2), mar = c(5, 7, 2, 2) + 0.2, oma = c(0, 0, 0, 3))

# Painel 1: mapa de calor dos resíduos padronizados
image(anos_ajuste, seq_len(X), Z, col = pal, zlim = c(-mapa_lim, mapa_lim), xlab = "", ylab = "", xaxt = "n", yaxt = "n", useRaster = TRUE)
axis(1, at = seq(1970, 2000, by = 10), cex.axis = 1.4)
axis(2, at = seq(1, X, by = 4), labels = idades[seq(1, X, by = 4)], las = 1, cex.axis = 1.4)
mtext("ano", side = 1, line = 3, cex = 1.4)
mtext("faixa etária", side = 2, line = 4.8, cex = 1.4)
box()

# Barra de cores: curta, à direita da figura
fields::image.plot(legend.only = TRUE, zlim = c(-mapa_lim, mapa_lim), col = pal,
                   smallplot = c(0.955, 0.965, 0.25, 0.75),
                   axis.args = list(at = -3:3, las = 1, cex.axis = 1.0),
                   legend.args = list(text = "", side = 4, line = 2, cex = 0.9))


# Painel 2: variância dos resíduos padronizados
plot(anos_ajuste, V_med, type = "n", ylim = range(c(V_li, V_ls, 1)), xaxt = "n", xlab = "", ylab = "", cex.axis = 1.4)
polygon(c(anos_ajuste, rev(anos_ajuste)), c(V_li, rev(V_ls)), col = grDevices::adjustcolor("gray70", alpha.f = 0.65), border = NA)
lines(anos_ajuste, V_med, col = "black", lwd = 2)
axis(1, at = seq(1970, 2000, by = 10), cex.axis = 1.4)
mtext("ano", side = 1, line = 3, cex = 1.4)
mtext("variância residual padronizada", side = 2, line = 4.8, cex = 1.4)
abline(h = 1, lty = 2, col = "gray60")
box()


par(old_par)
dev.off()



# Salvar resultados
write.csv(R_med, file.path(saida, "residuos_padronizados_mediana.csv"))
write.csv(data.frame(ano = anos_ajuste, mediana = V_med, li95 = V_li, ls95 = V_ls),
          file.path(saida, "variancia_residual_anual.csv"), row.names = FALSE)

# Resumo simples das previsoes, opcional para a tabela do artigo.
resultado_geral <- data.frame(
  MAE=mean(abs(previsao_csv$observado-previsao_csv$mediana)),
  RMSE=sqrt(mean((previsao_csv$observado-previsao_csv$mediana)^2)),
  Cobertura95=mean(previsao_csv$coberto),
  Largura_media95=mean(previsao_csv$ls95-previsao_csv$li95))
write.csv(resultado_geral,file.path(saida,"avaliacao_previsoes.csv"),row.names=FALSE)

cat("Previsoes e diagnosticos gerados em:",normalizePath(saida),"\n")
cat("Cobertura empirica dos intervalos preditivos 95%:",
    round(resultado_geral$Cobertura95,3),"\n")

############################################################
## ESPERANCA DE VIDA TEMPORARIA ATE 100 ANOS — 2001-2010
## Nao e a esperanca de vida completa; nao extrapola 100+.
## Suposicao: forca constante por faixa = taxa central mx.
############################################################

# duracoes das 21 faixas: 0, 1-4, 5-9, ..., 95-99
larguras <- c(1, 4, rep(5, 19))
stopifnot(length(larguras) == length(idades), sum(larguras) == 100)

# Calculo da integral da sobrevivencia dentro de cada faixa,
# evitando instabilidade numerica quando m esta muito proximo de zero.
vida_temporaria <- function(mx, inicio = 0L) {
  stopifnot(length(mx) == length(larguras), all(is.finite(mx)), all(mx > 0))
  if (inicio == 0L) primeiro <- 1L
  else if (inicio == 60L) primeiro <- match("60-64", idades)
  else stop("Esta rotina foi definida somente para inicio = 0 ou 60.")
  sobrevivencia <- 1
  total <- 0
  for (j in seq.int(primeiro, length(mx))) {
    d <- larguras[j]
    m <- mx[j]
    # Exposicao integrada sob forca constante m.
    total <- total + sobrevivencia * (-expm1(-m * d)) / m
    sobrevivencia <- sobrevivencia * exp(-m * d)
  }
  total
}

# Observado e distribuicao preditiva — e0^(100) e e60^(100)
EV <- array(NA_real_, dim=c(2,H,S),
            dimnames=list(c("e0_ate100","e60_ate100"),as.character(anos_previsao),NULL))
EV_OBS <- matrix(NA_real_,2,H,
                 dimnames=list(c("e0_ate100","e60_ate100"),as.character(anos_previsao)))
for (h in seq_len(H)) {
  EV_OBS[1,h] <- vida_temporaria(MX_TESTE[,h],0L)
  EV_OBS[2,h] <- vida_temporaria(MX_TESTE[,h],60L)
  for (s in seq_len(S)) {
    EV[1,h,s] <- vida_temporaria(MX_PRED[,h,s],0L)
    EV[2,h,s] <- vida_temporaria(MX_PRED[,h,s],60L)
  }
}
saveRDS(EV, file.path(saida,"amostras_preditivas_esperanca_temporaria.rds"))
EV_med <- apply(EV,c(1,2),median)
EV_li <- apply(EV,c(1,2),quantile,probs=.025)
EV_ls <- apply(EV,c(1,2),quantile,probs=.975)
EV_TAB <- expand.grid(indicador=rownames(EV_OBS), ano=anos_previsao,
                      KEEP.OUT.ATTRS=FALSE, stringsAsFactors=FALSE)
EV_TAB$observado <- as.vector(EV_OBS)
EV_TAB$mediana <- as.vector(EV_med)
EV_TAB$li95 <- as.vector(EV_li)
EV_TAB$ls95 <- as.vector(EV_ls)
write.csv(EV_TAB,file.path(saida,"esperanca_vida_temporaria_previsoes.csv"),row.names=FALSE)

pdf(file.path(saida,"Figura_6_Esperanca_Vida_Temporaria_USA.pdf"),
    width=9.2,height=3.45,bg="white")
par(mfrow=c(1,2),mar=c(3.8,4.0,0.75,0.6),mgp=c(2.2,.55,0),
    las=1,bty="o",bg="white",cex.axis=.85,cex.lab=.9)
for (j in 1:2) {
  limites <- range(EV_OBS[j,],EV_li[j,],EV_ls[j,])
  plot(anos_previsao,EV_OBS[j,],type="n",ylim=limites,xlab="ano",
       ylab=if (j==1) "anos esperados ate 100" else "anos esperados de 60 a 100")
  polygon(c(anos_previsao,rev(anos_previsao)),
          c(EV_li[j,],rev(EV_ls[j,])),col="gray85",border=NA)
  lines(anos_previsao,EV_med[j,],lwd=1.5,col="steelblue")
  lines(anos_previsao,EV_OBS[j,],lwd=1.1,col="black",lty=2)
  if (j==1) legend("topleft",c("observado","mediana","intervalo preditivo 95%"),
                  lty=c(2,1,NA),col=c("black","steelblue",NA),
                  fill=c(NA,NA,"gray85"),bty="n",cex=.68)
}
dev.off()
cat("Esperanca de vida temporaria ate 100 anos: tabela e Figura 6 geradas.\n")
