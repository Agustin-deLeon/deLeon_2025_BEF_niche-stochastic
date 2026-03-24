# ============================================================
# FIG (2 panels, 2 columnas) — BEF ~ R2.exp + lluvia_anual
# (mismo estilo que las anteriores: halo + color, puntos con alpha)
# Color = igual al de Figura 3B (col_bef = "#3A8EC1")
# ============================================================

# ---- data ----
# (evito attach; pero si querés mantenerlo, podés comentarlo)
# attach(nn.1)
nn.1<- read.csv("CATS_sin_aleatoria.csv"); nn.1<- nn.1[,-1]; nn.1<-nn.1[,c(1,6)]
read.csv("lluvias.csv")->lluvias
lluvias[-c(1:3),-1]->lluvias
cbind(nn.1, lluvias, prof_media=profundidad[,2], prof_CV=profundidad[,3],
      R2.exp = (1- (R2.unexp)), BEF=pendientes.covariables[,2])-> nn.1
#nn.1$R2.exp <- 1 - nn.1$R2.unexp
as.data.frame(nn.1)-> nn.1
par(mfrow=c(1,2))
bestglm(na.omit(nn.1[,-1]), family = gaussian, IC = "AIC", nvmax=2)$BestModels
m1<-(lm(BEF~R2.exp+lluvia_anual, data=nn.1))
summary(m1)
r.squaredGLMM(m1)

# ---- estética global (ajustá si querés igualar Fig 3) ----
par(mfrow = c(1,2), mar = c(9, 9, 4, 3), mgp = c(3.5, 1, 0),
    cex.lab = 2, cex.axis = 2)

ax_col <- "#2E2E2E"

# ---- color (igual Figura 3B) ----
col_bef <- "#3A8EC1"
pt_bef  <- adjustcolor(col_bef, alpha.f = 0.6)

# ---- helper: doble curva (halo negro abajo, color arriba) ----
double_curve <- function(fun, from, to, col, halo_col = ax_col,
                         halo_lwd = 6, col_lwd = 4.5, ...) {
  curve(fun, from = from, to = to, add = TRUE, col = halo_col, lwd = halo_lwd, ...)
  curve(fun, from = from, to = to, add = TRUE, col = col,      lwd = col_lwd,  ...)
}

# ---- modelo ----
m_bef <- lm(BEF ~ R2.exp + lluvia_anual, data = nn.1)
p <- coefficients(m_bef)

# ============================================================
# Panel A — BEF ~ R2.exp (ajustada por precipitación)
# ============================================================
BEF.neutral <- with(nn.1,
                    BEF - p["lluvia_anual"] * lluvia_anual +
                      mean(p["lluvia_anual"] * lluvia_anual, na.rm = TRUE))

plot(BEF.neutral ~ nn.1$R2.exp,
     xlab = "",
     ylab = "Biomass-Richness slope (BEF)",
     bty = "l",
     pch = 19, col = pt_bef, cex = 3,
     col.axis = ax_col, col.lab = ax_col)

double_curve(
  fun  = function(x) p["(Intercept)"] +
    mean(p["lluvia_anual"] * nn.1$lluvia_anual, na.rm = TRUE) +
    p["R2.exp"] * x,
  from = min(nn.1$R2.exp, na.rm = TRUE),
  to   = max(nn.1$R2.exp, na.rm = TRUE),
  col  = col_bef
)

mtext("Variance in species occurrences\nexplained by traits",
      side = 1, line = 6, cex = 2, col = ax_col)

# stats (tal cual tu texto, con estilo similar)
text(x = 0.42, y = 1.4, "F(2,16): 6.8\nP: 0.007\nR²: 0.46",
     cex = 2.0, pos = 1, col = ax_col)

# ============================================================
# Panel B — BEF ~ precipitación (ajustada por R2.exp)
# ============================================================
BEF.precipitation <- with(nn.1,
                          BEF - p["R2.exp"] * R2.exp +
                            mean(p["R2.exp"] * R2.exp, na.rm = TRUE))

plot(BEF.precipitation ~ nn.1$lluvia_anual,
     xlab = "Annual precipitation",
     ylab = "",
     bty = "l",
     pch = 19, col = pt_bef, cex = 3,
     col.axis = ax_col, col.lab = ax_col)

double_curve(
  fun  = function(x) p["(Intercept)"] +
    p["lluvia_anual"] * x +
    mean(p["R2.exp"] * nn.1$R2.exp, na.rm = TRUE),
  from = min(nn.1$lluvia_anual, na.rm = TRUE),
  to   = max(nn.1$lluvia_anual, na.rm = TRUE),
  col  = col_bef
)






library(ggplot2)
library(ggpubr)
library(ggbeeswarm)
library(reshape2)
library(dplyr)
pendientes_filtrado <- pendientes.covariables[, c(5:18)]
pendientes_largo <- reshape2::melt(pendientes_filtrado, na.rm = TRUE)
pendientes_largo$neutralidad <- NA
pendientes_largo$neutralidad[!is.na(pendientes_largo$value)] <- rep(nn.1$R2.exp, length.out = sum(!is.na(pendientes_largo$value)))
counts <- pendientes_largo %>%
  group_by(variable) %>%
  summarise(n = n())
ggplot(pendientes_largo, aes(x = variable, y = value)) +
  geom_violin(trim = FALSE, color = "black", alpha = 0.5, aes(fill = variable)) +
  geom_quasirandom(aes(color = neutralidad), alpha = 0.7, size = 3) +
  stat_summary(fun = "mean", geom = "point", shape = 18,
               color = "red", size = 3) +
  geom_text(data = counts,
            aes(x = variable, y = -10,
                label = paste("n =", n)),
            color = "black", size = 4) +
  geom_hline(yintercept = 0, color = "black",
             linetype = "dashed", size = 0.8) +
  labs(
    x = " ",
    y = "Estimated effect on biomass",
    color = expression(R^2 ~ "(Traits)")
  ) +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right",
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12)) +
  guides(fill = "none") +
  scale_color_gradient(low = "lightblue", high = "navy")









# ============================================================
# FIG S5 - FRIC, FDIV & RAO ~ R2(Traits)
# ============================================================
layout(matrix(c(1,3,2,4), nrow = 2))
par(mar = c(9, 9, 4, 3),
    mgp = c(3.5, 1, 0),
    cex.lab = 2.2, cex.axis = 1.8)

# ---- R2.exp (si no lo tenés ya creado) ----
resultados_año_completo$R2.exp <- 1 - resultados_año_completo$R2.unexp
Modelo_nulo$R2.exp            <- 1 - Modelo_nulo$R2.unexp

ax_col <- "#2E2E2E"

# ---- color rojizo ----
col_R  <- "#B2472F"
pt_R   <- adjustcolor(col_R, alpha.f = 0.6)

xlab_exp <- "Variance in species occurrences\nexplained by traits"

# ---- helper: doble curva (halo + color) ----
double_curve <- function(fun, from = NULL, to = NULL, col,
                         halo_col = ax_col, halo_lwd = 6, col_lwd = 4.5) {
  if (is.null(from) || is.null(to)) {
    curve(fun(x), add = TRUE, col = halo_col, lwd = halo_lwd)
    curve(fun(x), add = TRUE, col = col,      lwd = col_lwd)
  } else {
    curve(fun(x), from = from, to = to, add = TRUE, col = halo_col, lwd = halo_lwd)
    curve(fun(x), from = from, to = to, add = TRUE, col = col,      lwd = col_lwd)
  }
}

# ============================================================
# ---------------- A: FRic vs Explained-by-traits -------------
# ============================================================
mRic <- lm(data = resultados_año_completo, FRic ~ R2.exp)
p <- coefficients(mRic)

plot(FRic ~ R2.exp, data = resultados_año_completo,
     xlab = " ", ylab = "Metacommunity functional richness",
     bty = "l",
     pch = 19, col = pt_R, cex = 2.5,
     col.axis = ax_col, col.lab = ax_col,
     axes = FALSE)

box(col = ax_col)

double_curve(
  fun = function(x) p[1] + p[2]*x,
  col = col_R
)

mtext(xlab_exp, side = 1, line = 7, cex = 2, col = ax_col)

text(x = 0.425, y = 1.7e-06,
     "F(1,17): 9.3\nP: 0.013\nR²: 0.31",
     cex = 2, pos = 1, col = ax_col)

mtext("A", side = 3, line = 1.2, adj = 0,
      font = 2, cex = 2.2, col = ax_col)

# ============================================================
# ---------------- B: FDiv vs Explained-by-traits -------------
# ============================================================
mDiv <- lm(data = resultados_año_completo, FDiv ~ R2.exp)
p <- coefficients(mDiv)

plot(FDiv ~ R2.exp, data = resultados_año_completo,
     xlab = " ", ylab = "Metacommunity functional divergence",
     bty = "l",
     pch = 19, col = pt_R, cex = 2.5,
     col.axis = ax_col, col.lab = ax_col,
     axes = FALSE)

box(col = ax_col)

double_curve(
  fun = function(x) p[1] + p[2]*x,
  col = col_R
)

mtext(xlab_exp, side = 1, line = 7, cex = 2, col = ax_col)

text(x = 0.75, y = 0.84,
     "F(1,17): 9.7\nP: 0.006\nR²: 0.33",
     cex = 2, pos = 1, col = ax_col)

mtext("B", side = 3, line = 1.2, adj = 0,
      font = 2, cex = 2.2, col = ax_col)


# ============================================================
# ---------------- C: Rao vs Explained-by-traits --------------
# ============================================================
par(new = FALSE)

# Si tu variable observada es RaoQ (como en tu script), dejamos RaoQ:
mRao <- lm(data = resultados_año_completo, RaoQ ~ R2.exp)
p <- coefficients(mRao)

plot(RaoQ ~ R2.exp, data = resultados_año_completo,
     xlab = " ", ylab = "Metacomunity Rao's Q index",
     bty = "l",
     pch = 19, col = pt_R, cex = 2.5,
     col.axis = ax_col, col.lab = ax_col,
     axes = FALSE)

box(col = ax_col)

double_curve(
  fun = function(x) p[1] + p[2]*x,
  col = col_R
)

mtext(xlab_exp, side = 1, line = 7, cex = 2, col = ax_col)

text(x = 0.8, y = 0.055,
     "F(1,17): 19\nP < 0.001\nR²: 0.50",
     cex = 2, pos = 1, col = ax_col)

mtext("C", side = 3, line = 1.2, adj = 0,
      font = 2, cex = 2.2, col = ax_col)

# Panel D vacío
plot.new()

# ---------------- Insets -----

# ---------------- Inset A: z_FRic ----------------
mRic0 <- lm(data = Modelo_nulo, z_FRic ~ R2.exp)
p0 <- coefficients(mRic0)

par(fig = c(0.335, 0.465, 0.75, 0.92), new = TRUE)
par(mar = c(3, 3, 1, 1), mgp = c(2.2, 0.6, 0))

plot(z_FRic ~ R2.exp, data = Modelo_nulo,
     xlab = "", ylab = "", bty = "l",
     pch = 16, col = adjustcolor(col_R, alpha.f = 0.6), cex = 1.2,
     axes = FALSE)

double_curve(
  fun = function(x) p0[1] + p0[2]*x,
  col = col_R,
  halo_lwd = 5, col_lwd = 3.8
)

abline(h = c(-1.96, 1.96), lty = 2, col = ax_col)
abline(h = 0, lty = 1, col = ax_col)

box(col = ax_col)
axis(1, cex.axis = 0.9, tck = -0.04, col.axis = ax_col, col = ax_col)
axis(2, cex.axis = 0.9, tck = -0.04, col.axis = ax_col, col = ax_col)
mtext("Z-values\nfrom null model", side = 2, line = 1.6, cex = 1, col = ax_col)
text(x = 0.4, y = -0.6, paste0("R²: ", round(summary(mRic0)$r.squared, 2)),
     cex = 1, pos = 1, col = ax_col)

# ---------------- Inset B: z_FDiv ----------------
mDiv0 <- lm(data = Modelo_nulo, z_FDiv ~ R2.exp)
p0 <- coefficients(mDiv0)

par(fig = c(0.605, 0.725, 0.6, 0.77), new = TRUE)
par(mar = c(3, 3, 1, 1), mgp = c(2.2, 0.6, 0))

plot(z_FDiv ~ R2.exp, data = Modelo_nulo,
     xlab = "", ylab = "", bty = "l",
     pch = 16, col = adjustcolor(col_R, alpha.f = 0.6), cex = 1.2,
     axes = FALSE)

double_curve(
  fun = function(x) p0[1] + p0[2]*x,
  col = col_R,
  halo_lwd = 5, col_lwd = 3.8
)

abline(h = c(-1.96, 1.96), lty = 2, col = ax_col)
abline(h = 0, lty = 1, col = ax_col)

box(col = ax_col)
axis(1, cex.axis = 0.9, tck = -0.04, col.axis = ax_col, col = ax_col)
axis(2, cex.axis = 0.9, tck = -0.04, col.axis = ax_col, col = ax_col)
mtext("Z-values\nfrom null model", side = 2, line = 1.6, cex = 1, col = ax_col)
text(x = 0.4, y = -0.6, paste0("R²: ", round(summary(mDiv0)$r.squared, 2)),
     cex = 1, pos = 1, col = ax_col)

# ---------------- Inset C: z_RaoQ (o z_Rao) ----------------
mRao0 <- lm(data = Modelo_nulo, z_FRao ~ R2.exp)
p0 <- coefficients(mRao0)

par(fig = c(0.105, 0.235, 0.093, 0.263), new = TRUE)
par(mar = c(3, 3, 1, 1), mgp = c(2.2, 0.6, 0))

plot(Modelo_nulo$z_FRao ~ Modelo_nulo$R2.exp,
     xlab = "", ylab = "", bty = "l",
     pch = 16, col = adjustcolor(col_R, alpha.f = 0.6), cex = 1.2,
     axes = FALSE)

double_curve(
  fun = function(x) p0[1] + p0[2]*x,
  col = col_R,
  halo_lwd = 5, col_lwd = 3.8
)

abline(h = c(-1.96, 1.96), lty = 2, col = ax_col)
abline(h = 0, lty = 1, col = ax_col)

box(col = ax_col)
axis(1, cex.axis = 0.9, tck = -0.04, col.axis = ax_col, col = ax_col)
axis(2, cex.axis = 0.9, tck = -0.04, col.axis = ax_col, col = ax_col)
mtext("Z-values\nfrom null model", side = 2, line = 1.6, cex = 1, col = ax_col)
text(x = 0.4, y = -0.6, paste0("R²: ", round(summary(mRao0)$r.squared, 2)),
     cex = 1, pos = 1, col = ax_col)














plot_glm_visreg_R2 <- function(modelo,
                               xvar = "R2.exp",
                               xlab = "Variance in species occurrences\nexplained by traits",
                               ylab = "P(quadratic)",
                               col_line = "navy",
                               col_fill = adjustcolor("lightblue", 0.25),
                               ax_col = "#2E2E2E",
                               line_xlab = 6.5) {
  
  visreg::visreg(modelo, xvar,
                 scale = "response",
                 xlab = " ",
                 ylab = ylab,
                 cex.lab = 2.5,
                 line.par = list(col = col_line, lwd = 3),
                 fill.par = list(col = col_fill),
  )
  
  mtext(xlab,
        side = 1,
        line = line_xlab,
        cex = 2,
        col = ax_col)
}


plot_lm_partial_R2 <- function(mod, data, yvar,
                               ylab = "", main = "",
                               col_line = "#3A8EC1",
                               pt_col  = adjustcolor(col_line, alpha.f = 0.6),
                               cex_pt = 2.5,
                               ax_col = "#2E2E2E",
                               xlab_mtext = "Variance in species occurrences\nexplained by traits",
                               xlab_line = 7,
                               show_stats = TRUE,
                               stats_pos = "topleft",
                               stats_cex = 2) {
  
  # Residuales parciales (controla por lluvia_anual) si está en el modelo.
  # Si el modelo NO tiene lluvia_anual, plotea directo.
  has_rain <- "lluvia_anual" %in% all.vars(formula(mod))
  
  if (has_rain) {
    y_res <- resid(lm(data[[yvar]] ~ data$lluvia_anual, data = data))
    x_res <- resid(lm(data$R2.exp ~ data$lluvia_anual, data = data))
  } else {
    y_res <- data[[yvar]]
    x_res <- data$R2.exp
  }
  
  ok <- complete.cases(x_res, y_res)
  x_res <- x_res[ok]; y_res <- y_res[ok]
  
  plot(y_res ~ x_res,
       ylab = ylab, xlab = "",
       main = main,
       bty = "l",
       pch = 19, col = pt_col, cex = cex_pt,
       col.axis = ax_col, col.lab = ax_col)
  
  # Línea de regresión en el espacio ploteado
  ab <- coef(lm(y_res ~ x_res))
  abline(a = ab[1], b = ab[2], col = col_line, lwd = 3)
  
  # X label como mtext, estilo tuyo
  mtext(xlab_mtext, side = 1, line = xlab_line, cex = 2, col = ax_col)
  
  # Stats del modelo original
  if (show_stats) {
    sm <- summary(mod)
    if (!is.null(sm$fstatistic)) {
      f  <- sm$fstatistic[1]
      df1 <- sm$fstatistic[2]
      df2 <- sm$fstatistic[3]
      p   <- pf(f, df1, df2, lower.tail = FALSE)
      r2  <- sm$r.squared
      txt <- sprintf("F(%d,%d): %.2f\nP: %.3f\nR²: %.2f", df1, df2, f, p, r2)
      legend(stats_pos, legend = txt, bty = "n", text.col = ax_col, cex = stats_cex)
    }
  }
  
  invisible(list(x = x_res, y = y_res))
}



####################################################
###  ------------ Robustness ------------------   ##
####################################################

# ============================================================
# PANEL 1 (2x2) — FULL DATA
# A: logistic (m_logit)
# B: m_bef_S1 (hybrid der1/linear)   ~ R2.exp + lluvia_anual
# C: m_bef_S2 (hybrid der2/linear)   ~ R2.exp + lluvia_anual
# D: slot (e.g., base slope or ΔAIC vs R2.exp)
# ============================================================

par(mfrow = c(2,2), mar = c(9,8,4,2), cex.lab = 2.2, cex.axis = 1.6)

# A
plot_glm_visreg_R2(
  modelo = m_logit
)
text(x = 0.72, y = 0.9, "P: 0.23", cex = 2,
     pos = 1, col = ax_col)
mtext("A", side = 3, line = 1.0, adj = 0, font = 2, cex = 2.0, col = ax_col)

# B
plot_lm_partial_R2(
  mod = m_bef_S1,
  data = df_all,
  yvar = "dBdS_S1",
  ylab = "Biomass-Richness slope (BEF)\n(derivative S = 1)",
  col_line = col_S1, pt_col = pt_S1, ax_col = ax_col,
  show_stats = F
)
text(x = 0.42, y = 2.5, "F(1,17): 4.7\nP: 0.049\nR²: 0.21", cex = 2,
     pos = 1, col = ax_col)
mtext("B", side = 3, line = 1.0, adj = 0, font = 2, cex = 2.0, col = ax_col)

# C
plot_lm_partial_R2(
  mod = m_bef_S2,
  data = df_all,
  yvar = "dBdS_S2",
  ylab = "Biomass-Richness slope (BEF)\n(derivative S = 2)",
  main = NULL,
  col_line = col_S2, pt_col = pt_S2, ax_col = ax_col,
  show_stats = F
)
text(x = 0.42, y = 2, "F(1,17): 4.7\nP: 0.049\nR²: 0.21", cex = 2,
     pos = 1, col = ax_col)
mtext("C", side = 3, line = 1.0, adj = 0, font = 2, cex = 2.0, col = ax_col)



# ============================================================
# PANEL 2 (2x2) — TRIMMED 1–9 (SAME AESTHETIC AS PANEL 1)
# A: linear slope 1–9          (m_linear_base)  yvar = "linear"
# B: quadratic selection logistic (m_logit_trim9)
# C: hybrid der1 1–9           (m_linear_der1)  yvar = "linear_der1"
# D: hybrid der2 1–9           (m_linear_der2)  yvar = "linear_der2"
# ============================================================

df_panel2 <- if (exists("df_lm")) df_lm else df_trim9

par(mfrow = c(2,2), mar = c(9,8,4,2), cex.lab = 2.2, cex.axis = 1.6)

# ---- A: Linear (1–9) ----
p <- coefficients(m_linear_base)

linear.neutral <- with(df_lm,
                       linear - p["lluvia_anual"] * lluvia_anual +
                         mean(p["lluvia_anual"] * lluvia_anual, na.rm = TRUE))

plot(linear.neutral ~ df_lm$R2.exp,
     ylab = "Biomass-Richness slope (BEF)",
     xlab = "",
     bty = "l", pch = 19, col = pt_trim, cex = 2.5,
     col.axis = ax_col, col.lab = ax_col)

double_curve(
  fun  = function(x) p["(Intercept)"] +
    mean(p["lluvia_anual"] * df_lm$lluvia_anual, na.rm = TRUE) +
    p["R2.exp"] * x,
  from = min(df_lm$R2.exp, na.rm = TRUE),
  to   = max(df_lm$R2.exp, na.rm = TRUE),
  col  = col_trim
)

mtext("Variance in species occurrences\nexplained by traits",
      side = 1, line = 6, cex = 2, col = ax_col)

text(x = 0.35, y = 0.9, "F(2,16): 8.1\nP: 0.004\nR²: 0.5", cex = 2, pos = 1, col = ax_col)

mtext("A", side = 3, line = 1.0, adj = 0, font = 2, cex = 2.0, col = ax_col)

# ---- B: Logistic (1–9) ----
plot_glm_visreg_R2(
  modelo = m_logit_trim9,
  line_xlab = 6.5
)

text(x = 0.8, y = 0.8, "P: 0.81", cex = 2, pos = 1, col = ax_col)

mtext("B", side = 3, line = 1.0, adj = 0, font = 2, cex = 2.0, col = ax_col)


# ---- C: Hybrid der S=1 (1–9) ----
p <- coefficients(m_linear_der1)

plot(df_lm$linear_der1 ~ df_lm$R2.exp,
     ylab = "Richness–biomass relationship (BEF)\n(derivative S = 1)",
     xlab = "",
     bty = "l", pch = 19, col = pt_S1, cex = 2.5,
     col.axis = ax_col, col.lab = ax_col)

double_curve(
  fun  = function(x) p["(Intercept)"] + p["R2.exp"] * x,
  from = min(df_lm$R2.exp, na.rm = TRUE),
  to   = max(df_lm$R2.exp, na.rm = TRUE),
  col  = col_S1
)

mtext("Variance in species occurrences\nexplained by traits",
      side = 1, line = 6, cex = 2, col = ax_col)

text(x = 0.4, y = 1, "F(1,17): 3.7\nP: 0.07\nR²: 0.2", cex = 2, pos = 1, col = ax_col)

mtext("C", side = 3, line = 1.0, adj = 0, font = 2, cex = 2.0, col = ax_col)


# ---- D: Hybrid der S=2 (1–9) ----
p <- coefficients(m_linear_der2)

der2.neutral <- with(df_lm,
                     linear_der2 - p["lluvia_anual"] * lluvia_anual +
                       mean(p["lluvia_anual"] * lluvia_anual, na.rm = TRUE))

plot(der2.neutral ~ df_lm$R2.exp,
     ylab = "Biomass-Richness slope (BEF)\n(derivative S = 2)",
     xlab = "",
     bty = "l", pch = 19, col = pt_S2, cex = 2.5,
     col.axis = ax_col, col.lab = ax_col)

double_curve(
  fun  = function(x) p["(Intercept)"] +
    mean(p["lluvia_anual"] * df_lm$lluvia_anual, na.rm = TRUE) +
    p["R2.exp"] * x,
  from = min(df_lm$R2.exp, na.rm = TRUE),
  to   = max(df_lm$R2.exp, na.rm = TRUE),
  col  = col_S2
)

mtext("Variance in species occurrences\nexplained by traits",
      side = 1, line = 6, cex = 2, col = ax_col)

text(x = 0.4, y = 1, "F(2,16): 4\nP: 0.04\nR²: 0.33", cex = 2, pos = 1, col = ax_col)

mtext("D", side = 3, line = 1.0, adj = 0, font = 2, cex = 2.0, col = ax_col)



























library(glmmTMB)
library(MuMIn)
library(viridisLite)

par(mfrow=c(5,4), mar=c(5.6,6,3,2), cex.lab = 2)

slopes_linear <- list()
slopes_ponds  <- list()

for (ii in 2005:2023) {
  
  bm.temp <- bm.nuevo[bm.nuevo[, 1] == ii, ]
  
  df.um.ch <- bm.temp[, c(1, 3, 6, 5, 7:10, 14:26)]
  colnames(df.um.ch) <- c(
    "um.year", "pond.id", "um.biom", "um.rich",
    "rain_sampling", "rain_annual", "temp_sampling", "temp_annual",
    "DM", "ddmm", "Shape", "Islands", "log.Area", "log.Volume",
    "Mean.Depth", "Sd.Depth", "CV.Depth", "Hydroperiod", "Degree",
    "log.Betweenness", "Closenness"
  )
  
  # >>> NEW: drop richness == 0 (or lower) <<<
  df.um.ch <- df.um.ch[df.um.ch$um.rich > 0, ]
  
  df.um.ch <- na.omit(df.um.ch)
  df.um.ch$pond.id <- as.factor(df.um.ch$pond.id)
  
  title_txt <- paste("Year: ", ii)
  
  # Skip if too little data / only one pond
  if (nrow(df.um.ch) < 5 || length(unique(df.um.ch$pond.id)) < 2) {
    plot.new()
    title(main = title_txt)
    mtext("Insufficient data", side = 3, line = -2)
    next
  }
  
  m1 <- glmmTMB(
    um.biom ~ um.rich + (um.rich | pond.id),
    family = gaussian(),
    data = df.um.ch
  )
  
  fixef_m1 <- fixef(m1)$cond
  ranef_intercepts <- ranef(m1)$cond$pond.id[, "(Intercept)"]
  ranef_slopes <- ranef(m1)$cond$pond.id[, "um.rich"]
  
  slopes_ponds[[as.character(ii)]] <- ranef_slopes
  
  unique_ids <- levels(df.um.ch$pond.id)
  
  colors <- viridis(length(unique_ids), option = "D")
  names(colors) <- unique_ids
  
  plot(
    df.um.ch$um.rich, df.um.ch$um.biom,
    col = colors[as.character(df.um.ch$pond.id)],
    pch = 19,
    xlab = "Species richness",
    ylab = expression("Biomass (g / 0.04 m"^2*")"),
    main = title_txt,
    cex.lab = 2.2,
    cex.axis = 2,
    cex.main = 1.7
  )
  
  # ---- Pond-specific dashed lines (only within each pond's richness range) ----
  for (i in seq_along(unique_ids)) {
    
    pond_i <- unique_ids[i]
    current_data <- df.um.ch[df.um.ch$pond.id == pond_i, ]
    
    intercept <- fixef_m1["(Intercept)"] + ranef_intercepts[i]
    slope     <- fixef_m1["um.rich"]     + ranef_slopes[i]
    
    x_min <- min(current_data$um.rich)
    x_max <- max(current_data$um.rich)
    
    # draw dashed line only within that pond's richness range
    x_vals <- seq(x_min, x_max, length.out = 50)
    y_vals <- intercept + slope * x_vals
    
    lines(
      x_vals, y_vals,
      col = colors[pond_i],
      lwd = 1,
      lty = 2
    )
  }
  
  # ---- Overall fixed-effects solid line ----
  general_intercept <- fixef_m1["(Intercept)"]
  general_slope     <- fixef_m1["um.rich"]
  
  curve(
    general_intercept + general_slope * x,
    from = min(df.um.ch$um.rich, na.rm = TRUE),
    to   = max(df.um.ch$um.rich, na.rm = TRUE),
    add = TRUE,
    col = "black",
    lwd = 3,
    lty = 1
  )
  
  # ---- Save metrics ----
  s <- summary(m1)
  sigma.temp <- s$sigma^2
  sd_int <- attr(s$varcor$cond$pond.id, "stddev")["(Intercept)"]
  sd_slp <- attr(s$varcor$cond$pond.id, "stddev")["um.rich"]
  r2_vals <- r.squaredGLMM(m1)
  
  row_out <- cbind(
    t(as.matrix(fixef_m1)),
    sigma.temp,
    r2_vals[1],
    r2_vals[2],
    sd_int,
    sd_slp,
    s$coefficients$cond["um.rich", "Pr(>|z|)"]
  )
  
  colnames(row_out) <- c(
    "intercept", "slope", "sigma", "R2_m", "R2_c",
    "sd_intercept", "sd_slope", "p_value"
  )
  rownames(row_out) <- ii
  
  slopes_linear[[as.character(ii)]] <- row_out
}

slopes_linear <- do.call(rbind, slopes_linear)
pendientes.lineal <- do.call(rbind, pendientes.lineal)







library(glmmTMB)
library(MuMIn)
library(viridisLite)

par(mfrow=c(5,4), mar=c(5.5,6,3,2), cex.lab = 2)

slopes <- list()
sigma_vec  <- numeric()

for (ii in 2005:2023) {
  
  # ---- 1) Filter year ----
  bm.temp <- bm.nuevo[bm.nuevo[, 1] == ii, ]
  
  # ---- 2) Build df ----
  df.um.ch <- bm.temp[, c(1, 3, 6, 5, 7:10, 14:26)]
  colnames(df.um.ch) <- c(
    "um.year", "pond.id", "um.biom", "um.rich",
    "rain_sampling", "rain_annual", "temp_sampling", "temp_annual",
    "DM", "ddmm", "Shape", "Islands", "log.Area", "log.Volume",
    "Mean.Depth", "Sd.Depth", "CV.Depth", "Hydroperiod", "Degree",
    "log.Betweenness", "Closenness"
  )
  
  df.um.ch <- na.omit(df.um.ch)
  df.um.ch$pond.id <- as.factor(df.um.ch$pond.id)
  
  title_txt <- paste("Year:", ii)
  
  # If too little data, skip nicely
  if (nrow(df.um.ch) < 5 || length(unique(df.um.ch$pond.id)) < 2) {
    plot.new()
    title(main = title_txt)
    mtext("Insufficient data", side = 3, line = -2)
    next
  }
  
  # ---- 3) Fit candidate models ----
  m1  <- try(
    glmmTMB(um.biom ~ um.rich + (1 + um.rich | pond.id),
            family = gaussian(), data = df.um.ch),
    silent = TRUE
  )
  
  m12 <- try(
    glmmTMB(um.biom ~ um.rich + I(um.rich^2) + (1 + um.rich | pond.id),
            family = gaussian(), data = df.um.ch),
    silent = TRUE
  )
  
  # --- Individual AICs ---
  AIC_m1  <- if (!inherits(m1, "try-error"))  tryCatch(as.numeric(AIC(m1)),  error = function(e) NA_real_) else NA_real_
  AIC_m12 <- if (!inherits(m12, "try-error")) tryCatch(as.numeric(AIC(m12)), error = function(e) NA_real_) else NA_real_
  
  # Delta AIC (positive => quadratic better)
  deltaAIC_lin_minus_quad <- if (!is.na(AIC_m1) && !is.na(AIC_m12)) (AIC_m1 - AIC_m12) else NA_real_
  
  # ---- 4) Selection rule: keep quadratic only if ΔAIC > 2 ----
  best_model <- NULL
  best_name  <- NA_character_
  best_label <- NA_character_
  
  if (!inherits(m1, "try-error") && !inherits(m12, "try-error")) {
    if (!is.na(deltaAIC_lin_minus_quad) && deltaAIC_lin_minus_quad > 2) {
      best_model <- m12
      best_name  <- "m12"
      best_label <- "quadratic"
    } else {
      best_model <- m1
      best_name  <- "m1"
      best_label <- "linear"
    }
  } else if (!inherits(m1, "try-error")) {
    best_model <- m1
    best_name  <- "m1"
    best_label <- "linear"
  } else if (!inherits(m12, "try-error")) {
    best_model <- m12
    best_name  <- "m12"
    best_label <- "quadratic"
  } else {
    plot.new()
    title(main = title_txt)
    mtext("No model converged", side = 3, line = -2)
    next
  }
  
  s <- summary(best_model)
  
  # ---- 5) Coefs and random effects ----
  fixef_best <- fixef(best_model)$cond
  ran_int <- ranef(best_model)$cond$pond.id[, "(Intercept)"]
  ran_slp <- ranef(best_model)$cond$pond.id[, "um.rich"]
  
  unique_ids <- levels(df.um.ch$pond.id)
  
  colors <- viridis(length(unique_ids), option = "D")
  names(colors) <- unique_ids
  
  # ---- 6) Plot ----
  plot(df.um.ch$um.rich, df.um.ch$um.biom,
       col = colors[as.character(df.um.ch$pond.id)],
       pch = 19,
       xlab = "Species richness",
       ylab = expression("Biomass (g / 0.04 m"^2*")"),
       main = paste0(title_txt, " (", best_label, ")"))
  
  # Pond-specific lines
  for (k in seq_along(unique_ids)) {
    idk <- unique_ids[k]
    current_data <- df.um.ch[df.um.ch$pond.id == idk, ]
    
    x_min <- min(current_data$um.rich, na.rm = TRUE)
    x_max <- max(current_data$um.rich, na.rm = TRUE)
    
    intercept <- fixef_best["(Intercept)"] + ran_int[idk]
    slope1    <- fixef_best["um.rich"]     + ran_slp[idk]
    slope2    <- if ("I(um.rich^2)" %in% names(fixef_best)) fixef_best["I(um.rich^2)"] else 0
    
    curve(intercept + slope1 * x + slope2 * x^2,
          from = x_min, to = x_max, add = TRUE,
          col = colors[idk], lwd = 1, lty = 2)
  }
  
  # Overall fixed-effects line
  g0 <- fixef_best["(Intercept)"]
  g1 <- fixef_best["um.rich"]
  g2 <- if ("I(um.rich^2)" %in% names(fixef_best)) fixef_best["I(um.rich^2)"] else 0
  
  curve(g0 + g1 * x + g2 * x^2,
        from = min(df.um.ch$um.rich, na.rm = TRUE),
        to   = max(df.um.ch$um.rich, na.rm = TRUE),
        add = TRUE, col = "black", lwd = 3)
  
  # ---- 7) Save metrics ----
  p_umrich <- NA_real_
  if ("um.rich" %in% rownames(s$coefficients$cond)) {
    p_umrich <- s$coefficients$cond["um.rich", "Pr(>|z|)"]
  }
  
  sigma2 <- s$sigma^2
  sigma_vec <- c(sigma_vec, sigma2)
  
  r2m <- r2c <- NA_real_
  r2_try <- try(r.squaredGLMM(best_model), silent = TRUE)
  if (!inherits(r2_try, "try-error")) {
    r2m <- as.numeric(r2_try[1])
    r2c <- as.numeric(r2_try[2])
  }
  
  sd_int <- sd_slp <- NA_real_
  vc <- try(summary(best_model)$varcor$cond$pond.id, silent = TRUE)
  if (!inherits(vc, "try-error")) {
    sd_int <- attr(vc, "stddev")["(Intercept)"]
    sd_slp <- attr(vc, "stddev")["um.rich"]
  }
  
  row_out <- data.frame(
    year = ii,
    model = best_name,
    best_label = best_label,
    AIC_m1 = AIC_m1,
    AIC_m12 = AIC_m12,
    deltaAIC_lin_minus_quad = deltaAIC_lin_minus_quad, # >2 => quadratic
    intercept = unname(g0),
    slope     = unname(g1),
    quadratic = unname(g2),
    sigma     = sigma2,
    R2_m      = r2m,
    R2_c      = r2c,
    p_value   = p_umrich,
    stddev_int  = sd_int,
    stddev_slope = sd_slp
  )
  
  slopes[[as.character(ii)]] <- row_out
}

slopes_df <- do.call(rbind, slopes)





###############################
###############################
# Con covariables

find.x <- function(m, y_en, x_en, KK = 3) {
  out <- NULL
  require(bestglm)
  require(mgcv)
  
  X <- m[, x_en]
  y <- m[, y_en]
  Xy <- as.data.frame(cbind(X, y))
  
  # Model selection by AIC
  model <- bestglm(Xy, family = gaussian, IC = "AIC")
  model$BestModels -> models_tbl
  names(model$BestModel$coefficients)[-1] -> vars_best
  
  a.t <- NULL
  for (i in 1:nrow(models_tbl)) {
    vars_temp <- colnames(models_tbl)[which(models_tbl[i, ] == TRUE)]
    form <- as.formula(paste("y", "~", paste(vars_temp, collapse = "+")))
    fit <- lm(form, data = Xy)
    
    aa <- summary(fit)$coefficients
    a <- if ("um.rich" %in% rownames(aa)) aa["um.rich", 1] else NA_real_
    a.t <- rbind(a.t, cbind(row = i, a))
  }
  
  print(a.t)
  cbind(models_tbl, a = a.t[, 2]) -> models_tbl
  models_tbl[which(models_tbl$Criterion == median(models_tbl$Criterion)), ] -> temp
  
  vars_model <- colnames(temp)[which(temp[1, ] == TRUE)]
  out[[1]] <- vars_model
  out[[2]] <- models_tbl
  out[[3]] <- vars_best
  return(out)
}

library(visreg)

slopes_covariates <- NULL
par(mfrow = c(5, 4), mar = c(5.5, 6, 3, 2), cex.lab = 2)

for (ii in 2005:2023) {
  
  slopes_cov_temp <- as.data.frame(matrix(NA, ncol = 18, nrow = 1))
  colnames(slopes_cov_temp) <- c(
    "Intercept", "um.rich", "p.value", "R2", "Islands", "log.Area",
    "Mean.Depth", "CV.Depth", "Hydroperiod", "Degree", "log.Betweenness",
    "um.rich_Islands", "um.rich_log.Area",
    "um.rich_Mean.Depth", "um.rich_CV.Depth", "um.rich_Hydroperiod",
    "um.rich_Degree", "um.rich_log.Betweenness"
  )
  
  bm.temp <- bm.nuevo[bm.nuevo[, 1] == ii, ]
  
  df.um.ch <- bm.temp[, c(1, 3, 6, 5, 7:10, 14:16, 29, 34, 30:32, 23:25)]
  colnames(df.um.ch) <- c(
    "um.year", "pond.id", "um.biom", "um.rich", "rain_sampling", "rain_annual",
    "temp_sampling", "temp_annual", "DM", "ddmm", "Shape", "Islands", "log.Area",
    "Mean.Depth", "Sd.Depth", "CV.Depth", "Hydroperiod", "Degree", "log.Betweenness"
  )
  
  # >>> NEW: drop records with richness == 0 (or lower) <<<
  df.um.ch <- df.um.ch[df.um.ch$um.rich > 0, ]
  
  df.um.ch <- na.omit(df.um.ch)
  df.um.ch$pond.id <- as.factor(df.um.ch$pond.id)
  
  title_txt <- paste("Year:", ii)
  
  # If too little data, skip cleanly
  if (nrow(df.um.ch) < 5 || length(unique(df.um.ch$pond.id)) < 2) {
    plot.new()
    title(main = title_txt)
    mtext("Insufficient data", side = 3, line = -2)
    next
  }
  
  # Save richness mean/sd (for unscaled x-axis labels)
  rich_mean <- mean(df.um.ch$um.rich, na.rm = TRUE)
  rich_sd   <- sd(df.um.ch$um.rich,  na.rm = TRUE)
  
  # Scale richness + covariates (same as your original workflow)
  df.um.ch[, c(4, 12:19)] <- scale(df.um.ch[, c(4, 12:19)])
  
  # Interactions (scaled richness * scaled covariates)
  df.um.ch[, 20:27] <- df.um.ch[, 4] * df.um.ch[, c(12:19)]
  colnames(df.um.ch)[20:27] <- paste0("um.rich_", colnames(df.um.ch)[c(12:19)])
  
  # Covariate selection
  variables <- find.x(m = df.um.ch, y_en = 3, x_en = c(4, 12:27))
  
  # Defaults to avoid object-not-found edge cases
  variables_model2 <- character(0)
  variables_model  <- character(0)
  coefs <- matrix(NA_real_, nrow = 1, ncol = 4,
                  dimnames = list("(Intercept)", c("Estimate","Std. Error","t value","Pr(>|t|)")))
  
  if (ii == 2021) {
    
    variables_model2 <- variables[[1]]
    form <- as.formula(paste("um.biom", "~", paste(variables_model2, collapse = "+")))
    fit <- lm(form, data = df.um.ch)
    
    resid <- residuals(fit)
    df2021 <- cbind(resid, df.um.ch)
    fit_resid_2021 <- lm(resid ~ um.rich, data = df2021)
    
    visreg(fit_resid_2021, "um.rich",
           xlab = "",
           ylab = expression("Biomass (g / 0.04 m"^2*")"),
           main = "Year: 2021*",
           line.par = list(col = "darkgreen", lwd = 2),
           fill.par = list(col = rgb(0.5, 0.5, 1, 0.3)),
           xaxt = "n")
    
    ticks <- axTicks(1)
    axis(1, at = ticks, labels = round(ticks * rich_sd + rich_mean, 0))
    mtext("Species richness", side = 1, line = 3)
    
    summary(fit_resid_2021)
    mtext("Estimate: 0.01, p-value: 0.93", side = 3, line = 0, cex = 0.8, col = "black")
    
    slopes_cov_temp$um.rich  <- 0.01
    slopes_cov_temp$p.value  <- 0.93
    coefs <- summary(fit)$coefficients
    
  } else {
    
    if (length(variables[[1]]) == 0) {
      plot.new()
      title(main = title_txt)
      mtext("Empty selection vector (no model)", side = 3, line = -2, cex = 0.8)
      coefs <- matrix(NA_real_, nrow = 1, ncol = 4,
                      dimnames = list("(Intercept)", c("Estimate","Std. Error","t value","Pr(>|t|)")))
    } else {
      
      variables_model2 <- variables[[1]]
      variables_model  <- variables[[3]]
      
      form <- as.formula(paste("um.biom", "~", paste(variables_model2, collapse = "+")))
      fit <- lm(form, data = df.um.ch)
      coefs <- summary(fit)$coefficients
      
      # Robust check (prevents subscript out of bounds)
      if ("um.rich" %in% rownames(coefs)) {
        
        p_val <- coefs["um.rich", "Pr(>|t|)"]
        est   <- coefs["um.rich", "Estimate"]
        
        visreg(fit, "um.rich",
               xlab = "",
               ylab = expression("Biomass (g / 0.04 m"^2*")"),
               main = title_txt,
               line.par = list(col = "darkgreen", lwd = 2),
               fill.par = list(col = rgb(0.5, 0.5, 1, 0.3)),
               xaxt = "n")
        
        ticks <- axTicks(1)
        axis(1, at = ticks, labels = round(ticks * rich_sd + rich_mean, 0))
        mtext("Species richness", side = 1, line = 3)
        
        txt <- ifelse(p_val > 0.001,
                      sprintf("p-value: %.3f, Estimate: %.3f", p_val, est),
                      sprintf("p-value < 0.001, Estimate: %.3f", est))
        mtext(txt, side = 3, line = 0, cex = 0.8, col = "black")
        
        slopes_cov_temp$um.rich  <- est
        slopes_cov_temp$p.value  <- p_val
        
      } else {
        plot.new()
        title(main = title_txt)
        mtext("um.rich was not retained in the model", side = 3, line = -2, cex = 0.8)
      }
    }
  }
  
  # Intercept
  if ("(Intercept)" %in% rownames(coefs)) {
    slopes_cov_temp$Intercept <- coefs["(Intercept)", 1]
  }
  
  covariates <- c(
    "Islands", "log.Area", "log.Volumen", "Mean.Depth", "CV.Depth", "Hydroperiod",
    "Degree", "log.Betweenness", "Closenness", "um.rich_Islands", "um.rich_log.Area",
    "um.rich_log.Volumen", "um.rich_Mean.Depth", "um.rich_CV.Depth", "um.rich_Hydroperiod",
    "um.rich_Degree", "um.rich_log.Betweenness", "um.rich_Closenness"
  )
  
  for (v in covariates) {
    if ((v %in% variables_model2) && (v %in% rownames(coefs))) {
      slopes_cov_temp[[v]] <- coefs[v, 1]
    }
  }
  
  # R2 (MuMIn). If it fails, keep NA.
  r2_try <- try(r.squaredGLMM(fit)[1], silent = TRUE)
  if (!inherits(r2_try, "try-error")) slopes_cov_temp$R2 <- as.numeric(r2_try)
  
  slopes_covariates <- rbind(slopes_covariates, slopes_cov_temp)
}

pendientes.covariables<-slopes_covariates
