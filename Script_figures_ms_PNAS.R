# ============================================================
# (NUEVO) Crear la variable "explicada por rasgos" una sola vez
# ============================================================
nn.2$R2.exp <- 1 - nn.2$R2.unexp
df_model$R2.exp <- 1 - df_model$R2_unexp
resultados_año_completo$R2.exp <- 1 - resultados_año_completo$R2.unexp
Modelo_nulo$R2.exp <- 1 - Modelo_nulo$R2.unexp
biom_neutral$R2.exp <- 1 - biom_neutral$R2.unexp
nn.1$R2.exp <- 1 - nn.1$R2.unexp
R2.exp<- 1 - (R2.unexp)



# ============================================================
# FIGURA 1 —
# ============================================================
layout(matrix(1:2, nrow = 1))
par(mar = c(5, 8, 4, 3), mgp = c(3.5, 1, 0), cex.lab = 2.2, cex.axis = 1.8)

ax_col <- "#2E2E2E"

# Colores Fig 1 (los tuyos)
col_drift <- "darkgreen"
col_bef   <- "darkblue"

pt_drift <- adjustcolor(col_drift, alpha.f = 0.6)
pt_bef   <- adjustcolor(col_bef,   alpha.f = 0.6)

# ---- helper para línea doble (halo + color) ----
double_line <- function(x, y, col, halo_col = ax_col, halo_lwd = 4, col_lwd = 2.5) {
  lines(x, y, col = halo_col, lwd = halo_lwd)
  lines(x, y, col = col,      lwd = col_lwd)
}

# ===================== Plot A =====================
plot(R2.exp ~ year, data = nn.2,
     xlab = "Year",
     ylab = "Variance in species occurrences\nexplained by traits",
     pch = 19, col = pt_drift, cex = 2,
     bty = "l",
     type = "n",
     col.axis = ax_col, col.lab = ax_col)

# línea doble + puntos
double_line(nn.2$year, nn.2$R2.exp, col = col_drift)
points(nn.2$year, nn.2$R2.exp, pch = 19, col = pt_drift, cex = 2.5)

mtext("A", side = 3, line = 1.2, adj = 0, font = 2, cex = 2, col = ax_col)

# ===================== Plot B =====================
plot(pendientes.lineal[,2] ~ nn.2$year,
     xlab = "Year",
     ylab = "Biomass-Richness slope (BEF)",
     pch = 19, col = pt_bef, cex = 2,
     bty = "l",
     type = "n",
     col.axis = ax_col, col.lab = ax_col)

double_line(nn.2$year, pendientes.lineal[,2], col = col_bef)
points(nn.2$year, pendientes.lineal[,2], pch = 19, col = pt_bef, cex = 2.5)

mtext("B", side = 3, line = 1.2, adj = 0, font = 2, cex = 2, col = ax_col)









# ============================================================
# FIGURA 2 —
# ============================================================

layout(matrix(c(1,3,2,4), nrow = 2))
par(mar = c(9, 9, 4, 3),
    mgp = c(3.5, 1, 0),
    cex.lab = 2.2, cex.axis = 1.8)

# ---- Paleta por panel ----
col_A <- "#2E5EAA"  # azul
col_B <- "#B2472F"  # terracota
col_C <- "#1F8A70"  # teal
col_D <- "#6A78A8"  # azul violáceo

# puntos = color del panel con alpha 0.6
pt_A <- adjustcolor(col_A, alpha.f = 0.6)
pt_B <- adjustcolor(col_B, alpha.f = 0.6)
pt_C <- adjustcolor(col_C, alpha.f = 0.6)
pt_D <- adjustcolor(col_D, alpha.f = 0.6)

ax_col <- "#2E2E2E"

# ---- helper: doble curva (halo negro abajo, color arriba) ----
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
# ---------------- A: Explained-by-traits vs delta lluvia ------
# ============================================================
library(quantreg)
m1_mediana <- rq(R2.exp ~ delta_lluvia_anual + I(delta_lluvia_anual^2),
                 data = df_model, tau = 0.5)

plot(R2.exp ~ delta_lluvia_anual, data = df_model,
     xlab = "",
     ylab = "Variance in species occurrences\nexplained by traits",
     pch = 19, col = pt_A, cex = 2.5,
     bty = "l",
     col.axis = ax_col, col.lab = ax_col)

double_curve(
  fun = function(x)
    coef(m1_mediana)[1] +
    coef(m1_mediana)[2]*x +
    coef(m1_mediana)[3]*x^2,
  from = min(df_model$delta_lluvia_anual, na.rm = TRUE),
  to   = max(df_model$delta_lluvia_anual, na.rm = TRUE),
  col  = col_A
)

text(x = -75, y = 0.35, "P: 0.05", cex = 2, pos = 4, col = ax_col)
mtext("Inter-annual difference in\nprecipitation", side = 1, line = 7,
      cex = 2, col = ax_col)
mtext("A", side = 3, line = 1.2, adj = 0, font = 2, cex = 2.2, col = ax_col)

# ============================================================
# ---------------- B: Riqueza vs Explained-by-traits ----------
# ============================================================
p <- coefficients(lm(round(riqueza) ~ R2.exp))

plot(round(riqueza) ~ R2.exp,
     ylab = "Metacommunity species richness",
     xlab = "",
     bty = "l",
     pch = 19, col = pt_B, cex = 2.5,
     col.axis = ax_col, col.lab = ax_col)

double_curve(
  fun = function(x) p[1] + p[2]*x,
  col = col_B
)

text(x = 0.8, y = 100, "F(2,17): 4.3\nP: 0.05\nR²: 0.16",
     cex = 2, pos = 1, col = ax_col)

mtext("Variance in species occurrences\nexplained by traits",
      side = 1, line = 7, cex = 2, col = ax_col)
mtext("B", side = 3, line = 1.2, adj = 0, font = 2, cex = 2.2, col = ax_col)

# ============================================================
# ---------------- C: FDis vs Explained-by-traits -------------
# ============================================================
mDis <- lm(data = resultados_año_completo, FDis ~ R2.exp)
p <- coefficients(mDis)

par(mgp = c(3.2, 1, 0))
plot(FDis ~ R2.exp, data = resultados_año_completo,
     xlab = " ", ylab = "Functional dispersion",
     bty = "l",
     pch = 19, col = pt_C, cex = 2.5,
     col.axis = ax_col, col.lab = ax_col)

double_curve(
  fun = function(x) p[1] + p[2]*x,
  col = col_C
)

mtext("Variance in species occurrences\nexplained by traits",
      side = 1, line = 7, cex = 2, col = ax_col)

text(x = 0.4, y = 0.203, "F(1,17): 21\nP < 0.001\nR²: 0.53",
     cex = 2, pos = 1, col = ax_col)

mtext("C", side = 3, line = 1.2, adj = 0, font = 2, cex = 2.2, col = ax_col)

# ============================================================
# ---------------- D: FEve vs Explained-by-traits -------------
# ============================================================
mEve <- lm(data = resultados_año_completo, FEve ~ R2.exp)
p <- coefficients(mEve)

plot(FEve ~ R2.exp, data = resultados_año_completo,
     xlab = " ", ylab = "Functional evenness",
     bty = "l",
     pch = 19, col = pt_D, cex = 2.5,
     col.axis = ax_col, col.lab = ax_col)

double_curve(
  fun = function(x) p[1] + p[2]*x,
  col = col_D
)

text(x = 0.82, y = 0.505, "F(1,17): 0.07\nP: 0.8\nR² < 0.01",
     cex = 2, pos = 1, col = ax_col)

mtext("Variance in species occurrences\nexplained by traits",
      side = 1, line = 7, cex = 2, col = ax_col)
mtext("D", side = 3, line = 1.2, adj = 0, font = 2, cex = 2.2, col = ax_col)

# ============================================================
# ---------------- Inset C: Z_FDis vs Explained-by-traits -----
# ============================================================
mDis0 <- lm(data = Modelo_nulo, z_FDis ~ R2.exp)
p0 <- coefficients(mDis0)

par(fig = c(0.355, 0.475, 0.32, 0.47), new = TRUE)
par(mar = c(3, 3, 1, 1), mgp = c(2.2, 0.6, 0))

plot(z_FDis ~ R2.exp, data = Modelo_nulo,
     xlab = "", ylab = "", bty = "l",
     pch = 16, col = adjustcolor(col_C, alpha.f = 0.6), cex = 1.2,
     axes = FALSE)

double_curve(
  fun = function(x) p0[1] + p0[2]*x,
  col = col_C,
  halo_lwd = 5, col_lwd = 3.8
)

abline(h = c(-1.96, 1.96), lty = 2, col = ax_col)
abline(h = 0, lty = 1, col = ax_col)

box(col = ax_col)
axis(1, cex.axis = 0.9, tck = -0.04, col.axis = ax_col, col = ax_col)
axis(2, cex.axis = 0.9, tck = -0.04, col.axis = ax_col, col = ax_col)
mtext("Z-values\nfrom null model", side = 2, line = 1.6, cex = 1, col = ax_col)
text(x = 0.4, y = -0.6, "R²: 0.18", cex = 1, pos = 1, col = ax_col)

# ============================================================
# ---------------- Inset D: Z_FEve vs Explained-by-traits -----
# ============================================================
mEve0 <- lm(data = Modelo_nulo, z_FEve ~ R2.exp)
p0 <- coefficients(mEve0)

par(fig = c(0.855, 0.975, 0.28, 0.43), new = TRUE)
par(mar = c(3, 3, 1, 1), mgp = c(2.2, 0.6, 0))

plot(z_FEve ~ R2.exp, data = Modelo_nulo,
     xlab = "", ylab = "", bty = "l",
     pch = 16, col = adjustcolor(col_D, alpha.f = 0.6), cex = 1.2,
     axes = FALSE)

double_curve(
  fun = function(x) p0[1] + p0[2]*x,
  col = col_D,
  halo_lwd = 5, col_lwd = 3.8
)

abline(h = c(-1.96, 1.96), lty = 2, col = ax_col)
abline(h = 0, lty = 1, col = ax_col)

box(col = ax_col)
axis(1, cex.axis = 0.9, tck = -0.04, col.axis = ax_col, col = ax_col)
axis(2, cex.axis = 0.9, tck = -0.04, col.axis = ax_col, col = ax_col)
mtext("Z-values\nfrom null model", side = 2, line = 1.6, cex = 1, col = ax_col)
text(x = 0.4, y = -2.1, "R²: 0.21", cex = 1, pos = 1, col = ax_col)










# ============================================================
# Figura 3 — 
# ============================================================

# ---- Paleta (más viva) ----
col_bef <- "#3A8EC1"
col_biom  <- "#E69F00"
pt_bef <- adjustcolor(col_bef, alpha.f = 0.6)
pt_biom  <- adjustcolor(col_biom,  alpha.f = 0.6)
ax_col <- "#2E2E2E"

double_curve <- function(fun, from, to, col, halo_col = ax_col,
                         halo_lwd = 6, col_lwd = 4.5, ...) {
  curve(fun, from = from, to = to, add = TRUE, col = halo_col, lwd = halo_lwd, ...)
  curve(fun, from = from, to = to, add = TRUE, col = col,      lwd = col_lwd,  ...)
}

layout(matrix(1:4, nrow = 2, byrow = TRUE))
par(mar = c(9, 9, 4, 3), mgp = c(3.5, 1, 0), cex.lab = 2.2, cex.axis = 1.8)

# ============================================================
# A) Biomass model: biom ~ R2.exp + lluvia_anual
# ============================================================

# (mejor que attach)
# attach(biom_neutral)

m_biom <- lm(biom ~ R2.exp + lluvia_anual, data = biom_neutral)
p <- coefficients(m_biom)

# ---------- A: Biomass ~ explained-by-traits (ajustada por precip) ----------
Biom.neutral <- with(biom_neutral,
                     biom - p["lluvia_anual"] * lluvia_anual + mean(p["lluvia_anual"] * lluvia_anual, na.rm = TRUE))

plot(Biom.neutral ~ biom_neutral$R2.exp,
     ylab = expression("Mean biomass (g / 0.04 m"^2*")"), 
     xlab = "",
     bty = "l", pch = 19, col = pt_biom, cex = 2.5,
     col.axis = ax_col, col.lab = ax_col)

double_curve(
  fun  = function(x) p["(Intercept)"] + mean(p["lluvia_anual"] * biom_neutral$lluvia_anual, na.rm = TRUE) + p["R2.exp"] * x,
  from = min(biom_neutral$R2.exp, na.rm = TRUE),
  to   = max(biom_neutral$R2.exp, na.rm = TRUE),
  col  = col_biom
)

mtext("Variance in species occurrences\nexplained by traits", side = 1, line = 6, cex = 2, col = ax_col)
mtext("A", side = 3, line = 1.5, adj = 0, font = 2, cex = 2.5, col = ax_col)

text(x = 0.42, y = 8, "F(2,16): 9.8\nP: 0.002\nR²: 0.50", cex = 2,
     pos = 1, col = ax_col)


# ---------- Biomass ~ precip (ajustada por explained-by-traits) ----------
Biom.precipitation <- with(biom_neutral,
                           biom - p["R2.exp"] * R2.exp + mean(p["R2.exp"] * R2.exp, na.rm = TRUE))

plot(Biom.precipitation ~ biom_neutral$lluvia_anual,
     xlab = "Annual precipitation", ylab = "",
     bty = "l", pch = 19, col = pt_biom, cex = 2.5,
     col.axis = ax_col, col.lab = ax_col)

double_curve(
  fun  = function(x) p["(Intercept)"] + p["lluvia_anual"] * x + mean(p["R2.exp"] * biom_neutral$R2.exp, na.rm = TRUE),
  from = min(biom_neutral$lluvia_anual, na.rm = TRUE),
  to   = max(biom_neutral$lluvia_anual, na.rm = TRUE),
  col  = col_biom
)

#mtext("B", side = 3, line = 1.5, adj = 0, font = 2, cex = 2.5, col = ax_col)

# ============================================================
# B) BEF model: BEF ~ R2.exp + lluvia_anual  (data = nn.1)
# ============================================================

m_bef <- lm(BEF ~ R2.exp + lluvia_anual, data = nn.1)
p2 <- coefficients(m_bef)

# ---------- C: BEF ~ explained-by-traits (ajustada por precip) ----------
BEF.neutral <- with(nn.1,
                    BEF - p2["lluvia_anual"] * lluvia_anual + mean(p2["lluvia_anual"] * lluvia_anual, na.rm = TRUE))

plot(BEF.neutral ~ nn.1$R2.exp,
     xlab = "", ylab = "Biomass-Richness slope (BEF)",
     bty = "l", pch = 19, col = pt_bef, cex = 2.5,
     col.axis = ax_col, col.lab = ax_col)

double_curve(
  fun  = function(x) p2["(Intercept)"] + mean(p2["lluvia_anual"] * nn.1$lluvia_anual, na.rm = TRUE) + p2["R2.exp"] * x,
  from = min(nn.1$R2.exp, na.rm = TRUE),
  to   = max(nn.1$R2.exp, na.rm = TRUE),
  col  = col_bef
)

mtext("Variance in species occurrences\nexplained by traits", side = 1, line = 6, cex = 2, col = ax_col)
mtext("B", side = 3, line = 1.5, adj = 0, font = 2, cex = 2.5, col = ax_col)

text(x = 0.42, y = 0.8, "F(2,16): 8.3\nP: 0.003\nR²: 0.51", cex = 2, pos = 1, col = ax_col)

# ----------  BEF ~ precip (ajustada por explained-by-traits) ----------
BEF.precipitation <- with(nn.1,
                          BEF - p2["R2.exp"] * R2.exp + mean(p2["R2.exp"] * R2.exp, na.rm = TRUE))

plot(BEF.precipitation ~ nn.1$lluvia_anual,
     xlab = "Annual precipitation", ylab = "",
     bty = "l", pch = 19, col = pt_bef, cex = 2.5,
     col.axis = ax_col, col.lab = ax_col)

double_curve(
  fun  = function(x) p2["(Intercept)"] + p2["lluvia_anual"] * x + mean(p2["R2.exp"] * nn.1$R2.exp, na.rm = TRUE),
  from = min(nn.1$lluvia_anual, na.rm = TRUE),
  to   = max(nn.1$lluvia_anual, na.rm = TRUE),
  col  = col_bef
)

#mtext("D", side = 3, line = 1.5, adj = 0, font = 2, cex = 2.5, col = ax_col)