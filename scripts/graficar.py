"""Grafica los CSV que genera ./conduccion2d (requiere numpy y matplotlib)."""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

carpeta = Path(sys.argv[1] if len(sys.argv) > 1 else "resultados")
salida = Path(sys.argv[2] if len(sys.argv) > 2 else "figuras")
salida.mkdir(exist_ok=True)

x, y, T = np.loadtxt(carpeta / "campo_temperatura.csv", delimiter=",", skiprows=1).T
n = int(round(np.sqrt(len(T))))
X, Y, T = x.reshape(n, n), y.reshape(n, n), T.reshape(n, n)

w, it = np.loadtxt(carpeta / "barrido_omega.csv", delimiter=",", skiprows=1).T

fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.8))

cf = axs[0].contourf(X, Y, T, levels=np.linspace(0, 100, 21), cmap="inferno")
cs = axs[0].contour(X, Y, T, levels=[5, 10, 25, 50, 75], colors="white", linewidths=0.8)
axs[0].clabel(cs, fmt="%d °C", fontsize=8)
fig.colorbar(cf, ax=axs[0], label="Temperatura [°C]")
axs[0].set_aspect("equal")
axs[0].set_xlabel("x [m]")
axs[0].set_ylabel("y [m]")
axs[0].set_title("Campo de temperatura (borde superior a 100 °C)")

axs[1].semilogy(w, it, "o-", ms=3, color="#1f5fa8")
k = np.argmin(it)
axs[1].annotate(f"mínimo: ω = {w[k]:.2f}\n{int(it[k])} iteraciones", (w[k], it[k]),
                xytext=(w[k] - 0.45, it[k] * 3), arrowprops=dict(arrowstyle="->"))
axs[1].set_xlabel("Factor de relajación ω")
axs[1].set_ylabel("Iteraciones hasta converger")
axs[1].set_title("SOR: sensibilidad a ω (ω = 1 es Gauss-Seidel)")
axs[1].grid(True, which="both", alpha=0.3)

plt.tight_layout()
plt.savefig(salida / "resultados.png", dpi=150)
print(f"Figura guardada en {salida / 'resultados.png'}")
