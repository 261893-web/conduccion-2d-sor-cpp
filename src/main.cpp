// Demostración: placa cuadrada con el borde superior a 100 °C y los demás a 0 °C.
// 1) Compara Jacobi, Gauss-Seidel y SOR contra la solución por series de Fourier.
// 2) Barre ω para encontrar el factor de relajación que minimiza las iteraciones.

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>

#include "laplace2d.hpp"

using namespace calor;

int main(int argc, char** argv) {
    const std::string carpeta = argc > 1 ? argv[1] : "resultados";
    std::filesystem::create_directories(carpeta);

    Placa placa;
    placa.superior = [](double) { return 100.0; };
    const int n = 51;
    const double tol = 1e-8;

    // ---------- 1. comparación de métodos ----------
    std::printf("Placa 1 x 1 m, malla %d x %d, tolerancia %.0e\n\n", n, n, tol);
    std::printf("%-14s %12s %12s %16s\n", "Método", "Iteraciones", "Tiempo [ms]", "Error en (0.5, 0.25)");

    const double exacta = solucion_fourier(0.5, 0.25, 1.0, 1.0, 100.0);
    for (Metodo m : {Metodo::Jacobi, Metodo::GaussSeidel, Metodo::SOR}) {
        Malla malla(placa, n, n);
        Resultado r = malla.resolver(m, tol);
        double error = std::abs(malla.interpolar(0.5, 0.25) - exacta);
        std::printf("%-14s %12d %12.1f %16.4f\n", nombre(m).c_str(), r.iteraciones,
                    1e3 * r.segundos, error);
        if (m == Metodo::SOR) malla.guardar_csv(carpeta + "/campo_temperatura.csv");
    }
    std::printf("\nSolución exacta (Fourier) en (0.5, 0.25): %.4f °C\n",
                exacta);

    // ---------- 2. barrido de ω ----------
    Malla malla(placa, n, n);
    std::ofstream f(carpeta + "/barrido_omega.csv");
    f << "omega,iteraciones\n";
    int mejor_iter = 1 << 30;
    double mejor_omega = 1.0;
    for (int k = 0; k <= 49; ++k) {
        double w = 1.0 + 0.02 * k;
        malla.reiniciar_interior();
        Resultado r = malla.resolver(Metodo::SOR, tol, 50000, w);
        f << w << ',' << r.iteraciones << '\n';
        if (r.iteraciones < mejor_iter) {
            mejor_iter = r.iteraciones;
            mejor_omega = w;
        }
    }
    std::printf("ω óptimo teórico: %.4f | mejor ω del barrido: %.2f (%d iteraciones)\n",
                malla.omega_optimo(), mejor_omega, mejor_iter);
    std::printf("Resultados guardados en %s/\n", carpeta.c_str());
    return 0;
}
