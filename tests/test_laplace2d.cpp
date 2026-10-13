// Pruebas sin dependencias externas: cada prueba es una función y
// CHECK registra el fallo con el archivo y la línea.

#include <cmath>
#include <cstdio>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

#include "laplace2d.hpp"

using namespace calor;

static int fallos = 0;
#define CHECK(cond)                                                           \
    do {                                                                      \
        if (!(cond)) {                                                        \
            std::printf("  FALLA %s:%d  %s\n", __FILE__, __LINE__, #cond);    \
            ++fallos;                                                         \
        }                                                                     \
    } while (0)

static const double PI = 3.14159265358979323846;

static Placa placa_superior_caliente() {
    Placa p;
    p.superior = [](double) { return 100.0; };
    return p;
}

void serie_de_fourier_en_el_centro_vale_un_cuarto() {
    // Superponiendo los 4 bordes a 100 °C se obtiene 100 °C en todo el dominio,
    // así que por simetría cada borde aporta 25 °C en el centro.
    CHECK(std::abs(solucion_fourier(0.5, 0.5, 1, 1, 100) - 25.0) < 1e-6);
}

void los_tres_metodos_coinciden_con_fourier() {
    double exacta = solucion_fourier(0.5, 0.25, 1, 1, 100);
    for (Metodo m : {Metodo::Jacobi, Metodo::GaussSeidel, Metodo::SOR}) {
        Malla malla(placa_superior_caliente(), 41, 41);
        Resultado r = malla.resolver(m, 1e-10);
        CHECK(r.convergio);
        CHECK(std::abs(malla.interpolar(0.5, 0.25) - exacta) < 0.05);
    }
}

void campo_lineal_es_exacto() {
    // T = 3 + 2x − y cumple Laplace y el esquema de 5 puntos lo reproduce exacto
    Placa p;
    p.ancho = 2.0;
    p.alto = 1.0;
    auto T = [](double x, double y) { return 3 + 2 * x - y; };
    p.inferior = [&](double x) { return T(x, 0); };
    p.superior = [&](double x) { return T(x, 1); };
    p.izquierdo = [&](double y) { return T(0, y); };
    p.derecho = [&](double y) { return T(2, y); };
    Malla malla(p, 41, 21);
    malla.resolver(Metodo::SOR, 1e-12);
    double err = 0;
    for (int j = 0; j < malla.ny(); ++j)
        for (int i = 0; i < malla.nx(); ++i)
            err = std::max(err, std::abs(malla(i, j) - T(malla.x(i), malla.y(j))));
    CHECK(err < 1e-9);
}

void convergencia_de_segundo_orden() {
    // Borde superior senoidal: T = sin(πx)·sinh(πy)/sinh(π) es la solución exacta.
    Placa p;
    p.superior = [](double x) { return std::sin(PI * x); };
    auto exacta = [](double x, double y) {
        return std::sin(PI * x) * std::sinh(PI * y) / std::sinh(PI);
    };
    std::vector<double> errores;
    for (int n : {11, 21, 41}) {
        Malla malla(p, n, n);
        malla.resolver(Metodo::SOR, 1e-13);
        double err = 0;
        for (int j = 0; j < n; ++j)
            for (int i = 0; i < n; ++i)
                err = std::max(err, std::abs(malla(i, j) - exacta(malla.x(i), malla.y(j))));
        errores.push_back(err);
    }
    // al dividir h entre 2 el error debe caer ~4 veces
    for (size_t k = 1; k < errores.size(); ++k) {
        double orden = std::log2(errores[k - 1] / errores[k]);
        CHECK(orden > 1.9 && orden < 2.1);
    }
}

void sor_necesita_menos_iteraciones() {
    Malla a(placa_superior_caliente(), 31, 31), b(placa_superior_caliente(), 31, 31),
        c(placa_superior_caliente(), 31, 31);
    int jac = a.resolver(Metodo::Jacobi, 1e-8).iteraciones;
    int gs = b.resolver(Metodo::GaussSeidel, 1e-8).iteraciones;
    int sor = c.resolver(Metodo::SOR, 1e-8).iteraciones;
    CHECK(gs < jac);
    CHECK(sor * 5 < gs);
}

void entradas_invalidas() {
    bool lanzo = false;
    try {
        Malla m(placa_superior_caliente(), 21, 21);
        m.resolver(Metodo::SOR, 1e-8, 100, 2.5);
    } catch (const std::invalid_argument&) {
        lanzo = true;
    }
    CHECK(lanzo);

    lanzo = false;
    try {
        Malla m(placa_superior_caliente(), 21, 31);   // h distinto en x e y
    } catch (const std::invalid_argument&) {
        lanzo = true;
    }
    CHECK(lanzo);
}

int main() {
    std::vector<std::pair<std::string, std::function<void()>>> pruebas = {
        {"serie_de_fourier_en_el_centro_vale_un_cuarto", serie_de_fourier_en_el_centro_vale_un_cuarto},
        {"los_tres_metodos_coinciden_con_fourier", los_tres_metodos_coinciden_con_fourier},
        {"campo_lineal_es_exacto", campo_lineal_es_exacto},
        {"convergencia_de_segundo_orden", convergencia_de_segundo_orden},
        {"sor_necesita_menos_iteraciones", sor_necesita_menos_iteraciones},
        {"entradas_invalidas", entradas_invalidas},
    };
    for (auto& [nombre_prueba, f] : pruebas) {
        int antes = fallos;
        f();
        std::printf("%s %s\n", fallos == antes ? "[ OK ]" : "[FALLA]", nombre_prueba.c_str());
    }
    std::printf("\n%zu pruebas, %d verificaciones fallidas\n", pruebas.size(), fallos);
    return fallos == 0 ? 0 : 1;
}
