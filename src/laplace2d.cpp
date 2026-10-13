#include "laplace2d.hpp"

#include <chrono>
#include <cmath>
#include <fstream>
#include <stdexcept>

namespace calor {

namespace {
constexpr double PI = 3.14159265358979323846;
}

std::string nombre(Metodo m) {
    switch (m) {
        case Metodo::Jacobi: return "Jacobi";
        case Metodo::GaussSeidel: return "Gauss-Seidel";
        case Metodo::SOR: return "SOR";
    }
    return "?";
}

Malla::Malla(const Placa& placa, int nx, int ny) : nx_(nx), ny_(ny) {
    if (nx < 3 || ny < 3) {
        throw std::invalid_argument("La malla necesita al menos 3 nodos por dirección");
    }
    h_ = placa.ancho / (nx - 1);
    double hy = placa.alto / (ny - 1);
    if (std::abs(h_ - hy) > 1e-12 * h_) {
        throw std::invalid_argument("El paso debe ser igual en x e y (ajusta nx o ny)");
    }
    T_.assign(static_cast<size_t>(nx) * ny, 0.0);
    aplicar_bordes(placa);
}

void Malla::aplicar_bordes(const Placa& p) {
    for (int i = 0; i < nx_; ++i) {
        (*this)(i, 0) = p.inferior(x(i));
        (*this)(i, ny_ - 1) = p.superior(x(i));
    }
    // las esquinas quedan con el valor de los bordes horizontales;
    // el esquema de 5 puntos nunca las usa
    for (int j = 1; j < ny_ - 1; ++j) {
        (*this)(0, j) = p.izquierdo(y(j));
        (*this)(nx_ - 1, j) = p.derecho(y(j));
    }
}

void Malla::reiniciar_interior(double valor) {
    for (int j = 1; j < ny_ - 1; ++j)
        for (int i = 1; i < nx_ - 1; ++i) (*this)(i, j) = valor;
}

double Malla::omega_optimo() const {
    // radio espectral de Jacobi para el problema modelo
    double rho = 0.5 * (std::cos(PI / (nx_ - 1)) + std::cos(PI / (ny_ - 1)));
    return 2.0 / (1.0 + std::sqrt(1.0 - rho * rho));
}

Resultado Malla::resolver(Metodo metodo, double tol, int max_iter, double omega) {
    if (metodo == Metodo::SOR && omega <= 0) omega = omega_optimo();
    if (metodo == Metodo::SOR && (omega <= 0 || omega >= 2)) {
        throw std::invalid_argument("SOR solo converge con 0 < ω < 2");
    }
    if (metodo == Metodo::GaussSeidel) omega = 1.0;

    auto inicio = std::chrono::steady_clock::now();
    Resultado r;
    std::vector<double> anterior;

    for (r.iteraciones = 1; r.iteraciones <= max_iter; ++r.iteraciones) {
        double max_cambio = 0.0;

        if (metodo == Metodo::Jacobi) {
            anterior = T_;
            auto A = [&](int i, int j) { return anterior[j * nx_ + i]; };
            for (int j = 1; j < ny_ - 1; ++j) {
                for (int i = 1; i < nx_ - 1; ++i) {
                    double nuevo = 0.25 * (A(i + 1, j) + A(i - 1, j) + A(i, j + 1) + A(i, j - 1));
                    max_cambio = std::max(max_cambio, std::abs(nuevo - A(i, j)));
                    (*this)(i, j) = nuevo;
                }
            }
        } else {
            // Gauss-Seidel y SOR actualizan en el mismo arreglo
            double* T = T_.data();
            for (int j = 1; j < ny_ - 1; ++j) {
                for (int k = j * nx_ + 1; k < (j + 1) * nx_ - 1; ++k) {
                    double gs = 0.25 * (T[k + 1] + T[k - 1] + T[k + nx_] + T[k - nx_]);
                    double delta = omega * (gs - T[k]);
                    T[k] += delta;
                    max_cambio = std::max(max_cambio, std::abs(delta));
                }
            }
        }

        r.cambio_final = max_cambio;
        if (max_cambio < tol) {
            r.convergio = true;
            break;
        }
    }
    if (!r.convergio) r.iteraciones = max_iter;

    auto fin = std::chrono::steady_clock::now();
    r.segundos = std::chrono::duration<double>(fin - inicio).count();
    return r;
}

double Malla::interpolar(double xp, double yp) const {
    double fx = xp / h_, fy = yp / h_;
    int i = std::min(static_cast<int>(fx), nx_ - 2);
    int j = std::min(static_cast<int>(fy), ny_ - 2);
    double a = fx - i, b = fy - j;
    return (1 - a) * (1 - b) * (*this)(i, j) + a * (1 - b) * (*this)(i + 1, j) +
           (1 - a) * b * (*this)(i, j + 1) + a * b * (*this)(i + 1, j + 1);
}

void Malla::guardar_csv(const std::string& ruta) const {
    std::ofstream f(ruta);
    if (!f) throw std::runtime_error("No se pudo escribir " + ruta);
    f << "x,y,T\n";
    for (int j = 0; j < ny_; ++j)
        for (int i = 0; i < nx_; ++i) f << x(i) << ',' << y(j) << ',' << (*this)(i, j) << '\n';
}

double solucion_fourier(double x, double y, double W, double H, double T1, int terminos) {
    // θ = (2/π) Σ [(−1)^(n+1) + 1]/n · sin(nπx/W) · sinh(nπy/W) / sinh(nπH/W)
    // Solo sobreviven los n impares. El cociente de senos hiperbólicos se
    // escribe con exponenciales negativas para no desbordar con n grande.
    double suma = 0.0;
    for (int n = 1; n <= terminos; n += 2) {
        double k = n * PI / W;
        double cociente = std::exp(k * (y - H)) * (1 - std::exp(-2 * k * y)) /
                          (1 - std::exp(-2 * k * H));
        suma += 2.0 / n * std::sin(k * x) * cociente;
    }
    return T1 * 2.0 / PI * suma;
}

}  // namespace calor
