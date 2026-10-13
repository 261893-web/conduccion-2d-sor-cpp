#pragma once
// Conducción de calor 2D en estado estacionario (ecuación de Laplace)
// sobre una placa rectangular, resuelta con métodos iterativos.

#include <functional>
#include <string>
#include <vector>

namespace calor {

enum class Metodo { Jacobi, GaussSeidel, SOR };

std::string nombre(Metodo m);

// Temperatura impuesta en un borde, en función de la coordenada a lo largo del borde.
using Borde = std::function<double(double)>;

struct Placa {
    double ancho = 1.0;   // dirección x [m]
    double alto = 1.0;    // dirección y [m]
    Borde inferior = [](double) { return 0.0; };   // y = 0, función de x
    Borde superior = [](double) { return 0.0; };   // y = alto, función de x
    Borde izquierdo = [](double) { return 0.0; };  // x = 0, función de y
    Borde derecho = [](double) { return 0.0; };    // x = ancho, función de y
};

struct Resultado {
    int iteraciones = 0;
    double cambio_final = 0.0;   // máximo |ΔT| en la última iteración
    bool convergio = false;
    double segundos = 0.0;
};

class Malla {
public:
    // nx, ny: número de nodos en cada dirección (incluye bordes).
    // El paso h debe ser igual en x e y.
    Malla(const Placa& placa, int nx, int ny);

    Resultado resolver(Metodo metodo, double tolerancia = 1e-8,
                       int max_iter = 200000, double omega = -1.0);

    // ω óptimo teórico para SOR en un rectángulo con Dirichlet.
    double omega_optimo() const;

    double& operator()(int i, int j) { return T_[j * nx_ + i]; }
    double operator()(int i, int j) const { return T_[j * nx_ + i]; }

    int nx() const { return nx_; }
    int ny() const { return ny_; }
    double h() const { return h_; }
    double x(int i) const { return i * h_; }
    double y(int j) const { return j * h_; }

    // Temperatura en un punto cualquiera por interpolación bilineal.
    double interpolar(double x, double y) const;

    void reiniciar_interior(double valor = 0.0);
    void guardar_csv(const std::string& ruta) const;

private:
    void aplicar_bordes(const Placa& placa);

    int nx_, ny_;
    double h_;
    std::vector<double> T_;
};

// Solución analítica por series de Fourier para una placa con T = T1
// en el borde superior y T = 0 en los otros tres (Incropera, sección 4.2).
double solucion_fourier(double x, double y, double ancho, double alto,
                        double T1, int terminos = 201);

}  // namespace calor
