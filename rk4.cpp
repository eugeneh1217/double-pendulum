#include <array>
#include <iostream>
#include <string>

#define N 1

void _step(double t, std::array<double, N> const &x, std::array<double, N> &x_d)
{
    x_d[0] = x[0];
}

void runge_kutta(
        double t,
        std::array<double, N> const &x,
        double h,
        std::array<double, N> &x_next)
{
    std::array<double, N> k1, x_temp, k2, k3, k4;
    _step(t, x, k1);

    for (unsigned long i = 0; i < N; ++ i)
    {
        x_temp[i] = x[i] + h * k1[i] / 2;
    }
    _step(t + (h/2), x_temp, k2);

    for (unsigned long i = 0; i < N; ++ i)
    {
        x_temp[i] = x[i] + h * k2[i] / 2;
    }
    _step(t + (h/2), x_temp, k3);

    for (unsigned long i = 0; i < N; ++ i)
    {
        x_temp[i] = x[i] + h * k3[i];
    }
    _step(t + (h/2), x_temp, k4);

    for (unsigned long i = 0; i < N; ++ i)
    {
        x_next[i] = x[i] + (h / 6) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
    }
}

int main(int argc, char **argv)
{
    int limit = 10;
    if (argc == 2)
        limit = std::stoi(argv[1]);

    double t = 0;
    double h = 0.1;
    std::array<double, N> x = {1};
    std::cout << x[0] << std::endl;
    for (int i = 0; i < limit; ++ i)
    {
        runge_kutta(t, x, h, x);
        std::cout << x[0] << std::endl;
    }
    std::cout << "limit: " << limit << std::endl;
}