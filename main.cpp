#include <iostream>
#include <optional>
#include <map>
#include <fstream>
#include <vector>

#define map_value(X) {#X, &X}

#define checker(X) if(!X) {std::string error = #X; \
error += " is missing\n";                          \
throw std::runtime_error(error);                   \
}

auto input(std::ifstream& in){
    std::optional<double> a1, b1, g1, a2, b2, g2, A, B, N;
    std::map<std::string, std::optional<double>*> map{
            map_value(A),
            map_value(B),
            map_value(N),
            map_value(a1),
            map_value(b1),
            map_value(g1),
            map_value(a2),
            map_value(b2),
            map_value(g2)
    };
    std::string temp;
    while (std::getline(in, temp) && temp != "========"){
        if (temp.empty()){
            continue;
        }
        size_t index = std::find(temp.begin(), temp.end(), '=') - temp.begin();
        std::string key = temp.substr(0, index);
        std::string value = temp.substr(index + 1, temp.size());
        double d_value = std::stod(value);
        auto it = map.find(key);
        if (it == map.end()){
            continue;
        }
        *(it->second) = d_value;
    }
    checker(A);
    checker(B);
    checker(N);
    checker(a1);
    checker(b1);
    checker(g1);
    checker(a2);
    checker(b2);
    checker(g2);
    std::vector<double> x;
    x.resize(*N);
    for (auto& elem : x) {
        in >> elem;
    }
    return std::make_tuple(*A, *B, *a1, *b1, *g1, *a2, *b2, *g2, x);
}

using Value = decltype(input(std::declval<std::ifstream&>()));


#define LOG(X) std::cout << #X << " " << X << "\n"

#define LOG_RESULT(X)  \
LOG(X.value);           \
LOG(X.code);           \
LOG(X.points);         \
LOG(X.point_without_accuracity); \
LOG(X.point_with_min); \
LOG(X.point_with_max); \

double p(double x){
    return 42;
}

double q(double x){
    return x * x;
}

double f (double x){
    return x * 0.5;
}

double RungeKutta(double x, double y, double h){
    double k1 = h * (f(x) - q(x) * y) / p(x);
    double k2 = h * (f(x + h / 2.0) - q(x + h / 2.0) * (y + k1 / 2.0)) / p(x + h / 2.0);
    double k3 = h * (f(x + h / 2.0) - q(x + h / 2.0) * (y + k2 / 2.0)) / p(x + h / 2.0);
    double k4 = h * (f(x + h) - q(x + h) * (y + k3)) / p(x + h);
    double k5 = h * (f(x + h) - q(x + h) * (y + k4)) / p(x + h);

    return y + (k1 + 2.0 * k2 + 2.0 * k3 + k4 + k5) / 6.0;
}

double dif_run(Value value){
    auto&& [A, B, a1, b1, g1, a2, b2, g2, x] = value;
    auto N = x.size();
    double h = (B - A) / (N - 1);
    std::vector<double> x_points(N);
    std::vector<double> y_points(N);
    std::vector<double> y_prime_points(N);
    std::vector<double> A_coeffs(N);
    std::vector<double> B_coeffs(N);
    // Прямой ход метода прогонки
    // Устанавливаем начальные значения для граничных условий
    y_points[0] = b1;
    y_prime_points[0] = (g1 - a1 * y_points[0]) / a1;

    for (int i = 1; i < N; ++i) {
        x_points[i] = A + i * h;

        // Коэффициенты прогонки
        A_coeffs[i] = -q(x_points[i]) * h / p(x_points[i]) / (h + p(x_points[i - 1]) * A_coeffs[i - 1]);
        B_coeffs[i] = (f(x_points[i]) * h / p(x_points[i]) - B_coeffs[i - 1] * p(x_points[i - 1])) /
                      (h + p(x_points[i - 1]) * A_coeffs[i - 1]);

        // Используем метод Рунге-Кутта для вычисления y(x) и y'(x)
        y_points[i] = RungeKutta(x_points[i - 1], y_points[i - 1], h);
        y_prime_points[i] = (f(x_points[i]) - q(x_points[i]) * y_points[i]) / p(x_points[i]);
    }

    // Обратный ход метода прогонки
    for (int i = N - 2; i >= 0; --i) {
        y_points[i] = A_coeffs[i + 1] * y_points[i + 1] + B_coeffs[i + 1];
    }

    // Выводим результаты
    std::cout << "Point\t\t x\t\t y(x)\t\t y'(x)" << std::endl;
    for (int i = 0; i < N; ++i) {
        std::cout << i << "\t\t " << x_points[i] << "\t\t " << y_points[i] << "\t\t " << y_prime_points[i] << std::endl;
    }
}

int main() {
    std::ifstream in("in.txt");
    auto value = input(in);
    dif_run(value);
    return 0;
}
