#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <math.h>

class Vector {
public:
    double x, y;

    Vector(double x, double y) : x(x), y(y) {}

    Vector operator+(const Vector& rhs) const {
        return Vector(this->x + rhs.x, this->y + rhs.y);
    }

    Vector operator-(const Vector& rhs) const {
        return Vector(this->x - rhs.x, this->y - rhs.y);
    }

    friend Vector operator*(const Vector& v, double scalar) {
        return Vector(v.x * scalar, v.y * scalar);
    }

    friend Vector operator*(double scalar, const Vector& v) {
        return Vector(v.x * scalar, v.y * scalar);
    }

    Vector operator/(double scalar) const {
        return Vector(this->x / scalar, this->y / scalar);
    }

    friend std::ostream& operator<<(std::ostream& os, const Vector& v) {
        os << "(" << v.x << ", " << v.y << ")";
        return os;
    }
};


double objectiveFunction(const Vector& point) {
    return point.x * point.x + point.y * point.y;
}

Vector nelderMead(double alpha = 1, double beta = 0.5, double gamma = 2, int maxiter = 10) {
    std::vector<Vector> history_b;
    std::vector<Vector> history_g;
    std::vector<Vector> history_w;
    Vector v1(1.0, 1.0);
    Vector v2(2.0, 1.0);
    Vector v3(1.0, 2.0);

    for (int i = 0; i < maxiter; ++i) {
        std::vector<std::pair<Vector, double>> points = {
            {v1, objectiveFunction(v1)},
            {v2, objectiveFunction(v2)},
            {v3, objectiveFunction(v3)}
        };

        std::sort(points.begin(), points.end(),
            [](const std::pair<Vector, double>& a, const std::pair<Vector, double>& b) -> bool {
                return a.second < b.second;
            });

        Vector b = points[0].first;
        Vector g = points[1].first;
        Vector w = points[2].first;

        Vector mid = (g + b) * 0.5;

        Vector xr = mid + alpha * (mid - w);
        if (objectiveFunction(xr) < objectiveFunction(g)) {
            w = xr;
        }
        if (objectiveFunction(xr) < objectiveFunction(b)) {
            Vector xe = mid + gamma * (xr - mid);
            w = (objectiveFunction(xe) < objectiveFunction(xr)) ? xe : xr;
        }
        else if (objectiveFunction(xr) > objectiveFunction(g)) {
            Vector xc = mid + beta * (w - mid);
            w = (objectiveFunction(xc) < objectiveFunction(w)) ? xc : w;
        }
        history_b.push_back(b);
        history_g.push_back(g);
        history_w.push_back(w);

        v2 = b; v1 = g; v3 = w;
    }
    std::cout << "Iteration Best PointFunction Value" << std::endl;
    for (size_t i = 0; i < history_b.size(); ++i) {
        std::cout << std::setw(9) << i + 1 << " : " << history_b[i]
            << " = " << objectiveFunction(history_b[i]) << std::endl;
    }
    std::cout << "Iteration Good PointFunction Value" << std::endl;
    for (size_t i = 0; i < history_g.size(); ++i) {
        std::cout << std::setw(9) << i + 1 << " : " << history_g[i]
            << " = " << objectiveFunction(history_g[i]) << std::endl;
    }
    std::cout << "Iteration Worst PointFunction Value" << std::endl;
    for (size_t i = 0; i < history_w.size(); ++i) {
        std::cout << std::setw(9) << i + 1 << " : " << history_w[i]
            << " = " << objectiveFunction(history_w[i]) << std::endl;
    }


    return history_b.back();
}

int main() {
    std::cout << "Result of Nelder-Mead algorithm:" << std::endl;
    Vector xk = nelderMead();
    std::cout << "Best point is: " << xk << std::endl;
    return 0;
}
