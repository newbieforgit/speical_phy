#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <ctime>

const double g = 9.81;
const double l1 = 1.0;
const double l2 = 1.0;
const double m1 = 1.0;
const double m2 = 1.0;
const double dt = 0.001;
struct State {
    double theta1, omega1, theta2, omega2;
};

struct Derivatives {
    double dtheta1, domega1, dtheta2, domega2;
};

std::vector<std::vector<double>> inv_matrix(const std::vector<std::vector<double>>& matrix) {
    double a = matrix[0][0];
    double b = matrix[0][1];
    double c = matrix[1][0];
    double d = matrix[1][1];
    double det = a * d - b * c;
    return {{d / det, -b / det}, {-c / det, a / det}};
}

std::vector<std::vector<double>> matrix_product(const std::vector<std::vector<double>>& matrix1, const std::vector<std::vector<double>>& matrix2) {
    /*if (matrix1.size( ) != matrix2[0].size()){
        std::cout<< "wrong dimension" << std::endl;
        std::exit(1);
    }*/
    std::vector<std::vector<double>> result(matrix1.size(), std::vector<double>(matrix2[0].size(), 0));
    for (size_t i = 0; i < matrix1.size(); i++) {
        for (size_t j = 0; j < matrix2[0].size(); j++) {
            for (size_t k = 0; k < matrix2.size(); k++) {
                result[i][j] += matrix1[i][k] * matrix2[k][j];
            }
        }
    }
    return result;
}

std::vector<std::vector<double>> M_1(const State& state) {
    double mu = 1.0 + m1 / m2;
    double delta = state.theta1 - state.theta2;
    std::vector<std::vector<double>> om = {
        {mu * l1 * l1, l1 * l2 * std::cos(delta)},
        {l1 * l2 * std::cos(delta), l2 * l2}
    };
    return inv_matrix(om);
}

std::vector<std::vector<double>> M_2(const State& state) {
    double mu = 1.0 + m1/m2;
    double delta = state.theta1 - state.theta2;
    return {
        {-state.omega2 * state.omega2 * l1 * l2 * std::sin(delta) - mu * g * l1 * std::sin(state.theta1)},
        {state.omega1 * state.omega1 * l1 * l2 * std::sin(delta) - g * l2 * std::sin(state.theta2)}
    };
}

Derivatives calculate(const State& state) {
    Derivatives d;
    d.dtheta1 = state.omega1;
    d.dtheta2 = state.omega2;
    std::vector<std::vector<double>> new_omega = matrix_product(M_1(state), M_2(state));
    d.domega1 = new_omega[0][0];
    d.domega2 = new_omega[1][0];
    return d;
}

State rk4Step(const State& state) {
    Derivatives k1 = calculate(state);
    State temp1 = {state.theta1 + k1.dtheta1 * dt / 2, state.omega1 + k1.domega1 * dt / 2,
                  state.theta2 + k1.dtheta2 * dt / 2, state.omega2 + k1.domega2 * dt / 2};
    Derivatives k2 = calculate(temp1);
    
    State temp2 = {state.theta1 + k2.dtheta1 * dt / 2, state.omega1 + k2.domega1 * dt / 2,
            state.theta2 + k2.dtheta2 * dt / 2, state.omega2 + k2.domega2 * dt / 2};
    Derivatives k3 = calculate(temp2);
    
    State temp3 = {state.theta1 + k3.dtheta1 * dt, state.omega1 + k3.domega1 * dt,
            state.theta2 + k3.dtheta2 * dt, state.omega2 + k3.domega2 * dt};
    Derivatives k4 = calculate(temp3);
    
    State nextState;
    nextState.theta1 = state.theta1 + dt * (k1.dtheta1 + 2 * k2.dtheta1 + 2 * k3.dtheta1 + k4.dtheta1) / 6;
    nextState.omega1 = state.omega1 + dt * (k1.domega1 + 2 * k2.domega1 + 2 * k3.domega1 + k4.domega1) / 6;
    nextState.theta2 = state.theta2 + dt * (k1.dtheta2 + 2 * k2.dtheta2 + 2 * k3.dtheta2 + k4.dtheta2) / 6;
    nextState.omega2 = state.omega2 + dt * (k1.domega2 + 2 * k2.domega2 + 2 * k3.domega2 + k4.domega2) / 6;
     
    nextState.theta1 = std::fmod(nextState.theta1 + M_PI, 2 * M_PI)  - M_PI;
    nextState.theta2 = std::fmod(nextState.theta2 + M_PI, 2 * M_PI)  - M_PI;
    
    return nextState;
}

int main() {
    State state = { M_PI * 40 / 180, 0, M_PI * 55 / 180, 0 };
    time_t start, finish;
    double duration;
    start = clock();
    std::ofstream outFile("/Users/sylee(0221)/Downloads/special_phy/double_pendulum_x.csv");
    outFile << "Time,Theta1,Omega1,Theta2,Omega2\n";
    
    for (int i = 0; i < 5000; ++i) {
        state = rk4Step(state);
        outFile << i * dt << "," << state.theta1 * (180 / M_PI)  << "," << state.omega1 * (180 / M_PI) << "," << state.theta2 * (180 / M_PI) << "," << state.omega2 * (180 / M_PI) << "\n";
    }
    
    outFile.close();
    finish = clock();
    duration = double(finish -start)/CLOCKS_PER_SEC;
    std::cout << "Simulation complete. Data saved to double_pendulum_x.csv" <<"\n time passed "<< duration<<"sec" << std::endl;
    return 0;
}
