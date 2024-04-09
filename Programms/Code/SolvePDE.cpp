//
// Реализация функций для решения Уравнений в частных производных(УЧП)
//

#include <string>
#include "SolvePDE.h"
#include "algebra.cpp"
#include "FileIO.h"


//
//печать вектора
template<typename LT>
void out(vector<LT> vec)
{
    int n = vec.size();
    for (int i = 0; i < n; ++i)
    {
        cout << fixed << setprecision(2) << setw(8) << setfill(' ') << vec[i] << "  ";
    }
    cout << endl;

}
//TODO: вынести в отдельный файл
//решение СЛАУ методом правой прогонки
//A - B + C = -D
template<typename DT>
std::vector<DT> TridiagonalMatrixAlgorithm(
        std::vector<DT> a,
        std::vector<DT> b,
        std::vector<DT> c,
        std::vector<DT> d
        ){
    int n = b.size();
    std::vector<DT> alphas({ c[0] / b[0] });
    std::vector<DT> betas({d[0] / b[0]});
    for (int i = 1; i < n-1; ++i)
    {
        DT denom = b[i] - a[i] * alphas[i - 1];
        alphas.push_back(c[i] / denom);
        betas.push_back((d[i] + a[i] * betas[i-1]) / denom);
    }
    betas.push_back((d[n - 1] + a[n - 1] * betas[n - 2]) / (b[n-1] - a[n-1] * alphas[n - 2]));
    std::vector<DT> SolutionX({betas[n-1]});
    for (int i = n - 2; i >= 0; --i)
    {
        SolutionX.push_back(alphas[i] * SolutionX[n - i - 2] + betas[i]);
    }
    reverse(SolutionX.begin(), SolutionX.end());
    return SolutionX;
}

/*Инициализация начального состояния во всех точках стержня
 * arguments: n  - amount of points
 *            u0 - initial temperature
 * return:    initial state of the system
 * */

std::vector<double> init_state(int n, double u0)
{
    std::vector<double> result(n, u0);
    return result;
}

std::vector<double> init_state(int n, double h, PDE_data& test)
{
    std::vector<double> result(n,0);
    double x_i = 0;
    result[0] = test.G_left(x_i);
    for(int i = 1; i < n-1; ++i)
    {
        x_i += h;
        result[i] = test.initFunction(x_i);
    }
    result[n-1] = test.G_right(x_i+h);
    return result;
}

double left_point_state(double u0){
    return u0;
}

double right_point_state(double kappa, double mu, double y_im){
    return kappa * y_im + mu;
}

template<typename F>
double a(F K, double x_i, double x_im){
   // return 0.5 * (K(x_i)+K(x_im));
    //return K(x_i - 0.5*(x_i-x_im));
    //return sqrt(K(x_i)*K(x_im));
    return 2*K(x_i)*K(x_im) / (K(x_i)+K(x_im));
}

double w(double a, double u_i, double u_im, double h) {
    return a * (u_i-u_im)/h;
}

// То что написал Ваня
bool ExplicitScheme(double tau, double h, double sigma, PDE_data test, std::string filename="ExpScheme"){
    double c = test.c;
    double rho = test.rho;
    double t_0 = 0;
    double T = test.T;
    double x_0 = 0;
    double X = test.L;

    // Шаги по времени и пространству
    int num_time_steps = static_cast<int>((T-t_0) / tau);
    int num_space_steps = static_cast<int>((X - x_0)/h);

    // Инициализация начального состояния
    //std::vector<double> state_0 = init_state(num_space_steps, u_0); //TODO: расширить init_state
    std::vector<double> state_0 = init_state(num_space_steps+1, h, test);
    std::vector<double> As(num_space_steps+1, 0);
    std::vector<double> Cs(num_space_steps+1, 0);
    std::vector<double> Bs(num_space_steps+1, 0);
    std::vector<double> Fs(num_space_steps+1, 0);

    std::string path = "./OutputData/" + filename;
    std::ofstream fpoints(path);
    std::cout << "log[INFO]: Starting ExplicitScheme" << std::endl;
    std::cout << "log[INFO]: Opening a file \"" << filename << "\" to write..." << std::endl;
    if (fpoints.is_open())
    {
        double t_i = t_0;
        std::vector<double> state_i = state_0;
        int ind = 0;
        writeVectorToFile(fpoints, t_i, state_i);
        double x_i = x_0;
        for(int j = 0; j <= num_time_steps; ++j) {
            t_i += tau;
            Cs[0] = 1.;
            Bs[0] = 0.;
            As[0] = 0.;
            Fs[0] = state_0[0];
            Bs[num_space_steps] = 0.;
            As[num_space_steps] = 0.;
            Cs[num_space_steps] = 1.;
            Fs[num_space_steps] = state_0[num_space_steps];
            for (int i = 1; i < num_space_steps; ++i) {
                x_i += h;
                double a_i = a(test.K_ptr, x_i, x_i - h);
                double a_ip = a(test.K_ptr, x_i + h, x_i);
                As[i] = sigma / h * a_i;
                Bs[i] = sigma / h * a_ip;
                Cs[i] = (As[i] + Bs[i] + c * rho * h / tau);
                Fs[i] = (c * rho * h / tau * state_i[i] +
                        (1 - sigma) * (w(a_ip, state_i[i + 1], state_i[i], h) - w(a_i, state_i[i], state_i[i - 1], h)));

            }
            state_i = TridiagonalMatrixAlgorithm(As, Cs, Bs, Fs);
            writeVectorToFile(fpoints, t_i, state_i);
        }
        fpoints.close();
        return true;
    }
    else {
        std::cout << "log[ERROR]: Couldn't open or create a file" << std::endl;
        return false;
    }
};




// Моя переделка функций Вани:


std::vector<double> progonka(std::vector<double>& a, std::vector<double>& b, std::vector<double>& c, std::vector<double>& d) {

    int n = b.size();
    std::vector<double> alph(n, 0);
    std::vector<double> beth(n, 0);
    std::vector<double> solve(n, 0);
    double tmp;

    //vectinput();
    alph[0] = 0;
    beth[0] = 0;
    tmp = b[0];
    alph[1] = c[0] / tmp;
    beth[1] = d[0] / tmp;

    //std::cout << alph[1] << " " << beth[1] << "\n";
    for (int i = 1; i < n - 1; i++) {
        tmp = b[i] - a[i] * alph[i];
        alph[i + 1] = c[i] / tmp;//(b[i] - a[i] * alph[i]);
        beth[i + 1] = (d[i] + a[i] * beth[i]) / tmp; // (b[i] - a[i] * alph[i]);
    }

    solve[n - 1] = (d[n - 1] + a[n - 1] * beth[n - 1]) / (b[n - 1] - a[n - 1] * alph[n - 1]);

    for (int i = n - 2; i >= 0; i--) {
        solve[i] = alph[i + 1] * solve[i + 1] + beth[i + 1];
    }

    return solve;
}

/* Функция для решения PDE в случае 1 (по методичке) */
bool SolvePDE_1(PDE_data test, double tau, double h, double sigma, std::string filename){


    double t_0 = 0; // Начальное время
    double x_0 = 0; // Начальное условие?

    // Количество шагов по времени и пространству
    int num_time_steps = static_cast<int>((test.T-t_0) / tau);
    int num_space_steps = static_cast<int>((test.L - x_0) / h);

    // TODO: брать граничное условие из теста (здесь начальная температура)
    double u_0 = 10;

    // Инициализация начального состояния
    std::vector<double> state_0 = init_state(num_space_steps, u_0); //TODO: расширить init_state
    std::vector<double> As(num_space_steps, 0);
    std::vector<double> Cs(num_space_steps, 0);
    std::vector<double> Bs(num_space_steps, 0);
    std::vector<double> Fs(num_space_steps, 0);

    std::string path = "./OutputData/" + filename + ".txt";
    std::ofstream fpoints(path);
    std::cout << "log[INFO]: Starting ExplicitScheme" << std::endl;
    std::cout << "log[INFO]: Opening a file \"" << filename << "\" to write..." << std::endl;

    if (fpoints.is_open()) {

        double t_i = t_0;
        std::vector<double> state_i = state_0;
        int ind = 0;
        writeVectorToFile(fpoints, t_i, state_i);
        double x_i = x_0;
        Cs[0] = 1;
        Bs[0] = 0;
        As[0] = 0;
        Fs[0] = u_0;
        Bs[num_space_steps-1] = 0;
        As[num_space_steps-1] = 0;
        Cs[num_space_steps-1] = 1;
        Fs[num_space_steps-1] = u_0;

        for (int j = 0; j < num_time_steps; ++j) {
            t_i += tau;

            for (int i = 1; i < num_space_steps - 1; ++i) {
                x_i += h;

                //double a_i = a(test.K_ptr, x_i, x_i - h);
                double a_i = 0.5 * (test.K(x_i) + test.K(x_i - h));

                //double a_ip = a(test.K_ptr, x_i + h, x_i);
                double a_ip = 0.5 * (test.K(x_i + h) + test.K(x_i - h));

                As[i] = sigma / h * a_i;
                Bs[i] = sigma / h * a_ip;
                Cs[i] = As[i] + Bs[i] + test.c * test.rho * h / tau;
                Fs[i] = test.c * test.rho * h / tau * state_i[i]
                        + (1 - sigma) * (w(a_ip, state_i[i + 1],
                                           state_i[i], h) - w(a_i, state_i[i], state_i[i - 1], h));
            }

            Cs = (-1.) * Cs;
            Fs = (-1.) * Fs;
            state_i = progonka(As, Cs, Bs, Fs);
            writeVectorToFile(fpoints, t_i, state_i);
        }

        fpoints.close();
        return true;

    } else {
        std::cout << "log[ERROR]: Couldn't open or create a file" << std::endl;
        return false;
    }
}
