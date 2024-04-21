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
//A x_{i-1} - B x_i + C x_{i+1} = -D  (***)
// Если векторы диагоналей исходной системы заданы как A,B,C,D (A - диагональ опд главной, B - главная диагональ, C - диагональ над главной, D- правая часть)
// То для правильного расчёта необходимо передавать A, (-1.)*B, C, (-1.)*D
// Так как прогонка актуальная для системы (***)
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
    result[0] = test.initFunction(x_i);
    if(!test.G_left_type)
        result[0] = test.G_left(x_i);
    for(int i = 1; i < n-1; ++i)
    {
        x_i += h;
        result[i] = test.initFunction(x_i);
    }
    x_i+=h;
    result[n-1] = test.initFunction(x_i);
    if(!test.G_right_type)
        result[n-1] = test.G_right(x_i);
    return result;
}

// Если температура не передана (т.е. K не зависит т температуры u),
// то u присваивается фиктивное значение 0 (какая разница, чему равно u, если в формуле для K оно не используется в return)(аргумент есть, но он не участвует в вычислении - сделано для универсальности)
template<typename F>
double a(F K, double x_i, double x_im, double u_i=0, double u_im=0){
    //return 0.5 * (K(x_i, u_i)+K(x_im, u_im));
    //return K(x_i - 0.5*(x_i-x_im));
    //return sqrt(K(x_i)*K(x_im));

    // Предотвращаем деление на ноль
    if(K(x_i, u_i)+K(x_im, u_im) != 0)
        return 2*K(x_i, u_i)*K(x_im, u_im) / (K(x_i, u_i)+K(x_im, u_im));
    else
        return sqrt(K(x_i, u_i)*K(x_im, u_im));
}

double w(double a, double u_i, double u_im, double h) {
    return a * (u_i-u_im)/h;
}

//*******************IVAN's SANDBOX*************//

// Случай 1 (линейное ур-е)
bool FiniteScheme(double tau, double h, double sigma, PDE_data test, std::string filename="ExpScheme"){

    // Физические параметры
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

    // Создание файла
    std::string path = "./OutputData/" + filename;
    std::ofstream fpoints(path);
    std::cout << "log[INFO]: Starting ExplicitScheme" << std::endl;
    std::cout << "log[INFO]: Opening a file \"" << filename << "\" to write..." << std::endl;
    if (fpoints.is_open())
    {
        double t_i = t_0;
        std::vector<double> state_i = state_0;

        // Запишем шаг по сетке в файл
        std::vector<double> setka(num_space_steps+1);
        for (int i = 0; i < num_space_steps+1; i++) {
            setka[i] = h * i;
        }
        writeVectorToFile(fpoints, 0., setka);

        writeVectorToFile(fpoints, t_i, state_i);
        double x_i = x_0;

        // Эволюция системы во времени
        for(int j = 0; j <= num_time_steps; ++j) {
            t_i += tau;

            // Граничные условия слева

            // 1-го рода
            if(!test.G_left_type){
                Cs[0] = -1.;
                Bs[0] = 0.;
                As[0] = 0.;
                Fs[0] = -state_0[0];
            }

            // 2-го рода
            else {
                double a0 = a(test.K_ptr, x_0+h, x_0);
                double w0 = w(a0, state_i[1], state_i[0], h);
                double kappa = sigma*a0/h / (c*rho*h/(2*tau)+sigma*a0/h);
                double mu = (c*rho*state_i[0]*h/(2*tau)+sigma*test.G_left(t_i)+(1-sigma)*(test.G_left(t_i-tau)+w0))/(c*rho*h/(2*tau)+sigma*a0/h);
                Cs[0] = -1.;
                Bs[0] = -kappa;
                As[0] = 0;
                Fs[0] = -mu;
            }

            // Граничные условия справа
            // 1-го рода
            if(!test.G_right_type){
                Bs[num_space_steps] = 0.;
                As[num_space_steps] = 0.;
                Cs[num_space_steps] = -1.;
                Fs[num_space_steps] = -state_0[num_space_steps];
            }

            // 2-го рода
            else{
                double am = a(test.K_ptr, X, X-h);
                double wn = w(am, state_i[num_space_steps], state_i[num_space_steps-1], h);
                double denom = c * rho * h / (2 * tau) + sigma * am / h;
                double kappa = sigma * am /h / denom;
                double mu = (c * rho * state_i[num_space_steps] * h / (2 * tau) + sigma * test.G_right(t_i) + (1 - sigma) * (test.G_right(t_i-tau) - wn)) / denom;
                Cs[num_space_steps] = -1.;
                Bs[num_space_steps] = 0.;
                As[num_space_steps] = -kappa;
                Fs[num_space_steps] = -mu;
            }

            // Обход пространства
            for (int i = 1; i < num_space_steps; ++i) {
                x_i += h;
                double a_i = a(test.K_ptr, x_i, x_i - h);
                double a_ip = a(test.K_ptr, x_i + h, x_i);
                As[i] = sigma / h * a_i;
                Bs[i] = sigma / h * a_ip;
                Cs[i] = As[i] + Bs[i] + c * rho * h / tau;
                Fs[i] = c * rho * h / tau * state_i[i] +
                        (1 - sigma) * (w(a_ip, state_i[i + 1], state_i[i], h) - w(a_i, state_i[i], state_i[i - 1], h));
            }

            // Получение нового состояния системы
            // A - C + B = - F (не домножаем векторы на -1, так как уже считали домноженные)
            state_i = TridiagonalMatrixAlgorithm(As, Cs, Bs, Fs);

            // Запись в файл
            writeVectorToFile(fpoints, t_i, state_i);
        }
        fpoints.close();
        return true;

    } else {
        std::cout << "log[ERROR]: Couldn't open or create a file" << std::endl;
        return false;
    }
};

// Итерационный метод решения СЛАУ с трёхдиагональной матрицей
// A x_{i-1} + B x_i + C x_{i+1} = D
// Если исходная система задаётся диагоналями,
// То передавать векторы как они есть (не домножать на -1)
vector<double> TripleBigRelaxSolve(const vector<double>& a, const vector<double>& b,
                               const vector<double>& c, const vector<double>& d,
                               const vector<double>& x_0, double EPS=1e-6)
{
    int n = x_0.size();
    int max_iter = 10000;
    int iter = 0;
    double w = 1;
    //LT w = 1.1;
    vector<double> x_now(x_0);
    vector<double> x_prev;

    do {
        x_prev = x_now;
        x_now[0] = (d[0] - c[0] * x_prev[1]);
        x_now[0] *= w;
        x_now[0] /= b[0];
        x_now[0] += (1 - w) * x_prev[0];
        for (int i = 1; i < n-1; ++i)
        {
            x_now[i] = d[i];
            x_now[i] -= a[i] * x_now[i - 1];
            x_now[i] -= c[i] * x_prev[i + 1];
            x_now[i] *= w;
            x_now[i] /= b[i];
            x_now[i] += (1 - w) * x_prev[i];
        }
        x_now[n - 1] = d[n - 1] - a[n - 1] * x_now[n - 2];
        x_now[n-1] *= w;
        x_now[n-1] /= b[n-1];
        x_now[n-1] += (1 - w) * x_prev[n-1];
        ++iter;
    } while (norm(x_now - x_prev) > EPS && iter <= max_iter);

    vector<double> an_sol(n, 2);
    for (int i = 0; i < n; ++i)
        an_sol[i] -= (i+1) % 2;
    //for (int i = 0; i < n; ++i)
    //cout << x_now[i] << endl;

    //Relax_log_info.C_norm = C_norm;
    //Relax_log_info.aprior = aprior_iters;
    //cout << "Число итераций = " << iter << endl;
    //cout << "Достигнутая точность " << norm(x_now - an_sol) << endl;
    //Relax_log_info.error_vector = error_vector();
    //Relax_log_info.error = vec_norm(Relax_log_info.error_vector);

    return x_now;
}

// Случай 2 (квазилинейное уравнение)
bool IterationScheme(double tau, double h, double sigma, PDE_data test, std::string filename="ImpScheme"){

    // Физические параметры
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

    // Создание файла
    std::string path = "./OutputData/" + filename;
    std::ofstream fpoints(path);
    std::cout << "log[INFO]: Starting ExplicitScheme" << std::endl;
    std::cout << "log[INFO]: Opening a file \"" << filename << "\" to write..." << std::endl;
    if (fpoints.is_open()) {

        double t_i = t_0;
        std::vector<double> state_i = state_0;
        writeVectorToFile(fpoints, t_i, state_i);
        double x_i = x_0;

        // Эволюция системы во времени
        for(int j = 0; j <= num_time_steps; ++j) {
            t_i += tau;
            for(int s = 0; s<3; ++s) {
                // Граничные условия слева
                // 1-го рода
                if (!test.G_left_type) {
                    Cs[0] = -1.;
                    Bs[0] = 0.;
                    As[0] = 0.;
                    Fs[0] = -state_0[0];
                }
                    // 2-го рода
                else {
                    double a0 = a(test.K_ptr, x_0 + h, x_0, state_i[1], state_i[0]);
                    double w0 = w(a0, state_i[1], state_i[0], h);
                    double kappa = sigma * a0 / h / (c * rho * h / (2 * tau) + sigma * a0 / h);
                    double mu = (c * rho * state_i[0] * h / (2 * tau) + sigma * test.G_left(t_i) +
                                 (1 - sigma) * (test.G_left(t_i - tau) + w0)) /
                                (c * rho * h / (2 * tau) + sigma * a0 / h);
                    Cs[0] = -1.;
                    Bs[0] = -kappa;
                    As[0] = 0;
                    Fs[0] = -mu;
                }

                // Граничные условия справа
                // 1-го рода
                if (!test.G_right_type) {
                    Bs[num_space_steps] = 0.;
                    As[num_space_steps] = 0.;
                    Cs[num_space_steps] = -1.;
                    Fs[num_space_steps] = -state_0[num_space_steps];
                }
                    // 2-го рода
                else {
                    double am = a(test.K_ptr, X, X - h, state_i[num_space_steps], state_i[num_space_steps-1]);
                    double wn = w(am, state_i[num_space_steps], state_i[num_space_steps - 1], h);
                    double denom = c * rho * h / (2 * tau) + sigma * am / h;
                    double kappa = sigma * am / h / denom;
                    double mu = (c * rho * state_i[num_space_steps] * h / (2 * tau) + sigma * test.G_right(t_i) +
                                 (1 - sigma) * (test.G_right(t_i - tau) - wn)) / denom;
                    Cs[num_space_steps] = -1.;
                    Bs[num_space_steps] = 0.;
                    As[num_space_steps] = -kappa;
                    Fs[num_space_steps] = -mu;
                }

                // Обход пространства
                for (int i = 1; i < num_space_steps; ++i) {
                    x_i += h;
                    double a_i = a(test.K_ptr, x_i, x_i - h, state_i[i], state_i[i-1]);
                    double a_ip = a(test.K_ptr, x_i + h, x_i, state_i[i + 1], state_i[i]);
                    As[i] = sigma / h * a_i;
                    Bs[i] = sigma / h * a_ip;
                    Cs[i] = (As[i] + Bs[i] + c * rho * h / tau);
                    Fs[i] = (c * rho * h / tau * state_i[i] +
                             (1 - sigma) *
                             (w(a_ip, state_i[i + 1], state_i[i], h) - w(a_i, state_i[i], state_i[i - 1], h)));

                }
                // Получение нового состояния системы
                //state_i = TripleBigRelaxSolve(As, Cs, Bs, Fs, state_i);
                state_i = TridiagonalMatrixAlgorithm(As, Cs, Bs, Fs);
                // Запись в файл
            }

            writeVectorToFile(fpoints, t_i, state_i);
        }
        fpoints.close();
        return true;
    } else {
        std::cout << "log[ERROR]: Couldn't open or create a file" << std::endl;
        return false;
    }
};


//Для определения числа итераций "до сходимости"
bool infoIterationScheme(double tau, double h, double sigma, PDE_data test, std::string filename="ImpScheme", double EPS=1e-12){

    // Физические параметры
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

    // Создание файла
    std::string path = "./OutputData/" + filename;
    std::ofstream fpoints(path);
    std::cout << "log[INFO]: Starting ExplicitScheme" << std::endl;
    std::cout << "log[INFO]: Opening a file \"" << filename << "\" to write..." << std::endl;
    if (fpoints.is_open()) {

        double t_i = t_0;
        std::vector<double> state_i = state_0;
        std::vector<double> state_ipp = state_i;
        writeVectorToFile(fpoints, t_i, state_i);
        double x_i = x_0;
        int iter_counter = 0;
        // Эволюция системы во времени
        for(int j = 0; j <= num_time_steps; ++j) {
            t_i += tau;
            iter_counter = 0;
            do {
                state_i = state_ipp;
                // Граничные условия слева
                // 1-го рода
                if (!test.G_left_type) {
                    Cs[0] = -1.;
                    Bs[0] = 0.;
                    As[0] = 0.;
                    Fs[0] = -state_0[0];
                }
                    // 2-го рода
                else {
                    double a0 = a(test.K_ptr, x_0 + h, x_0, state_i[1], state_i[0]);
                    double w0 = w(a0, state_i[1], state_i[0], h);
                    double kappa = sigma * a0 / h / (c * rho * h / (2 * tau) + sigma * a0 / h);
                    double mu = (c * rho * state_i[0] * h / (2 * tau) + sigma * test.G_left(t_i) +
                                 (1 - sigma) * (test.G_left(t_i - tau) + w0)) /
                                (c * rho * h / (2 * tau) + sigma * a0 / h);
                    Cs[0] = -1.;
                    Bs[0] = -kappa;
                    As[0] = 0;
                    Fs[0] = -mu;
                }
                // Граничные условия справа
                // 1-го рода
                if (!test.G_right_type) {
                    Bs[num_space_steps] = 0.;
                    As[num_space_steps] = 0.;
                    Cs[num_space_steps] = -1.;
                    Fs[num_space_steps] = -state_0[num_space_steps];
                }
                    // 2-го рода
                else {
                    double am = a(test.K_ptr, X, X - h, state_i[num_space_steps],state_i[num_space_steps-1]);
                    double wn = w(am, state_i[num_space_steps], state_i[num_space_steps - 1], h);
                    double denom = c * rho * h / (2 * tau) + sigma * am / h;
                    double kappa = sigma * am / h / denom;
                    double mu = (c * rho * state_i[num_space_steps] * h / (2 * tau) + sigma * test.G_right(t_i) +
                                 (1 - sigma) * (test.G_right(t_i - tau) - wn)) / denom;
                    Cs[num_space_steps] = -1.;
                    Bs[num_space_steps] = 0.;
                    As[num_space_steps] = -kappa;
                    Fs[num_space_steps] = -mu;
                }

                // Обход пространства
                for (int i = 1; i < num_space_steps; ++i) {
                    x_i += h;
                    double a_i = a(test.K_ptr, x_i, x_i - h, state_i[i], state_i[i-1]);
                    double a_ip = a(test.K_ptr, x_i + h, x_i, state_i[i + 1],state_i[i]);
                    As[i] = sigma / h * a_i;
                    Bs[i] = sigma / h * a_ip;
                    Cs[i] = (As[i] + Bs[i] + c * rho * h / tau);
                    Fs[i] = (c * rho * h / tau * state_i[i] +
                             (1 - sigma) *
                             (w(a_ip, state_i[i + 1], state_i[i], h) - w(a_i, state_i[i], state_i[i - 1], h)));

                }
                // Получение нового состояния системы
                //state_ipp = TripleBigRelaxSolve(As, (-1.)*Cs, Bs, (-1.)*Fs, state_i);
                state_ipp = TridiagonalMatrixAlgorithm(As, Cs, Bs, Fs);
                ++iter_counter;
            } while(norm(state_ipp + (-1.)*state_i) >= EPS);
            std::cout << "Iterations on time-step" << t_i << " is " << iter_counter << std::endl;
            // Запись в файл
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

//*****************************************OLEG's sandbox**********************************//
// Моя переделка функций Вани:

//
//std::vector<double> progonka(std::vector<double>& a, std::vector<double>& b, std::vector<double>& c, std::vector<double>& d) {
//
//    int n = b.size();
//    std::vector<double> alph(n, 0);
//    std::vector<double> beth(n, 0);
//    std::vector<double> solve(n, 0);
//    double tmp;
//
//    //vectinput();
//    alph[0] = 0;
//    beth[0] = 0;
//    tmp = b[0];
//    alph[1] = c[0] / tmp;
//    beth[1] = d[0] / tmp;
//
//    //std::cout << alph[1] << " " << beth[1] << "\n";
//    for (int i = 1; i < n - 1; i++) {
//        tmp = b[i] - a[i] * alph[i];
//        alph[i + 1] = c[i] / tmp;//(b[i] - a[i] * alph[i]);
//        beth[i + 1] = (d[i] + a[i] * beth[i]) / tmp; // (b[i] - a[i] * alph[i]);
//    }
//
//    solve[n - 1] = (d[n - 1] + a[n - 1] * beth[n - 1]) / (b[n - 1] - a[n - 1] * alph[n - 1]);
//
//    for (int i = n - 2; i >= 0; i--) {
//        solve[i] = alph[i + 1] * solve[i + 1] + beth[i + 1];
//    }
//
//    return solve;
//}
//
///* Функция для решения PDE в случае 1 (по методичке) */
//bool SolvePDE_1(PDE_data test, double tau, double h, double sigma, std::string filename){
//
//
//    double t_0 = 0; // Начальное время
//    double x_0 = 0; // Начальное условие?
//
//    // Количество шагов по времени и пространству
//    int num_time_steps = static_cast<int>((test.T-t_0) / tau);
//    int num_space_steps = static_cast<int>((test.L - x_0) / h);
//
//    // TODO: брать граничное условие из теста (здесь начальная температура)
//    double u_0 = 10;
//
//    // Инициализация начального состояния
//    std::vector<double> state_0 = init_state(num_space_steps, u_0); //TODO: расширить init_state
//    std::vector<double> As(num_space_steps, 0);
//    std::vector<double> Cs(num_space_steps, 0);
//    std::vector<double> Bs(num_space_steps, 0);
//    std::vector<double> Fs(num_space_steps, 0);
//
//    std::string path = "./OutputData/" + filename + ".txt";
//    std::ofstream fpoints(path);
//    std::cout << "log[INFO]: Starting ExplicitScheme" << std::endl;
//    std::cout << "log[INFO]: Opening a file \"" << filename << "\" to write..." << std::endl;
//
//    if (fpoints.is_open()) {
//
//        double t_i = t_0;
//        std::vector<double> state_i = state_0;
//        int ind = 0;
//        writeVectorToFile(fpoints, t_i, state_i);
//        double x_i = x_0;
//        Cs[0] = 1;
//        Bs[0] = 0;
//        As[0] = 0;
//        Fs[0] = u_0;
//        Bs[num_space_steps-1] = 0;
//        As[num_space_steps-1] = 0;
//        Cs[num_space_steps-1] = 1;
//        Fs[num_space_steps-1] = u_0;
//
//        for (int j = 0; j < num_time_steps; ++j) {
//            t_i += tau;
//
//            for (int i = 1; i < num_space_steps - 1; ++i) {
//                x_i += h;
//
//                //double a_i = a(test.K_ptr, x_i, x_i - h);
//                double a_i = 0.5 * (test.K(x_i) + test.K(x_i - h));
//
//                //double a_ip = a(test.K_ptr, x_i + h, x_i);
//                double a_ip = 0.5 * (test.K(x_i + h) + test.K(x_i - h));
//
//                As[i] = sigma / h * a_i;
//                Bs[i] = sigma / h * a_ip;
//                Cs[i] = As[i] + Bs[i] + test.c * test.rho * h / tau;
//                Fs[i] = test.c * test.rho * h / tau * state_i[i]
//                        + (1 - sigma) * (w(a_ip, state_i[i + 1],
//                                           state_i[i], h) - w(a_i, state_i[i], state_i[i - 1], h));
//            }
//
//            Cs = (-1.) * Cs;
//            Fs = (-1.) * Fs;
//            state_i = progonka(As, Cs, Bs, Fs);
//            writeVectorToFile(fpoints, t_i, state_i);
//        }
//
//        fpoints.close();
//        return true;
//
//    } else {
//        std::cout << "log[ERROR]: Couldn't open or create a file" << std::endl;
//        return false;
//    }
//}
