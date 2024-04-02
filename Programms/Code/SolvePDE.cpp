//
// Реализация функций для решения Уравнений в частных производных(УЧП)
//

#include "SolvePDE.h"
#include "algebra.h"


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

std::vector<double> ExplicitScheme(double tau, double h, double sigma){
    std::vector<double> As;
    std::vector<double> Cs;
    std::vector<double> Bs;
    std::vector<double> Fs;
};

/* Функция для решения PDE в случае 1 (по методичке) */
void SolvePDE_1(PDE_data test, double h, double tau, double sigma, std::string filename){

}
