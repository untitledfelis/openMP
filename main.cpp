#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <omp.h>
#include <cmath>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <chrono>
#include <random>
#include <getopt.h>
#include <ctime>


using namespace std;

class Matrix {
public:
    Matrix operator*(const Matrix& matrix);
    Matrix(size_t rows, size_t cols)
      : m_rows(rows)
      , m_cols(cols)
      , m_data(rows * cols)
  {}

  size_t rows() const {
    return m_rows;
  }

  size_t cols() const {
    return m_cols;
  }

  double &operator()(size_t row, size_t col) {
    return m_data[row * m_cols + col];
  }

  const double &operator()(size_t row, size_t col) const  {
    return m_data[row * m_cols + col];
  }
private:
  size_t m_rows;
  size_t m_cols;
  vector<double> m_data;
};
Matrix Matrix::operator*(const Matrix& matrix){
    Matrix res(m_rows, matrix.m_cols);

    if (m_cols == matrix.m_rows) {

       // #pragma omp parallel for
        for (int i = 0; i < res.m_rows; ++i)
            for (int j = 0; j < res.m_cols; ++j)
                for (int k = 0; k < m_rows; ++k)
                    res.m_data[i * m_rows + j] += m_data[i * m_rows + k] * matrix.m_data[k * m_rows + j];
    }
    else
        cout << "\nWrong dimensions\n";
    return res;

}
double F(double x, double y){
    return 0;
}
double G(double x, double y){
    if (x==0) return 1-2*y;
    if (x==1) return -1+2*y;
    if (y==0) return 1-2*x;
    if (y==1) return -1+2*x;
}
void UOutput(Matrix &u) {
    cout << fixed << setprecision(3);
    for (size_t i = 0; i <11; i++) {
        for (int j = 0; j <11; j++)
        cout << setw(7) << u(i,j);
        cout << endl;
    }
}
void Init(Matrix &f,Matrix &u, double h){
    for (size_t i = 0; i < f.rows()-2; i++) {
        for (size_t j = 0; j < f.cols()-2; j++)
            f(i,j) = F((i+1)*h, (j+1)*h);
    }
    for (int i = 1; i < u.rows()-1; i++) {
        for (int j = 1; j < u.rows()-1; j++) u(i,j) = 0;
        u(i,0) = G(i*h, 0);
        u(i,u.rows()-1) = G(i*h, (u.rows()-1)*h);
    }
    for (int j = 0; j < u.rows(); j++) {
        u(0,j) = G(0,j*h);
        u(u.rows()-1,j) = G((u.rows()-1)*h,j*h);
    }
    UOutput(u);
}
void Calc(Matrix &f,Matrix &u, double h, double eps, int &IterCnt){
    double max;
    do{
        IterCnt++;
        max = 0;
        for (int i = 1; i < u.rows()-1; i++)
            for (int j = 1; j < u.cols()-1; j++){
                double u0 = u(i,j);
                u(i,j) = 0.25*(u(i-1,j) + u(i+1,j)+ u(i,j-1) + u(i,j+1) - h*h*f(i-1,j-1));
                double d = abs(u(i,j)-u0);
                if (d > max) max = d;
            }
    }
    while (max > eps);
}

void OMPCalc3(Matrix &f,Matrix &u, double h, double eps, int &IterCnt){
    double max;
    int N = f.rows();
    double *mx = new double[N];
    IterCnt = 0;
    do{
    IterCnt++;
    #pragma omp parallel for shared(u,N,max) num_threads(1)
    for (int k = 1; k < N-1; k++){
        mx[k] = 0;
        //schedule()// num_threads(1)
        for (int i = 1; i < k-1; i++){
            int j = k + 1 - i;
            double u0 = u(i,j);
            u(i,j) = 0.25*(u(i-1,j) + u(i+1,j)+ u(i,j-1) + u(i,j+1) - h*h*f(i-1,j-1));
            double d = abs(u(i,j)-u0);
            if (d > mx[i])
            mx[i] = d;
        }
    }
    #pragma omp parallel for shared(u,N,max) num_threads(1) //schedule(static, 1)
    for (int k = N-1; k > 0; k--) {
        //#pragma omp parallel for shared(u,N,max) num_threads(2)////
        for (int i = N-k+1; i < N-1; i++) {
            int j = 2*N - k - i + 1;
            double u0 = u(i,j);
            u(i,j) = 0.25*(u(i-1,j) + u(i+1,j)+ u(i,j-1) + u(i,j+1) - h*h*f(i-1,j-1));
            double d = abs(u(i,j)-u0);
            if (d > mx[i])
            mx[i] = d;
        }
    }
    max = 0;
   #pragma omp parallel for shared(N,max) num_threads(1)
    for (int i = 1; i < 3; i++) {
        double d = 0;
        for (int j = i; j < N+1; j+=2)
            if (d < mx[j])
                d = mx[j];
            if (d > max)
                #pragma omp critical
            if (d > max)
            max = d;
        }
    }
    while (max > eps);
}

int main(){
    size_t N = 1001;
    double eps = 0.001;
    double h = 0.001;
    int IterCnt1 = 0;
    int IterCnt2 = 0;

    size_t rows = N;
    size_t cols = N;

    Matrix f(rows, cols);
    Matrix u(rows, cols);


    Init(f, u, h);
    double tt= omp_get_wtime();
    //Calc(f, u, h, eps, IterCnt1);
    tt = omp_get_wtime() - tt;
    double tt1= omp_get_wtime();
    OMPCalc3(f, u, h, eps, IterCnt2);
    tt1 = omp_get_wtime() - tt1;
    cout << "Time1 = " << tt << " IterCnt = " << IterCnt1 << endl;
    cout << "Time2 = " << tt1 << " IterCnt = " << IterCnt2 << endl;
    UOutput(u);
    cin.get();
    return 0;
}
