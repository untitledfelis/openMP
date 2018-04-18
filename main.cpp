#include <iostream>
#include <sstream>
#include <vector>
#include <iostream>
#include <string>

#include <algorithm>
#include <chrono>
#include <math.h>
#include <random>
#include <getopt.h>
#include <ctime>
#include <omp.h>

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

//////////////////////////////////////////////
Matrix Matrix::operator*(const Matrix& matrix)
{
    Matrix res(m_rows, matrix.m_cols);

    if (m_cols == matrix.m_rows) {

        #pragma omp parallel for
        for (int i = 0; i < res.m_rows; ++i)
            for (int j = 0; j < res.m_cols; ++j)
                for (int k = 0; k < m_rows; ++k)
                    res.m_data[i * m_rows + j] += m_data[i * m_rows + k] * matrix.m_data[k * m_rows + j];
    }
    else
        cout << "\nWrong dimensions\n";
    return res;

}

//////////////////////////////////////////////
string toString(const Matrix& matrix)
{
	stringstream ss;
	for (size_t i = 0; i < matrix.rows(); ++i) {
		for (size_t j = 0; j < matrix.cols(); ++j) {
			ss << matrix(i, j) << " ";
		}
		ss << endl;
	}
	return ss.str();
}

//////////////////////////////////////
int main() {

  for (int th=4; th<5; th++) {
        for (int s=1000; s<2001; s+=500) {
                double sumTime=0;
                for (int b=0; b<10; b++) {
  double t1 = omp_get_wtime();
  omp_set_num_threads(th);

   size_t rows = s;
   size_t cols = s;

  Matrix A(rows, cols);
  Matrix B(rows, cols);
  Matrix C(rows, cols);

  random_device rd;
  #pragma omp parallel
 {
    uint32_t seed;
    #pragma omp critical
    {
      seed = rd();
    }
    mt19937 gen(seed);
    uniform_real_distribution<double> dist(0.0, 1.0);

  #pragma omp for
  for(size_t i = 0; i < A.rows(); ++i) {
    for(size_t j = 0; j < A.cols(); ++j) {

       A(i, j) = dist(gen);
       B(i, j) = dist(gen);
       C(i, j) = 0;
    }
  }
 }

  //cout << toString(A) << endl;
  //cout << toString(B) << endl;
  C = A * B;
  //cout << toString(C) << endl;
  double t2 = omp_get_wtime();
  sumTime += (t2 - t1);
  cout << "Threads: " << omp_get_max_threads();
  cout << " Size: "<<s;
  cout << " Worked time: " << (t2 - t1)<< "\n";
}
  double midTime = sumTime/10;
  cout<< "Middle time: "<<midTime<<endl;
}}
  return 0;
}

