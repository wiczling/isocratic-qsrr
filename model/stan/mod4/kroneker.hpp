#include <iostream>
#include <stan/math.hpp>
#include <Eigen/Dense>
#include <stdexcept>

namespace mod4_model_namespace {
// Kronecker product for matrices with potentially different scalar types
template <typename Derived1, typename Derived2>
Eigen::Matrix<typename stan::return_type<typename Derived1::Scalar, typename Derived2::Scalar>::type, Eigen::Dynamic, Eigen::Dynamic>
kronecker_prod(
  const Eigen::MatrixBase<Derived1>& A,
  const Eigen::MatrixBase<Derived2>& B,
  std::ostream* pstream__ = nullptr) {
  using T_return = typename stan::return_type<typename Derived1::Scalar, typename Derived2::Scalar>::type;
  
  int m = A.rows();
  int n = A.cols();
  int p = B.rows();
  int q = B.cols();
  
  // Check for valid dimensions
  if (m <= 0 || n <= 0 || p <= 0 || q <= 0) {
    throw std::invalid_argument("Kronecker product: Matrix dimensions must be positive");
  }
  
  // Initialize output matrix
  Eigen::Matrix<T_return, Eigen::Dynamic, Eigen::Dynamic> C(m * p, n * q);
  
  // Compute Kronecker product
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      C.block(i * p, j * q, p, q) = A(i, j) * B;
    }
  }
  
  return C;
}
}