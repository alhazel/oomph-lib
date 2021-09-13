// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2021 Matthias Heil and Andrew Hazel
// LIC//
// LIC// This library is free software; you can redistribute it and/or
// LIC// modify it under the terms of the GNU Lesser General Public
// LIC// License as published by the Free Software Foundation; either
// LIC// version 2.1 of the License, or (at your option) any later version.
// LIC//
// LIC// This library is distributed in the hope that it will be useful,
// LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
// LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// LIC// Lesser General Public License for more details.
// LIC//
// LIC// You should have received a copy of the GNU Lesser General Public
// LIC// License along with this library; if not, write to the Free Software
// LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// LIC// 02110-1301  USA.
// LIC//
// LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
// LIC//
// LIC//====================================================================

// Oomph-lib includes
#include "generic.h"

using namespace std;
using namespace oomph;

void print_complex_vector(Vector<complex<double>>& x)
{
  for (unsigned i = 0; i < x.size(); i++)
  {
    cout << x[i] << endl;
  }
}

void print_cr_complex_matrix(CRComplexMatrix& M)
{
  for (unsigned i = 0; i < M.nrow(); i++)
  {
    for (unsigned j = 0; j < M.ncol(); j++)
    {
      cout << M(i, j) << ", ";
    }
    cout << endl;
  }
}

int main()
{
  unsigned long n = 2;
  unsigned long m = 4;

  Vector<complex<double>> values(4);
  values[0] = complex<double>(1, 1);
  values[1] = complex<double>(2, 1);
  values[2] = complex<double>(3, 1);
  values[3] = complex<double>(4, 1);

  Vector<int> column_index(4);
  column_index[0] = 0;
  column_index[1] = 2;
  column_index[2] = 1;
  column_index[3] = 3;

  Vector<int> row_start(3); // length = n+1
  row_start[0] = 0;
  row_start[1] = 2;
  row_start[2] = 4; // Fictitous index = NNz

  Vector<complex<double>> values_square(5);
  values_square[0] = complex<double>(1, 1);
  values_square[1] = complex<double>(2, 1);
  values_square[2] = complex<double>(3, 1);
  values_square[3] = complex<double>(4, 1);
  values_square[4] = complex<double>(5, 1);

  Vector<int> column_index_square(5);
  column_index_square[0] = 0;
  column_index_square[1] = 2;
  column_index_square[2] = 1;
  column_index_square[3] = 0;
  column_index_square[4] = 3;

  Vector<int> row_start_square(5); // length = n+1
  row_start_square[0] = 0;
  row_start_square[1] = 2;
  row_start_square[2] = 3;
  row_start_square[3] = 4;
  row_start_square[4] = 5; // Fictitous index = NNz

  Vector<std::complex<double>> x(2);
  x[0] = complex<double>(2, 2);
  x[1] = complex<double>(-2, 3);
  x[2] = complex<double>(0.5, 2.2);
  x[3] = complex<double>(-1.4, 3);

  Vector<std::complex<double>> rhs(4);
  rhs[0] = complex<double>(1, 0);
  rhs[1] = complex<double>(0, 1);
  rhs[2] = complex<double>(2, 1);
  rhs[3] = complex<double>(-3, 1);

  Vector<std::complex<double>> soln(4);

  // test default constructor
  CRComplexMatrix matrix_default();
  // cout << matrix_default.nrow() << endl;
  // cout << matrix_default.ncol() << endl;

  // test full matrix constructor
  CRComplexMatrix matrix(values, column_index, row_start, n, m);
  cout << matrix.nrow() << endl;
  cout << matrix.ncol() << endl;
  print_cr_complex_matrix(matrix);

  CRComplexMatrix matrix_square(values_square, column_index_square, row_start_square, 4, 4);
  cout << matrix_square.nrow() << endl;
  cout << matrix_square.ncol() << endl;
  print_cr_complex_matrix(matrix_square);

  // test LU decomposition
  cout << matrix_square.ludecompose() << endl;

  // test residual
  matrix_square.lubksub(rhs);
  print_complex_vector(rhs);

  // test multiply
  matrix_square.multiply(x, soln);
  print_complex_vector(soln);

  // test multiply transposed
  matrix_square.multiply_transpose(x, soln);
  print_complex_vector(soln);

  return (EXIT_SUCCESS);
}