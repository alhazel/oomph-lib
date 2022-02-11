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

/// Base eigenproblem element class used to generate the Jacobian and mass
/// matrix which correspond to the eigenproblem (J - lambda M) v = 0
class BaseEigenElement : public GeneralisedElement
{
public:
  BaseEigenElement() {}

  void set_size(const unsigned& n)
  {
    N_value = n;

    Data_index = add_internal_data(new Data(N_value));
  }

  void fill_in_contribution_to_jacobian_and_mass_matrix(
    Vector<double>& residuals,
    DenseMatrix<double>& jacobian,
    DenseMatrix<double>& mass_matrix) = 0;

protected:
  unsigned N_value;
  unsigned Data_index;
};

/// IdentityEigenElement generates an eigenproblem whose Jacobian and mass
/// matrices are equal to the identiy matrix
class IdentityEigenElement : public BaseEigenElement
{
public:
  void fill_in_contribution_to_jacobian_and_mass_matrix(
    Vector<double>& residuals,
    DenseMatrix<double>& jacobian,
    DenseMatrix<double>& mass_matrix)
  {
    for (unsigned i = 0; i < N_value; i++)
    {
      unsigned local_eqn = internal_local_eqn(Data_index, i);
      for (unsigned j = 0; j < N_value; j++)
      {
        unsigned local_unknown = internal_local_eqn(Data_index, j);
        if (i == j)
        {
          jacobian(local_eqn, local_unknown) += 1;
          mass_matrix(local_eqn, local_unknown) += 1;
        }
      }
    }
  }
};

/// AsymmetricEigenElement generates an eigenproblem whose Jacobian and mass
/// matrices are simple asymmetric examples
class AsymmetricEigenElement : public BaseEigenElement
{
public:
  void fill_in_contribution_to_jacobian_and_mass_matrix(
    Vector<double>& residuals,
    DenseMatrix<double>& jacobian,
    DenseMatrix<double>& mass_matrix)
  {
    for (unsigned i = 0; i < N_value; i++)
    {
      unsigned local_eqn = internal_local_eqn(Data_index, i);
      for (unsigned j = 0; j < N_value; j++)
      {
        unsigned local_unknown = internal_local_eqn(Data_index, j);
        if (i == j)
        {
          jacobian(local_eqn, local_unknown) += i + 1.0;
          mass_matrix(local_eqn, local_unknown) += 1.0;
        }
        else if (i > j)
        {
          jacobian(local_eqn, local_unknown) += 1.0;
          mass_matrix(local_eqn, local_unknown) += 1.0;
        }
      }
    }
  }
};

/// RandomAsymmetricEigenElement generates an eigenproblem whose Jacobian and
/// mass matrices have random integer elements between -128 and 127. Seed = 0.
class RandomAsymmetricEigenElement : public BaseEigenElement
{
public:
  void fill_in_contribution_to_jacobian_and_mass_matrix(
    Vector<double>& residuals,
    DenseMatrix<double>& jacobian,
    DenseMatrix<double>& mass_matrix)
  {
    unsigned seed = 0;
    srand(seed);
    for (unsigned i = 0; i < N_value; i++)
    {
      unsigned local_eqn = internal_local_eqn(Data_index, i);
      for (unsigned j = 0; j < N_value; j++)
      {
        unsigned local_unknown = internal_local_eqn(Data_index, j);
        jacobian(local_eqn, local_unknown) += rand() % 256 - 128;
        mass_matrix(local_eqn, local_unknown) += rand() % 256 - 128;
      }
    }
  }
};

/// RosserSymmetricEigenElement generates the classic Rosser eigenproblem
/// matrices. Size 8x8.
class RosserSymmetricEigenElement : public BaseEigenElement
{
public:
  void set_size(const unsigned& n)
  {
    /// Override the input argument as the Rosser matrix is fixed at 8x8
    N_value = 8;

    Data_index = add_internal_data(new Data(N_value));
  }

  void fill_in_contribution_to_jacobian_and_mass_matrix(
    Vector<double>& residuals,
    DenseMatrix<double>& jacobian,
    DenseMatrix<double>& mass_matrix)
  {
    int A[8][8] = {{611, 196, -192, 407, -8, -52, -49, 29},
                   {196, 899, 113, -192, -71, -43, -8, -44},
                   {-192, 113, 899, 196, 61, 49, 8, 52},
                   {407, -192, 196, 611, 8, 44, 59, -23},
                   {-8, -71, 61, 8, 411, -599, 208, 208},
                   {-52, -43, 49, 44, -599, 411, 208, 208},
                   {-49, -8, 8, 59, 208, 208, 99, -911},
                   {29, -44, 52, -23, 208, 208, -911, 99}};


    for (unsigned i = 0; i < N_value; i++)
    {
      unsigned local_eqn = internal_local_eqn(Data_index, i);
      for (unsigned j = 0; j < N_value; j++)
      {
        unsigned local_unknown = internal_local_eqn(Data_index, j);
        jacobian(local_eqn, local_unknown) += double(A[i][j]);
        if (i == j)
        {
          mass_matrix(local_eqn, local_unknown) += 1;
        }
      }
    }
  }
};

/// Eigenproblem class. Creates a mesh with a single eigenelement and assigns
/// the appropriate equation numbers.
template<class ELEMENT>
class Eigenproblem : public Problem
{
public:
  Eigenproblem(const unsigned& size)
  {
    this->mesh_pt() = new Mesh;

    // First element
    ELEMENT* el_pt = new ELEMENT;
    el_pt->set_size(size);
    this->mesh_pt()->add_element_pt(el_pt);

    // Second element
    el_pt = new ELEMENT;
    el_pt->set_size(size);
    this->mesh_pt()->add_element_pt(el_pt);

    assign_eqn_numbers();
  }

  ~Eigenproblem()
  {
    delete this->mesh_pt();
  }
};


template<class ELEMENT>
class SolveEigenProblemTest
{
public:
  SolveEigenProblemTest(EigenSolver* const& eigen_solver_pt,
                        const unsigned& N,
                        const unsigned& n_timing_loops,
                        DocInfo* const& doc_info_pt,
                        bool do_adjoint_problem)
    : Eigen_solver_pt(eigen_solver_pt),
      Matrix_size(N),
      N_timing_loops(n_timing_loops),
      Doc_info_pt(doc_info_pt),
      Do_adjoint_problem(do_adjoint_problem)
  {
    Problem_pt = new Eigenproblem<ELEMENT>(Matrix_size);
#ifdef OOMPH_HAS_MPI
    Problem_pt->distribute();
#endif

    // Set up additional arguments
    // Output all eigenvalues
    N_eval = 8;

    // Store outputs
    Vector<complex<double>> eval(N_eval);
    Vector<DoubleVector> eigenvector_real(N_eval);
    Vector<DoubleVector> eigenvector_imag(N_eval);

    // Start clock
    clock_t t_start = clock();
    for (unsigned i = 0; i < N_timing_loops; i++)
    {
      // Call solve_eigenproblem
      Eigen_solver_pt->solve_eigenproblem(Problem_pt,
                                          N_eval,
                                          eval,
                                          eigenvector_real,
                                          eigenvector_imag,
                                          Do_adjoint_problem);
    }
    // Stop clock
    clock_t t_end = clock();

    // Document duration
    double t_length = (double)(t_end - t_start) / CLOCKS_PER_SEC;
    ofstream timing_stream;
    timing_stream.open("timing_distributed.dat", ios_base::app);
    timing_stream << "test" << Doc_info_pt->number()
                  << ", time: " << t_length / double(N_timing_loops) << endl;
    timing_stream.close();

    // Document solution
    string filename = Doc_info_pt->directory() + "test" +
                      to_string(Doc_info_pt->number()) + ".dat";

    ofstream output_stream;
    output_stream.open(filename);
    for (unsigned i = 0; i < N_eval; i++)
    {
      output_stream << eval[i].real() << " " << eval[i].imag() << endl;
    }
    output_stream.close();

    Doc_info_pt->number()++;
  }

private:
  EigenSolver* Eigen_solver_pt;
  unsigned Matrix_size;
  unsigned N_eval;
  unsigned N_timing_loops;
  Problem* Problem_pt;
  DocInfo* Doc_info_pt;
  bool Do_adjoint_problem;
};

template<class ELEMENT>
class SolveEigenProblemLegacyTest
{
public:
  SolveEigenProblemLegacyTest(EigenSolver* const& eigen_solver_pt,
                              const unsigned& N,
                              const unsigned& n_timing_loops,
                              DocInfo* const& doc_info_pt,
                              bool do_adjoint_problem)
    : Eigen_solver_pt(eigen_solver_pt),
      Matrix_size(N),
      N_timing_loops(n_timing_loops),
      Doc_info_pt(doc_info_pt),
      Do_adjoint_problem(do_adjoint_problem)
  {
    Problem_pt = new Eigenproblem<ELEMENT>(Matrix_size);

    // Set up additional arguments
    // Output all eigenvalues
    N_eval = 8;

    // Store outputs
    Vector<complex<double>> eval(N_eval);
    Vector<DoubleVector> evec(N_eval);

    // Start clock
    clock_t t_start = clock();
    for (unsigned i = 0; i < N_timing_loops; i++)
    {
      // Call solve_eigenproblem
      Eigen_solver_pt->solve_eigenproblem_legacy(
        Problem_pt, N_eval, eval, evec, Do_adjoint_problem);
    }
    // Stop clock
    clock_t t_end = clock();

    // Document duration
    double t_length = (double)(t_end - t_start) / CLOCKS_PER_SEC;
    ofstream timing_stream;
    timing_stream.open("timing_distributed.dat", ios_base::app);
    timing_stream << "test" << Doc_info_pt->number()
                  << ", time: " << t_length / double(N_timing_loops) << endl;
    timing_stream.close();

    // Document solution
    string filename = Doc_info_pt->directory() + "test" +
                      to_string(Doc_info_pt->number()) + ".dat";

    ofstream output_stream;
    output_stream.open(filename);
    for (unsigned i = 0; i < N_eval; i++)
    {
      output_stream << eval[i].real() << " " << eval[i].imag() << endl;
    }
    output_stream.close();

    Doc_info_pt->number()++;
  }

private:
  EigenSolver* Eigen_solver_pt;
  unsigned Matrix_size;
  unsigned N_eval;
  unsigned N_timing_loops;
  Problem* Problem_pt;
  DocInfo* Doc_info_pt;
  bool Do_adjoint_problem;
};

void test_anasazi(const unsigned N,
                  const unsigned n_timing_loops,
                  DocInfo* doc_info_pt)
{
  EigenSolver* eigen_solver_pt = new ANASAZI;

  const bool do_adjoint_problem = false;

  SolveEigenProblemTest<IdentityEigenElement>(
    eigen_solver_pt, N, n_timing_loops, doc_info_pt, do_adjoint_problem);
  SolveEigenProblemTest<AsymmetricEigenElement>(
    eigen_solver_pt, N, n_timing_loops, doc_info_pt, do_adjoint_problem);
  SolveEigenProblemTest<RandomAsymmetricEigenElement>(
    eigen_solver_pt, N, n_timing_loops, doc_info_pt, do_adjoint_problem);

  SolveEigenProblemLegacyTest<IdentityEigenElement>(
    eigen_solver_pt, N, n_timing_loops, doc_info_pt, do_adjoint_problem);
  SolveEigenProblemLegacyTest<AsymmetricEigenElement>(
    eigen_solver_pt, N, n_timing_loops, doc_info_pt, do_adjoint_problem);
  SolveEigenProblemLegacyTest<RandomAsymmetricEigenElement>(
    eigen_solver_pt, N, n_timing_loops, doc_info_pt, do_adjoint_problem);

  delete eigen_solver_pt;
}

/// Main function. Apply solver tests to each eigensolver.
int main(int argc, char** argv)
{
// Want to test Trilinos if we have it, so we must initialise MPI
// if we have compiled with it
#ifdef OOMPH_HAS_MPI
  MPI_Helpers::init(argc, argv);

  // Number of times to repeat the operation for better timings
  const unsigned n_timing_loops = 2;

  // Matrix dimensions
  const unsigned N = 64;

  DocInfo* doc_info_pt = new DocInfo;
#ifdef OOMPH_HAS_TRILINOS
  doc_info_pt->set_directory("RESLT_anasazi_distributed/");
  doc_info_pt->number() = 0;
  ofstream timing_stream;
  timing_stream.open("timing_distributed.dat", ios_base::app);
  timing_stream << "ANASAZI" << endl;
  timing_stream.close();

  Anasazi::Use_temporary_code_for_andrew_legacy_version = true;

  test_anasazi(N, n_timing_loops, doc_info_pt);
#endif

  delete doc_info_pt;

  MPI_Helpers::finalize();
#endif

  return 0;
}
