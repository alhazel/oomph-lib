// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2022 Matthias Heil and Andrew Hazel
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
#ifndef OOMPH_SOLID_PRECONDITIONERS_HEADER
#define OOMPH_SOLID_PRECONDITIONERS_HEADER


// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif


// oomphlib headers
#include "../generic/matrices.h"
#include "../generic/assembly_handler.h"
#include "../generic/problem.h"
#include "../generic/block_preconditioner.h"
#include "../generic/preconditioner.h"
#include "../generic/SuperLU_preconditioner.h"
#include "../generic/matrix_vector_product.h"


namespace oomph
{
  //===========================================================================
  /// The least-squares commutator (LSC; formerly BFBT)
  /// preconditioner. It uses blocks corresponding to the displacement/position
  /// and pressure unknowns, i.e. there are a total of 2x2 blocks,
  /// and all displacement/position components are treated as a
  /// single block of unknowns.
  ///
  /// Here are the details: An "ideal" preconditioner
  /// would solve the saddle point system
  /// \f[ \left( \begin{array}{cc} {\bf F} & {\bf G} \\ {\bf D} & {\bf 0} \end{array} \right) \left( \begin{array}{c} {\bf z}_u \\ {\bf z}_p \end{array} \right) = \left( \begin{array}{c} {\bf r}_u \\ {\bf r}_p \end{array} \right) \f]
  /// where \f$ {\bf F}\f$,  \f$ {\bf G} \f$, and \f$ {\bf D}\f$ are
  /// the blocks that arise in the Jacobian of the pressure-based
  /// equations of linear and nonlinear elasticity (with dofs in order
  /// of displacement/position and pressure).
  /// The use of this preconditioner would ensure the convergence
  /// of any iterative linear solver in a single iteration but its
  /// application is, of course, exactly as expensive as a direct solve.
  /// The LSC/BFBT preconditioner replaces the exact Jacobian by
  /// a block-triangular approximation
  /// \f[ \left( \begin{array}{cc} {\bf F} & {\bf G} \\ {\bf 0} & -{\bf M}_s \end{array} \right) \left( \begin{array}{c} {\bf z}_u \\ {\bf z}_p \end{array} \right) = \left( \begin{array}{c} {\bf r}_u \\ {\bf r}_p \end{array} \right), \f]
  /// where \f${\bf M}_s\f$ is an approximation to the pressure
  /// Schur-complement \f$ {\bf S} = {\bf D} {\bf F}^{-1}{\bf G}. \f$
  /// This system can be solved in two steps:
  /// -# Solve the second row for \f$ {\bf z}_p\f$ via
  /// \f[ {\bf z}_p = - {\bf M}_s^{-1} {\bf r}_p \f]
  /// -# Given \f$ {\bf z}_p \f$ , solve the first row for \f$ {\bf z}_u\f$ via
  /// \f[ {\bf z}_u = {\bf F}^{-1} \big( {\bf r}_u - {\bf G} {\bf z}_p \big) \f].
  /// In the LSC/BFBT preconditioner, the action of the inverse pressure
  /// Schur complement
  /// \f[ {\bf z}_p = - {\bf M}_s^{-1} {\bf r}_p \f]
  /// is approximated by
  /// \f[ {\bf z}_p = - \big({\bf D} \widehat{\bf Q}^{-1}{\bf G} \big)^{-1} \big({\bf D} \widehat{\bf Q}^{-1}{\bf F} \widehat{\bf Q}^{-1}{\bf G}\big) \big({\bf D} \widehat{\bf Q}^{-1}{\bf G} \big)^{-1} {\bf r}_p, \f]
  /// where  \f$ \widehat{\bf Q} \f$ is the diagonal of the
  /// displacement/position mass matrix. The evaluation of this expression
  /// involves two linear solves involving the matrix
  /// \f[ {\bf P} = \big({\bf D} \widehat{\bf Q}^{-1}{\bf G} \big) \f]
  /// which has the character of a matrix arising from the discretisation
  /// of a Poisson problem on the pressure space. We also have
  /// to evaluate matrix-vector products with the matrix
  /// \f[ {\bf E}={\bf D}\widehat{\bf Q}^{-1}{\bf F}\widehat{\bf Q}^{-1}{\bf G} \f]
  /// Details of the theory can be found in "Finite Elements and
  /// Fast Iterative Solvers with Applications in Incompressible Fluid
  /// Dynamics" by Howard C. Elman, David J. Silvester, and Andrew J. Wathen,
  /// published by Oxford University Press, 2006.
  ///
  /// In our implementation of the preconditioner, the linear systems
  /// can either be solved "exactly", using SuperLU (in its incarnation
  /// as an exact preconditioner; this is the default) or by any
  /// other Preconditioner (inexact solver) specified via the access functions
  /// \code
  /// PressureBasedSolidLSCPreconditioner::set_f_preconditioner(...)
  /// \endcode
  /// or
  /// \code
  /// PressureBasedSolidLSCPreconditioner::set_p_preconditioner(...)
  /// \endcode
  //===========================================================================
  class PressureBasedSolidLSCPreconditioner
    : public BlockPreconditioner<CRDoubleMatrix>
  {
  public:
    /// Constructor - sets defaults for control flags
    PressureBasedSolidLSCPreconditioner()
      : BlockPreconditioner<CRDoubleMatrix>()
    {
      // Flag to indicate that the preconditioner has been setup
      // previously -- if setup() is called again, data can
      // be wiped.
      Preconditioner_has_been_setup = false;

      // By default we use SuperLU for both p and f blocks
      Using_default_p_preconditioner = true;
      Using_default_f_preconditioner = true;

      // resize the mesh pt
      // note: meaningless if subsidiary preconditioner
      this->set_nmesh(1);
      Solid_mesh_pt = 0;

      // Set default preconditioners (inexact solvers) -- they are
      // members of this class!
      P_preconditioner_pt = 0;
      F_preconditioner_pt = 0;

      // Flag to determine if mass matrix diagonal Q^{-1}
      // is used for scaling.
      P_matrix_using_scaling = true;

      // set Doc_time to false
      Doc_time = false;

      // null the off diagonal Block matrix pt
      Bt_mat_vec_pt = 0;

      // null the F matrix vector product helper
      F_mat_vec_pt = 0;

      // null the QBt matrix vector product pt
      QBt_mat_vec_pt = 0;

      // null the E matrix vector product helper
      E_mat_vec_pt = 0;

      // by default we do not form the E matrix (BQFQBt)
      Form_BFBt_product = false;
    }

    /// Destructor
    ~PressureBasedSolidLSCPreconditioner()
    {
      clean_up_memory();
    }

    /// Broken copy constructor
    PressureBasedSolidLSCPreconditioner(
      const PressureBasedSolidLSCPreconditioner&) = delete;

    /// Broken assignment operator
    // Commented out broken assignment operator because this can lead to a
    // conflict warning when used in the virtual inheritence hierarchy.
    // Essentially the compiler doesn't realise that two separate
    // implementations of the broken function are the same and so, quite
    // rightly, it shouts.
    /*void operator=(const PressureBasedSolidLSCPreconditioner&) = delete;*/

    /// Setup the preconditioner
    void setup();

    /// Apply preconditioner to Vector r
    void preconditioner_solve(const DoubleVector& r, DoubleVector& z);

    /// specify the mesh containing the mesh containing the
    /// block-preconditionable solid elements. The dimension of the
    /// problem must also be specified.
    void set_solid_mesh(Mesh* mesh_pt)
    {
      Solid_mesh_pt = mesh_pt;
    }

    /// Enable mass matrix diagonal scaling in the
    /// Schur complement approximation
    void enable_p_matrix_scaling()
    {
      P_matrix_using_scaling = true;
    }

    /// Enable mass matrix diagonal scaling in the
    /// Schur complement approximation
    void disable_p_matrix_scaling()
    {
      P_matrix_using_scaling = false;
    }

    /// Return whether the mass matrix is using diagonal
    /// scaling or not
    bool is_p_matrix_using_scaling() const
    {
      return P_matrix_using_scaling;
    }

    /// Function to set a new pressure matrix preconditioner (inexact solver)
    void set_p_preconditioner(Preconditioner* new_p_preconditioner_pt)
    {
      // If the default preconditioner has been used
      // clean it up now...
      if (Using_default_p_preconditioner)
      {
        delete P_preconditioner_pt;
      }
      P_preconditioner_pt = new_p_preconditioner_pt;
      Using_default_p_preconditioner = false;
    }

    /// Function to (re-)set pressure matrix preconditioner  (inexact
    /// solver) to SuperLU
    void set_p_superlu_preconditioner()
    {
      if (!Using_default_p_preconditioner)
      {
        P_preconditioner_pt = new SuperLUPreconditioner;
        Using_default_p_preconditioner = true;
      }
    }

    /// Function to set a new momentum matrix preconditioner (inexact solver)
    void set_f_preconditioner(Preconditioner* new_f_preconditioner_pt)
    {
      // If the default preconditioner has been used
      // clean it up now...
      if (Using_default_f_preconditioner)
      {
        delete F_preconditioner_pt;
      }
      F_preconditioner_pt = new_f_preconditioner_pt;
      Using_default_f_preconditioner = false;
    }

    /// Function to (re-)set momentum matrix preconditioner (inexact
    /// solver) to SuperLU
    void set_f_superlu_preconditioner()
    {
      if (!Using_default_f_preconditioner)
      {
        F_preconditioner_pt = new SuperLUPreconditioner;
        Using_default_f_preconditioner = true;
      }
    }


    /// Enable documentation of time
    void enable_doc_time()
    {
      Doc_time = true;
    }

    /// Disable documentation of time
    void disable_doc_time()
    {
      Doc_time = false;
    }

    /// If this function is called then:
    /// in setup(...) : BFBt is computed.
    /// in preconditioner_solve(...) : a single matrix vector product with
    /// BFBt is performed.
    void enable_form_BFBt_product()
    {
      Form_BFBt_product = true;
    }

    /// if this function is called  then:
    /// in setup(...) : the matrices B, F are assembled and stored
    /// (the default behaviour) .
    /// in preconditioner_solve(...) : a sequence of matrix vector products
    /// with B, F, and Bt is performed.
    /// (Note: in this discussion no scaling was considered but B and Bt
    ///  are replaced with BQ and QBt with scaling)
    void disable_form_BFBt_product()
    {
      Form_BFBt_product = false;
    }

    /// Helper function to delete preconditioner data.
    void clean_up_memory();

  private:
    // oomph-lib objects
    // -----------------

    // Pointers to preconditioner (=inexact solver) objects
    // -----------------------------------------------------
    /// Pointer to the 'preconditioner' for the pressure matrix
    Preconditioner* P_preconditioner_pt;

    /// Pointer to the 'preconditioner' for the F matrix
    Preconditioner* F_preconditioner_pt;

    /// flag indicating whether the default F preconditioner is used
    bool Using_default_f_preconditioner;

    /// flag indicating whether the default P preconditioner is used
    bool Using_default_p_preconditioner;

    /// Control flag is true if the preconditioner has been setup
    /// (used so we can wipe the data when the preconditioner is
    /// called again)
    bool Preconditioner_has_been_setup;

    /// Control flag is true if mass matrix diagonal scaling
    /// is used in the Schur complement approximation
    bool P_matrix_using_scaling;

    /// Helper function to assemble the diagonal of the
    /// mass matrix from the elemental contributions defined in
    /// PressureBasedSolidEquations<DIM>::get_mass_matrix_diagonal(...).
    CRDoubleMatrix* assemble_mass_matrix_diagonal();

    /// Boolean indicating whether the momentum system preconditioner
    /// is a block preconditioner
    bool F_preconditioner_is_block_preconditioner;

    /// Set Doc_time to true for outputting results of timings
    bool Doc_time;

    /// MatrixVectorProduct operator for F if BFBt is not to be formed.
    MatrixVectorProduct* F_mat_vec_pt;

    /// MatrixVectorProduct operator for QBt if BFBt is not to be formed.
    MatrixVectorProduct* QBt_mat_vec_pt;

    /// MatrixVectorProduct operator for Bt;
    MatrixVectorProduct* Bt_mat_vec_pt;

    /// MatrixVectorProduct operator for E (BFBt) if BFBt is to be formed.
    MatrixVectorProduct* E_mat_vec_pt;

    /// indicates whether BFBt should be formed or the component matrices
    /// should be retained.
    /// If true then:
    /// in setup(...) : BFBt is computed.
    /// in preconditioner_solve(...) : a single matrix vector product with
    /// BFBt is performed.
    /// if false then:
    /// in setup(...) : the matrices B, F are assembled and stored.
    /// in preconditioner_solve(...) : a sequence of matrix vector products
    /// with B, F, and Bt is performed.
    /// (Note: in this discussion no scaling was considered but B and Bt
    ///  are replaced with BQ and QBt with scaling)
    bool Form_BFBt_product;

    /// the pointer to the mesh of block preconditionable solid
    /// elements.
    Mesh* Solid_mesh_pt;
  };


  /// ////////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////////


  //============================================================================
  /// The exact solid preconditioner. This extracts 2x2 blocks
  /// (corresponding to the displacement/position and pressure unknowns)
  /// and uses these to
  /// build a single preconditioner matrix for testing purposes.
  /// Iterative solvers should converge in a single step if this is used.
  /// If it doesn't something is wrong in the setup of the block matrices.
  //=============================================================================
  template<typename MATRIX>
  class PressureBasedSolidExactPreconditioner
    : public BlockPreconditioner<MATRIX>
  {
  public:
    /// Constructor - do nothing
    PressureBasedSolidExactPreconditioner() : BlockPreconditioner<MATRIX>() {}


    /// Destructor - do nothing
    ~PressureBasedSolidExactPreconditioner() {}


    /// Broken copy constructor
    PressureBasedSolidExactPreconditioner(
      const PressureBasedSolidExactPreconditioner&) = delete;


    /// Broken assignment operator
    /*void operator=(const PressureBasedSolidExactPreconditioner&) = delete;*/


    /// Setup the preconditioner
    void setup();

    /// Apply preconditioner to r
    void preconditioner_solve(const Vector<double>& r, Vector<double>& z);

  protected:
    /// Preconditioner matrix
    MATRIX P_matrix;
  };

} // namespace oomph
#endif
