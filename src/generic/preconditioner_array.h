// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2025 Matthias Heil and Andrew Hazel
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
// Include guards
#ifndef OOMPH_PRECONDITION_ARRAY_HEADER
#define OOMPH_PRECONDITION_ARRAY_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

// oomph-lib includes
#include "Vector.h"
#include "double_vector.h"
#include "matrices.h"
#include "preconditioner.h"


// Preconditioner array is only really possible and useful if we have MPI
// (and it uses a lot of MPI functions in the .cc file so changing that
// would be hard). So if we don't have MPI just define a dummy
// implementation that throws an error if you try to use it.
#ifndef OOMPH_HAS_MPI
namespace oomph
{
  class PreconditionerArray
  {
  public:
    void setup_preconditioners(Vector<CRDoubleMatrix*> matrix_pt,
                               Vector<Preconditioner*> prec_pt,
                               const OomphCommunicator* comm_pt)
    {
      throw OomphLibError("PreconditionerArray requires MPI",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }

    void solve_preconditioners(const Vector<DoubleVector>& r,
                               Vector<DoubleVector>& z)
    {
      throw OomphLibError("PreconditionerArray requires MPI",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
  };
} // namespace oomph

// Otherwise (if we have MPI do the real implementation)
#else

#include "mpi.h"


namespace oomph
{
  //=============================================================================
  /// PreconditionerArray -
  /// NOTE - first implementation, a number of assumptions / simplifications
  /// were made:
  /// 1. Only works with CRDoubleMatrices
  /// 2. The number of processors must be greater than the number of
  ///    preconditioners
  /// 3. Currently only very crude load balancing - each preconditioner will
  ///    be setup and applied with the same number of processors (or as near to
  ///    as possible to the same number of processors)
  /// 4. This class will, at the appropriate time, delete the all the
  ///    Preconditioners passed setup_preconditioners(...)
  /// 5. (but) Deletion of matrices passed to setup_preconditioners(...) is NOT
  ///    performed by this class
  /// 6. It is assumed that preconditioners do not require access to matrix
  ///    once setup(...) is called
  /// 7. The matrix on the subset of processors will be the same type
  ///    (distributed or global) as the matrix passed to
  ///    setup_preconditioners(...)
  /// 8. If the matrix is a distributed matrix - it will be assembled with
  ///    a uniform distribution on the subset of processors.
  //=============================================================================
  class PreconditionerArray
  {
  public:
    /// Constructor (empty)
    PreconditionerArray()
      : Preconditioner_pt(0),
        Global_communicator_pt(0),
        Local_communicator_pt(0)
    {
      Method = 0;
      Nprec = 0;
    };

    /// Broken copy constructor
    PreconditionerArray(const PreconditionerArray&) = delete;

    /// Broken assignment operator
    void operator=(const PreconditionerArray&) = delete;

    /// Destructor (empty)
    ~PreconditionerArray()
    {
      this->clean_up_memory();
    }

    /// Setup the preconditioners. Sets up each preconditioner in the
    /// array for the corresponding matrix in the vector matrix_pt.
    /// The number of preconditioners in the array is taken to be the length of
    /// prec_pt
    /// The preconditioners that are not used on this processor are deleted.
    void setup_preconditioners(Vector<CRDoubleMatrix*> matrix_pt,
                               Vector<Preconditioner*> prec_pt,
                               const OomphCommunicator* comm_pt);

    /// Applies each preconditioner to the corresponding vector in
    /// r and z
    void solve_preconditioners(const Vector<DoubleVector>& r,
                               Vector<DoubleVector>& z);

    /// Clean up memory.
    void clean_up_memory()
    {
      // delete the preconditioner pt
      delete Preconditioner_pt;
      Preconditioner_pt = 0;

      // delete the communicators
      delete Global_communicator_pt;
      Global_communicator_pt = 0;
      delete Local_communicator_pt;
      Local_communicator_pt = 0;

      // clear vectors
      First_row_for_proc.clear();
      Nrow_local_for_proc.clear();
      First_row_from_proc.clear();
      Nrow_local_from_proc.clear();
      First_proc_for_prec.clear();
      Nproc_for_prec.clear();

      // zero
      Color = 0;

#ifdef PARANOID
      // delete PARANOID check distribution pts
      for (unsigned i = 0; i < Nprec; i++)
      {
        delete Distribution_pt[i];
      }
      Distribution_pt.resize(0);
#endif
    }

    // access function to Method
    unsigned& method()
    {
      return Method;
    }

  private:
    /// helper method for computing the MPI_Isend and MPI_Irecv tags
    int compute_tag(const int& nproc,
                    const int& source,
                    const int& dest,
                    const int& type)
    {
      return source + (nproc * dest) + (nproc * nproc * type);
    }

    /// the number of preconditioner in the array
    unsigned Nprec;

    /// The pointer to the local preconditioner on this processor
    Preconditioner* Preconditioner_pt;

    /// The first_row component of the distribution of the processors
    /// over the preconditioners
    Vector<unsigned> First_proc_for_prec;

    /// The nrow_local component of the distribution of the processors
    /// over the preconditioners
    Vector<unsigned> Nproc_for_prec;

    /// Storage (indexed [i][j]) for the first row that will be sent
    /// from this processor to processor j for preconditioner i
    Vector<Vector<unsigned>> First_row_for_proc;

    /// Storage (indexed [i][j]) for the nrow_local that will be sent
    /// from this processor to processor j for preconditioner i
    Vector<Vector<unsigned>> Nrow_local_for_proc;

    /// Storage (indexed [i][j]) for the first row that will be received
    /// by this processor from processor j for preconditioner i
    Vector<Vector<unsigned>> First_row_from_proc;

    /// Storage (indexed [i][j]) for the nrow_local that will be
    /// received by this processor from processor j for preconditioner i
    Vector<Vector<unsigned>> Nrow_local_from_proc;

    /// the Color of this processor (or the preconditioner number)
    unsigned Color;

    /// pointer to the global communicator for this preconditioner array
    OomphCommunicator* Global_communicator_pt;

    /// Vector of communicators for the preconditioners
    OomphCommunicator* Local_communicator_pt;

#ifdef PARANOID
    // Vector of distribution of each preconditioner - for PARANOID checks only
    Vector<LinearAlgebraDistribution*> Distribution_pt;
#endif

    /// the communication method in the setup_preconditioners(...) method
    /// 1. Non-blocking Send with Blocking Recv
    /// 2. MPI_Datatypes with Non-blocking sends and receives
    unsigned Method;

  }; // PreconditionerArray
} // namespace oomph

// End of "if we have MPI"
#endif

// End of include guard
#endif
