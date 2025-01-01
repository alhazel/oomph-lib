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
#ifndef OOMPH_BIHARMONIC_PROBLEM_HEADER
#define OOMPH_BIHARMONIC_PROBLEM_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

#ifdef OOMPH_HAS_MPI
// mpi headers
#include "mpi.h"
#endif

// Generic C++ headers
#include <iostream>
#include <math.h>

// oomph-lib headers
#include "../generic/problem.h"
#include "../generic/hijacked_elements.h"
#include "../meshes/hermite_element_quad_mesh.template.h"
#include "../meshes/hermite_element_quad_mesh.template.cc"
#include "biharmonic_elements.h"
#include "biharmonic_flux_elements.h"


namespace oomph
{
  //=============================================================================
  /// Biharmonic Plate Problem Class - for problems where the load can be
  /// assumed to be acting normal to the surface of the plate and the
  /// deflections are small relative to the thickness of the plate. Developed
  /// for the topologically rectangular Hermite Element Mesh. Contains functions
  /// allowing the following boundary conditions to be applied (on a given
  /// edge):
  ///   + clamped :           u and du/dn imposed
  ///   + simply supported :  u and laplacian(u) imposed
  ///   + free :              laplacian(u) and dlaplacian(u)/dn imposed
  //=============================================================================
  template<unsigned DIM>
  class BiharmonicProblem : public Problem
  {
  public:
    /// Definition of a dirichlet boundary condition function pointer.
    /// Takes the position  along a boundary (s) in the macro element coordinate
    /// scheme and returns the value of the boundary condition at that point
    /// (u).
    typedef void (*DirichletBCFctPt)(const double& s, double& u);


    /// Definition of the Source Function.
    typedef void (*BiharmonicSourceFctPt)(const Vector<double>& x, double& f);

    /// Constructor
    BiharmonicProblem()
    {
      Bulk_element_mesh_pt = 0;
      Face_element_mesh_pt = 0;
    }

    /// Destructor. Delete the meshes
    virtual ~BiharmonicProblem()
    {
      delete Bulk_element_mesh_pt;
      delete Face_element_mesh_pt;
    };

    /// actions before solve, performs self test
    void actions_before_newton_solve()
    {
#ifdef PARANOID
      if (0 == self_test())
      {
        oomph_info << "self test passed" << std::endl;
      }
      else
      {
        oomph_info << "self test failed" << std::endl;
      }
#endif
    }

    /// action after solve
    void actions_after_newton_solve() {}

    /// documents the solution, and if an exact solution is provided,
    /// then the error between the numerical and exact solution is presented
    void doc_solution(
      DocInfo& doc_info,
      FiniteElement::SteadyExactSolutionFctPt exact_soln_pt = 0);

    /// Access function to the bulk element mesh pt
    Mesh* bulk_element_mesh_pt()
    {
      return Bulk_element_mesh_pt;
    }

  protected:
    /// builds the bulk mesh on a prescribed domain with a node spacing
    /// defined by spacing fn with n_x by n_y elements
    void build_bulk_mesh(
      const unsigned n_x,
      const unsigned n_y,
      TopologicallyRectangularDomain* domain_pt,
      HermiteQuadMesh<BiharmonicElement<2>>::MeshSpacingFnPtr spacing_fn = 0)
    {
      if (spacing_fn == 0)
      {
        Bulk_element_mesh_pt =
          new HermiteQuadMesh<Hijacked<BiharmonicElement<2>>>(
            n_x, n_y, domain_pt);
      }
      else
      {
        Bulk_element_mesh_pt =
          new HermiteQuadMesh<Hijacked<BiharmonicElement<2>>>(
            n_x, n_y, domain_pt, spacing_fn);
      }
      //   add_sub_mesh(Bulk_element_mesh_pt);
      //   build_global_mesh();
    }


    /// Build global mesh and assign equation numbers
    void build_global_mesh_and_assign_eqn_numbers()
    {
      add_sub_mesh(Bulk_element_mesh_pt);
      if (Face_element_mesh_pt != 0)
      {
        add_sub_mesh(Face_element_mesh_pt);
      }
      build_global_mesh();
      assign_eqn_numbers();
    }

    /// Impose a load to the surface of the plate.
    /// note : MUST be called before neumann boundary conditions are imposed,
    /// i.e. a free edge or a simply supported edge with laplacian(u) imposed
    void set_source_function(const BiharmonicSourceFctPt source_fct_pt)
    {
      // number of elements in mesh
      unsigned n_bulk_element = Bulk_element_mesh_pt->nelement();

      // loop over bulk elements
      for (unsigned i = 0; i < n_bulk_element; i++)
      {
        // upcast from generalised element to specific element
        BiharmonicElement<2>* element_pt = dynamic_cast<BiharmonicElement<2>*>(
          Bulk_element_mesh_pt->element_pt(i));

        // set the source function pointer
        element_pt->source_fct_pt() = source_fct_pt;
      }
    }

    /// Imposes the prescribed dirichlet BCs u (u_fn) and
    /// du/dn (dudn_fn) dirichlet BCs by 'pinning'
    void set_dirichlet_boundary_condition(const unsigned& b,
                                          DirichletBCFctPt u_fn = 0,
                                          DirichletBCFctPt dudn_fn = 0);

    /// Imposes the prescribed Neumann BCs laplacian(u)  (flux0_fct_pt)
    /// and dlaplacian(u)/dn (flux1_fct_pt) with flux edge elements
    void set_neumann_boundary_condition(
      const unsigned& b,
      BiharmonicFluxElement<2>::FluxFctPt flux0_fct_pt,
      BiharmonicFluxElement<2>::FluxFctPt flux1_fct_pt = 0);

  public:
    // NOTE: these two private meshes are required for the block
    // preconditioners.

    /// Mesh for BiharmonicElement<DIM> only - the block preconditioner
    /// assemble the global equation number to block number mapping from
    /// elements in this mesh only
    Mesh* Bulk_element_mesh_pt;

    /// mesh for face elements
    Mesh* Face_element_mesh_pt;
  };


  //=============================================================================
  /// Biharmonic Fluid Problem Class - describes stokes flow in 2D.
  /// Developed for the topologically rectangular Hermite Element Mesh. Contains
  /// functions allowing the following boundary conditions to be applied (on a
  /// given edge):
  ///   + wall :           v_n = 0 and v_t = 0 (psi must also be prescribed)
  ///   + traction free :  v_t = 0
  ///   + flow :           v_n and v_t are prescribed
  /// NOTE 1 : psi is the stream function
  ///            + fluid velocity normal to boundary v_n = +/- dpsi/dt
  ///            + fluid velocity tangential to boundary v_t = -/+ dpsi/dn
  /// NOTE 2 : when a solid wall boundary condition is applied to ensure that
  ///          v_n = 0 the the streamfunction psi must also be prescribed (and
  ///          constant)
  //=============================================================================
  template<unsigned DIM>
  class BiharmonicFluidProblem : public Problem
  {
  public:
    /// Definition of a dirichlet boundary condition function pointer.
    /// Takes the position  along a boundary (s) in the macro element coordinate
    /// scheme and returns the fluid velocity normal (dpsi/dt) to the boundary
    /// (u[0]) and the fluid velocity tangential (dpsidn) to the boundary
    /// (u[1]).
    typedef void (*FluidBCFctPt)(const double& s, Vector<double>& u);


    /// constructor
    BiharmonicFluidProblem()
    {
      // initialise the number of non bulk elements
      Npoint_element = 0;
    }


    /// actions before solve, performs self test
    void actions_before_newton_solve()
    {
#ifdef PARANOID
      if (0 == self_test())
      {
        oomph_info << "self test passed" << std::endl;
      }
      else
      {
        oomph_info << "self test failed" << std::endl;
      }
#endif
    }


    /// action after solve
    void actions_after_newton_solve() {}


    /// documents the solution, and if an exact solution is provided,
    /// then the error between the numerical and exact solution is presented
    void doc_solution(
      DocInfo& doc_info,
      FiniteElement::SteadyExactSolutionFctPt exact_soln_pt = 0);


  protected:
    /// Imposes a solid boundary on boundary b - no flow into boundary
    /// or along boundary v_n = 0 and v_t = 0. User must presribe the
    /// streamfunction psi to ensure dpsi/dt = 0 is imposed at all points on the
    /// boundary and not just at the nodes
    void impose_solid_boundary_on_edge(const unsigned& b,
                                       const double& psi = 0);

    /// Impose a traction free edge - i.e. v_t = 0 or dpsi/dn = 0. In
    /// general dpsi/dn = 0 can only be imposed using equation elements to set
    /// the DOFs dpsi/ds_n, however in the special case of  dt/ds_n = 0, then
    /// dpsi/ds_n = 0 and can be imposed using pinning - this is handled
    /// automatically in this function. For a more detailed description of the
    /// equations see the description of the class
    /// BiharmonicFluidBoundaryElement
    void impose_traction_free_edge(const unsigned& b);


    /// Impose a prescribed fluid flow comprising the velocity normal to
    /// the boundary (u_imposed_fn[0]) and the velocity tangential to the
    /// boundary (u_imposed_fn[1])
    void impose_fluid_flow_on_edge(const unsigned& b,
                                   FluidBCFctPt u_imposed_fn);


  private:
    // number of non-bulk elements - i.e. biharmonic fluid boundary elements
    unsigned Npoint_element;
  };


  //=============================================================================
  /// Point equation element used to impose the traction free edge (i.e.
  /// du/dn = 0) on the boundary when dt/ds_n != 0. The following equation is
  /// implemented :  du/ds_n = dt/ds_n * ds_t/dt * du/dt.
  /// The bulk biharmonic elements on the boundary must be hijackable and the
  /// du/ds_n and d2u/ds_nds_t boundary DOFs hijacked when these  elements are
  /// applied. At any node where dt/ds_n = 0 we can impose  du/ds_n = 0 and
  /// d2u/ds_nds_t = 0 using pinning - see
  /// BiharmonicFluidProblem::impose_traction_free_edge()
  //=============================================================================
  class BiharmonicFluidBoundaryElement : public virtual PointElement
  {
  public:
    // constructor
    BiharmonicFluidBoundaryElement(Node* node_pt, const unsigned s_fixed_index)
    {
      // set the node pt
      this->node_pt(0) = node_pt;

      // store fixed index on the boundary
      S_fixed_index = s_fixed_index;
    }

    /// Output function -- does nothing
    void output(std::ostream& outfile) {}


    /// Output function -- does nothing
    void output(std::ostream& outfile, const unsigned& n_plot) {}


    /// Output function -- does nothing
    void output_fluid_velocity(std::ostream& outfile, const unsigned& n_plot) {}


    /// C-style output function -- does nothing
    void output(FILE* file_pt) {}


    /// C-style output function -- does nothing
    void output(FILE* file_pt, const unsigned& n_plot) {}


    /// compute_error -- does nothing
    void compute_error(std::ostream& outfile,
                       FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
                       double& error,
                       double& norm)
    {
    }


    /// Compute the elemental residual vector - wrapper function called
    /// by get_residuals in GeneralisedElement
    inline void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      // create a dummy matrix
      DenseDoubleMatrix dummy(1);

      // call the generic residuals functions with flag set to zero
      fill_in_generic_residual_contribution_biharmonic_boundary(
        residuals, dummy, 0);
    }


    /// Compute the elemental residual vector and jacobian matrix -
    /// wrapper function called by get_jacobian in GeneralisedElement
    inline void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                                 DenseMatrix<double>& jacobian)
    {
      // call generic routine with flag set to 1
      fill_in_generic_residual_contribution_biharmonic_boundary(
        residuals, jacobian, 1);
    }


    /// Computes the elemental residual vector and the elemental jacobian
    /// matrix if JFLAG = 0
    /// Imposes the equations :  du/ds_n = dt/ds_n * ds_t/dt * du/dt
    virtual void fill_in_generic_residual_contribution_biharmonic_boundary(
      Vector<double>& residuals, DenseMatrix<double>& jacobian, unsigned JFLAG);

  private:
    // fixed local coordinate index on boundary
    unsigned S_fixed_index;
  };


} // namespace oomph
#endif
