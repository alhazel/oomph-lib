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
// Header file for KL beam elements
#ifndef OOMPH_KIRCHHOFF_LOVE_BEAM_ELEMENTS_HEADER
#define OOMPH_KIRCHHOFF_LOVE_BEAM_ELEMENTS_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

// OOMPH-LIB header files
#include "../generic/hermite_elements.h"
#include "../generic/geom_objects.h"
#include "../generic/fsi.h"
#include "../generic/block_preconditioner.h"

namespace oomph
{
  //=======================================================================
  /// A class for elements that solve the equations of Kirchhoff-Love
  /// large-displacement (but linearly-elastic) thin-beam theory.
  ///
  /// The variational principle has the form
  /// \f[ \int_0^{L} \left[ (\sigma_0 + \gamma) \ \delta \gamma + \frac{1}{12} \left(\frac{h}{R_0}\right)^2 \kappa \ \delta \kappa - \left( \left(\frac{R_0}{h}\right) {\bf f} - \Lambda^2 \frac{\partial^2 {\bf R}_w}{\partial t^2} \right) \cdot \delta {\bf R}_w \right] \ d\xi = 0, \f]
  /// where all lengths have been non-dimensionalised w.r.t. \f$ R_0 \f$.
  /// The strain and and bending "tensors" \f$\gamma\f$ and \f$\kappa\f$
  /// are computed relative to the shape of the beam's undeformed shape
  /// which is specified as a GeomObject.
  ///
  /// Time is scaled on the timescale \f$T\f$ and
  /// \f[ \Lambda = \frac{a}{T} \sqrt{\frac{\rho}{E_{eff}}}, \f]
  /// the ratio of the timescale used in the non-dimensionalisation of the
  /// equations to the natural timescale of the wall oscillations (in the
  /// wall's in-plane mode). \f$ \Lambda^2 \f$ can be interpreted as
  /// the non-dimensional wall density, therefore \f$ \Lambda=0\f$
  /// corresponds to the case without wall inertia.
  ///
  ///
  /// Note that:
  /// - the load vector \f$ {\bf f} \f$ is scaled
  ///   on the effective elastic modulus \f$ E_{eff}=E/(1-\nu^2)\f$
  ///   (rather than the
  ///   bending stiffness). Rescale the result yourself if you prefer
  ///   another non-dimensionalisation (the current version yields the
  ///   the most compact maths).
  /// - Poisson's ratio does not appear explicitly since it only occurs
  ///   in combination with Young's modulus \f$E\f$.
  ///
  /// Default values:
  /// - the 2nd Piola Kirchhoff pre-stress \f$ \sigma_0 \f$ is zero.
  /// - the wall thickness \f$ h/R_0\f$ is 1/20.
  /// - the timescale ratio \f$ \Lambda^2\f$ is 1.
  /// - the traction vector \f$ f \f$ evaluates to zero.
  ///
  /// Need to specify:
  /// - the undeformed wall shape (as a GeomObject).
  ///
  /// The governing equations can be switched from the principle of
  /// virtual displacements to a system of equations that forces the
  /// beam to deform into a shape specified by a SolidInitialCondition object.
  /// If \c SolidFiniteElement::solid_ic_pt()!=0 we solve the
  /// the equations
  /// \f[ \int_0^{L} \left( \frac{\partial^i {\bf R}_{IC}}{\partial t^i} - {\bf R}_w \right) \psi_{jk} \ d\xi = 0, \f]
  /// where  \f$ \partial^i {\bf R}_{IC}/\partial t^i\f$ is
  /// implemented by the SolidInitialCondition object, pointed to by
  /// \c SolidFiniteElement::shell_ic_pt().
  ///
  //=======================================================================
  class KirchhoffLoveBeamEquations : public virtual SolidFiniteElement
  {
  private:
    /// Static default value for 2nd Piola Kirchhoff prestress
    static double Default_sigma0_value;

    /// Static default value for timescale ratio (1.0 -- for natural scaling)
    static double Default_lambda_sq_value;

    /// Static default value for non-dim wall thickness
    // i.e. The reference value 'h_0'
    static double Default_h_value;

    /// Pointer to axial prestress
    double* Sigma0_pt;

    /// Pointer to wall thickness
    // i.e. The reference value 'h_0'
    double* H_pt;

    /// Pointer to Timescale ratio (non-dim. density)
    double* Lambda_sq_pt;

  protected:
    /// Default load function (zero traction)
    static void Zero_traction_fct(const Vector<double>& xi,
                                  const Vector<double>& x,
                                  const Vector<double>& N,
                                  Vector<double>& load);

    /// Pointer to load vector function: Its arguments are:
    /// Lagrangian coordinate, Eulerian coordinate, normal vector and
    /// load vector itself (not all of the input arguments will be
    /// required for all specific load functions but the list should
    /// cover all cases)
    void (*Load_vector_fct_pt)(const Vector<double>& xi,
                               const Vector<double>& x,
                               const Vector<double>& N,
                               Vector<double>& load);

    /// Default profile function (constant thickness 'h_0')
    static void Unit_profile_fct(const Vector<double>& xi,
                                 const Vector<double>& x,
                                 double& h_ratio);

    /// Pointer to wall profile function: Its arguments are:
    /// Lagrangian coordinate, Eulerian coordinate, and
    /// profile itself (not all of the input arguments will be
    /// required for all specific profile functions but the list should
    /// cover all cases)
    void (*Wall_profile_fct_pt)(const Vector<double>& xi,
                                const Vector<double>& x,
                                double& h_ratio);

    /// Pointer to the GeomObject that specifies the beam's
    /// undeformed midplane
    GeomObject* Undeformed_beam_pt;

  public:
    /// Constructor. Set default values for all physical parameters
    /// and zero traction.
    KirchhoffLoveBeamEquations() : Undeformed_beam_pt(0)
    {
      // Set physical parameter pointers to the default values
      Sigma0_pt = &Default_sigma0_value;
      Lambda_sq_pt = &Default_lambda_sq_value;
      // The reference thickness 'h_0'
      H_pt = &Default_h_value;
      // Zero traction
      Load_vector_fct_pt = &Zero_traction_fct;
      // Unit thickness profile
      Wall_profile_fct_pt = &Unit_profile_fct;
    }


    /// Reference to the load vector function pointer
    void (*&load_vector_fct_pt())(const Vector<double>& xi,
                                  const Vector<double>& x,
                                  const Vector<double>& N,
                                  Vector<double>& load)
    {
      return Load_vector_fct_pt;
    }


    /// Get the load vector: Pass number of integration point (dummy),
    /// Lagr. and Eulerian coordinate and normal vector and return the load
    /// vector (not all of the input arguments will be required for all specific
    /// load functions but the list should cover all cases). This function is
    /// virtual so it can be overloaded for FSI.
    virtual void load_vector(const unsigned& intpt,
                             const Vector<double>& xi,
                             const Vector<double>& x,
                             const Vector<double>& N,
                             Vector<double>& load)
    {
      Load_vector_fct_pt(xi, x, N, load);
    }

    /// Reference to the wall thickness ratio profile function pointer
    void (*&wall_profile_fct_pt())(const Vector<double>& xi,
                                   const Vector<double>& x,
                                   double& h_ratio)
    {
      return Wall_profile_fct_pt;
    }


    /// Get the wall profile: Pass Lagrangian & Eulerian coordinate
    /// and return the wall profile (not all of the input arguments will be
    /// required for all specific thickness functions but the list should cover
    /// all cases).
    void wall_profile(const Vector<double>& xi,
                      const Vector<double>& x,
                      double& h_ratio)
    {
      Wall_profile_fct_pt(xi, x, h_ratio);
    }


    /// Return the non-dimensional wall thickness
    // i.e. the reference value 'h_0'
    const double& h() const
    {
      return *H_pt;
    }

    /// Return the timescale ratio (non-dimensional density)
    const double& lambda_sq() const
    {
      return *Lambda_sq_pt;
    }

    /// Return the axial prestress
    const double& sigma0() const
    {
      return *Sigma0_pt;
    }

    /// Return a pointer to axial prestress
    double*& sigma0_pt()
    {
      return Sigma0_pt;
    }

    /// Return a pointer to non-dim. wall thickness
    //  i.e. the reference value 'h_0'
    double*& h_pt()
    {
      return H_pt;
    }

    /// Return a pointer to timescale ratio (nondim density)
    double*& lambda_sq_pt()
    {
      return Lambda_sq_pt;
    }

    /// Return a Pointer to geometric object that specifies the beam's
    /// undeformed geometry
    GeomObject*& undeformed_beam_pt()
    {
      return Undeformed_beam_pt;
    }

    /// Get normal vector on wall
    void get_normal(const Vector<double>& s, Vector<double>& N)
    {
      Vector<double> r(2);
      get_normal(s, r, N);
    }


    /// Get position vector to and normal vector on wall
    void get_normal(const Vector<double>& s,
                    Vector<double>& r,
                    Vector<double>& N);

    /// Get position vector to and non-unit tangent vector on wall:
    /// dr/ds
    void get_non_unit_tangent(const Vector<double>& s,
                              Vector<double>& r,
                              Vector<double>& drds);

    /// Return the residuals for the equations of Kirchhoff-Love beam
    /// theory with linear constitutive equations; if  Solid_ic_pt!=0, we
    /// assign residuals which force the assignement of an initial shape/
    /// veloc/accel to the dofs. This overloads the standard interface.
    void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      fill_in_contribution_to_residuals_beam(residuals);
    }


    /// Return the residuals for the equations of Kirchhoff-Love beam
    /// theory with linear constitutive equations; if  Solid_ic_pt!=0, we
    /// assign residuals which force the assignement of an initial shape/
    /// veloc/accel to the dofs.
    void fill_in_contribution_to_residuals_beam(Vector<double>& residuals);


    /// Get FE jacobian and residuals (Jacobian done by finite differences)
    virtual void fill_in_contribution_to_jacobian(
      Vector<double>& residuals, DenseMatrix<double>& jacobian);

    /// Get potential (strain) and kinetic energy of the element
    void get_energy(double& pot_en, double& kin_en);

    /// Get the potential energy due to stretching and bending and the
    /// kinetic energy of the element
    void get_energy(double& stretch, double& bend, double& kin_en);
  };


  //=========================================================================
  /// Hermite Kirchhoff Love beam. Implements KirchhoffLoveBeamEquations
  /// using 2-node Hermite elements as the underlying geometrical elements.
  //=========================================================================
  class HermiteBeamElement : public virtual SolidQHermiteElement<1>,
                             public KirchhoffLoveBeamEquations
  {
  public:
    /// Constructor (empty)
    HermiteBeamElement()
      : SolidQHermiteElement<1>(), KirchhoffLoveBeamEquations()
    {
      // Set the number of dimensions at each node (2D node on 1D surface)
      set_nodal_dimension(2);
    }

    /// Output function
    void output(std::ostream& outfile);

    /// Output function with specified number of plot points
    void output(std::ostream& outfile, const unsigned& n_plot);

    /// Output at previous time (t=0: present; t>0: previous)
    /// with specified number of plot points
    void output(const unsigned& t,
                std::ostream& outfile,
                const unsigned& n_plot) const;

    /// C-style output function
    void output(FILE* file_pt);

    /// C-style output function with specified number of plot points
    void output(FILE* file_pt, const unsigned& n_plot);

    /// C-style output at previous time (t=0: present; t>0: previous)
    /// with specified number of plot points
    void output(const unsigned& t, FILE* file_pt, const unsigned& n_plot) const;
  };

  //=========================================================================
  /// Hermite Kirchhoff Love beam "upgraded" to a FSIWallElement (and thus,
  /// by inheritance, a GeomObject), so it can be used in FSI.
  //=========================================================================
  class FSIHermiteBeamElement : public virtual HermiteBeamElement,
                                public virtual FSIWallElement
  {
  private:
    // Boolean flag to indicate whether the normal is directed into the fluid
    bool Normal_points_into_fluid;

  public:
    /// Constructor: Create beam element as FSIWallElement (and thus,
    /// by inheritance, a GeomObject). By default, we assume that the
    /// normal vector computed by KirchhoffLoveBeamEquations::get_normal(...)
    /// points into the fluid. If this is not the case, overwrite this
    /// with the access function
    /// FSIHermiteBeamElement::set_normal_pointing_out_of_fluid()
    FSIHermiteBeamElement()
      : HermiteBeamElement(), Normal_points_into_fluid(true)
    {
      unsigned n_lagr = 1;
      unsigned n_dim = 2;
      setup_fsi_wall_element(n_lagr, n_dim);
    }

    /// Destructor: empty
    ~FSIHermiteBeamElement() {}

    /// Set the normal computed by
    /// KirchhoffLoveBeamEquations::get_normal(...) to point into the fluid
    void set_normal_pointing_into_fluid()
    {
      Normal_points_into_fluid = true;
    }

    /// Set the normal computed by
    /// KirchhoffLoveBeamEquations::get_normal(...) to point out of the fluid
    void set_normal_pointing_out_of_fluid()
    {
      Normal_points_into_fluid = false;
    }


    /// Derivative of position vector w.r.t. the SolidFiniteElement's
    /// Lagrangian coordinates; evaluated at current time.
    void dposition_dlagrangian_at_local_coordinate(
      const Vector<double>& s, DenseMatrix<double>& drdxi) const;

    /// Get the load vector: Pass number of the integration point,
    /// Lagr. coordinate, Eulerian coordinate and normal vector
    /// and return the load vector. (Not all of the input arguments will be
    /// required for all specific load functions but the list should
    /// cover all cases). We first evaluate the load function defined via
    /// KirchhoffLoveBeamEquations::load_vector_fct_pt() -- this
    /// represents the non-FSI load on the beam, e.g. an external
    /// pressure load. Then we add to this the FSI load due to
    /// the traction exerted by the adjacent FSIFluidElements, taking
    /// the sign of the normal into account.
    void load_vector(const unsigned& intpt,
                     const Vector<double>& xi,
                     const Vector<double>& x,
                     const Vector<double>& N,
                     Vector<double>& load)
    {
      // Initially call the standard Load_vector_fct_pt
      Load_vector_fct_pt(xi, x, N, load);

      // Memory for the FSI load
      Vector<double> fsi_load(2);

      // Get the fluid load on the wall stress scale
      fluid_load_vector(intpt, N, fsi_load);

      // If the normal is outer to the fluid switch the direction
      double sign = 1.0;
      if (!Normal_points_into_fluid)
      {
        sign = -1.0;
      }

      // Add the FSI load to the load vector
      for (unsigned i = 0; i < 2; i++)
      {
        load[i] += sign * fsi_load[i];
      }
    }

    /// Get the Jacobian and residuals. Wrapper to generic FSI version;
    /// that catches the case when we replace the Jacobian by the
    /// mass matrix (for the consistent assignment of initial conditions).
    virtual void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                                  DenseMatrix<double>& jacobian)
    {
      // Call the standard beam element's jacobian function
      HermiteBeamElement::fill_in_contribution_to_jacobian(residuals, jacobian);
      // Now add the external interaction data by finite differences
      this->fill_in_jacobian_from_external_interaction_by_fd(jacobian);
    }

    /// Find the local coordinate s in this element
    /// that corresponds to the global "intrinsic" coordinate \f$ \zeta \f$
    /// (here identical to the Lagrangian coordinate \f$ \xi \f$).
    /// If the coordinate is contained within this element, the
    /// geom_object_pt points to "this" element; if the zeta coordinate
    /// is not contained in this element geom_object_pt=NULL.
    /// By default don't use any value passed in to the local coordinate s
    /// as the initial guess in the Newton method
    void locate_zeta(const Vector<double>& zeta,
                     GeomObject*& geom_object_pt,
                     Vector<double>& s,
                     const bool& use_coordinate_as_initial_guess = false);


    /// The number of "DOF types" that degrees of freedom in this element
    /// are sub-divided into: Just the solid degrees of freedom themselves.
    unsigned ndof_types() const
    {
      return 1;
    }

    /// Create a list of pairs for all unknowns in this element,
    /// so that the first entry in each pair contains the global equation
    /// number of the unknown, while the second one contains the number
    /// of the "DOF type" that this unknown is associated with.
    /// (Function can obviously only be called if the equation numbering
    /// scheme has been set up.)
    void get_dof_numbers_for_unknowns(
      std::list<std::pair<unsigned long, unsigned>>& dof_lookup_list) const;
  };


  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////


  //=======================================================================
  /// Face geometry for the HermiteBeam elements: Solid point element
  //=======================================================================
  template<>
  class FaceGeometry<HermiteBeamElement> : public virtual SolidPointElement
  {
  public:
    /// Constructor [this was only required explicitly
    /// from gcc 4.5.2 onwards...]
    FaceGeometry() {}
  };


  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////


  //======================================================================
  /// Element that allows the imposition of boundary
  /// conditions for a beam that is clamped but can slide
  /// along a line which is specified by a position vector
  /// to that line and the normal vector to it. The endpoint
  /// of the beam is forced to stay on that line and meet
  /// it at a right angle. This is achieved with Lagrange multipliers.
  //======================================================================
  class ClampedSlidingHermiteBeamBoundaryConditionElement
    : public virtual FaceGeometry<HermiteBeamElement>,
      public virtual SolidFaceElement
  {
  public:
    /// Constructor, takes the pointer to the "bulk" element, the
    /// index of the fixed local coordinate and its value represented
    /// by an integer (+/- 1), indicating that the face is located
    /// at the max. or min. value of the "fixed" local coordinate
    /// in the bulk element.
    ClampedSlidingHermiteBeamBoundaryConditionElement(
      FiniteElement* const& bulk_el_pt, const int& face_index);

    /// Broken empty constructor
    ClampedSlidingHermiteBeamBoundaryConditionElement()
    {
      throw OomphLibError("Don't call empty constructor for "
                          "ClampedSlidingHermiteBeamBoundaryConditionElement ",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }


    /// Broken copy constructor
    ClampedSlidingHermiteBeamBoundaryConditionElement(
      const ClampedSlidingHermiteBeamBoundaryConditionElement& dummy) = delete;

    /// Broken assignment operator
    // Commented out broken assignment operator because this can lead to a
    // conflict warning when used in the virtual inheritence hierarchy.
    // Essentially the compiler doesn't realise that two separate
    // implementations of the broken function are the same and so, quite
    // rightly, it shouts.
    /*void operator=(const ClampedSlidingHermiteBeamBoundaryConditionElement&) =
      delete;*/


    /// Set vectors to some point on the symmetry line, and
    /// normal to that line along which the end of the beam is sliding.
    void set_symmetry_line(const Vector<double>& vector_to_symmetry_line,
                           const Vector<double>& normal_to_symmetry_line)
    {
      Vector_to_symmetry_line[0] = vector_to_symmetry_line[0];
      Vector_to_symmetry_line[1] = vector_to_symmetry_line[1];
      Normal_to_symmetry_line[0] = normal_to_symmetry_line[0];
      Normal_to_symmetry_line[1] = normal_to_symmetry_line[1];
    }


    /// Fill in the element's contribution to its residual vector
    void fill_in_contribution_to_residuals(Vector<double>& residuals);


    /// Output function -- forward to broken version in FiniteElement
    /// until somebody decides what exactly they want to plot here...
    void output(std::ostream& outfile)
    {
      FiniteElement::output(outfile);
    }

    /// Output function -- forward to broken version in FiniteElement
    /// until somebody decides what exactly they want to plot here...
    void output(std::ostream& outfile, const unsigned& n_plot)
    {
      FiniteElement::output(outfile, n_plot);
    }

    /// C-style output function -- forward to broken version in FiniteElement
    /// until somebody decides what exactly they want to plot here...
    void output(FILE* file_pt)
    {
      FiniteElement::output(file_pt);
    }

    /// C-style output function -- forward to broken version in
    /// FiniteElement until somebody decides what exactly they want to plot
    /// here...
    void output(FILE* file_pt, const unsigned& n_plot)
    {
      FiniteElement::output(file_pt, n_plot);
    }

  private:
    /// Vector to some point on the symmetry line along which the
    /// end of the beam is sliding
    Vector<double> Vector_to_symmetry_line;

    /// Normal vector to the symmetry line along which the
    /// end of the beam is sliding
    Vector<double> Normal_to_symmetry_line;
  };


} // namespace oomph

#endif
