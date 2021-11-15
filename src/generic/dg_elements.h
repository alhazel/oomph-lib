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
// Header functions for classes that define Discontinuous Galerkin elements

// Include guards to prevent multiple inclusion of the header
#ifndef OOMPH_DG_ELEMENT_HEADER
#define OOMPH_DG_ELEMENT_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

// Oomph-lib header
#include "elements.h"
#include "mesh.h"

namespace oomph
{
  //=============================================================
  /// Base class for Discontinuous Galerkin Faces.
  /// These are responsible for calculating the normal fluxes
  /// that provide the communication between the discontinuous
  /// elements.
  //===============================================================
  class DGFaceElement : public virtual FaceElement
  {
    /// Vector of neighbouring face elements at the integration points
    Vector<FaceElement*> Neighbour_face_pt;

    /// Vector of neighbouring local coordinates at the integration points
    Vector<Vector<double>> Neighbour_local_coordinate;

    /// Vector of the vectors that will store the number of the
    /// bulk external data that correspond to the dofs in the neighbouring face
    /// This is only used if we are using implict timestepping and wish
    /// to assemble the jacobian matrix. In order that the data be correctly
    /// set up DGMesh::setup_face_neighbour_info() must be called with the
    /// boolean flag set to true.
    Vector<Vector<unsigned>> Neighbour_external_data;

  protected:
    /// Return the index at which the i-th unknown flux is stored.
    // The default return is suitable for single-physics problem
    virtual inline unsigned flux_index(const unsigned& i) const
    {
      return i;
    }

    /// Set the number of flux components
    virtual unsigned required_nflux()
    {
      return 0;
    }

  public:
    /// Empty Constructor
    DGFaceElement() : FaceElement() {}

    /// Empty Destructor
    virtual ~DGFaceElement() {}

    /// Access function for neighbouring face information
    FaceElement* neighbour_face_pt(const unsigned& i)
    {
      return Neighbour_face_pt[i];
    }

    /// Setup information from the neighbouring face, return a set
    /// of nodes (as data) in the neighour. The boolean flag
    /// is used to determine whether the data in the neighbour are
    /// added as external data to the bulk element --- required when
    /// computing the jacobian of the system
    void setup_neighbour_info(const bool& add_neighbour_data_to_bulk);

    /// Output information about the present element and its neighbour
    void report_info();

    // Get the value of the unknowns
    virtual void interpolated_u(const Vector<double>& s, Vector<double>& f);

    /// Get the data that are used to interpolate the unkowns
    /// in the element. These must be returned in order.
    virtual void get_interpolation_data(Vector<Data*>& interpolation_data);


    /// Calculate the normal numerical flux at the integration point.
    /// This is the most general interface that can be overloaded if desired
    virtual void numerical_flux_at_knot(const unsigned& ipt,
                                        const Shape& psi,
                                        Vector<double>& flux,
                                        DenseMatrix<double>& dflux_du_int,
                                        DenseMatrix<double>& dflux_du_ext,
                                        unsigned flag);

    /// Calculate the normal flux, which is the dot product of our
    /// approximation to the flux with the outer unit normal
    virtual void numerical_flux(const Vector<double>& n_out,
                                const Vector<double>& u_int,
                                const Vector<double>& u_ext,
                                Vector<double>& flux)
    {
      std::ostringstream error_stream;
      error_stream
        << "Empty numerical flux function called\n"
        << "This function should be overloaded with a specific flux\n"
        << "that is appropriate to the equations being solved.\n";
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    /// Calculate the derivative of the
    /// normal flux, which is the dot product of our
    /// approximation to the flux with the outer unit normal,
    /// with respect to the interior and exterior variables
    /// Default is to use finite differences
    virtual void dnumerical_flux_du(const Vector<double>& n_out,
                                    const Vector<double>& u_int,
                                    const Vector<double>& u_ext,
                                    DenseMatrix<double>& dflux_du_int,
                                    DenseMatrix<double>& dflux_du_ext);


    /// Add the contribution from integrating the numerical flux
    // over the face to the residuals
    void add_flux_contributions(Vector<double>& residuals,
                                DenseMatrix<double>& jacobian,
                                unsigned flag);
  };

  class DGMesh;
  class SlopeLimiter;

  //==================================================================
  /// A Base class for DGElements
  //=================================================================
  class DGElement : public virtual FiniteElement
  {
  protected:
    /// The DGFaceElement requires access to the nodal local equation
    /// information, so it's a friend
    friend class DGFaceElement;

    /// Vector of pointers to faces of the element
    Vector<FaceElement*> Face_element_pt;

    /// Pointer to Mesh, which will be responsible for the neighbour finding
    DGMesh* DG_mesh_pt;

    /// Pointer to storage for a mass matrix that can be recycled if
    /// desired
    DenseDoubleMatrix* M_pt;

    /// Pointer to storage for the average values of the of the
    /// variables over the element
    double* Average_value;

    /// Boolean flag to indicate whether to reuse the mass matrix
    bool Mass_matrix_reuse_is_enabled;

    /// Boolean flag to indicate whether the mass matrix has been
    /// computed
    bool Mass_matrix_has_been_computed;

    /// Boolean flag to indicate whether the mass matrix can be
    /// deleted (i.e. was it created by this element)
    bool Can_delete_mass_matrix;

    /// Set the number of flux components
    virtual unsigned required_nflux()
    {
      return 0;
    }

  public:
    /// Constructor, initialise the pointers to zero
    DGElement()
      : DG_mesh_pt(0),
        M_pt(0),
        Average_value(0),
        Mass_matrix_reuse_is_enabled(false),
        Mass_matrix_has_been_computed(false),
        Can_delete_mass_matrix(true)
    {
    }

    /// Virtual destructor, destroy the mass matrix, if we created it
    /// Clean-up storage associated with average values
    virtual ~DGElement()
    {
      if ((M_pt != 0) && Can_delete_mass_matrix)
      {
        delete M_pt;
      }
      if (this->Average_value != 0)
      {
        delete[] Average_value;
        Average_value = 0;
      }
    }

    /// Access function for the boolean to indicate whether the
    /// mass matrix has been computed
    bool mass_matrix_has_been_computed()
    {
      return Mass_matrix_has_been_computed;
    }

    /// Function that allows the reuse of the mass matrix
    void enable_mass_matrix_reuse()
    {
      Mass_matrix_reuse_is_enabled = true;
      Mass_matrix_has_been_computed = false;
    }

    /// Function that disables the reuse of the mass matrix
    void disable_mass_matrix_reuse()
    {
      // If we are using another element's mass matrix then reset our pointer
      // to zero
      if (!Can_delete_mass_matrix)
      {
        M_pt = 0;
      }
      // Otherwise we do not reuse the mass matrix
      Mass_matrix_reuse_is_enabled = false;
      // Recalculate the mass matrix
      Mass_matrix_has_been_computed = false;
    }


    /// Set the mass matrix to point to one in another element
    virtual void set_mass_matrix_from_element(DGElement* const& element_pt)
    {
      // If the element's mass matrix has not been computed, compute it!
      if (!element_pt->mass_matrix_has_been_computed())
      {
        element_pt->pre_compute_mass_matrix();
      }

      // Now set the mass matrix in this element to address that
      // of element_pt
      this->M_pt = element_pt->M_pt;
      // We must reuse the mass matrix, or there will be trouble
      // Because we will recalculate it in the original element
      Mass_matrix_reuse_is_enabled = true;
      Mass_matrix_has_been_computed = true;
      // We cannot delete the mass matrix
      Can_delete_mass_matrix = false;
    }

    /// Function that computes and stores the (inverse) mass matrix
    void pre_compute_mass_matrix();

    // Function that is used to construct all the faces of the DGElement
    virtual void build_all_faces() = 0;

    /// Function that returns the current value of the residuals
    /// multiplied by the inverse mass matrix (virtual so that it can be
    /// overloaded specific elements in which time saving tricks can be applied)
    virtual void get_inverse_mass_matrix_times_residuals(
      Vector<double>& minv_res);

    /// Construct all nodes and faces of the element.
    /// The vector of booleans boundary should be the same size
    /// as the number of nodes and if any entries are true
    /// that node will be constructed as a boundary node.
    void construct_boundary_nodes_and_faces(DGMesh* const& mesh_pt,
                                            std::vector<bool>& boundary_flag,
                                            TimeStepper* const& time_stepper_pt)
    {
      // Construct the nodes (This should not be used in a base class)
      const unsigned n_node = this->nnode();
#ifdef PARANOID
      if (boundary_flag.size() != n_node)
      {
        std::ostringstream error_stream;
        error_stream
          << "Size of boundary_flag vector is " << boundary_flag.size() << "\n"
          << "It must be the same as the number of nodes in the element "
          << n_node << "\n";

        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // Loop over the nodes
      for (unsigned n = 0; n < n_node; n++)
      {
        // If the node is marked on a boundary construct a boundary node
        if (boundary_flag[n])
        {
          (void)this->construct_boundary_node(n, time_stepper_pt);
        }
        // Otherwise, construct a normal node
        else
        {
          (void)this->construct_node(n, time_stepper_pt);
        }
      }

      // Make the faces
      this->build_all_faces();

      // Set the Mesh pointer
      DG_mesh_pt = mesh_pt;
    }


    /// Construct the nodes and faces of the element
    void construct_nodes_and_faces(DGMesh* const& mesh_pt,
                                   TimeStepper* const& time_stepper_pt)
    {
      // Loop over the nodes
      const unsigned n_node = this->nnode();
      for (unsigned n = 0; n < n_node; n++)
      {
        // Construct the node and ignore the return value
        (void)this->construct_node(n, time_stepper_pt);
      }

      // Make the faces
      this->build_all_faces();

      // Set the Mesh pointer
      DG_mesh_pt = mesh_pt;
    }

    // Set the mesh pointer of the element
    void set_mesh_pt(DGMesh*& mesh_pt)
    {
      DG_mesh_pt = mesh_pt;
    }

    /// Return the number of faces
    unsigned nface() const
    {
      return Face_element_pt.size();
    }

    /// Access function for the faces
    DGFaceElement* face_element_pt(const unsigned& i)
    {
      return dynamic_cast<DGFaceElement*>(Face_element_pt[i]);
    }

    /// Output the faces of the element
    void output_faces(std::ostream& outfile)
    {
      // Loop over the faces
      unsigned n_face = nface();
      for (unsigned f = 0; f < n_face; f++)
      {
        face_element_pt(f)->output(outfile);
      }
    }

    /// Return the neighbour info
    void get_neighbouring_face_and_local_coordinate(
      const int& face_index,
      const Vector<double>& s,
      FaceElement*& face_element_pt,
      Vector<double>& s_face);

    // Setup the face information
    /// The boolean flag determines whether the data from the neighbouring
    /// elements is added as external data to the element (required for
    /// correct computation of the jacobian)
    void setup_face_neighbour_info(const bool& add_face_data_as_external)
    {
      unsigned n_face = this->nface();
      for (unsigned f = 0; f < n_face; f++)
      {
        face_element_pt(f)->setup_neighbour_info(add_face_data_as_external);
        face_element_pt(f)->report_info();
      }
    }


    /// Loop over all faces and add their integrated numerical fluxes
    /// to the residuals
    void add_flux_contributions_to_residuals(Vector<double>& residuals,
                                             DenseMatrix<double>& jacobian,
                                             unsigned flag)
    {
      // Add up the contributions from each face
      unsigned n_face = this->nface();
      for (unsigned f = 0; f < n_face; f++)
      {
        face_element_pt(f)->add_flux_contributions(residuals, jacobian, flag);
      }
    }

    /// Limit the slope within the element
    void slope_limit(SlopeLimiter* const& slope_limiter_pt);

    /// Calculate the averages in the element
    virtual void calculate_element_averages(double*& average_values)
    {
      throw OomphLibError("Default (empty) version called",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }

    /// Calculate the elemental averages
    void calculate_averages()
    {
      this->calculate_element_averages(this->Average_value);
    }

    /// Return the average values
    double& average_value(const unsigned& i)
    {
      if (Average_value == 0)
      {
        throw OomphLibError("Averages not calculated yet",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
      return Average_value[i];
    }


    /// Return the average values
    const double& average_value(const unsigned& i) const
    {
      if (Average_value == 0)
      {
        throw OomphLibError("Averages not calculated yet",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
      return Average_value[i];
    }
  };

  class DGMesh : public Mesh
  {
  public:
    static double FaceTolerance;

    DGMesh() : Mesh() {}

    virtual ~DGMesh() {}

    virtual void neighbour_finder(FiniteElement* const& bulk_element_pt,
                                  const int& face_index,
                                  const Vector<double>& s_bulk,
                                  FaceElement*& face_element_pt,
                                  Vector<double>& s_face)
    {
      std::string error_message = "Empty neighbour_finder() has been called.\n";
      error_message +=
        "This function is implemented in the base class of a DGMesh.\n";
      error_message += "It must be overloaded in a specific DGMesh\n";

      throw OomphLibError(
        error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }


    // Setup the face information for all the elements
    /// If the boolean flag is set to true then the data from neighbouring
    /// faces will be added as external data to the bulk element.
    void setup_face_neighbour_info(
      const bool& add_face_data_as_external = false)
    {
      // Loop over all the elements and setup their face neighbour information
      const unsigned n_element = this->nelement();
      for (unsigned e = 0; e < n_element; e++)
      {
        dynamic_cast<DGElement*>(this->element_pt(e))
          ->setup_face_neighbour_info(add_face_data_as_external);
      }
    }

    // Limit the slopes on the entire mesh
    void limit_slopes(SlopeLimiter* const& slope_limiter_pt)
    {
      // Loop over all the elements and calculate the averages
      const unsigned n_element = this->nelement();
      for (unsigned e = 0; e < n_element; e++)
      {
        dynamic_cast<DGElement*>(this->element_pt(e))->calculate_averages();
      }

      // Now loop over again and limit the values
      for (unsigned e = 0; e < n_element; e++)
      {
        dynamic_cast<DGElement*>(this->element_pt(e))
          ->slope_limit(slope_limiter_pt);
      }
    }
  };


  //======================================================
  /// Base class for slope limiters
  //=====================================================
  class SlopeLimiter
  {
  public:
    /// Empty constructor
    SlopeLimiter() {}

    /// virtual destructor
    virtual ~SlopeLimiter() {}

    /// Basic function
    virtual void limit(const unsigned& i,
                       const Vector<DGElement*>& required_element_pt)
    {
      throw OomphLibError("Calling default empty limiter\n",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
  };

  class MinModLimiter : public SlopeLimiter
  {
    /// Constant that is used in the modified min-mod function to
    /// provide better properties at extrema
    double M;

    /// Boolean flag to indicate a MUSCL or straight min-mod limiter
    bool MUSCL;

  public:
    /// Constructor takes a value for the modification parameter M
    /// (default to zero --- classic min mod) and a flag to indicate whether
    /// we use MUSCL limiting or not --- default false
    MinModLimiter(const double& m = 0.0, const bool& muscl = false)
      : SlopeLimiter(), M(m), MUSCL(muscl)
    {
    }

    /// Empty destructor
    virtual ~MinModLimiter() {}

    /// The basic minmod function
    double minmod(Vector<double>& args);


    /// The modified  minmod function
    double minmodB(Vector<double>& args, const double& h);

    /// The limit function
    void limit(const unsigned& i,
               const Vector<DGElement*>& required_element_pt);
  };


} // namespace oomph
#endif
