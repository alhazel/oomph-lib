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
// Header file for specific surface elements

#ifndef OOMPH_NAVIER_STOKES_SURFACE_DRAG_TORQUE_ELEMENTS_HEADER
#define OOMPH_NAVIER_STOKES_SURFACE_DRAG_TORQUE_ELEMENTS_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

namespace oomph
{
  //======================================================================
  /// A class of elements that allow the determination of the
  /// drag and toque, relative to a given centre of rotation, along
  /// a domain boundary.
  /// The element operates as a FaceElement and attaches itself
  /// to a bulk element of the type specified by the template
  /// argument.
  //======================================================================
  template<class ELEMENT>
  class NavierStokesSurfaceDragTorqueElement
    : public virtual FaceGeometry<ELEMENT>,
      public virtual FaceElement,
      public virtual ElementWithDragFunction
  {
  public:
    /// Constructor, which takes a "bulk" element and the value of an
    /// index describing to which face the element should be attached.
    NavierStokesSurfaceDragTorqueElement(FiniteElement* const& element_pt,
                                         const int& face_index)
      : FaceGeometry<ELEMENT>(), FaceElement()
    {
      // Attach the geometrical information to the element. N.B. This function
      // also assigns nbulk_value from the required_nvalue of the bulk element
      element_pt->build_face_element(face_index, this);

      // Set the dimension from the dimension of the first node
      this->Dim = node_pt(0)->ndim();

      // Default centre of rotation is the origin
      this->Centre_of_rotation.resize(this->Dim, 0.0);
    }

    /// Set the translation and rotation of the rigid object
    /// as external data
    void set_translation_and_rotation(Data* const& object_data_pt)
    {
      this->Translation_index = this->add_external_data(object_data_pt);
    }


    /// Access function for the centre of rotation
    double& centre_of_rotation(const unsigned& i)
    {
      return this->Centre_of_rotation[i];
    }


    /// Function that specifies the drag force and the torque about
    /// the origin
    virtual void get_drag_and_torque(Vector<double>& drag_force,
                                     Vector<double>& drag_torque)
    {
      // Spatial dimension of element
      unsigned ndim = dim();

      // Initialise force
      for (unsigned i = 0; i < ndim + 1; i++)
      {
        drag_force[i] = 0.0;
      }

      // Now there will be one torque component in 2D (1D surface element)
      // and three in 3D (2D surface element)
      for (unsigned i = 0; i < 2 * ndim - 1; i++)
      {
        drag_torque[i] = 0.0;
      }

      // Vector of local coordinates in face element
      Vector<double> s(ndim);

      // Vector for global Eulerian coordinates
      Vector<double> x(ndim + 1);

      // Vector for local coordinates in bulk element
      Vector<double> s_bulk(ndim + 1);

      // Set the value of n_intpt
      unsigned n_intpt = integral_pt()->nweight();

      // Get pointer to assocated bulk element
      ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(bulk_element_pt());

      // Loop over the integration points
      for (unsigned ipt = 0; ipt < n_intpt; ipt++)
      {
        // Assign values of s in FaceElement and local coordinates in bulk
        // element
        for (unsigned i = 0; i < ndim; i++)
        {
          s[i] = integral_pt()->knot(ipt, i);
        }

        // Get the bulk coordinates
        this->get_local_coordinate_in_bulk(s, s_bulk);

        // Get the integral weight
        double w = integral_pt()->weight(ipt);

        // Get jacobian of mapping
        double J = J_eulerian(s);

        // Premultiply the weights and the Jacobian
        double W = w * J;

        // Get x position as Vector
        interpolated_x(s, x);

#ifdef PARANOID

        // Get x position as Vector from bulk element
        Vector<double> x_bulk(ndim + 1);
        bulk_el_pt->interpolated_x(s_bulk, x_bulk);

        double max_legal_error = 1.0e-5;
        double error = 0.0;
        for (unsigned i = 0; i < ndim + 1; i++)
        {
          error += fabs(x[i] - x_bulk[i]);
        }
        if (error > max_legal_error)
        {
          std::ostringstream error_stream;
          error_stream << "difference in Eulerian posn from bulk and face: "
                       << error << " exceeds threshold " << max_legal_error
                       << std::endl;
          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif

        // Outer unit normal of the fluid
        Vector<double> normal(ndim + 1);
        outer_unit_normal(s, normal);

        // Get velocity from bulk element
        Vector<double> veloc(ndim + 1);
        bulk_el_pt->interpolated_u_nst(s_bulk, veloc);

        // Get traction from bulk element this is directed into the fluid
        // so we need to reverse the sign
        Vector<double> traction(ndim + 1);
        bulk_el_pt->get_traction(s_bulk, normal, traction);

        // Integrate (note the minus sign)
        for (unsigned i = 0; i < ndim + 1; i++)
        {
          drag_force[i] -= traction[i] * W;
        }

        // Now Calculate the torque which is r x F
        // Scale X relative to the centre of rotation
        Vector<double> X(ndim + 1);
        for (unsigned i = 0; i < ndim + 1; i++)
        {
          X[i] = x[i] - Centre_of_rotation[i] -
                 this->external_data_pt(Translation_index)->value(i);
        }

        // In 2D it's just a single scalar
        // again note the minus sign
        if (ndim == 1)
        {
          drag_torque[0] -= (X[0] * traction[1] - X[1] * traction[0]) * W;
        }
        else if (ndim == 2)
        {
          drag_torque[0] -= (X[1] * traction[2] - X[2] * traction[1]) * W;
          drag_torque[1] -= (X[2] * traction[0] - X[0] * traction[2]) * W;
          drag_torque[2] -= (X[0] * traction[1] - X[1] * traction[0]) * W;
        }
        else
        {
          throw OomphLibError("Dimension of a surface element must be 1 or 2\n",
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      }

      // std::cout << "Forces " << drag_force[0] << " " << drag_force[1] << " "
      //          << drag_torque[0] << "\n";
    }


    /// Specify the value of nodal zeta from the face geometry
    /// The "global" intrinsic coordinate of the element when
    /// viewed as part of a geometric object should be given by
    /// the FaceElement representation, by default (needed to break
    /// indeterminacy if bulk element is SolidElement)
    double zeta_nodal(const unsigned& n,
                      const unsigned& k,
                      const unsigned& i) const
    {
      return FaceElement::zeta_nodal(n, k, i);
    }


    /// Output function
    void output(std::ostream& outfile, const unsigned& n_plot)
    {
      // Get pointer to assocated bulk element
      ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(bulk_element_pt());

      // Elemental dimension
      unsigned dim_el = dim();

      // Local coordinates
      Vector<double> s(dim_el);

      Vector<double> s_bulk(dim_el + 1);

      // Calculate the Eulerian coordinates and Lagrange multiplier
      Vector<double> x(dim_el + 1, 0.0);

      // Outer unit normal of the fluid
      Vector<double> normal(dim_el + 1);

      // Velocity from bulk element
      Vector<double> veloc(dim_el + 1);

      // Tractions (from bulk element)
      Vector<double> traction(dim_el + 1);

      // Tecplot header info
      outfile << this->tecplot_zone_string(n_plot);

      // Loop over plot points
      unsigned num_plot_points = this->nplot_points(n_plot);
      for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
      {
        // Get local coordinates of plot point
        this->get_s_plot(iplot, n_plot, s);

        this->get_local_coordinate_in_bulk(s, s_bulk);

        // Get x position from bulk
        bulk_el_pt->interpolated_x(s_bulk, x);
        // Get outer unit normal
        outer_unit_normal(s, normal);
        bulk_el_pt->interpolated_u_nst(s_bulk, veloc);
        // Get traction from bulk element this is directed into the fluid
        // so we need to reverse the sign
        bulk_el_pt->get_traction(s_bulk, normal, traction);

        // Now Calculate the torque which is r x F
        // Scale X relative to the centre of rotation
        Vector<double> X(dim_el + 1);
        for (unsigned i = 0; i < dim_el + 1; i++)
        {
          X[i] = x[i] - Centre_of_rotation[i] -
                 this->external_data_pt(Translation_index)->value(i);
        }

        for (unsigned i = 0; i < dim_el + 1; i++)
        {
          outfile << x[i] << "  ";
        }

        for (unsigned i = 0; i < dim_el + 1; i++)
        {
          outfile << normal[i] << " ";
        }

        for (unsigned i = 0; i < dim_el + 1; i++)
        {
          outfile << veloc[i] << " ";
        }

        for (unsigned i = 0; i < dim_el + 1; i++)
        {
          outfile << -1.0 * traction[i] << " ";
        }
        outfile << -(X[0] * traction[1] - X[1] * traction[0]);

        outfile << std::endl;
      }

      this->write_tecplot_zone_footer(outfile, n_plot);
    }


  private:
    /// The highest dimension of the problem
    unsigned Dim;

    /// The centre of rotation for the torque calculation
    Vector<double> Centre_of_rotation;

    /// The index of where the translation and rotation data is stored
    unsigned Translation_index;
  };


} // namespace oomph

#endif
