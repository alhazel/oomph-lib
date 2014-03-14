//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//           Version 0.90. August 3, 2009.
//LIC// 
//LIC// Copyright (C) 2006-2009 Matthias Heil and Andrew Hazel
//LIC// 
//LIC// This library is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU Lesser General Public
//LIC// License as published by the Free Software Foundation; either
//LIC// version 2.1 of the License, or (at your option) any later version.
//LIC// 
//LIC// This library is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//LIC// Lesser General Public License for more details.
//LIC// 
//LIC// You should have received a copy of the GNU Lesser General Public
//LIC// License along with this library; if not, write to the Free Software
//LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//LIC// 02110-1301  USA.
//LIC// 
//LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
//LIC// 
//LIC//====================================================================
//Header file for elements that are used to integrate fluid tractions
//This includes the guts (i.e. equations) because we want to inline them
//for faster operation, although it slows down the compilation!
#ifndef OOMPH_AXISYM_FLUID_TRACTION_ELEMENTS_HEADER
#define OOMPH_AXISYM_FLUID_TRACTION_ELEMENTS_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif


//OOMPH-LIB headers
#include "../generic/shape.h"
#include "../generic/elements.h"
#include "../generic/element_with_external_element.h"


namespace oomph
{
 

//=======================================================================
/// Namespace containing the zero traction function for axisymmetric
/// Navier Stokes traction elements
//=======================================================================
namespace AxisymmetricNavierStokesTractionElementHelper
 {

  //=======================================================================
  /// Default load function (zero traction)
  //=======================================================================
  extern void Zero_traction_fct(const double& time,
                                const Vector<double> &x,
                                const Vector<double>& N,
                                Vector<double>& load);
 
 }


//======================================================================
/// A class for elements that allow the imposition of an applied traction
/// in the axisym Navier Stokes eqns.
/// The geometrical information can be read from the FaceGeometry<ELEMENT> 
/// class and and thus, we can be generic enough without the need to have
/// a separate equations class.
//======================================================================
template <class ELEMENT>
class AxisymmetricNavierStokesTractionElement : 
public virtual FaceGeometry<ELEMENT>, 
 public virtual FaceElement
 {
   protected:
  
  /// Index at which the i-th velocity component is stored
  Vector<unsigned> U_index_axisymmetric_nst_traction;
  
 /// \short Pointer to an imposed traction function. Arguments:
 /// Eulerian coordinate; outer unit normal;
 /// applied traction. (Not all of the input arguments will be
 /// required for all specific load functions but the list should
 /// cover all cases)
  void (*Traction_fct_pt)(const double& time,
                          const Vector<double> &x, 
                          const Vector<double> &n,
                          Vector<double> &result);
 
 
 /// \short Get the traction vector: Pass number of integration point (dummy), 
 /// Eulerian coordinate and normal vector and return the load vector
 /// (not all of the input arguments will be
 /// required for all specific load functions but the list should
 /// cover all cases). This function is virtual so it can be 
 /// overloaded for FSI.
 virtual void get_traction(const double& time,
                           const unsigned& intpt,
                           const Vector<double>& x,
                           const Vector<double>& n,
                           Vector<double>& traction) const
 {
  Traction_fct_pt(time,x,n,traction);
 }
 
 
 /// \short Helper function that actually calculates the residuals
 // This small level of indirection is required to avoid calling
 // fill_in_contribution_to_residuals in fill_in_contribution_to_jacobian
 // which causes all kinds of pain if overloading later on
 void fill_in_contribution_to_residuals_axisymmetric_nst_traction(
  Vector<double> &residuals);
 
 
  public:
 
 /// \short Constructor, which takes a "bulk" element and the 
 /// value of the index and its limit
 AxisymmetricNavierStokesTractionElement
  (FiniteElement* const &element_pt, const int &face_index) : 
 FaceGeometry<ELEMENT>(), FaceElement()
  {    
   //Attach the geometrical information to the element. N.B. This function
   //also assigns nbulk_value from the required_nvalue of the bulk element
   element_pt->build_face_element(face_index,this);
   
   //Find the dimension of the problem
   unsigned n_dim = element_pt->nodal_dimension();
   
   //Find the index at which the velocity unknowns are stored
   ELEMENT* cast_element_pt = dynamic_cast<ELEMENT*>(element_pt);
   this->U_index_axisymmetric_nst_traction.resize(n_dim+1);
   for(unsigned i=0;i<n_dim+1;i++)
    {
     this->U_index_axisymmetric_nst_traction[i] = 
      cast_element_pt->u_index_axi_nst(i);
    }
   
   // Zero traction
   Traction_fct_pt=
    &AxisymmetricNavierStokesTractionElementHelper::
    Zero_traction_fct;
  }
 
 
 /// Reference to the traction function pointer
 void (* &traction_fct_pt())(const double& time,
                             const Vector<double>& x,
                             const Vector<double>& n,
                             Vector<double>& traction)
  {return Traction_fct_pt;}
 
 
 /// Return the residuals
 void fill_in_contribution_to_residuals(Vector<double> &residuals)
 {
  fill_in_contribution_to_residuals_axisymmetric_nst_traction(residuals);
 }
 
 
 
 /// Fill in contribution from Jacobian
 void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                       DenseMatrix<double> &jacobian)
 {
  //Call the residuals
  fill_in_contribution_to_residuals_axisymmetric_nst_traction(residuals);
 }
 
 /// Specify the value of nodal zeta from the face geometry
 /// \short The "global" intrinsic coordinate of the element when
 /// viewed as part of a geometric object should be given by
 /// the FaceElement representation, by default (needed to break
 /// indeterminacy if bulk element is SolidElement)
 double zeta_nodal(const unsigned &n, const unsigned &k,           
                   const unsigned &i) const 
 {return FaceElement::zeta_nodal(n,k,i);}     

 /// \short Output function
 void output(std::ostream &outfile)
  {
   unsigned nplot=5;
   output(outfile,nplot);
  }
 
 /// \short Number of scalars/fields output by this element. Reimplements
 /// broken virtual function in base class.
 unsigned nscalar_paraview() const
 {
  //Number of dimensions
  unsigned n_dim = this->nodal_dimension();

  return 2*(n_dim+1);
 }
 
 /// \short Write values of the k-th scalar field at the plot points. Needs 
 /// to be implemented for each new specific element type.
 void scalar_value_paraview(std::ofstream& file_out,
                            const unsigned& k,
                            const unsigned& nplot) const
 {
  //Number of dimensions
   unsigned n_dim = this->nodal_dimension();

   //Find out how many nodes there are
   const unsigned n_node = nnode();

   // Get continuous time from timestepper of first node
   double time=node_pt(0)->time_stepper_pt()->time_pt()->time();
  
   //Set up memory for the shape functions
   Shape psi(n_node);

   // Local and global coordinates
   Vector<double> s(n_dim-1);
   Vector<double> interpolated_x(n_dim);

   // Loop over plot points
   unsigned num_plot_points=this->nplot_points_paraview(nplot);
   for (unsigned iplot=0;iplot<num_plot_points;iplot++)
    {
     // Get local coordinates of plot point
     this->get_s_plot(iplot,nplot,s);

     // Outer unit normal
     Vector<double> unit_normal(n_dim);
     outer_unit_normal(s,unit_normal);

     //Find the shape functions
     shape(s,psi);

     //Initialise to zero
     for(unsigned i=0;i<n_dim;i++)
      {
       interpolated_x[i] = 0.0;
      }

     //Calculate stuff
     for(unsigned l=0;l<n_node;l++) 
      {
       //Loop over directions
       for(unsigned i=0;i<n_dim;i++)
        {
         interpolated_x[i] += this->nodal_position(l,i)*psi[l];
        }
      }

     //Get the imposed traction
     Vector<double> traction(3);

     //Dummy integration point
     unsigned ipt=0;
     get_traction(time,ipt,interpolated_x,unit_normal,traction);

     // Traction components
     if (k<n_dim+1) 
      {
       file_out << traction[k] << std::endl;
      }
     // Advection Diffusion 
     else if (k<2*n_dim+1 && k>=n_dim+1) 
      {
       file_out << unit_normal[k] << std::endl;
      }
     // Never get here
     else
      {
       std::stringstream error_stream;
       error_stream
        << "Axisymmetric Fluid Traction Navier-Stokes Elements only store "
        << 2*(n_dim+1) << " fields " << std::endl;
       throw OomphLibError(
        error_stream.str(),
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
      }
    }
 }

  /// \short Name of the i-th scalar field. Default implementation
 /// returns V1 for the first one, V2 for the second etc. Can (should!) be
 /// overloaded with more meaningful names in specific elements.
 std::string scalar_name_paraview(const unsigned& i) const
 {
  //Number of dimensions
  unsigned n_dim = this->nodal_dimension();
  
  // Traction components
  if (i<n_dim+1) 
   {
    return "Traction component "+StringConversion::to_string(i);
   }
  // Normals
  else if (i<2*n_dim+1 && i>=n_dim+1) 
   {
    return "Normal "+StringConversion::to_string(i%(n_dim+1));
   }
  // Never get here
  else
   {
    std::stringstream error_stream;
    error_stream
     << "Axisymmetric Fluid Traction Navier-Stokes Elements only store "
     << 2*(n_dim+1) << " fields " << std::endl;
    throw OomphLibError(
     error_stream.str(),
     OOMPH_CURRENT_FUNCTION,
     OOMPH_EXCEPTION_LOCATION);
   }
 }
 
 /// \short Output function
 void output(std::ostream &outfile, const unsigned &n_plot)
  {
   //Number of dimensions
   unsigned n_dim = this->nodal_dimension();

   //Find out how many nodes there are
   const unsigned n_node = nnode();

   // Get continuous time from timestepper of first node
   double time=node_pt(0)->time_stepper_pt()->time_pt()->time();
  
   //Set up memory for the shape functions
   Shape psi(n_node);

   // Local and global coordinates
   Vector<double> s(n_dim-1);
   Vector<double> interpolated_x(n_dim);

   // Tecplot header info
   outfile << this->tecplot_zone_string(n_plot);

   // Loop over plot points
   unsigned num_plot_points=this->nplot_points(n_plot);
   for (unsigned iplot=0;iplot<num_plot_points;iplot++)
    {
     // Get local coordinates of plot point
     this->get_s_plot(iplot,n_plot,s);

     // Outer unit normal
     Vector<double> unit_normal(n_dim);
     outer_unit_normal(s,unit_normal);

     //Find the shape functions
     shape(s,psi);

     //Initialise to zero
     for(unsigned i=0;i<n_dim;i++)
      {
       interpolated_x[i] = 0.0;
      }

     //Calculate stuff
     for(unsigned l=0;l<n_node;l++) 
      {
       //Loop over directions
       for(unsigned i=0;i<n_dim;i++)
        {
         interpolated_x[i] += this->nodal_position(l,i)*psi[l];
        }
      }

     //Get the imposed traction
     Vector<double> traction(3);

     //Dummy integration point
     unsigned ipt=0;
     get_traction(time,ipt,interpolated_x,unit_normal,traction);

     //Output the x,y,..
     for(unsigned i=0;i<n_dim;i++) 
      {
       outfile << interpolated_x[i] << " ";
      }

     //Output the traction components
     for(unsigned i=0;i<n_dim+1;i++)
      {
       outfile << traction[i] << " ";
      }

     // Output normal
     for(unsigned i=0;i<n_dim;i++) 
      {
       outfile << unit_normal[i] << " ";
      } 
     outfile << std::endl;

    }
  }
 
 /// \short C_style output function
 void output(FILE* file_pt)
 {FiniteElement::output(file_pt);}
 
 /// \short C-style output function
 void output(FILE* file_pt, const unsigned &n_plot)
 {FiniteElement::output(file_pt,n_plot);}
 
 
 /// \short Compute traction vector at specified local coordinate
 /// Should only be used for post-processing; ignores dependence
 /// on integration point!
 void traction(const double &time,
               const Vector<double>& s, 
               Vector<double>& traction);
 
 }; 

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

//=====================================================================
/// Compute traction vector at specified local coordinate
/// Should only be used for post-processing; ignores dependence
/// on integration point!
//=====================================================================
template<class ELEMENT>
void AxisymmetricNavierStokesTractionElement<ELEMENT>::
 traction(const double& time, const Vector<double>& s, Vector<double>& traction)
 {
  unsigned n_dim = this->nodal_dimension();
  
  // Position vector
  Vector<double> x(n_dim);
  interpolated_x(s,x);
  
  // Outer unit normal (only in r and z direction!)
  Vector<double> unit_normal(n_dim);
  outer_unit_normal(s,unit_normal);
  
  // Dummy
  unsigned ipt=0;
  
  // Traction vector
  get_traction(time,ipt,x,unit_normal,traction);
  
 }


//=====================================================================
/// Return the residuals for the 
/// AxisymmetricNavierStokesTractionElement equations
//=====================================================================
template<class ELEMENT>
 void AxisymmetricNavierStokesTractionElement<ELEMENT>::
 fill_in_contribution_to_residuals_axisymmetric_nst_traction(
  Vector<double> &residuals)
 {

  //Find out how many nodes there are
  unsigned n_node = nnode();
    
  // Get continuous time from timestepper of first node
  double time=node_pt(0)->time_stepper_pt()->time_pt()->time();
  
 #ifdef PARANOID
  //Find out how many positional dofs there are
  unsigned n_position_type = this->nnodal_position_type();  
  if(n_position_type != 1)
   {
    throw OomphLibError(
     "AxisymmetricNavierStokes is not yet implemented for more than one position type",
     OOMPH_CURRENT_FUNCTION,
     OOMPH_EXCEPTION_LOCATION);
   }
#endif
  
  //Find out the dimension of the node
  unsigned n_dim = this->nodal_dimension();
  
  //Cache the nodal indices at which the velocity components are stored
  unsigned u_nodal_index[n_dim+1];
  for(unsigned i=0;i<n_dim+1;i++)
   {
    u_nodal_index[i] = this->U_index_axisymmetric_nst_traction[i];
   }
  
  //Integer to hold the local equation number
  int local_eqn=0;
  
  //Set up memory for the shape functions
  //Note that in this case, the number of lagrangian coordinates is always
  //equal to the dimension of the nodes
  Shape psi(n_node);
  DShape dpsids(n_node,n_dim-1); 
  
  //Set the value of n_intpt
  unsigned n_intpt = integral_pt()->nweight();
  
  //Loop over the integration points
  for(unsigned ipt=0;ipt<n_intpt;ipt++)
   {
    //Get the integral weight
    double w = integral_pt()->weight(ipt);
    
    //Only need to call the local derivatives
    dshape_local_at_knot(ipt,psi,dpsids);
    
    //Calculate the Eulerian and Lagrangian coordinates 
    Vector<double> interpolated_x(n_dim,0.0);
    
    //Also calculate the surface Vectors (derivatives wrt local coordinates)
    DenseMatrix<double> interpolated_A(n_dim-1,n_dim,0.0);   
    
    //Calculate positions and derivatives
    for(unsigned l=0;l<n_node;l++) 
     {
      //Loop over directions
      for(unsigned i=0;i<n_dim;i++)
       {        
        //Calculate the Eulerian coords
        const double x_local = nodal_position(l,i);
        interpolated_x[i] += x_local*psi(l);
        
        //Loop over LOCAL derivative directions, to calculate the tangent(s)
        for(unsigned j=0;j<n_dim-1;j++)
         {
          interpolated_A(j,i) += x_local*dpsids(l,j);
         }
       }
     }
    
    //Now find the local metric tensor from the tangent Vectors
    DenseMatrix<double> A(n_dim-1);
    for(unsigned i=0;i<n_dim-1;i++)
     {
      for(unsigned j=0;j<n_dim-1;j++)
       {
        //Initialise surface metric tensor to zero
        A(i,j) = 0.0;
        
        //Take the dot product
        for(unsigned k=0;k<n_dim;k++)
         { 
          A(i,j) += interpolated_A(i,k)*interpolated_A(j,k);
         }
       }
     }
    
    //Get the outer unit normal
    Vector<double> interpolated_normal(n_dim);
    outer_unit_normal(ipt,interpolated_normal);
    
    //Find the determinant of the metric tensor
    double Adet =0.0;
    switch(n_dim)
     {
     case 2:
      Adet = A(0,0);
      break;
     case 3:
      Adet = A(0,0)*A(1,1) - A(0,1)*A(1,0);
      break;
     default:
      throw 
       OomphLibError(
        "Wrong dimension in AxisymmetricNavierStokesTractionElement",
        "AxisymmetricNavierStokesTractionElement::fill_in_contribution_to_residuals()",
        OOMPH_EXCEPTION_LOCATION);
     }
    
    //Premultiply the weights and the square-root of the determinant of 
    //the metric tensor
    double W = w*sqrt(Adet);

    //Now calculate the load
    Vector<double> traction(n_dim+1);
    get_traction(time,
                 ipt,
                 interpolated_x,
                 interpolated_normal,
                 traction);
    
    //Loop over the test functions, nodes of the element
    for(unsigned l=0;l<n_node;l++)
     {
      //Loop over the velocity components
      for(unsigned i=0;i<n_dim+1;i++)
       {
        // Equation number
        local_eqn = this->nodal_local_eqn(l,u_nodal_index[i]);
        /*IF it's not a boundary condition*/
        if(local_eqn >= 0)
         {
          //Add the loading terms to the residuals
          residuals[local_eqn] -= traction[i]*psi(l)*interpolated_x[0]*W;
         }
       }
     } //End of loop over shape functions
   } //End of loop over integration points
  
 }




/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////


//=======================================================================
/// Namespace containing the default Strouhal number of axisymmetric 
/// linearised FSI.
//=======================================================================
namespace LinearisedFSIAxisymmetricNStNoSlipBCHelper
 {
 
  /// Default for fluid Strouhal number
  extern double Default_strouhal_number;

 }



//======================================================================
/// \short A class for elements that allow the imposition of the linearised
/// FSI no slip condition from an adjacent linearly elastic axisymmetric 
/// solid. The element geometry is obtained from the FaceGeometry<ELEMENT> 
/// policy class.
//======================================================================
 template <class FLUID_BULK_ELEMENT, class SOLID_BULK_ELEMENT>
 class LinearisedFSIAxisymmetricNStNoSlipBCElementElement : 
  public virtual FaceGeometry<FLUID_BULK_ELEMENT>, 
  public virtual FaceElement, 
  public virtual ElementWithExternalElement
 {
  
 public:
  
  
  /// \short Constructor, takes the pointer to the "bulk" element and the 
  /// face index identifying the face to which the element is attached.
  /// The optional identifier can be used
  /// to distinguish the additional nodal values created by 
  /// this element from thos created by other FaceElements.
  LinearisedFSIAxisymmetricNStNoSlipBCElementElement(
   FiniteElement* const &bulk_el_pt,
   const int& face_index, 
   const unsigned &id=0); 
  
  /// Broken copy constructor
  LinearisedFSIAxisymmetricNStNoSlipBCElementElement(
   const LinearisedFSIAxisymmetricNStNoSlipBCElementElement& dummy) 
   { 
    BrokenCopy::broken_copy(
     "LinearisedFSIAxisymmetricNStNoSlipBCElementElement");
   } 
  
  /// Broken assignment operator
  void operator=(const LinearisedFSIAxisymmetricNStNoSlipBCElementElement&) 
   {
    BrokenCopy::broken_assign(
     "LinearisedFSIAxisymmetricNStNoSlipBCElementElement");
   }
  

  /// \short Access function for the pointer to the fluid Strouhal number
  /// (if not set, St defaults to 1)
  double* &st_pt() {return St_pt;}
  
  /// Access function for the fluid Strouhal number
  double st() const 
  {
   return *St_pt;
  }

  /// Add the element's contribution to its residual vector
  inline void fill_in_contribution_to_residuals(Vector<double> &residuals)
  {
   //Call the generic residuals function with flag set to 0
   //using a dummy matrix argument
   fill_in_generic_residual_contribution_fsi_no_slip_axisym(
    residuals,GeneralisedElement::Dummy_matrix,0);
  }
  

  /// \short Add the element's contribution to its residual vector and its
  /// Jacobian matrix
  inline void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                               DenseMatrix<double> &jacobian)
   {
    //Call the generic routine with the flag set to 1
    fill_in_generic_residual_contribution_fsi_no_slip_axisym
     (residuals,jacobian,1);
    
    //Derivatives w.r.t. external data
    fill_in_jacobian_from_external_interaction_by_fd(residuals,jacobian);
   }
  
  /// Output function
  void output(std::ostream &outfile) 
  {
   //Dummy
   unsigned nplot=0;
   output(outfile,nplot);
  }
  
  /// Output function: Output at Gauss points; n_plot is ignored.
  void output(std::ostream &outfile, const unsigned &n_plot)
  {
   outfile << "ZONE\n";
   
   //Get the value of Nintpt
   const unsigned n_intpt = integral_pt()->nweight();
   
   //Set the Vector to hold local coordinates
   Vector<double> s_int(Dim-1);
   
   //Loop over the integration points
   for(unsigned ipt=0;ipt<n_intpt;ipt++)
    {
     
     //Assign values of s
     for(unsigned i=0;i<(Dim-1);i++) 
      {
       s_int[i] = integral_pt()->knot(ipt,i);
      }
     
     // Get boundary coordinate
     Vector<double> zeta(1);
     interpolated_zeta(s_int,zeta);
     
     // Get velocity from adjacent solid
     SOLID_BULK_ELEMENT* ext_el_pt=dynamic_cast<SOLID_BULK_ELEMENT*>(
      external_element_pt(0,ipt));
     Vector<double> s_ext(external_element_local_coord(0,ipt));
     Vector<double> dudt(3);
     ext_el_pt->interpolated_du_dt_axisymmetric_linear_elasticity(s_ext,dudt);
     
     // Output
     outfile << ext_el_pt->interpolated_x(s_ext,0) << " "
             << ext_el_pt->interpolated_x(s_ext,1) << " "
             << dudt[0] << " " 
             << dudt[1] << " " 
             << dudt[2] << " " 
             << zeta[0] << std::endl;
    }
  }
 
  
  /// C-style output function
  void output(FILE* file_pt)
   {FaceGeometry<FLUID_BULK_ELEMENT>::output(file_pt);}

  /// C-style output function
  void output(FILE* file_pt, const unsigned &n_plot)
   {FaceGeometry<FLUID_BULK_ELEMENT>::output(file_pt,n_plot);}

  
  
 protected:
  
  /// \short Function to compute the shape and test functions and to return 
  /// the Jacobian of mapping between local and global (Eulerian)
  /// coordinates
  inline double shape_and_test(const Vector<double> &s, Shape &psi, Shape &test)
   const
   {
    //Find number of nodes
    unsigned n_node = nnode();
    
    //Get the shape functions
    shape(s,psi);
    
    //Set the test functions to be the same as the shape functions
    for(unsigned i=0;i<n_node;i++) {test[i] = psi[i];}
    
    //Return the value of the jacobian
    return J_eulerian(s);
   }
  
  
  /// \short Function to compute the shape and test functions and to return 
 /// the Jacobian of mapping between local and global (Eulerian)
 /// coordinates
  inline double shape_and_test_at_knot(const unsigned &ipt,
                                       Shape &psi, Shape &test)
   const
   {
    //Find number of nodes
    unsigned n_node = nnode();
    
    //Get the shape functions
    shape_at_knot(ipt,psi);
    
    //Set the test functions to be the same as the shape functions
    for(unsigned i=0;i<n_node;i++) {test[i] = psi[i];}
    
    //Return the value of the jacobian
    return J_eulerian_at_knot(ipt);
   }
  

  
 private:
  
  
  /// \short Add the element's contribution to its residual vector.
  /// flag=1(or 0): do (or don't) compute the contribution to the
  /// Jacobian as well. 
  void fill_in_generic_residual_contribution_fsi_no_slip_axisym(
   Vector<double> &residuals, DenseMatrix<double> &jacobian, 
   const unsigned& flag);
  
  ///The spatial dimension of the problem
  unsigned Dim;
  
  ///The index at which the unknowns are stored at the nodes
  Vector<unsigned> U_index_fsi_no_slip_axisym;
  
  /// Lagrange Id
  unsigned Id;

  /// Pointer to fluid Strouhal number
  double* St_pt;

 };
 
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
 
  
 
//===========================================================================
/// Constructor, takes the pointer to the "bulk" element, and the 
/// face index that identifies the face of the bulk element to which
/// this face element is to be attached.
/// The optional identifier can be used
/// to distinguish the additional nodal values created by 
/// this element from thos created by other FaceElements.
//===========================================================================
  template <class FLUID_BULK_ELEMENT, class SOLID_BULK_ELEMENT>
   LinearisedFSIAxisymmetricNStNoSlipBCElementElement<FLUID_BULK_ELEMENT, 
   SOLID_BULK_ELEMENT>::
   LinearisedFSIAxisymmetricNStNoSlipBCElementElement(
    FiniteElement* const &bulk_el_pt, 
    const int &face_index,
    const unsigned &id) : 
  FaceGeometry<FLUID_BULK_ELEMENT>(), FaceElement()
   { 
    // Set source element storage: one interaction with an external element
    // that provides the velocity of the adjacent linear elasticity 
    // element
    this->set_ninteraction(1); 
    
    //  Store the ID of the FaceElement -- this is used to distinguish
    // it from any others
    Id=id;

    // Initialise pointer to fluid Strouhal number. Defaults to 1
    St_pt=&LinearisedFSIAxisymmetricNStNoSlipBCHelper::Default_strouhal_number;
    
    // Let the bulk element build the FaceElement, i.e. setup the pointers 
    // to its nodes (by referring to the appropriate nodes in the bulk
    // element), etc.
    bulk_el_pt->build_face_element(face_index,this);
    
    // Extract the dimension of the problem from the dimension of 
    // the first node
    Dim = this->node_pt(0)->ndim();
    
    //Read the index from the (cast) bulk element.
    U_index_fsi_no_slip_axisym.resize(3);
    for (unsigned i=0;i<3;i++)
     {
      U_index_fsi_no_slip_axisym[i] = 
       dynamic_cast<FLUID_BULK_ELEMENT*>(bulk_el_pt)->u_index_axi_nst(i);
     }

    // We need Dim+1 additional values for each FaceElement node
    // to store the Lagrange multipliers.
    Vector<unsigned> n_additional_values(nnode(), Dim+1);
    
    // Now add storage for Lagrange multipliers and set the map containing 
    // the position of the first entry of this face element's 
    // additional values.
    add_additional_values(n_additional_values,id);
    
   }
 
 
//===========================================================================
/// \short Helper function to compute the element's residual vector and 
/// the Jacobian matrix.
//===========================================================================
template <class FLUID_BULK_ELEMENT, class SOLID_BULK_ELEMENT>
void LinearisedFSIAxisymmetricNStNoSlipBCElementElement<FLUID_BULK_ELEMENT,
 SOLID_BULK_ELEMENT>:: fill_in_generic_residual_contribution_fsi_no_slip_axisym(
  Vector<double> &residuals, DenseMatrix<double> &jacobian, 
  const unsigned& flag)
  {
   //Find out how many nodes there are
   const unsigned n_node = nnode();
   
   //Set up memory for the shape and test functions
   Shape psif(n_node), testf(n_node);
   
   //Set the value of Nintpt
   const unsigned n_intpt = integral_pt()->nweight();
   
   //Set the Vector to hold local coordinates
   Vector<double> s(Dim-1);

   // Cache the Strouhal number
   const double local_st=st();
   
   //Integers to hold the local equation and unknown numbers
   int local_eqn=0;
   
   //Loop over the integration points
   //--------------------------------
   for(unsigned ipt=0;ipt<n_intpt;ipt++)
    {
     
     //Assign values of s
    for(unsigned i=0;i<(Dim-1);i++) {s[i] = integral_pt()->knot(ipt,i);}
    
    //Get the integral weight
    double w = integral_pt()->weight(ipt);
    
    //Find the shape and test functions and return the Jacobian
    //of the mapping
    double J = shape_and_test(s,psif,testf);
 
    //Premultiply the weights and the Jacobian
    double W = w*J;
    
    //Calculate the Lagrange multiplier and the fluid veloc
    Vector<double> lambda(Dim+1,0.0);
    Vector<double> fluid_veloc(Dim+1,0.0);
    
    // Loop over nodes
    for(unsigned j=0;j<n_node;j++) 
     {
      Node* nod_pt=node_pt(j);
      
      // Cast to a boundary node
      BoundaryNodeBase *bnod_pt=dynamic_cast<BoundaryNodeBase*>(node_pt(j));
      
      // Get the index of the first nodal value associated with
      // this FaceElement
      unsigned first_index=
       bnod_pt->index_of_first_value_assigned_by_face_element(Id);
      
      //Assemble
      for(unsigned i=0;i<Dim+1;i++)
       {
        lambda[i]+=nod_pt->value(first_index+i)*psif(j);
        fluid_veloc[i]+=nod_pt->value(U_index_fsi_no_slip_axisym[i])*psif(j);
       }  
     }
    
    // Get velocity from adjacent solid
    SOLID_BULK_ELEMENT* ext_el_pt=dynamic_cast<SOLID_BULK_ELEMENT*>(
     external_element_pt(0,ipt));
    Vector<double> s_ext(external_element_local_coord(0,ipt));
    Vector<double> dudt(3);
    ext_el_pt->interpolated_du_dt_axisymmetric_linear_elasticity(s_ext,dudt);
    
    
    //Now add to the appropriate equations
    
    //Loop over the test functions
    for(unsigned l=0;l<n_node;l++)
     {
      
      // Loop over directions
      for (unsigned i=0;i<Dim+1;i++)
       {
        
        // Add contribution to bulk Navier Stokes equations where
        //-------------------------------------------------------
        // the Lagrange multiplier acts as a traction
        //-------------------------------------------
        local_eqn = nodal_local_eqn(l,U_index_fsi_no_slip_axisym[i]);
        
        /*IF it's not a boundary condition*/
        if(local_eqn >= 0)
         {
          //Add the Lagrange multiplier "traction" to the bulk
          residuals[local_eqn] -= lambda[i]*testf[l]*W;
          

          //Jacobian entries
          if(flag)
           {
            //Loop over the lagrange multiplier unknowns
            for(unsigned l2=0;l2<n_node;l2++)
             {
              // Cast to a boundary node
              BoundaryNodeBase *bnod_pt = 
               dynamic_cast<BoundaryNodeBase*>(node_pt(l2));
              
              // Local unknown
              int local_unknown=nodal_local_eqn
               (l2,bnod_pt->
                index_of_first_value_assigned_by_face_element(Id)+i); 
              
              // If it's not pinned
              if (local_unknown>=0)
               {
                jacobian(local_eqn,local_unknown) -=
                 psif[l2]*testf[l]*W;
               }
             }
           }
         }
         
        // Now do the Lagrange multiplier equations
        //-----------------------------------------
        // Cast to a boundary node
        BoundaryNodeBase *bnod_pt = 
         dynamic_cast<BoundaryNodeBase*>(node_pt(l));
        
        // Local eqn number:   
        int local_eqn=nodal_local_eqn
         (l,bnod_pt->index_of_first_value_assigned_by_face_element(Id)+i); 
        
        // If it's not pinned
        if (local_eqn>=0)
         {
          residuals[local_eqn]+=(local_st*dudt[i]-fluid_veloc[i])*testf(l)*W;
          
          //Jacobian entries
          if(flag)
           {
            // Loop over the velocity unknowns [derivs w.r.t. to 
            // wall velocity taken care of by fd-ing
            for(unsigned l2=0;l2<n_node;l2++)
             {
              int local_unknown = 
               nodal_local_eqn(l2,U_index_fsi_no_slip_axisym[i]);
              
              /*IF it's not a boundary condition*/
              if(local_unknown >= 0)
               {
                jacobian(local_eqn,local_unknown) -=
                 psif[l2]*testf[l]*W;
               }
             }
           }

         }
        
       }
     }
    }
  }
    
}



#endif
