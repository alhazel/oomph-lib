//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//           Version 0.85. June 9, 2008.
//LIC// 
//LIC// Copyright (C) 2006-2008 Matthias Heil and Andrew Hazel
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
//Header for a multi-physics problem that couples a Navier--Stokes
//mesh to an advection diffusion mesh, giving Boussinesq convection

//Oomph-lib headers, we require the generic, advection-diffusion
//and navier-stokes elements.
#include "generic.h"
#include "advection_diffusion.h"
#include "navier_stokes.h"

// Use the oomph and std namespaces 
using namespace oomph;
using namespace std;


/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////


//======================nst_bous_class=================================
/// Build a RefineableQCrouzeixRaviart element that inherits from 
/// ElementWithExternalElement so that it can "communicate" with 
/// an advection diffusion element that provides the temperature
/// in the body force term.
//=====================================================================
template<unsigned DIM>
class RefineableQCrouzeixRaviartBoussinesqElement : 
 public virtual RefineableQCrouzeixRaviartElement<DIM>,
 public virtual ElementWithExternalElement
{

public: 

 /// \short Constructor: call the underlying constructors and 
 /// initialise the pointer to the Rayleigh number to point
 /// to the default value of 0.0.
 RefineableQCrouzeixRaviartBoussinesqElement() : 
  RefineableQCrouzeixRaviartElement<DIM>(),
  ElementWithExternalElement()
  {
   Ra_pt = &Default_Physical_Constant_Value;

   //There is one interaction: The effect of the advection-diffusion
   //element onto the buoyancy term
   this->set_ninteraction(1);
  } 

 ///Access function for the Rayleigh number (const version)
 const double &ra() const {return *Ra_pt;}

 ///Access function for the pointer to the Rayleigh number
 double* &ra_pt() {return Ra_pt;}

 /// \short Call the underlying single-physics element's further_build()
 /// functions and make sure that the pointer to the Rayleigh number
 /// is passed to the sons
 void further_build()
  {
   RefineableQCrouzeixRaviartElement<DIM>::further_build();

   //Cast the pointer to the father element to the specific
   //element type
   RefineableQCrouzeixRaviartBoussinesqElement<DIM>* 
    cast_father_element_pt
    = dynamic_cast<RefineableQCrouzeixRaviartBoussinesqElement<DIM>*>(
     this->father_element_pt());

   //Set the pointer to the Rayleigh number to be the same as that in
   //the father
   this->Ra_pt = cast_father_element_pt->ra_pt();
  }

 /// \short Overload get_body_force_nst() to return the temperature-dependent
 /// buoyancy force, using the temperature computed by the 
 /// "external" advection diffusion element associated with 
 /// integration point \c ipt.
 void get_body_force_nst(const double& time, const unsigned& ipt, 
                         const Vector<double> &s, const Vector<double> &x, 
                         Vector<double> &body_force)
 {

  // Set interaction index -- there's only one interaction...
  const unsigned interaction=0;
  
  // Get a pointer to the external element that computes the
  // the temperature -- we know it's an advection diffusion element.
  const AdvectionDiffusionEquations<DIM>* adv_diff_el_pt=
   dynamic_cast<AdvectionDiffusionEquations<DIM>*>(
    external_element_pt(interaction,ipt));
  
  // Get the temperature interpolated from the external element
  const double interpolated_t =adv_diff_el_pt->
   interpolated_u_adv_diff(external_element_local_coord(interaction,ipt));
  
  // Get vector that indicates the direction of gravity from
  // the Navier-Stokes equations
  Vector<double> gravity(NavierStokesEquations<DIM>::g());
  
  // Set the temperature-dependent body force:
  for (unsigned i=0;i<DIM;i++)
   {
    body_force[i] = -gravity[i]*interpolated_t*ra();
   }

 } // end overloaded body force


 
// hierher

/*  /// Fill in the constituent elements' contribution to the residual vector. */
/*  void fill_in_contribution_to_residuals(Vector<double> &residuals) */
/*   { */
/*    //Call the residuals of the Navier-Stokes equations */
/*    RefineableNavierStokesEquations<DIM>::fill_in_contribution_to_residuals( */
/*     residuals); */
/*   } */

 /// \short Compute the element's residual vector and the Jacobian matrix.
 void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                       DenseMatrix<double> &jacobian)
 {

   //Get the analytical contribution from the basic Navier-Stokes element
   RefineableQCrouzeixRaviartElement<DIM>::
    fill_in_contribution_to_jacobian(residuals,jacobian);
   
#ifdef USE_FD_FOR_DERIVATIVES_WRT_EXTERNAL_DATA

   //Get the off-diagonal terms by finite differencing
   this->fill_in_jacobian_from_external_interaction_by_fd(residuals,jacobian);

#else

   //Get the off-diagonal terms analytically
   this->fill_in_off_diagonal_block_analytic(residuals,jacobian);

#endif

  }



 /// Add the element's contribution to its residuals vector,
 /// jacobian matrix and mass matrix
 void fill_in_contribution_to_jacobian_and_mass_matrix(
  Vector<double> &residuals, DenseMatrix<double> &jacobian, 
  DenseMatrix<double> &mass_matrix)
  {
   //Call the standard (Broken) function
   //which will prevent these elements from being used
   //in eigenproblems until replaced.
   FiniteElement::fill_in_contribution_to_jacobian_and_mass_matrix(
     residuals,jacobian,mass_matrix);
  }


 /// \short Fill in the derivatives of the body force with respect to the
 /// external unknowns
 void get_dbody_force_nst_dexternal_element_data(
  const unsigned& ipt, 
  DenseMatrix<double> &result, Vector<unsigned> &global_eqn_number);


 /// \short Compute the contribution of the external
 /// degrees of freedom (temperatures) on the Navier-Stokes equations
 void fill_in_off_diagonal_block_analytic(Vector<double> &residuals,
                                          DenseMatrix<double> &jacobian)
  {
   //Local storage for the index in the nodes at which the
   //Navier-Stokes velocities are stored (we know that this should be 0,1,2)
   unsigned u_nodal_nst[DIM];
   for(unsigned i=0;i<DIM;i++) 
    {u_nodal_nst[i] = this->u_index_nst(i);}

   //Find out how many nodes there are
   const unsigned n_node = this->nnode();
   
   //Set up memory for the shape and test functions and their derivatives
   Shape psif(n_node), testf(n_node);
   DShape dpsifdx(n_node,DIM), dtestfdx(n_node,DIM);
   
   //Number of integration points
   const unsigned n_intpt = this->integral_pt()->nweight();
   
   //Integers to store the local equations and unknowns
   int local_eqn=0, local_unknown=0;
   
   // Local storage for pointers to hang_info objects
   HangInfo *hang_info_pt=0;   

   //Loop over the integration points
   for(unsigned ipt=0;ipt<n_intpt;ipt++)
    {
     //Get the integral weight
     double w = this->integral_pt()->weight(ipt);
     
     //Call the derivatives of the shape and test functions
     double J = 
      this->dshape_and_dtest_eulerian_at_knot_nst(ipt,psif,dpsifdx,
                                                  testf,dtestfdx);
     
     //Premultiply the weights and the Jacobian
     double W = w*J;
     
     //Assemble the jacobian terms
     
     //Get the derivatives of the body force wrt the unknowns
     //of the external element
     DenseMatrix<double> dbody_dexternal_element_data;

     //Vector of global equation number corresponding to the external
     //element's data
     Vector<unsigned> global_eqn_number_of_external_element_data;

     //Get the appropriate derivatives
     this->get_dbody_force_nst_dexternal_element_data(
      ipt,dbody_dexternal_element_data,
      global_eqn_number_of_external_element_data);

     //Find out how many external data there are
     const unsigned n_external_element_data = 
      global_eqn_number_of_external_element_data.size();

     //Loop over the test functions
     for(unsigned l=0;l<n_node;l++)
      {
       //Assemble the contributions of the temperature to 
       //the Navier--Stokes equations (which arise through the buoyancy
       //body-force term)
       unsigned n_master = 1;
       double hang_weight = 1.0;
       
       //Local bool (is the node hanging)
       bool is_node_hanging = this->node_pt(l)->is_hanging();

       //If the node is hanging, get the number of master nodes
       if(is_node_hanging)
        {
         hang_info_pt = this->node_pt(l)->hanging_pt();
         n_master = hang_info_pt->nmaster();
        }
       //Otherwise there is just one master node, the node itself
       else 
        {
         n_master = 1;
        }
       
       //Loop over the master nodes
       for(unsigned m=0;m<n_master;m++)
        {
         //If the node is hanging get weight from master node
         if(is_node_hanging)
          {
           //Get the hang weight from the master node
           hang_weight = hang_info_pt->master_weight(m);
          }
         else
          {
           // Node contributes with full weight
           hang_weight = 1.0;
          }
         
         
         //Loop over the velocity components in the Navier--Stokes equtions
         for(unsigned i=0;i<DIM;i++)
          {
           //Get the equation number
           if(is_node_hanging)
            {
             //Get the equation number from the master node
             local_eqn = this->local_hang_eqn(hang_info_pt->master_node_pt(m),
                                              u_nodal_nst[i]);
            }
           else
            {
             // Local equation number
             local_eqn = this->nodal_local_eqn(l,u_nodal_nst[i]);
            }
           
           if(local_eqn >= 0)
            {
             //Loop over the external data
             for(unsigned l2=0;l2<n_external_element_data;l2++)
              { 
               //Find the local equation number corresponding to the global
               //unknown
               local_unknown = 
                this->local_eqn_number(
                 global_eqn_number_of_external_element_data[l2]);
               if(local_unknown >= 0)
                {
                 //Add contribution to jacobian matrix
                 jacobian(local_eqn,local_unknown) 
                  += dbody_dexternal_element_data(i,l2)*testf(l)*hang_weight*W;
                }
              }
            }
          }
        }
      }
    }
  } 


private:

 /// Pointer to a private data member, the Rayleigh number
 double* Ra_pt;

 /// The static default value of the Rayleigh number
 static double Default_Physical_Constant_Value;


};




///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////



//======================ad_bous_class==================================
/// Build an AdvectionDiffusionElement that inherits from 
/// ElementWithExternalElement so that it can "communicate" with the 
/// a NavierStokesElement that provides its wind.
//=====================================================================
template<unsigned DIM>
class RefineableQAdvectionDiffusionBoussinesqElement : 
 public virtual RefineableQAdvectionDiffusionElement<DIM,3>,
 public virtual ElementWithExternalElement
{

public:

 /// \short Constructor: call the underlying constructors
 RefineableQAdvectionDiffusionBoussinesqElement() : 
  RefineableQAdvectionDiffusionElement<DIM,3>(), ElementWithExternalElement()
  { 
   //There is one interaction
   this->set_ninteraction(1);
  }



// hierher why do we have these?

/*  /// \short Output function: */
/*  ///  Output x, y, theta at Nplot^DIM plot points */
/*  // Start of output function */
/*  void output(ostream &outfile, const unsigned &nplot) */
/*   { */
   
/*    //vector of local coordinates */
/*    Vector<double> s(DIM); */
   
/*    // Tecplot header info */
/*    outfile << this->tecplot_zone_string(nplot); */
   
/*    // Loop over plot points */
/*    unsigned num_plot_points=this->nplot_points(nplot); */
/*    for (unsigned iplot=0;iplot<num_plot_points;iplot++) */
/*     { */
/*      // Get local coordinates of plot point */
/*      this->get_s_plot(iplot,nplot,s); */
     
/*      // Output the position of the plot point */
/*      for(unsigned i=0;i<DIM;i++) */
/*       {outfile << this->interpolated_x(s,i) << " ";} */
     
/*      // Output the temperature (the advected variable) at the plot point */
/*      outfile << this->interpolated_u_adv_diff(s) << std::endl; */
/*     } */
/*    outfile << std::endl; */
   
/*    // Write tecplot footer (e.g. FE connectivity lists) */
/*    this->write_tecplot_zone_footer(outfile,nplot); */
   
/*   } //End of output function */


/*  ///  Overload the standard output function with the broken default */
/*  void output(ostream &outfile) {FiniteElement::output(outfile);} */

/*  /// \short C-style output function: Broken default */
/*  void output(FILE* file_pt) */
/*   {FiniteElement::output(file_pt);} */

/*  ///  \short C-style output function: Broken default */
/*  void output(FILE* file_pt, const unsigned &n_plot) */
/*   {FiniteElement::output(file_pt,n_plot);} */


 /// \short Overload the wind function in the advection-diffusion equations.
 /// This provides the coupling from the Navier--Stokes equations to the
 /// advection-diffusion equations because the wind is the fluid velocity,
 /// obtained from the "external" element
 void get_wind_adv_diff(const unsigned& ipt, const Vector<double> &s, 
                        const Vector<double>& x, Vector<double>& wind) const
 {
  // There is only one interaction
  unsigned interaction=0;
  
  // Dynamic cast "external" element to Navier Stokes element
  NavierStokesEquations<DIM>* nst_el_pt=
   dynamic_cast<NavierStokesEquations<DIM>*>
   (external_element_pt(interaction,ipt));
  
  //Wind is given by the velocity in the Navier Stokes element
  nst_el_pt->interpolated_u_nst
   (external_element_local_coord(interaction,ipt),wind);
  
 }  //end of get_wind_adv_diff
 
 
// hierher

/*  /// Just call the fill_in_residuals for AdvDiff */
/*  void fill_in_contribution_to_residuals(Vector<double> &residuals) */
/*   { */
/*    RefineableAdvectionDiffusionEquations<DIM>:: */
/*     fill_in_contribution_to_residuals(residuals); */
/*   } */


 ///\short Compute the element's residual vector and the Jacobian matrix.
 void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                       DenseMatrix<double> &jacobian)
  {
   //Get the contribution from the basic Navier--Stokes element
   RefineableQAdvectionDiffusionElement<DIM,3>::
    fill_in_contribution_to_jacobian(residuals,jacobian);
   
#ifdef USE_FD_FOR_DERIVATIVES_WRT_EXTERNAL_DATA

   //Get the off-diagonal terms by finite differencing
   this->fill_in_jacobian_from_external_interaction_by_fd(residuals,jacobian);

#else

   //Get the off-diagonal terms analytically
   this->fill_in_off_diagonal_block_analytic(residuals,jacobian);

#endif

  }


 /// \short Overload the function that must return all field data involved
 /// in the interaction with the external (Navier Stokes) element. 
 /// Only the velocity dofs in the Navier Stokes element affect the
 /// interaction with the current element. 
 void identify_all_field_data_for_external_interaction(
  Vector<std::set<FiniteElement*> > const &external_elements_pt,
  std::set<std::pair<Data*,unsigned> > &paired_interaction_data);

 /// \short Add the element's contribution to its residuals vector,
 /// jacobian matrix and mass matrix
 void fill_in_contribution_to_jacobian_and_mass_matrix(
  Vector<double> &residuals, DenseMatrix<double> &jacobian, 
  DenseMatrix<double> &mass_matrix)
  {
   //Call the standard (Broken) function
   //which will prevent these elements from being used
   //in eigenproblems until replaced.
   FiniteElement::fill_in_contribution_to_jacobian_and_mass_matrix(
     residuals,jacobian,mass_matrix);
  }

 
 /// \short Fill in the derivatives of the wind with respect to the
 /// external unknowns
 void get_dwind_adv_diff_dexternal_element_data(
  const unsigned& ipt, const unsigned &i,
  Vector<double> &result, Vector<unsigned> &global_eqn_number)
 {
  // The interaction index is 0 in this case
  unsigned interaction=0;
  
  // Dynamic cast "other" element to correct type
  RefineableQCrouzeixRaviartBoussinesqElement<DIM>* source_el_pt=
   dynamic_cast<RefineableQCrouzeixRaviartBoussinesqElement<DIM>*>
   (external_element_pt(interaction,ipt));
  
  // Get the external element's derivatives of the velocity with respect
  // to the data. The wind is just the Navier--Stokes velocity, so this
  // is all that's required
  source_el_pt->dinterpolated_u_nst_ddata(
   external_element_local_coord(interaction,ipt),i,result,
   global_eqn_number);
 } 

 
 /// \short Compute the contribution of the external
 /// degrees of freedom (velocities) on the advection-diffusion equations
 void fill_in_off_diagonal_block_analytic(Vector<double> &residuals,
                                          DenseMatrix<double> &jacobian)
  {
   //Local storage for the index in the nodes at which the temperature 
   //is stored
   const unsigned u_nodal_adv_diff = this->u_index_adv_diff();

   //Find out how many nodes there are
   const unsigned n_node = this->nnode();
   
   //Set up memory for the shape and test functions and their derivatives
   Shape psi(n_node), test(n_node);
   DShape dpsidx(n_node,DIM), dtestdx(n_node,DIM);
   
   //Number of integration points
   const unsigned n_intpt = this->integral_pt()->nweight();
   
   //Integers to store the local equations and unknowns
   int local_eqn=0, local_unknown=0;
   
   // Local storage for pointers to hang_info objects
   HangInfo *hang_info_pt=0;   

   //Get the peclet number
   const double peclet = this->pe();

   //Loop over the integration points
   for(unsigned ipt=0;ipt<n_intpt;ipt++)
    {
     //Get the integral weight
     double w = this->integral_pt()->weight(ipt);
     
     //Call the derivatives of the shape and test functions
     double J = 
      this->dshape_and_dtest_eulerian_at_knot_adv_diff(ipt,psi,dpsidx,
                                                       test,dtestdx);
     
     //Premultiply the weights and the Jacobian
     double W = w*J;
     
     //Calculate local values of the derivatives of the solution
     Vector<double> interpolated_dudx(DIM,0.0);
     // Loop over nodes
     for(unsigned l=0;l<n_node;l++) 
      {
       // Loop over directions
       for(unsigned j=0;j<DIM;j++)
        {
         interpolated_dudx[j] += 
          this->nodal_value(l,u_nodal_adv_diff)*dpsidx(l,j);
        }
      }
     
     //Get the derivatives of the wind wrt the unknowns
     //of the external element
     Vector<double> dwind_dexternal_element_data;
     //Vector of global equation number corresponding to the external
     //element's data
     Vector<unsigned> global_eqn_number_of_external_element_data;

     //Loop over the wind directions
     for(unsigned i2=0;i2<DIM;i2++)
      {
       //Get the appropriate derivatives
       this->get_dwind_adv_diff_dexternal_element_data(
        ipt,i2,dwind_dexternal_element_data,
        global_eqn_number_of_external_element_data);
       

       //Find out how many external data there are
       const unsigned n_external_element_data = 
        global_eqn_number_of_external_element_data.size();
       
       //Loop over the test functions
       for(unsigned l=0;l<n_node;l++)
        {
         //Assemble the contributions of the velocities to 
         //the advection-diffusion equations
         unsigned n_master = 1;
         double hang_weight = 1.0;
         
         //Local bool (is the node hanging)
         bool is_node_hanging = this->node_pt(l)->is_hanging();
         
         //If the node is hanging, get the number of master nodes
         if(is_node_hanging)
          {
           hang_info_pt = this->node_pt(l)->hanging_pt();
           n_master = hang_info_pt->nmaster();
          }
         //Otherwise there is just one master node, the node itself
         else 
          {
           n_master = 1;
          }
         
         //Loop over the master nodes
         for(unsigned m=0;m<n_master;m++)
          {
           //If the node is hanging get weight from master node
           if(is_node_hanging)
            {
             //Get the hang weight from the master node
             hang_weight = hang_info_pt->master_weight(m);
            }
           else
            {
             // Node contributes with full weight
             hang_weight = 1.0;
            }
           
           //Get the equation number
           if(is_node_hanging)
            {
             //Get the equation number from the master node
             local_eqn = this->local_hang_eqn(hang_info_pt->master_node_pt(m),
                                              u_nodal_adv_diff);
            }
           else
            {
             // Local equation number
             local_eqn = this->nodal_local_eqn(l,u_nodal_adv_diff);
            }
           
           if(local_eqn >= 0)
            {
             //Loop over the external data
             for(unsigned l2=0;l2<n_external_element_data;l2++)
              { 
               //Find the local equation number corresponding to the global
               //unknown
               local_unknown = 
                this->local_eqn_number(
                 global_eqn_number_of_external_element_data[l2]);
               if(local_unknown >= 0)
                {
                 //Add contribution to jacobian matrix
                 jacobian(local_eqn,local_unknown) 
                  -= peclet*dwind_dexternal_element_data[l2]*
                  interpolated_dudx[i2]*test(l)*hang_weight*W;
                }
              }
            }
          }
        }
      }
    }
  }

};



///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////


//================optimised_identification_of_field_data==================
/// Overload the function that must return all field data involved
/// in the interaction with the external (Navier Stokes) element. 
/// Only the velocity dofs in the Navier Stokes element affect the
/// interaction with the current element. 
//=======================================================================
template<unsigned DIM>
void RefineableQAdvectionDiffusionBoussinesqElement<DIM>::
identify_all_field_data_for_external_interaction(
 Vector<std::set<FiniteElement*> > const &external_elements_pt,
 std::set<std::pair<Data*,unsigned> > &paired_interaction_data)
{
 //There's only one interaction
 const unsigned interaction = 0;
 
 // Loop over each Navier Stokes element in the set of external elements that
 // affect the current element
 for(std::set<FiniteElement*>::iterator it=
      external_elements_pt[interaction].begin();
     it != external_elements_pt[interaction].end(); it++)
  {
   
   //Cast the external element to a fluid element
   NavierStokesEquations<DIM>* external_fluid_el_pt =
    dynamic_cast<NavierStokesEquations<DIM>*>(*it);
   
   // Loop over the nodes
   unsigned nnod=external_fluid_el_pt->nnode();
   for (unsigned j=0;j<nnod;j++)
    {
     // Pointer to node (in its incarnation as Data)
     Data* veloc_data_pt=external_fluid_el_pt->node_pt(j);
     
     // Get all velocity dofs
     for (unsigned i=0;i<DIM;i++)
      {
       // Which value corresponds to the i-th velocity?
       unsigned val=external_fluid_el_pt->u_index_nst(i);
       
       // Turn pointer to Data and index of value into pair
       // and add to the set
       paired_interaction_data.insert(std::make_pair(veloc_data_pt,val));
      }
    }
  }
} // done


/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////


//==========start_of_get_dbody_force=========== ===========================
/// Fill in the derivatives of the body force with respect to the external
/// unknowns in the Navier--Stokes equations
//=========================================================================
template<unsigned DIM>
void RefineableQCrouzeixRaviartBoussinesqElement<DIM>::
get_dbody_force_nst_dexternal_element_data(const unsigned &ipt,
                                           DenseMatrix<double> &result,
                                           Vector<unsigned> &global_eqn_number)
{
 // The interaction index is 0 in this case
 unsigned interaction=0;
 
 // Dynamic cast "other" element to correct type
 RefineableQAdvectionDiffusionBoussinesqElement<DIM>* source_el_pt=
  dynamic_cast<RefineableQAdvectionDiffusionBoussinesqElement<DIM>*>
  (external_element_pt(interaction,ipt));
 
 // Get vector that indicates the direction of gravity from
 // the Navier-Stokes equations
 Vector<double> gravity(NavierStokesEquations<DIM>::g());
 
 // Get the external element's derivatives
 Vector<double> du_adv_diff_ddata;
 source_el_pt->dinterpolated_u_adv_diff_ddata(
  external_element_local_coord(interaction,ipt),du_adv_diff_ddata,
  global_eqn_number);
  
 //Find the number of external data
 unsigned n_external_element_data = du_adv_diff_ddata.size();
 //Set the size of the matrix to be returned
 result.resize(DIM,n_external_element_data);

 // Temperature-dependent body force:
 for (unsigned i=0;i<DIM;i++)
  {
   //Loop over the external data
   for(unsigned n=0;n<n_external_element_data;n++)
    {
     result(i,n) = -gravity[i]*du_adv_diff_ddata[n]*ra();
    }
  }

} // end_of_get_dbody_force



/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////




//=========================================================================
/// Set the default physical value for the Rayleigh number to be zero 
/// in 2D and 3D
//=========================================================================
template<>
double RefineableQCrouzeixRaviartBoussinesqElement<2>::Default_Physical_Constant_Value=0.0;

template<>
double RefineableQCrouzeixRaviartBoussinesqElement<3>::Default_Physical_Constant_Value=0.0;

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

