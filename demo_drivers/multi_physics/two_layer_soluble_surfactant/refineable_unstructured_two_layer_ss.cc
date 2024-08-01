//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2023 Matthias Heil and Andrew Hazel
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
//Driver for a multi-physics problem that couples the free-surface
//Navier--Stokes equations to the surface transport equations for interfacial
//surfactant motion as well as the advection-diffusion equations in the bulk
//describing the transport of soluble surfactant and miceles. 

//At the moment the assumption that bulk surfactant is only present in
//the lower layer is handled by pinning all values in the upper layer
//and switching off the contributions to the transport equations in the upper
//elements to prevent the upper layer making any contribution at all.
//If the problem is to be generalised
//to allow surfactant transport in both layers then an additional field
//will have to be added with a new variable that is pinned (and hijacked)
//in the lower field. Wherever possible I've tried to stick to the notation
//used in Kalogriou and Blyth

//N.B. This all requires more careful validation
//Simple checks that pass are that mass of surfactant is conserved, for
//a suitably small timestep; and that the steady solution matches the
//analytic expression.

#include<iostream>

//Oomph-lib headers, 
//We require the generic header
#include "generic.h"
//Our custom coupling of advection-diffusion and Navier--Stokes
#include "double_buoyant_navier_stokes_elements.h"
//The fluid interface elements
#include "fluid_interface.h"
//The surfactant transport equations
#include "soluble_surfactant_transport_equations.h"

#include "constitutive.h"
#include "solid.h"

// The mesh is our standard rectangular quadmesh
#include "meshes/triangle_mesh.h"

// Use the oomph and std namespaces 
using namespace oomph;
using namespace std;

namespace oomph
{

namespace Control_Parameters
{
 bool Periodic_BCs = true;

 bool Pin_Micelle = true;
}

  
//======start_of_namespace============================================
/// Namespace for the physical parameters in the problem
//====================================================================
namespace Global_Physical_Variables
{
  /// Geometry
  //----------------
 double L = 1.0; //2.4;

  /// Fluid property Ratios
  //----------------------------
  
  //Density ratio:
  /// Ratio of density in upper fluid to density in lower
  /// fluid. Reynolds number etc. is based on density in lower fluid.
  double R = 1.0;

  //Vecosity ratio;
  /// Ratio of viscosity in upper fluid to viscosity in lower
  /// fluid. Reynolds number etc. is based on viscosity in lower fluid.
 double M = 17.0; //0.5

  //Dimensionless position of interface (relative to a total domain height of 1)
 double H0 = 0.2; //1.0/11.0;

  //Dimensionless size of perturbation to the interface
  double Ha = H0/5.0;
  
  /// Hydrodynamic Parameters
  //----------------------------
  
  /// Reynolds number
  double Re = 0.0; 

  /// We do not need the Weber number,
  /// because we have non-dimensionalised pressure
  /// on the viscous scale and so multiplying the normal stress balanced by the
  /// Reynolds number gives a term in the Capillary number only (Ca Re = We).
  
  /// Capillary number (of which the results are independent
  /// for a pinned surface)
  double Ca = 1.0; //0.1;
  
  /// In our non-dimensionalisation, we have a
  /// Reynolds number/ Froude number squared in
  /// (a material parameter that doesn't involve the velocity scale).
  /// 
  double ReInvFr = 0.0; // (Fr = 1)

  /// Surfactant Parameters
  //--------------------------
  
  /// Marangoni number
 double Ma = 8.0; //0.1;

  /// Surface Elasticity number (Capillary number x Marangoni number)
  double Beta_s =Ca*Ma;

  /// Surface Peclet number
 double Pe_s = 1.0e8; // 10.0; 
  
  /// Bulk Peclet number
  double Pe_b = 100.0;
  
  /// Micelle Pelect number
  double Pe_m = 100.0;

  /// Solubility Parameters
  //-------------------------
 
  /// Biot number
 double Biot = 1.0; //0.01; 
  
  /// The ratio of adsorption-desorption times
 double K_b = 1.0; //0.5; 

 // ratio of equilibrium concentrations
 double Beta_b = 1.0;

  // Reaction rate between bulk and micelle 
 double K_m = 0.0;

 /// Power of the concentration in bulk -> micelle flux expression
 double N = 10.0;
  
 /// The imposed pressure gradient
 double Delta_P = 1.0; 
  
 /// Timescales for transport equations (identically one from our
 /// non-dimensionalisation)
 Vector<double> Tau(2,1.0);
 
 /// Diffusivity  (should be 1/Pe_b, 1/Pe_m), which
 /// will be set in the main code
 Vector<double> D(2,1.0);
 
 /// Gravity vector, will be set in the main code
 Vector<double> Direction_of_gravity(2);

 /// Pseudo-solid Poisson ratio
 double Nu = 0.1;
 
 /// This next set of functions is only used if we do NOT have
 /// periodic conditions

 
 /// Function that prescribes the hydrostatic pressure field at the outlet
 /// Let's fix things so that the pressure at the top of the channel is zero.
 void hydrostatic_pressure_outlet_upper(const double& time,
                                        const Vector<double> &x, 
                                        const Vector<double> &n, 
                                        Vector<double> &traction)
 {
  traction[0] = ReInvFr*R*Direction_of_gravity[1]*(1.0 - x[1]);
  traction[1] = 0.0;
 }

 /// Function that prescribes hydrostatic pressure field at the inlet
 void hydrostatic_pressure_inlet_upper(const double& time, const Vector<double> &x, 
				       const Vector<double> &n,
				       Vector<double> &traction)
 {
   traction[0] = Delta_P + -ReInvFr*R*Direction_of_gravity[1]*(1.0 - x[1]);
   traction[1] = 0.0;
 }

 
 /// Function that prescribes the hydrostatic pressure field at the outlet
 /// Must match pressure in lower fluid --- This may be tricky if the
 /// interface is not pinned
 /// (i.e. we'll need to read out the interfacial position
 /// on the boundary). For now assume it's at H0.
  void hydrostatic_pressure_outlet_lower(const double& time,
                                         const Vector<double> &x, 
                                  const Vector<double> &n, 
                                  Vector<double> &traction)
 {
  traction[0] = ReInvFr*Direction_of_gravity[1]*(R*(1.0 - H0) + H0 - x[1]);
  traction[1] = 0.0;
 }

 /// Function that prescribes hydrostatic pressure field at the inlet
 void hydrostatic_pressure_inlet_lower(const double& time, const Vector<double> &x, 
				       const Vector<double> &n,
				       Vector<double> &traction)
 {
  traction[0] = Delta_P +
   -ReInvFr*Direction_of_gravity[1]*(R*(1.0 - H0) + H0 - x[1]);
   traction[1] = 0.0;
 }

  //end of traction functions

 //Set specificied angle if required
  double Inlet_Angle = 2.0*atan(1.0);

  
 /// Direction of the wall normal vector (at the inlet)
 Vector<double> Wall_normal;

 /// Function that specifies the wall unit normal at the inlet
 void wall_unit_normal_inlet_fct(const Vector<double> &x, 
                                 Vector<double> &normal)
 {
  normal=Wall_normal;
 }

 /// Function that specified the wall unit normal at the outlet
 void wall_unit_normal_outlet_fct(const Vector<double> &x, 
                                 Vector<double> &normal)
 {
  //Set the normal
  normal = Wall_normal;
  //and flip the sign
  unsigned n_dim = normal.size();
  for(unsigned i=0;i<n_dim;++i) {normal[i] *= -1.0;}
 }

  
} // end_of_namespace


} //end of oomph namespace



//A Comparison operator for the boundary nodes
class CompareNodeCoordinates
{
public:
///The actual comparison operator
 int operator() (Node* const &node1_pt,
                 Node* const &node2_pt)
  {
   unsigned n_dim = node1_pt->ndim();
   if(n_dim != node2_pt->ndim())
    {
     throw OomphLibError("Can't compare two nodes of different dimension",
                         "CompareNodeCoordinates::operator()",
                         OOMPH_EXCEPTION_LOCATION);
    }

   //Make sure to handle the finite precision problems
   //associated with coordinates!
   unsigned i=0;
   //const double tol = 3.0e-7;
   const double tol = 4.0e-7;
   while((i < n_dim) && (std::abs(node1_pt->x(i) - node2_pt->x(i)) < tol)) 
    {++i;}
   //If we've got to the end, they are the same so return false
   if(i==n_dim) {return false;}
   
   //Otherwise
   return node1_pt->x(i) < node2_pt->x(i);
  }
};


//==start_of_specific_element_class=============================
/// Element class used to insist that the vertical positions of
/// the periodic nodes coincide.
/// These are essentially point elements that are created
/// from existing periodic nodes.
//=============================================================
class DependentPositionPointElement : public virtual SolidPointElement,
                                      public virtual SolidFaceElement
{
  //Dependent node
  Node* Dependent_node_pt;

  //Face ID
  unsigned Id;

 //Boolean to indicate whether this is on the interface
 bool Is_on_interface;
 
   /// Fill in the residuals for the volume constraint
 void fill_in_generic_contribution_to_residuals_match_position(
  Vector<double> &residuals,
  DenseMatrix<double> &jacobian,
  unsigned flag)
  {
    //Cache the node
    Node* const nod_pt = this->node_pt(0);

    //Storage for the local equation number
    int local_eqn = 0;
    
    //Read out the index at which the Lagrange multiplier is stored
    const unsigned lm_index = dynamic_cast<BoundaryNodeBase*>(nod_pt)
     ->index_of_first_value_assigned_by_face_element(Id);
    
    //Read out the Lagrange multiplier
    const double lambda = nod_pt->value(lm_index);
    
    //Get the local equation for the dependent node,
    //This is the node that the Lagrange multiplier should be affecting
    local_eqn = this->external_local_eqn(0,1);
    if(local_eqn >=0)
     {
      residuals[local_eqn] = -lambda;
      
      if(flag)
       {
        int local_unknown = this->nodal_local_eqn(0,lm_index);
        if(local_unknown >= 0)
         {
          jacobian(local_eqn,local_unknown) = -1.0;
         }
       }
     }
    //End of Lagrange multiplier eqn
     
    //HARD-CODED DIRECTION THAT HAS TO MATCH IS Y (1)
    //If the master is pinned then we don't bother to do anything
    //Loop over the number of coordinate directions
    //Let's get the Lagrange multiplier
    local_eqn = this->nodal_local_eqn(0,lm_index);
    
    //Residuals for matching vertical positions between the nodes

    //If it's not a boundary condition (i.e. not pinned)
    //then add the lagrange multiplier
    if(local_eqn >=0)
     {
      residuals[local_eqn] = nod_pt->x(1) - Dependent_node_pt->x(1); 
      //Jacobian terms
      if(flag)
       {
        int local_unknown = this->position_local_eqn(0,0,1);
        //The entry in the jacobian is one
        if(local_unknown >= 0)
         {
          jacobian(local_eqn,local_unknown) = 1.0;
         }
	
        //The off diagonal term is minus one, but
        //we need to work out the corresponding external equation numbers
        //There is only one external datum so the first entry is the required number
        local_unknown = this->external_local_eqn(0,1);
        if(local_unknown >= 0)
         {
          jacobian(local_eqn,local_unknown) = -1.0;
         }
       } //End of Jacobian calculation
     }
    
  }

  
public:

  //Constructor: Make this from a single Node
 DependentPositionPointElement(Node* const &node_pt, Node *const &dependent_node_pt,
                               const unsigned &id=0, const bool &is_on_interface=false) :
  Dependent_node_pt(dependent_node_pt), Id(id), Is_on_interface(is_on_interface)
  {
    //There is only one node
    this->set_n_node(1);
    //The node must be the master
    this->node_pt(0) = node_pt;
    //THIS IS REALLY IMPORTANT.
    //IF YOU MISS IT THEN LOCAL EQUATION NUMBERING GOES WRONG
    this->set_nodal_dimension(node_pt->ndim());
    //Set the dependent node's variable position as external data
    SolidNode* solid_nod_pt = static_cast<SolidNode*>(Dependent_node_pt);
    this->add_external_data(solid_nod_pt->variable_position_pt());
    //There is one extra value (the Lagrange multiplier for each node)
    Vector<unsigned> n_additional_values(1,1);
    //We also need to add the lagrange multiplier
    this->add_additional_values(n_additional_values,id);
   }

  //Broken copy constructore
  DependentPositionPointElement(const DependentPositionPointElement&)
  {
    BrokenCopy::broken_copy("DependentPositionPointElement");
  }

  //Empty Destructor
  ~DependentPositionPointElement() {}
  
 /// Fill in the residuals for the volume constraint
 void fill_in_contribution_to_residuals( Vector<double> &residuals)
 {
  this->fill_in_generic_contribution_to_residuals_match_position(
   residuals,GeneralisedElement::Dummy_matrix,0);
 }
  
 /// Fill in the residuals and jacobian for the volume constraint
 void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                       DenseMatrix<double> &jacobian)
 {
  //No contribution to jacobian; see comment in that function
  this->fill_in_generic_contribution_to_residuals_match_position(
   residuals,jacobian,1);
 }
    
  
};


//Custom Mesh class to wrap up a few interfaces nicely so that the code looks
//as much like the unrefineable version as possible
template <class ELEMENT>
class RefineableSolidTwoLayerTriangleMesh :
 public virtual RefineableSolidTriangleMesh<ELEMENT>
{
public:
 //Not sure why I have to call all the constructors...
    RefineableSolidTwoLayerTriangleMesh(
      TriangleMeshParameters& triangle_mesh_parameters,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
     : TriangleMesh<ELEMENT>(triangle_mesh_parameters,time_stepper_pt),
     RefineableTriangleMesh<ELEMENT>(triangle_mesh_parameters,time_stepper_pt),
     RefineableSolidTriangleMesh<ELEMENT>(triangle_mesh_parameters,time_stepper_pt)
  {}

    /// Build mesh from specified triangulation and
    /// associated target areas for elements in it.
    RefineableSolidTwoLayerTriangleMesh(
      const Vector<double>& target_area,
      TriangulateIO& triangulate_io,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper,
      const bool& use_attributes = false,
      const bool& allow_automatic_creation_of_vertices_on_boundaries = true,
      OomphCommunicator* comm_pt = 0)
      : RefineableSolidTriangleMesh<ELEMENT>(
          target_area,
          triangulate_io,
          time_stepper_pt,
          use_attributes,
          allow_automatic_creation_of_vertices_on_boundaries,
          comm_pt) {}



 //Helper function for the custom update of the nodes on the p-th  PSLG of
 //the outer boundary polygon
 void add_common_nodes_to_boundary_vector(TriangleMeshPolygon*& polygon_pt,
                                          const unsigned &bound,
                                          const std::set<double> &vertical_coordinates,
                                          Vector<Vector<double> > &vector_vertex_node)
  {
   unsigned n_common = vertical_coordinates.size();
   std::cout << n_common << " nodes required to match both boundaries\n";

   //The requirements here are specific for this mesh, so I've not
   //attempted to write this as general code.
   
   //The idea is to loop over the existing nodes and add only the
   //missing nodes, but it's easier to create new vectors
   unsigned n_vertex = vector_vertex_node.size();
   
   //We only need to do something if the number of nodes differs.
   if(n_common!=n_vertex)
    {
     //Store the x value (assumed to be the same for all nodes)
     double x = vector_vertex_node[0][0];
     
     //Are we increasing or decreasing in y
     //If the last node has a greater y coordinate than the first
     //then ascending
     if(vector_vertex_node[n_vertex-1][1] > vector_vertex_node[0][1])
      {
       //Clear the vector
       vector_vertex_node.clear();
       vector_vertex_node.resize(n_common);
       unsigned i=0;
       //Loop over the nodes in the common set and add in
       for(std::set<double>::iterator it=vertical_coordinates.begin();
           it!=vertical_coordinates.end();++it)
        {
         //Create the coordinates
         Vector<double> new_vertex_coord(2);
         new_vertex_coord[0] = x; new_vertex_coord[1] = *it;
         vector_vertex_node[i] = new_vertex_coord;
         ++i;
        }
      }
     //Otherwise do it in reverse
     else
      {
       //Clear the vector
       vector_vertex_node.clear();
       vector_vertex_node.resize(n_common);
       unsigned i=0;
       //Loop over the nodes in the common set and add in
       for(std::set<double>::reverse_iterator it=vertical_coordinates.rbegin();
           it!=vertical_coordinates.rend();++it)
        {
         //Create the coordinates
         Vector<double> new_vertex_coord(2);
         new_vertex_coord[0] = x; new_vertex_coord[1] = *it;
         vector_vertex_node[i] = new_vertex_coord;
         ++i;
        }
      }
     
     
     // Now update the polyline according to the new vertices
     //Get the boundary id
     //Find the value of p corresponding to the boundary id
     const unsigned n_polyline = polygon_pt->npolyline();
     unsigned p=n_polyline+1;
     //Find the polyline in question
     for(unsigned n=0;n<n_polyline;++n)
      {
       if(polygon_pt->polyline_pt(n)->boundary_id() == bound)
        {
         p=n;
         break;
        }
      }

     //If we haven't found it, report error
     if(p==n_polyline+1)
      {
       throw OomphLibError("Boundary ID not found in polygon\n",
                           OOMPH_CURRENT_FUNCTION,
                           OOMPH_EXCEPTION_LOCATION);
      }

     //Create a new polyline that will have the new nodes on it
     TriangleMeshPolyLine* tmp_polyline_pt =
      new TriangleMeshPolyLine(vector_vertex_node, bound);

     // Create a temporal "curve section" version of the recently
     // created polyline
     TriangleMeshCurveSection* tmp_curve_section_pt = tmp_polyline_pt;

     // Tolerance below which the middle point can be deleted (ratio of
     // deflection to element length)
     double unrefinement_tolerance =
      polygon_pt->polyline_pt(p)->unrefinement_tolerance();

     // Tolerance to add points
     double refinement_tolerance =
      polygon_pt->polyline_pt(p)->refinement_tolerance();
     
     // Establish refinement and unrefinement tolerance on the new polyline
     tmp_polyline_pt->set_unrefinement_tolerance(unrefinement_tolerance);
     tmp_polyline_pt->set_refinement_tolerance(refinement_tolerance);

     // Establish the maximum length constraint
     double maximum_length = polygon_pt->polyline_pt(p)->maximum_length();
     tmp_polyline_pt->set_maximum_length(maximum_length);

#ifdef OOMPH_HAS_MPI
     // If the mesh is distributed check that the polyline still has
     // vertices
     if (this->is_mesh_distributed())
      {
       if (n_vertex >= 2)
        {
         // Pass the connection information from the old polyline to the
         // new one
         this->copy_connection_information(polygon_pt->polyline_pt(p),
                                           tmp_curve_section_pt);
        } // if (n_vertex >= 2)
      } // if (this->is_mesh_distributed())
     else
#endif
      {
       // Pass the connection information from the old polyline to the
       // new one
       this->copy_connection_information(polygon_pt->polyline_pt(p),
                                         tmp_curve_section_pt);
      }

     // Now update the polyline according to the new vertices but first
     // check if the object is allowed to delete the representation or
     // if it should be done by other object
     bool delete_it_on_destructor = false;
     
     std::set<TriangleMeshCurveSection*>::iterator it =
      this->Free_curve_section_pt.find(polygon_pt->curve_section_pt(p));

      if (it != this->Free_curve_section_pt.end())
      {
        this->Free_curve_section_pt.erase(it);
        delete polygon_pt->curve_section_pt(p);
        delete_it_on_destructor = true;
      }

      // -------------------------------------------------------
      // Copying the new representation
      polygon_pt->curve_section_pt(p) = tmp_polyline_pt;

      // Update the Boundary - Polyline map
      this->Boundary_curve_section_pt[bound] = polygon_pt->curve_section_pt(p);

      if (delete_it_on_destructor)
      {
        this->Free_curve_section_pt.insert(polygon_pt->curve_section_pt(p));
      }

    } //End of case if numbers of nodes are different
  } 

 //Snap the nodes onto the boundary
 void snap_nodes_onto_boundary(RefineableTriangleMesh<ELEMENT>*& new_mesh_pt,
                               const unsigned& b)
  {
   //Just call the default version
   RefineableSolidTriangleMesh<ELEMENT>::snap_nodes_onto_boundary(new_mesh_pt,b);
  }

 
 //Custom override
 void update_polygon_custom(TriangleMeshPolygon*& polygon_pt)
  {
   //Find the number of polylines
   const unsigned n_polyline = polygon_pt->npolyline();
   TriangleMeshPolyLine* polyline1_pt=0, *polyline3_pt=0;
   
   for(unsigned n=0;n<n_polyline;++n)
    {
     //Find the polyline
     TriangleMeshPolyLine* polyline_pt=polygon_pt->polyline_pt(n);
     unsigned boundary_id = polyline_pt->boundary_id();
     if(boundary_id==1) {polyline1_pt = polyline_pt;}
     if(boundary_id==3) {polyline3_pt = polyline_pt;}
    }

   //If the periodic boundaries do not have the same number
   //of nodes in the same places, then create the common set
   unsigned n_vertex1 = polyline1_pt->nvertex();
   unsigned n_vertex3 = polyline3_pt->nvertex();
   
   //Need to create the vector of vertex nodes
   Vector<Vector<double> > vector_vertex_node1(n_vertex1);
   Vector<Vector<double> > vector_vertex_node3(n_vertex3);
   
   //Let's find a common set of vertical coordinates
   std::set<double> vertical_coordinates;
   for(unsigned n=0;n<n_vertex1;++n)
    {
     //Push it back
     vector_vertex_node1[n] = polyline1_pt->vertex_coordinate(n);
     vertical_coordinates.insert(polyline1_pt->vertex_coordinate(n)[1]);
    }
   for(unsigned n=0;n<n_vertex3;++n)
    {
     vector_vertex_node3[n] = polyline3_pt->vertex_coordinate(n);
     vertical_coordinates.insert(polyline3_pt->vertex_coordinate(n)[1]);
    }

   {
    //Include a check for repeated nodes within a tolerance error
    //The idea is to collect nodes that are within the tolerance and then merge them
    //There should only ever be two nodes within a tolerance error
    std::set<double>::iterator it=vertical_coordinates.begin();
    while(it!=vertical_coordinates.end())
     {
      //Loop over the other elements in the set, starting at the current point
      for(std::set<double>::iterator it2=std::next(it,1);it2!=vertical_coordinates.end();
          ++it2)
       {
        const double y1 = *it;
        const double y2 = *it2;
        double dist = (y1-y2)*(y1-y2);
        dist = sqrt(dist);
        //If the two points are far enough apart then we're done
        if(dist > ToleranceForVertexMismatchInPolygons::Tolerable_error)
         {
          break;
         }
        //Otherwise there should be only one node, erase the second one
        else
         {
          vertical_coordinates.erase(it2);
          break;
         }
       }
      ++it;
     }
   }

   //Boolen to indicate whether we need to adjust the nodes
   bool adjust_nodes = false;
   //If the number of nodes is different we must do something
   if(n_vertex1 != n_vertex3)
    {
     std::cout << "Different numbers of vertices across the periodic boundary:\n";
     std::cout << n_vertex1 << " on boundary 1 and " << n_vertex3 << " on boundary 3\n";
     adjust_nodes = true;
    }
   else
    {
    //If the number of common nodes is bigger than the number on each boundary
    //then they are not in the same location, so we also need to so something
     if(vertical_coordinates.size() != n_vertex1)
      {
       std::cout << "Same number of vertices across the periodic boundary:\n";
       std::cout << "but they are in different locations!\n";
       adjust_nodes = true;
      }
    }

   if(adjust_nodes)
    {
     //Add the common nodes to each boundary
     add_common_nodes_to_boundary_vector(polygon_pt,1,vertical_coordinates,vector_vertex_node1);
     add_common_nodes_to_boundary_vector(polygon_pt,3,vertical_coordinates,vector_vertex_node3);
    }

  }

  
 
 //Destructor 
 virtual ~RefineableSolidTwoLayerTriangleMesh () {}

 //Return the numbers (Note hard-coded region)
 unsigned nlower() {return(this->nregion_element(1));} 

 unsigned nupper() {return(this->nregion_element(0));} 

 FiniteElement* lower_layer_element_pt(const unsigned &e)
 {return this->region_element_pt(1,e);}

 FiniteElement* upper_layer_element_pt(const unsigned &e)
 {return this->region_element_pt(0,e);}
 
}; // End of specific mesh class


/// ///////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////

//====== start_of_problem_class=======================================
/// 2D surfactant transport problem on rectangular domain, discretised 
/// with spine elements. The specific type
/// of element is specified via the template parameter.
//====================================================================
template<class ELEMENT,class INTERFACE_ELEMENT> 
class SurfactantProblem : public Problem
{

public:

 /// Constructor. The boolean indicates whether the free surface
 //should be pinned or not in the first instance
 SurfactantProblem(const bool &pin=true);

 /// Destructor. Empty
 ~SurfactantProblem() {}

 /// Release the free surface so that it can move
 void unpin_surface()
  {
   //Only bother if the surface is pinned
   if(Surface_pinned)
    {
     Surface_pinned = false;
     
     //Unpin the vertical positions of all nodes not on the top and bottom boundaries
     unsigned n_node = Bulk_mesh_pt->nnode();
     for(unsigned n=0;n<n_node;n++)
      {
       SolidNode* nod_pt = static_cast<SolidNode*>(Bulk_mesh_pt->node_pt(n));
       if(!(nod_pt->is_on_boundary(0) || nod_pt->is_on_boundary(2)))
        {
         nod_pt->unpin_position(1);
        }
      }

     // Loop over the elements to re-enable ALE
     unsigned n_element = Bulk_mesh_pt->nelement();
     for(unsigned i=0;i<n_element;i++)
      {
       // Upcast from GeneralsedElement to the present element
       ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(i));
       el_pt->enable_ALE();
      }

     //Unpin the interfacial surfactant concentrations and Lagrange multipliers
     unsigned n_interface = Surface_mesh_pt->nelement();
     for(unsigned i=0;i<n_interface;i++)
      {
       FiniteElement *el_pt = Surface_mesh_pt->finite_element_pt(i);
       
       //Need to unpin the values of surfactant concentration and Lagrange multiplier
       unsigned n_el_node = el_pt->nnode();
       for(unsigned n=0;n<n_el_node;++n)
        {
	  Node* nod_pt = el_pt->node_pt(n);
	  nod_pt->unpin(4);
	  nod_pt->unpin(5);
          
          //Also unpin the periodic LM
          if((Control_Parameters::Periodic_BCs) &&
             (nod_pt->is_on_boundary(1)))
           {
            unsigned n_value = dynamic_cast<BoundaryNodeBase*>(nod_pt)
             ->nvalue_assigned_by_face_element(Periodic_index);
            //Check that the Lagrange multiplier has been allocated
            if(n_value > 0)
             {
              unsigned periodic_lm_index =
               dynamic_cast<BoundaryNodeBase*>(nod_pt)
               ->index_of_first_value_assigned_by_face_element(Periodic_index);

               //Pin this lagrange multiplier
               nod_pt->unpin(periodic_lm_index);
             }
           }
        }
      }
     
     //Now unpin the bulk concentrations in the lower region
     const unsigned n_lower = Bulk_mesh_pt->nlower();
      // Loop over bulk elements in lower fluid
      for(unsigned e=0;e<n_lower;e++)
	{
	  // Upcast from GeneralisedElement to the present element
	  ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->
						  lower_layer_element_pt(e));
	  unsigned n_node = el_pt->nnode();
	  for(unsigned n=0;n<n_node;++n)
	    {
	      el_pt->node_pt(n)->unpin(2); //Bulk Concentration
              if(Control_Parameters::Pin_Micelle)
               {
                el_pt->node_pt(n)->pin(3); //Pin micelle concentration
               }
              else
               {
                el_pt->node_pt(n)->unpin(3); //Micelle Concentration
               }
	    }
	   }

      //We also need to unpin the Lagrange multipliers for periodicity
      //We must not clash with the lagrange multiplier used to drive
      //the interface so we pin that at the left-hand-side
      if(Control_Parameters::Periodic_BCs)
      {
	unsigned b=3;
	unsigned n_boundary_node = this->Bulk_mesh_pt->nboundary_node(b);
	for(unsigned n=0;n<n_boundary_node;++n)
	  {
	    Node* nod_pt = this->Bulk_mesh_pt->boundary_node_pt(b,n);
	    //If we are not on the top or bottom boundary unpin
	    if((!(nod_pt->is_on_boundary(0))) &&
	       (!(nod_pt->is_on_boundary(2))) && 
               (!(nod_pt->is_on_boundary(4))))
	      {
		unsigned lm_index = dynamic_cast<BoundaryNodeBase*>(nod_pt)
		  ->index_of_first_value_assigned_by_face_element(
								  Periodic_index);
		nod_pt->unpin(lm_index);
		
	      }
	    
	  }
      }


      
      
     //Reassign the equation number
     std::cout << "Surface unpinned to give " 
               << assign_eqn_numbers() << " equation numbers\n";
    }
  }

 /// Setup periodic boundaries
 void setup_periodic_boundaries()
  {
    //Now we need to make the boundaries periodic. Boundary 1 is made periodic from Boundary 3
    //Let's load in the boundary nodes
    Vector<Node*> boundary_nodes1, boundary_nodes3;
    unsigned n_boundary_node1 = Bulk_mesh_pt->nboundary_node(1);
    unsigned n_boundary_node3 = Bulk_mesh_pt->nboundary_node(3);

    if(n_boundary_node1 != n_boundary_node3)
     {
      throw OomphLibError("Different numbers of nodes on putative periodic boundary.",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
     }
    
    
    for(unsigned n=0;n<n_boundary_node1;++n)
     {boundary_nodes1.push_back(Bulk_mesh_pt->boundary_node_pt(1,n));}
    for(unsigned n=0;n<n_boundary_node3;++n)
     {boundary_nodes3.push_back(Bulk_mesh_pt->boundary_node_pt(3,n));}
    
    //Sort them
    std::sort(boundary_nodes1.begin(),boundary_nodes1.end(),CompareNodeCoordinates());
    std::sort(boundary_nodes3.begin(),boundary_nodes3.end(),CompareNodeCoordinates());
    

    //Check that the nodes are aligned
    /*for(unsigned n=0;n<n_boundary_node1;++n)
     {
      std::cout << n << " "  << boundary_nodes1[n]->x(1) << " " << boundary_nodes3[n]->x(1) << "\n";
     }

    //Let's output the nodal positions
    for(unsigned n=0;n<n_boundary_node1;++n)
     {          
      if(boundary_nodes1[n]->is_on_boundary(4))
       {std::cout << "1 " << n << " " << boundary_nodes1[n]->x(1) << "\n";}
      if(boundary_nodes3[n]->is_on_boundary(4))
       {std::cout << "3 " << n << " " << boundary_nodes3[n]->x(1) << "\n";}
       }*/

    //Now we can loop over the nodes and make them periodic
    for(unsigned n=0;n<n_boundary_node1;++n)
     {
      boundary_nodes1[n]->make_periodic(boundary_nodes3[n]);
     }
  }
 

 /// Update the problem specs before solve (empty)
  void actions_before_newton_solve() {}

 /// Update the problem after solve (empty)
  void actions_after_newton_solve(){}

 /// Remember to update the nodes if the surface is not pinned
 void actions_before_newton_convergence_check() {}

 /// Actions before the timestep (update the the time-dependent 
 /// boundary conditions)
 void actions_before_implicit_timestep() 
  {
    set_boundary_conditions(time_pt()->time());
    Bulk_mesh_pt->set_lagrangian_nodal_coordinates();
    
   //Also set lagrange multipliers to zero
    unsigned n_boundary_node = Bulk_mesh_pt->nboundary_node(4);
    for(unsigned n=0;n<n_boundary_node;++n)
      {
       //Read out the node
       Node* node_pt = Bulk_mesh_pt->boundary_node_pt(4,n);
       unsigned index = dynamic_cast<BoundaryNodeBase*>(node_pt)
        ->index_of_first_value_assigned_by_face_element(true,false,0);
       //The Lagrange multiplier is one greater than the index assigned by the face element
       //The first index is the interfacial surfactant concentration
       Bulk_mesh_pt->boundary_node_pt(4,n)->set_value(index+1,0.0);
      }

    //Lagrange multipliers associated with the periodicity constraint
    if(Control_Parameters::Periodic_BCs)
     {
      n_boundary_node = this->Bulk_mesh_pt->nboundary_node(3);
      for(unsigned n=0;n<n_boundary_node;++n)
       {
	Node* nod_pt = Bulk_mesh_pt->boundary_node_pt(3,n);
	if(!(nod_pt->is_on_boundary(0)) &&
	   !(nod_pt->is_on_boundary(2)) &&
           !(nod_pt->is_on_boundary(4)))
         {
          unsigned lm_index = dynamic_cast<BoundaryNodeBase*>(nod_pt)
           ->index_of_first_value_assigned_by_face_element(
            Periodic_index);
          nod_pt->set_value(lm_index,0.0);
         }
       }
     }
	 
  }
  
 
 /// Strip off the interface elements before adapting the bulk mesh
 void actions_before_adapt()
   {
     //Unfix the pressure to avoid double fixing
     this->unfix_pressure(0,0);
     
     if(Control_Parameters::Periodic_BCs)
       {
	 //Delete  the dependent elements and wipe the mesh
	 this->delete_dependent_position_elements();
       }
     
     //Backup the surface mesh
     Backed_up_surface_mesh_pt =
       new BackupMeshForProjection<TElement<1,3> >(Surface_mesh_pt,4,0); 
     
     // Delete the interface elements and wipe the surface mesh
     delete_interface_elements();
    
     // Rebuild the Problem's global mesh from its various sub-meshes
     rebuild_global_mesh();
   }

 /// Rebuild the mesh of interface elements after adapting the bulk mesh
 void actions_after_adapt()
  {
    // Determine number of bulk elements in lower/upper fluids
    const unsigned n_lower = Bulk_mesh_pt->nlower();
    const unsigned n_upper = Bulk_mesh_pt->nupper();
    
    
    // Loop over bulk elements in lower fluid
    for(unsigned e=0;e<n_lower;e++)
     {
      // Upcast from GeneralisedElement to the present element
      ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->
                                              lower_layer_element_pt(e));
      
      // Set the diffusivities number
      el_pt->diff_pt() = &Global_Physical_Variables::D;
      
      // Set the timescales
      el_pt->tau_pt() =&Global_Physical_Variables::Tau;
      
      // Set the Reynolds number
      el_pt->re_pt() = &Global_Physical_Variables::Re;
      
      // Set the Womersley number
      el_pt->re_st_pt() = &Global_Physical_Variables::Re;
      
      // Set the product of the Reynolds number and the inverse of the
      // Froude number
      el_pt->re_invfr_pt() = &Global_Physical_Variables::ReInvFr;
      
      // Set the km parameter
      el_pt->km_pt() = &Global_Physical_Variables::K_m;
      
      // Set the N parameter
      el_pt->n_pt() = &Global_Physical_Variables::N;
      
      
      // Set the direction of gravity
      el_pt->g_pt() = &Global_Physical_Variables::Direction_of_gravity;
      
      // Set the constitutive law
      el_pt->constitutive_law_pt() = this->Constitutive_law_pt;
      
     } // End of loop over bulk elements in lower fluid
    
    
    // Loop over bulk elements in upper fluid 
    for(unsigned e=0;e<n_upper;e++)
     {
      // Upcast from GeneralisedElement to the present element
      ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->
                                              upper_layer_element_pt(e));

      
      // Set the diffusivities number
      el_pt->diff_pt() = &Global_Physical_Variables::D;
      
      // Set the timescales
      el_pt->tau_pt() =&Global_Physical_Variables::Tau;
      
      // Set the Reynolds number
      el_pt->re_pt() = &Global_Physical_Variables::Re;
      
      // Set the Womersley number
      el_pt->re_st_pt() = &Global_Physical_Variables::Re;
      
      // Set the product of the Reynolds number and the inverse of the
      // Froude number
      el_pt->re_invfr_pt() = &Global_Physical_Variables::ReInvFr;
      
      // Set the viscosity ratio
      el_pt->viscosity_ratio_pt() = &Global_Physical_Variables::M;
      
      // Set the density ratio
      el_pt->density_ratio_pt() = &Global_Physical_Variables::R;
      
      // Set the km parameter
      el_pt->km_pt() = &Global_Physical_Variables::K_m;
      
      // Set the N parameter
      el_pt->n_pt() = &Global_Physical_Variables::N;
      
      // Set the direction of gravity
      el_pt->g_pt() = &Global_Physical_Variables::Direction_of_gravity;
      
      // Set the constitutive law
      el_pt->constitutive_law_pt() = Constitutive_law_pt;

      // Do not solve the transport equations (no surfactant here)
      el_pt->disable_advection_diffusion_reaction_equations();
      
     } // End of loop over bulk elements in upper fluid


    // Create the interface elements
    //Must be done *after* the physical parameters have been set, because we
    //use the viscosity ratio point to determine upper or lower elements!
    this->create_interface_elements();
    
    //Now we need to make sure that we set up the correct insoluble surfactant
    //concentration
    //We do this by interpolating from the previous mesh
    
    // Now project from backup of original contact mesh to new one
    Backed_up_surface_mesh_pt->project_onto_new_mesh(
						     Surface_mesh_pt);   
    //Now delete the backed up mesh
    delete Backed_up_surface_mesh_pt;
    Backed_up_surface_mesh_pt=0;
    
    if(Control_Parameters::Periodic_BCs)
      {
       std::cout << "Setting up periodic boundaries\n";
       //Need to reset periodic boundary conditions on the boundaries...
       this->setup_periodic_boundaries();
       
	//Create the dependent position elements
	this->create_dependent_position_elements();
      }
    
    // Rebuild the Problem's global mesh from its various sub-meshes
    this->rebuild_global_mesh();
    
    // Pin horizontal displacement of all nodes
    const unsigned n_node = Bulk_mesh_pt->nnode();
    for(unsigned n=0;n<n_node;n++) {Bulk_mesh_pt->node_pt(n)->pin_position(0); }
    

    //Need to create the monitor node again (use a different criterion in general)
    {
     //Choose one of the boundary nodes otherwise gets too tricky
     //Loop over nodes on the surface
     unsigned n_node = Bulk_mesh_pt->nboundary_node(4);
     for(unsigned n=0;n<n_node;++n)
      {
       Node* node_pt = Bulk_mesh_pt->boundary_node_pt(4,n);
       //Choose the right hand edge
       if(node_pt->is_on_boundary(1))
        {
         Monitor_node_pt = node_pt;
         break;
        }
      }
    }
  
 
    
    // Loop over bulk elements in upper fluid 
    for(unsigned e=0;e<n_upper;e++)
      {
     // Upcast from GeneralisedElement to the present element
     ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->
                                             upper_layer_element_pt(e));

     
     //Need to pin the values of concentration and micelle up here
     unsigned n_el_node = el_pt->nnode();
     for(unsigned n=0;n<n_el_node;++n)
       {
	 Node* const nod_pt = el_pt->node_pt(n);
         //For all nodes that are not on the boundary
	 if(nod_pt->is_on_boundary(4) == false)
	   {
	     //Set all the history values to zero as well
	     unsigned nprev_values = nod_pt->time_stepper_pt()->nprev_values();
	     for(unsigned t=0;t<=nprev_values;++t)
	       {
		 nod_pt->set_value(t,2,0.0);
		 nod_pt->set_value(t,3,0.0);
	       }
	     
	     nod_pt->pin(2);
	     nod_pt->pin(3);
	   }
       }
      } // End of loop over bulk elements in upper fluid
    
    
    // Loop over bulk elements in lower fluid and unpin
    for(unsigned e=0;e<n_lower;e++)
      {
	// Upcast from GeneralisedElement to the present element
	ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->
						lower_layer_element_pt(e));
	
	
	//No need to pin the values of concentration and micelle down here
	//But don't unpin if they are hanging
	unsigned n_el_node = el_pt->nnode();
	for(unsigned n=0;n<n_el_node;++n)
	  {
	    Node* const nod_pt = el_pt->node_pt(n);
	    if(!(nod_pt->is_constrained(2))) {nod_pt->unpin(2);}
	    if(!(nod_pt->is_constrained(3)))
             {
              if(Control_Parameters::Pin_Micelle)
               {
                //Pin micelle concentration
                nod_pt->pin(3);
               }
              else
               {
                //unpin micelle concentration
                nod_pt->unpin(3);
               }
             }
	  }
      }
   
   // Set the boundary conditions for this problem: All nodes are
   // free by default -- only need to pin the ones that have Dirichlet 
   // conditions here
   
   //Loop over the boundaries
   unsigned num_bound = Bulk_mesh_pt->nboundary();
   for(unsigned ibound=0;ibound<num_bound;ibound++)
    {
     //If we are on the side-walls, the concentrations
     //satisfy natural boundary conditions, so we only pin the
     //v-velocity for now if not periodic
     if(!Control_Parameters::Periodic_BCs)
      {
       if((ibound==1) || (ibound==3))
        {
         //Loop over the number of nodes on the boundary
         unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
         for(unsigned inod=0;inod<num_nod;inod++)
          {
           //Cast to a solid node, so that we can impose boundary conditions
           SolidNode* nod_pt = static_cast<SolidNode*>(
            Bulk_mesh_pt->boundary_node_pt(ibound,inod));
           nod_pt->pin(1);
           //Also pin the horizontal displacement of the nodes
           nod_pt->pin_position(0);
          }
        }
      }
     
     //If we on the top or bottom wall, velocity is pinned
     //as is the nodal position
     if((ibound==0) || (ibound==2))
      {
       //Loop over the number of nodes on the boundary
       unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
       for(unsigned inod=0;inod<num_nod;inod++)
	{
         SolidNode* nod_pt = static_cast<SolidNode*>(Bulk_mesh_pt->boundary_node_pt(ibound,inod));
         nod_pt->pin(0);
         nod_pt->pin(1);
         nod_pt->pin_position(0);
         nod_pt->pin_position(1);
	}
      }
    }
   
   // Now set the pressure in the first element at 'node' 0 to 0.0
   fix_pressure(0,0,0.0);

   // Reset the boundary conditions
   set_boundary_conditions(this->time_pt()->time());
  }

  /// Create the elements that dependent the position
  void create_dependent_position_elements()
  {
    //We know that boundary 1 is made periodic from boundary 3
    //We also know that the nodes have already been made periodic so we simply need
    //to create the dependent elements to match the vertical positions
    const unsigned n_boundary_node = this->Bulk_mesh_pt->nboundary_node(1);
    for(unsigned n=0;n<n_boundary_node;++n)
      {
	//Cache the boundary node
	Node* nod_pt = this->Bulk_mesh_pt->boundary_node_pt(1,n);
	//Cache the master node
	Node* master_node_pt = nod_pt->copied_node_pt();
	//If there is no master node, then the node has no
	//counterpart on the other side and is therefore
	//not a degree of freedom and does not need to be
	//dependent
	if(master_node_pt!=0)
	  {
	    
	   //Only create dependents for nodes that are not on the interface
           //The interface is driven by the kinematic condition, but we
           //also need to be able to make it periodic.
           //Thus, we need to have the kinematic condition on one "side"
           //and the periodicity condition on the other.
	    if(!(master_node_pt->is_on_boundary(0)) &&
	       !(master_node_pt->is_on_boundary(2)) &&
               !(master_node_pt->is_on_boundary(4)))
	      {
		//Create the dependent node with a position ID
		this->Dependent_position_mesh_pt->
		  add_element_pt(new DependentPositionPointElement(master_node_pt,
                                                                   nod_pt,Periodic_index));
	      }
            //If we are on the interface, create a custom version of the element to
            //avoid setting up a singular problem when combining Lagrange multipliers to enforce
            //both periodicity and the kinematic condition
            if(master_node_pt->is_on_boundary(4))
             {
              const bool is_on_interface = true;
              //Create the dependent node with a position ID
              this->Dependent_position_mesh_pt->
               add_element_pt(new DependentPositionPointElement(master_node_pt,
                                                                nod_pt,Periodic_index,
                                                              is_on_interface));
                                                              }
              
            
          }
        
      } //End of loop over boundary nodes
  }
 



 /// Delete the 1d interface elements
 void delete_dependent_position_elements()
  {
    // Determine number of interface elements
    const unsigned n_dependent_element = Dependent_position_mesh_pt->nelement();
    
    // Loop over the dependent position elements and delete
    for(unsigned e=0;e<n_dependent_element;e++)
      {
	delete Dependent_position_mesh_pt->element_pt(e);
      }
    
    // Wipe the mesh
    Dependent_position_mesh_pt->flush_element_and_node_storage();
  }

  
 /// Create the 1d interface elements
 void create_interface_elements()
  {
   //In the adaptive formulation the only way that we will know which elements
   //are on the lower or upper side is to use the viscosity ratio. This will work
   //even if the density ratio is set to 1 because we are distinguishing based on
   //the pointer.

   
   // Determine number of bulk elements adjacent to interface (boundary 4)
   const unsigned n_element = this->Bulk_mesh_pt->nboundary_element(4);
   
   // Loop over all those elements adjacent to the interface
   for(unsigned e=0;e<n_element;e++)
    {
     // Get pointer to the bulk element that is adjacent to the interface
     ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
      this->Bulk_mesh_pt->boundary_element_pt(4,e));

     // Find index of the face of element e that corresponds to the interface
     const int face_index = this->Bulk_mesh_pt->face_index_at_boundary(4,e);
       
     // We only want to attach interface elements to the bulk elements
     // which are BELOW the interface, and so we filter out those above by
     // referring to the viscosity_ratio_pt

     // If the adjacent element is on the lower side (viscosity ratio not changed)
     if(bulk_elem_pt->viscosity_ratio_pt() !=&Global_Physical_Variables::M)
      {
       // Create the interface element
       INTERFACE_ELEMENT* interface_element_element_pt =
        new INTERFACE_ELEMENT(bulk_elem_pt,face_index);
       
       // Add the interface element to the surface mesh
       this->Surface_mesh_pt->add_element_pt(interface_element_element_pt); 
      }
    }

   
 // --------------------------------------------------------
 // Complete the setup to make the elements fully functional
 // --------------------------------------------------------

 // Determine number of 1D interface elements in mesh
 const unsigned n_interface_element = this->Surface_mesh_pt->nelement();

 // Loop over the interface elements
 for(unsigned e=0;e<n_interface_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   INTERFACE_ELEMENT* el_pt = 
    dynamic_cast<INTERFACE_ELEMENT*>
    (Surface_mesh_pt->element_pt(e));

   // Set the Biot number
   el_pt->bi_pt() = &Global_Physical_Variables::Biot;

   // Set the Marangoni number
   el_pt->ma_pt() =&Global_Physical_Variables::Ma;

   // Set the Ca number
   el_pt->ca_pt() = &Global_Physical_Variables::Ca;

   // Set the surface elasticity number
   el_pt->beta_pt() = &Global_Physical_Variables::Beta_s;

   // Set the surface peclect number
   el_pt->peclet_s_pt() = &Global_Physical_Variables::Pe_s;

   // Set the surface peclect number multiplied by strouhal number
   el_pt->peclet_strouhal_s_pt() = &Global_Physical_Variables::Pe_s;

   // Set the reaction ratio
   el_pt->k_pt() = &Global_Physical_Variables::K_b;


   el_pt->beta_b_pt() = &Global_Physical_Variables::Beta_b;


  } // End of loop over interface elements

  }
     
     


 /// Delete the 1d interface elements
 void delete_interface_elements()
  {
    // Determine number of interface elements
    const unsigned n_interface_element = Surface_mesh_pt->nelement();
    
    // Loop over interface elements and delete
    for(unsigned e=0;e<n_interface_element;e++)
      {
	delete Surface_mesh_pt->element_pt(e);
      }
    
    // Wipe the mesh
    Surface_mesh_pt->flush_element_and_node_storage();
  }
  
 
 /// Fix pressure in element e at pressure dof pdof and set to pvalue
 void fix_pressure(const unsigned &e, const unsigned &pdof, 
                   const double &pvalue)
  {
   //Cast to specific element and fix pressure
   dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e))->
    fix_pressure(pdof,pvalue);
  } // end_of_fix_pressure


 /// UnFix pressure in element e at pressure dof pdof and set to pvalue
 void unfix_pressure(const unsigned &e, const unsigned &pdof)
  {
   //Cast to specific element and fix pressure
   dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e))->
    unfix_pressure(pdof);
  } // end_of_unfix_pressure


//Apply a prescribed deforamtion to the interface
void deform_interface(const double &epsilon,
		      const unsigned &n_periods)
{
 // Determine number of nodes in the "bulk" mesh
 const unsigned n_node = Bulk_mesh_pt->nnode();
 
 // Loop over all nodes in mesh
 for(unsigned n=0;n<n_node;n++)
  {
   // Determine eulerian position of node
   const double current_x_pos = Bulk_mesh_pt->node_pt(n)->x(0);
   const double current_y_pos = Bulk_mesh_pt->node_pt(n)->x(1);

   double y_scale = current_y_pos/Global_Physical_Variables::H0;
   if(current_y_pos > Global_Physical_Variables::H0)
    {
     y_scale = (1.0 - current_y_pos)/(1.0 - Global_Physical_Variables::H0);
    }
   
   // Determine new vertical position of node, *NEED TO THINK*
   const double new_y_pos = current_y_pos
     + y_scale*epsilon
    *(cos(n_periods*MathematicalConstants::Pi*
          (current_x_pos - Global_Physical_Variables::L)/Global_Physical_Variables::L));
   
   // Set new position
   Bulk_mesh_pt->node_pt(n)->x(1) = new_y_pos;
  }
} // End of deform_free_surface

 
  
 /// Doc the solution.
 void doc_solution(std::ofstream &trace);

 /// Set the boundary conditions
 void set_boundary_conditions(const double &time);

 //Return the global error norm to be used in adaptive timestepping
 double global_temporal_error_norm();
  
 /// Overloaded version of the problem's access function to 
 /// the mesh. Recasts the pointer to the base Mesh object to 
 /// the actual mesh type.
 RefineableSolidTwoLayerTriangleMesh<ELEMENT>* Bulk_mesh_pt;

 /// Storage to the mesh of Surface interface elements
 Mesh* Surface_mesh_pt;

 /// Storage for the back up the surface mesh
 BackupMeshForProjection<TElement<1,3> >* Backed_up_surface_mesh_pt;
  
 /// Storage to point elements if needed for non-periodic domains
 Mesh* Point_mesh_pt;

 /// Storage for any traction elements applied to the inlet
 Mesh* Inlet_traction_mesh_pt;

 /// Storage for any traction elements applied to the outlet
 Mesh* Outlet_traction_mesh_pt;

 /// Storage for elements that constraint the vertical positions
 /// of the periodic nodes on the boundaries
 Mesh* Dependent_position_mesh_pt;
 
 /// Pointer to the constitutive law used to determine the mesh deformation
 ConstitutiveLaw* Constitutive_law_pt;


 /// Triangle mesh polygon for outer boundary 
 TriangleMeshPolygon* Outer_boundary_polyline_pt; 

 
  //Calculate the minimum and maximum interface heights
  void interface_min_max(double &min, double &max)
  {
    //Loop over the interface and add each elemental contribution
    const unsigned n_interface  = Surface_mesh_pt->nelement();
    if(n_interface > 0)
      {
	//Set the initial values to the first nodal position
	max = Surface_mesh_pt->finite_element_pt(0)->node_pt(0)->x(1);
	min = max;

	//Loop over all elements and find the max
	for(unsigned i=0;i<n_interface;i++)
	  {
	    // Upcast from GeneralsedElement to the present element
	    INTERFACE_ELEMENT *el_pt =
	      dynamic_cast<INTERFACE_ELEMENT*>(
					       Surface_mesh_pt->element_pt(i));

	    const unsigned n_node = el_pt->nnode();
	    if(n_node != 3)
	      {
		std::cout << "Can't do max and min for elements that are not quadratic\n";
		return;
	      }


	    //Read out the y values from the nodes
	    Vector<double> y(3);
	    for(unsigned n=0;n<3;++n)
	      {
		y[n] = el_pt->node_pt(n)->x(1);
	      }

	    double local_max = y[0];
	    double local_min = y[0];
	    //Find maximum and minimum nodes
	    for(unsigned n=1;n<3;++n)
	      {
		if(y[n] > local_max) {local_max = y[n];}
		if(y[n] < local_min) {local_min = y[n];}
	      }

	    //If we have a linear case then we are done
	    //Check that it's not a degenerate (linear) case
	    if(std::abs(y[0] - 2*y[1] + y[2]) > 1.0e-10)
	      {
		//Calculate extreme value of the local coordinate based on known
		//quadratic basis functions
		//(This should really be inside the element class)
		Vector<double> extreme_s(1,0.5*(y[0] - y[2])/
					 (y[0] - 2.0*y[1] + y[2]));

		//In triangles rescale from zero to 1
		extreme_s[0] = 0.5*(1.0 + extreme_s[0]);

		//Find the extreme height if the local coordinate is within the
		//range of the element
		if((extreme_s[0] >= 0.0) && (extreme_s[0] <= 1.0))
		  {
		    double extreme_h = el_pt->interpolated_x(extreme_s,1);
		    //Check whether the extreme value is greater than any of the nodes.
		    if(extreme_h > local_max) {local_max = extreme_h;}
		    if(extreme_h < local_min) {local_min = extreme_h;}
		  }
		  }

	    //Now check whether local max and min are global
	    if(local_max > max) {max = local_max;}
	    if(local_min < min) {min = local_min;}
	  }
      }
  }
    

  
 //Return the l2 norm of height difference between the interface
 //and its undeformed value
 double l2_norm_of_height(const double &h0)
  {
   double norm = 0.0;
   //Loop over the interface and add each elemental contribution
   const unsigned n_interface  = Surface_mesh_pt->nelement();
   for(unsigned i=0;i<n_interface;i++)
    {
     // Upcast from GeneralsedElement to the present element
     INTERFACE_ELEMENT *el_pt =
      dynamic_cast<INTERFACE_ELEMENT*>(
       Surface_mesh_pt->element_pt(i));
     
     norm += el_pt->l2_norm_of_height(h0);
    }
   return norm;
  }


  /// Return the total concentrations of the surfactant
  /// integrated over the bulk or surface accordingly
 void compute_integrated_concentrations(double &surface,
                                        double &bulk,
                                        double &micelle)
  {
   //Initialise to zero
   surface = 0.0;
   //Loop over the interface and add each elemental contribution
   const unsigned n_interface  = Surface_mesh_pt->nelement();
   for(unsigned i=0;i<n_interface;i++)
    {
     // Upcast from GeneralsedElement to the present element
     INTERFACE_ELEMENT *el_pt =
      dynamic_cast<INTERFACE_ELEMENT*>(
       Surface_mesh_pt->element_pt(i));
     
     surface += el_pt->integrated_C();
    }
   
   //Initialise to zero
   bulk = 0.0; micelle = 0.0;

   // Loop over bulk elements in lower fluid and add each elemental
   // contribution
   const unsigned n_lower = Bulk_mesh_pt->nlower();
   for(unsigned e=0;e<n_lower;e++)
    {
     // Upcast from GeneralisedElement to the present element
     ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->
                                             lower_layer_element_pt(e));
     double int_M=0.0, int_C=0.0;
     el_pt->integrated_C_and_M(int_C,int_M);
     bulk += int_C;
     micelle += int_M;
    }
  }



  //Return the areas of the upper and lower fluid
 void compute_areas(double &lower_area,
		    double &upper_area)
  {
    // Loop over bulk elements in lower fluid and add each elemental
    // contribution
    const unsigned n_lower = Bulk_mesh_pt->nlower();
    for(unsigned e=0;e<n_lower;e++)
      {
	lower_area += Bulk_mesh_pt->lower_layer_element_pt(e)->size();
      }
    const unsigned n_upper = Bulk_mesh_pt->nupper();
    for(unsigned e=0;e<n_upper;e++)
      {
	upper_area += Bulk_mesh_pt->upper_layer_element_pt(e)->size();
      }
  }


 
private:
 
 /// DocInfo object
 DocInfo Doc_info;

 /// Boolean to indicate whether the surface is pinned
 bool Surface_pinned;

  /// Node used to monitor the interface height
  Node* Monitor_node_pt;

  /// Integer used to specify the Face ID of the Lagrange multipliers
  /// Used to enforce periodicity
  unsigned Periodic_index;
  
 /// Enumeration of channel boundaries
 enum 
 {
  Inflow_boundary_id=3,
  Upper_wall_boundary_id=2,
  Outflow_boundary_id=1,
  Bottom_wall_boundary_id=0,
  Interface_boundary_id=4,
 };



}; // end of problem class

//===========start_of_constructor=========================================
/// Constructor for convection problem
//========================================================================
template<class ELEMENT, class INTERFACE_ELEMENT>
SurfactantProblem<ELEMENT,INTERFACE_ELEMENT>::
SurfactantProblem(const bool &pin) : Surface_pinned(pin), Periodic_index(50)
{
 //Allocate an (adaptive) timestepper
  add_time_stepper_pt(new BDF<2>(true));
 //Don't worry if the estimated error is above the tolerance
  Keep_temporal_error_below_tolerance = false;
  
 // Set output directory
 Doc_info.set_directory("RESLT");

 //Domain length in x direction
 const double lx = 2.0*Global_Physical_Variables::L;
 
 // Interface height
 const double h=Global_Physical_Variables::H0;
 // Total height of the domain
 const double H=1.0; 

 

//Setup the triangle mesh
 // Build the boundary segments for outer boundary, consisting of
 //--------------------------------------------------------------
 // four separate polylines
 //------------------------
 Vector<TriangleMeshCurveSection*> boundary_polyline_pt(4);

 //Specify how many elements in the horizontal direction and vertical direction in each layer
 //This is rquired to ensure that we have the same number of nodes on each side
 //and so the boundaries can be made periodic
 unsigned n_x = 10;
 unsigned n_y1 = 8;
 unsigned n_y2 = 8;

  if(n_x < 2) {std::cout  << "n_x needs to be 2 or more\n";}
  if(n_y1 < 2) {std::cout << "n_y1 needs to be 2 or more\n";}
  if(n_y2 < 2) {std::cout << "n_y2 needs to be 2 or more\n";}
 
 
 // Each horizontal polyline has n_x +1 vertices -- provide storage for their
 // coordinates
 Vector<Vector<double> > vertex_coord_h(n_x+1);
 for(unsigned i=0;i<(n_x+1);i++)
  {
   vertex_coord_h[i].resize(2);
  }

 //Each vertical polyline has n_y1 + n_y2 +1 points
 unsigned n_vertical_points = n_y1 + n_y2 + 1;

 Vector<Vector<double> > vertex_coord(n_vertical_points);
 for(unsigned i=0;i<n_vertical_points;i++)
  {
   vertex_coord[i].resize(2);
  }

 
 // First polyline: Inflow
 double dy_1 = h/((double)(n_y1));
 double dy_2 = (H-h)/((double)(n_y2));

//Lower film - lower left corner
 vertex_coord[0][0]=0.0;
 vertex_coord[0][1]=0.0;
 //Add intermediate nodes
 for(unsigned i=1;i<n_y1;++i)
  {
   vertex_coord[i][0] = 0.0;
   vertex_coord[i][1] = i*dy_1;
  }
 //Interface height
 vertex_coord[n_y1][0] = 0.0;
 vertex_coord[n_y1][1] = h;

 
 //Upper film
 for(unsigned i=1;i<n_y2;++i)
  {
   vertex_coord[n_y1+i][0] = 0.0;
   vertex_coord[n_y1+i][1] = h + i*dy_2;
  }
 
 //Override to avoid FP error
 vertex_coord[n_y1+n_y2][0] = 0.0;
 vertex_coord[n_y1+n_y2][1] = H;


 /*
 vertex_coord_h[0][0]=0.0;
 vertex_coord_h[0][1]=0.0;
 vertex_coord_h[1][0]=0.0;
 vertex_coord_h[1][1]= h;
 vertex_coord_h[2][0]=0.0;
 vertex_coord_h[2][1]= H;*/

 
 
 // Build the 1st boundary polyline
 boundary_polyline_pt[0] = new TriangleMeshPolyLine(vertex_coord,
                                                   Inflow_boundary_id);

 // Second boundary polyline: Upper wall
 double dx = lx/((double)(n_x));
 //First vertex
 vertex_coord_h[0][0]=0.0;
 vertex_coord_h[0][1]=H;
 for(unsigned i=1;i<n_x;++i)
  {
   vertex_coord_h[i][0] = i*dx;
   vertex_coord_h[i][1] = H;
  }
 //Final vertex
 vertex_coord_h[n_x][0]=lx;
 vertex_coord_h[n_x][1]=H;


/* vertex_coord_h[0][0]=0.0;
 vertex_coord_h[0][1]=H;
 vertex_coord_h[1][0]=0.5*lx;
 vertex_coord_h[1][1]=H;
 vertex_coord_h[2][0]=lx;
 vertex_coord_h[2][1]=H;*/

 // Build the 2nd boundary polyline
 boundary_polyline_pt[1] = new TriangleMeshPolyLine(vertex_coord_h,
                                                   Upper_wall_boundary_id);

 // Third boundary polyline: Outflow
 /*vertex_coord_h[0][0]=lx;
 vertex_coord_h[0][1]=H;
 vertex_coord_h[1][0]=lx;
 vertex_coord_h[1][1]=h;
 vertex_coord_h[2][0]=lx;
 vertex_coord_h[2][1]=0.0;*/

 
  //Lower film - upper right corner
 vertex_coord[0][0]=lx;
 vertex_coord[0][1]=H;
 //Add intermediate nodes
 for(unsigned i=1;i<n_y2;++i)
  {
   vertex_coord[i][0] = lx;
   vertex_coord[i][1] = H - i*dy_2;
  }
 //Interface height
 vertex_coord[n_y2][0] = lx;
 vertex_coord[n_y2][1] = h;

 
 //Upper film
 for(unsigned i=1;i<n_y1;++i)
  {
   vertex_coord[n_y2+i][0] = lx;
   vertex_coord[n_y2+i][1] = h - i*dy_1;
  }
 
 //Override to avoid FP error
 vertex_coord[n_y1+n_y2][0] = lx;
 vertex_coord[n_y1+n_y2][1] = 0.0; 


 // Build the 3rd boundary polyline
 boundary_polyline_pt[2] = new TriangleMeshPolyLine(vertex_coord,
                                                   Outflow_boundary_id);

 // Fourth boundary polyline: Bottom wall
 //First vertex
 vertex_coord_h[0][0]=lx;
 vertex_coord_h[0][1]=0.0;
 for(unsigned i=1;i<n_x;++i)
  {
   vertex_coord_h[i][0] = lx - i*dx;
   vertex_coord_h[i][1] = 0.0;
  }
 //Final vertex
 vertex_coord_h[n_x][0]=0.0;
 vertex_coord_h[n_x][1]=0.0;


 /*vertex_coord_h[0][0]=lx;
 vertex_coord_h[0][1]=0.0;
 vertex_coord_h[1][0]=0.5*lx;
 vertex_coord_h[1][1]=0.0;
 vertex_coord_h[2][0]=0.0;
 vertex_coord_h[2][1]=0.0;*/

 // Build the 4th boundary polyline
 boundary_polyline_pt[3] = new TriangleMeshPolyLine(vertex_coord_h,
                                                    Bottom_wall_boundary_id);
 
 // Create the triangle mesh polygon for outer boundary
 Outer_boundary_polyline_pt = new TriangleMeshPolygon(boundary_polyline_pt);
 
 //Here we need to put the dividing internal line in
 Vector<TriangleMeshOpenCurve *> interface_pt(1);
 //Set the vertex coordinates
 //First vertex
 vertex_coord_h[0][0]=0.0;
 vertex_coord_h[0][1]=h;
 for(unsigned i=1;i<n_x;++i)
  {
   vertex_coord_h[i][0] = i*dx;
   vertex_coord_h[i][1] = h;
  }
 //Final vertex
 vertex_coord_h[n_x][0]=lx;
 vertex_coord_h[n_x][1]=h;


/* vertex_coord_h[0][0]=0.0;
 vertex_coord_h[0][1]=h;
 vertex_coord_h[1][0]=0.5*lx;
 vertex_coord_h[1][1]=h;
 vertex_coord_h[2][0]=lx;
 vertex_coord_h[2][1]=h;*/
 

//Create the internal line
  TriangleMeshPolyLine* interface_polyline_pt =
   new TriangleMeshPolyLine(vertex_coord_h,
                            Interface_boundary_id);

  // Do the connection with the destination boundary, in this case
  // the connection is done with the inflow boundary. The connection is to the
  // specified node.
  interface_polyline_pt->connect_initial_vertex_to_polyline(
   dynamic_cast<TriangleMeshPolyLine*>(boundary_polyline_pt[0]),n_y1/*1*/);

  // Do the connection with the destination boundary, in this case
  // the connection is done with the outflow boundary
  interface_polyline_pt->connect_final_vertex_to_polyline(
   dynamic_cast<TriangleMeshPolyLine*>(boundary_polyline_pt[2]),n_y2/*1*/);
 
  Vector<TriangleMeshCurveSection*> interface_curve_pt(1);
  interface_curve_pt[0] = interface_polyline_pt;
  
  interface_pt[0] = new TriangleMeshOpenCurve(interface_curve_pt);
  
 // Now build the mesh, based on the boundaries specified by
 //---------------------------------------------------------
 // polygons just created
 //----------------------

 // Convert to "closed curve" objects
 TriangleMeshClosedCurve* outer_closed_curve_pt=Outer_boundary_polyline_pt;

 unsigned n_internal_closed_boundaries = 0;
 Vector<TriangleMeshClosedCurve *>
  inner_boundaries_pt(n_internal_closed_boundaries);
 
 // Target area for initial mesh
 double uniform_element_area=0.01;

 // Use the TriangleMeshParameter object for gathering all
 // the necessary arguments for the TriangleMesh object
 TriangleMeshParameters triangle_mesh_parameters(
   outer_closed_curve_pt);

 //Define the inner boundaries
 triangle_mesh_parameters.internal_closed_curve_pt() = inner_boundaries_pt;
 
 // Define the holes on the boundary
 triangle_mesh_parameters.internal_open_curves_pt() = interface_pt;

 // Define the maximum element area
 triangle_mesh_parameters.element_area() =
   uniform_element_area;

 Vector<double> lower_region(2);
 lower_region[0] = 0.5*lx;
 lower_region[1] = 0.5*h;
 
 // Define the region
 triangle_mesh_parameters.add_region_coordinates(1, lower_region);

 //Don't allow addition of new nodes on boundaries
 triangle_mesh_parameters.disable_automatic_creation_of_vertices_on_boundaries();
 
 // Create the mesh
 Bulk_mesh_pt =
   new RefineableSolidTwoLayerTriangleMesh<ELEMENT>(
     triangle_mesh_parameters, this->time_stepper_pt());

 Bulk_mesh_pt->disable_automatic_creation_of_vertices_on_boundaries();
 

 // Create and set the error estimator for spatial adaptivity
 Bulk_mesh_pt->spatial_error_estimator_pt() = new Z2ErrorEstimator;

 // Set the maximum refinement level for the mesh to 4
 //Bulk_mesh_pt->max_refinement_level() = 4;

 // Set error targets for refinement
 Bulk_mesh_pt->max_permitted_error() = 1.0e-3;
 Bulk_mesh_pt->min_permitted_error() = 1.0e-5;

 //Set a smaller minimum element size
 Bulk_mesh_pt->min_element_size() = 0.001; //0.0001
 //Set a smaller maximum element size as well
 Bulk_mesh_pt->max_element_size() = 0.1;
 
 //Set the refinement tolerance for the interface
 /*interface_polyline_pt->set_refinement_tolerance(0.05);
 interface_polyline_pt->set_unrefinement_tolerance(0.001);
 //Need to be consistent on the side walls
 for(unsigned i=0;i<4;++i)
  {
   boundary_polyline_pt[i]->set_refinement_tolerance(0.05);
   boundary_polyline_pt[i]->set_unrefinement_tolerance(0.001);
   }*/

 
 // Define a constitutive law for the solid equations: generalised Hookean
 Constitutive_law_pt = new GeneralisedHookean(&Global_Physical_Variables::Nu);
 
 // Complete the build of all elements so they are fully functional 
 
 // Determine number of bulk elements in lower/upper fluids
 const unsigned n_upper = Bulk_mesh_pt->nupper();
 const unsigned n_lower = Bulk_mesh_pt->nlower();


 //std::ofstream f_lower("lower.dat");
 //std::ofstream f_upper("upper.dat");
 
 // Loop over bulk elements in lower fluid
 for(unsigned e=0;e<n_lower;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->
                                           lower_layer_element_pt(e));

   //Output element to the file, if required during debugging
   //el_pt->output(f_lower,2);
   
   // Set the diffusivities number
   el_pt->diff_pt() = &Global_Physical_Variables::D;

   // Set the timescales
   el_pt->tau_pt() =&Global_Physical_Variables::Tau;

   // Set the Reynolds number
   el_pt->re_pt() = &Global_Physical_Variables::Re;

   // Set the Womersley number
   el_pt->re_st_pt() = &Global_Physical_Variables::Re;

   // Set the product of the Reynolds number and the inverse of the
   // Froude number
   el_pt->re_invfr_pt() = &Global_Physical_Variables::ReInvFr;

   // Set the km parameter
   el_pt->km_pt() = &Global_Physical_Variables::K_m;

   // Set the N parameter
   el_pt->n_pt() = &Global_Physical_Variables::N;

   
   // Set the direction of gravity
   el_pt->g_pt() = &Global_Physical_Variables::Direction_of_gravity;

   // Set the constitutive law
   el_pt->constitutive_law_pt() = this->Constitutive_law_pt;
   
  } // End of loop over bulk elements in lower fluid

 
 // Loop over bulk elements in upper fluid 
 for(unsigned e=0;e<n_upper;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->
                                           upper_layer_element_pt(e));

   //Output element to file if required during debugging
   //el_pt->output(f_upper,2);
   
   // Set the diffusivities number
   el_pt->diff_pt() = &Global_Physical_Variables::D;

   // Set the timescales
   el_pt->tau_pt() =&Global_Physical_Variables::Tau;
   
   // Set the Reynolds number
   el_pt->re_pt() = &Global_Physical_Variables::Re;

   // Set the Womersley number
   el_pt->re_st_pt() = &Global_Physical_Variables::Re;

   // Set the product of the Reynolds number and the inverse of the
   // Froude number
   el_pt->re_invfr_pt() = &Global_Physical_Variables::ReInvFr;

   // Set the viscosity ratio
   el_pt->viscosity_ratio_pt() = &Global_Physical_Variables::M;

   // Set the density ratio
   el_pt->density_ratio_pt() = &Global_Physical_Variables::R;

   // Set the km parameter
   el_pt->km_pt() = &Global_Physical_Variables::K_m;

   // Set the N parameter
   el_pt->n_pt() = &Global_Physical_Variables::N;
   
   // Set the direction of gravity
   el_pt->g_pt() = &Global_Physical_Variables::Direction_of_gravity;

   // Set the constitutive law
   el_pt->constitutive_law_pt() = Constitutive_law_pt;

   //Turn off the solution of the advection-diffusion-reaction equations
   el_pt->disable_advection_diffusion_reaction_equations();
   
   //Need to pin the values of concentration and micelle up here
   unsigned n_el_node = el_pt->nnode();
   for(unsigned n=0;n<n_el_node;++n)
     {
       el_pt->node_pt(n)->set_value(2,0.0);
       el_pt->node_pt(n)->set_value(3,0.0);
       el_pt->node_pt(n)->pin(2);
       el_pt->node_pt(n)->pin(3);
     }
   
  } // End of loop over bulk elements in upper fluid

 
 //Create the surface mesh that will contain the interface elements
 //First create storage, but with no elements or nodes
 Surface_mesh_pt = new Mesh;
 //Make point elements at the end to compensate for
 //the unbalanced line tension terms 
 //if we DON'T have periodic boundaryc conditions
 if(!Control_Parameters::Periodic_BCs)
   {
     Point_mesh_pt = new Mesh;
   }


 //Must be done *after* the physical parameters have been set
 //because we are using the viscosity ratio pointer to determine
 //upper or lower elements
 create_interface_elements();
 //Setup the monitor node
 {
  //Choose one of the boundary nodes otherwise gets too tricky
  //Loop over nodes on the surface
  unsigned n_node = Bulk_mesh_pt->nboundary_node(4);
  for(unsigned n=0;n<n_node;++n)
   {
    Node* node_pt = Bulk_mesh_pt->boundary_node_pt(4,n);
    //Choose the right hand edge
    if(node_pt->is_on_boundary(1))
     {
      Monitor_node_pt = node_pt;
      break;
     }
   }
 }
  
  //Making periodic boundaries is a bit more involved for an unstructured mesh 
  if(Control_Parameters::Periodic_BCs)
   {
    this->setup_periodic_boundaries();
    
    Dependent_position_mesh_pt = new Mesh;
    //Create the dependent elements that ensure that the positions match
    create_dependent_position_elements();
   }
  
 // Add the two sub-meshes to the problem
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Surface_mesh_pt);
 if(!Control_Parameters::Periodic_BCs)
   {
    //add_sub_mesh(Point_mesh_pt);
    //add_sub_mesh(Inlet_traction_mesh_pt);
    //add_sub_mesh(Outlet_traction_mesh_pt);
   }
 else
   {
     add_sub_mesh(Dependent_position_mesh_pt);
     }
 
 // Combine all sub-meshes into a single mesh
 build_global_mesh();

 
 //Pin the positions of all nodes and the Lagrange multipliers if the surface is pinned
 if(Surface_pinned)
  {
   unsigned n_node = Bulk_mesh_pt->nnode();
   for(unsigned n=0;n<n_node;n++)
     {
       SolidNode* nod_pt = Bulk_mesh_pt->node_pt(n);
       nod_pt->pin_position(0);
       nod_pt->pin_position(1);
       //Don't forget to pin the Lagrange multipliers as well
       if(nod_pt->is_on_boundary(4))
        {
         //Lagrange multiplier is always value 5
         unsigned lagrange_multiplier_index = 5;
         nod_pt->pin(lagrange_multiplier_index);
         
         //There is one additional Lagrange multiplier for the Periodic
         //node on the interface that also needs to be pinned
         if((Control_Parameters::Periodic_BCs) &&
            (nod_pt->is_on_boundary(1)))
          {
           unsigned n_value = dynamic_cast<BoundaryNodeBase*>(nod_pt)
            ->nvalue_assigned_by_face_element(Periodic_index);
           //Check that the Lagrange multiplier has been allocated
           if(n_value > 0)
            {
             unsigned periodic_lm_index = dynamic_cast<BoundaryNodeBase*>(nod_pt)
             ->index_of_first_value_assigned_by_face_element(Periodic_index);
             //Pin this lagrange multiplier
             nod_pt->pin(periodic_lm_index);
            }
           }
	}
     }
  }


 // Set the boundary conditions for this problem: All nodes are
 // free by default -- only need to pin the ones that have Dirichlet 
 // conditions here
 
 //Loop over the boundaries
 unsigned num_bound = Bulk_mesh_pt->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   //If we are on the side-walls, the concentrations
   //satisfy natural boundary conditions, so we only pin the
   //v-velocity for now if not periodic
    if(!Control_Parameters::Periodic_BCs)
      {
	if((ibound==1) || (ibound==3))
	  {
	    //Loop over the number of nodes on the boundary
	    unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
	    for(unsigned inod=0;inod<num_nod;inod++)
	      {
               //Cast to a solid node, so that we can impose boundary conditions
               SolidNode* nod_pt = static_cast<SolidNode*>(
                Bulk_mesh_pt->boundary_node_pt(ibound,inod));
		nod_pt->pin(1);
                //Also pin the horizontal displacement of the nodes
                nod_pt->pin_position(0);
                
		/*if((ibound==4) || (ibound==5))
		  {
		  Bulk_mesh_pt->boundary_node_pt(ibound,inod)->pin(0);
		  }*/
	      }
	  }
      }

   //If we on the top or bottom wall, velocity is pinned
   //as is the nodal position
   if((ibound==0) || (ibound==2))
    {
      //Loop over the number of nodes on the boundary
      unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
      for(unsigned inod=0;inod<num_nod;inod++)
	{
         SolidNode* nod_pt = static_cast<SolidNode*>(Bulk_mesh_pt->boundary_node_pt(ibound,inod));
	  nod_pt->pin(0);
	  nod_pt->pin(1);
          nod_pt->pin_position(0);
          nod_pt->pin_position(1);
	}
    }
  }

 
 //Pin the zero-th pressure dof in element 0 and set its value to
 //zero:
 fix_pressure(0,0,0.0);


  // Loop over the interface elements to set up element-specific 
 // things that cannot be handled by the (argument-free!) ELEMENT 
 // constructor.
 unsigned n_interface  = Surface_mesh_pt->nelement();
 for(unsigned i=0;i<n_interface;i++)
  {
   // Upcast from GeneralsedElement to the present element
   INTERFACE_ELEMENT *el_pt = dynamic_cast<INTERFACE_ELEMENT*>(
    Surface_mesh_pt->element_pt(i));
   
   //Need to unpin the values of concentration and micelle on the interface
   unsigned n_el_node = el_pt->nnode();
   for(unsigned n=0;n<n_el_node;++n)
     {
       el_pt->node_pt(n)->unpin(2); //unpin surfactant
       if(Control_Parameters::Pin_Micelle)
        {
         el_pt->node_pt(n)->pin(3); //pin micelle concentration
        }
       else
        {
         el_pt->node_pt(n)->unpin(3); //unpin micelle concentration
        }
       el_pt->node_pt(n)->pin(4); //Pin Lagrange Multiplier initially
     }

  }

 //Pin all bulk concentrations ... to save
 //time and to avoid having to conserve global mass
 unsigned n_node = Bulk_mesh_pt->nnode();
 for(unsigned n=0;n<n_node;++n)
   {
    SolidNode* nod_pt = static_cast<SolidNode*>(Bulk_mesh_pt->node_pt(n));
    //Pin x position
     nod_pt->pin_position(0);
     nod_pt->pin_position(1);
     nod_pt->pin(2);
     nod_pt->pin(3);
     double y = nod_pt->x(1);
     nod_pt->set_value(0,y);
     nod_pt->set_value(1,0.0);
   }

 //We initialise need to pin the lagrange multipliers
 //associated with periodic boundary conditions on the side boundary
 if(Control_Parameters::Periodic_BCs)
 {
   unsigned b=3;
   unsigned n_boundary_node = this->Bulk_mesh_pt->nboundary_node(b);
   for(unsigned n=0;n<n_boundary_node;++n)
     {
       Node* nod_pt = this->Bulk_mesh_pt->boundary_node_pt(b,n);
       if(!(nod_pt->is_on_boundary(0)) &&
	  !(nod_pt->is_on_boundary(2)) &&
	  !(nod_pt->is_on_boundary(4)))
	 {
	   //Pin all
	   nod_pt->pin(dynamic_cast<BoundaryNodeBase*>(nod_pt)
		       ->index_of_first_value_assigned_by_face_element(
								       Periodic_index));
	 }
     }
 }

 
 //Pin one surface concentration at the surface, if only
 //solving for insoluble surfactant
 //Surface_mesh_pt->finite_element_pt(0)->node_pt(1)->pin(4);
  
 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << endl; 
    
} // end of constructor



//===========start_of_set_boundary_conditions================
/// Set the boundary conditions as a function of continuous 
/// time
//===========================================================
template<class ELEMENT,class INTERFACE_ELEMENT>
void SurfactantProblem<ELEMENT,INTERFACE_ELEMENT>::set_boundary_conditions(
 const double &time)
{
 //Set initial temperature profile
 if(time <= 0.0)
  {
   const double Gamma_init = 0.1; //0.5; //0.1; //0.75;
   const double C_init = Gamma_init/(Global_Physical_Variables::K_b*(1.0-Gamma_init));
   double M_init = 0.0; 
   //If not pinning the micelle concentration then set the initial value 
   if(!Control_Parameters::Pin_Micelle)
    {
     M_init = pow(C_init,Global_Physical_Variables::N);
    }
    
    unsigned n_lower = Bulk_mesh_pt->nlower();
    for(unsigned e=0;e<n_lower;e++)
      {
	FiniteElement* el_pt =
	  Bulk_mesh_pt->lower_layer_element_pt(e);
	unsigned n_node = el_pt->nnode();
	for(unsigned n=0;n<n_node;n++)
	  {
	    Node* nod_pt = el_pt->node_pt(n);
	    //And in uniformly distributed surfactant
	    //Be careful about upper and lower layers
	    //If these are not set 
	    nod_pt->set_value(2,C_init); 
	    nod_pt->set_value(3,M_init); 
	    
	    //Set the velocityq
	    /*double y = nod_pt->x(1);
	      nod_pt->set_value(0,y);
	      nod_pt->set_value(1,0.0);*/
	  }
      }
  
   //Set the initial surface concentration to be one
   unsigned n_surface_element = Surface_mesh_pt->nelement();
   for(unsigned e=0;e<n_surface_element;++e)
     {
       unsigned n_el_node = Surface_mesh_pt->finite_element_pt(e)->nnode();
       for(unsigned n=0;n<n_el_node;++n)
	 {
	   Surface_mesh_pt->finite_element_pt(e)->node_pt(n)->set_value(4,Gamma_init);
	 }
     }
  } //End of initial conditions

 // Loop over the boundaries
 unsigned num_bound = Bulk_mesh_pt->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   // Loop over the nodes on boundary 
   unsigned num_nod=Bulk_mesh_pt->nboundary_node(ibound);
   for(unsigned inod=0;inod<num_nod;inod++)
    {
     // Get pointer to node
     Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);

     if(!Control_Parameters::Periodic_BCs)
       {
	 //If we are on the side walls we only set the v-velocity.
        if((ibound==1) || (ibound==3))
	   {
	     nod_pt->set_value(1,0.0); 
	   }
       }

     //If we are on the top boundary, do not set the velocities
     //(yet)
     if(ibound==2)
       {
	 nod_pt->set_value(0,1.0); nod_pt->set_value(1,0.0);
       }

     //If we are on the bottom boundary
     if(ibound==0)
       {
	 nod_pt->set_value(0,0.0); nod_pt->set_value(1,0.0);
       }
    }
  }
} // end_of_set_boundary_conditions

//===============start_doc_solution=======================================
/// Doc the solution
//========================================================================
template<class ELEMENT,class INTERFACE_ELEMENT>
void SurfactantProblem<ELEMENT,INTERFACE_ELEMENT>::doc_solution(
 ofstream &trace)
{ 
 //Declare an output stream and filename
 ofstream some_file;
 char filename[100];

 // Number of plot points: npts x npts
 unsigned npts=5;

 // Output solution 
 //-----------------
 sprintf(filename,"%s/soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 unsigned n_element = Bulk_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   Bulk_mesh_pt->finite_element_pt(e)->output(some_file,npts);
  }
 some_file.close();

 //Output the interface
 sprintf(filename,"%s/int%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);

 unsigned n_interface  = Surface_mesh_pt->nelement();
 for(unsigned i=0;i<n_interface;i++)
  {
   Surface_mesh_pt->finite_element_pt(i)->output(some_file,npts);
  }
 some_file.close();

 //Let's get the mases
 double surface=0.0, bulk=0.0, micelle=0.0;
 this->compute_integrated_concentrations(surface,bulk,micelle);

 //double upper_area = 0.0, lower_area = 0.0;
 //this->compute_areas(lower_area,upper_area);

 //Compute the max and min
 double max = 0.0; double min = 0.0;
 this->interface_min_max(min,max);
 
 trace << time_pt()->time() << " " 
       << this->Monitor_node_pt->x(1) << " " 
       << this->Monitor_node_pt->value(2) << " "
       << this->Monitor_node_pt->value(3) << " "
       << this->Monitor_node_pt->value(4) << " "
       << std::sqrt(this->l2_norm_of_height(Global_Physical_Variables::H0)/(2.0*Global_Physical_Variables::L)) << " "
       << surface << " " << bulk << " " << micelle << " "
       << surface + (bulk + micelle)/(Global_Physical_Variables::Beta_b)
       << " " << min << " " << max << " " 
   //<< " " << upper_area << " " << lower_area << " " << upper_area + lower_area 
       << std::endl;


 Doc_info.number()++;
} // end of doc



//Specify the global error norm
template<class ELEMENT,class INTERFACE_ELEMENT>
double SurfactantProblem<ELEMENT,INTERFACE_ELEMENT>::global_temporal_error_norm()
{
 //Temp
 double global_error = 0.0;

 //Base the error norm on the interface motion
 //Number of interface nodes
 unsigned n_node = Bulk_mesh_pt->nboundary_node(4);

 //Loop over the nodes and calculate the errors in the positions
 for(unsigned long n=0;n<n_node;n++)
  {
   Node* node_pt = Bulk_mesh_pt->boundary_node_pt(4,n);
   
   //Find number of dimensions of the node
   unsigned n_dim = node_pt->ndim();
   //Set the position error to zero
   double node_position_error = 0.0;
   //Loop over the dimensions 
   for(unsigned j=0;j<n_dim;j++)
   {
     //Get position error
     double error = 
      node_pt->position_time_stepper_pt()->
      temporal_error_in_position(node_pt,j);

     //Add the square of the individual error to the position error
     node_position_error += error*error;
    }
    
   //Divide the position error by the number of dimensions
   node_position_error /= n_dim;
   //Now add to the global error
   global_error += node_position_error;
  }
   
   //Now the global error must be divided by the number of nodes
 global_error /= n_node;

 //Return the square root of the errr
 return sqrt(global_error);
}
 

//=======start_of_main================================================
/// Driver code for 2D Boussinesq convection problem
//====================================================================
int main(int argc, char **argv)
{
 ofstream trace("RESLT/trace.dat");

 // Set the direction of gravity
 Global_Physical_Variables::Direction_of_gravity[0] = 0.0;
 Global_Physical_Variables::Direction_of_gravity[1] = -1.0;

 //Set the diffusivities (inverse peclect numbers)
 Global_Physical_Variables::D[0] = 1.0/Global_Physical_Variables::Pe_b;
 Global_Physical_Variables::D[1] = 1.0/Global_Physical_Variables::Pe_m;

 Global_Physical_Variables::Wall_normal.resize(2);
 Global_Physical_Variables::Wall_normal[0] = -1.0;
 Global_Physical_Variables::Wall_normal[1] = 0.0;
 
 
 //Construct our problem
 SurfactantProblem<Hijacked<
  ProjectableDoubleBuoyantElement<PseudoSolidNodeUpdateElement<DoubleBuoyantTCrouzeixRaviartElement<2>, TPVDElement<2,3> > > >,
		   ElasticLineSolubleSurfactantTransportInterfaceElement<
                    ProjectableDoubleBuoyantElement<PseudoSolidNodeUpdateElement<DoubleBuoyantTCrouzeixRaviartElement<2> , TPVDElement<2,3> > >  > > 
                   problem;

 
 // Apply the boundary condition at time zero
 problem.set_boundary_conditions(0.0);
 
 //Perform a single steady Newton solve
 problem.steady_newton_solve();

 //Document the solution
 problem.doc_solution(trace);

 //Now release the interface for real fun
 problem.unpin_surface();

 //Set the timestep
 double dt = 0.01; //0.1;

 //Initialise the value of the timestep and set initial values 
 //of previous time levels assuming an impulsive start.
 problem.deform_interface(Global_Physical_Variables::Ha,1);
 problem.doc_solution(trace);
 problem.assign_initial_values_impulsive(dt);

 //Set the number of timesteps to our default value
 unsigned n_steps = 1000;

 //Set the freqency of refinement
 unsigned refine_after_n_steps = 10;
 
 //If we have a command line argument, perform fewer steps 
 //(used for self-test runs)
 if(argc > 1) {n_steps = 5; refine_after_n_steps=3;}


 
 //This is the code for reloading a previously dumped solution
 //If you want to reload then uncomment this 
 /* {
   //The dumpfile must be renamed to restart_file.dat
   std::ifstream restart_file("restart_file.dat");
   bool unsteady_restart;
   //Load the restart file
   problem.read(restart_file,unsteady_restart);
 
   //Note that paramater values are not saved in the restart file, so this will
   //not work if the parameter values have been changed in the driver code.
   //This can be modified if you want to dump all the parameters as well.

   //The next value of dt may need to be adjusted for a restart.

  }*/
 
  
 
 //Perform n_steps timesteps
 for(unsigned i=0;i<n_steps;++i)
  {
   double dt_next = dt;
   //Spatial adapt every 5 steps with a fixed timestep
   if((i>0) && (i%refine_after_n_steps==0))
    {
      problem.unsteady_newton_solve(dt,1,false);
    }
    else
    {
     //Unsteady timestep (for fixed mesh)
     dt_next = problem.adaptive_unsteady_newton_solve(dt,1.0e-5);
    }
   dt = dt_next;
   problem.doc_solution(trace);
   //Dump the most recently converged solution in a form suitable for restarts.
   //This will be overwritten after each successful solution.
   problem.dump("dumpfile.dat");
  }

} // end of main









