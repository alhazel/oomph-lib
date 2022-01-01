//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2022 Matthias Heil and Andrew Hazel
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
// Driver code for a simple unstructured solid problem using a mesh
// generated from an input file generated by the triangle mesh generator
// Triangle.

//Generic routines
#include "generic.h"
#include "solid.h"
#include "constitutive.h"


// The mesh
#include "meshes/triangle_mesh.h"

using namespace std;
using namespace oomph;

/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////



//=======start_namespace==========================================
/// Global variables
//================================================================
namespace Global_Physical_Variables
{
 /// Poisson's ratio
 double Nu=0.3;

 /// Pointer to constitutive law
 ConstitutiveLaw* Constitutive_law_pt=0;
 
 /// Uniform pressure
 double P = 0.0;

 /// Constant pressure load. The arguments to this function are imposed
 /// on us by the SolidTractionElements which allow the traction to 
 /// depend on the Lagrangian and Eulerian coordinates x and xi, and on the 
 /// outer unit normal to the surface. Here we only need the outer unit
 /// normal.
 void constant_pressure(const Vector<double> &xi, const Vector<double> &x,
                        const Vector<double> &n, Vector<double> &traction)
 {
  unsigned dim = traction.size();
  for(unsigned i=0;i<dim;i++)
   {
    traction[i] = -P*n[i];
   }
 } 

} //end namespace



//==============start_problem=========================================
/// Unstructured solid problem
//====================================================================
template<class ELEMENT> 
class UnstructuredSolidProblem : public Problem
{

public:

 /// Constructor: 
 UnstructuredSolidProblem();

 /// Destructor (empty)
 ~UnstructuredSolidProblem(){}
 
 /// Set the problem to be incompressible
 void set_incompressible() {Incompressible=true;}

 /// Doc the solution
 void doc_solution(DocInfo& doc_info);

 /// Calculate the strain energy
 double get_strain_energy();

 /// Remove Traction Mesh
 void actions_before_adapt();

 /// Add on the traction elements after adaptation
 void actions_after_adapt();

private:
 
 /// Bulk mesh
 RefineableSolidTriangleMesh<ELEMENT>* Solid_mesh_pt;
 
 /// Pointer to mesh of traction elements
 SolidMesh* Traction_mesh_pt;

 /// Triangle mesh polygon for outer boundary 
 TriangleMeshPolygon* Outer_boundary_polyline_pt; 

 /// Boolean flag used in an incompressible problem
 bool Incompressible;

};



//===============start_constructor========================================
/// Constructor for unstructured solid problem
//========================================================================
template<class ELEMENT>
UnstructuredSolidProblem<ELEMENT>::UnstructuredSolidProblem() :
 Incompressible(false)
{  
 // Build the boundary segments for outer boundary, consisting of
 //--------------------------------------------------------------
 // four separeate polyline segments
 //---------------------------------
 Vector<TriangleMeshCurveSection*> boundary_segment_pt(4);
 
 // Initialize boundary segment
 Vector<Vector<double> > bound_seg(2);
 for(unsigned i=0;i<2;i++) {bound_seg[i].resize(2);}
 
 // First boundary segment
 bound_seg[0][0]=0.0;
 bound_seg[0][1]=0.0;
 bound_seg[1][0]=0.0;
 bound_seg[1][1]=5.0;
 
 // Specify 1st boundary id
 unsigned bound_id = 0;

 // Build the 1st boundary segment
 boundary_segment_pt[0] = new TriangleMeshPolyLine(bound_seg,bound_id);
 
 // Second boundary segment
 bound_seg[0][0]=0.0;
 bound_seg[0][1]=5.0;
 bound_seg[1][0]=1.0;
 bound_seg[1][1]=5.0;

 // Specify 2nd boundary id
 bound_id = 1;

 // Build the 2nd boundary segment
 boundary_segment_pt[1] = new TriangleMeshPolyLine(bound_seg,bound_id);

 // Third boundary segment
 bound_seg[0][0]=1.0;
 bound_seg[0][1]=5.0;
 bound_seg[1][0]=1.0;
 bound_seg[1][1]=0.0;

 // Specify 3rd boundary id
 bound_id = 2;

 // Build the 3rd boundary segment
 boundary_segment_pt[2] = new TriangleMeshPolyLine(bound_seg,bound_id);

 // Fourth boundary segment
 bound_seg[0][0]=1.0;
 bound_seg[0][1]=0.0;
 bound_seg[1][0]=0.0;
 bound_seg[1][1]=0.0;

 // Specify 4th boundary id
 bound_id = 3;

 // Build the 4th boundary segment
 boundary_segment_pt[3] = new TriangleMeshPolyLine(bound_seg,bound_id);
  
 // Create the triangle mesh polygon for outer boundary using boundary segment
 Outer_boundary_polyline_pt = new TriangleMeshPolygon(boundary_segment_pt);


 // There are no holes
 //-------------------------------
 
 // Now build the mesh, based on the boundaries specified by
 //---------------------------------------------------------
 // polygons just created
 //----------------------
 double uniform_element_area=0.2;

 TriangleMeshClosedCurve* closed_curve_pt=Outer_boundary_polyline_pt;

 // Use the TriangleMeshParameters object for gathering all
 // the necessary arguments for the TriangleMesh object
 TriangleMeshParameters triangle_mesh_parameters(
   closed_curve_pt);

 // Define the maximum element area
 triangle_mesh_parameters.element_area() =
   uniform_element_area;

 // Create the mesh
 Solid_mesh_pt =
   new RefineableSolidTriangleMesh<ELEMENT>(
     triangle_mesh_parameters);
 
 //hierher
 // Disable the use of an iterative solver for the projection
 // stage during mesh adaptation
 Solid_mesh_pt->disable_iterative_solver_for_projection();
 
 // Set error estimator for bulk mesh
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 Solid_mesh_pt->spatial_error_estimator_pt()=error_estimator_pt;


 // Set targets for spatial adaptivity
 Solid_mesh_pt->max_permitted_error()=0.0001;
 Solid_mesh_pt->min_permitted_error()=0.001; 
 Solid_mesh_pt->max_element_size()=0.2;
 Solid_mesh_pt->min_element_size()=0.001; 
   
 // Output mesh boundaries
 this->Solid_mesh_pt->output_boundaries("boundaries.dat");

 // Make the traction mesh
 Traction_mesh_pt=new SolidMesh;
 
 // Add sub meshes
 add_sub_mesh(Solid_mesh_pt);
 add_sub_mesh(Traction_mesh_pt);
 
 // Build the global mesh
 build_global_mesh();

 //Call actions after adapt:
 // 1) to build the traction elements
 // 2) to pin the nodes on the lower boundary (boundary 3)
 // 3) to complete the build of the elements
 // Note there is slight duplication here because we rebuild the global mesh
 // twice.
 this->actions_after_adapt();
   
 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 
 
} //end constructor


//========================================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void UnstructuredSolidProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{ 

 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5;

 // Output solution
 //----------------
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Solid_mesh_pt->output(some_file,npts);
 some_file.close();

 // Output traction
 //----------------
 sprintf(filename,"%s/traction%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Traction_mesh_pt->output(some_file,npts);
 some_file.close();

 // Output boundaries
 //------------------
 sprintf(filename,"%s/boundaries.dat",doc_info.directory().c_str());
 some_file.open(filename);
 Solid_mesh_pt->output_boundaries(some_file);
 some_file.close();

}

//================start_get_strain_energy================================
/// Calculate the strain energy in the entire elastic solid
//=======================================================================
template<class ELEMENT>
double UnstructuredSolidProblem<ELEMENT>::get_strain_energy()
{
 double strain_energy=0.0;
 const unsigned n_element = Solid_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   //Cast to a solid element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Solid_mesh_pt->element_pt(e));
   
   double pot_en, kin_en;
   el_pt->get_energy(pot_en,kin_en);
   strain_energy += pot_en;
  }
 
 return strain_energy;
} // end_get_strain_energy


//==============start_actions_before_adapt================================
/// Actions before adapt: remove the traction elements in the surface mesh
//========================================================================
template<class ELEMENT>
void UnstructuredSolidProblem<ELEMENT>::actions_before_adapt()
{
 // How many surface elements are in the surface mesh
 unsigned n_element = Traction_mesh_pt->nelement();
 
 // Loop over the surface elements and kill them
 for(unsigned e=0;e<n_element;e++) {delete Traction_mesh_pt->element_pt(e);}
 
 // Wipe the mesh
 Traction_mesh_pt->flush_element_and_node_storage();

} // end_actions_before_adapt

//=================start_actions_after_adapt=============================
 /// Need to add on the traction elements after adaptation
//=======================================================================
template<class ELEMENT>
void UnstructuredSolidProblem<ELEMENT>::actions_after_adapt()
{
 //The boundary in question is boundary 0
 unsigned b=0;
 
 // How many bulk elements are adjacent to boundary b?
 unsigned n_element = Solid_mesh_pt->nboundary_element(b);
 // Loop over the bulk elements adjacent to boundary b
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk element that is adjacent to boundary b
   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
    Solid_mesh_pt->boundary_element_pt(b,e));
   
   //Find the index of the face of element e along boundary b
   int face_index = Solid_mesh_pt->face_index_at_boundary(b,e);
   
   //Create solid traction element
   SolidTractionElement<ELEMENT> *el_pt = 
    new SolidTractionElement<ELEMENT>(bulk_elem_pt,face_index);   
   
   // Add to mesh
   Traction_mesh_pt->add_element_pt(el_pt);
   
   //Set the traction function
   el_pt->traction_fct_pt() = Global_Physical_Variables::constant_pressure;
  }  
 
 //Now rebuild the global mesh
 this->rebuild_global_mesh();
 
 //(Re)set the boundary conditions
 //Pin both positions at lower boundary (boundary 3)
 unsigned ibound=3;
 unsigned num_nod= mesh_pt()->nboundary_node(ibound);
 for (unsigned inod=0;inod<num_nod;inod++)
  {  
   // Get node
   SolidNode* nod_pt=Solid_mesh_pt->boundary_node_pt(ibound,inod);
   
   // Pin both directions
   for (unsigned i=0;i<2;i++) {nod_pt->pin_position(i);}
  }
 //End of set boundary conditions 
 
 // Complete the build of all elements so they are fully functional
 n_element = Solid_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   //Cast to a solid element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Solid_mesh_pt->element_pt(e));
   
   // Set the constitutive law
   el_pt->constitutive_law_pt() =
    Global_Physical_Variables::Constitutive_law_pt;
   
   //Set the incompressibility flag if required
   if(Incompressible) 
    {
     //Need another dynamic cast
     dynamic_cast<TPVDElementWithContinuousPressure<2>*>(el_pt)
      ->set_incompressible();
    }
  }

} // end_actions_after_adapt


//===========start_main===================================================
/// Demonstrate how to solve an unstructured solid problem
//========================================================================
int main(int argc, char **argv)
{

 //Doc info object
 DocInfo doc_info;

 // Output directory
 doc_info.set_directory("RESLT");
 
 // Create generalised Hookean constitutive equations
 Global_Physical_Variables::Constitutive_law_pt = 
  new GeneralisedHookean(&Global_Physical_Variables::Nu);
 
 {
  std::ofstream strain("RESLT/s_energy.dat");
  std::cout << "Running with pure displacement formulation\n";

  //Set up the problem
  UnstructuredSolidProblem<ProjectablePVDElement<TPVDElement<2,3> > > problem;
  
  //Output initial configuration
  problem.doc_solution(doc_info);
  doc_info.number()++;
  
  // Parameter study
  Global_Physical_Variables::P=0.0;
  double pressure_increment=0.1e-2;
  
  unsigned nstep=5;

  for (unsigned istep=0;istep<nstep;istep++)
   {
    // Solve the problem with one round of adaptivity
    problem.newton_solve(1);

    double strain_energy = problem.get_strain_energy();
    std::cout << "Strain energy is " << strain_energy << "\n";
    //Output strain energy to file
    strain << Global_Physical_Variables::P << " " << strain_energy << std::endl;

    //Output solution
    problem.doc_solution(doc_info);
    doc_info.number()++;
    
    //Reverse direction of increment 
    if(istep==2) {pressure_increment *= -1.0;}

    // Increase (or decrease) load
    Global_Physical_Variables::P+=pressure_increment;
   }

  strain.close();
 } //end_displacement_formulation


 //Repeat for displacement/pressure formulation 
 {
  std::ofstream strain("RESLT_pres_disp/s_energy.dat");
  std::cout << "Running with pressure/displacement formulation\n";

  // Change output directory
  doc_info.set_directory("RESLT_pres_disp");
  //Reset doc_info number
  doc_info.number() = 0;
  
  //Set up the problem
  UnstructuredSolidProblem<
  ProjectablePVDElementWithContinuousPressure<
  TPVDElementWithContinuousPressure<2> > > problem;
  
  //Output initial configuration
  problem.doc_solution(doc_info);
  doc_info.number()++;
  
  // Parameter study
  Global_Physical_Variables::P=0.0;
  double pressure_increment=0.1e-2;
  
  unsigned nstep=5;
  for (unsigned istep=0;istep<nstep;istep++)
   {
    // Solve the problem
    problem.newton_solve(1);

    double strain_energy = problem.get_strain_energy();
    std::cout << "Strain energy is "<< strain_energy << "\n";
    //Output strain energy to file
    strain << Global_Physical_Variables::P << " " << strain_energy << std::endl;
    
    //Output solution
    problem.doc_solution(doc_info);
    doc_info.number()++;

    if(istep==2) {pressure_increment *= -1.0;}    
    // Increase (or decrease) pressure load
    Global_Physical_Variables::P+=pressure_increment;
   }

  strain.close();
 }


 //Repeat for displacement/pressure formulation 
 //enforcing incompressibility
 {
  std::ofstream strain("RESLT_pres_disp_incomp/s_energy.dat");
  std::cout << 
   "Running with pressure/displacement formulation (incompressible) \n";

  // Change output directory
  doc_info.set_directory("RESLT_pres_disp_incomp");
  //Reset doc_info number
  doc_info.number() = 0;
  
  //Set up the problem
  UnstructuredSolidProblem<
  ProjectablePVDElementWithContinuousPressure<
  TPVDElementWithContinuousPressure<2> > > problem;
  
  //Loop over all elements and set incompressibility flag
  {
   const unsigned n_element = problem.mesh_pt()->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     //Cast the element to the equation base of our 2D elastiticy elements
     PVDEquationsWithPressure<2> *cast_el_pt =
      dynamic_cast<PVDEquationsWithPressure<2>*>(
       problem.mesh_pt()->element_pt(e));
     
     //If the cast was successful, it's a bulk element, 
     //so set the incompressibilty flag
     if(cast_el_pt) {cast_el_pt->set_incompressible();}
    }
  }

  //Turn on the incompressibity flag so that elements stay incompressible
  //after refinement
  problem.set_incompressible();

  //Output initial configuration
  problem.doc_solution(doc_info);
  doc_info.number()++;
  
  // Parameter study
  Global_Physical_Variables::P=0.0;
  double pressure_increment=0.1e-2;
  
  unsigned nstep=5;
  for (unsigned istep=0;istep<nstep;istep++)
   {
    // Solve the problem
    problem.newton_solve(1);
    
    double strain_energy = problem.get_strain_energy();
    std::cout << "Strain energy is " << strain_energy << "\n";
    //Output strain energy to file
    strain << Global_Physical_Variables::P << " " << strain_energy << std::endl;
    
    //Output solution
    problem.doc_solution(doc_info);
    doc_info.number()++;

    if(istep==2) {pressure_increment *= -1.0;}        
    // Increase (or decrease) pressure load
    Global_Physical_Variables::P+=pressure_increment;
   }

  strain.close();
 }

} // end main



