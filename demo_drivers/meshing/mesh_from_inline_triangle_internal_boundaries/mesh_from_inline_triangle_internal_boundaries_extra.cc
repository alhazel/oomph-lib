//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2025 Matthias Heil and Andrew Hazel
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
// Driver code for a simple test poisson problem using a mesh
// generated from an input file generated by the triangle mesh generator
// Triangle.

//Generic routines
#include "generic.h"

// Poisson equations
#include "poisson.h"

// The mesh
#include "meshes/triangle_mesh.h"

using namespace std;

using namespace oomph;

//====================================================================
/// Namespace for exact solution for Poisson equation with sharp step 
//====================================================================
namespace TanhSolnForPoisson
{

  /// Parameter for steepness of step
  double Alpha;

  /// Parameter for angle of step
  double Beta;


  /// Exact solution as a Vector
  void get_exact_u(const Vector<double>& x, Vector<double>& u)
  {
    u[0]=tanh(1.0-Alpha*(Beta*x[0]-x[1]));
  }


  /// Exact solution as a scalar
  void get_exact_u(const Vector<double>& x, double& u)
  {
    u=tanh(1.0-Alpha*(Beta*x[0]-x[1]));
  }


  /// Source function to make it an exact solution
  void get_source(const Vector<double>& x, double& source)
  {
    source = 2.0*tanh(-1.0+Alpha*(Beta*x[0]-x[1]))*
        (1.0-pow(tanh(-1.0+Alpha*(Beta*x[0]-x[1])),2.0))*
        Alpha*Alpha*Beta*Beta+2.0*tanh(-1.0+Alpha*(Beta*x[0]-x[1]))*
        (1.0-pow(tanh(-1.0+Alpha*(Beta*x[0]-x[1])),2.0))*Alpha*Alpha;
  }

}







//====================================================================
/// Poisson problem.
//====================================================================

// Poisson problem
template<class ELEMENT> 
class PoissonProblem : public Problem
{

public:


  /// Constructor: Pass pointer to source function and names of
  /// two triangle input files
  PoissonProblem(PoissonEquations<2>::PoissonSourceFctPt source_fct_pt);

  /// Destructor (empty)
  ~PoissonProblem(){}

  /// Update the problem specs before solve: (Re)set boundary conditions
  void actions_before_newton_solve()
  {
    //Loop over the boundaries
    unsigned num_bound = mesh_pt()->nboundary();
    for(unsigned ibound=0;ibound<num_bound;ibound++)
      {
        // Loop over the nodes on boundary
        unsigned num_nod=mesh_pt()->nboundary_node(ibound);
        for (unsigned inod=0;inod<num_nod;inod++)
          {
            Node* nod_pt=mesh_pt()->boundary_node_pt(ibound,inod);
            double u;
            Vector<double> x(2);
            x[0]=nod_pt->x(0);
            x[1]=nod_pt->x(1);
            TanhSolnForPoisson::get_exact_u(x,u);
            nod_pt->set_value(0,u);
          }
      }
  }

  /// Update the problem specs before solve (empty)
  void actions_after_newton_solve()
  {}


  /// Access function for the specific mesh
  TriangleMesh<ELEMENT>* mesh_pt()
      {
    return dynamic_cast<TriangleMesh<ELEMENT>*>(Problem::mesh_pt());
      }

  /// Doc the solution
  void doc_solution(DocInfo& doc_info);

private:

  /// Pointer to source function
  PoissonEquations<2>::PoissonSourceFctPt Source_fct_pt;

  /// Trace file to document norm of solution
  ofstream Trace_file;

};



//========================================================================
/// Constructor for Poisson problem
//========================================================================

template<class ELEMENT>
PoissonProblem<ELEMENT>::
PoissonProblem(PoissonEquations<2>::PoissonSourceFctPt source_fct_pt)
: Source_fct_pt(source_fct_pt)
  {

  // Setup parameters for exact tanh solution

  // Steepness of step
  TanhSolnForPoisson::Alpha=1.0;

  // Orientation of step
  TanhSolnForPoisson::Beta=1.4;

  unsigned boundary_id;

  // *********************************************************************
  // Begin - Outer boundary

  // Create storage for the vertices coordinates that define each boundary
  const unsigned num_vertices_b0 = 4;
  Vector<Vector <double> > vertices(num_vertices_b0);

  vertices[0].resize(2);
  vertices[1].resize(2);
  vertices[2].resize(2);
  vertices[3].resize(2);

  // Building the outer boundary (vertices)
  vertices[0][0] = 0;
  vertices[0][1] = 0;

  vertices[1][0] = 0;
  vertices[1][1] = 1;

  vertices[2][0] = 0;
  vertices[2][1] = 3;

  vertices[3][0] = 3;
  vertices[3][1] = 3;

  // The outer boundary is represented by two TriangleMeshPolyLine objects,
  // this is the first one
  boundary_id = 0;
  TriangleMeshPolyLine *boundary0_pt =
      new TriangleMeshPolyLine(vertices, boundary_id);

  // Create storage for the vertices coordinates that define each boundary
  const unsigned num_vertices_b1 = 5;
  vertices.resize(num_vertices_b1);

  vertices[0].resize(2);
  vertices[1].resize(2);
  vertices[2].resize(2);
  vertices[3].resize(2);
  vertices[4].resize(2);

  // More vertices for the outer boundary
  vertices[0][0] = 3;
  vertices[0][1] = 3;

  vertices[1][0] = 3;
  vertices[1][1] = 0.5;

  vertices[2][0] = 3;
  vertices[2][1] = 0;

  vertices[3][0] = 1;
  vertices[3][1] = 0;

  vertices[4][0] = 0;
  vertices[4][1] = 0;

  // The outer boundary is represented by two TriangleMeshPolyLine objects,
  // this is the second one
  boundary_id = 1;
  TriangleMeshPolyLine *boundary1_pt =
      new TriangleMeshPolyLine(vertices, boundary_id);

  // A vector for storing the outer boundary representation
  Vector<TriangleMeshCurveSection*> outer_boundary_polyLines_pt(2);

  outer_boundary_polyLines_pt[0] = boundary0_pt;
  outer_boundary_polyLines_pt[1] = boundary1_pt;

  TriangleMeshClosedCurve *outer_boundary_pt =
      new TriangleMeshPolygon(outer_boundary_polyLines_pt);
  // End - Outer boundary
  // *********************************************************************

  // *********************************************************************
  // Begin - Internal closed boundaries
  // We define no internal closed boundaries
  unsigned n_internal_closed_boundaries = 0;
  Vector<TriangleMeshClosedCurve *>
  inner_boundaries_pt(n_internal_closed_boundaries);
  // End - Internal closed boundaries
  // *********************************************************************

  // *********************************************************************
  // Internal open boundaries
  // Total number of open curves in the domain
  unsigned n_open_curves = 5;

  // We want internal open curves
  Vector<TriangleMeshOpenCurve *> inner_open_boundaries_pt(n_open_curves);

  // *********************************************************************
  // We start by creating the internal boundaries
  // The boundary 2 is defined by its two vertices
  // Open curve 1
  vertices.resize(2);
  vertices[0].resize(2);
  vertices[1].resize(2);

  vertices[0][0] = 0.5;
  vertices[0][1] = 2.0;

  vertices[1][0] = 0.5;
  vertices[1][1] = 2.5;

  boundary_id = 2;
  TriangleMeshPolyLine *boundary2_pt =
    new TriangleMeshPolyLine(vertices, boundary_id);

  // Each internal open curve is defined by a vector of TriangleMeshCurveSection,
  // on this example we only need one curve section for each internal boundary
  Vector<TriangleMeshCurveSection *> internal_curve_section1_pt(1);

  internal_curve_section1_pt[0] = boundary2_pt;

  // The open curve that define this boundary is composed of just one
  // curve section
  inner_open_boundaries_pt[0] =
      new TriangleMeshOpenCurve(internal_curve_section1_pt);

  // *********************************************************************
  // We define the curved boundary as a TriangleMeshCurviline
  // Open curve 2
  double x_centre = 2.0;
  double y_centre = 2.0;
  double r_circle = 0.5;

  Circle * boundary_circle1_pt = new Circle(x_centre, y_centre, r_circle);

  // Number of segments used for representing the curve boundary
  unsigned n_segments = 20;

  // The intrinsic coordinates for the beginning and end of the curve
  double s_start = 0.0;
  double s_end = MathematicalConstants::Pi;

  boundary_id = 3;
  TriangleMeshCurviLine *boundary3_pt =
      new TriangleMeshCurviLine(boundary_circle1_pt,
          s_start,
          s_end,
          n_segments,
          boundary_id);

  Vector<TriangleMeshCurveSection *> internal_curve_section2_pt(1);
  internal_curve_section2_pt[0] = boundary3_pt;

  // The open curve that define this boundary is composed of just one
  // curve section
  inner_open_boundaries_pt[1] =
      new TriangleMeshOpenCurve(internal_curve_section2_pt);

  // *********************************************************************
  // We define the curved boundary as a TriangleMeshCurviline
  // Open curve 3
  x_centre = 0.0;
  y_centre = 0.0;
  r_circle = 1.0;

  Circle * boundary_circle2_pt = new Circle(x_centre, y_centre, r_circle);

  // Number of segments used for representing the curve boundary
  n_segments = 20;

  // The intrinsic coordinates for the beginning and end of the curve
  s_start = 0.0;
  s_end = MathematicalConstants::Pi / 2.0;

  // State the vertex number for connection on the destination
  // boundaries
  unsigned vertex_to_connect_initial = 3;
  unsigned vertex_to_connect_final = 1;

  boundary_id = 4;
  TriangleMeshCurviLine *boundary4_pt =
      new TriangleMeshCurviLine(boundary_circle2_pt,
          s_start,
          s_end,
          n_segments,
          boundary_id);

  // Do the connection with the destination boundary, in this case
  // the connection is done with the outer boundary
  boundary4_pt->connect_initial_vertex_to_polyline(
                                boundary1_pt,
                                vertex_to_connect_initial);

  // Do the connection with the destination boundary, in this case
  // the connection is done with the outer boundary
  boundary4_pt->connect_final_vertex_to_polyline(
                                boundary0_pt,
                                vertex_to_connect_final);

  Vector<TriangleMeshCurveSection *> internal_curve_section3_pt(1);
  internal_curve_section3_pt[0] = boundary4_pt;

  // The open curve that define this boundary is composed of just one
  // curve section
  inner_open_boundaries_pt[2] =
      new TriangleMeshOpenCurve(internal_curve_section3_pt);

  // *********************************************************************
  // This boundary is connected to the outer boundary on the initial end
  // and to an internal boundary on the final end
  // Open curve 4
  vertices.resize(3);
  vertices[0].resize(2);
  vertices[1].resize(2);
  vertices[2].resize(2);

  // We need to specify the vertices for the boundary on the first end
  // The connection method performs a checking on the matching of the
  // vertices used for connection
  vertices[0][0] = 3.0;
  vertices[0][1] = 0.5;

  // These values are necessary when we subsequently perform the connection
  // of the last boundary
  vertices[1][0] = (3.0 + sqrt(3.0)) * 0.5;
  vertices[1][1] = 0.5;

  // We need to specify the vertices for the boundary on the last end
  // The connection method performs a checking on the matching of the
  // vertices used for connection
  vertices[2][0] = sqrt(3.0) * 0.5;
  vertices[2][1] = 0.5;

  // State the vertex number for connection with the destination
  // boundary
  vertex_to_connect_initial = 1;

  // State the s value for connection with the destination boundary
  double s_connection_final = atan2(0.5, sqrt(3.0) * 0.5);

  boundary_id = 5;
  TriangleMeshPolyLine *boundary5_pt =
    new TriangleMeshPolyLine(vertices, boundary_id);

  // Do the connection with the destination boundary, in this case
  // the connection is done with the outer boundary
  boundary5_pt->connect_initial_vertex_to_polyline(
                                boundary1_pt,
                                vertex_to_connect_initial);

  // Do the connection with the destination boundary, in this case
  // the connection is done with an internal boundary
  boundary5_pt->connect_final_vertex_to_curviline(
                                boundary4_pt,
                                s_connection_final);

  Vector<TriangleMeshCurveSection *> internal_curve_section4_pt(1);
  internal_curve_section4_pt[0] = boundary5_pt;

  // The open curve that define this boundary is composed of just one
  // curve section
  inner_open_boundaries_pt[3] =
      new TriangleMeshOpenCurve(internal_curve_section4_pt);

  // *********************************************************************
  // We define the curved boundary as a TriangleMeshCurviline
  // Open curve 5
  x_centre = 1.5;
  y_centre = 0.0;
  r_circle = 1.0;

  Circle * boundary_circle3_pt = new Circle(x_centre, y_centre, r_circle);

  // Number of segments used for representing the curve boundary
  n_segments = 20;

  // These numbers can be easily obtained by computing the
  // intersection of the circle (circle3) and the other boundaries
  // (in this case, boundaries 4 and 5)
  s_start = atan2(0.5, (3.0 + sqrt(3.0)) * 0.5 - 1.5);
  s_end = atan2(0.25*sqrt(7.0), -0.75);

  // State the vertex number for connection on the destination
  // boundaries
  vertex_to_connect_initial = 1;

  // State the s value for connection with the destination boundary
  s_connection_final = atan2(0.25*sqrt(7.0), 0.75);

  boundary_id = 6;
  TriangleMeshCurviLine *boundary6_pt =
    new TriangleMeshCurviLine(boundary_circle3_pt,
          s_start,
          s_end,
          n_segments,
          boundary_id);

  // Do the connection with the destination boundary, in this case
  // the connection is done with an internal boundary
  boundary6_pt->connect_initial_vertex_to_polyline(
                                boundary5_pt,
                                vertex_to_connect_initial);

  // Do the connection with the destination boundary, in this case
  // the connection is done with an internal boundary
  boundary6_pt->connect_final_vertex_to_curviline(
                                boundary4_pt,
                                s_connection_final);

  Vector<TriangleMeshCurveSection *> internal_curve_section5_pt(1);
  internal_curve_section5_pt[0] = boundary6_pt;

  // The open curve that define this boundary is composed of just one
  // curve section
  inner_open_boundaries_pt[4] =
      new TriangleMeshOpenCurve(internal_curve_section5_pt);

  // *********************************************************************
  // *********************************************************************

  //Create mesh

  // Use the TriangleMeshParameters object for helping on the manage of the
  // TriangleMesh parameters. The only parameter that needs to take is the
  // outer boundary.
  TriangleMeshParameters triangle_mesh_parameters(outer_boundary_pt);

  // Specify the internal closed boundaries
  triangle_mesh_parameters.internal_closed_curve_pt() = inner_boundaries_pt;

  // Specify the internal open boundaries
  triangle_mesh_parameters.internal_open_curves_pt() = inner_open_boundaries_pt;

  // One can define an specific maximum element area
  double element_area = 0.2;
  triangle_mesh_parameters.element_area() = element_area;

  // Adding a hole on the domain
  Vector<Vector <double> > additional_hole(1);

  // Define the coordinates on the domain
  additional_hole[0].resize(2);
  additional_hole[0][0] = 1.5;
  additional_hole[0][1] = 0.75;

  // Pass information about the additional holes coordinates
  triangle_mesh_parameters.extra_holes_coordinates() = additional_hole;

  // Adding a region on the domain
  Vector<double> region(2);

  // Define the coordinates of the regions on the domain
  region[0] = 0.5;
  region[1] = 0.5;

  // Pass information about the defined regions
  triangle_mesh_parameters.add_region_coordinates(1, region);
  
  // Specify different target area of region 1:
  triangle_mesh_parameters.set_target_area_for_region(1,0.01);

  // Pass the TriangleMeshParameters object to the TriangleMesh one
  Problem::mesh_pt() = new TriangleMesh<ELEMENT>(triangle_mesh_parameters);

  // Set the boundary conditions for this problem: All nodes are
  // free by default -- just pin the ones that have Dirichlet conditions
  // here.
  unsigned num_bound = mesh_pt()->nboundary();
  for(unsigned ibound=0;ibound<num_bound;ibound++)
    {
      unsigned num_nod= mesh_pt()->nboundary_node(ibound);
      for (unsigned inod=0;inod<num_nod;inod++)
        {
          mesh_pt()->boundary_node_pt(ibound,inod)->pin(0);
        }
    }

  // Complete the build of all elements so they are fully functional

  //Find number of elements in mesh
  unsigned n_element = mesh_pt()->nelement();

  // Loop over the elements to set up element-specific
  // things that cannot be handled by constructor
  for(unsigned i=0;i<n_element;i++)
    {
      // Upcast from GeneralElement to the present element
      ELEMENT* el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));

      //Set the source function pointer
      el_pt->source_fct_pt() = Source_fct_pt;
    }

  // Open trace file
  char filename[100];
  sprintf(filename,"RESLT/trace.dat");
  Trace_file.open(filename);

  // Setup equation numbering scheme
  cout <<"Number of equations: " << assign_eqn_numbers() << std::endl;

  }



//========================================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void PoissonProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{ 

  ofstream some_file;
  char filename[100];

  // Number of plot points
  unsigned npts;
  npts=10;


  // Output boundaries
  //------------------
  sprintf(filename,"%s/boundaries.dat",doc_info.directory().c_str());
  some_file.open(filename);
  mesh_pt()->output_boundaries(some_file);
  some_file.close();

  // Output solution
  //----------------
  sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
      doc_info.number());
  some_file.open(filename);
  mesh_pt()->output(some_file,npts);
  some_file.close();


  // Output exact solution
  //----------------------
  sprintf(filename,"%s/exact_soln%i.dat",doc_info.directory().c_str(),
      doc_info.number());
  some_file.open(filename);
  mesh_pt()->output_fct(some_file,npts,TanhSolnForPoisson::get_exact_u);
  some_file.close();


  // Output regions
  //---------------
  {
   unsigned nr=mesh_pt()->nregion();
   for (unsigned r=0;r<nr;r++)
    {
     ofstream some_file;
     sprintf(filename,"%s/soln_in_region%i_%i.dat",doc_info.directory().c_str(),
             r,doc_info.number());
     some_file.open(filename);
     unsigned nel=mesh_pt()->nregion_element(r);
     for(unsigned e=0;e<nel;e++)
      {
       mesh_pt()->region_element_pt(r,e)->output(some_file,npts);
      }
     some_file.close();
    }
         }

  // Doc error
  //----------
  double error,norm;
  sprintf(filename,"%s/error%i.dat",doc_info.directory().c_str(),
      doc_info.number());
  some_file.open(filename);
  mesh_pt()->compute_error(some_file,TanhSolnForPoisson::get_exact_u,
      error,norm);
  some_file.close();
  cout << "error: " << sqrt(error) << std::endl;
  cout << "norm : " << sqrt(norm) << std::endl << std::endl;

  Trace_file << sqrt(norm) << " " << sqrt(error) << std::endl;

}





//========================================================================
/// Demonstrate how to solve Poisson problem
//========================================================================
int main(int argc, char* argv[])
{

  // Store command line arguments
  CommandLineArgs::setup(argc,argv);

  // Label for output
  DocInfo doc_info;

  // Output directory
  doc_info.set_directory("RESLT");


  // Do the problem with cubic elements
  //-----------------------------------
  {
    cout << std::endl << "Cubic elements" << std::endl;
    cout <<         "==============" << std::endl << std::endl;

    //Set up the problem
    PoissonProblem<TPoissonElement<2,4> >
    problem(&TanhSolnForPoisson::get_source);

    // Solve the problem
    problem.newton_solve();

    //Output solution
    problem.doc_solution(doc_info);

    //Increment counter for solutions
    doc_info.number()++;
  }


  // Do the problem with quadratic elements
  //---------------------------------------
  {
    cout << std::endl  << "Quadratic elements" << std::endl;
    cout <<               "===================" << std::endl << std::endl;

    //Set up the problem
    PoissonProblem<TPoissonElement<2,3> >
    problem(&TanhSolnForPoisson::get_source);

    // Solve the problem
    problem.newton_solve();

    //Output solution
    problem.doc_solution(doc_info);

    //Increment counter for solutions
    doc_info.number()++;
  }



  // Do the problem with linear elements
  //------------------------------------
  {
    cout << std::endl << "Linear elements" << std::endl;
    cout <<              "===============" << std::endl << std::endl;

    //Set up the problem
    PoissonProblem<TPoissonElement<2,2> >
    problem(&TanhSolnForPoisson::get_source);

    // Solve the problem
    problem.newton_solve();

    //Output solution
    problem.doc_solution(doc_info);

    //Increment counter for solutions
    doc_info.number()++;
  }

}



