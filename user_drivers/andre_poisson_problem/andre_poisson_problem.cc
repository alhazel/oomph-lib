//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2024 Matthias Heil and Andrew Hazel
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
#include <fenv.h> 

//Generic routines
#include "generic.h" 

// The equations
#include "poisson.h"

// The mesh
#include "meshes/triangle_mesh.h"

using namespace std;
using namespace oomph;




//===== start_of_namespace=============================================
/// Namespace for exact solution for Poisson equation with "sharp step" 
//=====================================================================
namespace TanhSolnForPoisson
{

 /// Parameter for steepness of "step"
 double Alpha=5.0;

 /// Parameter for angle Phi of "step"
 double TanPhi=0.0;

 /// Exact solution as a Vector
 void get_exact_u(const Vector<double>& x, Vector<double>& u)
 {
  u[0]=tanh(1.0-Alpha*(TanPhi*x[0]-x[1]));
 }

 /// Source function required to make the solution above an exact solution 
 void get_source(const Vector<double>& x, double& source)
 {
   source = -1.0;
 }
 

 /// Zero function -- used to compute norm of the computed solution by 
 /// computing the norm of the error when compared against this.
 void zero(const Vector<double>& x, Vector<double>& u)
 {
  u[0]=0.0;
 }

} // end of namespace



/// ////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////
/// This geometric object represents the boundary shape that
/// is a homotopy between a triangle and an ellipse.
class HomotopyTriangleEllipse : public GeomObject
  {
  public:
    /// Constructor: 1 Lagrangian coordinate, 2 Eulerian coords. Pass
    /// half axes and homotopy parameter as Data:
    /// \code
    /// Geom_data_pt[0]->value(0) = A
    /// Geom_data_pt[0]->value(1) = B
    /// Geom_data_pt[0]->value(2) = m
    /// \endcode
    HomotopyTriangleEllipse(const Vector<Data*>& geom_data_pt) :
      GeomObject(1, 2)
    {
#ifdef PARANOID
      if (geom_data_pt.size() != 1)
      {
        std::ostringstream error_message;
        error_message << "geom_data_pt should have size 1, not "
                      << geom_data_pt.size() << std::endl;

        if (geom_data_pt[0]->nvalue() != 3)
        {
          error_message << "geom_data_pt[0] should have 3 values, not "
                        << geom_data_pt[0]->nvalue() << std::endl;
        }

        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
     }
#endif
      Geom_data_pt.resize(1);
      Geom_data_pt[0] = geom_data_pt[0];

      // Data has been created externally: Must not clean up
      Must_clean_up = false;
    }


    /// Constructor: 1 Lagrangian coordinate, 2 Eulerian coords. Pass
    /// half axes A and B, and homotopy parameter m; both pinned.
    HomotopyTriangleEllipse(const double& A, const double& B,
			    const double& m) : GeomObject(1, 2)
    {
      // Resize Data for ellipse object:
      Geom_data_pt.resize(1);

      // Create data: Two values, no timedependence, free by default
      Geom_data_pt[0] = new Data(3);

      // I've created the data, I need to clean up
      Must_clean_up = true;

      // Pin the data
      Geom_data_pt[0]->pin(0);
      Geom_data_pt[0]->pin(1);
      Geom_data_pt[0]->pin(2);

      // Set half axes
      Geom_data_pt[0]->set_value(0, A);
      Geom_data_pt[0]->set_value(1, B);
      Geom_data_pt[0]->set_value(2, m);
    }

    /// Broken copy constructor
    HomotopyTriangleEllipse(const HomotopyTriangleEllipse& dummy) = delete;

    /// Broken assignment operator
    void operator=(const HomotopyTriangleEllipse&) = delete;

    /// Destructor:  Clean up if necessary
    ~HomotopyTriangleEllipse()
    {
      // Do I need to clean up?
      if (Must_clean_up)
      {
        delete Geom_data_pt[0];
        Geom_data_pt[0] = 0;
      }
    }

    /// Set horizontal half axis
    void set_A_ellips(const double& a)
    {
      Geom_data_pt[0]->set_value(0, a);
    }

    /// Set vertical half axis
    void set_B_ellips(const double& b)
    {
      Geom_data_pt[0]->set_value(1, b);
    }

    //Set homotopy parameter
    void set_m_ellips(const double &m)
    {
      Geom_data_pt[0]->set_value(2,m);
    }

    /// Access function for horizontal half axis
    double a_ellips()
    {
      return Geom_data_pt[0]->value(0);
    }

    /// Access function for vertical half axis
    double b_ellips()
    {
      return Geom_data_pt[0]->value(1);
    }

    /// Access function for homotopy parameter
    double m_ellips()
    {
      return Geom_data_pt[0]->value(2);
    }

    /// Position Vector at Lagrangian coordinate zeta
    void position(const Vector<double>& zeta, Vector<double>& r) const
    {
      const double pi = 4.0*atan(1.0);
      //Radius of the ellipse
      Vector<double> r_ellipse(2);
      const double a = Geom_data_pt[0]->value(0);
      const double b = Geom_data_pt[0]->value(1);
      const double m = Geom_data_pt[0]->value(2);
      const double theta = zeta[0];
      
      // Position Vector
      r_ellipse[0] = a * cos(theta);
      r_ellipse[1] = b * sin(theta);

      //IMPORTANT NOTE: CENTRE OF ELLIPSE IS NOT THE SAME AS
      //CENTRE OF TRIANGLE. DOES THIS MATTER FOR HOMOTOPY?
      
      //Radius of the exscribed triangle in polar coordinates
      Vector<double> r_triangle(2);
      //Aspect ratio
      double lambda = b/a;
      //Half-length of triangle side
      double L = a*(1.0 + sqrt(1.0 + 3.0*lambda*lambda))/sqrt(3.0);

      //X-coordinate of the triangle centre
      double tri_centre = L*tan(pi/6.0) - a;
      
      //Find angle mod 2*pi/3.0
      double theta_mod = std::fmod(theta,2.0*pi/3.0);
      //Then this is the formula for the radius of the triangle
      //from its centre
      double r_polar = L*tan(pi/6.0)/
	(sin(5.0*pi/6.0 - theta_mod));

      //Now we simply convert to Cartesians
      r_triangle[0] = tri_centre + r_polar * cos(theta);
      r_triangle[1] = r_polar * sin(theta);
      
      //The actual position is the weighted combination of the ellipse and
      //triangle
      for(unsigned i=0;i<2;i++)
	{
	  r[i] = (1.0-m)*r_ellipse[i] + m*r_triangle[i];
	}
	  
    }


    /// Parametrised position on object: r(zeta). Evaluated at
    /// previous timestep. t=0: current time; t>0: previous
    /// timestep.
    void position(const unsigned& t,
                  const Vector<double>& zeta,
                  Vector<double>& r) const
    {
      // If we have done the construction, it's a Steady HomotopyTriangleEllipse,
      // so all time-history values of the position are equal to the position
      if (Must_clean_up)
      {
        position(zeta, r);
        return;
      }

      // Otherwise check that the value of t is within range
#ifdef PARANOID
      if (t > Geom_data_pt[0]->time_stepper_pt()->nprev_values())
      {
        std::ostringstream error_message;
        error_message << "t > nprev_values() " << t << " "
                      << Geom_data_pt[0]->time_stepper_pt()->nprev_values()
                      << std::endl;

        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // Position Vector
      r[0] = Geom_data_pt[0]->value(t, 0) * cos(zeta[0]);
      r[1] = Geom_data_pt[0]->value(t, 1) * sin(zeta[0]);
    }


    /// Derivative of position Vector w.r.t. to coordinates:
    /// \f$ \frac{dR_i}{d \zeta_\alpha}\f$ = drdzeta(alpha,i).
    void dposition(const Vector<double>& zeta,
                   DenseMatrix<double>& drdzeta) const
    {
      exit(1);
      // Components of the single tangent Vector
      drdzeta(0, 0) = -Geom_data_pt[0]->value(0) * sin(zeta[0]);
      drdzeta(0, 1) = Geom_data_pt[0]->value(1) * cos(zeta[0]);
    }


    /// 2nd derivative of position Vector w.r.t. to coordinates:
    /// \f$ \frac{d^2R_i}{d \zeta_\alpha d \zeta_\beta}\f$ =
    /// ddrdzeta(alpha,beta,i).
    /// Evaluated at current time.
    void d2position(const Vector<double>& zeta,
                    RankThreeTensor<double>& ddrdzeta) const
    {
      exit(1);
      // Components of the derivative of the tangent Vector
      ddrdzeta(0, 0, 0) = -Geom_data_pt[0]->value(0) * cos(zeta[0]);
      ddrdzeta(0, 0, 1) = -Geom_data_pt[0]->value(1) * sin(zeta[0]);
    }

    /// Position Vector and 1st and 2nd derivs to coordinates:
    /// \f$ \frac{dR_i}{d \zeta_\alpha}\f$ = drdzeta(alpha,i).
    /// \f$ \frac{d^2R_i}{d \zeta_\alpha d \zeta_\beta}\f$ =
    /// ddrdzeta(alpha,beta,i).
    /// Evaluated at current time.
    void d2position(const Vector<double>& zeta,
                    Vector<double>& r,
                    DenseMatrix<double>& drdzeta,
                    RankThreeTensor<double>& ddrdzeta) const
    {
      exit(1);
      double a = Geom_data_pt[0]->value(0);
      double b = Geom_data_pt[0]->value(1);
      // Position Vector
      r[0] = a * cos(zeta[0]);
      r[1] = b * sin(zeta[0]);

      // Components of the single tangent Vector
      drdzeta(0, 0) = -a * sin(zeta[0]);
      drdzeta(0, 1) = b * cos(zeta[0]);

      // Components of the derivative of the tangent Vector
      ddrdzeta(0, 0, 0) = -a * cos(zeta[0]);
      ddrdzeta(0, 0, 1) = -b * sin(zeta[0]);
    }


    /// How many items of Data does the shape of the object depend on?
    unsigned ngeom_data() const
    {
      return Geom_data_pt.size();
    }

    /// Return pointer to the j-th Data item that the object's
    /// shape depends on
    Data* geom_data_pt(const unsigned& j)
    {
      return Geom_data_pt[j];
    }

  private:
    /// Vector of pointers to Data items that affects the object's shape
    Vector<Data*> Geom_data_pt;

    /// Do I need to clean up?
    bool Must_clean_up;
  };


//==start_of_problem_class============================================
/// Class definition
//====================================================================
template<class ELEMENT>
class UnstructuredPoissonProblem : public virtual Problem
{

public:

 /// Constructor
 UnstructuredPoissonProblem();
    
 /// Destructor
 ~UnstructuredPoissonProblem(){};

 /// Actions before adapt. Empty
 void actions_before_adapt() {}
 
 /// Actions after adapt: 
 /// Setup the problem again -- remember that the mesh has been
 /// completely rebuilt and its element's don't have any
 /// pointers to source fcts etc. yet
 void actions_after_adapt()
  {
   complete_problem_setup();
  }
 
 /// Update after solve (empty)
 void actions_after_newton_solve(){}

 /// Update the problem specs before solve: Re-apply boundary conditons
 void actions_before_newton_solve()
  {
   apply_boundary_conditions();
  }
  
 /// Doc the solution
 void doc_solution(const std::string& comment="");
 

private:

 /// Doc info object for labeling output
 DocInfo Doc_info;

 /// Helper function to apply boundary conditions
 void apply_boundary_conditions();

 /// Helper function to (re-)set boundary condition
 /// and complete the build of  all elements
 void complete_problem_setup();

 /// Pointers to specific mesh
 RefineableTriangleMesh<ELEMENT>* My_mesh_pt;

 /// Trace file to document norm of solution
 ofstream Trace_file;

}; // end_of_problem_class





//==start_constructor=====================================================
/// Constructor
//========================================================================
template<class ELEMENT>
UnstructuredPoissonProblem<ELEMENT>::UnstructuredPoissonProblem()
{          
 // Intrinsic coordinate along GeomObject
 Vector<double> zeta(1);

 // Position vector on GeomObject
 Vector<double> posn(2);
 
 // HomotopyTriangleEllipse defining the outer boundary
 double x_center = 0.0;
 double y_center = 0.0;
 double A = 1.0;
 double B = 0.6;
 double M = 0.1;
 HomotopyTriangleEllipse * outer_boundary_ellipse_pt =
   new HomotopyTriangleEllipse(A,B,M);

 //Let's see the shape
 {
   std::ofstream shape("bound.dat");
   Vector<double> pos(2);
   Vector<double> zeta(1);
   zeta[0] = 0.0;
   unsigned n_incr=900;
   double theta_incr = 8.0*atan(1.0)/n_incr;
   for(unsigned i=0;i<n_incr;++i)
     {
       outer_boundary_ellipse_pt->position(zeta,pos);
       shape << pos[0] << " " << pos[1] << std::endl;
       zeta[0] += theta_incr;
     }
   shape.close();
   exit(1);
   }
       
 // Pointer to the closed curve that defines the outer boundary
 TriangleMeshClosedCurve* closed_curve_pt=0;

 // Build outer boundary as Polygon?
 //---------------------------------
 bool polygon_for_outer_boundary=false;
#ifdef OUTER_POLYGON
 polygon_for_outer_boundary=true;
#endif
 if (polygon_for_outer_boundary)
  { 
   // Number of segments that make up the boundary
   unsigned n_seg = 5; 
   double unit_zeta = 0.5*MathematicalConstants::Pi/double(n_seg);
   
   // The boundary is bounded by two distinct boundaries, each
   // represented by its own polyline
   Vector<TriangleMeshCurveSection*> boundary_polyline_pt(2);
   
   // Vertex coordinates on boundary
   Vector<Vector<double> > bound_coords(n_seg+1);
   
   // First part of the boundary 
   //---------------------------
   for(unsigned ipoint=0; ipoint<n_seg+1;ipoint++)
    {
     // Resize the vector 
     bound_coords[ipoint].resize(2);
     
     // Get the coordinates
     zeta[0]=unit_zeta*double(ipoint);
     outer_boundary_ellipse_pt->position(zeta,posn);
     bound_coords[ipoint][0]=posn[0]+x_center;
     bound_coords[ipoint][1]=posn[1]+y_center;
    }
   
   // Build the 1st boundary polyline
   unsigned boundary_id=0;
   boundary_polyline_pt[0]=new TriangleMeshPolyLine(bound_coords,boundary_id);
   
   // Second part of the boundary
   //----------------------------
   unit_zeta*=3.0;
   for(unsigned ipoint=0; ipoint<n_seg+1;ipoint++)
    {
     // Resize the vector 
     bound_coords[ipoint].resize(2);
     
     // Get the coordinates
     zeta[0]=(unit_zeta*double(ipoint))+0.5*MathematicalConstants::Pi;
     outer_boundary_ellipse_pt->position(zeta,posn);
     bound_coords[ipoint][0]=posn[0]+x_center;
     bound_coords[ipoint][1]=posn[1]+y_center;
    }
   
   // Build the 2nd boundary polyline
   boundary_id=1;
   boundary_polyline_pt[1]=new TriangleMeshPolyLine(bound_coords,boundary_id);
   

   // Create the triangle mesh polygon for outer boundary
   //----------------------------------------------------
   TriangleMeshPolygon *outer_polygon =
    new TriangleMeshPolygon(boundary_polyline_pt);

   // Enable redistribution of polylines
   outer_polygon->
    enable_redistribution_of_segments_between_polylines();

   // Set the pointer
   closed_curve_pt = outer_polygon;

  }
 // Build outer boundary as curvilinear
 //------------------------------------
 else
  {   

   // Provide storage for pointers to the two parts of the curvilinear boundary
   Vector<TriangleMeshCurveSection*> outer_curvilinear_boundary_pt(2);
   
   // First bit
   //----------
   double zeta_start=0.0;
   double zeta_end=MathematicalConstants::Pi;
   unsigned nsegment=5;
   unsigned boundary_id=0;
   outer_curvilinear_boundary_pt[0]=new TriangleMeshCurviLine(
    outer_boundary_ellipse_pt,zeta_start,zeta_end,nsegment,boundary_id);
   
   // Second bit
   //-----------
   zeta_start=MathematicalConstants::Pi;
   zeta_end=2.0*MathematicalConstants::Pi;
   nsegment=8;
   boundary_id=1;
   outer_curvilinear_boundary_pt[1]=new TriangleMeshCurviLine(
    outer_boundary_ellipse_pt,zeta_start,zeta_end,nsegment,boundary_id);
   
   // Combine to curvilinear boundary and define the
   //--------------------------------
   // outer boundary
   //--------------------------------
   closed_curve_pt=
     new TriangleMeshClosedCurve(outer_curvilinear_boundary_pt);
   
  }
 

 
 // Uncomment this as an exercise to observe how a
 // layer of fine elements get left behind near the boundary
 // once the tanh step has swept past: 

 // closed_curve_pt->disable_polyline_refinement();
 // closed_curve_pt->disable_polyline_unrefinement();
 
 // Now build the mesh
 //===================

 // Use the TriangleMeshParameters object for helping on the manage of the
 // TriangleMesh parameters
 TriangleMeshParameters triangle_mesh_parameters(closed_curve_pt);

 // Specify the maximum area element
 double uniform_element_area=0.2;
 triangle_mesh_parameters.element_area() = uniform_element_area;
 
 // Create the mesh
 My_mesh_pt=new 
  RefineableTriangleMesh<ELEMENT>(triangle_mesh_parameters);
 
 // Store as the problem's one and only mesh
 Problem::mesh_pt()=My_mesh_pt;

 // Set error estimator for bulk mesh
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 My_mesh_pt->spatial_error_estimator_pt()=error_estimator_pt;

 // Set element size limits
 My_mesh_pt->max_element_size()=0.2;
 My_mesh_pt->min_element_size()=0.002; 
 
 // Set boundary condition and complete the build of all elements
 complete_problem_setup();

 // Open trace file
 char filename[100];
 sprintf(filename,"RESLT/trace.dat");
 Trace_file.open(filename);

 // Setup equation numbering scheme
 oomph_info <<"Number of equations: " 
            << this->assign_eqn_numbers() << std::endl;
 
} // end_of_constructor




//==start_of_complete======================================================
 /// Set boundary condition exactly, and complete the build of 
 /// all elements
//========================================================================
template<class ELEMENT>
void UnstructuredPoissonProblem<ELEMENT>::complete_problem_setup()
{   

 // Set the boundary conditions for problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here. 
 unsigned nbound=My_mesh_pt->nboundary();
 for(unsigned ibound=0;ibound<nbound;ibound++)
  {
   unsigned num_nod=My_mesh_pt->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Get node
     Node* nod_pt=My_mesh_pt->boundary_node_pt(ibound,inod);
     
     // Pin one-and-only unknown value
     nod_pt->pin(0);
    }   
  } // end loop over boundaries
 
 
 // Complete the build of all elements so they are fully functional
 unsigned n_element = My_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(My_mesh_pt->element_pt(e));
   
   //Set the source function pointer
   el_pt->source_fct_pt() = &TanhSolnForPoisson::get_source;
  }
 
 // Re-apply Dirichlet boundary conditions (projection ignores
 // boundary conditions!)
 apply_boundary_conditions();
}




//==start_of_apply_bc=====================================================
 /// Helper function to apply boundary conditions
//========================================================================
template<class ELEMENT>
void UnstructuredPoissonProblem<ELEMENT>::apply_boundary_conditions()
{
 
 // Loop over all boundary nodes
 unsigned nbound=this->My_mesh_pt->nboundary();
 for(unsigned ibound=0;ibound<nbound;ibound++)
  {
   unsigned num_nod=this->My_mesh_pt->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Get node
     Node* nod_pt=this->My_mesh_pt->boundary_node_pt(ibound,inod);
     
     //Boundaries are at zero
     nod_pt->set_value(0,0.0);
    }
  } 

} // end set bc


//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void UnstructuredPoissonProblem<ELEMENT>::doc_solution(const 
                                                       std::string& comment)
{ 
 ofstream some_file;
 char filename[100];
 
 // Number of plot points
 unsigned npts;
 npts=5; 
 
 sprintf(filename,"RESLT/soln%i.dat",Doc_info.number());
 some_file.open(filename);
 this->My_mesh_pt->output(some_file,npts); 
 some_file << "TEXT X = 22, Y = 92, CS=FRAME T = \"" 
           << comment << "\"\n";
 some_file.close();
 
 // Output exact solution 
 //----------------------
 sprintf(filename,"RESLT/exact_soln%i.dat",Doc_info.number());
 some_file.open(filename);
 My_mesh_pt->output_fct(some_file,npts,TanhSolnForPoisson::get_exact_u); 
 some_file.close();
 
 // Output boundaries
 //------------------
 sprintf(filename,"RESLT/boundaries%i.dat",Doc_info.number());
 some_file.open(filename);
 My_mesh_pt->output_boundaries(some_file);
 some_file.close();


 // Doc error and return of the square of the L2 error
 //---------------------------------------------------
 double error,norm,dummy_error,zero_norm;
 sprintf(filename,"RESLT/error%i.dat",Doc_info.number());
 some_file.open(filename);
 My_mesh_pt->compute_error(some_file,TanhSolnForPoisson::get_exact_u,
                           error,norm); 
 
 My_mesh_pt->compute_error(some_file,TanhSolnForPoisson::zero,
                           dummy_error,zero_norm); 
 some_file.close();

 // Doc L2 error and norm of solution
 oomph_info << "\nNorm of error   : " << sqrt(error) << std::endl; 
 oomph_info << "Norm of exact solution: " << sqrt(norm) << std::endl;
 oomph_info << "Norm of computed solution: " << sqrt(dummy_error) << std::endl;
 Trace_file << sqrt(norm) << " " << sqrt(dummy_error) << std::endl;

 // Increment the doc_info number
 Doc_info.number()++;

} // end of doc


//=======start_of_main========================================
/// Driver code for demo of inline triangle mesh generation
//============================================================
int main(int argc, char **argv)
{
 // Store command line arguments
 CommandLineArgs::setup(argc,argv);
 
 // Define possible command line arguments and parse the ones that
 // were actually specified
 
 // Validation?
 CommandLineArgs::specify_command_line_flag("--validation");

 // Parse command line
 CommandLineArgs::parse_and_assign(); 
 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();
 
 // Create problem
 UnstructuredPoissonProblem<ProjectablePoissonElement<TPoissonElement<2,3> > >
  problem;
 
//  // Solve, adapt and doc manually
//  unsigned nadapt=4;
//  for (unsigned i=0;i<nadapt;i++)
//   {
//    problem.newton_solve();   
//    std::stringstream comment_stream;
//    comment_stream << "Solution after " << i << " manual adaptations";
//    problem.doc_solution(comment_stream.str());
//    if (i!=(nadapt-1))  problem.adapt();
//   }
 
 
 // Doc the initial mesh
 //=====================
 {
   std::stringstream comment_stream;
   comment_stream << "Initial mesh ";
   problem.doc_solution(comment_stream.str());
 } 

 // Solve with spatial adaptation
 //==============================
 unsigned max_adapt=3;
 if (CommandLineArgs::command_line_flag_has_been_set("--validation"))
   {
     max_adapt=1;
   }
 problem.newton_solve(max_adapt);
 
 // Doc the solution
 //=================
 std::stringstream comment_stream;
 comment_stream << "Solution for tan(phi) = " << TanhSolnForPoisson::TanPhi; 
 problem.doc_solution(comment_stream.str());
 
} //End of main
