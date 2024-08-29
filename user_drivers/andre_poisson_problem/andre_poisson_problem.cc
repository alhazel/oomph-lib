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


/// Namespace for source function (constant -1) 
//=====================================================================
namespace ConstantSourceFunction
{
 /// Source function required to make the solution above an exact solution 
  void get_source(const Vector<double>& x, double& source)
  {
    source = -1.0;
  }
} // end of namespace




/// ////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////
/// This geometric object represents the boundary shape that
/// is a homotopy between a triangle and an ellipse by making
/// a transform between the boundaries. The difficulty with
/// this approach is that it does not lead to smooth boundaries.
/////////////////////////////////////////////////////////////////
class HomotopyTriangleEllipse : public GeomObject
  {
  public:
    /// Constructor: 1 Lagrangian coordinate, 2 Eulerian coords. Pass
    /// half axes A and B, and homotopy parameter m; both pinned.
    HomotopyTriangleEllipse(const double& A, const double& B,
			    const double& m) : GeomObject(1, 2)
    {
      // Resize Data for ellipse object:
      Geom_data_pt.resize(1);

      // Create data: Two values, no timedependence, free by default
      Geom_data_pt[0] = new Data(3);

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
      //Delete the allocated data
        delete Geom_data_pt[0];
        Geom_data_pt[0] = 0;
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
      //No time dependence in boundary
      position(zeta, r);
      return;
    }

    /// Derivative of position Vector w.r.t. to coordinates:
    /// \f$ \frac{dR_i}{d \zeta_\alpha}\f$ = drdzeta(alpha,i).
    void dposition(const Vector<double>& zeta,
                   DenseMatrix<double>& drdzeta) const
    {
      //Not implemented so break
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
      //Not implemented so break
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
      //Not implemented so break
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
  };

/// ////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////
/// This geometric object represents the boundary shape that
/// is a homotopy between a triangle and an ellipse based on
/// the zero contour of a function. The formula was provided
/// by Andre von Bories Lopes.
/// ////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////
class HomotopyTriangleEllipseCartesian : public GeomObject
  {
  public:
    /// Constructor: 1 Lagrangian coordinate, 2 Eulerian coords. Pass
    /// half axes A and B, and homotopy parameter m; both pinned.
    HomotopyTriangleEllipseCartesian(const double& A, const double& B,
			    const double& m) : GeomObject(1, 2)
    {
      // Resize Data for ellipse object:
      Geom_data_pt.resize(1);

      // Create data: Two values, no timedependence, free by default
      Geom_data_pt[0] = new Data(3);

      // Pin the data
      Geom_data_pt[0]->pin(0);
      Geom_data_pt[0]->pin(1);
      Geom_data_pt[0]->pin(2);

      // Set half axes
      Geom_data_pt[0]->set_value(0, A);
      Geom_data_pt[0]->set_value(1, B);
      Geom_data_pt[0]->set_value(2, m);


      //Now we find and store the x_limits of the shape
      X_limit.resize(2);
      X_limit[0] = this->find_x_limit(-1);
      X_limit[1] = this->find_x_limit(1);
    }

    /// Broken copy constructor
    HomotopyTriangleEllipseCartesian(
     const HomotopyTriangleEllipseCartesian& dummy) = delete;

    /// Broken assignment operator
    void operator=(const HomotopyTriangleEllipseCartesian&) = delete;

    /// Destructor:  Clean up if necessary
    ~HomotopyTriangleEllipseCartesian()
    {
      //Clean up the allocated storage
        delete Geom_data_pt[0];
        Geom_data_pt[0] = 0;
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

    /// Return the square of the height y as a function of x
    // using the formula given by Andre. Note that this could be negative
    // which is why we use it in the root-finding to calculate the x limits.
    double y_square_value(const double &x) const
    {
      //Read out the geometric information
      const double a = Geom_data_pt[0]->value(0);
      const double b = Geom_data_pt[0]->value(1);
      const double m = Geom_data_pt[0]->value(2);
      const double lambda = b/a;

      //Formula provided by Andre
      return ((((2.0 - 2.0*pow(x, 2.0))*pow(lambda, 2.0) 
		- m*(2.0*pow(x, 2.0) + 2.0*(lambda*x) + 2.0*pow(lambda, 2.0) 
		     + 2.0*(x*pow(lambda, 3.0))))*sqrt(3.0*pow(lambda, 2.0)
						       + 1.0) 
	       + m*(3.0*pow(lambda, 5.0) + 3.0*(x*pow(lambda, 4.0)) 
		    + (pow(x, 2.0) + 4.0)*pow(lambda, 3.0) 
		    + (pow(x, 3.0) + 2.0*pow(x, 2.0) + 4.0*x - 2.0)*
		    pow(lambda, 2.0) + lambda*(pow(x, 2.0) + 1.0) 
		    + pow(x, 3.0) + x) + (2.0 - 2.0*pow(x, 2.0))*
	       pow(lambda, 2.0)) 
	      / ((2.0 - 2.0*m)*sqrt(3.0*pow(lambda, 2.0) + 1) 
		 + m*(3.0*pow(lambda, 3.0) + 3.0*(x*pow(lambda, 2.0)) + 3.0*lambda
		      + 3.0*x - 2.0) + 2.0));
    }


    //Return the sign of the value of the y_square value
    //to within a tolerance default 1.0e-14
    int sign_y_square(const double &x,
		      const double &tolerance=1.0e-14) const
    {
      //Get y
      double y = y_square_value(x);
      //If it's close enough to zero return zero
      if(std::abs(y) < tolerance) {return 0;}
      //Otherwise just return the sign
      else
	{
	  if(y > 0) {return 1;}
	  else {return -1;}
	}
    }


    //Find the minimum or maximum value of x using bisection
    //based on sign change of the y_square function
    double find_x_limit(const int &sign_incr) const
    {
      //Increment in x
      double x_incr = sign_incr*0.1;
      //Initial value
      double x = 0.0;
      //Storage for previous sign
      int old_sign = sign_y_square(x);

      int counter=0;
      int max_count=200;
      
      //Bisection loop
      do
	{
	  //Take a step
	  x += x_incr;
	  //What's the new sign
	  int sign = sign_y_square(x);
	  //If we're close enough to zero return the value of x
	  if(sign == 0) {return x;}
	  //Otherwise check for a sign change
	  if(sign != old_sign)
	    {
	      //Reverse direction and half the increment if there
	      //is a sign change
	      x_incr *= -0.5;
	    }
	  //Update the old sign
	  old_sign = sign;
	  ++counter;
	  if(counter > max_count)
	    {
	      throw
		OomphLibError("Bisection method failed.",
			      OOMPH_CURRENT_FUNCTION,
			      OOMPH_EXCEPTION_LOCATION);
	    }
	    }
      while(true); //Can loop forever because the loop will break
                   //when the root is found.
      
    }

    

    
    /// Position Vector at Lagrangian coordinate zeta
    void position(const Vector<double>& zeta, Vector<double>& r) const
    {
      const double pi = 4.0*atan(1.0);
      const double theta = zeta[0];
      
      //The upper and lower values of x have been found, so let's create a
      //simple linear mapping between theta and these values:
      double x_range = X_limit[1] - X_limit[0];

      //Convert theta to x
      double x = 0.0;
      //Upper half-plane
      if(theta < pi)
	{
	  x = X_limit[1] - (theta/pi)*x_range;
	}
      //Lower half plane
      else
	{
	  x = X_limit[1] + (theta/pi - 2.0)*x_range;
	}

      //Now we get y
      double y = y_square_value(x);
      //Deal with finite precision
      if(std::abs(y) < 1.0e-14) {y = 0.0;}

      //This shoudln't happen
      if(y < 0.0)
	{
	  std::cout << "Y is " << y << "\n";
	  throw OomphLibError("Negative value in argument to square-root\n",
			      OOMPH_CURRENT_FUNCTION,
			      OOMPH_EXCEPTION_LOCATION);
	}
      else
	{
	  //Upper half
	  if(theta < pi)
	    {
	      y = sqrt(y);
	    }
	  //Lower half
	  else
	    {
	      y = -sqrt(y);
	    }
	}

      //Fill in the return values
      r[0] = x; r[1] = y;
    }


    /// Parametrised position on object: r(zeta). Evaluated at
    /// previous timestep. t=0: current time; t>0: previous
    /// timestep.
    void position(const unsigned& t,
                  const Vector<double>& zeta,
                  Vector<double>& r) const
    {
      //Boundary is all steady
      position(zeta, r);
      return;
    }
    

    /// Derivative of position Vector w.r.t. to coordinates:
    /// \f$ \frac{dR_i}{d \zeta_\alpha}\f$ = drdzeta(alpha,i).
    void dposition(const Vector<double>& zeta,
                   DenseMatrix<double>& drdzeta) const
    {
      //Not implemented so die
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
      //Not implemented so die
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
      //Not implemented so die
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

    /// Storage for x_limits
    Vector<double> X_limit;

};


namespace oomph
{

  ///Overload Poisson Elements to allow specific output calculations
  class MyPoissonElement:
    public ProjectablePoissonElement<TPoissonElement<2,3> >
  {
  public:
    
    ///Constructor
    MyPoissonElement() : ProjectablePoissonElement<TPoissonElement<2,3> >()
    {
    }
    
    ///Destructor
    ~MyPoissonElement()
    {
    }
    
    
    //Calculate the integrated value of u over the element
    double integrate_u()
    {
      // Initialise
      double integrated_u=0.0;
      
      // Vector of local coordinates
      Vector<double> s(2);
      
      // Set the value of n_intpt
      unsigned n_intpt = integral_pt()->nweight();
      
      // Loop over the integration points
      for (unsigned ipt = 0; ipt < n_intpt; ipt++)
	{
	  // Assign values of s
	  for (unsigned i = 0; i < 2; i++)
	    {
	      s[i] = integral_pt()->knot(ipt, i);
	    }
	  
	  // Get the integral weight
	  double w = integral_pt()->weight(ipt);
	  
	  // Get jacobian of mapping
	  double J = J_eulerian(s);
	  
	  // Get FE function value
	  double u_fe = interpolated_u_poisson(s);
	  
	  //Add to the integral
	  integrated_u += u_fe*w*J;
	}
      
      return integrated_u;
    }
    
  };
  
  
  //=======================================================================
  /// Face geometry for element is the same as that for the underlying
  /// wrapped element
  //=======================================================================
  template<>
  class FaceGeometry<MyPoissonElement>
    : public virtual
  FaceGeometry<TPoissonElement<2,3> >
  {
  public:
    FaceGeometry() :
      FaceGeometry<TPoissonElement<2,3> >() {}
  };
  
  
  //=======================================================================
  /// Face geometry of the Face Geometry for element is the same as
  /// that for the underlying wrapped element
  //=======================================================================
  template<>
  class FaceGeometry<FaceGeometry<MyPoissonElement> >
    : public virtual
  FaceGeometry<FaceGeometry<TPoissonElement<2,3> > > 
  {
  public:
    FaceGeometry() :
      FaceGeometry<FaceGeometry<TPoissonElement<2,3> > > () {}
  };


} // namespace oomph



//==start_of_problem_class============================================
/// Class definition
//====================================================================
template<class ELEMENT>
class UnstructuredPoissonProblem : public virtual Problem
{

public:

 /// Constructor
  UnstructuredPoissonProblem(const double &a, const double &b, const double &m);
    
 /// Destructor
 ~UnstructuredPoissonProblem()
  {
    Trace_file.close();
  };

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
  void doc_solution(const std::string& comment="",
		    const bool &write_trace=true);
 

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

  /// Storage for geometric parameters
  double A, B, M;
  
}; // end_of_problem_class





//==start_constructor=====================================================
/// Constructor
//========================================================================
template<class ELEMENT>
UnstructuredPoissonProblem<ELEMENT>::UnstructuredPoissonProblem
(const double &geometry_a, const double &geometry_b,
 const double &geometry_m) : A(geometry_a), B(geometry_b), M(geometry_m)
{          
 // Intrinsic coordinate along GeomObject
 Vector<double> zeta(1);

 // Position vector on GeomObject
 Vector<double> posn(2);
 
 // HomotopyTriangleEllipse defining the outer boundary
 HomotopyTriangleEllipseCartesian * outer_boundary_ellipse_pt =
   new HomotopyTriangleEllipseCartesian(geometry_a,
					geometry_b,
				        geometry_m);

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
   }
       
 // Pointer to the closed curve that defines the outer boundary
 TriangleMeshClosedCurve* closed_curve_pt=0;

 // Build outer boundary as curvilinear
 //------------------------------------
 
 // Provide storage for pointers to the two parts of the curvilinear boundary
 Vector<TriangleMeshCurveSection*> outer_curvilinear_boundary_pt(2);
 
 // First bit (upper half of the curve)
 //------------------------------------
 double zeta_start=0.0;
 double zeta_end=MathematicalConstants::Pi;
 unsigned nsegment=20;
 unsigned boundary_id=0;
 outer_curvilinear_boundary_pt[0]=
   new TriangleMeshCurviLine(
			     outer_boundary_ellipse_pt,
			     zeta_start,zeta_end,nsegment,boundary_id);
 
 // Second bit (lower half of the curve)
 //--------------------------------------
 zeta_start=MathematicalConstants::Pi;
 zeta_end=2.0*MathematicalConstants::Pi;
 nsegment=20;
 boundary_id=1;
 outer_curvilinear_boundary_pt[1]=
   new TriangleMeshCurviLine(
			     outer_boundary_ellipse_pt,
			     zeta_start,zeta_end,nsegment,boundary_id);
   
 // Combine to curvilinear boundary and define the
 //--------------------------------
 // outer boundary
 //--------------------------------
 closed_curve_pt=
   new TriangleMeshClosedCurve(outer_curvilinear_boundary_pt);
 
  
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
 My_mesh_pt->max_element_size()=0.1;
 My_mesh_pt->min_element_size()=0.001; 
 
 // Set boundary condition and complete the build of all elements
 complete_problem_setup();

 // Open trace file
 char filename[100];
 sprintf(filename,"RESLT/trace.dat");
 Trace_file.open(filename,std::ios_base::app);

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
   el_pt->source_fct_pt() = &ConstantSourceFunction::get_source;
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
void UnstructuredPoissonProblem<ELEMENT>::
doc_solution(const std::string& comment, const bool &write_trace)
{ 
 ofstream some_file;
 char filename[100];
 
 // Number of plot points
 unsigned npts;
 npts=5; 
 
 sprintf(filename,"RESLT/soln_A%g_B%g_M%g_%i.dat",A,B,M,Doc_info.number());
 some_file.open(filename);
 this->My_mesh_pt->output(some_file,npts); 
 some_file << "TEXT X = 22, Y = 92, CS=FRAME T = \"" 
           << comment << "\"\n";
 some_file.close();
  
 // Output boundaries
 //------------------
 sprintf(filename,"RESLT/boundaries_A%g_B%g_M%g_%i.dat",
	 A,B,M,Doc_info.number());
 some_file.open(filename);
 My_mesh_pt->output_boundaries(some_file);
 some_file.close();

 if(write_trace)
   {
     //Calculate the area and flux
     double area=0.0;
     double flux = 0.0;
     unsigned n_element = this->My_mesh_pt->nelement();
     for(unsigned e=0;e<n_element;++e)
       {
	 MyPoissonElement* element_pt=dynamic_cast<MyPoissonElement*>
	   (this->My_mesh_pt->finite_element_pt(e));
	 area += element_pt->size();
	 flux += element_pt->integrate_u();
       }
     
     //Calculate the perimeter
     double perimeter=0.0;
     //Loop over the two boundaries
     for(unsigned b=0;b<2;++b)
       {
	 unsigned n_bound_element = this->My_mesh_pt->nboundary_element(b);
	 for(unsigned e=0;e<n_bound_element;++e)
	   {
	     //Create a temporary flux element adjacent to the boundary
	     FiniteElement* boundary_el_pt =
	       new PoissonFluxElement<MyPoissonElement>
	       (this->My_mesh_pt->boundary_element_pt(b,e),
		this->My_mesh_pt->face_index_at_boundary(b,e));
	     
	     //Add the size of the element to the perimeter
	     perimeter += boundary_el_pt->size();
	     //Delete the element
	     delete boundary_el_pt;
	   }
       }
     
     
     double fRe = 8.0*area*area*area/(perimeter*perimeter*flux);
     
     Trace_file << A  << " " << B << " "
		<< M<< " " << area << " " << perimeter
		<< " " << flux << " " << fRe << "\n";
   }
 
 // Increment the doc_info number
 Doc_info.number()++;

} // end of doc


//=======start_of_main========================================
/// Driver code for demo of inline triangle mesh generation
//============================================================
int main(int argc, char **argv)
{
  //Sanity check (circle)
  {
    //Set the geometric parameters
    double A = 1.0;
    double B = 1.0;
    double M = 0.0;
    
    // Create problem
    UnstructuredPoissonProblem<MyPoissonElement> problem(A,B,M);
    
    // Doc the initial mesh
    //=====================
    {
      std::stringstream comment_stream;
      comment_stream << "Initial mesh ";
      problem.doc_solution(comment_stream.str(),false);
    } 
    
    // Solve with spatial adaptation
    //==============================
    unsigned max_adapt=3;
    problem.newton_solve(max_adapt);
    
    // Doc the solution
    //=================
    problem.doc_solution();
  }


  //First geometry
  {
    //Set the geometric parameters
    double A = 1.0;
    double B = 0.5766;//0.6;
    double M = 0.3529;//0.1;
    
    // Create problem
    UnstructuredPoissonProblem<MyPoissonElement> problem(A,B,M);
    
    // Doc the initial mesh
    //=====================
    {
      std::stringstream comment_stream;
      comment_stream << "Initial mesh ";
      problem.doc_solution(comment_stream.str(),false);
    } 
    
    // Solve with spatial adaptation
    //==============================
    unsigned max_adapt=3;
    problem.newton_solve(max_adapt);
    
    // Doc the solution
    //=================
    problem.doc_solution();
  }


  //Second Geometry
  {
    double A = 1.0;
    double B = 1.4223;//0.6;
    double M = 0.3336;//0.1;
    
    // Create problem
    UnstructuredPoissonProblem<MyPoissonElement> problem(A,B,M);
    
    // Doc the initial mesh
    //=====================
    {
      std::stringstream comment_stream;
      comment_stream << "Initial mesh ";
      problem.doc_solution(comment_stream.str(),false);
    } 
    
    // Solve with spatial adaptation
    //==============================
    unsigned max_adapt=3;
    problem.newton_solve(max_adapt);
    
    // Doc the solution
    //=================
    problem.doc_solution();
  }
 

  
} //End of main
