
// This header file
// contains a bunch of derived classes which are easier to find when separated,
// we have derived elements whcih output some specified quantities, a transvsersely
// isotropic strain energy function constitutive law, several strain energy functions
// and the definition of a displaced boundary
void contour_aligned_fibres_superposition(const Vector<double> &xi, const double &s
 										,const double &t1, const double &t2, Vector<double> &M);

double phi_spatial(const double &z, const double &z0, const double &t1, const double &t2, const double &s, const double &phi);

//====================================================================
//Continuous Pressure QPVDElements
//====================================================================

namespace oomph
{
 template<unsigned DIM>
 class MyQPVDElementWithContinuousPressure : public virtual QPVDElementWithContinuousPressure<DIM>
{
public:

 /// Constructor: Call constructor of underlying element
  MyQPVDElementWithContinuousPressure() :  QPVDElementWithContinuousPressure<DIM>() {}

  /// Overload output function
  void output(std::ostream &outfile, const unsigned &n_plot)
   {
    Vector<double> s(3);
    Vector<double> x(3);
    Vector<double> xi(3);
    DenseMatrix<double> sigma(3,3);
    DenseMatrix<double> deformation_gradient(3,3);
    DenseMatrix<double> cauchy_stress(3,3);
    DenseMatrix<double> right_cauchy_green(3,3);
    Vector<double> M(3);
    Vector<double> FM(3);

   // initialise to zero
   for(unsigned i = 0; i < 3; i++)
   {
 	  M[i] = 0;
 	  FM[i] = 0;
   }

    //Tecplot header info
    outfile << "ZONE I=" << n_plot << ", J=" << n_plot << ", K=" << n_plot << std::endl;

      //Loop over element nodes
      for(unsigned l3=0;l3<n_plot;l3++)
       {
        s[2] = -1.0 + l3*2.0/(n_plot-1);
        for(unsigned l2=0;l2<n_plot;l2++)
         {
          s[1] = -1.0 + l2*2.0/(n_plot-1);
          for(unsigned l1=0;l1<n_plot;l1++)
           {
            s[0] = -1.0 + l1*2.0/(n_plot-1);

        // Get Eulerian and Lagrangian coordinates and the stress
        this->interpolated_x(s,x);
        this->interpolated_xi(s,xi);
        this->get_stress(s,sigma);

        //Output the x,y,z coordinates
        for(unsigned i=0;i<3;i++)
   	 { outfile << x[i] << " ";}

 	 // number of local nodes for the element
 	 unsigned n_node = this->nnode();
 	 // number of coordinate types required to interpolate between the nodes
 	 unsigned n_position_type = this->nnodal_position_type();
 	 // create the shape functions
 	 Shape psi(n_node,n_position_type);
 	 // and their derivatives
 	 DShape dpsi(n_node,n_position_type,3);
 	 // calculate the shape functions
 	 this->shape(s,psi);
 	 // and their derivatives
 	 this->dshape_lagrangian(s,psi,dpsi);

 	 // initialise to zero
 	 for(unsigned i=0;i<3;i++)
 	 {
 		 for(unsigned j=0;j<3;j++)
 		 {
 			 deformation_gradient(i,j) = 0.0;
 			 cauchy_stress(i,j) = 0.0;
 			 right_cauchy_green(i,j) = 0.0;
 		 }
 	 }

 	 // loop over the number of nodes
 	 for(unsigned j=0;j<n_node;j++)
 	 {
 		 // loop over the dimensions
 		 for(unsigned i=0;i<3;i++)
 		 {
 			 // loop over the derivative directions
 			 for(unsigned k=0;k<3;k++)
 			 {
 				 deformation_gradient(i,k) += this->nodal_position(j,i)*dpsi(j,k);
 			 }
 		 }
 	 }

 	 // calculate the cauchy stress from the second piola kirchhoff stress (sigma)
 	 for(unsigned i=0;i<3;i++)
 	 {
 		 for(unsigned j=0;j<3;j++)
 		 {
 			 for(unsigned k=0;k<3;k++)
 			 {
 				 for(unsigned l=0;l<3;l++)
 				 {
 					 cauchy_stress(i,j) +=
 					 deformation_gradient(i,k)*sigma(k,l)*deformation_gradient(j,l);
 				 }
 			 }
 		 }
 	 }


 	 for(unsigned i=0;i<3;i++)
 	 {
 		 for(unsigned j=0;j<3;j++)
 		 {
 			 for(unsigned k=0;k<3;k++)
 			 {
 				 right_cauchy_green(i,j) += deformation_gradient(k,i)*deformation_gradient(k,j);
 			 }
 		 }
 	 }

         //Output xi[0], xi[1], xi[2]
         for(unsigned i=0;i<3;i++)
         {outfile << xi[i] << " ";}


		 // find the lagrangian fibre direction
		 if (Global_Parameters::contour_fibres)
		 {
			 contour_aligned_fibres_superposition(xi, Global_Parameters::s , Global_Parameters::t1, Global_Parameters::t2 , M);
		 }
		 else
		 {
			 M[0] = 0.0;
			 M[1] = 0.0;
			 M[2] = 1.0;
	 	}

 		// compute the deformed fibre direction
 		for(unsigned i = 0; i < 3; i++)
 		{
 			for(unsigned j = 0; j < 3; j++)
 			{
 				FM[i] = deformation_gradient(i,j)*M[j];
 			}
 		}
 		// turn into unit vector
 		double FM_magnitude = sqrt(FM[0]*FM[0] + FM[1]*FM[1] + FM[2]*FM[2]);
 		for(unsigned i = 0; i < 3; i++)
 		{
 			FM[i] = FM[i]/FM_magnitude;
 		}


 		// calculate I4
 		double I4 = 0.0;
 		for(unsigned i=0;i<3;i++)
 		  {
 			for(unsigned j=0;j<3;j++)
 			  {
 				I4 += M[i]*right_cauchy_green(i,j)*M[j];
 			  }
 		  }

        //Output stress
        outfile    << sigma(0,0) << " "
        	       << sigma(1,1) << " "
        	       << sigma(2,2) << " "
        	       << sigma(0,1) << " "
        	       << sigma(0,2) << " "
        	       << sigma(1,2) << " "
			   	   // << cauchy_stress(0,0) << " "
				   // << cauchy_stress(1,1) << " "
				   // << cauchy_stress(2,2) << " "
				   // << cauchy_stress(0,1) << " "
				   // << cauchy_stress(0,2) << " "
				   // << cauchy_stress(1,2) << " "
				   << deformation_gradient(0,0) << " "
				   << deformation_gradient(1,1) << " "
				   << deformation_gradient(2,2) << " "
				   << deformation_gradient(0,1) << " "
				   << deformation_gradient(1,0) << " "
				   << deformation_gradient(0,2) << " "
				   << deformation_gradient(2,0) << " "
				   << deformation_gradient(1,2) << " "
				   << deformation_gradient(2,1) << " "
 			       << FM[0] << " "
 			       << FM[1] << " "
 			       << FM[2] << " "
 			       << I4 << " "
                   << std::endl;


	}
}
}
} // end of output



};

//=========================================================================
///FaceGeometry of the 3D QPVDElementWithContinuousPressure
//=========================================================================
template<>
class FaceGeometry<MyQPVDElementWithContinuousPressure<3> >:
 public virtual SolidQElement<2,3>
{
  public:
 //Make sure that we call the constructor of the SolidQElement
 //Only the Intel compiler seems to need this!
 FaceGeometry() : SolidQElement<2,3>() {}
};

//===============================================================
/// FaceGeometry of FaceGeometry
/// for 3D MyQPVDElementWithContinuousPressure element
//===============================================================
template<>
class FaceGeometry<FaceGeometry<MyQPVDElementWithContinuousPressure<3> > >:
public virtual SolidQElement<1,3>
{
  public:
 /// Constructor must call constructor of the underlying element
  FaceGeometry() : SolidQElement<1,3>() {}
};

//////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

  template<unsigned DIM>
  class MyRefineableQPVDElementWithContinuousPressure : public virtual RefineableQPVDElementWithContinuousPressure<DIM>
{

public:

 /// Constructor: Call constructor of underlying element
  MyRefineableQPVDElementWithContinuousPressure() :  RefineableQPVDElementWithContinuousPressure<DIM>() {}

  /// Overload output function
  void output(std::ostream &outfile, const unsigned &n_plot)
   {
	Vector<double> s(3);
	Vector<double> x(3);
	Vector<double> xi(3);
	DenseMatrix<double> sigma(3,3);
	DenseMatrix<double> deformation_gradient(3,3);
	DenseMatrix<double> cauchy_stress(3,3);
	DenseMatrix<double> right_cauchy_green(3,3);
	Vector<double> M(3);
	Vector<double> FM(3);



	//Tecplot header info
	outfile << "ZONE I=" << n_plot << ", J=" << n_plot << ", K=" << n_plot << std::endl;

	  //Loop over element nodes
	  for(unsigned l3=0;l3<n_plot;l3++)
	   {
		s[2] = -1.0 + l3*2.0/(n_plot-1);
		for(unsigned l2=0;l2<n_plot;l2++)
		 {
		  s[1] = -1.0 + l2*2.0/(n_plot-1);
		  for(unsigned l1=0;l1<n_plot;l1++)
		   {
			s[0] = -1.0 + l1*2.0/(n_plot-1);

		// Get Eulerian and Lagrangian coordinates and the stress
		this->interpolated_x(s,x);
		this->interpolated_xi(s,xi);
		this->get_stress(s,sigma);

		//Output the x,y,z coordinates
	    for(unsigned i=0;i<3;i++)
	   { outfile << x[i] << " ";}

   // number of local nodes for the element
   unsigned n_node = this->nnode();
   // number of coordinate types required to interpolate between the nodes
   unsigned n_position_type = this->nnodal_position_type();
   // create the shape functions
   Shape psi(n_node,n_position_type);
   // and their derivatives
   DShape dpsi(n_node,n_position_type,3);
   // calculate the shape functions
   this->shape(s,psi);


   // and their derivatives
   this->dshape_lagrangian(s,psi,dpsi);

   // initialise to zero
   for(unsigned i=0;i<3;i++)
   {
	   for(unsigned j=0;j<3;j++)
	   {
		   deformation_gradient(i,j) = 0.0;
		   cauchy_stress(i,j) = 0.0;
		   right_cauchy_green(i,j) = 0.0;
	   }
   }


   // initialise to zero
   for(unsigned i = 0; i < 3; i++)
   {
	M[i] = 0;
	FM[i] = 0;
   }

   // loop over the number of nodes
   for(unsigned j=0;j<n_node;j++)
   {
	   // loop over the dimensions
	   for(unsigned i=0;i<3;i++)
	   {
		   // loop over the derivative directions
		   for(unsigned k=0;k<3;k++)
		   {
			   deformation_gradient(i,k) += this->nodal_position(j,i)*dpsi(j,k);
		   }
	   }
   }

   // calculate the cauchy stress from the second piola kirchhoff stress (sigma)
   for(unsigned i=0;i<3;i++)
   {
	   for(unsigned j=0;j<3;j++)
	   {
		   for(unsigned k=0;k<3;k++)
		   {
			   for(unsigned l=0;l<3;l++)
			   {
				   cauchy_stress(i,j) +=
				   deformation_gradient(i,k)*sigma(k,l)*deformation_gradient(j,l);
			   }
		   }
	   }
   }


   for(unsigned i=0;i<3;i++)
   {
	   for(unsigned j=0;j<3;j++)
	   {
		   for(unsigned k=0;k<3;k++)
		   {
			   right_cauchy_green(i,j) += deformation_gradient(k,i)*deformation_gradient(k,j);
		   }
	   }
   }

		 //Output xi[0], xi[1], xi[2]
		 for(unsigned i=0;i<3;i++)
		 {outfile << xi[i] << " ";}


	   // find the lagrangian fibre direction
	   if (Global_Parameters::contour_fibres)
	   {
		   contour_aligned_fibres_superposition(xi,Global_Parameters::s , Global_Parameters::t1, Global_Parameters::t2 , M);
	   }
	   else
	   {
		   M[0] = 0.0;
		   M[1] = 0.0;
		   M[2] = 1.0;
	  }

	  // compute the deformed fibre direction
	  for(unsigned i = 0; i < 3; i++)
	  {
		  for(unsigned j = 0; j < 3; j++)
		  {
			  FM[i] = deformation_gradient(i,j)*M[j];
		  }
	  }
	  // turn into unit vector
	  double FM_magnitude = sqrt(FM[0]*FM[0] + FM[1]*FM[1] + FM[2]*FM[2]);
	  for(unsigned i = 0; i < 3; i++)
	  {
		  FM[i] = FM[i]/FM_magnitude;
	  }


	  // calculate I4
	  double I4 = 0.0;
	  for(unsigned i=0;i<3;i++)
		{
		  for(unsigned j=0;j<3;j++)
			{
			  I4 += M[i]*right_cauchy_green(i,j)*M[j];
			}
		}

		//Output stress
        outfile    << sigma(0,0) << " "
        	       << sigma(1,1) << " "
        	       << sigma(2,2) << " "
        	       << sigma(0,1) << " "
        	       << sigma(0,2) << " "
        	       << sigma(1,2) << " "
				   // << cauchy_stress(0,0) << " "
           	       // << cauchy_stress(1,1) << " "
           	       // << cauchy_stress(2,2) << " "
           	       // << cauchy_stress(0,1) << " "
           	       // << cauchy_stress(0,2) << " "
           	       // << cauchy_stress(1,2) << " "
				   << deformation_gradient(0,0) << " "
				   << deformation_gradient(1,1) << " "
				   << deformation_gradient(2,2) << " "
				   << deformation_gradient(0,1) << " "
				   << deformation_gradient(1,0) << " "
				   << deformation_gradient(0,2) << " "
				   << deformation_gradient(2,0) << " "
				   << deformation_gradient(1,2) << " "
				   << deformation_gradient(2,1) << " "
 			       << FM[0] << " "
 			       << FM[1] << " "
 			       << FM[2] << " "
 			       << I4 << " "
                   << std::endl;

	}
}
}
} // end of output



  }; //end MyRefineableQPVDElementWithContinuousPressure

//=========================================================================
///FaceGeometry of the 3D MyRefineableQPVDElementWithContinuousPressure
//=========================================================================
template<>
class FaceGeometry<MyRefineableQPVDElementWithContinuousPressure<3> >:
 public virtual SolidQElement<2,3>
{
  public:
 //Make sure that we call the constructor of the SolidQElement
 //Only the Intel compiler seems to need this!
 FaceGeometry() : SolidQElement<2,3>() {}
};

//===============================================================
/// FaceGeometry of FaceGeometry
/// for 3D MyRefineableQPVDElementWithContinuousPressure element
//===============================================================
template<>
class FaceGeometry<FaceGeometry<MyRefineableQPVDElementWithContinuousPressure<3> > >:
public virtual SolidQElement<1,3>
{
  public:
 /// Constructor must call constructor of the underlying element
  FaceGeometry() : SolidQElement<1,3>() {}
};


//===============================================================
// A class derived from ImposeDisplacementByLagrangeMultiplierElements
// that include a new output function designed to output the
// resultant force over the top of the cylinder
//===============================================================
template<class ELEMENT>
class MyLagrangeMultiplierElement : public ImposeDisplacementByLagrangeMultiplierElement<ELEMENT>
{
public:


  MyLagrangeMultiplierElement(FiniteElement* const &element_pt, const int &face_index):
    ImposeDisplacementByLagrangeMultiplierElement<ELEMENT>(element_pt,face_index) { }






		/// \short Output function
	  void output(std::ostream &outfile, const unsigned &n_plot)
	   {
	    // Elemental dimension
	    unsigned dim_el=this->dim();

	    //Find the number of positional types
	    unsigned n_position_type = this->nnodal_position_type();

	 #ifdef PARANOID
	    if(n_position_type!=1)
	     {
	      throw OomphLibError(
	       "ImposeDisplacementByLagrangeMultiplierElement cannot (currently) be used with elements that have generalised positional dofs",
	       OOMPH_CURRENT_FUNCTION,
	       OOMPH_EXCEPTION_LOCATION);
	     }
	 #endif


	    //Local coord
	    Vector<double> s(dim_el);

	    // # of nodes,
	    unsigned n_node=this->nnode();
	    Shape psi(n_node,n_position_type);

	    // Tecplot header info
	    //outfile << this->tecplot_zone_string(n_plot);

	    // Loop over plot points
	    unsigned num_plot_points=this->nplot_points(n_plot);
	    for (unsigned iplot=0;iplot<num_plot_points;iplot++)
	     {
	      // Get local coordinates of plot point
	      this->get_s_plot(iplot,n_plot,s);

	      // Get shape function
	      this->shape(s,psi);

	      //Calculate the Eulerian coordinates and Lagrange multiplier
	      Vector<double> x(dim_el+1,0.0);
	      Vector<double> lambda(dim_el+1,0.0);
	      Vector<double> zeta(dim_el,0.0);
//				Vector<double> resultant_force(dim_el+1,0.0);
	      for(unsigned j=0;j<n_node;j++)
	       {
	        // Cast to a boundary node
	        BoundaryNodeBase *bnod_pt =
	         dynamic_cast<BoundaryNodeBase*>(this->node_pt(j));

	        // get the node pt
	        Node* nod_pt = this->node_pt(j);

	        // Get the index of the first nodal value associated with
	        // this FaceElement
	        unsigned first_index=
	         bnod_pt->index_of_first_value_assigned_by_face_element(this->Id);

	        // higher dimensional quantities
	        for(unsigned i=0;i<dim_el+1;i++)
	         {
	          x[i]+=this->nodal_position(j,i)*psi(j,0); // need to sort
	                                              // this out properly
	                                              // for generalised dofs
	          lambda[i]+=nod_pt->value
	           (first_index+i)*psi(j,0);
	         }
	        //In-element quantities
	        for(unsigned i=0;i<dim_el;i++)
	         {
	          //Loop over positional types
	          for (unsigned k=0;k<n_position_type;k++)
	           {
	            zeta[i]+=this->zeta_nodal(j,k,i)*psi(j,k);
	           }
	         }
	       }

	      // Get prescribed wall shape
	      Vector<double> r_prescribed(dim_el+1);
	      this->Boundary_shape_geom_object_pt->position(zeta,r_prescribed);


			}

				Vector<double> resultant_force(dim_el+1,0.0);
				this->compute_resultant_force(resultant_force);


			 	outfile << this->size() << " ";
			  outfile << resultant_force[0] << " ";
				outfile << resultant_force[1] << " ";
				outfile << resultant_force[2] << std::endl;


	   }


	 /// \short Compute the resultant force
 void compute_resultant_force(Vector<double> &resultant_force)
  {
   //Find out how many positional dofs there are
   unsigned n_position_type = this->nnodal_position_type();

#ifdef PARANOID
   if(n_position_type!=1)
    {
     throw OomphLibError(
      "ImposeDisplacementByLagrangeMultiplierElement cannot (currently) be used with elements that have generalised positional dofs",
      OOMPH_CURRENT_FUNCTION,
      OOMPH_EXCEPTION_LOCATION);
    }
#endif

   //Find out how many nodes there are
   unsigned n_node = this->nnode();

   // Dimension of element
   unsigned dim_el=this->dim();

   //Set up memory for the shape functions
   Shape psi(n_node);
   DShape dpsids(n_node,dim_el);

   //Set the value of n_intpt
   unsigned n_intpt = this->integral_pt()->nweight();


   // Initialise resultant force
   //Vector<double> resultant_force(dim_el+1,0.0);

   //Loop over the integration points
   for(unsigned ipt=0;ipt<n_intpt;ipt++)
    {
     //Get the integral weight
     double w = this->integral_pt()->weight(ipt);

     //Only need to call the local derivatives
     this->dshape_local_at_knot(ipt,psi,dpsids);

     //Calculate the Eulerian coordinates and Lagrange multiplier
     Vector<double> x(dim_el+1,0.0);
     Vector<double> lambda(dim_el+1,0.0);
     Vector<double> zeta(dim_el,0.0);
     DenseMatrix<double> interpolated_a(dim_el,dim_el+1,0.0);

     // Loop over nodes
     for(unsigned j=0;j<n_node;j++)
      {
       Node* nod_pt=this->node_pt(j);

       // Cast to a boundary node
       BoundaryNodeBase *bnod_pt =
        dynamic_cast<BoundaryNodeBase*>(this->node_pt(j));

       // Get the index of the first nodal value associated with
       // this FaceElement
       unsigned first_index=
        bnod_pt->index_of_first_value_assigned_by_face_element(this->Id);

       //Assemble higher-dimensional quantities
       for(unsigned i=0;i<dim_el+1;i++)
        {
         x[i]+=this->nodal_position(j,i)*psi(j);
         lambda[i]+=nod_pt->value(first_index+i)*psi(j);
         for(unsigned ii=0;ii<dim_el;ii++)
          {
           interpolated_a(ii,i) +=
            this->lagrangian_position(j,i)*dpsids(j,ii);
          }
        }
       if (!this->Sparsify)
        {
         for(unsigned k=0;k<n_position_type;k++)
          {
           //Assemble in-element quantities: boundary coordinate
           for(unsigned i=0;i<dim_el;i++)
            {
             zeta[i]+=this->zeta_nodal(j,k,i)*psi(j,k);
            }
          }
        }
      }

     if (this->Sparsify) zeta=this->Zeta_sub_geom_object[ipt];


     //Now find the local undeformed metric tensor from the tangent Vectors
     DenseMatrix<double> a(dim_el);
     for(unsigned i=0;i<dim_el;i++)
      {
       for(unsigned j=0;j<dim_el;j++)
        {
         //Initialise surface metric tensor to zero
         a(i,j) = 0.0;
         //Take the dot product
         for(unsigned k=0;k<dim_el+1;k++)
          {
           a(i,j) += interpolated_a(i,k)*interpolated_a(j,k);
          }
        }
      }


     //Find the determinant of the metric tensor
     double adet =0.0;
     switch(dim_el+1)
      {

      case 2:
       adet = a(0,0);
       break;

      case 3:
       adet = a(0,0)*a(1,1) - a(0,1)*a(1,0);
       break;

      default:
       throw
        OomphLibError(
         "Wrong dimension fill_in_generic_contribution_to_residuals_displ_lagr_multiplier",
         "ImposeDisplacementByLagrangeMultiplierElement::fill_in_generic_contribution_to_residuals_displ_lagr_multiplier()",
         OOMPH_EXCEPTION_LOCATION);
      }

     // Get prescribed wall shape
     Vector<double> r_prescribed(dim_el+1);
     if (!this->Sparsify)
      {
       this->Boundary_shape_geom_object_pt->position(zeta,r_prescribed);
      }
     else
      {
       this->Sub_geom_object_pt[ipt]->position(zeta,r_prescribed);
      }

     //Premultiply the weights and the square-root of the determinant of
     //the metric tensor
     double W = w*sqrt(adet);

     // Assemble error

     //Loop over directions
     for(unsigned i=0;i<dim_el+1;i++)
      {
       resultant_force[i] += -lambda[i]*W;
      }


    } //End of loop over the integration points




   //return resultant_force;

  }


};





  //===============================================================
  // A class for transversely isotropic constitutive laws -
  // inherits from isotropic constitutive class
  //===============================================================
  class TransverselyIsotropicStrainEnergyFunctionConstitutiveLaw : public
  IsotropicStrainEnergyFunctionConstitutiveLaw
  {
  public:

    /// constructor takes a pointer to the strain energy function
    TransverselyIsotropicStrainEnergyFunctionConstitutiveLaw(
							     StrainEnergyFunction* const &strain_energy_function_pt) :
      IsotropicStrainEnergyFunctionConstitutiveLaw(strain_energy_function_pt) {}

    //===========================================================================
    /// Calculate the deviatoric part
    /// \f$ \overline{ \sigma^{ij}}\f$  of the contravariant
    /// 2nd Piola Kirchhoff stress tensor \f$ \sigma^{ij}\f$.
    /// Also return the contravariant deformed metric
    /// tensor and the determinant of the deformed metric tensor.
    /// Uses correct 3D invariants for 2D (plane strain) problems.
    /// This is the version for the pure incompressible formulation.
    //============================================================================
    void calculate_second_piola_kirchhoff_stress(
						 const DenseMatrix<double> &g,
						 const DenseMatrix<double> &G,
						 DenseMatrix<double> &sigma_dev,
						 DenseMatrix<double> &Gup,
						 double &detG,
						 const Vector<double> &xi);

  };

//===============================================================
// member function which calculates the 2nd piola-kirchoff stress
//===============================================================
void TransverselyIsotropicStrainEnergyFunctionConstitutiveLaw::
calculate_second_piola_kirchhoff_stress(const DenseMatrix<double> &g,
					const DenseMatrix<double> &G,DenseMatrix<double> &sigma_dev,
					DenseMatrix<double> &Gup, double &detG, const Vector<double> &xi)
{
  // Error checking
#ifdef PARANOID
  error_checking_in_input(g,G,sigma_dev);
#endif

  // find the dimension of the problem
  unsigned dim = g.nrow();

#ifdef PARANOID
  if (dim==1 || dim==2)
    {
      std::string function_name =
	"TransverselyIsotropicStrainEnergyFunctionConstitutiveLaw::";
      function_name += "calculate_second_piola_kirchhoff_stress()";

   throw OomphLibError(
		       "Check constitutive equations carefully when dim=1",
		       OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

  // calculate the contravariant undeformed and deformed metric tensors
  // and get the determinants of the metric tensors
  DenseMatrix<double> gup(dim);
  //double  detg = calculate_contravariant(g,gup);
  (void)calculate_contravariant(g,gup);
  detG = calculate_contravariant(G,Gup);

  // Define the direction of the fibres in this transversely isotropic model
  // the fibres point in the direction of a vector M - alignment
  Vector<double> M(3,0.0);
  if (Global_Parameters::contour_fibres)
  {
	  contour_aligned_fibres_superposition(xi,Global_Parameters::s , Global_Parameters::t1, Global_Parameters::t2 , M);
	  //std::cout << M[0] << std::endl;
  }
  else
  {
	  M[0] = 0.0;
	  M[1] = 0.0;
	  M[2] = 1.0;
 }

  // calculate the strain invariants
  Vector<double> I(5,0.0);
  /// The third strain invariant is the volume change = detG/detg, set as 1
  I[2] = 1.0;
  // the first strain invariant is equal to g^{ij}G_{ij}
  // the second is equal to G^{rs}g_{rs}I3
  // the third is related to the direction of the fibres
  //  - all can be put within one set of loops
  for(unsigned i=0;i<dim;i++)
    {
      for(unsigned j=0;j<dim;j++)
	{
	  I[0] += gup(i,j)*G(i,j);
	  I[1] += Gup(i,j)*g(i,j);
	  I[3] += M[i]*M[j]*G(i,j);
	}
    }

  // the fifth strain invariant needs four loops
  // I5 = M_i * M_j * g^{rs}G_{is}G_{sj}
  for(unsigned i=0;i<dim;i++)
    {
      for(unsigned j=0;j<dim;j++)
	{
	  for(unsigned r=0;r<dim;r++)
	     {
	       for(unsigned s=0;s<dim;s++)
		 {
		   I[4] += M[i]*G(i,r)*gup(r,s)*G(s,j)*M[j];
		 }
	     }
	}
    }

  //Second strain invariant is multiplied by the third.
  //I[1] *= I[2];

  //Calculate the derivatives of the strain energy function wrt the
  //strain invariants
  Vector<double> dWdI(5,0.0);
   Strain_energy_function_pt->derivatives(I,dWdI);

	 //Only bother to compute the tensor B^{ij} (Green & Zerna notation)
   //if the derivative wrt the second strain invariant is non-zero
   DenseMatrix<double> Bup(dim,dim,0.0);
   if(std::fabs(dWdI[1]) > 0.0)
     {
       for(unsigned i=0;i<dim;i++)
	 {
	   for(unsigned j=0;j<dim;j++)
	     {
	       Bup(i,j) = I[0]*gup(i,j);
	       for(unsigned r=0;r<dim;r++)
		 {
		   for(unsigned s=0;s<dim;s++)
		     {
		       Bup(i,j) -= gup(i,r)*gup(j,s)*G(r,s);
		     }
		 }
	     }
	 }
     }

   //Only bother to compute the tensor B^{ij} (Green & Zerna notation)
   //if the derivative wrt the second strain invariant is non-zero
   DenseMatrix<double> Mup(dim,dim,0.0);
   if(std::fabs(dWdI[3]) > 0.0)
     {
       for(unsigned i=0;i<dim;i++)
	 {
	   for(unsigned j=0;j<dim;j++)
	     {
	       Mup(i,j) = M[i]*M[j];
	     }
	 }
     }

   //Only bother to compute the tensor C^{ij} (Green & Zerna notation)
   //if the derivative wrt the fifth strain invariant is non-zero
   DenseMatrix<double> Cup(dim,dim,0.0);
   if(std::fabs(dWdI[4]) > 0.0)
     {
       for(unsigned i=0;i<dim;i++)
	 {
	   for(unsigned j=0;j<dim;j++)
	     {
	       for(unsigned r=0;r<dim;r++)
		 {
		   for(unsigned s=0;s<dim;s++)
		     {
		       Cup(i,j) += M[i]*G(r,s)*gup(r,j)*M[s] + G(r,s)*gup(r,i)*M[s]*M[j];
		     }
		 }
	     }
	 }
     }

   // Now set the values of the functions phi and psi (Green and Zerna)
   // and chi and omega (derived due to anisotropy)
   double phi = 2.0*dWdI[0]*(Global_Parameters::phi_sv ? (1.0 -
   							phi_spatial(xi[2], 0.5*Global_Parameters::length, Global_Parameters::t1, Global_Parameters::t2,
							Global_Parameters::s, Global_Parameters::Phi)) : (1.0 - Global_Parameters::Phi));
   double psi = 2.0*dWdI[1];
   double chi = 2.0*dWdI[3]*(Global_Parameters::phi_sv ? (phi_spatial(xi[2],
   							0.5*Global_Parameters::length, Global_Parameters::t1, Global_Parameters::t2,
   							Global_Parameters::s, Global_Parameters::Phi)) : Global_Parameters::Phi);
   double omega = 2.0*dWdI[4];
   // std::cout << phi_spatial(xi[2],
   // 							0.5*Global_Parameters::length, Global_Parameters::t1, Global_Parameters::t2,
   // 							Global_Parameters::s, Global_Parameters::Phi) << std::endl;
   //Calculate the trace/dim of the first two terms of the stress tensor
   // this gives the mechanical pressure
   double K;
   K = (I[0]*phi + 2.0*psi*I[1] + chi*I[3] + 2.0*omega*I[4])/3.0;

   // Now take the mechanical pressure P*G^{ij} away from sigma to
   // get the deviatoric part of the 2nd piola-kirchoff stress
   for(unsigned i=0;i<dim;i++)
     {
       for(unsigned j=0;j<dim;j++)
       {
	 sigma_dev(i,j) = phi*gup(i,j) + psi*Bup(i,j) + chi*Mup(i,j)
	   + omega*Cup(i,j) - K*Gup(i,j);
       }
     }


}


//===============================================================
// HGO (Holzapfel-Gasser-Ogden) strain energy function - this
// inherits from the standard strain energy function.
// W = (C1/2)*(I1-3) + k1/(2*k2)(exp(k2(I4-1)^2)-1)
//===============================================================
class HolzapfelGasserOgdenSEF : public StrainEnergyFunction
{

  public:

 /// Constructor takes the pointer to the value of the constants
 HolzapfelGasserOgdenSEF(double* c1_pt, double* k1_pt, double* k2_pt) : StrainEnergyFunction(),
   C1_pt(c1_pt), K1_pt(k1_pt), K2_pt(k2_pt)  {}


 /// Empty Virtual destructor
 virtual ~HolzapfelGasserOgdenSEF(){}

 /// Return the strain energy in terms of strain tensor
 double W(const DenseMatrix<double> &gamma)
  {return StrainEnergyFunction::W(gamma);}

 /// Return the strain energy in terms of the strain invariants
 double W(const Vector<double> &I)
  {std::cout << "computing strain energy" << std::endl;
		return (*C1_pt)*(I[0]-3.0)/2.0 + ((*K1_pt)/(2*(*K2_pt)))*(exp((*K2_pt)*(pow(I[3]-1.0,2.0)))-1.0);}


 /// \short Return the derivatives of the strain energy function with
 /// respect to the strain invariants
 void derivatives(Vector<double> &I, Vector<double> &dWdI)
  {
    dWdI[0] =(*C1_pt)/2.0;
   dWdI[1] = 0.0;
   dWdI[2] = 0.0;
   dWdI[3] = (*K1_pt)*(I[3]-1.0)*(exp((*K2_pt)*(pow(I[3]-1.0,2.0))));
   dWdI[4] = 0.0;
	 //std::cout << "lambda = " << sqrt(I[3]) << ", W = " << W(I) << std::endl;
  }

 /// \short Pure virtual function in which the user must declare if the
 /// constitutive equation requires an incompressible formulation
 /// in which the volume constraint is enforced explicitly.
 /// Used as a sanity check in PARANOID mode. True
  bool requires_incompressibility_constraint(){return true;}


  private:

 /// Pointer to constants
 double* C1_pt;
 double* K1_pt;
 double* K2_pt;

};



//====== Displaced Boundary ===============================
// Longitudinally displaced boundary in 3D
//===========================================================
class DisplacedBoundary : public GeomObject
{

public:

  /// constructor: specify amplitude of deflection from straight  horizontal line
  DisplacedBoundary(const double& ampl) : GeomObject(2,3)
  {
    Ampl=ampl;
  }

  /// Broken copy constructor
  DisplacedBoundary(const DisplacedBoundary& dummy)
  {
    BrokenCopy::broken_copy("DisplacedBoundary");
  }

  /// Broken assignment operator
  void operator=(const DisplacedBoundary&)
  {
    BrokenCopy::broken_assign("DisplacedBoundary");
  }

  /// Empty destructor
  ~DisplacedBoundary(){}


  // top surface of the cylinder shrinks
  void position(const Vector<double>& zeta, Vector<double>& r) const
  {
		if (Global_Parameters::pin_top_surface) {
			//top surface only changes in the z-direction
			r[0] = zeta[0];
			r[1] = zeta[1];
			r[2] =  Global_Parameters::length + Ampl;
		}
		else {
 	//change from (x,y) -> (x,y)/sqrt(lambda)
		r[0] = zeta[0]/pow(1.0+Ampl/Global_Parameters::length, 0.5);
	    r[1] = zeta[1]/pow(1.0+Ampl/Global_Parameters::length, 0.5);
	    r[2] = Global_Parameters::length + Ampl;
		if (Global_Parameters::print_height){
			std::cout << "boundary is at height: " << r[2] << std::endl;
			Global_Parameters::print_height = false;
		}
		}
		//std::cout << "z on surface = " << r[2] << std::endl;
		// if (r[0]*r[0] + r[1]*r[1] >= zeta[0]*zeta[0] + zeta[1]*zeta[1])
		// {
		// 	std::cout << "r increased" << std::endl;
		// }
  }

  /// Short parameterised position on object: r(zeta). Evaluated at
  /// previous timestep.
  void position(const unsigned& t, const Vector<double>& zeta,
		Vector<double>& r) const
  {
    position(zeta,r);
  }

  /// Access to amplitude
  double& ampl() {return Ampl;}

  /// how many items of data does the shape of the object depend on?
  unsigned ngeom_data() const
  {
    return 0;
  }

private:

  /// amplitude of perturbation
  double Ampl;

};
}

// triangular distribution
double triangular_distribution(const double &x, const double &a, const double &b, const double &c) {

	if ((a <= x) && (x < c)) {
		return 2.0*(x-a)/((b-a)*(c-a));
	}
	else if ((c <= x) && (x < b)) {
		return 2.0*(b-x)/((b-a)*(b-c));
	}
	else {return 0.0;}


}

// fibril stress when there is rupturing but no yielding, required for ER model
double fibril_stress(const double &x, const double &colE, const double &lambdaC) {

	if (x >= lambdaC) {
		return colE*(x/lambdaC - 1.0);
	}
	else {return 0.0;}

}


// testing integral for the ER model
double integrand_test(const double &x, const double &t, const double &upper_lim, const double &lower_lim,
											const double &a, const double &b, const double &c,
	 										const double &colE) {

	//double lambdaC = (std::min(x,b)-std::max(a,x/lambdaR))*t/2.0 + (std::min(x,b)+std::max(a,x/lambdaR))/2.0;
	double lambdaC = 0.5*(upper_lim - lower_lim)*t + 0.5*(upper_lim + lower_lim);

	return fibril_stress(x, colE, lambdaC)*triangular_distribution(lambdaC, a, b, c)*0.5*(upper_lim - lower_lim);

}

//===============================================================
// distribution of fibril lengths in tendon
//===============================================================
class CrimpDistributionModel : public StrainEnergyFunction
{
public:

  CrimpDistributionModel(double* mu_pt, double* e_pt, double* phi_pt, double* a_pt, double* b_pt, double* c_pt) :
    StrainEnergyFunction(),
    Mu_pt(mu_pt), E_pt(e_pt), Phi_pt(phi_pt), A_pt(a_pt), B_pt(b_pt), C_pt(c_pt) {}

  /// Empty virtual constructor
  virtual ~CrimpDistributionModel(){}

  /// Return the strain energy in terms of the strain tensor
  double W(const DenseMatrix<double> &gamma)
  {return StrainEnergyFunction::W(gamma);}

  /// return the strain energy in terms of the strain invariants
  double W(const Vector<double> &I)
  {
	  return 0;
  }

  // return the derivatives of W wrt strain invariants
  void derivatives(Vector<double> &I, Vector<double> &dWdI)
  {
    dWdI[0] = (*Mu_pt)/2.0; //(1.0-(*Phi_pt))*(*Mu_pt)/2.0;
    dWdI[1] = 0.0;
    dWdI[2] = 0.0;

	// initialise dWdI[3]
	double temp = 0.0;
	// set the integration scheme
	Gauss<1,4> integration_scheme;

	for (unsigned i = 0; i < 4; i++){
		double t = integration_scheme.knot(i,0);
		if (sqrt(I[3]) <= *C_pt) {
			temp += integrand_test(sqrt(I[3]), t, sqrt(I[3]), 1, *A_pt, *B_pt, *C_pt, *E_pt)*integration_scheme.weight(i);
		}
		else if ((sqrt(I[3]) > *C_pt) && (sqrt(I[3]) <= *B_pt)) {
			temp += integrand_test(sqrt(I[3]), t, *C_pt, *A_pt, *A_pt, *B_pt, *C_pt, *E_pt)*integration_scheme.weight(i);
			temp += integrand_test(sqrt(I[3]), t, sqrt(I[3]), *C_pt, *A_pt, *B_pt, *C_pt, *E_pt)*integration_scheme.weight(i);
		}
		else if (sqrt(I[3]) > *B_pt){
			temp += integrand_test(sqrt(I[3]), t, *C_pt, *A_pt, *A_pt, *B_pt, *C_pt, *E_pt)*integration_scheme.weight(i);
			temp += integrand_test(sqrt(I[3]), t, *B_pt, *C_pt, *A_pt, *B_pt, *C_pt, *E_pt)*integration_scheme.weight(i);
		}
		else {
			temp += 0;
		}
	}

	dWdI[3] = temp/(2.0*I[3]); //*Phi_pt*temp/(2.0*I[3]);

   }



    bool requires_incompressibility_constraint(){return true;}

  private:

    double* Mu_pt;
    double* E_pt;
    double* Phi_pt;
    double* A_pt;
	double* B_pt;
	double* C_pt;

};

void transform_mesh_superposition(const Vector<double> &x_initial,
					 const double &s, const double &t1, const double &t2,
					 Vector<double> &x_new)
{
  // takes an initial coordinate in cartesian, transforms to cylindrical coordinates
  // apply shape transformation and then convert back to cartesian

  // transform initial position to cylindrical polar coordinates
  Vector<double> r_initial(3,0.0);
  // r = sqrt(x^2 + y^2)
  r_initial[0] = sqrt(pow(x_initial[0],2.0) + pow(x_initial[1],2.0));
  // theta = arctan(y/x)
  r_initial[1] = atan2(x_initial[1],x_initial[0]);
  // z remains unchanged
  r_initial[2] = x_initial[2];

  // apply the transformation to the position vector in cylindrical polars
  Vector<double> r_new(3,0.0);
  r_new[0] = r_initial[0]*((1.0 - t1)*r_initial[2]/Global_Parameters::length + t1)*(1.0 - (1.0 - t2)
  				*sin(MathematicalConstants::Pi*r_initial[2]/Global_Parameters::length));
  r_new[1] = r_initial[1];
  r_new[2] = r_initial[2];

  // transform back to cartesian
  x_new[0] = r_new[0]*cos(r_new[1]);
  if (Global_Parameters::constant_thickness){
	   x_new[1] = x_initial[1];//s*r_new[0]*sin(r_new[1]);
  }
  else{
	  x_new[1] = s*r_new[0]*sin(r_new[1]);
  }
  x_new[2] = r_new[2];


}



// function to work out fibre alignment for contour-aligned fibres
void contour_aligned_fibres_superposition(const Vector<double> &xi,
				 const double &s, const double &t1, const double &t2,
				 Vector<double> &M)
{
  // if the fibres are near the centre, they are aligned with the z-axis
  if (xi[0] == 0.0 && xi[1] == 0.0)
	{
  M[0] = 0.0;
  M[1] = 0.0;
  M[2] = 1.0;
	}
  else
	{
  double pi = 4.0*atan(1.0);
  double rp = sqrt(pow(xi[0],2.0) + pow(xi[1],2.0));
  double thetap = atan2(xi[1],xi[0]);

  // find ap, the semi major axis (at z=0) of the elliptic, tapered cylinder that our
  // point lies on

  double ap;
  double Mx;
  double My;

  if (Global_Parameters::constant_thickness){
	  ap = rp / sqrt(pow(cos(thetap)*(((1.0 - t1)*(xi[2]/Global_Parameters::length) + t1)*(1.0 - (1.0 - t2)
					  *sin(MathematicalConstants::Pi*xi[2]/Global_Parameters::length))),2.0)+ pow(s*sin(thetap),2.0));
					  // find the components of M, let Mz = 1
	  Mx = ap*cos(thetap)*((-pi/Global_Parameters::length)*(1.0-t2)*(t1 + (1.0-t1)*xi[2]/Global_Parameters::length)*cos(pi*xi[2]/Global_Parameters::length)
					  + (1.0/Global_Parameters::length)*(1.0-t1)*(1.0-(1.0-t2)*sin(pi*xi[2]/Global_Parameters::length)));
	  My = 0.0;
  }
  else{
	  ap = rp / ((((1.0 - t1)*(xi[2]/Global_Parameters::length) + t1)*(1.0 - (1.0 - t2)
					  *sin(MathematicalConstants::Pi*xi[2]/Global_Parameters::length)))
				  *sqrt(cos(thetap)*cos(thetap) + s*s*sin(thetap)*sin(thetap)));
      Mx = ap*cos(thetap)*((-pi/Global_Parameters::length)*(1.0-t2)*(t1 + (1.0-t1)*xi[2]/Global_Parameters::length)*cos(pi*xi[2]/Global_Parameters::length)
					 + (1.0/Global_Parameters::length)*(1.0-t1)*(1.0-(1.0-t2)*sin(pi*xi[2]/Global_Parameters::length)));
	  My = ap*s*sin(thetap)*((-pi/Global_Parameters::length)*(1.0-t2)*(t1 + (1.0-t1)*xi[2]/Global_Parameters::length)*cos(pi*xi[2]/Global_Parameters::length)
					  + (1.0/Global_Parameters::length)*(1.0-t1)*(1.0-(1.0-t2)*sin(pi*xi[2]/Global_Parameters::length)));
  }

  double Mz = 1.0;

  // find the unit vector
  double M2  = sqrt(pow(Mx,2.0) + pow(My,2.0) + pow(Mz, 2.0));

  M[0] = Mx/M2;
  M[1] = My/M2;
  M[2] = Mz/M2;
	}

}

double get_new_radius(const double &initial_radius, const double &s, const double &t1, const double &t2)
{
	double pi = 4.0*atan(1.0);
	// double v = (-8.0/pow(pi,3.0))*(t1 - 1.0)*(t1 - 1.0)*(t2 - 1.0)
	// 			+ (2.0/pi)*(1.0 + t1*t1)*(t2 - 1.0)
	// 			- (1.0/(4.0*pi*pi))*(t1 - 1.0)*(t1 - 1.0)*(t2 - 1.0)*(t2 - 1.0)
	// 			+ (1.0/6.0)*(1.0 + t1 + t1*t1)*(3.0 - 2.0*t2 + t2*t2);
    double v;
	if (Global_Parameters::constant_thickness){
		  v = Global_Parameters::length*(t1 + 1.0)*(2*t2 + pi - 2.0)/(2.0 * pi);
	  }
    else{
		  v = Global_Parameters::length*(-96.0*(t1 - 1.0)*(t1 - 1.0)*(t2 - 1.0)
	  				+ 24.0*pi*pi*(1.0 + t1*t1)*(t2 - 1.0)
	  				- 3.0*pi*(t1 - 1.0)*(t1 - 1.0)*(t2 - 1.0)*(t2 - 1.0)
	  				+ 2*pow(pi,3.0)*(1.0 + t1 + t1*t1)*(3.0 - 2.0*t2 + t2*t2))/(12.0*pow(pi,3.0));
	  }
	double new_radius = initial_radius/sqrt(s * v);
	return new_radius;
}

double transform_radius(const double &z, const double &t1, const double &t2, const double &s)
{
	return ((1.0 - t1)*z/Global_Parameters::length + t1)*(1.0 - (1.0 - t2)*sin(MathematicalConstants::Pi*z/Global_Parameters::length));
}

double phi_spatial(const double &z, const double &z0, const double &t1, const double &t2, const double &s, const double &phi)
{
	if (Global_Parameters::constant_thickness){
		return phi*transform_radius(z0, t1, t2, s)/transform_radius(z, t1, t2, s);
	}
	else{
		return phi*pow(transform_radius(z0, t1, t2, s), 2.0)/pow(transform_radius(z, t1, t2, s), 2.0);
	}

}
