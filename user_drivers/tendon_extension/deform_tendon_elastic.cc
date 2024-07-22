// code to ensure that the basic transvsersely isotropic cylinder
// is working as expected

//Generic routines

//Generic routines
#include "generic.h"
// Edited source code to pass xi through as a parameter
#include "solid.h"
// edited so that all five strain invariants are looped over
#include "constitutive.h"
// The mesh
#include "meshes/tube_mesh.h"
// a separate header file containing only material parameters
#include "parameters.h"


using namespace oomph;

// a header file which includes derived classes and functions which can be used anywhere
#include "MyDerivedClasses.h"
// include header file with the elastic-rupture constants in



// Global Objects
//=================================================================
namespace Global_Objects
{
  /// pointer to strain energy function
  StrainEnergyFunction* Strain_energy_function_pt;

  /// Pointer to constitutive law
  ConstitutiveLaw* Constitutive_law_pt;

  /// GeomObject specifying the shape of the boundary. (flat)
  //class DisplacedBoundary;//(const double& ampl);
  DisplacedBoundary Boundary_geom_object(0.0);

} // end of namespace


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


//===============================================================
/// Tube upgraded to become a solid mesh
//===============================================================
template<class ELEMENT>
class ElasticTubeMesh : public virtual TubeMesh<ELEMENT>,
			public virtual SolidMesh
{

 public:

  /// short constructor
  ElasticTubeMesh(GeomObject* volume_pt,
		  const Vector<double>& centreline_limits,
		  const Vector<double>& theta_positions,
		  const Vector<double>& radius_box,
		  const unsigned& nlayer,
		  TimeStepper* time_stepper_pt=
		  &Mesh::Default_TimeStepper) :
    TubeMesh<ELEMENT>(volume_pt,centreline_limits,theta_positions,radius_box,
		      nlayer,time_stepper_pt),
   SolidMesh()

  {
    // Assign the initial lagrangian coordinates
    set_lagrangian_nodal_coordinates();

    Vector<double> zeta(2);
    /// Loop over the top boundary
    {
      unsigned b=2;
      unsigned n_node = this->nboundary_node(b);
      for (unsigned n=0;n<n_node;n++)
    	{
	  SolidNode* nod_pt=static_cast<SolidNode*>(this->boundary_node_pt(b,n));
	  zeta[0] = nod_pt->x(0);
    	  zeta[1] = nod_pt->x(1);
    	  nod_pt->set_coordinates_on_boundary(b,zeta);
    	}
      this->Boundary_coordinate_exists[b]=true;
    }


  }


  /// Empty Destructor
  virtual ~ElasticTubeMesh() { }

};

//===============================================================
/// Tube upgraded to become a refineable solid mesh
//===============================================================
template<class ELEMENT>
class RefineableElasticTubeMesh : public virtual  RefineableTubeMesh<ELEMENT>,
				  public virtual  SolidMesh
{

 public:

  /// short constructor
  RefineableElasticTubeMesh(GeomObject* volume_pt,
			    const Vector<double>& centreline_limits,
			    const Vector<double>& theta_positions,
			    const Vector<double>& radius_box,
			    const unsigned& nlayer,
			    TimeStepper* time_stepper_pt=
			    &Mesh::Default_TimeStepper) :
    RefineableTubeMesh<ELEMENT>(volume_pt,centreline_limits,theta_positions,radius_box,
				nlayer,time_stepper_pt),
   SolidMesh()
  {

    this->setup_octree_forest();

    /// Assign the initial lagrangian coordinates
    set_lagrangian_nodal_coordinates();

    Vector<double> zeta(2);
    /// Loop over the top boundary
    {
      unsigned b=2;
      unsigned n_node = this->nboundary_node(b);
      for (unsigned n=0;n<n_node;n++)
    	{
	  Node* nod_pt=this->boundary_node_pt(b,n);
	  zeta[0] = nod_pt->x(0);
    	  zeta[1] = nod_pt->x(1);
    	  nod_pt->set_coordinates_on_boundary(b,zeta);
    	}
      this->Boundary_coordinate_exists[b]=true;
    }

  }

  /// Empty Destructor
  virtual ~RefineableElasticTubeMesh() { }

};


//=start_of_MyCylinder===================================
//A geometric object that represents the geometry of the domain
//================================================================
class MyEllipticalCylinder : public GeomObject
{
public:

 /// Constructor that takes the radius of the tube as its argument
 MyEllipticalCylinder(const double &radius1, const double &radius2) :
  GeomObject(3,3), Radius1(radius1), Radius2(radius2) { }

/// Destructor
virtual~MyEllipticalCylinder(){}

///Lagrangian coordinate xi
void position (const Vector<double>& xi, Vector<double>& r) const
{
 r[0] = xi[2]*Radius1*cos(xi[1]);
 r[1] = xi[2]*Radius2*sin(xi[1]);
 r[2] = xi[0];}


/// Return the position of the tube as a function of time
/// (doesn't move as a function of time)
void position(const unsigned& t,
              const Vector<double>& xi, Vector<double>& r) const
  {
   position(xi,r);
  }

private:

 ///Storage for the radius of the tube
 double Radius1;
 double Radius2;

};






//==================================================================
// Deformation of elastic tube
//==================================================================
template<class ELEMENT>
class CylinderExtensionProblem : public Problem
{

public:

  /// Constructor:
  CylinderExtensionProblem();

  // Run simulation
  void run(const std::string &dirname);

  /// doc the solution
  void doc_solution(DocInfo& doc_info);

	// dump solution
 void dump_it(std::ofstream& dump_file);

	// restart fron dumped solution
 void restart(std::ifstream& restart_file);

  /// update function (empty)
  void actions_after_newton_solve() {}

 #ifdef NOREFINE

  /// Access function for the solid mesh
  ElasticTubeMesh<ELEMENT>*& solid_mesh_pt()
  {return Solid_mesh_pt;}

 #else

  /// Access function for the refineable solid mesh
  RefineableElasticTubeMesh<ELEMENT>*& solid_mesh_pt()
  {return Solid_mesh_pt;}

 #endif

 #ifndef NOREFINE

  /// Actions before adapt
  void actions_before_adapt()
  {
    /// kill the elements and wipe the surface mesh
    delete_lagrange_multiplier_elements();
    /// Rebuild the problems global mesh from its various sub meshes
    rebuild_global_mesh();

  } // end of actions_before_adapt

  /// Actions after adapt
  void actions_after_adapt()
  {
	  std::cout << "setting BCS!!!!!!!\n" << std::endl;
		/// setup boundary conditions
	  {
	    unsigned b=0;
	    unsigned n_node = Solid_mesh_pt->nboundary_node(b);


			if (Global_Parameters::pin_bottom_surface){
				for(unsigned n=0;n<n_node;n++)
				{
					// Pin all nodes:
					//std::cout << "pinning boundary nodes in all directions" << std::endl;
					Solid_mesh_pt->boundary_node_pt(b,n)->pin_position(0);
					Solid_mesh_pt->boundary_node_pt(b,n)->pin_position(1);
					Solid_mesh_pt->boundary_node_pt(b,n)->pin_position(2);
				}
			}
			else {
				for(unsigned n=0;n<n_node;n++){
					// Pin nodes in z-direction and to prevent rotations
					// pin only in z direction and to prevent rotations and translations:
					//std::cout << "only pinning boundary nodes in z-direction" << std::endl;
					const double tolerance = 1e-10;
					// pin in the z-direction
					Solid_mesh_pt->boundary_node_pt(b,n)->pin_position(2);
					// if the node lies on the x axis, pin in y direction
					if (((Solid_mesh_pt->boundary_node_pt(b,n)->x(0)) < tolerance)
					&& ((Solid_mesh_pt->boundary_node_pt(b,n)->x(0)) > (-tolerance)))
					{
						Solid_mesh_pt->boundary_node_pt(b,n)->pin_position(0);
					}
					// if the node lies on the y axis, pin in the x direction
					if (((Solid_mesh_pt->boundary_node_pt(b,n)->x(1)) <  tolerance)
					&& ((Solid_mesh_pt->boundary_node_pt(b,n)->x(1)) > (-tolerance)))
					{
						Solid_mesh_pt->boundary_node_pt(b,n)->pin_position(1);
					}
				}

			}

	  } // end of setup boundary conditions

    //Sort out the boundary coordinate
    /// Loop over the top boundary
    {
      Vector<double> zeta(2);
      unsigned b=2;
      unsigned n_node = solid_mesh_pt()->nboundary_node(b);
      for (unsigned n=0;n<n_node;n++)
    	{
	  SolidNode* nod_pt=static_cast<SolidNode*>(solid_mesh_pt()->boundary_node_pt(b,n));
	  zeta[0] = nod_pt->x(0);
    	  zeta[1] = nod_pt->x(1);
    	  nod_pt->set_coordinates_on_boundary(b,zeta);
    	}
    }

    /// Create the elements that impose the displaced boundary constraint
    /// and attach them to the bulk elements that are adjacent to boundary 2
    create_lagrange_multiplier_elements();

    /// Rebuild the problem's global mesh from its various sub-meshes
    rebuild_global_mesh();

    /// Pin the redundant solid pressures // PVDEquationsBase
    PVDEquationsBase<3>::pin_redundant_nodal_solid_pressures(Solid_mesh_pt->element_pt());


  } // end of  after adapt

 #endif



  void refine_boundary_elements()
  {
		// number of elements in the mesh
  unsigned no_of_element = Solid_mesh_pt->nelement();

	// initialise two vectors of elements to be refined
  Vector<unsigned> top_element_vector(0);
	Vector<unsigned> bottom_element_vector(0);

  // loop  over all the elements
  for(unsigned e=0;e<no_of_element;e++)
    {
      // loop over the nodes of each element
      unsigned no_nodes = Solid_mesh_pt->finite_element_pt(e)->nnode();
      for (unsigned i=0; i<no_nodes; i++)
			{
					// if element is on the top boundary, add to vector
				  if(Solid_mesh_pt->finite_element_pt(e)->node_pt(i)->is_on_boundary(2))
				    {
				      //std::cout << "element number " << e << " is on boundary 2" << std::endl;
				      top_element_vector.push_back(e);
				    }
					if(Solid_mesh_pt->finite_element_pt(e)->node_pt(i)->is_on_boundary(0))
				    {
				      //std::cout << "element number " << e << " is on boundary 2" << std::endl;
				      bottom_element_vector.push_back(e);
				    }
			}
    }
  // vector should already be sorted, delete repeated elements
 top_element_vector.erase(unique(top_element_vector.begin(),top_element_vector.end()),top_element_vector.end());
 bottom_element_vector.erase(unique(bottom_element_vector.begin(),bottom_element_vector.end()),bottom_element_vector.end());

// refine elements in both vectors
 Solid_mesh_pt->refine_selected_elements(top_element_vector);
 Solid_mesh_pt->refine_selected_elements(bottom_element_vector);
 std::cout << "no of elements: " << Solid_mesh_pt->nelement() << std::endl;
  }





private:

  /// Create elements that enforce prescribed boundary motion
  void create_lagrange_multiplier_elements();

  /// delete elements that enforce prescribed boundary motion
  void delete_lagrange_multiplier_elements();

  /// updated before solve: empty
  void actions_before_newton_solve(){}

	void actions_before_newton_convergence_check(){

		Solid_mesh_pt->node_update();
	}





 #ifdef NOREFINE

  /// pointer to solid mesh
  ElasticTubeMesh<ELEMENT>* Solid_mesh_pt;

 #else

  /// pointer to solid mesh
  RefineableElasticTubeMesh<ELEMENT>* Solid_mesh_pt;

 #endif

  /// pointers to meshes of lagrange multiplier elements
  SolidMesh* Lagrange_multiplier_mesh_pt;

  /// pointer to GeomObject that specifies the domain volume
  GeomObject *Volume_pt;

  /// Docinfo object for output
  DocInfo Doc_info;

};

//==============================================================
// Creates Lagrange multiplier elements
//==============================================================
template<class ELEMENT>
void CylinderExtensionProblem<ELEMENT>::
create_lagrange_multiplier_elements()
{
  /// Lagrange multiplier elements are located on boundary 2
  unsigned b=2;

  // How many bulk elements are adjacent to boundary b?
  unsigned n_element = solid_mesh_pt()->nboundary_element(b);

  // loop over the bulk elements adjacent to boundary b
  for(unsigned e=0;e<n_element;e++)
    {
      // Get pointer to the bulk element that is adjacent to boundary b
      ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(solid_mesh_pt()->boundary_element_pt(b,e));

      // find the index of the face of element e along boundary b
      int face_index = solid_mesh_pt()->face_index_at_boundary(b,e);


      /// create new element and add to mesh
      Lagrange_multiplier_mesh_pt->add_element_pt( new MyLagrangeMultiplierElement<ELEMENT>(bulk_elem_pt,face_index));

    }

  // Loop over the elements in the lagrange multiplier element mesh

  n_element = Lagrange_multiplier_mesh_pt->nelement();
  for(unsigned i=0;i<n_element;i++)
    {

      ///cast to lagrange multiplier element
      MyLagrangeMultiplierElement<ELEMENT> *el_pt =
	dynamic_cast<MyLagrangeMultiplierElement<ELEMENT>*>
	(Lagrange_multiplier_mesh_pt->element_pt(i));

      /// Set the GeomObject that defines the boundary shape and
      /// specify which bulk boundary we are attached to ( need to extract
      /// the boundary coordinate from the bulk nodes)
      el_pt->set_boundary_shape_geom_object_pt(&Global_Objects::Boundary_geom_object,b);

		}




}// end of create_lagrane_multiplier_elements

///===============================================================
/// Deletes lagrange multiplier elements
///===============================================================
template<class ELEMENT>
void CylinderExtensionProblem<ELEMENT>::delete_lagrange_multiplier_elements()
{
  /// How many surface elements are in the surface mesh?
  unsigned n_element = Lagrange_multiplier_mesh_pt->nelement();

  /// Loop over the surface elements
  for(unsigned e=0;e<n_element;e++)
    {
      /// kill the surface element
      delete Lagrange_multiplier_mesh_pt->element_pt(e);
    }

  /// wipe the mesh
  Lagrange_multiplier_mesh_pt->flush_element_and_node_storage();

} // end of delete_lagrange_multiplier_elements

///===============================================================
/// Constructor
///===============================================================
template<class ELEMENT>
CylinderExtensionProblem<ELEMENT>::CylinderExtensionProblem()
{
  /// Create GeomObject that specifies the domain geometry
  /// the radius of the reference cylinder

  //Volume_pt = new MyCylinder(Global_Parameters::A0);
  //Volume_pt = new EllipticalTube(Global_Parameters::A0, Global_Parameters::s*Global_Parameters::A0);
  double A;
  if (Global_Parameters::match_volume)
  {
	  A = get_new_radius(Global_Parameters::A0, Global_Parameters::s, Global_Parameters::t1, Global_Parameters::t2);
  }
  else
  {
	  A = Global_Parameters::A0;
  }
  std::cout << "top radius = " << A << std::endl;
  Volume_pt = new MyEllipticalCylinder(A, Global_Parameters::s*A);

  /// Define pi
  const double pi = 4.0*atan(1.0);


  /// Set the centreline coordinates spanning the mesh
  Vector<double> centreline_limits(2);
  centreline_limits[0] = 0.0;
  centreline_limits[1] = Global_Parameters::length;

  /// Set the positions of the angles that divide the outer ring
  Vector<double> theta_positions(4);
  theta_positions[0] = -0.75*pi;
  theta_positions[1] = -0.25*pi;
  theta_positions[2] = 0.25*pi;
  theta_positions[3] = 0.75*pi;






  /// Define the radial fraction of the central box
  Vector<double> radial_frac(4,0.5);

#ifdef NOREFINE


  /// now create the mesh
  Solid_mesh_pt = new ElasticTubeMesh<ELEMENT>
    (Volume_pt,centreline_limits,theta_positions,radial_frac,Global_Parameters::nlayer);


#else
  /// now create the mesh
  Solid_mesh_pt = new RefineableElasticTubeMesh<ELEMENT>
    (Volume_pt,centreline_limits,theta_positions,radial_frac,Global_Parameters::nlayer);

  // /// Set error estimator
  Solid_mesh_pt->spatial_error_estimator_pt() = new Z2ErrorEstimator;
  // Solid_mesh_pt->max_permitted_error() = 0.002;
  // Solid_mesh_pt->min_permitted_error() = 0.0001;
  Solid_mesh_pt->max_permitted_error() = 0.0;
  Solid_mesh_pt->min_permitted_error() = 0.0;

#endif

// transform cylinder into new shape
// loop over the nodes
unsigned n_node_cylinder=Solid_mesh_pt->nnode();
std::cout << "number of nodes in the cylinder" << n_node_cylinder << std::endl;
for(unsigned i=0;i<n_node_cylinder;i++)
  {
	// initial position vector of the node
	Vector<double> x_initial(3,0.0);
	// transformed coordinates
	Vector<double> x_new(3,0.0);
	x_initial[0] = Solid_mesh_pt->node_pt(i)->x(0);
	x_initial[1] = Solid_mesh_pt->node_pt(i)->x(1);
	x_initial[2] = Solid_mesh_pt->node_pt(i)->x(2);

	// use the global function to convert the coordinates as required
	transform_mesh_superposition(x_initial, 1.0,
	 				Global_Parameters::t1, Global_Parameters::t2, x_new);
  // set the coordinates of the mesh to the new transformed coordinates
  Solid_mesh_pt->node_pt(i)->x(0) = x_new[0];
  Solid_mesh_pt->node_pt(i)->x(1) = x_new[1];
  Solid_mesh_pt->node_pt(i)->x(2) = x_new[2];
  }

  // set these positions as the unstrained node positions -- important when
  // applying transformations to the mesh
  Solid_mesh_pt->set_lagrangian_nodal_coordinates();
  /// loop over the elements in the mesh to set parameters/functions pointers
  unsigned n_element=Solid_mesh_pt->nelement();
  for(unsigned i=0;i<n_element;i++)
    {
      /// cast to a solid element
      ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Solid_mesh_pt->element_pt(i));

      /// set the constitutive law
      el_pt->constitutive_law_pt() =
	Global_Objects::Constitutive_law_pt;

      /// Material is incompressible
      PVDEquationsWithPressure<3>* cast_el_pt =
	dynamic_cast<PVDEquationsWithPressure<3>*>(Solid_mesh_pt->element_pt(i));
      if (cast_el_pt!=0)
	{
	  cast_el_pt->set_incompressible();
	}

    }


  // REFINE ELEMENTS
  // refine all or just the top elements a specified number of times
  for(int i=0; i<Global_Parameters::no_uniform; i++)
    {
      Solid_mesh_pt->refine_uniformly();
    }
  for(int i=0; i<Global_Parameters::no_boundary; i++)
    {
      refine_boundary_elements();
    }

  std::cout << "number of elements: "<< Solid_mesh_pt->nelement() << std::endl;
	std::cout << "number of nodes: "<< Solid_mesh_pt->nnode() << std::endl;

  /// Construct the mesh of elements that enforce prescribed boundary motion
  Lagrange_multiplier_mesh_pt = new SolidMesh;
  create_lagrange_multiplier_elements();


  /// Solid mesh is the first sub-mesh
  add_sub_mesh(Solid_mesh_pt);

  /// add lagrange multiplier sub-mesh
  add_sub_mesh(Lagrange_multiplier_mesh_pt);

  /// build combined global mesh
  build_global_mesh();

  Solid_mesh_pt->disable_adaptation();
  /// setup boundary conditions
  {
		unsigned b=0;
		unsigned n_node = Solid_mesh_pt->nboundary_node(b);


		if (Global_Parameters::pin_bottom_surface){
			for(unsigned n=0;n<n_node;n++)
			{
				// Pin all nodes:
				//std::cout << "pinning boundary nodes in all directions" << std::endl;
				Solid_mesh_pt->boundary_node_pt(b,n)->pin_position(0);
				Solid_mesh_pt->boundary_node_pt(b,n)->pin_position(1);
				Solid_mesh_pt->boundary_node_pt(b,n)->pin_position(2);

 			}
		}
		else {
			for(unsigned n=0;n<n_node;n++){
				// Pin nodes in z-direction and to prevent rotations
				// pin only in z direction and to prevent rotations and translations:
				//std::cout << "only pinning boundary nodes in z-direction" << std::endl;
				const double tolerance = 1e-10;
				// pin in the z-direction
				Solid_mesh_pt->boundary_node_pt(b,n)->pin_position(2);
				// if the node lies on the x axis, pin in y direction
				if (((Solid_mesh_pt->boundary_node_pt(b,n)->x(0)) < tolerance)
				&& ((Solid_mesh_pt->boundary_node_pt(b,n)->x(0)) > (-tolerance)))
				{
					Solid_mesh_pt->boundary_node_pt(b,n)->pin_position(0);
				}
				// if the node lies on the y axis, pin in the x direction
				if (((Solid_mesh_pt->boundary_node_pt(b,n)->x(1)) <  tolerance)
				&& ((Solid_mesh_pt->boundary_node_pt(b,n)->x(1)) > (-tolerance)))
				{
					Solid_mesh_pt->boundary_node_pt(b,n)->pin_position(1);
				}
			}

		}
		// for(unsigned n=0;n<n_node;n++){
		// 	// Pin nodes in z-direction and to prevent rotations
		// 	// pin only in z direction and to prevent rotations and translations:
		// 	//std::cout << "only pinning boundary nodes in z-direction" << std::endl;
		// 	const double tolerance = 1e-10;
		// 	// pin in the z-direction
		// 	Solid_mesh_pt->boundary_node_pt(b,n)->pin_position(2);
		// 	if the node lies on the x axis, pin in y direction
		// 	if (((Solid_mesh_pt->boundary_node_pt(b,n)->x(0)) < tolerance)
		// 	&& ((Solid_mesh_pt->boundary_node_pt(b,n)->x(0)) > (-tolerance))
		// 	&& ((Solid_mesh_pt->boundary_node_pt(b,n)->x(1)) <  tolerance)
		// 	&& ((Solid_mesh_pt->boundary_node_pt(b,n)->x(1)) > (-tolerance)))
		// 	{
		// 		Solid_mesh_pt->boundary_node_pt(b,n)->pin_position(0);
		// 		Solid_mesh_pt->boundary_node_pt(b,n)->pin_position(1);
		// 	}
		// }

	} // end of setup boundary conditions

  /// pin the redundant solid pressures
  PVDEquationsBase<3>::pin_redundant_nodal_solid_pressures(Solid_mesh_pt->element_pt());

  /// attach the boundary conditions to the mesh
  std::cout << "Number of equations: " << assign_eqn_numbers() << std::endl;


}

///==========================================================
/// Doc the solution
///==========================================================
template<class ELEMENT>
void CylinderExtensionProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{
 std::ofstream some_file;
  char filename[150];

  /// number of plot points
  unsigned npts = 5;

  ///  the shape  of the deformed boday
  sprintf(filename, "%s/soln%s.dat", doc_info.directory().c_str(),
	  doc_info.label().c_str());
  some_file.open(filename);
  Solid_mesh_pt->output(some_file,npts);
  some_file.close();

  // output information about the surface elements on the top of the tendon
  sprintf(filename, "%s/resultant_force%s.dat",doc_info.directory().c_str(),
	  doc_info.label().c_str());
  some_file.open(filename);
  Lagrange_multiplier_mesh_pt->output(some_file,npts);
  some_file.close();

	// // Write restart file
 // sprintf(filename,"%s/restart%i.dat",doc_info.directory().c_str(),
 //         doc_info.number());
 // some_file.open(filename);
 // dump_it(some_file);
 // some_file.close();



}

//=====start_of_dump_it===================================================
/// Dump the solution to disk to allow for restart
//========================================================================
template<class ELEMENT>
void CylinderExtensionProblem<ELEMENT>::dump_it(std::ofstream& dump_file)
{
 // Call generic dump()
 Problem::dump(dump_file);
} // end of dump_it


//=================start_of_restart=======================================
/// Read solution from disk for restart
//========================================================================
template<class ELEMENT>
void CylinderExtensionProblem<ELEMENT>::restart(std::ifstream& restart_file)
{
 // Read the generic problem data from restart file
 Problem::read(restart_file);
} // end of restart


//======================================================
/// Driver for simple elastic problem
//======================================================
int main(int argc, char **argv)
{

	// store command line arguments
	CommandLineArgs::setup(argc,argv);
	// t1
	CommandLineArgs::specify_command_line_flag("--t1",&Global_Parameters::t1);
	// t2
	CommandLineArgs::specify_command_line_flag("--t2",&Global_Parameters::t2);
	// ellipticity
	CommandLineArgs::specify_command_line_flag("--s",&Global_Parameters::s);
	// counter for the naming of files
	CommandLineArgs::specify_command_line_flag("--counter",&Global_Parameters::counter);
	// directory name
	CommandLineArgs::specify_command_line_flag("--dir",&Global_Parameters::dir);
	// Parse command line
	CommandLineArgs::parse_and_assign();
	// Doc what has actually been specified on the command line
	CommandLineArgs::doc_specified_flags();





	// create an output file so that data is more organised
        std::ofstream metadata;
	metadata.open(Global_Parameters::dir + "/parameters.txt");
	if (Global_Parameters::match_volume)
    {
  	  metadata << "top radius: " << get_new_radius(Global_Parameters::A0, Global_Parameters::s, Global_Parameters::t1, Global_Parameters::t2);
    }
    else
    {
      metadata << "top radius: " << Global_Parameters::A0 << std::endl;
    }
	metadata << "s: " << Global_Parameters::s << std::endl;
	metadata << "t1: " << Global_Parameters::t1 << std::endl;
	metadata << "t2: " << Global_Parameters::t2 << std::endl;
	metadata << "nlayer: " << Global_Parameters::nlayer << std::endl;
	metadata << "initial step size: " << Global_Parameters::initial_step_size << std::endl;
	metadata << "no x uniform refinement: " << Global_Parameters::no_uniform << std::endl;
	metadata << "no x boundary refinement: " << Global_Parameters::no_boundary << "\n" << std::endl;

	if (Global_Parameters::model_selection == 1)
	{
		// transversely isotropic HGO model
		Global_Objects::Strain_energy_function_pt =
		new HolzapfelGasserOgdenSEF(&Global_Parameters::C1,
		&Global_Parameters::K1,
		&Global_Parameters::K2);

		// output pars
		metadata << "C1: " << Global_Parameters::C1 << std::endl;
		metadata << "K1: " << Global_Parameters::K1 << std::endl;
		metadata << "K2: " << Global_Parameters::K2 << std::endl;
	}
	else if (Global_Parameters::model_selection == 2)
	{
		// crimped fibril model
		Global_Objects::Strain_energy_function_pt =
		new CrimpDistributionModel(&Global_Parameters::mu,
		&Global_Parameters::Ecol,
		&Global_Parameters::Phi,
		&Global_Parameters::a,
		&Global_Parameters::b,
		&Global_Parameters::c);

		// output pars
		metadata << "mu: " << Global_Parameters::mu << std::endl;
		metadata << "E: " << Global_Parameters::Ecol << std::endl;
		metadata << "phi: " << Global_Parameters::Phi << std::endl;
		metadata << "a: " << Global_Parameters::a << std::endl;
		metadata << "b: " << Global_Parameters::b << std::endl;
		metadata << "c: " << Global_Parameters::c << std::endl;
	}

	metadata.close();

	Global_Objects::Constitutive_law_pt =
	new TransverselyIsotropicStrainEnergyFunctionConstitutiveLaw(
	Global_Objects::Strain_energy_function_pt);

	// create a file name and label using the ellipticity and taper
        std::string filename = Global_Parameters::dir + "/data";
	std::cout << filename << std::endl;
	#ifdef NOREFINE
	///set up the problem with pressure/displacement formulation
	CylinderExtensionProblem<MyQPVDElementWithContinuousPressure<3> > problem;
	//problem.run(filename);

	#else
	/// set up the problem with pressure/displacement formulation
	CylinderExtensionProblem<MyRefineableQPVDElementWithContinuousPressure<3> > problem;
	//problem.run(filename);

	#endif

	// output
	DocInfo doc_info;
	// set output director
	doc_info.set_directory(filename);


	const unsigned n_dof = problem.ndof();
	Vector<double> solution_backup(n_dof);
	Vector<double> base_solution_backup(n_dof);
	for(unsigned n=0;n<n_dof;n++) {solution_backup[n] = problem.dof(n);}


	// output initial nodal positions
        std::ofstream nodal_position;
	nodal_position.open(Global_Parameters::dir + "/initial_nodal_positions.txt");
	int nnodes = problem.solid_mesh_pt()->nnode();
	for(int i = 0; i < nnodes; i++){
		nodal_position << problem.solid_mesh_pt()->node_pt(i)->x(0) << " "
		<< problem.solid_mesh_pt()->node_pt(i)->x(1) << " "
		<< problem.solid_mesh_pt()->node_pt(i)->x(2) << " "
		<< problem.solid_mesh_pt()->node_pt(i)->position_is_pinned(0) << " "
		<< problem.solid_mesh_pt()->node_pt(i)->position_is_pinned(1) << " "
		<< problem.solid_mesh_pt()->node_pt(i)->position_is_pinned(2) << std::endl;
	}

	/// Initial parameter value
	Global_Objects::Boundary_geom_object.ampl()=0.0;

	double step_size = Global_Parameters::initial_step_size;


	// integers for labelling output
	int stretch_no = 0;


	/// doc initial configuration
	doc_info.label() = to_string(stretch_no);
	problem.doc_solution(doc_info);

	while(Global_Objects::Boundary_geom_object.ampl() < Global_Parameters::max_strain)
	{
		// print the output label
		std::cout << "DOC LABEL = " << doc_info.label() << std::endl;

		Global_Objects::Boundary_geom_object.ampl() += step_size;
		std::cout << "BOUNDARY DISPLACEMENT = " << Global_Objects::Boundary_geom_object.ampl() << std::endl;
		stretch_no ++;

		// int nnodes = problem.solid_mesh_pt()->nnode();
		// for(int i = 0; i < nnodes; i++){
		// 	if (problem.solid_mesh_pt()->node_pt(i)->is_on_boundary(2)){
		// 		std::cout << "height of boundary = " << problem.solid_mesh_pt()->node_pt(i)->x(2);
		// 	}
		// }


		// break if step size becomes too small
		if(abs(step_size) < 1e-5)
		{
			break;
		}

		try{
			problem.newton_solve(0);
			// if(step_size < Global_Parameters::initial_step_size){
			// 	step_size *= 2.0;
			// }
			std::cout << "solve successful. boundary displacement = " << Global_Objects::Boundary_geom_object.ampl() << std::endl;
		}
		catch (OomphLibError& error){
			std::cout << "ERROR! DECREASING STEP SIZE. NEW STEP SIZE = " << 0.5*step_size << std::endl;
			// restore parameter value
			Global_Objects::Boundary_geom_object.ampl() -= step_size;
			stretch_no --;
			step_size *= 0.5;
			//restore from backup
			for(unsigned n=0;n<n_dof;n++){problem.dof(n) = base_solution_backup[n];}
			continue;
		}

		// save the base solution and set the solution backup to the base solution
		for(unsigned n=0;n<n_dof;n++)
		{
			base_solution_backup[n] = problem.dof(n);
		}

		doc_info.label() = to_string(stretch_no);
		std::cout << doc_info.label() << std::endl;
		std::cout << "DOCUMENTING SOLUTION" << std::endl;
		problem.doc_solution(doc_info);
		Global_Parameters::print_height = true;
	}




}
