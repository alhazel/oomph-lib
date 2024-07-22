
namespace Global_Parameters
{

  // output directory
 std::string dir;

  int model_selection = 2;

  // Holzapfel-Gasser-Ogden
  double C1 = 1.0;
  double K1 = 10.0;
  double K2 = 10.0;

  // crimped fibril model
  double mu = 1.0;
  double Ecol = 50.0;
  double a = 1.0;
  double b = 1.05;
  double c = 1.025;


  double s;
  double t1;
  double t2;


  // achilles:
  double A0 = 1.64;
  double length = 17;
  //double Phi = 0.65; // for constant phi
  double Phi = 0.85; // for non-constant phi

  // ACL
  // double A0 = 7.3;
  // double length = 24;
  // double Phi = 0.6; // for constant phi
  // //double Phi = 0.84; // for non-constant phi

  // number of layers in the initial mesh
  unsigned nlayer = 5;
  // no of times the bulk elements will be refined uniformly
  int no_uniform = 0;
  // no of times the top elements will be refined
  int no_boundary = 0;

  double initial_step_size = 0.005 * length;

  double max_strain = 0.1 * length;
  int counter;

  bool pin_top_surface = false;
  bool pin_bottom_surface = false;
  bool contour_fibres = true;
  bool match_volume = false;
  bool phi_sv = false;
  bool print_height = true;
  bool constant_thickness = false;

}
