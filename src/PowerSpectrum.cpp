#include"PowerSpectrum.h"

//====================================================
// Constructors
//====================================================

PowerSpectrum::PowerSpectrum(
    BackgroundCosmology *cosmo,
    RecombinationHistory *rec,
    Perturbations *pert,
    double A_s,
    double n_s,
    double kpivot_mpc) :
  cosmo(cosmo),
  rec(rec),
  pert(pert),
  A_s(A_s),
  n_s(n_s),
  kpivot_mpc(kpivot_mpc)
{}

//====================================================
// Do all the solving
//====================================================
void PowerSpectrum::solve(){

  Vector k_array = Utils::linspace(k_min, k_max, n_k);
  Vector log_k_array = log(k_array);

  generate_bessel_function_splines();

  //=========================================================================
  // TODO: Line of sight integration to get Theta_ell(k)
  // Implement line_of_sight_integration
  //=========================================================================
  line_of_sight_integration(k_array);

  //=========================================================================
  // TODO: Integration to get Cell by solving dCell^f/dlogk = Delta(k) * f_ell(k)^2
  // Implement solve_for_cell
  //=========================================================================
  auto cell_TT = solve_for_cell(log_k_array, thetaT_ell_of_k_spline, thetaT_ell_of_k_spline);
  cell_TT_spline.create(ells, cell_TT, "Cell_TT_of_ell");

  //=========================================================================
  // TODO: Do the same for polarization...
  //=========================================================================
  // ...
  // ...
  // ...
  // ...
}

//====================================================
// Generate splines of j_ell(z) needed for LOS integration
//====================================================

void PowerSpectrum::generate_bessel_function_splines(){
  Utils::StartTiming("besselspline");

  // Make storage for the splines
  j_ell_splines = std::vector<Spline>(ells.size());

  //=============================================================================
  //Callin has an algorithm for choosing an appropriate x range based on ell
  //but for now I will use the one from page 5 of the powerpoint in the milestone description
  //and just see how it works

  //on second thought I'm not sure the algorithm in Callin actually refers to this
  //and not the integration to find C_ell
  //=============================================================================

  double eta_0 = cosmo->eta_of_x(0.0);
  double x_max = k_max*eta_0;
  int n_x = x_max*16/(2*M_PI);
  Vector x_bessel = Utils::linspace(0.0, x_max, n_x);

  //find the Bessel function
  Vector j_ell(n_x);

  for(size_t i = 0; i < ells.size(); i++){
    const int ell = ells[i];
    for(int ix = 0; ix < n_x; ix++){
      const double x = x_bessel[ix];
      j_ell[ix] = Utils::j_ell(ell, x);
    }

    // Make the j_ell_splines[i] spline
    j_ell_splines[i].create(x_bessel, j_ell, "j_ell_spline");
  }

  Utils::EndTiming("besselspline");
}

//====================================================
// Do the line of sight integration for a single
// source function
//====================================================

Vector2D PowerSpectrum::line_of_sight_integration_single(
    Vector & k_array,
    std::function<double(double,double)> &source_function){
  Utils::StartTiming("lineofsight");

  // Make storage for the results
  Vector2D result = Vector2D(ells.size(), Vector(k_array.size()));

  // Set up x array to integrate over
  double delta_x = 2.*M_PI/6.;
  int nx         = abs(0.0 - Constants.x_start)/delta_x;
  Vector x_array = Utils::linspace(Constants.x_start, 0.0, nx);

  // For the Bessel functions
  double eta_start = cosmo->eta_of_x(Constants.x_start);
  double eta_0     = cosmo->eta_of_x(0.0);

  double Theta;
  for(size_t ik = 0; ik < k_array.size(); ik++){
    const double k = k_array[ik];
    //=============================================================================
    // Implement to solve for the general line of sight integral
    // Theta_ell(k) = Int dx jell(k(eta0-eta)) * S(x,k) for all ell for given k
    //=============================================================================
    for(size_t il = 0; il < ells.size(); il++){
      // Trapezoidal integral over x:
      Theta = source_function(x_array[0], k)*j_ell_splines[il](k*(eta_0 - eta_start))
              + source_function(x_array[nx-1], k)*j_ell_splines[il](0.0);

      for(int ix = 1; ix < nx-2; ix++){
        const double x = x_array[ix];
        double eta_x   = cosmo->eta_of_x(x);

        Theta += 2.*source_function(x_array[ix], k)*j_ell_splines[il](k*(eta_0 - eta_x));
      }

      Theta = Theta*delta_x/2.;

      // Store the result for Theta_ell(k)
      result[il][ik] = Theta;
    }
  }

  Utils::EndTiming("lineofsight");
  return result;
}

//====================================================
// Do the line of sight integration
//====================================================
void PowerSpectrum::line_of_sight_integration(Vector & k_array){
  const int n_k        = k_array.size();
  const int n          = 100;
  const int nells      = ells.size();

  // Make storage for the splines we are to create
  thetaT_ell_of_k_spline = std::vector<Spline>(nells);

  //============================================================================
  // TODO: Solve for Theta_ell(k) and spline the result
  //============================================================================

  // Make a function returning the source function
  std::function<double(double,double)> source_function_T = [&](double x, double k){
    return pert->get_Source_T(x,k);
  };

  // Do the line of sight integration
  Vector2D thetaT_ell_of_k = line_of_sight_integration_single(k_array, source_function_T);

  // Spline the result and store it in thetaT_ell_of_k_spline
  for (size_t il = 0; il < ells.size(); il++){
    thetaT_ell_of_k_spline[il].create(k_array, thetaT_ell_of_k[il], "ThetaT_ell_of_k");
  }

  if(Constants.polarization){
    // not included
  }
}

//====================================================
// Compute Cell (could be TT or TE or EE)
// Cell = Int_0^inf 4 * pi * P(k) f_ell g_ell dk/k
//====================================================
Vector PowerSpectrum::solve_for_cell(
    Vector & log_k_array,
    std::vector<Spline> & f_ell_spline,
    std::vector<Spline> & g_ell_spline){
  const int nells      = ells.size();

  //============================================================================
  // TODO: Integrate Cell = Int 4 * pi * P(k) f_ell g_ell dk/k
  // or equivalently solve the ODE system dCell/dlogk = 4 * pi * P(k) * f_ell * g_ell
  //============================================================================

  // ...
  // ...
  // ...
  // ...

  Vector result;

  return result;
}

//====================================================
// The primordial power-spectrum
//====================================================

double PowerSpectrum::primordial_power_spectrum(const double k) const{
  return A_s * pow( Constants.Mpc * k / kpivot_mpc , n_s - 1.0);
}

//====================================================
// P(k) in units of (Mpc)^3
//====================================================

double PowerSpectrum::get_matter_power_spectrum(const double x, const double k_mpc) const{
  double pofk = 0.0;

  //=============================================================================
  // TODO: Compute the matter power spectrum
  //=============================================================================

  // ...
  // ...
  // ...

  return pofk;
}

//====================================================
// Get methods
//====================================================
double PowerSpectrum::get_cell_TT(const double ell) const{
  return cell_TT_spline(ell);
}
double PowerSpectrum::get_cell_TE(const double ell) const{
  return cell_TE_spline(ell);
}
double PowerSpectrum::get_cell_EE(const double ell) const{
  return cell_EE_spline(ell);
}

//====================================================
// Output the cells to file
//====================================================

void PowerSpectrum::output(std::string filename) const{
  // Output in standard units of muK^2
  std::ofstream fp(filename.c_str());
  const int ellmax = int(ells[ells.size()-1]);
  auto ellvalues = Utils::linspace(2, ellmax, ellmax-1);
  auto print_data = [&] (const double ell) {
    double normfactor  = (ell * (ell+1)) / (2.0 * M_PI) * pow(1e6 * cosmo->get_TCMB(), 2);
    double normfactorN = (ell * (ell+1)) / (2.0 * M_PI)
      * pow(1e6 * cosmo->get_TCMB() *  pow(4.0/11.0, 1.0/3.0), 2);
    double normfactorL = (ell * (ell+1)) * (ell * (ell+1)) / (2.0 * M_PI);
    fp << ell                                 << " ";
    fp << cell_TT_spline( ell ) * normfactor  << " ";
    if(Constants.polarization){
      fp << cell_EE_spline( ell ) * normfactor  << " ";
      fp << cell_TE_spline( ell ) * normfactor  << " ";
    }
    fp << "\n";
  };
  std::for_each(ellvalues.begin(), ellvalues.end(), print_data);
}
