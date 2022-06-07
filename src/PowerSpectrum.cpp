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

  generate_bessel_function_splines();

  // Implement line_of_sight_integration
  // make sure you have enough points for the bessel functions
  Vector k_array = Utils::linspace(k_min, k_max, n_k);
  line_of_sight_integration(k_array);

  // Integration to get Cell by solving dCell^f/dlogk = Delta(k) * f_ell(k)^2
  double eta_0 = cosmo->eta_of_x(0.);
  int nk = (k_max - k_min)*16*eta_0/M_PI;
  Vector log_k_array = Utils::linspace(log(k_min), log(k_max), nk);

  auto cell_TT = solve_for_cell(log_k_array, thetaT_ell_of_k_spline, thetaT_ell_of_k_spline);
  cell_TT_spline.create(ells, cell_TT, "Cell_TT_of_ell");
}

//====================================================
// Generate splines of j_ell(z) needed for LOS integration
//====================================================

void PowerSpectrum::generate_bessel_function_splines(){
  Utils::StartTiming("besselspline");

  // Make storage for the splines
  j_ell_splines = std::vector<Spline>(ells.size());

  double eta_0 = cosmo->eta_of_x(0.0);
  double x_max = k_max*eta_0;
  int n_x = x_max*16/(2*M_PI);
  Vector x_bessel = Utils::linspace(0.0, x_max, n_x);

  //find the Bessel function
  Vector j_ell(n_x);

  for(size_t il = 0; il < ells.size(); il++){
    const int ell = ells[il];

    for(int ix = 0; ix < n_x; ix++){
      const double x = x_bessel[ix];
      j_ell[ix] = Utils::j_ell(ell, x);
    }

    // Make the j_ell_splines[il] spline
    j_ell_splines[il].create(x_bessel, j_ell, "j_ell_spline");
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
  //double delta_x = 2.*M_PI/6.;
  double delta_x = 0.05;
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

      Theta = source_function(Constants.x_start, k)*j_ell_splines[il](k*(eta_0 - eta_start))
              + source_function(0.0, k)*j_ell_splines[il](k*(eta_0 - eta_0));

      for(int ix = 1; ix < nx-1; ix++){
        const double x = x_array[ix];
        double eta_x   = cosmo->eta_of_x(x);

        Theta += 2.*source_function(x, k)*j_ell_splines[il](k*(eta_0 - eta_x));
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
  const int nells      = ells.size();

  // Make storage for the splines we are to create
  thetaT_ell_of_k_spline = std::vector<Spline>(nells);

  //============================================================================
  // TODO: Solve for Theta_ell(k) and spline the result
  //============================================================================

  // Make a function returning the source function
  // to test source function: Set S(k, x) = g
  // CMB power spectrum should be C(l)(l+1) = constant
  // or set S(k, x) = SW or ISW, etc
  std::function<double(double,double)> source_function_T = [&](double x, double k){
    return pert->get_Source_T(x,k);
    //return rec->g_tilde_of_x(x);
  };

  // Do the line of sight integration
  Vector2D thetaT_ell_of_k = line_of_sight_integration_single(k_array, source_function_T);

  // Spline the result and store it in thetaT_ell_of_k_spline
  for (size_t il = 0; il < nells; il++){
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
  const int nk         = log_k_array.size();

  //============================================================================
  // TODO: Integrate Cell = Int 4 * pi * P(k) f_ell g_ell dk/k
  // or equivalently solve the ODE system dCell/dlogk = 4 * pi * P(k) * f_ell * g_ell
  //============================================================================

  Vector result(nells);

  double C_ell;
  const double delta_log_k = abs(log(k_max) - log(k_min))/nk;

  for(int il = 0; il < nells; il++){
    C_ell = 4.*M_PI*(primordial_power_spectrum(k_min)*f_ell_spline[il](k_min)*g_ell_spline[il](k_min)
                      + primordial_power_spectrum(k_max)*f_ell_spline[il](k_max)*g_ell_spline[il](k_max));

    for(size_t ik = 1; ik < nk-1; ik++){
      const double k = exp(log_k_array[ik]);
      C_ell += 8.*M_PI*primordial_power_spectrum(k)*f_ell_spline[il](k)*g_ell_spline[il](k);
    }

    C_ell = C_ell*delta_log_k/2.;

    // Store the result
    result[il] = C_ell;
  }

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

double PowerSpectrum::get_matter_power_spectrum(const double x, const double k) const{
  const double c      = Constants.c;
  const double Phi    = pert->get_Phi(x, k);
  const double OmegaM = cosmo->get_OmegaM(0.);
  const double a      = exp(x);
  const double H0     = cosmo->get_H0();

  double DeltaM = pow(c*k/H0, 2.)*2.*a*Phi/(3.*OmegaM);
  double pofk   = pow(abs(DeltaM), 2.)*primordial_power_spectrum(k)*2.*pow(M_PI, 2.)/pow(k, 3.);

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

void PowerSpectrum::output(std::string filename1, std::string filename2, std::string filename3) const{
  // Output in standard units of muK^2
  std::ofstream fp(filename1.c_str());
  const int ellmax = int(ells[ells.size()-1]);
  auto ellvalues = Utils::linspace(2, ellmax, ellmax-1);

  auto print_data = [&] (const double ell) {
    double normfactor  = (ell * (ell+1.)) / (2.0 * M_PI) * pow(1e6 * cosmo->get_TCMB(0.), 2.);
    double normfactorN = (ell * (ell+1.)) / (2.0 * M_PI)
      * pow(1e6 * cosmo->get_TCMB() *  pow(4.0/11.0, 1.0/3.0), 2);
    double normfactorL = (ell * (ell+1)) * (ell * (ell+1)) / (2.0 * M_PI);

    fp << ell << " ";
    fp << cell_TT_spline( ell ) * normfactor;

    if(Constants.polarization){
      fp << " ";
      fp << cell_EE_spline( ell ) * normfactor << " ";
      fp << cell_TE_spline( ell ) * normfactor;
    }

    fp << "\n";
  };

  std::for_each(ellvalues.begin(), ellvalues.end(), print_data);

  std::ofstream fp2(filename2.c_str());
  Vector k_array = Utils::linspace(k_min, k_max, n_k);

  auto print_data2 = [&] (const double k) {
    fp2 << k << " ";
    fp2 << get_matter_power_spectrum(0.0, k) << " ";
    fp2 << thetaT_ell_of_k_spline[4](k)      << " ";
    fp2 << thetaT_ell_of_k_spline[12](k)     << " ";
    fp2 << thetaT_ell_of_k_spline[19](k)     << " ";
    fp2 << thetaT_ell_of_k_spline[24](k)     << " ";
    fp2 << thetaT_ell_of_k_spline[32](k)     << " ";
    fp2 << thetaT_ell_of_k_spline[42](k)     << "\n";
  };
  std::for_each(k_array.begin(), k_array.end(), print_data2);

  std::ofstream fp3(filename3.c_str());
  Vector x_array = Utils::linspace(0, 150, 1000);

  auto print_data3 = [&] (const double x) {
    fp3 << x                     << " ";
    fp3 << j_ell_splines[0](x)  << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data3);
}
