#include"RecombinationHistory.h"

//====================================================
// Constructors
//====================================================

RecombinationHistory::RecombinationHistory(
    BackgroundCosmology *cosmo,
    double Yp) :
  cosmo(cosmo),
  Yp(Yp)
{}

//====================================================
// Do all the solving we need to do
//====================================================

void RecombinationHistory::solve(){

  // Compute and spline Xe, ne
  solve_number_density_electrons();

  // Compute and spline tau, dtaudx, ddtauddx, g, dgdx, ddgddx, ...
  solve_for_optical_depth_tau();
}

//====================================================
// Solve for X_e and n_e using Saha and Peebles and spline the result
//====================================================

void RecombinationHistory::solve_number_density_electrons(){
  Utils::StartTiming("Xe");

  Vector x_array = Utils::linspace(x_start, x_end, npts_rec_arrays);
  Vector Xe_arr(npts_rec_arrays);
  Vector ne_arr(npts_rec_arrays);

  // Calculate recombination history
  bool saha_regime = true;
  for(int i = 0; i < npts_rec_arrays; i++){
    //Get X_e from solving the Saha equation
    auto Xe_ne_data = electron_fraction_from_saha_equation(x_array[i]);

    // Electron fraction and number density
    const double Xe_current = Xe_ne_data.first;
    const double ne_current = Xe_ne_data.second;

    // Are we still in the Saha regime?
    if(Xe_current < Xe_saha_limit){
      saha_regime = false;
    }
    if(saha_regime){
      Xe_arr[i] = Xe_current;
      ne_arr[i] = ne_current;
    }
    else {
      // The Peebles ODE equation
      // debugging hint: if(rand() % 1000 == 0) Cout << x << “\n”;
      auto x_array_peebles = Vector({x_array[i-1], x_array[i]});
      ODESolver peebles_Xe_ode;
      ODEFunction dXedx = [&](double x, const double *Xe, double *dXedx){
        return rhs_peebles_ode(x, Xe, dXedx);
      };

      //get initial conditions
      double peeblesini = Xe_arr[i-1];
      Vector peebles_ic{peeblesini};

      //solve and retrieve Xe
      peebles_Xe_ode.solve(dXedx, x_array_peebles, peebles_ic);
      Xe_arr[i] = peebles_Xe_ode.get_data_by_component(0)[1];
      if(Xe_arr[i] != Xe_arr[i]) exit(1); //check for nans

      //get ne from Xe
      double OmegaB   = cosmo->get_OmegaB(0);
      double H0       = cosmo->get_H0();
      double rho_crit = 3.*pow(H0, 2.)/(8.*M_PI*Constants.G);
      double a = exp(x_array[i]);

      ne_arr[i] = Xe_arr[i] * OmegaB*rho_crit/(Constants.m_H*pow(a, 3.));
    }
    //std::cout << Xe_arr[i] << std::endl;
  }

  //create splines
  Vector log_Xe = log(Xe_arr);
  Vector log_ne = log(ne_arr);
  log_Xe_of_x_spline.create(x_array, log_Xe, "Xe");
  log_ne_of_x_spline.create(x_array, log_ne, "ne");

  Utils::EndTiming("Xe");
}

//====================================================
// Solve the Saha equation to get ne and Xe
//====================================================
std::pair<double,double> RecombinationHistory::electron_fraction_from_saha_equation(double x) const{
  const double a           = exp(x);

  // Physical constants
  const double k_b         = Constants.k_b;
  const double G           = Constants.G;
  const double m_e         = Constants.m_e;
  const double hbar        = Constants.hbar;
  const double c           = Constants.c;
  const double m_H         = Constants.m_H;
  const double epsilon_0   = Constants.epsilon_0; //units J

  // Fetch cosmological parameters
  double OmegaB   = cosmo->get_OmegaB(0);
  double Tb       = cosmo->get_TCMB(x)*k_b; //units J
  double H0       = cosmo->get_H0(); //units 1/s
  double rho_crit = 3.*pow(H0, 2.)/(8.*M_PI*Constants.G);

  // Electron fraction and number density
  double Xe = 0.0;
  double ne = 0.0;

  //some things to make the Saha equation easier to read
  double nb = OmegaB*rho_crit/(m_H*pow(a, 3.)); //units 1/m^3
  double b = 1./nb*pow(m_e*Tb/(2.*M_PI*hbar*hbar), 3./2.)*exp(-epsilon_0/Tb); //dimensionless

  //find Xe and nb using the Saha equation
  if(4./b < 1e-8){
    //large b approximation
    Xe = 1.;
  }
  else{
    Xe = (-b + b*sqrt(1+4./b))/2.;
  }
  ne = Xe*nb;

  return std::pair<double,double>(Xe, ne);
}

//====================================================
// The right hand side of the dXedx Peebles ODE
//====================================================
int RecombinationHistory::rhs_peebles_ode(double x, const double *Xe, double *dXedx){

  // Current value of a and X_e
  const double X_e         = Xe[0];
  const double a           = exp(x);

  // Physical constants in SI units
  const double k_b          = Constants.k_b;
  const double G            = Constants.G;
  const double c            = Constants.c;
  const double m_e          = Constants.m_e;
  const double hbar         = Constants.hbar;
  const double m_H          = Constants.m_H;
  const double sigma_T      = Constants.sigma_T;
  const double lambda_2s1s  = Constants.lambda_2s1s; //units 1/s
  const double epsilon_0    = Constants.epsilon_0; //units J

  // Cosmological parameters
  double OmegaB      = cosmo->get_OmegaB(0);
  double H           = cosmo->H_of_x(x); //units 1/s
  double H0          = cosmo->get_H0();
  double Tb          = cosmo->get_TCMB(x)*Constants.k_b; //units J

  //Write expression for dXedx
  double n_H = 3.*pow(H0, 2.)*OmegaB/(8*M_PI*G*m_H*pow(a, 3.)); //units 1/m^3
  double n1s = (1-X_e)*n_H; //units 1/m^3

  double lambda_alpha = H*pow(3.*epsilon_0/(hbar*c), 3.)/(pow(8.*M_PI, 2.)*n1s); //units 1/s
  double phi2 = 0.448*log(epsilon_0/Tb); //dimensionless
  double alpha2 = 24.*c/sqrt(27.*M_PI)*sigma_T*sqrt(epsilon_0/Tb)*phi2; //units m^3/s
  double beta2 = alpha2*pow(m_e*Tb/(2*M_PI*hbar*hbar), 3./2.)*exp(-epsilon_0/(4.*Tb)); //units 1/s
  double beta = alpha2*pow(m_e*Tb/(2*M_PI*hbar*hbar), 3./2.)*exp(-epsilon_0/Tb); //units 1/s

  double Cr = (lambda_2s1s + lambda_alpha)/(lambda_2s1s + lambda_alpha + beta2); //dimensionless

  dXedx[0] = Cr/H*(beta*(1-X_e) - n_H*alpha2*pow(X_e, 2.));
  //std::cout << x << " " << X_e << " " << dXedx[0] << "\n";

  return GSL_SUCCESS;
}

//====================================================
// Solve for the optical depth tau, compute the
// visibility function and spline the result
//====================================================

void RecombinationHistory::solve_for_optical_depth_tau(){
  Utils::StartTiming("opticaldepth");

  // Set up x-arrays to integrate over.
  const int npts = 1000;
  Vector x_array = Utils::linspace(x_end, x_start, npts); //integrate backwards

  // The ODE
  ODEFunction dtaudx = [&](double x, const double *tau, double *dtaudx){

    const double c       = Constants.c;
    const double sigma_T = Constants.sigma_T;

    double ne = ne_of_x(x);
    double H  = cosmo->H_of_x(x); //units 1/s

    // Set the derivative for photon optical depth
    dtaudx[0] = -c*ne*sigma_T/H;

    return GSL_SUCCESS;
  };

  //set up and solve tau(x)
  double tauini = 0.0;
  Vector tau_ic{tauini};

  ODESolver tau_ode;
  tau_ode.solve(dtaudx, x_array, tau_ic);

  auto tau_array = tau_ode.get_data_by_component(0);

  //reverse order of x and tau so that x is strictly increasing
  std::reverse(x_array.begin(), x_array.end());
  std::reverse(tau_array.begin(), tau_array.end());
  //create splines
  tau_of_x_spline.create(x_array, tau_array, "tau");

  //compute visibility function
  Vector g_array(npts);
  for(int i = 0; i < npts; i++){
    double x    = x_array[i];
    double dtau = tau_of_x_spline.deriv_x(x);
    double tau  = tau_of_x_spline(x);

    g_array[i] = -dtau*exp(-tau);
  }

  //spline g_tilde
  g_tilde_of_x_spline.create(x_array, g_array, "g");

  Utils::EndTiming("opticaldepth");
}

//====================================================
// Get methods
//====================================================

double RecombinationHistory::tau_of_x(double x) const{
  return tau_of_x_spline(x);
}

double RecombinationHistory::dtaudx_of_x(double x) const{
  return tau_of_x_spline.deriv_x(x);
}

double RecombinationHistory::ddtauddx_of_x(double x) const{
  return tau_of_x_spline.deriv_xx(x);
}

double RecombinationHistory::g_tilde_of_x(double x) const{
  return g_tilde_of_x_spline(x);
}

double RecombinationHistory::dgdx_tilde_of_x(double x) const{
  return g_tilde_of_x_spline.deriv_x(x);
}

double RecombinationHistory::ddgddx_tilde_of_x(double x) const{
  return g_tilde_of_x_spline.deriv_xx(x);
}

double RecombinationHistory::Xe_of_x(double x) const{
  double logXe = log_Xe_of_x_spline(x);
  return exp(logXe);
}

double RecombinationHistory::ne_of_x(double x) const{
  double logne = log_ne_of_x_spline(x);
  return exp(logne);
}

double RecombinationHistory::get_Yp() const{
  return Yp;
}

//====================================================
// Print some useful info about the class
//====================================================
void RecombinationHistory::info() const{
  //x_decouple = 


  std::cout << "\n";
  std::cout << "Info about recombination/reionization history class:\n";
  std::cout << "Yp:            " << Yp                << "\n";
  std::cout << "x(decoupling): " << -6.986841         << "\n"; //you still need to automate this instead of trial and error
  std::cout << "z(decoupling): " << exp(6.986841) - 1 << "\n";
  std::cout << "x(rec):        " << -7.163781         << "\n"; //also should be automated
  std::cout << "z(rec):        " << exp(7.163781) - 1 << "\n";
  std::cout << "Freezeout Xe:  " << Xe_of_x(0.)       << "\n";
  std::cout << std::endl;
}

//====================================================
// Output the data computed to file
//====================================================
void RecombinationHistory::output(const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts       = 5000;
  const double x_min   = x_start;
  const double x_max   = x_end;

  Vector x_array = Utils::linspace(x_min, x_max, npts);
  auto print_data = [&] (const double x) {
    fp << x                    << " ";
    fp << Xe_of_x(x)           << " ";
    fp << ne_of_x(x)           << " ";
    fp << tau_of_x(x)          << " ";
    fp << dtaudx_of_x(x)       << " ";
    fp << ddtauddx_of_x(x)     << " ";
    fp << g_tilde_of_x(x)      << " ";
    fp << dgdx_tilde_of_x(x)   << " ";
    fp << ddgddx_tilde_of_x(x) << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}
