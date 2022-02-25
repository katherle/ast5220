#include "BackgroundCosmology.h"

//====================================================
// Constructors
//====================================================

BackgroundCosmology::BackgroundCosmology(
    double h,
    double OmegaB,
    double OmegaCDM,
    double OmegaK,
    double Neff,
    double TCMB) :
  h(h),
  OmegaB(OmegaB),
  OmegaCDM(OmegaCDM),
  OmegaK(OmegaK),
  Neff(Neff),
  TCMB(TCMB)
{
  //define derived variables
  OmegaNu = 0.0;

  H0 = h * Constants.H0_over_h;
  OmegaR = 2.*pow(M_PI, 2)/30.*pow(Constants.k_b * TCMB, 4)/(pow(Constants.hbar, 3) * pow(Constants.c, 5))*8*M_PI*Constants.G/(3*pow(H0, 2));

  OmegaLambda = 1. - OmegaK - OmegaB - OmegaCDM - OmegaR - OmegaNu;
}
//====================================================
// Do all the solving. Compute eta(x)
//====================================================

// Solve the background
void BackgroundCosmology::solve(){
  Utils::StartTiming("Eta");

  //set up x points
  int npts = 100;
  Vector x_array;
  x_array = Utils::linspace(Constants.x_start, Constants.x_end, npts);

  // The ODE for deta/dx
  ODEFunction detadx = [&](double x, const double *eta, double *detadx){
    detadx[0] = Constants.c/Hp_of_x(x);
    return GSL_SUCCESS;
  };

  //set initial conditions
  //use analytical approximation
  double etaini = Constants.c/Hp_of_x(Constants.x_start);
  Vector eta_ic{etaini};

  //solve ODE
  ODESolver ode;
  ode.solve(detadx, x_array, eta_ic);

  //retrieve solution
  auto eta_array = ode.get_data_by_component(0);

  //create spline of solution
  eta_of_x_spline.create(x_array, eta_array, "Eta(x)");

  Utils::EndTiming("Eta");

  //repeat for cosmic time
  Utils::StartTiming("t");

  // The ODE for dt/dx
  ODEFunction dtdx = [&](double x, const double *t, double *dtdx){
    dtdx[0] = 1/H_of_x(x);
    return GSL_SUCCESS;
  };

  //set initial conditions
  //use analytical approximation
  double tini = 1/(2.*H_of_x(Constants.x_start));
  Vector t_ic{tini};

  //solve ODE
  ODESolver ode_t;
  ode_t.solve(dtdx, x_array, t_ic);

  //retrieve solution
  auto t_array = ode_t.get_data_by_component(0);

  //create spline of solution
  t_of_x_spline.create(x_array, t_array, "t(x)");

  Utils::EndTiming("t");
}

//====================================================
// Get methods
//====================================================

double BackgroundCosmology::eta_of_x(double x) const{
  return eta_of_x_spline(x);
}

double BackgroundCosmology::get_detadx(double x) const{
  return eta_of_x_spline.deriv_x(x);
}

double BackgroundCosmology::t_of_x(double x) const{
  return t_of_x_spline(x);
}

double BackgroundCosmology::get_dtdx(double x) const{
  return t_of_x_spline.deriv_x(x);
}

double BackgroundCosmology::comoving_distance(double x) const{
  return(eta_of_x(0.0) - eta_of_x(x));
}

double BackgroundCosmology::proper_distance(double x) const{
  double X = comoving_distance(x);
  double paren = sqrt(abs(OmegaK))*H0*X/Constants.c;

  if(OmegaK < 0.0){
    return X*sin(paren)/paren;
  }
  else if(OmegaK == 0.0){
    return X;
  }
  else {
    return X*sinh(paren)/paren;
  }
}

double BackgroundCosmology::angular_distance(double x) const{
  double r = proper_distance(x);
  return pow(M_E, x)*r;
}

double BackgroundCosmology::luminosity_distance(double x) const{
  double da = angular_distance(x);
  return da/pow(M_E, 2.*x);
}

double BackgroundCosmology::H_of_x(double x) const{
  //implement Friedmann eq
  double OmegaM = OmegaB + OmegaCDM;
  double OmegaRad = OmegaR + OmegaNu;
  double H = H0 * sqrt(OmegaM/pow(M_E, 3.*x) + OmegaRad/pow(M_E, 4.*x) + OmegaK/pow(M_E, 2.*x) + OmegaLambda);

  return H;
}

double BackgroundCosmology::Hp_of_x(double x) const{
  //Hp = aH
  return pow(M_E, x)*H_of_x(x);
}

double BackgroundCosmology::dHdx_of_x(double x) const{
  double OmegaM = OmegaB + OmegaCDM;
  double OmegaRad = OmegaR + OmegaNu;

  return H0*H0*(-3.*OmegaM/pow(M_E, 3.*x) - 4.*OmegaRad/pow(M_E, 4.*x) - 2.*OmegaK/pow(M_E, 2.*x))/(2.*H_of_x(x));
}

double BackgroundCosmology::dHpdx_of_x(double x) const{
  return pow(M_E, x)*(H_of_x(x) + dHdx_of_x(x));
}

double BackgroundCosmology::ddHddx_of_x(double x) const{
  double OmegaM = OmegaB + OmegaCDM;
  double OmegaRad = OmegaR + OmegaNu;

  double t1 = (9.+3.*dHdx_of_x(x)/H_of_x(x))*OmegaM/pow(M_E, 3.*x);
  double t2 = (16.+4.*dHdx_of_x(x)/H_of_x(x))*OmegaRad/pow(M_E, 4.*x);
  double t3 = (4.+2.*dHdx_of_x(x)/H_of_x(x))*OmegaK/pow(M_E, 2.*x);

  return H0*H0/(2.*H_of_x(x))*(t1 + t2 + t3);
}

double BackgroundCosmology::ddHpddx_of_x(double x) const{
  return pow(M_E, x)*(H_of_x(x) + 2.*dHdx_of_x(x) + ddHddx_of_x(x));
}

double BackgroundCosmology::get_OmegaB(double x) const{
  if(x == 0.0){
    return OmegaB;
  }
  else {
    return pow(H0, 2.)*OmegaB/(pow(M_E, 3.*x)*pow(H_of_x(x), 2.));
  }
}

double BackgroundCosmology::get_OmegaM(double x) const{

    return get_OmegaB(x)+get_OmegaCDM(x);
}

double BackgroundCosmology::get_OmegaR(double x) const{
  if(x == 0.0) {
    return OmegaR;
  }
  else {
    return pow(H0, 2.)*OmegaR/(pow(M_E, 4.*x)*pow(H_of_x(x), 2.));
  }
}

double BackgroundCosmology::get_OmegaRtot(double x) const{
  return get_OmegaR(x)+get_OmegaNu(x);
}

double BackgroundCosmology::get_OmegaNu(double x) const{
  if(x == 0.0) {
    return OmegaNu;
  }
  else {
    return pow(H0, 2.)*OmegaNu/(pow(M_E, 4.*x)*pow(H_of_x(x), 2.));
  }
}

double BackgroundCosmology::get_OmegaCDM(double x) const{
  if(x == 0.0) {
    return OmegaCDM;
  }
  else {
    return pow(H0, 2.)*OmegaCDM/(pow(M_E, 3.*x)*pow(H_of_x(x), 2.));
  }
}

double BackgroundCosmology::get_OmegaLambda(double x) const{
  if(x == 0.0) {
    return OmegaLambda;
  }
  else {
    return pow(H0, 2.)*OmegaLambda/pow(H_of_x(x), 2.);
  }
}

double BackgroundCosmology::get_OmegaK(double x) const{
  if(x == 0.0) {
    return OmegaK;
  }
  else {
    return pow(H0, 2.)*OmegaK/(pow(M_E, 2.)*pow(H_of_x(x), 2.));
  }
}

double BackgroundCosmology::get_H0() const{
  return H0;
}

double BackgroundCosmology::get_h() const{
  return h;
}

double BackgroundCosmology::get_Neff() const{
  return Neff;
}

double BackgroundCosmology::get_TCMB(double x) const{
  if(x == 0.0) return TCMB;
  return TCMB * exp(-x);
}

//====================================================
// Print out info about the class
//====================================================
void BackgroundCosmology::info() const{
  std::cout << "\n";
  std::cout << "Info about cosmology class:\n";
  std::cout << "OmegaB:      " << OmegaB      << "\n";
  std::cout << "OmegaCDM:    " << OmegaCDM    << "\n";
  std::cout << "OmegaLambda: " << OmegaLambda << "\n";
  std::cout << "OmegaK:      " << OmegaK      << "\n";
  std::cout << "OmegaNu:     " << OmegaNu     << "\n";
  std::cout << "OmegaR:      " << OmegaR      << "\n";
  std::cout << "Neff:        " << Neff        << "\n";
  std::cout << "h:           " << h           << "\n";
  std::cout << "TCMB:        " << TCMB        << "\n";
  std::cout << std::endl;
}

//====================================================
// Output some data to file
//====================================================
void BackgroundCosmology::output(const std::string filename) const{
  const double x_min = -10.0;
  const double x_max =  0.0;
  const int    n_pts =  100;

  Vector x_array = Utils::linspace(x_min, x_max, n_pts);

  std::ofstream fp(filename.c_str());
  auto print_data = [&] (const double x) {
    fp << x                      << " ";
    fp << eta_of_x(x)            << " ";
    fp << t_of_x(x)              << " ";
    fp << H_of_x(x)              << " ";
    fp << Hp_of_x(x)             << " ";
    fp << dHpdx_of_x(x)          << " ";
    fp << ddHpddx_of_x(x)        << " ";
    fp << get_OmegaB(x)          << " ";
    fp << get_OmegaCDM(x)        << " ";
    fp << get_OmegaLambda(x)     << " ";
    fp << get_OmegaR(x)          << " ";
    fp << get_OmegaNu(x)         << " ";
    fp << luminosity_distance(x) << " ";
    fp <<"\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}
