#include"Perturbations.h"

//====================================================
// Constructors
//====================================================

Perturbations::Perturbations(
    BackgroundCosmology *cosmo,
    RecombinationHistory *rec) :
  cosmo(cosmo),
  rec(rec)
{}

//====================================================
// Do all the solving
//====================================================

void Perturbations::solve(){

  // Integrate all the perturbation equation and spline the result
  integrate_perturbations();

  // Compute source functions and spline the result
  compute_source_functions();
}

//====================================================
// The main work: integrate all the perturbations
// and spline the results
//====================================================

void Perturbations::integrate_perturbations(){
  Utils::StartTiming("integrateperturbation");

  Vector k_array = Utils::linspace(log(k_min), log(k_max), n_k);
  //log spacing
  for (int k = 0; k < n_k; k++){
    k_array[k] = exp(k_array[k]);
  }
  Vector x_array = Utils::linspace(x_start, 0.0, n_x);

  //declare vectors to store data in later
  Vector delta_cdm(n_x*n_k, 0.);
  Vector v_cdm(n_x*n_k, 0.);
  Vector delta_b(n_x*n_k, 0.);
  Vector v_b(n_x*n_k, 0.);
  Vector Theta0(n_x*n_k, 0.);
  Vector Theta1(n_x*n_k, 0.);
  Vector Theta2(n_x*n_k, 0.);
  Vector Theta3(n_x*n_k, 0.);
  Vector Theta4(n_x*n_k, 0.);
  Vector Phi(n_x*n_k, 0.);
  Vector Psi(n_x*n_k, 0.);
  Vector Pi(n_x*n_k, 0.);

  // Loop over all wavenumbers
  for(int ik = 0; ik < n_k; ik++){

    // Progress bar...
    if( (10*ik) / n_k != (10*ik+10) / n_k ) {
      std::cout << (100*ik+100)/n_k << "% " << std::flush;
      if(ik == n_k-1) std::cout << std::endl;
    }

    double k = k_array[ik];

    // Find value to integrate to
    double x_end_tight = get_tight_coupling_time(k);

    //create x array that only goes to x_end_tight:
    Vector x_dummy(n_x);
    for(int i = 0; i<n_x; i++){
      x_dummy[i] = abs(x_array[i] - x_end_tight);
    }
    auto it = std::min_element(x_dummy.begin(), x_dummy.end());
    int index_tc = std::distance(x_dummy.begin(), it);
    Vector x_tc = {x_array.begin(), x_array.begin()+index_tc};

    // Set up initial conditions in the tight coupling regime
    Vector y_tight_coupling_ini = set_ic(x_start, k);

    // The tight coupling ODE system
    ODEFunction dydx_tight_coupling = [&](double x, const double *y, double *dydx){
      return rhs_tight_coupling_ode(x, k, y, dydx);
    };

    // Integrate from x_start -> x_end_tight
    ODESolver y_tight_coupling;
    y_tight_coupling.solve(dydx_tight_coupling, x_tc, y_tight_coupling_ini);

    auto y_tc = y_tight_coupling.get_data();
    auto y_tc_final = y_tight_coupling.get_final_data();

    //Full equation integration

    // Set up initial conditions
    auto y_full_ini = set_ic_after_tight_coupling(y_tc_final, x_end_tight, k);

    // The full ODE system
    ODEFunction dydx_full = [&](double x, const double *y, double *dydx){
      return rhs_full_ode(x, k, y, dydx);
    };

    // Integrate from x_end_tight -> x_end
    Vector x_full = {x_array.begin() + index_tc, x_array.end()};
    ODESolver y_full_ode;
    y_full_ode.solve(dydx_full, x_full, y_full_ini);

    auto y_full = y_full_ode.get_data();

    //store the data in vectors:
    const double c = Constants.c;
    const double H0 = cosmo->get_H0();
    double x;
    double Hp;
    double OmegaR;
    double dtau;

    for (int i = 0; i<index_tc; i++){
      x = x_array[i];
      Hp = cosmo->Hp_of_x(x);
      OmegaR = cosmo->get_OmegaR(x);
      dtau = rec->dtaudx_of_x(x);

      delta_cdm[i + n_x*ik] = y_tc[i][Constants.ind_deltacdm_tc];
      v_cdm[i + n_x*ik]     = y_tc[i][Constants.ind_vcdm_tc];
      delta_b[i + n_x*ik]   = y_tc[i][Constants.ind_deltab_tc];
      v_b[i + n_x*ik]       = y_tc[i][Constants.ind_vb_tc];
      Theta0[i + n_x*ik]    = y_tc[i][Constants.ind_start_theta_tc];
      Theta1[i + n_x*ik]    = y_tc[i][Constants.ind_start_theta_tc + 1];
      Theta2[i + n_x*ik]    = -20.*c*k/(45.*Hp*dtau)*Theta1[i + n_x*ik];
      Theta3[i + n_x*ik]    = -3.*c*k/(7.*Hp*dtau)*Theta2[i + n_x*ik];
      Theta4[i + n_x*ik]    = -4.*c*k/(9.*Hp*dtau)*Theta3[i + n_x*ik];
      Phi[i + n_x*ik]       = y_tc[i][Constants.ind_Phi_tc];
      Psi[i + n_x*ik]       = -Phi[i + n_x*ik] - 12.*pow(H0/(c*k), 2.)*(OmegaR*Theta2[i + n_x*ik] + 0.);
      Pi[i + n_x*ik]        = Theta2[i + n_x*ik]; //no polarization
    }
    for (int i = index_tc; i < n_x; i++){
      x = x_array[i];
      OmegaR = cosmo->get_OmegaR(x);

      delta_cdm[i + n_x*ik] = y_full[i-index_tc][Constants.ind_deltacdm];
      v_cdm[i + n_x*ik]     = y_full[i-index_tc][Constants.ind_vcdm];
      delta_b[i + n_x*ik]   = y_full[i-index_tc][Constants.ind_deltab];
      v_b[i + n_x*ik]       = y_full[i-index_tc][Constants.ind_vb];
      Theta0[i + n_x*ik]    = y_full[i-index_tc][Constants.ind_start_theta];
      Theta1[i + n_x*ik]    = y_full[i-index_tc][Constants.ind_start_theta + 1];
      Theta2[i + n_x*ik]    = y_full[i-index_tc][Constants.ind_start_theta + 2];
      Theta3[i + n_x*ik]    = y_full[i-index_tc][Constants.ind_start_theta + 3];
      Theta4[i + n_x*ik]    = y_full[i-index_tc][Constants.ind_start_theta + 4];
      Phi[i + n_x*ik]       = y_full[i-index_tc][Constants.ind_Phi];
      Psi[i + n_x*ik]       = -Phi[i + n_x*ik] - 12.*pow(H0/(c*k), 2.)*(OmegaR*Theta2[i + n_x*ik] + 0.);
      Pi[i + n_x*ik]        = Theta2[i + n_x*ik]; //no polarization
    }
  }
  Utils::EndTiming("integrateperturbation");

  delta_cdm_spline.create(x_array, k_array, delta_cdm, "delta_cdm_spline");
  v_cdm_spline.create(x_array, k_array, v_cdm, "v_cdm_spline");
  delta_b_spline.create(x_array, k_array, delta_b, "delta_b_spline");
  v_b_spline.create(x_array, k_array, v_b, "v_b_spline");

  Theta_spline = std::vector<Spline2D>(5);
  Theta_spline[0].create(x_array, k_array, Theta0, "Theta0_spline");
  Theta_spline[1].create(x_array, k_array, Theta1, "Theta1_spline");
  Theta_spline[2].create(x_array, k_array, Theta2, "Theta2_spline");
  Theta_spline[3].create(x_array, k_array, Theta3, "Theta3_spline");
  Theta_spline[4].create(x_array, k_array, Theta4, "Theta4_spline");

  Phi_spline.create(x_array, k_array, Phi, "Phi_spline");
  Psi_spline.create(x_array, k_array, Psi, "Psi_spline");
  Pi_spline.create(x_array, k_array, Pi, "Pi_spline");
}

//====================================================
// Set IC at the start of the run (this is in the
// tight coupling regime)
//====================================================
Vector Perturbations::set_ic(const double x, const double k) const{

  // For integration of perturbations in tight coupling regime
  // (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  const int n_ell_tot_tc        = Constants.n_ell_tot_tc;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // useful quantities
  double f_v    = 0.0;
  if(neutrinos){
    double OmegaNu    = cosmo->get_OmegaNu();
    double OmegaRtot  = cosmo->get_OmegaRtot();
    f_v               = OmegaNu/(OmegaRtot);
  }
  const double c      = Constants.c;
  const double Hp     = cosmo->Hp_of_x(x);
  const double dtau   = rec->dtaudx_of_x(x);

  // The vector we are going to fill
  Vector y_tc(n_ell_tot_tc);

  // References to the tight coupling quantities
  double &delta_cdm    =  y_tc[Constants.ind_deltacdm_tc];
  double &delta_b      =  y_tc[Constants.ind_deltab_tc];
  double &v_cdm        =  y_tc[Constants.ind_vcdm_tc];
  double &v_b          =  y_tc[Constants.ind_vb_tc];
  double &Phi          =  y_tc[Constants.ind_Phi_tc];
  double *Theta        = &y_tc[Constants.ind_start_theta_tc];
  double *Nu           = &y_tc[Constants.ind_start_nu_tc];

  // SET: Scalar quantities (Gravitational potential, baryons and CDM)
  double Psi = -1./(3./2. + 2.*f_v/5.);
  Phi       = -(1. + 2.*f_v/5.)*Psi;
  delta_cdm = -3.*Psi/2.;
  delta_b   = delta_cdm;
  v_cdm     = -c*k*Psi/(2.*Hp);
  v_b       = v_cdm;

  // SET: Photon temperature perturbations (Theta_ell)
  //theta0
  Theta[0] = -Psi/2.;

  //theta1
  Theta[1] = c*k*Psi/(6.*Hp);


  // SET: Neutrino perturbations (N_ell)
  if(neutrinos){
    // Not included
  }

  return y_tc;
}

//====================================================
// Set IC for the full ODE system after tight coupling
// regime ends
//====================================================

Vector Perturbations::set_ic_after_tight_coupling(
    const Vector &y_tc,
    const double x,
    const double k) const{

  // Make the vector we are going to fill
  Vector y(Constants.n_ell_tot_full);

  //=============================================================================
  // Compute where in the y array each component belongs and where corresponding
  // components are located in the y_tc array
  //=============================================================================

  // Number of multipoles we have in the full regime
  const int n_ell_theta         = Constants.n_ell_theta;
  const int n_ell_thetap        = Constants.n_ell_thetap;
  const int n_ell_neutrinos     = Constants.n_ell_neutrinos;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // Number of multipoles we have in the tight coupling regime
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;

  // References to the tight coupling quantities
  const double &delta_cdm_tc    =  y_tc[Constants.ind_deltacdm_tc];
  const double &delta_b_tc      =  y_tc[Constants.ind_deltab_tc];
  const double &v_cdm_tc        =  y_tc[Constants.ind_vcdm_tc];
  const double &v_b_tc          =  y_tc[Constants.ind_vb_tc];
  const double &Phi_tc          =  y_tc[Constants.ind_Phi_tc];
  const double *Theta_tc        = &y_tc[Constants.ind_start_theta_tc];
  const double *Nu_tc           = &y_tc[Constants.ind_start_nu_tc];

  // References to the quantities we are going to set
  double &delta_cdm       =  y[Constants.ind_deltacdm];
  double &delta_b         =  y[Constants.ind_deltab];
  double &v_cdm           =  y[Constants.ind_vcdm];
  double &v_b             =  y[Constants.ind_vb];
  double &Phi             =  y[Constants.ind_Phi];
  double *Theta           = &y[Constants.ind_start_theta];
  double *Theta_p         = &y[Constants.ind_start_thetap];
  double *Nu              = &y[Constants.ind_start_nu];

  // Useful quantities
  double f_v    = 0.0;
  if(neutrinos){
    double OmegaNu    = cosmo->get_OmegaNu();
    double OmegaRtot  = cosmo->get_OmegaRtot();
    f_v               = OmegaNu/(OmegaRtot);
  }
  const double c = Constants.c;
  const double Hp     = cosmo->Hp_of_x(x);
  const double dtau   = rec->dtaudx_of_x(x);

  // SET: Scalar quantities (Gravitational potental, baryons and CDM)
  Phi = Phi_tc;
  delta_cdm = delta_cdm_tc;
  delta_b = delta_b_tc;
  v_cdm = v_cdm_tc;
  v_b = v_b_tc;

  // SET: Photon temperature perturbations (Theta_ell)
  Theta[0] = Theta_tc[0]; //theta0

  Theta[1] = Theta_tc[1]; //theta1

  Theta[2] = -20.*c*k/(45.*Hp*dtau); //theta2
  if(polarization){
    Theta[2] = -8.*c*k/(15.*Hp*dtau);
  }

  //theta_ell
  for(int l = 3; l < n_ell_theta; l++){
    Theta[l]= -l/(2*l + 1) * c*k/(Hp*dtau)*Theta[l-1];
  }

  // SET: Photon polarization perturbations (Theta_p_ell)
  if(polarization){
    //not included
  }

  // SET: Neutrino perturbations (N_ell)
  if(neutrinos){
    //not included
  }

  return y;
}

//====================================================
// The time when tight coupling ends
//====================================================

double Perturbations::get_tight_coupling_time(const double k) const{
  const double c     = Constants.c;
  const double x_rec = rec->get_x_rec();

  // Remember all the three conditions in Callin:
  // |tau'| > 10
  // |k/(Hp tau')| < 1/10
  // x < x_rec

  Vector x = {x_rec, x_start};
  double x1;
  int i = 0, max_i  = 100;
  do{
    x1 = x[1];
    double dtaudx = rec->dtaudx_of_x(x1);
    double ckH = c*k/(cosmo->Hp_of_x(x1));

    if (abs(dtaudx) < 10.*std::max(1., ckH)){
      x[1] -= abs(x[0] - x[1])/2.;
      x[0] = x1;
    }
    if (abs(dtaudx) > 10.*std::max(1., ckH)){
      x[1] += abs(x[0] - x[1])/2.;
      x[0] = x1;
    }
    i++;
  }
  while(x[1] != x[0] && i < max_i);

  // if(i == max_i){
  //   std::cout << "get_tight_coupling_time Warning: max number of iterations reached (k = " << k*Constants.Mpc << "/Mpc)" << std::endl;
  // }

  return(x[1]);
}

//====================================================
// After integrsating the perturbation compute the
// source function(s)
//====================================================
void Perturbations::compute_source_functions(){
  Utils::StartTiming("source");

  Vector k_array = Utils::linspace(k_min, k_max, n_x);
  Vector x_array = Utils::linspace(x_start, 0.0, n_x);

  // Make storage for the source functions (in 1D array to be able to pass it to the spline)
  Vector ST_array(k_array.size() * x_array.size());
  //in principle there should be another source function for polarization
  //however I have neglected polarization and set pi=theta2

  // Compute source functions
  for(auto ix = 0; ix < x_array.size(); ix++){
    const double x = x_array[ix];
    const double a = exp(x);
    #pragma omp parallel for schedule(dynamic, 1)
    for(auto ik = 0; ik < k_array.size(); ik++){
      const double k = k_array[ik];

      // NB: This is the format the data needs to be stored
      // in a 1D array for the 2D spline routine source(ix,ik) -> S_array[ix + nx * ik]
      const int index = ix + n_x * ik;

      // Fetch all the things we need
      const double c         = Constants.c;
      const double Hp        = cosmo->Hp_of_x(x);
      const double dHp       = cosmo->dHpdx_of_x(x);
      const double ddHp      = cosmo->ddHpddx_of_x(x);
      const double H0        = cosmo->get_H0();
      const double OmegaCDM  = cosmo->get_OmegaCDM();
      const double OmegaB    = cosmo->get_OmegaB();
      const double OmegaR    = cosmo->get_OmegaR();
      const double R         = 4.*OmegaR/(3.*OmegaB*a);

      const double g_tilde   = rec->g_tilde_of_x(x);
      const double dg_tilde  = rec->dgdx_tilde_of_x(x);
      const double ddg_tilde = rec->ddgddx_tilde_of_x(x);
      const double tau       = rec->tau_of_x(x);
      const double dtau      = rec->dtaudx_of_x(x);
      const double ddtau     = rec->ddtauddx_of_x(x);

      const double theta0    = get_Theta(x, k, 0);
      const double theta1    = get_Theta(x, k, 1);
      const double pi        = get_Pi(x, k); // = theta2
      const double theta3    = get_Theta(x, k, 3);
      const double theta4    = get_Theta(x, k, 4);
      const double psi       = get_Psi(x, k);
      const double dpsi      = Psi_spline.deriv_x(x, k);
      const double phi       = get_Phi(x, k);
      const double dphi      = Phi_spline.deriv_x(x, k);
      const double delta_cdm = get_delta_cdm(x, k);
      const double delta_b   = get_delta_b(x, k);
      const double vb        = get_v_b(x, k);
      const double dvb       = v_b_spline.deriv_x(x, k);

      const double x_tc      = get_tight_coupling_time(k);

      double dtheta1 = c*k/(3.*Hp)*theta0 -2.*c*k/(3.*Hp)*pi + c*k/(3.*Hp)*psi + dtau*(theta1 + vb/3.);
      double dpi     = 2.*c*k/(5.*Hp)*theta1 - 3.*c*k/(5.*Hp)*theta3 + 9./10.*dtau*pi;
      double dtheta3 = 3.*c*k/(7.*Hp)*pi - 4.*c*k/(7.*Hp)*theta4 + dtau*theta3;
      double ddpi    = 2.*c*k/(5.*Hp)*(dtheta1 - dHp*theta1/Hp) - 3.*c*k/(5.*Hp)*(dtheta3 - dHp*theta3/Hp) + 9./10.*(ddtau*pi + dtau*dpi);

      if (x < x_tc) {
        double dtheta0 = -c*k*theta1/Hp - dphi;
        double q       = (-((1.-R)*dtau + (1.+R)*ddtau)*(3.*theta1 + vb) - c*k*psi/Hp + (1.-dHp/Hp)*c*k/Hp*(-theta0 + 2.*pi) - c*k/Hp*dtheta0)/((1.+R)*dtau + dHp/Hp - 1.);
        dtheta1        = (q - dvb)/3.;

        dpi            = -20.*c*k/45.*(dtheta1/(Hp*dtau) - theta1/pow(Hp*dtau, 2.)*(Hp*ddtau + dHp*dtau));
        dtheta3        = -3.*c*k/7.*(dpi/(Hp*dtau) - pi/pow(Hp*dtau, 2.)*(Hp*ddtau + dHp*dtau));
        ddpi           = 2.*c*k/(5.*Hp)*(dtheta1 - dHp*theta1/Hp) - 3.*c*k/(5.*Hp)*(dtheta3 - dHp*theta3/Hp) + 9./10.*(ddtau*pi + dtau*dpi);
      }

      // Temperature source
      double sw = g_tilde*(theta0 + psi + pi/4.);
      double isw = exp(-tau)*(dpsi-dphi);
      double doppler = 1./(c*k)*(dHp*g_tilde*vb + Hp*dg_tilde*vb + Hp*g_tilde*dvb);
      double quadrupole = 3./pow(2.*c*k, 2.)*(g_tilde*pi*(dHp*dHp + Hp*ddHp)
                                              + 3.*Hp*dHp*(dg_tilde*pi + g_tilde*dpi)
                                              + Hp*Hp*(ddg_tilde*pi + 2.*dg_tilde*dpi + g_tilde*ddpi));

      ST_array[index] = sw + isw - doppler + quadrupole;
    }
  }

  // Spline the source functions
  ST_spline.create (x_array, k_array, ST_array, "Source_Temp_x_k");
  if(Constants.polarization){
    //SE_spline.create (x_array, k_array, SE_array, "Source_Pol_x_k");
  }

  Utils::EndTiming("source");
}

//====================================================
// The right hand side of the perturbations ODE
// in the tight coupling regime
//====================================================

// Derivatives in the tight coupling regime
int Perturbations::rhs_tight_coupling_ode(double x, double k, const double *y, double *dydx){

  //=============================================================================
  // Compute where in the y / dydx array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  const bool neutrinos          = Constants.neutrinos;
  const bool polarization       = Constants.polarization;

  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm_tc];
  const double &delta_b         =  y[Constants.ind_deltab_tc];
  const double &v_cdm           =  y[Constants.ind_vcdm_tc];
  const double &v_b             =  y[Constants.ind_vb_tc];
  const double &Phi             =  y[Constants.ind_Phi_tc];
  const double *Theta           = &y[Constants.ind_start_theta_tc];
  const double *Nu              = &y[Constants.ind_start_nu_tc];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm_tc];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab_tc];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm_tc];
  double &dv_bdx          =  dydx[Constants.ind_vb_tc];
  double &dPhidx          =  dydx[Constants.ind_Phi_tc];
  double *dThetadx        = &dydx[Constants.ind_start_theta_tc];
  double *dNudx           = &dydx[Constants.ind_start_nu_tc];

  // Useful quantities
  const double c          = Constants.c;
  const double Hp         = cosmo->Hp_of_x(x);
  const double dHp        = cosmo->dHpdx_of_x(x);
  const double H0         = cosmo->get_H0();
  const double a          = exp(x);
  const double OmegaR     = cosmo->get_OmegaR();
  const double OmegaNu    = cosmo->get_OmegaNu();
  const double OmegaCDM   = cosmo->get_OmegaCDM();
  const double OmegaB     = cosmo->get_OmegaB();
  const double dtau       = rec->dtaudx_of_x(x);
  const double ddtau      = rec->ddtauddx_of_x(x);
  const double R          = 4.*OmegaR/(3.*OmegaB*a);

  double Theta0           = Theta[0];
  double Theta1           = Theta[1];
  double Theta2           = -20.*c*k/(45.*Hp*dtau)*Theta1;
  if(polarization){
    Theta2                = -8.*c*k/(15.*Hp*dtau)*Theta1;
  }

  // SET: Scalar quantities (Phi, delta, v, ...)
  double Psi     = -Phi - 12.*pow(H0/(c*k*a), 2.)*(OmegaR*Theta2 + 0.);
  dPhidx         = Psi - pow(c*k/Hp, 2.)*Phi/3. + pow(H0/Hp, 2.)/2.*(OmegaCDM*delta_cdm/a + OmegaB*delta_b/a + 4.*OmegaR*Theta0/(a*a) + 0.);
  ddelta_cdmdx   = c*k*v_cdm/Hp - 3.*dPhidx;
  dv_cdmdx       = -v_cdm - c*k*Psi/Hp;
  ddelta_bdx     = c*k*v_b/Hp - 3.*dPhidx;

  // SET: Photon multipoles (Theta_ell)
  dThetadx[0]    = -c*k*Theta1/Hp - dPhidx; //Theta0'
  double q       = (-((1.-R)*dtau + (1.+R)*ddtau)*(3.*Theta1 + v_b) - c*k*Psi/Hp + (1.-dHp/Hp)*c*k/Hp*(-Theta0 + 2.*Theta2) - c*k/Hp*dThetadx[0])/((1.+R)*dtau + dHp/Hp - 1.);
  dv_bdx         = (-v_b - c*k*Psi/Hp + R*(q + c*k/Hp*(-Theta0 + 2.*Theta2) - c*k*Psi/Hp))/(1.+R);
  dThetadx[1]    = (q - dv_bdx)/3.; //Theta1'

  // SET: Neutrino mutlipoles (Nu_ell)
  if(neutrinos){
    // Not included
  }

  return GSL_SUCCESS;
}

//====================================================
// The right hand side of the full ODE
//====================================================

int Perturbations::rhs_full_ode(double x, double k, const double *y, double *dydx){

  //=============================================================================
  // Compute where in the y / dydx array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

  // Index and number of the different quantities
  const int n_ell_theta         = Constants.n_ell_theta;
  const int n_ell_thetap        = Constants.n_ell_thetap;
  const int n_ell_neutrinos     = Constants.n_ell_neutrinos;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm];
  const double &delta_b         =  y[Constants.ind_deltab];
  const double &v_cdm           =  y[Constants.ind_vcdm];
  const double &v_b             =  y[Constants.ind_vb];
  const double &Phi             =  y[Constants.ind_Phi];
  const double *Theta           = &y[Constants.ind_start_theta];
  const double *Theta_p         = &y[Constants.ind_start_thetap];
  const double *Nu              = &y[Constants.ind_start_nu];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm];
  double &dv_bdx          =  dydx[Constants.ind_vb];
  double &dPhidx          =  dydx[Constants.ind_Phi];
  double *dThetadx        = &dydx[Constants.ind_start_theta];
  double *dTheta_pdx      = &dydx[Constants.ind_start_thetap];
  double *dNudx           = &dydx[Constants.ind_start_nu];

  // Useful quantities
  const double c          = Constants.c;
  const double Hp         = cosmo->Hp_of_x(x);
  const double H0         = cosmo->get_H0();
  const double a          = exp(x);
  const double OmegaR     = cosmo->get_OmegaR();
  const double OmegaCDM   = cosmo->get_OmegaCDM();
  const double OmegaB     = cosmo->get_OmegaB();
  const double dtau       = rec->dtaudx_of_x(x);
  const double R          = 4.*OmegaR/(3.*OmegaB*a);
  const double eta        = cosmo->eta_of_x(x);

  double Theta0     = Theta[0];
  double Theta1     = Theta[1];
  double Theta2     = Theta[2];

  // SET: Scalar quantities (Phi, delta, v, ...)
  double Psi    = -Phi - 12.*pow(H0/(c*k*a), 2.)*(OmegaR*Theta2 + 0.);
  dPhidx        = Psi - pow(c*k/Hp, 2.)*Phi/3. + pow(H0/Hp, 2.)/2.*(OmegaCDM*delta_cdm/a + OmegaB*delta_b/a + 4.*OmegaR*Theta0/(a*a) + 0.);
  ddelta_cdmdx  = c*k/Hp*v_cdm - 3.*dPhidx;
  dv_cdmdx      = -v_cdm - c*k/Hp*Psi;
  ddelta_bdx    = c*k/Hp*v_b -3.*dPhidx;
  dv_bdx        = -v_b - c*k/Hp*Psi + dtau*R*(3.*Theta1 + v_b);

  // SET: Photon multipoles (Theta_ell)
  dThetadx[0]   = -c*k/Hp*Theta1 - dPhidx; //theta'0
  dThetadx[1]   = c*k/(3.*Hp)*Theta0 -2.*c*k/(3.*Hp)*Theta2 + c*k/(3.*Hp)*Psi + dtau*(Theta1 + v_b/3.);

  double Theta_prev;
  double Theta_curr;
  double Theta_next;
  for(int l = 2; l < n_ell_theta-1; l++){
    Theta_prev  = Theta[l-1];
    Theta_curr  = Theta[l];
    Theta_next  = Theta[l+1];
    dThetadx[l] = l*c*k/((2.*l+1.)*Hp) * Theta_prev - (l+1.)*c*k/((2.*l+1.)*Hp) * Theta_next + dtau*Theta_curr;
    if(l == 2){
      dThetadx[l] = l*c*k/((2.*l+1.)*Hp) * Theta_prev - (l+1.)*c*k/((2.*l+1.)*Hp) * Theta_next + 9./10.*dtau*Theta_curr;
    }
  }
  Theta_prev    = Theta[n_ell_theta-2];
  Theta_curr    = Theta[n_ell_theta-1];
  dThetadx[n_ell_theta-1]   = c*k*Theta_prev/Hp - c*(n_ell_theta)*Theta_curr/(Hp*eta) + dtau*Theta_curr;


  // SET: Photon polarization multipoles (Theta_p_ell)
  if(polarization){
    // Not included
  }

  // SET: Neutrino mutlipoles (Nu_ell)
  if(neutrinos){
    // Not included
  }

  return GSL_SUCCESS;
}

//====================================================
// Get methods
//====================================================

double Perturbations::get_delta_cdm(const double x, const double k) const{
  return delta_cdm_spline(x,k);
}
double Perturbations::get_delta_b(const double x, const double k) const{
  return delta_b_spline(x,k);
}
double Perturbations::get_v_cdm(const double x, const double k) const{
  return v_cdm_spline(x,k);
}
double Perturbations::get_v_b(const double x, const double k) const{
  return v_b_spline(x,k);
}
double Perturbations::get_Phi(const double x, const double k) const{
  return Phi_spline(x,k);
}
double Perturbations::get_Psi(const double x, const double k) const{
  return Psi_spline(x,k);
}
double Perturbations::get_Pi(const double x, const double k) const{
  return Pi_spline(x,k);
}
double Perturbations::get_Source_T(const double x, const double k) const{
  return ST_spline(x,k);
}
double Perturbations::get_Source_E(const double x, const double k) const{
  return SE_spline(x,k);
}
double Perturbations::get_Theta(const double x, const double k, const int ell) const{
  return Theta_spline[ell](x,k);
}
double Perturbations::get_Theta_p(const double x, const double k, const int ell) const{
  return Theta_p_spline[ell](x,k);
}
double Perturbations::get_Nu(const double x, const double k, const int ell) const{
  return Nu_spline[ell](x,k);
}

//====================================================
// Print some useful info about the class
//====================================================

void Perturbations::info() const{
  std::cout << "\n";
  std::cout << "Info about perturbations class:\n";
  std::cout << "x_start:       " << x_start                << "\n";
  std::cout << "x_end:         " << x_end                  << "\n";
  std::cout << "n_x:           " << n_x              << "\n";
  std::cout << "k_min (1/Mpc): " << k_min * Constants.Mpc  << "\n";
  std::cout << "k_max (1/Mpc): " << k_max * Constants.Mpc  << "\n";
  std::cout << "n_k:           " << n_k              << "\n";
  if(Constants.polarization)
    std::cout << "We include polarization\n";
  else
    std::cout << "We do not include polarization\n";
  if(Constants.neutrinos)
    std::cout << "We include neutrinos\n";
  else
    std::cout << "We do not include neutrinos\n";

  std::cout << "Information about the perturbation system:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm         << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab           << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm             << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb               << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi              << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta      << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta          << "\n";
  if(Constants.polarization){
    std::cout << "ind_start_thetap:   " << Constants.ind_start_thetap   << "\n";
    std::cout << "n_ell_thetap:       " << Constants.n_ell_thetap       << "\n";
  }
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu       << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos    << "\n";
  }
  std::cout << "n_ell_tot_full:     " << Constants.n_ell_tot_full       << "\n";

  std::cout << "Information about the perturbation system in tight coupling:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm_tc      << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab_tc        << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm_tc          << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb_tc            << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi_tc           << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta_tc   << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta_tc       << "\n";
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu_tc    << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos_tc << "\n";
  }
  std::cout << "n_ell_tot_tc:       " << Constants.n_ell_tot_tc         << "\n";
  std::cout << std::endl;
}

//====================================================
// Output some results to file for a given value of k
//====================================================

void Perturbations::output(const double k, const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts = 5000;
  auto x_array = Utils::linspace(x_start, 0.0, npts);

  std::cout << "Tight coupling end (k = " << k*Constants.Mpc << "/Mpc): " << get_tight_coupling_time(k) << std::endl;

  auto print_data = [&] (const double x) {
    double arg = k * (cosmo->eta_of_x(0.0) - cosmo->eta_of_x(x));
    fp << x                  << " ";
    fp << get_delta_cdm(x,k) << " ";
    fp << get_v_cdm(x,k)     << " ";
    fp << get_delta_b(x,k)   << " ";
    fp << get_v_b(x, k)      << " ";
    fp << get_Theta(x,k,0)   << " ";
    fp << get_Theta(x,k,1)   << " ";
    fp << get_Phi(x,k)       << " ";
    fp << get_Psi(x,k)       << " ";
    fp << get_Pi(x,k)        << " ";
    fp << get_Source_T(x,k)  << " ";
    fp << get_Source_T(x,k) * Utils::j_ell(5,   arg)           << " ";
    fp << get_Source_T(x,k) * Utils::j_ell(100,  arg)          << " ";
    fp << get_Source_T(x,k) * Utils::j_ell(500, arg)           << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}
