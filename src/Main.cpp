#include "Utils.h"
#include "BackgroundCosmology.h"
#include "RecombinationHistory.h"
#include "Perturbations.h"
#include "PowerSpectrum.h"

int main(int argc, char **argv){
  Utils::StartTiming("Everything");

  //=========================================================================
  // Parameters
  //=========================================================================

  // Background parameters
  double h           = 0.67;
  double OmegaB      = 0.05;
  double OmegaCDM    = 0.267;
  double OmegaK      = 0.0;
  double Neff        = 0.0;
  double TCMB        = 2.7255;

  // Recombination parameters
  double Yp          = 0.245;

  // Power-spectrum parameters
  double A_s         = 2.1e-9;
  double n_s         = 0.965;
  double kpivot_mpc  = 0.05;

  //=========================================================================
  // Module I
  //=========================================================================

  // Set up and solve the background
  BackgroundCosmology cosmo(h, OmegaB, OmegaCDM, OmegaK, Neff, TCMB);
  cosmo.solve();
  cosmo.info();

  // Output background evolution quantities
  cosmo.output("cosmology.txt");

  //test
  double x = Constants.x_end;
  std::cout << "Testing:" << std::endl;
  std::cout << "Should equal 1: " << cosmo.Hp_of_x(x)/(pow(M_E, x)*cosmo.H_of_x(x)) << std::endl;
  std::cout << "Should equal 1: " << cosmo.get_detadx(0.0)*cosmo.Hp_of_x(0.0)/Constants.c << std::endl;
  std::cout << "Should equal 1: " << cosmo.dHpdx_of_x(x)/cosmo.Hp_of_x(x) << std::endl;
  std::cout << "Should equal -1/2: " << cosmo.dHpdx_of_x(-3.)/cosmo.Hp_of_x(-3.) << std::endl;
  std::cout << "Should equal -1: " << cosmo.dHpdx_of_x(Constants.x_start)/cosmo.Hp_of_x(Constants.x_start) << std::endl;
  std::cout << "Should equal 1: " << cosmo.ddHpddx_of_x(x)/cosmo.Hp_of_x(x) << std::endl;
  std::cout << "Should equal 1/4: " << cosmo.ddHpddx_of_x(-3.)/cosmo.Hp_of_x(-3.) << std::endl;
  std::cout << "Should equal 1: " << cosmo.ddHpddx_of_x(Constants.x_start)/cosmo.Hp_of_x(Constants.x_start) << std::endl;

  // Remove when module is completed
  return 0;

  //=========================================================================
  // Module II
  //=========================================================================

  // Solve the recombination history
  RecombinationHistory rec(&cosmo, Yp);
  rec.solve();
  rec.info();

  // Output recombination quantities
  rec.output("recombination.txt");

  // Remove when module is completed
  return 0;

  //=========================================================================
  // Module III
  //=========================================================================

  // Solve the perturbations
  Perturbations pert(&cosmo, &rec);
  pert.solve();
  pert.info();

  // Output perturbation quantities
  double kvalue = 0.01 / Constants.Mpc;
  pert.output(kvalue, "perturbations_k0.01.txt");

  // Remove when module is completed
  return 0;

  //=========================================================================
  // Module IV
  //=========================================================================

  PowerSpectrum power(&cosmo, &rec, &pert, A_s, n_s, kpivot_mpc);
  power.solve();
  power.output("cells.txt");

  // Remove when module is completed
  return 0;

  Utils::EndTiming("Everything");
}
