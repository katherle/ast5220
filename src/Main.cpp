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

  // test with toy cosmology:
  // double h           = 0.7;
  // double OmegaB      = 0.05;
  // double OmegaCDM    = 0.45;

  //Recombination parameters
  double Yp          = 0.0;

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

  //=========================================================================
  // Module II
  //=========================================================================

  // Solve the recombination history
  RecombinationHistory rec(&cosmo, Yp);
  rec.solve();
  rec.info();

  // Output recombination quantities
  rec.output("recombination.txt");

  //=========================================================================
  // Module III
  //=========================================================================

  // Solve the perturbations
  Perturbations pert(&cosmo, &rec);
  pert.solve();
  pert.info();

  // Output perturbation quantities
  double kvalue = 0.1 / Constants.Mpc;
  pert.output(kvalue, "perturbations_k0.1.txt");

  kvalue = 0.01 / Constants.Mpc;
  pert.output(kvalue, "perturbations_k0.01.txt");

  kvalue = 0.001 / Constants.Mpc;
  pert.output(kvalue, "perturbations_k0.001.txt");

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
