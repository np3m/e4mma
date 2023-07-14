#include "fore.h"

using namespace std;
using namespace o2scl;

int main(int argc, char *argv[]) {
  fore f;

  double Y_p = 0.105;
  double n_B = 0.1;
  double T = 29.286445646252375;

  double mu_n = 1314900.817136293;
  double mu_p = 966.3503814841417;
  double mu_pi = 818.3244888867657;
  double meff_n = 575.450474823058;
  double meff_p = 715.9852675623772;
  double U_n = -57.04326039625805;
  double U_p = -93.55858031069548;
  double m_n = 9.395654e+02;
  double m_p = 9.382721e+02;

  f.nucleon_mod.resize(2,2);
  f.nucleon_mod(0,0)=meff_n;
  f.nucleon_mod(0,1)=mu_n-U_n-m_n;
  f.nucleon_mod(1,0)=meff_p;
  f.nucleon_mod(1,1)=mu_p-U_p-m_p;

  f.load_pion();
  cout <<"loading pions done" << endl;
  f.verbose=1;

  f.single_point_data(Y_p, T, n_B, mu_n, mu_p, meff_n, meff_p);

  return 0;
}