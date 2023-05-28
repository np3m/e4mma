#include "fore.h"


// heaviside function
int fore::heaviside(double x) {
  if (x>=0) {
    return 1;
  }
  return 0;
}

funct fore::NR_fermion_distro(double m, double T, double mu) {
  funct f;
  if (T==0) {
    f=[this, m,T,mu] (double p) -> double { return heaviside((mu-p*p)/(2*m)); };
  } else {
    f=[this, m,T,mu] (double p) -> double { 
      double v = p*p/(2*m)-mu;
      if (-v/T>0){
        return 1/(exp(v/T)+1);
      } else {
        return exp(-v/T)/(1+exp(-v/T));
      } };
  }
  return f;
}
  
void fore::get_phase_shifts() {

  vector<string> fnames={"p11_central.o2", "p13_central.o2",
    "p31_central.o2","p33_central.o2", "s11_central.o2",
    "s31_central.o2"};

  double MN=m_p;
  double a=MN*MN+m_pi*m_pi;
  double b=4*m_pi*m_pi*MN*MN;

  for(size_t i=0;i<6;i++) {

    if (verbose>-1){
      cout << "Reading file: " << fnames[i] << endl;
    }
      
    hdf_file hf;
    hf.open(((string)"piN_phase_shifts/")+fnames[i]);
    table_units<> t;
    hdf_input(hf,t);
    hf.close();

    for(size_t j=0;j<100;j++) {
      t.set("c1",j,t.get("c1",j)*1.0e3);
      t.set("c2",j,t.get("c2",j)*o2scl_const::pi/180.0);
    }
      
    if (i==0) {
      com_energy.resize(100);
      for(size_t j=0;j<100;j++) {
        com_energy(j) = (t.get("c1",j));
      }
      com_momentum.resize(100);
      for(size_t j=0;j<100;j++) {
        if (j==0) {
          com_momentum(j) = 0.0;
        } else {
          com_momentum(j) = sqrt(fabs(pow(com_energy[j]*
                                               com_energy[j]-a,2.0)-b))/
                                 (2*com_energy[j]);
        }
      }
    }

    for(size_t j=0;j<100;j++) {
      if (i==0) d11m.push_back(t.get("c2",j));
      else if (i==1) d11.push_back(t.get("c2",j));
      else if (i==2) d31m.push_back(t.get("c2",j));
      else if (i==3) d31.push_back(t.get("c2",j));
      else if (i==4) d10.push_back(t.get("c2",j));
      else d30.push_back(t.get("c2",j));
    }
  }

  // Also save everything in a table to check later or plot using o2graph
  if(false) {
    hdf_file hf1;
    hf1.open_or_create("phase_shifts.o2");
    table_units<> t1;
    t1.set_nlines(com_energy.size());
    t1.line_of_names((std::string)"com_mom com_e d10 d11 d11m d30 d31 d31m");

    for(size_t i=0; i<100; i++){
      t1.set("com_e", i, com_energy[i]); t1.set("com_mom", i, com_momentum[i]);
      t1.set("d10", i, d10[i]); t1.set("d11", i, d11[i]);
      t1.set("d11m", i, d11m[i]); t1.set("d30", i, d30[i]);
      t1.set("d31", i, d31[i]); t1.set("d31m", i, d31m[i]);
    }
    hdf_output(hf1,t1,"phase_shifts");
    hf1.close();
  }
  
  return;
}

void fore::interp_phase_shift_sum(std::vector<double> params) {
  double a10n=params[0]; double a11n=params[1]; double a30n=params[2];
  double a31n=params[3]; double a10p=params[0]; double a11p=params[1];
  double a30p=params[2]; double a31p=params[3];

  int p_wave_gate=0;
  if (include_p_wave) {
    p_wave_gate=1;
  }

  d3_sum.resize(100); d3d1_sum.resize(100);

  // Neutrons
  for(size_t i=0;i<100;i++) {
    d3_sum(i) = a30n*2*d30[i]+p_wave_gate*a31n*3*(d31[i]+d31m[i]);
  }
  neutron_sum.set(100,com_momentum,d3_sum,itp_linear);

  // Protons
  for(size_t i=0;i<100;i++) {
    d3d1_sum(i) = 2*(a30p*d30[i]+a10p*d10[i])+
                       p_wave_gate*3*(a31p*(d31[i]+d31m[i])+
                                      a11p*(d11[i]*d11m[i]));
  }
  proton_sum.set(100,com_momentum,d3d1_sum,itp_linear);

  neutron_phase_shift_sum = [this] (double p_cm)-> double {
    double values;
    
    if (p_cm<com_momentum[0]) { p_cm=com_momentum[0]; }
    if (p_cm>com_momentum[com_momentum.size()-1]) {
      p_cm=com_momentum[com_momentum.size()-1]; }
    values=neutron_sum.eval(p_cm);
    
    if (!const_after_data) {
      if (p_cm>com_momentum[com_momentum.size()-1]) {
        values=0.0; }
    }
    return values; };

  proton_phase_shift_sum = [this] (double p_cm)-> double {
    double values;
    
    if (p_cm<com_momentum[0]) { p_cm=com_momentum[0]; }
    if (p_cm>com_momentum[com_momentum.size()-1]) {
      p_cm=com_momentum[com_momentum.size()-1]; }
    values=proton_sum.eval(p_cm);
    
    if (!const_after_data) {
      if (p_cm>com_momentum[com_momentum.size()-1]) {
        values=0.0; }
    }
    return values; };
  return;
}
  
double fore::inner_integral(double p, double k, funct phase_shift_sum,
                      double m_N, double last_val) {
  // Reduced mass
  //cout << "ii: m_pi, m_N: " << m_pi << " " << m_N << endl;
  double m_bar=(m_pi*m_N)/(m_pi+m_N); 
  // Fore flips the bounds to keep the larger number as the upper bound.
  // Minus sign has been included in return value.
  double lower_bound=m_bar*fabs(p/m_pi-k/m_N);
  double upper_bound=m_bar*(p/m_pi+k/m_N);
  
  double int_val=0.0;
  if (lower_bound>=com_momentum[com_momentum.size()-1]) {
    int_val=upper_bound*last_val-lower_bound*last_val;
  } else if (upper_bound>com_momentum[com_momentum.size()-1]) {
    double integral = fixed_quad.integ(phase_shift_sum, lower_bound,
                                    com_momentum[com_momentum.size()-1]);
    int_val = integral;
    int_val += upper_bound*last_val - com_momentum[com_momentum.size()-1]*last_val;
      
  } else {
    double integral = fixed_quad.integ(phase_shift_sum, lower_bound,
                                    upper_bound);
    int_val = integral;
  }
  //cout << "ii: " << -1*int_val << endl;
  return -1*int_val;
}

double fore::self_energy(double p, std::vector<double> params, ubmatrix nuc_mod, double T) {

  double last_val_n = neutron_phase_shift_sum(com_momentum[com_momentum.size()-1]);
  funct fn = NR_fermion_distro(nuc_mod(0,0),T,nuc_mod(0,1));

  double last_val_p = proton_phase_shift_sum(com_momentum[com_momentum.size()-1]);
  funct fp = NR_fermion_distro(nuc_mod(1,0),T,nuc_mod(1,1));

  if(verbose>1){
    cout << "se: fn, fp: " << fn(p) << " " << fp(p) << endl;
  }
  
  funct fiin=[this,p,fn,last_val_n] (double k) -> double { return k*fn(k)*
                    inner_integral(p,k,neutron_phase_shift_sum,m_n,last_val_n); };
  funct fiip=[this,p,fp,last_val_p] (double k) -> double { return k*fp(k)*
                    inner_integral(p,k,proton_phase_shift_sum,m_p,last_val_p); };

  if(verbose>1){
    cout << "se: fiin, fiip: " << fiin(p) << " " << fiip(p) << endl;
  }
  // The result and the uncertainty
  double res_n, err_n, res_p, err_p, tot, res_n1, err_n1, res_p1, err_p1;
  gu.tol_rel=1.49e-04;
  gu.tol_abs=0.01;
  //kron.tol_rel=1.49e-04;
  //kron.tol_abs=0.01;

  int ret_n = gu.integ_err(fiin, 0.0, 0.0, res_n, err_n);
  //int ret_n1 = kron.integ_err(fiin, 0.0, 1.0e8, res_n1, err_n1);
  int ret_p = gu.integ_err(fiip, 0.0, 0.0, res_p, err_p);
  //int ret_p1 = kron.integ_err(fiip, 0.0, 1.0e8, res_p1, err_p1);

  if(verbose>0){
    cout << "Result(gu): " << res_n << " Uncertainty: " << err_n << endl;
    //cout << "Result(kron): " << res_n1 << " Uncertainty: " << err_n1 << endl;
    cout << "Result(gu): " << res_p << " Uncertainty: " << err_p << endl;
    //cout << "Result(kron): " << res_p1 << " Uncertainty: " << err_p1 << endl;
    //cout << "Number of iterations: " << gu.last_iter << endl;
  }

  double m_bar = (m_pi*m_n)/(m_pi+m_n);  // Reduced mass
  tot = ((m_n+m_pi)/(2*pi*p*m_bar*m_bar))*res_n;
  double m_bar_p = (m_pi*m_p)/(m_pi+m_p);  // Reduced mass
  tot += ((m_p+m_pi)/(2*pi*p*m_bar_p*m_bar_p))*res_p;

  if(verbose>0){  
    cout << "se: p, se: " << p << " " << tot << endl;
  }
  return tot;
}

funct fore::self_energy_interp(vector<double> params, ubmatrix nuc_mod, double T) {

  se_list.resize(p_list.size()); p_list_New.resize(p_list.size());

  for (size_t i=0; i<p_list.size(); i++) {
    p_list_New[i] = p_list[i];
    se_list[i] = self_energy(p_list[i], params, nuc_mod, T);
    //cout << "se: " << p_list[i] << " " << se_list[i] << endl;
  }

  if(false) {
    hdf_file hf2;
    hf2.open_or_create("self_energy.o2");
    table_units<> t2;
    t2.set_nlines(p_list.size());
    t2.line_of_names((std::string)"p se");

    for(size_t i=0; i<p_list.size(); i++){
      t2.set("p", i, p_list[i]); t2.set("se", i, se_list[i]);
    }
    hdf_output(hf2,t2,"self_energy");
    hf2.close();
  }

  interp_se.set(p_list.size(), p_list_New, se_list, itp_cspline);

  funct se_func=[this] (double p) -> double {
    if (p < p_list[0]) {
      p = p_list[0];
    } else if (p > p_list[p_list.size()-1]) {
      return 0;
    }
    return interp_se.eval(p);
  };
  return se_func;
}

funct fore::rel_boson_distro(double m, double T, double mu, funct sigma) {
  funct distro = [this,m,T,mu,sigma](double p) -> double { 
    double v = sqrt(p*p+m*m)+sigma(p)-mu;
    if (-v/T>0){
      return 1/(exp(v/T)-1);
    } else {
      return exp(-v/T)/(1-exp(-v/T));
    } };
  return distro;
}

bool fore::condensation_exists(funct sigma_pi, double mu) {
  funct pi_en = [this,sigma_pi](double p) -> double {return sqrt(p*p+m_pi*m_pi) + sigma_pi(p); };
  double low_p_min, sigma_pi_min, high_p_min, x3; double x1=0.0; double x2=200.0;

  mcn.min_bkt(x1, 0, 1.0e4, low_p_min, pi_en);
  mcn.min_bkt(x2, 0, 1.0e4, sigma_pi_min, sigma_pi);
  if (verbose>0) {
    std::cout << "pion energy min: " << low_p_min << " at momentum: " << x1 << std::endl;
    std::cout << "self energy min: " << sigma_pi_min << " at momentum: " << x2 << std::endl;
    std::cout << "self energy at p=200: " << sigma_pi(200) << std::endl;
  }

  if (sigma_pi_min==0.0){ 
    std::cout << "condensation finding error. bad assumptions" << std::endl;
    for (size_t i=0;i<p_list.size();i++) {
      cout << p_list[i] << sigma_pi(p_list[i]) << endl;
    }
  }

  mcn.min_bkt(x3,x1,x2, high_p_min, pi_en);
  if (verbose>0) { 
    std::cout << "high p min: " << high_p_min << " at " << x3 << std::endl;
    std::cout << "high p min: " << pi_en(0) << " at 0" << std::endl;
    std::cout << "high p min: " << pi_en(x2) << " at " << x2 << std::endl;
  }
  double global_min = min(low_p_min, high_p_min);
  if (global_min<=mu) { return true;}
  else {return false;}
}

double fore::rel_pion_number_density(double T, double mu, vector<double> params, ubmatrix nuc_mod) {
    
  funct sigma_pi = self_energy_interp(params, nuc_mod, T);

  if (condensation_exists(sigma_pi,mu)) {return 0.0;}

  funct integrand = [this,T,mu,sigma_pi](double p) -> double { 
    return p*p*rel_boson_distro(m_pi,T,mu,sigma_pi)(p); };

  gu.tol_abs=0.0; gu.tol_rel=1.49e-08;  
  double res, err;
  int ret = gu.integ_err(integrand,0.0,0.0,res,err);

  if (fabs(res)<fabs(err)) {
    std::cout << "Large integral error in rel_pion_number_density: " << std::endl;
  }

  return res/(2*pi*pi);
}

double fore::non_int_rel_pion_number_density(double T, double mu) {
  if (m_pi<mu){ return 0.0;}
  funct temp = [this](double p) -> double {return 0; };

  funct integrand = [this,T,mu,temp](double p) -> double { return p*p*rel_boson_distro(m_pi,T,mu,temp)(p); };

  gu.tol_abs=0.0; gu.tol_rel=1.49e-08;
  double res, err;
  int ret = gu.integ_err(integrand,0.0,0.0,res,err);

  if (fabs(res)<fabs(err)) {
    std::cout << "Large integral error in non_int_rel_pion_number_density: " << std::endl;
  }

  return res/(2*pi*pi);
}

double fore::pion_entropy(double T, double mu, vector<double> params, ubmatrix nuc_mod) {

  funct sigma_pi = self_energy_interp(params, nuc_mod, T);
  funct distro = rel_boson_distro(m_pi, T, mu, sigma_pi);

  funct integrand = [this,distro](double p) -> double {
    double d =distro(p);
    if(0<d){
      return p*p*(d*log(d)-(1+d)*log(1+d));
    }
    return 0;
  };

  gu.tol_abs=0.0; gu.tol_rel=1.49e-08;
  double res,err;
  int ret = gu.integ_err(integrand,0.0,0.0,res,err);

  if (fabs(res)<fabs(err)) {
    std::cout << "Large integral error in pion_entropy: " << std::endl;
  }

  return -res/(2*pi*pi);
}

double fore::pion_energy(double T, double mu, vector<double> params, ubmatrix nuc_mod) {

  funct sigma_pi = self_energy_interp(params, nuc_mod, T);

  funct integrand = [this,T,mu,sigma_pi](double p) -> double {
    return (p*p*(sqrt(p*p+m_pi*m_pi)+sigma_pi(p))*
            rel_boson_distro(m_pi,T,mu,sigma_pi)(p));
  };

  gu.tol_abs=0.0; gu.tol_rel=1.49e-08;
  double res,err;
  int ret = gu.integ_err(integrand,0.0,0.0,res,err);

  if (fabs(res)<fabs(err)) {
    std::cout << "Large integral error in pion_energy: " << std::endl;
  }

  return res/(2*pi*pi);
}

double fore::boson_entropy(funct distro) {

  funct integrand = [this,distro](double p) -> double {
    double d =distro(p);
    if(0<d){
      return p*p*(d*log(d)-(1+d)*log(1+d));
      }
    return 0;
  };

  gu.tol_abs=0.0; gu.tol_rel=1.49e-08;
  double res,err;
  int ret = gu.integ_err(integrand,0.0,0.0,res,err);

  if (fabs(res)<fabs(err)) {
    std::cout << "Large integral error in boson_entropy: " << std::endl;
  }

  return -res/(2*pi*pi);
}

void fore::single_point_data(double Y_p, double T, double n_B, double mu_n, 
                double mu_p, double meff_n, double meff_p, double U_n, double U_p){

  Y_pi=0.0, e_pi=0.0, s_pi=0.0, press_pi=0.0;

  funct sigma_pi = self_energy_interp(pseudo_pot_params, nucleon_mod, T);

  if (mu_pi>(m_pi+sigma_pi(0))){
    flag = 0;
  } else {
    Y_pi = rel_pion_number_density(T, mu_pi, pseudo_pot_params, nucleon_mod);

    if (Y_pi!=0) {
      Y_pi = Y_pi/n_B;
    }

    if (Y_pi==0){
      flag = 0;
    } else if (Y_pi>0) {
      e_pi = pion_energy(T, mu_pi, pseudo_pot_params, nucleon_mod);
      s_pi = pion_entropy(T, mu_pi, pseudo_pot_params, nucleon_mod)/n_B;

      press_pi = T*s_pi*n_B - e_pi + Y_pi*n_B*mu_pi;

      e_pi = e_pi/(hc_mev_fm*hc_mev_fm*hc_mev_fm);
      press_pi = press_pi/(hc_mev_fm*hc_mev_fm*hc_mev_fm);
    } else {
      flag = -1;
    }
    cout << "Y_pi, e, s, P: " << Y_pi << " " << e_pi << " " << s_pi << " " << press_pi << endl;
  }
  return;
}

void fore::calc_mu(boson &b, fermion &n, fermion &p, double T, double n_B) {

  // Convert units from 1/fm to MeV to use here..
  T = T*hc_mev_fm;
  mu_n = n.mu*hc_mev_fm;       // Neutron Chemical potential
  mu_p = p.mu*hc_mev_fm;       // Proton chemical potential
  mu_pi = b.mu*hc_mev_fm;      // Pion chemical potential
  meff_n = n.ms*hc_mev_fm;     // Neutron effective mass
  meff_p = p.ms*hc_mev_fm;     // Proton effective mass
  U_n=(n.mu+n.nu)*hc_mev_fm;   // Neutron single particle potential
  U_p=(p.mu+p.nu)*hc_mev_fm;   // Proton single particle potential

  m_n = n.m*hc_mev_fm;         // Neutron rest mass
  m_p = p.m*hc_mev_fm;         // Proton rest mass
  m_pi = b.m*hc_mev_fm;        // Pion rest mass

  if (verbose>0) {
    cout << "n_B, T: " << n_B << " " << T << endl;
    cout << "mu_n, mu_p, mu_pi: " << mu_n << " " << mu_p << " " << mu_pi << endl;
    cout << "meff_n, meff_p: " << meff_n << " " << meff_p << endl;
    cout << "U_n, U_p: " << U_n << " " << U_p << endl;
    cout << "m_n, m_p, m_pi:  " << m_n << " " << m_p << " " << m_pi << endl;
  }

  double n_0 = .16;
  double n_0_MeV3 = n_0*(hc_mev_fm*hc_mev_fm*hc_mev_fm);

  nucleon_mod.resize(2,2);
  nucleon_mod(0,0)=meff_n;
  nucleon_mod(0,1)=mu_n-U_n-m_n;
  nucleon_mod(1,0)=meff_p;
  nucleon_mod(1,1)=mu_p-U_p-m_p;

if (verbose>0) {
  cout << "nuc_mod(0,1): " << nucleon_mod(0,0) << " " << nucleon_mod(0,1) << endl;
  cout << "nuc_mod(2,3): " << nucleon_mod(1,0) << " " << nucleon_mod(1,1) << endl;
  cout << "p_list.size(): " << p_list.size() << endl;
  cout << "p_list: " << p_list[0] << " - " << p_list[p_list.size()-1] << std::endl;
}

  //con_chk(T,n_B);
  //exit(-1);

  single_point_data(Y_p, T, n_B, mu_n, mu_p, meff_n, meff_p, U_n, U_p);
  b.n=Y_pi;
  b.pr=press_pi;
  b.ed=e_pi;
  b.en=s_pi;
  return;
}

void fore::load_pion(){
  verbose=0;
  pseudo_pot_params = {1, 0.3081465, 1, 0.3081465};

  include_p_wave=true;
  const_after_data=true;
  p_list.clear();

  for(double x=0.01;x<500.01;x+=(500.0-0.01)/29.0) {
    p_list.push_back(x);
  }
  for(double x=520.0;x<2000.01;x+=(2000-520)/9.0) {
    p_list.push_back(x);
  }
  for(double x=2150.0;x<10000.01;x+=(1e4-2150)/9.0) {
    p_list.push_back(x);
  }
  get_phase_shifts();
  interp_phase_shift_sum(pseudo_pot_params);
  return;
}

void fore::con_chk(double T, double n_B){

  cout << "starting condensation check." << endl;

  hdf_file hf;
  hf.open("data/fid_3_5_22.o2");
  size_t n_nB, n_Ye, n_T;

  if (verbose>2) cout << "Reading n_nB." << endl;
  hf.get_szt("n_nB",n_nB);
  if (verbose>2) cout << "Reading n_Ye." << endl;
  hf.get_szt("n_Ye",n_Ye);
  if (verbose>2) cout << "Reading n_T." << endl;
  hf.get_szt("n_T",n_T);
  if (n_nB==0 || n_Ye==0 || n_T==0) {
    O2SCL_ERR("One of the grid counts is zero.",o2scl::exc_efailed);
  }
  cout << "Done with numbers" << endl;

  std::vector<double> nB_grid, Ye_grid, T_grid;

  if (verbose>2) cout << "Reading nB_grid." << endl;
  hf.getd_vec("nB_grid",nB_grid);
  if (verbose>2) cout << "Reading Ye_grid." << endl;
  hf.getd_vec("Ye_grid",Ye_grid);
  if (verbose>2) cout << "Reading T_grid." << endl;
  hf.getd_vec("T_grid",T_grid);

  cout << "Done with grids" << endl;

  o2scl::tensor_grid<> pi_con;

  size_t st[3]={n_nB,n_Ye,n_T};
  vector<vector<double> > grid={nB_grid,Ye_grid,T_grid};
  pi_con.resize(3,st);
  pi_con.set_grid(grid);

  for (size_t i=0;i<n_nB;i+=10) {
    for (size_t j=0;j<n_Ye;j+=7) {
      for (size_t k=0;k<n_T;k+=15) {
        vector<size_t> ix={i,j,k};
        bool status = condensation_exists(self_energy_interp(pseudo_pot_params, nucleon_mod, T), mu_pi);
        if (status==true) { pi_con.set(ix,2); }
        else pi_con.set(ix,1);
        cout << "index: " << i << " " << j << " " << k << " " << endl;
      }
    }
  }
  
  hf.open_or_create("con_chk");
  hdf_output(hf,pi_con,"pi_con");
  hf.close();

  return;
}
         