#include "fore.h"


  // heaviside function
  int fore::heaviside(double x) {
    if (x>=0) {
      return 1;
    }
    return 0;
  }

  // Nucleon distribution function for a non-rel fermion.
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

      cout << "Reading file: " << fnames[i] << endl;
      
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
    return -1*int_val;
  }
  
  // Computes self energy for a pion momentum p

  double fore::self_energy(double p, std::vector<double> params, ubmatrix nuc_mod, double T) {

    double last_val_n = neutron_phase_shift_sum(com_momentum[com_momentum.size()-1]);
    funct fn = NR_fermion_distro(nuc_mod(0,0),T,nuc_mod(0,1));

    double last_val_p = proton_phase_shift_sum(com_momentum[com_momentum.size()-1]);
    funct fp = NR_fermion_distro(nuc_mod(1,0),T,nuc_mod(1,1));

    funct fiin=[this,p,fn,last_val_n] (double k) -> double { return k*fn(k)*
                      inner_integral(p,k,neutron_phase_shift_sum,m_n,last_val_n); };
    funct fiip=[this,p,fp,last_val_p] (double k) -> double { return k*fp(k)*
                      inner_integral(p,k,proton_phase_shift_sum,m_p,last_val_p); };
    
    // The result and the uncertainty
    double res_n, err_n, res_p, err_p, tot;
    quad.tol_rel=1.0e-4;

    int ret1 = quad.integ_err(fiin, 0.0, 0.0, res_n, err_n);
    int ret2 = quad.integ_err(fiip, 0.0, 0.0, res_p, err_p);

    double m_bar = (m_pi*m_n)/(m_pi+m_n);  // Reduced mass
    tot = ((m_n+m_pi)/(2*pi*p*m_bar*m_bar))*res_n;
    double m_bar_p = (m_pi*m_p)/(m_pi+m_p);  // Reduced mass
    tot += ((m_p+m_pi)/(2*pi*p*m_bar_p*m_bar_p))*res_p;
    
    return tot;
  }

  // Sets up the interpolation vector and returns self energy for a given pion 
  // momentum through interpolation

  funct fore::self_energy_interp(vector<double> params, ubmatrix nuc_mod, double T) {

    se_list.resize(p_list.size()); p_list_New.resize(p_list.size());

    for (size_t i=0; i<p_list.size(); i++) {
      p_list_New[i] = p_list[i];
      se_list[i] = self_energy(p_list[i], params, nuc_mod, T);
      cout << "se: " << p_list[i] << " " << se_list[i] << endl;
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
  //---------------------------------------------------------------------------------------
  //---------------------------------------------------------------------------------------
  // Nucleon distribution function for a non-rel boson.

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

  //Return True if the given self energy has a minimum less than $\mu$.

  //The current self energy tends to have a single minimum in the
  //hundreds of MeV. This translates to a minimum in the pion energy
  //either at p=0 or near the minimum of $\Sigma_\pi$. So we will search
  //for the location of the potential minimum at finite momentum
  //starting the minimum function near the minimum of $\Sigma_\pi$ which
  //should only have one minima and another to search near zero.

  bool fore::condensation_exists(funct sigma_pi, double mu) {

    funct pion_energy = [this,sigma_pi](double p) -> double {return sqrt(p*p+m_pi*m_pi) + sigma_pi(p); };
    
    double low_p_min, sigma_pi_min, high_p_min; double x1=0.0; double x2=200.0;

    mcn.min_bkt(x1, 0.0, numeric_limits<double>::infinity(), low_p_min, pion_energy);
    mcn.min_bkt(x2, 0.0, numeric_limits<double>::infinity(), sigma_pi_min, sigma_pi);

    if (sigma_pi_min==0.0){ 
      std::cout << "condensation finding error. bad assumptions" << std::endl; }

    mcn.min_bkt(sigma_pi_min,0.0,numeric_limits<double>::infinity(), high_p_min, pion_energy); 
    double global_min = min(low_p_min, high_p_min);
    if (global_min<=mu) { return true;}
    else {return false;}
  }

  // Computes the number density of interacting pions

  double fore::rel_pion_number_density(double T, double mu, vector<double> params, ubmatrix nuc_mod) {
    
    funct sigma_pi = self_energy_interp(params, nuc_mod, T);

    if (condensation_exists(sigma_pi,mu)) {return 0.0;}

    funct integrand = [this,T,mu,sigma_pi](double p) -> double { 
      return p*p*rel_boson_distro(m_pi,T,mu,sigma_pi)(p); };
    
    double res, err;
    int ret = quad.integ_err(integrand,0.0,0.0,res,err);

    if (fabs(res)<fabs(err)) {
      std::cout << "Large integral error in rel_pion_number_density: " << std::endl;
    }

    return res/(2*pi*pi);
  }

  // Computes the number density of non-interacting pions

  double fore::non_int_rel_pion_number_density(double T, double mu) {
    if (m_pi<mu){ return 0.0;}
    funct temp = [this](double p) -> double {return 0; };

    funct integrand = [this,T,mu,temp](double p) -> double { return p*p*rel_boson_distro(m_pi,T,mu,temp)(p); };

    double res, err;
    int ret = quad.integ_err(integrand,0.0,0.0,res,err);

    if (fabs(res)<fabs(err)) {
      std::cout << "Large integral error in non_int_rel_pion_number_density: " << std::endl;
    }

    return res/(2*pi*pi);
  }

  // Calculate entropy of bosons given their distribution function.

  //  Input is unitless distribution function which takes in momentum in
  //  MeV Output is MeV^3

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

    double res,err;
    int ret = quad.integ_err(integrand,0.0,0.0,res,err);

    if (fabs(res)<fabs(err)) {
      std::cout << "Large integral error in pion_entropy: " << std::endl;
    }

    return -res/(2*pi*pi);
  }

  // Energy of pions in medium given by nucleon_mod using 
  // pseudo-potential.

  //  Function includes interaction energy.

  double fore::pion_energy(double T, double mu, vector<double> params, ubmatrix nuc_mod) {

    funct sigma_pi = self_energy_interp(params, nuc_mod, T);

    funct integrand = [this,T,mu,sigma_pi](double p) -> double {
      return (p*p*(sqrt(p*p+m_pi*m_pi)+sigma_pi(p))*
              rel_boson_distro(m_pi,T,mu,sigma_pi)(p));
    };

    double res,err;
    int ret = quad.integ_err(integrand,0.0,0.0,res,err);

    if (fabs(res)<fabs(err)) {
      std::cout << "Large integral error in pion_energy: " << std::endl;
    }

    return res/(2*pi*pi);
  }

  // Computes non-interacting boson entropy without self-energy contributions
  // given their distribution function.

  double fore::boson_entropy(funct distro) {

    funct integrand = [this,distro](double p) -> double {
      double d =distro(p);
      if(0<d){
        return p*p*(d*log(d)-(1+d)*log(1+d));
        }
      return 0;
    };

    double res,err;
    int ret = quad.integ_err(integrand,0.0,0.0,res,err);

    if (fabs(res)<fabs(err)) {
      std::cout << "Large integral error in boson_entropy: " << std::endl;
    }

    return -res/(2*pi*pi);
  }

  void fore::single_point_data(double Y_p, double T, double n_B, double mu_n, 
                  double mu_p, double meff_n, double meff_p, double U_n, double U_p){
    funct sigma_pi = self_energy_interp(pseudo_pot_params, nucleon_mod, T);
    if (mu_pi>(m_pi+sigma_pi(0))){
      flag = 0;
      cout << "flag 0: "  << mu_pi << " " << m_pi << " " << sigma_pi(0) << " " << flag << endl;
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

        e_pi = e_pi/(hbar_c*hbar_c*hbar_c);
        press_pi = press_pi/(hbar_c*hbar_c*hbar_c);
      } else {
        flag = -1;
      }
      cout << "e s P: " << e_pi << " " << s_pi << " " << press_pi << " " << flag << endl;
    }
    return;
  }

  void fore::calc_mu(boson &b, fermion &n, fermion &p, double T, double n_B) {

    // Convert units from 1/fm to MeV to use here..

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

    cout << "n_B, T: " << n_B << " " << T << endl;
    cout << "mu_n, mu_p, mu_pi: " << mu_n << " " << mu_p << " " << mu_pi << endl;
    cout << "meff_n, meff_p: " << meff_n << " " << meff_p << endl;
    cout << "U_n, U_p: " << U_n << " " << U_p << endl;
    cout << "m_n, m_p, m_pi:  " << m_n << " " << m_p << " " << m_pi << endl;

    hbar = 6.582119514e-22;

    pseudo_pot_params = {1, 0.3081465, 1, 0.3081465};

    include_p_wave=true;
    const_after_data=true;

    for(double x=0.01;x<500.01;x+=(500.0-0.01)/29.0) {
        p_list.push_back(x);
    }
    for(double x=520.0;x<2000.01;x+=(2000-520)/9.0) {
        p_list.push_back(x);
    }
    for(double x=2150.0;x<10000.01;x+=(1e4-2150)/9.0) {
        p_list.push_back(x);
    }

    double n_0 = .16;
    double n_0_MeV3 = n_0*(hc_mev_fm*hc_mev_fm*hc_mev_fm);
    
    
    nucleon_mod.resize(2,2);
    nucleon_mod(0,0)=n.ms*hc_mev_fm;
    nucleon_mod(0,1)=(n.mu-U_n-n.m)*hc_mev_fm;
    nucleon_mod(1,0)=p.ms*hc_mev_fm;
    nucleon_mod(1,1)=(p.mu-U_p-p.m)*hc_mev_fm;

    cout << "p_list.size(): " << p_list.size() << endl;

    get_phase_shifts();
    
    interp_phase_shift_sum(pseudo_pot_params);

    single_point_data(Y_p, T, n_B, mu_n, mu_p, meff_n, meff_p, U_n, U_p);
    b.n=Y_pi;
    b.pr=press_pi;
    b.ed=e_pi;
    b.en=s_pi;
    return;
  }
         