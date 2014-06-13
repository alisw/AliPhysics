{
  if (!Setup::stage == 0) { //generation step configurables
    //look at higgs decays
    Setup::decay_particle=25;
    Setup::mass_power=1;
    Setup::mass_scale_on=true;
    
    // Setup histograms
    int n_bins=60;
    double default_min_bin=0.0;
    double default_max_bin=1.1;
    Setup::SetHistogramDefaults(n_bins,default_min_bin,default_max_bin);
    
    // Setup User Histograms
    Setup::UserTreeAnalysis = "RhoRhoPHOTOSUserTreeAnalysis";
    
    // Description
    Setup::gen1_desc_1=" Pythia + Tauola + Photos Interface Test";
    Setup::gen1_desc_2=" $H \\rightarrow 2 \\pi^0 \\pi^+ \\pi^- \\nu_{\\tau} \\bar{\\nu_{\\tau}} $";
    Setup::gen1_desc_3=" No photon symmetrization";
    
    Setup::SuppressDecay(111); // suppress pi0 decays
  }
  else{ //Setup for analysis step
    Setup::user_analysis=MCTest01;
    //Setup::rebin_factor=4; // to reduce no of bins by rebin_factor
  }
};

