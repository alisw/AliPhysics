{
  Setup::decay_particle=23;
  Setup::debug_mode=false; //verbose output from MC-Tester?
  // Setup histograms 
  int n_bins=120;
  double default_min_bin=0.0; // just to try svn write privilege
  double default_max_bin=1.1;
  Setup::SetHistogramDefaults(n_bins,default_min_bin,default_max_bin);

  Setup::mass_scale_on=true;

  Setup::UserTreeAnalysis = "UserTreeAnalysis";
  Setup::UTA_params[0]=0.01/91.187; //0.010966475/.977; 
  // p_t threshold as fraction of particle energy in 
  // mothers frame 
  Setup::UTA_params[1]=2;
  Setup::UTA_params[2]=0.0;
  Setup::UTA_params[3]=1.0;
  Setup::UTA_params[4]=22;

  Setup::UTA_nparams=5;
  
  Setup::SuppressDecay(22);
  Setup::SuppressDecay(23);
  //Setup::SuppressDecay(11);
  //Setup::SuppressDecay(-11);
  //Setup::UserTreeAnalysis="UserTreeAnalysis";
};
