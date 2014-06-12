{

  if (!Setup::stage == 0) { //generation step configurables

    Setup::decay_particle=23;
    Setup::debug_mode=false;

    // Setup histograms 
    int n_bins=120;
    double default_min_bin=0.0;
    double default_max_bin=1.1;
    Setup::SetHistogramDefaults(n_bins,default_min_bin,default_max_bin);
    Setup::mass_scale_on=true;
    
  
    // Description
    Setup::gen1_desc_1=" Pythia + Photos Interface Test";
    Setup::gen1_desc_2=" $Z \\rightarrow \\mu^+ \\mu^-$. Photons filtered below 10 MeV";
    Setup::gen1_desc_3=" No photon symmetrization";
    
    //Filter photons
    Setup::UserTreeAnalysis = "UserTreeAnalysis";
    Setup::UTA_params[0]=0.01/91.187; //10 MeV
    // p_t threshold as fraction of particle energy in 
    // mothers frame 
    Setup::UTA_params[1]=2;
    Setup::UTA_params[2]=0.0;
    Setup::UTA_params[3]=1.0;
    Setup::UTA_params[4]=22;
    
    Setup::UTA_nparams=5;
    
    Setup::SuppressDecay(22);
    Setup::SuppressDecay(23);
     /**************************************************************************
                          Settings for old FORTRAN tests
               Uncomment when generating comparison with these files
    ***************************************************************************/
/*
    n_bins=1200;
    default_min_bin=0.0;
    default_max_bin=120.0;
    Setup::SetHistogramDefaults(n_bins,default_min_bin,default_max_bin);
    Setup::mass_scale_on=false;
    Setup::mass_power=1;

    Setup::UserTreeAnalysis = "UserTreeAnalysis";
    Setup::UTA_params[0]=1./91.187;
    Setup::UTA_params[1]=1;
*/
    /**************************************************************************/
   
  }
  else{ //Setup for analysis step
    Setup::user_analysis=MCTest01;
    //Setup::rebin_factor=4; // to reduce no of bins by rebin_factor
    Setup::use_log_y=true;
  }
};
