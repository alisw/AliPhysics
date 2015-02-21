void loadlibPWGJE()
{
  // this macro should load PWGJE with all dependencies from ROOT
  // (it should not require to be run with AliRoot)
  //
  // please use the macro before a commit to check
  // that the PWGJE libraries not only compile but can also be loaded
  // with the known dependencies
  //
  // if you have to add dependendencies make sure that people get informed
  // (especially the operators of trains should know what to change)

  gSystem->Load("libCore");
  gSystem->Load("libPhysics");
  gSystem->Load("libMinuit");

  gSystem->Load("libVMC");

  gSystem->Load("libNet");
  gSystem->Load("libTree");

  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");

  gSystem->Load("libCGAL");
  gSystem->Load("libfastjet");

  gSystem->Load("libCORRFW");
  gSystem->Load("libPWGTools");
  gSystem->Load("libOADB");
  gSystem->Load("libPWGEMCAL");

  gSystem->Load("libCGAL");
  gSystem->Load("libsiscone");
  gSystem->Load("libsiscone_spherical");
  gSystem->Load("libfastjetplugins");
  gSystem->Load("libfastjettools");
  gSystem->Load("libfastjetcontribfragile");

  gSystem->Load("libPWGJE");
  gSystem->Load("libPWGJEStrangenessInJets");

  //
  // EMCAL jet framework
  //  if you only need the EMCAL jet framework, use $ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/loadlibPWGJEEMCAL.C
  gSystem->Load("libTender");

 
  gSystem->Load("libPWGJEEMCALJetTasks");
}
