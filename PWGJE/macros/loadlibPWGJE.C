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

  gSystem->Load("libPWGJE");
  gSystem->Load("libPWGJEStrangenessInJets");

  // EMCAL jet framework
  gSystem->Load("libTENDER.so");
  gSystem->Load("libCORRFW.so");
  gSystem->Load("libPWGTools.so");
  gSystem->Load("libCDB.so");
  gSystem->Load("libProof.so");
  gSystem->Load("libfastjetplugins.so");
  gSystem->Load("libfastjettools.so");
  gSystem->Load("libfastjetcontribfragile.so");
  gSystem->Load("libRAWDatabase.so");
  gSystem->Load("libSTEER.so");
  gSystem->Load("libPWGEMCAL.so");
  gSystem->Load("libPWGJEEMCALJetTasks.so");
}
