//
//  simple test macro for FASTJET in JETAN
//
//  This macro only does a technical test of FASTJETAN. This is
//    NOT a good starting point for developing analysis code.

//  More complete examples can be found at:Somewhat outdated; for a 
//   more recent set of macros, look at:
//    $ALICE_ROOT/PWGJE//EMCALJetTasks/macros/runEMCalJetAnalysis.C
//    $ALICE_ROOT/PWGJE/macros/examples
//  and documentation at: https://twiki.cern.ch/twiki/bin/view/ALICE/EMCalJet
//
void testFastJet(const char* file="testdata.dat")

{
      gSystem->Load("libTree");
      gSystem->Load("libEG");
      gSystem->Load("libVMC");
      gSystem->Load("libPhysics");

      gSystem->Load("libCGAL");

      gSystem->Load("libfastjet");
      gSystem->Load("libsiscone");
      gSystem->Load("libsiscone_spherical");
      gSystem->Load("libfastjetplugins");
     
      gSystem->Load("libANALYSIS");
      gSystem->Load("libSTEERBase");
      gSystem->Load("libAOD");
      gSystem->Load("libESD");
      gSystem->Load("libANALYSISalice");
      
      gSystem->Load("libJETAN");
      gSystem->Load("libFASTJETAN");

      AliFastJetFinder* jetFinder = new AliFastJetFinder();
      jetFinder->RunTest(file);
      
      cout << "bye bye " << endl;
}
