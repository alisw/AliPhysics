
void testFastJet(const char* file="testdata.dat")

{
      gSystem->Load("libTree.so");
      gSystem->Load("libEG.so");
      gSystem->Load("libVMC.so");
      gSystem->Load("libPhysics.so");

      gSystem->Load("libCGAL.dylib");

      gSystem->Load("libfastjet.so");
      gSystem->Load("libSISConePlugin.so");
      
     
      gSystem->Load("libANALYSIS.so");
      gSystem->Load("libSTEERBase.so");
      gSystem->Load("libAOD.so");
      gSystem->Load("libESD.so");
      gSystem->Load("libANALYSISalice.so");
      
      gSystem->Load("libJETAN.so");
      gSystem->Load("libFASTJETAN.so");

      AliFastJetFinder* jetFinder = new AliFastJetFinder();
      jetFinder->RunTest(file);
      
      cout << "bye bye " << endl;
}
