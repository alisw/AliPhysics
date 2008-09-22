
void testFastJet(const char* file="testdata.dat")
{
      gSystem->Load("libCGAL.so");
      gSystem->Load("libfastjet.so");
     
      gSystem->Load("libANALYSIS.so");
      gSystem->Load("libSTEERBase.so");
      gSystem->Load("libAOD.so");
      gSystem->Load("libESD.so");
      gSystem->Load("libANALYSISalice.so");
      gSystem->Load("libJETAN.so");
      
      AliFastJetFinder* jetFinder = new AliFastJetFinder();
      jetFinder->RunTest(file);
      
      cout << "bye bye " << endl;
}