void TestReadTable() {

  LoadLibs();
  TClonesArray * arr = AliParticleYield::ReadFromASCIIFile("./PbPb_2760.txt");
  std::cout << "------------------------------ All Part ------------------------------" << std::endl;
  arr->Print();
  AliParticleYield::WriteThermusFile(arr, "thermus.txt" );

  // Get it as tree
  TTree * tree = AliParticleYield::ReadFromASCIIFileAsTree("./PbPb_2760.txt");

  // examples on how to extract sub arrays;
  delete arr;
  arr =AliParticleYield::GetEntriesMatchingSelection(tree, "fCentr == \"V0M0010\" && fStatus == 0");
  std::cout << "------------------------------ CENTR = 0-10%, Status = 0 ------------------------------" << std::endl;
  arr->Print();

  delete arr;
  arr =AliParticleYield::GetEntriesMatchingSelection(tree, "fCentr == \"V0M0020\" && !IsTypeRatio()");
  std::cout << "------------------------------ CENTR = 0-20%, no ratios ------------------------------" << std::endl;
  arr->Print();


						    

  AliParticleYield::SaveAsASCIIFile(arr,"pippo.txt");

}

void LoadLibs() {

  gSystem->Load("libCore.so");  
  gSystem->Load("libGeom.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libVMC");
  gSystem->Load("libTree");
  gSystem->Load("libProof");
  gSystem->Load("libMatrix");
  gSystem->Load("libMinuit");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libOADB");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libTENDER");
  gSystem->Load("libCORRFW");
  gSystem->Load("libMinuit");
  gSystem->Load("libPWGTools");
  //  gROOT->LoadMacro("AliParticleYield.cxx+");

  gSystem->Load("libPWGLFthermalfits");
  

}
