void TestReadTable() {

  gROOT->LoadMacro("AliParticleYield.cxx+");
  TClonesArray * arr = AliParticleYield::ReadFromASCIIFile("./PbPb_2760.txt");
  std::cout << "------------------------------ All Part ------------------------------" << std::endl;
  arr->Print();

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
