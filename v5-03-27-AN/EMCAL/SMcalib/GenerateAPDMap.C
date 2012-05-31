void GenerateAPDMap()
{
  gSystem->Load("AliEMCALMapAPD_cxx");
  AliEMCALMapAPD *mapAPD = new AliEMCALMapAPD();

  Int_t iSM[2] = {0,1}; // allow for two SuperModules
  mapAPD->GenerateDummyAPDInfo(1, iSM); // space for one SuperModule, with number iSM[0] = 0

  mapAPD->ReadMapAPDInfoSingleStripBasis(0, 0, "APD/APDStripModW15.txt"); 
  mapAPD->ReadMapAPDInfoSingleStripBasis(0, 1, "APD/APDStripModW21.txt"); 
  mapAPD->ReadMapAPDInfoSingleStripBasis(0, 2, "APD/APDStripModW3.txt"); 
  mapAPD->ReadMapAPDInfoSingleStripBasis(0, 3, "APD/APDStripModW5.txt"); 
  mapAPD->ReadMapAPDInfoSingleStripBasis(0, 4, "APD/APDStripModW18.txt"); 
  mapAPD->ReadMapAPDInfoSingleStripBasis(0, 5, "APD/APDStripModW8.txt"); 
  mapAPD->ReadMapAPDInfoSingleStripBasis(0, 6, "APD/APDStripModW16.txt"); 
  mapAPD->ReadMapAPDInfoSingleStripBasis(0, 7, "APD/APDStripModW17.txt");
  mapAPD->ReadMapAPDInfoSingleStripBasis(0, 8, "APD/APDStripModW14.txt");  

  // not yet installed: use dummy values for now..
  mapAPD->ReadMapAPDInfoSingleStripBasis(0, 9, "APD/APDStripModW0.txt"); 
  mapAPD->ReadMapAPDInfoSingleStripBasis(0,10, "APD/APDStripModW7.txt"); 
  mapAPD->ReadMapAPDInfoSingleStripBasis(0,11, "APD/APDStripModW10.txt"); 
  mapAPD->ReadMapAPDInfoSingleStripBasis(0,12, "APD/APDStripModW11.txt"); 
  mapAPD->ReadMapAPDInfoSingleStripBasis(0,13, "APD/APDStripModW12.txt"); 
  mapAPD->ReadMapAPDInfoSingleStripBasis(0,14, "APD/APDStripModW13.txt"); 
  mapAPD->ReadMapAPDInfoSingleStripBasis(0,15, "APD/APDStripModW19.txt"); 
  mapAPD->ReadMapAPDInfoSingleStripBasis(0,16, "APD/APDStripModW1.txt"); 
  mapAPD->ReadMapAPDInfoSingleStripBasis(0,17, "APD/APDStripModW2.txt"); 
  mapAPD->ReadMapAPDInfoSingleStripBasis(0,18, "APD/APDStripModW4.txt"); 
  mapAPD->ReadMapAPDInfoSingleStripBasis(0,19, "APD/APDStripModW6.txt"); 

  mapAPD->WriteMapAPDInfo("APD/APDSuperModW1.txt");

  int nProblems = mapAPD->CheckForDuplicates();
}
