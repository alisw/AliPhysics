//
// This macro is a part of "Alice PPR fast flow analysis package"
// 
// The macro creates new pico DST structure in the file 
// "flowPiocEvent.root" and create structure "flowData" inside
// 
// Once Ntuple is created it can be filled with data using  
// macro AliFlowReconstruction.C
//
// Sylwester Radomski, GSI
// email: S.Radomski@gsi.de
// 22. Oct. 2002
//

AliFlowCreatePico() {

  const char* fileName = "flowPicoEvent.root";
  const char* structure = "runNumber:evNumber:Mult:truePsi:trueV2:Psi:PsiA:PsiB:V2";
  
  TFile *file = new TFile(fileName, "recreate");
  TNtuple *data = new TNtuple("flowData","Flow Data", structure);
  
  data->Print();

  data->Write();
  file->Close();

  gSystem->Exit(0);
}
