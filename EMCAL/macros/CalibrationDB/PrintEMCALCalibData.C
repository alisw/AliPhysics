// Macro to print the values stored in the OCDB with AliEMCALCalibData
// either local file or alien file

// Author: Gustavo Conesa (LPSC-IN2P3)


void PrintEMCALCalibData(TString file = 
			 /*"alien:///alice/data/2014/OCDB/EMCAL/Calib/Data/Run0_999999999_v1_s0.root"*/
			 "$ALICE_ROOT/OCDB/EMCAL/Calib/Data/Run0_999999999_v0_s0.root"){
  
  if(file.Contains("alien:///"))
    TGrid::Connect("alien://");

  TFile * f = TFile::Open(file,"READ");
  
  AliCDBEntry * cdb = (AliCDBEntry*) f->Get("AliCDBEntry");    
  AliEMCALCalibData * cparam =  cdb->GetObject();
    
  //cparam->Print("gain");
  //cparam->Print("ped");  
  cparam->Print("all");

}
