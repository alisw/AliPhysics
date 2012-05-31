// Macro to print the values stored in the OCDB with AliEMCALCalibData
// These parameters are used during simulation and reconstruction

// Author: Gustavo Conesa (LPSC-IN2P3)


void PrintEMCALCalibData(char * file = "$ALICE_ROOT/OCDB/EMCAL/Calib/Data/Run0_999999999_v0_s0.root"){
  
  
  TFile * f = new TFile(file,"READ");
  
  AliCDBEntry * cdb = (AliCDBEntry*) f->Get("AliCDBEntry");    
  AliEMCALCalibData * cparam =  cdb->GetObject();
  
  
  cout<<"============== "<<cparam->GetName()<<" ==============="<<endl;
  
  cparam->Print("gain");
  //cparam->Print("ped");  
}
