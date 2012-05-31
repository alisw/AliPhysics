// Macro to print the values stored in the OCDB with AliEMCALSimParam
// These parameters are used during simulation

// Author: Gustavo Conesa (IN2P3-LPSC)


void PrintEMCALSimParam(char * file = "$ALICE_ROOT/OCDB/EMCAL/Calib/SimParam/Run0_999999999_v0_s0.root"){
  
  
  TFile * f = new TFile(file,"READ");
  
  AliCDBEntry * cdb = (AliCDBEntry*) f->Get("AliCDBEntry");    
  AliEMCALSimParam * sparam =  cdb->GetObject();
  
  
  cout<<"============== "<<sparam->GetName()<<" ==============="<<endl;
  
  sparam->Print("");
  
}
