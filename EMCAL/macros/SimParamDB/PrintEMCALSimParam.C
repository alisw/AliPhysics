// Macro to print the values stored in the OCDB with AliEMCALRecParam
// These parameters are used during reconstruction

// Author: Gustavo Conesa (INFN-LNF)


void PrintEMCALSimParam(char * file = "$ALICE_ROOT/OCDB/EMCAL/Calib/SimParam/Run0_999999999_v0_s0.root"){
  
  
  TFile * f = new TFile(file,"READ");
  
  AliCDBEntry * cdb = (AliCDBEntry*) f->Get("AliCDBEntry");    
  AliEMCALSimParam * sparam =  cdb->GetObject();
  
  
  cout<<"============== "<<sparam->GetName()<<" ==============="<<endl;
  
  sparam->Print("");
  
}
