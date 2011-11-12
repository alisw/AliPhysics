// Macro to print the values stored in the OCDB with AliEMCALRecParam
// These parameters are used during reconstruction

// Author: Gustavo Conesa (INFN-LNF)


void PrintEMCALRecParam(char * file = "$ALICE_ROOT/OCDB/EMCAL/Calib/RecoParam/Run0_999999999_v0_s0.root")
{

TFile * f = new TFile(file,"READ");

AliCDBEntry * cdb = (AliCDBEntry*) f->Get("AliCDBEntry");
TObjArray * array = (TObjArray *) cdb->GetObject();

//Loop on the different event species and print parameters.
Int_t nSpecies = array->GetEntriesFast();
for(Int_t i = 0; i < nSpecies ; i++){

AliEMCALRecParam * rparam = array->At(i);

cout<<"================================================"<<endl;

cout<<"============== "<<rparam->GetName()<<" ==============="<<endl;

cout<<"================================================"<<endl;

//rparam->Print("reco");//Print only clusterizer parameters
//rparam->Print("pid");//Print only pid parameters
//rparam->Print("raw");//Print only raw digitization parameters
  rparam->Print("");// Print all


}




}
