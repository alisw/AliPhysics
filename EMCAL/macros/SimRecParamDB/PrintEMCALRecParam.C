/// \file PrintEMCALRecParam.C
/// \brief Print OCDB reconstruction parameters
///
/// Macro to print the values stored in the OCDB with AliEMCALRecParam
/// These parameters are used during reconstruction
///
/// \author : Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)
///

///
/// Main method, access local or grid file, print its content.
///
/// \param file: string with location of input file
///
void PrintEMCALRecParam(TString file = 
                        /*"alien:///alice/data/2015/OCDB/EMCAL/Calib/RecoParam/Run0_999999999_v1_s0.root"*/
                       "$ALICE_ROOT/OCDB/EMCAL/Calib/RecoParam/Run0_999999999_v0_s0.root"
)
{
  
  if(file.Contains("alien://"))
    TGrid::Connect("alien://");
  
  TFile * f = TFile::Open(file);
  
  AliCDBEntry * cdb = (AliCDBEntry*) f->Get("AliCDBEntry");
  TObjArray * array = (TObjArray *) cdb->GetObject();
  
  //Loop on the different event species and print parameters.
  Int_t nSpecies = array->GetEntriesFast();
  for(Int_t i = 0; i < nSpecies ; i++){
    
    AliEMCALRecParam * rparam = array->At(i);
    
    cout<<"================================================"<<endl;
    
    cout<<"============== "<<rparam->GetName()<<" ==============="<<endl;
    
    cout<<"================================================"<<endl;
    
    rparam->Print("reco");//Print only clusterizer parameters
    //rparam->Print("pid");//Print only pid parameters
    rparam->Print("raw");//Print only raw digitization parameters
                         //rparam->Print("");// Print all
    
  }
}
