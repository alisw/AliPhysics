/// \file PrintEMCALSimParam.C
/// \brief Print OCDB simulation parameters
///
/// Macro to print the values stored in the OCDB with AliEMCALSimParam
/// These parameters are used during simulation
///
/// \author : Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)
///

///
/// Main method, access local or grid file, print its content.
///
/// \param file: string with location of input file
///
void PrintEMCALSimParam(TString file = 
                        "alien:///alice/data/2015/OCDB/EMCAL/Calib/SimParam/Run0_999999999_v3_s0.root"
//"$ALICE_ROOT/OCDB/EMCAL/Calib/SimParam/Run0_999999999_v0_s0.root"
                        )
{
  if(file.Contains("alien://"))
    TGrid::Connect("alien://");
  
  TFile * f = TFile::Open(file);
  
  AliCDBEntry * cdb = (AliCDBEntry*) f->Get("AliCDBEntry");    
  AliEMCALSimParam * sparam =  cdb->GetObject();
  
  
  cout<<"============== "<<sparam->GetName()<<" ==============="<<endl;
  
  sparam->Print("");
}
