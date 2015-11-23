/// \file PrintOCDBObject.C
/// \brief Print OCDB alignment matrices
///
/// Macro to print the aligment matrices stored in the OCDB 
/// These parameters are used during reconstruction and simulation
///
/// \author : Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)
///

///
/// Main method, access local or grid file, print its content.
///
/// \param file: string with location of input file
///
void PrintOCDBObject(TString file = 
                     "$ALICE_ROOT/OCDB/EMCAL/Align/Data/Run0_999999999_v0_s0.root"
                     //"alien:///alice/data/2015/OCDB/EMCAL/Align/Data/Run0_999999999_v1_s0.root"
                     )
{
  if(file.Contains("alien://"))
    TGrid::Connect("alien://");

  TFile *f = TFile::Open(file);

  AliCDBEntry *entry = (AliCDBEntry*)f->Get("AliCDBEntry");
  
  TClonesArray *aligndata = dynamic_cast<TClonesArray*>  (entry->GetObject());
  
  aligndata->Print("");
}
