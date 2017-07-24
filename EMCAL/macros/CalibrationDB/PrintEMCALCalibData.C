///
/// \file PrintEMCALCalibData.C
/// \ingroup EMCAL_CalibDB
/// \brief Print energy calibration parameters in OCDB
///
/// Macro to print the values stored in the OCDB with AliEMCALCalibData, energy calibration,
/// either local file or alien file
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)
///

///
/// Main method
///
/// \param file: full file path
///
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
