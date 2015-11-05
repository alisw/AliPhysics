/// \file PrintEMCALCalibTime.C
/// \brief Print time parameters in OCDB
///
/// Macro to print the values stored in the OCDB with AliEMCALCalibTime
/// either local file or alien file
///
///
/// \author : Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)
///

///
/// Main method
///
/// \param file: full file path
///
//______________________________________
void PrintEMCALCalibTime(TString file = 
			 /*"alien:///alice/data/2015/OCDB/EMCAL/Calib/Time/Run0_999999999_v1_s0.root"*/
			 "$ALICE_ROOT/OCDB/EMCAL/Calib/Time/Run0_999999999_v0_s0.root"
/*"/AliceSoft/aliroot/master/src/EMCAL/macros/CalibrationDB/DeCalibGausDB/EMCAL/Calib/Time/Run0_999999999_v0_s0.root"*/
)
{  
  if(file.Contains("alien:///"))
    TGrid::Connect("alien://");

  TFile * f = TFile::Open(file,"READ");
  
  AliCDBEntry * cdb = (AliCDBEntry*) f->Get("AliCDBEntry");    
  AliEMCALCalibTime * cparam =  cdb->GetObject();
    
  cparam->Print("");
}
