#include "EMCAL/macros/CalibrationDB/AliEMCALSetCDB.C"
#include "EMCAL/macros/mapping/MakeEMCALAltroMapping.C"
#include "EMCAL/macros/PeakFinder/MakeEMCALPF.C"
#include "EMCAL/macros/PedestalDB/AliEMCALPedestalCDB.C"
#include "EMCAL/macros/RecParamDB/AliEMCALSetRecParamCDB.C"
#include "EMCAL/macros/SimParamDB/AliEMCALSetSimParamCDB.C"
#include "EMCAL/macros/Shuttle/MakeOCDBConfigPreprocessor.C"
#include "EMCAL/macros/Shuttle/MakeOCDBTempTree.C"

void MakeEMCALCDBObjects()
{
  //AliEMCALSetCDB();
  //MakeEMCALAltroMapping();
  //MakeEMCALPF();
  AliEMCALPedestalCDB();
  AliEMCALSetRecParamCDB();
  AliEMCALSetSimParamCDB();
  MakeOCDBConfigPreprocessor();
  MakeOCDBTempTree();

}
