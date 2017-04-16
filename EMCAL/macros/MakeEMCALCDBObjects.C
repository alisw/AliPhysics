////////////////////////////////////////////////////////////////////////////////
///
/// \file MakeEMCALCDBObjects.C
/// \ingroup EMCAL_macros
/// \brief Generate OCDB objects
///
///  Call the different single macros generating OCDB files with default values
///
/// \author Raffaele Grosso, Raffaele.Grosso@cern.ch (CERN)
///
////////////////////////////////////////////////////////////////////////////////

#include "EMCAL/macros/CalibrationDB/AliEMCALSetCDB.C"
#include "EMCAL/macros/mapping/MakeEMCALAltroMapping.C"
#include "EMCAL/macros/PeakFinder/MakeEMCALPF.C"
#include "EMCAL/macros/PedestalDB/AliEMCALPedestalCDB.C"
#include "EMCAL/macros/SimRecParamDB/AliEMCALSetRecParamCDB.C"
#include "EMCAL/macros/SimRecParamDB/AliEMCALSetSimParamCDB.C"
#include "EMCAL/macros/Shuttle/MakeOCDBConfigPreprocessor.C"
#include "EMCAL/macros/Shuttle/MakeOCDBTempTree.C"

///
/// Main method, call all the macros
///
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
