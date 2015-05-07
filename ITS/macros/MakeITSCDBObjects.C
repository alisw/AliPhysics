#include "ITS/defaultObjSPD.C"
#include "ITS/macrosSDD/StoreCalibSDD.C"
#include "ITS/macrosSDD/StoreDriftSpeedSDD.C"
#include "ITS/macrosSDD/StoreRespSDD.C"
#include "ITS/macrosSDD/StoreDDLMapSDD.C"
#include "ITS/macrosSDD/StoreMapsSDD.C"
#include "ITS/MakeCalibrationSSD.C"

void MakeITSCDBObjects()
{
  defaultObjSPD(0, AliCDBRunRange::Infinity(), "local://$ALICE_ROOT/OCDB");      // all SPD objects
  StoreCalibSDD();        // ITS/Calib/CalibSDD
  StoreDriftSpeedSDD();   // ITS/Calib/DriftSpeedSDD
  StoreRespSDD();         // ITS/Calib/RespSDD
  StoreDDLMapSDD();       // ITS/Calib/DDLMapSDD
  StoreMapsSDD();         // ITS/Calib/MapsTimeSDD
  MakeCalibrationSSD(); // all SSD objects
}

