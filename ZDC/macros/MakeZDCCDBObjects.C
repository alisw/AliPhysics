#include "ZDC/MakeZDCZeroMisAlignment.C"
#include "ZDC/MakeZDCChMapCalib.C"
#include "ZDC/MakeZDCEnergyCalib.C"
#include "ZDC/MakeZDCPedestalCalib.C"
#include "ZDC/MakeZDCRecoParam.C"
#include "ZDC/MakeZDCSaturationCalib.C"
#include "ZDC/MakeZDCTDCCalib.C"
#include "ZDC/MakeZDCTowerCalib.C"

void MakeZDCCDBObjects()
{
  // ZDC/Calib/Calib
  // ZDC/Calib/MBCalib
  MakeZDCZeroMisAlignment();    // ZDC/Align/Data                       
  MakeZDCChMapCalib();          // ZDC/Calib/ChMap                      
  MakeZDCEnergyCalib();         // ZDC/Calib/EnergyCalib                
  MakeZDCPedestalCalib();       // ZDC/Calib/Pedestals                  
  //MakeZDCRecoParam();           // ZDC/Calib/RecoParam                  
  MakeZDCSaturationCalib();     // ZDC/Calib/SaturationCalib            
  MakeZDCTDCCalib();            // ZDC/Calib/TDCCalib                   
  MakeZDCTowerCalib();          // ZDC/Calib/TowerCalib                 
}
