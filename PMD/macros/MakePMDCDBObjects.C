#include "PMD/ocdbmacro/MakePMDZeroMisAlignment.C"
#include "PMD/ocdbmacro/MakePMDDDLinfoCDB.C"
#include "PMD/ocdbmacro/MakePMDGainCDB.C"
#include "PMD/ocdbmacro/MakePMDHotCDB.C"
#include "PMD/ocdbmacro/MakePMDMappingCDB.C"
#include "PMD/ocdbmacro/MakePMDNoiseCutCDB.C"
#include "PMD/ocdbmacro/MakePMDPedCDB.C"
#include "PMD/ocdbmacro/MakePMDRecoParam.C"

void MakePMDCDBObjects()
{
  MakePMDZeroMisAlignment();  // PMD/Align/Data                           
  MakePMDDDLinfoCDB();        // PMD/Calib/Ddlinfo                        
  MakePMDGainCDB();           // PMD/Calib/Gain                           
  MakePMDHotCDB();            // PMD/Calib/Hot                            
  MakePMDMappingCDB();        // PMD/Calib/Mapping                        
  MakePMDNoiseCutCDB();       // PMD/Calib/NoiseCut                       
  MakePMDPedCDB();            // PMD/Calib/Ped                            
  MakePMDRecoParam();         // PMD/Calib/RecoParam                      
}

