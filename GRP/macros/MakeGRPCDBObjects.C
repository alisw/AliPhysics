#include "GRP/MakeCosmicTriggersEntry.C"
#include "GRP/MakeLHCClockPhaseEntry.C"
#include "GRP/UpdateCDBVertexDiamond.C"
#include "GRP/MakeQAThresholdsEntry.C"
#include "GRP/MakeGRPRecoParam.C"
#include "GRP/UpdateCDBCTPConfig.C"
#include "GRP/MakeCTPEntries.C"
#include "GRP/MakeCTPDummyEntries.C"
#include "GRP/MakeCTPLTUConfigEntry.C"
#include "GRP/MakeCTPTimeAlignEntry.C"
#include "GRP/UpdateCDBIdealGeom.C"
#include "GRP/UpdateCDBGRPEntry.C"
//#include "GRP/MakeLHCDataEntry.C"

void MakeGRPCDBObjects(const char* *outputCDB = "local://$ALICE_ROOT/../AliRoot/OCDB")
{
//GRP/Align/Data
// GRP/CTP/Aliases
  MakeCosmicTriggersEntry("CosmicTriggers.txt", "local://");    // GRP/Calib/CosmicTriggers          
  MakeLHCClockPhaseEntry();     // GRP/Calib/LHCClockPhase           
  UpdateCDBVertexDiamond("MeanVertex", 0.0, 0.0, 0.0, 0.006103, 0.006103, 5.3);     // GRP/Calib/MeanVertex
  UpdateCDBVertexDiamond("MeanVertexSPD", 0.0, 0.0, 0.0, 0.006, 0.006, 3.8);     // GRP/Calib/MeanVertexSPD
  UpdateCDBVertexDiamond("MeanVertexTPC", 0.0, 0.0, 0.0, 0.006103, 0.006103, 5.3);     // GRP/Calib/MeanVertexTPC           
  MakeQAThresholdsEntry();      // GRP/Calib/QAThresholds            
  //MakeGRPRecoParam();           // GRP/Calib/RecoParam               
  //UpdateCDBCTPConfig(const char *CTPcfg, const char* cdbUri, const char* cfgFile);         // GRP/CTP/Config                    
  MakeCTPEntries();             // GRP/CTP/CTPtiming  GRP/CTP/Scalers
  MakeCTPDummyEntries();        // GRP/CTP/DummyConfig GRP/CTP/DummyCTPtime GRP/CTP/DummyLTUConfig GRP/CTP/DummyScalers              
  MakeCTPLTUConfigEntry("local://");      // GRP/CTP/LTUConfig                 
  MakeCTPTimeAlignEntry();      // GRP/CTP/TimeAlign                 
  UpdateCDBIdealGeom("local://", "$ALICE_ROOT/../AliRoot/macros/ConfigRaw2015.C");         // GRP/Geometry/Data                 
  UpdateCDBGRPEntry();          // GRP/GRP/Data                      
  //MakeLHCDataEntry();           // GRP/GRP/LHCData                   
}
