#include "TRD/MakeTRDZeroMisAlignment.C"
#include "TRD/TRDbase/AliTRDCreateDummyCDB.C"
#include "TRD/TRDbase/AliTRDCreateDummyCDB_DP.C"
#include "TRD/TRDbase/AliTRDCreateDummyCDB_DCS.C"
#include "TRD/TRDbase/AliTRDCreateLocalGain.C"
#include "TRD/TRDbase/AliTRDCreateLocalGain.C"
#include "TRD/TRDbase/AliTRDCreateOCDBPIDLQ.C"
#include "TRD/TRDbase/AliTRDmakePIDThresholds.C"
//#include "TRD/TRDbase/AliTRDmakeRecoParam.C"
#include "TRD/Macros/AliTRDcreateTrapConfigCDB.C"
#include "TRD/TRDbase/AliTRDmakeTrkDB.C"

void MakeTRDCDBObjects()
{
  MakeTRDZeroMisAlignment(); //TRD/Align/Data
  AliTRDCreateDummyCDB();    //ChamberExB ChamberGainFactor ChamberStatus ChamberT0 ChamberVdrift DetNoise LocalGainFactor LocalT0 LocalVdrift MonitoringData PadNoise PadStatus PRFWidth
  // TRD/Calib/RecoParam ... could be created also with AliTRDmakeRecoParam.C
  AliTRDCreateDummyCDB_DP(); // trd_chamberStatus trd_envTemp trd_gasCO2 trd_gasH2O trd_gasO2 trd_gasOverpressure trd_goofieCO2 trd_goofieGain trd_goofieHv trd_goofieN2 trd_goofiePeakArea trd_goofiePeakPos trd_goofiePressure trd_goofieTemp trd_goofieVelocity trd_hvAnodeImon trd_hvAnodeUmon trd_hvDriftImon trd_hvDriftUmon trd_gaschromatographXe trd_gaschromatographCO2 trd_gaschromatographN2
  AliTRDCreateDummyCDB_DCS(); //TRD/Calib/DCS 
  AliTRDCreateLocalGain();    // Gaintbl_Uniform_FGAN0_2011-01 Gaintbl_Uniform_FGAN0_2012-01 Gaintbl_Uniform_FGAN8_2011-01 Gaintbl_Uniform_FGAN8_2012-01 Krypton_2011-01 Krypton_2011-02 Krypton_2011-03 Krypton_2012-01 Krypton_2015-01
  AliTRDCreateOCDBPIDLQ();      // TRD/Calib/PIDLQ1D
  AliTRDmakePIDThresholds();    // TRD/Calib/PIDThresholds
  AliTRDcreateTrapConfigCDB();   // TRD/Calib/TrapConfig
  AliTRDmakeTrkDB();            // TRD/Calib/TrkAttach
}

/*  Macro not yet defined for:
TRD/Calib/PHQ
TRD/Calib/PIDLQ
TRD/Calib/PIDNN
*/
