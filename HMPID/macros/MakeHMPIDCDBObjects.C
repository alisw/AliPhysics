#include "HMPID/MakeHMPIDAlignmentObjs.C"
#include "HMPID/Hshuttle.C"
#include "HMPID/Hshuttle.C"
#include "HMPID/Hshuttle.C"
#include "HMPID/Hshuttle.C"
#include "HMPID/MakeHMPIDQeEffMaps.C"
#include "HMPID/Hshuttle.C"
#include "HMPID/MakeHMPIDRecoParamV1.C"

void MakeHMPIDCDBObjects()
{
  MakeHMPIDAlignmentObjs(kTRUE, "ideal");   //  HMPID/Align/Data                                    
  Hshuttle();                // HMPID/Calib/DaqSig HMPID/Calib/Masked HMPID/Calib/Nmean HMPID/Calib/NoiseMap HMPID/Calib/Qthre
  MakeHMPIDQeEffMaps();      // HMPID/Calib/QeMap                                  
  MakeHMPIDRecoParamV1();    // HMPID/Calib/RecoParam                              
}
