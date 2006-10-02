#include "AliITSEventHeader.h"
////////////////////////////////////////////////////
//  Base class to define                          //
//  ITS Event Header                              //
//                                                //
//  Origin: E. Crescio crescio@to.infn.it         //
//                                                //
////////////////////////////////////////////////////

ClassImp(AliITSEventHeader)

  
//_____________________________________________________________
  AliITSEventHeader::AliITSEventHeader():AliDetectorEventHeader(),
fEventTypeSDD(),
fJitterSDD(0)
{
  //
  // Defaulst Constructor
  //
  
  SetEventTypeSDD();
  
  for(Int_t idet=0;idet<3;idet++){
    fL1TriggerType[idet]=0;
    fOrbitNumber[idet]=0;
    fBunchCross[idet]=0;
    fBlockAttr[idet]=0;
    fTriggerClass[idet]=0;
    fStatusBits[idet]=0;
    fMiniEvId[idet]=0;
    fSubDet[idet]=0;
    fVersion[idet]=0;
  }
  
  
}

//_____________________________________________________________
AliITSEventHeader::AliITSEventHeader(const char* name):AliDetectorEventHeader(name),
fEventTypeSDD(),
fJitterSDD(-123)
{
  //
  // Constructor
  //
  SetEventTypeSDD();
  for(Int_t idet=0;idet<3;idet++){
    fL1TriggerType[idet]=0;
    fOrbitNumber[idet]=0;
    fBunchCross[idet]=0;
    fBlockAttr[idet]=0;
    fTriggerClass[idet]=0;
    fStatusBits[idet]=0;
    fMiniEvId[idet]=0;
    fSubDet[idet]=0;
    fVersion[idet]=0;
  }
  

}





