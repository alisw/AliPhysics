//////////////////////////////////////////////////////////////////////////////////
// Class AliHLTT0CalibObject for HLT calibration
// vector fT0CalibParams:
// 0-23 mean CFD
//                                       24-47 diff CFD
//                                       48 T0AC shift
//                                       49 T0A shift  
//                                       50 T0C shift
/////////////////////////////////////////////////////////////////////////////////////
#include "AliHLTT0CalibObject.h"

ClassImp(AliHLTT0CalibObject)

//------------------------------------
AliHLTT0CalibObject::AliHLTT0CalibObject() : TObject()
{
 
  for (int i=0; i<52; i++) fT0CalibParams[i] = 0; 
}
//------------------------------------
AliHLTT0CalibObject::AliHLTT0CalibObject(const AliHLTT0CalibObject &r) : TObject()
{
  //
  // Copy constructor
  ((AliHLTT0CalibObject &) r).Copy(*this);
}
