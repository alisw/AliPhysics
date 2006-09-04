/**************************************************************************
 * Author: Panos Christakoglou.                                           *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

//-----------------------------------------------------------------
//                   AliRunTagCuts class
//   This is the class to deal with the run tag level cuts
//   Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-----------------------------------------------------------------

class AliLog;
class AliESD;

#include "AliRunTag.h"
#include "AliRunTagCuts.h"

ClassImp(AliRunTagCuts)


//___________________________________________________________________________
AliRunTagCuts::AliRunTagCuts() :
  TObject(),
  fAliceRunId(-1),                  
  fAliceRunIdFlag(kFALSE),              
  fAliceMagneticField(-1.),          
  fAliceMagneticFieldFlag(kFALSE),      
  fAliceRunStartTimeMin(-1),        
  fAliceRunStartTimeMax(-1),        
  fAliceRunStartTimeFlag(kFALSE),       
  fAliceRunStopTimeMin(-1),         
  fAliceRunStopTimeMax(-1),         
  fAliceRunStopTimeFlag(kFALSE),        
  fAlirootVersion(""),              
  fAlirootVersionFlag(kFALSE),          
  fRootVersion(""),                 
  fRootVersionFlag(kFALSE),             
  fGeant3Version(""),               
  fGeant3VersionFlag(kFALSE),           
  fAliceRunQuality(0),             
  fAliceRunQualityFlag(kFALSE),         
  fAliceBeamEnergy(-1),             
  fAliceBeamEnergyFlag(kFALSE),         
  fAliceBeamType(""),               
  fAliceBeamTypeFlag(kFALSE),           
  fAliceCalibrationVersion(-1),    
  fAliceCalibrationVersionFlag(kFALSE),
  fAliceDataType(-1),                
  fAliceDataTypeFlag(kFALSE)
{
  //Default constructor which calls the Reset method.
}

//___________________________________________________________________________
AliRunTagCuts::~AliRunTagCuts() {  
  //Defaut destructor.
}

//___________________________________________________________________________
void AliRunTagCuts::Reset() {
  //Sets dummy values to every private member.
  fAliceRunId = -1;                  
  fAliceRunIdFlag = kFALSE;              
  fAliceMagneticField = -1.;          
  fAliceMagneticFieldFlag = kFALSE;      
  fAliceRunStartTimeMin = -1;        
  fAliceRunStartTimeMax = -1;        
  fAliceRunStartTimeFlag = kFALSE;       
  fAliceRunStopTimeMin = -1;         
  fAliceRunStopTimeMax = -1;         
  fAliceRunStopTimeFlag = kFALSE;        
  fAlirootVersion = "";              
  fAlirootVersionFlag = kFALSE;          
  fRootVersion = "";                 
  fRootVersionFlag = kFALSE;             
  fGeant3Version = "";               
  fGeant3VersionFlag = kFALSE;           
  fAliceRunQuality = 0;             
  fAliceRunQualityFlag = kFALSE;         
  fAliceBeamEnergy = -1;             
  fAliceBeamEnergyFlag = kFALSE;         
  fAliceBeamType = "";               
  fAliceBeamTypeFlag = kFALSE;           
  fAliceCalibrationVersion = -1;    
  fAliceCalibrationVersionFlag = kFALSE;
  fAliceDataType = -1;                
  fAliceDataTypeFlag = kFALSE;           
}

//___________________________________________________________________________
Bool_t AliRunTagCuts::IsAccepted(AliRunTag *RunTag) const {
  //Returns true if the event is accepted otherwise false.
  if(fAliceRunIdFlag)
    if((RunTag->GetRunId() != fAliceRunId))
      return kFALSE;
  if(fAliceMagneticFieldFlag)
    if((RunTag->GetMagneticField() != fAliceMagneticField))
      return kFALSE;
  if(fAliceRunStartTimeFlag)
    if((RunTag->GetRunStartTime() < fAliceRunStartTimeMin) || (RunTag->GetRunStartTime() > fAliceRunStartTimeMax))
      return kFALSE;
  if(fAliceRunStopTimeFlag)
    if((RunTag->GetRunStopTime() < fAliceRunStopTimeMin) || (RunTag->GetRunStopTime() > fAliceRunStopTimeMax))
      return kFALSE;
  if(fAlirootVersionFlag)
    if((RunTag->GetAlirootVersion() != fAlirootVersion))
      return kFALSE;
  if(fRootVersionFlag)
    if((RunTag->GetRootVersion() != fRootVersion))
      return kFALSE;
  if(fGeant3VersionFlag)
    if((RunTag->GetGeant3Version() != fGeant3Version))
      return kFALSE;
  if(fAliceRunQualityFlag)
    if(RunTag->GetRunQuality())
      return kFALSE;
  if(fAliceBeamEnergyFlag)
    if(RunTag->GetBeamEnergy() != fAliceBeamEnergy)
      return kFALSE;
  if(fAliceBeamTypeFlag)
    if(RunTag->GetBeamType() != fAliceBeamType)
      return kFALSE;
  if(fAliceCalibrationVersionFlag)
    if(RunTag->GetBeamEnergy() != fAliceBeamEnergy)
      return kFALSE;
  if(fAliceDataTypeFlag)
    if(RunTag->GetDataType() != fAliceDataType)
      return kFALSE;
 
  return kTRUE;
}
