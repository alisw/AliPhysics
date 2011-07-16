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
  fAliceDipoleField(-1.),
  fAliceDipoleFieldFlag(kFALSE),
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
  fLHCPeriod(""),
  fLHCPeriodFlag(kFALSE),
  fRecPass(""),
  fRecPassFlag(kFALSE),
  fProdName(""),
  fProdNameFlag(kFALSE),
  fAliceRunValidation(0),             
  fAliceRunValidationFlag(kFALSE),         
  fAliceRunQualities(""),
  fAliceRunQualitiesFlag(kFALSE),
  fAliceBeamEnergy(-1),             
  fAliceBeamEnergyFlag(kFALSE),         
  fAliceBeamType(""),               
  fAliceBeamTypeFlag(kFALSE),           
  fAliceCalibrationVersion(-1),    
  fAliceCalibrationVersionFlag(kFALSE),
  fAliceDataType(-1),                
  fAliceDataTypeFlag(kFALSE),
  fBeamTriggerRange(),
  fBeamTriggerFlag(kFALSE),
  fCollisionTriggerRange(),
  fCollisionTriggerFlag(kFALSE),
  fEmptyTriggerRange(),
  fEmptyTriggerFlag(kFALSE),
  fASideTriggerRange(),
  fASideTriggerFlag(kFALSE),
  fCSideTriggerRange(),
  fCSideTriggerFlag(kFALSE),
  fHMTriggerRange(),
  fHMTriggerFlag(kFALSE),
  fMuonTriggerRange(),
  fMuonTriggerFlag(kFALSE),
  fCollisionRateRange(),
  fCollisionRateFlag(kFALSE),
  fMeanVertexRange(),
  fMeanVertexFlag(kFALSE),
  fVertexQualityRange(),
  fVertexQualityFlag(kFALSE)
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
  fAliceDipoleField = -1.;
  fAliceDipoleFieldFlag = kFALSE;
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
  fLHCPeriod = "";
  fLHCPeriodFlag = kFALSE;
  fRecPass = "";
  fRecPassFlag = kFALSE;
  fProdName = "";
  fProdNameFlag = kFALSE;
  fAliceRunValidation = 0;             
  fAliceRunValidationFlag = kFALSE;         
  fAliceRunQualities = "";
  fAliceRunQualitiesFlag = kFALSE;
  fAliceBeamEnergy = -1;             
  fAliceBeamEnergyFlag = kFALSE;         
  fAliceBeamType = "";               
  fAliceBeamTypeFlag = kFALSE;           
  fAliceCalibrationVersion = -1;    
  fAliceCalibrationVersionFlag = kFALSE;
  fAliceDataType = -1;                
  fAliceDataTypeFlag = kFALSE;           
  fBeamTriggerRange[0] = 0;
  fBeamTriggerRange[1] = 0;
  fBeamTriggerFlag = kFALSE;
  fCollisionTriggerRange[0] = 0;
  fCollisionTriggerRange[1] = 0;
  fCollisionTriggerFlag = kFALSE;
  fEmptyTriggerRange[0] = 0;
  fEmptyTriggerRange[1] = 0;
  fEmptyTriggerFlag = kFALSE;
  fASideTriggerRange[0] = 0;
  fASideTriggerRange[1] = 0;
  fASideTriggerFlag = kFALSE;
  fCSideTriggerRange[0] = 0;
  fCSideTriggerRange[1] = 0;
  fCSideTriggerFlag = kFALSE;
  fHMTriggerRange[0] = 0;
  fHMTriggerRange[1] = 0;
  fHMTriggerFlag = kFALSE;
  fMuonTriggerRange[0] = 0;
  fMuonTriggerRange[1] = 0;
  fMuonTriggerFlag = kFALSE;
  fCollisionRateRange[0] = 0;
  fCollisionRateRange[1] = 0;
  fCollisionRateFlag = kFALSE;
  fMeanVertexRange[0] = 0;
  fMeanVertexRange[1] = 0;
  fMeanVertexFlag = kFALSE;
  fVertexQualityRange[0] = 0;
  fVertexQualityRange[1] = 0;
  fVertexQualityFlag = kFALSE;
}

void AliRunTagCuts::AddRunQualityValue(Int_t qval)
{
  // Adds to the list of selected run qualities
  fAliceRunQualities += qval;
  fAliceRunQualities += " ";
  fAliceRunQualitiesFlag = kTRUE;
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
  if(fAliceDipoleFieldFlag)
    if((RunTag->GetDipoleField() != fAliceDipoleField))
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
  if(fLHCPeriodFlag)
    if((RunTag->GetLHCPeriod() != fLHCPeriod))
      return kFALSE;
  if(fRecPassFlag)
    if((RunTag->GetReconstructionPass() != fRecPass))
      return kFALSE;
  if(fProdNameFlag)
    if((RunTag->GetProductionName() != fProdName))
      return kFALSE;
  if(fAliceRunValidationFlag)
    if(RunTag->GetRunValidation())
      return kFALSE;
  if (fAliceRunQualitiesFlag) {
    TObjArray *tQualities = fAliceRunQualities.Tokenize(" ");
    Bool_t tQual = kFALSE;

    TString tRQual = "";
    tRQual += RunTag->GetRunQuality();

    for (int iqual=0; iqual<tQualities->GetEntries(); iqual++)
      if (((TObjString *) tQualities->At(iqual))->GetString().Contains(tRQual))
	tQual = kTRUE;
	//      if (EvTag->GetFiredTriggerClasses().Contains(((TObjString *) tClasses->At(iqual))->GetString()))
  
    tQualities->Delete();
    delete tQualities;

    if (!tQual)
      return kFALSE;
  }
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
  if (fBeamTriggerFlag)
    if ((RunTag->GetBeamTriggers() < fBeamTriggerRange[0]) || (RunTag->GetBeamTriggers() > fBeamTriggerRange[1]))
      return kFALSE;
  if (fCollisionTriggerFlag)
    if ((RunTag->GetCollisionTriggers() < fCollisionTriggerRange[0]) || (RunTag->GetCollisionTriggers() > fCollisionTriggerRange[1]))
      return kFALSE;
  if (fEmptyTriggerFlag)
    if ((RunTag->GetEmptyTriggers() < fEmptyTriggerRange[0]) || (RunTag->GetEmptyTriggers() > fEmptyTriggerRange[1]))
      return kFALSE;
  if (fCSideTriggerFlag)
    if ((RunTag->GetCSideTriggers() < fCSideTriggerRange[0]) || (RunTag->GetCSideTriggers() > fCSideTriggerRange[1]))
      return kFALSE;
  if (fASideTriggerFlag)
    if ((RunTag->GetASideTriggers() < fASideTriggerRange[0]) || (RunTag->GetASideTriggers() > fASideTriggerRange[1]))
      return kFALSE;
  if (fHMTriggerFlag)
    if ((RunTag->GetHMTriggers() < fHMTriggerRange[0]) || (RunTag->GetHMTriggers() > fHMTriggerRange[1]))
      return kFALSE;
  if (fMuonTriggerFlag)
    if ((RunTag->GetMuonTriggers() < fMuonTriggerRange[0]) || (RunTag->GetMuonTriggers() > fMuonTriggerRange[1]))
      return kFALSE;
  if (fCollisionRateFlag)
    if ((RunTag->GetCollisionRate() < fCollisionRateRange[0]) || (RunTag->GetCollisionRate() > fCollisionRateRange[1]))
      return kFALSE;
  if (fMeanVertexFlag)
    if ((RunTag->GetMeanVertex() < fMeanVertexRange[0]) || (RunTag->GetMeanVertex() > fMeanVertexRange[1]))
      return kFALSE;
  if (fVertexQualityFlag)
    if ((RunTag->GetVertexQuality() < fVertexQualityRange[0]) || (RunTag->GetVertexQuality() > fVertexQualityRange[1]))
      return kFALSE;

  return kTRUE;
}
