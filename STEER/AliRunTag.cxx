/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
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
//           Implementation of the RunTag class
//   This is the class to deal with the tags in the run level
//   Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-----------------------------------------------------------------

#include <stdlib.h>
#include "AliRunTag.h"
#include "AliDetectorTag.h"
#include "AliEventTag.h"

ClassImp(AliRunTag)

//___________________________________________________________________________
  AliRunTag::AliRunTag() :
    TObject(),
    fAliceRunId(-1),
    fAliceMagneticField(0.0),
    fAliceDipoleField(0.0),
    fAliceRunStartTime(0),
    fAliceRunStopTime(0),
    fAlirootVersion(0),
    fRootVersion(0),
    fGeant3Version(0),
    fLHCPeriod(0),
    fRecPass(0),
    fProductionName(0),
    fAliceRunValidated(0),
    fAliceRunGlobalQuality(0),
    fAliceBeamEnergy(0.0),
    fAliceBeamType(0),
    fAliceCalibrationVersion(0),
    fAliceDataType(0),
    //    fNumEvents(0),
    fNumFiles(0),
    fBeamTriggers(0),
    fCollisionTriggers(0),
    fEmptyTriggers(0),
    fASideTriggers(0),
    fCSideTriggers(0),
    fHMTriggers(0),
    fMuonTriggers(0),
    fCollisionRate(0.0),
    fMeanVertex(0.0),
    fVertexQuality(0.0),
    fNumDetectors(0),
    fFileTags(100),
    fDetectorTag(),
    fLHCTag(), 
    fActiveTriggerClasses(""),
    fQA(),  
    fQALength(0), 
    fQAArray(NULL), 
    fESLength(0), 
    fEventSpecies(NULL)
{
  //Default constructor
}

//___________________________________________________________________________
AliRunTag::~AliRunTag() {
  //Destructor
  //  fEventTag.Delete();
  if ( fQAArray ) 
    delete [] fQAArray ; 
  if ( fEventSpecies )
    delete [] fEventSpecies ; 
  fFileTags.Delete();
}

//___________________________________________________________________________
AliRunTag::AliRunTag(const AliRunTag& tag):
TObject(),
fAliceRunId(tag.fAliceRunId),
fAliceMagneticField(tag.fAliceMagneticField),
fAliceDipoleField(tag.fAliceDipoleField),
fAliceRunStartTime(tag.fAliceRunStartTime),
fAliceRunStopTime(tag.fAliceRunStopTime),
fAlirootVersion(tag.fAlirootVersion),
fRootVersion(tag.fRootVersion),
fGeant3Version(tag.fGeant3Version),
fLHCPeriod(tag.fLHCPeriod),
fRecPass(tag.fRecPass),
fProductionName(tag.fProductionName),
fAliceRunValidated(tag.fAliceRunValidated),
fAliceRunGlobalQuality(tag.fAliceRunGlobalQuality),
fAliceBeamEnergy(tag.fAliceBeamEnergy),
fAliceBeamType(tag.fAliceBeamType),
fAliceCalibrationVersion(tag.fAliceCalibrationVersion),
fAliceDataType(tag.fAliceDataType),
//fNumEvents(tag.fNumEvents),
fNumFiles(0),
fBeamTriggers(tag.fBeamTriggers),
fCollisionTriggers(tag.fCollisionTriggers),
fEmptyTriggers(tag.fEmptyTriggers),
fASideTriggers(tag.fASideTriggers),
fCSideTriggers(tag.fCSideTriggers),
fHMTriggers(tag.fHMTriggers),
fMuonTriggers(tag.fMuonTriggers),
fCollisionRate(tag.fCollisionRate),
fMeanVertex(tag.fMeanVertex),
fVertexQuality(tag.fVertexQuality),
fNumDetectors(tag.fNumDetectors),
fFileTags(100),
fDetectorTag(tag.fDetectorTag),
fLHCTag(tag.fLHCTag), 
fActiveTriggerClasses(tag.fActiveTriggerClasses),
fQA(tag.fQA),
fQALength(tag.fQALength),
fQAArray(NULL), 
fESLength(tag.fESLength),
fEventSpecies(NULL)
{
  //copy constructor
  if (fQALength == 0 ) 
    fQAArray = NULL ; 
  else {
    fQAArray = new ULong_t[fQALength] ; 
    memcpy(fQAArray, tag.fQAArray, fQALength*sizeof(ULong_t)) ;
  }
  if (fESLength == 0 ) 
    fEventSpecies = NULL ; 
  else {
    fEventSpecies = new Bool_t[fESLength] ; 
    memcpy(fEventSpecies, tag.fEventSpecies, fESLength*sizeof(Bool_t)) ;
  }
  for (int ifl=0; ifl<tag.fNumFiles; ifl++) {
    AddFileTag(new AliFileTag(*tag.GetFileTag(ifl)));
  }
    

}

//___________________________________________________________________________
AliRunTag& AliRunTag::operator = (const AliRunTag& tag) {
//assignment operator
  if(&tag != this) {
    fAliceRunId               = tag.fAliceRunId ; 
    fAliceMagneticField       = tag.fAliceMagneticField ;
    fAliceDipoleField         = tag.fAliceDipoleField ;
    fAliceRunStartTime        = tag.fAliceRunStartTime ; 
    fAliceRunStopTime         = tag.fAliceRunStopTime ; 
    fAlirootVersion           = tag.fAlirootVersion ; 
    fRootVersion              = tag.fRootVersion ;
    fGeant3Version            = tag.fGeant3Version ; 
    fLHCPeriod                = tag.fLHCPeriod ; 
    fRecPass                  = tag.fRecPass ; 
    fProductionName           = tag.fProductionName ; 
    fAliceRunValidated        = tag.fAliceRunValidated ; 
    fAliceRunGlobalQuality    = tag.fAliceRunGlobalQuality ; 
    fAliceBeamEnergy          = tag.fAliceBeamEnergy ;
    fAliceBeamType            = tag.fAliceBeamType ; 
    fAliceCalibrationVersion  = tag.fAliceCalibrationVersion ; 
    fAliceDataType            = tag.fAliceDataType ; 
    //    fNumEvents                = tag.fNumEvents ;
    fNumFiles                 = tag.fNumFiles;
    fBeamTriggers             = tag.fBeamTriggers;
    fCollisionTriggers	      = tag.fCollisionTriggers;
    fEmptyTriggers	      = tag.fEmptyTriggers;
    fASideTriggers	      = tag.fASideTriggers;
    fCSideTriggers	      = tag.fCSideTriggers;
    fHMTriggers               = tag.fHMTriggers;
    fMuonTriggers             = tag.fMuonTriggers;
    fCollisionRate	      = tag.fCollisionRate;
    fMeanVertex		      = tag.fMeanVertex;
    fVertexQuality	      = tag.fVertexQuality;
    fNumDetectors             = tag.fNumDetectors ; 
    fDetectorTag              = tag.fDetectorTag ;
    fLHCTag                   = tag.fLHCTag ;  
    fActiveTriggerClasses     = tag.fActiveTriggerClasses;
    fQA                       = tag.fQA ;      
    fQALength                 = tag.fQALength ; 
    if (fQAArray) 
      delete [] fQAArray ; 
    if (fQALength == 0 ) 
      fQAArray = NULL ; 
    else {
      fQAArray = new ULong_t[fQALength] ; 
      memcpy(fQAArray, tag.fQAArray, fQALength*sizeof(ULong_t)) ;
    }
    fESLength                 = tag.fESLength ; 
    if (fEventSpecies)
      delete [] fEventSpecies ; 
    if (fESLength == 0 ) 
      fEventSpecies = NULL ; 
    else {
      fEventSpecies = new Bool_t[fESLength] ; 
      memcpy(fEventSpecies, tag.fEventSpecies, fESLength*sizeof(Bool_t)) ;
    }
    for (int ifl=0; ifl<fNumFiles; ifl++) {
      AddFileTag(new AliFileTag(*tag.GetFileTag(ifl)));
    }
//     for (int ifile=0; ifile<tag.GetFileTags()->GetEntries(); ifile++)
//       AddFileTag(*((AliFileTag *) tag.GetFileTags()->At(ifile)));
  }
  return *this ; 
}

//___________________________________________________________________________
void AliRunTag::CopyStandardContent(AliRunTag *oldtag) {
  //function that copies the run, lhc and detector levels
  SetRunId(oldtag->GetRunId());
  SetMagneticField(oldtag->GetMagneticField());
  SetDipoleField(oldtag->GetDipoleField());
  SetRunStartTime(oldtag->GetRunStartTime());
  SetRunStopTime(oldtag->GetRunStopTime());
  SetAlirootVersion(oldtag->GetAlirootVersion());
  SetRootVersion(oldtag->GetRootVersion());
  SetGeant3Version(oldtag->GetGeant3Version());
  SetLHCPeriod(oldtag->GetLHCPeriod());
  SetReconstructionPass(oldtag->GetReconstructionPass());
  SetProductionName(oldtag->GetProductionName());
  SetRunValidation(oldtag->GetRunValidation());
  SetRunQuality(oldtag->GetRunQuality());
  SetBeamEnergy(oldtag->GetBeamEnergy());
  SetBeamType(oldtag->GetBeamType());
  SetCalibVersion(oldtag->GetCalibVersion());
  SetDataType(oldtag->GetDataType());
  SetBeamTriggers(oldtag->GetBeamTriggers());
  SetCollisionTriggers(oldtag->GetCollisionTriggers());
  SetEmptyTriggers(oldtag->GetEmptyTriggers());
  SetASideTriggers(oldtag->GetASideTriggers());
  SetCSideTriggers(oldtag->GetCSideTriggers());
  SetHMTriggers(oldtag->GetHMTriggers());
  SetMuonTriggers(oldtag->GetMuonTriggers());
  SetCollisionRate(oldtag->GetCollisionRate());
  SetMeanVertex(oldtag->GetMeanVertex());
  SetVertexQuality(oldtag->GetVertexQuality());
  SetLHCTag(oldtag->GetLHCTag()->GetLuminosity(),oldtag->GetLHCTag()->GetLHCState());
  SetDetectorTag(oldtag->GetDetectorTags()->GetIntDetectorMaskDAQ(), oldtag->GetDetectorTags()->GetIntDetectorMaskReco());
  SetActiveTriggerClasses(oldtag->GetActiveTriggerClasses());
  SetQA(*(oldtag->GetQA())) ;  	
  SetQAArray(oldtag->GetQAArray(), oldtag->GetQALength()) ;  
  SetEventSpecies(oldtag->GetEventSpecies(), oldtag->GetESLength()) ;  
  for (int ifile=0; ifile<oldtag->GetNFiles(); ifile++) {
    AliFileTag *ntag = new AliFileTag();
    ntag->CopyFileInfo((const AliFileTag &) *(oldtag->GetFileTag(ifile)));
    AddFileTag(ntag);
  }
}

void AliRunTag::UpdateFromRunTable(AliRunTag *tabtag)
{
  SetBeamTriggers(tabtag->GetBeamTriggers());
  SetCollisionTriggers(tabtag->GetCollisionTriggers());
  SetEmptyTriggers(tabtag->GetEmptyTriggers());
  SetASideTriggers(tabtag->GetASideTriggers());
  SetCSideTriggers(tabtag->GetCSideTriggers());
  SetHMTriggers(tabtag->GetHMTriggers());
  SetMuonTriggers(tabtag->GetMuonTriggers());
  SetCollisionRate(tabtag->GetCollisionRate());
  SetMeanVertex(tabtag->GetMeanVertex());
  SetVertexQuality(tabtag->GetVertexQuality());
  SetRunQuality(tabtag->GetRunQuality());
  fLHCTag.UpdateFromRunTable(*tabtag->GetLHCTag());
  fDetectorTag.UpdateFromRunTable(*tabtag->GetDetectorTags());
}

//___________________________________________________________________________
void AliRunTag::SetQAArray(ULong_t * qa, Int_t qalength) {
  //Setter for the qa bits 
  if (qa && qalength > 0) {
    fQALength = qalength ; 
    if (fQAArray) 
      delete [] fQAArray ; 
    fQAArray = new ULong_t[qalength] ; 
    memcpy(fQAArray, qa, qalength*sizeof(ULong_t)) ;
  }
}

//___________________________________________________________________________
void AliRunTag::SetEventSpecies(Bool_t * es, Int_t eslength) {
  //setter for the eventspecices 
  if (es && eslength >0 ) {
    fESLength = eslength ; 
    if (fEventSpecies) 
      delete [] fEventSpecies ; 
    fEventSpecies = new Bool_t[eslength] ;
    memcpy(fEventSpecies, es, eslength*sizeof(Bool_t)) ; 
  }
}


//___________________________________________________________________________
void AliRunTag::SetLHCTag(Float_t lumin, TString type) {
  //Setter for the LHC tags
  fLHCTag.SetLuminosity(lumin);
  fLHCTag.SetLHCState(type);
}

//___________________________________________________________________________
void AliRunTag::SetDetectorTag(UInt_t mask, UInt_t maskReco) {
  //Setter for the detector tags
  fDetectorTag.SetDetectorMaskDAQ(mask);
  if (maskReco == 0)
    fDetectorTag.SetDetectorMaskReco(mask);
  else
    fDetectorTag.SetDetectorMaskReco(maskReco);

  int ndet = 0;
  for (int iter=0; iter<32; iter++)  
    ndet += (mask & (1 << iter)) > 0;
  
  fNumDetectors = ndet;
}

//___________________________________________________________________________
void AliRunTag::AddEventTag(const AliEventTag & EvTag) {
  //Adds an entry to the event tag TClonesArray
  ((AliFileTag *) fFileTags[fNumFiles-1])->AddEventTag(EvTag);
  //  new(fEventTag[fNumEvents++]) AliEventTag(EvTag);
}

void AliRunTag::AddFileTag(AliFileTag *t) {
  //Adds an entry for each file tag
  if (fNumFiles == fFileTags.GetSize()-1) fFileTags.Expand(fFileTags.GetSize()*2);
  //  new(fFileTags[fNumFiles++]) AliFileTag(t);
  fFileTags[fNumFiles++] = t;
}

//___________________________________________________________________________
void AliRunTag::Clear(const char *) {
  //Resets the number of events and detectors
  //  fEventTag.Delete();
  //  fNumEvents = 0;
  fFileTags.Delete();
  fNumFiles = 0;
  if ( fQAArray ) {
    delete [] fQAArray ;
    fQAArray = 0x0;
  } 
  fQALength=0;
  if ( fEventSpecies ) {
    delete [] fEventSpecies ;
    fEventSpecies = 0x0;
  } 
  fESLength=0;
}

const AliEventTag* AliRunTag::GetEventTag(int evt) const
{
  int curev = evt;
  int curf = 0;

  if (evt >= GetNEvents()) return 0;
 
  while (curev >= ((AliFileTag *) fFileTags[curf])->GetNEvents()) {
    curev -= ((AliFileTag *) fFileTags[curf])->GetNEvents();
    curf++;
  }
  return ((AliFileTag *) fFileTags[curf])->GetEventTag(curev);
}

AliFileTag *AliRunTag::GetFileTagForEvent(int evt) 
{
  // Returns FileTag in which the given event is
  int curev = evt;
  int curf = 0;

  if (evt >= GetNEvents()) return 0;
 
  while (curev >= ((AliFileTag *) fFileTags[curf])->GetNEvents()) {
    curev -= ((AliFileTag *) fFileTags[curf])->GetNEvents();
    curf++;
  }
  return (AliFileTag *) fFileTags[curf];
}

Int_t       AliRunTag::GetNEvents() const
{
  Int_t evtot = 0;
  for (int iter=0; iter<fNumFiles; iter++)
    evtot += ((AliFileTag *) fFileTags[iter])->GetNEvents();

  return evtot;
}

Int_t      AliRunTag::GetFileId(const char *guid)
{
  for (int iter=0; iter<fNumFiles; iter++) {
    if (!strcmp(((AliFileTag *) fFileTags[iter])->GetGUID(), guid))
      return iter;
  }
  return -1;
}

