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
    fAliceRunQuality(0),
    fAliceBeamEnergy(0.0),
    fAliceBeamType(0),
    fAliceCalibrationVersion(0),
    fAliceDataType(0),
    fNumEvents(0),
    fNumDetectors(0),
    fEventTag("AliEventTag", 1000),
    fDetectorTag(),
    fLHCTag(), 
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
  fEventTag.Delete();
  if ( fQAArray ) 
    delete [] fQAArray ; 
  if ( fEventSpecies )
    delete [] fEventSpecies ; 
}

//___________________________________________________________________________
AliRunTag::AliRunTag(const AliRunTag& tag):
TObject(),
fAliceRunId(tag.fAliceRunId),
fAliceMagneticField(tag.fAliceMagneticField),
fAliceDipoleField(tag.fAliceDipoleField),
fAliceRunStartTime(tag.fAliceRunStartTime),
fAliceRunStopTime(fAliceRunStopTime),
fAlirootVersion(tag.fAlirootVersion),
fRootVersion(tag.fRootVersion),
fGeant3Version(tag.fGeant3Version),
fLHCPeriod(tag.fLHCPeriod),
fRecPass(tag.fRecPass),
fProductionName(tag.fProductionName),
fAliceRunQuality(tag.fAliceRunQuality),
fAliceBeamEnergy(tag.fAliceBeamEnergy),
fAliceBeamType(tag.fAliceBeamType),
fAliceCalibrationVersion(tag.fAliceCalibrationVersion),
fAliceDataType(tag.fAliceDataType),
fNumEvents(tag.fNumEvents),
fNumDetectors(tag.fNumDetectors),
fEventTag(tag.fEventTag),
fDetectorTag(tag.fDetectorTag),
fLHCTag(tag.fLHCTag), 
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
    fAliceRunQuality          = tag.fAliceRunQuality ; 
    fAliceBeamEnergy          = tag.fAliceBeamEnergy ;
    fAliceBeamType            = tag.fAliceBeamType ; 
    fAliceCalibrationVersion  = tag.fAliceCalibrationVersion ; 
    fAliceDataType            = tag.fAliceDataType ; 
    fNumEvents                = tag.fNumEvents ;
    fNumDetectors             = tag.fNumDetectors ; 
    fEventTag                 = tag.fEventTag ;
    fDetectorTag              = tag.fDetectorTag ;
    fLHCTag                   = tag.fLHCTag ;  
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
  SetRunQuality(oldtag->GetRunQuality());
  SetBeamEnergy(oldtag->GetBeamEnergy());
  SetBeamType(oldtag->GetBeamType());
  SetCalibVersion(oldtag->GetCalibVersion());
  SetDataType(oldtag->GetDataType());
  SetLHCTag(oldtag->GetLHCTag()->GetLuminosity(),oldtag->GetLHCTag()->GetLHCState());
  SetDetectorTag(oldtag->GetDetectorTags()->GetIntDetectorMaskDAQ(), oldtag->GetDetectorTags()->GetIntDetectorMaskReco());
  SetQA(*(oldtag->GetQA())) ;  	
  SetQAArray(oldtag->GetQAArray(), oldtag->GetQALength()) ;  
  SetEventSpecies(oldtag->GetEventSpecies(), oldtag->GetESLength()) ;  
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
  fLHCTag.SetLHCTag(lumin,type);
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
  new(fEventTag[fNumEvents++]) AliEventTag(EvTag);
}

//___________________________________________________________________________
void AliRunTag::Clear(const char *) {
  //Resets the number of events and detectors
  fEventTag.Delete();
  fNumEvents = 0;
  fDetectorTag.Clear();
  fNumDetectors = 0;
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
