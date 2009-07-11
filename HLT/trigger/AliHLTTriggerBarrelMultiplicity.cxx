// $Id$
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/// @file   AliHLTTriggerBarrelMultiplicity.cxx
/// @author Matthias Richter
/// @date   2009-06-30
/// @brief  HLT trigger component for charged particle multiplicity in
///         the central barrel.

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTTriggerBarrelMultiplicity.h"
#include "AliESDEvent.h"
#include "AliHLTTriggerDecision.h"
#include "AliHLTDomainEntry.h"
#include "AliHLTGlobalBarrelTrack.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTriggerBarrelMultiplicity)

AliHLTTriggerBarrelMultiplicity::AliHLTTriggerBarrelMultiplicity()
  : AliHLTTrigger()
  , fPtMin(0.0)
  , fPtMax(0.0)
  , fMinTracks(1)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

const char* AliHLTTriggerBarrelMultiplicity::fgkOCDBEntry="HLT/ConfigHLT/BarrelMultiplicityTrigger";

AliHLTTriggerBarrelMultiplicity::~AliHLTTriggerBarrelMultiplicity()
{
  // see header file for class documentation
}

const char* AliHLTTriggerBarrelMultiplicity::GetTriggerName() const
{
  // see header file for class documentation
  return "BarrelMultiplicityTrigger";
}

AliHLTComponent* AliHLTTriggerBarrelMultiplicity::Spawn()
{
  // see header file for class documentation
  return new AliHLTTriggerBarrelMultiplicity;
}

int AliHLTTriggerBarrelMultiplicity::DoTrigger()
{
  // see header file for class documentation
  int iResult=0;
  int numberOfTracks=-1;

  // try the ESD as input
  const TObject* obj = GetFirstInputObject(kAliHLTAllDataTypes, "AliESDEvent");
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(const_cast<TObject*>(obj));
  TString description;
  TString ptcut;
  if (esd != NULL) {
    numberOfTracks=0;
    esd->GetStdContent();
    
    for (Int_t i = 0; i < esd->GetNumberOfTracks(); i++) {
      if (CheckCondition(esd->GetTrack(i))) numberOfTracks++;
    }
  }

  // try the AliHLTExternal track data as input
  if (iResult>=0 && numberOfTracks<0) {
    for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(kAliHLTDataTypeTrack);
	 pBlock!=NULL; pBlock=GetNextInputBlock()) {
      if (numberOfTracks<0) numberOfTracks=0;
      vector<AliHLTGlobalBarrelTrack> tracks;
      if ((iResult=AliHLTGlobalBarrelTrack::ConvertTrackDataArray(reinterpret_cast<const AliHLTTracksData*>(pBlock->fPtr), pBlock->fSize, tracks))>0) {
	for (vector<AliHLTGlobalBarrelTrack>::iterator element=tracks.begin();
	     element!=tracks.end(); element++) {
	  if (CheckCondition(&(*element))) numberOfTracks++;
	}
      } else if (iResult<0) {
	HLTError("can not extract tracks from data block of type %s (specification %08x) of size %d: error %d", 
		 DataType2Text(pBlock->fDataType).c_str(), pBlock->fSpecification, pBlock->fSize, iResult);
      }
    }
  }

  if (iResult>=0 && numberOfTracks>=0) {
    if (fPtMax>fPtMin) {
      ptcut.Form(" %.02f GeV/c <= pt < %.02f GeV/c", fPtMin, fPtMax);
    } else {
      ptcut.Form(" pt >= %.02f GeV/c", fPtMin);
    }
    if (numberOfTracks>=fMinTracks) {
      description.Form("Event contains %d track(s) with ", numberOfTracks);
      description+=ptcut;
      SetDescription(description.Data());
      // Enable the central detectors for readout.
      GetReadoutList().Enable(
			      AliHLTReadoutList::kITSSPD |
			      AliHLTReadoutList::kITSSDD |
			      AliHLTReadoutList::kITSSSD |
			      AliHLTReadoutList::kTPC |
			      AliHLTReadoutList::kTRD |
			      AliHLTReadoutList::kTOF |
			      AliHLTReadoutList::kHMPID |
			      AliHLTReadoutList::kPHOS
			      );
      // Add the available HLT information for readout too.
      GetTriggerDomain().Add("CLUSTERS", "TPC ");
      TriggerEvent(true);
      return 0;
    }
    description.Form("No tracks matching the tresholds found in the central barrel (min tracks %d, %s)",
		     fMinTracks, ptcut.Data());
  } else {
    description.Form("No input blocks found");
  }
  SetDescription(description.Data());
  TriggerEvent(false);
  return iResult;
}

template<class T>
bool AliHLTTriggerBarrelMultiplicity::CheckCondition(T* track)
{
  // see header file for class documentation
  if (track && track->Pt() >= fPtMin &&
      (fPtMax<=fPtMin || track->Pt() < fPtMax)) {
    return true;
  }
  return false;
}

int AliHLTTriggerBarrelMultiplicity::DoInit(int argc, const char** argv)
{
  // see header file for class documentation

  // first configure the default
  int iResult=ConfigureFromCDBTObjString(fgkOCDBEntry);

  // configure from the command line parameters if specified
  if (iResult>=0 && argc>0)
    iResult=ConfigureFromArgumentString(argc, argv);
  return iResult;
}

int AliHLTTriggerBarrelMultiplicity::DoDeinit()
{
  // see header file for class documentation
  return 0;
}

int AliHLTTriggerBarrelMultiplicity::Reconfigure(const char* cdbEntry, const char* /*chainId*/)
{
  // see header file for class documentation

  // configure from the specified antry or the default one
  const char* entry=cdbEntry;
  if (!entry || entry[0]==0) entry=fgkOCDBEntry;

  return ConfigureFromCDBTObjString(entry);
}

int AliHLTTriggerBarrelMultiplicity::ScanConfigurationArgument(int argc, const char** argv)
{
  // see header file for class documentation
  if (argc<=0) return 0;
  int i=0;
  TString argument=argv[i];

  // -maxpt
  if (argument.CompareTo("-maxpt")==0) {
    if (++i>=argc) return -EPROTO;
    argument=argv[i];
    fPtMax=argument.Atof();
    return 2;
  }    

  // -minpt
  if (argument.CompareTo("-minpt")==0) {
    if (++i>=argc) return -EPROTO;
    argument=argv[i];
    fPtMin=argument.Atof();
    return 2;
  }    

  // -mintracks
  if (argument.CompareTo("-mintracks")==0) {
    if (++i>=argc) return -EPROTO;
    argument=argv[i];
    fMinTracks=argument.Atoi();
    return 2;
  }    
  
  // unknown argument
  return -EINVAL;
}
