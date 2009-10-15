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

/// @file   AliHLTTriggerBarrelGeomMultiplicity.cxx
/// @author Oystein Djuvsland
/// @date   2009-10-08
/// @brief  HLT trigger component for charged particle multiplicity in
///         the central barrel.

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTTriggerBarrelGeomMultiplicity.h"
#include "AliESDEvent.h"
#include "AliHLTTriggerDecision.h"
#include "AliHLTDomainEntry.h"
#include "AliHLTGlobalBarrelTrack.h"
#include "TObjArray.h"
#include "TObjString.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTriggerBarrelGeomMultiplicity)

AliHLTTriggerBarrelGeomMultiplicity::AliHLTTriggerBarrelGeomMultiplicity()
  : AliHLTTrigger()
  , fMinTracks(1)
  , fSolenoidBz(0.0)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

}

const char* AliHLTTriggerBarrelGeomMultiplicity::fgkOCDBEntry="HLT/ConfigHLT/BarrelGeomMultiplicityTrigger";

AliHLTTriggerBarrelGeomMultiplicity::~AliHLTTriggerBarrelGeomMultiplicity()
{
  // see header file for class documentation
}

const char* AliHLTTriggerBarrelGeomMultiplicity::GetTriggerName() const
{
  // see header file for class documentation
  return "BarrelGeomMultiplicityTrigger";
}

AliHLTComponent* AliHLTTriggerBarrelGeomMultiplicity::Spawn()
{
  // see header file for class documentation
  return new AliHLTTriggerBarrelGeomMultiplicity;
}

int AliHLTTriggerBarrelGeomMultiplicity::DoTrigger()
{
  // see header file for class documentation
  int iResult=0;
  int numberOfTracks=-1;

  // try the ESD as input
  const TObject* obj = GetFirstInputObject(kAliHLTAllDataTypes, "AliESDEvent");
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(const_cast<TObject*>(obj));
  TString description;
  TString ptcut,tdca,ldca,dcaref,op1st,op2nd;
  if (esd != NULL) {
    numberOfTracks=0;
    esd->GetStdContent();
    
    for (Int_t i = 0; i < esd->GetNumberOfTracks(); i++) {
      if (CheckCondition(esd->GetTrack(i), esd->GetMagneticField())) numberOfTracks++;
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
	  if (CheckCondition(&(*element), fSolenoidBz)) numberOfTracks++;
	}
      } else if (iResult<0) {
	HLTError("can not extract tracks from data block of type %s (specification %08x) of size %d: error %d", 
		 DataType2Text(pBlock->fDataType).c_str(), pBlock->fSpecification, pBlock->fSize, iResult);
      }
    }
  }
  if (numberOfTracks>=fMinTracks) {

    ApplyTrigger();

    description.Form("Event contains %d track(s) with ", numberOfTracks);
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
    description.Form("No tracks matching the tresholds found in the central barrel (min tracks %d, ", fMinTracks);
    description+=ptcut;
    description+=op1st;
    description+=ldca;
    description+=op2nd;
    description+=tdca;
    description+=dcaref;
    description+=")";
  } else {
    description.Form("No input blocks found");
  }
  SetDescription(description.Data());
  TriggerEvent(false);
  return iResult;
}

template<class T>
bool AliHLTTriggerBarrelGeomMultiplicity::CheckCondition(T* track, float b)
{
  // see header file for class documentation
  if (!track) return false;

  ret = IsInDetectors(track, b);

  return ret;

}

template<class T>
bool AliHLTTriggerBarrelGeomMultiplicity::IsInDetector(T* track, b)
{

  for(Int_t i = 0; i < fDetectorList->GetEntries(); i++)
    {
      AliHLTTriggerDetectorGeom *det = static_cast<AliHLTTriggerDetectorGeom*>(fDetectorList.At(i))
      Double_t trackPoint[3];

      det->GetInitialPoint(trackPoint);

      bool ret = track->Intersect(trackPoint, det->NormVector(), b);

      if(track->Eta() >= det->EtaMin() && 
	 track->Eta() <= det->EtaMax() &&
	 track->Phi() >= det->PhiMin() &&
	 track->Phi() <= det->PhiMax())
	{
	  return true;
	}
    }
  return false;
}

int AliHLTTriggerBarrelGeomMultiplicity::DoInit(int argc, const char** argv)
{
  // see header file for class documentation

  // first configure the default
  int iResult=0;
  iResult=ConfigureFromCDBTObjString(kAliHLTCDBSolenoidBz);
  if (iResult>=0) iResult=ConfigureFromCDBTObjString(fgkOCDBEntry);

  // configure from the command line parameters if specified
  if (iResult>=0 && argc>0)
    iResult=ConfigureFromArgumentString(argc, argv);
  return iResult;
}

int AliHLTTriggerBarrelGeomMultiplicity::DoDeinit()
{
  // see header file for class documentation
  return 0;
}

int AliHLTTriggerBarrelGeomMultiplicity::Reconfigure(const char* cdbEntry, const char* /*chainId*/)
{
  // see header file for class documentation

  // configure from the specified antry or the default one
  const char* entry=cdbEntry;
  if (!entry || entry[0]==0) {
    ConfigureFromCDBTObjString(kAliHLTCDBSolenoidBz);
    entry=fgkOCDBEntry;
  }

  return ConfigureFromCDBTObjString(entry);
}

int AliHLTTriggerBarrelGeomMultiplicity::ReadPreprocessorValues(const char* /*modules*/)
{
  // see header file for class documentation

  // TODO 2009-09-10: implementation
  // for the moment very quick, just reload the magnetic field
  return ConfigureFromCDBTObjString(kAliHLTCDBSolenoidBz);
}

int AliHLTTriggerBarrelGeomMultiplicity::ScanConfigurationArgument(int argc, const char** argv)
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

  // -dca-reference
  // reference point for the transverse and longitudinal dca cut
  if (argument.CompareTo("-dca-reference")==0) {
    if (++i>=argc) return -EPROTO;
    argument=argv[i];
    // scan x,y,z
    TObjArray* pTokens=argument.Tokenize("'");
    if (pTokens) {
      for (int c=0; c<pTokens->GetEntriesFast() && c<fgkDCAReferenceSize; c++) {
	argument=((TObjString*)pTokens->At(c))->GetString();
	fDCAReference[i]=argument.Atof();
      }
      delete pTokens;
    }
    return 2;
  }

  // -min-ldca
  // minimum longitudinal dca to reference point
  if (argument.CompareTo("-min-ldca")==0) {
    if (++i>=argc) return -EPROTO;
    argument=argv[i];
    fMinLDca=argument.Atof();
    return 2;
  }
  
  // -max-ldca
  // maximum longitudinal dca to reference point
  if (argument.CompareTo("-max-ldca")==0) {
    if (++i>=argc) return -EPROTO;
    argument=argv[i];
    fMaxLDca=argument.Atof();
    return 2;
  }

  // -min-tdca
  // minimum transverse dca to reference point
  if (argument.CompareTo("-min-tdca")==0) {
    if (++i>=argc) return -EPROTO;
    argument=argv[i];
    fMinTDca=argument.Atof();
    return 2;
  }
  
  // -max-tdca
  // maximum transverse dca to reference point
  if (argument.CompareTo("-max-tdca")==0) {
    if (++i>=argc) return -EPROTO;
    argument=argv[i];
    fMaxTDca=argument.Atof();
    return 2;
  }

  // -solenoidBz
  if (argument.CompareTo("-solenoidBz")==0) {
    if (++i>=argc) return -EPROTO;
    argument=argv[i];
    fSolenoidBz=argument.Atof();
    return 2;
  }

  // unknown argument
  return -EINVAL;
}
