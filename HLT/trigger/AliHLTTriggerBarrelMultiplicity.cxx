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
#include "TObjArray.h"
#include "TObjString.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTriggerBarrelMultiplicity)

AliHLTTriggerBarrelMultiplicity::AliHLTTriggerBarrelMultiplicity()
  : AliHLTTrigger()
  , fPtMin(0.0)
  , fPtMax(0.0)
  , fMinTracks(1)
  , fDCAReference()
  , fMinLDca(-1.)
  , fMaxLDca(-1.)
  , fMinTDca(-1.)
  , fMaxTDca(-1.)
  , fSolenoidBz(0.0)
  , fName()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

  for (int i=0; i<fgkDCAReferenceSize; i++) fDCAReference[i]=0.0;
}

const char* AliHLTTriggerBarrelMultiplicity::fgkDefaultOCDBEntry="HLT/ConfigHLT/BarrelMultiplicityTrigger";

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

  bool condition=false;
  if (iResult>=0 && numberOfTracks>=0) {
    if (fPtMax>fPtMin) {
      ptcut.Form(" %.02f GeV/c <= pt < %.02f GeV/c", fPtMin, fPtMax);
    } else {
      ptcut.Form(" pt >= %.02f GeV/c", fPtMin);
    }

    if (fMinTDca>=0.0) {
      if (fMaxTDca>=0.0) {
	tdca.Form(", %.02f<=transverse_dca<=%.02f", fMinTDca, fMaxTDca);
      } else {
	tdca.Form(" transverse_dca >= %.02f", fMinTDca);
      }
    } else if (fMaxTDca>=0.0) {
	tdca.Form(" transverse_dca<=%.02f", fMaxTDca);
    }
    if (!tdca.IsNull()) {
      if (op1st.IsNull()) op1st=" && ";
      else op2nd=" && ";
    }

    if (fMinLDca>=0.0) {
      if (fMaxLDca>=0.0) {
	ldca.Form(" %.02f<=longitudinal_dca<=%.02f", fMinLDca, fMaxLDca);
      } else {
	ldca.Form(" longitudinal_dca >= %.02f", fMinLDca);
      }
    } else if (fMaxLDca>=0.0) {
	ldca.Form(" longitudinal_dca<=%.02f", fMaxLDca);
    }
    if (!ldca.IsNull()) {
      if (op1st.IsNull()) op1st=" && ";
      else op2nd=" && ";
    }

    if (fMinTDca>=0.0 || fMaxTDca>=0 || fMinLDca>=0.0 || fMaxLDca>=0) {
      dcaref.Form(" (%.01f,%.01f,%.01f)", fDCAReference[0], fDCAReference[1], fDCAReference[2]);
    }

    if (numberOfTracks>=fMinTracks) {
      description.Form("Event contains %d track(s) with ", numberOfTracks);
      description+=ptcut;
      description+=op1st;
      description+=ldca;
      description+=op2nd;
      description+=tdca;
      description+=dcaref;
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
			      AliHLTReadoutList::kPHOS |
			      AliHLTReadoutList::kEMCAL
			      );
      condition=true;
    } else {
    description.Form("No tracks matching the tresholds found in the central barrel (min tracks %d, ", fMinTracks);
    description+=ptcut;
    description+=op1st;
    description+=ldca;
    description+=op2nd;
    description+=tdca;
    description+=dcaref;
    description+=")";
    }
  } else {
    description.Form("No input blocks found");
  }

  // add a specific trigger decision object with initialized name
  // the readout list however is fixed 
  AliHLTTriggerDecision decision(
				 condition,
				 fName.IsNull()?GetTriggerName():fName.Data(),
				 GetReadoutList(),
				 description.Data()
				 );
  TriggerEvent(&decision, kAliHLTDataTypeTObject|kAliHLTDataOriginOut);

  return iResult;
}

template<class T>
bool AliHLTTriggerBarrelMultiplicity::CheckCondition(T* track, float b)
{
  // see header file for class documentation
  if (!track) return false;

  // check on ptransverse momentum
  if (TMath::Abs(track->Pt()) < fPtMin || (fPtMax>fPtMin && TMath::Abs(track->Pt()) > fPtMax)) {
    return false;
  }

  // check on transverse and longitudinal DCA
  if (fMinTDca>=0.0 || fMaxTDca>=0 || fMinLDca>=0.0 || fMaxLDca>=0) {
    Float_t dz[2]={0.0,0.0};
    track->GetDZ(fDCAReference[0], fDCAReference[1], fDCAReference[2], b, dz);
    HLTDebug("checking dca condition: transversal %f logitudinal %f", dz[0], dz[1]);
    if (fMinTDca>=0 && TMath::Abs(dz[0])<fMinTDca) return false;
    if (fMaxTDca>=0 && TMath::Abs(dz[0])>fMaxTDca) return false;
    if (fMinLDca>=0 && TMath::Abs(dz[1])<fMinLDca) return false;
    if (fMaxLDca>=0 && TMath::Abs(dz[1])>fMaxLDca) return false;
  }

  return true;
}

int AliHLTTriggerBarrelMultiplicity::DoInit(int argc, const char** argv)
{
  // see header file for class documentation
  int iResult=0;

  // check if the -triggername argument is used
  // the name of the trigger determines the following initialization
  vector<const char*> remainingArgs;
  for (int i=0; i<argc; i++) {
    if (strcmp(argv[i], "-triggername")==0) {
      if (++i<argc) fName=argv[i];
      else {
	HLTError("invalid parameter for argument '-triggername', string expected");
	return -EINVAL;
      }
      continue;
    }
    remainingArgs.push_back(argv[i]);
  }

  // first configure the default
  TString cdbPath;
  if (!fName.IsNull()) {
    cdbPath="HLT/ConfigHLT/";
    cdbPath+=fName;
  } else {
    cdbPath=fgkDefaultOCDBEntry;
  }
  iResult=ConfigureFromCDBTObjString(cdbPath);

  // configure from the command line parameters if specified
  if (iResult>=0 && argc>0)
    iResult=ConfigureFromArgumentString(remainingArgs.size(), &(remainingArgs[0]));

  fSolenoidBz=GetBz();
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
  TString cdbPath;
  if (!cdbEntry || cdbEntry[0]==0) {
    if (!fName.IsNull()) {
      cdbPath="HLT/ConfigHLT/";
      cdbPath+=fName;
    } else {
      cdbPath=fgkDefaultOCDBEntry;
    }
  } else {
    cdbPath=cdbEntry;
  }

  return ConfigureFromCDBTObjString(cdbPath);
}

int AliHLTTriggerBarrelMultiplicity::ReadPreprocessorValues(const char* /*modules*/)
{
  // see header file for class documentation

  // nothing to do for the moment
  return 0;
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
    HLTWarning("argument -solenoidBz is deprecated, magnetic field set up globally (%f)", GetBz());
    return 2;
  }

  // unknown argument
  return -EINVAL;
}
