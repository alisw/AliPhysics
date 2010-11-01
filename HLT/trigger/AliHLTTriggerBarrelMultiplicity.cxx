// $Id$
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                  Jochen Thaeder <jochen@thaeder.de>                    *
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
/// @author Matthias Richter, Jochen Thaeder
/// @date   2009-06-30
/// @brief  HLT trigger component for charged particle multiplicity in
///         the central barrel.

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTTriggerBarrelMultiplicity.h"
#include "AliHLTESDTrackCuts.h"
#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "AliHLTTriggerDecision.h"
#include "AliHLTDomainEntry.h"
#include "AliHLTGlobalBarrelTrack.h"
#include "AliHLTErrorGuard.h"
#include "TObjArray.h"
#include "TObjString.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTriggerBarrelMultiplicity)

AliHLTTriggerBarrelMultiplicity::AliHLTTriggerBarrelMultiplicity()
  : AliHLTTrigger()
  , fHLTESDTrackCuts(NULL)
  , fMinTracks(1)
  , fName()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

const char* AliHLTTriggerBarrelMultiplicity::fgkDefaultOCDBEntry="HLT/ConfigHLT/BarrelMultiplicityTrigger";

AliHLTTriggerBarrelMultiplicity::~AliHLTTriggerBarrelMultiplicity()
{
  // see header file for class documentation
}

const char* AliHLTTriggerBarrelMultiplicity::GetTriggerName() const
{
  // see header file for class documentation

  if (!fName.IsNull())
    return fName.Data();
  else
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

  if (!IsDataEvent()) {
    IgnoreEvent();  // dont generate any trigger decision.
  }

  int iResult=0;
  int numberOfTracks=-1;

  // try the ESD as input
  const TObject* obj = GetFirstInputObject(kAliHLTAllDataTypes, "AliESDEvent");
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(const_cast<TObject*>(obj));
  
  if (esd != NULL) {
    numberOfTracks=0;
    esd->GetStdContent();
    
    for (Int_t i = 0; i < esd->GetNumberOfTracks(); i++) {
      AliESDtrack *esdTrack = esd->GetTrack(i);
      if ( !esdTrack )
	continue;

      if ( fHLTESDTrackCuts->IsSelected(esdTrack) )
	numberOfTracks++;
    }
  }

  // try the AliHLTExternal track data as input
  // TODO: 2010-08-27
  // AliHLTTrackCuts needs an AliESDtrack object and not just AliExternalTrackParam
  // this part needs to be revised to work correctly with the track array as input
  // - think about specific conversion method in AliHLTGlobalBarrelTrack
  // - make sure that all necessary parameters are set
  // - clarify what to do about the track flags
  if (iResult>=0 && numberOfTracks<0) {
    for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(kAliHLTDataTypeTrack);
	 pBlock!=NULL; pBlock=GetNextInputBlock()) {
      if (numberOfTracks<0) numberOfTracks=0;
      vector<AliHLTGlobalBarrelTrack> tracks;
      if ((iResult=AliHLTGlobalBarrelTrack::ConvertTrackDataArray(reinterpret_cast<const AliHLTTracksData*>(pBlock->fPtr), pBlock->fSize, tracks))>0) {
	for (vector<AliHLTGlobalBarrelTrack>::iterator element=tracks.begin();
	     element!=tracks.end(); element++) {
	  ALIHLTERRORGUARD(1, "component needs to be revised to work with track array as input");
	  // TODO CHECK CONDITION HERE
	}
      } else if (iResult<0) {
	HLTError("can not extract tracks from data block of type %s (specification %08x) of size %d: error %d", 
		 DataType2Text(pBlock->fDataType).c_str(), pBlock->fSpecification, pBlock->fSize, iResult);
      }
    }
  }

  bool condition=false;
  TString description;

  if (iResult>=0 && numberOfTracks>=0) {
    if (numberOfTracks>=fMinTracks) {
      description.Form("Event contains %d track(s) with : ", numberOfTracks);
      description += fHLTESDTrackCuts->GetTitle();
      condition=true;
    } else {
      description.Form("No tracks matching the tresholds found in the central barrel (min tracks %d) with : ", fMinTracks);
      description += fHLTESDTrackCuts->GetTitle();
    }
  } else {
    if(IsDataEvent()) {
      description.Form("No input blocks found");
    } else {
      description.Form("No DataEvent found");
    }
  }
  
  SetDescription(description.Data());

  // add a specific trigger decision object with initialized name
  // the readout list however is fixed 
  AliHLTTriggerDecision decision(
				 condition,
				 GetTriggerName(),
				 GetReadoutList(),
				 GetDescription()
				 );
  TriggerEvent(&decision, kAliHLTDataTypeTObject|kAliHLTDataOriginOut);

  return iResult;
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

  // get path from triggername, use default object otherwise
  TString cdbPath;
  if (!fName.IsNull()) {
    cdbPath="HLT/ConfigHLT/";
    cdbPath+=fName;
  } else {
    cdbPath=fgkDefaultOCDBEntry;
  }

  // -- Check if CDB object is AliHLTESDTrackCuts or TObjString 
  //    and configure from it. Replace "-" by "_._" if needed in the cdbPath
  iResult = ConfigureFromCDBObject(cdbPath);

  // -- Configure from the command line parameters if specified
  if (iResult>=0 && argc>0)
    iResult=ConfigureFromArgumentString(remainingArgs.size(), &(remainingArgs[0]));

  // -- Check if we have the track cuts for triggering
  if (!fHLTESDTrackCuts) {
    HLTError("No AliHLTESDTrackCuts object has been created as basis for triggering.");
    iResult=-ENOENT;
  }
  else {
    if (!fName.IsNull()) {
      if (fName.Contains("Barrel_pT_Single"))
	fMinTracks = 1;
    }
  }

  return iResult;
}

int AliHLTTriggerBarrelMultiplicity::DoDeinit()
{
  // see header file for class documentation

  if (fHLTESDTrackCuts)
    delete fHLTESDTrackCuts;
  fHLTESDTrackCuts = NULL;

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

  return ConfigureFromCDBObject(cdbPath);
}

int AliHLTTriggerBarrelMultiplicity::ReadPreprocessorValues(const char* /*modules*/)
{
  // see header file for class documentation

  // nothing to do for the moment
  return 0;
}

Int_t AliHLTTriggerBarrelMultiplicity::ConfigureFromCDBObject(TString cdbPath)
{
  // see header file for class documentation

  Int_t iResult = 0;
  TString arguments;

  // -- check for "-" and replace by "_._" in the path name
  cdbPath.ReplaceAll("-",1,"_._",3);

  TObject* pCDBObject = LoadAndExtractOCDBObject(cdbPath);
  if (pCDBObject) {
    AliHLTESDTrackCuts *pCuts = dynamic_cast<AliHLTESDTrackCuts*>(pCDBObject);
    if (pCuts) {
      HLTInfo("Received AliHLTESDTrackCuts configuration object : \'%s\'", pCuts->GetTitle());
      if (fHLTESDTrackCuts)
	delete fHLTESDTrackCuts;
      fHLTESDTrackCuts = pCuts;
    }
    else {
      TObjString* pString = dynamic_cast<TObjString*>(pCDBObject);
      if (pString) {
	HLTInfo("Received configuration object string: \'%s\'", pString->GetString().Data());
	arguments+=pString->GetString().Data();
      } 
      else {
	HLTError("Configuration object \"%s\" has wrong type, required AliHLTESDTrackCuts or TObjString", cdbPath.Data());
	iResult=-EINVAL;
      }
    }
  } 
  else {
    HLTError("Can not fetch object \"%s\" from CDB", cdbPath.Data());
    iResult=-ENOENT;
  }
  
  if ( iResult>=0 && !arguments.IsNull() ) {
    const Char_t* array = arguments.Data();
    iResult = ConfigureFromArgumentString(1, &array);
  }

  return iResult;
}

int AliHLTTriggerBarrelMultiplicity::ScanConfigurationArgument(int argc, const char** argv)
{
  // see header file for class documentation
  if (argc<=0) return 0;
  int i=0;
  TString argument=argv[i];

  if (!fHLTESDTrackCuts)
    fHLTESDTrackCuts = new AliHLTESDTrackCuts("AliHLTESDTrackCuts","No track cuts");

  // -maxpt
  if (argument.CompareTo("-maxpt")==0) {
    if (++i>=argc) return -EPROTO;
    argument=argv[i];

    Float_t minPt, maxPt;
    fHLTESDTrackCuts->GetPtRange(minPt,maxPt);
    maxPt = argument.Atof(); 
    fHLTESDTrackCuts->SetPtRange(minPt,maxPt);

    TString title = fHLTESDTrackCuts->GetTitle();
    if (!title.CompareTo("No track cuts")) title = "";
    else title += " && ";
    title += Form("p_t < %f", maxPt);
    fHLTESDTrackCuts->SetTitle(title);
    return 2;
  }    

  // -minpt
  if (argument.CompareTo("-minpt")==0) {
    if (++i>=argc) return -EPROTO;
    argument=argv[i];

    Float_t minPt, maxPt;
    fHLTESDTrackCuts->GetPtRange(minPt,maxPt);
    minPt = argument.Atof(); 
    fHLTESDTrackCuts->SetPtRange(minPt,maxPt);

    TString title = fHLTESDTrackCuts->GetTitle();
    if (!title.CompareTo("No track cuts")) title = "";
    else title += " && ";
    title += Form("p_t > %f", minPt);
    fHLTESDTrackCuts->SetTitle(title);
    return 2;
  }    

  // -mintracks
  if (argument.CompareTo("-mintracks")==0) {
    if (++i>=argc) return -EPROTO;
    argument=argv[i];
    fMinTracks=argument.Atoi();
    return 2;
  }    

  // -min-ldca
  // minimum longitudinal dca to vertex
  if (argument.CompareTo("-min-ldca")==0) {
    if (++i>=argc) return -EPROTO;
    argument=argv[i];

    fHLTESDTrackCuts->SetMinDCAToVertexZ(argument.Atof());
    TString title = fHLTESDTrackCuts->GetTitle();
    if (!title.CompareTo("No track cuts")) title = "";
    else title += " && ";
    title += Form("DCAz > %f", argument.Atof());
    fHLTESDTrackCuts->SetTitle(title);
    return 2;
  }
  
  // -max-ldca
  // maximum longitudinal dca to vertex
  if (argument.CompareTo("-max-ldca")==0) {
    if (++i>=argc) return -EPROTO;
    argument=argv[i];

    fHLTESDTrackCuts->SetMaxDCAToVertexZ(argument.Atof());
    TString title = fHLTESDTrackCuts->GetTitle();
    if (!title.CompareTo("No track cuts")) title = "";
    else title += " && ";
    title += Form("DCAz < %f", argument.Atof());
    fHLTESDTrackCuts->SetTitle(title);
    return 2;
  }

  // -min-tdca
  // minimum transverse dca to vertex
  if (argument.CompareTo("-min-tdca")==0) {
    if (++i>=argc) return -EPROTO;
    argument=argv[i];

    fHLTESDTrackCuts->SetMinDCAToVertexXY(argument.Atof());
    TString title = fHLTESDTrackCuts->GetTitle();
    if (!title.CompareTo("No track cuts")) title = "";
    else title += " && ";
    title += Form("DCAr > %f", argument.Atof());
    fHLTESDTrackCuts->SetTitle(title);
    return 2;
  }
  
  // -max-tdca
  // maximum transverse dca to vertex
  if (argument.CompareTo("-max-tdca")==0) {
    if (++i>=argc) return -EPROTO;
    argument=argv[i];

    fHLTESDTrackCuts->SetMaxDCAToVertexXY(argument.Atof());
    TString title = fHLTESDTrackCuts->GetTitle();
    if (!title.CompareTo("No track cuts")) title = "";
    else title += " && ";
    title += Form("DCAr < %f", argument.Atof());
    fHLTESDTrackCuts->SetTitle(title);
    return 2;
  }

  // -- deprecated

  // -dca-reference
  // reference point for the transverse and longitudinal dca cut
  if (argument.CompareTo("-dca-reference")==0) {
    if (++i>=argc) return -EPROTO;
    HLTWarning("argument -dca-reference deprecated, ESDTrackCuts only allow for DCA to vertex");
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
