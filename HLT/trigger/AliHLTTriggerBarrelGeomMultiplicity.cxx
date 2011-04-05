// $Id$
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Oystein Djuvsland                                     *
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
/// @brief  HLT trigger component for charged particle multiplicity 
///         within a geometrical acceptance in the central barrel.

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTTriggerBarrelGeomMultiplicity.h"
#include "AliHLTTriggerDetectorGeom.h"
#include "AliHLTTriggerDecisionParameters.h"
#include "AliESDEvent.h"
#include "AliHLTTriggerDecision.h"
#include "AliHLTDomainEntry.h"
#include "AliHLTGlobalBarrelTrack.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "TFile.h"
#include "AliHLTTrigger.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTriggerBarrelGeomMultiplicity)

AliHLTTriggerBarrelGeomMultiplicity::AliHLTTriggerBarrelGeomMultiplicity()
  : AliHLTTrigger()
  , fSolenoidBz(0)
  , fMinTracks(1)
  , fDetectorArray(0)
  , fTriggerDecisionPars(0)
  , fTriggerName(0)
  , fOCDBEntry(0)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

  fDetectorArray = new TObjArray;

}

AliHLTTriggerBarrelGeomMultiplicity::~AliHLTTriggerBarrelGeomMultiplicity()
{
  // see header file for class documentation

  if (fDetectorArray != NULL) delete fDetectorArray;
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

int AliHLTTriggerBarrelGeomMultiplicity::Reconfigure(const char *cdbEntry, const char *chainId)
{
  // see header file for class documentation

  // configure from the specified entry or the default
  const char* entry=cdbEntry;

  if (!entry)
    {
      HLTDebug("No CDB path specified");
      entry = fOCDBEntry; 
    }

  return GetDetectorGeomsFromCDBObject(entry, chainId);
} 

int AliHLTTriggerBarrelGeomMultiplicity::DoTrigger()
{
  // see header file for class documentation
  int iResult=0;
  int numberOfTracks=-1;

  if (!fTriggerDecisionPars) {
    iResult=-ENODEV;
  }

  // try the ESD as input
  const TObject* obj = GetFirstInputObject(kAliHLTAllDataTypes, "AliESDEvent");
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(const_cast<TObject*>(obj));
  TString description;

  if (esd != NULL) 
    {
      numberOfTracks=0;
      esd->GetStdContent();
      for (Int_t i = 0; i < esd->GetNumberOfTracks(); i++) 
	{
	  if (CheckCondition(esd->GetTrack(i), esd->GetMagneticField())) numberOfTracks++;
	}
    }

  // try the AliHLTExternal track data as input
  if (iResult>=0 && numberOfTracks<0) 
    {
      for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(kAliHLTDataTypeTrack);
	   pBlock!=NULL; pBlock=GetNextInputBlock()) 
	{
	  if (numberOfTracks<0) numberOfTracks=0;
	  vector<AliHLTGlobalBarrelTrack> tracks;
	  if ((iResult=AliHLTGlobalBarrelTrack::ConvertTrackDataArray(reinterpret_cast<const AliHLTTracksData*>(pBlock->fPtr), pBlock->fSize, tracks))>0) 
	    {
	      for (vector<AliHLTGlobalBarrelTrack>::iterator element=tracks.begin();
		   element!=tracks.end(); element++) 
		{
		  if (CheckCondition(&(*element), fSolenoidBz)) numberOfTracks++;
		}
	    } 
	  else if (iResult<0) 
	    {
	      HLTError("can not extract tracks from data block of type %s (specification %08x) of size %d: error %d", 
		       DataType2Text(pBlock->fDataType).c_str(), pBlock->fSpecification, pBlock->fSize, iResult);
	    }
	}
    }

  bool condition=false;
  description="Geometrical conditions not matched";
  AliHLTReadoutList readout;

  if (numberOfTracks>=fMinTracks) 
    {
      condition=true;
      description=fTriggerDecisionPars->GetDescription();
      readout=fTriggerDecisionPars->GetReadoutListParameter();
      HLTDebug("Geometrical acceptance trigger %s triggered", fTriggerDecisionPars->GetTriggerName().Data());
    }

  AliHLTTriggerDecision decision(
				 condition,
				 fTriggerDecisionPars->GetTriggerName().Data(),
				 AliHLTTriggerDomain(readout),
				 description.Data()
				 );
  TriggerEvent(&decision, kAliHLTDataTypeTObject|kAliHLTDataOriginOut);

  return iResult;

}


template<class T>
bool AliHLTTriggerBarrelGeomMultiplicity::CheckCondition(T* track, float b)
{

  bool ret = false;

  // see header file for class documentation
  if (!track) return false;

  ret = IsInDetectors(track, b);

  return ret;

}

template<class T>
bool AliHLTTriggerBarrelGeomMultiplicity::IsInDetectors(T* track, float b)
{
  // See header file for class documentation  
  for(Int_t i = 0; i < fDetectorArray->GetEntries(); i++)
    {
      AliHLTTriggerDetectorGeom *det = static_cast<AliHLTTriggerDetectorGeom*>(fDetectorArray->At(i));

      Double_t trackPoint[3];
      Double_t normVector[3];

      det->GetInitialPoint(trackPoint);
      det->GetNormVector(normVector);

      bool ret = track->Intersect(trackPoint, normVector, b);

      if(ret)
	{
	  if(det->IsInDetector(trackPoint)) return true;
	}
    }
  return false;
}

int AliHLTTriggerBarrelGeomMultiplicity::DoInit(int argc, const char** argv)
{
  // see header file for class documentation

  // first configure the default
  int iResult=0;

  // Matthias 05.04.2011 code audit
  // looks like somebody has to commission this component
  HLTWarning("this component is not tested and needs most likely a major revision!");

  if (iResult>=0 && argc>0)
    iResult=ConfigureFromArgumentString(argc, argv);

  if (!fTriggerDecisionPars) {
    HLTError("decision parameter not initialized");
    iResult=-ENODEV;
  }
  fSolenoidBz=GetBz();

  return iResult;
}

int AliHLTTriggerBarrelGeomMultiplicity::DoDeinit()
 {
  // see header file for class documentation
   if (fTriggerName) delete fTriggerName;
   fTriggerName=NULL;
  return 0;
}

int AliHLTTriggerBarrelGeomMultiplicity::ReadPreprocessorValues(const char* /*modules*/)
{
    // see header file for function documentation

  // nothing to do for the moment
  return 0;
}

int AliHLTTriggerBarrelGeomMultiplicity::GetDetectorGeomsFromCDBObject(const char *cdbEntry, const char* chainId)
{
    // see header file for function documentation
  int nDetectorGeoms=0;

  if(fDetectorArray)
    {
      fDetectorArray->Clear();
    }
  else
    {
      fDetectorArray = new TObjArray();
    }
  
  const char *path = cdbEntry;

  if(!path) path = fOCDBEntry;
  
  if(path)
    {
      //     const char* chainId=GetChainId();
      HLTInfo("configure from entry %s, chain id %s", path, (chainId!=NULL && chainId[0]!=0)?chainId:"<none>");
      AliCDBEntry *pEntry = AliCDBManager::Instance()->Get(path/*,GetRunNo()*/);
      if (pEntry) 
	{
	  TObjArray* pArr=dynamic_cast<TObjArray*>(pEntry->GetObject());
	  if (pArr) 
	    {

	      for(int i = 0; i < pArr->GetEntries(); i++)
		{
		  if(!strcmp(pArr->At(i)->ClassName(), "AliHLTTriggerDecisionParameters"))
		    {
		      fTriggerDecisionPars = dynamic_cast<AliHLTTriggerDecisionParameters*>(pArr->At(i));
		    }
		  else if(pArr->At(i)->InheritsFrom("AliHLTTriggerDetectorGeom"))
		    {
		      fDetectorArray->AddLast(dynamic_cast<AliHLTTriggerDetectorGeom*>(pArr->At(i)));
		      nDetectorGeoms++;
		      HLTDebug("received detector geometry of type %s", pArr->At(i)->ClassName());
		    }
		  else
		    {
		      HLTWarning("Unknown object of type %s in configuration object", pArr->At(i)->ClassName());
		    }
		}
	    } 
	  else 
	    {
	      HLTError("configuration object \"%s\" has wrong type, required TObjArray", path);
	      nDetectorGeoms=-EINVAL;
	    }
	}
      else 
	{
	  HLTError("can not fetch object \"%s\" from OCDB", path);
	  nDetectorGeoms=-ENOENT;
	}
    }

  HLTInfo("received %d detector geometries", nDetectorGeoms);

  return nDetectorGeoms;
}

int AliHLTTriggerBarrelGeomMultiplicity::GetDetectorGeomsFromFile(const char *filename)
{
    // see header file for function documentation
  int nDetectorGeoms=0;

  if(fDetectorArray)
    {
      fDetectorArray->Clear();
    }
  else
    {
      fDetectorArray = new TObjArray();
    }
  

  if (filename) 
    {
      TFile *geomfile = TFile::Open(filename, "READ");
      
      if(geomfile)
	{
	  HLTInfo("configure from file \"%s\"", filename);
	  TObjArray* pArr=dynamic_cast<TObjArray*>(geomfile->Get("GeomConf"));
	  if (pArr) 
	    {

	      for(int i = 0; i < pArr->GetEntries(); i++)
		{
		  if(!strcmp(pArr->At(i)->ClassName(), "AliHLTTriggerDecisionParameters"))
		    {
		      fTriggerDecisionPars = dynamic_cast<AliHLTTriggerDecisionParameters*>(pArr->At(i));
		    }
		  else if(pArr->At(i)->InheritsFrom("AliHLTTriggerDetectorGeom"))
		    {
		      fDetectorArray->AddLast(dynamic_cast<AliHLTTriggerDetectorGeom*>(pArr->At(i)));
		      nDetectorGeoms++;
		      HLTDebug("received detector geometry of type %s", pArr->At(i)->ClassName());
		    }
		  else
		    {
		      HLTWarning("Unknown object of type %s in configuration object", pArr->At(i)->ClassName());
		    }
		}
	    } 
	  else 
	    {
	      HLTError("configuration object has wrong type, required TObjArray");
	      nDetectorGeoms=-EINVAL;
	    }
	  } 
      else 
	{
	  HLTError("could not open file \"%s\"", filename);
	  nDetectorGeoms=-ENOENT;
	}
    }
  else
    {
      HLTError("ROOT file name not specified");
    }
  HLTInfo("received %d detector geometries", nDetectorGeoms);

  return nDetectorGeoms;
}

int AliHLTTriggerBarrelGeomMultiplicity::ScanConfigurationArgument(int argc, const char** argv)
{
  // See header file for class documentation
  if (argc<=0) return 0;
  int i=0;
  TString argument=argv[i];

  if (argument.CompareTo("-geomfile")==0) 
    {
      if (++i>=argc) return -EPROTO;
    
      GetDetectorGeomsFromFile(argv[i]);
    
      return 2;
    }    

  if (argument.CompareTo("-triggername")==0) 
    {
      if (++i>=argc || argv[i]==NULL) return -EPROTO;

      int namelen=strlen(argv[i])+1;
      fTriggerName = new char[namelen];
      if (!fTriggerName) return -ENOMEM;
      snprintf(fTriggerName, namelen, "%s", argv[i]);
      
      fOCDBEntry = fTriggerName;

      return 2;
  }    
  return 0;
}

