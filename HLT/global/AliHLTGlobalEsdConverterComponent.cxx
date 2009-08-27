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

/** @file   AliHLTGlobalEsdConverterComponent.cxx
    @author Matthias Richter
    @date   
    @brief  Global ESD converter component.
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include <cassert>
#include "AliHLTGlobalEsdConverterComponent.h"
#include "AliHLTGlobalBarrelTrack.h"
#include "AliHLTExternalTrackParam.h"
#include "AliHLTTrackMCLabel.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliPID.h"
#include "TTree.h"
#include "TList.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTGlobalEsdConverterComponent)

AliHLTGlobalEsdConverterComponent::AliHLTGlobalEsdConverterComponent()
  : AliHLTProcessor()
  , fWriteTree(0)
  , fVerbosity(0)
  , fESD(NULL)
  , fSolenoidBz(5)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTGlobalEsdConverterComponent::~AliHLTGlobalEsdConverterComponent()
{
  // see header file for class documentation
  if (fESD) delete fESD;
  fESD=NULL;
}

int AliHLTGlobalEsdConverterComponent::Configure(const char* arguments)
{
  // see header file for class documentation
  int iResult=0;
  if (!arguments) return iResult;

  TString allArgs=arguments;
  TString argument;
  int bMissingParam=0;

  TObjArray* pTokens=allArgs.Tokenize(" ");
  if (pTokens) {
    for (int i=0; i<pTokens->GetEntries() && iResult>=0; i++) {
      argument=((TObjString*)pTokens->At(i))->GetString();
      if (argument.IsNull()) continue;
      
      if (argument.CompareTo("-solenoidBz")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	HLTInfo("Magnetic Field set to: %s", ((TObjString*)pTokens->At(i))->GetString().Data());
	fSolenoidBz=((TObjString*)pTokens->At(i))->GetString().Atof();
	continue;
      } else {
	HLTError("unknown argument %s", argument.Data());
	iResult=-EINVAL;
	break;
      }
    }
    delete pTokens;
  }
  if (bMissingParam) {
    HLTError("missing parameter for argument %s", argument.Data());
    iResult=-EINVAL;
  }

  return iResult;
}

int AliHLTGlobalEsdConverterComponent::Reconfigure(const char* cdbEntry, const char* chainId)
{
  // see header file for class documentation
  int iResult=0;
  const char* path=kAliHLTCDBSolenoidBz;
  const char* defaultNotify="";
  if (cdbEntry) {
    path=cdbEntry;
    defaultNotify=" (default)";
  }
  if (path) {
    HLTInfo("reconfigure from entry %s%s, chain id %s", path, defaultNotify,(chainId!=NULL && chainId[0]!=0)?chainId:"<none>");
    AliCDBEntry *pEntry = AliCDBManager::Instance()->Get(path/*,GetRunNo()*/);
    if (pEntry) {
      TObjString* pString=dynamic_cast<TObjString*>(pEntry->GetObject());
      if (pString) {
	HLTInfo("received configuration object string: \'%s\'", pString->GetString().Data());
	iResult=Configure(pString->GetString().Data());
      } else {
	HLTError("configuration object \"%s\" has wrong type, required TObjString", path);
      }
    } else {
      HLTError("can not fetch object \"%s\" from CDB", path);
    }
  }
  
  return iResult;
}

void AliHLTGlobalEsdConverterComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list)
{
  // see header file for class documentation
  list.push_back(kAliHLTDataTypeTrack);
  list.push_back(kAliHLTDataTypeTrackMC);
}

AliHLTComponentDataType AliHLTGlobalEsdConverterComponent::GetOutputDataType()
{
  // see header file for class documentation
  return kAliHLTDataTypeESDObject|kAliHLTDataOriginOut;
}

void AliHLTGlobalEsdConverterComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
{
  // see header file for class documentation
  constBase=2000000;
  inputMultiplier=10.0;
}

int AliHLTGlobalEsdConverterComponent::DoInit(int argc, const char** argv)
{
  // see header file for class documentation
  int iResult=0;
  TString argument="";
  int bMissingParam=0;
  for (int i=0; i<argc && iResult>=0; i++) {
    argument=argv[i];
    if (argument.IsNull()) continue;

    // -notree
    if (argument.CompareTo("-notree")==0) {
      fWriteTree=0;

      // -tree
    } else if (argument.CompareTo("-tree")==0) {
      fWriteTree=1;

    } else {
      HLTError("unknown argument %s", argument.Data());
      break;
    }
  }
  if (bMissingParam) {
    HLTError("missing parameter for argument %s", argument.Data());
    iResult=-EINVAL;
  }

  if (iResult>=0) {
    iResult=Reconfigure(NULL, NULL);
  }

  if (iResult>=0) {
    fESD = new AliESDEvent;
    if (fESD) {
      fESD->CreateStdContent();
    } else {
      iResult=-ENOMEM;
    }
  }

  return iResult;
}

int AliHLTGlobalEsdConverterComponent::DoDeinit()
{
  // see header file for class documentation
  if (fESD) delete fESD;
  fESD=NULL;

  return 0;
}

int AliHLTGlobalEsdConverterComponent::DoEvent(const AliHLTComponentEventData& /*evtData*/, 
					       AliHLTComponentTriggerData& /*trigData*/)
{
  // see header file for class documentation
  int iResult=0;
  if (!fESD) return -ENODEV;

  AliESDEvent* pESD = fESD;

  pESD->Reset(); 
  pESD->SetMagneticField(fSolenoidBz);

  TTree* pTree = NULL;
  if (fWriteTree)
    pTree = new TTree("esdTree", "Tree with HLT ESD objects");
 
  if (pTree) {
    pTree->SetDirectory(0);
  }

  if ((iResult=ProcessBlocks(pTree, pESD))>=0) {
    // TODO: set the specification correctly
    if (pTree) {
      // the esd structure is written to the user info and is
      // needed in te ReadFromTree method to read all objects correctly
      pTree->GetUserInfo()->Add(pESD);
      pESD->WriteToTree(pTree);
      iResult=PushBack(pTree, kAliHLTDataTypeESDTree|kAliHLTDataOriginOut, 0);
    } else {
      iResult=PushBack(pESD, kAliHLTDataTypeESDObject|kAliHLTDataOriginOut, 0);
    }
  }
  if (pTree) {
    // clear user info list to prevent objects from being deleted
    pTree->GetUserInfo()->Clear();
    delete pTree;
  }
  return iResult;
}

int AliHLTGlobalEsdConverterComponent::ProcessBlocks(TTree* pTree, AliESDEvent* pESD)
{
  // see header file for class documentation

  int iResult=0;
  int iAddedDataBlocks=0;

  // Barrel tracking

  // in the first attempt this component reads the TPC tracks and updates in the
  // second step from the ITS tracks


  // first read MC information (if present)
  std::map<int,int> mcLabels;

  for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(kAliHLTDataTypeTrackMC|kAliHLTDataOriginTPC);
       pBlock!=NULL; pBlock=GetNextInputBlock()) {
    AliHLTTrackMCData* dataPtr = reinterpret_cast<AliHLTTrackMCData*>( pBlock->fPtr );
    if (sizeof(AliHLTTrackMCData)+dataPtr->fCount*sizeof(AliHLTTrackMCLabel)==pBlock->fSize) {
      for( unsigned int il=0; il<dataPtr->fCount; il++ ){
	AliHLTTrackMCLabel &lab = dataPtr->fLabels[il];
	mcLabels[lab.fTrackID] = lab.fMCLabel;
      }
    } else {
      HLTWarning("data mismatch in block %s (0x%08x): count %d, size %d -> ignoring track MC information", 
		 DataType2Text(pBlock->fDataType).c_str(), pBlock->fSpecification, 
		 dataPtr->fCount, pBlock->fSize);
    }
  }

  // convert the TPC tracks to ESD tracks
  for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(kAliHLTDataTypeTrack|kAliHLTDataOriginTPC);
       pBlock!=NULL; pBlock=GetNextInputBlock()) {
    vector<AliHLTGlobalBarrelTrack> tracks;
    if ((iResult=AliHLTGlobalBarrelTrack::ConvertTrackDataArray(reinterpret_cast<const AliHLTTracksData*>(pBlock->fPtr), pBlock->fSize, tracks))>=0) {
      for (vector<AliHLTGlobalBarrelTrack>::iterator element=tracks.begin();
	   element!=tracks.end(); element++) {
	Float_t points[4] = {
	  element->GetX(),
	  element->GetY(),
	  element->GetLastPointX(),
	  element->GetLastPointY()
	};

	Int_t mcLabel = -1;
	if( mcLabels.find(element->TrackID())!=mcLabels.end() )
	  mcLabel = mcLabels[element->TrackID()];
	element->SetLabel( mcLabel );

	AliESDtrack iotrack;
	iotrack.UpdateTrackParams(&(*element),AliESDtrack::kTPCin);
	iotrack.SetTPCPoints(points);

	pESD->AddTrack(&iotrack);
	if (fVerbosity>0) element->Print();
      }
      HLTInfo("converted %d track(s) to AliESDtrack and added to ESD", tracks.size());
      iAddedDataBlocks++;
    } else if (iResult<0) {
      HLTError("can not extract tracks from data block of type %s (specification %08x) of size %d: error %d", 
	       DataType2Text(pBlock->fDataType).c_str(), pBlock->fSpecification, pBlock->fSize, iResult);
    }
  }

  // now update ESD tracks with the ITS info
  for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(kAliHLTDataTypeTrack|kAliHLTDataOriginITS);
       pBlock!=NULL; pBlock=GetNextInputBlock()) {
    vector<AliHLTGlobalBarrelTrack> tracks;
    if ((iResult=AliHLTGlobalBarrelTrack::ConvertTrackDataArray(reinterpret_cast<const AliHLTTracksData*>(pBlock->fPtr), pBlock->fSize, tracks))>0) {
      for (vector<AliHLTGlobalBarrelTrack>::iterator element=tracks.begin();
	   element!=tracks.end(); element++) {
	int ncl=0;
	const UInt_t* pointsArray=element->GetPoints();
	for( unsigned il=0; il<element->GetNumberOfPoints(); il++ ){
	  // TODO: check what needs to be done with the clusters
	  if( pointsArray[il]<~(UInt_t)0 ) {/*tITS.SetClusterIndex(ncl,  tr.fClusterIds[il]);*/}
	  ncl++;
	}
	//tITS.SetNumberOfClusters( ncl );
	int tpcID=element->TrackID();
	// the ITS tracker assigns the TPC track used as seed for a certain track to
	// the trackID
	if( tpcID<0 || tpcID>=pESD->GetNumberOfTracks()) continue;

	AliESDtrack *tESD = pESD->GetTrack( tpcID );
	if( tESD ) tESD->UpdateTrackParams( &(*element), AliESDtrack::kITSin );
      }
    }
  }

  // convert the HLT TRD tracks to ESD tracks                        
  for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(kAliHLTDataTypeTrack | kAliHLTDataOriginTRD);
       pBlock!=NULL; pBlock=GetNextInputBlock()) {
    vector<AliHLTGlobalBarrelTrack> tracks;
    if ((iResult=AliHLTGlobalBarrelTrack::ConvertTrackDataArray(reinterpret_cast<const AliHLTTracksData*>(pBlock->fPtr), pBlock->fSize, tracks))>0) {
      for (vector<AliHLTGlobalBarrelTrack>::iterator element=tracks.begin();
	   element!=tracks.end(); element++) {
	
	Double_t TRDpid[AliPID::kSPECIES], eProb(0.2), restProb((1-eProb)/(AliPID::kSPECIES-1)); //eprob(element->GetTRDpid...);
	for(Int_t i=0; i<AliPID::kSPECIES; i++){
	  switch(i){
	  case AliPID::kElectron: TRDpid[AliPID::kElectron]=eProb; break;
	  default: TRDpid[i]=restProb; break;
	  }
	}
	
	AliESDtrack iotrack;
	iotrack.UpdateTrackParams(&(*element),AliESDtrack::kTRDin);
	iotrack.SetTRDpid(TRDpid);
	
	pESD->AddTrack(&iotrack);
	if (fVerbosity>0) element->Print();
      }
      HLTInfo("converted %d track(s) to AliESDtrack and added to ESD", tracks.size());
      iAddedDataBlocks++;
    } else if (iResult<0) {
      HLTError("can not extract tracks from data block of type %s (specification %08x) of size %d: error %d", 
	       DataType2Text(pBlock->fDataType).c_str(), pBlock->fSpecification, pBlock->fSize, iResult);
    }
  }
      
      // primary vertex & V0's 
      /*
      //AliHLTVertexer vertexer;
      //vertexer.SetESD( pESD );
      //vertexer.FindPrimaryVertex();
      //vertexer.FindV0s();
      */
      if (iAddedDataBlocks>0 && pTree) {
       pTree->Fill();
      }
  
  if (iResult>=0) iResult=iAddedDataBlocks;
  return iResult;
}
