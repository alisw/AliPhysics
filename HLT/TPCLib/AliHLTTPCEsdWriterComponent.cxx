// @(#) $Id$

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

/** @file   AliHLTTPCEsdWriterComponent.cxx
    @author Matthias Richter
    @date   
    @brief  Writer component to store tracks of the HLT TPC conformal
            mapping tracker in the AliESD format
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include <cassert>
#include "AliHLTTPCEsdWriterComponent.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "TTree.h"
#include "TList.h"
#include "AliHLTTPCTrack.h"
#include "AliHLTTPCTrackArray.h"
#include "AliHLTTPCTrackletDataFormat.h"
#include "AliHLTTPCDefinitions.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTExternalTrackParam.h"
#include "AliHLTGlobalBarrelTrack.h"
#include "AliHLTTrackMCLabel.h"
#include "TGeoGlobalMagField.h"

#include <vector>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCEsdWriterComponent)

AliHLTTPCEsdWriterComponent::AliHLTTPCEsdWriterComponent()
  :
  fSolenoidBz(0)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTPCEsdWriterComponent::~AliHLTTPCEsdWriterComponent()
{
  // see header file for class documentation
}

AliHLTTPCEsdWriterComponent::AliWriter::AliWriter()
  :
  fTree(NULL),
  fESD(NULL),
  fBase(new AliHLTTPCEsdWriterComponent)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTPCEsdWriterComponent::AliWriter::~AliWriter()
{
  // see header file for class documentation
  if (fBase) delete fBase;
  fBase=NULL;
}

void AliHLTTPCEsdWriterComponent::AliWriter::GetInputDataTypes(AliHLTComponentDataTypeList& list)
{
  // see header file for class documentation
  list.push_back(AliHLTTPCDefinitions::fgkTrackSegmentsDataType);
  list.push_back(AliHLTTPCDefinitions::fgkTracksDataType);
  list.push_back( kAliHLTDataTypeTrack|kAliHLTDataOriginTPC );
}

int AliHLTTPCEsdWriterComponent::AliWriter::InitWriter()
{
  // see header file for class documentation
  int iResult=0;
  fESD = new AliESDEvent;
  if (fESD) {
    fESD->CreateStdContent();
    fTree = new TTree("esdTree", "Tree with HLT ESD objects");
    if (fTree) {
      fTree->SetDirectory(0);
      fESD->WriteToTree(fTree);
    }
  }
  if (fTree==NULL) {
    iResult=-ENOMEM;
  }

  if (iResult>=0) {
    iResult=fBase->Reconfigure(NULL, NULL);
  }

  return iResult;
}

int AliHLTTPCEsdWriterComponent::AliWriter::CloseWriter()
{
  // see header file for class documentation
  int iResult=0;
  if (fTree) {
    // the esd structure is written to the user info and is
    // needed in te ReadFromTree method to read all objects correctly
    if (fESD) fTree->GetUserInfo()->Add(fESD);
    WriteObject(kAliHLTVoidEventID, fTree);
    fTree->GetUserInfo()->Clear();
    TTree* pTree=fTree;
    fTree=NULL;
    delete pTree;
  } else {
    HLTWarning("not initialized");
  }

  if (fESD) {
    delete fESD;
  }
  iResult=AliHLTRootFileWriterComponent::CloseWriter();
  return iResult;
}

int AliHLTTPCEsdWriterComponent::AliWriter::DumpEvent( const AliHLTComponentEventData& evtData,
					    const AliHLTComponentBlockData* blocks, 
					    AliHLTComponentTriggerData& /*trigData*/ )
{
  // see header file for class documentation
  int iResult=0;
  TTree* pTree=fTree;
  assert(fBase);
  if (pTree && fBase) {
    if (fESD) {
      AliESDEvent* pESD=fESD;

      pESD->Reset(); 
      pESD->SetMagneticField(GetBz());
      iResult=fBase->ProcessBlocks(pTree, pESD, blocks, (int)evtData.fBlockCnt);

    } else {
      iResult=-ENOMEM;
    }
  }
  return iResult;
}

int AliHLTTPCEsdWriterComponent::AliWriter::ScanArgument(int argc, const char** argv)
{
  // see header file for class documentation
  int iResult=AliHLTRootFileWriterComponent::ScanArgument(argc, argv);
  return iResult;
}

int AliHLTTPCEsdWriterComponent::ProcessBlocks(TTree* pTree, AliESDEvent* pESD,
					       const AliHLTComponentBlockData* blocks,
					       int nBlocks, int* pMinSlice,
					       int* pMaxSlice)
{
  // see header file for class documentation

  int iResult=0;
  int iAddedDataBlocks=0;
  
  if (pESD && blocks) {
      const AliHLTComponentBlockData* iter = NULL;
      int bIsTrackSegs=0;

      // first read MC information (if present)

      std::map<int,int> mcLabels;

      for (int ndx=0; ndx<nBlocks && iResult>=0; ndx++) {
	iter = blocks+ndx;
	if(iter->fDataType == (kAliHLTDataTypeTrackMC|kAliHLTDataOriginTPC) ) {
	  AliHLTTrackMCData* dataPtr = ( AliHLTTrackMCData* )( iter->fPtr );
	  for( unsigned int il=0; il<dataPtr->fCount; il++ ){
	    AliHLTTrackMCLabel &lab = dataPtr->fLabels[il];
	    mcLabels[lab.fTrackID] = lab.fMCLabel;
	  }
	}
      }


      // do the conversion of tracks

      for (int ndx=0; ndx<nBlocks && iResult>=0; ndx++) {
	iter = blocks+ndx;
	
	if( iter->fDataType == ( kAliHLTDataTypeTrack|kAliHLTDataOriginTPC ) ){	  
	  AliHLTTracksData* dataPtr = ( AliHLTTracksData* ) iter->fPtr;
	  int nTracks = dataPtr->fCount;
	  AliHLTExternalTrackParam* currOutTrack = dataPtr->fTracklets;

	  for( int itr=0; itr<nTracks; itr++ ){
	    AliHLTGlobalBarrelTrack t(*currOutTrack);
	    Float_t points[4] = {currOutTrack->fX, currOutTrack->fY, currOutTrack->fLastX, currOutTrack->fLastY };
	    
	    Int_t mcLabel = -1;
	    if( mcLabels.find(currOutTrack->fTrackID)!=mcLabels.end() )
	      mcLabel = mcLabels[currOutTrack->fTrackID];
	    
	    t.SetLabel( mcLabel );
	    
	    AliESDtrack iotrack;
	    iotrack.UpdateTrackParams( &t,AliESDtrack::kTPCin);
	    iotrack.SetTPCPoints(points);
	    pESD->AddTrack(&iotrack);
	    unsigned int dSize = sizeof( AliHLTExternalTrackParam ) + currOutTrack->fNPoints * sizeof( unsigned int );
	    currOutTrack = ( AliHLTExternalTrackParam* )( (( Byte_t * )currOutTrack) + dSize );
	  }
	  iAddedDataBlocks++;
	}

	
	if ( (bIsTrackSegs=(iter->fDataType == AliHLTTPCDefinitions::fgkTrackSegmentsDataType))==1 ||
	     iter->fDataType == AliHLTTPCDefinitions::fgkTracksDataType ) {
	  Int_t minslice=AliHLTTPCDefinitions::GetMinSliceNr(iter->fSpecification);
	  Int_t maxslice=AliHLTTPCDefinitions::GetMaxSliceNr(iter->fSpecification);
	  if (bIsTrackSegs==0) {
	    // slice parameter and data specification ignored, tracks already in global coordinates
	    minslice=-1;
	    maxslice=-1;
	    if (pMinSlice) *pMinSlice=0;
	    if (pMaxSlice) *pMaxSlice=AliHLTTPCTransform::GetNSlice()-1;
	  } else {
	    if (pMinSlice && (*pMinSlice==-1 || *pMinSlice>minslice)) *pMinSlice=minslice;
	    if (pMaxSlice && (*pMaxSlice==-1 || *pMaxSlice<maxslice)) *pMaxSlice=maxslice;
	  }
	  //HLTDebug("dataspec %#x minslice %d", iter->fSpecification, minslice);
	  if (minslice >=-1 && minslice<AliHLTTPCTransform::GetNSlice()) {
	    if (minslice!=maxslice) {
	      HLTWarning("data from multiple sectors in one block: "
			 "possible mismatch in treatment of local coordinate system");
	    }
	    AliHLTTPCTrackArray tracks;
	    AliHLTTPCTrackletData* inPtr = (AliHLTTPCTrackletData*) iter->fPtr;	    
	    HLTDebug("reading block %d (slice %d): %d tracklets", ndx, minslice, inPtr->fTrackletCnt);
	    if ((iResult=tracks.FillTracksChecked(inPtr->fTracklets, inPtr->fTrackletCnt, iter->fSize, minslice, 0/*don't rotate*/))>=0) {
	      if ((iResult=Tracks2ESD(&tracks, pESD ))>=0) {
		iAddedDataBlocks++;
	      }
	    }
	  } else {
	    HLTError("invalid sector number");
	    iResult=-EBADF;
	  }
	}      
      }

      if (iAddedDataBlocks>0 && pTree) {
	pTree->Fill();
      }
  
  } else {
    iResult=-EINVAL;
  }
  if (iResult>=0) iResult=iAddedDataBlocks;
  return iResult;
}


int AliHLTTPCEsdWriterComponent::Tracks2ESD(AliHLTTPCTrackArray* pTracks, AliESDEvent* pESD )
{
  // see header file for class documentation
  int iResult=0;

  if (pTracks && pESD) {    
 
    for (int i=0; i<pTracks->GetNTracks() && iResult>=0; i++) {
      AliHLTTPCTrack* pTrack=(*pTracks)[i];
      if (pTrack) {
	
	if( pTrack->Convert2AliKalmanTrack() ){	  
	  HLTError("conversion to AliKalmanTrack failed for track %d of %d", i, pTracks->GetNTracks());	
	  continue;
	}

	Float_t points[4] = {pTrack->GetFirstPointX(), pTrack->GetFirstPointY(), pTrack->GetLastPointX(), pTrack->GetLastPointY() };

	if(pTrack->GetSector() == -1){ // Set first and last points for global tracks
	  Double_t s = TMath::Sin( pTrack->GetAlpha() );
	  Double_t c = TMath::Cos( pTrack->GetAlpha() );
	  points[0] =  pTrack->GetFirstPointX()*c + pTrack->GetFirstPointY()*s;
	  points[1] = -pTrack->GetFirstPointX()*s + pTrack->GetFirstPointY()*c;	  
	  points[2] =  pTrack->GetLastPointX() *c + pTrack->GetLastPointY() *s;
	  points[3] = -pTrack->GetLastPointX() *s + pTrack->GetLastPointY() *c;	  
	}

	AliESDtrack iotrack;
	iotrack.UpdateTrackParams(pTrack,AliESDtrack::kTPCin);
	iotrack.SetTPCPoints(points);

	pESD->AddTrack(&iotrack);
      } else {
	HLTError("internal mismatch in array");
	iResult=-EFAULT;
      }
    }
    
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTTPCEsdWriterComponent::Configure(const char* arguments)
{
  // see header file for class documentation
  int iResult=0;
  if (!arguments) return iResult;

  // Check field
  if (!TGeoGlobalMagField::Instance()) {
    HLTError("magnetic field not initialized, please set up TGeoGlobalMagField and AliMagF");
    return -ENODEV;
  }

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
	HLTWarning("parameter -solenoidBz deprecated, magnetic field handled by global AliMagF object and TGeoGlobalMagField");
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

int AliHLTTPCEsdWriterComponent::Reconfigure(const char* cdbEntry, const char* chainId)
{
  // see header file for class documentation
  int iResult=0;
  const char* path=NULL;
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

AliHLTTPCEsdWriterComponent::AliConverter::AliConverter()
  :
  fESD(NULL),
  fBase(new AliHLTTPCEsdWriterComponent),
  fWriteTree(0)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTPCEsdWriterComponent::AliConverter::~AliConverter()
{
  // see header file for class documentation
  if (fBase) delete fBase;
  fBase=NULL;

  if (fESD) delete fESD;
  fESD=NULL;
}

void AliHLTTPCEsdWriterComponent::AliConverter::GetInputDataTypes(AliHLTComponentDataTypeList& list)
{
  // see header file for class documentation
  list.push_back(AliHLTTPCDefinitions::fgkTrackSegmentsDataType);
  list.push_back(AliHLTTPCDefinitions::fgkTracksDataType);
  list.push_back( kAliHLTDataTypeTrack|kAliHLTDataOriginTPC );
}

AliHLTComponentDataType AliHLTTPCEsdWriterComponent::AliConverter::GetOutputDataType()
{
  // see header file for class documentation
  return kAliHLTDataTypeESDTree;
}

void AliHLTTPCEsdWriterComponent::AliConverter::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
{
  // see header file for class documentation
  constBase=2000000;
  inputMultiplier=10.0;
}

int AliHLTTPCEsdWriterComponent::AliConverter::DoInit(int argc, const char** argv)
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

      // -solenoidBz
    } else if (argument.CompareTo("-solenoidBz")==0) {
      HLTWarning("parameter -solenoidBz deprecated, magnetic field handled by global AliMagF object and TGeoGlobalMagField");
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
    iResult=fBase->Reconfigure(NULL, NULL);
  }

  return iResult;
}

int AliHLTTPCEsdWriterComponent::AliConverter::DoDeinit()
{
  // see header file for class documentation
  return 0;
}

int AliHLTTPCEsdWriterComponent::AliConverter::DoEvent(const AliHLTComponentEventData& evtData, 
						       const AliHLTComponentBlockData* blocks, 
						       AliHLTComponentTriggerData& /*trigData*/,
						       AliHLTUInt8_t* /*outputPtr*/, 
						       AliHLTUInt32_t& size,
						       AliHLTComponentBlockDataList& /*outputBlocks*/ )
{
  // see header file for class documentation
  int iResult=0;
  // no direct writing to the output buffer
  size=0;

  assert(fBase);
  if (!fESD) {
    fESD = new AliESDEvent;
    if (fESD) {
      fESD->CreateStdContent();
    } else {
      iResult=-ENOMEM;
    }
  }

  AliESDEvent* pESD = fESD;

  if (pESD && fBase) {
  
    TTree* pTree = NULL;
    // TODO: Matthias 06.12.2007
    // Tried to write the ESD directly instead to a tree, but this did not work
    // out. Information in the ESD is different, needs investigation.
    
    if (fWriteTree)
      pTree = new TTree("esdTree", "Tree with HLT ESD objects");
 
    if (pTree) {
      pTree->SetDirectory(0);
    }

    pESD->Reset(); 
    pESD->SetMagneticField(GetBz());
    if ((iResult=fBase->ProcessBlocks(pTree, pESD, blocks, (int)evtData.fBlockCnt))>0) {
      // TODO: set the specification correctly
      if (pTree) {
	// the esd structure is written to the user info and is
	// needed in te ReadFromTree method to read all objects correctly
	pTree->GetUserInfo()->Add(pESD);
	pESD->WriteToTree(pTree);
	iResult=PushBack(pTree, kAliHLTDataTypeESDTree|kAliHLTDataOriginTPC, 0);
      } else {
	iResult=PushBack(pESD, kAliHLTDataTypeESDObject|kAliHLTDataOriginTPC, 0);
      }
    }
    if (pTree) {
      // clear user info list to prevent objects from being deleted
      pTree->GetUserInfo()->Clear();
      delete pTree;
    }
  }
  return iResult;
}

