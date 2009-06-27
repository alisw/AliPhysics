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
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "TTree.h"
#include "TList.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTGlobalEsdConverterComponent)

AliHLTGlobalEsdConverterComponent::AliHLTGlobalEsdConverterComponent()
  : AliHLTProcessor()
  , fESD(NULL)
  , fSolenoidBz(5)
  , fWriteTree(0)

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
  return kAliHLTDataTypeESDObject|kAliHLTDataOriginHLT;
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

  if ((iResult=ProcessBlocks(pTree, pESD))>0) {
    // TODO: set the specification correctly
    if (pTree) {
      // the esd structure is written to the user info and is
      // needed in te ReadFromTree method to read all objects correctly
      pTree->GetUserInfo()->Add(pESD);
      pESD->WriteToTree(pTree);
      iResult=PushBack(pTree, kAliHLTDataTypeESDTree|kAliHLTDataOriginHLT, 0);
    } else {
      iResult=PushBack(pESD, kAliHLTDataTypeESDObject|kAliHLTDataOriginHLT, 0);
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
  std::vector<int> trackIdESD2TPCmap; // map esd index -> tpc index

  for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(kAliHLTDataTypeTrack|kAliHLTDataOriginTPC);
       pBlock!=NULL; pBlock=GetNextInputBlock()) {
    vector<AliHLTGlobalBarrelTrack> tracks;
    if ((iResult=AliHLTGlobalBarrelTrack::ConvertTrackDataArray(reinterpret_cast<const AliHLTTracksData*>(pBlock->fPtr), pBlock->fSize, tracks))>0) {
      for (vector<AliHLTGlobalBarrelTrack>::iterator element=tracks.begin();
	   element!=tracks.end(); element++) {
	Float_t points[4] = {
	  element->GetX(),
	  element->GetY(),
	  element->GetLastPointX(),
	  element->GetLastPointY()
	};
	AliESDtrack iotrack;
	iotrack.UpdateTrackParams(&(*element),AliESDtrack::kTPCin);
	iotrack.SetTPCPoints(points);

	pESD->AddTrack(&iotrack);
      }
      iAddedDataBlocks++;
    } else if (iResult<0) {
      HLTError("can not extract tracks from data block of type %s (specification %08x) of size $d", 
	       DataType2Text(pBlock->fDataType).c_str(), pBlock->fSpecification, pBlock->fSize);
    }
  }

      // ITS updated tracks
      /*
      int nESDTracks = pESD->GetNumberOfTracks();

      // create map of tpc->esd track indices

      int *trackIdTPC2ESDmap = new int[ nTPCTracks ];
      {
	for( int i=0; i<nTPCTracks; i++ ) trackIdTPC2ESDmap[i] = -1;
	for( unsigned int i=0; i<trackIdESD2TPCmap.size(); i++ ){
	  int tpcId = trackIdESD2TPCmap[i];
	  if( tpcId>=0 && tpcId<nTPCTracks ) trackIdTPC2ESDmap[tpcId] = i;
	}
      }

      for (int ndx=0; ndx<nBlocks && iResult>=0; ndx++) {
	iter = blocks+ndx;
	if(iter->fDataType == fgkITSTracksDataType ) {
	  AliHLTITSTrackDataHeader *inPtr = reinterpret_cast<AliHLTITSTrackDataHeader*>( iter->fPtr );
	  int nTracks = inPtr->fTrackletCnt;	  
	  for( int itr=0; itr<nTracks; itr++ ){
	    AliHLTITSTrackData &tr = inPtr->fTracks[itr];
	    AliHLTITSTrack tITS( tr.fTrackParam );
	    int ncl=0;
	    for( int il=0; il<6; il++ ){
	      if( tr.fClusterIds[il]>=0 ) tITS.SetClusterIndex(ncl++,  tr.fClusterIds[il]);
	    }
	    tITS.SetNumberOfClusters( ncl );
	    int tpcId = tr.fTPCId;
	    if( tpcId<0 || tpcId>=nTPCTracks ) continue;
	    int esdID = trackIdTPC2ESDmap[tpcId];
	    if( esdID<0 || esdID>=nESDTracks ) continue;
	    AliESDtrack *tESD = pESD->GetTrack( esdID );
	    if( tESD ) tESD->UpdateTrackParams( &tITS, AliESDtrack::kITSin );
	  }
	}
      }   
      delete[] trackIdTPC2ESDmap;
      
      // primary vertex & V0's 

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
