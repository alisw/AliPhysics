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

//  @file   AliHLTGlobalEsdConverterComponent.cxx
//  @author Matthias Richter
//  @date   
//  @brief  Global ESD converter component.
// 

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
#include "AliHLTCTPData.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDMuonTrack.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliPID.h"
#include "TTree.h"
#include "TList.h"
#include "TClonesArray.h"
#include "AliHLTESDCaloClusterMaker.h"
#include "AliHLTCaloClusterDataStruct.h"
#include "AliHLTCaloClusterReader.h"
#include "AliESDCaloCluster.h"
#include "AliHLTGlobalVertexerComponent.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTGlobalEsdConverterComponent)

AliHLTGlobalEsdConverterComponent::AliHLTGlobalEsdConverterComponent()
  : AliHLTProcessor()
  , fWriteTree(0)
  , fVerbosity(0)
  , fESD(NULL)
  , fSolenoidBz(-5.00668)
  , fBenchmark("EsdConverter")
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
	HLTWarning("argument -solenoidBz is deprecated, magnetic field set up globally (%f)", GetBz());
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

void AliHLTGlobalEsdConverterComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list)
{
  // see header file for class documentation
  list.push_back(kAliHLTDataTypeTrack);
  list.push_back(kAliHLTDataTypeTrackMC);
  list.push_back(kAliHLTDataTypeCaloCluster);
  list.push_back(kAliHLTDataTypedEdx );
  list.push_back(kAliHLTDataTypeESDVertex );
  list.push_back(kAliHLTDataTypeESDObject);
  list.push_back(kAliHLTDataTypeTObject);
  list.push_back(kAliHLTDataTypeGlobalVertexer);
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

  // default list of skiped ESD objects
  TString skipObjects=
    // "AliESDRun,"
    // "AliESDHeader,"
    "AliESDZDC,"
    "AliESDFMD,"
    // "AliESDVZERO,"
    // "AliESDTZERO,"
    // "TPCVertex,"
    // "SPDVertex,"
    // "PrimaryVertex,"
    // "AliMultiplicity,"
    // "PHOSTrigger,"
    // "EMCALTrigger,"
    // "SPDPileupVertices,"
    // "TrkPileupVertices,"
    "Cascades,"
    "Kinks,"
    "AliRawDataErrorLogs,"
    "AliESDACORDE";

  iResult=Reconfigure(NULL, NULL);
  TString allArgs = "";
  for ( int i = 0; i < argc; i++ ) {
    if ( !allArgs.IsNull() ) allArgs += " ";
    allArgs += argv[i];
  }

  TObjArray* pTokens=allArgs.Tokenize(" ");
  if (pTokens) {
    for (int i=0; i<pTokens->GetEntries() && iResult>=0; i++) {
      argument=((TObjString*)pTokens->At(i))->GetString();	
      if (argument.IsNull()) continue;

      // -notree
      if (argument.CompareTo("-notree")==0) {
	fWriteTree=0;
	
	// -tree
      } else if (argument.CompareTo("-tree")==0) {
	fWriteTree=1;
      } else if (argument.CompareTo("-solenoidBz")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	HLTInfo("Magnetic Field set to: %s", ((TObjString*)pTokens->At(i))->GetString().Data());
	HLTWarning("argument '-solenoidBz' is deprecated, solenoid field initiaized from CDB settings");
	continue;
      } else if (argument.Contains("-skipobject=")) {
	argument.ReplaceAll("-skipobject=", "");
	skipObjects=argument;
      } else {
	HLTError("unknown argument %s", argument.Data());
	iResult=-EINVAL;
	break;
      }
    }
  }
  if (bMissingParam) {
    HLTError("missing parameter for argument %s", argument.Data());
    iResult=-EINVAL;
  }

  fSolenoidBz=GetBz();

  if (iResult>=0) {
    fESD = new AliESDEvent;
    if (fESD) {
      fESD->CreateStdContent();

      // remove some of the objects which are not needed
      if (fESD->GetList() && !skipObjects.IsNull()) {
	pTokens=skipObjects.Tokenize(",");
	if (pTokens) {
	  const char* id=NULL;
	  TIter next(pTokens);
	  TObject* pObject=NULL;
	  while ((pObject=next())!=NULL) {
	    id=((TObjString*)pObject)->GetString().Data();
	    TObject* pObj=fESD->GetList()->FindObject(id);
	    if (pObj) {
	      HLTDebug("removing object %s", id);
	      fESD->GetList()->Remove(pObj);
	      delete pObj;
	    } else {
	      HLTWarning("failed to remove object '%s' from ESD", id);
	    }
	  }
	  fESD->GetStdContent();
	  delete pTokens;
	}
      }
    } else {
      iResult=-ENOMEM;
    }

    SetupCTPData();
  }

  fBenchmark.SetTimer(0,"total");

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
					       AliHLTComponentTriggerData& trigData)
{
  // see header file for class documentation
  int iResult=0;
  if (!fESD) return -ENODEV;

  if (IsDataEvent()) fBenchmark.StartNewEvent();
  fBenchmark.Start(0);

  AliESDEvent* pESD = fESD;

  pESD->Reset(); 
  pESD->SetMagneticField(fSolenoidBz);
  pESD->SetRunNumber(GetRunNo());
  pESD->SetPeriodNumber(GetPeriodNumber());
  pESD->SetOrbitNumber(GetOrbitNumber());
  pESD->SetBunchCrossNumber(GetBunchCrossNumber());
  pESD->SetTimeStamp(GetTimeStamp());

  const AliHLTCTPData* pCTPData=CTPData();
  if (pCTPData) {
    AliHLTUInt64_t mask=pCTPData->ActiveTriggers(trigData);
    for (int index=0; index<gkNCTPTriggerClasses; index++) {
      if ((mask&((AliHLTUInt64_t)0x1<<index)) == 0) continue;
      pESD->SetTriggerClass(pCTPData->Name(index), index);
    }
    pESD->SetTriggerMask(mask);
  }

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
    fBenchmark.AddOutput(GetLastObjectSize());
  }
  if (pTree) {
    // clear user info list to prevent objects from being deleted
    pTree->GetUserInfo()->Clear();
    delete pTree;
  }

  fBenchmark.Stop(0);
  HLTInfo( fBenchmark.GetStatistics() );

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
  std::map<int,int> mcLabelsTPC;
  std::map<int,int> mcLabelsITS;

  for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(kAliHLTDataTypeTrackMC|kAliHLTDataOriginTPC);
       pBlock!=NULL; pBlock=GetNextInputBlock()) {

    fBenchmark.AddInput(pBlock->fSize);

    AliHLTTrackMCData* dataPtr = reinterpret_cast<AliHLTTrackMCData*>( pBlock->fPtr );
    if (sizeof(AliHLTTrackMCData)+dataPtr->fCount*sizeof(AliHLTTrackMCLabel)==pBlock->fSize) {
      for( unsigned int il=0; il<dataPtr->fCount; il++ ){
	AliHLTTrackMCLabel &lab = dataPtr->fLabels[il];
	mcLabelsTPC[lab.fTrackID] = lab.fMCLabel;
      }
    } else {
      HLTWarning("data mismatch in block %s (0x%08x): count %d, size %d -> ignoring track MC information", 
		 DataType2Text(pBlock->fDataType).c_str(), pBlock->fSpecification, 
		 dataPtr->fCount, pBlock->fSize);
    }
  }
 
  for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(kAliHLTDataTypeTrackMC|kAliHLTDataOriginITS);
       pBlock!=NULL; pBlock=GetNextInputBlock()) {

    fBenchmark.AddInput(pBlock->fSize);

    AliHLTTrackMCData* dataPtr = reinterpret_cast<AliHLTTrackMCData*>( pBlock->fPtr );
    if (sizeof(AliHLTTrackMCData)+dataPtr->fCount*sizeof(AliHLTTrackMCLabel)==pBlock->fSize) {
      for( unsigned int il=0; il<dataPtr->fCount; il++ ){
	AliHLTTrackMCLabel &lab = dataPtr->fLabels[il];
	mcLabelsITS[lab.fTrackID] = lab.fMCLabel;
      }
    } else {
      HLTWarning("data mismatch in block %s (0x%08x): count %d, size %d -> ignoring track MC information", 
		 DataType2Text(pBlock->fDataType).c_str(), pBlock->fSpecification, 
		 dataPtr->fCount, pBlock->fSize);
    }
  }


  // read dEdx information (if present)

  AliHLTFloat32_t *dEdxTPC = 0; 
  Int_t ndEdxTPC = 0;
  for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(kAliHLTDataTypedEdx|kAliHLTDataOriginTPC);
       pBlock!=NULL; pBlock=GetNextInputBlock()) {
    fBenchmark.AddInput(pBlock->fSize);
    dEdxTPC = reinterpret_cast<AliHLTFloat32_t*>( pBlock->fPtr );
    ndEdxTPC = pBlock->fSize / sizeof(AliHLTFloat32_t);
    break;
  }

  // convert the TPC tracks to ESD tracks
  for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(kAliHLTDataTypeTrack|kAliHLTDataOriginTPC);
       pBlock!=NULL; pBlock=GetNextInputBlock()) {
    fBenchmark.AddInput(pBlock->fSize);
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
	if( mcLabelsTPC.find(element->TrackID())!=mcLabelsTPC.end() )
	  mcLabel = mcLabelsTPC[element->TrackID()];
	element->SetLabel( mcLabel );

	AliESDtrack iotrack;

	// for the moment, the number of clusters are not set when processing the
	// kTPCin update, only at kTPCout
	// there ar emainly three parameters updated for kTPCout
	//   number of clusters
	//   chi2
	//   pid signal
	// The first one can be updated already at that stage here, while the two others
	// eventually require to update from the ITS tracks before. The exact scheme
	// needs to be checked
	// 2010-05-12 TODO: the outer parameter is set when updating with kTPCout but
	// also when updated with kITSout. So the value from here is overwritten further
	// down. Comes along with the necessity to check the full sequence.
	iotrack.SetID( element->TrackID() );
	{
	  AliHLTGlobalBarrelTrack outPar(*element);	  
	  outPar.AliExternalTrackParam::PropagateTo( element->GetLastPointX(), fSolenoidBz );
	  outPar.SetLabel(element->GetLabel());
	  iotrack.UpdateTrackParams(&outPar,AliESDtrack::kTPCout);
	}
	iotrack.UpdateTrackParams(&(*element),AliESDtrack::kTPCin);
	iotrack.SetTPCPoints(points);
	if( element->TrackID()<ndEdxTPC ){
	  iotrack.SetTPCsignal( dEdxTPC[element->TrackID()], 0, 0 ); 
	  //AliTPCseed s;
	  //s.Set( element->GetX(), element->GetAlpha(),
	  //element->GetParameter(), element->GetCovariance() );
	  //s.SetdEdx( dEdxTPC[element->TrackID()] );
	  //s.CookPID();
	  //iotrack.SetTPCpid(s.TPCrPIDs() );
	} else {
	  if( dEdxTPC ) HLTWarning("Wrong number of dEdx TPC labels");
	}
	iotrack.SetLabel(mcLabel);
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


  // Get ITS SPD vertex
  for( const AliHLTComponentBlockData *i= GetFirstInputBlock(kAliHLTDataTypeESDVertex|kAliHLTDataOriginITS); i!=NULL; i=GetNextInputBlock() ){
    fBenchmark.AddInput(i->fSize);
  }

  for ( const TObject *iter = GetFirstInputObject(kAliHLTDataTypeESDVertex|kAliHLTDataOriginITS); iter != NULL; iter = GetNextInputObject() ) {
    AliESDVertex *vtx = dynamic_cast<AliESDVertex*>(const_cast<TObject*>( iter ) );
    pESD->SetPrimaryVertexSPD( vtx );
  }

  // now update ESD tracks with the ITSOut info
  for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(kAliHLTDataTypeTrack|kAliHLTDataOriginITSOut);
       pBlock!=NULL; pBlock=GetNextInputBlock()) {
    fBenchmark.AddInput(pBlock->fSize);
    vector<AliHLTGlobalBarrelTrack> tracks;
    if ((iResult=AliHLTGlobalBarrelTrack::ConvertTrackDataArray(reinterpret_cast<const AliHLTTracksData*>(pBlock->fPtr), pBlock->fSize, tracks))>0) {
      for (vector<AliHLTGlobalBarrelTrack>::iterator element=tracks.begin();
	   element!=tracks.end(); element++) {
	int tpcID=element->TrackID();
	// the ITS tracker assigns the TPC track used as seed for a certain track to
	// the trackID
	if( tpcID<0 || tpcID>=pESD->GetNumberOfTracks()) continue;
	AliESDtrack *tESD = pESD->GetTrack( tpcID );
	element->SetLabel(tESD->GetLabel());
	if( tESD ) tESD->UpdateTrackParams( &(*element), AliESDtrack::kITSout );
      }
    }
  }

  // now update ESD tracks with the ITS info
  for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(kAliHLTDataTypeTrack|kAliHLTDataOriginITS);
       pBlock!=NULL; pBlock=GetNextInputBlock()) {
    fBenchmark.AddInput(pBlock->fSize);
    vector<AliHLTGlobalBarrelTrack> tracks;
    if ((iResult=AliHLTGlobalBarrelTrack::ConvertTrackDataArray(reinterpret_cast<const AliHLTTracksData*>(pBlock->fPtr), pBlock->fSize, tracks))>0) {
      for (vector<AliHLTGlobalBarrelTrack>::iterator element=tracks.begin();
	   element!=tracks.end(); element++) {
	int tpcID=element->TrackID();
	// the ITS tracker assigns the TPC track used as seed for a certain track to
	// the trackID
	if( tpcID<0 || tpcID>=pESD->GetNumberOfTracks()) continue;
	Int_t mcLabel = -1;
	if( mcLabelsITS.find(element->TrackID())!=mcLabelsITS.end() )
	  mcLabel = mcLabelsITS[element->TrackID()];
	AliESDtrack *tESD = pESD->GetTrack( tpcID );
	if (!tESD) continue;
	// the labels for the TPC and ITS tracking params can be different, e.g.
	// there can be a decay. The ITS label should then be the better one, the
	// TPC label is saved in a member of AliESDtrack
	if (mcLabel>=0) {
	  // upadte only if the ITS label is available, otherwise keep TPC label
	  element->SetLabel( mcLabel );
	}
	tESD->UpdateTrackParams( &(*element), AliESDtrack::kITSin );

	// TODO: add a proper refit
	//tESD->UpdateTrackParams( &(*element), AliESDtrack::kTPCrefit );
      }
    }
  }

  // update with  vertices and vertex-fitted tracks

  for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(kAliHLTDataTypeGlobalVertexer);
       pBlock!=NULL; pBlock=GetNextInputBlock()) {
    fBenchmark.AddInput(pBlock->fSize);   
    AliHLTGlobalVertexerComponent::FillESD( pESD, reinterpret_cast<AliHLTGlobalVertexerComponent::AliHLTGlobalVertexerData* >(pBlock->fPtr) );
  }

  // loop over all tracks and set the TPC refit flag by updating with the
  // original TPC inner parameter if not yet set
  // TODO: replace this by a proper refit
  // code is comented for the moment as it does not fully solve the problems with
  // the display
  /*
  for (int i=0; i<pESD->GetNumberOfTracks(); i++) {
    if (!pESD->GetTrack(i) || 
	!pESD->GetTrack(i)->GetTPCInnerParam() ||
	pESD->GetTrack(i)->IsOn(AliESDtrack::kTPCrefit)) continue;
    AliESDtrack* tESD=pESD->GetTrack(i);
    AliHLTGlobalBarrelTrack inner(*tESD->GetTPCInnerParam());
    inner.SetLabel(tESD->GetLabel());
    tESD->UpdateTrackParams(&inner, AliESDtrack::kTPCrefit);
  }
  */

  // convert the HLT TRD tracks to ESD tracks                        
  for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(kAliHLTDataTypeTrack | kAliHLTDataOriginTRD);
       pBlock!=NULL; pBlock=GetNextInputBlock()) {
    fBenchmark.AddInput(pBlock->fSize);
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
	iotrack.UpdateTrackParams(&(*element),AliESDtrack::kTRDout);
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
  for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(kAliHLTDataTypeCaloCluster | kAliHLTDataOriginAny); pBlock!=NULL; pBlock=GetNextInputBlock()) 
    {
      fBenchmark.AddInput(pBlock->fSize);
      AliHLTCaloClusterHeaderStruct *caloClusterHeaderPtr = reinterpret_cast<AliHLTCaloClusterHeaderStruct*>(pBlock->fPtr);

      HLTDebug("%d HLT clusters from spec: 0x%X", caloClusterHeaderPtr->fNClusters, pBlock->fSpecification);

      //AliHLTCaloClusterReader reader;
      //reader.SetMemory(caloClusterHeaderPtr);

      AliHLTESDCaloClusterMaker clusterMaker;

      int nClusters = clusterMaker.FillESD(pESD, caloClusterHeaderPtr);
     
      HLTInfo("converted %d cluster(s) to AliESDCaloCluster and added to ESD", nClusters);
      iAddedDataBlocks++;
    }
  
  // Add tracks from MUON.
  for( const AliHLTComponentBlockData *i= GetFirstInputBlock(kAliHLTAnyDataType | kAliHLTDataOriginMUON); i!=NULL; i=GetNextInputBlock() ){
    fBenchmark.AddInput(i->fSize);
  }

  for (const TObject* obj = GetFirstInputObject(kAliHLTAnyDataType | kAliHLTDataOriginMUON);
       obj != NULL;
       obj = GetNextInputObject()
      )
  {
    const TClonesArray* tracklist = NULL;
    if (obj->IsA() == AliESDEvent::Class())
    {
      const AliESDEvent* event = static_cast<const AliESDEvent*>(obj);
      HLTDebug("Received a MUON ESD with specification: 0x%X", GetSpecification(obj));
      if (event->GetList() == NULL) continue;
      tracklist = dynamic_cast<const TClonesArray*>(event->GetList()->FindObject("MuonTracks"));
      if (tracklist == NULL) continue;
    }
    else if (obj->IsA() == TClonesArray::Class())
    {
      tracklist = static_cast<const TClonesArray*>(obj);
      HLTDebug("Received a MUON TClonesArray of tracks with specification: 0x%X", GetSpecification(obj));
    }
    else
    {
      // Cannot handle this object type.
      continue;
    }
    HLTDebug("Received %d MUON tracks.", tracklist->GetEntriesFast());
    if (tracklist->GetEntriesFast() > 0)
    {
      const AliESDMuonTrack* track = dynamic_cast<const AliESDMuonTrack*>(tracklist->UncheckedAt(0));
      if (track == NULL)
      {
        HLTError(Form("%s from MUON does not contain AliESDMuonTrack objects.", obj->ClassName()));
        continue;
      }
    }
    for (Int_t i = 0; i < tracklist->GetEntriesFast(); ++i)
    {
      const AliESDMuonTrack* track = static_cast<const AliESDMuonTrack*>(tracklist->UncheckedAt(i));
      pESD->AddMuonTrack(track);
    }
  }
  
  if (iAddedDataBlocks>0 && pTree) {
    pTree->Fill();
  }
  
  if (iResult>=0) iResult=iAddedDataBlocks;
  return iResult;
}
