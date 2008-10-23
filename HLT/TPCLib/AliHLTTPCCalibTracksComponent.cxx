// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        *
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Kalliopi Kanaki <Kalliopi.Kanaki@ift.uib.no>          *
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

/** @file   AliHLTTPCCalibTracksComponent.cxx
    @author Kalliopi Kanaki
    @date   
    @brief  A calibration component for the TPC.
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt   

#if __GNUC__>= 3
using namespace std;
#endif

#include "AliHLTTPCCalibTracksComponent.h"

#include "AliHLTTPCLogging.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCClusterDataFormat.h"
#include "AliHLTTPCTrackletDataFormat.h"
#include "AliHLTTPCOfflineCluster.h"

#include "AliTPCcalibDB.h"
#include "AliTPCClusterParam.h"
#include "AliTPCcalibTracksCuts.h"
#include "AliTPCcalibTracksGain.h"
#include "AliTPCcalibTracks.h"
#include "AliTPCcalibAlign.h"

#include "AliTPCclusterMI.h"
#include "AliTPCseed.h"

#include "AliCDBEntry.h"
#include "AliCDBManager.h"

#include <stdlib.h>
#include <errno.h>
#include "TString.h"

//AliHLTTPCCalibTracksComponent gAliHLTTPCCalibTracksComponent;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCCalibTracksComponent)

AliHLTTPCCalibTracksComponent::AliHLTTPCCalibTracksComponent()
  :
  fClustParam(0),
  fTrackCuts(0),
  fCalibTracksGain(NULL),
  fCalibAlign(NULL),
  fCalibTracks(NULL),
  fMinPatch(5),
  fMaxPatch(0),
  fSpecification(0),
  fEnableAnalysis(kFALSE),
  fSeed(0),
  fClusters(),
  fTracks(),
  pConv(NULL),
  fOffArray(),
  fReadMergedTracks(0),
  fReadSliceTracks(0)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTPCCalibTracksComponent::~AliHLTTPCCalibTracksComponent(){
// see header file for class documentation
}

const char* AliHLTTPCCalibTracksComponent::GetComponentID(){
// see header file for class documentation
  return "TPCCalibTracks";
}

void AliHLTTPCCalibTracksComponent::GetInputDataTypes(vector<AliHLTComponentDataType>& list){
// see header file for class documentation

  list.clear(); 
  list.push_back(AliHLTTPCDefinitions::fgkClustersDataType);
  list.push_back(AliHLTTPCDefinitions::fgkTrackSegmentsDataType);
  list.push_back(AliHLTTPCDefinitions::fgkTracksDataType);
}

AliHLTComponentDataType AliHLTTPCCalibTracksComponent::GetOutputDataType(){
// see header file for class documentation

  return kAliHLTMultipleDataType;
}

int AliHLTTPCCalibTracksComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList){
// create output data type

  tgtList.clear();
  tgtList.push_back(AliHLTTPCDefinitions::fgkOfflineCalibAlignDataType);
  tgtList.push_back(AliHLTTPCDefinitions::fgkOfflineCalibTracksDataType);
  tgtList.push_back(AliHLTTPCDefinitions::fgkOfflineCalibTracksGainDataType);

  return tgtList.size();
}

void AliHLTTPCCalibTracksComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier){
// see header file for class documentation

  constBase = 30000000;
  inputMultiplier = 1;
}

AliHLTComponent* AliHLTTPCCalibTracksComponent::Spawn(){
// see header file for class documentation
  
  return new AliHLTTPCCalibTracksComponent();
}  


Int_t AliHLTTPCCalibTracksComponent::ScanArgument(Int_t argc, const char** argv){
// see header file for class documentation

  Int_t iResult = 0;
  TString argument = "";
  TString parameter = "";

  if(!argc) return -EINVAL;
  argument = argv[iResult];
  
  if(argument.IsNull()) return -EINVAL;

//   // -rcuformat
//   if ( argument.CompareTo("-rcuformat") == 0 ) {
// 
//     if ( ++iResult >= argc  ) {
//       iResult = -EPROTO;
//     }
//     else {
//       parameter = argv[1];
//       if ( parameter.CompareTo("old") == 0 ) {
//         HLTInfo( "RCU Format is set to old." );
//       }
//       else if ( parameter.CompareTo("new") == 0 ) {
//         HLTInfo( "RCU Format is set to new." );
//       }
//       else {
// 	HLTError( "Cannot convert rcu format specifier '%s'.", argv[1] );
// 	iResult = -EPROTO;
//       }
//     } 
//   }


  if(argument.CompareTo("-enable-analysis") == 0) {   
  
    HLTInfo( "Enable analysis before shipping data to FXS." );
    fEnableAnalysis = kTRUE;
    //cout << "fEnableAnalysis: " << fEnableAnalysis << endl;
  }
  else if(argument.CompareTo("-read-slice-tracks") == 0){
    HLTInfo( "Reading slice tracks..." );
    fReadSliceTracks = kTRUE;  
    //cout << "fReadSliceTracks: " << fReadSliceTracks << endl;
  }
  else if(argument.CompareTo("-read-merged-tracks") == 0){
    HLTInfo( "Reading merged tracks..." );
    fReadMergedTracks = kTRUE; 
    //cout << "fReadMergedTracks: " << fReadMergedTracks << endl;
  }
  else {
    iResult = -EINVAL;
  }
  
  return iResult;



//---------------------------------


 
//   Int_t i = 0;
//   Char_t* cpErr;
//   
//   int iResult = 0;
//   
//   TString configuration = "";
//   TString argument = "";
//   for(int j=0; j<argc && iResult>=0; j++){
//     
//     argument=argv[j];
//     if (!configuration.IsNull()) configuration+=" ";
//     configuration+=argument;    
//   }
//    
//   if (!configuration.IsNull()) {
//     iResult = Configure(configuration.Data());
//   } else {
//     iResult = Reconfigure(NULL, NULL);
//   }
//  
//   while ( i < argc ) {      
//     if (!strcmp( argv[i], "-enable-analysis")) {
//         fEnableAnalysis = strtoul( argv[i+1], &cpErr ,0);
//             
//     if ( *cpErr ) {
//         HLTError("Cannot convert enable-analysis specifier '%s'.", argv[i+1]);
//         return EINVAL;
//     }
//       i+=2;
//       continue;
//     }
//     
//     if (!strcmp( argv[i], "-read-merged-tracks")) {
//         fReadMergedTracks = strtoul( argv[i+1], &cpErr ,0);
//             
//     if ( *cpErr ) {
//         HLTError("Cannot convert read-merged-tracks specifier '%s'.", argv[i+1]);
//         return EINVAL;
//     }
//       i+=2;
//       continue;
//     }
//     
//     if (!strcmp( argv[i], "-read-slice-tracks")) {
//         fReadSliceTracks = strtoul( argv[i+1], &cpErr ,0);
//     
//     if ( *cpErr ) {
//         HLTError("Cannot convert read-slice-tracks specifier '%s'.", argv[i+1]);
//         return EINVAL;
//     }
//       i+=2;
//       continue;
//     }
// 
//     Logging(kHLTLogError, "HLT::TPCCalibTracks::ScanArgument", "Unknown Option", "Unknown option '%s'", argv[i] );
//     return EINVAL;
// 
//   } // end while
//   
//   return 0;
}

Int_t AliHLTTPCCalibTracksComponent::InitCalibration(){
// see header file for class documentation
  
  HLTInfo("init calibration component");
  int iResult=0;

  if(fCalibTracksGain) return EINPROGRESS; 
  fCalibTracksGain = new AliTPCcalibTracksGain();
  
  if(fCalibAlign) return EINPROGRESS;  
  fCalibAlign = new AliTPCcalibAlign();
 
  if(fCalibTracks) return EINPROGRESS;  
  fCalibTracks = new AliTPCcalibTracks();

  
//   // Init parameters and cuts
//   fClustParam = AliTPCcalibDB::Instance()->GetClusterParam();
//   fTrackCuts = new AliTPCcalibTracksCuts(20, 0.4, 0.5, 0.13, 0.018);
// 
//   // Init calibration componenets
//   fCalibAlign	   = new AliTPCcalibAlign("TPCcalibAlign","TPCcalibAlign");
//   fCalibTracksGain = new AliTPCcalibTracksGain("TPCcalibTracksGain","TPCcalibTracksGain",fTrackCuts);
//   fCalibTracks     = new AliTPCcalibTracks("TPCcalibTracks","TPCcalibTracks",fClustParam,fTrackCuts);
// 
//   if (!fTrackCuts || !fClustParam ||  !fCalibAlign || !fCalibTracksGain || !fCalibTracks) {
//     HLTError("failed creating internal objects");
//     iResult=-ENOMEM;
//   }
  return iResult;  

}

Int_t AliHLTTPCCalibTracksComponent::DeinitCalibration(){
// see header file for class documentation
   
   if(fClustParam) delete fClustParam; fClustParam = 0; 
   if(fTrackCuts)  delete fTrackCuts;  fTrackCuts  = 0; 

   if(fCalibTracksGain) delete fCalibTracksGain; fCalibTracksGain = NULL;
   if(fCalibAlign)      delete fCalibAlign;	 fCalibAlign	  = NULL;
   //if(fCalibTracks)     delete fCalibTracks;	 fCalibTracks     = NULL;
   
  return 0;
}

Int_t AliHLTTPCCalibTracksComponent::ProcessCalibration(const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/){
// see header file for class documentation
  
  
  const AliHLTComponentBlockData *iter = NULL;
  //AliHLTUInt8_t slice = 0;
  //AliHLTUInt8_t patch = 0;
      
  if(GetFirstInputBlock( kAliHLTDataTypeSOR ) || GetFirstInputBlock( kAliHLTDataTypeEOR)) return 0;

  Int_t TotalTrack = 0;

  //if(fReadMergedTracks){ //Reading Merged Tracks 
     for(iter = GetFirstInputBlock(AliHLTTPCDefinitions::fgkTracksDataType); iter != NULL; iter = GetNextInputBlock()){
         if(iter->fDataType != AliHLTTPCDefinitions::fgkTracksDataType) continue; 
         ReadTracks(iter,TotalTrack);  
     }  
  //} // end if reading merged tracks

  //else 
  //if(fReadSliceTracks){ //Reading Tracks from slice 
     for(iter = GetFirstInputBlock(AliHLTTPCDefinitions::fgkTrackSegmentsDataType); iter != NULL; iter = GetNextInputBlock()){
         if(iter->fDataType != AliHLTTPCDefinitions::fgkTrackSegmentsDataType) continue; 
         ReadTracks(iter,TotalTrack);
     }
  //} // end if reading slice tracks
	  
  //else HLTError("No input preference (merged/slice tracks) specified!");
 
  //int TotalSpacePoint = 0;
  //Int_t nClustersUsed = 0;
  //fOffArray.Clear();
  
  for(iter = GetFirstInputBlock(AliHLTTPCDefinitions::fgkClustersDataType); iter != NULL; iter = GetNextInputBlock()){
      if(iter->fDataType!=AliHLTTPCDefinitions::fgkClustersDataType) continue;

      //AliHLTUInt8_t slice = AliHLTTPCDefinitions::GetMinSliceNr(*iter);
      //AliHLTUInt8_t patch = AliHLTTPCDefinitions::GetMinPatchNr(*iter);

      const AliHLTTPCClusterData *clusterData = (const AliHLTTPCClusterData*)iter->fPtr;
      Int_t nSpacepoint = (Int_t)clusterData->fSpacePointCnt;	
      
      //cout << "nSpacepoint: " << nSpacepoint << endl;      
      //TotalSpacePoint += nSpacepoint;
      
      AliHLTTPCSpacePointData *clusters = (AliHLTTPCSpacePointData*)clusterData->fSpacePoints;
     
      Int_t thissector, thisrow = -99;
      pConv = new AliHLTTPCOfflineCluster();
      
      for(Int_t i=0; i<nSpacepoint; i++){
       
          UInt_t idCluster = clusters[i].fID;
          Int_t  sliceCl   = (idCluster>>25) & 0x7f;
          Int_t  patchCl   = (idCluster>>22) & 0x7;
          UInt_t pos       = idCluster&0x3fffff;
	  
	  //cout << idCluster <<"\t"<< patchCl << "\t"<< pos << endl;
	  //cout << fTrackClusterID[sliceCl][patchCl].size()<< endl;
              
          for(UInt_t id=0; id<fTrackClusterID[sliceCl][patchCl].size(); id++){
	    //cout << "HELLOOOOOOOOOOOOO" << endl;
   	      
	      if(fTrackClusterID[sliceCl][patchCl][id]==pos){
    	         clusters[i].fUsed = kTRUE;
	         AliHLTTPCTransform::Slice2Sector(sliceCl,(Int_t)clusters[i].fPadRow,thissector,thisrow);
	         //AliHLTTPCTransform::Slice2Sector(fCurrentslice,fCurrentRow,thissector,thisrow);
    	         //nClustersUsed++;     
		 
		 //cout <<  (Int_t)clusters[i].fPadRow << endl; 
	      
                 
		 AliTPCclusterMI *offClus = pConv->ConvertHLTToOffline(clusters[i]);
		 offClus->SetDetector(thissector);
		 //offClus->IsUsed(1);
		 fOffArray.Add(pConv->ConvertHLTToOffline(clusters[i]));
		 //delete pConv;
		 //delete offClus;

    	     } // end if setting used clusters
          } //end for loop over cluster id       
      } // end for loop over spacepoints
      
      //----- turn TObjArray of AliTPCclusterMI to AliTPCseed -----
     
      delete pConv;
      for(Int_t k=0; k<fOffArray.GetSize(); k++) delete fOffArray.At(k);
      
      fOffArray.Clear();
     
      //fCalibTracksGain->Process(fSeed);
      //fCalibAlign->Process(fSeed);
      //fCalibTracks->Process(fSeed);
  } // end for loop over cluster data type
    
  //fSpecification = AliHLTTPCDefinitions::EncodeDataSpecification( slice, slice, fMinPatch, fMaxPatch );

  //PushBack((TObject*)fCalibTracksGain, AliHLTTPCDefinitions::fgkOfflineCalibTracksGainDataType, fSpecification);
  //PushBack((TObject*)fCalibAlign,      AliHLTTPCDefinitions::fgkOfflineCalibAlignDataType,	fSpecification);
  //PushBack((TObject*)fCalibTracks,     AliHLTTPCDefinitions::fgkOfflineCalibTracksDataType,	fSpecification);
  
  //fClusters.clear();
  //fTracks.clear();

  return 0;
} // end ProcessCalibration


Int_t AliHLTTPCCalibTracksComponent::ShipDataToFXS(const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/){
// see header file for class documentation
    
  if(fEnableAnalysis){
     fCalibTracksGain->Analyze();
     fCalibAlign->Analyze();
     fCalibTracks->Analyze();
  }
  
  PushToFXS((TObject*)fCalibAlign,	"TPC", "TPCcalibAlign"); 
  PushToFXS((TObject*)fCalibTracksGain, "TPC", "TPCcalibTracksGain");
  PushToFXS((TObject*)fCalibTracks,	"TPC", "TPCcalibTracks");
  
  return 0;
} // end ShipDataToFXS


void AliHLTTPCCalibTracksComponent::ReadTracks(const AliHLTComponentBlockData *iter, Int_t &tt){

  HLTDebug("Input Data - TPC cluster - Slice/Patch: %d/%d.", 
	   AliHLTTPCDefinitions::GetMinSliceNr( *iter ), 
	   AliHLTTPCDefinitions::GetMinPatchNr( *iter ));
 

  const AliHLTTPCTrackletData *trackData = (const AliHLTTPCTrackletData*)iter->fPtr;
  AliHLTUInt32_t nTracks = trackData->fTrackletCnt;

  tt += nTracks;
 
  //cout << "=============================================================================" << endl;
  //cout << "=============================================================================" << endl;
  //cout << nTracks << endl;
 
  AliHLTTPCTrackSegmentData *tracks = (AliHLTTPCTrackSegmentData*)trackData->fTracklets;
  
  for(AliHLTUInt32_t i=0;i<nTracks;i++){
     
      fTracks.push_back(tracks[i]);
      UInt_t nHits = tracks->fNPoints;
      const UInt_t *hitnum = tracks->fPointIDs;
      
      for(UInt_t h=0; h<nHits; h++){
          
	  UInt_t idTrack = hitnum[h];
          Int_t sliceTrack = (idTrack>>25) & 0x7f;
          Int_t patchTrack = (idTrack>>22) & 0x7;
          UInt_t pos = idTrack&0x3fffff;
          fTrackClusterID[sliceTrack][patchTrack].push_back(pos);
      } // end for loop over hits
       
      UChar_t *tmpP = (UChar_t*)tracks;
      tmpP += sizeof(AliHLTTPCTrackSegmentData)+tracks->fNPoints*sizeof(UInt_t);
      tracks = (AliHLTTPCTrackSegmentData*)tmpP;
  } // end for loop over tracks
} // end ReadTracks()


Int_t AliHLTTPCCalibTracksComponent::Configure(const char* arguments){ 
// see header file for class documentation
  
  Int_t iResult=0;
  if (!arguments) return iResult;
  HLTInfo("parsing configuration string \'%s\'", arguments);

  TString allArgs = arguments;
  TString argument;
  int bMissingParam = 0;

  TObjArray* pTokens=allArgs.Tokenize(" ");
  if (pTokens) {
    for (int i=0; i<pTokens->GetEntries() && iResult>=0; i++) {
      argument=((TObjString*)pTokens->At(i))->GetString();
      if (argument.IsNull()) continue;
     
      if (argument.CompareTo("-enable-analysis")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	HLTInfo("got \'-enable-analysis\': %s", ((TObjString*)pTokens->At(i))->GetString().Data());
	
      } 
      else if (argument.CompareTo("-read-merged-tracks")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	HLTInfo("got \'-read-merged-tracks\': %s", ((TObjString*)pTokens->At(i))->GetString().Data());
	
      } 
      else if (argument.CompareTo("-read-slice-tracks")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	HLTInfo("got \'-read-slice-tracks\': %s", ((TObjString*)pTokens->At(i))->GetString().Data());
	
      } 
      else {
	HLTError("unknown argument %s", argument.Data());
	iResult=-EINVAL;
	break;
      }
    } // end for
  
    delete pTokens;
  
  } // end if pTokens
  
  if (bMissingParam) {
    HLTError("missing parameter for argument %s", argument.Data());
    iResult=-EINVAL;
  }
  return iResult;  
}

int AliHLTTPCCalibTracksComponent::Reconfigure(const char* cdbEntry, const char* chainId){ 
// see header file for class documentation
 
  int iResult=0;
  const char* path = "HLT/ConfigTPC/TPCCalibTracksComponent";
  const char* defaultNotify="";
  if (cdbEntry) {
      path=cdbEntry;
      defaultNotify = " (default)";
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
      HLTError("cannot fetch object \"%s\" from CDB", path);
    }
  }
  return iResult;
}

