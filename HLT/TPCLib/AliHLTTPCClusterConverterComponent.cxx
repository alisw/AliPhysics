// $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Kalliopi Kanaki <Kalliopi.Kanaki@ift.uib.no>          *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   AliHLTTPCClusterConverterComponent.cxx
    @author Kalliopi Kanaki
    @date   
    @brief  The TPC cluster format conversion component.
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#if __GNUC__>= 3
using namespace std;
#endif

#include "AliHLTTPCClusterConverterComponent.h"
#include "AliHLTTPCOfflineCluster.h"
#include "AliTPCclusterMI.h"

#include "AliHLTTPCClusterDataFormat.h"

#include "AliHLTTPCTrackSegmentData.h"
#include "AliHLTTPCTrackletDataFormat.h"

#include "AliCDBEntry.h"
#include "AliCDBManager.h"

#include "AliHLTTPCMemHandler.h"
#include "AliHLTTPCDefinitions.h"

#include "TString.h"
#include "TObjString.h"
#include "TObjArray.h"
//#include "AliHLTTPC.h"
//#include <stdlib.h>
//#include <cerrno>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCClusterConverterComponent)

AliHLTTPCClusterConverterComponent::AliHLTTPCClusterConverterComponent()
:
fClusters(0),
fTracks(0),
fOffArray(0)
{  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTPCClusterConverterComponent::~AliHLTTPCClusterConverterComponent(){
 // see header file for class documentation
}

// Public functions to implement AliHLTComponent's interface.
// These functions are required for the registration process

const char* AliHLTTPCClusterConverterComponent::GetComponentID(){
// see header file for class documentation
  return "TPCClusterConverter";
}

void AliHLTTPCClusterConverterComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list){
// see header file for class documentation

  list.clear();
  list.push_back(AliHLTTPCDefinitions::fgkTrackSegmentsDataType);
}

AliHLTComponentDataType AliHLTTPCClusterConverterComponent::GetOutputDataType(){
// see header file for class documentation
  return kAliHLTDataTypeTObject;
}

int AliHLTTPCClusterConverterComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList){
// see header file for class documentation

  tgtList.clear();
  tgtList.push_back(kAliHLTDataTypeTObject);
  return tgtList.size();
}

void AliHLTTPCClusterConverterComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ){
// see header file for class documentation

  constBase = 0;
  inputMultiplier = 1;
}

AliHLTComponent* AliHLTTPCClusterConverterComponent::Spawn(){
// see header file for class documentation
  return new AliHLTTPCClusterConverterComponent;
}


int AliHLTTPCClusterConverterComponent::DoInit( int argc, const char** argv ){
// see header file for class documentation 
    
  //Int_t i = 0;
  //Char_t* cpErr;  
  int iResult = 0;
  
  TString configuration = "";
  TString argument = "";

  for (int i=0; i<argc && iResult>=0; i++) {    
    argument=argv[i];
    if (!configuration.IsNull()) configuration+=" ";
    configuration+=argument;    
  }
   
  if (!configuration.IsNull()) {
    iResult=Configure(configuration.Data());
  } else {
    iResult=Reconfigure(NULL, NULL);
  }

  
//   while ( i < argc ) {      
//     if (!strcmp( argv[i], "-apply-noisemap")) {
//         fApplyNoiseMap = strtoul( argv[i+1], &cpErr ,0);
//             
//     if ( *cpErr ) {
//         HLTError("Cannot convert apply-noisemap specifier '%s'.", argv[i+1]);
//         return EINVAL;
//     }
//       i+=2;
//       continue;
//     }
//     
//     if (!strcmp( argv[i], "-plot-side-a")) {
//         fPlotSideA = strtoul( argv[i+1], &cpErr ,0);
//             
//     if ( *cpErr ) {
//         HLTError("Cannot convert plot-side-a specifier '%s'.", argv[i+1]);
//         return EINVAL;
//     }
//       i+=2;
//       continue;
//     }
//     
//     if (!strcmp( argv[i], "-plot-side-c")) {
//         fPlotSideC = strtoul( argv[i+1], &cpErr ,0);
//     
//     if ( *cpErr ) {
//         HLTError("Cannot convert plot-side-c specifier '%s'.", argv[i+1]);
//         return EINVAL;
//     }
//       i+=2;
//       continue;
//     }
// 
//     if (!strcmp( argv[i], "-reset-histograms")) {
//         fResetHistograms = strtoul( argv[i+1], &cpErr ,0);
//     
//     if ( *cpErr ) {
//         HLTError("Cannot convert reset-histograms specifier '%s'.", argv[i+1]);
//         return EINVAL;
//     }
//       i+=2;
//       continue;
//     }
//                     
//     Logging(kHLTLogError, "HLT::TPCNoiseMap::DoInit", "Unknown Option", "Unknown option '%s'", argv[i] );
//     return EINVAL;
// 
//   } // end while
  
  return 0;
} // end DoInit()

int AliHLTTPCClusterConverterComponent::DoDeinit(){
// see header file for class documentation
  return 0;
}

int AliHLTTPCClusterConverterComponent::DoEvent(const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& /*trigData*/){
// see header file for class documentation

  //Logging(kHLTLogDebug, "HLT::TPCClusterConverter::DoEvent", "DoEvent", "DoEvent()");
  HLTInfo("--- Entering DoEvent() in TPCClusterConverter ---");

  if(GetFirstInputBlock(kAliHLTDataTypeSOR) || GetFirstInputBlock(kAliHLTDataTypeEOR)) return 0;
  
  Int_t numOfTotalTracks = 0;
  Int_t numOfTotalSpacePoints = 0 ;
  const AliHLTComponentBlockData *iter = NULL;

  
  // ========== LOOP OVER CLUSTER DATA =====================//
 
  for(iter = GetFirstInputBlock(AliHLTTPCDefinitions::fgkClustersDataType); iter != NULL; iter = GetNextInputBlock()){
      //AliHLTUInt8_t slice = AliHLTTPCDefinitions::GetMinSliceNr( *iter );
      //AliHLTUInt8_t patch = AliHLTTPCDefinitions::GetMinPatchNr( *iter );
  
      const AliHLTTPCClusterData *clusterData = (const AliHLTTPCClusterData*)iter->fPtr;
      numOfTotalSpacePoints+= (Int_t)clusterData->fSpacePointCnt;

      AliHLTTPCSpacePointData *clusters = (AliHLTTPCSpacePointData*)clusterData->fSpacePoints;

      for(Int_t i=0; i<(Int_t)clusterData->fSpacePointCnt; i++){ fClusters.push_back(clusters[i]); }
  } 
  
  
  
  
//   // ========== LOOP OVER TRACKS =====================//
//   
//   for(iter = GetFirstInputBlock(AliHLTTPCDefinitions::fgkTracksDataType); iter != NULL; iter = GetNextInputBlock()){
//       AliHLTUInt8_t slice = AliHLTTPCDefinitions::GetMinSliceNr( *iter );
//       AliHLTUInt8_t patch = AliHLTTPCDefinitions::GetMinPatchNr( *iter );
// 
//       const AliHLTTPCTrackletData* trackData = (const AliHLTTPCTrackletData*)iter->fPtr;
//       numOfTotalTracks += (Int_t)trackData->fTrackletCnt;
// 
//       AliHLTTPCTrackSegmentData *tracks = (AliHLTTPCTrackSegmentData*)trackData->fTracklets;
// 
//       for(int i=0;i<(Int_t)trackData->fTrackletCnt;i++){ fTracks.push_back(tracks[i]); }
//   }
 
 
 
  // ========== LOOP OVER TRACK SEGMENTS =====================//
 
  for(iter = GetFirstInputBlock(AliHLTTPCDefinitions::fgkTrackSegmentsDataType); iter != NULL; iter = GetNextInputBlock()){

      HLTInfo("Event 0x%08LX (%Lu) received datatype: %s - required datatype: %s", 
     	       evtData.fEventID, evtData.fEventID,
     	       DataType2Text(iter->fDataType).c_str(), 
     	       DataType2Text(AliHLTTPCDefinitions::fgkTrackSegmentsDataType | kAliHLTDataOriginTPC).c_str());

      if(iter->fDataType == AliHLTTPCDefinitions::fgkTrackSegmentsDataType && GetEventCount()<2){
     	 HLTWarning("data type %s is depricated, use %s (fgkTrackSegmentsDataType)!", 
     	 DataType2Text(AliHLTTPCDefinitions::fgkTrackSegmentsDataType).c_str(),
     	 DataType2Text(AliHLTTPCDefinitions::fgkTrackSegmentsDataType | kAliHLTDataOriginTPC).c_str());
      }      
    
      if(iter->fDataType != (AliHLTTPCDefinitions::fgkTrackSegmentsDataType | kAliHLTDataOriginTPC)) continue;
      
      if(iter->fDataType!=AliHLTTPCDefinitions::fgkTrackSegmentsDataType){
         HLTDebug("Data block type is not of type AliHLTTPCDefinitions::fgkTrackSegmentsDataType"); continue;
      } // endif
     
      const AliHLTTPCTrackletData *trackData = (const AliHLTTPCTrackletData*)iter->fPtr;
      numOfTotalTracks += (Int_t)trackData->fTrackletCnt;

      AliHLTTPCTrackSegmentData *tracks = (AliHLTTPCTrackSegmentData*)trackData->fTracklets;
      for(Int_t i=0; i<(Int_t)trackData->fTrackletCnt; i++){ fTracks.push_back(tracks[i]); } 

  } // end of loop over track segments
    

 
  // ========== TRIPLE LOOP FOR SETTING THE fUsed CLUSTERS =====================//

      Int_t nClustersUsed = 0;
      for(Int_t tr=0; tr<numOfTotalTracks; tr++){
	 Int_t nHits = fTracks[tr].fNPoints;
         UInt_t *hitnum = fTracks[tr].fPointIDs;
         //HLTInfo("Hits %d ", nHits);
       
           for(Int_t h=0; h<nHits; h++){
               UInt_t  idTrack        = hitnum[h];
               Int_t   sliceTrack     = AliHLTTPCSpacePointData::GetSlice(idTrack);
               Int_t   partitionTrack = AliHLTTPCSpacePointData::GetPatch(idTrack);
               UInt_t  posTrack       = AliHLTTPCSpacePointData::GetNumber(idTrack);
                 
		 fOffArray->Clear();
                 for(Int_t cl=0; cl<numOfTotalSpacePoints; cl++){       
                     UInt_t  idCluster        = fClusters[cl].fID;
                     Int_t   sliceCluster     = AliHLTTPCSpacePointData::GetSlice(idCluster);
                     Int_t   partitionCluster = AliHLTTPCSpacePointData::GetPatch(idCluster);
                     UInt_t  posCluster       = AliHLTTPCSpacePointData::GetNumber(idCluster);
		      		     
                     
		     if(sliceCluster==sliceTrack && partitionCluster==partitionTrack && posCluster==posTrack){
         	        fClusters[cl].fUsed = kTRUE;
         	        fClusters[cl].fTrackN = tr;
         	        nClustersUsed++;
		        fOffArray->Add(new AliHLTTPCOfflineCluster(fClusters[cl]));
                     } // end if
                 } // end for clusters
		 PushBack((TObject*)fOffArray, kAliHLTDataTypeTObject, 0);
           } // end for hits      
  } // end for loop over track segments
  return 0;
} // end DoEvent()

int AliHLTTPCClusterConverterComponent::Configure(const char* arguments){
// see header file for class documentation

  int iResult=0;
  if (!arguments) return iResult;
  HLTInfo("parsing configuration string \'%s\'", arguments);

  TString allArgs=arguments;
  TString argument;
  int bMissingParam=0;

  TObjArray* pTokens = allArgs.Tokenize(" ");
  if (pTokens) {
    for (int i=0; i<pTokens->GetEntries() && iResult>=0; i++) {
      argument=((TObjString*)pTokens->At(i))->GetString();
      if (argument.IsNull()) continue;
     
      if (argument.CompareTo("-apply-noisemap")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	HLTInfo("got \'-apply-noisemap\': %s", ((TObjString*)pTokens->At(i))->GetString().Data());
	
      } 
      else if (argument.CompareTo("-plot-side-c")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	HLTInfo("got \'-plot-side-c\': %s", ((TObjString*)pTokens->At(i))->GetString().Data());
	
      } 
      else if (argument.CompareTo("-plot-side-a")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	HLTInfo("got \'-plot-side-a\': %s", ((TObjString*)pTokens->At(i))->GetString().Data());
	
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

int AliHLTTPCClusterConverterComponent::Reconfigure(const char* cdbEntry, const char* chainId){
// see header file for class documentation
 
   int iResult=0;
  const char* path="HLT/ConfigTPC/ClusterConverterComponent";
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
