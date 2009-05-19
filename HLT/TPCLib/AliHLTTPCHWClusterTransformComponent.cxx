
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

/** @file   AliHLTTPCHWClusterTransformComponent.cxx
    @author Kalliopi Kanaki
    @date   
    @brief 
*/

#if __GNUC__>= 3
using namespace std;
#endif
#include "AliHLTTPCHWClusterTransformComponent.h"
#include "AliHLTTPCDefinitions.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCSpacePointData.h"
#include "AliHLTTPCClusterDataFormat.h"

#include <cstdlib>
#include <cerrno>
#include <sys/time.h>

ClassImp(AliHLTTPCHWClusterTransformComponent) //ROOT macro for the implementation of ROOT specific class methods

AliHLTTPCHWClusterTransformComponent::AliHLTTPCHWClusterTransformComponent()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTPCHWClusterTransformComponent::~AliHLTTPCHWClusterTransformComponent() { 
// see header file for class documentation
}

const char* AliHLTTPCHWClusterTransformComponent::GetComponentID() { 
// see header file for class documentation

  return "TPCHWClusterTransform";
}

void AliHLTTPCHWClusterTransformComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list) { 
  // see header file for class documentation

  list.clear(); 
  list.push_back( AliHLTTPCDefinitions::fgkHWClustersDataType );
}

AliHLTComponentDataType AliHLTTPCHWClusterTransformComponent::GetOutputDataType() { 
  // see header file for class documentation

  return AliHLTTPCDefinitions::fgkClustersDataType;
}

int AliHLTTPCHWClusterTransformComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList) { 
  // see header file for class documentation

  tgtList.clear();
  tgtList.push_back(AliHLTTPCDefinitions::fgkClustersDataType);
  return tgtList.size();
}

void AliHLTTPCHWClusterTransformComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ) { 
  // see header file for class documentation

  constBase = 0;
  inputMultiplier = 2.0;
}

AliHLTComponent* AliHLTTPCHWClusterTransformComponent::Spawn() { 
  // see header file for class documentation

  return new AliHLTTPCHWClusterTransformComponent();
}
	
int AliHLTTPCHWClusterTransformComponent::DoInit( int /*argc*/, const char** /*argv*/ ) { 
// see header file for class documentation

//   Int_t i = 0;
//   Char_t* cpErr;
//   
//   int iResult=0;
//   
//   TString configuration="";
//   TString argument="";
//   for (int j=0; j<argc && iResult>=0; j++) {
//     
//     argument=argv[j];
//     if (!configuration.IsNull()) configuration+=" ";
//     configuration+=argument;    
//   }
//    
//   if (!configuration.IsNull()) {
//     iResult=Configure(configuration.Data());
//   } else {
//     iResult=Reconfigure(NULL, NULL);
//   }

  return 0;
} // end DoInit()

int AliHLTTPCHWClusterTransformComponent::DoDeinit() { 
  // see header file for class documentation 
   
  return 0;
}

int AliHLTTPCHWClusterTransformComponent::DoEvent(const AliHLTComponentEventData& evtData, 
					      const AliHLTComponentBlockData* blocks, 
					      AliHLTComponentTriggerData& /*trigData*/, AliHLTUInt8_t* outputPtr, 
					      AliHLTUInt32_t& size, 
					      vector<AliHLTComponentBlockData>& outputBlocks ){
  // see header file for class documentation
    
  if(GetFirstInputBlock( kAliHLTDataTypeSOR ) || GetFirstInputBlock( kAliHLTDataTypeEOR )){
     size = 0;
     return 0;
  }

  const AliHLTComponentBlockData *iter = NULL;    
  unsigned long ndx;

  AliHLTTPCClusterData* outPtr;

  AliHLTUInt8_t* outBPtr;
  UInt_t offset, mysize, nSize, tSize = 0;

  outBPtr = outputPtr;
  outPtr  = (AliHLTTPCClusterData*)outBPtr;
  
  AliHLTTPCSpacePointData *spacePoints = outPtr->fSpacePoints;

  unsigned long maxPoints = 0;

  for(ndx=0; ndx<evtData.fBlockCnt; ndx++){
     
     iter   = blocks+ndx;
     mysize = 0;
     offset = tSize;
 
     HLTDebug("Event 0x%08LX (%Lu) received datatype: %s - required datatype: %s",
     	      evtData.fEventID, evtData.fEventID, 
     	      DataType2Text( iter->fDataType).c_str(), 
	      DataType2Text(AliHLTTPCDefinitions::fgkHWClustersDataType).c_str());                       
 
     if(iter->fDataType != AliHLTTPCDefinitions::fgkHWClustersDataType) continue;                        
  
     UInt_t minSlice     = AliHLTTPCDefinitions::GetMinSliceNr(*iter); 
     UInt_t minPartition = AliHLTTPCDefinitions::GetMinPatchNr(*iter);
     //UInt_t maxSlice     = AliHLTTPCDefinitions::GetMaxSliceNr(*iter); 
     //UInt_t maxPartition = AliHLTTPCDefinitions::GetMaxPatchNr(*iter);
    
     outPtr = (AliHLTTPCClusterData*)outBPtr;
     maxPoints = (size-tSize-sizeof(AliHLTTPCClusterData))/sizeof(AliHLTTPCSpacePointData);
     
     AliHLTUInt32_t *buffer;     
     buffer = (AliHLTUInt32_t*)iter->fPtr;  
     
     /*  
     //cluster fabrication
     buffer = new AliHLTUInt32_t[14];
     //header
     buffer[0]=0xffffffff;
     buffer[1]=0xffffffff;
     buffer[2]=0xffffffff;
     buffer[3]=0xffffffff;
     buffer[4]=0xffffffff;
     buffer[5]=0xffffffff;
     buffer[6]=0xffffffff;
     buffer[7]=0xffffffff;
     //cluster 1
     buffer[8]=0xC60002EF;
     buffer[9]=0x0;
     buffer[10]=0x0;
     buffer[11]=0x0;
     buffer[12]=0x0;
     //RCU trailer
     buffer[13]=0x80000000;
     */

     Int_t sector=-99, thisrow=-99;

     // PrintDebug(buffer, 14);
     
     // skip the first 8 32-bit CDH words
     buffer += 8;

     //PrintDebug(buffer, (Int_t)iter->fSize/sizeof(AliHLTUInt32_t));

     unsigned long nAddedClusters = 0;
     
     for(UInt_t nWords=0; nWords<(iter->fSize/sizeof(AliHLTUInt32_t)); nWords+=5){
     //     for(UInt_t nWords=0; nWords<5; nWords+=5){
         
    	// check if bit 31 and 30 of the 32-bit word is 11 -> cluster (10 is RCU trailer)
	AliHLTUInt32_t bit3130 = (buffer[nWords]>>30); // shift 30 to the right
       
	if(bit3130 == 0x3){ //beginning of a cluster

	   //PrintDebug(&buffer[nWords], 5);
	  
	   if(nAddedClusters>=maxPoints){
	      HLTWarning("No more space to add clusters, exiting!");
              break;
           }
	   
	   AliHLTTPCSpacePointData cluster = { 0.,0.,0.,0,0,0.,0.,0,0,kFALSE,0 };
              
	   //get the first word
    	   AliHLTUInt32_t  rowCharge = buffer[nWords];
    	   AliHLTUInt8_t  *rowPtr    = reinterpret_cast<AliHLTUInt8_t*>(&rowCharge);	  
	   rowPtr+=3; // this is to run for little endian architecture, the word is read from right to left
	
    	   cluster.fPadRow  = (UChar_t)((*rowPtr)&0x3f);
  	   cluster.fCharge  = (UInt_t)rowCharge&0xFFFFFF; //24-bit mask to get out the charge
    	   cluster.fSigmaY2 = (Float_t)buffer[nWords+3];
    	   cluster.fSigmaZ2 = (Float_t)buffer[nWords+4];
	       	   
    	   Float_t xyz[3]; xyz[0] = xyz[1] = xyz[2] = -99.;
	   	   
	   cluster.fPadRow += AliHLTTPCTransform::GetFirstRow(minPartition);             	   
	   AliHLTTPCTransform::Slice2Sector(minSlice, cluster.fPadRow, sector, thisrow);	      
    	   AliHLTTPCTransform::Raw2Local(xyz, sector, thisrow, (Float_t)buffer[nWords+1], (Float_t)buffer[nWords+2]);    	   
    	   
	   cluster.fX = xyz[0];
    	   cluster.fY = xyz[1];
    	   cluster.fZ = xyz[2]; 
	   
	   HLTDebug("X: %f, Y: %f, Z: %f, Charge: %d", cluster.fX,cluster.fY,cluster.fZ, cluster.fCharge);           
	   
	   spacePoints[nAddedClusters] = cluster;
	   	   
           nAddedClusters++; 
	} // end of clusters starting with 11=0x3
 	else if(bit3130 == 0x2){ // we have reached the beginning of the RCU trailer - 10=0x2
	  break;
  	}
     } // end of loop over clusters
     
     HLTDebug("Number of found clusters: %d", nAddedClusters);
     
     outPtr->fSpacePointCnt = nAddedClusters;
     nSize = sizeof(AliHLTTPCSpacePointData)*outPtr->fSpacePointCnt;
     mysize += nSize+sizeof(AliHLTTPCClusterData);
 
     AliHLTComponentBlockData bd;
     FillBlockData( bd );
     bd.fOffset = offset;
     bd.fSize = mysize;
     bd.fSpecification = iter->fSpecification;
     bd.fDataType = AliHLTTPCDefinitions::fgkClustersDataType;
     outputBlocks.push_back( bd );
     
     tSize   += mysize;
     outBPtr += mysize;
     outPtr   = (AliHLTTPCClusterData*)outBPtr;
  
   } // end of loop over data blocks
   size = tSize;
   return 0;
} // end DoEvent()
  
// int AliHLTTPCHWClusterTransformComponent::Configure(const char* arguments) { 
//   // see header file for class documentation
//   
//   int iResult=0;
//   if (!arguments) return iResult;
//   HLTInfo("parsing configuration string \'%s\'", arguments);
// 
//   TString allArgs=arguments;
//   TString argument;
//   int bMissingParam=0;
// 
//   TObjArray* pTokens=allArgs.Tokenize(" ");
//   if (pTokens) {
//     for (int i=0; i<pTokens->GetEntries() && iResult>=0; i++) {
//       argument=((TObjString*)pTokens->At(i))->GetString();
//       if (argument.IsNull()) continue;
//      
//       if (argument.CompareTo("-sum-noise-histograms")==0) {
// 	fNoiseHistograms = kTRUE;
// 	HLTInfo("got \'-sum-noise-histograms\': %s", ((TObjString*)pTokens->At(i))->GetString().Data());
// 	
//       } else if (argument.CompareTo("-sum-krypton-histograms")==0) {
// 	fKryptonHistograms = kTRUE;
// 	HLTInfo("got \'-sum-krypton-histograms\': %s", ((TObjString*)pTokens->At(i))->GetString().Data());
// 	
//       } else if (argument.CompareTo("-use-general")==0) {
// 	fUseGeneral = kTRUE;
// 	HLTInfo("got \'-use-general\': %s", ((TObjString*)pTokens->At(i))->GetString().Data());
// 	
//       } else if (argument.CompareTo("-ignore-specification")==0) {
// 	fIgnoreSpecification = kTRUE;
// 	HLTInfo("got \'-ignore-specification\': %s", ((TObjString*)pTokens->At(i))->GetString().Data());
//       }
//       else {
// 	HLTError("unknown argument %s", argument.Data());
// 	iResult=-EINVAL;
// 	break;
//       }
//     } // end for
//     
//   
//     delete pTokens;
//   
//   } // end if pTokens
//   
//   if (bMissingParam) {
//     HLTError("missing parameter for argument %s", argument.Data());
//     iResult=-EINVAL;
//   }
//   return iResult;
// }

// int AliHLTTPCHWClusterTransformComponent::Reconfigure(const char* /*cdbEntry*/, const char* /*chainId*/) { 
//   // see header file for class documentation
// 
//   int iResult=0;
//   const char* path="HLT/ConfigTPC/TPCHistogramHandlerComponent";
//   const char* defaultNotify="";
//   if (cdbEntry) {
//     path=cdbEntry;
//     defaultNotify=" (default)";
//   }
//   
//   if (path) {
//     HLTInfo("reconfigure from entry %s%s, chain id %s", path, defaultNotify,(chainId!=NULL && chainId[0]!=0)?chainId:"<none>");
//     AliCDBEntry *pEntry = AliCDBManager::Instance()->Get(path/*,GetRunNo()*/);
//     if (pEntry) {
//       TObjString* pString=dynamic_cast<TObjString*>(pEntry->GetObject());
//       if (pString) {
// 	HLTInfo("received configuration object string: \'%s\'", pString->GetString().Data());
// 	iResult=Configure(pString->GetString().Data());
//       } else {
// 	HLTError("configuration object \"%s\" has wrong type, required TObjString", path);
//       }
//     } else {
//       HLTError("cannot fetch object \"%s\" from CDB", path);
//     }
//   }
//   return iResult;
// }


void AliHLTTPCHWClusterTransformComponent::PrintDebug(AliHLTUInt32_t *buffer, Int_t size){
// see header file for class documentation 

  HLTInfo("The size is: %d", size);
  for(Int_t n32bit=0; n32bit<size; n32bit++){
    
    AliHLTUInt8_t *wordPtr = reinterpret_cast<AliHLTUInt8_t*>(&buffer[n32bit]);
    //    cout << "word ptr initialized"<<endl;
    for(Int_t w=3;w>=0;w--){
      //     cout <<"accessing word"<<endl;
      AliHLTUInt8_t word = wordPtr[w];
      //     cout<< "word was accessed"<<endl; 
      for(int n=7; n>=0; n--){
	//print the byte values
	if((((word>>n)<<7)&0x80) != 0){
	  printf("1");
	}
	else{
	  printf("0");
	}
      }
      printf("  ");
    }
    printf("\n");
  }
} // end of PrintDebug
