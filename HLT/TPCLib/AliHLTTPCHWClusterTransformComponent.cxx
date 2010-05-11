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

#include "AliTPCcalibDB.h"
#include "AliTPCTransform.h"
#include "AliTPCCalPad.h"

#include "AliCDBManager.h"
#include "AliCDBEntry.h"

#include "TMath.h"
#include "TObjString.h" 
#include <cstdlib>
#include <cerrno>
#include <sys/time.h>

ClassImp(AliHLTTPCHWClusterTransformComponent) //ROOT macro for the implementation of ROOT specific class methods

const char* AliHLTTPCHWClusterTransformComponent::fgkOCDBEntryHWTransform="HLT/ConfigTPC/TPCHWClusterTransform";

AliHLTTPCHWClusterTransformComponent::AliHLTTPCHWClusterTransformComponent()
:
fOfflineTransform(NULL),
fDataId(kFALSE),
fChargeThreshold(10),
fOfflineTPCRecoParam()
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
  inputMultiplier = 3.0;
}

AliHLTComponent* AliHLTTPCHWClusterTransformComponent::Spawn() { 
  // see header file for class documentation

  return new AliHLTTPCHWClusterTransformComponent();
}
	
int AliHLTTPCHWClusterTransformComponent::DoInit( int argc, const char** argv ) { 
// see header file for class documentation

  AliTPCcalibDB* pCalib=AliTPCcalibDB::Instance();
  if(!pCalib ||
     !(fOfflineTransform = AliTPCcalibDB::Instance()->GetTransform())){
    HLTError("Cannot retrieve offline transform from AliTPCcalibDB (%p)", pCalib);
    return -ENOENT;
  }
  // set the flags in the reco param
  fOfflineTPCRecoParam.SetUseExBCorrection(1);
  fOfflineTPCRecoParam.SetUseTOFCorrection(1);
  fOfflineTransform->SetCurrentRecoParam(&fOfflineTPCRecoParam);
  
  pCalib->SetExBField(GetBz());

  int iResult=0;
  iResult = ConfigureFromCDBTObjString(fgkOCDBEntryHWTransform);

  if (iResult>=0 && argc>0)
    iResult=ConfigureFromArgumentString(argc, argv);

  return iResult;
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
     
     HLTDebug("minSlice: %d, minPartition: %d", minSlice, minPartition);
    
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
  	   cluster.fCharge  = ((UInt_t)rowCharge&0xFFFFFF)>>6; //24-bit mask to get out the charge and division with 64(>>6) for the gain correction

	   Float_t tmpPad   = *((Float_t*)&buffer[nWords+1]);
	   Float_t tmpTime  = *((Float_t*)&buffer[nWords+2]);
	   cluster.fSigmaY2 = *((Float_t*)&buffer[nWords+3]);
    	   cluster.fSigmaZ2 = *((Float_t*)&buffer[nWords+4]);
	   
	   
	   if(cluster.fCharge<fChargeThreshold) continue;
	   
	   // correct expressions for the error calculation
	   // Kenneth: 12.11.2009 I'm not sure if this is a correct calculation. Leave it out for now since it is anyway not used later since it caused segfaults.
	   // cluster.fSigmaY2 = TMath::Sqrt( *((Float_t*)&buffer[nWords+3]) - *((Float_t*)&buffer[nWords+1])* (*((Float_t*)&buffer[nWords+1])) );
	   // cluster.fSigmaZ2 = TMath::Sqrt( *((Float_t*)&buffer[nWords+3]) - *((Float_t*)&buffer[nWords+1])* (*((Float_t*)&buffer[nWords+1])) );

    	   Float_t xyz[3]; xyz[0] = xyz[1] = xyz[2] = -99.;
    	  	   
	   HLTDebug("padrow: %d, charge: %d, pad: %f, time: %f, errY: %f, errZ: %f \n", cluster.fPadRow, (UInt_t)cluster.fCharge, tmpPad, tmpTime, cluster.fSigmaY2, cluster.fSigmaZ2);        	   
	   
	   //fOfflineTransform=NULL;
	   
	   if(fOfflineTransform == NULL){	   	   
	      cluster.fPadRow += AliHLTTPCTransform::GetFirstRow(minPartition);             	   
	      AliHLTTPCTransform::Slice2Sector(minSlice, cluster.fPadRow, sector, thisrow);	      
    	      AliHLTTPCTransform::Raw2Local(xyz, sector, thisrow, tmpPad, tmpTime); 
	      if(minSlice>17) xyz[1]=(-1)*xyz[1];	   
	      cluster.fX = xyz[0];
    	      cluster.fY = xyz[1];
    	      cluster.fZ = xyz[2]; 
     	   } else {	
	   
	         
	    
	     cluster.fPadRow += AliHLTTPCTransform::GetFirstRow(minPartition);
	     
	     AliHLTTPCTransform::Slice2Sector(minSlice, (UInt_t)cluster.fPadRow, sector, thisrow);	     
	     
	      
	     Double_t x[3] = {thisrow,tmpPad+.5,tmpTime}; 
	     Int_t iSector[1]= {sector};
	     fOfflineTransform->Transform(x,iSector,0,1);
	     cluster.fX = x[0];
	     cluster.fY = x[1];
	     cluster.fZ = x[2];		     
 	   }	   

	   HLTDebug("cluster X: %f, Y: %f, Z: %f, charge: %d \n", cluster.fX, cluster.fY, cluster.fZ, (UInt_t)cluster.fCharge);
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
     
     if(fDataId==kFALSE) bd.fDataType = AliHLTTPCDefinitions::fgkClustersDataType;
     else                bd.fDataType = AliHLTTPCDefinitions::fgkAlterClustersDataType;
     
     //HLTDebug("datatype: %s", DataType2Text(bd.fDataType).c_str());
     
     outputBlocks.push_back( bd );
     
     tSize   += mysize;
     outBPtr += mysize;
     outPtr   = (AliHLTTPCClusterData*)outBPtr;
  
   } // end of loop over data blocks
   size = tSize;
   return 0;
} // end DoEvent()

int AliHLTTPCHWClusterTransformComponent::ScanConfigurationArgument(int argc, const char** argv){

  // see header file for class documentation

  if (argc<=0) return 0;
  int i=0;
  TString argument=argv[i];

  if (argument.CompareTo("-solenoidBz")==0){
    if (++i>=argc) return -EPROTO;
    argument=argv[i];
    AliTPCcalibDB*  calib=AliTPCcalibDB::Instance();
    if(!calib){
      HLTError("CalibDB not available");
    }
    Float_t magneticField = argument.Atof();
    calib->SetExBField(magneticField);
    HLTInfo("SolenoidBz is set to %f in the calibDB",magneticField);
    return 2;
  }

  if (argument.CompareTo("-change-dataId")==0){
    HLTDebug("Change data ID received.");
    fDataId = kTRUE;
    return 1;
  }
  
  if (argument.CompareTo("-charge-threshold")==0) {
    if (++i>=argc) return -EPROTO;
    argument=argv[i];
    fChargeThreshold=argument.Atof();
    HLTInfo("The charge threshold has been set to %d.", fChargeThreshold);
    return 2;
  }    

  // unknown argument
  return -EINVAL;
}

int AliHLTTPCHWClusterTransformComponent::Reconfigure(const char* /*cdbEntry*/, const char* /*chainId*/) { 
  // see header file for class documentation
  return ConfigureFromCDBTObjString(fgkOCDBEntryHWTransform);
}

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

void AliHLTTPCHWClusterTransformComponent::GetOCDBObjectDescription( TMap* const targetMap){
// Get a list of OCDB object description needed for the particular component
  if (!targetMap) return;
  targetMap->Add(new TObjString("HLT/ConfigTPC/TPCHWClusterTransform"), new TObjString("component argument for the charge threshold"));
}
