// $Id$
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Jochen Thaeder <thaeder@kip.uni-heidelberg.de>        *
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

/** @file   AliHLTTPCCalibTimeComponent.cxx
    @author Kalliopi Kanaki
    @date   2009-07-08
    @brief  A calibration component for interfacing the offline calculation of TPC drift velocity correction 
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt   

#if __GNUC__>= 3
using namespace std;
#endif


#include "AliHLTTPCCalibTimeComponent.h"
#include "AliHLTTPCDefinitions.h"
#include "AliHLTMisc.h"

#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDfriend.h"

#include "AliTPCcalibTime.h"
#include "AliTPCcalibCalib.h"
#include "AliTPCseed.h"
#include "AliTPCcalibDB.h"
#include "AliTPCClusterParam.h"

#include "AliHLTTPCOfflineCluster.h"
#include "AliHLTTPCSpacePointData.h"
#include "AliHLTTPCTrackletDataFormat.h"
#include "AliHLTExternalTrackParam.h"
#include "AliHLTGlobalBarrelTrack.h"
#include "AliHLTTPCTransform.h"

#include "TObjArray.h"
#include "TString.h"
#include "TFile.h"

#include "THnSparse.h"
#include "TGraphErrors.h"

#include <cstdlib>
#include <cerrno>

#include "AliHLTReadoutList.h"

ClassImp(AliHLTTPCCalibTimeComponent) // ROOT macro for the implementation of ROOT specific class methods

AliHLTTPCCalibTimeComponent::AliHLTTPCCalibTimeComponent()
  :
   fCalibTime(NULL)
  ,fCal(NULL)
  ,fESDevent(NULL)
  ,fESDtrack(NULL)
  ,fESDfriend(NULL)
  ,fSeedArray(NULL)
  ,fOutputSize(50000)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt  
  
  for(int i=0; i<fkNPartition; i++){
      fPartitionClusters[i]  = 0;    
      fNPartitionClusters[i] = 0;    
  }
}

const char* AliHLTTPCCalibTimeComponent::fgkOCDBEntry="HLT/ConfigTPC/TPCCalibTime";

AliHLTTPCCalibTimeComponent::~AliHLTTPCCalibTimeComponent(){
// see header file for class documentation 
 
  for(int i=0; i<fkNPartition; i++){
      delete[] fPartitionClusters[i];
  }
}

const char* AliHLTTPCCalibTimeComponent::GetComponentID() {
// see header file for class documentation
  return "TPCCalibTime";
}

void AliHLTTPCCalibTimeComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list){
// see header file for class documentation

  list.clear(); 
  list.push_back( kAliHLTDataTypeTObjArray|kAliHLTDataOriginTPC ); // output of TPCCalibSeedMaker
  list.push_back( AliHLTTPCDefinitions::fgkClustersDataType );     // output of the TPC CF
  list.push_back( kAliHLTDataTypeTrack|kAliHLTDataOriginTPC );     // output of the global merger
  list.push_back( kAliHLTDataTypeESDObject|kAliHLTDataOriginOut ); // output of global esd converter
}

AliHLTComponentDataType AliHLTTPCCalibTimeComponent::GetOutputDataType() {
// see header file for class documentation

  return AliHLTTPCDefinitions::fgkCalibCEDataType|kAliHLTDataOriginOut;
}

void AliHLTTPCCalibTimeComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ){
// see header file for class documentation

  constBase = fOutputSize;
  inputMultiplier = 0; // to be estimated
}

AliHLTComponent* AliHLTTPCCalibTimeComponent::Spawn(){
// see header file for class documentation

  return new AliHLTTPCCalibTimeComponent();
}  


Int_t AliHLTTPCCalibTimeComponent::ScanConfigurationArgument( Int_t argc, const char** argv ){
// see header file for class documentation
 
  if (argc<=0) return 0;
  int i=0;
  TString argument=argv[i];

  // -output-size
  if (argument.CompareTo("-output-size")==0) {
    if (++i>=argc) return -EPROTO;
    argument=argv[i];
    fOutputSize=argument.Atoi();
    return 2;
  }
  return -EINVAL;
}

Int_t AliHLTTPCCalibTimeComponent::InitCalibration(){
// see header file for class documentation
  
  AliTPCcalibDB::Instance()->SetRun(AliHLTMisc::Instance().GetCDBRunNo());
  AliTPCcalibDB::Instance()->GetClusterParam()->SetInstance(AliTPCcalibDB::Instance()->GetClusterParam());
  

//   AliTPCcalibDB *calib = AliTPCcalibDB::Instance();
// 
//   if(!calib){
//     HLTError("AliTPCcalibDB does not exist");
//     return -ENOENT;
//   }
//   
//   AliTPCClusterParam *clusPar = calib->GetClusterParam();
//   if(!clusPar){
//     HLTError("OCDB entry TPC/Calib/ClusterParam (AliTPCcalibDB::GetClusterParam()) is not available.");
//     return -ENOENT;
//   }

  // first configure the default
  int iResult=0;
  if (iResult>=0) iResult=ConfigureFromCDBTObjString(fgkOCDBEntry);
    
  if(fCalibTime) return EINPROGRESS;
  fCal = new AliTPCcalibCalib();
  
  fSeedArray = new TObjArray();
  
  return iResult;
 
}

Int_t AliHLTTPCCalibTimeComponent::DeinitCalibration() {
// see header file for class documentation

  if(fCalibTime) delete fCalibTime; fCalibTime = NULL;
  if(fCal)       delete fCal;	    fCal       = NULL;
  if(fSeedArray) delete fSeedArray; fSeedArray = NULL;

  //if(fESDfriend) delete fESDfriend; fESDfriend = NULL;
  
  //if(arr) delete arr; arr = NULL;
  
  return 0;
}

int AliHLTTPCCalibTimeComponent::Reconfigure(const char* cdbEntry, const char* /*chainId*/){
// see header file for class documentation

  // configure from the specified antry or the default one
  const char* entry=cdbEntry;
  if(!entry || entry[0]==0){
     entry=fgkOCDBEntry;
  }

  return ConfigureFromCDBTObjString(entry);
}

Int_t AliHLTTPCCalibTimeComponent::ProcessCalibration( const AliHLTComponentEventData& /*evtData*/,  AliHLTComponentTriggerData& /*trigData*/ ){
// see header file for class documentation

  if(GetFirstInputBlock( kAliHLTDataTypeSOR ) || GetFirstInputBlock( kAliHLTDataTypeEOR )) return 0;

  const AliHLTComponentBlockData *iter = NULL;      

  //--------------- output over TObjArray of AliTPCseed objects (output of TPCSeedMaker) -------------------//
  
  // A previous component in the chain (TPCSeedMaker) has processed the TPC clusters and tracks and created a TClonesArray of AliTPCseed objects
  // In this loop the iterator accesses this array stored in memory, in order to use it in the next loop over the AliESDevent of the HLT
  
//   for(iter = (TObject*)GetFirstInputObject(kAliHLTDataTypeTObjArray|kAliHLTDataOriginTPC); iter != NULL; iter = (TObject*)GetNextInputObject()){  
//               
//       if(GetDataType(iter) != (kAliHLTDataTypeTObjArray | kAliHLTDataOriginTPC)) continue;      
//       fSeedArray = dynamic_cast<TClonesArray*>(iter);      
//    }
 
 // int nInputClusters = 0;
 // int nInputTracks = 0;

  //TObjArray *arr = new TObjArray(1000);
  //arr->SetOwner(kTRUE);
  fSeedArray->Clear();

  
  for( int i=0; i<fkNPartition; i++ ){
    delete[] fPartitionClusters[i];    
    fPartitionClusters[i] = 0;
    fNPartitionClusters[i] = 0;    
  }

  
  
  //------------------- loop over clusters -------------//
  
  for(iter = GetFirstInputBlock(AliHLTTPCDefinitions::fgkClustersDataType); iter != NULL; iter = GetNextInputBlock()){
      
      if ( iter->fDataType != AliHLTTPCDefinitions::fgkClustersDataType ) continue;
    
      Int_t slice     = AliHLTTPCDefinitions::GetMinSliceNr(iter->fSpecification);
      Int_t partition = AliHLTTPCDefinitions::GetMinPatchNr(iter->fSpecification);
       
      Int_t slicepartition = slice*6+partition;
      
      if(slicepartition > fkNPartition){
         HLTWarning("Wrong header of TPC cluster data, slice %d, partition %d", slice, partition );
         continue;
      }
      
      AliHLTTPCClusterData* inPtrSP = ( AliHLTTPCClusterData* )( iter->fPtr );
      // nInputClusters += inPtrSP->fSpacePointCnt;

      delete[] fPartitionClusters[slicepartition];
      fPartitionClusters[slicepartition]  = new AliTPCclusterMI[inPtrSP->fSpacePointCnt];
      fNPartitionClusters[slicepartition] = inPtrSP->fSpacePointCnt;
    
      // create  offline clusters out of the HLT clusters
      // todo: check which cluster information is really needed for the dEdx
      for(unsigned int i = 0; i < inPtrSP->fSpacePointCnt; i++){          
	  AliHLTTPCSpacePointData *chlt = &( inPtrSP->fSpacePoints[i] );
          AliTPCclusterMI *c = fPartitionClusters[slicepartition]+i;
          c->SetX(chlt->fX);
          c->SetY(chlt->fY);
          c->SetZ(chlt->fZ);
          c->SetSigmaY2(chlt->fSigmaY2);
          c->SetSigmaYZ( 0 );
          c->SetSigmaZ2(chlt->fSigmaZ2);
          c->SetQ( chlt->fCharge );
          c->SetMax( chlt->fQMax );
          Int_t sector, row;
          Float_t padtime[3] = {0,chlt->fY,chlt->fZ};
          AliHLTTPCTransform::Slice2Sector(slice,chlt->fPadRow, sector, row);
          AliHLTTPCTransform::Local2Raw( padtime, sector, row);
          c->SetDetector( sector );
          c->SetRow( row );
          c->SetPad( (Int_t) padtime[1] );
          c->SetTimeBin( (Int_t) padtime[2] );
      }      
  } // end of loop over blocks of clusters
  
  
  
  
  //---------- loop over merged tracks ------------------ //
  int nTracks = 0;
  for(const AliHLTComponentBlockData *pBlock = GetFirstInputBlock(kAliHLTDataTypeTrack|kAliHLTDataOriginTPC); pBlock != NULL; pBlock = GetNextInputBlock()){
 
      AliHLTTracksData *dataPtr = (AliHLTTracksData*) pBlock->fPtr;
      nTracks = dataPtr->fCount;
    
      AliHLTExternalTrackParam *currTrack = dataPtr->fTracklets;
      
      //nInputTracks += nTracks;
    
      for(int itr=0;  itr<nTracks && ( (AliHLTUInt8_t *)currTrack < ((AliHLTUInt8_t *) pBlock->fPtr)+pBlock->fSize); itr++){    
          
	  // create an offline track
          AliHLTGlobalBarrelTrack gb(*currTrack);
          AliTPCseed tTPC;
          tTPC.Set(gb.GetX(), gb.GetAlpha(), gb.GetParameter(), gb.GetCovariance());
            
          // set the cluster pointers     
          for(UInt_t ic=0; ic<currTrack->fNPoints; ic++){	
      
              tTPC.SetNumberOfClusters(currTrack->fNPoints);
          
	      UInt_t id = currTrack->fPointIDs[ic];
	      int iSlice = AliHLTTPCSpacePointData::GetSlice(id);
	      int iPartition = AliHLTTPCSpacePointData::GetPatch(id);
	      int iCluster = AliHLTTPCSpacePointData::GetNumber(id);	
	
	      if(iSlice<0 || iSlice>36 || iPartition<0 || iPartition>5){
	         HLTError("Corrupted TPC cluster Id: slice %d, partition %d, cluster %d", iSlice, iPartition, iCluster);
	         continue;
	      }
	
	      AliTPCclusterMI *partitionClusters = fPartitionClusters[iSlice*6 + iPartition];
	      
	      if(!partitionClusters){
	          HLTError("Clusters are missed for slice %d, partition %d", iSlice, iPartition );
	          continue;
	      }
	
	      if(iCluster >= fNPartitionClusters[iSlice*6 + iPartition]){
	         HLTError("TPC slice %d, partition %d: ClusterID==%d >= N Clusters==%d ", iSlice, iPartition,iCluster, fNPartitionClusters[iSlice*6 + iPartition]);
	         continue;
	      }
	
	      AliTPCclusterMI *c = &(partitionClusters[iCluster]);	  	
    	      int sec = c->GetDetector();
	      int row = c->GetRow();
	      if(sec >= 36) row = row + AliHLTTPCTransform::GetNRowLow();
	
	      tTPC.SetClusterPointer(row, c);	
	
	      AliTPCTrackerPoint &point = *( tTPC.GetTrackPoint( row ) );
	      //tTPC.Propagate( TMath::DegToRad()*(sec%18*20.+10.), c->GetX(), fSolenoidBz );
	      Double_t angle2 = tTPC.GetSnp()*tTPC.GetSnp();
	      angle2 = (angle2<1) ?TMath::Sqrt(angle2/(1-angle2)) :10.; 
	      point.SetAngleY( angle2 );
	      point.SetAngleZ( tTPC.GetTgl() );
          } // end of associated cluster loop
       
      AliTPCseed *seed = &(tTPC);
      fSeedArray->Add(seed);
     
      unsigned int step = sizeof( AliHLTExternalTrackParam ) + currTrack->fNPoints * sizeof( unsigned int );
      currTrack = ( AliHLTExternalTrackParam* )( (( Byte_t * )currTrack) + step );  

      }// end of vector track loop           
  } // end of loop over blocks of merged tracks  
  
  HLTInfo("Number of reconstructed tracks %d, number of produced seeds %d\n", nTracks, fSeedArray->GetEntries());

 
  //----------- loop over output of global esd converter ----------------//
  
  // In this loop we access the AliESDevent that was produced by the HLT and is stored in memory. There should exist 1 object 
  // of type kAliHLTDataTypeESDObject per event.
 
  TObject *iterOb = NULL; 
  for(iterOb = (TObject*)GetFirstInputObject(kAliHLTDataTypeESDObject | kAliHLTDataOriginOut); iterOb != NULL; iterOb = (TObject*)GetNextInputObject()){   
      
      if(GetDataType(iterOb) != (kAliHLTDataTypeESDObject | kAliHLTDataOriginOut)) continue;
            
      fESDevent = dynamic_cast<AliESDEvent*>(iterOb);
      fESDevent->GetStdContent();    
                 
      HLTInfo("Number of seeds: %i\n", fSeedArray->GetEntriesFast()); // access of the info from the previous loop over the AliTPCseed array        
   
      fCal->UpdateEventInfo(fESDevent);     
      for(Int_t i=0; i<fSeedArray->GetEntriesFast(); i++){  // loop over TObjArray with seeds
	 
 	  AliTPCseed *seed = (AliTPCseed*)fSeedArray->UncheckedAt(i);                    
          fESDtrack = fESDevent->GetTrack(i);
   	  if(!fESDtrack || !seed) continue; 

          //if(fESDtrack->GetID() != seed->GetLabel()) { 
          //   HLTWarning("Mismatch of track id between seed and ESD track: %i, %i\n", fESDtrack->GetID(), seed->GetLabel());
          //   continue;	    
          //}
	  
	  if(seed->GetNumberOfClusters()==0) continue;      
	  fCal->RefitTrack(fESDtrack, seed, GetBz()); // update AliESDtrack and AliTPCseed info, acccording to Marian's request
	       
	  AliTPCseed *seedCopy = new AliTPCseed(*seed, kTRUE); 
	  fESDtrack->AddCalibObject(seedCopy);  // add the AliTPCseed as a friend track to the AliESDtrack (to be accessed in TPC/AliTPCcalibTime.cxx)              
	
	  //fESDfriendTrack = const_cast<AliESDfriendTrack*>(fESDtrack->GetFriendTrack());        
      }
  } 
  
  if(!fCalibTime){ // create the calibration object that will call the offline functions
  
     Int_t startTime = fESDevent->GetTimeStamp()-60*60*1;  //Start time one hour before first event, will make precise cuts later.
     Int_t   endTime = fESDevent->GetTimeStamp()+60*60*23; //End time 23 hours after first event.
     fCalibTime = new AliTPCcalibTime("calibTime","time dependent Vdrift calibration", startTime, endTime, 20*60);
     fCalibTime->SetStreamLevel(20);
     fCalibTime->SetDebugLevel(20);
     printf("fCalibTime = %i, startTime = %i, endTime = %i \n", fCalibTime!=0, startTime, endTime);
  }
  
  fESDfriend = new AliESDfriend();
  fESDevent->GetESDfriend(fESDfriend);
  fESDevent->SetESDfriend(fESDfriend);
  fESDevent->AddObject(fESDfriend); 
  // create the AliESDfriend and add it to the event, now both the friend tracks and the friends are available for the offline functions to be called

  fCalibTime->UpdateEventInfo(fESDevent); // needed for getting the run number and time stamp information correct on the offline side
  fCalibTime->Process(fESDevent);         // first offline function called
  
  // delete fESDfriend;
  
  //PushBack( (TObject*)fCalibTime, AliHLTTPCDefinitions::fgkCalibCEDataType | kAliHLTDataOriginOut, 0x0);
  
  return 0;
}

Int_t AliHLTTPCCalibTimeComponent::ShipDataToFXS( const AliHLTComponentEventData& /*evtData*/,  AliHLTComponentTriggerData& /*trigData*/ ){
// see header file for class documentation

  HLTInfo("Shipping data to FXS...\n");
 
  fCalibTime->Analyze(); // called at the end of the run or event modulo
  
  // the rest of the histogram and graph declarations were copied by Dag as a first attempt to get the start/end time bin "automatically". Perhaps we need some more
  // thinking here to avoid copying all these lines that might chain in offline without HLT realizing.
  
  THnSparse* addHist = fCalibTime->GetHistoDrift("all");
  if(!addHist) return -1;

  //Identifying used range of histogram
 
  Int_t startTimeBin = 0;
  Int_t endTimeBin   = 0;

  TH1D *histoTime = addHist->Projection(0);
  if(histoTime){
    startTimeBin = histoTime->FindFirstBinAbove(0);
    endTimeBin   = histoTime->FindLastBinAbove(0);
    printf("startTimeBin       = %i endTimeBin       = %i\n", startTimeBin, endTimeBin);
    printf("startTimeBinCentre = %f endTimeBinCentre = %f\n", histoTime->GetBinCenter(startTimeBin), histoTime->GetBinCenter(endTimeBin));
    printf("startTimeBinWidth  = %f endTimeBinWidth  = %f\n", histoTime->GetBinWidth(startTimeBin),  histoTime->GetBinWidth(endTimeBin));
    delete histoTime; histoTime = 0;
  }

  Int_t startPtBin = 0;
  Int_t endPtBin   = 0;
  TH1D *histoPt = addHist->Projection(1);
  if(histoPt){
    startPtBin = histoPt->FindFirstBinAbove(0);
    endPtBin   = histoPt->FindLastBinAbove(0);
    printf("startPtBin       = %i endPtBin       = %i\n", startPtBin, endPtBin);
    printf("startPtBinCentre = %f endPtBinCentre = %f\n", histoPt->GetBinCenter(startPtBin), histoPt->GetBinCenter(endPtBin));
    printf("startPtinWidth   = %f endPtBinWidth  = %f\n", histoPt->GetBinWidth(startPtBin),  histoPt->GetBinWidth(endPtBin));
    delete histoPt; histoPt = 0;
  }

  Int_t startVdBin = 0;
  Int_t endVdBin   = 0;
  TH1D *histoVd = addHist->Projection(2);
  if(histoVd){
    startVdBin = histoVd->FindFirstBinAbove(0);
    endVdBin   = histoVd->FindLastBinAbove(0);
    printf("startVdBin       = %i endVdBin       = %i\n", startVdBin, endVdBin);
    printf("startVdBinCentre = %f endVdBinCentre = %f\n", histoVd->GetBinCenter(startVdBin), histoVd->GetBinCenter(endVdBin));
    printf("startVdBinWidth  = %f endVdBinWidth  = %f\n", histoVd->GetBinWidth(startVdBin),  histoVd->GetBinWidth(endVdBin));
    delete histoVd; histoVd = 0;
  }

  Int_t startRunBin = 0;
  Int_t endRunBin   = 0;
  TH1D *histoRun = addHist->Projection(3);
  if(histoRun){
    startRunBin = histoRun->FindFirstBinAbove(0);
    endRunBin   = histoRun->FindLastBinAbove(0);
    printf("startRunBin       = %i endRunBin       = %i\n", startRunBin, endRunBin);
    printf("startRunBinCentre = %f endRunBinCentre = %f\n", histoRun->GetBinCenter(startRunBin), histoRun->GetBinCenter(endRunBin));
    printf("startRunBinWidth  = %f endRunBinWidth  = %f\n", histoRun->GetBinWidth(startRunBin),  histoRun->GetBinWidth(endRunBin));
    delete histoRun; histoRun = 0;
  }

  TObjArray *vdriftArray = new TObjArray();
  if(!vdriftArray) return -2;

  TObjArray *array = fCalibTime->GetHistoDrift();
  if(!array) return -3;

  TIterator *iterator = array->MakeIterator();
  if(!iterator) return -4;

  iterator->Reset();
  THnSparse *hist = NULL;
  while((hist = (THnSparseF*)iterator->Next())){
       
         //if(!hist) continue;
         hist->Print();
         hist->GetAxis(0)->SetRange(startTimeBin, endTimeBin);
         hist->GetAxis(1)->SetRange(startPtBin,   endPtBin);
         hist->GetAxis(0)->SetRange(startVdBin,   endVdBin);
         hist->GetAxis(3)->SetRange(startRunBin,  endRunBin);
         
	 TString name = hist->GetName();
         Int_t dim[4] = {0,1,2,3};
         THnSparse *newHist = hist->Projection(4,dim);
         newHist->SetName(name);
         vdriftArray->Add(newHist);
         
	 TGraphErrors *graph = AliTPCcalibBase::FitSlices(newHist,2,0,400,100,0.05,0.95, kTRUE);
         printf("name = %s graph = %i\n", name.Data(), graph==0);
         if(!graph || !graph->GetN()) continue;
         printf("name = %s graph = %i, N = %i\n", name.Data(), graph==0, graph->GetN());
         Int_t pos = name.Index("_");
         name = name(pos,name.Capacity()-pos);
         TString graphName = graph->ClassName();
         graphName+=name;
         graphName.ToUpper();
         graph->SetName(graphName);
         printf("name = %s\n", graphName.Data());
         vdriftArray->Add(graph);

         //Currently, AliSplineFits can not be given names...
         //AliSplineFit* fit=new AliSplineFit();
         //fit->SetGraph(graph);
         //fit->SetMinPoints(graph->GetN()+1);
         //fit->InitKnots(graph,2,0,0.001);
         //fit->SplineFit(0);
         //TString fiName=fit->ClassName();
         //fiName+=type;
         //fiName+=trigger;
         //fiName.ToUpper();
         //fit->SetName(fiName.Data());
         //printf("name=%s\n", fiName.Data());
         //vdriftArray->Add(fit);
  }
  
  THnSparse    *laserHist  = NULL;
  TGraphErrors *laserGraph = NULL;
  TString laserName = "";

  //Histograms and graphs for A side lasers
  laserHist = fCalibTime->GetHistVdriftLaserA(1);
  if(laserHist){
    
     laserName=laserHist->ClassName();
     laserName+="_MEAN_DRIFT_LASER_ALL_A";
     laserName.ToUpper();
     laserHist->SetName(laserName);
     vdriftArray->Add(laserHist);
     laserGraph=AliTPCcalibBase::FitSlices(laserHist,2,0,400,100,0.05,0.95, kTRUE);
     if(laserGraph && laserGraph->GetN()){
       laserName=laserGraph->GetName();
       laserName+="_MEAN_DRIFT_LASER_ALL_A";
       laserName.ToUpper();
       laserGraph->SetName(laserName);
       vdriftArray->Add(laserGraph);
     }
  }

  //Histograms and graphs for C side lasers
  laserHist=fCalibTime->GetHistVdriftLaserC(1);
  if(laserHist){
     laserName=laserHist->ClassName();
     laserName+="_MEAN_DRIFT_LASER_ALL_C";
     laserName.ToUpper();
     laserHist->SetName(laserName);
     vdriftArray->Add(laserHist);
     laserGraph=AliTPCcalibBase::FitSlices(laserHist,2,0,400,100,0.05,0.95, kTRUE);
     if(laserGraph && laserGraph->GetN()){
       laserName=laserGraph->GetName();
       laserName+="_MEAN_DRIFT_LASER_ALL_C";
       laserName.ToUpper();
       laserGraph->SetName(laserName);
       vdriftArray->Add(laserGraph);
     }
  }

  //Meatdata set in off-line...
  //AliCDBMetaData *metaData= new AliCDBMetaData();
  //metaData->SetObjectClassName("TObjArray");
  //metaData->SetResponsible("Dag Toppe Larsen");
  //metaData->SetBeamPeriod(1);
  //metaData->SetAliRootVersion("05-25-01"); //root version
  //metaData->SetComment("Calibration of the time dependence of the drift velocity due to pressure and temperature changes");
  //AliCDBId* id1=NULL;
  //if(end) id1=new AliCDBId("TPC/Calib/TimeDrift", runNumber, end);
  //else    id1=new AliCDBId("TPC/Calib/TimeDrift", runNumber, runNumber);
  //AliCDBStorage* gStorage = AliCDBManager::Instance()->GetStorage("local://$ALICE_ROOT/OCDB");
  //gStorage->Put(vdriftArray, (*id1), metaData);
  //printf("done runNumber=%i, end=%i\n", runNumber, end);

  static AliHLTReadoutList rdList(AliHLTReadoutList::kTPC);
  

  TFile *file = TFile::Open("vdrift.root", "RECREATE");
  vdriftArray->Write();
  file->Close();
  delete file;

  file = TFile::Open("calibTime.root", "RECREATE");
  fCalibTime->Write();
  file->Close();
  delete file;
 
  // the vdriftArray is pushed to the HLT-FXSsubscriber 
  PushToFXS( (TObject*)vdriftArray, "TPC", "TIMEDRIFT", &rdList );
 
  //Should array be deleted now?
  //  if(vdriftArray){
  //      vdriftArray.Clear();
  //    delete vdriftArray;
  //    vdriftArray=0;
  //  }
  
  return 0;
} 

