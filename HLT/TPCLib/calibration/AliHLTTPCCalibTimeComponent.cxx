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

#include "AliESDEvent.h"
#include "AliESDtrack.h"

#include "AliTPCcalibTime.h"
#include "AliTPCseed.h"
#include "AliTPCcalibCalib.h"

#include "TObjArray.h"
#include "TString.h"

#include "THnSparse.h"
#include "TGraphErrors.h"

#include <cstdlib>
#include <cerrno>

#include "AliHLTReadoutList.h"

ClassImp(AliHLTTPCCalibTimeComponent) // ROOT macro for the implementation of ROOT specific class methods

AliHLTTPCCalibTimeComponent::AliHLTTPCCalibTimeComponent()
  :
  fCalibTime(NULL),
  fESDevent(NULL),
  fESDtrack(NULL),
  fESDfriendTrack(NULL),
  fSeedArray(NULL),
  fMinPartition(5),
  fMaxPartition(0),
  fMinSlice(35),
  fMaxSlice(0),
  fSpecification(0) ,
  fEnableAnalysis(kTRUE)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTPCCalibTimeComponent::~AliHLTTPCCalibTimeComponent() {
// see header file for class documentation
}


const char* AliHLTTPCCalibTimeComponent::GetComponentID() {
// see header file for class documentation

  return "TPCCalibTime";
}

void AliHLTTPCCalibTimeComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list) {
// see header file for class documentation

  list.clear(); 
  list.push_back( kAliHLTDataTypeTObjArray|kAliHLTDataOriginTPC ); // output of TPCCalibSeedMaker
  list.push_back( kAliHLTDataTypeESDObject|kAliHLTDataOriginOut ); // output of global esd converter
}

AliHLTComponentDataType AliHLTTPCCalibTimeComponent::GetOutputDataType() {
// see header file for class documentation

  return AliHLTTPCDefinitions::fgkCalibCEDataType|kAliHLTDataOriginOut;
}

void AliHLTTPCCalibTimeComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ) {
// see header file for class documentation

  constBase = 20000;
  inputMultiplier = (2.0); // to be estimated
}

AliHLTComponent* AliHLTTPCCalibTimeComponent::Spawn() {
// see header file for class documentation

  return new AliHLTTPCCalibTimeComponent();
}  


Int_t AliHLTTPCCalibTimeComponent::ScanArgument( Int_t argc, const char** argv ) {
// see header file for class documentation

  Int_t iResult = 0;
  TString argument = "";
  TString parameter = "";

  if(!argc) return -EINVAL;

  argument = argv[iResult];
  
  if(argument.IsNull()) return -EINVAL;

  if( argument.CompareTo("-enable-analysis") == 0 ){
    HLTInfo( "Analysis before shipping data to FXS enabled." );
    fEnableAnalysis = kTRUE;
  }
  else {
    iResult = -EINVAL;
  }      
  return iResult;
}

Int_t AliHLTTPCCalibTimeComponent::InitCalibration() {
// see header file for class documentation
    
  if(fCalibTime) return EINPROGRESS;
  //fCalibTime = new AliTPCcalibTime();
  //fCalibTime = new AliTPCcalibTime("calibTime","time dependent Vdrift calibration",-2, 2, 1);

  return 0;
}

Int_t AliHLTTPCCalibTimeComponent::DeinitCalibration() {
// see header file for class documentation

  if(fCalibTime) delete fCalibTime; fCalibTime = NULL;

  return 0;
}

Int_t AliHLTTPCCalibTimeComponent::ProcessCalibration( const AliHLTComponentEventData& /*evtData*/,  AliHLTComponentTriggerData& /*trigData*/ ){
// see header file for class documentation

  if(GetFirstInputBlock( kAliHLTDataTypeSOR ) || GetFirstInputBlock( kAliHLTDataTypeEOR )) return 0;

  TObject *iter = NULL;

  //--------------- output over TObjArray of AliTPCseed objects (output of TPCSeedMaker) -------------------//
  
  for(iter = (TObject*)GetFirstInputObject(kAliHLTDataTypeTObjArray|kAliHLTDataOriginTPC); iter != NULL; iter = (TObject*)GetNextInputObject()){  
              
      if(GetDataType(iter) != (kAliHLTDataTypeTObjArray | kAliHLTDataOriginTPC)) continue;      
      fSeedArray = dynamic_cast<TObjArray*>(iter);      
   }

 
  //----------- loop over output of global esd converter ----------------//
 
  for(iter = (TObject*)GetFirstInputObject(kAliHLTDataTypeESDObject | kAliHLTDataOriginOut); iter != NULL; iter = (TObject*)GetNextInputObject()){   
      
    if(GetDataType(iter) != (kAliHLTDataTypeESDObject | kAliHLTDataOriginOut)) continue;
  
    AliHLTUInt8_t slice     = AliHLTTPCDefinitions::GetMinSliceNr(GetSpecification(iter)); 
    AliHLTUInt8_t partition = AliHLTTPCDefinitions::GetMinPatchNr(GetSpecification(iter));
      
    if( partition < fMinPartition ) fMinPartition = partition;
    if( partition > fMaxPartition ) fMaxPartition = partition;
    if( slice < fMinSlice ) fMinSlice = slice;
    if( slice > fMaxSlice ) fMaxSlice = slice;
            
    fESDevent = dynamic_cast<AliESDEvent*>(iter);
    fESDevent->CreateStdContent();
    
    //fESDevent->SetTimeStamp(1256910155);
    //fESDevent->SetRunNumber(84714);
    //HLTInfo("time stamp and event number -----: %d, %d\n", fESDevent->GetTimeStamp(), fESDevent->GetRunNumber());
    //fCalibTime->UpdateEventInfo(fESDevent);
      
    HLTDebug("# Seeds: %i\n", fSeedArray->GetEntriesFast());
    
    AliTPCcalibCalib *cal = new AliTPCcalibCalib();
      
    for(Int_t i=0; i<fSeedArray->GetEntriesFast(); i++){        
        
	AliTPCseed *seed = (AliTPCseed*)fSeedArray->UncheckedAt(i);
        if(!seed) continue; 
                    
        fESDtrack = fESDevent->GetTrack(i);
        fESDfriendTrack = const_cast<AliESDfriendTrack*>(fESDtrack->GetFriendTrack());
      
        cal->RefitTrack(fESDtrack, seed, GetBz());

        AliTPCseed *seedCopy = new AliTPCseed(*seed, kTRUE); 
        fESDtrack->AddCalibObject(seedCopy);      
        //fESDfriendTrack->AddCalibObject(seed);
      
        fCalibTime->ProcessAlignITS(fESDtrack,fESDfriendTrack);
    }
  } 
    
  if(!fCalibTime) {
    Int_t startTime=fESDevent->GetTimeStamp()-60*60*1;  //Start time one hour before first event, will make precise cuts later.
    Int_t   endTime=fESDevent->GetTimeStamp()+60*60*23; //End time 23 hours after first event.
    fCalibTime = new AliTPCcalibTime("calibTime","time dependent Vdrift calibration", startTime, endTime, 20*60);
    printf("fCalibTime=%i, startTime=%i, endTime=%i\n", fCalibTime!=0, startTime, endTime);
  }

  fCalibTime->UpdateEventInfo(fESDevent);
  fCalibTime->Process(fESDevent);
  
  
  fSpecification = AliHLTTPCDefinitions::EncodeDataSpecification( fMinSlice, fMaxSlice, fMinPartition, fMaxPartition );
  PushBack( (TObject*)fCalibTime, AliHLTTPCDefinitions::fgkCalibCEDataType | kAliHLTDataOriginOut, fSpecification);
  
  return 0;
}

Int_t AliHLTTPCCalibTimeComponent::ShipDataToFXS( const AliHLTComponentEventData& /*evtData*/,  AliHLTComponentTriggerData& /*trigData*/ ){
// see header file for class documentation

  HLTInfo("Shipping data to FXS...\n");

  if(fEnableAnalysis) fCalibTime->Analyze();

  THnSparse* addHist=fCalibTime->GetHistoDrift("all");
  if(!addHist) return -1;

  //Identifying used range of histogram
  Int_t startTimeBin=0;
  Int_t endTimeBin=0;
  TH1D* histoTime=addHist->Projection(0);
  if(histoTime) {
    startTimeBin=histoTime->FindFirstBinAbove(0);
    endTimeBin  =histoTime->FindLastBinAbove(0);
    printf("startTimeBin=%i endTimeBin=%i\n", startTimeBin, endTimeBin);
    printf("startTimeBinCentre=%f endTimeBinCentre=%f\n", histoTime->GetBinCenter(startTimeBin), histoTime->GetBinCenter(endTimeBin));
    printf("startTimeBinWidth=%f endTimeBinWidth=%f\n", histoTime->GetBinWidth(startTimeBin), histoTime->GetBinWidth(endTimeBin));
    delete histoTime;
    histoTime=0;
  }

  Int_t startPtBin=0;
  Int_t endPtBin=0;
  TH1D* histoPt=addHist->Projection(1);
  if(histoPt) {
    startPtBin=histoPt->FindFirstBinAbove(0);
    endPtBin  =histoPt->FindLastBinAbove(0);
    printf("startPtBin=%i endPtBin=%i\n", startPtBin, endPtBin);
    printf("startPtBinCentre=%f endPtBinCentre=%f\n", histoPt->GetBinCenter(startPtBin), histoPt->GetBinCenter(endPtBin));
    printf("startPtinWidth=%f endPtBinWidth=%f\n", histoPt->GetBinWidth(startPtBin), histoPt->GetBinWidth(endPtBin));
    delete histoPt;
    histoPt=0;
  }

  Int_t startVdBin=0;
  Int_t endVdBin=0;
  TH1D* histoVd=addHist->Projection(2);
  if(histoVd) {
    startVdBin=histoVd->FindFirstBinAbove(0);
    endVdBin  =histoVd->FindLastBinAbove(0);
    printf("startVdBin=%i endVdBin=%i\n", startVdBin, endVdBin);
    printf("startVdBinCentre=%f endVdBinCentre=%f\n", histoVd->GetBinCenter(startVdBin), histoVd->GetBinCenter(endVdBin));
    printf("startVdBinWidth=%f endVdBinWidth=%f\n", histoVd->GetBinWidth(startVdBin), histoVd->GetBinWidth(endVdBin));
    delete histoVd;
    histoVd=0;
  }

  TH1D* histoRun=addHist->Projection(3);
  Int_t startRunBin=0;
  Int_t endRunBin=0;
  if(histoRun) {
    startRunBin=histoRun->FindFirstBinAbove(0);
    endRunBin  =histoRun->FindLastBinAbove(0);
    printf("startRunBin=%i endRunBin=%i\n", startRunBin, endRunBin);
    printf("startRunBinCentre=%f endRunBinCentre=%f\n", histoRun->GetBinCenter(startRunBin), histoRun->GetBinCenter(endRunBin));
    printf("startRunBinWidth=%f endRunBinWidth=%f\n", histoRun->GetBinWidth(startRunBin), histoRun->GetBinWidth(endRunBin));
    delete histoRun;
    histoRun=0;
  }

  TObjArray* vdriftArray = new TObjArray();
  if(!vdriftArray) return -2;

  TObjArray* array=fCalibTime->GetHistoDrift();
  if(!array) return -3;

  TIterator* iterator = array->MakeIterator();
  if(!iterator) return -4;

  iterator->Reset();
  THnSparse* hist=NULL;
  while((hist=(THnSparseF*)iterator->Next())){
    if(!hist) continue;
    hist->Print();
    hist->GetAxis(0)->SetRange(startTimeBin, endTimeBin);
    hist->GetAxis(1)->SetRange(startPtBin,   endPtBin);
    hist->GetAxis(0)->SetRange(startVdBin,   endVdBin);
    hist->GetAxis(3)->SetRange(startRunBin,  endRunBin);
    TString name=hist->GetName();
    Int_t dim[4]={0,1,2,3};
    THnSparse* newHist=hist->Projection(4,dim);
    newHist->SetName(name);
    vdriftArray->Add(newHist);
    TGraphErrors* graph=AliTPCcalibBase::FitSlices(newHist,2,0,400,100,0.05,0.95, kTRUE);
    printf("name=%s graph=%i\n", name.Data(), graph==0);
    if(!graph || !graph->GetN()) continue;
    printf("name=%s graph=%i, N=%i\n", name.Data(), graph==0, graph->GetN());
    Int_t pos=name.Index("_");
    name=name(pos,name.Capacity()-pos);
    TString graphName=graph->ClassName();
    graphName+=name;
    graphName.ToUpper();
    graph->SetName(graphName);
    printf("name=%s\n", graphName.Data());
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
  THnSparse* laserHist=NULL;
  TGraphErrors* laserGraph=NULL;
  TString laserName="";

  //Histograms and graphs for A side lasers
  laserHist=fCalibTime->GetHistVdriftLaserA(1);
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
  PushToFXS( (TObject*)vdriftArray, "TPC", "Time", rdList.Buffer() );
  //PushToFXS( (TObject*)vdriftArray, "TPC", "Time");

  //Should array be deleted now?
  //  if(vdriftArray){
  //      vdriftArray.Clear();
  //    delete vdriftArray;
  //    vdriftArray=0;
  //  }
  
  return 0;
} 

