/*
***********************************************************
  Manager for event plane corrections framework
  Contact: Jaap Onderwaater, j.onderwaater@gsi.de, jacobus.onderwaater@cern.ch
  Instructions in AddTask_EPcorrectionsExample.C
  2014/12/10
  *********************************************************
*/

#include "AliTrigger.h"
//#include "AliSysInfo.h"
#include <iostream>

#include <TROOT.h>
#include <TTimeStamp.h>
#include <TStopwatch.h>
#include <AliInputEventHandler.h>
#include <AliESDInputHandler.h>
#include <AliAODInputHandler.h>
#include <AliAnalysisManager.h>
#include <AliCentrality.h>
#include <AliESDEvent.h>
#include "AliEventPlaneESDVarManager.h"
#include "AliEventPlaneCuts.h"
#include "AliEventPlaneVarManager.h"
#include "AliEventPlaneManager.h"
#include "AliEventPlaneHelper.h"
#include "AliEventPlaneHistos.h"
//#include "AliEventPlaneReducedVarManager.h"
//#include "AliReducedEvent.h"
//#include "trainsimulator/localHelper.h"
// make a change
#include "AliAnalysisTaskEventPlaneCalibration.h"

//#ifdef ALIEVENTPLANEREDUCEDVARMANAGER_H
//#define FILL AliEventPlaneReducedVarManager
//#endif

//#ifdef ALIEVENTPLANEESDVARMANAGER_H
#define FILL AliEventPlaneESDVarManager
//#endif

#ifdef ALIEVENTPLANEVARMANAGER_H
#define VAR AliEventPlaneVarManager
#endif

#ifdef ALIEVENTPLANEHISTOS_H
#define HIST AliEventPlaneHistos
#endif

#ifdef ALIEVENTPLANEHELPER_H
#define HELPER AliEventPlaneHelper
#endif

// change 2

using std::cout;
using std::endl;


ClassImp(AliAnalysisTaskEventPlaneCalibration)


//_________________________________________________________________________________
AliAnalysisTaskEventPlaneCalibration::AliAnalysisTaskEventPlaneCalibration() :
 AliAnalysisTaskSE(),
 fRunLightWeight(kFALSE),
 fCalibrateByRun(kTRUE),
 fUseFriendEvent(kFALSE),
 fTriggerMask(0),
 fFillTPC(kFALSE),
 fFillVZERO(kFALSE),
 fFillTZERO(kFALSE),
 fFillZDC(kFALSE),
 fFillFMD(kFALSE),
 fChannelEqualization(kFALSE),
 fRecenterQvec(kFALSE),
 fRotateQvec(kFALSE),
 fTwistQvec(kFALSE),
 fScaleQvec(kFALSE),
 fInitialized(kFALSE),
 fTree(0x0),
 fHistosFile(0x0),
 fQvectorFile(0x0),
 fListHistos(),
 fListHistosQA(),
 fFriendEventNo(0),
 fEventPlaneManager(0x0),
 fEventCuts(0x0),
 fOffset(0),
 fNfiles(1),
 fIsAOD(kFALSE),
 fIsESD(kFALSE),
 fIsReduced(kFALSE)
{
  //
  // Default constructor
  //
}

//_________________________________________________________________________________
AliAnalysisTaskEventPlaneCalibration::AliAnalysisTaskEventPlaneCalibration(const char* name) :
 AliAnalysisTaskSE(name),
 fRunLightWeight(kFALSE),
 fCalibrateByRun(kTRUE),
 fUseFriendEvent(kFALSE),
 fTriggerMask(0),
 fFillTPC(kFALSE),
 fFillVZERO(kFALSE),
 fFillTZERO(kFALSE),
 fFillZDC(kFALSE),
 fFillFMD(kFALSE),
 fChannelEqualization(kFALSE),
 fRecenterQvec(kFALSE),
 fRotateQvec(kFALSE),
 fTwistQvec(kFALSE),
 fScaleQvec(kFALSE),
 fInitialized(kFALSE),
 fTree(0x0),
 fHistosFile(0x0),
 fQvectorFile(0x0),
 fListHistos(),
 fListHistosQA(),
 fFriendEventNo(0),
 fEventPlaneManager(0x0),
 fEventCuts(0x0),
 fOffset(0),
 fNfiles(1),
 fIsAOD(kFALSE),
 fIsESD(kFALSE),
 fIsReduced(kFALSE)
{
  //
  // Constructor
  //
  

  //localHelper::GetInstance()->DefineInput(0,AliReducedEvent::Class());
  //localHelper::GetInstance()->DefineOutput(1, TList::Class());   // Calibration histograms
  //localHelper::GetInstance()->DefineOutput(2, TTree::Class());   // Calibrated qvector tree
  //localHelper::GetInstance()->DefineOutput(3, TList::Class());   // QA histograms
  //if(IsReduced()) DefineInput(0,AliReducedEvent::Class());
  //else DefineInput(0,TChain::Class());
  DefineInput(0,TChain::Class());
  DefineOutput(1, TList::Class());   // Calibration histograms
  DefineOutput(2, TTree::Class());   // Calibrated qvector tree
  DefineOutput(3, TList::Class());   // QA histograms

  fListHistos.SetName("CalibrationHistograms");
  fListHistos.SetOwner(kFALSE);
  fListHistosQA.SetName("CalibrationQA");
  fListHistosQA.SetOwner(kFALSE);
}


//_________________________________________________________________________________
void AliAnalysisTaskEventPlaneCalibration::UserCreateOutputObjects()
{
  //
  // Add all histogram manager histogram lists to the output TList
  //

  Float_t values[AliEventPlaneVarManager::kNVars]={0};
  if (!fListHistos.IsEmpty()) return; //already initialised

  fTree = new TTree("Qvectors","Qvector values");
  fTree->SetDirectory(0);

  //TObjArray* histObj = HIST::Instance()->HistList();

  //for(Int_t i=0; i<histObj->GetEntries(); ++i) {
  //  cout<<"create stuff"<<endl;
  //    fListHistos.Add(static_cast<THashList*>(histObj->At(i)));
  //}

  AliEventPlaneConfiguration* EPconf = 0x0;
  //TClonesArray* QvecList[fEventPlaneManager->NEventPlaneDetectors()];
  TClonesArray* QvecList[50];
  for(Int_t idet=0; idet<AliEventPlaneManager::kNdetectors; idet++){
    TClonesArray* epConfList=fEventPlaneManager->GetEventPlaneConfigurations(idet);
    TIter nextEPconf(epConfList);
    while((EPconf=static_cast<AliEventPlaneConfiguration*>(nextEPconf()))) {
      if(!EPconf) continue;   
      //QvecList[EPconf->GlobalIndex()]  = fEventPlaneManager->GetQvectors(EPconf->GlobalIndex());
      QvecList[EPconf->GlobalIndex()]  = EPconf->Qvectors();
      QvecList[EPconf->GlobalIndex()]->SetName(EPconf->EventPlaneDetectorName());
      fTree->Branch(EPconf->EventPlaneDetectorName(),&QvecList[EPconf->GlobalIndex()],256000,1);
    }
  }



  fListHistos.SetName("CalibrationHistograms");
  fListHistos.SetOwner();
  fListHistosQA.SetName("CalibrationQA");
  fListHistosQA.SetOwner();

  PostData(1, &fListHistos);
  PostData(2, fTree);
  PostData(3, &fListHistosQA);
  //localHelper::GetInstance()->PostData(1, &fListHistos);
  //localHelper::GetInstance()->PostData(2, fTree);
  //localHelper::GetInstance()->PostData(3, &fListHistosQA);


}


  
//________________________________________________________________________________________________________
void AliAnalysisTaskEventPlaneCalibration::UserExec(Option_t *){
  //
  // Main loop. Called for every event
  //

  //AliReducedEvent* reducedEvent=0x0;
  //AliReducedEvent* event = dynamic_cast<AliReducedEvent*>(localHelper::GetInstance()->GetInputData(0));
  AliVEvent* event = InputEvent();
  //TObject* event = (TObject*) Vevent;
  //TObject* event = (TObject*) reducedEvent;

   //Was event selected ?
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  Bool_t isESD = (event->IsA()==AliESDEvent::Class());
  Bool_t isAOD = (event->IsA()==AliAODEvent::Class());
  if(isESD) SetInputESD();
  if(isAOD) SetInputAOD();
  
  AliESDEvent* esdEvent = 0x0;
  if(isESD) esdEvent = static_cast<AliESDEvent*>(event);
  AliAODEvent* aodEvent = 0x0;
  if(isAOD) aodEvent = static_cast<AliAODEvent*>(event);
  
  UInt_t isSelected = AliTrigger::kAny;
  if(fTriggerMask==0) fTriggerMask=AliTrigger::kAny;
  if(!IsReduced()){
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  if(inputHandler){
    if((isESD && inputHandler->GetEventSelection()) || isAOD){
      isSelected = inputHandler->IsEventSelected();
      isSelected&=fTriggerMask;
    }
  }}
  
  AliEventPlaneManager* EPmanager = EventPlaneManager();

  EPmanager->ClearEvent();
  fRunLightWeight = EPmanager->RunLightWeight();



  Float_t values[VAR::kNVars]={0};
  for(Int_t i=0; i<VAR::kNVars; ++i) values[i]=-9999.;
  
  Int_t currentRunNumber = 0;
  FILL::FillEventInfo(event, values);
  //values[VAR::kIsPhysicsSelection] = ( isSelected==0 ? 0 : 1 );
  EPmanager->SetEventStatus( isSelected==0 ? 0 : 1 );
  //std::cout<<"!!!!  "<<values[VAR::kIsPhysicsSelection]<<std::endl;
  currentRunNumber = values[VAR::kRunNo];
  //if(IsReduced()) currentRunNumber = reducedEvent->RunNo();
  //else currentRunNumber = Vevent->GetRunNumber();

  if(!fInitialized){
  

  // define histograms -----------------------------------------------
  TString histClasses = "";

  if(!fRunLightWeight) histClasses += "Event_NoCuts;Event_AnalysisEvents;Event_CalibrationEvents;Event_MultCutOuts;";
  if(!fRunLightWeight) histClasses += "OfflineTriggers_NoCuts;OfflineTriggers_WithCuts;";


  //      /// Set the correlation detector indices (this prevents having to do many string operations)
  Int_t cor1=-1, cor2=-1;
  
 
 AliEventPlaneConfiguration* EPconf = 0x0;
 for(Int_t idet=0; idet<AliEventPlaneManager::kNdetectors; idet++){
   AliEventPlaneConfiguration* EPconf2 = 0x0;
   TClonesArray* epConfList=EPmanager->GetEventPlaneConfigurations(idet);
   for(Int_t iconf=0; iconf<epConfList->GetEntriesFast(); iconf++){
   EPconf = (AliEventPlaneConfiguration*) epConfList->At(iconf);
    if(!EPconf) continue;        


       TString detStrCor1 = EPconf->CorrelationDetectorName(0);
       TString detStrCor2 = EPconf->CorrelationDetectorName(1);

          for(Int_t idet2=0; idet2<AliEventPlaneManager::kNdetectors; idet2++){
               TClonesArray* epConfList2=EPmanager->GetEventPlaneConfigurations(idet2);
               for(Int_t iconf2=0; iconf2<epConfList2->GetEntriesFast(); iconf2++){
               EPconf2 = (AliEventPlaneConfiguration*) epConfList2->At(iconf2);


                 if(detStrCor1.EqualTo(EPconf2->EventPlaneDetectorName())) cor1 = EPconf2->GlobalIndex();
                 if(detStrCor2.EqualTo(EPconf2->EventPlaneDetectorName())) cor2 = EPconf2->GlobalIndex();
       }    
  
     }


     EPconf->SetCorrelationDetectorIndices(cor1, cor2);
   }

  }
  



  //cout<<"RUN  "<<currentRunNumber<<endl;
  
  InitializeCalibrationHistograms(currentRunNumber, EPmanager);

  for(Int_t idet=0; idet<AliEventPlaneManager::kNdetectors; idet++){ 
    TClonesArray* epConfList=EPmanager->GetEventPlaneConfigurations(idet);
    for(Int_t iconf=0; iconf<epConfList->GetEntriesFast(); iconf++){
    EPconf = (AliEventPlaneConfiguration*) epConfList->At(iconf);
     if(!EPconf) continue;
        if(idet==AliEventPlaneManager::kTPC)   fFillTPC = kTRUE;
        if(idet==AliEventPlaneManager::kVZERO) fFillVZERO = kTRUE;
        if(idet==AliEventPlaneManager::kTZERO) fFillTZERO = kTRUE;
        if(idet==AliEventPlaneManager::kZDC)   fFillZDC = kTRUE;
        if(idet==AliEventPlaneManager::kFMD)   fFillFMD = kTRUE;
        EPconf->CreateCalibrationHistograms();
        if(idet==AliEventPlaneManager::kTPC){
          //histClasses += "Tracks_"+EPconf->EventPlaneDetectorName()+";";
          histClasses += "TrackQA_"+EPconf->EventPlaneDetectorName()+";";
        }
    }
  }


  DefineHistograms(histClasses.Data(), currentRunNumber, EPmanager);

fInitialized = kTRUE;
cout<<"initialized"<<endl;
}


  
  //if(IsReduced()) VAR::FillEventInfo(reducedEvent, values);
  //else VAR::FillEventInfo(Vevent, values);
    
    
      HIST::Instance()->FillHistClass("Event_NoCuts", values);
      for(UShort_t ibit=0; ibit<64; ++ibit) {
        //VAR::FillEventOfflineTriggers(ibit, reducedEvent, values);
        HIST::Instance()->FillHistClass("OfflineTriggers_NoCuts", values);
      }


    ///// FILL REDUCED DETECTOR INFO FOR ALL EVENTFRIEND DETECTORS ////

//cout<<"RUNNING CALIBRATION"<<endl;
 
    //VZERO  
    if(fFillVZERO) FILL::FillVZERO(EPmanager,event);

    //cout<<"crash 1"<<endl;
    //TPC
    if(fFillTPC)     FILL::FillTPC(EPmanager, event, values);
    //if(fFillTPC&&IsESD())     EPmanager->FillTPC(*esdEvent, values);
    //if(fFillTPC&&IsReduced()) EPmanager->FillTPC(reducedEvent, values);
    
    ////cout<<"crash 2"<<endl;
    ////TZERO  
    if(fFillTZERO) FILL::FillTZERO(EPmanager,event);
    //  
    ////cout<<"crash 3"<<endl;
    //ZDC  
    if(fFillZDC) FILL::FillZDC(EPmanager,event);
    //  
    ////cout<<"crash 4"<<endl;
    ////FMD  
    //if(fFillFMD) FILL::FillFMD(EPmanager,event);

    //cout<<"crash 5"<<endl;
    EPmanager->ClearUnusedDetectors();

    // use only selected triggers for event plane calibration averages
    if(IsEventSelected(values)){//&&TriggerSelected(event)) {
        HIST::Instance()->FillHistClass("Event_AnalysisEvents", values);
        if(isSelected) HIST::Instance()->FillHistClass("Event_CalibrationEvents", values);
        //for(UShort_t ibit=0; ibit<64; ++ibit) {
          //VAR::FillEventOfflineTriggers(ibit, event, values);
          //HIST::Instance()->FillHistClass("OfflineTriggers_WithCuts", values);
        //}

    //cout<<"crash 6"<<endl;

    // Get raw Qvectors
    EPmanager->GetQvector(0,0,values);

    //cout<<"crash 6a"<<endl;
    // Fill raw Qvector histograms
    //FillQvecHistograms(values, 0);
    //cout<<"crash 6b"<<endl;
    //FillRPcorrelationHistograms(values,  0);
    EPmanager->FillCorrelationHistograms(values[VAR::kCentVZERO],  0);
    //cout<<"crash 6c"<<endl;

    // Calibrate channels
    EPmanager->CalibrateChannels(values);
          
    //cout<<"crash 7"<<endl;
    // Get Qvectors with calibrated channels
    EPmanager->GetQvector(kFALSE,kTRUE,values);
    EPmanager->FillCorrelationHistograms(values[VAR::kCentVZERO],  1);

    // Fill equalized Qvector histograms
    //FillQvecHistograms(values, 1);
    //FillRPcorrelationHistograms(values,  1);

    // Recenter Qvectors
    EPmanager->RecenterQvec(values);
    EPmanager->FillCorrelationHistograms(values[VAR::kCentVZERO],  2);

    EPmanager->U2nTwistQvec(values, VAR::kCentVZERO);
    // Fill recentered Qvector histograms
    //FillQvecHistograms(values, 2);
    //FillRPcorrelationHistograms(values,  2);

    //cout<<"crash 8"<<endl;
      //Fill channel multiplicity histogram
      //FillMultHistograms(values);
      //cout<<"crash 9"<<endl;

      // Fill rotated Qvector histograms
 //--    EPmanager->RotateQvec();
      //FillQvecHistograms(values, 3);
      //FillRPcorrelationHistograms(values,  3);

      //cout<<"ENTRIES   "<<EPmanager->EventPlaneConfiguration("FMDAcor")->CalibrationHistogramQ(0,2,1)->GetEntries()<<endl;

 //--    // Twist Qvectors
 //--    EPmanager->U2nTwistQvec(values, VAR::kCentVZERO);
      //EPmanager->TwoDetectorCorrelationTwistQvec(values, VAR::kCentVZERO);
 //--    EPmanager->ThreeDetectorCorrelationTPCTwistQvec(values, VAR::kCentVZERO);

      //cout<<"crash 10"<<endl;
      // Fill twisted Qvector histograms
      //FillQvecHistograms(values, 4);
      //cout<<"crash 11"<<endl;
      //FillRPcorrelationHistograms(values,  4);
      //cout<<"crash 12"<<endl;


      // Fill rescaled Qvector histograms
 //--    EPmanager->U2nRescalingQvec(values, VAR::kCentVZERO);
 //--    EPmanager->ThreeDetectorCorrelationTPCRescalingQvec(values,VAR::kCentVZERO);
      //FillQvecHistograms(values, 5);
      //FillRPcorrelationHistograms(values,  5);
      

    }  // end if trigger selection
      

    //if(gkCreateUpdatedEP) {
    fTree->Fill();
      //fHistosFile->cd();
    //}

    EPmanager->ClearEvent();
      
  }  // end loop over events

  
//__________________________________________________________________
void AliAnalysisTaskEventPlaneCalibration::FinishTaskOutput()
{
  //
  // Finish Task 
  //
  TObjArray* histObj = HIST::Instance()->HistList();

  for(Int_t i=0; i<histObj->GetEntries(); ++i) {
    TString out = histObj->At(i)->GetName();

    fListHistosQA.Add(static_cast<THashList*>(histObj->At(i)));
  }


 for(Int_t idet=0; idet<AliEventPlaneManager::kNdetectors; idet++){
   AliEventPlaneConfiguration* EPconf = 0x0;
   TClonesArray* epConfList=EventPlaneManager()->GetEventPlaneConfigurations(idet);
   if(epConfList->GetEntriesFast()>0){
     for(Int_t iconf=0; iconf<epConfList->GetEntriesFast(); iconf++){
       THashList* detector = new THashList();
       THashList* detectorM = new THashList();
       THashList* detectorC = new THashList();
       THashList* detectorMqa = new THashList();
       THashList* detectorCqa = new THashList();
       EPconf = (AliEventPlaneConfiguration*) epConfList->At(iconf);
       if(!EPconf) continue;
       for(Int_t is=0; is<(EPconf->CalibrationStep()+1); ++is){
         for(Int_t ih=fgkEPMinHarmonics; ih<=fgkEPMaxHarmonics; ++ih) {
           for(Int_t ic=0; ic<2; ++ic){
           detector->Add(EPconf->CalibrationHistogramQ(is,ih,ic));
           }
         }
         detector->Add(EPconf->CalibrationHistogramE(is));
       }

       for(Int_t ih=fgkEPMinHarmonics; ih<=fgkEPMaxHarmonics; ++ih){
         if(EPconf->U2nHistogram(ih,0)){
           detector->Add(EPconf->U2nHistogram(ih,0));
           detector->Add(EPconf->U2nHistogram(ih,1));
         }
       }
       if(EPconf->U2nHistogramE()) detector->Add(EPconf->U2nHistogramE());

       detector->SetName(Form("Qvec%s",EPconf->EventPlaneDetectorName().Data()));
       fListHistos.Add(detector);


       if(EPconf->doChannelEqualization()){
         for(Int_t is=0; is<3; ++is){
           if(EPconf->EqualizationHistogramM(is)&&is==0) detectorM->Add(EPconf->EqualizationHistogramM(is));
           if(EPconf->EqualizationHistogramE(is)&&is==0) detectorM->Add(EPconf->EqualizationHistogramE(is));
           if(EPconf->EqualizationHistogramM(is)&&is!=0) detectorMqa->Add(EPconf->EqualizationHistogramM(is));
           if(EPconf->EqualizationHistogramE(is)&&is!=0) detectorMqa->Add(EPconf->EqualizationHistogramE(is));
         }
         detectorM->SetName(Form("Mult%s",EPconf->EventPlaneDetectorName().Data()));
         detectorMqa->SetName(Form("Mult%s",EPconf->EventPlaneDetectorName().Data()));
         fListHistos.Add(detectorM);
         fListHistosQA.Add(detectorMqa);
       }



      for(Int_t is=0; is<(EPconf->CalibrationStep()+1); ++is){
        for(Int_t icomb=0; icomb<3; ++icomb){ 
          for(Int_t icomp=0; icomp<4; ++icomp){ 
            for(Int_t ih=fgkEPMinHarmonics; ih<=fgkEPMaxHarmonics; ++ih){ 
              if(EPconf->CorrelationProf(is,icomb,ih,icomp)&&is<=2) detectorC->Add(EPconf->CorrelationProf(is,icomb,ih,icomp));
              if(EPconf->CorrelationProf(is,icomb,ih,icomp)&&is>2) detectorCqa->Add(EPconf->CorrelationProf(is,icomb,ih,icomp));
              if(EPconf->CorrelationEpProf(is,icomb,ih,icomp)) detectorCqa->Add(EPconf->CorrelationEpProf(is,icomb,ih,icomp));
      }}}}
      detectorC->SetName(Form("Correlations%s",EPconf->EventPlaneDetectorName().Data()));
      detectorCqa->SetName(Form("Correlations%s",EPconf->EventPlaneDetectorName().Data()));
      fListHistos.Add(detectorC);
      fListHistosQA.Add(detectorCqa);

     }
   }
 }



}



//__________________________________________________________________
Bool_t AliAnalysisTaskEventPlaneCalibration::IsEventSelected(Float_t* values) {
  if(!fEventCuts) return kTRUE;
  return fEventCuts->IsSelected(values);
}


//__________________________________________________________________
void AliAnalysisTaskEventPlaneCalibration::InitializeCalibrationHistograms(Int_t currentRunNumber, AliEventPlaneManager* EPmanager){


  AliEventPlaneConfiguration* EPconf=0x0;
  for(Int_t idet=0; idet<AliEventPlaneManager::kNdetectors; idet++){
    TClonesArray* epConfList=EPmanager->GetEventPlaneConfigurations(idet);
    TIter nextEPconf(epConfList);
    while((EPconf=static_cast<AliEventPlaneConfiguration*>(nextEPconf()))) {
      if(!EPconf) continue;

      // Initialize equalization histograms
      //if(EPconf->doChannelEqualization()&&EPconf->CalibrationStep()>0)
        InitEqualizationHistograms(currentRunNumber, EPconf); 
  
      // Initialize recentering histograms
      //if(EPconf->CalibrationStep()>1&&EPconf->doRecentering())
        InitQvecCalibrationHistograms(currentRunNumber, EPconf); 
  
      // Initialize correlation histograms
      //if(EPconf->CalibrationStep()>2&&EPconf->TwistAndScalingMethod()==2)
        //InitQvecCorrelationHistograms(currentRunNumber, EPconf); 
      //
      // Initialize double harmonic histograms
      //if(EPconf->CalibrationStep()>2&&(EPconf->TwistAndScalingMethod()==0||EPconf->TwistAndScalingMethod()==1))
        //InitQvec2nHistograms(currentRunNumber, EPconf); 
    }
  }


}

//__________________________________________________________________
void AliAnalysisTaskEventPlaneCalibration::InitEqualizationHistograms(Int_t runNo, AliEventPlaneConfiguration* epConf) {
  //
  //  Initialize the calibration file and set the calibration histograms for a given detector
  //
  TString path = epConf->EqualizationHistPath();

  if(path!="") {
   cout << "Loading the "<<epConf->EventPlaneDetectorName()<<" channel multiplicities for run " << runNo << " from path " ;
   if(fCalibrateByRun) path=Form(path.Data(), runNo);
   else path+="/allCalibrationHistogramsMerged.root";
   cout << Form(path.Data(),runNo) << endl;
   epConf->ConnectInputMultiplicityHistograms(path);
  }
}



//_____________________________________________________________________
void AliAnalysisTaskEventPlaneCalibration::InitQvecCalibrationHistograms(Int_t runNo, AliEventPlaneConfiguration* epConf) {

  TString path = epConf->RecenteringHistPath();

  if(path!="") {
   TString det = epConf->EventPlaneDetectorName();
   cout << "Loading the "<<det<<" calibration for run " << runNo << " from path ";
   if(fCalibrateByRun) path=Form(path.Data(), runNo);
   else path+="/allCalibrationHistogramsMerged.root";
   cout << Form(path.Data(),runNo) << endl;
   epConf->ConnectInputCalibrationHistograms(path);
 }
}



//_______________________



////__________________________________________________________________
//Bool_t AliAnalysisTaskEventPlaneCalibration::TriggerSelected(AliReducedEvent* event) {
//  //
//  // Trigger selection for the events used to calculate calibration and recentering averages
//  //
//  //if(!(event->TriggerMask() & (ULong64_t(1)<<1))) return kFALSE;  // min bias for calibration
//  return kTRUE;
//}
////__________________________________________________________________
//Bool_t AliAnalysisTaskEventPlaneCalibration::TriggerSelected2(AliReducedEvent* event) {
//  //
//  // 
//  //
//  if(!(event->TriggerMask() & ((ULong64_t(1)<<4) | (ULong64_t(1)<<7)))) return kFALSE;  // Semi-central OR central trigger 
//	  //if(!(event->TriggerMask() & (ULong64_t(1)<<4))) return kFALSE;  // central trigger
//  return kTRUE;
//}




//__________________________________________________________________
void AliAnalysisTaskEventPlaneCalibration::DefineHistograms(const Char_t* histClasses, Int_t runNumber, AliEventPlaneManager* EPmanager) {
  //
  // define the histograms
  //
  cout << "Defining histograms ..." << endl;
  cout << "histogram classes: " << histClasses << endl;
  
  //fHistosFile=new TFile(output,"RECREATE");
  
  TString classesStr(histClasses);
  TObjArray* arr=classesStr.Tokenize(";");
  
  const Int_t kNRunBins = 3000;
  Double_t runHistRange[2] = {137000.,140000.};
  
  for(Int_t iclass=0; iclass<arr->GetEntries(); ++iclass) {
    TString classStr = arr->At(iclass)->GetName();
    cout << "hist class: " << classStr.Data() << endl;
    
    // Event wise histograms
    if(classStr.Contains("Event")) {
      HIST::Instance()->AddHistClass(classStr.Data());
      HIST::Instance()->AddHistogram(classStr.Data(),"RunNo","Run numbers;Run", kFALSE, kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo);
      HIST::Instance()->AddHistogram(classStr.Data(),"BC","Bunch crossing;BC", kFALSE,3000,0.,3000.,VAR::kBC);
      HIST::Instance()->AddHistogram(classStr.Data(),"IsPhysicsSelection","Physics selection flag;;", kFALSE,
		   2,-0.5,1.5,VAR::kIsPhysicsSelection, 0,0.0,0.0,VAR::kNothing, 0,0.0,0.0,VAR::kNothing, "off;on");
      
      HIST::Instance()->AddHistogram(classStr.Data(),"VtxZ","Vtx Z;vtx Z (cm)", kFALSE,300,0.0,0.0,VAR::kVtxZ);
      //HIST::Instance()->AddHistogram(classStr.Data(),"VtxZ","Vtx Z;vtx Z (cm)", kFALSE,300,-15.,15.,VAR::kVtxZ);
      HIST::Instance()->AddHistogram(classStr.Data(),"VtxX","Vtx X;vtx X (cm)", kFALSE,300,-1.,1.,VAR::kVtxX);
      HIST::Instance()->AddHistogram(classStr.Data(),"VtxY","Vtx Y;vtx Y (cm)", kFALSE,300,-1.,1.,VAR::kVtxY);


      if(!fRunLightWeight) HIST::Instance()->AddHistogram(classStr.Data(),"CentVZEROvsCentSPD","Centrality(VZERO);centrality VZERO (percents);centrality SPD (percents)", kFALSE,
                   100, 0.0, 100.0, VAR::kCentVZERO, 100, 0.0, 100.0, VAR::kCentSPD);
      if(!fRunLightWeight) HIST::Instance()->AddHistogram(classStr.Data(),"CentTPCvsCentSPD","Centrality(TPC);centrality TPC (percents);centrality SPD (percents)", kFALSE,
                   100, 0.0, 100.0, VAR::kCentTPC, 100, 0.0, 100.0, VAR::kCentSPD);
      if(!fRunLightWeight) HIST::Instance()->AddHistogram(classStr.Data(),"CentTPCvsCentVZERO","Centrality(TPC);centrality TPC (percents);centrality VZERO (percents)", kFALSE,
                   100, 0.0, 100.0, VAR::kCentTPC, 100, 0.0, 100.0, VAR::kCentVZERO);
      if(!fRunLightWeight) HIST::Instance()->AddHistogram(classStr.Data(),"CentTPCvsCentZDC","Centrality(TPC);centrality TPC (percents);centrality ZDC (percents)", kFALSE,
                   100, 0.0, 100.0, VAR::kCentTPC, 100, 0.0, 100.0, VAR::kCentZDC);
      if(!fRunLightWeight) HIST::Instance()->AddHistogram(classStr.Data(),"CentZDCvsCentVZERO","Centrality(ZDC);centrality ZDC (percents);centrality VZERO (percents)", kFALSE,
                   100, 0.0, 100.0, VAR::kCentZDC, 100, 0.0, 100.0, VAR::kCentVZERO);
      if(!fRunLightWeight) HIST::Instance()->AddHistogram(classStr.Data(),"CentZDCvsCentSPD","Centrality(ZDC);centrality ZDC (percents);centrality SPD (percents)", kFALSE,
                   100, 0.0, 100.0, VAR::kCentZDC, 100, 0.0, 100.0, VAR::kCentSPD);

      if(!fRunLightWeight) HIST::Instance()->AddHistogram(classStr.Data(),"MultVZEROvsCentVZERO","Multiplicity;multiplicity VZERO;VZERO centrality", kFALSE,
                   100, 0.0, 32000.0, VAR::kVZEROTotalMult, 100,0.,100., VAR::kCentVZERO);
      if(!fRunLightWeight) HIST::Instance()->AddHistogram(classStr.Data(),"MultSPDvsCentPSD","Multiplicity;SPD tracklets;SPD centrality", kFALSE,
                   100, 0.0, 3000.0, VAR::kSPDntracklets, 100,0.,100., VAR::kCentSPD);
      if(!fRunLightWeight) HIST::Instance()->AddHistogram(classStr.Data(),"MultTPCvsCentSPD","Multiplicity;TPC selected tracks;TPC centrality", kFALSE,
                   100, 0.0, 3500.0, VAR::kNtracksSelected, 100,0.,100., VAR::kCentTPC);
      if(!fRunLightWeight) HIST::Instance()->AddHistogram(classStr.Data(),"MultZDCvsCentZDC","Multiplicity;multiplicity ZDC;ZDC centrality", kFALSE,
                   100, 0.0, 300000.0, VAR::kZDCTotalEnergy, 100,0.,100., VAR::kCentZDC);


      if(!fRunLightWeight) HIST::Instance()->AddHistogram(classStr.Data(),"MultTPCvsMultVZERO","Multiplicity;tracks TPC;multiplicity VZERO", kFALSE,
                   100, 0.0, 3500.0, VAR::kNtracksSelected, 100, 0.0, 32000.0, VAR::kVZEROTotalMult);
      if(!fRunLightWeight) HIST::Instance()->AddHistogram(classStr.Data(),"MultTPCvsMultSPD","Multiplicity;tracklets SPD;tracks TPC", kFALSE,
                   100, 0.0, 3500.0, VAR::kNtracksSelected, 100, 0.0, 3000.0, VAR::kSPDntracklets);
      if(!fRunLightWeight) HIST::Instance()->AddHistogram(classStr.Data(),"MultSPDvsMultVZERO","Multiplicity;tracklets SPD;multiplicity VZERO", kFALSE,
                   100, 0.0, 32000.0, VAR::kVZEROTotalMult, 100, 0.0, 3000.0, VAR::kSPDntracklets);
      if(!fRunLightWeight) HIST::Instance()->AddHistogram(classStr.Data(),"MultTPCvsMultZDC","Multiplicity;tracks TPC;energy ZDC", kFALSE,
                   100, 0.0, 3500.0, VAR::kNtracksSelected, 100, 0.0, 300000.0, VAR::kZDCTotalEnergy);
      if(!fRunLightWeight) HIST::Instance()->AddHistogram(classStr.Data(),"MultVZEROvsMultZDC","Multiplicity;multiplicity VZERO;energy ZDC", kFALSE,
                   100, 0.0, 32000.0, VAR::kVZEROTotalMult, 100, 0.0, 300000.0, VAR::kZDCTotalEnergy);
      if(!fRunLightWeight) HIST::Instance()->AddHistogram(classStr.Data(),"MultSPDvsMultZDC","Multiplicity;tracklets SPD;energy ZDC", kFALSE,
                   100, 0.0, 3000.0, VAR::kSPDntracklets, 100, 0.0, 300000.0, VAR::kZDCTotalEnergy);
      if(!fRunLightWeight) HIST::Instance()->AddHistogram(classStr.Data(),"MultTPCvsMultFMD1","Multiplicity;tracks TPC;multiplicity FMD1", kFALSE,
                   100, 0.0, 3500.0, VAR::kNtracksSelected, 100, 0.0, 10000.0, VAR::kFMD1TotalMult);
      if(!fRunLightWeight) HIST::Instance()->AddHistogram(classStr.Data(),"MultTPCvsMultFMD2I","Multiplicity;tracks TPC;multiplicity FMD2I", kFALSE,
                   100, 0.0, 3500.0, VAR::kNtracksSelected, 100, 0.0, 10000.0, VAR::kFMD2ITotalMult);
      if(!fRunLightWeight) HIST::Instance()->AddHistogram(classStr.Data(),"MultTPCvsMultFMD2O","Multiplicity;tracks TPC;multiplicity FMD2O", kFALSE,
                   100, 0.0, 3500.0, VAR::kNtracksSelected, 100, 0.0, 10000.0, VAR::kFMD2OTotalMult);
      if(!fRunLightWeight) HIST::Instance()->AddHistogram(classStr.Data(),"MultTPCvsMultFMD3I","Multiplicity;tracks TPC;multiplicity FMD3I", kFALSE,
                   100, 0.0, 3500.0, VAR::kNtracksSelected, 100, 0.0, 10000.0, VAR::kFMD3ITotalMult);
      if(!fRunLightWeight) HIST::Instance()->AddHistogram(classStr.Data(),"MultTPCvsMultFMD3O","Multiplicity;tracks TPC;multiplicity FMD3O", kFALSE,
                   100, 0.0, 3500.0, VAR::kNtracksSelected, 100, 0.0, 10000.0, VAR::kFMD3OTotalMult);

      

      if(!fRunLightWeight) HIST::Instance()->AddHistogram(classStr.Data(),"MultVZERO","Multiplicity;multiplicity VZERO", kFALSE,
                   320, 0.0, 25000.0, VAR::kVZEROTotalMult);
      if(!fRunLightWeight) HIST::Instance()->AddHistogram(classStr.Data(),"MultVZEROA","Multiplicity;multiplicity VZEROA", kFALSE,
                   250, 0.0, 10000.0, VAR::kVZEROATotalMult);
      if(!fRunLightWeight) HIST::Instance()->AddHistogram(classStr.Data(),"MultVZEROC","Multiplicity;multiplicity VZEROC", kFALSE,
                   250, 0.0, 15000.0, VAR::kVZEROCTotalMult);
      if(!fRunLightWeight) HIST::Instance()->AddHistogram(classStr.Data(),"MultZDC","Multiplicity;multiplicity ZDC", kFALSE,
                   200, 0.0, 300000.0, VAR::kZDCTotalEnergy);
      if(!fRunLightWeight) HIST::Instance()->AddHistogram(classStr.Data(),"MultZDCA","Multiplicity;multiplicity ZDCA", kFALSE,
                   200, 0.0, 150000.0, VAR::kZDCATotalEnergy);
      if(!fRunLightWeight) HIST::Instance()->AddHistogram(classStr.Data(),"MultZDCC","Multiplicity;multiplicity ZDCC", kFALSE,
                   200, 0.0, 150000.0, VAR::kZDCCTotalEnergy);
      if(!fRunLightWeight) HIST::Instance()->AddHistogram(classStr.Data(),"MultFMD1","Multiplicity;multiplicity FMD1", kFALSE,
                   300, 0.0, 10000.0, VAR::kFMD1TotalMult);
      if(!fRunLightWeight) HIST::Instance()->AddHistogram(classStr.Data(),"MultFMD2I","Multiplicity;multiplicity FMD2I", kFALSE,
                   300, 0.0, 10000.0, VAR::kFMD2ITotalMult);
      if(!fRunLightWeight) HIST::Instance()->AddHistogram(classStr.Data(),"MultFMD2O","Multiplicity;multiplicity FMD2O", kFALSE,
                   300, 0.0, 10000.0, VAR::kFMD2OTotalMult);
      if(!fRunLightWeight) HIST::Instance()->AddHistogram(classStr.Data(),"MultFMD3I","Multiplicity;multiplicity FMD3I", kFALSE,
                   300, 0.0, 10000.0, VAR::kFMD3ITotalMult);
      if(!fRunLightWeight) HIST::Instance()->AddHistogram(classStr.Data(),"MultFMD3O","Multiplicity;multiplicity FMD3O", kFALSE,
                   300, 0.0, 10000.0, VAR::kFMD3OTotalMult);
      if(!fRunLightWeight) HIST::Instance()->AddHistogram(classStr.Data(),"MultTZEROA","Multiplicity;multiplicity TZEROA", kFALSE,
                   300, 0.0, 3000.0, VAR::kTZEROATotalMult);
      if(!fRunLightWeight) HIST::Instance()->AddHistogram(classStr.Data(),"MultTZEROC","Multiplicity;multiplicity TZEROC", kFALSE,
                   300, 0.0, 3000.0, VAR::kTZEROCTotalMult);



            
      if(!fRunLightWeight) HIST::Instance()->AddHistogram(classStr.Data(),"CentVZERO","Centrality(VZERO);centrality VZERO (percents)", kFALSE,
                   100, 0.0, 100.0, VAR::kCentVZERO);
      if(!fRunLightWeight) HIST::Instance()->AddHistogram(classStr.Data(),"CentSPD","Centrality(SPD);centrality SPD (percents)", kFALSE,
                   100, 0.0, 100.0, VAR::kCentSPD);
      if(!fRunLightWeight) HIST::Instance()->AddHistogram(classStr.Data(),"CentTPC","Centrality(TPC);centrality TPC (percents)", kFALSE,
                   100, 0.0, 100.0, VAR::kCentTPC);
      if(!fRunLightWeight) HIST::Instance()->AddHistogram(classStr.Data(),"CentZDC","Centrality(ZDC);centrality ZDC (percents)", kFALSE,
                   100, 0.0, 100.0, VAR::kCentZDC);

      HIST::Instance()->AddHistogram(classStr.Data(),"CentQuality","Centrality quality;centrality quality", kFALSE,
                   100, -50.5, 49.5, VAR::kCentQuality);
      HIST::Instance()->AddHistogram(classStr.Data(),"CentVZERO_Run_prof","<Centrality(VZERO)> vs run;Run; centrality VZERO (%)", kTRUE,
                   kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 100, 0.0, 100.0, VAR::kCentVZERO);
      HIST::Instance()->AddHistogram(classStr.Data(),"CentSPD_Run_prof","<Centrality(SPD)> vs run;Run; centrality SPD (%)", kTRUE,
                   kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 100, 0.0, 100.0, VAR::kCentSPD);
      HIST::Instance()->AddHistogram(classStr.Data(),"CentTPC_Run_prof","<Centrality(TPC)> vs run;Run; centrality TPC (%)", kTRUE,
                   kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 100, 0.0, 100.0, VAR::kCentTPC);
      HIST::Instance()->AddHistogram(classStr.Data(),"CentZDC_Run_prof","<Centrality(ZDC)> vs run;Run; centrality ZDC (%)", kTRUE,
                   kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 100, 0.0, 100.0, VAR::kCentZDC);
      
      
      HIST::Instance()->AddHistogram(classStr.Data(),"NV0sTotal","Number of V0 candidates per event;# pairs", kFALSE,
                   1000,0.,30000.,VAR::kNV0total);
      HIST::Instance()->AddHistogram(classStr.Data(),"NV0sSelected","Number of selected V0 candidates per event;# pairs", kFALSE,
                   1000,0.,10000.,VAR::kNV0selected);
      HIST::Instance()->AddHistogram(classStr.Data(),"NPairs","Number of candidates per event;# pairs", kFALSE,
                   5000,0.,5000.,VAR::kNdielectrons);
      HIST::Instance()->AddHistogram(classStr.Data(),"NPairsSelected", "Number of selected pairs per event; #pairs", kFALSE,
	           5000,0.,5000.,VAR::kNpairsSelected);
      HIST::Instance()->AddHistogram(classStr.Data(),"NTracksTotal","Number of total tracks per event;# tracks", kFALSE,
                   1000,0.,30000.,VAR::kNtracksTotal);
      HIST::Instance()->AddHistogram(classStr.Data(),"NTracksSelected","Number of selected tracks per event;# tracks", kFALSE,
                   1000,0.,30000.,VAR::kNtracksSelected);
      HIST::Instance()->AddHistogram(classStr.Data(),"SPDntracklets", "SPD #tracklets in |#eta|<1.0; tracklets", kFALSE,
                   3000, -0.5, 2999.5, VAR::kSPDntracklets);
      
      HIST::Instance()->AddHistogram(classStr.Data(),"NV0total_Run_prof", "<Number of total V0s> per run; Run; #tracks", kTRUE,
	           kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 100, 0., 10000., VAR::kNV0total);
      HIST::Instance()->AddHistogram(classStr.Data(),"NV0selected_Run_prof", "<Number of selected V0s> per run; Run; #tracks", kTRUE,
	           kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 100, 0., 10000., VAR::kNV0selected);
      HIST::Instance()->AddHistogram(classStr.Data(),"Ndielectrons_Run_prof", "<Number of dielectrons> per run; Run; #tracks", kTRUE,
	           kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 100, 0., 10000., VAR::kNdielectrons);
      HIST::Instance()->AddHistogram(classStr.Data(),"NpairsSelected_Run_prof", "<Number of selected pairs> per run; Run; #tracks", kTRUE,
	           kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 100, 0., 10000., VAR::kNpairsSelected);
      HIST::Instance()->AddHistogram(classStr.Data(),"NTracksTotal_Run_prof", "<Number of tracks> per run; Run; #tracks", kTRUE,
	           kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 100, 0., 10000., VAR::kNtracksTotal);
      HIST::Instance()->AddHistogram(classStr.Data(),"NTracksSelected_Run_prof", "<Number of selected tracks> per run; Run; #tracks", kTRUE,
	           kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 100, 0., 10000., VAR::kNtracksSelected);
      HIST::Instance()->AddHistogram(classStr.Data(),"SPDntracklets_Run_prof", "<SPD ntracklets> per run; Run; #tracks", kTRUE,
	           kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 100, 0., 10000., VAR::kSPDntracklets);
      
      HIST::Instance()->AddHistogram(classStr.Data(),"VtxZ_CentVZERO","Centrality(VZERO) vs vtx. Z;vtx Z (cm); centrality VZERO (%)", kFALSE,
                   300,-15.,15.,VAR::kVtxZ, 100, 0.0, 100.0, VAR::kCentVZERO);
      HIST::Instance()->AddHistogram(classStr.Data(),"VtxZ_CentSPD","Centrality(SPD) vs vtx. Z;vtx Z (cm); centrality SPD (%)", kFALSE,
                   300,-15.,15.,VAR::kVtxZ, 100, 0.0, 100.0, VAR::kCentSPD);
      HIST::Instance()->AddHistogram(classStr.Data(),"VtxZ_CentTPC","Centrality(TPC) vs vtx. Z;vtx Z (cm); centrality TPC (%)", kFALSE,
                   300,-15.,15.,VAR::kVtxZ, 100, 0.0, 100.0, VAR::kCentTPC);
      continue;
    }  // end if className contains "Event"    
      
    
    // Offline trigger histograms
    if(classStr.Contains("OfflineTriggers")) {
      HIST::Instance()->AddHistClass(classStr.Data());

      TString triggerNames = "";
      for(Int_t i=0; i<64; ++i) {triggerNames += Form("%s",VAR::fOfflineTriggerNames[i]); triggerNames+=";";}
      
      HIST::Instance()->AddHistogram(classStr.Data(), "Triggers", "Offline triggers fired; ; ;", kFALSE,
	           64, -0.5, 63.5, VAR::kOfflineTrigger, 2, -0.5, 1.5, VAR::kOfflineTriggerFired, 0, 0.0, 0.0, VAR::kNothing, triggerNames.Data(), "off;on");
      HIST::Instance()->AddHistogram(classStr.Data(), "Triggers2", "Offline triggers fired; ; ;", kFALSE,
	           64, -0.5, 63.5, VAR::kOfflineTriggerFired2, 0, 0.0, 0.0, VAR::kNothing, 0, 0.0, 0.0, VAR::kNothing, triggerNames.Data());
      HIST::Instance()->AddHistogram(classStr.Data(), "CentVZERO_Triggers2", "Offline triggers fired vs centrality VZERO; ; centrality VZERO;", kFALSE,
	           64, -0.5, 63.5, VAR::kOfflineTriggerFired2, 20, 0.0, 100.0, VAR::kCentVZERO, 0, 0.0, 0.0, VAR::kNothing, triggerNames.Data());
      HIST::Instance()->AddHistogram(classStr.Data(), "CentTPC_Triggers2", "Offline triggers fired vs centrality TPC; ; centrality TPC;", kFALSE,
	           64, -0.5, 63.5, VAR::kOfflineTriggerFired2, 20, 0.0, 100.0, VAR::kCentTPC, 0, 0.0, 0.0, VAR::kNothing, triggerNames.Data());
      HIST::Instance()->AddHistogram(classStr.Data(), "CentSPD_Triggers2", "Offline triggers fired vs centrality SPD; ; centrality SPD;", kFALSE,
	           64, -0.5, 63.5, VAR::kOfflineTriggerFired2, 20, 0.0, 100.0, VAR::kCentSPD, 0, 0.0, 0.0, VAR::kNothing, triggerNames.Data());
      HIST::Instance()->AddHistogram(classStr.Data(), "CentZDC_Triggers2", "Offline triggers fired vs centrality ZDC; ; centrality ZDC;", kFALSE,
	           64, -0.5, 63.5, VAR::kOfflineTriggerFired2, 20, 0.0, 100.0, VAR::kCentZDC, 0, 0.0, 0.0, VAR::kNothing, triggerNames.Data());
      HIST::Instance()->AddHistogram(classStr.Data(), "VtxZ_Triggers2", "Offline triggers fired vs vtxZ; ; vtx Z (cm.);", kFALSE,
	           64, -0.5, 63.5, VAR::kOfflineTriggerFired2, 200, -20.0, 20.0, VAR::kVtxZ, 0, 0.0, 0.0, VAR::kNothing, triggerNames.Data());
      continue;
    }

    // Track histograms
    if(classStr.Contains("Tracks")) {
      HIST::Instance()->AddHistClass(classStr.Data());
            for(Int_t ih=0; ih<fgkEPMaxHarmonics; ++ih) {
	HIST::Instance()->AddHistogram(classStr.Data(), Form("Cos%dPhi_CentVZERO",ih+1), Form("<cos%d #varphi> vs (CentVZERO); centrality VZERO; <cos%d #varphi>", ih+1, ih+1), kTRUE,
                     20, 0.0, 100.0, VAR::kCentVZERO, 500, -1., 1., VAR::kCosNPhi+ih);
	HIST::Instance()->AddHistogram(classStr.Data(), Form("Sin%dPhi_CentVZERO",ih+1), Form("<sin%d #varphi> vs (CentVZERO); centrality VZERO; <sin%d #varphi>", ih+1, ih+1), kTRUE,
                     20, 0.0, 100.0, VAR::kCentVZERO, 500, -1.0, 1.0, VAR::kSinNPhi+ih);
            }
    }

     
    // Track histograms
    if(classStr.Contains("TrackQA")) {
      HIST::Instance()->AddHistClass(classStr.Data());
    
      HIST::Instance()->AddHistogram(classStr.Data(), "Pt", "p_{T} distribution; p_{T} (GeV/c^{2});", kFALSE,
                   1000, 0.0, 50.0, VAR::kPt);
      HIST::Instance()->AddHistogram(classStr.Data(), "Eta", "#eta illumination; #eta;", kFALSE,
                   1000, -1.5, 1.5, VAR::kEta);
      HIST::Instance()->AddHistogram(classStr.Data(), "Phi", "#varphi illumination; #varphi;", kFALSE,
                   1000, 0.0, 6.3, VAR::kPhi);
      HIST::Instance()->AddHistogram(classStr.Data(), "DCAxy", "DCAxy; DCAxy (cm.)", kFALSE,
                   1000, -10.0, 10.0, VAR::kDcaXY);
      HIST::Instance()->AddHistogram(classStr.Data(), "DCAz", "DCAz; DCAz (cm.)", kFALSE,
                   1000, -10.0, 10.0, VAR::kDcaZ);
      HIST::Instance()->AddHistogram(classStr.Data(), "TPCncls", "TPCncls; TPCncls", kFALSE,
                   160, 0.0, 160.0, VAR::kTPCncls);

      // run dependence
      HIST::Instance()->AddHistogram(classStr.Data(), "Pt_Run", "<p_{T}> vs run; run;", kTRUE,
                   kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 1000, 0.0, 50.0, VAR::kPt);
      HIST::Instance()->AddHistogram(classStr.Data(), "Eta_Run", "<#eta> vs run; run;", kTRUE,
                   kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 1000, -1.5, 1.5, VAR::kEta);      
      HIST::Instance()->AddHistogram(classStr.Data(), "Phi_Run", "<#varphi> vs run; run;", kTRUE,
                   kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 1000, 0.0, 6.3, VAR::kPhi);      
      HIST::Instance()->AddHistogram(classStr.Data(), "DCAxy_Run", "<DCAxy> vs run; run;", kTRUE,
                   kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 1000, -10.0, 10.0, VAR::kDcaXY);
      HIST::Instance()->AddHistogram(classStr.Data(), "DCAz_Run", "<DCAz> vs run; run;", kTRUE,
                   kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 1000, -10.0, 10.0, VAR::kDcaZ);

      // correlations between parameters
      if(!fRunLightWeight) HIST::Instance()->AddHistogram(classStr.Data(), "Eta_Pt_prof", "<p_{T}> vs #eta; #eta; p_{T} (GeV/c);", kTRUE,
                   300, -1.5, +1.5, VAR::kEta, 100, 0.0, 10.0, VAR::kPt);
      if(!fRunLightWeight) HIST::Instance()->AddHistogram(classStr.Data(), "Phi_Pt", "p_{T} vs #varphi; #varphi (rad.); p_{T} (GeV/c)", kFALSE,
                   300, -0.01, 6.3, VAR::kPhi, 100, 0.0, 2.2, VAR::kPt);
      if(!fRunLightWeight) HIST::Instance()->AddHistogram(classStr.Data(), "Phi_Pt_prof", "<p_{T}> vs #varphi; #varphi (rad.); p_{T} (GeV/c)", kTRUE,
                   300, 0.0, 6.3, VAR::kPhi, 100, 0.0, 10.0, VAR::kPt);
      if(!fRunLightWeight) HIST::Instance()->AddHistogram(classStr.Data(), "Eta_Phi", "#varphi vs #eta; #eta; #varphi (rad.);", kFALSE,
                   200, -1.0, +1.0, VAR::kEta, 100, 0.0, 6.3, VAR::kPhi);
      if(!fRunLightWeight) HIST::Instance()->AddHistogram(classStr.Data(), "TPCncls_Eta_Phi_prof", "<TPC ncls> vs #varphi vs #eta; #eta; #varphi (rad.);TPC ncls", kTRUE,
                   200, -1.0, +1.0, VAR::kEta, 100, 0.0, 6.3, VAR::kPhi, 10, 0.0, 200., VAR::kTPCncls);
      if(!fRunLightWeight) HIST::Instance()->AddHistogram(classStr.Data(), "DCAxy_Eta_Phi_prof", "<DCAxy> vs #varphi vs #eta; #eta; #varphi (rad.);DCAxy (cm)", kTRUE,
                   200, -1.0, +1.0, VAR::kEta, 100, 0.0, 6.3, VAR::kPhi, 10, 0.0, 200., VAR::kDcaXY);
      if(!fRunLightWeight) HIST::Instance()->AddHistogram(classStr.Data(), "DCAz_Eta_Phi_prof", "<DCAz> vs #varphi vs #eta; #eta; #varphi (rad.);DCAz (cm)", kTRUE,
                   200, -1.0, +1.0, VAR::kEta, 100, 0.0, 6.3, VAR::kPhi, 10, 0.0, 200., VAR::kDcaZ);
      if(!fRunLightWeight) HIST::Instance()->AddHistogram(classStr.Data(), "Pt_DCAxy", "DCAxy vs p_{T}; p_{T} (GeV/c); DCA_{xy} (cm)", kFALSE,
                   100, 0.0, 10.0, VAR::kPt, 500, -2.0, 2.0, VAR::kDcaXY);
      if(!fRunLightWeight) HIST::Instance()->AddHistogram(classStr.Data(), "Pt_DCAz", "DCAz vs p_{T}; p_{T} (GeV/c); DCA_{z} (cm)", kFALSE,
                   100, 0.0, 10.0, VAR::kPt, 500, -2.0, 2.0, VAR::kDcaZ);
      if(!fRunLightWeight) HIST::Instance()->AddHistogram(classStr.Data(), "Eta_DCAxy", "DCAxy vs #eta; #eta; DCA_{xy} (cm)", kFALSE,
                   100, -1.0, 1.0, VAR::kEta, 500, -2.0, 2.0, VAR::kDcaXY);
      if(!fRunLightWeight) HIST::Instance()->AddHistogram(classStr.Data(), "Eta_DCAz", "DCAz vs #eta; #eta; DCA_{z} (cm)", kFALSE,
                   100, -1.0, 1.0, VAR::kEta, 500, -2.0, 2.0, VAR::kDcaZ);
      
      for(Int_t ih=0; ih<fgkEPMaxHarmonics; ++ih) {
	//HIST::Instance()->AddHistogram(classStr.Data(), Form("Cos%dPhi_CentVZERO",ih+1), Form("<cos%d #varphi> vs (CentVZERO); centrality VZERO; <cos%d #varphi>", ih+1, ih+1), kTRUE,
  //                   20, 0.0, 100.0, VAR::kCentVZERO, 500, -1., 1., VAR::kCosNPhi+ih);
	//HIST::Instance()->AddHistogram(classStr.Data(), Form("Sin%dPhi_CentVZERO",ih+1), Form("<sin%d #varphi> vs (CentVZERO); centrality VZERO; <sin%d #varphi>", ih+1, ih+1), kTRUE,
  //                   20, 0.0, 100.0, VAR::kCentVZERO, 500, -1.0, 1.0, VAR::kSinNPhi+ih);
	if(!fRunLightWeight) HIST::Instance()->AddHistogram(classStr.Data(), Form("Cos%dPhi_Pt_Eta",ih+1), Form("<cos%d #varphi> vs (#eta,p_{T}); #eta; p_{T} (GeV/c); <cos%d #varphi>", ih+1, ih+1), kTRUE,
                     20, -1.0, 1.0, VAR::kEta, 30, 0.0, 3.0, VAR::kPt, 500, -1.0, 1.0, VAR::kCosNPhi+ih);
	if(!fRunLightWeight) HIST::Instance()->AddHistogram(classStr.Data(), Form("Sin%dPhi_Pt_Eta",ih+1), Form("<sin%d #varphi> vs (#eta,p_{T}); #eta; p_{T} (GeV/c); <sin%d #varphi>", ih+1, ih+1), kTRUE,
                     20, -1.0, 1.0, VAR::kEta, 30, 0.0, 3.0, VAR::kPt, 500, -1.0, 1.0, VAR::kSinNPhi+ih);
	if(!fRunLightWeight) HIST::Instance()->AddHistogram(classStr.Data(), Form("Cos%dPhi_CentVZERO_VtxZ",ih+1), Form("<cos%d #varphi> vs (CentVZERO,VtxZ); Z (cm.); centrality VZERO; <cos%d #varphi>", ih+1, ih+1), kTRUE,
                     30, -15.0, 15.0, VAR::kVtxZ, 20, 0.0, 100.0, VAR::kCentVZERO, 500, -1., 1., VAR::kCosNPhi+ih);
	if(!fRunLightWeight) HIST::Instance()->AddHistogram(classStr.Data(), Form("Sin%dPhi_CentVZERO_VtxZ",ih+1), Form("<sin%d #varphi> vs (CentVZERO,VtxZ); Z (cm.); centrality VZERO; <sin%d #varphi>", ih+1, ih+1), kTRUE,
                     30, -15.0, 15.0, VAR::kVtxZ, 20, 0.0, 100.0, VAR::kCentVZERO, 500, -1.0, 1.0, VAR::kSinNPhi+ih);
      }
    }
    

}

  cout << " done" << endl;
}







