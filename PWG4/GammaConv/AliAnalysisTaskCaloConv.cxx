/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Dmitri Peressounko (RRC KI)                                    *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

////////////////////////////////////////////////
//--------------------------------------------- 
// Class used to do analysis on conversion + calorimeter pairs
// A lot of code cut-and-pasted from AliV0Reader and GammaConversion classes
// 
//---------------------------------------------
////////////////////////////////////////////////

// root
#include "TChain.h"
#include "TH3.h"
#include "TH2.h"
#include "TDirectory.h"
#include "TLorentzVector.h"
#include "TRandom.h"

// analysis
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisTaskCaloConv.h"
#include "AliStack.h"
#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliESDpid.h"
#include "AliESDtrackCuts.h"
#include "AliESDtrackCuts.h"
#include "AliCFContainer.h"   // for CF
#include "AliESDCaloCluster.h" 
#include "AliPHOSGeoUtils.h" 
#include "AliEMCALGeoUtils.h" 
#include "AliCFContainer.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliESDv0.h"
#include "AliKFParticle.h"
#include "AliKFVertex.h"

class Riostream;
class TFile;

ClassImp(AliAnalysisTaskCaloConv)


AliAnalysisTaskCaloConv::AliAnalysisTaskCaloConv():
AliAnalysisTaskSE(),
  fESDEvent(NULL),	
  fESDpid(NULL),
  fESDtrackCuts(NULL),
  fStack(NULL),
  fOutputContainer(NULL),
  fCFOutputContainer(NULL),
  fConvCFCont(0x0),
  fPHOSCFCont(0x0),
  fEMCALCFCont(0x0),
  fPi0CFCont(0x0),
  fCentr(0.),
  fTriggerCINT1B(kFALSE),
  fToUseCF(kFALSE),
  fMinOpeningAngleGhostCut(0.),
  fPHOSgeom(0x0),
  fEMCALgeom(0x0),
  fPi0Thresh1(0.5),
  fPi0Thresh2(1.),
  fBadDistCutPHOS(3.3),
  fBadDistCutEMCAL(6.),
  fGammaV0s(),
  fGammaPHOS(),
  fGammaEMCAL(),
  fConvEvent(NULL) ,
  fPHOSEvent(NULL),
  fEMCALEvent(NULL),
  fnSigmaAboveElectronLine(5.),
  fnSigmaBelowElectronLine(-3.),
  fnSigmaAbovePionLine(0.),
  fpnSigmaAbovePionLine(1.),
  fprobCut(0.),
  fmaxR(180.),
  fmaxZ(240.),
  fetaCut(0.9),
  fptCut(0.02),
  fchi2CutConversion(30.)
{
  // Default constructor
  Int_t nBin=10 ;
  for(Int_t i=0;i<nBin;i++){
    fPHOSEvents[i]=0 ;
    fEMCALEvents[i]=0;
    fConvEvents[i]=0;
  }
  char key[55] ;
  for(Int_t i=0; i<6; i++){
    snprintf(key,55,"PHOS_BadMap_mod%d",i) ;
    fPHOSBadMap[i]=new TH2I(key,"Bad Modules map",64,0.,64.,56,0.,56.) ;
  }
  for(Int_t i=0; i<10; i++){
    snprintf(key,55,"EMCAL_BadMap_mod%d",i) ;
    fEMCALBadMap[i] = new TH2I(key,"Bad Modules map",24,0.,24.,48,0.,48.) ;
  }
}
AliAnalysisTaskCaloConv::AliAnalysisTaskCaloConv(const char* name):
  AliAnalysisTaskSE(name),
  fESDEvent(NULL),
  fESDpid(NULL),
  fESDtrackCuts(NULL),
  fStack(NULL),
  fOutputContainer(NULL),
  fCFOutputContainer(NULL),
  fConvCFCont(0x0),
  fPHOSCFCont(0x0),
  fEMCALCFCont(0x0),
  fPi0CFCont(0x0),
  fCentr(0.),
  fTriggerCINT1B(kFALSE),
  fToUseCF(kFALSE),
  fMinOpeningAngleGhostCut(0.),
  fPHOSgeom(0x0),
  fEMCALgeom(0x0),
  fPi0Thresh1(0.5),
  fPi0Thresh2(1.),
  fBadDistCutPHOS(3.3),
  fBadDistCutEMCAL(6.),
  fGammaV0s(),
  fGammaPHOS(),
  fGammaEMCAL(),
  fConvEvent(NULL) ,
  fPHOSEvent(NULL),
  fEMCALEvent(NULL),
  fnSigmaAboveElectronLine(5.),
  fnSigmaBelowElectronLine(-3.),
  fnSigmaAbovePionLine(0.),
  fpnSigmaAbovePionLine(1.),
  fprobCut(0.),
  fmaxR(180.),
  fmaxZ(240.),
  fetaCut(0.9),
  fptCut(0.02),
  fchi2CutConversion(30.)
{
  // Common I/O in slot 0
  DefineInput (0, TChain::Class());
  DefineOutput(0, TTree::Class());
	
  // Your private output
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());  // for CF
	
  Int_t nBin=10 ;
  for(Int_t i=0;i<nBin;i++){
    fPHOSEvents[i]=0 ;
    fEMCALEvents[i]=0;
    fConvEvents[i]=0;
  }
  char key[55] ;
  for(Int_t i=0; i<6; i++){
    snprintf(key,55,"PHOS_BadMap_mod%d",i) ;
    fPHOSBadMap[i]=new TH2I(key,"Bad Modules map",64,0.,64.,56,0.,56.) ;
  }
  for(Int_t i=0; i<10; i++){
    snprintf(key,55,"EMCAL_BadMap_mod%d",i) ;
    fEMCALBadMap[i] = new TH2I(key,"Bad Modules map",24,0.,24.,48,0.,48.) ;
  }
//  fESDpid = new AliESDpid;
}
//_____________________________________________________
AliAnalysisTaskCaloConv::~AliAnalysisTaskCaloConv() 
{
  // Remove all pointers
	
  if (AliAnalysisManager::GetAnalysisManager()->GetAnalysisType() !=
      AliAnalysisManager::kProofAnalysis) {

    if(fOutputContainer){
      fOutputContainer->Clear() ; 
      delete fOutputContainer ;
    }
    if(fCFOutputContainer){
      fCFOutputContainer->Clear() ; 
      delete fCFOutputContainer ;
    }

    if(fPHOSgeom){
      delete fPHOSgeom ;
      fPHOSgeom=0x0 ;
    }

    if(fEMCALgeom){
      delete fEMCALgeom ;
      fEMCALgeom=0x0;
    }

    for(Int_t ivtx=0; ivtx<10; ivtx++){
      if(fPHOSEvents[ivtx]){
        delete fPHOSEvents[ivtx] ;
          fPHOSEvents[ivtx]=0x0 ;
    }
      if(fEMCALEvents[ivtx]){
        delete fEMCALEvents[ivtx] ;
        fEMCALEvents[ivtx]=0x0 ;
      }
      if(fConvEvents[ivtx]){
        delete fConvEvents[ivtx] ;
        fConvEvents[ivtx]=0x0 ;
      }
    }
    for(Int_t i=0; i<6; i++)
      if(fPHOSBadMap[i]){
        delete fPHOSBadMap[i] ;
          fPHOSBadMap[i]=0 ;
      }
    for(Int_t i=0; i<10; i++)
      if(fEMCALBadMap[i]){
       delete fEMCALBadMap[i];
       fEMCALBadMap[i]=0 ;
    }
  } 
}
//_____________________________________________________
void AliAnalysisTaskCaloConv::Init()
{
  // Initialization
  // AliLog::SetGlobalLogLevel(AliLog::kError);
}
//_____________________________________________________
void AliAnalysisTaskCaloConv::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event
  // First select conversion and calorimeter photons 
  // then construct inv. mass distributions
  //First try to find Stack information.
  if(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()){
    if(static_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler())->MCEvent())
      fStack = static_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler())->MCEvent()->Stack();
  }


  AliESDInputHandler *esdHandler=dynamic_cast<AliESDInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if(!fESDpid){
    if( esdHandler && esdHandler->GetESDpid()){
      fESDpid=new AliESDpid(*(esdHandler->GetESDpid())) ;
    } 
    else {
      fESDpid=new AliESDpid;
      Double_t alephParameters[5];
      if(fStack){// simulation
        alephParameters[0] = 2.15898e+00/50.;
        alephParameters[1] = 1.75295e+01;
        alephParameters[2] = 3.40030e-09;
        alephParameters[3] = 1.96178e+00;
        alephParameters[4] = 3.91720e+00;
        fESDpid->GetTOFResponse().SetTimeResolution(80.);
      }
      else{// data
        alephParameters[0] = 0.0283086;
        alephParameters[1] = 2.63394e+01;
        alephParameters[2] = 5.04114e-11;
        alephParameters[3] = 2.12543e+00;
        alephParameters[4] = 4.88663e+00;
        fESDpid->GetTOFResponse().SetTimeResolution(130.);
        fESDpid->GetTPCResponse().SetMip(47.9);
      }

      fESDpid->GetTPCResponse().SetBetheBlochParameters(
        alephParameters[0],alephParameters[1],alephParameters[2],
        alephParameters[3],alephParameters[4]);
      fESDpid->GetTPCResponse().SetSigma(3.79301e-03, 2.21280e+04);
    }
  }

  if(!fESDtrackCuts){
//    fESDtrackCuts= AliESDtrackCuts::GetStandardITSTPCTrackCuts2010();
    fESDtrackCuts = new AliESDtrackCuts;

    // TPC
    fESDtrackCuts->SetMinNClustersTPC(70);
    fESDtrackCuts->SetMaxChi2PerClusterTPC(4);
    fESDtrackCuts->SetAcceptKinkDaughters(kFALSE);
    fESDtrackCuts->SetRequireTPCRefit(kTRUE);
    // ITS
    fESDtrackCuts->SetRequireITSRefit(kTRUE);
    fESDtrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                                           AliESDtrackCuts::kAny);
//    if(selPrimaries) {
      // 7*(0.0026+0.0050/pt^1.01)
      fESDtrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
//    }
    fESDtrackCuts->SetMaxDCAToVertexZ(2);
    fESDtrackCuts->SetDCAToVertex2D(kFALSE);
    fESDtrackCuts->SetRequireSigmaToVertex(kFALSE);



/*
    fESDtrackCuts = new AliESDtrackCuts("AliESDtrackCuts","AliESDtrackCuts");

    fESDtrackCuts->SetAcceptKinkDaughters(kFALSE); 
    fESDtrackCuts->SetMinNClustersTPC(70); 
    fESDtrackCuts->SetMaxChi2PerClusterTPC(4); 
    fESDtrackCuts->SetRequireTPCRefit(kTRUE); 
    fESDtrackCuts->SetRequireITSRefit(kTRUE); 
//    fESDtrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny); //TEMPORARY <-> REMOVE 
*/
  } 


  fESDEvent=(AliESDEvent*)InputEvent();
  FillHistogram("hEventsTrig",0.5) ;
  Bool_t isSelected =(esdHandler->IsEventSelected()& AliVEvent::kMB) == AliVEvent::kMB;
  if(!isSelected){
    printf("Not selected !!!!! \n") ;
    PostData(1, fOutputContainer);
    return ;
  }
  FillHistogram("hEventsTrig",1.5) ;

  //Take Only events with proper trigger
  //No trigger in MC data => no check
//  if(!fStack && !fESDEvent->IsTriggerClassFired("CINT1B-ABCE-NOPF-ALL")){
  //for LHC10e
  if(!fStack && !fESDEvent->IsTriggerClassFired("CINT1-B-NOPF-ALLNOTRD")){
printf("Trigger failed \n") ;
    PostData(1, fOutputContainer);
    return ;
  } 
  FillHistogram("hEventsTrig",2.5) ;

  //checks if we have a prim vertex
  if(fESDEvent->GetPrimaryVertex()->GetNContributors()<=0) {
    PostData(1, fOutputContainer);
    return ;
  }
  FillHistogram("hEventsTrig",3.5) ;

  if(TMath::Abs(fESDEvent->GetPrimaryVertex()->GetZ())>10.){
    PostData(1, fOutputContainer);
    return ;
  }
  FillHistogram("hEventsTrig",4.5) ;


  //Calculate charged multiplicity
  Int_t trackCounter = 0; 
  for (Int_t i=0;i<fESDEvent->GetNumberOfTracks();++i) { 
    AliESDtrack *track = new AliESDtrack(*fESDEvent->GetTrack(i)) ;
    if(fESDtrackCuts->AcceptTrack(track) &&  TMath::Abs(track->Eta())< 0.9) 
      trackCounter++; 
    delete track; 
  } 
  fCentr=trackCounter+0.5 ;
  FillHistogram("hMult",fCentr) ;
	
  //Init geometry if not done yet
  InitGeometry();

  //Select conversion and calorimeter photons
  //Conversion photons should go first since they are used in calibration in SelectCALOPhotons()
  SelectConvPhotons() ;
  SelectPHOSPhotons() ;
  SelectEMCALPhotons() ;
  //Fill MC histograms if MC is present
  ProcessMC();
  FillRealMixed() ;

  PostData(1, fOutputContainer);
  if(fToUseCF)
    PostData(2, fCFOutputContainer);  // for CF


}
//____________________________________________________________
void AliAnalysisTaskCaloConv::ConnectInputData(Option_t *option){
  // see header file for documentation

  AliAnalysisTaskSE::ConnectInputData(option);

}
//____________________________________________________________
void AliAnalysisTaskCaloConv::UserCreateOutputObjects()
{
  //UserCreateOutputObjects
  if(fDebug)gDirectory->Print() ;
  // Create the output container
  if(fOutputContainer != NULL){
    delete fOutputContainer;
    fOutputContainer = NULL;
  }
  fOutputContainer = new TList();
  fOutputContainer->SetOwner(kTRUE);

  if(fCFOutputContainer != NULL){
    delete fCFOutputContainer;
    fCFOutputContainer = NULL;
  }
  //===========Correction Framework ======================
  if(fToUseCF){
    fCFOutputContainer = new TList();
    fCFOutputContainer->SetOwner(kTRUE);

    //bins: pt,eta,mass
    Int_t iBin[3]={500,40,100};
    fConvCFCont = new AliCFContainer("ConvContainer","container for converted photons", 23,3,iBin);
    fConvCFCont->SetBinLimits(0,0.,50.);
    fConvCFCont->SetBinLimits(1,-2.,2.) ;
    fConvCFCont->SetBinLimits(2,0.,1.);
    fCFOutputContainer->Add(fConvCFCont) ;
  
    fPHOSCFCont = new AliCFContainer("PHOSContainer","container for PHOS photons", 10,2,iBin);
    fPHOSCFCont->SetBinLimits(0,0.,50.);
    fPHOSCFCont->SetBinLimits(1,-2.,2.) ;
    fCFOutputContainer->Add(fPHOSCFCont) ;
  
    fEMCALCFCont = new AliCFContainer("EMCALContainer","container for EMCAL photons", 10,2,iBin);
    fEMCALCFCont->SetBinLimits(0,0.,50.);
    fEMCALCFCont->SetBinLimits(1,-2.,2.) ;
    fCFOutputContainer->Add(fEMCALCFCont) ;
  
    fPi0CFCont = new AliCFContainer("Pi0Container","container for EMCAL photons", 10,2,iBin);
    fPi0CFCont->SetBinLimits(0,0.,50.);
    fPi0CFCont->SetBinLimits(1,-2.,2.) ;
    fCFOutputContainer->Add(fPi0CFCont) ;
  
  }
  //========================================================

  //Adding the histograms to the output container
  Int_t firstRun= 125000 ;
  Int_t lastRun = 135000 ;
  Int_t nRuns =lastRun-firstRun+1 ;

  //Run QA histigrams
  fOutputContainer->Add(new TH2F("hRunTrigger","Triggers fired",nRuns,float(firstRun),float(lastRun),2,0.,2.)) ;
  fOutputContainer->Add(new TH1F("hRunEvents","Events per run",nRuns,float(firstRun),float(lastRun))) ;
  fOutputContainer->Add(new TH1F("hRunConvs","Conversion photons per run",nRuns,float(firstRun),float(lastRun))) ;
  fOutputContainer->Add(new TH1F("hRunPHOS","PHOS photons per run",nRuns,float(firstRun),float(lastRun))) ;
  fOutputContainer->Add(new TH1F("hRunEMCAL","EMCAL photons per run",nRuns,float(firstRun),float(lastRun))) ;
  fOutputContainer->Add(new TH1F("hVtxBin","Vtx distribution",10,0.,10.)) ;
  fOutputContainer->Add(new TH1F("hEvents","Events processed",1,0.,1.)) ;
  fOutputContainer->Add(new TH1F("hEventsTrig","Events processed",10,0.,10.)) ;

  fOutputContainer->Add(new TH2F("hQA_PHOS_mod1_soft","number of clusters per cell",64,0.,64.,56,0.,56.)) ;
  fOutputContainer->Add(new TH2F("hQA_PHOS_mod2_soft","number of clusters per cell",64,0.,64.,56,0.,56.)) ;
  fOutputContainer->Add(new TH2F("hQA_PHOS_mod3_soft","number of clusters per cell",64,0.,64.,56,0.,56.)) ;
  fOutputContainer->Add(new TH2F("hQA_PHOS_mod4_soft","number of clusters per cell",64,0.,64.,56,0.,56.)) ;
  fOutputContainer->Add(new TH2F("hQA_PHOS_mod5_soft","number of clusters per cell",64,0.,64.,56,0.,56.)) ;
  fOutputContainer->Add(new TH2F("hQA_PHOS_mod1_hard","number of clusters per cell",64,0.,64.,56,0.,56.)) ;
  fOutputContainer->Add(new TH2F("hQA_PHOS_mod2_hard","number of clusters per cell",64,0.,64.,56,0.,56.)) ;
  fOutputContainer->Add(new TH2F("hQA_PHOS_mod3_hard","number of clusters per cell",64,0.,64.,56,0.,56.)) ;
  fOutputContainer->Add(new TH2F("hQA_PHOS_mod4_hard","number of clusters per cell",64,0.,64.,56,0.,56.)) ;
  fOutputContainer->Add(new TH2F("hQA_PHOS_mod5_hard","number of clusters per cell",64,0.,64.,56,0.,56.)) ;

  fOutputContainer->Add(new TH2F("hQA_EMCAL_SM0_soft","number of clusters per cell",24,0.,24,48,0.,48.)) ;
  fOutputContainer->Add(new TH2F("hQA_EMCAL_SM1_soft","number of clusters per cell",24,0.,24,48,0.,48.)) ;
  fOutputContainer->Add(new TH2F("hQA_EMCAL_SM2_soft","number of clusters per cell",24,0.,24,48,0.,48.)) ;
  fOutputContainer->Add(new TH2F("hQA_EMCAL_SM3_soft","number of clusters per cell",24,0.,24,48,0.,48.)) ;
  fOutputContainer->Add(new TH2F("hQA_EMCAL_SM4_soft","number of clusters per cell",24,0.,24,48,0.,48.)) ;
  fOutputContainer->Add(new TH2F("hQA_EMCAL_SM5_soft","number of clusters per cell",24,0.,24,48,0.,48.)) ;
  fOutputContainer->Add(new TH2F("hQA_EMCAL_SM6_soft","number of clusters per cell",24,0.,24,48,0.,48.)) ;
  fOutputContainer->Add(new TH2F("hQA_EMCAL_SM7_soft","number of clusters per cell",24,0.,24,48,0.,48.)) ;
  fOutputContainer->Add(new TH2F("hQA_EMCAL_SM8_soft","number of clusters per cell",24,0.,24,48,0.,48.)) ;
  fOutputContainer->Add(new TH2F("hQA_EMCAL_SM9_soft","number of clusters per cell",24,0.,24,48,0.,48.)) ;
  fOutputContainer->Add(new TH2F("hQA_EMCAL_SM0_hard","number of clusters per cell",24,0.,24,48,0.,48.)) ;
  fOutputContainer->Add(new TH2F("hQA_EMCAL_SM1_hard","number of clusters per cell",24,0.,24,48,0.,48.)) ;
  fOutputContainer->Add(new TH2F("hQA_EMCAL_SM2_hard","number of clusters per cell",24,0.,24,48,0.,48.)) ;
  fOutputContainer->Add(new TH2F("hQA_EMCAL_SM3_hard","number of clusters per cell",24,0.,24,48,0.,48.)) ;
  fOutputContainer->Add(new TH2F("hQA_EMCAL_SM4_hard","number of clusters per cell",24,0.,24,48,0.,48.)) ;
  fOutputContainer->Add(new TH2F("hQA_EMCAL_SM5_hard","number of clusters per cell",24,0.,24,48,0.,48.)) ;
  fOutputContainer->Add(new TH2F("hQA_EMCAL_SM6_hard","number of clusters per cell",24,0.,24,48,0.,48.)) ;
  fOutputContainer->Add(new TH2F("hQA_EMCAL_SM7_hard","number of clusters per cell",24,0.,24,48,0.,48.)) ;
  fOutputContainer->Add(new TH2F("hQA_EMCAL_SM8_hard","number of clusters per cell",24,0.,24,48,0.,48.)) ;
  fOutputContainer->Add(new TH2F("hQA_EMCAL_SM9_hard","number of clusters per cell",24,0.,24,48,0.,48.)) ;

  fOutputContainer->Add(new TH2F("hQA_ConvPhiEta","Number of V0s phi eta",100,0.,TMath::TwoPi(),40,-1.5,1.5)) ;

  fOutputContainer->Add(new TH2F("hdEdx","dEdx of acceptaed electrons",1000,0.,10.,150,0.,150.)) ;

  fOutputContainer->Add(new TH1F("hMult","Multiplicity",200,0.,200.)) ;

  Int_t npt=200 ;
  Double_t ptmax=20. ;
  //Calibration of PHOS
  fOutputContainer->Add(new TH3F("PHOS_mod1_th1","Inv.Mass distr. per channel",64,0.,64,56,0.,56,100,0.,0.5)) ;
  fOutputContainer->Add(new TH3F("PHOS_mod2_th1","Inv.Mass distr. per channel",64,0.,64,56,0.,56,100,0.,0.5)) ;
  fOutputContainer->Add(new TH3F("PHOS_mod3_th1","Inv.Mass distr. per channel",64,0.,64,56,0.,56,100,0.,0.5)) ;
  fOutputContainer->Add(new TH3F("PHOS_mod4_th1","Inv.Mass distr. per channel",64,0.,64,56,0.,56,100,0.,0.5)) ;
  fOutputContainer->Add(new TH3F("PHOS_mod5_th1","Inv.Mass distr. per channel",64,0.,64,56,0.,56,100,0.,0.5)) ;

  fOutputContainer->Add(new TH3F("PHOS_mod1_th2","Inv.Mass distr. per channel",64,0.,64,56,0.,56,100,0.,0.5)) ;
  fOutputContainer->Add(new TH3F("PHOS_mod2_th2","Inv.Mass distr. per channel",64,0.,64,56,0.,56,100,0.,0.5)) ;
  fOutputContainer->Add(new TH3F("PHOS_mod3_th2","Inv.Mass distr. per channel",64,0.,64,56,0.,56,100,0.,0.5)) ;
  fOutputContainer->Add(new TH3F("PHOS_mod4_th2","Inv.Mass distr. per channel",64,0.,64,56,0.,56,100,0.,0.5)) ;
  fOutputContainer->Add(new TH3F("PHOS_mod5_th2","Inv.Mass distr. per channel",64,0.,64,56,0.,56,100,0.,0.5)) ;

  //Pi0 histograms
  //Vary Conversion cuts
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_OnFly","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_Offline","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH3F("PHOS_Re_mvsPt_OnFly_mult","Mass vs pt",400,0.,1.,npt,0.,ptmax,150,0.,150.)) ;
  fOutputContainer->Add(new TH3F("PHOS_Re_mvsPt_Offline_mult","Mass vs pt",400,0.,1.,npt,0.,ptmax,150,0.,150.)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_On_dEdx","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_Off_dEdx","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
//  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_On_Prob","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
//  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_Off_Prob","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_On_R120","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_Off_R120","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_On_Z","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_Off_Z","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_On_chi","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_Off_chi","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_On_Eta","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_Off_Eta","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_On_Wcut","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_On_Wcut_Neu","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_Off_Wcut","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_Off_Wcut_Neu","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_On_ArmQt","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_Off_ArmQt","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;

  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_OnFly","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_Offline","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_On_dEdx","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_Off_dEdx","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
//  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_On_Prob","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
//  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_Off_Prob","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_On_R120","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_Off_R120","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_On_Z","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_Off_Z","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_On_chi","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_Off_chi","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_On_Eta","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_Off_Eta","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_On_Wcut","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_Off_Wcut","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_On_ArmQt","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_Off_ArmQt","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;

  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_OnFly","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_Offline","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH3F("EMCAL_Re_mvsPt_OnFly_mult","Mass vs pt",400,0.,1.,npt,0.,ptmax,30,0.,60.)) ;
  fOutputContainer->Add(new TH3F("EMCAL_Re_mvsPt_Offline_mult","Mass vs pt",400,0.,1.,npt,0.,ptmax,30,0.,60.)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_On_ArmQt","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_Off_ArmQt","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_On_dEdx","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_Off_dEdx","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
//  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_On_Prob","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
//  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_Off_Prob","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_On_R120","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_Off_R120","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_On_Z","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_Off_Z","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_On_chi","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_Off_chi","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_On_Eta","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_Off_Eta","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_On_Wcut","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_Off_Wcut","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;

  fOutputContainer->Add(new TH2F("EMCAL_Mi_mvsPt_OnFly","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Mi_mvsPt_Offline","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Mi_mvsPt_On_ArmQt","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Mi_mvsPt_Off_ArmQt","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Mi_mvsPt_On_dEdx","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Mi_mvsPt_Off_dEdx","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
//  fOutputContainer->Add(new TH2F("EMCAL_Mi_mvsPt_On_Prob","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
//  fOutputContainer->Add(new TH2F("EMCAL_Mi_mvsPt_Off_Prob","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Mi_mvsPt_On_R120","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Mi_mvsPt_Off_R120","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Mi_mvsPt_On_Z","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Mi_mvsPt_Off_Z","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Mi_mvsPt_On_chi","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Mi_mvsPt_Off_chi","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Mi_mvsPt_On_Eta","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Mi_mvsPt_Off_Eta","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Mi_mvsPt_On_Wcut","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Mi_mvsPt_Off_Wcut","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;

  //PHOS PID variations
  fOutputContainer->Add(new TH3F("PHOS_Re_mvsPt_alpha","Mass vs pt vs PHOS E",400,0.,1.,npt,0.,ptmax,10,0.,1.)) ;
  fOutputContainer->Add(new TH3F("PHOS_Re_mvsPt_E","Mass vs pt vs PHOS E",400,0.,1.,npt,0.,ptmax,100,0.,10.)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_all","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_all_dist","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_Disp","Mass vs pt, disp cut",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_TOF","Mass vs pt, TOF cut",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_Neutral","Mass vs pt, Neutral cut",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_DispNeutral","Mass vs pt, Disp & neutral cut",400,0.,1.,npt,0.,ptmax)) ;

  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_all","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_all_dist","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_Disp","Mass vs pt, disp cut",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_TOF","Mass vs pt, TOF cut",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_Neutral","Mass vs pt, Neutral cut",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_DispNeutral","Mass vs pt, Disp & neutral cut",400,0.,1.,npt,0.,ptmax)) ;

  char key[155] ;
  for(Int_t mod=1; mod<=5;mod++){
    snprintf(key,155,"PHOS_Re_mvsPt_mod%d_single",mod) ;
    fOutputContainer->Add(new TH2F(key,"Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
    snprintf(key,155,"PHOS_Re_mvsPt_mod%d_all",mod) ;
    fOutputContainer->Add(new TH2F(key,"Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  }

  //Single photon spectrum
  //Conversion
  fOutputContainer->Add(new TH1F("Single_conv_OnFly","Single photon spectrum",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("Single_conv_Offline","Single photon spectrum",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("Single_conv_On_ArmQt","Single photon spectrum",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("Single_conv_Off_ArmQt","Single photon spectrum",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("Single_conv_On_dEdx","Single photon spectrum",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("Single_conv_Off_dEdx","Single photon spectrum",npt,0.,ptmax)) ;
//  fOutputContainer->Add(new TH1F("Single_conv_On_Prob","Single photon spectrum",npt,0.,ptmax)) ;
//  fOutputContainer->Add(new TH1F("Single_conv_Off_Prob","Single photon spectrum",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("Single_conv_On_R120","Single photon spectrum",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("Single_conv_Off_R120","Single photon spectrum",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("Single_conv_On_Z","Single photon spectrum",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("Single_conv_Off_Z","Single photon spectrum",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("Single_conv_On_chi","Single photon spectrum",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("Single_conv_Off_chi","Single photon spectrum",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("Single_conv_On_Eta","Single photon spectrum",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("Single_conv_Off_Eta","Single photon spectrum",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("Single_conv_On_Wcut","Single photon spectrum",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("Single_conv_Off_Wcut","Single photon spectrum",npt,0.,ptmax)) ;

  //PHOS
  fOutputContainer->Add(new TH2F("PHOS_single_all_mult","Single photon spectrum",npt,0.,ptmax,150,0.,150.)) ;
  fOutputContainer->Add(new TH2F("PHOS_single_disp_mult","Single photon spectrum",npt,0.,ptmax,150,0.,150.)) ;
  fOutputContainer->Add(new TH2F("PHOS_single_neu_mult","Single photon spectrum",npt,0.,ptmax,150,0.,150.)) ;
   
  for(Int_t mod=1; mod<=5;mod++){
    snprintf(key,155,"PHOS_single_mod%d_all",mod) ;
    fOutputContainer->Add(new TH1F(key,"Single photon spectrum",npt,0.,ptmax)) ;
    snprintf(key,155,"PHOS_single_mod%d_disp",mod) ;
    fOutputContainer->Add(new TH1F(key,"Single photon spectrum",npt,0.,ptmax)) ;
    snprintf(key,155,"PHOS_single_mod%d_neutral",mod) ;
    fOutputContainer->Add(new TH1F(key,"Single photon spectrum",npt,0.,ptmax)) ;
    snprintf(key,155,"PHOS_single_mod%d_dispneutral",mod) ;
    fOutputContainer->Add(new TH1F(key,"Single photon spectrum",npt,0.,ptmax)) ;
    snprintf(key,155,"PHOS_single_mod%d_dist1",mod) ;
    fOutputContainer->Add(new TH1F(key,"Single photon spectrum",npt,0.,ptmax)) ;
    snprintf(key,155,"PHOS_single_mod%d_dist2",mod) ;
    fOutputContainer->Add(new TH1F(key,"Single photon spectrum",npt,0.,ptmax)) ;
  }

  for(Int_t mod=0; mod<4;mod++){
    snprintf(key,155,"EMCAL_mod%d_th1",mod) ;
    fOutputContainer->Add(new TH3F(key,"Inv.Mass distr. per channel",24,0.,24,48,0.,48,100,0.,0.5)) ;
    snprintf(key,155,"EMCAL_mod%d_th2",mod) ;
    fOutputContainer->Add(new TH3F(key,"Inv.Mass distr. per channel",24,0.,24,48,0.,48,100,0.,0.5)) ;

    snprintf(key,155,"EMCAL_Re_mvsPt_mod%d_single",mod) ;
    fOutputContainer->Add(new TH2F(key,"Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;

    snprintf(key,155,"EMCAL_Re_mvsPt_mod%d_all",mod) ;
    fOutputContainer->Add(new TH2F(key,"Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
    snprintf(key,155,"EMCAL_Re_mvsPt_mod%d_Disp",mod) ;
    fOutputContainer->Add(new TH2F(key,"Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
    snprintf(key,155,"EMCAL_Re_mvsPt_mod%d_TOF",mod) ;
    fOutputContainer->Add(new TH2F(key,"Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
    snprintf(key,155,"EMCAL_Re_mvsPt_mod%d_Neutral",mod) ;
    fOutputContainer->Add(new TH2F(key,"Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
    snprintf(key,155,"EMCAL_Re_mvsPt_mod%d_DispNeutral",mod) ;
    fOutputContainer->Add(new TH2F(key,"Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;

    snprintf(key,155,"EMCAL_Mi_mvsPt_mod%d_all",mod) ;
    fOutputContainer->Add(new TH2F(key,"Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
    snprintf(key,155,"EMCAL_Mi_mvsPt_mod%d_Disp",mod) ;
    fOutputContainer->Add(new TH2F(key,"Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
    snprintf(key,155,"EMCAL_Mi_mvsPt_mod%d_TOF",mod) ;
    fOutputContainer->Add(new TH2F(key,"Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
    snprintf(key,155,"EMCAL_Mi_mvsPt_mod%d_Neutral",mod) ;
    fOutputContainer->Add(new TH2F(key,"Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
    snprintf(key,155,"EMCAL_Mi_mvsPt_mod%d_DispNeutral",mod) ;
    fOutputContainer->Add(new TH2F(key,"Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  }

  //MC info
  fOutputContainer->Add(new TH2F("hMC_CaloConv_pi0_unitEta","Primary #pi^{0}",npt,0.,ptmax,150,0.,150.)) ;
  fOutputContainer->Add(new TH2F("hMC_CaloConv_eta_unitEta","Primary #pi^{0}",npt,0.,ptmax,150,0.,150.)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_allpi0","Primary #pi^{0}",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_alleta","Primary #pi^{0}",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_pi0PHOSacc","#pi^{0} decayed in PHOS acc",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_etaPHOSacc","#pi^{0} decayed in PHOS acc",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_pi0EMCALacc","#pi^{0} decayed in EMCAL acc",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_etaEMCALacc","#pi^{0} decayed in EMCAL acc",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_pi0_PHOS_conv","#pi^{0} decayed in PHOS acc asnd conv. photon",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_eta_PHOS_conv","#pi^{0} decayed in PHOS acc asnd conv. photon",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_pi0_EMCAL_conv","#pi^{0} decayed in EMCAL acc asnd conv. photon",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_eta_EMCAL_conv","#pi^{0} decayed in EMCAL acc asnd conv. photon",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_pi0_bothphot_conv","#pi^{0} both photons converted",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_eta_bothphot_conv","#pi^{0} both photons converted",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_pi0_convPhotInCalo","#pi^{0} photon in calo converted",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_eta_convPhotInCalo","#pi^{0} photon in calo converted",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_pi0_v0onfly_PHOSacc","#pi^{0} photon converted and V0 found",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_eta_v0onfly_PHOSacc","#pi^{0} photon converted and V0 found",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_pi0_v0offline_PHOSacc","#pi^{0} photon converted and V0 found",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_eta_v0offline_PHOSacc","#pi^{0} photon converted and V0 found",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_pi0_v0onfly_EMCALacc","#pi^{0} photon converted and V0 found",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_eta_v0onfly_EMCALacc","#pi^{0} photon converted and V0 found",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_pi0_v0offline_EMCALacc","#pi^{0} photon converted and V0 found",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_eta_v0offline_EMCALacc","#pi^{0} photon converted and V0 found",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_pi0_v0onfly_PHOSclu","#pi^{0} V0 and cluster in PHOS found",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_eta_v0onfly_PHOSclu","#pi^{0} V0 and cluster in PHOS found",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_pi0_v0onfly_PHOSclu_pid","#pi^{0} V0 and cluster in PHOS found",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_eta_v0onfly_PHOSclu_pid","#pi^{0} V0 and cluster in PHOS found",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_pi0_v0offline_PHOSclu","#pi^{0} V0 and cluster in PHOS found",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_eta_v0offline_PHOSclu","#pi^{0} V0 and cluster in PHOS found",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_pi0_v0offline_PHOSclu_pid","#pi^{0} V0 and cluster in PHOS found",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_eta_v0offline_PHOSclu_pid","#pi^{0} V0 and cluster in PHOS found",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_pi0_v0onfly_PHOSclu_good","#pi^{0} V0 and cluster in PHOS found",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_eta_v0onfly_PHOSclu_good","#pi^{0} V0 and cluster in PHOS found",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_pi0_v0offline_PHOSclu_good","#pi^{0} V0 and cluster in PHOS found",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_eta_v0offline_PHOSclu_good","#pi^{0} V0 and cluster in PHOS found",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_pi0_v0onfly_PHOSclu_mod1","#pi^{0} V0 and cluster in PHOS found",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_eta_v0onfly_PHOSclu_mod1","#pi^{0} V0 and cluster in PHOS found",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_pi0_v0onfly_PHOSclu_mod2","#pi^{0} V0 and cluster in PHOS found",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_eta_v0onfly_PHOSclu_mod2","#pi^{0} V0 and cluster in PHOS found",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_pi0_v0onfly_PHOSclu_mod3","#pi^{0} V0 and cluster in PHOS found",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_eta_v0onfly_PHOSclu_mod3","#pi^{0} V0 and cluster in PHOS found",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_pi0_v0onfly_PHOSclu_mod4","#pi^{0} V0 and cluster in PHOS found",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_eta_v0onfly_PHOSclu_mod4","#pi^{0} V0 and cluster in PHOS found",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_pi0_v0onfly_PHOSclu_mod5","#pi^{0} V0 and cluster in PHOS found",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_eta_v0onfly_PHOSclu_mod5","#pi^{0} V0 and cluster in PHOS found",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_pi0_v0offline_PHOSclu_mod1","#pi^{0} V0 and cluster in PHOS found",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_eta_v0offline_PHOSclu_mod1","#pi^{0} V0 and cluster in PHOS found",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_pi0_v0offline_PHOSclu_mod2","#pi^{0} V0 and cluster in PHOS found",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_eta_v0offline_PHOSclu_mod2","#pi^{0} V0 and cluster in PHOS found",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_pi0_v0offline_PHOSclu_mod3","#pi^{0} V0 and cluster in PHOS found",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_eta_v0offline_PHOSclu_mod3","#pi^{0} V0 and cluster in PHOS found",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_pi0_v0offline_PHOSclu_mod4","#pi^{0} V0 and cluster in PHOS found",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_eta_v0offline_PHOSclu_mod4","#pi^{0} V0 and cluster in PHOS found",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_pi0_v0offline_PHOSclu_mod5","#pi^{0} V0 and cluster in PHOS found",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_eta_v0offline_PHOSclu_mod5","#pi^{0} V0 and cluster in PHOS found",npt,0.,ptmax)) ;

  fOutputContainer->Add(new TH1F("hMC_CaloConv_pi0_v0on_PHOSclu_ptRec","#pi^{0} V0 and cluster in PHOS found(rec pt)",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_eta_v0on_PHOSclu_ptRec","#pi^{0} V0 and cluster in PHOS found(rec pt)",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_pi0_v0off_PHOSclu_ptRec","#pi^{0} V0 and cluster in PHOS found(rec pt)",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_eta_v0off_PHOSclu_ptRec","#pi^{0} V0 and cluster in PHOS found(rec pt)",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("hMC_CaloConv_pi0_v0on_PHOSclu_mvsPt","m vs pt for rec pi0s",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("hMC_CaloConv_eta_v0on_PHOSclu_mvsPt","m vs pt for rec pi0s",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("hMC_CaloConv_pi0_v0off_PHOSclu_mvsPt","m vs pt for rec pi0s",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("hMC_CaloConv_eta_v0off_PHOSclu_mvsPt","m vs pt for rec pi0s",400,0.,1.,npt,0.,ptmax)) ;

  fOutputContainer->Add(new TH1F("hMC_CaloConv_pi0_v0off_EMCALclu_ptRec","#pi^{0} V0 and cluster in EMCAL found (rec pt)",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_pi0_v0on_EMCALclu_ptRec","#pi^{0} V0 and cluster in EMCAL found (rec pt)",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_pi0_v0onfly_EMCALclu","#pi^{0} V0 and cluster in EMCAL found",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_pi0_v0onfly_EMCALclu_pid","#pi^{0} V0 and cluster in EMCAL found",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_pi0_v0offline_EMCALclu","#pi^{0} V0 and cluster in EMCAL found",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_pi0_v0offline_EMCALclu_pid","#pi^{0} V0 and cluster in EMCAL found",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_pi0_v0onfly_EMCALclu_good","#pi^{0} V0 and cluster in EMCAL found",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_pi0_v0offline_EMCALclu_good","#pi^{0} V0 and cluster in EMCAL found",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("hMC_CaloConv_pi0_v0on_EMCALclu_mvsPt","m vs pt for rec pi0s",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("hMC_CaloConv_pi0_v0off_EMCALclu_mvsPt","m vs pt for rec pi0s",400,0.,1.,npt,0.,ptmax)) ;

  fOutputContainer->Add(new TH3F("hMC_Resid_PHOS_Phot_mvsPt","m vs pt for rec pi0s",400,0.,1.,npt,0.,ptmax,10,0.,1.)) ;
  fOutputContainer->Add(new TH3F("hMC_Resid_PHOS_Pi0_mvsPt","m vs pt for rec pi0s",400,0.,1.,npt,0.,ptmax,10,0.,1.)) ;
  fOutputContainer->Add(new TH3F("hMC_Resid_PHOS_eta_mvsPt","m vs pt for rec pi0s",400,0.,1.,npt,0.,ptmax,10,0.,1.)) ;
  fOutputContainer->Add(new TH3F("hMC_Resid_PHOS_K_mvsPt","m vs pt for rec pi0s",400,0.,1.,npt,0.,ptmax,10,0.,1.)) ;
  fOutputContainer->Add(new TH3F("hMC_Resid_PHOS_pi_mvsPt","m vs pt for rec pi0s",400,0.,1.,npt,0.,ptmax,10,0.,1.)) ;
  fOutputContainer->Add(new TH3F("hMC_Resid_PHOS_pbar_mvsPt","m vs pt for rec pi0s",400,0.,1.,npt,0.,ptmax,10,0.,1.)) ;
  fOutputContainer->Add(new TH3F("hMC_Resid_PHOS_other_mvsPt","m vs pt for rec pi0s",400,0.,1.,npt,0.,ptmax,10,0.,1.)) ;
  fOutputContainer->Add(new TH3F("hMC_Resid_EMCAL_Phot_mvsPt","m vs pt for rec pi0s",400,0.,1.,npt,0.,ptmax,10,0.,1.)) ;
  fOutputContainer->Add(new TH3F("hMC_Resid_EMCAL_Pi0_mvsPt","m vs pt for rec pi0s",400,0.,1.,npt,0.,ptmax,10,0.,1.)) ;
  fOutputContainer->Add(new TH3F("hMC_Resid_EMCAL_eta_mvsPt","m vs pt for rec pi0s",400,0.,1.,npt,0.,ptmax,10,0.,1.)) ;
  fOutputContainer->Add(new TH3F("hMC_Resid_EMCAL_K_mvsPt","m vs pt for rec pi0s",400,0.,1.,npt,0.,ptmax,10,0.,1.)) ;
  fOutputContainer->Add(new TH3F("hMC_Resid_EMCAL_pi_mvsPt","m vs pt for rec pi0s",400,0.,1.,npt,0.,ptmax,10,0.,1.)) ;
  fOutputContainer->Add(new TH3F("hMC_Resid_EMCAL_pbar_mvsPt","m vs pt for rec pi0s",400,0.,1.,npt,0.,ptmax,10,0.,1.)) ;
  fOutputContainer->Add(new TH3F("hMC_Resid_EMCAL_other_mvsPt","m vs pt for rec pi0s",400,0.,1.,npt,0.,ptmax,10,0.,1.)) ;


  fOutputContainer->Add(new TH1F("hMC_CaloConv_phot","Primary photons",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_gammaPHOSacc","Photons in PHOS acc",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_gammaEMCALacc","Photons in EMCAL acc",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_gamma_conv","Converted photons",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_gamma_v0","Converted photons with V0",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("hMC_CaloConv_gamma_v0_devsE","Converted photons with V0",200,-1.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_gamma_PHOSclu","Photons with cluster in PHOS",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_gamma_PHOSclu_dist1","Photons with cluster in PHOS",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_gamma_PHOSclu_dist2","Photons with cluster in PHOS",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_gamma_EMCALclu","Photons with cluster in EMCAL",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_gamma_EMCALclu_dist1","Photons with cluster in EMCAL",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_gamma_EMCALclu_dist2","Photons with cluster in EMCAL",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_gamma_PHOSclu_recE","Photons with cluster in PHOS",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_gamma_EMCALclu_recE","Photons with cluster in EMCAL",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("hMC_CaloConv_gamma_PHOSclu_devsE","Photons with cluster in PHOS",200,-1.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("hMC_CaloConv_gamma_EMCALclu_devsE","Photons with cluster in EMCAL",200,-1.,1.,npt,0.,ptmax)) ;
	

  //Non-linearity test
  char keym[55] ;
  for(Int_t iw=0;iw<10;iw++){ //resolution
    for(Int_t in=0;in<10;in++){
      snprintf(keym,55,"hMC_nonlinearity_w%d_n%d",iw,in) ;
      fOutputContainer->Add(new TH2F(keym,"m vs pt, nonlinearity test" ,200,0.,0.5,npt,0.,ptmax)) ;
      snprintf(keym,55,"hMC_nonlinearity_ConvPHOS_w%d_n%d",iw,in) ;
      fOutputContainer->Add(new TH2F(keym,"m vs pt, nonlinearity test" ,200,0.,0.5,npt,0.,ptmax)) ;
      snprintf(keym,55,"hMC_nonlinearity_EMCAL_w%d_n%d",iw,in) ;
      fOutputContainer->Add(new TH2F(keym,"m vs pt, nonlinearity test" ,200,0.,0.5,npt,0.,ptmax)) ;
      snprintf(keym,55,"hMC_CaloConv_pi0_v0onfly_PHOSclu_w%d_n%d",iw,in) ;
      fOutputContainer->Add(new TH1F(keym,"#pi^{0} V0 and cluster in PHOS found",npt,0.,ptmax)) ;
      snprintf(keym,55,"hMC_CaloConv_pi0_v0onfly_ConvPHOSclu_w%d_n%d",iw,in) ;
      fOutputContainer->Add(new TH1F(keym,"#pi^{0} V0 and cluster in PHOS found",npt,0.,ptmax)) ;
      snprintf(keym,55,"hMC_CaloConv_pi0_v0onfly_EMCALclu_w%d_n%d",iw,in) ;
      fOutputContainer->Add(new TH1F(keym,"#pi^{0} V0 and cluster in PHOS found",npt,0.,ptmax)) ;
    }
  }


  fOutputContainer->Add(new TH3F("All_chi2_eta_pt","MC chi2 vs eta vs phi",100,0.,100.,200,-2.,2.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("All_w_vs_m","MC w vs m",300,0.,TMath::Pi(),400,0.,1.)) ;
  fOutputContainer->Add(new TH3F("MC_V0_pt_eta_phi","MC pt vs eta vs phi",npt,0.,ptmax,200,-2.,2.,200,0.,TMath::TwoPi())) ;
  fOutputContainer->Add(new TH3F("MC_V0_m_eta_pt","MC m vs eta vs phi",400,0.,1.,200,-2.,2.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH3F("MC_V0_chi2_eta_pt","MC chi2 vs eta vs phi",100,0.,100.,200,-2.,2.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("MC_V0_w_vs_m","MC w vs m",300,0.,TMath::Pi(),400,0.,1.)) ;
 
  fOutputContainer->SetName(GetName());
}
//______________________________________________________________________
 void AliAnalysisTaskCaloConv::InitGeometry()
{
  //If not done yet, create Geometry for PHOS and EMCAL
  //and read misalignment matrixes from ESD/AOD (AOD not implemented yet)
  //

   if(fPHOSgeom && fEMCALgeom){ //already initialized
     return ;
   }

   AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent()) ;
   if(!esd )
     AliFatal("Can not read geometry matrixes from ESD/AOD: NO ESD") ;
   if(!fPHOSgeom){//reading PHOS matrixes
     fPHOSgeom = new AliPHOSGeoUtils("IHEP","");
     for(Int_t mod=0; mod<5; mod++){
       if(esd){
         const TGeoHMatrix* m=esd->GetPHOSMatrix(mod) ;
         if(m)
           fPHOSgeom->SetMisalMatrix(m, mod) ;
       }
     }
   }
   if(!fEMCALgeom){
     fEMCALgeom = new AliEMCALGeoUtils("EMCAL_FIRSTYEAR");
     for(Int_t mod=0; mod < (fEMCALgeom->GetEMCGeometry())->GetNumberOfSuperModules(); mod++){ //<---Gustavo, could you check???
       if(esd){
         const TGeoHMatrix* m=esd->GetEMCALMatrix(mod) ;
         if(m)
           fEMCALgeom->SetMisalMatrix(m, mod) ;
       }
     }
   }
}
//________________________________________________________________
void AliAnalysisTaskCaloConv::SelectPHOSPhotons(){
  //SelectPHOSPhotons
  // Loop over all CaloClusters 
  if(fPHOSEvent)
    fPHOSEvent->Clear() ;
  else
    fPHOSEvent = new TClonesArray("TLorentzVector",10) ;
  Int_t inPHOS = 0;
  TLorentzVector pi0 ;

  //vertex
  Double_t vtx[3];
  vtx[0] = fESDEvent->GetPrimaryVertex()->GetX();
  vtx[1] = fESDEvent->GetPrimaryVertex()->GetY();
  vtx[2] = fESDEvent->GetPrimaryVertex()->GetZ();
  for (Int_t i=0; i<fESDEvent->GetNumberOfCaloClusters(); i++) {
    AliESDCaloCluster * clu = fESDEvent->GetCaloCluster(i);
    if(!clu->IsPHOS())
      continue ;
    TLorentzVector p ;
    clu ->GetMomentum(p ,vtx);
    if(p.Energy()<0.25){
      continue ;
    }
    if(clu->GetNCells()<=2){
      continue ;
    }

    Bool_t isNeutral = kTRUE ;
    Bool_t isDispOK = kTRUE ;
    Bool_t isTOFOK = kTRUE ;
    Int_t iMod,iX,iZ ;

    isNeutral = clu->GetEmcCpvDistance()>5. ;  //To be improved
    isDispOK = kFALSE ;
    Double_t l0=clu->GetM02(),l1=clu->GetM20() ;
    if(l1>= 0   && l0>= 0   && l1 < 0.1 && l0 < 0.1) isDispOK=kFALSE ;
    if(l1>= 0   && l0 > 0.5 && l1 < 0.1 && l0 < 1.5) isDispOK=kTRUE ;
    if(l1>= 0   && l0 > 2.0 && l1 < 0.1 && l0 < 2.7) isDispOK=kFALSE ;
    if(l1>= 0   && l0 > 2.7 && l1 < 0.1 && l0 < 4.0) isDispOK=kFALSE ;
    if(l1 > 0.1 && l1 < 0.7 && l0 > 0.7 && l0 < 2.1) isDispOK=kTRUE ;
    if(l1 > 0.1 && l1 < 0.3 && l0 > 3.0 && l0 < 5.0) isDispOK=kFALSE  ;
    if(l1 > 0.3 && l1 < 0.7 && l0 > 2.5 && l0 < 4.0) isDispOK=kFALSE ;
    if(l1 > 0.7 && l1 < 1.3 && l0 > 1.0 && l0 < 1.6) isDispOK=kTRUE ;
    if(l1 > 0.7 && l1 < 1.3 && l0 > 1.6 && l0 < 3.5) isDispOK=kTRUE ;
    if(l1 > 1.3 && l1 < 3.5 && l0 > 1.3 && l0 < 3.5) isDispOK=kTRUE ;

    Float_t xyz[3] = {0,0,0};
    clu->GetPosition(xyz);   //Global position in ALICE system
    TVector3 global(xyz) ;
    Int_t relid[4] ;
    if(!fPHOSgeom->GlobalPos2RelId(global,relid)){
      printf("PHOS_beyond: x=%f, y=%f, z=%f \n",xyz[0],xyz[1],xyz[2]) ;
      continue ;
    }
    iMod=relid[0] ;
    iX=relid[2];
    iZ=relid[3] ;
    if(!IsGoodChannel("PHOS",iMod,iX,iZ))
      continue ;

    Bool_t closeToBad=(clu->GetDistanceToBadChannel()>fBadDistCutPHOS) ;

    p.SetBit(kCaloPIDdisp,isDispOK) ;
    p.SetBit(kCaloPIDtof,isTOFOK) ;
    p.SetBit(kCaloPIDneutral,isNeutral) ;
    p.SetBit(BIT(17+iMod),kTRUE) ;
    p.SetBit(kCaloDistBad,closeToBad) ;
    new((*fPHOSEvent)[inPHOS]) TLorentzVector(p) ;
    fGammaPHOS[inPHOS] = i ;
    inPHOS++ ;


    //Single photon spectra
    Double_t pt= p.Pt() ;
    TString skey="PHOS_single_"; skey+="mod" ; skey+=iMod ; skey+="_all" ;
    FillHistogram(skey,pt) ;
    FillHistogram("PHOS_single_all_mult",pt,fCentr) ;
    if(isDispOK){
      skey="PHOS_single_"; skey+="mod" ; skey+=iMod ; skey+="_disp" ;
      FillHistogram("PHOS_single_disp_mult",pt,fCentr) ;
      FillHistogram(skey,pt) ;
    }
    if(isNeutral){
      skey="PHOS_single_"; skey+="mod" ; skey+=iMod ; skey+="_neutral" ;
      FillHistogram(skey,pt) ;
      FillHistogram("PHOS_single_neu_mult",pt,fCentr) ;
    }
    if(isNeutral && isDispOK){
      skey="PHOS_single_"; skey+="mod" ; skey+=iMod ; skey+="_dispneutral" ;
      FillHistogram(skey,pt) ;
    }
    //Distance to bad channel
    if(clu->GetDistanceToBadChannel()>fBadDistCutPHOS){
      skey="PHOS_single_"; skey+="mod" ; skey+=iMod ; skey+="_dist1" ;
      FillHistogram(skey,pt) ;
    }
    if(clu->GetDistanceToBadChannel()>2.*fBadDistCutPHOS){
      skey="PHOS_single_"; skey+="mod" ; skey+=iMod ; skey+="_dist2" ;
      FillHistogram(skey,pt) ;
    }
 
    //Fill QA
    if(clu->E()>0.5){
      skey="hQA_PHOS_mod"; skey+=iMod; skey+="_soft" ; 
      FillHistogram(skey,iX-0.5, iZ-0.5,1.) ;
      if(clu->E()>1.5){
        skey="hQA_PHOS_mod"; skey+=iMod; skey+="_hard" ;
        FillHistogram(skey,iX-0.5, iZ-0.5,1.) ;
      }
    }

    //Fill histogams for calibration
    if(clu->E()<fPi0Thresh1 ) continue;
    for(Int_t iconv=0;iconv<fConvEvent->GetEntriesFast();iconv++){
      TLorentzVector *gammaConv=static_cast<TLorentzVector*>(fConvEvent->At(iconv)) ;
      pi0=*gammaConv+p;
      skey="PHOS_"; skey+="mod" ; skey+=iMod ; skey+="_th1" ;
      FillHistogram(skey,iX-0.5, iZ-0.5,pi0.M());
      if(isNeutral){
        skey="PHOS_"; skey+="mod"; skey+=iMod ; skey+="_th2" ;
        FillHistogram(skey,iX-0.5, iZ-0.5,pi0.M());
      }
    }
  }
}
//____________________________________________________________
void AliAnalysisTaskCaloConv::SelectEMCALPhotons(){
  //SelectEMCALPhotons
  // Loop over all CaloClusters
  if(fEMCALEvent)
    fEMCALEvent->Clear() ;
  else
    fEMCALEvent = new TClonesArray("TLorentzVector",10) ;
  Int_t inEMCAL = 0 ; //, inEMCALRecal=0;
  TLorentzVector pi0 ;
return ; //Untill EMCAL geometry will be fixed

  //vertex
  Double_t vtx[3];
  vtx[0] = fESDEvent->GetPrimaryVertex()->GetX();
  vtx[1] = fESDEvent->GetPrimaryVertex()->GetY();
  vtx[2] = fESDEvent->GetPrimaryVertex()->GetZ();
  Int_t nEMCAL=0 ; //For QA
  for (Int_t i=0; i<fESDEvent->GetNumberOfCaloClusters(); i++) {
    AliESDCaloCluster * clu = fESDEvent->GetCaloCluster(i);
    if(clu->IsPHOS())
      continue ;
    TLorentzVector p ;
    TLorentzVector pRecal ;
    clu ->GetMomentum(p ,vtx);
    Bool_t isNeutral = kTRUE ;
    Bool_t isDispOK = kTRUE ;
    Bool_t isTOFOK = kTRUE ;
    Int_t iMod,iX,iZ ;

    if(clu->E()>0.1 ){ 
      nEMCAL++ ;
      isNeutral = clu->GetEmcCpvDistance()>10. ;  //To be improved
      if(clu->GetTOF()>550.e-9 && clu->GetTOF()<750.e-9)
        isTOFOK=kTRUE ;
      else
        isTOFOK=kFALSE ;
      Float_t phi = p.Phi();
      if(phi < 0) phi+=TMath::TwoPi();
      Int_t absId = -1;
      fEMCALgeom->GetAbsCellIdFromEtaPhi(p.Eta(),phi, absId);
      iMod=fEMCALgeom->GetSuperModuleNumber(absId) ;
      Int_t imod = -1, iTower = -1, iIphi = -1, iIeta = -1, iphi = -1, ieta = -1;
      fEMCALgeom->GetCellIndex(absId,imod,iTower,iIphi,iIeta);
      fEMCALgeom->GetCellPhiEtaIndexInSModule(imod,iTower, iIphi, iIeta,iphi,ieta);
      if(imod<0 || imod>5){
        printf("EMCAL: Beyond the geometry!\n") ;
        printf("phi=%f, eta=%f, absId=%d, SM=%d \n",p.Eta(),phi, absId, imod) ;
        continue ;
      }
   
      iX=iphi+1 ;
      iZ=ieta+1 ;
      if(!IsGoodChannel("EMCAL",iMod,iX,iZ))
        continue ;
      p.SetBit(kCaloPIDdisp,isDispOK) ;
      p.SetBit(kCaloPIDtof,isTOFOK) ;
      p.SetBit(kCaloPIDneutral,isNeutral) ;
      p.SetBit(BIT(17+imod),kTRUE) ;
      new((*fEMCALEvent)[inEMCAL]) TLorentzVector(p) ;
      fGammaEMCAL[inEMCAL] = i ;
      inEMCAL++ ;
 

      //Fill QA histograms
      if(clu->E()>0.5 && iMod>=0){ //Sometimes modules is negative not found??
        TString skey="hQA_EMCAL_SM";skey+=iMod ; skey+="_soft" ;
        FillHistogram(skey,iX-0.5, iZ-0.5,1.) ;
        if(clu->E()>1.5){
          skey="hQA_EMCAL_SM";skey+=iMod ; skey+="_hard" ;
          FillHistogram(skey,iX-0.5, iZ-0.5,1.) ;
        }
      }

      //Fill histograms for recalibration
      if(clu->E()<fPi0Thresh1) continue ;
      for(Int_t iconv=0;iconv<fConvEvent->GetEntriesFast();iconv++){
        TLorentzVector *gammaConv=static_cast<TLorentzVector*>(fConvEvent->At(iconv)) ;
        pi0=*gammaConv+p;
        TString skey="EMCAL_"; skey+="mod" ; skey+=iMod ; skey+="_th1" ;
        FillHistogram(skey,iX-0.5, iZ-0.5,pi0.M());
        if(clu->E()>fPi0Thresh2){
          skey="EMCAL_"; skey+="mod"; skey+=iMod ; skey+="_th2" ;
          FillHistogram(skey,iX-0.5, iZ-0.5,pi0.M());
        }
      }
    }
  }
}
//______________________________________________________________________
void AliAnalysisTaskCaloConv::SelectConvPhotons(){
  //Fill list of conversion photons
  //that is scan v0s and select photon-like

  //set some constants
  const Double_t cutSigmaMass=0.0001;  //Constraint on photon mass
  const Bool_t useImprovedVertex=kTRUE ; //Use verted with converted photon?
//  const Double_t zrSlope = TMath::Tan(2*TMath::ATan(TMath::Exp(-fetaCut)));
  const Double_t zrSlope12 = TMath::Tan(2*TMath::ATan(TMath::Exp(-1.2)));
  const Double_t zrSlope09 = TMath::Tan(2*TMath::ATan(TMath::Exp(-0.9)));
  const Double_t zOffset = 7.;

  if(!fConvEvent)
    fConvEvent = new TClonesArray("TLorentzVector",10) ;
  else
    fConvEvent->Clear() ;

  //No primary vertex in event => scip
  if(fESDEvent->GetPrimaryVertex()->GetNContributors()<=0) {
    return;
  }

  Int_t inConv=0 ;
  for(Int_t iv0=0; iv0<fESDEvent->GetNumberOfV0s();iv0++){
    AliESDv0 * v0 = fESDEvent->GetV0(iv0) ;

    AliESDtrack * pos = fESDEvent->GetTrack(v0->GetPindex()) ;
    AliESDtrack * neg = fESDEvent->GetTrack(v0->GetNindex()) ;
    const AliExternalTrackParam * paramPos = v0->GetParamP() ;
    const AliExternalTrackParam * paramNeg = v0->GetParamN() ;
    if(pos->GetSign() <0){//change tracks
      pos=neg ;
      neg=fESDEvent->GetTrack(v0->GetPindex()) ;
      paramPos=paramNeg ;
      paramNeg=v0->GetParamP() ;
    }
    AliKFParticle negKF(*paramNeg,11);
    AliKFParticle posKF(*paramPos,-11);
    AliKFParticle photKF(negKF,posKF) ;
    photKF.SetMassConstraint(0,cutSigmaMass);

    if(useImprovedVertex){
      AliKFVertex primaryVertexImproved(*(fESDEvent->GetPrimaryVertex()));
      //if Vtx do created
      if(primaryVertexImproved.GetNContributors()>1){
        primaryVertexImproved+=photKF;
        photKF.SetProductionVertex(primaryVertexImproved);
      }
    }
    Double_t m=0., width=0. ;
    photKF.GetMass(m,width);

    TLorentzVector photLV;
//    photLV.SetXYZM(negKF.Px()+posKF.Px(),negKF.Py()+posKF.Py(),negKF.Pz()+negKF.Pz(),0.) ;
//    photLV.SetXYZT(photKF.GetPx(),photKF.GetPy(),photKF.GetPz(),photKF.GetE()) ;
    photLV.SetXYZM(photKF.GetPx(),photKF.GetPy(),photKF.GetPz(),0.) ;  //Produces slightly better pi0 width

    //Parameters for correction function
    Double_t a[3]={photLV.Pt(),photLV.Eta(),m} ;
    if(fToUseCF)
      fConvCFCont->Fill(a,0) ;

    //select V0 finder
    Bool_t isOnFly=kTRUE ;
    //select V0Finder
    if (v0->GetOnFlyStatus()){
      if(fToUseCF)
        fConvCFCont->Fill(a,1) ;
    }
    else{
      isOnFly=kFALSE ;
      if(fToUseCF)
        fConvCFCont->Fill(a,2) ;
    }

    //Number of TPC clusters
    if(neg->GetNcls(1) <2 || pos->GetNcls(1) <2){
      continue ; 
    }

    //remove like sign pairs 
    if(pos->GetSign() == neg->GetSign()){ 
      continue ;
    }
    if(fToUseCF){
      if(isOnFly)
        fConvCFCont->Fill(a,3) ;
      else
        fConvCFCont->Fill(a,4) ;
    }

    if( !(pos->GetStatus() & AliESDtrack::kTPCrefit) ||
        !(neg->GetStatus() & AliESDtrack::kTPCrefit) ){
      continue;
    }
    if(fToUseCF){
      if(isOnFly)
        fConvCFCont->Fill(a,5) ;
      else
        fConvCFCont->Fill(a,6) ;
    }
 
    if( neg->GetKinkIndex(0) > 0 ||
        pos->GetKinkIndex(0) > 0) {
      continue ;
    }

    //First rough PID
    if( fESDpid->NumberOfSigmasTPC(pos,AliPID::kElectron)<fnSigmaBelowElectronLine ||
        fESDpid->NumberOfSigmasTPC(pos,AliPID::kElectron)>fnSigmaAboveElectronLine ||
        fESDpid->NumberOfSigmasTPC(neg,AliPID::kElectron)<fnSigmaBelowElectronLine ||
        fESDpid->NumberOfSigmasTPC(neg,AliPID::kElectron)>fnSigmaAboveElectronLine ){
        continue ;
    }
    const Double_t minPnSigmaAbovePionLine = 1. ;
    const Double_t maxPnSigmaAbovePionLine = 3. ;
    const Double_t nSigmaAbovePionLine = 0 ;
    if(pos->P()>minPnSigmaAbovePionLine && pos->P()<maxPnSigmaAbovePionLine ){
      if(fESDpid->NumberOfSigmasTPC(pos,AliPID::kPion)<nSigmaAbovePionLine){
          continue ;
        }
    }
    if(neg->P()>minPnSigmaAbovePionLine && neg->P()<maxPnSigmaAbovePionLine){
      if(fESDpid->NumberOfSigmasTPC(neg,AliPID::kPion)<nSigmaAbovePionLine){
          continue ;
      }
    }
    //Strict dEdx
    Bool_t isdEdx=kTRUE;
    if(pos->P()>minPnSigmaAbovePionLine && pos->P()<maxPnSigmaAbovePionLine ){
      if(fESDpid->NumberOfSigmasTPC(pos,AliPID::kPion)<2.){
        isdEdx=kFALSE;
      }
    }
    if(neg->P()>minPnSigmaAbovePionLine && neg->P()<maxPnSigmaAbovePionLine){
      if(fESDpid->NumberOfSigmasTPC(neg,AliPID::kPion)<2.){
        isdEdx=kFALSE;
      }
    }
 

    //Kaon rejection
    const Double_t minPKaonRejection=1.5 ;
    const Double_t sigmaAroundLine=1. ;
    if(neg->P()<minPKaonRejection ){
      if(TMath::Abs(fESDpid->NumberOfSigmasTPC(neg,AliPID::kKaon))<sigmaAroundLine){
        isdEdx=kFALSE;
      }
    }
    if(pos->P()<minPKaonRejection ){
      if(TMath::Abs(fESDpid->NumberOfSigmasTPC(pos,AliPID::kKaon))<sigmaAroundLine){
        isdEdx=kFALSE;
      }
    }

    //Proton rejection
    const Double_t minPProtonRejection=2. ;
    if(neg->P()<minPProtonRejection){
      if(TMath::Abs(fESDpid->NumberOfSigmasTPC(neg,AliPID::kProton))<sigmaAroundLine){
        isdEdx=kFALSE;
      }
    }
    if(pos->P()<minPProtonRejection ){
      if(TMath::Abs(fESDpid->NumberOfSigmasTPC(pos,AliPID::kProton))<sigmaAroundLine){
        isdEdx=kFALSE;
      }
    }

    const Double_t minPPionRejection=0.5 ;
    if(neg->P()<minPPionRejection ){
      if(TMath::Abs(fESDpid->NumberOfSigmasTPC(neg,AliPID::kPion))<sigmaAroundLine){
        isdEdx=kFALSE;
      }
    }
    if(pos->P()<minPPionRejection ){
      if( TMath::Abs(fESDpid->NumberOfSigmasTPC(pos,AliPID::kPion))<sigmaAroundLine){
        isdEdx=kFALSE;
      }
    }


    if(isdEdx){
      FillHistogram("hdEdx",paramPos->GetP(),pos->GetTPCsignal()) ;
      FillHistogram("hdEdx",paramNeg->GetP(),neg->GetTPCsignal()) ;
    }


    //Check the pid probability
    Bool_t isProb=kTRUE ;
    Double_t posProbArray[10];
    Double_t negProbArray[10];
    neg->GetTPCpid(negProbArray);
    pos->GetTPCpid(posProbArray);
    if(negProbArray[AliPID::kElectron]<fprobCut || posProbArray[AliPID::kElectron]<fprobCut){
      isProb=kFALSE ;
    }
    if(fToUseCF){
      if(isOnFly)
        fConvCFCont->Fill(a,9) ;
      else
        fConvCFCont->Fill(a,10) ;
    }

    Double_t v0x=0.,v0y=0.,v0z=0.;
    v0->GetXYZ(v0x,v0y,v0z) ;
    Double_t r=TMath::Sqrt(v0x*v0x + v0y*v0y) ;
    //Remove Dalitz
    const Double_t rMin=2.8 ;
    if(r<rMin)
      continue ;
    if(r>fmaxR){ // cuts on distance from collision point
      continue;
    }
    Bool_t isStrictR=kFALSE ;
    if(r<120.)
      isStrictR=kTRUE ;
    if(fToUseCF){
      if(isOnFly)
        fConvCFCont->Fill(a,11) ;
      else
        fConvCFCont->Fill(a,12) ;
    }


    if((TMath::Abs(v0z)*zrSlope12)-zOffset > r ){ // cuts out regions where we do not reconstruct
      continue;
    }
    if(fToUseCF){
      if(isOnFly)
        fConvCFCont->Fill(a,13) ;
      else
        fConvCFCont->Fill(a,14) ;
    }

    if(TMath::Abs(v0z) > fmaxZ ){ // cuts out regions where we do not reconstruct
      continue;
    }
    Bool_t isStrictZ=kFALSE ;
    if((TMath::Abs(v0z)*zrSlope09)-zOffset < r )
      isStrictZ=kTRUE ;

    if(fToUseCF){
      if(isOnFly)
        fConvCFCont->Fill(a,15) ;
      else
        fConvCFCont->Fill(a,16) ;
    }
 
    if(photKF.GetNDF()<=0){
      continue;
    }
    if(fToUseCF){
      if(isOnFly)
        fConvCFCont->Fill(a,17) ;
      else
        fConvCFCont->Fill(a,18) ;
    }

    Double_t chi2V0 = photKF.GetChi2()/photKF.GetNDF();
    FillHistogram("All_chi2_eta_pt",chi2V0,photLV.Eta(),photLV.Pt()) ;

    if(chi2V0 > fchi2CutConversion || chi2V0 <=0){
      continue;
    }
    Bool_t isStrictChi=kFALSE ;
    if(chi2V0 < 0.7*fchi2CutConversion && chi2V0 >0){
      isStrictChi=kTRUE;
    }
 
    if(fToUseCF){
      if(isOnFly)
        fConvCFCont->Fill(a,19) ;
      else
        fConvCFCont->Fill(a,20) ;
    }

    const Double_t wideEtaCut=1.2 ;
    if(TMath::Abs(photLV.Eta())> wideEtaCut){
      continue;
    }
    if(TMath::Abs(paramPos->Eta())> wideEtaCut ||  
       TMath::Abs(paramNeg->Eta())> wideEtaCut ){
      continue ;
    }

    Bool_t isWideEta=kTRUE ;
    if(TMath::Abs(photLV.Eta())< fetaCut && 
       TMath::Abs(paramPos->Eta())<fetaCut  && 
       TMath::Abs(paramNeg->Eta()) < fetaCut){
      isWideEta=kFALSE;
    }
    

    if(photLV.Pt()<fptCut){
      continue;
    }
    if(fToUseCF){
      if(isOnFly)
        fConvCFCont->Fill(a,21) ;
      else
        fConvCFCont->Fill(a,22) ;
    }


    //Just QA plot
    if(photLV.Pt()>0.5){
       Double_t phi=photLV.Phi() ;
       while(phi<0.)phi+=TMath::TwoPi() ;
       while(phi>TMath::TwoPi())phi-=TMath::TwoPi() ;
       FillHistogram("hQA_ConvPhiEta",phi,photLV.Eta()) ; 
    }
    
    Double_t w=PlanarityAngle(paramPos,paramNeg) ;
    Bool_t isPlanarityCut = (0.08-0.22*w > m || 0.15*(w-2.4)>m) ;
    FillHistogram("All_w_vs_m",w,m) ;

    const Double_t armenterosAlphaCut=0.05 ;
    Double_t armenterosQtAlfa[2]={0.,0.}  ;
    GetArmenterosQtAlfa(&posKF, &negKF, &photKF, armenterosQtAlfa ) ;
    Bool_t isArmQt=(armenterosQtAlfa[1]<armenterosAlphaCut) ;

    photLV.SetBit(kConvOnFly,isOnFly) ;
    photLV.SetBit(kConvArmQt,isArmQt) ;
    photLV.SetBit(kConvdEdx,isdEdx) ;
    photLV.SetBit(kConvProb,isProb) ;
    photLV.SetBit(kConvR,isStrictR) ;
    photLV.SetBit(kConvZR,isStrictZ) ;
    photLV.SetBit(kConvNDF,isStrictChi) ;
    photLV.SetBit(kConvEta,isWideEta) ;
    photLV.SetBit(kConvPlan,isPlanarityCut) ;

    new((*fConvEvent)[inConv]) TLorentzVector(photLV) ;
    fGammaV0s[inConv] = iv0 ;
    inConv++ ;

    //Single photon spectrum
    Double_t pt=photLV.Pt() ;
    if(isOnFly){
      //Default
      if(!isWideEta)
        FillHistogram("Single_conv_OnFly",pt) ;
      if(isdEdx && !isWideEta)
        FillHistogram("Single_conv_On_dEdx",pt) ;
      if(!isWideEta && isStrictR)
        FillHistogram("Single_conv_On_R120",pt) ; 
      if( !isWideEta && isStrictZ)
        FillHistogram("Single_conv_On_Z",pt) ;
      if(!isWideEta && isStrictChi)
        FillHistogram("Single_conv_On_chi",pt) ;
      if(1)
        FillHistogram("Single_conv_On_Eta",pt) ;
      if(!isWideEta && isPlanarityCut)
        FillHistogram("Single_conv_On_Wcut",pt) ;
      if(!isWideEta && isArmQt)
        FillHistogram("Single_conv_On_ArmQt",pt) ;
    }
    else{
      if(!isWideEta)
        FillHistogram("Single_conv_Offline",pt) ;
      if(isdEdx && !isWideEta)
        FillHistogram("Single_conv_Off_dEdx",pt) ;
      if(!isWideEta && isStrictR)
        FillHistogram("Single_conv_Off_R120",pt) ; 
      if(!isWideEta && isStrictZ)
        FillHistogram("Single_conv_Off_Z",pt) ;
      if(!isWideEta && isStrictChi)
        FillHistogram("Single_conv_Off_chi",pt) ;
      if(1)
        FillHistogram("Single_conv_Off_Eta",pt) ;
      if(!isWideEta && isPlanarityCut)
        FillHistogram("Single_conv_Off_Wcut",pt) ;
      if(!isWideEta && isArmQt)
        FillHistogram("Single_conv_Off_ArmQt",pt) ;
    }

    //Fill MC information
    if(fStack){
      TParticle * negativeMC = fStack->Particle(TMath::Abs(neg->GetLabel()));
      TParticle * positiveMC = fStack->Particle(TMath::Abs(pos->GetLabel()));

      if(negativeMC && positiveMC){
        if(negativeMC->GetMother(0) != positiveMC->GetMother(0))
          continue ;
      }

      if(TMath::Abs(negativeMC->GetPdgCode())!=11 || TMath::Abs(positiveMC->GetPdgCode())!=11){
        continue;
      }
      if(negativeMC->GetPdgCode()==positiveMC->GetPdgCode()){
        continue;
      }

      TParticle * v0Gamma = fStack->Particle(negativeMC->GetMother(0));
 
      if(negativeMC->GetUniqueID() != 5 || positiveMC->GetUniqueID() !=5){ // id 5 is conversion
        continue;
      }
      if(v0Gamma->GetPdgCode() == 22){
        FillHistogram("MC_V0_pt_eta_phi",v0Gamma->Pt(),v0Gamma->Eta(),v0Gamma->Phi()) ;  
        FillHistogram("MC_V0_m_eta_pt",m,v0Gamma->Eta(),v0Gamma->Pt()) ;
        FillHistogram("MC_V0_chi2_eta_pt",chi2V0,v0Gamma->Eta(),v0Gamma->Pt()) ;
        FillHistogram("MC_V0_w_vs_m",w,m) ;
      }
    }
  }
}
//______________________________________________________________________
void AliAnalysisTaskCaloConv::FillRealMixed(){
  // Fills real (same event) and Mixed (different events) inv.mass dsitributions
  // Moves current event to the list of stored events at the end

  Double_t vtx[3];
  vtx[0] = fESDEvent->GetPrimaryVertex()->GetX();
  vtx[1] = fESDEvent->GetPrimaryVertex()->GetY();
  vtx[2] = fESDEvent->GetPrimaryVertex()->GetZ();

  //Vtx class z-bin
  Int_t zvtx = (Int_t)((vtx[2]+10.)/2.) ;
  if(zvtx<0)zvtx=0 ;
  if(zvtx>9)zvtx=9 ;

  Double_t run = fESDEvent->GetRunNumber()+0.5;
  TString trigClasses = fESDEvent->GetFiredTriggerClasses();
  if(trigClasses.Contains("CINT1B-ABCE-NOPF-ALL"))
    FillHistogram("hRunTrigger",run,0.5) ;
  else
    FillHistogram("hRunTrigger",run,1.5) ;

  FillHistogram("hEvents",0.5) ;
  FillHistogram("hVtxBin",zvtx-0.5) ;
  FillHistogram("hRunEvents",run) ;

  //check if containers for mixed events exist
  //create new if necessary
  if(!fPHOSEvents[zvtx]) fPHOSEvents[zvtx]=new TList() ;
  if(!fEMCALEvents[zvtx]) fEMCALEvents[zvtx]=new TList() ;
  if(!fConvEvents[zvtx]) fConvEvents[zvtx]=new TList() ;

  Int_t nPHOS=fPHOSEvent->GetEntriesFast() ;
  Int_t nEMCAL=fEMCALEvent->GetEntriesFast() ;
  Int_t nConv = fConvEvent->GetEntriesFast() ;
  //Some QA histograms
  //Calculate number of good converion photons
  Int_t nConvGood=0 ;
  for(Int_t iConv = 0; iConv<nConv; iConv++){
    TLorentzVector * cnv = static_cast<TLorentzVector*>(fConvEvent->At(iConv)) ;
    if(cnv->TestBit(kConvOnFly) && !cnv->TestBit(kConvEta)){
      nConvGood++ ;
    }
  }
  FillHistogram("hRunConvs",run,double(nConvGood)) ;
  FillHistogram("hRunPHOS", run,double(nPHOS)) ;
  FillHistogram("hRunEMCAL",run,double(nEMCAL)) ;

  //Fill Real distributions
  for(Int_t iPHOS=0; iPHOS<nPHOS;iPHOS++){
    TLorentzVector * cal = static_cast<TLorentzVector*>(fPHOSEvent->At(iPHOS)) ;
    for(Int_t iConv = 0; iConv<nConv; iConv++){
      TLorentzVector * cnv = static_cast<TLorentzVector*>(fConvEvent->At(iConv)) ;
      TLorentzVector pi=*cal + *cnv ;
      Double_t alpha=TMath::Abs(cal->Energy()-cnv->Energy())/(cal->Energy()+cnv->Energy()) ;
      if(cnv->TestBit(kConvOnFly) && !cnv->TestBit(kConvEta)){
        FillHistogram("PHOS_Re_mvsPt_all",pi.M(),pi.Pt()) ;
        char keym[55] ;
        for(Int_t iw=0;iw<10;iw++){ //resolution
          for(Int_t in=0;in<10;in++){
            snprintf(keym,55,"hMC_nonlinearity_w%d_n%d",iw,in) ;
            Double_t mMod=0.,ptMod=0. ;
            Recalibrate(mMod, ptMod, cal, cnv, iw, in) ;
            FillHistogram(keym,mMod,ptMod) ;
            snprintf(keym,55,"hMC_nonlinearity_ConvPHOS_w%d_n%d",iw,in) ;
            RecalibrateConvPHOS(mMod, ptMod, cal, cnv, iw, in) ;
            FillHistogram(keym,mMod,ptMod) ;
 
          }
        }
        if(cal->TestBit(kCaloDistBad))
          FillHistogram("PHOS_Re_mvsPt_all_dist",pi.M(),pi.Pt()) ;
        FillHistogram("PHOS_Re_mvsPt_E",pi.M(),pi.Pt(),cal->Energy()) ;
        FillHistogram("PHOS_Re_mvsPt_alpha",pi.M(),pi.Pt(),alpha) ;
        if(cal->TestBit(kCaloPIDdisp))
          FillHistogram("PHOS_Re_mvsPt_Disp",pi.M(),pi.Pt()) ;
        if(cal->TestBit(kCaloPIDtof))
          FillHistogram("PHOS_Re_mvsPt_TOF",pi.M(),pi.Pt()) ;
        if(cal->TestBit(kCaloPIDneutral))
          FillHistogram("PHOS_Re_mvsPt_Neutral",pi.M(),pi.Pt()) ;
        if(cal->TestBit(kCaloPIDneutral) && cal->TestBit(kCaloPIDdisp))
          FillHistogram("PHOS_Re_mvsPt_DispNeutral",pi.M(),pi.Pt()) ;
      }
      //Vary Conversion cuts
      if(cnv->TestBit(kConvOnFly)){
        //Default
        if(!cnv->TestBit(kConvEta)){
          FillHistogram("PHOS_Re_mvsPt_OnFly",pi.M(),pi.Pt()) ;
          FillHistogram("PHOS_Re_mvsPt_OnFly_mult",pi.M(),pi.Pt(),fCentr) ;
        }
        if(cnv->TestBit(kConvdEdx) && !cnv->TestBit(kConvEta))
          FillHistogram("PHOS_Re_mvsPt_On_dEdx",pi.M(),pi.Pt()) ;
        if(!cnv->TestBit(kConvEta) && cnv->TestBit(kConvR))
          FillHistogram("PHOS_Re_mvsPt_On_R120",pi.M(),pi.Pt()) ;
        if(!cnv->TestBit(kConvEta) && cnv->TestBit(kConvZR))
          FillHistogram("PHOS_Re_mvsPt_On_Z",pi.M(),pi.Pt()) ;
        if( !cnv->TestBit(kConvEta) && cnv->TestBit(kConvNDF)) 
          FillHistogram("PHOS_Re_mvsPt_On_chi",pi.M(),pi.Pt()) ;
        if(1)
          FillHistogram("PHOS_Re_mvsPt_On_Eta",pi.M(),pi.Pt()) ;
        if( !cnv->TestBit(kConvEta) && cnv->TestBit(kConvPlan)){ 
          FillHistogram("PHOS_Re_mvsPt_On_Wcut",pi.M(),pi.Pt()) ;
          if(cal->TestBit(kCaloPIDneutral))
            FillHistogram("PHOS_Re_mvsPt_On_Wcut_Neu",pi.M(),pi.Pt()) ;
        }
        if( !cnv->TestBit(kConvEta) && cnv->TestBit(kConvArmQt)) 
          FillHistogram("PHOS_Re_mvsPt_On_ArmQt",pi.M(),pi.Pt()) ;
      }
      else{
        //Default
        if(!cnv->TestBit(kConvEta)){
          FillHistogram("PHOS_Re_mvsPt_Offline",pi.M(),pi.Pt()) ;
          FillHistogram("PHOS_Re_mvsPt_Offline_mult",pi.M(),pi.Pt(),fCentr) ;
        }
        if(cnv->TestBit(kConvdEdx) && !cnv->TestBit(kConvEta))
          FillHistogram("PHOS_Re_mvsPt_Off_dEdx",pi.M(),pi.Pt()) ;
        if(!cnv->TestBit(kConvEta) && cnv->TestBit(kConvR))
          FillHistogram("PHOS_Re_mvsPt_Off_R120",pi.M(),pi.Pt()) ;
        if(!cnv->TestBit(kConvEta) && cnv->TestBit(kConvZR))
          FillHistogram("PHOS_Re_mvsPt_Off_Z",pi.M(),pi.Pt()) ;
        if( !cnv->TestBit(kConvEta) && cnv->TestBit(kConvNDF))
          FillHistogram("PHOS_Re_mvsPt_Off_chi",pi.M(),pi.Pt()) ;
        if(1)
          FillHistogram("PHOS_Re_mvsPt_Off_Eta",pi.M(),pi.Pt()) ;
        if(!cnv->TestBit(kConvEta) && cnv->TestBit(kConvPlan)){
          FillHistogram("PHOS_Re_mvsPt_Off_Wcut",pi.M(),pi.Pt()) ;
          if(cal->TestBit(kCaloPIDneutral))
            FillHistogram("PHOS_Re_mvsPt_Off_Wcut_Neu",pi.M(),pi.Pt()) ;
        }
        if( !cnv->TestBit(kConvEta) && cnv->TestBit(kConvArmQt)) 
          FillHistogram("PHOS_Re_mvsPt_Off_ArmQt",pi.M(),pi.Pt()) ;
      }
    }
  }
  //PHOS module-dependent histograms
  for(Int_t iPHOS=0; iPHOS<nPHOS;iPHOS++){
    TLorentzVector * cal = static_cast<TLorentzVector*>(fPHOSEvent->At(iPHOS)) ;
    Int_t mod=1;
    while(!cal->TestBit(BIT(17+mod)) && mod<5)mod++ ;
    TString base("PHOS_Re_mvsPt_mod") ; base+=mod ;
    TString full ;
    for(Int_t iConv = 0; iConv<nConv; iConv++){
      TLorentzVector * cnv = static_cast<TLorentzVector*>(fConvEvent->At(iConv)) ;
      TLorentzVector pi=*cal + *cnv ;
      full=base ; full+="_single" ;
      if(cnv->TestBit(kConvOnFly) && !cnv->TestBit(kConvEta)){
        FillHistogram(full,pi.M(),cal->Pt()) ;
        full=base ; full+="_all" ;
        FillHistogram(full,pi.M(),pi.Pt()) ;
      }
    }
  }
  for(Int_t iEMCAL=0; iEMCAL<nEMCAL;iEMCAL++){
    TLorentzVector * cal = static_cast<TLorentzVector*>(fEMCALEvent->At(iEMCAL)) ;
    Int_t mod=0;
    while(!cal->TestBit(BIT(17+mod)) && mod<6)mod++ ;
    TString base("EMCAL_Re_mvsPt_mod") ; base+=mod ;
    TString full ;
    for(Int_t iConv = 0; iConv<nConv; iConv++){
      TLorentzVector * cnv = static_cast<TLorentzVector*>(fConvEvent->At(iConv)) ;
      TLorentzVector pi=*cal + *cnv ;
//      Double_t alpha=TMath::Abs(cal->Energy()-cnv->Energy())/(cal->Energy()+cnv->Energy()) ;
      full=base+"_single" ;
      if(cnv->TestBit(kConvOnFly) && !cnv->TestBit(kConvEta)){
        FillHistogram(full,pi.M(),cal->Pt()) ;
        full=base+"_all" ;
        FillHistogram(full,pi.M(),pi.Pt()) ;
        char keym[55] ;
        for(Int_t iw=0;iw<10;iw++){ //resolution
          for(Int_t in=0;in<10;in++){
            snprintf(keym,55,"hMC_nonlinearity_EMCAL_w%d_n%d",iw,in) ;
            Double_t mMod=0.,ptMod=0. ;
            RecalibrateEMCAL(mMod, ptMod, cal, cnv, iw, in) ;
            FillHistogram(keym,mMod,ptMod) ;
          }
        }
      }
      //Vary Conversion cuts
      if(cnv->TestBit(kConvOnFly)){
        //Default
        if(!cnv->TestBit(kConvEta)){
          FillHistogram("EMCAL_Re_mvsPt_OnFly",pi.M(),pi.Pt()) ;
          FillHistogram("EMCAL_Re_mvsPt_OnFly_mult",pi.M(),pi.Pt(),fCentr) ;
        }
        if(cnv->TestBit(kConvdEdx) && !cnv->TestBit(kConvEta))
          FillHistogram("EMCAL_Re_mvsPt_On_dEdx",pi.M(),pi.Pt()) ;
        if(!cnv->TestBit(kConvEta) && cnv->TestBit(kConvR))
          FillHistogram("EMCAL_Re_mvsPt_On_R120",pi.M(),pi.Pt()) ;
        if(!cnv->TestBit(kConvEta) && cnv->TestBit(kConvZR))
          FillHistogram("EMCAL_Re_mvsPt_On_Z",pi.M(),pi.Pt()) ;
        if(!cnv->TestBit(kConvEta) && cnv->TestBit(kConvNDF))
          FillHistogram("EMCAL_Re_mvsPt_On_chi",pi.M(),pi.Pt()) ;
        if(1)
          FillHistogram("EMCAL_Re_mvsPt_On_Eta",pi.M(),pi.Pt()) ;
        if(!cnv->TestBit(kConvEta) && cnv->TestBit(kConvPlan))
          FillHistogram("EMCAL_Re_mvsPt_On_Wcut",pi.M(),pi.Pt()) ;
        if(!cnv->TestBit(kConvEta) && cnv->TestBit(kConvArmQt))
          FillHistogram("EMCAL_Re_mvsPt_On_ArmQt",pi.M(),pi.Pt()) ;
      }
      else{
        //Default
        if(!cnv->TestBit(kConvEta)){
          FillHistogram("EMCAL_Re_mvsPt_Offline",pi.M(),pi.Pt()) ;
          FillHistogram("EMCAL_Re_mvsPt_Offline_mult",pi.M(),pi.Pt(),fCentr) ;
        }
        if(cnv->TestBit(kConvdEdx) && !cnv->TestBit(kConvEta))
          FillHistogram("EMCAL_Re_mvsPt_Off_dEdx",pi.M(),pi.Pt()) ;
        if(!cnv->TestBit(kConvEta) && cnv->TestBit(kConvR))
          FillHistogram("EMCAL_Re_mvsPt_Off_R120",pi.M(),pi.Pt()) ;
        if(!cnv->TestBit(kConvEta) && cnv->TestBit(kConvZR))
          FillHistogram("EMCAL_Re_mvsPt_Off_Z",pi.M(),pi.Pt()) ;
        if(!cnv->TestBit(kConvEta) && cnv->TestBit(kConvNDF))
          FillHistogram("EMCAL_Re_mvsPt_Off_chi",pi.M(),pi.Pt()) ;
        if(1)
          FillHistogram("EMCAL_Re_mvsPt_Off_Eta",pi.M(),pi.Pt()) ;
        if(!cnv->TestBit(kConvEta) && cnv->TestBit(kConvPlan))
          FillHistogram("EMCAL_Re_mvsPt_Off_Wcut",pi.M(),pi.Pt()) ;
        if(!cnv->TestBit(kConvEta) && cnv->TestBit(kConvArmQt))
          FillHistogram("EMCAL_Re_mvsPt_Off_ArmQt",pi.M(),pi.Pt()) ;
      }
    }
  }
  //Now fill mixed
  TList * prevPHOS = fPHOSEvents[zvtx] ;
  TList * prevEMCAL = fEMCALEvents[zvtx] ;
  TList * prevConv = fConvEvents[zvtx] ;
 
  //PHOS
  for(Int_t iPHOS=0; iPHOS<nPHOS;iPHOS++){
    TLorentzVector * cal = static_cast<TLorentzVector*>(fPHOSEvent->At(iPHOS)) ;
    for(Int_t ev=0; ev<prevConv->GetSize();ev++){
      TClonesArray * mixConv = static_cast<TClonesArray*>(prevConv->At(ev)) ;
      for(Int_t iConv = 0; iConv<mixConv->GetEntriesFast(); iConv++){
        TLorentzVector * cnv = static_cast<TLorentzVector*>(mixConv->At(iConv)) ;
        TLorentzVector pi=*cal + *cnv ;
        if(cnv->TestBit(kConvOnFly) && !cnv->TestBit(kConvEta)){
          FillHistogram("PHOS_Mi_mvsPt_all",pi.M(),pi.Pt()) ;
          if(cal->TestBit(kCaloPIDdisp))
            FillHistogram("PHOS_Mi_mvsPt_Disp",pi.M(),pi.Pt()) ;
          if(cal->TestBit(kCaloPIDtof))
            FillHistogram("PHOS_Mi_mvsPt_TOF",pi.M(),pi.Pt()) ;
            if(cal->TestBit(kCaloPIDneutral))
            FillHistogram("PHOS_Mi_mvsPt_Neutral",pi.M(),pi.Pt()) ;
          if(cal->TestBit(kCaloPIDneutral) && cal->TestBit(kCaloPIDdisp))
            FillHistogram("PHOS_Mi_mvsPt_DispNeutral",pi.M(),pi.Pt()) ;
        }
        //Vary Conversion cuts
        if(cnv->TestBit(kConvOnFly)){
          //Default
          if(!cnv->TestBit(kConvEta))
            FillHistogram("PHOS_Mi_mvsPt_OnFly",pi.M(),pi.Pt()) ;
          if( cnv->TestBit(kConvdEdx) && !cnv->TestBit(kConvEta))
            FillHistogram("PHOS_Mi_mvsPt_On_dEdx",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvEta) && cnv->TestBit(kConvR))
            FillHistogram("PHOS_Mi_mvsPt_On_R120",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvEta) && cnv->TestBit(kConvZR))
            FillHistogram("PHOS_Mi_mvsPt_On_Z",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvEta) && cnv->TestBit(kConvNDF))
            FillHistogram("PHOS_Mi_mvsPt_On_chi",pi.M(),pi.Pt()) ;
          if(1)
            FillHistogram("PHOS_Mi_mvsPt_On_Eta",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvEta) && cnv->TestBit(kConvPlan))
            FillHistogram("PHOS_Mi_mvsPt_On_Wcut",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvEta) && cnv->TestBit(kConvArmQt))
            FillHistogram("PHOS_Mi_mvsPt_On_ArmQt",pi.M(),pi.Pt()) ;
        }
        else{
          //Default
          if(!cnv->TestBit(kConvEta))
            FillHistogram("PHOS_Mi_mvsPt_Offline",pi.M(),pi.Pt()) ;
          if(cnv->TestBit(kConvdEdx) && !cnv->TestBit(kConvEta))
            FillHistogram("PHOS_Mi_mvsPt_Off_dEdx",pi.M(),pi.Pt()) ;
          if( !cnv->TestBit(kConvEta) && cnv->TestBit(kConvR))
            FillHistogram("PHOS_Mi_mvsPt_Off_R120",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvEta) && cnv->TestBit(kConvZR))
            FillHistogram("PHOS_Mi_mvsPt_Off_Z",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvEta) && cnv->TestBit(kConvNDF))
            FillHistogram("PHOS_Mi_mvsPt_Off_chi",pi.M(),pi.Pt()) ;
          if(1)
            FillHistogram("PHOS_Mi_mvsPt_Off_Eta",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvEta) && cnv->TestBit(kConvPlan))
            FillHistogram("PHOS_Mi_mvsPt_Off_Wcut",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvEta) && cnv->TestBit(kConvArmQt))
            FillHistogram("PHOS_Mi_mvsPt_Off_ArmQt",pi.M(),pi.Pt()) ;
        }
      }
    }
  }
  for(Int_t iConv = 0; iConv<nConv; iConv++){
    TLorentzVector * cnv = static_cast<TLorentzVector*>(fConvEvent->At(iConv)) ;
    for(Int_t ev=0; ev<prevPHOS->GetSize();ev++){
      TClonesArray * mixPHOS = static_cast<TClonesArray*>(prevPHOS->At(ev)) ;
      for(Int_t iPHOS=0; iPHOS<mixPHOS->GetEntriesFast();iPHOS++){
        TLorentzVector * cal = static_cast<TLorentzVector*>(mixPHOS->At(iPHOS)) ;
        TLorentzVector pi=*cal + *cnv ;
        if(cnv->TestBit(kConvOnFly) && !cnv->TestBit(kConvEta)){
          FillHistogram("PHOS_Mi_mvsPt_all",pi.M(),pi.Pt()) ;
          if(cal->TestBit(kCaloPIDdisp))
            FillHistogram("PHOS_Mi_mvsPt_Disp",pi.M(),pi.Pt()) ;
          if(cal->TestBit(kCaloPIDtof))
            FillHistogram("PHOS_Mi_mvsPt_TOF",pi.M(),pi.Pt()) ;
          if(cal->TestBit(kCaloPIDneutral))
            FillHistogram("PHOS_Mi_mvsPt_Neutral",pi.M(),pi.Pt()) ;
          if(cal->TestBit(kCaloPIDneutral) && cal->TestBit(kCaloPIDdisp))
            FillHistogram("PHOS_Mi_mvsPt_DispNeutral",pi.M(),pi.Pt()) ;
        }
        //Vary Conversion cuts
        if(cnv->TestBit(kConvOnFly)){
          //Default
          if(!cnv->TestBit(kConvEta))
            FillHistogram("PHOS_Mi_mvsPt_OnFly",pi.M(),pi.Pt()) ;
          if(cnv->TestBit(kConvdEdx) && !cnv->TestBit(kConvEta))
            FillHistogram("PHOS_Mi_mvsPt_On_dEdx",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvEta) && cnv->TestBit(kConvR))
            FillHistogram("PHOS_Mi_mvsPt_On_R120",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvEta) && cnv->TestBit(kConvZR))
            FillHistogram("PHOS_Mi_mvsPt_On_Z",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvEta) && cnv->TestBit(kConvNDF))
            FillHistogram("PHOS_Mi_mvsPt_On_chi",pi.M(),pi.Pt()) ;
          if(1)
            FillHistogram("PHOS_Mi_mvsPt_On_Eta",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvEta) && cnv->TestBit(kConvPlan))
            FillHistogram("PHOS_Mi_mvsPt_On_Wcut",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvEta) && cnv->TestBit(kConvArmQt))
            FillHistogram("PHOS_Mi_mvsPt_On_ArmQt",pi.M(),pi.Pt()) ;
        }
        else{
          //Default
          if(!cnv->TestBit(kConvEta))
            FillHistogram("PHOS_Mi_mvsPt_Offline",pi.M(),pi.Pt()) ;
          if(cnv->TestBit(kConvdEdx) && !cnv->TestBit(kConvEta))
            FillHistogram("PHOS_Mi_mvsPt_Off_dEdx",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvEta) && cnv->TestBit(kConvR))
            FillHistogram("PHOS_Mi_mvsPt_Off_R120",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvEta) && cnv->TestBit(kConvZR))
            FillHistogram("PHOS_Mi_mvsPt_Off_Z",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvEta) && cnv->TestBit(kConvNDF))
            FillHistogram("PHOS_Mi_mvsPt_Off_chi",pi.M(),pi.Pt()) ;
          if(1)
            FillHistogram("PHOS_Mi_mvsPt_Off_Eta",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvEta) && cnv->TestBit(kConvPlan))
            FillHistogram("PHOS_Mi_mvsPt_Off_Wcut",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvEta) && cnv->TestBit(kConvArmQt))
            FillHistogram("PHOS_Mi_mvsPt_Off_ArmQt",pi.M(),pi.Pt()) ;
        }
      }
    }
  }
 
/*
  //PHOS module dependent
  for(Int_t iPHOS=0; iPHOS<nPHOS;iPHOS++){
    TLorentzVector * cal = static_cast<TLorentzVector*>(fPHOSEvent->At(iPHOS)) ;
    Int_t mod=1;
    while(!cal->TestBit(BIT(17+mod)) && mod<5)mod++ ;
    TString base("PHOS_Mi_mvsPt_mod") ; base+=mod ;
    TString full ;
    for(Int_t ev=0; ev<prevConv->GetSize();ev++){
      TClonesArray * mixConv = static_cast<TClonesArray*>(prevConv->At(ev)) ;
      for(Int_t iConv = 0; iConv<mixConv->GetEntriesFast(); iConv++){
        TLorentzVector * cnv = static_cast<TLorentzVector*>(mixConv->At(iConv)) ;
        TLorentzVector pi=*cal + *cnv ;
        full=base+"_all" ;
        if(cnv->TestBit(kConvOnFly) && !cnv->TestBit(kConvEta)){
          FillHistogram(full,pi.M(),pi.Pt()) ;
        }
      }
    }
  }
  for(Int_t iConv = 0; iConv<nConv; iConv++){
    TLorentzVector * cnv = static_cast<TLorentzVector*>(fConvEvent->At(iConv)) ;
    for(Int_t ev=0; ev<prevPHOS->GetSize();ev++){
      TClonesArray * mixPHOS = static_cast<TClonesArray*>(prevPHOS->At(ev)) ;
      for(Int_t iPHOS=0; iPHOS<mixPHOS->GetEntriesFast();iPHOS++){
        TLorentzVector * cal = static_cast<TLorentzVector*>(mixPHOS->At(iPHOS)) ;
        Int_t mod=1;
        while(!cal->TestBit(BIT(17+mod)) && mod<5)mod++ ;
        TString base("PHOS_Mi_mvsPt_mod") ; base+=mod ;
        TString full ;
        TLorentzVector pi=*cal + *cnv ;
        full=base+"_all" ;
        if(cnv->TestBit(kConvOnFly) && !cnv->TestBit(kConvEta)){
          FillHistogram(full,pi.M(),pi.Pt()) ;
        }
      }
    }
  }
*/

  //EMCAL
  for(Int_t iEMCAL=0; iEMCAL<nEMCAL;iEMCAL++){
    TLorentzVector * cal = static_cast<TLorentzVector*>(fEMCALEvent->At(iEMCAL)) ;
    Int_t mod=0;
    while(!cal->TestBit(BIT(17+mod)) && mod<6)mod++ ;
    TString base("EMCAL_Mi_mvsPt_mod") ; base+=mod ;
    TString full ;
    for(Int_t ev=0; ev<prevConv->GetSize();ev++){
      TClonesArray * mixConv = static_cast<TClonesArray*>(prevConv->At(ev)) ;
      for(Int_t iConv = 0; iConv<mixConv->GetEntriesFast(); iConv++){
        TLorentzVector * cnv = static_cast<TLorentzVector*>(mixConv->At(iConv)) ;
        TLorentzVector pi=*cal + *cnv ;
        full=base+"_all" ;
        if(cnv->TestBit(kConvOnFly) && !cnv->TestBit(kConvEta)){
          FillHistogram(full,pi.M(),pi.Pt()) ;
        } 
        if(cnv->TestBit(kConvOnFly)){
          //Default
          if(!cnv->TestBit(kConvEta))
            FillHistogram("EMCAL_Mi_mvsPt_OnFly",pi.M(),pi.Pt()) ;
          if(cnv->TestBit(kConvdEdx) && !cnv->TestBit(kConvEta))
            FillHistogram("EMCAL_Mi_mvsPt_On_dEdx",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvEta) && cnv->TestBit(kConvR))
            FillHistogram("EMCAL_Mi_mvsPt_On_R120",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvEta) && cnv->TestBit(kConvZR))
            FillHistogram("EMCAL_Mi_mvsPt_On_Z",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvEta) && cnv->TestBit(kConvNDF))
            FillHistogram("EMCAL_Mi_mvsPt_On_chi",pi.M(),pi.Pt()) ;
          if(1)
            FillHistogram("EMCAL_Mi_mvsPt_On_Eta",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvEta) && cnv->TestBit(kConvPlan))
            FillHistogram("EMCAL_Mi_mvsPt_On_Wcut",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvEta) && cnv->TestBit(kConvArmQt))
            FillHistogram("EMCAL_Mi_mvsPt_On_ArmQt",pi.M(),pi.Pt()) ;
        }
        else{
          //Default
          if(!cnv->TestBit(kConvEta))
            FillHistogram("EMCAL_Mi_mvsPt_Offline",pi.M(),pi.Pt()) ;
          if( cnv->TestBit(kConvdEdx) && !cnv->TestBit(kConvEta))
            FillHistogram("EMCAL_Mi_mvsPt_Off_dEdx",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvEta) && cnv->TestBit(kConvR))
            FillHistogram("EMCAL_Mi_mvsPt_Off_R120",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvEta) && cnv->TestBit(kConvZR))
            FillHistogram("EMCAL_Mi_mvsPt_Off_Z",pi.M(),pi.Pt()) ;
          if( !cnv->TestBit(kConvEta) && cnv->TestBit(kConvNDF))
            FillHistogram("EMCAL_Mi_mvsPt_Off_chi",pi.M(),pi.Pt()) ;
          if(1)
            FillHistogram("EMCAL_Mi_mvsPt_Off_Eta",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvEta) && cnv->TestBit(kConvPlan))
            FillHistogram("EMCAL_Mi_mvsPt_Off_Wcut",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvEta) && cnv->TestBit(kConvArmQt))
            FillHistogram("EMCAL_Mi_mvsPt_Off_ArmQt",pi.M(),pi.Pt()) ;
        }
      }
    }
  }
/*
  for(Int_t iConv = 0; iConv<nConv; iConv++){
    TLorentzVector * cnv = static_cast<TLorentzVector*>(fConvEvent->At(iConv)) ;
    for(Int_t ev=0; ev<prevEMCAL->GetSize();ev++){
      TClonesArray * mixEMCAL = static_cast<TClonesArray*>(prevEMCAL->At(ev)) ;
      for(Int_t iEMCAL=0; iEMCAL<mixEMCAL->GetEntriesFast();iEMCAL++){
        TLorentzVector * cal = static_cast<TLorentzVector*>(mixEMCAL->At(iEMCAL)) ;
        Int_t mod=0;
        while(!cal->TestBit(BIT(17+mod)) && mod<6)mod++ ;
        TString base("EMCAL_Mi_mvsPt_mod") ; base+=mod ;
        TString full ;
        TLorentzVector pi=*cal + *cnv ;
        full=base+"_all" ;
        if(cnv->TestBit(kConvOnFly) && !cnv->TestBit(kConvEta)){
          FillHistogram(full,pi.M(),pi.Pt()) ;
          if(cal->TestBit(kCaloPIDdisp)){
            full=base+"_Disp" ;
            FillHistogram(full,pi.M(),pi.Pt()) ;
          }
          if(cal->TestBit(kCaloPIDtof)){
            full=base+"_TOF" ;
            FillHistogram(full,pi.M(),pi.Pt()) ;
          }
          if(cal->TestBit(kCaloPIDneutral)){
            full=base+"_Neutral" ;
              FillHistogram(full,pi.M(),pi.Pt()) ;
          }
          if(cal->TestBit(kCaloPIDneutral) && cal->TestBit(kCaloPIDdisp)){
            full=base+"_DispNeutral" ;
            FillHistogram(full,pi.M(),pi.Pt()) ;
          }
        }
        if(cnv->TestBit(kConvOnFly)){
          //Default
          if(!cnv->TestBit(kConvEta))
            FillHistogram("EMCAL_Mi_mvsPt_OnFly",pi.M(),pi.Pt()) ;
          if(cnv->TestBit(kConvdEdx) && !cnv->TestBit(kConvEta))
            FillHistogram("EMCAL_Mi_mvsPt_On_dEdx",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvEta) && cnv->TestBit(kConvR))
            FillHistogram("EMCAL_Mi_mvsPt_On_R120",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvEta) && cnv->TestBit(kConvZR))
            FillHistogram("EMCAL_Mi_mvsPt_On_Z",pi.M(),pi.Pt()) ;
          if( !cnv->TestBit(kConvEta) && cnv->TestBit(kConvNDF))
            FillHistogram("EMCAL_Mi_mvsPt_On_chi",pi.M(),pi.Pt()) ;
          if(1)
            FillHistogram("EMCAL_Mi_mvsPt_On_Eta",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvEta) && cnv->TestBit(kConvPlan)) 
            FillHistogram("EMCAL_Mi_mvsPt_On_Wcut",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvEta) && cnv->TestBit(kConvArmQt)) 
            FillHistogram("EMCAL_Mi_mvsPt_On_ArmQt",pi.M(),pi.Pt()) ;
        }
        else{
          //Default
          if(!cnv->TestBit(kConvEta))
            FillHistogram("EMCAL_Mi_mvsPt_Offline",pi.M(),pi.Pt()) ;
          if(cnv->TestBit(kConvdEdx) && !cnv->TestBit(kConvEta))
            FillHistogram("EMCAL_Mi_mvsPt_Off_dEdx",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvEta) && cnv->TestBit(kConvR))
            FillHistogram("EMCAL_Mi_mvsPt_Off_R120",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvEta) && cnv->TestBit(kConvZR)) 
            FillHistogram("EMCAL_Mi_mvsPt_Off_Z",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvEta) && cnv->TestBit(kConvNDF))
            FillHistogram("EMCAL_Mi_mvsPt_Off_chi",pi.M(),pi.Pt()) ;
          if(1)
            FillHistogram("EMCAL_Mi_mvsPt_Off_Eta",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvEta) && cnv->TestBit(kConvPlan))
            FillHistogram("EMCAL_Mi_mvsPt_Off_Wcut",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvEta) && cnv->TestBit(kConvArmQt)) 
            FillHistogram("EMCAL_Mi_mvsPt_Off_ArmQt",pi.M(),pi.Pt()) ;
        }
      }
    }
  }
*/


      
  //Now we either add current events to stack or remove
  //Here we have some difficulty: conversion and calorimeter photons have different spectra.
  //So to correctly reproduce combinatorial background we have to preserve average number of 
  //photons of both kinds per event. Therefore we should not reject empty PHOS/EMCAL events
  //though it will cost some memory. Reject only those events where no photons anywhere

  if((fPHOSEvent->GetEntriesFast()>0 || fEMCALEvent->GetEntriesFast()>0) && fConvEvent->GetEntriesFast()>0){
    prevPHOS->AddFirst(fPHOSEvent) ;
    fPHOSEvent=0; 
    prevEMCAL->AddFirst(fEMCALEvent) ;
    fEMCALEvent=0 ;
    prevConv->AddFirst(fConvEvent) ;
    fConvEvent=0 ;
    if(prevPHOS->GetSize()>100){//Remove redundant events
      TClonesArray * tmp = static_cast<TClonesArray*>(prevPHOS->Last()) ;
      prevPHOS->RemoveLast() ;
      delete tmp ;
      tmp = static_cast<TClonesArray*>(prevEMCAL->Last()) ;
      prevEMCAL->RemoveLast() ;
      delete tmp ;
      tmp = static_cast<TClonesArray*>(prevConv->Last()) ;
      prevConv->RemoveLast() ;
      delete tmp ;
    }
  }

}
//___________________________________________________________________________
void AliAnalysisTaskCaloConv::ProcessMC(){
  //ProcessMC
  //fill histograms for efficiensy etc. calculation
  if(!fStack) return ;
  
  const Double_t rcut = 1. ; //cut for primary particles
  Double_t vtx[3];
  vtx[0] = fESDEvent->GetPrimaryVertex()->GetX();
  vtx[1] = fESDEvent->GetPrimaryVertex()->GetY();
  vtx[2] = fESDEvent->GetPrimaryVertex()->GetZ();

  Int_t nPHOS=fPHOSEvent->GetEntriesFast() ;
  Int_t nEMCAL=fEMCALEvent->GetEntriesFast() ;
  Int_t nConv = fConvEvent->GetEntriesFast() ;
 
  //---------First pi0/eta-----------------------------
  char partName[10] ;
  char hkey[55] ;
  for (Int_t iTracks = 0; iTracks < fStack->GetNtrack(); iTracks++) {
    TParticle* particle = (TParticle *)fStack->Particle(iTracks);
    if(particle->GetPdgCode() == 111)
      snprintf(partName,10,"pi0") ;
    else
      if(particle->GetPdgCode() == 221)
        snprintf(partName,10,"eta") ;
      else
        continue ;
     
    //Primary particle
    if(particle->R() >rcut)
      continue ;

    Double_t pt = particle->Pt() ;
    //Total number of pi0 with creation radius <1 cm
    snprintf(hkey,55,"hMC_CaloConv_all%s",partName) ;
    FillHistogram(hkey,pt) ;
    if(TMath::Abs(particle->Y())<1.){
      snprintf(hkey,55,"hMC_CaloConv_%s_unitEta",partName) ;
      FillHistogram(hkey,pt,fCentr) ;
    }

    //Check if one of photons converted
    if(particle->GetNDaughters()!=2)
     continue ; //Do not account Dalitz decays

    TParticle * gamma1 = fStack->Particle(particle->GetFirstDaughter());
    TParticle * gamma2 = fStack->Particle(particle->GetLastDaughter());
    //Number of pi0s decayed into acceptance
    Bool_t inAcc1 = (TMath::Abs(gamma1->Eta())<0.9) ;
    Bool_t inAcc2 = (TMath::Abs(gamma2->Eta())<0.9) ;
    Int_t mod1,mod2 ;
    Double_t x=0.,z=0. ;
    Bool_t hitPHOS1 = fPHOSgeom->ImpactOnEmc(gamma1, mod1, z,x) ;
    Bool_t hitPHOS2 = fPHOSgeom->ImpactOnEmc(gamma2, mod2, z,x) ;
    Bool_t hitEMCAL1= fEMCALgeom->Impact(gamma1) ;
    Bool_t hitEMCAL2= fEMCALgeom->Impact(gamma2) ;
 
    Bool_t goodPair=kFALSE ;
    if((inAcc1 && hitPHOS2) || (inAcc2 && hitPHOS1)){
      snprintf(hkey,55,"hMC_CaloConv_%sPHOSacc",partName) ;
      FillHistogram(hkey,pt) ;
      goodPair=kTRUE ;
    } 
    if((inAcc1 && hitEMCAL2) || (inAcc2 && hitEMCAL1)){
      snprintf(hkey,55,"hMC_CaloConv_%sEMCALacc",partName) ;
      FillHistogram(hkey,pt) ;
       goodPair=kTRUE ;
    }
    if(!goodPair){
      continue ;
    }
 
    Bool_t converted1 = kFALSE ;
    if(gamma1->GetNDaughters()==2){
      TParticle * e1=fStack->Particle(gamma1->GetFirstDaughter()) ;
      TParticle * e2=fStack->Particle(gamma1->GetLastDaughter()) ;
      if(TMath::Abs(e1->GetPdgCode())==11 && TMath::Abs(e2->GetPdgCode())==11){ //conversion
        if(e1->R()<180.)
          converted1 = kTRUE ;
      }
    }
    Bool_t converted2 = kFALSE ;
    if(gamma2->GetNDaughters()==2){
      TParticle * e1=fStack->Particle(gamma2->GetFirstDaughter()) ;
      TParticle * e2=fStack->Particle(gamma2->GetLastDaughter()) ;
      if(TMath::Abs(e1->GetPdgCode())==11 && TMath::Abs(e2->GetPdgCode())==11){ //conversion
        if(e1->R()<180.)
          converted2 = kTRUE ;
      }
    }
  
    //Number of pi0s with one photon converted
    if((converted1 && !converted2 && hitPHOS2) || (!converted1 && hitPHOS1 && converted2)) {
      snprintf(hkey,55,"hMC_CaloConv_%s_PHOS_conv",partName) ;
      FillHistogram(hkey,pt) ;
    }
 
    if((converted1 && !converted2 && hitEMCAL2) || (!converted1 && hitEMCAL1 && converted2)) {
      snprintf(hkey,55,"hMC_CaloConv_%s_EMCAL_conv",partName) ;
      FillHistogram(hkey,pt) ;
    }
 
    //Both converted
    if(converted1 && converted2) {
      snprintf(hkey,55,"hMC_CaloConv_%s_bothphot_conv",partName) ;
      FillHistogram(hkey,pt) ;
        continue ;
    }
 
    //photon pointing calorimeter converted
    if((converted1 && hitPHOS1 && !hitEMCAL2) || (converted2 && hitPHOS2 && !hitEMCAL1) || 
       (converted1 && hitEMCAL1 && !hitPHOS2) || (converted2 && hitEMCAL2 && !hitPHOS1)){
      snprintf(hkey,55,"hMC_CaloConv_%s_convPhotInCalo",partName) ;
      FillHistogram(hkey,pt) ;
       continue ;
    }
   
    //Converted pi0 with v0 and photon PHOS or EMCAL
    Bool_t foundV01onfly=kFALSE, foundV01offline=kFALSE, foundV02onfly=kFALSE, foundV02offline=kFALSE ;
    Bool_t foundV01onflyPID=kFALSE, foundV01offlinePID=kFALSE, foundV02onflyPID=kFALSE, foundV02offlinePID=kFALSE ;
    TLorentzVector pConvOn,pConvOff ;
    for(Int_t iv0=0; iv0<fESDEvent->GetNumberOfV0s();iv0++){
      AliESDv0 * v0 = fESDEvent->GetV0(iv0) ;
 
      TParticle * negativeMC = fStack->Particle(TMath::Abs(fESDEvent->GetTrack(v0->GetNindex())->GetLabel()));
      TParticle * positiveMC = fStack->Particle(TMath::Abs(fESDEvent->GetTrack(v0->GetPindex())->GetLabel()));
 
      if(negativeMC->GetMother(0) != positiveMC->GetMother(0))
        continue ;
 
      if(TMath::Abs(negativeMC->GetPdgCode())!=11 || TMath::Abs(positiveMC->GetPdgCode())!=11){
        continue;
       }
      if(negativeMC->GetPdgCode()==positiveMC->GetPdgCode()){
        continue;
      }
  
      TParticle * v0Gamma = fStack->Particle(negativeMC->GetMother(0));
      Bool_t same = (v0Gamma == gamma1) ;
      TParticle * tmp = v0Gamma ;
      while(!same && tmp->GetFirstMother()>=0){
        tmp = fStack->Particle(tmp->GetFirstMother());
       same = (tmp == gamma1) ;
     }
     if(same){
       if(v0->GetOnFlyStatus())
         foundV01onfly = kTRUE ;
       else
         foundV01offline= kTRUE ;
       for(Int_t iconv=0; iconv<nConv;iconv++){
         if(fGammaV0s[iconv] == iv0){
           TLorentzVector * cnv = static_cast<TLorentzVector*>(fConvEvent->At(iconv)) ;
           //default cuts
           if(!cnv->TestBit(kConvEta)){
             if(v0->GetOnFlyStatus()){
               pConvOn= *cnv ; 
               foundV01onflyPID = kTRUE ;
             }
             else{
               pConvOff= *cnv ;
               foundV01offlinePID = kTRUE ;
             }
           }
           break ;
         }
       }
       continue ;
     } 
     same = (v0Gamma == gamma2) ;
     tmp = v0Gamma ;
     while(!same && tmp->GetFirstMother()>=0){
       tmp = fStack->Particle(tmp->GetFirstMother());
       same = (tmp == gamma2) ;
     }
     if(same){
       if(v0->GetOnFlyStatus())
         foundV02onfly = kTRUE ;
       else
         foundV02offline = kTRUE ;
       for(Int_t iconv=0; iconv<nConv;iconv++){
         if(fGammaV0s[iconv] == iv0){
           TLorentzVector * cnv = static_cast<TLorentzVector*>(fConvEvent->At(iconv)) ;
           //default cuts
           if(!cnv->TestBit(kConvEta)){
             if(v0->GetOnFlyStatus()){
               pConvOn= *cnv ;
               foundV02onflyPID = kTRUE ;
             }
             else{
               pConvOff= *cnv ;
               foundV02offlinePID = kTRUE ;
             }
           }
           break ;
         }
       }
     } 
   }

   goodPair=kFALSE ;
   if((foundV01onfly && hitPHOS2) || (foundV02onfly && hitPHOS1)){
     snprintf(hkey,55,"hMC_CaloConv_%s_v0onfly_PHOSacc",partName) ;
     FillHistogram(hkey,pt) ;
     goodPair=kTRUE;
   }
   if((foundV01offline && hitPHOS2) || (foundV02offline && hitPHOS1)){
     snprintf(hkey,55,"hMC_CaloConv_%s_v0offline_PHOSacc",partName) ;
     FillHistogram(hkey,pt) ;
     goodPair=kTRUE;
   }
   if((foundV01onfly && hitEMCAL2) || (foundV02onfly && hitEMCAL1)){
     snprintf(hkey,55,"hMC_CaloConv_%s_v0onfly_EMCALacc",partName) ;
     FillHistogram(hkey,pt) ;
     goodPair=kTRUE;
   }
   if((foundV01offline && hitEMCAL2) || (foundV02offline && hitEMCAL1)){
     snprintf(hkey,55,"hMC_CaloConv_%s_v0offline_EMCALacc",partName) ;
     FillHistogram(hkey,pt) ;
     goodPair=kTRUE;
   }
   if(!goodPair){
     continue ;
   }

   //Converted pi0 with v0 and cluster in PHOS/EMCAL
   Bool_t cluInPHOS = kFALSE,cluInEMCAL=kFALSE ;
   Bool_t cluInPHOSpid = kFALSE,cluInEMCALpid=kFALSE ;
   Bool_t closeToBad= kFALSE ;
   TLorentzVector pCalo ;
   for (Int_t i=0; i<fESDEvent->GetNumberOfCaloClusters(); i++) {
     AliESDCaloCluster * clu = fESDEvent->GetCaloCluster(i);
     Int_t iprim = clu->GetLabel() ; //# of particle hit PHOS/EMCAL
     Bool_t matched = kFALSE ;
     while(iprim>=0 ) {
       if(iprim==particle->GetFirstDaughter() || iprim==particle->GetLastDaughter()){
         matched=kTRUE ;
         break ;
       }
       else{
         iprim=fStack->Particle(iprim)->GetFirstMother() ;
       }
     }
     if(!matched)
       continue ;
     if(clu->IsPHOS() && (hitPHOS1 || hitPHOS2)){
       cluInPHOS=kTRUE ;
       //Check if cluster passed PID
       for(Int_t inPHOS=0; inPHOS<nPHOS;inPHOS++){
         if(fGammaPHOS[inPHOS] == i){
           cluInPHOSpid=kTRUE ;
           break ;
         }
       }
       clu->GetMomentum(pCalo ,vtx);
       if(clu->GetDistanceToBadChannel()<fBadDistCutPHOS)
         closeToBad=kTRUE ;
       break ;
     }
     if(!clu->IsPHOS() && (hitEMCAL1 || hitEMCAL2)){
       cluInEMCAL=kTRUE ;
       //Check if cluster passed PID
       for(Int_t inEMCAL=0; inEMCAL<nEMCAL;inEMCAL++){
         if(fGammaEMCAL[inEMCAL] == i){
           cluInPHOSpid=kTRUE ;
           break ;
         }
       }
       clu->GetMomentum(pCalo ,vtx);
       if(clu->GetDistanceToBadChannel()<fBadDistCutEMCAL)
         closeToBad=kTRUE ;
       break ;
     }
   }

   if(cluInPHOS){
     //OnFly
     if(foundV01onfly ||foundV02onfly){
       snprintf(hkey,55,"hMC_CaloConv_%s_v0onfly_PHOSclu",partName) ;
       FillHistogram(hkey,pt) ;

       if((foundV01onflyPID ||foundV02onflyPID) && cluInPHOSpid){
       snprintf(hkey,55,"hMC_CaloConv_%s_v0onfly_PHOSclu_pid",partName) ;
         FillHistogram(hkey,pt) ;
         for(Int_t iw=0;iw<10;iw++){ //resolution
           for(Int_t in=0;in<10;in++){
             char keym[55] ;
             snprintf(keym,55,"hMC_CaloConv_%s_v0onfly_PHOSclu_w%d_n%d",partName,iw,in) ;
             Double_t mMod=0.,ptMod=0. ;
             Recalibrate(mMod, ptMod, &pCalo, &pConvOn, iw, in) ;
             FillHistogram(keym,ptMod) ;
             snprintf(keym,55,"hMC_CaloConv_%s_v0onfly_ConvPHOSclu_w%d_n%d",partName,iw,in) ;
             RecalibrateConvPHOS(mMod, ptMod, &pCalo, &pConvOn, iw, in) ;
             FillHistogram(keym,ptMod) ;
           }
         }
         if(!closeToBad){
           snprintf(hkey,55,"hMC_CaloConv_%s_v0onfly_PHOSclu_good",partName) ;
           FillHistogram(hkey,pt) ;
         }
         Double_t m=(pCalo+pConvOn).M() ;
         Double_t ptm=(pCalo+pConvOn).Pt() ;
         snprintf(hkey,55,"hMC_CaloConv_%s_v0on_PHOSclu_ptRec",partName) ;
         FillHistogram(hkey,ptm) ;
         snprintf(hkey,55,"hMC_CaloConv_%s_v0on_PHOSclu_mvsPt",partName) ;
         FillHistogram(hkey,m,ptm) ;
       }
     }

     //Offline
     if(foundV01offline ||foundV02offline){
       snprintf(hkey,55,"hMC_CaloConv_%s_v0offline_PHOSclu",partName) ;
       FillHistogram(hkey,pt) ;
       if((foundV01offlinePID ||foundV02offlinePID) && cluInPHOSpid){
         snprintf(hkey,55,"hMC_CaloConv_%s_v0offline_PHOSclu_pid",partName) ;
         FillHistogram(hkey,pt) ;
         Double_t m=(pCalo+pConvOff).M() ;
         Double_t ptm=(pCalo+pConvOff).Pt() ;
         snprintf(hkey,55,"hMC_CaloConv_%s_v0off_PHOSclu_ptRec",partName) ;
         FillHistogram(hkey,ptm) ;
         snprintf(hkey,55,"hMC_CaloConv_%s_v0off_PHOSclu_mvsPt",partName) ;
         FillHistogram(hkey,m,ptm) ;
         if(!closeToBad){
           snprintf(hkey,55,"hMC_CaloConv_%s_v0offline_PHOSclu_good",partName) ;
           FillHistogram(hkey,pt) ;
         }
       }
     } 

     if((foundV01onflyPID ||foundV02onflyPID) && cluInPHOSpid){
       TString base("hMC_CaloConv_") ; base+=partName; base+="_v0onfly_PHOSclu_mod" ;
       if(hitPHOS1)
         base+=mod1 ;
       else
         base+=mod2 ;
       FillHistogram(base.Data(),pt) ;
     }

     if((foundV01offlinePID ||foundV02offlinePID) && cluInPHOSpid){
       TString base("hMC_CaloConv_") ; base+=partName; base+="_v0offline_PHOSclu_mod" ;
       if(hitPHOS1)
         base+=mod1 ;
       else
         base+=mod2 ;
       FillHistogram(base.Data(),pt) ;
     }
   }
   if(cluInEMCAL && strcmp(partName,"pi0")==0){
     //OnFly
     if(foundV01onfly ||foundV02onfly){
       FillHistogram("hMC_CaloConv_pi0_v0onfly_EMCALclu",pt) ;

       if((foundV01onflyPID ||foundV02onflyPID) && cluInEMCALpid){
         FillHistogram("hMC_CaloConv_pi0_v0onfly_EMCALclu_pid",pt) ;
         for(Int_t iw=0;iw<10;iw++){ //resolution
           for(Int_t in=0;in<10;in++){
             char keym[55] ;
             snprintf(keym,55,"hMC_CaloConv_pi0_v0onfly_EMCALclu_w%d_n%d",iw,in) ;
             Double_t mMod=0.,ptMod=0. ;
             RecalibrateEMCAL(mMod, ptMod, &pCalo, &pConvOn, iw, in) ;
             FillHistogram(keym,ptMod) ;
           }
         }
         if(!closeToBad)
           FillHistogram("hMC_CaloConv_pi0_v0onfly_EMCALclu_good",pt) ;
         Double_t m=(pCalo+pConvOn).M() ;
         Double_t ptm=(pCalo+pConvOn).Pt() ;
         FillHistogram("hMC_CaloConv_pi0_v0on_EMCALclu_ptRec",ptm) ;
         FillHistogram("hMC_CaloConv_pi0_v0on_EMCALclu_mvsPt",m,ptm) ;
       }
     }

     //Offline
     if(foundV01offline ||foundV02offline){
       FillHistogram("hMC_CaloConv_pi0_v0offline_EMCALclu",pt) ;
       if((foundV01offlinePID ||foundV02offlinePID) && cluInEMCALpid){
         FillHistogram("hMC_CaloConv_pi0_v0offline_EMCALclu_pid",pt) ;
         Double_t m=(pCalo+pConvOff).M() ;
         Double_t ptm=(pCalo+pConvOff).Pt() ;
         FillHistogram("hMC_CaloConv_pi0_v0off_EMCALclu_ptRec",ptm) ;
         FillHistogram("hMC_CaloConv_pi0_v0off_EMCALclu_mvsPt",m,ptm) ;
         if(!closeToBad)
           FillHistogram("hMC_CaloConv_pi0_v0offline_EMCALclu_good",pt) ;
       }
     }
   }
  }

  //Construct Inv mass distributions for residual correlations
  if(fPHOSEvent && fConvEvent){
    for(Int_t iPHOS=0; iPHOS<fPHOSEvent->GetEntriesFast();iPHOS++){
      TLorentzVector * cal = static_cast<TLorentzVector*>(fPHOSEvent->At(iPHOS)) ;
      Int_t iclu=fGammaPHOS[iPHOS] ;
      AliESDCaloCluster * clu = fESDEvent->GetCaloCluster(iclu);
      Int_t iprimPHOS = clu->GetLabel() ; //# of particle hit PHOS/EMCAL
      for(Int_t iConv = 0; iConv<fConvEvent->GetEntriesFast(); iConv++){
        TLorentzVector * cnv = static_cast<TLorentzVector*>(fConvEvent->At(iConv)) ;
        if(!cnv->TestBit(kConvOnFly) || cnv->TestBit(kConvEta)) 
          continue;

        Int_t iv0=fGammaV0s[iConv] ;
        AliESDv0 * v0 = fESDEvent->GetV0(iv0) ;
        Int_t iprimNeg = TMath::Abs(fESDEvent->GetTrack(v0->GetNindex())->GetLabel()) ;
        Int_t iprimPos = TMath::Abs(fESDEvent->GetTrack(v0->GetPindex())->GetLabel()) ;

        //Check if there was a common ancistor
        Bool_t found = kFALSE ;
        Int_t curPHOS=iprimPHOS ;
        Int_t commonA=-1 ; 
        while(!found && curPHOS>-1){
          Int_t curNeg=iprimNeg ;
          while(!found && curNeg>-1){
            if(curNeg==curPHOS){
              found=kTRUE ;
              commonA=curPHOS ;
            }
            else{
              curNeg=fStack->Particle(curNeg)->GetFirstMother() ;
            }
          }
          curPHOS=fStack->Particle(curPHOS)->GetFirstMother() ;
        }
        found = kFALSE ;
        curPHOS=iprimPHOS ;
        Int_t commonB=-1 ;
        while(!found && curPHOS>-1){
          Int_t curPos=iprimPos ;
          while(!found && curPos>-1){
            if(curPos==curPHOS){
              found=kTRUE ;
              commonB=curPHOS ;
            }
            else{
              curPos=fStack->Particle(curPos)->GetFirstMother() ;
            }
          }
          curPHOS=fStack->Particle(curPHOS)->GetFirstMother() ;
        }
        if(commonA != commonB){
           //Strange
           AliInfo(Form("CommonA=%d, commonB=%d",commonA,commonB)) ; 
        }
        if(commonA>-1){//There was common particles
          Int_t pdg = fStack->Particle(commonA)->GetPdgCode() ;
          TLorentzVector pi=*cal + *cnv ;
          Double_t m=pi.M() ;
          Double_t pt=pi.Pt() ;
          Double_t alpha=TMath::Abs(cal->Energy()-cnv->Energy())/(cal->Energy()+cnv->Energy()) ;
          switch(pdg){
          case 11:
          case -11:
          case 22: //conversion
            FillHistogram("hMC_Resid_PHOS_Phot_mvsPt",m,pt,alpha) ;
            break ;
          case 111: //pi0
            FillHistogram("hMC_Resid_PHOS_Pi0_mvsPt",m,pt,alpha) ;
            break ;
          case 221: //eta
              FillHistogram("hMC_Resid_PHOS_eta_mvsPt",m,pt,alpha) ;
            break ;
          case 321: //K+
          case -321: //K-
          case 310:  //K0s
          case 130:  //K0L
            FillHistogram("hMC_Resid_PHOS_K_mvsPt",m,pt,alpha) ;
            break ;
          case 211:
          case -211: 
            FillHistogram("hMC_Resid_PHOS_pi_mvsPt",m,pt,alpha) ;
            break ;
          case -2212:  //pbar
          case -2112:  //nbar
            FillHistogram("hMC_Resid_PHOS_pbar_mvsPt",m,pt,alpha) ;
            break ;
          default: //else
            FillHistogram("hMC_Resid_PHOS_other_mvsPt",m,pt,alpha) ;
            break ;
          }
        }
       }
     }
   }
  

  if(fEMCALEvent && fConvEvent){
    for(Int_t iEMCAL=0; iEMCAL<fEMCALEvent->GetEntriesFast();iEMCAL++){
      TLorentzVector * cal = static_cast<TLorentzVector*>(fEMCALEvent->At(iEMCAL)) ;
      Int_t iclu=fGammaEMCAL[iEMCAL] ;
      AliESDCaloCluster * clu = fESDEvent->GetCaloCluster(iclu);
      Int_t iprimEMCAL = clu->GetLabel() ; //# of particle hit EMCAL
      for(Int_t iConv = 0; iConv<fConvEvent->GetEntriesFast(); iConv++){
        TLorentzVector * cnv = static_cast<TLorentzVector*>(fConvEvent->At(iConv)) ;
        if(!cnv->TestBit(kConvOnFly) || cnv->TestBit(kConvEta)) 
          continue;
        Int_t iv0=fGammaV0s[iConv] ;
        AliESDv0 * v0 = fESDEvent->GetV0(iv0) ;
        Int_t iprimNeg = TMath::Abs(fESDEvent->GetTrack(v0->GetNindex())->GetLabel()) ;
        Int_t iprimPos = TMath::Abs(fESDEvent->GetTrack(v0->GetPindex())->GetLabel()) ;
  
        //Check if there was a common ancistor
        Bool_t found = kFALSE ;
        Int_t curEMCAL=iprimEMCAL ;
        Int_t commonA=-1 ;
        while(!found && curEMCAL>-1){
          Int_t curNeg=iprimNeg ;
          while(!found && curNeg>-1){
            if(curNeg==curEMCAL){
              found=kTRUE ;
              commonA=curEMCAL ;
            }
            else{
              curNeg=fStack->Particle(curNeg)->GetFirstMother() ;
            }
          }
          curEMCAL=fStack->Particle(curEMCAL)->GetFirstMother() ;
        }
        found = kFALSE ;
        curEMCAL=iprimEMCAL ;
        Int_t commonB=-1 ;
        while(!found && curEMCAL>-1){
          Int_t curPos=iprimPos ;
          while(!found && curPos>-1){
            if(curPos==curEMCAL){
              found=kTRUE ;
              commonB=curEMCAL ;
            }
            else{
              curPos=fStack->Particle(curPos)->GetFirstMother() ;
            }
          }
          curEMCAL=fStack->Particle(curEMCAL)->GetFirstMother() ;
        }
  
        if(commonA != commonB){
           //Strange
           AliInfo(Form("CommonA=%d, commonB=%d",commonA,commonB)) ;
        }
        if(commonA>-1){//There was common particles
          Int_t pdg = fStack->Particle(commonA)->GetPdgCode() ;
          TLorentzVector pi=*cal + *cnv ;
          Double_t m=pi.M() ;
          Double_t pt=pi.Pt() ;
          Double_t alpha=TMath::Abs(cal->Energy()-cnv->Energy())/(cal->Energy()+cnv->Energy()) ;
          switch(pdg){
          case 11:
          case -11:
          case 22: //conversion
            FillHistogram("hMC_Resid_EMCAL_Phot_mvsPt",m,pt,alpha) ;
            break ;
          case 111: //pi0
            FillHistogram("hMC_Resid_EMCAL_Pi0_mvsPt",m,pt,alpha) ;
            break ;
          case 221: //eta
            FillHistogram("hMC_Resid_EMCAL_eta_mvsPt",m,pt,alpha) ;
            break ;
          case 321: //K+
          case -321: //K-
          case 310:  //K0s
          case 130:  //K0L
            FillHistogram("hMC_Resid_EMCAL_K_mvsPt",m,pt,alpha) ;
            break ;
          case 211:
          case -211:
            FillHistogram("hMC_Resid_EMCAL_pi_mvsPt",m,pt,alpha) ;
            break ;
          case -2212:  //pbar
          case -2112:  //nbar
            FillHistogram("hMC_Resid_EMCAL_pbar_mvsPt",m,pt,alpha) ;
            break ;
          default: //else
            FillHistogram("hMC_Resid_EMCAL_other_mvsPt",m,pt,alpha) ;
            break ;
          }
        }
       }
     }
   } 
   

   //------------- now photons ----------------
   for (Int_t iTracks = 0; iTracks < fStack->GetNtrack(); iTracks++) {
     TParticle* particle = (TParticle *)fStack->Particle(iTracks);
     if(particle->GetPdgCode() != 22)
       continue ;

     if(particle->R() >rcut)
       continue ;

     if(TMath::Abs(particle->Eta())>0.9)
       continue ;

     Double_t pt = particle->Pt() ;
     //Total number of pi0 with creation radius <1 cm
     FillHistogram("hMC_CaloConv_phot",pt) ;

     Int_t mod ;
     Double_t x=0.,z=0. ;
     Bool_t hitPHOS = fPHOSgeom->ImpactOnEmc(particle, mod, z,x) ;
     Bool_t hitEMCAL= fEMCALgeom->Impact(particle) ;

     //Photons in PHOS/EMCAL acceptance
     if(hitPHOS)
       FillHistogram("hMC_CaloConv_gammaPHOSacc",pt) ;
     if(hitEMCAL)
       FillHistogram("hMC_CaloConv_gammaEMCALacc",pt) ;

     //number of photons converted
     Bool_t converted = kFALSE ;
     if(particle->GetNDaughters()==2){
       TParticle * e1=fStack->Particle(particle->GetFirstDaughter()) ;
       TParticle * e2=fStack->Particle(particle->GetLastDaughter()) ;
       if(TMath::Abs(e1->GetPdgCode())==11 && TMath::Abs(e2->GetPdgCode())==11){ //conversion
         if(e1->R()<180.)
           converted = kTRUE ;
       }
     }
     if(converted) 
       FillHistogram("hMC_CaloConv_gamma_conv",pt) ;

     //Converted photons with V0
     TLorentzVector pConv ;
     Bool_t foundV0=kFALSE ;
     for(Int_t iv0=0; iv0<fESDEvent->GetNumberOfV0s();iv0++){
       AliESDv0 * v0 = fESDEvent->GetV0(iv0) ;

       TParticle * negativeMC = fStack->Particle(TMath::Abs(fESDEvent->GetTrack(v0->GetNindex())->GetLabel()));
       TParticle * positiveMC = fStack->Particle(TMath::Abs(fESDEvent->GetTrack(v0->GetPindex())->GetLabel()));

       if(negativeMC && positiveMC){
         if(negativeMC->GetMother(0) != positiveMC->GetMother(0))
           continue ;
       }

       if(TMath::Abs(negativeMC->GetPdgCode())!=11 || TMath::Abs(positiveMC->GetPdgCode())!=11){
         continue;
       }
       if(negativeMC->GetPdgCode()==positiveMC->GetPdgCode()){
         continue;
       }

       TParticle * v0Gamma = fStack->Particle(negativeMC->GetMother(0));
       Bool_t same = (v0Gamma == particle) ;
       TParticle * tmp = v0Gamma ;
       while(!same && tmp->GetFirstMother()>=0){
         tmp = fStack->Particle(tmp->GetFirstMother());
         same = (tmp == particle) ;
       }
       if(same){
         foundV0 = kTRUE ;
         const AliExternalTrackParam * paramPos = v0->GetParamP() ;
         const AliExternalTrackParam * paramNeg = v0->GetParamN() ;
         AliKFParticle negKF(*paramNeg,11);
         AliKFParticle posKF(*paramPos,-11);
         pConv.SetXYZM(negKF.Px()+posKF.Px(),negKF.Py()+posKF.Py(),negKF.Pz()+negKF.Pz(),0.) ;
         break ;
       }
     }
     if(foundV0){
       FillHistogram("hMC_CaloConv_gamma_v0",pt) ;
       FillHistogram("hMC_CaloConv_gamma_v0_devsE",(particle->Energy()-pConv.E())/particle->Energy(),particle->Energy()) ;
     }

      //Registered in PHOS/EMCAL
     Bool_t cluInPHOS = kFALSE,cluInEMCAL=kFALSE ;
     TLorentzVector pCalo ;
     Bool_t dist1=kFALSE, dist2=kFALSE ;
     for (Int_t i=0; i<fESDEvent->GetNumberOfCaloClusters(); i++) {
       AliESDCaloCluster * clu = fESDEvent->GetCaloCluster(i);
       Int_t iprim = clu->GetLabel() ; //# of particle hit PHOS/EMCAL
       Bool_t matched = kFALSE ;
       while(iprim>=0 ) {
         if(iprim==iTracks){
           matched=kTRUE ;
           break ;
         }
         else{
           iprim=fStack->Particle(iprim)->GetFirstMother() ;
         }
       }
       if(!matched)
         continue ;
       if(clu->IsPHOS() && hitPHOS){
         cluInPHOS=kTRUE ;
         clu->GetMomentum(pCalo ,vtx);
         if(clu->GetDistanceToBadChannel()<fBadDistCutPHOS)
           dist1=kTRUE ;
         if(clu->GetDistanceToBadChannel()<2.*fBadDistCutPHOS)
           dist2=kTRUE ;
         break ;
       }
       if(!clu->IsPHOS() && hitEMCAL){
         cluInEMCAL=kTRUE ;
         clu->GetMomentum(pCalo ,vtx);
         if(clu->GetDistanceToBadChannel()<fBadDistCutEMCAL)
           dist1=kTRUE ;
         if(clu->GetDistanceToBadChannel()<2.*fBadDistCutEMCAL)
           dist2=kTRUE ;
         break ;
       }
     }

     if(cluInPHOS){
       FillHistogram("hMC_CaloConv_gamma_PHOSclu",pt) ;
       if(!dist1)
         FillHistogram("hMC_CaloConv_gamma_PHOSclu_dist1",pt) ;
       if(!dist2)
         FillHistogram("hMC_CaloConv_gamma_PHOSclu_dist2",pt) ;
       FillHistogram("hMC_CaloConv_gamma_PHOSclu_recE",pCalo.E()) ;
       FillHistogram("hMC_CaloConv_gamma_PHOSclu_devsE",(particle->Energy()-pCalo.E())/particle->Energy(),particle->Energy()) ;
     }
     if(cluInEMCAL){
       FillHistogram("hMC_CaloConv_gamma_EMCALclu",pt) ;
       if(!dist1)
         FillHistogram("hMC_CaloConv_gamma_EMCALclu_dist1",pt) ;
       if(!dist2)
         FillHistogram("hMC_CaloConv_gamma_EMCALclu_dist2",pt) ;
       FillHistogram("hMC_CaloConv_gamma_EMCALclu_recE",pCalo.E()) ;
       FillHistogram("hMC_CaloConv_gamma_EMCALclu_devsE",(particle->Energy()-pCalo.E())/particle->Energy(),particle->Energy()) ;
     }
  }
}
//_____________________________________________________________________________
void AliAnalysisTaskCaloConv::FillHistogram(const char * key,Double_t x)const{
  //FillHistogram
  TH1F * tmp = dynamic_cast<TH1F*>(fOutputContainer->FindObject(key)) ;
  if(!tmp){
    AliInfo(Form("can not find histogram <%s> ",key)) ;
    return ;
  }
  tmp->Fill(x) ;
}
//_____________________________________________________________________________
void AliAnalysisTaskCaloConv::FillHistogram(const char * key,Double_t x,Double_t y)const{
  //FillHistogram
  TObject * tmp = fOutputContainer->FindObject(key) ;
  if(!tmp){
    AliInfo(Form("can not find histogram <%s> ",key)) ;
    return ;
  }
  if(tmp->IsA() == TClass::GetClass("TH1F")){
    ((TH1F*)tmp)->Fill(x,y) ;
    return ;
  }
  if(tmp->IsA() == TClass::GetClass("TH2F")){
    ((TH2F*)tmp)->Fill(x,y) ;
    return ;
  }
  AliError(Form("Calling FillHistogram with 2 parameters for histo <%s> of type %s",key,tmp->IsA()->GetName())) ;
}

//_____________________________________________________________________________
void AliAnalysisTaskCaloConv::FillHistogram(const char * key,Double_t x,Double_t y, Double_t z) const{
  //Fills 1D histograms with key
  TObject * tmp = fOutputContainer->FindObject(key) ;
  if(!tmp){
    AliInfo(Form("can not find histogram <%s> ",key)) ;
    return ;
  }
  if(tmp->IsA() == TClass::GetClass("TH2F")){
    ((TH2F*)tmp)->Fill(x,y,z) ;
    return ;
  }
  if(tmp->IsA() == TClass::GetClass("TH3F")){
    ((TH3F*)tmp)->Fill(x,y,z) ;
    return ;
  }
}
//______________________________________________________________________________
Double_t AliAnalysisTaskCaloConv::PlanarityAngle(const AliExternalTrackParam * paramPos,const AliExternalTrackParam * paramNeg)const {
  //calculate angle between e+e- plain and perpendicular to MF
  //We need sign of MagField to calculate orienation

  TVector3 u(paramPos->Px()+paramNeg->Px(),paramPos->Py()+paramNeg->Py(),paramPos->Pz()+paramNeg->Pz()) ;
  u.Unit() ;
  TVector3 vPos(paramPos->Px(),paramPos->Py(),paramPos->Pz()) ;
  TVector3 vNeg(paramNeg->Px(),paramNeg->Py(),paramNeg->Pz()) ;
  TVector3 v=vPos.Cross(vNeg) ;
  TVector3 w = u.Cross(v);
  TVector3 z(0,0,1.);
  TVector3 ua=u.Cross(z);
  Double_t wa = w.Angle(ua);
  Double_t mfield=fESDEvent->GetMagneticField() ;
  if(mfield>0.)
    return wa;      //normal field
  else
    return TMath::Pi()-wa ; //reverse field

}
//______________________________________________________________________________
Bool_t AliAnalysisTaskCaloConv::IsGoodChannel(const char * det, Int_t mod, Int_t ix, Int_t iz){
//Check if this channel belogs to the good ones

  if(strcmp(det,"PHOS")==0){
    if(mod>5 || mod<1){
      AliError(Form("No bad map for PHOS module %d ",mod)) ;
      return kTRUE ;
    } 
    if(!fPHOSBadMap[mod]){
      AliError(Form("No Bad map for PHOS module %d",mod)) ;
      return kTRUE ;
    }
    if(fPHOSBadMap[mod]->GetBinContent(ix,iz)>0)
      return kFALSE ;
    else
      return kTRUE ;
  }
  else{
    if(strcmp(det,"EMCAL")==0){
      if(mod>9 || mod<0){
        AliError(Form("No bad map for EMCAL module %d ",mod)) ;
        return kTRUE ;
      }
      if(!fEMCALBadMap[mod]){
        AliError(Form("No bad map for EMCAL module %d ",mod)) ;
        return kTRUE ;
      }
      if(fEMCALBadMap[mod]->GetBinContent(ix,iz)>0)
        return kFALSE ;
      else
        return kTRUE ;
    }
    else{
      AliError(Form("Can not find bad channels for detector %s ",det)) ;
    }
  } 
   
  return kTRUE ;
}
//______________________________________________________________________________
void AliAnalysisTaskCaloConv::Recalibrate(Double_t &m, Double_t &pt, const TLorentzVector *calo, const TLorentzVector * conv, Int_t iw, Int_t in) {
  //Apply decalibration and non-linearity
  TLorentzVector calo2(*calo) ;
  Double_t en=calo2.E() ;

  Double_t sigma=0.06+0.005*iw ; //additional smearing
  //Nonlinearity
  Double_t a=0.02*(in%6-2.5) ;
  Double_t b=0.5+1.*((Int_t)in/6) ;
  Double_t enNew=1.-a*TMath::Exp(-en/b) ;
  Double_t corr=gRandom->Gaus(enNew,sigma) ;
  calo2*=corr ;

  m=(calo2+ *conv).M() ;
  pt=(calo2+ *conv).Pt() ;

}
//______________________________________________________________________________
void AliAnalysisTaskCaloConv::RecalibrateEMCAL(Double_t &m, Double_t &pt, const TLorentzVector *calo, const TLorentzVector * conv, Int_t iw, Int_t in) {
  //Apply decalibration and non-linearity
  TLorentzVector calo2(*calo) ;
  Double_t en=calo2.E() ;

  Double_t sigma=0.04+0.005*iw ; //additional smearing
  //Nonlinearity
  Double_t a=0.02*(in%6-2.5) ;
  Double_t b=0.25+0.5*((Int_t)in/6) ;
  Double_t enNew=1.-a*TMath::Exp(-en/b) ;
  Double_t corr=gRandom->Gaus(enNew,sigma) ;
  calo2*=corr ;

  m=(calo2+ *conv).M() ;
  pt=(calo2+ *conv).Pt() ;

}
//______________________________________________________________________________
void AliAnalysisTaskCaloConv::RecalibrateConvPHOS(Double_t &m, Double_t &pt, const TLorentzVector *calo, const TLorentzVector * conv, Int_t iw, Int_t in) {
  //Apply decalibration and non-linearity

  //First default PHOS smearing
  TLorentzVector calo2(*calo) ;
  Double_t en=calo2.E() ;

  Double_t sigma=0.065 ; //additional smearing
  //Nonlinearity
  Double_t a=0.15 ;
  Double_t b=0.45 ;
  Double_t enNew=1.+a*TMath::Exp(-en/b) ;
  Double_t corr=gRandom->Gaus(enNew,sigma) ;
  calo2*=corr ;

  //Now conversion photon
  TLorentzVector conv2(*conv) ;
  //linear offset in z:
  Double_t eta=conv2.Eta() ;
  Double_t c=1.e-3*iw ;
  eta+= c *TMath::Sign(0.9-TMath::Abs(eta),eta) ;

  //Smear energy and add nonlinearity
  //Nonlinearity
  Double_t enConv=conv2.E() ;
  Double_t ac=0.02*(in%5) ;
  Double_t bc=0.25+0.5*((Int_t)in/5) ;
  Double_t enNewc=1.+ac*TMath::Exp(-enConv/bc) ;
  corr=gRandom->Gaus(enNewc,0.01) ;
  Double_t ptc=conv2.Pt()*corr ;
  conv2.SetPtEtaPhiM(ptc,eta,conv2.Phi(),0.) ;

  m =(calo2 + conv2).M() ;
  pt=(calo2 + conv2).Pt() ;

}
//______________________________________________________________________________
void AliAnalysisTaskCaloConv::GetArmenterosQtAlfa(AliKFParticle* positiveKFParticle, AliKFParticle * negativeKFParticle, AliKFParticle * gammaKFCandidate, Double_t armenterosQtAlfa[2] ){
  //see header file for documentation

  TVector3 momentumVectorPositiveKF(positiveKFParticle->GetPx(),positiveKFParticle->GetPy(),positiveKFParticle->GetPz());
  TVector3 momentumVectorNegativeKF(negativeKFParticle->GetPx(),negativeKFParticle->GetPy(),negativeKFParticle->GetPz());
  TVector3 vecV0(gammaKFCandidate->GetPx(),gammaKFCandidate->GetPy(),gammaKFCandidate->GetPz());

  Float_t thetaV0pos=TMath::ACos(( momentumVectorPositiveKF* vecV0)/(momentumVectorPositiveKF.Mag() * vecV0.Mag()));
  Float_t thetaV0neg=TMath::ACos(( momentumVectorNegativeKF* vecV0)/(momentumVectorNegativeKF.Mag() * vecV0.Mag()));

  Float_t alfa =((momentumVectorPositiveKF.Mag())*TMath::Cos(thetaV0pos)-(momentumVectorNegativeKF.Mag())*TMath::Cos(thetaV0neg))/
    ((momentumVectorPositiveKF.Mag())*TMath::Cos(thetaV0pos)+(momentumVectorNegativeKF.Mag())*TMath::Cos(thetaV0neg)) ;


  Float_t qt = momentumVectorPositiveKF.Mag()*TMath::Sin(thetaV0pos);

  armenterosQtAlfa[0]=qt;
  armenterosQtAlfa[1]=alfa;

}




 



