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

// analysis
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisTaskCaloConv.h"
#include "AliStack.h"
#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliESDpid.h"
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
  fStack(NULL),
  fOutputContainer(NULL),
  fCFOutputContainer(NULL),
  fConvCFCont(0x0),
  fPHOSCFCont(0x0),
  fEMCALCFCont(0x0),
  fPi0CFCont(0x0),
  fTriggerCINT1B(kFALSE),
  fToUseCF(kFALSE),
  fMinOpeningAngleGhostCut(0.),
  fPHOSgeom(0x0),
  fEMCALgeom(0x0),
  fPi0Thresh1(0.5),
  fPi0Thresh2(1.),
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
    sprintf(key,"PHOS_BadMap_mod%d",i) ;
    fPHOSBadMap[i]=new TH2I(key,"Bad Modules map",64,0.,64.,56,0.,56.) ;
  }
  for(Int_t i=0; i<10; i++){
    sprintf(key,"EMCAL_BadMap_mod%d",i) ;
    fEMCALBadMap[i] = new TH2I(key,"Bad Modules map",24,0.,24.,48,0.,48.) ;
  }
}
AliAnalysisTaskCaloConv::AliAnalysisTaskCaloConv(const char* name):
  AliAnalysisTaskSE(name),
  fESDEvent(NULL),
  fESDpid(NULL),
  fStack(NULL),
  fOutputContainer(NULL),
  fCFOutputContainer(NULL),
  fConvCFCont(0x0),
  fPHOSCFCont(0x0),
  fEMCALCFCont(0x0),
  fPi0CFCont(0x0),
  fTriggerCINT1B(kFALSE),
  fToUseCF(kFALSE),
  fMinOpeningAngleGhostCut(0.),
  fPHOSgeom(0x0),
  fEMCALgeom(0x0),
  fPi0Thresh1(0.5),
  fPi0Thresh2(1.),
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
    sprintf(key,"PHOS_BadMap_mod%d",i) ;
    fPHOSBadMap[i]=new TH2I(key,"Bad Modules map",64,0.,64.,56,0.,56.) ;
  }
  for(Int_t i=0; i<10; i++){
    sprintf(key,"EMCAL_BadMap_mod%d",i) ;
    fEMCALBadMap[i] = new TH2I(key,"Bad Modules map",24,0.,24.,48,0.,48.) ;
  }
//  fESDpid = new AliESDpid;
}
//_____________________________________________________
AliAnalysisTaskCaloConv::~AliAnalysisTaskCaloConv() 
{
  // Remove all pointers
	
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


  if(!fESDpid){
    AliESDInputHandler *esdHandler=dynamic_cast<AliESDInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
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


 

  fESDEvent=(AliESDEvent*)InputEvent();
  //Take Only events with proper trigger
  //No trigger in MC data => no check
  if(!fStack && !fESDEvent->IsTriggerClassFired("CINT1B-ABCE-NOPF-ALL")){
    return ;
  } 
	
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
  Int_t firstRun= 114700 ;
  Int_t lastRun = 128000 ;
  Int_t nRuns =lastRun-firstRun+1 ;

  //Run QA histigrams
  fOutputContainer->Add(new TH2F("hRunTrigger","Triggers fired",nRuns,float(firstRun),float(lastRun),2,0.,2.)) ;
  fOutputContainer->Add(new TH1F("hRunEvents","Events per run",nRuns,float(firstRun),float(lastRun))) ;
  fOutputContainer->Add(new TH1F("hRunConvs","Conversion photons per run",nRuns,float(firstRun),float(lastRun))) ;
  fOutputContainer->Add(new TH1F("hRunPHOS","PHOS photons per run",nRuns,float(firstRun),float(lastRun))) ;
  fOutputContainer->Add(new TH1F("hRunEMCAL","EMCAL photons per run",nRuns,float(firstRun),float(lastRun))) ;
  fOutputContainer->Add(new TH1F("hVtxBin","Vtx distribution",10,0.,10.)) ;
  fOutputContainer->Add(new TH1F("hEvents","Events processed",1,0.,1.)) ;

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


  //Check of geometry
  fOutputContainer->Add(new TH3F("PHOS_beyond","PHOS clusters not in PHOS",200,-300.,300.,200,-500.,0.,200,-100.,100.)) ;

  fOutputContainer->Add(new TH2F("hdEdx","dEdx of acceptaed electrons",1000,0.,10.,150,0.,150.)) ;

  Int_t npt=200 ;
  Double_t ptmax=20. ;
  //Calibration of PHOS
  fOutputContainer->Add(new TH3F("PHOS_mod1_th1","Inv.Mass distr. per channel",64,0.,64,56,0.,56,200,0.,1.)) ;
  fOutputContainer->Add(new TH3F("PHOS_mod2_th1","Inv.Mass distr. per channel",64,0.,64,56,0.,56,200,0.,1.)) ;
  fOutputContainer->Add(new TH3F("PHOS_mod3_th1","Inv.Mass distr. per channel",64,0.,64,56,0.,56,200,0.,1.)) ;
  fOutputContainer->Add(new TH3F("PHOS_mod4_th1","Inv.Mass distr. per channel",64,0.,64,56,0.,56,200,0.,1.)) ;
  fOutputContainer->Add(new TH3F("PHOS_mod5_th1","Inv.Mass distr. per channel",64,0.,64,56,0.,56,200,0.,1.)) ;

  fOutputContainer->Add(new TH3F("PHOS_mod1_th2","Inv.Mass distr. per channel",64,0.,64,56,0.,56,200,0.,1.)) ;
  fOutputContainer->Add(new TH3F("PHOS_mod2_th2","Inv.Mass distr. per channel",64,0.,64,56,0.,56,200,0.,1.)) ;
  fOutputContainer->Add(new TH3F("PHOS_mod3_th2","Inv.Mass distr. per channel",64,0.,64,56,0.,56,200,0.,1.)) ;
  fOutputContainer->Add(new TH3F("PHOS_mod4_th2","Inv.Mass distr. per channel",64,0.,64,56,0.,56,200,0.,1.)) ;
  fOutputContainer->Add(new TH3F("PHOS_mod5_th2","Inv.Mass distr. per channel",64,0.,64,56,0.,56,200,0.,1.)) ;

  //Pi0 histograms
  //Vary Conversion cuts
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_OnFly","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_Offline","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_On_Kink","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_Off_Kink","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_On_dEdx","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_Off_dEdx","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_On_Prob","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_Off_Prob","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_On_R120","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_Off_R120","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_On_Z","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_Off_Z","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_On_chi","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_Off_chi","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_On_Eta","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_Off_Eta","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_On_Wcut","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_Off_Wcut","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;

  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_OnFly","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_Offline","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_On_Kink","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_Off_Kink","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_On_dEdx","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_Off_dEdx","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_On_Prob","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_Off_Prob","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
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

  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_OnFly","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_Offline","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_On_Kink","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_Off_Kink","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_On_dEdx","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_Off_dEdx","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_On_Prob","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_Off_Prob","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
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
  fOutputContainer->Add(new TH2F("EMCAL_Mi_mvsPt_On_Kink","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Mi_mvsPt_Off_Kink","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Mi_mvsPt_On_dEdx","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Mi_mvsPt_Off_dEdx","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Mi_mvsPt_On_Prob","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Mi_mvsPt_Off_Prob","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
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
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_Disp","Mass vs pt, disp cut",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_TOF","Mass vs pt, TOF cut",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_Neutral","Mass vs pt, Neutral cut",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_DispNeutral","Mass vs pt, Disp & neutral cut",400,0.,1.,npt,0.,ptmax)) ;

  //PHOS PID module-by-module
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_mod1_all","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_mod2_all","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_mod3_all","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_mod4_all","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_mod5_all","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_mod1_Disp","Mass vs pt, disp cut",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_mod2_Disp","Mass vs pt, disp cut",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_mod3_Disp","Mass vs pt, disp cut",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_mod4_Disp","Mass vs pt, disp cut",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_mod5_Disp","Mass vs pt, disp cut",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_mod1_TOF","Mass vs pt, TOF cut",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_mod2_TOF","Mass vs pt, TOF cut",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_mod3_TOF","Mass vs pt, TOF cut",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_mod4_TOF","Mass vs pt, TOF cut",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_mod5_TOF","Mass vs pt, TOF cut",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_mod1_Neutral","Mass vs pt, Neutral cut",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_mod2_Neutral","Mass vs pt, Neutral cut",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_mod3_Neutral","Mass vs pt, Neutral cut",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_mod4_Neutral","Mass vs pt, Neutral cut",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_mod5_Neutral","Mass vs pt, Neutral cut",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_mod1_DispNeutral","Mass vs pt, Disp & neutral cut",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_mod2_DispNeutral","Mass vs pt, Disp & neutral cut",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_mod3_DispNeutral","Mass vs pt, Disp & neutral cut",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_mod4_DispNeutral","Mass vs pt, Disp & neutral cut",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_mod5_DispNeutral","Mass vs pt, Disp & neutral cut",400,0.,1.,npt,0.,ptmax)) ;

  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_all","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_Disp","Mass vs pt, disp cut",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_TOF","Mass vs pt, TOF cut",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_Neutral","Mass vs pt, Neutral cut",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_DispNeutral","Mass vs pt, Disp & neutral cut",400,0.,1.,npt,0.,ptmax)) ;

  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_mod1_single","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_mod2_single","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_mod3_single","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_mod4_single","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Re_mvsPt_mod5_single","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;

  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_mod1_all","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_mod2_all","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_mod3_all","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_mod4_all","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_mod5_all","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_mod1_Disp","Mass vs pt, disp cut",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_mod2_Disp","Mass vs pt, disp cut",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_mod3_Disp","Mass vs pt, disp cut",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_mod4_Disp","Mass vs pt, disp cut",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_mod5_Disp","Mass vs pt, disp cut",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_mod1_TOF","Mass vs pt, TOF cut",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_mod2_TOF","Mass vs pt, TOF cut",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_mod3_TOF","Mass vs pt, TOF cut",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_mod4_TOF","Mass vs pt, TOF cut",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_mod5_TOF","Mass vs pt, TOF cut",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_mod1_Neutral","Mass vs pt, Neutral cut",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_mod2_Neutral","Mass vs pt, Neutral cut",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_mod3_Neutral","Mass vs pt, Neutral cut",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_mod4_Neutral","Mass vs pt, Neutral cut",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_mod5_Neutral","Mass vs pt, Neutral cut",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_mod1_DispNeutral","Mass vs pt, Disp & neutral cut",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_mod2_DispNeutral","Mass vs pt, Disp & neutral cut",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_mod3_DispNeutral","Mass vs pt, Disp & neutral cut",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_mod4_DispNeutral","Mass vs pt, Disp & neutral cut",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("PHOS_Mi_mvsPt_mod5_DispNeutral","Mass vs pt, Disp & neutral cut",400,0.,1.,npt,0.,ptmax)) ;

  //Single photon spectrum
  //Conversion
  fOutputContainer->Add(new TH1F("Single_conv","Single photon spectrum",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("Single_conv_OnFly","Single photon spectrum",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("Single_conv_Offline","Single photon spectrum",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("Single_conv_On_Kink","Single photon spectrum",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("Single_conv_Off_Kink","Single photon spectrum",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("Single_conv_On_dEdx","Single photon spectrum",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("Single_conv_Off_dEdx","Single photon spectrum",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("Single_conv_On_Prob","Single photon spectrum",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("Single_conv_Off_Prob","Single photon spectrum",npt,0.,ptmax)) ;
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
  fOutputContainer->Add(new TH1F("PHOS_single_mod1_all","Single photon spectrum",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("PHOS_single_mod2_all","Single photon spectrum",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("PHOS_single_mod3_all","Single photon spectrum",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("PHOS_single_mod4_all","Single photon spectrum",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("PHOS_single_mod5_all","Single photon spectrum",npt,0.,ptmax)) ;

  fOutputContainer->Add(new TH1F("PHOS_single_mod1_disp","Single photon spectrum",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("PHOS_single_mod2_disp","Single photon spectrum",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("PHOS_single_mod3_disp","Single photon spectrum",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("PHOS_single_mod4_disp","Single photon spectrum",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("PHOS_single_mod5_disp","Single photon spectrum",npt,0.,ptmax)) ;

  fOutputContainer->Add(new TH1F("PHOS_single_mod1_neutral","Single photon spectrum",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("PHOS_single_mod2_neutral","Single photon spectrum",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("PHOS_single_mod3_neutral","Single photon spectrum",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("PHOS_single_mod4_neutral","Single photon spectrum",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("PHOS_single_mod5_neutral","Single photon spectrum",npt,0.,ptmax)) ;

  fOutputContainer->Add(new TH1F("PHOS_single_mod1_dispneutral","Single photon spectrum",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("PHOS_single_mod2_dispneutral","Single photon spectrum",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("PHOS_single_mod3_dispneutral","Single photon spectrum",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("PHOS_single_mod4_dispneutral","Single photon spectrum",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("PHOS_single_mod5_dispneutral","Single photon spectrum",npt,0.,ptmax)) ;

  fOutputContainer->Add(new TH1F("PHOS_single_mod1_dist1","Single photon spectrum",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("PHOS_single_mod2_dist1","Single photon spectrum",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("PHOS_single_mod3_dist1","Single photon spectrum",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("PHOS_single_mod4_dist1","Single photon spectrum",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("PHOS_single_mod5_dist1","Single photon spectrum",npt,0.,ptmax)) ;

  fOutputContainer->Add(new TH1F("PHOS_single_mod1_dist2","Single photon spectrum",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("PHOS_single_mod2_dist2","Single photon spectrum",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("PHOS_single_mod3_dist2","Single photon spectrum",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("PHOS_single_mod4_dist2","Single photon spectrum",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("PHOS_single_mod5_dist2","Single photon spectrum",npt,0.,ptmax)) ;

  fOutputContainer->Add(new TH3F("EMCAL_Ep_0","PHOS E/p ratio",24,0.,24,48,0.,48,200,0.,2.)) ;
  fOutputContainer->Add(new TH3F("EMCAL_Ep_1","PHOS E/p ratio",24,0.,24,48,0.,48,200,0.,2.)) ;
  fOutputContainer->Add(new TH3F("EMCAL_Ep_2","PHOS E/p ratio",24,0.,24,48,0.,48,200,0.,2.)) ;
  fOutputContainer->Add(new TH3F("EMCAL_Ep_3","PHOS E/p ratio",24,0.,24,48,0.,48,200,0.,2.)) ;
  fOutputContainer->Add(new TH3F("EMCAL_EpPE","PHOS E/p vs E vs p",200,0.,2.,npt,0.,ptmax,npt,0.,ptmax)) ;

  fOutputContainer->Add(new TH3F("EMCAL_beyond","EMCAL clusters not in EMCAL",200,-300.,300.,200,-500.,0.,200,-100.,100.)) ;

  fOutputContainer->Add(new TH3F("EMCAL_mod0_th1","Inv.Mass distr. per channel",24,0.,24,48,0.,48,200,0.,1.)) ;
  fOutputContainer->Add(new TH3F("EMCAL_mod1_th1","Inv.Mass distr. per channel",24,0.,24,48,0.,48,200,0.,1.)) ;
  fOutputContainer->Add(new TH3F("EMCAL_mod2_th1","Inv.Mass distr. per channel",24,0.,24,48,0.,48,200,0.,1.)) ;
  fOutputContainer->Add(new TH3F("EMCAL_mod3_th1","Inv.Mass distr. per channel",24,0.,24,48,0.,48,200,0.,1.)) ;

  fOutputContainer->Add(new TH3F("EMCAL_mod0_th2","Inv.Mass distr. per channel",24,0.,24,48,0.,48,200,0.,1.)) ;
  fOutputContainer->Add(new TH3F("EMCAL_mod1_th2","Inv.Mass distr. per channel",24,0.,24,48,0.,48,200,0.,1.)) ;
  fOutputContainer->Add(new TH3F("EMCAL_mod2_th2","Inv.Mass distr. per channel",24,0.,24,48,0.,48,200,0.,1.)) ;
  fOutputContainer->Add(new TH3F("EMCAL_mod3_th2","Inv.Mass distr. per channel",24,0.,24,48,0.,48,200,0.,1.)) ;

  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_mod0_single","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_mod1_single","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_mod2_single","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_mod3_single","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;

  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_mod0_all","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_mod1_all","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_mod2_all","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_mod3_all","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_mod4_all","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_mod5_all","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_mod0_Disp","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_mod1_Disp","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_mod2_Disp","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_mod3_Disp","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_mod4_Disp","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_mod5_Disp","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_mod0_TOF","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_mod1_TOF","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_mod2_TOF","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_mod3_TOF","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_mod4_TOF","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_mod5_TOF","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_mod0_Neutral","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_mod1_Neutral","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_mod2_Neutral","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_mod3_Neutral","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_mod4_Neutral","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_mod5_Neutral","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_mod0_DispNeutral","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_mod1_DispNeutral","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_mod2_DispNeutral","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_mod3_DispNeutral","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_mod4_DispNeutral","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Re_mvsPt_mod5_DispNeutral","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;

  fOutputContainer->Add(new TH2F("EMCAL_Mi_mvsPt_mod0_all","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Mi_mvsPt_mod1_all","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Mi_mvsPt_mod2_all","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Mi_mvsPt_mod3_all","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Mi_mvsPt_mod4_all","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Mi_mvsPt_mod5_all","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Mi_mvsPt_mod0_Disp","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Mi_mvsPt_mod1_Disp","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Mi_mvsPt_mod2_Disp","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Mi_mvsPt_mod3_Disp","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Mi_mvsPt_mod4_Disp","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Mi_mvsPt_mod5_Disp","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Mi_mvsPt_mod0_TOF","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Mi_mvsPt_mod1_TOF","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Mi_mvsPt_mod2_TOF","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Mi_mvsPt_mod3_TOF","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Mi_mvsPt_mod4_TOF","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Mi_mvsPt_mod5_TOF","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Mi_mvsPt_mod0_Neutral","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Mi_mvsPt_mod1_Neutral","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Mi_mvsPt_mod2_Neutral","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Mi_mvsPt_mod3_Neutral","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Mi_mvsPt_mod4_Neutral","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Mi_mvsPt_mod5_Neutral","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Mi_mvsPt_mod0_DispNeutral","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Mi_mvsPt_mod1_DispNeutral","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Mi_mvsPt_mod2_DispNeutral","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Mi_mvsPt_mod3_DispNeutral","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Mi_mvsPt_mod4_DispNeutral","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("EMCAL_Mi_mvsPt_mod5_DispNeutral","Mass vs pt",400,0.,1.,npt,0.,ptmax)) ;

  //MC info
  fOutputContainer->Add(new TH1F("hMC_CaloConv_pi0_unitEta","Primary #pi^{0}",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_allPi0","Primary #pi^{0}",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_pi0PHOSacc","#pi^{0} decayed in PHOS acc",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_pi0EMCALacc","#pi^{0} decayed in EMCAL acc",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_pi0_PHOS_conv","#pi^{0} decayed in PHOS acc asnd conv. photon",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_pi0_EMCAL_conv","#pi^{0} decayed in EMCAL acc asnd conv. photon",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_pi0_bothphot_conv","#pi^{0} both photons converted",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_pi0__convPhotInCalo","#pi^{0} photon in calo converted",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_pi0_v0_PHOSacc","#pi^{0} photon converted and V0 found",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_pi0_v0_EMCALacc","#pi^{0} photon converted and V0 found",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_pi0_v0_PHOSclu","#pi^{0} V0 and cluster in PHOS found",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_pi0_v0_EMCALclu_ptRec","#pi^{0} V0 and cluster in EMCAL found (rec pt)",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_pi0_v0_PHOSclu_ptRec","#pi^{0} V0 and cluster in PHOS found(rec pt)",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_pi0_v0_EMCALclu","#pi^{0} V0 and cluster in EMCAL found",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("hMC_CaloConv_pi0_v0_PHOSclu_mvsPt","m vs pt for rec pi0s",400,0.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("hMC_CaloConv_pi0_v0_EMCALclu_mvsPt","m vs pt for rec pi0s",400,0.,1.,npt,0.,ptmax)) ;


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
  fOutputContainer->Add(new TH1F("hMC_CaloConv_gamma_EMCALclu","Photons with cluster in EMCAL",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_gamma_PHOSclu_recE","Photons with cluster in PHOS",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH1F("hMC_CaloConv_gamma_EMCALclu_recE","Photons with cluster in EMCAL",npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("hMC_CaloConv_gamma_PHOSclu_devsE","Photons with cluster in PHOS",200,-1.,1.,npt,0.,ptmax)) ;
  fOutputContainer->Add(new TH2F("hMC_CaloConv_gamma_EMCALclu_devsE","Photons with cluster in EMCAL",200,-1.,1.,npt,0.,ptmax)) ;
	

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
    if(p.Energy()<0.25)
      continue ;
    if(clu->GetNCells()<=2)
      continue ;

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
      FillHistogram("PHOS_beyond",xyz[0],xyz[1],xyz[2]) ;
      printf("PHOS_beyond: x=%f, y=%f, z=%f \n",xyz[0],xyz[1],xyz[2]) ;
      continue ;
    }
    iMod=relid[0] ;
    iX=relid[2];
    iZ=relid[3] ;
    if(!IsGoodChannel("PHOS",iMod,iX,iZ))
      continue ;

    p.SetBit(kCaloPIDdisp,isDispOK) ;
    p.SetBit(kCaloPIDtof,isTOFOK) ;
    p.SetBit(kCaloPIDneutral,isNeutral) ;
    p.SetBit(BIT(16+iMod),kTRUE) ;
    new((*fPHOSEvent)[inPHOS]) TLorentzVector(p) ;
    fGammaPHOS[inPHOS] = i ;
    inPHOS++ ;


    //Single photon spectra
    Double_t pt= p.Pt() ;
    TString skey="PHOS_single_"; skey+="mod" ; skey+=iMod ; skey+="_all" ;
    FillHistogram(skey,pt) ;
    if(isDispOK){
      skey="PHOS_single_"; skey+="mod" ; skey+=iMod ; skey+="_disp" ;
      FillHistogram(skey,pt) ;
    }
    if(isNeutral){
      skey="PHOS_single_"; skey+="mod" ; skey+=iMod ; skey+="_neutral" ;
      FillHistogram(skey,pt) ;
    }
    if(isNeutral && isDispOK){
      skey="PHOS_single_"; skey+="mod" ; skey+=iMod ; skey+="_dispneutral" ;
      FillHistogram(skey,pt) ;
    }
    //Distance to bad channel
    if(clu->GetDistanceToBadChannel()>2.2){
      skey="PHOS_single_"; skey+="mod" ; skey+=iMod ; skey+="_dist1" ;
      FillHistogram(skey,pt) ;
    }
    if(clu->GetDistanceToBadChannel()>4.4){
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

  // Loop over all CaloClusters
  if(fEMCALEvent)
    fEMCALEvent->Clear() ;
  else
    fEMCALEvent = new TClonesArray("TLorentzVector",10) ;
  Int_t inEMCAL = 0 ; //, inEMCALRecal=0;
  TLorentzVector pi0 ;

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
      fGammaPHOS[inEMCAL] = i ;
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
//printf("st 1: px=%f, py=%f, pz=%f, E=%f \n",photKF.GetPx(),photKF.GetPy(),photKF.GetPz(),photKF.GetE()) ;
    photKF.SetMassConstraint(0,cutSigmaMass);
//printf("st 2: px=%f, py=%f, pz=%f, E=%f \n",photKF.GetPx(),photKF.GetPy(),photKF.GetPz(),photKF.GetE()) ;

    if(useImprovedVertex){
      AliKFVertex primaryVertexImproved(*(fESDEvent->GetPrimaryVertex()));
      primaryVertexImproved+=photKF;
      photKF.SetProductionVertex(primaryVertexImproved);
    }
//printf("st 3: px=%f, py=%f, pz=%f, E=%f \n",photKF.GetPx(),photKF.GetPy(),photKF.GetPz(),photKF.GetE()) ;
    Double_t m, width ;
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
//printf("... likesign \n") ;
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
//printf("... status \n") ;
      continue;
    }
    if(fToUseCF){
      if(isOnFly)
        fConvCFCont->Fill(a,5) ;
      else
        fConvCFCont->Fill(a,6) ;
    }
 
    Bool_t isKink=kFALSE ;
    if( neg->GetKinkIndex(0) > 0 ||
        pos->GetKinkIndex(0) > 0) {
      isKink=kTRUE;
    }
    if(!isKink && fToUseCF){
      if(isOnFly)
        fConvCFCont->Fill(a,7) ;
      else
        fConvCFCont->Fill(a,8) ;
    }

    //First rough PID
    if( fESDpid->NumberOfSigmasTPC(pos,AliPID::kElectron)<-4. ||
        fESDpid->NumberOfSigmasTPC(pos,AliPID::kElectron)>6. ||
        fESDpid->NumberOfSigmasTPC(neg,AliPID::kElectron)<-4. ||
        fESDpid->NumberOfSigmasTPC(neg,AliPID::kElectron)>6. ){
        continue ;
    }

    Bool_t isdEdx=kTRUE;
    if( fESDpid->NumberOfSigmasTPC(pos,AliPID::kElectron)<fnSigmaBelowElectronLine ||
        fESDpid->NumberOfSigmasTPC(pos,AliPID::kElectron)>fnSigmaAboveElectronLine ||
        fESDpid->NumberOfSigmasTPC(neg,AliPID::kElectron)<fnSigmaBelowElectronLine ||
        fESDpid->NumberOfSigmasTPC(neg,AliPID::kElectron)>fnSigmaAboveElectronLine ){
//printf("... dEdx 1 \n") ;
         isdEdx=kFALSE;
    }
    const Double_t minPnSigmaAbovePionLine = 0.5 ;
    const Double_t maxPnSigmaAbovePionLine = 100. ;
    const Double_t nSigmaAbovePionLine = 0 ;
    if(pos->P()>minPnSigmaAbovePionLine && pos->P()<maxPnSigmaAbovePionLine ){
      if(fESDpid->NumberOfSigmasTPC(pos,AliPID::kPion)<nSigmaAbovePionLine){
//printf("... dEdx 2 \n") ;
          isdEdx=kFALSE;
        }
    }
    if(neg->P()>minPnSigmaAbovePionLine && neg->P()<maxPnSigmaAbovePionLine){
      if(fESDpid->NumberOfSigmasTPC(neg,AliPID::kPion)<nSigmaAbovePionLine){
//printf("... dEdx 3 \n") ;
          isdEdx=kFALSE;
      }
    }

    //Kaon rejection
    const Double_t minPKaonRejection=1.5 ;
    const Double_t sigmaAroundLine=1. ;
    if(neg->P()<minPKaonRejection ){
      if(TMath::Abs(fESDpid->NumberOfSigmasTPC(neg,AliPID::kKaon))<sigmaAroundLine){
//printf("... dEdx 4 \n") ;
        isdEdx=kFALSE;
      }
    }
    if(pos->P()<minPKaonRejection ){
      if(TMath::Abs(fESDpid->NumberOfSigmasTPC(pos,AliPID::kKaon))<sigmaAroundLine){
//printf("... dEdx 5 \n") ;
        isdEdx=kFALSE;
      }
    }

    //Proton rejection
    const Double_t minPProtonRejection=2. ;
    if(neg->P()<minPProtonRejection){
      if(TMath::Abs(fESDpid->NumberOfSigmasTPC(neg,AliPID::kProton))<sigmaAroundLine){
//printf("... dEdx 6 \n") ;
        isdEdx=kFALSE;
      }
    }
    if(pos->P()<minPProtonRejection ){
      if(TMath::Abs(fESDpid->NumberOfSigmasTPC(pos,AliPID::kProton))<sigmaAroundLine){
//printf("... dEdx 7 \n") ;
        isdEdx=kFALSE;
      }
    }

    const Double_t minPPionRejection=0.5 ;
    if(neg->P()<minPPionRejection ){
      if(TMath::Abs(fESDpid->NumberOfSigmasTPC(neg,AliPID::kPion))<sigmaAroundLine){
//printf("... dEdx 8 \n") ;
        isdEdx=kFALSE;
      }
    }
    if(pos->P()<minPPionRejection ){
      if( TMath::Abs(fESDpid->NumberOfSigmasTPC(pos,AliPID::kPion))<sigmaAroundLine){
//printf("... dEdx 9 \n") ;
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
    if(!isKink && isProb && fToUseCF){
      if(isOnFly)
        fConvCFCont->Fill(a,9) ;
      else
        fConvCFCont->Fill(a,10) ;
    }

    Double_t v0x,v0y,v0z;
    v0->GetXYZ(v0x,v0y,v0z) ;
    Double_t r=TMath::Sqrt(v0x*v0x + v0y*v0y) ;
    if(r>fmaxR){ // cuts on distance from collision point
//printf("... maxR \n") ;
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
//printf("... ZR slope=%f, offset=%f, z=%f, zs=%f, r=%f \n",zrSlope,zOffset,v0z,TMath::Abs(v0z)*zrSlope-zOffset,r) ;
      continue;
    }
    if(fToUseCF){
      if(isOnFly)
        fConvCFCont->Fill(a,13) ;
      else
        fConvCFCont->Fill(a,14) ;
    }

    if(TMath::Abs(v0z) > fmaxZ ){ // cuts out regions where we do not reconstruct
//printf("... maxZ \n") ;
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
//printf("... NDF \n") ;
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
//printf("... chi2 \n") ;
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
//printf("... ETA \n") ;
      continue;
    }
    if(TMath::Abs(paramPos->Eta())> wideEtaCut ||  
       TMath::Abs(paramNeg->Eta())> wideEtaCut ){
//printf("... ETA pls mns \n") ;
      continue ;
    }

    Bool_t isWideEta=kTRUE ;
    if(TMath::Abs(photLV.Eta())< fetaCut && 
       TMath::Abs(paramPos->Eta())<fetaCut  && 
       TMath::Abs(paramNeg->Eta()) < fetaCut){
      isWideEta=kFALSE;
    }
    

    if(photLV.Pt()<fptCut){
//printf("... pt \n") ;
      continue;
    }
    if(fToUseCF){
      if(isOnFly)
        fConvCFCont->Fill(a,21) ;
      else
        fConvCFCont->Fill(a,22) ;
    }


    //Just QA plot
/*
    if(photLV.Pt()>0.5){
       Double_t phi=photLV.Phi() ;
       while(phi<0.)phi+=TMath::TwoPi() ;
       while(phi>TMath::TwoPi())phi-=TMath::TwoPi() ;
       FillHistogram("ConvPhiEta",phi,photLV.Eta()) ; 
    }
*/
    
    Double_t w=PlanarityAngle(paramPos,paramNeg) ;
    Bool_t isPlanarityCut = (0.08-0.22*w > m || 0.15*(w-2.4)>m) ;
    FillHistogram("All_w_vs_m",w,m) ;

    photLV.SetBit(kConvOnFly,isOnFly) ;
    photLV.SetBit(kConvKink,isKink) ;
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

//    if(isOnFly)
//      printf("CaloConv: v0(%d): onFly \n",inConv) ;
//    else
//      printf("CaloConv: v0(%d): Offline \n",inConv) ;

    //Single photon spectrum
    Double_t pt=photLV.Pt() ;
    if(isOnFly){
      //Default
      if(!isKink && isdEdx && isProb && !isWideEta)
        FillHistogram("Single_conv_OnFly",pt) ;
      if(isdEdx && isProb && !isWideEta)
        FillHistogram("Single_conv_On_Kink",pt) ;
      if(!isKink && isProb && !isWideEta)
        FillHistogram("Single_conv_On_dEdx",pt) ;
      if(!isKink && isdEdx && !isWideEta)
        FillHistogram("Single_conv_On_Prob",pt) ;
      if(!isKink && isdEdx && isProb && !isWideEta && isStrictR)
        FillHistogram("Single_conv_On_R120",pt) ; 
      if(!isKink && isdEdx && isProb && !isWideEta && isStrictZ)
        FillHistogram("Single_conv_On_Z",pt) ;
      if(!isKink && isdEdx && isProb && !isWideEta && isStrictChi)
        FillHistogram("Single_conv_On_chi",pt) ;
      if(!isKink && isdEdx && isProb)
        FillHistogram("Single_conv_On_Eta",pt) ;
      if(!isKink && isdEdx && isProb && !isWideEta && isPlanarityCut)
        FillHistogram("Single_conv_On_Wcut",pt) ;
    }
    else{
      if(!isKink && isdEdx && isProb && !isWideEta)
        FillHistogram("Single_conv_Offline",pt) ;
      if(isdEdx && isProb && !isWideEta)
        FillHistogram("Single_conv_Off_Kink",pt) ;
      if(!isKink && isProb && !isWideEta)
        FillHistogram("Single_conv_Off_dEdx",pt) ;
      if(!isKink && isdEdx && !isWideEta)
        FillHistogram("Single_conv_Off_Prob",pt) ;
      if(!isKink && isdEdx && isProb && !isWideEta && isStrictR)
        FillHistogram("Single_conv_Off_R120",pt) ; 
      if(!isKink && isdEdx && isProb && !isWideEta && isStrictZ)
        FillHistogram("Single_conv_Off_Z",pt) ;
      if(!isKink && isdEdx && isProb && !isWideEta && isStrictChi)
        FillHistogram("Single_conv_Off_chi",pt) ;
      if(!isKink && isdEdx && isProb)
        FillHistogram("Single_conv_Off_Eta",pt) ;
      if(!isKink && isdEdx && isProb && !isWideEta && isPlanarityCut)
        FillHistogram("Single_conv_Off_Wcut",pt) ;
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
    if(cnv->TestBit(kConvOnFly) && !cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta)){
      nConvGood++ ;
    }
  }
  FillHistogram("hRunConvs",run,double(nConvGood)) ;
  FillHistogram("hRunPHOS", run,double(nPHOS)) ;
  FillHistogram("hRunEMCAL",run,double(nEMCAL)) ;
//printf("CaloConv::FillRe:  nConv=%d \n",nConvGood) ;

  //Fill Real distributions
  for(Int_t iPHOS=0; iPHOS<nPHOS;iPHOS++){
    TLorentzVector * cal = static_cast<TLorentzVector*>(fPHOSEvent->At(iPHOS)) ;
    for(Int_t iConv = 0; iConv<nConv; iConv++){
      TLorentzVector * cnv = static_cast<TLorentzVector*>(fConvEvent->At(iConv)) ;
      TLorentzVector pi=*cal + *cnv ;
      Double_t alpha=TMath::Abs(cal->Energy()-cnv->Energy())/(cal->Energy()+cnv->Energy()) ;
      if(cnv->TestBit(kConvOnFly) && !cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta)){
        //Non-linearity check
        FillHistogram("PHOS_Re_mvsPt_all",pi.M(),pi.Pt()) ;
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
        if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta))
          FillHistogram("PHOS_Re_mvsPt_OnFly",pi.M(),pi.Pt()) ;
        if(cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta))
          FillHistogram("PHOS_Re_mvsPt_On_Kink",pi.M(),pi.Pt()) ;
        if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta))
          FillHistogram("PHOS_Re_mvsPt_On_dEdx",pi.M(),pi.Pt()) ;
        if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && !cnv->TestBit(kConvEta))
          FillHistogram("PHOS_Re_mvsPt_On_Prob",pi.M(),pi.Pt()) ;
        if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta) && cnv->TestBit(kConvR))
          FillHistogram("PHOS_Re_mvsPt_On_R120",pi.M(),pi.Pt()) ;
        if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta) && cnv->TestBit(kConvZR))
          FillHistogram("PHOS_Re_mvsPt_On_Z",pi.M(),pi.Pt()) ;
        if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta) && cnv->TestBit(kConvNDF)) 
          FillHistogram("PHOS_Re_mvsPt_On_chi",pi.M(),pi.Pt()) ;
        if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb))
          FillHistogram("PHOS_Re_mvsPt_On_Eta",pi.M(),pi.Pt()) ;
        if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta) && cnv->TestBit(kConvPlan)) 
          FillHistogram("PHOS_Re_mvsPt_On_Wcut",pi.M(),pi.Pt()) ;
      }
      else{
        //Default
        if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta))
          FillHistogram("PHOS_Re_mvsPt_Offline",pi.M(),pi.Pt()) ;
        if(cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta))
          FillHistogram("PHOS_Re_mvsPt_Off_Kink",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta))
          FillHistogram("PHOS_Re_mvsPt_Off_dEdx",pi.M(),pi.Pt()) ;
        if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && !cnv->TestBit(kConvEta))
          FillHistogram("PHOS_Re_mvsPt_Off_Prob",pi.M(),pi.Pt()) ;
        if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta) && cnv->TestBit(kConvR))
          FillHistogram("PHOS_Re_mvsPt_Off_R120",pi.M(),pi.Pt()) ;
        if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta) && cnv->TestBit(kConvZR))
          FillHistogram("PHOS_Re_mvsPt_Off_Z",pi.M(),pi.Pt()) ;
        if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta) && cnv->TestBit(kConvNDF))
          FillHistogram("PHOS_Re_mvsPt_Off_chi",pi.M(),pi.Pt()) ;
        if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb))
          FillHistogram("PHOS_Re_mvsPt_Off_Eta",pi.M(),pi.Pt()) ;
        if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta) && cnv->TestBit(kConvPlan))
          FillHistogram("PHOS_Re_mvsPt_Off_Wcut",pi.M(),pi.Pt()) ;
      }
    }
  }
  //PHOS module-dependent histograms
  for(Int_t iPHOS=0; iPHOS<nPHOS;iPHOS++){
    TLorentzVector * cal = static_cast<TLorentzVector*>(fPHOSEvent->At(iPHOS)) ;
    Int_t mod=1;
    while(!cal->TestBit(BIT(16+mod)) && mod<5)mod++ ;
    TString base("PHOS_Re_mvsPt_mod") ; base+=mod ;
    TString full ;
    for(Int_t iConv = 0; iConv<nConv; iConv++){
      TLorentzVector * cnv = static_cast<TLorentzVector*>(fConvEvent->At(iConv)) ;
      TLorentzVector pi=*cal + *cnv ;
      full=base ; full+="_single" ;
      if(cnv->TestBit(kConvOnFly) && !cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta)){
        FillHistogram(full,pi.M(),cal->Pt()) ;
        full=base ; full+="_all" ;
        FillHistogram(full,pi.M(),pi.Pt()) ;
        if(cal->TestBit(kCaloPIDdisp)){
          full=base ; full+="_Disp" ;
          FillHistogram(full,pi.M(),pi.Pt()) ;
        }
        if(cal->TestBit(kCaloPIDtof)){
          full=base ; full+="_TOF" ;
          FillHistogram(full,pi.M(),pi.Pt()) ;
        }
          if(cal->TestBit(kCaloPIDneutral)){
          full=base ; full+="_Neutral" ;
          FillHistogram(full,pi.M(),pi.Pt()) ;
        }
        if(cal->TestBit(kCaloPIDneutral) && cal->TestBit(kCaloPIDdisp)){
          full=base ; full+="_DispNeutral" ;
            FillHistogram(full,pi.M(),pi.Pt()) ;
        }
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
      if(cnv->TestBit(kConvOnFly) && !cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta)){
        FillHistogram(full,pi.M(),cal->Pt()) ;
        full=base+"_all" ;
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
      //Vary Conversion cuts
      if(cnv->TestBit(kConvOnFly)){
        //Default
        if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta))
          FillHistogram("EMCAL_Re_mvsPt_OnFly",pi.M(),pi.Pt()) ;
        if(cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta))
          FillHistogram("EMCAL_Re_mvsPt_On_Kink",pi.M(),pi.Pt()) ;
        if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta))
          FillHistogram("EMCAL_Re_mvsPt_On_dEdx",pi.M(),pi.Pt()) ;
        if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && !cnv->TestBit(kConvEta))
          FillHistogram("EMCAL_Re_mvsPt_On_Prob",pi.M(),pi.Pt()) ;
        if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta) && cnv->TestBit(kConvR))
          FillHistogram("EMCAL_Re_mvsPt_On_R120",pi.M(),pi.Pt()) ;
        if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta) && cnv->TestBit(kConvZR))
          FillHistogram("EMCAL_Re_mvsPt_On_Z",pi.M(),pi.Pt()) ;
        if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta) && cnv->TestBit(kConvNDF))
          FillHistogram("EMCAL_Re_mvsPt_On_chi",pi.M(),pi.Pt()) ;
        if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb))
          FillHistogram("EMCAL_Re_mvsPt_On_Eta",pi.M(),pi.Pt()) ;
        if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta) && cnv->TestBit(kConvPlan))
          FillHistogram("EMCAL_Re_mvsPt_On_Wcut",pi.M(),pi.Pt()) ;
      }
      else{
        //Default
        if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta))
          FillHistogram("EMCAL_Re_mvsPt_Offline",pi.M(),pi.Pt()) ;
        if(cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta))
          FillHistogram("EMCAL_Re_mvsPt_Off_Kink",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta))
          FillHistogram("EMCAL_Re_mvsPt_Off_dEdx",pi.M(),pi.Pt()) ;
        if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && !cnv->TestBit(kConvEta))
          FillHistogram("EMCAL_Re_mvsPt_Off_Prob",pi.M(),pi.Pt()) ;
        if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta) && cnv->TestBit(kConvR))
          FillHistogram("EMCAL_Re_mvsPt_Off_R120",pi.M(),pi.Pt()) ;
        if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta) && cnv->TestBit(kConvZR))
          FillHistogram("EMCAL_Re_mvsPt_Off_Z",pi.M(),pi.Pt()) ;
        if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta) && cnv->TestBit(kConvNDF))
          FillHistogram("EMCAL_Re_mvsPt_Off_chi",pi.M(),pi.Pt()) ;
        if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb))
          FillHistogram("EMCAL_Re_mvsPt_Off_Eta",pi.M(),pi.Pt()) ;
        if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta) && cnv->TestBit(kConvPlan))
          FillHistogram("EMCAL_Re_mvsPt_Off_Wcut",pi.M(),pi.Pt()) ;
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
        if(cnv->TestBit(kConvOnFly) && !cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta)){
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
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta))
            FillHistogram("PHOS_Mi_mvsPt_OnFly",pi.M(),pi.Pt()) ;
          if(cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta))
            FillHistogram("PHOS_Mi_mvsPt_On_Kink",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta))
            FillHistogram("PHOS_Mi_mvsPt_On_dEdx",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && !cnv->TestBit(kConvEta))
            FillHistogram("PHOS_Mi_mvsPt_On_Prob",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta) && cnv->TestBit(kConvR))
            FillHistogram("PHOS_Mi_mvsPt_On_R120",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta) && cnv->TestBit(kConvZR))
            FillHistogram("PHOS_Mi_mvsPt_On_Z",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta) && cnv->TestBit(kConvNDF))
            FillHistogram("PHOS_Mi_mvsPt_On_chi",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb))
            FillHistogram("PHOS_Mi_mvsPt_On_Eta",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta) && cnv->TestBit(kConvPlan))
            FillHistogram("PHOS_Mi_mvsPt_On_Wcut",pi.M(),pi.Pt()) ;
        }
        else{
          //Default
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta))
            FillHistogram("PHOS_Mi_mvsPt_Offline",pi.M(),pi.Pt()) ;
          if(cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta))
            FillHistogram("PHOS_Mi_mvsPt_Off_Kink",pi.M(),pi.Pt()) ;
            if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta))
            FillHistogram("PHOS_Mi_mvsPt_Off_dEdx",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && !cnv->TestBit(kConvEta))
            FillHistogram("PHOS_Mi_mvsPt_Off_Prob",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta) && cnv->TestBit(kConvR))
            FillHistogram("PHOS_Mi_mvsPt_Off_R120",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta) && cnv->TestBit(kConvZR))
            FillHistogram("PHOS_Mi_mvsPt_Off_Z",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta) && cnv->TestBit(kConvNDF))
            FillHistogram("PHOS_Mi_mvsPt_Off_chi",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb))
            FillHistogram("PHOS_Mi_mvsPt_Off_Eta",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta) && cnv->TestBit(kConvPlan))
            FillHistogram("PHOS_Mi_mvsPt_Off_Wcut",pi.M(),pi.Pt()) ;
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
        if(cnv->TestBit(kConvOnFly) && !cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta)){
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
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta))
            FillHistogram("PHOS_Mi_mvsPt_OnFly",pi.M(),pi.Pt()) ;
          if(cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta))
            FillHistogram("PHOS_Mi_mvsPt_On_Kink",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta))
            FillHistogram("PHOS_Mi_mvsPt_On_dEdx",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && !cnv->TestBit(kConvEta))
            FillHistogram("PHOS_Mi_mvsPt_On_Prob",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta) && cnv->TestBit(kConvR))
            FillHistogram("PHOS_Mi_mvsPt_On_R120",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta) && cnv->TestBit(kConvZR))
            FillHistogram("PHOS_Mi_mvsPt_On_Z",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta) && cnv->TestBit(kConvNDF))
            FillHistogram("PHOS_Mi_mvsPt_On_chi",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb))
            FillHistogram("PHOS_Mi_mvsPt_On_Eta",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta) && cnv->TestBit(kConvPlan))
            FillHistogram("PHOS_Mi_mvsPt_On_Wcut",pi.M(),pi.Pt()) ;
        }
        else{
          //Default
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta))
            FillHistogram("PHOS_Mi_mvsPt_Offline",pi.M(),pi.Pt()) ;
          if(cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta))
            FillHistogram("PHOS_Mi_mvsPt_Off_Kink",pi.M(),pi.Pt()) ;
            if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta))
            FillHistogram("PHOS_Mi_mvsPt_Off_dEdx",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && !cnv->TestBit(kConvEta))
            FillHistogram("PHOS_Mi_mvsPt_Off_Prob",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta) && cnv->TestBit(kConvR))
            FillHistogram("PHOS_Mi_mvsPt_Off_R120",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta) && cnv->TestBit(kConvZR))
            FillHistogram("PHOS_Mi_mvsPt_Off_Z",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta) && cnv->TestBit(kConvNDF))
            FillHistogram("PHOS_Mi_mvsPt_Off_chi",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb))
            FillHistogram("PHOS_Mi_mvsPt_Off_Eta",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta) && cnv->TestBit(kConvPlan))
            FillHistogram("PHOS_Mi_mvsPt_Off_Wcut",pi.M(),pi.Pt()) ;
        }
      }
    }
  }
 
  //PHOS module dependent
  for(Int_t iPHOS=0; iPHOS<nPHOS;iPHOS++){
    TLorentzVector * cal = static_cast<TLorentzVector*>(fPHOSEvent->At(iPHOS)) ;
    Int_t mod=1;
    while(!cal->TestBit(BIT(16+mod)) && mod<5)mod++ ;
    TString base("PHOS_Mi_mvsPt_mod") ; base+=mod ;
    TString full ;
    for(Int_t ev=0; ev<prevConv->GetSize();ev++){
      TClonesArray * mixConv = static_cast<TClonesArray*>(prevConv->At(ev)) ;
      for(Int_t iConv = 0; iConv<mixConv->GetEntriesFast(); iConv++){
        TLorentzVector * cnv = static_cast<TLorentzVector*>(mixConv->At(iConv)) ;
        TLorentzVector pi=*cal + *cnv ;
        full=base+"_all" ;
        if(cnv->TestBit(kConvOnFly) && !cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta)){
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
        while(!cal->TestBit(BIT(16+mod)) && mod<5)mod++ ;
        TString base("PHOS_Mi_mvsPt_mod") ; base+=mod ;
        TString full ;
        TLorentzVector pi=*cal + *cnv ;
        full=base+"_all" ;
        if(cnv->TestBit(kConvOnFly) && !cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta)){
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
      }
    }
  }


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
        if(cnv->TestBit(kConvOnFly) && !cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta)){
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
            full=base+"_Neutral" ;
            FillHistogram(full,pi.M(),pi.Pt()) ;
          }
        } 
        if(cnv->TestBit(kConvOnFly)){
          //Default
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta))
            FillHistogram("EMCAL_Mi_mvsPt_OnFly",pi.M(),pi.Pt()) ;
          if(cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta))
            FillHistogram("EMCAL_Mi_mvsPt_On_Kink",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta))
            FillHistogram("EMCAL_Mi_mvsPt_On_dEdx",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && !cnv->TestBit(kConvEta))
            FillHistogram("EMCAL_Mi_mvsPt_On_Prob",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta) && cnv->TestBit(kConvR))
            FillHistogram("EMCAL_Mi_mvsPt_On_R120",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta) && cnv->TestBit(kConvZR))
            FillHistogram("EMCAL_Mi_mvsPt_On_Z",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta) && cnv->TestBit(kConvNDF))
            FillHistogram("EMCAL_Mi_mvsPt_On_chi",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb))
            FillHistogram("EMCAL_Mi_mvsPt_On_Eta",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta) && cnv->TestBit(kConvPlan))
            FillHistogram("EMCAL_Mi_mvsPt_On_Wcut",pi.M(),pi.Pt()) ;
        }
        else{
          //Default
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta))
            FillHistogram("EMCAL_Mi_mvsPt_Offline",pi.M(),pi.Pt()) ;
          if(cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta))
            FillHistogram("EMCAL_Mi_mvsPt_Off_Kink",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta))
            FillHistogram("EMCAL_Mi_mvsPt_Off_dEdx",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && !cnv->TestBit(kConvEta))
            FillHistogram("EMCAL_Mi_mvsPt_Off_Prob",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta) && cnv->TestBit(kConvR))
            FillHistogram("EMCAL_Mi_mvsPt_Off_R120",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta) && cnv->TestBit(kConvZR))
            FillHistogram("EMCAL_Mi_mvsPt_Off_Z",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta) && cnv->TestBit(kConvNDF))
            FillHistogram("EMCAL_Mi_mvsPt_Off_chi",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb))
            FillHistogram("EMCAL_Mi_mvsPt_Off_Eta",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta) && cnv->TestBit(kConvPlan))
            FillHistogram("EMCAL_Mi_mvsPt_Off_Wcut",pi.M(),pi.Pt()) ;
        }
      }
    }
  }
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
        if(cnv->TestBit(kConvOnFly) && !cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta)){
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
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta))
            FillHistogram("EMCAL_Mi_mvsPt_OnFly",pi.M(),pi.Pt()) ;
          if(cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta))
            FillHistogram("EMCAL_Mi_mvsPt_On_Kink",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta))
            FillHistogram("EMCAL_Mi_mvsPt_On_dEdx",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && !cnv->TestBit(kConvEta))
            FillHistogram("EMCAL_Mi_mvsPt_On_Prob",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta) && cnv->TestBit(kConvR))
            FillHistogram("EMCAL_Mi_mvsPt_On_R120",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta) && cnv->TestBit(kConvZR))
            FillHistogram("EMCAL_Mi_mvsPt_On_Z",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta) && cnv->TestBit(kConvNDF))
            FillHistogram("EMCAL_Mi_mvsPt_On_chi",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb))
            FillHistogram("EMCAL_Mi_mvsPt_On_Eta",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta) && cnv->TestBit(kConvPlan)) 
            FillHistogram("EMCAL_Mi_mvsPt_On_Wcut",pi.M(),pi.Pt()) ;
        }
        else{
          //Default
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta))
            FillHistogram("EMCAL_Mi_mvsPt_Offline",pi.M(),pi.Pt()) ;
          if(cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta))
            FillHistogram("EMCAL_Mi_mvsPt_Off_Kink",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta))
            FillHistogram("EMCAL_Mi_mvsPt_Off_dEdx",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && !cnv->TestBit(kConvEta))
            FillHistogram("EMCAL_Mi_mvsPt_Off_Prob",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta) && cnv->TestBit(kConvR))
            FillHistogram("EMCAL_Mi_mvsPt_Off_R120",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta) && cnv->TestBit(kConvZR)) 
            FillHistogram("EMCAL_Mi_mvsPt_Off_Z",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta) && cnv->TestBit(kConvNDF))
            FillHistogram("EMCAL_Mi_mvsPt_Off_chi",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb))
            FillHistogram("EMCAL_Mi_mvsPt_Off_Eta",pi.M(),pi.Pt()) ;
          if(!cnv->TestBit(kConvKink) && cnv->TestBit(kConvdEdx) && cnv->TestBit(kConvProb) && !cnv->TestBit(kConvEta) && cnv->TestBit(kConvPlan))
            FillHistogram("EMCAL_Mi_mvsPt_Off_Wcut",pi.M(),pi.Pt()) ;
        }
      }
    }
  }


      
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

  //fill histograms for efficiensy etc. calculation
  if(!fStack) return ;
  
  const Double_t rcut = 1. ; //cut for primary particles
  Double_t vtx[3];
  vtx[0] = fESDEvent->GetPrimaryVertex()->GetX();
  vtx[1] = fESDEvent->GetPrimaryVertex()->GetY();
  vtx[2] = fESDEvent->GetPrimaryVertex()->GetZ();

  //---------First pi0s-----------------------------
  for (Int_t iTracks = 0; iTracks < fStack->GetNtrack(); iTracks++) {
    TParticle* particle = (TParticle *)fStack->Particle(iTracks);
    if(particle->GetPdgCode() != 111)
      continue ;
     
    //Primary particle
    if(particle->R() >rcut)
      continue ;

    Double_t pt = particle->Pt() ;
    //Total number of pi0 with creation radius <1 cm
    FillHistogram("hMC_CaloConv_allPi0",pt) ;
    if(TMath::Abs(particle->Y())<1.)
      FillHistogram("hMC_CaloConv_pi0_unitEta",pt) ;

    //Check if one of photons converted
    if(particle->GetNDaughters()!=2)
     continue ; //Do not account Dalitz decays

    TParticle * gamma1 = fStack->Particle(particle->GetFirstDaughter());
    TParticle * gamma2 = fStack->Particle(particle->GetLastDaughter());
    //Number of pi0s decayed into acceptance
    Bool_t inAcc1 = (TMath::Abs(gamma1->Eta())<0.9) ;
    Bool_t inAcc2 = (TMath::Abs(gamma2->Eta())<0.9) ;
    Int_t mod ;
    Double_t x,z ;
    Bool_t hitPHOS1 = fPHOSgeom->ImpactOnEmc(gamma1, mod, z,x) ;
    Bool_t hitPHOS2 = fPHOSgeom->ImpactOnEmc(gamma2, mod, z,x) ;
    Bool_t hitEMCAL1= fEMCALgeom->Impact(gamma1) ;
    Bool_t hitEMCAL2= fEMCALgeom->Impact(gamma2) ;
 
    Bool_t goodPair=kFALSE ;
    if((inAcc1 && hitPHOS2) || (inAcc2 && hitPHOS1)){
      FillHistogram("hMC_CaloConv_pi0PHOSacc",pt) ;
      goodPair=kTRUE ;
    } 
    if((inAcc1 && hitEMCAL2) || (inAcc2 && hitEMCAL1)){
      FillHistogram("hMC_CaloConv_pi0EMCALacc",pt) ;
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
    if((converted1 && !converted2 && hitPHOS2) || (!converted1 && hitPHOS1 && converted2)) 
       FillHistogram("hMC_CaloConv_pi0_PHOS_conv",pt) ;
 
    if((converted1 && !converted2 && hitEMCAL2) || (!converted1 && hitEMCAL1 && converted2)) 
       FillHistogram("hMC_CaloConv_pi0_EMCAL_conv",pt) ;
 
    //Both converted
    if(converted1 && converted2) {
      FillHistogram("hMC_CaloConv_pi0_bothphot_conv",pt) ;
        continue ;
    }
 
    //photon pointing calorimeter converted
    if((converted1 && hitPHOS1 && !hitEMCAL2) || (converted2 && hitPHOS2 && !hitEMCAL1) || 
       (converted1 && hitEMCAL1 && !hitPHOS2) || (converted2 && hitEMCAL2 && !hitPHOS1)){
       FillHistogram("hMC_CaloConv_pi0__convPhotInCalo",pt) ;
       continue ;
    }
   
    //Converted pi0 with v0 and photon PHOS or EMCAL
    Bool_t foundV01=kFALSE, foundV02=kFALSE ;
    TLorentzVector pConv ;
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
      Bool_t same = (v0Gamma == gamma1) ;
      TParticle * tmp = v0Gamma ;
      while(!same && tmp->GetFirstMother()>=0){
        tmp = fStack->Particle(tmp->GetFirstMother());
       same = (tmp == gamma1) ;
     }
     if(same){
       foundV01 = kTRUE ;
       const AliExternalTrackParam * paramPos = v0->GetParamP() ;
       const AliExternalTrackParam * paramNeg = v0->GetParamN() ;
       AliKFParticle negKF(*paramNeg,11);
       AliKFParticle posKF(*paramPos,-11);
       pConv.SetXYZM(negKF.Px()+posKF.Px(),negKF.Py()+posKF.Py(),negKF.Pz()+negKF.Pz(),0.) ;
       break ;
     } 
     same = (v0Gamma == gamma2) ;
     tmp = v0Gamma ;
     while(!same && tmp->GetFirstMother()>=0){
       tmp = fStack->Particle(tmp->GetFirstMother());
       same = (tmp == gamma2) ;
     }
     if(same){
       foundV02 = kTRUE ;
       const AliExternalTrackParam * paramPos = v0->GetParamP() ;
       const AliExternalTrackParam * paramNeg = v0->GetParamN() ;
       AliKFParticle negKF(*paramNeg,11);
       AliKFParticle posKF(*paramPos,-11);
       pConv.SetXYZM(negKF.Px()+posKF.Px(),negKF.Py()+posKF.Py(),negKF.Pz()+negKF.Pz(),0.) ;
       break ;
     } 
   }

   goodPair=kFALSE ;
   if((foundV01 && hitPHOS2) || (foundV02 && hitPHOS1)){
     FillHistogram("hMC_CaloConv_pi0_v0_PHOSacc",pt) ;
     goodPair=kTRUE;
   }
   if((foundV01 && hitEMCAL2) || (foundV02 && hitEMCAL1)){
     FillHistogram("hMC_CaloConv_pi0_v0_EMCALacc",pt) ;
     goodPair=kTRUE;
   }
   if(!goodPair){
     continue ;
   }

   //Converted pi0 with v0 and cluster in PHOS/EMCAL
   Bool_t cluInPHOS = kFALSE,cluInEMCAL=kFALSE ;
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
       clu->GetMomentum(pCalo ,vtx);
       break ;
     }
     if(!clu->IsPHOS() && (hitEMCAL1 || hitEMCAL2)){
       cluInEMCAL=kTRUE ;
       clu->GetMomentum(pCalo ,vtx);
       break ;
     }
   }

   if(cluInPHOS){
     FillHistogram("hMC_CaloConv_pi0_v0_PHOSclu",pt) ;
     Double_t m=(pCalo+pConv).M() ;
     Double_t ptm=(pCalo+pConv).Pt() ;
     FillHistogram("hMC_CaloConv_pi0_v0_PHOSclu_ptRec",ptm) ;
     FillHistogram("hMC_CaloConv_pi0_v0_PHOSclu_mvsPt",m,ptm) ;
   }
   if(cluInEMCAL){
     FillHistogram("hMC_CaloConv_pi0_v0_EMCALclu",pt) ;
     Double_t m=(pCalo+pConv).M() ;
     Double_t ptm=(pCalo+pConv).Pt() ;
     FillHistogram("hMC_CaloConv_pi0_v0_EMCALclu_ptRec",ptm) ;
     FillHistogram("hMC_CaloConv_pi0_v0_EMCALclu_mvsPt",m,ptm) ;
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
        if(!cnv->TestBit(kConvOnFly) || cnv->TestBit(kConvKink) || !cnv->TestBit(kConvdEdx) || !cnv->TestBit(kConvProb) || cnv->TestBit(kConvEta)) 
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
        if(!cnv->TestBit(kConvOnFly) || cnv->TestBit(kConvKink) || !cnv->TestBit(kConvdEdx) || !cnv->TestBit(kConvProb) || cnv->TestBit(kConvEta)) 
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
     Double_t x,z ;
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
         break ;
       }
       if(!clu->IsPHOS() && hitEMCAL){
         cluInEMCAL=kTRUE ;
         clu->GetMomentum(pCalo ,vtx);
         break ;
       }
     }

     if(cluInPHOS){
       FillHistogram("hMC_CaloConv_gamma_PHOSclu",pt) ;
       FillHistogram("hMC_CaloConv_gamma_PHOSclu_recE",pCalo.E()) ;
       FillHistogram("hMC_CaloConv_gamma_PHOSclu_devsE",(particle->Energy()-pCalo.E())/particle->Energy(),particle->Energy()) ;
     }
     if(cluInEMCAL){
       FillHistogram("hMC_CaloConv_gamma_EMCALclu",pt) ;
       FillHistogram("hMC_CaloConv_gamma_EMCALclu_recE",pCalo.E()) ;
       FillHistogram("hMC_CaloConv_gamma_EMCALclu_devsE",(particle->Energy()-pCalo.E())/particle->Energy(),particle->Energy()) ;
     }
  }
}
//_____________________________________________________________________________
void AliAnalysisTaskCaloConv::FillHistogram(const char * key,Double_t x)const{
  TH1F * tmp = dynamic_cast<TH1F*>(fOutputContainer->FindObject(key)) ;
  if(!tmp)
    AliError(Form("can not find histogram <%s> ",key)) ;
  tmp->Fill(x) ;
}
//_____________________________________________________________________________
void AliAnalysisTaskCaloConv::FillHistogram(const char * key,Double_t x,Double_t y)const{
  TObject * tmp = fOutputContainer->FindObject(key) ;
  if(!tmp)
    AliError(Form("can not find histogram <%s> ",key)) ;
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
  if(!tmp)
    AliError(Form("can not find histogram <%s> ",key)) ;
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

 



