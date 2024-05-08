#include "TChain.h"
#include "TTree.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TKey.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1I.h"
#include "TH2I.h"
#include "TH3F.h"
#include "TProfile.h"
#include "TParticle.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TRandom.h"
#include "TROOT.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "THashList.h"
#include "AliMagF.h"
#include "TGeoGlobalMagField.h"

#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskGammaPHOSPbPbRun2.h"
#include "AliCaloPhoton.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSEsdCluster.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliESDCaloCells.h"
#include "AliAODVertex.h"
#include "AliESDtrackCuts.h"
#include "AliLog.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
#include "AliCDBManager.h"
#include "AliCentrality.h" 
#include "AliEventplane.h"
#include "AliOADBContainer.h"
#include "AliEPFlattener.h"
#include "AliFlowTrackCuts.h"
#include "AliFlowEvent.h"
#include "AliFlowVector.h"
#include "AliMultSelection.h"
#include "AliAODMCParticle.h"

// Analysis task to fill histograms with PHOS ESD clusters and cells
// Authors: Dmitri Peressounko
// Date   : 28.05.2011

ClassImp(AliAnalysisTaskGammaPHOSPbPbRun2)

//________________________________________________________________________
AliAnalysisTaskGammaPHOSPbPbRun2::AliAnalysisTaskGammaPHOSPbPbRun2(const char *name) 
: AliAnalysisTaskSE(name),
  fOutputContainer(0x0),
  fOutputContainer2(0x0),
  fOutputContainer3(0x0),
  fOutputContainer4(0x0),
  fOutputContainer5(0x0),
  fOutputContainer6(0x0),
  fEvent(0x0),
  fPHOSEvent(0x0),
  fPIDResponse(0x0), 
  fQV0A(0.),
  fQV0C(0.),
  fQTPC(0.),
  fRP(0.),
  fRPV0A(0.),
  fRPV0C(0.),
  fRPQ(0.),
  fRPQV0A(0.),
  fRPQV0C(0.),
  fV0Ares(1.),
  fV0Cres(1.),
  fTPCres(1.),
  fV0AQres(1.),
  fV0CQres(1.),
  fTPCQres(1.),
  fHaveTPCRP(0),
  fPhiDist(0x0),
  fV0AFlat(0x0),
  fV0CFlat(0x0),
  fTPCFlat(0x0),
  fV0AQFlat(0x0),
  fV0CQFlat(0x0),
  fTPCQFlat(0x0),
  fHarmonics(3),
  fDistCut(0),
  fRunNumber(0),
  fCentrality(0.),
  fCenBin(0),
  fPHOSGeo(0x0),
  fEventCounter(0),
  fInPHOS(0),

  fTPCfinalC2(0x0),   //HIstos with flattening parameters
  fTPCfinalS2(0x0),
  fTPCfinalC4(0x0),
  fTPCfinalS4(0x0),

  fV0AfinalC2(0x0),
  fV0AfinalS2(0x0),
  fV0AfinalC4(0x0),
  fV0AfinalS4(0x0),
    
  fV0CfinalC2(0x0),
  fV0CfinalS2(0x0),
  fV0CfinalC4(0x0),
  fV0CfinalS4(0x0), 
  fTPCfinalQC2(0x0),   //HIstos with flattening parameters
  fTPCfinalQS2(0x0),
  fTPCfinalQC4(0x0),
  fTPCfinalQS4(0x0),

  fV0AfinalQC2(0x0),
  fV0AfinalQS2(0x0),
  fV0AfinalQC4(0x0),
  fV0AfinalQS4(0x0),
    
  fV0CfinalQC2(0x0),
  fV0CfinalQS2(0x0),
  fV0CfinalQC4(0x0),
  fV0CfinalQS4(0x0),
  fCutsV0(0x0),
  fCutsTPC(0x0),
  fFlowEvent(0x0),
  
  fMCArray(0x0),
  fIsMC(kFALSE),
  fPidCuts(0),
  fCenBinEdges(0),
  fTOF(30.e-9),
  fNCenBins(0),
  fCurrFileName(0), 
  fCheckMCCrossSection(kFALSE),
  fh1Xsec(0), fh1Trials(0), fAvgTrials(-1)

{
  // Constructor
  for(Int_t i=0; i<1;i++){
    for(Int_t j=0; j<10;j++)
      for(Int_t k=0; k<11;k++)
	fPHOSEvents[i][j][k]=0;
  }

  // define cuts 
  fPidCuts.emplace_back("All", "no cuts");
  fPidCuts.emplace_back("Allcore", "no cuts, core energy");
  fPidCuts.emplace_back("CPV", "cpv");
  fPidCuts.emplace_back("CPVcore", "cpv, core energy");
  fPidCuts.emplace_back("Disp2", "disp");
  fPidCuts.emplace_back("Disp2core", "disp, core energy");
  fPidCuts.emplace_back("Both2", "cpv+disp");
  fPidCuts.emplace_back("Both2core", "cpv+disp, core energy");

  //fCenBinEdges = {0., 10., 20., 30., 40., 50., 80., 100};

  //fNCenBins = fCenBinEdges.size() - 1;

  // Output slots #0 write into a TH1 container

  // Initialize the PHOS geometry
  fPHOSGeo = AliPHOSGeometry::GetInstance("Run2");

  //We have to apply re-calibration for pass1 LCH10h
  // Initialize decalibration factors in the form of the OCDB object
  fCutsV0 = new AliFlowTrackCuts(Form("V0%d",fHarmonics));
  fCutsV0 = fCutsV0->GetStandardVZEROOnlyTrackCuts(); // select vzero tracks
  fCutsV0->SetVZEROgainEqualizationPerRing(kFALSE);
  fCutsV0->SetApplyRecentering(kTRUE);
//  fCutsV0A->SetEtaRange(2.,10.); 
  fFlowEvent = new AliFlowEvent(10000);

  fCutsTPC= new AliFlowTrackCuts(Form("TPC%d",fHarmonics));
  fCutsTPC=fCutsTPC->GetStandardTPCStandaloneTrackCuts();
  fCutsTPC->SetEtaMin(-0.9);
  fCutsTPC->SetEtaMax(0.9);

  DefineOutput(1,THashList::Class());
  DefineOutput(2,THashList::Class());
  DefineOutput(3,THashList::Class());
  DefineOutput(4,THashList::Class());
  DefineOutput(5,THashList::Class());
  DefineOutput(6,THashList::Class());
 }

//________________________________________________________________________
void AliAnalysisTaskGammaPHOSPbPbRun2::UserCreateOutputObjects()
{

  // Create histograms
  // Called once
  
  if (fOutputContainer != NULL) {
    delete fOutputContainer;
  }

  if (fOutputContainer2 != NULL) {
    delete fOutputContainer2;
  }

  if (fOutputContainer3 != NULL) {
    delete fOutputContainer3;
  }

  if (fOutputContainer4 != NULL) {
    delete fOutputContainer4;
  }

  if (fOutputContainer5 != NULL) {
    delete fOutputContainer5;
  }

  if (fOutputContainer6 != NULL) {
    delete fOutputContainer6;
  }

  fOutputContainer   = new THashList();
  fOutputContainer2  = new THashList();
  fOutputContainer3  = new THashList();
  fOutputContainer4  = new THashList();
  fOutputContainer5  = new THashList();
  fOutputContainer6  = new THashList();

  fOutputContainer ->SetOwner(kTRUE);
  fOutputContainer2->SetOwner(kTRUE);
  fOutputContainer3->SetOwner(kTRUE);
  fOutputContainer4->SetOwner(kTRUE);
  fOutputContainer5->SetOwner(kTRUE);
  fOutputContainer6->SetOwner(kTRUE);
  
  //========QA histograms=======

  fh1Xsec = new TH1F("hXsec","xsec from pyxsec.root", 1, 0, 1);
  fh1Xsec->GetXaxis()->SetBinLabel(1,"<#sigma>");
  fOutputContainer->Add(fh1Xsec);
  
  fh1Trials = new TH1F("hTrials","trials root file", 1, 0, 1);
  fh1Trials->GetXaxis()->SetBinLabel(1,"#sum{ntrials}");
  fOutputContainer->Add(fh1Trials);

  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  if (!man) {
      AliFatal("Could not find manager");
  }	

  AliInputEventHandler* inputHandler = dynamic_cast<AliInputEventHandler*> (man->GetInputEventHandler());
  if (!inputHandler) {
      AliFatal("No input event handler");
  }	 

  fPIDResponse = dynamic_cast<AliPIDResponse *>(inputHandler->GetPIDResponse());
  if (!fPIDResponse) {
      AliFatal("PIDResponse object was not created"); // Escalated to fatal. This task is unusable without PID response.
  }	

  fPIDResponse -> SetUseTPCMultiplicityCorrection(kFALSE);
  fPIDResponse->SetCurrentMCEvent(MCEvent()); //!!

  AddQAHistograms(); 
  AddSinglePhotonHistograms();
  AddEventPlaneHistograms();
  AddMCHistograms();
  AddDistBadHistograms();
  AddPhiTitleHistograms();
 
  PostData(1, fOutputContainer);
  PostData(2, fOutputContainer2);
  PostData(3, fOutputContainer3);
  PostData(4, fOutputContainer4);
  PostData(5, fOutputContainer5);
  PostData(6, fOutputContainer6);

}

//________________________________________________________________________
void AliAnalysisTaskGammaPHOSPbPbRun2::UserExec(Option_t *) 
{
  // Main loop, called for each event
  // Analyze ESD/AOD
  
  fEvent = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fEvent) {
    Printf("ERROR: Could not retrieve event");
    return;
  }

  fRunNumber = ConvertRunNumber(fEvent->GetRunNumber());
  FillHistogram("hSelEvents", 0.5, fRunNumber-0.5);
  FillHistogram("hTotSelEvents", 0.5);

  // PID
  fPIDResponse -> SetUseTPCMultiplicityCorrection(kFALSE);
  AliPIDCombined *pidcomb=new AliPIDCombined();
  pidcomb->SetDefaultTPCPriors();
  pidcomb->SetSelectedSpecies(AliPID::kSPECIESC);
  pidcomb->SetDetectorMask(AliPIDResponse::kDetTPC|AliPIDResponse::kDetTOF|AliPIDResponse::kDetITS|AliPIDResponse::kDetTRD);
  
  // Get PHOS rotation matrices from ESD and set them to the PHOS geometry
  if(fEventCounter == 0) {
//   OpenInfoCalbration(fEvent->GetRunNumber());
   for(Int_t mod=0; mod<5; mod++) {
      if(!fEvent->GetPHOSMatrix(mod)) continue;
      fPHOSGeo->SetMisalMatrix(fEvent->GetPHOSMatrix(mod),mod);
      Printf("PHOS geo matrix %p for module # %d is set\n", fEvent->GetPHOSMatrix(mod), mod);
    }
    Int_t run = fEvent->GetRunNumber();
    Int_t runTemporary = 170593;  //temporary fix for flattening
  
//    TFile * fflatFine = TFile::Open("EP_final.root");
//    TFile * fflatFineQ = TFile::Open("EPq_final.root");
  
    AliOADBContainer flatContainer("phosFlat");
    flatContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSflat.root","phosFlat");
    TObjArray *arr = (TObjArray*)flatContainer.GetObject(/*run*/ runTemporary,"phosFlat");
    if(!arr){
      AliError(Form("Can not read Flattening for run %d. \n From file $ALICE_PHYSICS/OADB/PHOS/PHOSflat.root",/*run*/ runTemporary));    
      arr = (TObjArray*)flatContainer.GetObject(1,"phosFlat"); //default
    }
    AliOADBContainer flatContainerQ("phosQFlat");
    flatContainerQ.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSflat.root","phosQFlat");
    TObjArray *arrQ = (TObjArray*)flatContainerQ.GetObject(/*run*/ runTemporary,"phosQFlat");
    if(!arrQ){
      AliError(Form("Can not read QFlattening for run %d. \n From file $ALICE_PHYSICS/OADB/PHOS/PHOSflat.root",/*run*/ runTemporary));    
    }
        
    AliInfo(Form("Setting PHOS flattening with name %s \n",arr->GetName()));
    if(fHarmonics==2){
      AliEPFlattener * h = (AliEPFlattener*)arr->At(0);  
      if(fTPCFlat) delete fTPCFlat;
      fTPCFlat = new AliEPFlattener();
      fTPCFlat = h;
      h = (AliEPFlattener*)arr->At(1);  
      if(fV0AFlat) delete fV0AFlat;
      fV0AFlat = new AliEPFlattener();
      fV0AFlat = h;
      h = (AliEPFlattener*)arr->At(2);  
      if(fV0CFlat) delete fV0CFlat;
      fV0CFlat = new AliEPFlattener();
      fV0CFlat = h;
      
      h = (AliEPFlattener*)arrQ->At(0);  
      if(fTPCQFlat) delete fTPCQFlat;
      fTPCQFlat = new AliEPFlattener();
      fTPCQFlat = h;
      h = (AliEPFlattener*)arrQ->At(1);  
      if(fV0AQFlat) delete fV0AQFlat;
      fV0AQFlat = new AliEPFlattener();
      fV0AQFlat = h;
      h = (AliEPFlattener*)arrQ->At(2);  
      if(fV0CQFlat) delete fV0CQFlat;
      fV0CQFlat = new AliEPFlattener();
      fV0CQFlat = h;
    }    
    if(fHarmonics==3){
      AliEPFlattener * h = (AliEPFlattener*)arr->At(3);  
      if(fTPCFlat) delete fTPCFlat;
      fTPCFlat = new AliEPFlattener();
      fTPCFlat = h;
      h = (AliEPFlattener*)arr->At(4);  
      if(fV0AFlat) delete fV0AFlat;
      fV0AFlat = new AliEPFlattener();
      fV0AFlat = h;
      h = (AliEPFlattener*)arr->At(5);  
      if(fV0CFlat) delete fV0CFlat;
      fV0CFlat = new AliEPFlattener();
      fV0CFlat = h;
      
      h = (AliEPFlattener*)arrQ->At(3);  
      if(fTPCQFlat) delete fTPCQFlat;
      fTPCQFlat = new AliEPFlattener();
      fTPCQFlat = h;
      h = (AliEPFlattener*)arrQ->At(4);  
      if(fV0AQFlat) delete fV0AQFlat;
      fV0AQFlat = new AliEPFlattener();
      fV0AQFlat = h;
      h = (AliEPFlattener*)arrQ->At(5);  
      if(fV0CQFlat) delete fV0CQFlat;
      fV0CQFlat = new AliEPFlattener();
      fV0CQFlat = h;
    }
      
    // TPC Event Plane Weights
    AliOADBContainer *fEPContainer=NULL;
    TString oadbfilename = (Form("%s/COMMON/EVENTPLANE/data/epphidist.root", AliAnalysisManager::GetOADBPath()));

    TFile foadb(oadbfilename);
    if(!foadb.IsOpen()) AliFatal(Form("Cannot open OADB file %s", oadbfilename.Data()));

    AliInfo("Using Standard OADB");
    fEPContainer = (AliOADBContainer*) foadb.Get("epphidist");
    if (!fEPContainer) AliFatal("Cannot fetch OADB container for EP selection");
    foadb.Close();

    fPhiDist = (TH1F*) fEPContainer->GetObject(run, "Default");

    Bool_t emptybins;
    int iter = 0;
    while (iter<3){
      emptybins = kFALSE;
      for (int i=1; i<fPhiDist->GetNbinsX(); i++){
         if (!((fPhiDist->GetBinContent(i))>0)) {
            emptybins = kTRUE;
         }
      }
      if (emptybins) {
        fPhiDist->Rebin();
        iter++;
      }
      else iter = 3;
   }

    fEventCounter++;
  }

  // Checks if we have a primary vertex
  // Get primary vertices form ESD
  const AliAODVertex *esdVertex5 = fEvent->GetPrimaryVertex();

  fVtx5[0] = esdVertex5->GetX();
  fVtx5[1] = esdVertex5->GetY();
  fVtx5[2] = esdVertex5->GetZ();
  
  
  FillHistogram("hZvertex", esdVertex5->GetZ(), fRunNumber - 0.5);
  if (TMath::Abs(esdVertex5->GetZ()) > 10. ){
    return;
  }
  FillHistogram("hSelEvents",1.5,fRunNumber-0.5);
  FillHistogram("hTotSelEvents",1.5);

  FillHistogram("hSelEvents",2.5,fRunNumber-0.5);
  FillHistogram("hTotSelEvents",2.5);  
  
  
  //No dependence on zVtx observed, save memory
  Int_t zvtx=0;

  AliMultSelection *multSelection =static_cast<AliMultSelection*>(fEvent->FindListObject("MultSelection"));
  if(multSelection) fCentrality = multSelection->GetMultiplicityPercentile("V0M");

  if( fCentrality < 0. || fCentrality > 100. ){
    return;
  }

  FillHistogram("hSelEvents",3.5, fRunNumber-0.5);
  FillHistogram("hTotSelEvents", 3.5);

  fCenBin = 0;

  for (Int_t i = 0; i < fNCenBins; i ++) {
    if (fCentrality < fCenBinEdges[i + 1]) break;
    else 
    fCenBin++;
  }

  FillHistogram("hCentralityBins", fCenBin + 0.5);

  Double_t cWeight=CentralityWeight(fCentrality);

  //Calculate EP resolutions from centrality
  EvalResolution();
  EvalQResolution();
  
  EvalV0ReactionPlane(fEvent);

  //TPC evals second order RP need sign of V0 RP to find direction
  //reaction plane
  Double_t rpFull=0.,dPsi=0.;
  
  fHaveTPCRP=GetTPCEventPlane(rpFull,dPsi);
   
  if(fHaveTPCRP){
    //rpFull is second order EP, in the range 0-pi
    //Set correct direction
    while(rpFull<0)rpFull+=TMath::TwoPi()/fHarmonics;
    while(rpFull>TMath::TwoPi()/fHarmonics)rpFull-=(TMath::TwoPi()/fHarmonics);
    FillHistogram("phiRP",rpFull,fCentrality);    
    //Apply flattening
    fRP = fTPCFlat->MakeFlat(rpFull,fCentrality);
    while(fRP<0)fRP+=TMath::TwoPi()/fHarmonics;
    while(fRP>TMath::TwoPi()/fHarmonics)fRP-=(TMath::TwoPi()/fHarmonics);
    FillHistogram("phiRPQ",fRP,fCentrality,fQTPC); //Yes, there is no difference between RP and RPQ so far 
    FillHistogram("phiRPflat",fRP,fCentrality);      
  }

  ApplyFinalQFlattening();
  
  fCutsV0->SetEvent(fEvent,0x0);
  fFlowEvent->ClearFast();
  fFlowEvent->Fill(fCutsV0, fCutsTPC);
  fFlowEvent->TagSubeventsInEta(-10.,-1.,1.,10.);
  fFlowEvent->SetReferenceMultiplicity(fEvent->GetNumberOfTracks());
  fFlowEvent->DefineDeadZone(0., 0, 0, 0);
  AliFlowVector qArray[2];
  fFlowEvent->Get2Qsub(qArray,fHarmonics);
   
  
  Double_t x= qArray[1].Phi()/Double_t(fHarmonics);
  while(x<0)x+=TMath::TwoPi()/fHarmonics;
  while(x>TMath::TwoPi()/fHarmonics)x-=TMath::TwoPi()/fHarmonics;
  Double_t y= qArray[0].Phi()/Double_t(fHarmonics);
  while(y<0)y+=TMath::TwoPi()/fHarmonics;
  while(y>TMath::TwoPi()/fHarmonics)y-=TMath::TwoPi()/fHarmonics;
  FillHistogram("phiRPV0AFlow",x,fCentrality);
  FillHistogram("phiRPV0CFlow",y,fCentrality);
  
  FillHistogram("phiV0ACorrel",x,fRPV0A,fCentrality);
  FillHistogram("phiV0CCorrel",y,fRPV0C,fCentrality);
    
  FillHistogram("phiRPV0Aflat",fRPV0A,fCentrality);
  FillHistogram("phiRPV0Cflat",fRPV0C,fCentrality);
  
  FillHistogram("phiRPV0AQflat",fRPQV0A,fCentrality,fQV0A);
  FillHistogram("phiRPV0CQflat",fRPQV0C,fCentrality,fQV0C);

  FillHistogram("phiRPV0AC",fRPV0A,fRPV0C,fCentrality);

  FillHistogram("cos2V0ARP",TMath::Cos(fHarmonics*fRPQV0A),fCentrality,fQV0A);
  FillHistogram("sin2V0ARP",TMath::Sin(fHarmonics*fRPQV0A),fCentrality,fQV0A);
  FillHistogram("cos2V0CRP",TMath::Cos(fHarmonics*fRPQV0C),fCentrality,fQV0C);
  FillHistogram("sin2V0CRP",TMath::Sin(fHarmonics*fRPQV0C),fCentrality,fQV0C);
  FillHistogram("QV0A",fCentrality,fQV0A);
  FillHistogram("QV0C",fCentrality,fQV0C);
  FillHistogram("QTPC",fCentrality,fQTPC);
  
  if(fHaveTPCRP){
    FillHistogram("phiRPQflat",fRPQ,fCentrality,fQTPC);  
    FillHistogram("cos2TPCAC",TMath::Cos(fHarmonics*dPsi),fCentrality);
    FillHistogram("qcos2TPCAC",TMath::Cos(fHarmonics*dPsi),fCentrality,fQTPC);
    FillHistogram("cos2TPCRP",TMath::Cos(fHarmonics*fRP),fCentrality,fQTPC);
    FillHistogram("sin2TPCRP",TMath::Sin(fHarmonics*fRP),fCentrality,fQTPC);    

    Double_t cosAC = TMath::Cos(fHarmonics*(fRPV0A-fRPV0C));
    Double_t cosAT = TMath::Cos(fHarmonics*(fRP-fRPV0A));
    Double_t cosCT = TMath::Cos(fHarmonics*(fRP-fRPV0C));
    Double_t qcosAC = TMath::Cos(fHarmonics*(fRPQV0A-fRPQV0C));
    Double_t qcosAT = TMath::Cos(fHarmonics*(fRPQ-fRPQV0A));
    Double_t qcosCT = TMath::Cos(fHarmonics*(fRPQ-fRPQV0C));

    FillHistogram("phiRPV0ATPC",fRP,fRPV0A,fCentrality,fQV0A*fQTPC);

    FillHistogram("cos2V0AC",fCentrality,cosAC);
    FillHistogram("qcos2V0AC",fCentrality,qcosAC,fQV0A*fQV0C);
    
    FillHistogram("cos2V0ATPC",fCentrality,cosAT);
    FillHistogram("qcos2V0ATPC",fCentrality,qcosAT,fQV0A*fQTPC);

    FillHistogram("cos2V0CTPC",fCentrality,cosCT);
    FillHistogram("qcos2V0CTPC",fCentrality,qcosCT,fQV0C*fQTPC);

    FillHistogram("phiRPV0CTPC",fRP,fRPV0C,fCentrality,fQV0C*fQTPC);
    
    //calculate event-by-event resolutions
    if(cosAC!=0. && cosAT!=0. && cosCT!=0.){
       Double_t res = TMath::Sqrt(TMath::Abs(cosAC*cosCT/cosAT));
       FillHistogram("ResV0C",fCentrality,res); 
       FillHistogram("InvResV0C",fCentrality,1./res); 
        
       res = TMath::Sqrt(TMath::Abs(cosAC*cosAT/cosCT));
       FillHistogram("ResV0A",fCentrality,res); 
       FillHistogram("InvResV0A",fCentrality,1./res); 

       res = TMath::Sqrt(TMath::Abs(cosAT*cosCT/cosAC));
       FillHistogram("ResTPC",fCentrality,res); 
       FillHistogram("InvResTPC",fCentrality,1./res); 
    }
  }     
 
  FillHistogram("hSelEvents",4.5,fRunNumber-0.5);
  FillHistogram("hTotSelEvents",4.5);
  //All event selections done
  FillHistogram("hCentrality",fCentrality,fRunNumber-0.5);
  FillHistogram("hCentralityCorr",fCentrality,cWeight);
  //Reaction plane is defined in the range (0;pi)
  //We have 10 bins
  Double_t averageRP = fRPV0A+fRPV0C+fRP;
  if(fHaveTPCRP)
    averageRP/=3.;
  else
    averageRP/=2.;
  Int_t irp=Int_t(10.*averageRP*fHarmonics/TMath::TwoPi());
  if(irp>9)irp=9;
  
  if(!fPHOSEvents[zvtx][fCenBin][irp]) 
    fPHOSEvents[zvtx][fCenBin][irp]=new TList();
  TList * prevPHOS = fPHOSEvents[zvtx][fCenBin][irp];

  fMCArray = (TClonesArray*)fEvent->FindListObject(AliAODMCParticle::StdBranchName());

  if (fMCArray) fIsMC = kTRUE;
  
  ProcessMC();

  if(fPHOSEvent)
    fPHOSEvent->Clear();
  else
    fPHOSEvent = new TClonesArray("AliCaloPhoton",200);

  TVector3 vertex(fVtx5);
  
  Int_t multClust = fEvent->GetNumberOfCaloClusters();
  Int_t fInPHOS=0;  //inMod1=0, inMod2=0, inMod3=0;
  // Double_t avrgEm1=0,avrgEm2=0,avrgEm3=0; //average cluster energy

  AliAODCaloCells * cells = fEvent->GetPHOSCells();
  FillHistogram("hCenTrack",fCentrality, fEvent->GetNumberOfTracks());
  
  //QA PHOS cells
  Int_t nCellModule[4] = {0, 0, 0, 0};

  Int_t a[10]={0}; //left
  
  Int_t nCells=cells->GetNumberOfCells();

  for (Int_t iCell=0; iCell<nCells; iCell++) {
    Int_t cellAbsId = cells->GetCellNumber(iCell);
    Int_t relId[4];
    fPHOSGeo->AbsToRelNumbering(cellAbsId,relId);
    Int_t mod1  = relId[0];
    Int_t cellX = relId[2];
    Int_t cellZ = relId[3];
    // Test if event is complete 
    if(cellX<29)
      a[2*mod1]++;
    else
      a[2*mod1 + 1]++;

    Float_t energy = cells->GetAmplitude(iCell);
    FillHistogram("hCellEnergy",energy);
    if(mod1==1) {
      nCellModule[0]++;
      FillHistogram("hCellEnergyM1",cells->GetAmplitude(iCell));
      FillHistogram("hCellNXZM1",cellX,cellZ,1.);
      FillHistogram("hCellEXZM1",cellX,cellZ,energy);
    }
    else if (mod1==2) {
      nCellModule[1]++;
      FillHistogram("hCellEnergyM2",cells->GetAmplitude(iCell));
      FillHistogram("hCellNXZM2",cellX,cellZ,1.);
      FillHistogram("hCellEXZM2",cellX,cellZ,energy);
    }
    else if (mod1==3) {
      nCellModule[2]++;
      FillHistogram("hCellEnergyM3",cells->GetAmplitude(iCell));
      FillHistogram("hCellNXZM3",cellX,cellZ,1.);
      FillHistogram("hCellEXZM3",cellX,cellZ,energy);
    }
    else if (mod1==4) {
      nCellModule[3]++;
      FillHistogram("hCellEnergyM4",cells->GetAmplitude(iCell));
      FillHistogram("hCellNXZM4",cellX,cellZ,1.);
      FillHistogram("hCellEXZM4",cellX,cellZ,energy);
    }
  }

  for(Int_t iMod = 2; iMod < 8; iMod ++){
    if(a[iMod] == 0){
       FillHistogram("hBadMod",float(iMod));
       FillHistogram("hBadCentrality",fCentrality);
       AliError("Bad centrality!");
       return;
    }
  }

  FillHistogram("hCellMultEventM1",nCellModule[0]);
  FillHistogram("hCellMultEventM2",nCellModule[1]);
  FillHistogram("hCellMultEventM3",nCellModule[2]);
  FillHistogram("hCellMultEventM4",nCellModule[3]);
  FillHistogram("hCenPHOSCells",fCentrality,nCells);
  FillHistogram("hCenPHOSCellsM12",fCentrality,nCellModule[0],nCellModule[1]);
  FillHistogram("hCenPHOSCellsM13",fCentrality,nCellModule[0],nCellModule[2]);
  FillHistogram("hCenPHOSCellsM14",fCentrality,nCellModule[0],nCellModule[3]);
  FillHistogram("hCenPHOSCellsM23",fCentrality,nCellModule[1],nCellModule[2]);
  FillHistogram("hCenPHOSCellsM24",fCentrality,nCellModule[1],nCellModule[3]);
  FillHistogram("hCenPHOSCellsM34",fCentrality,nCellModule[2],nCellModule[3]);
  
  
  TVector3 localPos;
  Int_t nAll=0,nCPV=0, nDisp=0, nBoth=0; 

  for (Int_t i=0; i < multClust; i++) {
    AliAODCaloCluster *clu = fEvent->GetCaloCluster(i);

    //cluster selection
    if ( !clu->IsPHOS() || clu->E()<0.3) 
      continue;
    if (!fMCArray && (TMath::Abs(clu->GetTOF()) > fTOF)) 
      continue; // TOF cut for real data only!
    if ((fRunNumber > 209122) && (clu->GetType() != AliVCluster::kPHOSNeutral)) 
       continue; //Run2 cut
    if ((clu->GetM02()) < 0.2 && (clu->E() > 1.))
       continue; // exotic clusters
    if(fDistCut && (clu->GetDistanceToBadChannel() < 2.5))
       continue;

    Float_t  position[3];
    clu->GetPosition(position);
    TVector3 global(position);
    Int_t relId[4];
    fPHOSGeo->GlobalPos2RelId(global,relId);
    Int_t mod  = relId[0];
    Int_t cellX = relId[2];
    Int_t cellZ = relId[3];
    
    
    FillHistogram("hCluEvsClu",clu->E(),clu->GetNCells());
    FillHistogram("hCluEvsCluM",clu->E(),clu->GetM02());

    FillHistogram(Form("hTofM%d",mod),clu->E(),clu->GetTOF());
    if(clu->E()>1.)
      FillHistogram("hCenTOF",fCentrality,clu->GetTOF());
    if((clu->GetTOF()>1.5e-7) || (clu->GetTOF() <-2.5e-7) )
      continue;
    
    //Apply re-Calibration
    if(clu->E() < 0.3) continue;
    if(clu->GetNCells() < 3) continue;
    if(clu->GetM02() < 0.2)   continue;
      
    TLorentzVector pv1;
    clu->GetMomentum(pv1 ,fVtx5);

    Double_t ecore=clu->GetCoreEnergy(); 
    
    FillHistogram(Form("hCluLowM%d", mod), cellX, cellZ, 1.);
    if(clu->E() > 1.5){
      FillHistogram(Form("hCluHighM%d", mod), cellX, cellZ, 1.);
    }
    
    if(fInPHOS>=fPHOSEvent->GetSize()){
      fPHOSEvent->Expand(fInPHOS+50);
    }
    new((*fPHOSEvent)[fInPHOS]) AliCaloPhoton(pv1.X(),pv1.Py(),pv1.Z(),pv1.E());
    AliCaloPhoton * ph = (AliCaloPhoton*)fPHOSEvent->At(fInPHOS);
    ph->SetModule(mod);
    pv1*= ecore/pv1.E();
    ph->SetMomV2(&pv1);
    ph->SetNCells(clu->GetNCells());    

    ph->SetPrimary(GetPrimaryLabel(clu));
    ph->SetPrimaryAtVertex(GetPrimaryLabelAtVertex(clu));
    
//  PID cuts
    Bool_t dispBit1  = clu->GetDispersion() < 2.5 * 2.5;
    Bool_t dispBit2 = clu->Chi2() < 2.5 * 2.5;
    ph->SetDispBit(dispBit1);
    ph->SetDisp2Bit(dispBit2);
    ph->SetTOFBit(fMCArray ? kTRUE : clu->GetTOF() < fTOF);
    ph->SetTime(clu->GetTOF());

    //Track matching
    Bool_t cpvBit=(clu->GetEmcCpvDistance()>2.5);
    ph->SetCPVBit(cpvBit);   

    if(ph->IsDisp2OK()){
      FillHistogram(Form("hCluDispM%d",mod),cellX,cellZ,1.);
    }
    
    if(cpvBit){
      FillHistogram(Form("hCluVetoM%d",mod),cellX,cellZ,1.);
    }

    Double_t distBC = clu->GetDistanceToBadChannel();
    if(distBC > 3.){
      FillHistogram(Form("hPhotAll_DistBad2_cen%d",fCenBin),ph->Pt(),cWeight);
      FillHistogram(Form("hPhotAll_DistBad2core_cen%d",fCenBin),ph->GetMomV2()->Pt(),cWeight);
      if(ph->IsDisp2OK()){
        FillHistogram(Form("hPhotAll_DistBad2Disp2_cen%d",fCenBin),ph->Pt(),cWeight);
        FillHistogram(Form("hPhotAll_DistBad2Disp2core_cen%d",fCenBin),ph->GetMomV2()->Pt(),cWeight);	
      }
      if(distBC>4.){
        FillHistogram(Form("hPhotAll_DistBad4_cen%d",fCenBin),ph->Pt(),cWeight);
        FillHistogram(Form("hPhotAll_DistBad4core_cen%d",fCenBin),ph->GetMomV2()->Pt(),cWeight);
        if(ph->IsDisp2OK()){
          FillHistogram(Form("hPhotAll_DistBad4Disp2_cen%d",fCenBin),ph->Pt(),cWeight);
          FillHistogram(Form("hPhotAll_DistBad4Disp2core_cen%d",fCenBin),ph->GetMomV2()->Pt(),cWeight);	
        }
        if(distBC>6.){
          FillHistogram(Form("hPhotAll_DistBad6_cen%d",fCenBin),ph->Pt(),cWeight);
          FillHistogram(Form("hPhotAll_DistBad6core_cen%d",fCenBin),ph->GetMomV2()->Pt(),cWeight);
          if(ph->IsDisp2OK()){
            FillHistogram(Form("hPhotAll_DistBad6Disp2_cen%d",fCenBin),ph->Pt(),cWeight);
            FillHistogram(Form("hPhotAll_DistBad6Disp2core_cen%d",fCenBin),ph->GetMomV2()->Pt(),cWeight);	
          }	  
	}    
      }
    }
    
    
    ph->SetEMCx(float(cellX));
    ph->SetEMCz(float(cellZ));
    ph->SetLambdas(clu->GetM20(), clu->GetM02());
    ph->SetUnfolded(clu->GetNExMax() < 2); // Remember, if it is unfolded          
    nAll++;
    if(cpvBit)nCPV++;
    if(ph->IsDisp2OK())nDisp++;
    if(cpvBit && ph->IsDisp2OK()) nBoth++;
    fInPHOS++;
  }
  
  FillHistogram("hCenPHOS",fCentrality,fInPHOS,cWeight);
  if(fInPHOS==0){
    fEventCounter++;
    return; 
  }
  
  //Parameterization of the mean number of PHOS clusters vs centrality
  //This is 1/<Nclu> parameterization
  Double_t wA=1./fV0Ares;
  Double_t wC=1./fV0Cres;
  Double_t wT=1./fTPCres;
  Double_t wQA=1./fV0AQres;
  Double_t wQC=1./fV0CQres;
  Double_t wQT=1./fTPCQres;
  Double_t nCorr=1./PHOSMultiplicity(); //multiplicity vs. centrality
  
   FillHistogram("resV0A",fCentrality,fV0Ares,cWeight);
   FillHistogram("resV0C",fCentrality,fV0Cres,cWeight);
   FillHistogram("resTPC",fCentrality,fTPCres,cWeight);
   FillHistogram("qresV0A",fCentrality,fV0AQres,cWeight);
   FillHistogram("qresV0C",fCentrality,fV0CQres,cWeight);
   FillHistogram("qresTPC",fCentrality,fTPCQres,cWeight);
   
   TVector2 vQmA(fQV0A*TMath::Cos(fHarmonics*fRPQV0A),fQV0A*TMath::Sin(fHarmonics*fRPQV0A));
   TVector2 vQmB(fQV0C*TMath::Cos(fHarmonics*fRPQV0C),fQV0C*TMath::Sin(fHarmonics*fRPQV0C));
   
   FillHistogram("VaVbcen",fCentrality,vQmA*vQmB,cWeight);
   
  for (Int_t i1=0; i1 < fInPHOS; i1++) {
    AliCaloPhoton * ph1=(AliCaloPhoton*)fPHOSEvent->At(i1);
    
    Double_t dphiA=ph1->Phi()-fRPV0A;
    Double_t cosA=wA*TMath::Cos(fHarmonics*dphiA);
    Double_t qcosA=wQA*TMath::Cos(fHarmonics*(ph1->Phi()-fRPQV0A));
    
    Double_t dphiC=ph1->Phi()-fRPV0C;
    Double_t cosC=wC*TMath::Cos(fHarmonics*dphiC);
    Double_t qcosC=wQC*TMath::Cos(fHarmonics*(ph1->Phi()-fRPQV0C));

    Double_t dphiT=ph1->Phi()-fRP;
    Double_t cosT=wT*TMath::Cos(fHarmonics*dphiT);
    Double_t qcosT=wQT*TMath::Cos(fHarmonics*(ph1->Phi()-fRPQ));
         
    Double_t pt  = ph1->Pt();
    Double_t ptV = ph1->GetMomV2()->Pt();
    
    //Official SP
    //calculate vU
    TVector2 vU;
    Double_t dUX = TMath::Cos(fHarmonics*ph1->Phi());
    Double_t dUY = TMath::Sin(fHarmonics*ph1->Phi());
    vU.Set(dUX,dUY);

    Double_t dUQA = vU*vQmA;
    Double_t dUQB = vU*vQmB;

    Int_t pdg = 0, pdg_naive = 0;

    if (fMCArray) {
       pdg = ((AliAODMCParticle*)fMCArray->At(ph1->GetPrimaryAtVertex()))->GetPdgCode();
       pdg_naive = ((AliAODMCParticle*)fMCArray->At(ph1->GetPrimary()))->GetPdgCode();
    }


    std::vector<TString> passed_cuts = {"All", "Allcore"};

    if (ph1->IsCPVOK()) {
       passed_cuts.emplace_back("CPV");
       passed_cuts.emplace_back("CPVcore");
     }
    if (ph1->IsDisp2OK()) {
       passed_cuts.emplace_back("Disp2");
       passed_cuts.emplace_back("Disp2core");
     }
    if (ph1->IsCPVOK() && ph1->IsDisp2OK()) { 
       passed_cuts.emplace_back("Both2");
       passed_cuts.emplace_back("Both2core");
     }

    for (auto  cut : passed_cuts) {
     
      Double_t pT = cut.Contains("core") ?  ptV : pt;
      
      FillHistogram(Form("hPhotPt%s_TOF_cen%d", cut.Data(), fCenBin), pT, 1.e9*ph1->GetTime());

      if (!ph1->IsTOFOK()) continue;

      FillHistogram(Form("hPhotPt%s_cen%d", cut.Data(), fCenBin),pT,cWeight);
      FillHistogram(Form("hPhotPt%s_M%d_cen%d", cut.Data(), ph1->Module(),fCenBin),pT,cWeight);
      FillHistogram(Form("hPhotPtPdgAtVertex%s_cen%d", cut.Data(), fCenBin), pT, pdg, cWeight);
      FillHistogram(Form("hPhotPtPdg%s_cen%d", cut.Data(), fCenBin), pT, pdg_naive, cWeight);

      FillHistogram(Form("hPhotPhiSPV0A%s", cut.Data()),pT,fCentrality,dUQA,cWeight);	
      FillHistogram(Form("hPhotPhiSPV0C%s", cut.Data()),pT,fCentrality,dUQB,cWeight);	
      FillHistogram(Form("hPhotPhiR2V0A%s_cen%d", cut.Data(), fCenBin),pT,qcosA,fQV0A*cWeight*nCorr);
      FillHistogram(Form("hPhotPhiRV0A%s_cen%d", cut.Data(), fCenBin),pT,qcosA,fQV0A*cWeight);
      FillHistogram(Form("hPhotPhiR2V0C%s_cen%d", cut.Data(), fCenBin),pT,qcosC,fQV0C*cWeight*nCorr);
      FillHistogram(Form("hPhotPhiRV0C%s_cen%d", cut.Data(), fCenBin),pT,qcosC,fQV0C*cWeight);
      FillHistogram(Form("hPhotcosV0A%s_cen%d", cut.Data(), fCenBin),pT,cosA,cWeight);
      FillHistogram(Form("hPhotcosV0C%s_cen%d", cut.Data(), fCenBin),pT,cosC,cWeight);
      FillHistogram(Form("hPhotcosNV0A%s_cen%d", cut.Data(), fCenBin),pT,cosA,cWeight*nCorr);
      FillHistogram(Form("hPhotcosNV0C%s_cen%d", cut.Data(), fCenBin),pT,cosC,cWeight*nCorr);
      FillHistogram(Form("hPhotRa%s_cen%d", cut.Data(), fCenBin),pT,cWeight*wA);
      FillHistogram(Form("hPhotRc%s_cen%d", cut.Data(), fCenBin),pT,cWeight*wC);
      FillHistogram(Form("hPhotRNa%s_cen%d", cut.Data(), fCenBin),pT,cWeight*wA*nCorr);
      FillHistogram(Form("hPhotRNc%s_cen%d", cut.Data(), fCenBin),pT,cWeight*wC*nCorr);

      if (fHaveTPCRP) {
        FillHistogram(Form("hPhotRt%s_cen%d", cut.Data(), fCenBin),pT,cWeight*wT);
        FillHistogram(Form("hPhotRNt%s_cen%d", cut.Data(), fCenBin),pT,cWeight*wT*nCorr);
        FillHistogram(Form("hPhotPhiR2TPC%s_cen%d", cut.Data(), fCenBin),pT,qcosT,fQTPC*cWeight*nCorr);
        FillHistogram(Form("hPhotPhiRTPC%s_cen%d", cut.Data(), fCenBin), pT,qcosT,fQTPC*cWeight);
        FillHistogram(Form("hPhotcosTPC%s_cen%d", cut.Data(), fCenBin),pT,cosT,cWeight);
        FillHistogram(Form("hPhotcosNTPC%s_cen%d", cut.Data(), fCenBin),pT,cosT,cWeight*nCorr);
      }
    }
  }

  //tracks
  for (Int_t i1 = 0; i1 < fInPHOS - 1; i1 ++) {
    AliCaloPhoton * ph1=(AliCaloPhoton*)fPHOSEvent->At(i1);
    TestMatchingTrackPID(ph1, kFALSE);

    for(Int_t ev = 0; ev < prevPHOS->GetSize(); ev++){
      TClonesArray * mixPHOS = static_cast<TClonesArray*>(prevPHOS->At(ev));
        
      for (Int_t i2 = 0; i2 < fInPHOS - 1; i2 ++) {
         auto ph2 = (AliCaloPhoton*)mixPHOS->At(i2);
         TestMatchingTrackPID(ph2, kTRUE);
      }
    }
  }

  //pi0
  for (Int_t i1 = 0; i1 < fInPHOS -1; i1 ++) {
   
    AliCaloPhoton * ph1=(AliCaloPhoton*)fPHOSEvent->At(i1);
    Double_t pt1  = ph1->Pt();
    Double_t ptV1 = ph1->GetMomV2()->Pt();

    for (Int_t i2 = i1 + 1; i2 < fInPHOS; i2++) {
      AliCaloPhoton * ph2=(AliCaloPhoton*)fPHOSEvent->At(i2);
      Double_t pt2  = ph2->Pt();
      Double_t ptV2 = ph2->GetMomV2()->Pt();

      TLorentzVector p12   = *ph1  + *ph2;
      TLorentzVector p12V  = *ph1->GetMomV2()  + *ph2->GetMomV2();

      Bool_t cpvBitPi01 = ph1->IsCPVOK() &&  ph2->IsCPVOK();
      Bool_t dispBitPi01  = ph1->IsDisp2OK() &&  ph2->IsDisp2OK();

      std::vector<TString> passed_cuts_Pi0 = {"All", "Allcore"};

      if (cpvBitPi01) {
          passed_cuts_Pi0.emplace_back("CPV");
	  passed_cuts_Pi0.emplace_back("CPVcore");
      }

      if (dispBitPi01) {
          passed_cuts_Pi0.emplace_back("Disp2");
	  passed_cuts_Pi0.emplace_back("Disp2core");
      }

      if (cpvBitPi01 &&  dispBitPi01) {
          passed_cuts_Pi0.emplace_back("Both2");
	  passed_cuts_Pi0.emplace_back("Both2core");
      }


      for (auto cut : passed_cuts_Pi0) {
	 Double_t pT1 = cut.Contains("core") ?  ptV1 : pt1;
	 Double_t pT2 = cut.Contains("core") ?  ptV2 : pt2;
         Double_t pT   = cut.Contains("core") ?  p12V.Pt() : p12.Pt();
         Double_t mInv = cut.Contains("core") ?  p12V.M()  : p12.M();

         FillHistogram(Form("hPi0%sSingle_NoTOF_cen%d", cut.Data(), fCenBin), mInv, pT1);
    
         if (ph1->IsTOFOK()) 
            FillHistogram(Form("hPi0%sSingle_FirstTOF_cen%d", cut.Data(), fCenBin), mInv, pT1);

	 if (!ph1->IsTOFOK() || !ph2->IsTOFOK()) continue;
	
	 FillHistogram(Form("hPi0%sPt_cen%d", cut.Data(), fCenBin), mInv, pT);
	 FillHistogram(Form("hPi0%sSingle_cen%d", cut.Data(), fCenBin), mInv, pT1);
	 FillHistogram(Form("hPi0%sSingle_cen%d", cut.Data(), fCenBin), mInv, pT2);
      }
      
      if (!ph1->IsTOFOK() || !ph2->IsTOFOK()) continue;

      Bool_t cpvBitPi0First  = ph1->IsCPVOK();
      Bool_t dispBitPi0First  = ph1->IsDisp2OK();

      std::vector<TString> passed_cuts_Pi0First = {};
      if (cpvBitPi01)  {
         passed_cuts_Pi0First.emplace_back("CPV");
	 passed_cuts_Pi0First.emplace_back("CPVcore");
      } else {
         passed_cuts_Pi0First.emplace_back("AntiCPV");
	 passed_cuts_Pi0First.emplace_back("AntiCPVcore");
      }
      
      if (dispBitPi0First) {
         passed_cuts_Pi0First.emplace_back("Disp2");
	 passed_cuts_Pi0First.emplace_back("Disp2core");
      } else {
         passed_cuts_Pi0First.emplace_back("AntiDisp2");
	 passed_cuts_Pi0First.emplace_back("AntiDisp2core");
      }

      if (cpvBitPi0First && dispBitPi0First) {
         passed_cuts_Pi0First.emplace_back("Both2");
	 passed_cuts_Pi0First.emplace_back("Both2core");
      } else {
         passed_cuts_Pi0First.emplace_back("AntiBoth2");
	 passed_cuts_Pi0First.emplace_back("AntiBoth2core");
      }

      for (auto cut : passed_cuts_Pi0First) {
	 Double_t pT1  = cut.Contains("core") ?  ptV1 : pt1;
         Double_t mInv = cut.Contains("core") ?  p12V.M()  : p12.M();

	 FillHistogram(Form("hPi0%sFirst_cen%d", cut.Data(), fCenBin), mInv, pT1);
      }
    } // end of loop i2
  } // end of loop i1
  
  //now mixed

  for (Int_t i1 = 0; i1 < fInPHOS - 1; i1 ++) {
   
    AliCaloPhoton * ph1=(AliCaloPhoton*)fPHOSEvent->At(i1);
    Double_t pt1  = ph1->Pt();
    Double_t ptV1 = ph1->GetMomV2()->Pt();

    for(Int_t ev = 0; ev < prevPHOS->GetSize(); ev++){
      TClonesArray * mixPHOS = static_cast<TClonesArray*>(prevPHOS->At(ev));

      for(Int_t i2 = 0; i2 < mixPHOS->GetEntriesFast();i2++){
       //for (Int_t i2 = i1+1; i2 < inPHOS; i2++) {
        AliCaloPhoton * ph2=(AliCaloPhoton*)mixPHOS->At(i2);
        Double_t pt2  = ph2->Pt();
        Double_t ptV2 = ph2->GetMomV2()->Pt();

        TLorentzVector p12   = *ph1  + *ph2;
        TLorentzVector p12V  = *ph1->GetMomV2()  + *ph2->GetMomV2();

        Bool_t cpvBitPi01   = ph1->IsCPVOK() &&  ph2->IsCPVOK();
        Bool_t dispBitPi01  = ph1->IsDisp2OK() &&  ph2->IsDisp2OK();

        std::vector<TString> passed_cuts_Pi0 = {"All", "Allcore"};

        if (cpvBitPi01) {
            passed_cuts_Pi0.emplace_back("CPV");
            passed_cuts_Pi0.emplace_back("CPVcore");
        }

        if (dispBitPi01) {
            passed_cuts_Pi0.emplace_back("Disp2");
            passed_cuts_Pi0.emplace_back("Disp2core");
        }

        if (cpvBitPi01 &&  dispBitPi01) {
            passed_cuts_Pi0.emplace_back("Both2");
            passed_cuts_Pi0.emplace_back("Both2core");
        }


        for (auto cut : passed_cuts_Pi0) {
           Double_t pT1 = cut.Contains("core") ?  ptV1 : pt1;
           Double_t pT2 = cut.Contains("core") ?  ptV2 : pt2;
           Double_t pT   = cut.Contains("core") ?  p12V.Pt() : p12.Pt();
           Double_t mInv = cut.Contains("core") ?  p12V.M()  : p12.M();


           FillHistogram(Form("hMiPi0%sSingle_NoTOF_cen%d", cut.Data(), fCenBin), mInv, pT1);
    
           if (ph1->IsTOFOK()) 
              FillHistogram(Form("hMiPi0%sSingle_FirstTOF_cen%d", cut.Data(), fCenBin), mInv, pT1);

    	   if (!ph1->IsTOFOK() || !ph2->IsTOFOK()) continue;
	
           FillHistogram(Form("hMiPi0%sPt_cen%d", cut.Data(), fCenBin), mInv, pT);
           FillHistogram(Form("hMiPi0%sSingle_cen%d", cut.Data(), fCenBin), mInv, pT1);
           FillHistogram(Form("hMiPi0%sSingle_cen%d", cut.Data(), fCenBin), mInv, pT2);
        }
        
    	if (!ph1->IsTOFOK() || !ph2->IsTOFOK()) continue;

        Bool_t cpvBitPi0First  = ph1->IsCPVOK();
        Bool_t dispBitPi0First  = ph1->IsDisp2OK();

        std::vector<TString> passed_cuts_Pi0First = {};
        if (cpvBitPi01)  {
           passed_cuts_Pi0First.emplace_back("CPV");
           passed_cuts_Pi0First.emplace_back("CPVcore");
        } else {
           passed_cuts_Pi0First.emplace_back("AntiCPV");
           passed_cuts_Pi0First.emplace_back("AntiCPVcore");
        }
        
        if (dispBitPi0First) {
           passed_cuts_Pi0First.emplace_back("Disp2");
           passed_cuts_Pi0First.emplace_back("Disp2core");
        } else {
           passed_cuts_Pi0First.emplace_back("AntiDisp2");
           passed_cuts_Pi0First.emplace_back("AntiDisp2core");
        }

        if (cpvBitPi0First && dispBitPi0First) {
           passed_cuts_Pi0First.emplace_back("Both2");
           passed_cuts_Pi0First.emplace_back("Both2core");
        } else {
           passed_cuts_Pi0First.emplace_back("AntiBoth2");
           passed_cuts_Pi0First.emplace_back("AntiBoth2core");
        }

        for (auto cut : passed_cuts_Pi0First) {
           Double_t pT1 = cut.Contains("core") ?  ptV1 : pt1;
           Double_t mInv = cut.Contains("core") ?  p12V.M()  : p12.M();

           FillHistogram(Form("hMiPi0%sFirst_cen%d", cut.Data(), fCenBin), mInv, pT1);
        }
      } // end of loop i2
    } // prevPHOS loop
  } // end of loop i1

  //Now we either add current events to stack or remove
  //If no photons in current event - no need to add it to mixed
  //const Int_t kMixEvents[6]={5,5,5,10,10,30};
  if(fPHOSEvent->GetEntriesFast()>0){
    prevPHOS->AddFirst(fPHOSEvent);
    fPHOSEvent=0;
    //if(prevPHOS->GetSize()>kMixEvents[fCenBin]){//Remove redundant events
    if(prevPHOS->GetSize()>100){//Remove redundant events
      TClonesArray * tmp = static_cast<TClonesArray*>(prevPHOS->Last());
      prevPHOS->RemoveLast();
      delete tmp;
    }
  }
  // Post output data.
  fEventCounter++;
}

//________________________________________________________________________
void AliAnalysisTaskGammaPHOSPbPbRun2::Terminate(Option_t *)
{
  // Draw result to the screen
  // Called once at the end of the query
  
}

//_____________________________________________________________________________
void AliAnalysisTaskGammaPHOSPbPbRun2::FillHistogram(const char * key,Double_t x)const
{
  //FillHistogram
  // TObject * tmp = (fOutputContainer->FindObject(key))?(fOutputContainer->FindObject(key)):(fOutputContainer2->FindObject(key));

  TObject *tmp=0x0;
  TObject *tmp1 = fOutputContainer->FindObject(key);
  TObject *tmp2 = fOutputContainer2->FindObject(key);
  TObject *tmp3 = fOutputContainer3->FindObject(key);
  TObject *tmp4 = fOutputContainer4->FindObject(key);
  TObject *tmp5 = fOutputContainer5->FindObject(key);
  TObject *tmp6 = fOutputContainer6->FindObject(key);

  if (tmp1) tmp = tmp1;
  else if (tmp2) tmp = tmp2;
  else if (tmp3) tmp = tmp3;
  else if (tmp4) tmp = tmp4;
  else if (tmp5) tmp = tmp5;
  else if (tmp6) tmp = tmp6;

  if(!tmp){
    AliInfo(Form("can not find histogram <%s> ",key));
    return;
  }
  if(tmp->IsA() == TClass::GetClass("TH1I")){
    ((TH1I*)tmp)->Fill(x);
    return;
  }
  if(tmp->IsA() == TClass::GetClass("TH1F")){
    ((TH1F*)tmp)->Fill(x);
    return;
  }
  if(tmp->IsA() == TClass::GetClass("TProfile")){
    ((TH1F*)tmp)->Fill(x,1.);
    return;
  }
  if(tmp->IsA() == TClass::GetClass("TH1D")){
    ((TH1D*)tmp)->Fill(x);
    return;
  }  
  AliInfo(Form("can not find 1D histogram <%s> ",key));
}
//_____________________________________________________________________________
void AliAnalysisTaskGammaPHOSPbPbRun2::FillHistogram(const char * key,Double_t x,Double_t y)const
{
  //FillHistogram
  //TObject * tmp = (fOutputContainer->FindObject(key))?(fOutputContainer->FindObject(key)):(fOutputContainer2->FindObject(key));
  //TObject * tmp = fOutputContainer->FindObject(key);

  TObject *tmp=0x0;
  TObject *tmp1 = fOutputContainer->FindObject(key);
  TObject *tmp2 = fOutputContainer2->FindObject(key);
  TObject *tmp3 = fOutputContainer3->FindObject(key);
  TObject *tmp4 = fOutputContainer4->FindObject(key);
  TObject *tmp5 = fOutputContainer5->FindObject(key);
  TObject *tmp6 = fOutputContainer6->FindObject(key);

  if (tmp1) tmp = tmp1;
  else if (tmp2) tmp = tmp2;
  else if (tmp3) tmp = tmp3;
  else if (tmp4) tmp = tmp4;
  else if (tmp5) tmp = tmp5;
  else if (tmp6) tmp = tmp6;

  if(!tmp){
    AliInfo(Form("can not find histogram <%s> ",key));
    return;
  }
  if(tmp->IsA() == TClass::GetClass("TH1F")){
    ((TH1F*)tmp)->Fill(x,y);
    return;
  }
  if(tmp->IsA() == TClass::GetClass("TProfile")){
    ((TProfile*)tmp)->Fill(x,y);
    return;
  }
  if(tmp->IsA() == TClass::GetClass("TH2F")){
    ((TH2F*)tmp)->Fill(x,y);
    return;
  }
  AliError(Form("Calling FillHistogram with 2 parameters for histo <%s> of type %s",key,tmp->IsA()->GetName()));
}

//_____________________________________________________________________________
void AliAnalysisTaskGammaPHOSPbPbRun2::FillHistogram(const char * key,Double_t x,Double_t y, Double_t z) const
{

  //TObject * tmp = (fOutputContainer->FindObject(key))?(fOutputContainer->FindObject(key)):(fOutputContainer2->FindObject(key));

  TObject *tmp=0x0;
  TObject *tmp1 = fOutputContainer->FindObject(key);
  TObject *tmp2 = fOutputContainer2->FindObject(key);
  TObject *tmp3 = fOutputContainer3->FindObject(key);
  TObject *tmp4 = fOutputContainer4->FindObject(key);
  TObject *tmp5 = fOutputContainer5->FindObject(key);
  TObject *tmp6 = fOutputContainer6->FindObject(key);

  if (tmp1) tmp = tmp1;
  else if (tmp2) tmp = tmp2;
  else if (tmp3) tmp = tmp3;
  else if (tmp4) tmp = tmp4;
  else if (tmp5) tmp = tmp5;
  else if (tmp6) tmp = tmp6;

  if(!tmp){
    AliInfo(Form("can not find histogram <%s> ",key));
    return;
  }
  if(tmp->IsA() == TClass::GetClass("TH2F")){
    ((TH2F*)tmp)->Fill(x,y,z);
    return;
  }
  if(tmp->IsA() == TClass::GetClass("TProfile")){
    ((TProfile*)tmp)->Fill(x,y,z);
    return;
  }
  if(tmp->IsA() == TClass::GetClass("TProfile2D")){
    ((TProfile2D*)tmp)->Fill(x,y,z);
    return;
  }
  if(tmp->IsA() == TClass::GetClass("TH3F")){
    ((TH3F*)tmp)->Fill(x,y,z);
    return;
  }
  AliInfo(Form("can not find 2D/3D histogram <%s> ",key));
}
//_____________________________________________________________________________
void AliAnalysisTaskGammaPHOSPbPbRun2::FillHistogram(const char * key,Double_t x,Double_t y, Double_t z, Double_t w) const
{

  //TObject * tmp = (fOutputContainer->FindObject(key))?(fOutputContainer->FindObject(key)):(fOutputContainer2->FindObject(key));

  TObject *tmp=0x0;
  TObject *tmp1 = fOutputContainer->FindObject(key);
  TObject *tmp2 = fOutputContainer2->FindObject(key);
  TObject *tmp3 = fOutputContainer3->FindObject(key);
  TObject *tmp4 = fOutputContainer4->FindObject(key);
  TObject *tmp5 = fOutputContainer5->FindObject(key);
  TObject *tmp6 = fOutputContainer6->FindObject(key);

  if (tmp1) tmp = tmp1;
  else if (tmp2) tmp = tmp2;
  else if (tmp3) tmp = tmp3;
  else if (tmp4) tmp = tmp4;
  else if (tmp5) tmp = tmp5;
  else if (tmp6) tmp = tmp6;

  if(!tmp){
    AliInfo(Form("can not find histogram <%s> ",key));
    return;
  }
  if(tmp->IsA() == TClass::GetClass("TH3F")){
    ((TH3F*)tmp)->Fill(x,y,z,w);
    return;
  }
  if(tmp->IsA() == TClass::GetClass("TProfile2D")){
    ((TProfile2D*)tmp)->Fill(x,y,z,w);
    return;
  }
  AliInfo(Form("can not find 3D histogram <%s> ",key));
}

//___________________________________________________________________________
Int_t AliAnalysisTaskGammaPHOSPbPbRun2::ConvertRunNumber(Int_t run){

  switch(run){	
    case 295584: return 1;
    case 295585: return 2;
    case 295586: return 3;
    case 295587: return 4;
    case 295588: return 5;
    case 295589: return 6;
    case 295610: return 7;
    case 295611: return 8;
    case 295612: return 9;
    case 295615: return 10;
    case 295666: return 11;
    case 295667: return 12;
    case 295668: return 13;
    case 295673: return 14;
    case 295675: return 15;
    case 295676: return 16;
    case 295677: return 17;
    case 295712: return 18;
    case 295717: return 19;
    case 295718: return 20;
    case 295719: return 21;
    case 295721: return 22;
    case 295723: return 23;
    case 295725: return 24;
    case 295753: return 25;
    case 295754: return 26;
    case 295755: return 27;
    case 295756: return 28;
    case 295758: return 29;
    case 295759: return 30;
    case 295762: return 31;
    case 295763: return 32;
    case 295788: return 33;
    case 295791: return 34;
    case 295818: return 35;
    case 295819: return 36;
    case 295822: return 37;
    case 295825: return 38;
    case 295826: return 39;
    case 295829: return 40;
    case 295831: return 41;
    case 295853: return 42;
    case 295854: return 43;
    case 295855: return 44;
    case 295856: return 45;
    case 295859: return 46;
    case 295860: return 47;
    case 295861: return 48;
    case 295881: return 49;
    case 295908: return 50;
    case 295910: return 51;
    case 295913: return 52;
    case 295936: return 53;
    case 295937: return 54;
    case 295941: return 55;
    case 295942: return 56;
    case 296016: return 57;
    case 296060: return 58;
    case 296061: return 59;
    case 296062: return 60;
    case 296063: return 61;
    case 296065: return 62;
    case 296066: return 63;
    case 296068: return 64;
    case 296074: return 65;
    case 296123: return 66;
    case 296132: return 67;
    case 296133: return 68;
    case 296134: return 69;
    case 296135: return 70;
    case 296142: return 71;
    case 296143: return 72;
    case 296191: return 73;
    case 296192: return 74;
    case 296194: return 75;
    case 296195: return 76;
    case 296196: return 77;
    case 296197: return 78;
    case 296198: return 79;
    case 296240: return 80;
    case 296241: return 81;
    case 296242: return 82;
    case 296243: return 83;
    case 296244: return 84;
    case 296246: return 85;
    case 296247: return 86;
    case 296269: return 87;
    case 296270: return 88;
    case 296273: return 89;
    case 296279: return 90;
    case 296280: return 91;
    case 296303: return 92;
    case 296304: return 93;
    case 296309: return 94;
    case 296312: return 95;
    case 296376: return 96;
    case 296377: return 97;
    case 296378: return 98;
    case 296379: return 99;
    case 296380: return 100;
    case 296381: return 101;
    case 296383: return 102;
    case 296414: return 103;
    case 296415: return 104;
    case 296419: return 105;
    case 296420: return 106;
    case 296423: return 107;
    case 296424: return 108;
    case 296433: return 109;
    case 296472: return 110;
    case 296509: return 111;
    case 296510: return 112;
    case 296511: return 113;
    case 296512: return 114;
    case 296514: return 115;
    case 296516: return 116;
    case 296547: return 117;
    case 296548: return 118;
    case 296549: return 119;
    case 296550: return 120;
    case 296551: return 121;
    case 296552: return 122;
    case 296553: return 123;
    case 296594: return 124;
    case 296615: return 125;
    case 296616: return 126;
    case 296618: return 127;
    case 296619: return 128;
    case 296621: return 129;
    case 296622: return 130;
    case 296623: return 131;
    case 296690: return 132;
    case 296691: return 133;
    case 296693: return 134;
    case 296694: return 135;
    case 296749: return 136;
    case 296750: return 137;
    case 296752: return 138;
    case 296781: return 139;
    case 296784: return 140;
    case 296785: return 141;
    case 296786: return 142;
    case 296787: return 143;
    case 296790: return 144;
    case 296794: return 145;
    case 296799: return 146;
    case 296835: return 147;
    case 296836: return 148;
    case 296848: return 149;
    case 296849: return 150;
    case 296850: return 151;
    case 296851: return 152;
    case 296852: return 153;
    case 296890: return 154;
    case 296894: return 155;
    case 296899: return 156;
    case 296900: return 157;
    case 296903: return 158;
    case 296930: return 159;
    case 296931: return 160;
    case 296932: return 161;
    case 296934: return 162;
    case 296935: return 163;
    case 296938: return 164;
    case 296941: return 165;
    case 296966: return 166;
    case 297029: return 167;
    case 297035: return 168;
    case 297085: return 169;
    case 297117: return 170;
    case 297118: return 171;
    case 297119: return 172;
    case 297123: return 173;
    case 297124: return 174;
    case 297128: return 175;
    case 297129: return 176;
    case 297132: return 177;
    case 297133: return 178;
    case 297193: return 179;
    case 297194: return 180;
    case 297196: return 181;
    case 297218: return 182;
    case 297219: return 183;
    case 297221: return 184;
    case 297222: return 185;
    case 297310: return 186;
    case 297311: return 187;
    case 297312: return 188;
    case 297315: return 189;

    default : return 199;
  } 

}

//____________________________________________________________________________
void  AliAnalysisTaskGammaPHOSPbPbRun2::EvalV0ReactionPlane(AliAODEvent * event){

  AliEventplane *eventPlane = event->GetEventplane();
  if( ! eventPlane ) { AliError("Event has no event plane"); return; }
  Double_t qx = 0, qy = 0;
  //V0A
  fRPV0A = eventPlane->CalculateVZEROEventPlane(event,8, fHarmonics, qx, qy);
  fQV0A=TMath::Sqrt(qx*qx+qy*qy);
  //V0C
  fRPV0C = eventPlane->CalculateVZEROEventPlane(event,9, fHarmonics, qx, qy);
  fQV0C=TMath::Sqrt(qx*qx+qy*qy);  

  
  while(fRPV0A<0)fRPV0A+=TMath::TwoPi()/fHarmonics;
  while(fRPV0A>TMath::TwoPi()/fHarmonics)fRPV0A-=TMath::TwoPi()/fHarmonics;
//  fRPV0A=ApplyFlatteningV0A(fRPV0A,fCentrality);
  while(fRPV0A<0)fRPV0A+=TMath::TwoPi()/fHarmonics;
  while(fRPV0A>TMath::TwoPi()/fHarmonics)fRPV0A-=TMath::TwoPi()/fHarmonics;
  FillHistogram("phiRPV0A",fRPV0A,fCentrality);
  fRPV0A = fV0AFlat->MakeFlat(fRPV0A,fCentrality);
  
  while(fRPV0C<0)fRPV0C+=TMath::TwoPi()/fHarmonics;
  while(fRPV0C>TMath::TwoPi()/fHarmonics)fRPV0C-=TMath::TwoPi()/fHarmonics;
//  fRPV0C=ApplyFlatteningV0C(fRPV0C,fCentrality);
  while(fRPV0C<0)fRPV0C+=TMath::TwoPi()/fHarmonics;
  while(fRPV0C>TMath::TwoPi()/fHarmonics)fRPV0C-=TMath::TwoPi()/fHarmonics;
  FillHistogram("phiRPV0C",fRPV0C,fCentrality);
  fRPV0C = fV0CFlat->MakeFlat(fRPV0C,fCentrality);
 
  //So far no difference between fRPV0A and fRPV0AQ
  FillHistogram("phiRPV0AQ",fRPV0A,fCentrality,fQV0A);
  FillHistogram("phiRPV0CQ",fRPV0C,fCentrality,fQV0C);
  
}

//____________________________________________________________________________
Bool_t AliAnalysisTaskGammaPHOSPbPbRun2::GetTPCEventPlane(Double_t &epAngle, Double_t &qsubRes)
{

    AliEventplane* ep=new AliEventplane();

    float mQx=0, mQy=0;
    float mQx1=0, mQy1=0, mQx2=0, mQy2=0;
    AliAODTrack* track;
    Double_t weight;
    Int_t idtemp = -1;
    int trackcounter1=0, trackcounter2=0;

    Int_t maxID=0;

    TObjArray *tracklist=GetEventPlaneTracks(maxID);
    ep->GetQContributionXArray()->Set(maxID);
    ep->GetQContributionYArray()->Set(maxID);
    ep->GetQContributionXArraysub1()->Set(maxID);
    ep->GetQContributionYArraysub1()->Set(maxID);
    ep->GetQContributionXArraysub2()->Set(maxID);
    ep->GetQContributionYArraysub2()->Set(maxID);

    for (int i=0; i<maxID; i++){
        weight = 1;
        track = dynamic_cast<AliAODTrack*> (tracklist->At(i));
        if (track) {

            weight=GetWeight(track);
            idtemp = i; // track->GetID();
            // TPC only tracks have negative id ((-1)*IDESD - 1) in AOD
            //if (fIsAOD && (fUseTPCOnlyTracks)) idtemp = idtemp*(-1) - 1;

            Double_t qx=weight*cos(Double_t(fHarmonics)*track->Phi());
            Double_t qy=weight*sin(Double_t(fHarmonics)*track->Phi());
            ep->GetQContributionXArray()->AddAt(qx,idtemp);
            ep->GetQContributionYArray()->AddAt(qy,idtemp);

            mQx += (qx);
            mQy += (qy);

            // This loop splits the track set into 2 random subsets
            if( trackcounter1 < int(maxID/2.) && trackcounter2 < int(maxID/2.)){
                float random = gRandom->Rndm();
                if(random < .5){
                    mQx1 += (qx);
                    mQy1 += (qy);
                    ep->GetQContributionXArraysub1()->AddAt(qx,idtemp);
                    ep->GetQContributionYArraysub1()->AddAt(qy,idtemp);
                    trackcounter1++;
                    }
                else {
                    mQx2 += (qx);
                    mQy2 += (qy);
                    ep->GetQContributionXArraysub2()->AddAt(qx,idtemp);
                    ep->GetQContributionYArraysub2()->AddAt(qy,idtemp);
                    trackcounter2++;
                }
            }
            else{
                if( trackcounter1 >= int(maxID/2.)){
                    mQx2 += (qx);
                    mQy2 += (qy);
                    ep->GetQContributionXArraysub2()->AddAt(qx,idtemp);
                    ep->GetQContributionYArraysub2()->AddAt(qy,idtemp);
                    trackcounter2++;
                }
                else {
                    mQx1 += (qx);
                    mQy1 += (qy);
                    ep->GetQContributionXArraysub1()->AddAt(qx,idtemp);
                    ep->GetQContributionYArraysub1()->AddAt(qy,idtemp);
                    trackcounter1++;
                }
            }
        }
    }

    tracklist->Clear();
    delete tracklist;
    tracklist = NULL;

    TVector2 *mQ=new TVector2();
    mQ->Set(mQx,mQy);
    epAngle=mQ->Phi()/Double_t(fHarmonics);
    
    fQTPC = mQ->Mod();
   
    TVector2 *qsub1=new TVector2();
    TVector2 *qsub2=new TVector2();
    qsub1->Set(mQx1,mQy1);
    qsub2->Set(mQx2,mQy2);

    ep->SetQVector(mQ);
    ep->SetEventplaneQ(epAngle);
    ep->SetQsub(qsub1,qsub2);
    qsubRes=qsub1->Phi()/Double_t(fHarmonics) - qsub2->Phi()/Double_t(fHarmonics);

    Int_t ntracks=trackcounter1+trackcounter2;

    delete ep;
    
    if(ntracks<3)return kFALSE;// <3 -> no subevents
    return kTRUE;
}
//_________________________________________________________________________
TObjArray* AliAnalysisTaskGammaPHOSPbPbRun2::GetEventPlaneTracks(Int_t &maxID)
{
  
    AliAODEvent * event = (AliAODEvent*)InputEvent();
        
    TObjArray *tracklist1=new TObjArray();
    Int_t nt = event->GetNumberOfTracks();
    for (Int_t i=0; i<nt; i++) {
      AliAODTrack *aodTrack=(AliAODTrack*)event->GetTrack(i);
      if(!aodTrack) continue;

      if(!aodTrack->IsHybridGlobalConstrainedGlobal())
        continue;    
      if( aodTrack->Pt()<0.15 ||  aodTrack->Pt()>20.)
	continue;
      if( TMath::Abs(aodTrack->Eta())<0.5 || TMath::Abs(aodTrack->Eta())>0.8)
	continue;
      tracklist1->Add(new AliAODTrack(*(aodTrack)));
    }  
    tracklist1->SetOwner(kTRUE);
    maxID=tracklist1->GetEntries();
    if(!tracklist1)AliError("No tracklist");
    return tracklist1;
}
//_________________________________________________________________________
Double_t AliAnalysisTaskGammaPHOSPbPbRun2::GetWeight(TObject* track1)
{
    Double_t ptweight=1;
    AliVTrack* track = dynamic_cast<AliVTrack*>(track1);
    if (track) {
        if (track->Pt()<2) ptweight=track->Pt();
        else ptweight=2;
    }

    return ptweight*GetPhiWeight(track);
}
//_________________________________________________________________________
Double_t AliAnalysisTaskGammaPHOSPbPbRun2::GetPhiWeight(TObject* track1)
{
  Double_t phiweight=1;
  AliVTrack* track = dynamic_cast<AliVTrack*>(track1);

  if (fPhiDist && track) {
      Double_t nParticles = fPhiDist->Integral();
      Double_t nPhibins = fPhiDist->GetNbinsX();

      Double_t Phi = track->Phi();

      while (Phi<0) Phi += TMath::TwoPi();
      while (Phi>TMath::TwoPi()) Phi -= TMath::TwoPi();

      Double_t PhiDistValue = fPhiDist->GetBinContent(1+TMath::FloorNint((track->Phi())*nPhibins/TMath::TwoPi()));

      if (PhiDistValue > 0) phiweight = nParticles/nPhibins/PhiDistValue;
  }
  return phiweight;
}
//_________________________________________________________________________
void AliAnalysisTaskGammaPHOSPbPbRun2::EvalResolution()
{
  //Estimate resolutions for each of EP estimators
  
  Double_t x= fCentrality;
  if(x<0.5)x=0.5;
  if(fHarmonics==2){  
    fV0Cres=3.294365e-01*(1.+x*1.507363e-01+x*x*1.272448e-03+x*x*x*4.229882e-05-x*x*x*x*6.939672e-07)/(1.-x*8.165966e-04+x*x*3.989182e-03-x*x*x*7.727254e-05+x*x*x*x*9.914099e-07);
    fV0Ares=2.905439e-01*(1.+x*1.372979e+00+x*x*6.166033e-01-x*x*x*7.905871e-03+x*x*x*x*1.507491e-06)/(1.+x*1.902084e+00+x*x*1.735136e-01-x*x*x*1.689358e-03+x*x*x*x*1.364915e-05);
    fTPCres=3.808781e-01*(1.+x*2.386826e-01+x*x*5.192646e-02+x*x*x*4.098382e-04-x*x*x*x*1.453986e-05)/(1.+x*1.259547e-01+x*x*2.285279e-02+x*x*x*1.030309e-04-x*x*x*x*2.666835e-06);
    
    //LHC01h    
//    fV0Cres=3.325009e-01*(1+x*1.939315e-01-x*x*9.652051e-04+x*x*x*1.720295e-05-x*x*x*x*3.731234e-07 )/(1+x*2.802502e-02+x*x*1.874919e-03-x*x*x*3.980602e-05+x*x*x*x*4.447131e-07);
//    fV0Ares=2.403080e-01*(1+x*2.078031e-01+x*x*3.030761e-04+x*x*x*6.242597e-05+x*x*x*x*-1.072181e-06)/(1+x*9.135497e-03+x*x*3.208661e-03-x*x*x*5.901547e-05+x*x*x*x*8.083149e-07);
//    fTPCres=3.907111e-01*(1+x*2.036863e-01+x*x*7.435000e-04+x*x*x*6.967697e-05+x*x*x*x*-1.132069e-06)/(1+x*2.751258e-02+x*x*4.000528e-03-x*x*x*6.829516e-05+x*x*x*x*7.677720e-07);
  }
  if(fHarmonics==3){
    fV0Cres=2.836819e-01*(1+x*3.017023e-01-x*x*1.083738e-02+x*x*x*9.760039e-05)/(1.+x*2.490451e-01-x*x*7.528651e-03+x*x*x*6.073070e-05);
    fV0Ares=2.002936e-01*(1+x*2.258327e-01+x*x*3.836017e-04-x*x*x*6.525109e-05)/(1.+x*1.739020e-01+x*x*1.889962e-03-x*x*x*2.786371e-05);
    fTPCres=4.357285e-01*(1+x*1.959005e-01+x*x*2.670202e-03-x*x*x*7.921636e-05)/(1.+x*1.582993e-01+x*x*2.707513e-03-x*x*x*2.773870e-05);

    //LHC10h
//    fV0Cres=TMath::Max(0.02,2.687135e-01*(1.+x*8.117353e-01+x*x*2.271778e-03-x*x*x*1.798932e-04)/(1.+x*6.721703e-01+x*x*3.829202e-03+x*x*x*1.202772e-05));
//    fV0Ares=TMath::Max(0.02,2.030253e-01*(1.+x*2.118380e-01+x*x*2.894326e-02-x*x*x*4.375165e-04)/(1.+x*1.908986e-01+x*x*2.324712e-02-x*x*x*7.957708e-05));
//    fTPCres=TMath::Max(0.02,4.829372e-01*(1.+x*1.300387e-01+x*x*1.579848e-04+x*x*x*3.856068e-05-x*x*x*x*8.861878e-07)/(1.+x*1.013947e-01+x*x*1.112352e-03+x*x*x*1.766346e-05+x*x*x*x*2.207416e-07));
  }
}
//_________________________________________________________________________
void AliAnalysisTaskGammaPHOSPbPbRun2::EvalQResolution()
{
  //Estimate resolutions for each of EP estimators
  
  Double_t x= fCentrality;
  if(x<0.5)x=0.5;
  if(fHarmonics==2){
    fTPCQres=4.857694e-01*(1+x*2.920904e-01-x*x*6.514101e-03+x*x*x*1.018369e-04-x*x*x*x*6.986208e-07)/(1+x*9.235900e-02+x*x*6.325890e-05-x*x*x*3.163401e-05+x*x*x*x*4.747334e-07);
    fV0CQres=4.387220e-01*(1+x*1.357633e-01+x*x*1.995734e-02-x*x*x*3.375979e-04+x*x*x*x*1.024155e-06)/(1+x*2.594425e-02+x*x*1.252949e-02-x*x*x*2.449652e-04+x*x*x*x*1.850385e-06);
    fV0AQres=3.242703e-01*(1+x*2.721827e-01+x*x*9.379002e-02-x*x*x*2.090242e-03+x*x*x*x*1.154136e-05)/(1+x*1.793834e-01+x*x*3.256506e-02-x*x*x*7.118371e-04+x*x*x*x*5.327807e-06);
  }
  if(fHarmonics==3){
    if(x>70)x=70;//bad statistics above
    fTPCQres=6.054041e-01*(1+x*2.871979e-01-x*x*8.501828e-03+x*x*x*6.042861e-05)/(1+x*2.422995e-01-x*x*6.217416e-03+x*x*x*4.229649e-05);
    fV0CQres=3.666094e-01*(1+x*6.818115e+00+x*x*2.199777e-01-x*x*x*3.688469e-03)/(1+x*6.666333e+00+x*x*1.403492e-01+x*x*x*1.229574e-03);
    fV0AQres=1.927797e-01*(1+x*1.306917e+01-x*x*1.261409e-01-x*x*x*2.002856e-04)/(1+x*9.054694e+00-x*x*1.144981e-01+x*x*x*2.580067e-03);
  }
}
//_________________________________________________________________________
void AliAnalysisTaskGammaPHOSPbPbRun2::ApplyFinalFlattening()
{//apply final fine flattening
  
  Int_t ibin = fTPCfinalC2->FindBin(fCentrality);
  
  //TPC
  Double_t v2c =2.*fTPCfinalC2->GetBinContent(ibin)/fHarmonics;
  Double_t v2s =2.*fTPCfinalS2->GetBinContent(ibin)/fHarmonics; 
  Double_t v3c =2.*fTPCfinalC4->GetBinContent(ibin)/fHarmonics/2.;
  Double_t v3s =2.*fTPCfinalS4->GetBinContent(ibin)/fHarmonics/2.; 
  fRP = fRP +
  v2c*TMath::Sin(fHarmonics*fRP)-v2s*TMath::Cos(fHarmonics*fRP)+
  v3c*TMath::Sin(2.*fHarmonics*fRP)-v3s*TMath::Cos(2.*fHarmonics*fRP);
  
  while(fRP<0)fRP+=TMath::TwoPi()/fHarmonics;
  while(fRP>TMath::TwoPi()/fHarmonics)fRP-=(TMath::TwoPi()/fHarmonics);

  //V0A
  v2c =2.*fV0AfinalC2->GetBinContent(ibin)/fHarmonics;
  v2s =2.*fV0AfinalS2->GetBinContent(ibin)/fHarmonics; 
  v3c =2.*fV0AfinalC4->GetBinContent(ibin)/fHarmonics/2.;
  v3s =2.*fV0AfinalS4->GetBinContent(ibin)/fHarmonics/2.; 
  fRPV0A = fRPV0A +
  v2c*TMath::Sin(fHarmonics*fRPV0A)-v2s*TMath::Cos(fHarmonics*fRPV0A)+
  v3c*TMath::Sin(2.*fHarmonics*fRPV0A)-v3s*TMath::Cos(2.*fHarmonics*fRPV0A);
 
  while(fRPV0A<0)fRPV0A+=TMath::TwoPi()/fHarmonics;
  while(fRPV0A>TMath::TwoPi()/fHarmonics)fRPV0A-=(TMath::TwoPi()/fHarmonics);
  
  //V0C
  v2c =2.*fV0CfinalC2->GetBinContent(ibin)/fHarmonics;
  v2s =2.*fV0CfinalS2->GetBinContent(ibin)/fHarmonics; 
  v3c =2.*fV0CfinalC4->GetBinContent(ibin)/fHarmonics/2.;
  v3s =2.*fV0CfinalS4->GetBinContent(ibin)/fHarmonics/2.; 
  fRPV0C = fRPV0C +
  v2c*TMath::Sin(fHarmonics*fRPV0C)-v2s*TMath::Cos(fHarmonics*fRPV0C)+
  v3c*TMath::Sin(2.*fHarmonics*fRPV0C)-v3s*TMath::Cos(2.*fHarmonics*fRPV0C);
  
  while(fRPV0C<0)fRPV0C+=TMath::TwoPi()/fHarmonics;
  while(fRPV0C>TMath::TwoPi()/fHarmonics)fRPV0C-=(TMath::TwoPi()/fHarmonics);
}  
//_________________________________________________________________________
void AliAnalysisTaskGammaPHOSPbPbRun2::ApplyFinalQFlattening()
{
  //apply final fine flattening
  //Apply Q-flattening
  fRPQ    = fTPCQFlat->MakeFlat(fRP,fCentrality);
  fRPQV0A = fV0AQFlat->MakeFlat(fRPV0A,fCentrality);
  fRPQV0C = fV0CQFlat->MakeFlat(fRPV0C,fCentrality);
    
  while(fRPQ<0)fRPQ+=TMath::TwoPi()/fHarmonics;
  while(fRPQ>TMath::TwoPi()/fHarmonics)fRPQ-=(TMath::TwoPi()/fHarmonics);
  
  while(fRPQV0A<0)fRPQV0A+=TMath::TwoPi()/fHarmonics;
  while(fRPQV0A>TMath::TwoPi()/fHarmonics)fRPQV0A-=(TMath::TwoPi()/fHarmonics);
  
  while(fRPQV0C<0)fRPQV0C+=TMath::TwoPi()/fHarmonics;
  while(fRPQV0C>TMath::TwoPi()/fHarmonics)fRPQV0C-=(TMath::TwoPi()/fHarmonics);
}  
//_________________________________________________________________________
Double_t AliAnalysisTaskGammaPHOSPbPbRun2::PHOSMultiplicity()
{
   Double_t x=fCentrality;
   return TMath::Max(2.,3.65177e+01-x*1.56822+x*x*3.06817e-02-x*x*x*3.14334e-04+x*x*x*x*1.27239e-06);
}
//_________________________________________________________________________
Double_t AliAnalysisTaskGammaPHOSPbPbRun2::CentralityWeight(Double_t c)
{
  //Weight to make flat centrality distribution
  Double_t weight=1.;
  return weight;//!!!
  //Central
  if(c<10.)
    weight = (4.81061e+05-1.18325e+04*c+2.33302e+04*c*c-8.46768e+03*c*c*c+1.16784e+03*c*c*c*c-5.60438e+01*c*c*c*c*c)/500552.;

  //SemiCentral
  if(c>14. && c <50.) //flat region
    weight = 1.;
  else if(c>10. && c<=14.)
    weight = 1.+7.56198e-03*TMath::Power(c-14.,2)+1.25110e-02*TMath::Power(c-14.,3); 
    else if(c>50. && c<=56.)
       weight = 1.-2.50575e-02*TMath::Power(c-50.,2)-1.29059e-03*TMath::Power(c-50.,3); 
     
  if(weight>0.01)
    return 1./weight;
  else
    return 0.;
}

//_____________________________________________________________________________
Int_t AliAnalysisTaskGammaPHOSPbPbRun2::GetPrimaryLabelAtVertex(AliVCluster *clu) //Returns label at vertex
{
   if (!fMCArray) 
     return 0;
      
   Int_t iPrimaryAtVertex = clu->GetLabel();
   AliAODMCParticle *particle0 =  (AliAODMCParticle*) fMCArray->At(iPrimaryAtVertex);

   while (getR(particle0) > 1.0) {
      iPrimaryAtVertex = particle0->GetMother();
      particle0 = (AliAODMCParticle*) fMCArray->At(particle0->GetMother());
   }

   Int_t nn = iPrimaryAtVertex;

   if (!particle0) return 0;

   if (particle0->GetPdgCode() == 22 || particle0->GetPdgCode() == 11) {
     //iPrimaryAtVertex = particle0->GetMother();
     for (Int_t i = 0; i < nn; i++) {
       AliAODMCParticle* particle =  (AliAODMCParticle*) fMCArray->At(i);
       if(particle->GetPdgCode() != 310 && particle->GetPdgCode() != 130) 
         continue;
       Int_t iSecondDaughter = particle->GetDaughterLabel(1); 
       if(iSecondDaughter != iPrimaryAtVertex) 
         continue;
       else
         iPrimaryAtVertex = i;
     }
   } else 
       iPrimaryAtVertex = nn;
  
   return iPrimaryAtVertex;   
}

//_____________________________________________________________________________
Int_t AliAnalysisTaskGammaPHOSPbPbRun2::GetPrimaryLabel(AliVCluster *clu) //Returns label of the impinging particle 
{
   if (!fMCArray) 
     return 0;
      
   return clu->GetLabel();
}

//__________________________________________________________________________
void AliAnalysisTaskGammaPHOSPbPbRun2::ProcessMC()
{
   if (!fMCArray)
     return;

   for (Int_t iPart  = 0; iPart < fMCArray->GetEntriesFast(); iPart++) {
     AliAODMCParticle* particle =  (AliAODMCParticle*) fMCArray->At(iPart);

     if (getR(particle) > 1.) 
	continue;
     if (TMath::Abs(particle->Pt() < 1.e-6))
        continue;

     Double_t w = 1; // primary particle weight
     Double_t pt = particle->Pt();
     Double_t y = particle->Y();
     Int_t pdg = particle->GetPdgCode();

     Double_t vz = TMath::Abs(particle->Zv() - fVtx5[2]);
     Double_t vx = particle->Xv();
     Double_t vy = particle->Yv();
     Double_t vr = TMath::Hypot(vx, vy);

     switch(pdg) {
     case kGamma: 
       FillHistogram(Form("hPrimPhot_cen%d", fCenBin), pt, y, w);
       FillHistogram(Form("hMC_gamma_vertex_cen%d", fCenBin), vz, vr);
       if (vr < 1) FillHistogram(Form("hPrimPhotAtVertex_cen%d", fCenBin), pt, y, w);
       break;
     case kElectron:
     case -kElectron:
       FillHistogram(Form("hPrimEl_cen%d", fCenBin), pt, y, w);
       break;
     case kPi0:
       FillHistogram(Form("hPrimPi0_cen%d", fCenBin), pt, y, w);
       FillHistogram(Form("hMC_pi0_vertex_cen%d", fCenBin), vz, vr);
       if (vr < 1) FillHistogram(Form("hPrimPi0AtVertex_cen%d", fCenBin), pt, y, w);
       break;
     case 221:
       FillHistogram(Form("hPrimEta_cen%d", fCenBin), pt, y, w);
       FillHistogram(Form("hMC_eta_vertex_cen%d", fCenBin), vz, vr);
       break;
     case kPiPlus:
     case kPiMinus:
       FillHistogram(Form("hPrimPipm_cen%d", fCenBin), pt, y, w);
     case kProton:
       FillHistogram(Form("hPrimP_cen%d", fCenBin), pt, y, w);
       break;
     case kProtonBar:
       FillHistogram(Form("hPrimPbar_cen%d", fCenBin), pt, y, w);
       break;
     case kNeutron:
       FillHistogram(Form("hPrimN_cen%d", fCenBin), pt, y, w);
       break;
     case kNeutronBar:
       FillHistogram(Form("hPrimNbar_cen%d", fCenBin), pt, y, w);
       break;
     case 310:
       FillHistogram(Form("hPrimK0S_cen%d", fCenBin), pt, y, w);
       break;
     case 130:
       FillHistogram(Form("hPrimK0L_cen%d", fCenBin), pt, y, w);
       break;
     case  321:
     case -321:
       FillHistogram(Form("hPrimKpm_cen%d", fCenBin), pt, y, w);
       break;
     default:
       FillHistogram(Form("hPrimOther_cen%d", fCenBin), pt, y, w);
     }
   }

}

//__________________________________________________________________________
Double_t AliAnalysisTaskGammaPHOSPbPbRun2::getR(AliAODMCParticle *particle)
{
    if (!particle) {
      printf("no particle!\n");
      return 0.;
    }
  
    Double_t vx = particle->Xv();
    Double_t vy = particle->Yv();
    Double_t vz = particle->Zv() - fVtx5[2];
 
    return TMath::Sqrt(vx*vx + vy*vy + vz*vz);
}

//__________________________________________________________________________
void AliAnalysisTaskGammaPHOSPbPbRun2::AddQAHistograms()
{
  const Int_t nRuns=200;
  Int_t nPtPhot = 400;
  Double_t ptPhotMax = 40.;

  Int_t nM       = 500;
  Double_t mMin  = 0.0;
  Double_t mMax  = 1.0;

  fOutputContainer->Add(new TH2F("hSelEvents","Event selection", 10,0.,10.,nRuns,0.,float(nRuns)));
  fOutputContainer->Add(new TH1F("hTotSelEvents","Event selection", 10,0.,10.));
  
  //vertex distribution
  fOutputContainer->Add(new TH2F("hZvertex","Z vertex position", 50, -25., 25., nRuns,0.,float(nRuns)));
  
  //Centrality
  fOutputContainer->Add(new TH1F("hCentralityBins", "Centrality bins", fNCenBins, 0, 1*fNCenBins));
  auto hh = (TH1F*)fOutputContainer->Last();
  for (Int_t i = 0; i < fNCenBins; i ++) {
    Int_t blow =  Int_t(fCenBinEdges[i]);
    Int_t bup =  Int_t(fCenBinEdges[i+1]);
    hh->GetXaxis()->SetBinLabel(i+1, Form("cen%d: %d-%d", i, blow, bup));
  }
  
  fOutputContainer->Add(new TH2F("hCentrality","Event centrality", 100,0.,100.,nRuns,0.,float(nRuns)));
  fOutputContainer->Add(new TH1F("hCentralityCorr","Event centrality", 100,0.,100.));
  fOutputContainer->Add(new TH1F("hBadCentrality","Bad PHOS event centrality", 100,0.,100.));
  fOutputContainer->Add(new TH2F("hCenPHOS","Centrality vs PHOSclusters", 100,0.,100.,200,0.,200.));
  fOutputContainer->Add(new TH2F("hCenPHOSCells","Centrality vs PHOS cells", 100,0.,100.,100,0.,1000.));
  fOutputContainer->Add(new TH3F("hCenPHOSCellsM12","Centrality vs PHOS cells", 100,0.,100.,100,0.,1000.,100,0.,1000.));
  fOutputContainer->Add(new TH3F("hCenPHOSCellsM23","Centrality vs PHOS cells", 100,0.,100.,100,0.,1000.,100,0.,1000.));
  fOutputContainer->Add(new TH3F("hCenPHOSCellsM13","Centrality vs PHOS cells", 100,0.,100.,100,0.,1000.,100,0.,1000.));
  fOutputContainer->Add(new TH3F("hCenPHOSCellsM14","Centrality vs PHOS cells", 100,0.,100.,100,0.,1000.,100,0.,1000.));
  fOutputContainer->Add(new TH3F("hCenPHOSCellsM24","Centrality vs PHOS cells", 100,0.,100.,100,0.,1000.,100,0.,1000.));
  fOutputContainer->Add(new TH3F("hCenPHOSCellsM34","Centrality vs PHOS cells", 100,0.,100.,100,0.,1000.,100,0.,1000.));
  fOutputContainer->Add(new TH2F("hCenTrack","Centrality vs tracks", 100,0.,100.,100,0.,15000.));  
  fOutputContainer->Add(new TH2F("hCluEvsClu","ClusterMult vs E",200,0.,20.,100,0.,100.));
  fOutputContainer->Add(new TH2F("hCluEvsCluM","ClusterMult vs E",200,0.,20.,100,0.,20.));
  fOutputContainer->Add(new TH2F("hCenTOF","Centrality vs PHOS TOF", 100,0.,100.,600,-6.e-6,6.e-6));
  
  //PHOS QA
  fOutputContainer->Add(new TH1F("hBadMod","PHOS module without cells",6,2.,8.));
  fOutputContainer->Add(new TH1I("hCellMultEvent"  ,"PHOS cell multiplicity per event"    ,2000,0,2000));
  fOutputContainer->Add(new TH1F("hCellEnergy"  ,"Cell energy"            ,3000,0.,30.));

  //Pi0 masses
  TString histTitle;

  for (Int_t i = 0; i < fNCenBins; i ++) {
    for  (auto cut : fPidCuts) {
       TString cutName  = cut.first;
       TString cutTitle = cut.second;
       
       //Pt
       histTitle = "(M,p_{T})_{#gamma#gamma}, " + cutTitle +";M^{#gamma#gamma}_{inv};p_{1T}+p_{2T}, GeV/c";
       fOutputContainer6->Add(new TH2F(Form("hPi0%sPt_cen%d", cutName.Data(), i), histTitle,  nM, mMin, mMax, nPtPhot, 0., ptPhotMax));
       fOutputContainer6->Add(new TH2F(Form("hMiPi0%sPt_cen%d", cutName.Data(), i), histTitle , nM, mMin, mMax, nPtPhot, 0., ptPhotMax));
       histTitle = histTitle + ";asym";
       
       histTitle = "(M,p_{T})_{#gamma#gamma}, " + cutTitle +";M^{#gamma#gamma}_{inv};p_{T}, GeV/c";
       fOutputContainer6->Add(new TH2F(Form("hPi0%sSingle_cen%d", cutName.Data(), i), histTitle,  nM, mMin, mMax, nPtPhot, 0., ptPhotMax));
       fOutputContainer6->Add(new TH2F(Form("hMiPi0%sSingle_cen%d", cutName.Data(), i), histTitle , nM, mMin, mMax, nPtPhot, 0., ptPhotMax));

       fOutputContainer6->Add(new TH2F(Form("hPi0%sSingle_FirstTOF_cen%d", cutName.Data(), i), histTitle,  nM, mMin, mMax, nPtPhot, 0., ptPhotMax));
       fOutputContainer6->Add(new TH2F(Form("hMiPi0%sSingle_FirstTOF_cen%d", cutName.Data(), i), histTitle , nM, mMin, mMax, nPtPhot, 0., ptPhotMax));

       fOutputContainer6->Add(new TH2F(Form("hPi0%sSingle_NoTOF_cen%d", cutName.Data(), i), histTitle,  nM, mMin, mMax, nPtPhot, 0., ptPhotMax));
       fOutputContainer6->Add(new TH2F(Form("hMiPi0%sSingle_NoTOF_cen%d", cutName.Data(), i), histTitle , nM, mMin, mMax, nPtPhot, 0., ptPhotMax));

       //cuts on first photon only
       histTitle = "(M,p_{T})_{#gamma#gamma}, " + cutTitle +";M^{#gamma#gamma}_{inv};p_{1T}, GeV/c";
       fOutputContainer6->Add(new TH2F(Form("hPi0%sFirst_cen%d", cutName.Data(), i), histTitle,  nM, mMin, mMax, nPtPhot, 0., ptPhotMax));
       fOutputContainer6->Add(new TH2F(Form("hMiPi0%sFirst_cen%d", cutName.Data(), i), histTitle , nM, mMin, mMax, nPtPhot, 0., ptPhotMax));

       if (cutName.Contains("All")) continue;

       fOutputContainer6->Add(new TH2F(Form("hPi0Anti%sFirst_cen%d", cutName.Data(), i), histTitle,  nM, mMin, mMax, nPtPhot, 0., ptPhotMax));
       fOutputContainer6->Add(new TH2F(Form("hMiPi0Anti%sFirst_cen%d", cutName.Data(), i), histTitle , nM, mMin, mMax, nPtPhot, 0., ptPhotMax));

    }
  }

  //Bad Map
  for (Int_t imod = 1; imod <= 4; imod ++) {
     fOutputContainer->Add(new TH1I(Form("hCellMultEventM%d", imod),Form("PHOS cell multiplicity per event, M%d", imod),2000,0,2000));
     fOutputContainer->Add(new TH1F(Form("hCellEnergyM%d", imod),Form("Cell energy in module %d", imod),3000,0.,30.));
     fOutputContainer->Add(new TH2F(Form("hCellNXZM%d", imod), Form("Cell (X,Z), M%d", imod) ,64,0.5,64.5, 56,0.5,56.5));
     fOutputContainer->Add(new TH2F(Form("hCellEXZM%d", imod), Form("Cell E(X,Z), M%d", imod) ,64,0.5,64.5, 56,0.5,56.5));
     fOutputContainer->Add(new TH2F(Form("hCluLowM%d", imod),  Form("Cell (X,Z), M%d", imod), 64,0.5,64.5, 56,0.5,56.5));
     fOutputContainer->Add(new TH2F(Form("hCluHighM%d", imod), Form("Cell (X,Z), M%d", imod), 64,0.5,64.5, 56,0.5,56.5));
     fOutputContainer->Add(new TH2F(Form("hCluVetoM%d", imod) , Form("Cell (X,Z), M%d", imod), 64,0.5,64.5, 56,0.5,56.5));
     fOutputContainer->Add(new TH2F(Form("hCluDispM%d", imod), Form("Cell (X,Z), M%d", imod), 64,0.5,64.5, 56,0.5,56.5));
     fOutputContainer->Add(new TH2F(Form("hTofM%d", imod)     ,Form("TOF in M%d", imod),     100,0.,20.,400,-4.e-6,4.e-6));

     for (Int_t imod2 = 1; imod2 <= imod; imod2 ++) {
       fOutputContainer->Add(new TH2F(Form("hPi0M%d%d", imod2, imod), Form("Pairs in modules %d and %d", imod2, imod),  nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
       fOutputContainer->Add(new TH2F(Form("hMiPi0M%d%d", imod2, imod), Form("Pairs in modules %d and %d", imod2, imod) ,nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
     }
  }
 
}

//__________________________________________________________________________
void AliAnalysisTaskGammaPHOSPbPbRun2::AddEventPlaneHistograms()
{
  //Reaction plane
//  fOutputContainer->Add(new TH3F("hPHOSphi","cos" ,10,0.,100.,20,0.,10.,100,-TMath::Pi(),TMath::Pi()));

   fOutputContainer4->Add(new TProfile("ResV0A","Resolution V0A", 100,0.,100.));
   fOutputContainer4->Add(new TProfile("ResV0C","Resolution V0C", 100,0.,100.));
   fOutputContainer4->Add(new TProfile("ResTPC","Resolution TPC", 100,0.,100.));
   fOutputContainer4->Add(new TProfile("InvResV0A","Inverse Resolution V0A", 100,0.,100.));
   fOutputContainer4->Add(new TProfile("InvResV0C","Inverse Resolution V0C", 100,0.,100.));
   fOutputContainer4->Add(new TProfile("InvResTPC","Inverse Resolution TPC", 100,0.,100.));  
   
   fOutputContainer4->Add(new TProfile("cos2TPCAC","RP correlation between TPC subs",100,0.,100.));
   fOutputContainer4->Add(new TProfile("cos2V0AC","RP correlation between VO A and C sides",100,0.,100.));
   fOutputContainer4->Add(new TProfile("cos2V0ATPC","RP correlation between TPC and V0A",100,0.,100.));
   fOutputContainer4->Add(new TProfile("cos2V0CTPC","RP correlation between TPC and V0C",100,0.,100.));
   fOutputContainer4->Add(new TProfile("qcos2TPCAC","RP correlation between TPC subs",100,0.,100.));
   fOutputContainer4->Add(new TProfile("qcos2V0AC","RP correlation between VO A and C sides",100,0.,100.));
   fOutputContainer4->Add(new TProfile("qcos2V0ATPC","RP correlation between TPC and V0A",100,0.,100.));
   fOutputContainer4->Add(new TProfile("qcos2V0CTPC","RP correlation between TPC and V0C",100,0.,100.));
   
   fOutputContainer4->Add(new TH2F("QV0A","Q_{V0A}",100,0.,100.,200,0.,200.));
   fOutputContainer4->Add(new TH2F("QV0C","Q_{V0C}",100,0.,100.,200,0.,300.));
   fOutputContainer4->Add(new TH2F("QTPC","Q_{TPC}",100,0.,100.,200,0.,100.));
   

   fOutputContainer4->Add(new TProfile("resV0A","Estimated resolution of V0A",100,0.,100.));
   fOutputContainer4->Add(new TProfile("resV0C","Estimated resolution of V0C",100,0.,100.));
   fOutputContainer4->Add(new TProfile("resTPC","Estimated resolution of TPC",100,0.,100.));
   fOutputContainer4->Add(new TProfile("qresV0A","Estimated resolution of V0A",100,0.,100.));
   fOutputContainer4->Add(new TProfile("qresV0C","Estimated resolution of V0C",100,0.,100.));
   fOutputContainer4->Add(new TProfile("qresTPC","Estimated resolution of TPC",100,0.,100.));
   
   fOutputContainer4->Add(new TProfile("VaVbcen","Estimated resolution of TPC",200,0.,100.));
   
   fOutputContainer4->Add(new TH2F("phiRP","RP distribution with TPC", 100,0.,TMath::TwoPi()/fHarmonics,100,0.,100.));
   fOutputContainer4->Add(new TH2F("phiRPflat","RP distribution with TPC flat", 100,0.,TMath::TwoPi()/fHarmonics,100,0.,100.));
   fOutputContainer4->Add(new TH2F("phiRPV0A","RP distribution with V0A", 100,0.,TMath::TwoPi()/fHarmonics,100,0.,100.));
   fOutputContainer4->Add(new TH2F("phiRPV0C","RP distribution with V0C", 100,0.,TMath::TwoPi()/fHarmonics,100,0.,100.));
   fOutputContainer4->Add(new TH2F("phiRPV0Aflat","RP distribution with V0 flat", 100,0.,TMath::TwoPi()/fHarmonics,100,0.,100.));
   fOutputContainer4->Add(new TH2F("phiRPV0Cflat","RP distribution with V0 flat", 100,0.,TMath::TwoPi()/fHarmonics,100,0.,100.));
   fOutputContainer4->Add(new TH2F("phiRPV0AFlow","RP V0A, Flow Package", 100,0.,TMath::TwoPi()/fHarmonics,100,0.,100.));
   fOutputContainer4->Add(new TH2F("phiRPV0CFlow","RP V0C, Flow Package", 100,0.,TMath::TwoPi()/fHarmonics,100,0.,100.));
                                   
   fOutputContainer4->Add(new TProfile2D("phiRPQ","RP distribution with TPC flat", 100,0.,TMath::TwoPi()/fHarmonics,100,0.,100.));
   fOutputContainer4->Add(new TProfile2D("phiRPV0AQ","RP distribution with V0 flat", 100,0.,TMath::TwoPi()/fHarmonics,100,0.,100.));
   fOutputContainer4->Add(new TProfile2D("phiRPV0CQ","RP distribution with V0 flat", 100,0.,TMath::TwoPi()/fHarmonics,100,0.,100.));
   fOutputContainer4->Add(new TProfile2D("phiRPQflat","RP distribution with TPC flat", 100,0.,TMath::TwoPi()/fHarmonics,100,0.,100.));
   fOutputContainer4->Add(new TProfile2D("phiRPV0AQflat","RP distribution with V0 flat", 100,0.,TMath::TwoPi()/fHarmonics,100,0.,100.));
   fOutputContainer4->Add(new TProfile2D("phiRPV0CQflat","RP distribution with V0 flat", 100,0.,TMath::TwoPi()/fHarmonics,100,0.,100.));
   fOutputContainer4->Add(new TH3F("phiRPV0AC","RP distribution with V0A + TPC", 100,0.,TMath::TwoPi()/fHarmonics,100,0.,TMath::TwoPi()/fHarmonics,100,0.,100.));
   fOutputContainer4->Add(new TH3F("phiRPV0ATPC","RP distribution with V0A + TPC", 100,0.,TMath::TwoPi()/fHarmonics,100,0.,TMath::TwoPi()/fHarmonics,100,0.,100.));
   fOutputContainer4->Add(new TH3F("phiRPV0CTPC","RP distribution with V0C + TPC", 100,0.,TMath::TwoPi()/fHarmonics,100,0.,TMath::TwoPi()/fHarmonics,100,0.,100.));

   fOutputContainer4->Add(new TH3F("phiV0ACorrel","V0A my vs Flow", 100,0.,TMath::TwoPi()/fHarmonics,100,0.,TMath::TwoPi()/fHarmonics,100,0.,100.));
   fOutputContainer4->Add(new TH3F("phiV0CCorrel","V0C my vs Flow", 100,0.,TMath::TwoPi()/fHarmonics,100,0.,TMath::TwoPi()/fHarmonics,100,0.,100.));    
   
   fOutputContainer4->Add(new TProfile2D("cos2TPCRP","cos(2Psi_{EP}^{TPC}", 200,-1.,1.,100,0.,100.));
   fOutputContainer4->Add(new TProfile2D("cos2V0ARP","cos(2Psi_{EP}^{V0A}", 200,-1.,1.,100,0.,100.));
   fOutputContainer4->Add(new TProfile2D("cos2V0CRP","cos(2Psi_{EP}^{V0C}", 200,-1.,1.,100,0.,100.));
   fOutputContainer4->Add(new TProfile2D("sin2TPCRP","sin(2Psi_{EP}^{TPC}", 200,-1.,1.,100,0.,100.));
   fOutputContainer4->Add(new TProfile2D("sin2V0ARP","sin(2Psi_{EP}^{V0A}", 200,-1.,1.,100,0.,100.));
   fOutputContainer4->Add(new TProfile2D("sin2V0CRP","sin(2Psi_{EP}^{V0C}", 200,-1.,1.,100,0.,100.));  
}

//__________________________________________________________________________
void AliAnalysisTaskGammaPHOSPbPbRun2::AddMCHistograms()
{
  Int_t nPtPhot = 400;
  Double_t ptPhotMax = 40.;

  // Add MC histograms
  for(Int_t cent = 0; cent < fNCenBins; cent ++){
    fOutputContainer5->Add(new TH2F(Form("hPrimPhot_cen%d", cent), "Primary spectrum", nPtPhot, 0., ptPhotMax, 240, -1.2, 1.2));
    fOutputContainer5->Add(new TH2F(Form("hPrimEl_cen%d", cent), "Primary spectrum", nPtPhot, 0., ptPhotMax, 240, -1.2, 1.2));
    fOutputContainer5->Add(new TH2F(Form("hPrimPi0_cen%d", cent), "Primary spectrum", nPtPhot, 0., ptPhotMax, 240, -1.2, 1.2));
    fOutputContainer5->Add(new TH2F(Form("hPrimEta_cen%d", cent), "Primary spectrum", nPtPhot, 0., ptPhotMax, 240, -1.2, 1.2));
    fOutputContainer5->Add(new TH2F(Form("hPrimPipm_cen%d", cent), "Primary spectrum", nPtPhot, 0., ptPhotMax, 240, -1.2, 1.2));
    fOutputContainer5->Add(new TH2F(Form("hPrimP_cen%d", cent), "Primary spectrum", nPtPhot, 0., ptPhotMax, 240, -1.2, 1.2));
    fOutputContainer5->Add(new TH2F(Form("hPrimPbar_cen%d", cent), "Primary spectrum", nPtPhot, 0., ptPhotMax, 240, -1.2, 1.2));
    fOutputContainer5->Add(new TH2F(Form("hPrimN_cen%d", cent), "Primary spectrum", nPtPhot, 0., ptPhotMax, 240, -1.2, 1.2));
    fOutputContainer5->Add(new TH2F(Form("hPrimNbar_cen%d", cent), "Primary spectrum", nPtPhot, 0., ptPhotMax, 240, -1.2, 1.2));
    fOutputContainer5->Add(new TH2F(Form("hPrimK0S_cen%d", cent), "Primary spectrum", nPtPhot, 0., ptPhotMax, 240, -1.2, 1.2));
    fOutputContainer5->Add(new TH2F(Form("hPrimK0L_cen%d", cent), "Primary spectrum", nPtPhot, 0., ptPhotMax, 240, -1.2, 1.2));
    fOutputContainer5->Add(new TH2F(Form("hPrimKpm_cen%d", cent), "Primary spectrum", nPtPhot, 0., ptPhotMax, 240, -1.2, 1.2));
    fOutputContainer5->Add(new TH2F(Form("hPrimOther_cen%d", cent), "Primary spectrum", nPtPhot, 0., ptPhotMax, 240, -1.2, 1.2));

    fOutputContainer5->Add(new TH2F(Form("hMC_gamma_vertex_cen%d", cent), "Creation vertex",25,0.,25.,1000,0.,500.));
    fOutputContainer5->Add(new TH2F(Form("hMC_pi0_vertex_cen%d", cent), "Creation vertex",25,0.,25.,1000,0.,500.));
    fOutputContainer5->Add(new TH2F(Form("hMC_eta_vertex_cen%d", cent), "Creation vertex",25,0.,25.,1000,0.,500.));

    fOutputContainer5->Add(new TH2F(Form("hPrimPhotAtVertex_cen%d", cent), "MC vertex #gamma {pt,y}", nPtPhot, 0., ptPhotMax, 240, -1.2, 1.2));
    fOutputContainer5->Add(new TH2F(Form("hPrimPi0AtVertex_cen%d", cent), "MC vertex #pi^{0} {pt, y}", nPtPhot, 0., ptPhotMax, 240, -1.2, 1.2));

  }
}

//__________________________________________________________________________
void AliAnalysisTaskGammaPHOSPbPbRun2::AddDistBadHistograms()
{
  Int_t nPtPhot = 400;
  Double_t ptPhotMax = 40.;

  for(Int_t cent = 0; cent < fNCenBins; cent ++){

    fOutputContainer->Add(new TH1F(Form("hPhotAll_DistBad2_cen%d",cent),"All clusters",nPtPhot,0.,ptPhotMax));
    fOutputContainer->Add(new TH1F(Form("hPhotAll_DistBad4_cen%d",cent),"All clusters",nPtPhot,0.,ptPhotMax));
    fOutputContainer->Add(new TH1F(Form("hPhotAll_DistBad6_cen%d",cent),"All clusters",nPtPhot,0.,ptPhotMax));
    fOutputContainer->Add(new TH1F(Form("hPhotAll_DistBad2core_cen%d",cent),"All clusters",nPtPhot,0.,ptPhotMax));
    fOutputContainer->Add(new TH1F(Form("hPhotAll_DistBad4core_cen%d",cent),"All clusters",nPtPhot,0.,ptPhotMax));
    fOutputContainer->Add(new TH1F(Form("hPhotAll_DistBad6core_cen%d",cent),"All clusters",nPtPhot,0.,ptPhotMax));
    fOutputContainer->Add(new TH1F(Form("hPhotAll_DistBad2Disp2_cen%d",cent),"All clusters",nPtPhot,0.,ptPhotMax));
    fOutputContainer->Add(new TH1F(Form("hPhotAll_DistBad4Disp2_cen%d",cent),"All clusters",nPtPhot,0.,ptPhotMax));
    fOutputContainer->Add(new TH1F(Form("hPhotAll_DistBad6Disp2_cen%d",cent),"All clusters",nPtPhot,0.,ptPhotMax));
    fOutputContainer->Add(new TH1F(Form("hPhotAll_DistBad2Disp2core_cen%d",cent),"All clusters",nPtPhot,0.,ptPhotMax));
    fOutputContainer->Add(new TH1F(Form("hPhotAll_DistBad4Disp2core_cen%d",cent),"All clusters",nPtPhot,0.,ptPhotMax));
    fOutputContainer->Add(new TH1F(Form("hPhotAll_DistBad6Disp2core_cen%d",cent),"All clusters",nPtPhot,0.,ptPhotMax));
  }
}

//__________________________________________________________________________
void AliAnalysisTaskGammaPHOSPbPbRun2::AddPhiTitleHistograms()
{
  Int_t nPtPhot = 400;
  Double_t ptPhotMax = 40.;

  for(Int_t cent = 0; cent < 7; cent ++){
    char phiTitle[15];
    Int_t nPt  = 150;
    Double_t xPt[151];
    for(Int_t i=0; i <= 150; i++) 
       xPt[i]=0.1*i;
     const Int_t nPhi=10;
     Double_t xPhi[nPhi];
     for(Int_t i=0; i <= nPhi; i++) 
       xPhi[i] = i * 0.1*TMath::TwoPi()/fHarmonics;
      
     for(Int_t iRP=0; iRP<3; iRP++){
       if(iRP==0)
         snprintf(phiTitle, 15,"TPC");
       if(iRP==1)
         snprintf(phiTitle, 15,"V0A");
       if(iRP==2)
         snprintf(phiTitle, 15,"V0C");

       for  (auto cut : fPidCuts) {
          TString cutName  = cut.first;
          TString cutTitle = cut.second;

          fOutputContainer4->Add(new TProfile2D(Form("hPhotPhiSP%s%s",phiTitle,cutName.Data()),"(M,p_{T},d#phi)_{#gamma#gamma}" ,nPt,xPt,200,0.,100.));
          ((TProfile2D*)fOutputContainer4->Last())->Sumw2();
          
          for(Int_t cent=0; cent < fNCenBins; cent++){
            
            
            fOutputContainer4->Add(new TProfile(Form("hPhotPhiR2%s%s_cen%d",phiTitle,cutName.Data(),cent),"(M,p_{T},d#phi)_{#gamma#gamma}" ,nPt,xPt));
            ((TProfile*)fOutputContainer4->Last())->Sumw2();

            fOutputContainer4->Add(new TProfile(Form("hPhotPhiR%s%s_cen%d",phiTitle,cutName.Data(),cent),"(M,p_{T},d#phi)_{#gamma#gamma}" ,nPt,xPt));
            ((TProfile*)fOutputContainer4->Last())->Sumw2();

            fOutputContainer4->Add(new TProfile(Form("hPhotcos%s%s_cen%d",phiTitle,cutName.Data(),cent),"All clusters",nPtPhot,0.,ptPhotMax));
            ((TProfile*)fOutputContainer4->Last())->Sumw2();
          
            fOutputContainer4->Add(new TProfile(Form("hPhotcosN%s%s_cen%d",phiTitle,cutName.Data(),cent),"All clusters",nPtPhot,0.,ptPhotMax));
            ((TProfile*)fOutputContainer4->Last())->Sumw2();
          }
       }
     }
  }
}

//__________________________________________________________________________
void AliAnalysisTaskGammaPHOSPbPbRun2::AddSinglePhotonHistograms()
{
  //Single photon and pi0 spectrum
  Int_t nPtPhot = 400;
  Double_t ptPhotMax = 40.;
    
  for(Int_t cent=0; cent < fNCenBins; cent++){

    for  (auto cut : fPidCuts) {
      
      TString cutName  = cut.first;
      TString cutTitle = cut.second;

      fOutputContainer2->Add(new TH1F(Form("hPhotPt%s_cen%d",  cutName.Data(), cent), cutTitle.Data(),nPtPhot,0.,ptPhotMax));

      fOutputContainer->Add(new TH2F(Form("hPhotPt%s_TOF_cen%d", cutName.Data(), cent), "Cluster TOF vs pt;p_{T}, GeV/c;TOF, ns", nPtPhot, 0., ptPhotMax, 100, 0, 500.));

      for(Int_t m=1; m<=4; m++){
        fOutputContainer2->Add(new TH1F(Form("hPhotPt%s_M%d_cen%d",cutName.Data(),m,cent), cutTitle.Data(),nPtPhot,0.,ptPhotMax));
      }

      fOutputContainer2->Add(new TH2F(Form("hPhotPtPdg%s_cen%d", cutName.Data(), cent), Form("MC cluster pdg vs pt, %s;p_{T};particle pdg", cutTitle.Data()), nPtPhot, 0., ptPhotMax, 8000, -4000, 4000));
      fOutputContainer2->Add(new TH2F(Form("hPhotPtPdgAtVertex%s_cen%d", cutName.Data(), cent), Form("MC cluster pdg at vertex vs pt, %s;p_{T};particle pdg", cutTitle.Data()), nPtPhot, 0., ptPhotMax, 8000, -4000, 4000));

      fOutputContainer2->Add(new TH1F(Form("hPhotRa%s_cen%d",cutName.Data(),cent), cutTitle.Data(),nPtPhot,0.,ptPhotMax));
      ((TH1F*)fOutputContainer2->Last())->Sumw2();
      fOutputContainer2->Add(new TH1F(Form("hPhotRc%s_cen%d",cutName.Data(),cent), cutTitle.Data(),nPtPhot,0.,ptPhotMax));
      ((TH1F*)fOutputContainer2->Last())->Sumw2();
      fOutputContainer2->Add(new TH1F(Form("hPhotRt%s_cen%d",cutName.Data(),cent), cutTitle.Data(),nPtPhot,0.,ptPhotMax));
      ((TH1F*)fOutputContainer2->Last())->Sumw2();
      fOutputContainer2->Add(new TH1F(Form("hPhotRNa%s_cen%d",cutName.Data(),cent), cutTitle.Data(),nPtPhot,0.,ptPhotMax));
      ((TH1F*)fOutputContainer2->Last())->Sumw2();
      fOutputContainer2->Add(new TH1F(Form("hPhotRNc%s_cen%d",cutName.Data(),cent), cutTitle.Data(),nPtPhot,0.,ptPhotMax));
      ((TH1F*)fOutputContainer2->Last())->Sumw2();
      fOutputContainer2->Add(new TH1F(Form("hPhotRNt%s_cen%d",cutName.Data(),cent), cutTitle.Data(),nPtPhot,0.,ptPhotMax));
      ((TH1F*)fOutputContainer2->Last())->Sumw2();

      fOutputContainer3->Add(new TH2F(Form("hTracksBetaPlus%s_cen%d",  cutName.Data(),cent), cutTitle.Data(),nPtPhot,0.,ptPhotMax, 100, 0,10.));
      fOutputContainer3->Add(new TH2F(Form("hTracksBetaMinus%s_cen%d",  cutName.Data(),cent), cutTitle.Data(),nPtPhot,0.,ptPhotMax, 100,0.,10.));
      fOutputContainer3->Add(new TH2F(Form("hTracksPiPlus%s_cen%d",  cutName.Data(),cent), cutTitle.Data(),nPtPhot,0.,ptPhotMax, 100, 0, 10.));
      fOutputContainer3->Add(new TH2F(Form("hTracksPiMinus%s_cen%d",  cutName.Data(),cent), cutTitle.Data(),nPtPhot,0.,ptPhotMax, 100, 0, 10.));
      fOutputContainer3->Add(new TH2F(Form("hTracksKaPlus%s_cen%d",  cutName.Data(),cent), cutTitle.Data(),nPtPhot,0.,ptPhotMax, 100, 0, 10.));
      fOutputContainer3->Add(new TH2F(Form("hTracksKaMinus%s_cen%d",  cutName.Data(),cent), cutTitle.Data(),nPtPhot,0.,ptPhotMax, 100, 0, 10.));
      fOutputContainer3->Add(new TH2F(Form("hTracksPPlus%s_cen%d",  cutName.Data(),cent), cutTitle.Data(),nPtPhot,0.,ptPhotMax, 100, 0, 10.));
      fOutputContainer3->Add(new TH2F(Form("hTracksPMinus%s_cen%d",  cutName.Data(),cent), cutTitle.Data(),nPtPhot,0.,ptPhotMax, 100, 0, 10.));

      fOutputContainer3->Add(new TH2F(Form("hTracksUndef%s_cen%d",  cutName.Data(),cent), cutTitle.Data(),nPtPhot,0.,ptPhotMax, 100, 0, 10.));

      fOutputContainer3->Add(new TH2F(Form("hMiTracksBetaPlus%s_cen%d",  cutName.Data(),cent), cutTitle.Data(),nPtPhot,0.,ptPhotMax, 100, 0,10.));
      fOutputContainer3->Add(new TH2F(Form("hMiTracksBetaMinus%s_cen%d",  cutName.Data(),cent), cutTitle.Data(),nPtPhot,0.,ptPhotMax, 100,0.,10.));
      fOutputContainer3->Add(new TH2F(Form("hMiTracksPiPlus%s_cen%d",  cutName.Data(),cent), cutTitle.Data(),nPtPhot,0.,ptPhotMax, 100, 0, 10.));
      fOutputContainer3->Add(new TH2F(Form("hMiTracksPiMinus%s_cen%d",  cutName.Data(),cent), cutTitle.Data(),nPtPhot,0.,ptPhotMax, 100, 0, 10.));
      fOutputContainer3->Add(new TH2F(Form("hMiTracksKaPlus%s_cen%d",  cutName.Data(),cent), cutTitle.Data(),nPtPhot,0.,ptPhotMax, 100, 0, 10.));
      fOutputContainer3->Add(new TH2F(Form("hMiTracksKaMinus%s_cen%d",  cutName.Data(),cent), cutTitle.Data(),nPtPhot,0.,ptPhotMax, 100, 0, 10.));
      fOutputContainer3->Add(new TH2F(Form("hMiTracksPPlus%s_cen%d",  cutName.Data(),cent), cutTitle.Data(),nPtPhot,0.,ptPhotMax, 100, 0, 10.));
      fOutputContainer3->Add(new TH2F(Form("hMiTracksPMinus%s_cen%d",  cutName.Data(),cent), cutTitle.Data(),nPtPhot,0.,ptPhotMax, 100, 0, 10.));

      fOutputContainer3->Add(new TH2F(Form("hMiTracksUndef%s_cen%d",  cutName.Data(),cent), cutTitle.Data(),nPtPhot,0.,ptPhotMax, 100, 0, 10.));

      fOutputContainer3->Add(new TH1F(Form("hPhotEmcDistance_cen%d", cent), "Emc distance", 200, 0, 20.));
    }        
    
  }
}

//-----------------------------------------------------------------------------
void AliAnalysisTaskGammaPHOSPbPbRun2::TestMatchingTrackPID(AliCaloPhoton *ph, Bool_t mix)
{
    if (!ph) return;
    Int_t mod = ph->Module();
    TVector3 locpos(ph->EMCx(), 0.,  ph->EMCz());

    Double_t dx = 999.,dz = 999., pttrack = 0.;
    Int_t charge = 0;
    Int_t itr = FindTrackMatching(mod, &locpos, dx,dz, pttrack, charge);
    if (itr < 0) 
       return;
        
    Double_t dist = TestCPV(dx, dz, pttrack, charge);

    FillHistogram(Form("hPhotEmcDistance_cen%d", fCenBin), dist);

    auto trackMatched = (AliAODTrack*)fEvent->GetTrack(itr);

    AliPIDCombined *pidcomb=new AliPIDCombined();

    pidcomb->SetDefaultTPCPriors();
    //pidcomb->SetEnablePriors(kFALSE);
    pidcomb->SetSelectedSpecies(AliPID::kSPECIESC);
    pidcomb->SetDetectorMask(AliPIDResponse::kDetTPC|AliPIDResponse::kDetTOF|AliPIDResponse::kDetITS|AliPIDResponse::kDetTRD);
      
    if (!trackMatched->TestFilterBit(32) ||  !(trackMatched->GetTPCsignal() > 0.) ) 
	return;

    Bool_t pidPion = kFALSE , pidKaon = kFALSE , pidProton = kFALSE , pidElectron = kFALSE, pidUndef = kTRUE;
    const Float_t nsigmaMax = 3.;
      
    Float_t  nsigmaElectron =   TMath::Abs( fPIDResponse->NumberOfSigmasTPC( trackMatched, AliPID::kElectron ));
    Float_t  nsigmaPion =       TMath::Abs( fPIDResponse->NumberOfSigmasTPC( trackMatched, AliPID::kPion ));
    Float_t  nsigmaKaon =       TMath::Abs( fPIDResponse->NumberOfSigmasTPC( trackMatched, AliPID::kKaon ) );
    Float_t  nsigmaProton =     TMath::Abs( fPIDResponse->NumberOfSigmasTPC( trackMatched, AliPID::kProton ));
   
    //find min. sigma
    if ((nsigmaPion < nsigmaMax) && ( nsigmaPion < nsigmaKaon) && (nsigmaPion < nsigmaProton) && (nsigmaPion < nsigmaElectron)) {
	pidUndef = kFALSE;    
	pidPion  = kTRUE;
    }
    if ((nsigmaProton < nsigmaMax) && ( nsigmaProton < nsigmaKaon) && (nsigmaProton < nsigmaPion) && (nsigmaProton < nsigmaElectron)) {
       pidUndef  = kFALSE;  
       pidProton = kTRUE;
    }
    if ((nsigmaKaon < nsigmaMax) && ( nsigmaKaon < nsigmaPion) && (nsigmaKaon < nsigmaProton) && (nsigmaKaon < nsigmaElectron)) {
       pidUndef  = kFALSE;  
       pidKaon= kTRUE;
    }
    if ((nsigmaElectron < nsigmaMax) && ( nsigmaElectron < nsigmaKaon) && (nsigmaElectron < nsigmaPion) && (nsigmaElectron < nsigmaProton)) {
       pidUndef  = kFALSE;  
       pidElectron= kTRUE;
    }

    TString trackName;
    TString chStr = charge > 0 ? "Plus" : "Minus";
   
    if (pidPion)     trackName = "Pi"   + chStr;
    else if (pidElectron) trackName = "Beta" + chStr;
    else if (pidProton)   trackName = "P"    + chStr;
    else if (pidKaon)     trackName = "Ka"   + chStr;
    else if (pidUndef)    trackName = "Undef";

    const Bool_t CPVBit   = ph->IsCPVOK();
   // const Bool_t CPV2Bit   = ph->IsCPVOK();
   // const Bool_t DispBit = ph->IsDispOK();
    const Bool_t Disp2Bit = ph->IsDisp2OK();

    std::vector<TString> passed_cuts = {"All", "Allcore"};
    if (CPVBit) {
      passed_cuts.emplace_back("CPV");
      passed_cuts.emplace_back("CPVcore");
    }
    if (Disp2Bit) {
      passed_cuts.emplace_back("Disp2");
      passed_cuts.emplace_back("Disp2core");
    }
    if (CPVBit && Disp2Bit) {
      passed_cuts.emplace_back("Both2");
      passed_cuts.emplace_back("Both2core");
    }

    for (auto cut : passed_cuts) { 
      FillHistogram(Form("h%sTracks%s%s_cen%d", mix ? "Mi" : "", trackName.Data(), cut.Data(), fCenBin), ph->Pt(), dist);
    }
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskGammaPHOSPbPbRun2::Notify()
{
  //
  // Implemented Notify() to read the cross sections
  // and number of trials from pyxsec.root
  //
  
  //if(!fCheckMCCrossSection) return kTRUE;

  // Fetch the aod also from the input in,
  // have todo it in notify
  
  Float_t xsection = 0;
  Float_t trials   = 1;
  fAvgTrials = -1;
  
  TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
  if(!tree) return kFALSE;
  
  TFile *curfile = tree->GetCurrentFile();
  
  if(!curfile) return kFALSE;
  
  if(fCurrFileName == curfile->GetName()) return kFALSE;
  
  fCurrFileName = TString(curfile->GetName());
  
  if (!fh1Xsec||!fh1Trials) {
  //  Printf("%s%d No Histogram fh1Xsec",(char*)__FILE__,__LINE__);
    return kFALSE;
  }
  
  Bool_t ok = PythiaInfoFromFile(fCurrFileName, xsection, trials);
  
  if (!ok) 
    return kFALSE;
  
  fh1Xsec->Fill("<#sigma>",xsection);
  
  // construct a poor man average trials
  Float_t nEntries = (Float_t)tree->GetTree()->GetEntries();
  
  if(trials >= nEntries && nEntries > 0.) fAvgTrials = trials/nEntries;
  
  fh1Trials->Fill("#sum{ntrials}",trials);
  
  //printf("AliAnalysisTaskGammaPHOSPP::Notify() - xs %f, trial %f, avg trials %f\n",xsection,trials, fAvgTrials);
  
  if (fDebug) 
    Printf("Reading File %s",fInputHandler->GetTree()->GetCurrentFile()->GetName());
  
  return kTRUE;
}

//_____________________________________________________________________________________________________
Bool_t AliAnalysisTaskGammaPHOSPbPbRun2::PythiaInfoFromFile(TString file,Float_t & xsec,Float_t & trials)
{
  //
  // get the cross section and the trails either from pyxsec.root or from pysec_hists.root
  // This is to called in Notify and should provide the path to the AOD/ESD file
    
  xsec   = 0;
  trials = 1;
  
  if (file.Contains("root_archive.zip#")) {
    Ssiz_t pos1 = file.Index("root_archive",12,0,TString::kExact);
    Ssiz_t pos  = file.Index("#",1,pos1,TString::kExact);
    Ssiz_t pos2 = file.Index(".root",5,TString::kExact);
    file.Replace(pos+1,pos2-pos1,"");
  } else {
    // not an archive take the basename....
    file.ReplaceAll(gSystem->BaseName(file.Data()),"");
  }
  
  //Printf("%s",file.Data());
  
  TFile *fxsec = TFile::Open(Form("%s%s",file.Data(),"pyxsec.root")); // problem that we cannot really test the existance of a file in a archive so we have to lvie with open error message from root
  if (!fxsec) {
    // next trial fetch the histgram file
    fxsec = TFile::Open(Form("%s%s",file.Data(),"pyxsec_hists.root"));
    if (!fxsec) {
      // not a severe condition but inciate that we have no information
      return kFALSE;
    } else {
      // find the tlist we want to be independtent of the name so use the Tkey
      TKey* key = (TKey*)fxsec->GetListOfKeys()->At(0);
      if (!key) {
        fxsec->Close();
        return kFALSE;
      }
      
      TList *list = dynamic_cast<TList*>(key->ReadObj());
      if (!list) {
        fxsec->Close();
        return kFALSE;
      }
      
      xsec    = ((TProfile*)list->FindObject("h1Xsec"))  ->GetBinContent(1);
      trials  = ((TH1F*)    list->FindObject("h1Trials"))->GetBinContent(1);
      fxsec->Close();
    }
  } else {
      TTree *xtree = (TTree*)fxsec->Get("Xsection");
      if(!xtree) {
        fxsec->Close();
        return kFALSE;
      }
      
      UInt_t   ntrials  = 0;
      Double_t  xsection  = 0;
      xtree->SetBranchAddress("xsection",&xsection);
      xtree->SetBranchAddress("ntrials",&ntrials);
      xtree->GetEntry(0);
      trials = ntrials;
      xsec = xsection;
      fxsec->Close();
  }
  return kTRUE;
}
//______________________________________________
void AliAnalysisTaskGammaPHOSPbPbRun2::SetCentralityIntervals(TString mode) {
   fCenBinEdges = {0.};   
   if (mode.Contains("0")) fCenBinEdges.emplace_back(1.);
   if (mode.Contains("1")) fCenBinEdges.emplace_back(5.);
   if (mode.Contains("2")) fCenBinEdges.emplace_back(10.);
   if (mode.Contains("3")) fCenBinEdges.emplace_back(15.);
   if (mode.Contains("4")) fCenBinEdges.emplace_back(20.);
   if (mode.Contains("5")) fCenBinEdges.emplace_back(30.);
   if (mode.Contains("6")) fCenBinEdges.emplace_back(40.);
   if (mode.Contains("7")) fCenBinEdges.emplace_back(50.);
   if (mode.Contains("8")) fCenBinEdges.emplace_back(60.);
   if (mode.Contains("9")) fCenBinEdges.emplace_back(70.);
   if (mode.Contains("a")) fCenBinEdges.emplace_back(80.);
   if (mode.Contains("b")) fCenBinEdges.emplace_back(90.);
   fCenBinEdges.emplace_back(100.);

   SetNCenBins(fCenBinEdges.size() - 1);
   
   for (Int_t i = 0; i < fNCenBins; i ++) {
       printf("Centrality bin nr %d: %.f < cent < %.f\n", i, fCenBinEdges[i], fCenBinEdges[i+1]);
   }
   return;
}
//_______________________________________________
Int_t AliAnalysisTaskGammaPHOSPbPbRun2::FindTrackMatching(Int_t mod,TVector3 *locpos,
					    Double_t &dx, Double_t &dz,
					    Double_t &pt,Int_t &charge){
  //Find track with closest extrapolation to cluster
  AliESDEvent *esd = 0x0;
  AliAODEvent *aod = fEvent;
  //AliAODEvent *aod = dynamic_cast<AliAODEvent*>(InputEvent());
  /*
  if(fTender){
    esd= dynamic_cast<AliESDEvent*>(fTender->GetEvent());
    aod= dynamic_cast<AliAODEvent*>(fTender->GetEvent());
  }else{ 
    esd=dynamic_cast<AliESDEvent*>(fTask->InputEvent());
    aod=dynamic_cast<AliAODEvent*>(fTask->InputEvent());
  }
  */
  
  if(!esd && !aod){
    AliError("Neither AOD nor ESD was found");
    return -1;
  }
  Double_t  magF =0.;
   if(esd)
     magF = esd->GetMagneticField();
   if(aod)
     magF = aod->GetMagneticField();
 
  Double_t magSign = 1.0;
  if(magF<0)magSign = -1.0;
  
  if (!TGeoGlobalMagField::Instance()->GetField()) {
    AliError("Margnetic filed was not initialized, use default");
    AliMagF* field = new AliMagF("Maps","Maps", magSign, magSign, AliMagF::k5kG);
    TGeoGlobalMagField::Instance()->SetField(field);
  }

  // *** Start the matching
  Int_t nt = 0;
  if(esd)
    nt = esd->GetNumberOfTracks();
  else
    nt = aod->GetNumberOfTracks();
      
  //Calculate actual distance to PHOS module
  TVector3 globaPos;
  fPHOSGeo->Local2Global(mod, 0.,0., globaPos);
  const Double_t rPHOS = globaPos.Pt(); //Distance to center of  PHOS module
  const Double_t kYmax = 72.+10.; //Size of the module (with some reserve) in phi direction
  const Double_t kZmax = 64.+10.; //Size of the module (with some reserve) in z direction
  const Double_t kAlpha0=330./180.*TMath::Pi(); //First PHOS module angular direction
  const Double_t kAlpha= 20./180.*TMath::Pi(); //PHOS module angular size
  Double_t minDistance = 1.e6;


  Double_t gposTrack[3]; 

  Double_t bz = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->SolenoidField();
  bz = TMath::Sign(0.5*kAlmost0Field,bz) + bz;

  Double_t b[3]; 
  Int_t itr=-1;
  AliESDtrack *esdTrack=0x0;
  AliAODTrack *aodTrack=0x0;
  Double_t xyz[3] = {0}, pxpypz[3] = {0}, cv[21] = {0};
  for (Int_t i=0; i<nt; i++) {
      if(esd)
        esdTrack=esd->GetTrack(i);
      else
        aodTrack=(AliAODTrack*)aod->GetTrack(i);

      // Skip the tracks having "wrong" status (has to be checked/tuned)
      if(esd){
        ULong_t status = esdTrack->GetStatus();
        if((status & AliESDtrack::kTPCout) == 0) continue;
      }
      else{
          
          
      }
      
      
      //Continue extrapolation from TPC outer surface
      AliExternalTrackParam outerParam;
      if(esdTrack){
          outerParam = *(esdTrack->GetOuterParam());
      }
      if(aodTrack){            
        aodTrack->GetPxPyPz(pxpypz);
        aodTrack->GetXYZ(xyz);
        aodTrack->GetCovarianceXYZPxPyPz(cv);
        outerParam.Set(xyz,pxpypz,cv,aodTrack->Charge());
      }
      
      Double_t z; 
      if(!outerParam.GetZAt(rPHOS,bz,z))
        continue;

      if (TMath::Abs(z) > kZmax) 
        continue; // Some tracks miss the PHOS in Z

     
      //Direction to the current PHOS module
      Double_t phiMod=kAlpha0-kAlpha*mod;
      if(!outerParam.RotateParamOnly(phiMod)) continue; //RS use faster rotation if errors are not needed 
    
      Double_t y;                       // Some tracks do not reach the PHOS
      if (!outerParam.GetYAt(rPHOS,bz,y)) continue; //    because of the bending
      
      if(TMath::Abs(y) < kYmax){
        outerParam.GetBxByBz(b);
        outerParam.PropagateToBxByBz(rPHOS,b);        // Propagate to the matching module
        //outerParam.CorrectForMaterial(...); // Correct for the TOF material, if needed
        outerParam.GetXYZ(gposTrack);
        TVector3 globalPositionTr(gposTrack);
        TVector3 localPositionTr;
        fPHOSGeo->Global2Local(localPositionTr,globalPositionTr,mod);
        Double_t ddx = locpos->X()-localPositionTr.X();
        Double_t ddz = locpos->Z()-localPositionTr.Z();
        Double_t d2 = ddx*ddx + ddz*ddz;
        if(d2 < minDistance) {
	  dx = ddx;
  	  dz = ddz;
	  minDistance=d2;
	  itr=i;
          if(esdTrack){
            pt=esdTrack->Pt();
	    charge=esdTrack->Charge();
          }
          else{            
           pt=aodTrack->Pt();
           charge=aodTrack->Charge();
          }
        }
      }
    }//Scanned all tracks
 
   return itr;
}

//____________________________________________________________________________
Double_t AliAnalysisTaskGammaPHOSPbPbRun2::TestCPV(Double_t dx, Double_t dz, Double_t pt, Int_t charge){
  //Parameterization of LHC10h period
  //_true if neutral_
  
 // AliAODEvent *aod= dynamic_cast<AliAODEvent*>(InputEvent());
  
  Double_t mf = fEvent->GetMagneticField();
  
  Double_t meanX=0;
  Double_t meanZ=0.;
  Double_t sx=0.; 
  Double_t sz=0.; 

  if(fRunNumber<209122){ //Run1
    sx=TMath::Min(5.4,2.59719e+02*TMath::Exp(-pt/1.02053e-01)+
              6.58365e-01*5.91917e-01*5.91917e-01/((pt-9.61306e-01)*(pt-9.61306e-01)+5.91917e-01*5.91917e-01)+1.59219);
    sz=TMath::Min(2.75,4.90341e+02*1.91456e-02*1.91456e-02/(pt*pt+1.91456e-02*1.91456e-02)+1.60);
  
    if(mf<0.){ //field --
      meanZ = -0.468318;
      if(charge>0)
        meanX=TMath::Min(7.3, 3.89994*1.20679*1.20679/(pt*pt+1.20679*1.20679)+0.249029+2.49088e+07*TMath::Exp(-pt*3.33650e+01));
      else
        meanX=-TMath::Min(7.7,3.86040*0.912499*0.912499/(pt*pt+0.912499*0.912499)+1.23114+4.48277e+05*TMath::Exp(-pt*2.57070e+01));
    }
    else{ //Field ++
      meanZ= -0.468318;
      if(charge>0)
        meanX=-TMath::Min(8.0,3.86040*1.31357*1.31357/(pt*pt+1.31357*1.31357)+0.880579+7.56199e+06*TMath::Exp(-pt*3.08451e+01));
      else
        meanX= TMath::Min(6.85, 3.89994*1.16240*1.16240/(pt*pt+1.16240*1.16240)-0.120787+2.20275e+05*TMath::Exp(-pt*2.40913e+01));     
    }

  }
  else{//Run2
  
    sx = TMath::Min(5.2, 1.111 + 0.56 * TMath::Exp(-0.031 * pt*pt) + 4.8 /TMath::Power(pt+0.61,3));
    sz = TMath::Min(3.3, 1.12  + 0.35 * TMath::Exp(-0.032 * pt*pt) + 0.75/TMath::Power(pt+0.24,3));

    if(mf<0.){ //field --
      meanZ = 0.102;
      if(charge>0)
        meanX =  TMath::Min(5.8, 0.42 + 0.70 * TMath::Exp(-0.015 * pt*pt) + 35.8/TMath::Power(pt+1.41,3));
      else
        meanX = -TMath::Min(5.8, 0.17 + 0.64 * TMath::Exp(-0.019 * pt*pt) + 26.1/TMath::Power(pt+1.21,3));
    }
    else{ //Field ++
      meanZ= 0.102;
      if(charge>0)
        meanX = -TMath::Min(5.8, 0.58 + 0.68 * TMath::Exp(-0.027 * pt*pt) + 28.0/TMath::Power(pt+1.28,3));
      else
        meanX =  TMath::Min(5.8, 0.11 + 0.67 * TMath::Exp(-0.015 * pt*pt) + 29.9/TMath::Power(pt+1.29,3));
    }

  }
  Double_t rz=(dz-meanZ)/sz;
  Double_t rx=(dx-meanX)/sx;
  return TMath::Sqrt(rx*rx+rz*rz);
}
