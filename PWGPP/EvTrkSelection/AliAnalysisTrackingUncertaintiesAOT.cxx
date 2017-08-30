/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, proviyaded that the above copyright notice appears in all *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purapose. It is         *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// Analysis task for the systematic study of the uncertainties related to   //
// the tracking and ITS-TPC matching efficiency for different particle      //
// species.                                                                 //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////


#include "Riostream.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TH3F.h"
#include "THn.h"
#include "TList.h"
#include "TMath.h"
#include "TParticlePDG.h"
//
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliPID.h"
#include "AliESDtrackCuts.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliESDpid.h"
#include "AliESDUtils.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMultSelection.h"
#include "AliStack.h"
#include "AliLog.h"
//
#include "AliAnalysisTrackingUncertaintiesAOT.h"


ClassImp(AliAnalysisTrackingUncertaintiesAOT)


//________________________________________________________________________
AliAnalysisTrackingUncertaintiesAOT::AliAnalysisTrackingUncertaintiesAOT()
: AliAnalysisTaskSE("TaskTestPA"),
  fUseCentrality(kCentOff),
  fMaxDCAxy(2.4),
  fMaxDCAz(3.2),
  fMaxEta(0.8),
  fSPDlayerReq(AliESDtrackCuts::kAny),
  fCrossRowsOverFndCltTPC(0.8),
  fminCent(0.),
  fmaxCent(100.),
  fTriggerClass("CINT1B"),
  fTriggerMask(AliVEvent::kMB),
  fESD(0),
  fESDpid(0),
  fspecie(0),
  fHistNEvents(0x0),
  fHistCent(0x0),
  fHistMC(0x0),
  fHistMCTPConly(0x0),
  fHistData(0x0),
  fHistAllV0multNTPCout(0),
  fHistSelV0multNTPCout(0),
  fMC(0),
  fRequireVtxTracks(kTRUE),
  fUsePtLogAxis(kFALSE),
  fUseGenPt(kFALSE),
  fDoCutV0multTPCout(kFALSE),
  fListHist(0x0),
  fESDtrackCuts(0x0),
  fVertex(0x0),
  fMultSelectionObjectName("MultSelection")
{
    
}
//________________________________________________________________________
AliAnalysisTrackingUncertaintiesAOT::AliAnalysisTrackingUncertaintiesAOT(const char *name)
  : AliAnalysisTaskSE(name),
    fUseCentrality(kCentOff),
    fMaxDCAxy(2.4),
    fMaxDCAz(3.2),
    fMaxEta(0.8),
    fSPDlayerReq(AliESDtrackCuts::kAny),
    fCrossRowsOverFndCltTPC(0.8),
    fminCent(0.),
    fmaxCent(100.),
    fTriggerClass("CINT1B"),
    fTriggerMask(AliVEvent::kMB),
    fESD(0),
    fESDpid(0),
    fspecie(0),
    fHistNEvents(0x0),
    fHistCent(0x0),
    fHistMC(0x0),
    fHistMCTPConly(0x0),
    fHistData(0x0),
    fHistAllV0multNTPCout(0),
    fHistSelV0multNTPCout(0),
    fMC(0),
    fRequireVtxTracks(kTRUE),
    fUsePtLogAxis(kFALSE),
    fUseGenPt(kFALSE),
    fDoCutV0multTPCout(kFALSE),
    fListHist(0x0),
    fESDtrackCuts(0x0),
    fVertex(0x0),
    fMultSelectionObjectName("MultSelection")
{
  //
  // standard constructur
  //
    
  DefineOutput(1, TList::Class());
    
}
//________________________________________________________________________
AliAnalysisTrackingUncertaintiesAOT::~AliAnalysisTrackingUncertaintiesAOT()
{
  // Destructor
  if (fListHist) {
        
    delete fHistNEvents;
    delete fHistCent;
      
    delete fHistMC;
    delete fHistMCTPConly;
    delete fHistData;
      
    delete fHistAllV0multNTPCout;
    delete fHistSelV0multNTPCout;

    delete fListHist;
    fListHist = 0;
  }
}


//________________________________________________________________________
void AliAnalysisTrackingUncertaintiesAOT::UserCreateOutputObjects()
{
  // create track cuts
  //reproduce filtering cuts
  fESDtrackCuts = new AliESDtrackCuts("AliESDtrackCuts","AliESDtrackCuts");
  fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kFALSE,0);
  fESDtrackCuts->SetEtaRange(-fMaxEta, fMaxEta);
  fESDtrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(fCrossRowsOverFndCltTPC);
  
    
  //
  // Create histograms
  //
  fListHist = new TList();
  fListHist->SetOwner(kTRUE);
  //
  // (1.) basic QA and statistics histograms
  //
  fHistNEvents = new TH1F("fHistNEvents", "Number of processed events",7,-1.5,5.5);
  fHistNEvents->Sumw2();
  fHistNEvents->SetMinimum(0);
  fHistNEvents->GetXaxis()->SetBinLabel(1,"Read from ESD");
  fHistNEvents->GetXaxis()->SetBinLabel(2,"Pass Phys. Sel. + Trig");
  fHistNEvents->GetXaxis()->SetBinLabel(3,"With SPD vertex");
  fHistNEvents->GetXaxis()->SetBinLabel(4,"Vertex contributors >0");
  fHistNEvents->GetXaxis()->SetBinLabel(5,"|zVertex|<10");
  fHistNEvents->GetXaxis()->SetBinLabel(6,"Error on zVertex<0.5");
  fHistNEvents->GetXaxis()->SetBinLabel(7,"Good Z vertex");
    
  fListHist->Add(fHistNEvents);
    
  fHistCent = new TH1F("histCent","Selected centrality; Percentile",100,0.,100.);
  fListHist->Add(fHistCent);
    
  TH2F * histVertexSelection = new TH2F("histVertexSelection", "vertex selection; vertex z (cm); accepted/rejected", 100, -50., 50., 2, -0.5, 1.5);
  fListHist->Add(histVertexSelection);
  TH2F * histTPCITS = new TH2F("histTPCITS", "TPC vs ITS clusters", 100, 0., 500., 1000, 0, 15000);
  TH2F * histTPCCL1 = new TH2F("histTPCCL1", "TPC vs CL1 ", 100, 0., 200., 1000, 0, 15000);
  TH2F * histTPCntrkl = new TH2F("histTPCntrkl", "TPC vs n tracklets", 100, 0., 200., 1000, 0, 15000);
  fListHist->Add(histTPCITS);
  fListHist->Add(histTPCCL1);
  fListHist->Add(histTPCntrkl);
    
    
  if(fDoCutV0multTPCout) {
    fHistAllV0multNTPCout = new TH2F("HistAllV0multNTPCout", "V0mult vs # TPCout (all) ;V0mult ;# TPCout", 1000, 0., 40000, 1000, 0, 30000);
    fHistSelV0multNTPCout = new TH2F("HistSelV0multNTPCout", "V0mult vs # TPCout (sel) ;V0mult ;# TPCout", 1000, 0., 40000, 1000, 0, 30000);
    fHistAllV0multNTPCout->Sumw2();
    fHistSelV0multNTPCout->Sumw2();
    fListHist->Add(fHistAllV0multNTPCout);
    fListHist->Add(fHistSelV0multNTPCout);
  }

  //
  // (2.) track cut variation histograms
  //
  InitializeTrackCutHistograms();
    
  //add default track cuts in the output list
  fListHist->Add(fESDtrackCuts);
    
  //THnSparses to store DCA distributions for primary/secondary
  //fraction extraction and associated correction factor
  Int_t nEtaBins = 2*fMaxEta/0.1;
  if(fMC) {
    const Int_t nvars = 10;
    Int_t nBins[nvars]   = {300,   64,   29,    29,   18,   nEtaBins,    3,   2,    5,    2};
    Double_t xmin[nvars] = {-3., -3.2,  0.5,   0.5,   0.,   -fMaxEta, -0.5, -2., -0.5, -0.5};
    Double_t xmax[nvars] = {3.,   3.2, 15.0,  15.0, 6.28,    fMaxEta,  2.5,  2.,  4.5,  1.5};
    TString axis[nvars]  = {"DCAxy","DCAz","track p_{T}","particle p_{T}","phi","eta","type (0=prim,1=sec,2=mat)","track label (1=lab>0,-1=lab<0)","species (0=e,1=#pi,2=k,3=p) - MC truth","TOFbc"};
        
    fHistMC  = new THnSparseF("fHistMC","fHistMC", nvars, nBins, xmin, xmax);
    fHistMCTPConly  = new THnSparseF("fHistMCTPConly","fHistMCTPConly", nvars, nBins, xmin, xmax);
    for (Int_t j=0; j<nvars; j++) {
      fHistMC->GetAxis(j)->SetTitle(Form("%s",axis[j].Data()));
      fHistMCTPConly->GetAxis(j)->SetTitle(Form("%s",axis[j].Data()));
    }
    fListHist->Add(fHistMC);
    fListHist->Add(fHistMCTPConly);
  }
  else {
    const Int_t nvars = 6;
    Int_t nBins[nvars]   = {300,   64,    29,    18,  nEtaBins,    2};
    Double_t xmin[nvars] = {-3., -3.2,   0.5,    0.,  -fMaxEta, -0.5};
    Double_t xmax[nvars] = {3.,   3.2,  15.0,  6.28,   fMaxEta,  1.5};
    TString axis[nvars]  = {"DCAxy","DCAz","track p_{T}","phi","eta","TOFbc"};
        
    fHistData  = new THnSparseF("fHistData","fHistData", nvars, nBins, xmin, xmax);
        
    for (Int_t j=0; j<nvars; j++) fHistData->GetAxis(j)->SetTitle(Form("%s",axis[j].Data()));
    fListHist->Add(fHistData);
  }
    
  TH1F * histTOFBC = new TH1F("histTOFBC", "TOF  BC!=0;pt (GeV/c)", 100, 0., 25.);
  fListHist->Add(histTOFBC);
  TH1F * histTOFBC0 = new TH1F("histTOFBC0", "TOF  BC==0;pt (GeV/c)", 100, 0., 25);
  fListHist->Add(histTOFBC0);
    
  //
  // post data
  //
  PostData(1, fListHist);
}

//________________________________________________________________________
void AliAnalysisTrackingUncertaintiesAOT::UserExec(Option_t *)
{
  //
  // main event loop
  //
  Printf("Main event loop");
  fESD = dynamic_cast<AliESDEvent*>( InputEvent() );
  if (!fESD) {
    AliWarning("AliAnalysisTrackingUncertaintiesAOT::Exec(): bad ESD");
    PostData(1, fListHist);
    return;
  }
  fHistNEvents->Fill(-1);
    
  if (!fESDtrackCuts) {
    PostData(1, fListHist);
    return;
  }
    
  if (!fESDpid) fESDpid = ((AliESDInputHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetESDpid();
    
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler = (AliInputEventHandler*)man->GetInputEventHandler();
  if(!inputHandler) {
    AliWarning("No AliInputEventHandler!");
    return;
  }
    
  //
  // Physics Selection
  //
  UInt_t maskPhysSel = inputHandler->IsEventSelected();
  TString firedTriggerClasses = fESD->GetFiredTriggerClasses();
  if(!fMC && (fESD->GetRunNumber()<136851 || fESD->GetRunNumber()>139517)) {
    if(!(firedTriggerClasses.Contains(fTriggerClass.Data()))) return;
  }
  if((maskPhysSel & fTriggerMask)==0) return;
    
  fHistNEvents->Fill(0);
  AliStack* stack = 0x0;
  if(fMC){
    AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if (!eventHandler) {
      AliWarning("ERROR: Could not retrieve MC event handler");
      PostData(1, fListHist);
      return;
    }
    AliMCEvent* mcEvent = eventHandler->MCEvent();
    if (!mcEvent) {
      AliWarning("ERROR: Could not retrieve MC event");
      PostData(1, fListHist);
      return;
    }
    stack = mcEvent->Stack();
    if (!stack) {
      AliWarning("ERROR: stack not available\n");
      PostData(1, fListHist);
      return;
    }
  }
  //
  // Monitor vertex position and selection
  //
  TH2F * histVertexSelection = (TH2F *) fListHist->FindObject("histVertexSelection");
  //
  if (IsVertexAccepted(fESD)) {
    histVertexSelection->Fill(fVertex->GetZ(), 0);
  } else {
    if(fVertex) histVertexSelection->Fill(fVertex->GetZ(), 1);
    PostData(1, fListHist);
    return;
  }
    
  if (fUseCentrality!=kCentOff) {
    if(!IsEventSelectedInCentrality(fESD)) return;
  }
  //
  // Fill track cut variation histograms
  //
  ProcessTracks(stack);
  //
  // Post output data
  //
  PostData(1, fListHist);
}
//---------------------------------------------------------------------------
void AliAnalysisTrackingUncertaintiesAOT::SetUseCentrality(AliAnalysisTrackingUncertaintiesAOT::ECentrality flag) {
  //
  // set centrality estimator
  //
  fUseCentrality=flag;
  if(fUseCentrality<kCentOff||fUseCentrality>=kCentInvalid) Printf("Centrality estimator not valid");
    
  return;
}

//________________________________________________________________________
Bool_t AliAnalysisTrackingUncertaintiesAOT::IsEventSelectedInCentrality(AliESDEvent *ESDevent) {
    
  if(fUseCentrality<kCentOff||fUseCentrality>=kCentInvalid){
    Printf("Centrality estimator not valid");
    return kFALSE;
  }else{
    if(ESDevent->GetRunNumber()<244824) {
      Printf("use OLD centrality fw -- not yet implemented!");
      return kFALSE;
    }
    else {
      Double_t cent=-999;
      AliMultSelection *multSelection = (AliMultSelection*)ESDevent->FindListObject(fMultSelectionObjectName);
      if(!multSelection){
	Printf("AliMultSelection could not be found in the esd event list of objects");
	return kFALSE;
      }
      if(fUseCentrality==kCentV0M){
	cent=multSelection->GetMultiplicityPercentile("V0M");
      }else if(fUseCentrality==kCentV0A){
	cent=multSelection->GetMultiplicityPercentile("V0A");
      }else if(fUseCentrality==kCentZNA){
	cent=multSelection->GetMultiplicityPercentile("ZNA");
      }else if(fUseCentrality==kCentCL1){
	cent=multSelection->GetMultiplicityPercentile("CL1");
      }else {
	Printf("CENTRALITY ESTIMATE WITH ESTIMATOR %d NOT YET IMPLEMENTED FOR NEW FRAMEWORK",(Int_t)fUseCentrality);
	return kFALSE;
      }
      Int_t qual = multSelection->GetEvSelCode();
      if(qual == 199 ) cent=-999;
            
      if(cent>=fminCent && cent<fmaxCent) {
	fHistCent->Fill(cent);
	return kTRUE;
      }
      else return kFALSE;
    }
  }
}
//________________________________________________________________________
void AliAnalysisTrackingUncertaintiesAOT::ProcessTracks(AliStack *stack) {
    
  //
  // fill track cut variation histograms - undo cuts step-by-step and fill histograms
  //
  // initialize histograms
  //
  THnSparseF * histTpcItsMatch = (THnSparseF*) fListHist->FindObject("histTpcItsMatch");
    
  Float_t dca[2], cov[3]; // dca_xy, dca_z, sigma_xy, sigma_xy_z, sigma_z for the vertex cut
    
  TParticle *part=0;
  TParticlePDG *pdgPart=0;
  Int_t code=-999, isph=-1,mfl=-999,uniqueID=-999;
  Float_t partType=-1;
  Float_t label = 1;
  Float_t specie = -10;
  Int_t nTPC=0;
  Int_t nITS=0;
  Int_t ntracklets=0;
  fESD->InitMagneticField();
  Int_t ncl1=0;
  const AliMultiplicity* mult=fESD->GetMultiplicity();
  if(mult){
    ntracklets = mult->GetNumberOfTracklets();
    ncl1 = mult->GetNumberOfITSClusters(1);
  }
    
  if(fDoCutV0multTPCout) { //cut on #tracks kTPCout and V0mult
    Float_t V0mult=0.;
    AliESDVZERO *esdVZERO = (AliESDVZERO*)fESD->GetVZEROData();
    if(esdVZERO) {
      for(int ich=0; ich<64;ich++) V0mult += esdVZERO->GetMultiplicity(ich);
    }
    Int_t nTPCout=0;
    for (Int_t i=0;i<fESD->GetNumberOfTracks();++i) {
      AliESDtrack *track =fESD->GetTrack(i);
      if (!track) continue;
      track->SetESDEvent(fESD);
      if(!track->RelateToVertex(fVertex,fESD->GetMagneticField(),100)) continue;
      if((track->GetStatus() & AliESDtrack::kTPCout)) nTPCout++;
    }
    fHistAllV0multNTPCout->Fill(V0mult,nTPCout);
    if(nTPCout > (0.32*V0mult+750)) return;
    else fHistSelV0multNTPCout->Fill(V0mult,nTPCout);
  }
    
  for (Int_t i=0;i<fESD->GetNumberOfTracks();++i) {
        
    isph=-1;
    code=-999;
    mfl=-999;
    uniqueID=-999;
        
    AliESDtrack *track =fESD->GetTrack(i);
    if (!track) continue;
    track->SetESDEvent(fESD);
    if(!track->RelateToVertex(fVertex,fESD->GetMagneticField(),100)) continue;
    //fill TPCcls histo
    nTPC+=track->GetTPCncls();
    nITS+=track->GetITSNcls();
        
    track->GetImpactParameters(dca, cov);
    if(fMC){
      part    = (TParticle*)stack->Particle(TMath::Abs(track->GetLabel()));
      pdgPart = part->GetPDG();
      code    = pdgPart->PdgCode();
      if(stack->IsPhysicalPrimary(TMath::Abs(track->GetLabel()))) isph=1;
      else {
	isph = 0;
	uniqueID = part->GetUniqueID();
      }
      Int_t indexMoth=part->GetFirstMother();
      if(indexMoth>=0){
	TParticle* moth = stack->Particle(indexMoth);
	Float_t codemoth = TMath::Abs(moth->GetPdgCode());
	mfl = Int_t (codemoth/ TMath::Power(10, Int_t(TMath::Log10(codemoth))));
      }
      if(track->GetLabel()<0.) label = -1;
      if(isph==1) {
	partType = 0; //primaries in MC
      }
      else if(isph==0) {
	if(mfl==3 && uniqueID == kPDecay) partType = 1; //secondaries from strangeness
	else partType = 2;  //from material
      }
      if(TMath::Abs(code)==11)   specie = 0;
      if(TMath::Abs(code)==211)  specie = 1;
      if(TMath::Abs(code)==321)  specie = 2;
      if(TMath::Abs(code)==2212) specie = 3;
    }
    //
    // relevant variables
    //
    Int_t nclsTPC       = track->GetTPCncls();
    Float_t pT          = track->Pt();
    Float_t eta         = track->Eta();
    Float_t phi         = track->Phi();
    Float_t chi2TPC     = track->GetTPCchi2();
    Float_t ncrTPC      = track->GetTPCCrossedRows();
    Int_t nclsTPCF      = track->GetTPCNclsF();
    Float_t nCRoverFC   = track->GetTPCCrossedRows();
    Double_t chi2ITS    = track->GetITSchi2();
    Int_t nclsITS       = track->GetITSclusters(0);
        
    if (nclsTPC != 0) {
      chi2TPC /= nclsTPC;
    } else {
      chi2TPC = 999.;
    }
        
    if (nclsTPCF !=0) {
      nCRoverFC /= nclsTPCF;
    } else {
      nCRoverFC = 999.;
    }
        
    if (nclsITS != 0){
      chi2ITS /= nclsITS;
    }else {
      chi2ITS = 999.;
    }
        
    //
    //  fill TPC->ITS matching efficiency histogram
    //
    Bool_t isMatched = kFALSE;
    //  -> if MC is available: fill it only for true primaries,
    //  -> Postprocessing: plot histogram with 1 divided by histogram with 0 as a function of pT/eta/phi
    // remove all ITS requirements
    Bool_t refit    = fESDtrackCuts->GetRequireITSRefit();
    Float_t chi2tpc = fESDtrackCuts->GetMaxChi2TPCConstrainedGlobal();
    Float_t chi2its = fESDtrackCuts->GetMaxChi2PerClusterITS();
        
    fESDtrackCuts->SetRequireITSRefit(kFALSE);
    fESDtrackCuts->SetMaxChi2TPCConstrainedGlobal(99999.);
    fESDtrackCuts->SetMaxChi2PerClusterITS(999999.);
    fESDtrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kOff);
    //DCA
    fESDtrackCuts->SetMaxDCAToVertexXY(fMaxDCAxy);
    fESDtrackCuts->SetMaxDCAToVertexZ(fMaxDCAz);
        
    //TPC refit
    Bool_t tpcrefit=kFALSE;
    Int_t bcTOF_d=0;
    if((track->GetStatus() & AliESDtrack::kTPCrefit)) tpcrefit=kTRUE;
    if(tpcrefit){
      if(fESDtrackCuts->AcceptTrack(track)) {
	TH1F * histTOFBC = (TH1F *) fListHist->FindObject("histTOFBC");
	TH1F * histTOFBC0 = (TH1F *) fListHist->FindObject("histTOFBC0");
	if(track->GetTOFBunchCrossing()!=0){
	  histTOFBC->Fill(track->Pt());
	  bcTOF_d=1;
	}
	histTOFBC0->Fill(track->Pt());
	for(int iSpec=0; iSpec<5; iSpec++) {
	  if(fspecie&BIT(iSpec)) {
	    if(IsConsistentWithPid(iSpec, track)) {
	      Double_t vecHistTpcItsMatch[kNumberOfAxes] = {static_cast<Double_t>(isMatched), pT, eta, phi, (Double_t)iSpec, (Double_t)isph,(Double_t)bcTOF_d,dca[0]};
	      if(fMC && fUseGenPt) vecHistTpcItsMatch[1] = part->Pt();
	      histTpcItsMatch->Fill(vecHistTpcItsMatch);
	      if(fMC){
		Double_t vec4Sparse[10] = {dca[0],dca[1],pT,part->Pt(),phi,eta,partType,label,specie,(Double_t)bcTOF_d};
		fHistMCTPConly->Fill(vec4Sparse);
	      }
	      break;
	    }
	  }
	}
      }
    }
        
    //apply back the cuts and go for ITS-TPC request
    fESDtrackCuts->SetRequireITSRefit(refit);
    fESDtrackCuts->SetMaxChi2TPCConstrainedGlobal(chi2tpc);
    fESDtrackCuts->SetMaxChi2PerClusterITS(chi2its);
    //set the SPD cluster requirement
    fESDtrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,fSPDlayerReq);
    //set is matched
    isMatched=kTRUE;
    Int_t bcTOF_n=0;
    if (fESDtrackCuts->AcceptTrack(track)) {
      if(track->GetTOFBunchCrossing()!=0){
	bcTOF_n=1;
      }
            
      for(int iSpec=0; iSpec<5; iSpec++) {
	if(fspecie&BIT(iSpec)) {
	  if(IsConsistentWithPid(iSpec, track)) {
	    Double_t vecHistTpcItsMatch[kNumberOfAxes] = {static_cast<Double_t>(isMatched), pT, eta, phi, (Double_t)iSpec, (Double_t)isph,(Double_t)bcTOF_n,dca[0]};
	    if(fMC && fUseGenPt) vecHistTpcItsMatch[1] = part->Pt();
	    histTpcItsMatch->Fill(vecHistTpcItsMatch);
	    if(fMC){
	      Double_t vec4Sparse[10] = {dca[0],dca[1],pT,part->Pt(),phi,eta,partType,label,specie,(Double_t)bcTOF_n};
	      fHistMC->Fill(vec4Sparse);
	    }
	    else {
	      Double_t vec4Sparse[6] = {dca[0],dca[1],pT,phi,eta,(Double_t)bcTOF_n};
	      fHistData->Fill(vec4Sparse);
	    }
	    break;
	  }
	}
      }
    }
  } // end of track loop
  TH2F * histTPCITS = (TH2F *) fListHist->FindObject("histTPCITS");
  TH2F * histTPCCL1 = (TH2F *) fListHist->FindObject("histTPCCL1");
  TH2F * histTPCntrkl = (TH2F *) fListHist->FindObject("histTPCntrkl");
  histTPCITS->Fill(nITS,nTPC);
  histTPCCL1->Fill(ncl1,nTPC);
  histTPCntrkl->Fill(ntracklets,nTPC);
}



//________________________________________________________________________
void AliAnalysisTrackingUncertaintiesAOT::Terminate(Option_t *)
{
  // Draw result to the screen
  // Called once at the end of the query
    
    
}


//________________________________________________________________________
void AliAnalysisTrackingUncertaintiesAOT::BinLogAxis(const THnSparseF *h, Int_t axisNumber) {
  //
  // Method for the correct logarithmic binning of histograms
  //
  TAxis *axis = h->GetAxis(axisNumber);
  int bins = axis->GetNbins();
    
  Double_t from = axis->GetXmin();
  Double_t to = axis->GetXmax();
  Double_t *newBins = new Double_t[bins + 1];
    
  newBins[0] = from;
  Double_t factor = pow(to/from, 1./bins);
    
  for (int i = 1; i <= bins; i++) {
    newBins[i] = factor * newBins[i-1];
  }
  axis->Set(bins, newBins);
  delete [] newBins;
    
}


//________________________________________________________________________
Bool_t AliAnalysisTrackingUncertaintiesAOT::IsVertexAccepted(AliESDEvent * esd) {
  //
  // function to check if a proper vertex is reconstructed and write z-positionin vertexZ
  //
  Bool_t vertexOkay = kFALSE;
  const AliESDVertex *vtx = esd->GetPrimaryVertexTracks();
  //fVertex = new AliESDVertex(*esd->GetPrimaryVertexTracks());
  if(!vtx) return vertexOkay;
  fHistNEvents->Fill(1);
    
  if (vtx->GetNContributors() < 1) {
    if(fRequireVtxTracks) {
      return vertexOkay;
    }
    vtx = esd->GetPrimaryVertexSPD();
    if (vtx->GetNContributors() < 1) vertexOkay = kFALSE;
    else vertexOkay = kTRUE;
    Double_t cov[6]={0};
    vtx->GetCovarianceMatrix(cov);
    Double_t zRes = TMath::Sqrt(cov[5]);
    if (vtx->IsFromVertexerZ() && (zRes>0.5)) vertexOkay = kFALSE;
    fVertex = (AliESDVertex*)vtx;
  } else {
    vertexOkay = kTRUE;
    fVertex = (AliESDVertex*)vtx;
  }
  if(TMath::Abs(vtx->GetZ())>10) vertexOkay=kFALSE;
  return vertexOkay;
}


//________________________________________________________________________
Bool_t AliAnalysisTrackingUncertaintiesAOT::IsConsistentWithPid(Int_t type, const AliESDtrack * const tr) {
  //
  // just check if the PID is consistent with a given hypothesis in order to
  // investigate effects which are only dependent on the energy loss.
  //
  if (type == 0) return IsElectron(tr);
  if (type == 1)     return IsPion(tr);
  if (type == 2)     return IsKaon(tr);
  if (type == 3)   return IsProton(tr);
  if (type == 4)          return kTRUE;
    
  return kFALSE;
    
}

//________________________________________________________________________
Bool_t AliAnalysisTrackingUncertaintiesAOT::IsElectron(const AliESDtrack * const tr, Bool_t useTPCTOF) const {
  //
  // Selection of electron candidates using the upper half of the TPC sigma band, starting at
  // the mean ignoring its shift, and going up to 3 sigma above the mean. In case TOF information
  // is available, tracks which are incompatible with electrons within 3 sigma are rejected. If
  // no TOF information is used, the momentum regions where the kaon and the proton line cross
  // the electron line are cut out using a 3 sigma cut around the kaon or proton line.
  //
    
  Float_t nsigmaElectronTPC = fESDpid->NumberOfSigmasTPC(tr, AliPID::kElectron);
  if(nsigmaElectronTPC < 0 || nsigmaElectronTPC  > 3) return kFALSE;
    
  if(useTPCTOF){
    Float_t nsigmaElectronTOF = fESDpid->NumberOfSigmasTOF(tr, AliPID::kElectron);
    if(TMath::Abs(nsigmaElectronTOF) > 3) return kFALSE;
    else return kTRUE;
  } else {
    Float_t nsigmaKaonTPC = fESDpid->NumberOfSigmasTPC(tr, AliPID::kKaon),
      nsigmaProtonTPC =fESDpid->NumberOfSigmasTPC(tr, AliPID::kProton);
    if(TMath::Abs(nsigmaKaonTPC < 3) || TMath::Abs(nsigmaProtonTPC < 3)) return kFALSE;
    else return kTRUE;
  }
}

//________________________________________________________________________
Bool_t AliAnalysisTrackingUncertaintiesAOT::IsPion(const AliESDtrack * const tr,Bool_t useTPCTOF) const{
  //
  // Selectron of pion candidates
  //
  Float_t nsigmaPionTPC = fESDpid->NumberOfSigmasTPC(tr, AliPID::kPion);
  if (TMath::Abs(nsigmaPionTPC) < 3) return kTRUE;
  if(useTPCTOF){
    Float_t nsigmaPionTOF = fESDpid->NumberOfSigmasTOF(tr, AliPID::kPion);
    if(TMath::Abs(nsigmaPionTOF) > 3) return kFALSE;
    else return kTRUE;
  }
  return kFALSE;
    
}

//________________________________________________________________________
Bool_t AliAnalysisTrackingUncertaintiesAOT::IsKaon(const AliESDtrack * const tr,Bool_t useTPCTOF) const {
  //
  // Selection of kaon candidates
  //
  Float_t nsigmaKaonTPC = fESDpid->NumberOfSigmasTPC(tr, AliPID::kKaon);
  if (TMath::Abs(nsigmaKaonTPC) < 3) return kTRUE;
  if(useTPCTOF){
    Float_t nsigmaKaonTOF = fESDpid->NumberOfSigmasTOF(tr, AliPID::kKaon);
    if(TMath::Abs(nsigmaKaonTOF) > 3) return kFALSE;
    else return kTRUE;
  }
    
  return kFALSE;
    
}

//________________________________________________________________________
Bool_t AliAnalysisTrackingUncertaintiesAOT::IsProton(const AliESDtrack * const tr,Bool_t useTPCTOF) const{
  //
  // Selection of proton candidates
  //
  Float_t nsigmaProtonTPC = fESDpid->NumberOfSigmasTPC(tr, AliPID::kProton);
  if (TMath::Abs(nsigmaProtonTPC) < 3) return kTRUE;
  if(useTPCTOF){
    Float_t nsigmaProtonTOF = fESDpid->NumberOfSigmasTOF(tr, AliPID::kProton);
    if(TMath::Abs(nsigmaProtonTOF) > 3) return kFALSE;
    else return kTRUE;
  }
    
  return kFALSE;
    
}


//________________________________________________________________________
void AliAnalysisTrackingUncertaintiesAOT::InitializeTrackCutHistograms() {
  //
  // create histograms for the track cut studies
  //
  //  match TPC->ITS
  //                                  0-is matched, 1-pt, 2-eta,   3-phi,   4-pid(0-3 -> electron-proton, 4 -> undef, 5 -> all) 6-bcTOF 7-DCAxy
  Int_t nEtaBins = 2*fMaxEta/0.1;
  Int_t    binsTpcItsMatch[kNumberOfAxes] = {    2,   29, nEtaBins,            18,  6,      3,    2,  30};
  Double_t minTpcItsMatch[kNumberOfAxes]  = { -0.5,  0.5, -fMaxEta,             0, -0.5, -1.5, -0.5, -3.};
  Double_t maxTpcItsMatch[kNumberOfAxes]  = {  1.5, 15.0,  fMaxEta, 2*TMath::Pi(),  5.5,  1.5,  1.5,  3.};
  //
  TString axisNameTpcItsMatch[kNumberOfAxes]  = {"isMatched","pT","eta","phi","pid","primSec","bcTOF","dcaxy"};
  TString axisTitleTpcItsMatch[kNumberOfAxes] = {"isMatched","pT","eta","phi","pid","primSec","bcTOF","dcaxy"};
  //
  THnSparseF * histTpcItsMatch = new THnSparseF("histTpcItsMatch","TPC -> ITS matching",kNumberOfAxes, binsTpcItsMatch, minTpcItsMatch, maxTpcItsMatch);
  if(fUsePtLogAxis)BinLogAxis(histTpcItsMatch, 1);
  fListHist->Add(histTpcItsMatch);
  //
  for (Int_t iaxis=0; iaxis<kNumberOfAxes;iaxis++){
    histTpcItsMatch->GetAxis(iaxis)->SetName(axisNameTpcItsMatch[iaxis]);
    histTpcItsMatch->GetAxis(iaxis)->SetTitle(axisTitleTpcItsMatch[iaxis]);
  }
}
