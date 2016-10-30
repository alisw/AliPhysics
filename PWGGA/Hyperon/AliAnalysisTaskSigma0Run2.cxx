#include <TChain.h>
#include <TH1.h>
#include <TH2.h>

#include <AliVEvent.h>
#include <AliVTrack.h>
#include <AliAODEvent.h>
#include <AliAODTrack.h>
#include <AliESDtrackCuts.h>
#include <AliPIDResponse.h>
#include <AliAnalysisManager.h>
#include <AliInputEventHandler.h>
#include <AliMCEvent.h>
#include "AliAODConversionMother.h"
#include "AliV0ParticleStrange.h"

#include "AliAnalysisTaskSigma0Run2.h"


ClassImp(AliAnalysisTaskSigma0Run2)

AliAnalysisTaskSigma0Run2::AliAnalysisTaskSigma0Run2() :
AliAnalysisTaskSE("AliAnalysisTaskSigma0Run2"),
fOutputContainer(),
fIsHeavyIon(kFALSE),
fIsMC(kFALSE),
fIsQA(kFALSE),
fMCEvent(NULL),
fMCStack(NULL),
fReco(),
fQA(),
fMC(),
fMCtruth(),
fESDCuts(0x0),
fESDEvent(NULL),
fPIDResponse(0x0),
fV0Reader(NULL),
fV0ReaderName("V0ReaderV1"),
fV0ReaderStrange(NULL),
fV0ReaderStrangeName("V0ReaderStrange"),
fReaderGammas(NULL),
fReaderV0s(NULL),
fHistLambdaInvMass(NULL),
fHistLambdaInvMassPt(NULL),
fHistLambdaInvMassEta(NULL),
fHistLambdaAngle(NULL),
fHistLambdaR(NULL),
fHistPhotonPt(NULL),
fHistPhotonInvMassPt(NULL),
fHistPhotonInvMassEta(NULL),
fHistPhotonAngle(NULL),
fHistPhotonR(NULL),
fHistSigmaInvMass(NULL),
fHistSigmaInvMassPt(NULL),
fHistSigmaInvMassEta(NULL),
fHistSigmaAngle(NULL),
fHistSigmaR(NULL),
fHistNevents(NULL),
fHistNTracksPt(NULL),
fHistNTracks(NULL),
fHistNV0(NULL),
fHistPt(NULL),
fHistZvertex(NULL),
fHistPtMC(NULL)

{
  // default constructor
}

//___________________________________________________
AliAnalysisTaskSigma0Run2::AliAnalysisTaskSigma0Run2(const char* name) :
AliAnalysisTaskSE(name),
fIsHeavyIon(kFALSE),
fIsMC(kFALSE),
fIsQA(kFALSE),
fOutputContainer(),
fMCEvent(NULL),
fMCStack(NULL),
fReco(),
fQA(),
fMC(),
fMCtruth(),
fESDCuts(0x0),
fESDEvent(NULL),
fPIDResponse(0x0),
fV0Reader(NULL),
fV0ReaderName("V0ReaderV1"),
fV0ReaderStrange(NULL),
fV0ReaderStrangeName("V0ReaderStrange"),
fReaderGammas(NULL),
fReaderV0s(NULL),
fHistLambdaInvMass(NULL),
fHistLambdaInvMassPt(NULL),
fHistLambdaInvMassEta(NULL),
fHistLambdaAngle(NULL),
fHistLambdaR(NULL),
fHistPhotonPt(NULL),
fHistPhotonInvMassPt(NULL),
fHistPhotonInvMassEta(NULL),
fHistPhotonAngle(NULL),
fHistPhotonR(NULL),
fHistSigmaInvMass(NULL),
fHistSigmaInvMassPt(NULL),
fHistSigmaInvMassEta(NULL),
fHistSigmaAngle(NULL),
fHistSigmaR(NULL),
fHistNevents(NULL),
fHistNTracks(NULL),
fHistNTracksPt(NULL),
fHistNV0(NULL),
fHistPt(NULL),
fHistZvertex(NULL),
fHistPtMC(NULL)
{
  DefineInput(0,  TChain::Class());
  DefineOutput(1, TList::Class());
}

//___________________________________________________
AliAnalysisTaskSigma0Run2::~AliAnalysisTaskSigma0Run2()
{  
  delete fESDCuts;
  delete fV0Reader;
  delete fV0ReaderName;
  delete fV0ReaderStrange;
  delete fReaderGammas;
  delete fReaderV0s;
}

//___________________________________________________
void AliAnalysisTaskSigma0Run2::UserCreateOutputObjects()
{
  if(fOutputContainer != NULL){
    delete fOutputContainer;
    fOutputContainer = NULL;
  }
  if(fOutputContainer == NULL){
    fOutputContainer = new TList();
    fOutputContainer->SetOwner(kTRUE);
  }
  
  fReco = new TList();
  fReco->SetName("Reco");
  fReco->SetOwner(kTRUE);
  fOutputContainer->Add(fReco);
  
  //RECO Histos
  fHistLambdaInvMass = new TH1F("fHistLambdaInvMass", "#Lambda candidates ; Invariant mass [GeV/#it{c}^{2}]; # of entries",2000,0,2);
  fReco->Add(fHistLambdaInvMass);
  fHistLambdaInvMassPt = new TH2F("fHistLambdaInvMassPt", "#Lambda candidates ; #it{p}_{T}; Invariant mass [GeV/#it{c}^{2}]",100, 0, 20, 2000,0,2);
  fReco->Add(fHistLambdaInvMassPt);
  fHistLambdaInvMassEta = new TH2F("fHistLambdaInvMassEta", "#Lambda candidates ; #eta; Invariant mass [GeV/#it{c}^{2}]",100, 0, 20, 2000,0,2);
  fReco->Add(fHistLambdaInvMassEta);
  fHistLambdaAngle = new TH2F("fHistLambdaAngle", "#Lambda candidates ; #it{p}_{T}; Opening angle",100, 0, 20, 2000,0,2);
  fReco->Add(fHistLambdaAngle);
  fHistLambdaR = new TH2F("fHistLambdaR", "#Lambda candidates ; #it{p}_{T}; Opening angle",100, 0, 20, 2000,0,2);
  fReco->Add(fHistLambdaR);

  fHistPhotonPt = new TH1F("fHistPhotonPt", "Photon candidates ; #it{p}_{T}; # of entries",5000,0,5);
  fReco->Add(fHistPhotonPt);
  fHistPhotonInvMassPt = new TH2F("fHistPhotonInvMassPt", "Photon candidates ; #it{p}_{T}; Invariant mass [GeV/#it{c}^{2}]",100, 0, 5, 2000,0,2);
  fReco->Add(fHistPhotonInvMassPt);
  fHistPhotonInvMassEta = new TH2F("fHistPhotonInvMassEta", "Photon candidates ; #it{p}_{T}; #eta",100, 0, 5, 2000,0,2);
  fReco->Add(fHistPhotonInvMassEta);
  fHistPhotonAngle = new TH2F("fHistPhotonAngle", "Photon candidates ; #it{p}_{T}; Opening angle",100, 0, 5, 2000,0,2);
  fReco->Add(fHistPhotonAngle);
  fHistPhotonR = new TH2F("fHistPhotonR", "Photon candidates ; #it{p}_{T}; Conversion radius",100, 0, 5, 2000,0,200);
  fReco->Add(fHistPhotonR);

  
  fHistSigmaInvMass = new TH1F("fHistSigmaInvMass", "#Sigma candidates ; Invariant mass [GeV/#it{c}^{2}]; # of entries",2000,0,2);
  fReco->Add(fHistSigmaInvMass);
  fHistSigmaInvMassPt = new TH2F("fHistSigmaInvMassPt", "#Sigma candidates ; #it{p}_{T}; Invariant mass [GeV/#it{c}^{2}]",100, 0, 20, 2000,0,2);
  fReco->Add(fHistSigmaInvMassPt);
  fHistSigmaInvMassEta = new TH2F("fHistSigmaInvMassEta", "#Sigma candidates ; #eta; Invariant mass [GeV/#it{c}^{2}]",100, 0, 20, 2000,0,2);
  fReco->Add(fHistSigmaInvMassEta);
  fHistSigmaAngle = new TH2F("fHistSigmaAngle", "#Sigma candidates ; #it{p}_{T}; Opening angle",100, 0, 20, 2000,0,2);
  fReco->Add(fHistSigmaAngle);
  fHistSigmaR = new TH2F("fHistSigmaR", "#Sigma candidates ; #it{p}_{T}; Opening angle",100, 0, 20, 2000,0,2);
  fReco->Add(fHistSigmaR);
  
  //MC Histos  
  if(fIsMC){
    fMC = new TList();
    fMC->SetName("MC");
    fMC->SetOwner(kTRUE);
    fOutputContainer->Add(fMC);
    fMCtruth = new TList();
    fMCtruth->SetName("MC truth");
    fMCtruth->SetOwner(kTRUE);
    fOutputContainer->Add(fMCtruth);
    
    fHistPtMC = new TH1F("hPtMC","MC Pt distribution; #it{p}_{T} (GeV/#it{c}); #tracks", 100, 0., 5.);
    fMC->Add(fHistPtMC);
  
  
  //MC truth histos
  
  }
  
  
  //QA histos
  if(fIsQA==kTRUE){
    fQA= new TList();
    fQA->SetName("QA");
    fQA->SetOwner(kTRUE);
    fOutputContainer->Add(fQA);
    
    fHistNevents = new TH1F("hNevents","Number of processed events;;#events",1,0,1);
    fQA->Add(fHistNevents);
    fHistNTracks = new TH1F("fHistNTracks","Number of tracks ;#tracks;#entries",500,0,500);
    fQA->Add(fHistNTracks);
    fHistNTracksPt = new TH2F("fHistNTracksPt","Number of tracks ; #it{p}_{T};#tracks",100,0,20,500,0,500);
    fQA->Add(fHistNTracksPt);
    fHistPt = new TH1F("hPt","Pt distribution; #it{p}_{T} (GeV/#it{c}); #entries", 1000, 0., 100);
    fQA->Add(fHistPt);
    fHistNV0 = new TH1F("fHistNV0","Number of V0s per event ; #V0;#entries",100,0,100);
    fQA->Add(fHistNV0);
    fHistZvertex = new TH1F("fHistZvertex", "Position z vertex; z vertex pos; #entries", 200, -25, 25);
    fQA->Add(fHistZvertex);
  }
  
  // Get the PID manager from the input event handler
  AliAnalysisManager*   man          = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = static_cast<AliInputEventHandler*>(man->GetInputEventHandler());
  fPIDResponse  = inputHandler->GetPIDResponse();
  
  fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data());
  if(!fV0Reader){printf("Error: No V0 Reader");return;} // GetV0Reader
  
  if(fV0Reader){
    if((AliConvEventCuts*)fV0Reader->GetEventCuts()){
      if(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetCutHistograms()){
        fOutputContainer->Add(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetCutHistograms());
      }
    }
    
    if((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts()){
      if(((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms()){
        fOutputContainer->Add(((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms());
      }
    }
  }
   
  fV0ReaderStrange=(AliV0ReaderStrange*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderStrangeName.Data());
  if(!fV0ReaderStrange){printf("Error: No V0 Reader Strange");return;} // GetV0Reader
  
  if(fV0ReaderStrange){  
    if((AliV0CutsStrange*)fV0ReaderStrange->GetV0Cuts()){
      if((AliV0CutsStrange*)fV0ReaderStrange->GetV0Cuts()->GetCutHistograms()){
        fOutputContainer->Add((AliV0CutsStrange*)fV0ReaderStrange->GetV0Cuts()->GetCutHistograms());
      }
    }
  }
  
  PostData(1, fOutputContainer);
}

//___________________________________________________
void AliAnalysisTaskSigma0Run2::UserExec(Option_t */*option*/)
{    
  AliVEvent *event=InputEvent();
  
  fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data());
  if(!fV0Reader){printf("Error: No V0 Reader");return;} // GetV0Reader
  
  fV0ReaderStrange=(AliV0ReaderStrange*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderStrangeName.Data());
  if(!fV0ReaderStrange){printf("Error: No V0 Reader Strange");return;} // GetV0Reader
  
  //   Int_t eventQuality = ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEventQuality();
  //   if(event->IsIncompleteDAQ()==kTRUE) eventQuality = 2;  // incomplete event
  //   if(eventQuality == 2 || eventQuality == 3){// Event Not Accepted due to MC event missing or wrong trigger for V0ReaderV1 or because it is incomplete
  //     for(Int_t iCut = 0; iCut<fnCuts; iCut++){
  //       fHistoNEvents[iCut]->Fill(eventQuality);
  //     }
  //     return;
  //   }
  
  const AliVVertex *vertex = event->GetPrimaryVertex();
  fHistZvertex->Fill(vertex->GetZ());
  
  fReaderGammas    = fV0Reader->GetReconstructedGammas(); // Gammas from default Cut
  fReaderV0s       = fV0ReaderStrange->GetReconstructedV0s();
  
  // check the analysis type
  Bool_t isAOD = event->IsA() == AliAODEvent::Class();
  
  // in case of ESD analysis create the track cut object
  if (!isAOD && !fESDCuts) fESDCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010();
  
  Int_t nV0 = event->GetNumberOfV0s();
  fHistNV0->Fill(nV0);
  
  fHistNevents->Fill(0);
  
  // loop over all tracks
  const Int_t ntracks = event->GetNumberOfTracks();  
  fHistNTracks->Fill(ntracks);
  
  
  if(fMCEvent){ // Process MC Particle
    fMCStack = fMCEvent->Stack();
    ProcessMCParticles();
  }
  
  
  
  
  for(Int_t iLambda = 0; iLambda < fReaderV0s->GetEntriesFast(); ++iLambda){
    AliV0ParticleStrange* LambdaCandidate =dynamic_cast<AliV0ParticleStrange*>(fReaderV0s->At(iLambda));
//     AliV0ParticleStrange* LambdaCandidate = (AliV0ParticleStrange*) fReaderV0s->At(iLambda);
    if(!LambdaCandidate) continue;
    fHistLambdaInvMass->Fill(LambdaCandidate->GetMass());
    fHistLambdaInvMassPt->Fill(LambdaCandidate->GetMass(), LambdaCandidate->GetPhotonPt());
    fHistLambdaInvMassEta->Fill(LambdaCandidate->GetMass(), LambdaCandidate->GetPhotonEta());
    //   fHistLambdaAngle(NULL),
    fHistLambdaR->Fill(LambdaCandidate->GetMass(), LambdaCandidate->GetConversionRadius());
    if(LambdaCandidate->GetMass()>1.23 || LambdaCandidate->GetMass()<1.) continue;
    
    for(Int_t iGamma=0; iGamma<fReaderGammas->GetEntriesFast(); ++iGamma){
      AliAODConversionPhoton* PhotonCandidate =dynamic_cast<AliAODConversionPhoton*>(fReaderGammas->At(iGamma));
      if(!PhotonCandidate) continue;
      AliAODConversionPhoton *vPhotonCandidate = new AliAODConversionPhoton(PhotonCandidate);
      fHistPhotonPt->Fill(PhotonCandidate->GetPhotonPt());
      fHistPhotonInvMassPt->Fill(PhotonCandidate->GetPhotonPt(), PhotonCandidate->M());
      fHistPhotonInvMassEta->Fill(PhotonCandidate->GetPhotonPt(), PhotonCandidate->GetPhotonEta());
      //   fHistPhotonAngle(NULL),
      fHistPhotonR->Fill(PhotonCandidate->GetPhotonPt(), PhotonCandidate->GetConversionRadius());
      
      AliAODConversionMother *SigmaCandidate = new AliAODConversionMother(LambdaCandidate, vPhotonCandidate); 
      fHistSigmaInvMass->Fill(SigmaCandidate->M());
      fHistSigmaInvMass->Fill(SigmaCandidate->M());
      fHistSigmaInvMassPt->Fill(SigmaCandidate->M(), SigmaCandidate->Pt());
      fHistSigmaInvMassEta->Fill(SigmaCandidate->M(), SigmaCandidate->Eta());
      //       fGammaRadius =  PhotonCandidate->GetConversionRadius();  // fRConvPhoton  
      //     fGammaDCAzToPrimVtx = PhotonCandidate->GetDCAzToPrimVtx();
      
      //   fHistSigmaAngle(NULL),
      //   fHistSigmaR->Fill(SigmaCandidate->GetMass(), SigmaCandidate->GetConversionRadius());
      
      if( fMCEvent ) {
//         LambdaCandidate->GetMCLabelPositive();
//         LambdaCandidate->GetMCLabelNegative();
        
        
      }
      
      
      delete SigmaCandidate;
    }
    
  }
  
  Float_t highest_pT=0;
  
  for (Int_t itrack=0 ;itrack<ntracks; ++itrack) {
    AliVTrack *vtrack = static_cast<AliVTrack*>(event->GetTrack(itrack));
    if(!vtrack) continue;
    if (isAOD) {
      // see also documentation for specific AOD productions
      // https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PWGPPAODTrackCuts
      AliAODTrack *aodTrack = static_cast<AliAODTrack*>(vtrack);
      if (!aodTrack->TestFilterBit(AliAODTrack::kTrkGlobal)) continue;
    } 
    else {
      if (!fESDCuts->IsSelected(vtrack)) continue;
    }
    
    Double_t pt = vtrack->Pt();
    if(pt>highest_pT) highest_pT=pt;
    
    // fill the Pt histogram
    fHistPt->Fill(pt);
    
    // get the MC track and fill the MC hist if it exists
//     AliVParticle *mcPart = GetMCTrack(vtrack);
//     if (mcPart) fHistPtMC->Fill(mcPart->Pt());
    
  }
  
  fHistNTracksPt->Fill(highest_pT, ntracks);
  
  // flush the data
  PostData(1, fOutputContainer);
}


//____________________________________________________________
void AliAnalysisTaskSigma0Run2::ProcessMCParticles()
{  
  // Loop over all primary MC particle
  const AliVVertex* primVtxMC 	= fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX = primVtxMC->GetX();
  Double_t mcProdVtxY = primVtxMC->GetY();
  Double_t mcProdVtxZ = primVtxMC->GetZ();
  
  for(Int_t i = 0; i < fMCStack->GetNtrack(); i++) {
    
    TParticle* particle = (TParticle *)fMCStack->Particle(i);
    if (!particle) continue;
    
  }
}