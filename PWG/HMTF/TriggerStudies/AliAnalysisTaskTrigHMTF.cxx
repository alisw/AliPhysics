#include "AliAnalysisTaskTrigHMTF.h"
#include "AliESDEvent.h"
#include "AliESDVZERO.h"
#include "AliESDVZEROfriend.h"
#include "AliESDTZERO.h"
#include "AliVVZERO.h"
#include "AliVAD.h"


// root
#include "TChain.h"
#include "TList.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TDirectory.h"

ClassImp(AliAnalysisTaskTrigHMTF)


//-----------------------------------------------------------------------------
AliAnalysisTaskTrigHMTF::AliAnalysisTaskTrigHMTF(const char* name) :
  AliAnalysisTaskSE(name),
  fTree(0x0),
  fRunNumber(0),
  fBC(0),
  fClassesFired(),
  fInputsL0(),
  fV0ATime(),
  fV0CTime(),
  fBBFlag(),
  fBGFlag(),
  fIsFriend(0),
  fBBFlagPF(),
  fBGFlagPF(),
  fV0ADecision(),
  fV0CDecision(),
  fMTotV0A(0),
  fMTotV0C(0),
  fMRingV0A(),
  fMRingV0C(),
  fTriggerChargeA(0),
  fTriggerChargeC(0),
  fIR1(),
  fIR2(),
  fNofTracklets(0),
  fNofITSClusters(),
  fFOmap(),
  fFiredChipMap(),
  fVertexContributors(0),
  fPileupContributors(0),
  fIsPileupSPD(0),
  fIsIncomplete(0),
  fT0A(),
  fT0C(),
  fTVX(),
  fADATime(),
  fADCTime(),
  fADADecision(),
  fADCDecision(),
  fMTotADA(0),
  fMTotADC(0),
  fTriggerChargeADA(0),
  fTriggerChargeADC(0),
  fVz(0),
  fNITSsaTracks(0),
  fNITSsaTracksHits()

{
  DefineInput(0,TChain::Class());
  DefineOutput(1,TTree::Class());
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void AliAnalysisTaskTrigHMTF::UserCreateOutputObjects(){
  TDirectory *owd = gDirectory;
  OpenFile(1);
  fTree = new TTree("events","events");
  owd->cd();
  fTree->Branch("fIsIncomplete",&fIsIncomplete);
  fTree->Branch("fRunNumber",&fRunNumber);
  fTree->Branch("fBC",&fBC);
  fTree->Branch("fClassesFired",&fClassesFired);
  fTree->Branch("fInputsL0",&fInputsL0);
  fTree->Branch("fVertexContributors",&fVertexContributors);
  fTree->Branch("fPileupContributors",&fPileupContributors);
  fTree->Branch("fIsPileupSPD",&fIsPileupSPD);
  fTree->Branch("fV0ATime",&fV0ATime);
  fTree->Branch("fV0CTime",&fV0CTime);
  fTree->Branch("fBBFlag",&fBBFlag,"fBBFlag[64]/O");
  fTree->Branch("fBGFlag",&fBGFlag,"fBGFlag[64]/O");
  fTree->Branch("fIsFriend",&fIsFriend);
  fTree->Branch("fBBFlagPF",&fBBFlagPF,"fBBFlagPF[21]/l");
  fTree->Branch("fBGFlagPF",&fBGFlagPF,"fBGFlagPF[21]/l");
  fTree->Branch("fV0ADecision",&fV0ADecision);
  fTree->Branch("fV0CDecision",&fV0CDecision);
  fTree->Branch("fMTotV0A",&fMTotV0A);
  fTree->Branch("fMTotV0C",&fMTotV0C);
  fTree->Branch("fMRingV0A",&fMRingV0A,"fMRingV0A[4]/F");
  fTree->Branch("fMRingV0C",&fMRingV0C,"fMRingV0C[4]/F");
  fTree->Branch("fTriggerChargeA",&fTriggerChargeA);
  fTree->Branch("fTriggerChargeC",&fTriggerChargeC);
  fTree->Branch("fIR1",&fIR1);
  fTree->Branch("fIR2",&fIR2);
  fTree->Branch("fNofTracklets",&fNofTracklets);
  fTree->Branch("fNofITSClusters",&fNofITSClusters,"fNofITSClusters[6]/I");
  fTree->Branch("fFOmap",&fFOmap);
  fTree->Branch("fFiredChipMap",&fFiredChipMap);
  fTree->Branch("fT0A",&fT0A,"fT0A[5]/F");
  fTree->Branch("fT0C",&fT0C,"fT0C[5]/F");
  fTree->Branch("fTVX",&fTVX,"fTVX[5]/F");
  fTree->Branch("fADATime",&fADATime);
  fTree->Branch("fADCTime",&fADCTime);
  fTree->Branch("fADADecision",&fADADecision);
  fTree->Branch("fADCDecision",&fADCDecision);
  fTree->Branch("fMTotADA",&fMTotADA);
  fTree->Branch("fMTotADC",&fMTotADC);
  fTree->Branch("fTriggerChargeADA",&fTriggerChargeADA);
  fTree->Branch("fTriggerChargeADC",&fTriggerChargeADC);
  fTree->Branch("fVz", &fVz);
  fTree->Branch("fNITSsaTracks", &fNITSsaTracks);
  fTree->Branch("fNITSsaTracksHits",&fNITSsaTracksHits, "fNITSsaTracksHits[6]/I");
  PostData(1,fTree);
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void AliAnalysisTaskTrigHMTF::UserExec(Option_t *){
  fClassesFired.SetString(fInputEvent->GetFiredTriggerClasses());
  if (!fClassesFired.String().Contains("CINT7-B-NOPF-CENT") &&
      !fClassesFired.String().Contains("CINT7-B-NOPF-ALL") &&
      !fClassesFired.String().Contains("CINT7-I-NOPF-CENT") &&
      !fClassesFired.String().Contains("CVHMV0M-B") &&
      !fClassesFired.String().Contains("CVHMSH1-B") &&
      !fClassesFired.String().Contains("CVHMSH2-B")
  ) return;
  fIsIncomplete = fInputEvent->IsIncompleteDAQ();
  fRunNumber    = fInputEvent->GetRunNumber();
  fBC           = fInputEvent->GetHeader()->GetBunchCrossNumber();
  fInputsL0     = fInputEvent->GetHeader()->GetL0TriggerInputs();
  fIR1          = fInputEvent->GetHeader()->GetIRInt1InteractionMap();
  fIR2          = fInputEvent->GetHeader()->GetIRInt2InteractionMap();
  fNofTracklets = fInputEvent->GetMultiplicity()->GetNumberOfTracklets();
  fFOmap        = fInputEvent->GetMultiplicity()->GetFastOrFiredChipMap();
  fFiredChipMap = fInputEvent->GetMultiplicity()->GetFiredChipMap();
  for (Int_t i=0;i<6;i++) fNofITSClusters[i] = fInputEvent->GetNumberOfITSClusters(i);

  AliVVZERO* vzero = fInputEvent->GetVZEROData();
  fV0ADecision    = vzero->GetV0ADecision();
  fV0CDecision    = vzero->GetV0CDecision();
  fMTotV0A        = vzero->GetMTotV0A();
  fMTotV0C        = vzero->GetMTotV0C();
  fTriggerChargeA = vzero->GetTriggerChargeA();
  fTriggerChargeC = vzero->GetTriggerChargeC();
  fV0ATime        = vzero->GetV0ATime();
  fV0CTime        = vzero->GetV0CTime();
  for (Int_t i=0;i<64;i++) fBBFlag[i] = vzero->GetBBFlag(i);
  for (Int_t i=0;i<64;i++) fBGFlag[i] = vzero->GetBGFlag(i);
  for (Int_t r=0;r<4;r++) fMRingV0A[r] = vzero->GetMRingV0A(r);
  for (Int_t r=0;r<4;r++) fMRingV0C[r] = vzero->GetMRingV0C(r);

  fVz = fInputEvent->GetPrimaryVertexSPD()->GetZ();
  fVertexContributors = fInputEvent->GetPrimaryVertexSPD()->GetNContributors();
  fIsPileupSPD = fInputEvent->IsPileupFromSPD(3,0.8,3.,2.,5.);

  
  AliESDEvent* esd = (AliESDEvent*) fInputEvent;

  AliVAD* ad = (AliVAD*) esd->GetADData();
  fADADecision      = ad->GetADADecision();
  fADCDecision      = ad->GetADCDecision();
  fMTotADA          = ad->GetMTotADA();
  fMTotADC          = ad->GetMTotADC();
  fTriggerChargeADA = ad->GetTriggerChargeA();
  fTriggerChargeADC = ad->GetTriggerChargeC();
  fADATime          = ad->GetADATime();
  fADCTime          = ad->GetADCTime();
  
  fPileupContributors = 0;
  TClonesArray* vertices = esd->GetPileupVerticesSPD();
  for (Int_t i=0;i<vertices->GetEntriesFast();i++){
    AliESDVertex* vertex = (AliESDVertex*) vertices->At(i);
    fPileupContributors += vertex->GetNContributors();
  }

  AliESDTZERO* tzero = (AliESDTZERO*) esd->GetESDTZERO();
  for (Int_t i=0; i<5; i++){
    fT0A[i] = tzero->GetOrA(i);
    fT0C[i] = tzero->GetOrC(i);
    fTVX[i] = tzero->GetTVDC(i);
  }

  for (Int_t ev=0;ev<21;ev++){
    fBBFlagPF[ev] = 0;
    fBGFlagPF[ev] = 0;
  }

  for (Int_t i=0;i<64;i++){
    for (Int_t ev=0;ev<21;ev++){
      fBBFlagPF[ev]|=ULong64_t(vzero->GetPFBBFlag(i,ev)) << i;
      fBGFlagPF[ev]|=ULong64_t(vzero->GetPFBGFlag(i,ev)) << i;
    }
  }

  // if the PF BB/BG flag information is not stored in the ESDs/AODs, use ESD friends
  fIsFriend=0;
  if(!(vzero->TestBit(AliESDVZERO::kPastFutureFlagsFilled)))
    {
      AliESDfriend* esdFriend = ESDfriend();
      if (esdFriend) {
	AliESDVZEROfriend* esdV0friend = esdFriend->GetVZEROfriend();
	if (esdV0friend) {
	  fIsFriend=1;
	  for (Int_t i=0;i<64;i++){
	    for (Int_t ev=0;ev<21;ev++){
	      fBBFlagPF[ev]|=ULong64_t(esdV0friend->GetBBFlag(i,ev)) << i;
	      fBGFlagPF[ev]|=ULong64_t(esdV0friend->GetBGFlag(i,ev)) << i;
	    }
	  }
	}
      }  
    }

  // Loop over tracks
  UInt_t ntrack = esd->GetNumberOfTracks();
  // Reset hits
  for(Int_t ilayer = 0; ilayer < 6; ilayer++){    
     fNITSsaTracksHits[ilayer]=0;
  }
  fNITSsaTracks=0;
  for(UInt_t itrack = 0; itrack < ntrack; itrack++){
    AliESDtrack * trk = esd->GetTrack(itrack);
    UInt_t status = trk->GetStatus();
    if (!(status & AliESDtrack::kITSpureSA )) continue;
    if (!(status & AliESDtrack::kITSrefit  )) continue;
    fNITSsaTracks++;
    UChar_t clumap=trk->GetITSClusterMap();
    for(Int_t ilayer = 0; ilayer < 6; ilayer++){
      if(clumap&(1<<ilayer)) fNITSsaTracksHits[ilayer]++;      
    }    
  }
  
  fTree->Fill();
  PostData(1,fTree);
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void AliAnalysisTaskTrigHMTF::Terminate(Option_t *){
}
