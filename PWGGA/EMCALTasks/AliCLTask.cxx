// Author: C.Loizides

#include "AliCLTask.h"
#include <complex>
#include <Riostream.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TNtuple.h>
#include <TNtupleD.h>
#include <TProfile.h>
#include <TTree.h>
#include <AliAODEvent.h>
#include <AliAODMCParticle.h>
#include <AliAODTrack.h>
#include <AliAnalysisManager.h>
#include <AliCaloPID.h>
#include <AliCalorimeterUtils.h>
#include <AliEventCuts.h>
#include <AliInputEventHandler.h>
#include <AliPIDResponse.h>
#include <AliPhysicsSelection.h>
#include <AliPhysicsSelectionTask.h>
#include <AliTriggerAnalysis.h>

ClassImp(AliCLTask)

//________________________________________________________________________
AliCLTask::AliCLTask(const char *name, Bool_t dolist) : 
  AliAnalysisTaskSE(name), fDoList(dolist), fMcMode(0), fDoClust(0), fDoSkip(0), fDoFilter(0), fPtCut(3), fECut(5), fM2Cut(2), fCounter(0), fCounterC(0), 
  fPidRes(0), fPidCalo(0), fCaloUtils(0), fMyTree(0), fMyInfo(0), fMyTracks(0), fMyClus(0), fMyParts(0), fEventCuts(0), fOutput(0), fTana(0)
{
  // Standard constructor.

  if (name) {
    DefineOutput(1, TTree::Class());
    if (fDoList)
      DefineOutput(2, TList::Class());
  }
}

//________________________________________________________________________
AliCLTask::~AliCLTask()
{
  // Destructor
}

//________________________________________________________________________
void AliCLTask::UserExec(Option_t *option)
{
  // Run various functions.

  fMyTracks->Clear();
  fMyClus->Clear();
  if (fMyParts)
    fMyParts->Clear();

  AliAODEvent *aodev = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!aodev) 
    return;

  const AliVVertex *pvertex = InputEvent()->GetPrimaryVertex();
  if (pvertex->GetNContributors()<1) {
    fCounter++;
    return;
  }

  Double_t vertex[3];
  pvertex->GetXYZ(vertex);
  if (TMath::Abs(vertex[2])>15) {
    fCounter++;
    return;
  }

  if (!fPidRes) {
    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
    fPidRes = inputHandler->GetPIDResponse();
  }

  if (!fTana) {
    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
    TObjArray *l = man->GetTopTasks();
    const Int_t n = l->GetEntries();
    for (Int_t i=0;i<n;++i) {
      AliPhysicsSelectionTask *task = dynamic_cast<AliPhysicsSelectionTask*>(l->At(i));
      if (task) {
	fTana=task->GetPhysicsSelection()->GetTriggerAnalysis();
	break;
      }
    }
  } else {
    fTana = new AliTriggerAnalysis("myown");
  }


  AliAODHeader *h = static_cast<AliAODHeader*>(InputEvent()->GetHeader());
  fMyInfo->fTrig     = h->GetOfflineTrigger();
  if (fTana) {
    fMyInfo->fSPDClsVsTrkBG   = fTana->IsSPDClusterVsTrackletBG(aodev);
    fMyInfo->fV0MOnVsOfPileup = fTana->IsV0MOnVsOfPileup(aodev);
    fMyInfo->fSPDOnVsOfPileup = fTana->IsSPDOnVsOfPileup(aodev);
    fMyInfo->fV0PFPileup      = fTana->IsV0PFPileup(aodev);
    fMyInfo->fSPDVtxPileup    = fTana->IsSPDVtxPileup(aodev);
  }
  fMyInfo->fEvAcc    = fEventCuts->AcceptEvent(aodev); 
  fMyInfo->fSpdPu1   = aodev->IsPileupFromSPD(3,0.8);
  fMyInfo->fSpdPu2   = aodev->IsPileupFromSPD(5,0.8);;
  fMyInfo->fSpdPu    = aodev->IsPileupFromSPDInMultBins();
  fMyInfo->fRun      = h->GetRunNumber();
  fMyInfo->fVz       = vertex[2];
  fMyInfo->fVcon     = pvertex->GetNContributors();
  fMyInfo->fNtracks  = aodev->GetNumberOfTracks();
  fMyInfo->fMult     = h->GetRefMultiplicityComb05();
  fMyInfo->fSkippedV = fCounter;
  fMyInfo->fSkippedC = fCounterC;

  TClonesArray *mcarr = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("mcparticles"));
  Int_t labels[99999] = {-1};

  for(Int_t i=0,newcl=0; i<fMyInfo->fNtracks; ++i) {
    AliAODTrack* track = static_cast<AliAODTrack*>(aodev->GetTrack(i));
    if (!track)
      continue;
    Double_t pt = track->Pt();
    if (pt<0.8)
      continue;
    if (pt>100)
      continue;
    Double_t eta = track->Eta();
    if (TMath::Abs(eta)>1.0) 
      continue;
    Double_t fraction = track->GetTPCFoundFraction();
    if (fraction<=0)
      continue;
    if (fDoFilter) {
      if (!track->IsGlobalConstrained() && !track->IsHybridGlobalConstrainedGlobal() && !track->IsHybridTPCConstrainedGlobal())
	continue;
    }
    Double_t nsigmaTPC = 999.0;	
    Double_t nsigmaTOF = 999.0;
    AliPIDResponse::EDetPidStatus statusTPC = AliPIDResponse::kDetNoSignal;
    AliPIDResponse::EDetPidStatus statusTOF = AliPIDResponse::kDetNoSignal;
    if (fPidRes) { 
      statusTPC = fPidRes->NumberOfSigmas(AliPIDResponse::kTPC, track, (AliPID::EParticleType) 0, nsigmaTPC);
      statusTOF = fPidRes->NumberOfSigmas(AliPIDResponse::kTOF, track, (AliPID::EParticleType) 0, nsigmaTOF);
    }
    Bool_t tpcIsOk = (statusTPC == AliPIDResponse::kDetPidOk);
    Bool_t tofIsOk = (statusTOF == AliPIDResponse::kDetPidOk);
    Double_t m2tof = -1;
    if (tofIsOk)
      m2tof = GetM2tof(track);
    if ((pt<fPtCut) && (m2tof<fM2Cut))
      continue;
    Double_t dedx = track->GetTPCsignal();
    Double_t nSigmaTPCDeut = -1;
    if (tpcIsOk && fPidRes) 
      nSigmaTPCDeut = fPidRes->NumberOfSigmasTPC(track,(AliPID::EParticleType)5);
    Double_t dca[2] = {0,0};   // 0 - dcar, 1 - dcaz -> output parameter
    Double_t cov[3] = {0,0,0}; // covariance matrix -> output parameter
    Double_t kmaxd = 10000.0;  // something big
    AliExternalTrackParam param;
    param.CopyFromVTrack(track);
    if (!param.PropagateToDCA(aodev->GetPrimaryVertex(), aodev->GetMagneticField(), kmaxd, dca, cov)) {
      dca[0] = -1;
      dca[1] = -1;
    }
    CLTrack *myt = (CLTrack*)fMyTracks->ConstructedAt(newcl);
    myt->fPt        = pt;
    myt->fPhi       = TVector2::Phi_0_2pi(track->Phi());
    myt->fEta       = track->Eta();
    myt->fCh        = track->Charge();
    for (Int_t i=14;i<=20;++i) {
      if (track->TestBit(BIT(i)))
	myt->SetBit(BIT(i));
    }
    myt->fIsTpcOk   = tpcIsOk;
    myt->fIsTofOk   = tofIsOk;
    myt->fTsig      = dedx;
    myt->fNsd       = nSigmaTPCDeut;
    myt->fM2tof     = m2tof;
    myt->fTfrac     = fraction;
    myt->fDcar      = dca[0];
    myt->fDcaz      = dca[1];
    myt->fMap       = track->GetFilterMap();
    myt->fChi2pdf   = track->GetTPCchi2perNDF();
    if (mcarr) {
      Int_t lab     = TMath::Abs(track->GetLabel());
      AliAODMCParticle *p = dynamic_cast<AliAODMCParticle*>(mcarr->At(lab));
      if (p) { 
	myt->fId    = p->PdgCode();
	if (lab>99999) {
	  AliFatal(Form("Increase array for labels: %d",lab));
	}
	labels[lab] = newcl;
      }
    }
    ++newcl;
    if (fOutput) {
      TH2F *h1 = (TH2F*)fOutput->At(0);
      h1->Fill(pt,m2tof);
    }
  }

  const Int_t nclus = InputEvent()->GetNumberOfCaloClusters();
  for (Int_t ii=0,newcl=0;ii<nclus;++ii) {
    AliVCluster *clus = InputEvent()->GetCaloCluster(ii);
    if (!clus->IsEMCAL())
      continue;
    if (clus->GetM02()<0.1) 
      continue;
    Double_t e=clus->E();
    if (e<fECut) 
      continue;
    if (e>200) 
      continue;
    TLorentzVector clusterVec;
    clus->GetMomentum(clusterVec,vertex);
    CLClus *myc = (CLClus*)fMyClus->ConstructedAt(newcl++);
    myc->fPt = clusterVec.Pt();
    myc->fEta = clusterVec.Eta();
    myc->fPhi = TVector2::Phi_0_2pi(clusterVec.Phi());
    myc->fM02 = clus->GetM02();
    if (fDoClust) {
      Int_t nMax     = -1;
      Double_t mass  = -1;   
      Double_t angle = -1;
      TLorentzVector l1,l2;
      Int_t absId1     =-1, absId2   =-1;
      Float_t distbad1 =-1, distbad2 =1;
      Bool_t fidcut1   = 0, fidcut2 = 0;
      AliVCaloCells *cells = InputEvent()->GetEMCALCells();
      fPidCalo->GetIdentifiedParticleTypeFromClusterSplitting(clus,cells,fCaloUtils,vertex,
							      nMax,mass,angle,l1,l2,absId1,absId2,distbad1,distbad2,fidcut1,fidcut2);
      myc->fNmax = nMax;
      myc->fM12  = mass;
      myc->fE1   = l1.E();
      myc->fE2   = l2.E();
    }
  }

  if (mcarr) {
    if (!fMyParts) {
      fMyParts = new TClonesArray("CLMCpart",100);
      fMyTree->Branch("mc", &fMyParts, 32*1024, 99);
    }
    const Int_t n = mcarr->GetEntries();
    for (Int_t i=0,newcl=0;i<n;++i) {
      AliAODMCParticle *p = static_cast<AliAODMCParticle*>(mcarr->At(i));
      Double_t eta = p->Eta();
      if (TMath::Abs(eta)>1)
	continue;
      Int_t apg = TMath::Abs(p->PdgCode()); 
      if (apg==1000010020) {
	CLMCpart *myp = (CLMCpart*)fMyParts->ConstructedAt(newcl);
	myp->fPt     = p->Pt();
	myp->fEta    = eta;
	myp->fPhi    = TVector2::Phi_0_2pi(p->Phi());
	myp->fId     = p->PdgCode();
        myp->fIsPrim = p->IsPhysicalPrimary();
	Int_t trkl   = labels[i];
	if (trkl>-1) {
	  CLTrack *myt = dynamic_cast<CLTrack*>(fMyTracks->At(trkl));
	  if (myt) {
	    myp->fLab    = trkl+1; //+1 so that zero means no label
	    myt->fLab    = newcl+1;
	  }
	}
	++newcl;
      }
    } 
  }

  if (fMcMode) {
    if (mcarr && (fMyParts->GetEntries()==0)) {
      return; // only accept event if at least one MC particle
    }
  }
  if (fDoSkip) {
    if ((fMyTracks->GetEntries()==0) && (fMyClus->GetEntries()==0)) {
      ++fCounterC;
      return;
    }
  }
  fCounter  = 0;
  fCounterC = 0;
  fMyTree->Fill();
}

//________________________________________________________________________
void AliCLTask::UserCreateOutputObjects()
{
  // Create histograms

  fMyTree = new TTree("CLTree", "CLTree created from AliCLTask");
  if (1) {
    fMyTree->SetDirectory(0);
  } else {
    TFile *f = OpenFile(1); 
    if (f) {
      f->SetCompressionLevel(2);
      fMyTree->SetDirectory(f);
      fMyTree->SetAutoFlush(-4*1024*1024);
      fMyTree->SetAutoSave(0);
    }
  }
  fMyInfo = new CLInfo;
  fMyTree->Branch("info", &fMyInfo, 1024, 99);
  fMyTracks = new TClonesArray("CLTrack",100);
  fMyTree->Branch("tracks", &fMyTracks, 32*1024, 99);
  fMyClus = new TClonesArray("CLClus",100);
  fMyTree->Branch("clus", &fMyClus, 32*1024, 99);

  fPidCalo = new AliCaloPID;
  fPidCalo->Print("");
  fCaloUtils = new AliCalorimeterUtils;
  fCaloUtils->SetEMCALGeometryName("EMCAL_COMPLETE12SMV1_DCAL_8SM");
  fCaloUtils->SetImportGeometryFromFile(0);
  fCaloUtils->InitEMCALGeometry();
  fCaloUtils->Print("");

  fEventCuts = new AliEventCuts;
  fEventCuts->Print();

  PostData(1, fMyTree);

  if (fDoList) {
    fOutput = new TList;
    TH2F *h1 = new TH2F("hM2vsPt",";p_{T} (GeV); m^{2} (GeV)",120,0,6,200,0,10);
    h1->Sumw2();
    fOutput->Add(h1);
    PostData(2, fOutput);
  }
}

//________________________________________________________________________
Double_t AliCLTask::GetM2tof(AliAODTrack *track) const
{
  if (!fPidRes) 
    return 0;
  const Double_t start_time = fPidRes->GetTOFResponse().GetStartTime(track->P());                 // in ps
  Double_t stop_time  = track->GetTOFsignal();                                                    // in ps
  if (fMyParts&&fPidRes->IsTunedOnData()) {
    stop_time = fPidRes->GetTOFsignalTunedOnData(track);
  }
  const Double_t tof        = (stop_time - start_time);                                           // in ps
  const Double_t length     = fPidRes->GetTOFResponse().GetExpectedSignal(track,AliPID::kPion);   // in ps
  if ((length<=1)||(tof<=1)) 
    return -1;
  const Double_t mom        = track->P();
  Double_t m2               = mom*mom * (tof*tof/length/length - 1);
  //cout << m2 << " " << length << " " << tof << " " << c << " " << mom << endl;
  return m2;
}
