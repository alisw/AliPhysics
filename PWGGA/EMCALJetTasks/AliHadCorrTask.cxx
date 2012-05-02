// $Id$
//
// Hadronic correction task.
//
//

#include <TChain.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TParticle.h>
#include "AliAODCaloCluster.h"
#include "AliAODEvent.h"
#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliESDCaloCluster.h"
#include "AliESDtrack.h"
#include "AliFJWrapper.h"
#include "AliPicoTrack.h"
#include "AliVEventHandler.h"
#include "AliHadCorrTask.h"

ClassImp(AliHadCorrTask)

//________________________________________________________________________
AliHadCorrTask::AliHadCorrTask() : 
  AliAnalysisTaskSE("AliHadCorrTask"),
  fTracksName(),
  fCaloName(),
  fOutCaloName(),
  fHadCorr(0),
  fMinPt(0.15),
  fOutClusters(0),
  fOutputList(0),
  fHistNclusvsCent(0),
  fHistNclusMatchvsCent(0),
  fHistEbefore(0),
  fHistEafter(0),
  fHistEoPCent(0),
  fHistNMatchCent(0)
{
  // Default constructor.

  for(Int_t i=0; i<4; ++i) {
    fHistMatchEvsP[i]   = 0;
    fHistMatchdRvsEP[i] = 0;
    for(Int_t j=0; j<5; ++j) 
      fHistMatchEtaPhi[i][j] = 0;
  }
}

//________________________________________________________________________
AliHadCorrTask::AliHadCorrTask(const char *name) : 
  AliAnalysisTaskSE(name),
  fTracksName("Tracks"),
  fCaloName("CaloClusters"),
  fOutCaloName("CaloClustersCorr"),
  fHadCorr(0),
  fMinPt(0.15),
  fOutClusters(0),
  fOutputList(0),
  fHistNclusvsCent(0),
  fHistNclusMatchvsCent(0),
  fHistEbefore(0),
  fHistEafter(0),
  fHistEoPCent(0),
  fHistNMatchCent(0)
{
  // Standard constructor.

  for(Int_t i=0; i<4; ++i) {
    fHistMatchEvsP[i]   = 0;
    fHistMatchdRvsEP[i] = 0;
    for(Int_t j=0; j<5; ++j) 
      fHistMatchEtaPhi[i][j] = 0;
  }

  fBranchNames="ESD:AliESDRun.,AliESDHeader.,PrimaryVertex.";

  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());
}

//________________________________________________________________________
AliHadCorrTask::~AliHadCorrTask()
{
  // Destructor
}

//________________________________________________________________________
Int_t AliHadCorrTask::GetCentBin(Double_t cent) const
{
  // Get centrality bin.

  Int_t centbin = -1;
  if (cent>=0 && cent<10) 
    centbin=0;
  else if (cent>=10 && cent<30) 
    centbin=1;
  else if (cent>=30 && cent<50) 
    centbin=2;
  else if (cent>=50 && cent<=100) 
    centbin=3;

  return centbin;
}

//________________________________________________________________________
Int_t AliHadCorrTask::GetMomBin(Double_t p) const
{
  // Get momenum bin.

  Int_t pbin=-1;
  if (p>=0 && p<0.5) 
    pbin=0;
  else if (p>=0.5 && p<2.) 
    pbin=1;
  else if (p>=2. && p<3.) 
    pbin=2;
  else if (p>=3. && p<5.) 
    pbin=3;
  else if (p>=5.) 
    pbin=4;

  return pbin;
}

//________________________________________________________________________
void AliHadCorrTask::UserCreateOutputObjects()
{
  // Create my user objects.

  AliVEventHandler* handler = AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
  if (!handler) {
    AliError("Input handler not available!");
    return;
  }
 
  if (handler->InheritsFrom("AliESDInputHandler")) {
    fOutClusters = new TClonesArray("AliESDCaloCluster");
  } else if (handler->InheritsFrom("AliAODInputHandler")) {
    fOutClusters = new TClonesArray("AliAODCaloCluster");
  } else {
    AliError("Input handler not recognized!");
    return;
  }
  fOutClusters->SetName(fOutCaloName);
 
  OpenFile(1);
  fOutputList = new TList();
  fOutputList->SetOwner();

  for(Int_t icent=0; icent<4; ++icent) {
    for(Int_t ipt=0; ipt<5; ++ipt){
      TString name(Form("fHistMatchEtaPhi_%i_%i",icent,ipt));
      fHistMatchEtaPhi[icent][ipt] = new TH2F(name, name, 400, -0.2, 0.2, 1600, -0.8, 0.8);
      fOutputList->Add(fHistMatchEtaPhi[icent][ipt]);
    }

    TString name(Form("fHistMatchEvsP_%i",icent));
    fHistMatchEvsP[icent] = new TH2F(name, name, 400, 0., 200., 1000, 0., 10.);
    fOutputList->Add(fHistMatchEvsP[icent]);

    name = Form("fHistMatchdRvsEP_%i",icent);
    fHistMatchdRvsEP[icent] = new TH2F(name, name, 1000, 0., 1., 1000, 0., 10.);
    fOutputList->Add(fHistMatchdRvsEP[icent]);
  }

  fHistNclusvsCent      = new TH1F("Nclusvscent",      "NclusVsCent",      100, 0, 100);
  fHistNclusMatchvsCent = new TH1F("NclusMatchvscent", "NclusMatchVsCent", 100, 0, 100);
  fHistEbefore          = new TH1F("Ebefore",          "Ebefore",          100, 0, 100);
  fHistEafter           = new TH1F("Eafter",           "Eafter",           100, 0, 100);
  fHistEoPCent          = new TH2F("EoPCent",          "EoPCent",          100, 0, 100, 1000, 0,   10);
  fHistNMatchCent       = new TH2F("NMatchesCent",     "NMatchesCent",     100, 0, 100, 101, -0.5, 100.5);
  fOutputList->Add(fHistNclusMatchvsCent);
  fOutputList->Add(fHistNclusvsCent);
  fOutputList->Add(fHistEbefore);
  fOutputList->Add(fHistEafter);
  fOutputList->Add(fHistEoPCent);
  fOutputList->Add(fHistNMatchCent);

  PostData(1, fOutputList);
}

//________________________________________________________________________
void AliHadCorrTask::UserExec(Option_t *) 
{
  // Execute per event.

  // post output in event if not yet present
  if (!(InputEvent()->FindListObject(fOutCaloName)))
    InputEvent()->AddObject(fOutClusters);
  
  // delete output
  fOutClusters->Delete();

  // esd or aod mode
  Bool_t esdMode = kTRUE;
  if (dynamic_cast<AliAODEvent*>(InputEvent()))
      esdMode = kFALSE;

  // optimization in case autobranch loading is off
  AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
  if (fCaloName == "CaloClusters")
    am->LoadBranch("CaloClusters");
  if (fTracksName == "Tracks")
    am->LoadBranch("Tracks");

  // get centrality 
  Double_t cent = -1; 
  AliCentrality *centrality = InputEvent()->GetCentrality() ;
  if (centrality)
    cent = centrality->GetCentralityPercentile("V0M");
  else
    cent=99; // probably pp data
  if (cent<0) {
    AliError(Form("Centrality negative: %f", cent));
    return;
  }
  Int_t centbin = GetCentBin(cent);

  // get input collections
  TClonesArray *tracks = 0;
  TClonesArray *clus   = 0;
  TList *l = InputEvent()->GetList();
  tracks = dynamic_cast<TClonesArray*>(l->FindObject(fTracksName));
  if (!tracks) {
    AliError(Form("Pointer to tracks %s == 0", fTracksName.Data() ));
    return;
  }
  clus = dynamic_cast<TClonesArray*>(l->FindObject(fCaloName));
  if (!clus) {
    AliError(Form("Pointer to clus %s == 0", fCaloName.Data() ));
    return;
  }

  // get primary vertex
  Double_t vertex[3] = {0, 0, 0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);

  // loop over clusters and tracks
  const Int_t Nclus = clus->GetEntries();
  const Int_t Ntrks = tracks->GetEntries();
  for (Int_t iClus = 0, clusCount=0; iClus < Nclus; ++iClus) {
    AliVCluster *c = dynamic_cast<AliVCluster*>(clus->At(iClus));
    if (!c)
      continue;
    if (!c->IsEMCAL())
      continue;

    // make primary particle out of cluster
    TLorentzVector nPart;
    c->GetMomentum(nPart, vertex);

    Double_t etclus = nPart.Pt();
    if (etclus<fMinPt) 
      continue;

    Double_t energyclus = nPart.P();

    // reset cluster/track matching
    c->SetEmcCpvDistance(-1);
    c->SetTrackDistance(999,999);

    // run cluster-track matching
    Int_t    imin       = -1;
    Double_t dEtaMin    = 1e9;
    Double_t dPhiMin    = 1e9;
    Double_t dRmin      = 1e9;
    Double_t totalTrkP  = 0.0; // count total track momentum
    Int_t    Nmatches   = 0;   // count total number of matches
    for (Int_t t = 0; t<Ntrks; ++t) {
      AliVTrack *track = dynamic_cast<AliVTrack*>(tracks->At(t));
      if (!track)
        continue;
      if (!track->IsEMCAL())
        continue;
      if (track->Pt()<fMinPt)
        continue;

      Double_t etadiff=999;
      Double_t phidiff=999;
      AliPicoTrack::GetEtaPhiDiff(track,c,phidiff,etadiff);
      Double_t dR = TMath::Sqrt(etadiff*etadiff+phidiff*phidiff);
      if(dR<dRmin){
        dEtaMin = etadiff;
        dPhiMin = phidiff;
        dRmin = dR;
        imin = t;
      }
      if (fHadCorr>1) {
        Double_t mom = track->P();
        Int_t mombin = GetMomBin(mom);
        if (mombin>-1) {
          fHistMatchEtaPhi[centbin][mombin]->Fill(etadiff,phidiff);
          fHistMatchdRvsEP[centbin]->Fill(dR,energyclus/mom);
        }
      }
      if (TMath::Abs(phidiff)<0.05 && TMath::Abs(etadiff)<0.025) { // pp cuts!!!
        ++Nmatches;
        totalTrkP += track->P();
      }
    }

    // store closest track
    c->SetEmcCpvDistance(imin);
    c->SetTrackDistance(dPhiMin, dEtaMin);

    fHistNclusvsCent->Fill(cent);
    fHistEbefore->Fill(cent,energyclus);
    fHistNMatchCent->Fill(cent,Nmatches);
    if(Nmatches>0) 
      fHistNclusMatchvsCent->Fill(cent);

    // apply correction / subtraction
    if (fHadCorr>0) {
      // to subtract only the closest track set fHadCor to a %
      // to subtract all tracks within the cut set fHadCor to %+1
      if (fHadCorr>1) {
        if (totalTrkP>0) {
          Double_t EoP = energyclus/totalTrkP;
          fHistEoPCent->Fill(cent,EoP);
          fHistMatchEvsP[centbin]->Fill(energyclus,EoP);
        }
        energyclus -= (fHadCorr-1)*totalTrkP;
      } else if (imin>=0) {
        AliVTrack *t = dynamic_cast<AliVTrack*>(tracks->At(imin));
        if (t) {
          Double_t mom = t->P();
          Int_t mombin = GetMomBin(mom);
          fHistMatchEtaPhi[centbin][mombin]->Fill(dEtaMin,dPhiMin);
          if (mom>0){
            fHistMatchEvsP[centbin]->Fill(energyclus,energyclus/mom);
            fHistEoPCent->Fill(cent,energyclus/mom);
            fHistMatchdRvsEP[centbin]->Fill(dRmin,energyclus/mom);
          }
          if (TMath::Abs(dPhiMin)<0.05 && TMath::Abs(dEtaMin)<0.025) { // pp cuts!!!
            energyclus -= fHadCorr*t->P();
          }
        }
      }
    }
    if (energyclus<0) 
      energyclus = 0;
    fHistEafter->Fill(cent,energyclus);

    if (energyclus>0) { // create corrected cluster
      AliVCluster *oc;
      if (esdMode) {
        oc = new ((*fOutClusters)[clusCount]) AliESDCaloCluster(*(dynamic_cast<AliESDCaloCluster*>(c)));
      } else {
        oc = new ((*fOutClusters)[clusCount]) AliAODCaloCluster(*(dynamic_cast<AliAODCaloCluster*>(c)));
      }
      oc->SetE(energyclus);
      ++clusCount;
    }
  }
}

//________________________________________________________________________
void AliHadCorrTask::Terminate(Option_t *) 
{
  // Nothing to be done.
}

