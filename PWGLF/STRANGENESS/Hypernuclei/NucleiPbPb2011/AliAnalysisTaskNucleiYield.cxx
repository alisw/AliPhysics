#include "AliAnalysisTaskNucleiYield.h"

// ROOT includes
#include <TAxis.h>
#include <TChain.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TList.h>
#include <TMath.h>
#include <TParticle.h>
#include <TClonesArray.h>
#include <TTree.h>
#include <TRandom3.h>

// ALIROOT includes
#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliPID.h"
#include "AliTPCPIDResponse.h"
#include "AliTOFPIDResponse.h"
#include "AliVTrack.h"
#include "AliVVertex.h"
#include "AliVEvent.h"
#include "AliPIDResponse.h"
#include "AliVParticle.h"
#include "AliMCEvent.h"
#include "AliInputEventHandler.h"
#include "AliVEventHandler.h"
#include "AliStack.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliAODVertex.h"

#define LIGHT_SPEED 2.99792457999999984e-02 // in the units that TOF like
#define EPS 1.e-15

using TMath::TwoPi;

///\cond CLASSIMP
ClassImp(AliAnalysisTaskNucleiYield);
///\endcond

/// Standard and default constructor of the class.
///
/// \param taskname Name of the task
/// \param partname Name of the analysed particle
///
AliAnalysisTaskNucleiYield::AliAnalysisTaskNucleiYield(TString taskname)
:AliAnalysisTaskSE(taskname.Data())						//
,fList(0x0)
,fPDG(0x0)
,fPDGMass(0)
,fPDGMassOverZ(0)
,fIsMC(kFALSE)
,fPID(0x0)
,fMagField(0.f)
,fPrimaryVertex(0x0)
,fMmc()
,fAmc()
,fDCAzLimit(10.)
,fDCAzNbins(400)
,fPtCorrectionA(3)
,fPtCorrectionM(3)
,fTOFlowBoundary(-2.4)
,fTOFhighBoundary(3.6)
,fTOFnBins(75)
,fRequireITSrefit(kTRUE)
,fRequireTPCrefit(kTRUE)
,fRequireNoKinks(kTRUE)
,fRequireITSrecPoints(2u)
,fRequireITSsignal(0u)
,fRequireSPDrecPoints(1u)
,fRequireTPCrecPoints(70u)
,fRequireTPCsignal(70u)
,fRequireEtaMin(-0.8f)
,fRequireEtaMax(0.8f)
,fRequireYmin(-0.5f)
,fRequireYmax(0.5f)
,fRequireMaxChi2(4.f)
,fRequireMaxDCAxy(0.5f)
,fRequireMaxDCAz(1.f)
,fRequireTPCpidSigmas(3.f)
,fRequireITSpidSigmas(-1.f)
,fRequireMinEnergyLoss(0.)
,fRequireMagneticField(0)
,fParticle(AliPID::kUnknown)
,fCentBins(0x0)
,fDCABins(0x0)
,fPtBins(0x0)
,fCustomTPCpid(0)
,fFlatteningProbs(0)
,fCentrality(0x0)
,fFlattenedCentrality(0x0)
,fCentralityClasses(0x0)
,fProduction(0x0)
,fAITS_TPC(0x0)
,fAITS_TPC_TOF(0x0)
,fATotal(0x0)
,fAPtCorrection(0x0)
,fMITS_TPC(0x0)
,fMITS_TPC_TOF(0x0)
,fMTotal(0x0)
,fMPtCorrection(0x0)
,fMDCAPrimaryTPC(0x0)
,fMDCASecondaryTPC(0x0)
,fMDCAPrimaryTOF(0x0)
,fMDCASecondaryTOF(0x0)
,fATOFsignal(0x0)
,fATPCcounts(0x0)
,fMDCAxyTPC(0x0)
,fMDCAzTPC(0x0)
,fMDCAxyTOF(0x0)
,fMDCAzTOF(0x0)
,fMTOFsignal(0x0)
,fMTPCcounts(0x0) {
  gRandom->SetSeed(0);
  Float_t aCorrection[3] = {-2.10154e-03,-4.53472e-01,-3.01246e+00};
  Float_t mCorrection[3] = {-2.00277e-03,-4.93461e-01,-3.05463e+00};
  fPtCorrectionA.Set(3, aCorrection);
  fPtCorrectionM.Set(3,mCorrection);
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}

/// Standard destructor
///
AliAnalysisTaskNucleiYield::~AliAnalysisTaskNucleiYield(){
  // Destructor
  
  if (AliAnalysisManager::GetAnalysisManager()->IsProofMode()) return;
  if (fList) delete fList;
}

/// This function creates all the histograms and all the objects in general used during the analysis
/// \return void
///
void AliAnalysisTaskNucleiYield::UserCreateOutputObjects(){
  fList = new TList();
  fList->SetOwner(kTRUE);

  const Int_t nPtBins = fPtBins.GetSize() - 1;
  const Int_t nCentBins = fCentBins.GetSize() - 1;
  const Int_t nDCAbins = fDCABins.GetSize() - 1;
  const Float_t *pTbins = fPtBins.GetArray();
  const Float_t *centBins = fCentBins.GetArray();
  const Float_t *dcaBins = fDCABins.GetArray();
  
  fCentrality = new TH1F("fCentrality",";Centrality (%);Events / 1%;",100,0.,100.);
  fCentralityClasses = new TH1F("fCentralityClasses",";Centrality classes(%);Events / Class;",nCentBins,centBins);
  fFlattenedCentrality = new TH1F("fFlattenCentrality","Centrality distribution after the flattening;Centrality (%);Events / 1%;",100,0.,100.);
  fList->Add(fCentrality);
  fList->Add(fCentralityClasses);
  fList->Add(fFlattenedCentrality);
  
  if (fIsMC) {
    fMTotal = new TH2F("fMTotal",";Centrality (%);p_{T} (GeV/c); Counts",
                       nCentBins,centBins,nPtBins,pTbins);
    fATotal = new TH2F("fATotal",";Centrality (%);p_{T} (GeV/c); Counts",
                       nCentBins,centBins,nPtBins,pTbins);
    fMITS_TPC = new TH2F("fMITS_TPC",";Centrality (%);p_{T} (GeV/c); Counts",
                         nCentBins,centBins,nPtBins,pTbins);
    fAITS_TPC = new TH2F("fAITS_TPC",";Centrality (%);p_{T} (GeV/c); Counts",
                         nCentBins,centBins,nPtBins,pTbins);
    fMITS_TPC_TOF = new TH2F("fMITS_TPC_TOF",";Centrality (%);p_{T} (GeV/c); Counts",
                             nCentBins,centBins,nPtBins,pTbins);
    fAITS_TPC_TOF = new TH2F("fAITS_TPC_TOF",";Centrality (%);p_{T} (GeV/c); Counts",
                             nCentBins,centBins,nPtBins,pTbins);
    fMDCAPrimaryTPC = new TH3F("fMDCAPrimaryTPC",";Centrality (%);p_{T} (GeV/c); DCA_{xy} (cm)",
                            nCentBins,centBins,nPtBins,pTbins,nDCAbins,dcaBins);
    fMDCASecondaryTPC = new TH3F("fMDCASecondaryTPC",";Centrality (%);p_{T} (GeV/c); DCA_{xy} (cm)",
                              nCentBins,centBins,nPtBins,pTbins,nDCAbins,dcaBins);
    fMDCAPrimaryTOF = new TH3F("fMDCAPrimaryTOF",";Centrality (%);p_{T} (GeV/c); DCA_{xy} (cm)",
                            nCentBins,centBins,nPtBins,pTbins,nDCAbins,dcaBins);
    fMDCASecondaryTOF = new TH3F("fMDCASecondaryTOF",";Centrality (%);p_{T} (GeV/c); DCA_{xy} (cm)",
                              nCentBins,centBins,nPtBins,pTbins,nDCAbins,dcaBins);
    fAPtCorrection = new TH2F("fAPtCorrection",
                              ";p_{T}^{rec} (GeV/c);p_{T}^{MC}-p_{T}^{rec} (GeV/c);Entries",
                              160,0.4,6.,80,-1.,1.);
    fMPtCorrection = new TH2F("fMPtCorrection",
                              ";p_{T}^{rec} (GeV/c);p_{T}^{MC}-p_{T}^{rec} (GeV/c);Entries",
                              160,0.4,6.,80,-1.,1.);
    fProduction = new TH1F("fProduction",";p (GeV/c);Entries",100,-10,10);
    
    fList->Add(fProduction);
    fList->Add(fMTotal);
    fList->Add(fATotal);
    fList->Add(fMITS_TPC);
    fList->Add(fAITS_TPC);
    fList->Add(fMITS_TPC_TOF);
    fList->Add(fAITS_TPC_TOF);
    fList->Add(fMDCAPrimaryTPC);
    fList->Add(fMDCASecondaryTPC);
    fList->Add(fMDCAPrimaryTOF);
    fList->Add(fMDCASecondaryTOF);
    fList->Add(fAPtCorrection);
    fList->Add(fMPtCorrection);
  
  } else {
  
    float tofBins[fTOFnBins + 1];
    const float deltaTOF = (fTOFhighBoundary - fTOFlowBoundary) / fTOFnBins;
    for (int i = 0; i <= fTOFnBins; ++i)
      tofBins[i] = i * deltaTOF + fTOFlowBoundary;
    float dcazBins[fDCAzNbins + 1];
    const float deltaDCAz = 2.f * fDCAzLimit / fDCAzNbins;
    for (int i = 0; i <= fDCAzNbins; ++i)
      dcazBins[i] = i * deltaDCAz - fDCAzLimit;
    const int nTpcBins = 40;
    float tpcBins[nTpcBins + 1];
    for (int i = 0; i <= nTpcBins; ++i) {
      tpcBins[i] = -4.f + i * 0.2f;
    }
    
    fATOFsignal = new TH3F("fATOFsignal",
                           ";Centrality (%);p_{T} (GeV/c);m_{TOF}^{2}-m_{PDG}^{2} (GeV/c^{2})^{2}",
                           nCentBins,centBins,nPtBins,pTbins,fTOFnBins,tofBins);
    fATPCcounts = new TH2F("fATPCcounts",";Centrality (%);p_{T} (GeV/c); Entries",
                           nCentBins,centBins,nPtBins,pTbins);
    fMDCAxyTPC = new TH3F("fMDCAxyTPC",";Centrality (%);p_{T} (GeV/c); DCA_{xy} (cm)",
                       nCentBins,centBins,nPtBins,pTbins,nDCAbins,dcaBins);
    fMDCAzTPC = new TH3F("fMDCAzTPC",";Centrality (%);p_{T} (GeV/c); DCA_{z} (cm)",
                      nCentBins,centBins,nPtBins,pTbins,fDCAzNbins,dcazBins);
    fMDCAxyTOF = new TH3F("fMDCAxyTOF",";Centrality (%);p_{T} (GeV/c); DCA_{xy} (cm)",
                       nCentBins,centBins,nPtBins,pTbins,nDCAbins,dcaBins);
    fMDCAzTOF = new TH3F("fMDCAzTOF",";Centrality (%);p_{T} (GeV/c); DCA_{z} (cm)",
                      nCentBins,centBins,nPtBins,pTbins,fDCAzNbins,dcazBins);
    fMTOFsignal = new TH3F("fMTOFsignal",
                           ";Centrality (%);p_{T} (GeV/c);m_{TOF}^{2}-m_{PDG}^{2} (GeV/c^{2})^{2}",
                           nCentBins,centBins,nPtBins,pTbins,fTOFnBins,tofBins);
    fMTPCcounts = new TH2F("fMTPCcounts",";Centrality (%);p_{T} (GeV/c); Entries",
                           nCentBins,centBins,nPtBins,pTbins);
    
    fMTPCphiCounts = new TH2F("fMTPCphiCounts",";#phi;p_{T} (GeV/c);Counts",64,0.,TMath::TwoPi(),
                              36,0.2,2.);
    fMTOFphiSignal = new TH3F("fMTOFphiSignal",
                             ";#phi;p_{T} (GeV/c);m_{TOF}^{2}-m_{PDG}^{2} (GeV/c^{2})^{2}",
                              64,0.,TMath::TwoPi(),36,0.2,2.,
                              fTOFnBins,fTOFlowBoundary,fTOFhighBoundary);
    fATPCphiCounts = new TH2F("fATPCphiCounts",";#phi;p_{T} (GeV/c);Counts",64,0.,TMath::TwoPi(),
                              36,0.2,2.);
    fATOFphiSignal = new TH3F("fATOFphiSignal",
                              ";#phi;p_{T} (GeV/c);m_{TOF}^{2}-m_{PDG}^{2} (GeV/c^{2})^{2}",
                              64,0.,TMath::TwoPi(),36,0.2,2.,
                              fTOFnBins,fTOFlowBoundary,fTOFhighBoundary);
    fList->Add(fATOFsignal);
    fList->Add(fATPCcounts);
    fList->Add(fMDCAxyTPC);
    fList->Add(fMDCAzTPC);
    fList->Add(fMDCAxyTOF);
    fList->Add(fMDCAzTOF);
    fList->Add(fMTOFsignal);
    fList->Add(fMTPCcounts);
    fList->Add(fMTPCphiCounts);
    fList->Add(fMTOFphiSignal);
    fList->Add(fATPCphiCounts);
    fList->Add(fATOFphiSignal);
  }
  
  PostData(1,fList);
}

/// This is the function that is evaluated for each event. The analysis code stays here.
///
/// \param options Deprecated parameter
/// \return void
///
void AliAnalysisTaskNucleiYield::UserExec(Option_t *){
  /// What is inside this function:
  AliVEvent *ev = dynamic_cast<AliVEvent*> (InputEvent());
  
  /// Check event selection mask
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* handl = (AliInputEventHandler*)mgr->GetInputEventHandler();
  UInt_t mask = handl->IsEventSelected();
  if (!(mask & 0xffffffff)) {
    PostData(1, fList);
    return;
  }
  if (!((mask & AliVEvent::kMB) || (mask & AliVEvent::kCentral) || (mask & AliVEvent::kSemiCentral))) {
    PostData(1, fList);
    return;
  }
  
  fPID = handl->GetPIDResponse();
  
  AliCentrality *centr = ev->GetCentrality();
  
  /// Centrality selection in PbPb, percentile determined with V0
  float centrality = centr->GetCentralityPercentile("V0M");
  /// Select only events with centralities between 0 and 80 %
  if (centrality < 0. || centrality > 80.) {
    PostData(1, fList);
    return;
  }

  /// Primary vertex displacement cut
  fPrimaryVertex = (AliVVertex*)ev->GetPrimaryVertex();
  if(TMath::Abs(fPrimaryVertex->GetZ()) > 10.) {
    PostData(1, fList);
    return;
  }

  /// Magnetic field cut 
  fMagField = ev->GetMagneticField();
  if (fRequireMagneticField != 0 && fRequireMagneticField * fMagField < 0.) {
    PostData(1, fList);
    return;
  }
  
  fCentrality->Fill(centrality);
  if (Flatten(centrality)) {
    PostData(1, fList);
    return;
  }
  fCentralityClasses->Fill(centrality);
  fFlattenedCentrality->Fill(centrality);
  
  if (fParticle == AliPID::kUnknown) {
    ::Error("AliAnalysisTaskNucleiYield::UserExec", "No particle type set");
    PostData(1, fList);
    return;
  }
  
  TClonesArray *stack = 0x0;
  if (fIsMC) {
    // get branch "mcparticles"
    stack = (TClonesArray*)ev->GetList()->FindObject(AliAODMCParticle::StdBranchName());
    if (!stack) {
      PostData(1, fList);
      return;
    }
    
    /// Making the list of deuterons in acceptance
    fMmc.SetOwner(kFALSE);
    fAmc.SetOwner(kFALSE);
    for (int iMC = 0; iMC < stack->GetEntriesFast(); ++iMC) {
      AliAODMCParticle *part = (AliAODMCParticle*)stack->UncheckedAt(iMC);
      const int pdg = part->GetPdgCode();
      if (pdg == fPDG) fProduction->Fill(part->P());
      else if (pdg == -fPDG) fProduction->Fill(-part->P());
      if (part->Y() > fRequireYmax || part->Y() < fRequireYmin) continue;
      if (pdg == fPDG) {
        if (part->IsPhysicalPrimary()) fMTotal->Fill(centrality,part->Pt());
        fMmc.Add(part);
      } else if (pdg == -fPDG) {
        if (part->IsPhysicalPrimary()) fATotal->Fill(centrality,part->Pt());
        fAmc.Add(part);
      } 
    }
  }
  
  /// Checking how many deuterons in acceptance are reconstructed well
  for (Int_t iT = 0; iT < (Int_t)ev->GetNumberOfTracks(); ++iT) {
    AliAODTrack *track = dynamic_cast<AliAODTrack*>(ev->GetTrack(iT));
    if (track->GetID() <= 0) continue;
    Double_t dca[2];
    if (!AcceptTrack(track,dca)) continue;
    const float beta = HasTOF(track);
    float pT = track->Pt();
    PtCorrection(pT,track->Charge() > 0);
    if (fIsMC) {
      AliAODMCParticle *part = (AliAODMCParticle*)stack->At(TMath::Abs(track->GetLabel()));
      if (!part) continue;
      if (part->GetPdgCode() == fPDG) {
        if (part->IsPhysicalPrimary()) {
          fMITS_TPC->Fill(centrality,pT);
          fMDCAPrimaryTPC->Fill(centrality,part->Pt(),dca[0]);
        } else fMDCASecondaryTPC->Fill(centrality,part->Pt(),dca[0]);
        fMPtCorrection->Fill(pT,part->Pt() - pT);
        if (beta > 0) {
          if (part->IsPhysicalPrimary()) {
            fMDCAPrimaryTOF->Fill(centrality,part->Pt(),dca[0]);
            fMITS_TPC_TOF->Fill(centrality,part->Pt());
          } else {
            fMDCASecondaryTOF->Fill(centrality,part->Pt(),dca[0]);
          }
        }
      } else if (part->GetPdgCode() == -fPDG) {
        fAPtCorrection->Fill(pT,part->Pt() - pT);
        if (part->IsPhysicalPrimary()) {
          fAITS_TPC->Fill(centrality,part->Pt());
          if (beta > 0) fAITS_TPC_TOF->Fill(centrality,part->Pt());
        }
      }
    } else {
      if (!PassesPIDSelection(track)) continue;
      if (track->Charge() > 0) {
        fMDCAxyTPC->Fill(centrality, pT, dca[0]);
        fMDCAzTPC->Fill(centrality, pT, dca[1]);
        fMTPCcounts->Fill(centrality, pT);
        fMTPCphiCounts->Fill(track->Phi(),pT);
      } else {
        fATPCphiCounts->Fill(track->Phi(),pT);
        fATPCcounts->Fill(centrality, pT);
      }
      if (beta < 0) continue;
      /// \f$ m = \frac{p}{\beta\gamma} \f$
      const float m = track->GetTPCmomentum() * (TMath::Sqrt(1.f - beta * beta) / beta);
      if (track->Charge() > 0) {
        fMDCAxyTOF->Fill(centrality, pT, dca[0]);
        fMDCAzTOF->Fill(centrality, pT, dca[1]);
        fMTOFsignal->Fill(centrality, pT, m * m - fPDGMassOverZ * fPDGMassOverZ);
        fMTOFphiSignal->Fill(track->Phi(),pT,m * m - fPDGMassOverZ * fPDGMassOverZ);
      } else {
        fATOFsignal->Fill(centrality, pT, m * m - fPDGMassOverZ * fPDGMassOverZ);
        fATOFphiSignal->Fill(track->Phi(),pT,m * m - fPDGMassOverZ * fPDGMassOverZ);
      }
    }
    
  } // End AOD track loop
  
  fAmc.Clear();
  fMmc.Clear();
  
  //  Post output data.
  PostData(1,fList);
}

/// Merge the output. Called once at the end of the query.
///
/// \return void
///
void AliAnalysisTaskNucleiYield::Terminate(Option_t *) {
  return;
}

/// This function checks whether a track passes the cuts required in this task
///
/// \param track Track that is going to be checked
/// \param dca[2] Projections on the transverse plane and on z of the distance of closest approach
///               of the track to the primary vertex
/// \return Boolean value: true means that the track has passed all the cuts.
///
Bool_t AliAnalysisTaskNucleiYield::AcceptTrack(AliAODTrack *track, Double_t dca[2]) {
  ULong_t status = track->GetStatus();
  if (!(status & AliVTrack::kITSrefit) && fRequireITSrefit) return kFALSE;
  if (!(status & AliVTrack::kTPCrefit) && fRequireTPCrefit) return kFALSE;
  if (track->Eta() < fRequireEtaMin || track->Eta() > fRequireEtaMax) return kFALSE;
  if (track->Y(fPDGMass) < fRequireYmin || track->Y(fPDGMass) > fRequireYmax) return kFALSE;
  AliAODVertex *vtx1 = (AliAODVertex*)track->GetProdVertex();
  if(Int_t(vtx1->GetType()) == AliAODVertex::kKink && fRequireNoKinks) return kFALSE;
  unsigned int nSPD = 0, nITS = 0;
  for (int i = 0; i < 6; ++i) {
    if (track->HasPointOnITSLayer(i)) {
      if(i < 2) nSPD++;
      nITS++;
    }
  }
  if (track->GetTPCsignal() < fRequireMinEnergyLoss) return kFALSE;
  if (nITS < fRequireITSrecPoints) return kFALSE;
  if (nSPD < fRequireSPDrecPoints) return kFALSE;
  if (track->Chi2perNDF() > fRequireMaxChi2) return kFALSE;
  Double_t cov[3];
  if (!track->PropagateToDCA(fPrimaryVertex, fMagField, 100, dca, cov)) return kFALSE;
  if (TMath::Abs(dca[0]) > fRequireMaxDCAxy) return kFALSE;
  if (TMath::Abs(dca[1]) > fRequireMaxDCAz) return kFALSE;
  return kTRUE;
}

/// This function checks whether a track has or has not a prolongation in TOF.
///
/// \param track Track that has to be checked
/// \return \f$\beta\f$ of the particle, -1 means that there is no correct prolongation in TOF.
///
Float_t AliAnalysisTaskNucleiYield::HasTOF(AliAODTrack *track) {
  Bool_t hasTOFout  = track->GetStatus() & AliVTrack::kTOFout;
  Bool_t hasTOFtime = track->GetStatus() & AliVTrack::kTIME;
  const float len = track->GetIntegratedLength();
  Bool_t hasTOF = Bool_t(hasTOFout & hasTOFtime) && len > 350.f;
  
  if (!hasTOF) return -1.;
  const float p = track->GetTPCmomentum();
  const float tim = track->GetTOFsignal() - fPID->GetTOFResponse().GetStartTime(p);
  if (tim < len / LIGHT_SPEED) return -1.;
  else {
    const float beta = len / (tim * LIGHT_SPEED);
    if (beta < EPS || beta > (1. - EPS))
      return -1.f;
    else
      return beta;
  }
}

/// This functions sets the centrality bins used in the analysis
///
/// \param nbins Number of centrality bins
/// \param bins Array with nbins + 1 elements contanining the edges of the bins
/// \return void
///
void AliAnalysisTaskNucleiYield::SetCentBins(Int_t nbins, Float_t *bins) {
  fCentBins.Set(nbins + 1, bins);
}

/// This functions sets the \f$\mathrm{DCA}_{xy}\f$ bins used in the analysis
///
/// \param nbins Number of \f$\mathrm{DCA}_{xy}\f$ bins
/// \param bins Array with nbins + 1 elements contanining the edges of the bins
/// \return void
///
void AliAnalysisTaskNucleiYield::SetDCABins(Int_t nbins, Float_t min, Float_t max) {
  const float delta = (max - min) / nbins;
  fDCABins.Set(nbins + 1);
  for (int iB = 0; iB < nbins; ++iB) {
    fDCABins[iB] = min + iB * delta;
  }
  fDCABins[nbins] = max;
}

/// This functions sets the \f$p_{\mathrm{T}}\f$ bins used in the analysis
///
/// \param nbins Number of \f$p_{\mathrm{T}}\f$ bins
/// \param bins Array with nbins + 1 elements contanining the edges of the bins
/// \return void
///
void AliAnalysisTaskNucleiYield::SetPtBins(Int_t nbins, Float_t *bins) {
  fPtBins.Set(nbins + 1, bins);
}

/// This function allows to set a custom parametrisation for the TPC response function
///
/// \param par Array of 5 values corresponding to the Bethe Bloch parametrisation
/// \param sigma Sigma of the parametrisation
/// \return void
///
void AliAnalysisTaskNucleiYield::SetCustomTPCpid(Float_t *par, Float_t sigma) {
  if (par == 0x0 && sigma <= 0) {
    fCustomTPCpid.Set(1);
  } else {
    fCustomTPCpid.Set(6);
    for (int i = 0; i < 5; ++i)
      fCustomTPCpid.AddAt(par[i],i);
    fCustomTPCpid.AddAt(sigma, 5);
  }
}

/// This function checks if the track passes the PID selection
///
/// \param t Track to be tested
/// \param sigmas Number of sigmas
/// \return Boolean value: true means that the track passes the PID selection
///
Bool_t AliAnalysisTaskNucleiYield::PassesPIDSelection(AliAODTrack *t) {
  if (fCustomTPCpid.GetSize() < 6 || fIsMC) {
    bool itsPID = kTRUE;
    AliITSPIDResponse &itsPidResp = fPID->GetITSResponse();
    AliTPCPIDResponse &tpcPidResp = fPID->GetTPCResponse();
    if (fRequireITSpidSigmas > 0) {
      itsPID = TMath::Abs(itsPidResp.GetNumberOfSigmas(t, fParticle)) < fRequireITSpidSigmas;
    }
    return itsPID && TMath::Abs(tpcPidResp.GetNumberOfSigmas(t, fParticle)) < fRequireTPCpidSigmas;
  } else {
    const float p = t->GetTPCmomentum() / fPDGMassOverZ;
    const float r = AliExternalTrackParam::BetheBlochAleph(p, fCustomTPCpid[0], fCustomTPCpid[1],
                                                           fCustomTPCpid[2], fCustomTPCpid[3],
                                                           fCustomTPCpid[4]);
    return TMath::Abs(t->GetTPCsignal() - r) < fRequireTPCpidSigmas * fCustomTPCpid[5] * r;
  }
  
}

/// This function sets the number of TOF bins and the boundaries of the histograms
///
/// \param nbins Number of bins
/// \param min Lower boundary of the histogram
/// \param max Higher boundary of the histogram
/// \return void
///
void AliAnalysisTaskNucleiYield::SetTOFBins(Int_t nbins, Float_t min, Float_t max) {
  fTOFnBins = nbins;
  fTOFlowBoundary = min;
  fTOFhighBoundary = max;
}

/// This function sets the number of DCA\f$_{z}\f$ bins and the boundaries of the histogram
///
/// \param nbins Number of bins
/// \param limit Boundaries of the histogram (symmetrical with respect to zero)
/// \return void
///
void AliAnalysisTaskNucleiYield::SetDCAzBins(Int_t nbins, Float_t limit) {
  fDCAzNbins = nbins;
  fDCAzLimit = limit;
}

/// This function sets the particle type to be analysed
///
/// \param part Particle type
/// \return void
///
void AliAnalysisTaskNucleiYield::SetParticleType(AliPID::EParticleType part) {
  fParticle = part;
  fPDGMass = AliPID::ParticleMass(part);
  fPDGMassOverZ = AliPID::ParticleMassZ(part);
}

/// This function provides the flattening of the centrality distribution.
/// Please check hardcoded values! It is better to provide those number by yourself: the probability
/// is computed as \f[\mathrm{Probability}=\frac{C_{i}}{C_{ref}} \f] where \f$C_{i}\f$ is the
/// centrality in the bin _i_ and \f$C_{ref}\f$ is the centrality of the reference bin (i.e. the
/// value around you want the centrality to fluctuate).
///
/// \param cent Event centrality
/// \return Boolean value: true means that the event must be skipped
///
Bool_t AliAnalysisTaskNucleiYield::Flatten(float cent) {
  if (fFlatteningProbs.GetSize() <= 0) {
    Float_t prob[13] = {
      0.839266,0.822364,0.807522,0.804727,0.806675,
      0.828297,0.820842,0.834088,0.861455,1.,
      0.38112,0.661154,0.953928
    };
    fFlatteningProbs.Set(13, prob);
  }
  if (!fIsMC) {
    if (cent >= fFlatteningProbs.GetSize()) return kFALSE;
    else return gRandom->Rndm() > fFlatteningProbs[int(cent)];
  } else {
    // This flattening is required since a strange peak in VOM distribution is observed in MC
    return (cent < 1) && gRandom->Rndm() > 0.35;
  }
}

/// This function provides the correction for wrongly calculated \f$p_{\mathrm{T}}\f$.
///
/// \param pt \f$p_{\mathrm{T}}\f$ of the track
/// \param positiveCharge True if the track has positive sign.
/// \return void
///
void AliAnalysisTaskNucleiYield::PtCorrection(float &pt, bool positiveCharge) {
  Float_t *par = (positiveCharge) ? fPtCorrectionM.GetArray() : fPtCorrectionA.GetArray();
  const Float_t correction = par[0] + par[1] * TMath::Exp(par[2] * pt);
  pt -= correction;
}
