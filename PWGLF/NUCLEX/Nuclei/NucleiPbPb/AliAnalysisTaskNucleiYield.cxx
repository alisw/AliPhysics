#include "AliAnalysisTaskNucleiYield.h"

// ROOT includes
#include <TAxis.h>
#include <TChain.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TF1.h>
#include <TList.h>
#include <TMath.h>
#include <TParticle.h>
#include <TClonesArray.h>
#include <TTree.h>
#include <TRandom3.h>

// ALIROOT includes
#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliPDG.h"
#include "AliMultSelection.h"
#include "AliTPCPIDResponse.h"
#include "AliTOFPIDResponse.h"
#include "AliVTrack.h"
#include "AliVVertex.h"
#include "AliVEvent.h"
#include "AliVParticle.h"
#include "AliMCEvent.h"
#include "AliInputEventHandler.h"
#include "AliVEventHandler.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliAODVertex.h"

#define LIGHT_SPEED 2.99792457999999984e-02 // in the units that TOF likes
#define EPS 1.e-16

using TMath::TwoPi;

///\cond CLASSIMP
ClassImp(AliAnalysisTaskNucleiYield);
///\endcond

const TString kNames[5]= {"pion","kaon","proton","deuteron","triton"};
const AliPID::EParticleType kSpecies[5] = {AliPID::kPion, AliPID::kKaon, AliPID::kProton, AliPID::kDeuteron, AliPID::kTriton};

static double TOFsignal(double *x, double *par) {
  double &norm = par[0];
  double &mean = par[1];
  double &sigma = par[2];
  double &tail = par[3];

  if (x[0] <= (tail + mean))
    return norm * TMath::Gaus(x[0], mean, sigma);
  else
    return norm * TMath::Gaus(tail + mean, mean, sigma) * TMath::Exp(-tail * (x[0] - tail - mean) / (sigma * sigma));
}

/// Method for the correct logarithmic binning of histograms.
///
/// \param h Histogram that has to be correctly binned
///
static void BinLogAxis(const TH1 *h) {
  TAxis *axis = const_cast<TAxis*>(h->GetXaxis());
  const Int_t bins = axis->GetNbins();

  const Double_t from = axis->GetXmin();
  const Double_t to = axis->GetXmax();
  Double_t *newBins = new Double_t[bins + 1];

  newBins[0] = from;
  Double_t factor = pow(to / from, 1. / bins);

  for (Int_t i = 1; i <= bins; i++) {
    newBins[i] = factor * newBins[i - 1];
  }
  axis->Set(bins, newBins);
  delete [] newBins;
}

/// Standard and default constructor of the class.
///
/// \param taskname Name of the task
/// \param partname Name of the analysed particle
///
AliAnalysisTaskNucleiYield::AliAnalysisTaskNucleiYield(TString taskname)
  :AliAnalysisTaskSE(taskname.Data())
   ,fEventCut{false}
   ,fFilterBit{BIT(4)}
   ,fPropagateTracks{true}
   ,fList{nullptr}
   ,fCutVec{}
   ,fPDG{0}
   ,fPDGMass{0}
   ,fPDGMassOverZ{0}
   ,fCharge{1.f}
   ,fIsMC{false}
   ,fPID{nullptr}
   ,fMagField{0.f}
   ,fDCAzLimit{10.}
   ,fDCAzNbins{400}
   ,fPtCorrectionA{3}
   ,fPtCorrectionM{3}
   ,fTOFlowBoundary{-2.4}
   ,fTOFhighBoundary{3.6}
   ,fTOFnBins{75}
   ,fDisableITSatHighPt{100.f}
   ,fDisableTPCpidAtHighPt{100.f}
   ,fEnablePtCorrection{false}
   ,fRequireITSrefit{true}
   ,fRequireTPCrefit{true}
   ,fRequireNoKinks{true}
   ,fRequireITSrecPoints{2u}
   ,fRequireITSsignal{0u}
   ,fRequireSDDrecPoints{0u}
   ,fRequireSPDrecPoints{1u}
   ,fRequireTPCsignal{70u}
   ,fRequireEtaMin{-0.8f}
   ,fRequireEtaMax{0.8f}
   ,fRequireYmin{-0.5f}
   ,fRequireYmax{0.5f}
   ,fRequireMaxChi2{4.f}
   ,fRequireMaxDCAxy{0.12f}
   ,fRequireMaxDCAz{1.f}
   ,fRequireTPCpidSigmas{3.f}
   ,fRequireITSpidSigmas{-1.f}
   ,fRequireTOFpidSigmas{-1.f}
   ,fRequireMinEnergyLoss{0.}
   ,fRequireVetoSPD{false}
   ,fRequireMaxMomentum{-1.}
   ,fFixForLHC14a6{true}
   ,fParticle{AliPID::kUnknown}
   ,fCentBins{0}
   ,fDCABins{0}
   ,fPtBins{0}
   ,fCustomTPCpid{0}
   ,fFlatteningProbs{0}
   ,fFlattenedCentrality{nullptr}
   ,fCentralityClasses{nullptr}
   ,fProduction{nullptr}
   ,fReconstructed{{nullptr}}
   ,fTotal{nullptr}
   ,fPtCorrection{nullptr}
   ,fDCAPrimary{{nullptr}}
   ,fDCASecondary{{nullptr}}
   ,fDCASecondaryWeak{{nullptr}}
   ,fTOFsignal{nullptr}
   ,fTPCcounts{nullptr}
   ,fDCAxy{{nullptr}}
   ,fDCAz{{nullptr}}
   ,fTOFtemplates{nullptr}
   ,fEnableFlattening{true} {
     gRandom->SetSeed(0); //TODO: provide a simple method to avoid "complete randomness"
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
  if (fList) delete fList;
  if (fTOFfunction) delete fTOFfunction;
}

/// This function creates all the histograms and all the objects in general used during the analysis
/// \return void
///
void AliAnalysisTaskNucleiYield::UserCreateOutputObjects() {

  fList = new TList();
  fList->SetOwner(true);

  const Int_t nPtBins = fPtBins.GetSize() - 1;
  const Int_t nCentBins = fCentBins.GetSize() - 1;
  const Int_t nDCAbins = fDCABins.GetSize() - 1;
  const Float_t *pTbins = fPtBins.GetArray();
  const Float_t *centBins = fCentBins.GetArray();
  const Float_t *dcaBins = fDCABins.GetArray();

  fCentralityClasses = new TH1F("fCentralityClasses",";Centrality classes(%);Events / Class;",
      nCentBins,centBins);
  fFlattenedCentrality = new TH1F("fFlattenCentrality","After the flattening;Centrality (%); \
      Events / 1%;",100,0.,100.);
  fList->Add(fCentralityClasses);
  fList->Add(fFlattenedCentrality);

  char   letter[2] = {'A','M'};
  string tpctof[2] = {"TPC","TOF"};
  string tpctofMC[2] = {"TPC","TPC_TOF"};

  if (fIsMC) {
    fProduction = new TH1F("fProduction",";#it{p} (GeV/#it{c});Entries",100,-10,10);
    fList->Add(fProduction);

    for (int iC = 0; iC < 2; ++iC) {
      fTotal[iC] = new TH2F(Form("f%cTotal",letter[iC]),";Centrality (%);#it{p}_{T} (GeV/#it{c}); Counts",
          nCentBins,centBins,nPtBins,pTbins);
      fPtCorrection[iC] = new TH2F(Form("f%cPtCorrection",letter[iC]),
          ";#it{p}_{T}^{rec} (GeV/#it{c});#it{p}_{T}^{MC}-#it{p}_{T}^{rec} (GeV/#it{c});Entries",
          160,0.4,6.,80,-1.,1.);
      fList->Add(fTotal[iC]);
      fList->Add(fPtCorrection[iC]);
      for (int iT = 0; iT < 2; ++iT) {
        fReconstructed[iT][iC] = new TH2F(Form("f%cITS_%s",letter[iC],tpctofMC[iT].data()),";Centrality (%);#it{p}_{T} (GeV/#it{c}); Counts",
            nCentBins,centBins,nPtBins,pTbins);
        fDCAPrimary[iT][iC] = new TH3F(Form("f%cDCAPrimary%s",letter[iC],tpctof[iT].data()),";Centrality (%);#it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)",
            nCentBins,centBins,nPtBins,pTbins,nDCAbins,dcaBins);
        fDCASecondary[iT][iC] = new TH3F(Form("f%cDCASecondary%s",letter[iC],tpctof[iT].data()),";Centrality (%);#it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)",
            nCentBins,centBins,nPtBins,pTbins,nDCAbins,dcaBins);
        fDCASecondaryWeak[iT][iC] = new TH3F(Form("f%cDCASecondaryWeak%s",letter[iC],tpctof[iT].data()),";Centrality (%);#it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)",
            nCentBins,centBins,nPtBins,pTbins,nDCAbins,dcaBins);
        fList->Add(fReconstructed[iT][iC]);
        fList->Add(fDCAPrimary[iT][iC]);
        fList->Add(fDCASecondary[iT][iC]);
        fList->Add(fDCASecondaryWeak[iT][iC]);
      }
    }
  } else {

    float tofBins[fTOFnBins + 1];
    const float deltaTOF = (fTOFhighBoundary - fTOFlowBoundary) / fTOFnBins;
    for (int i = 0; i <= fTOFnBins; ++i)
      tofBins[i] = i * deltaTOF + fTOFlowBoundary;
    float dcazBins[fDCAzNbins + 1];
    const float deltaDCAz = 2.f * fDCAzLimit / fDCAzNbins;
    for (int i = 0; i <= fDCAzNbins; ++i)
      dcazBins[i] = i * deltaDCAz - fDCAzLimit;
    const int nSigmaBins = 240;
    float sigmaBins[nSigmaBins + 1];
    for (int i = 0; i <= nSigmaBins; ++i)
      sigmaBins[i] = -6.f + i * 0.05;

    for (int iC = 0; iC < 2; ++iC) {
      fTOFsignal[iC] = new TH3F(Form("f%cTOFsignal",letter[iC]),
          ";Centrality (%);#it{p}_{T} (GeV/#it{c});#it{m}^{2}-m_{PDG}^{2} (GeV/#it{c}^{2})^{2}",
          nCentBins,centBins,nPtBins,pTbins,fTOFnBins,tofBins);
      fTPCcounts[iC] = new TH3F(Form("f%cTPCcounts",letter[iC]),";Centrality (%);#it{p}_{T} (GeV/#it{c}); n_{#sigma} d",
          nCentBins,centBins,nPtBins,pTbins,nSigmaBins,sigmaBins);

      fList->Add(fTOFsignal[iC]);
      fList->Add(fTPCcounts[iC]);

      for (int iT = 0; iT < 2; ++iT) {
        fDCAxy[iT][iC] = new TH3F(Form("f%cDCAxy%s",letter[iC],tpctof[iT].data()),";Centrality (%);#it{p}_{T} (GeV/#it[c}); DCA_{xy} (cm)",
            nCentBins,centBins,nPtBins,pTbins,nDCAbins,dcaBins);
        fDCAz[iT][iC] = new TH3F(Form("f%cDCAz%s",letter[iC],tpctof[iT].data()),";Centrality (%);#it{p}_{T} (GeV/#it{c}); DCA_{z} (cm)",
            nCentBins,centBins,nPtBins,pTbins,fDCAzNbins,dcazBins);
        fList->Add(fDCAxy[iT][iC]);
        fList->Add(fDCAz[iT][iC]);
      }
    }

    for (int iS = 0; iS < 5; ++iS) {
      fTOFtemplates[iS] = new TH3F(Form("fTOFtemplates%i",iS),
          Form("%s;Centrality (%%);#it{p}_{T} (GeV/#it{c});#it{m}^{2}-m_{PDG}^{2} (GeV/#it{c}^{2})^{2}",kNames[iS].Data()),
          nCentBins,centBins,nPtBins,pTbins,fTOFnBins,tofBins);
      fList->Add(fTOFtemplates[iS]);
    }
  }

  fTOFfunction = new TF1("fTOFfunction", TOFsignal, -2440., 2440., 4);
  if (fTOFfunctionPars.GetSize() == 4)
    fTOFfunction->SetParameters(fTOFfunctionPars.GetArray());

  AliPDG::AddParticlesToPdgDataBase();
  PostData(1,fList);
}

/// This is the function that is evaluated for each event. The analysis code stays here.
///
/// \param options Deprecated parameter
/// \return void
///
void AliAnalysisTaskNucleiYield::UserExec(Option_t *){
  /// The first check performed is on the particle type requested. If this type is AliPID::Unknown -
  /// the default one - the task does not process the information.
  if (fParticle == AliPID::kUnknown) {
    ::Error("AliAnalysisTaskNucleiYield::UserExec", "No particle type set");
    PostData(1, fList);
    return;
  }

  AliVEvent *ev = InputEvent();
  if (!fEventCut.AcceptEvent(ev)) {
    PostData(1, fList);
    return;
  }

  /// To perform the majority of the analysis - and also this one - the standard PID handler is
  /// required.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* handl = (AliInputEventHandler*)mgr->GetInputEventHandler();
  fPID = handl->GetPIDResponse();


  /// The centrality selection in PbPb uses the percentile determined with V0.
  float centrality = fEventCut.GetCentrality();

  /// The magnetic field
  fMagField = ev->GetMagneticField();

  /// At the stage of event selection the Flattening is applied. This technique makes flat the
  /// centrality distribution using a pseudo-random selection based on prior probabilities.
  /// A complete description of this technique is present in the documentation of the Flatten
  /// function.

  if (Flatten(centrality) && fEnableFlattening) {
    PostData(1, fList);
    return;
  }
  fCentralityClasses->Fill(centrality);
  fFlattenedCentrality->Fill(centrality);

  TClonesArray *stack = nullptr;
  if (fIsMC) {
    // get branch "mcparticles"
    stack = (TClonesArray*)ev->GetList()->FindObject(AliAODMCParticle::StdBranchName());
    if (!stack)
      ::Fatal("AliAnalysisTaskNucleiYield::UserExec","MC analysis requested on a sample without the MC particle array.");

    /// Making the list of the deuterons we want to measure
    for (int iMC = 0; iMC < stack->GetEntriesFast(); ++iMC) {
      AliAODMCParticle *part = (AliAODMCParticle*)stack->UncheckedAt(iMC);
      const int pdg = std::abs(part->GetPdgCode());
      const int iC = part->Charge() > 0 ? 1 : 0;
      const int mult = -1 + 2 * iC;
      if (pdg != fPDG) continue;
      fProduction->Fill(mult * part->P());
      if (part->Y() > fRequireYmax || part->Y() < fRequireYmin) continue;
      if (part->IsPhysicalPrimary()) fTotal[iC]->Fill(centrality,part->Pt());
    }
  }

  /// Checking how many deuterons in acceptance are reconstructed well
  for (Int_t iT = 0; iT < (Int_t)ev->GetNumberOfTracks(); ++iT) {
    AliAODTrack *track = dynamic_cast<AliAODTrack*>(ev->GetTrack(iT));

    if (track->GetID() <= 0) continue;
    Double_t dca[2] = {0.};
    if (!track->TestFilterBit(fFilterBit) && fFilterBit) continue;
    if (!AcceptTrack(track,dca)) continue;
    const float beta = HasTOF(track,fPID);
    const int iTof = beta > EPS ? 1 : 0;
    float pT = track->Pt() * fCharge;
    if (fEnablePtCorrection) PtCorrection(pT,track->Charge() > 0);

    if (fIsMC) {
      AliAODMCParticle *part = (AliAODMCParticle*)stack->At(TMath::Abs(track->GetLabel()));
      /// Workaround: if the AOD are filtered with an AliRoot tag before v5-08-18, hyper-nuclei prongs
      /// are marked as SecondaryFromMaterial.
      const int mother_id = part->GetMother();
      AliAODMCParticle* mother = (mother_id >= 0) ? (AliAODMCParticle*)stack->At(mother_id) : 0x0;
      const int mother_pdg = mother ? TMath::Abs(mother->GetPdgCode()) : 0;
      const bool isFromHyperNucleus = (mother_pdg > 1000000000 && (mother_pdg / 10000000) % 10 != 0);
      if (!part) continue;
      const int iC = part->Charge() > 0 ? 1 : 0;
      if (std::abs(part->GetPdgCode()) == fPDG) {
        for (int iR = iTof; iR >= 0; iR--) {
          if (part->IsPhysicalPrimary()) {
            if (TMath::Abs(dca[0]) <= fRequireMaxDCAxy) fReconstructed[iR][iC]->Fill(centrality,pT);
            fDCAPrimary[iR][iC]->Fill(centrality,pT,dca[0]);
            if (!iR) fPtCorrection[iC]->Fill(pT,part->Pt()-pT); // Fill it only once.
          } else if (part->IsSecondaryFromMaterial() && !isFromHyperNucleus)
            fDCASecondary[iR][iC]->Fill(centrality,pT,dca[0]);
          else
            fDCASecondaryWeak[iR][iC]->Fill(centrality,pT,dca[0]);
        }
      }
    } else {
      bool pid_check = PassesPIDSelection(track);
      const int iC = (track->Charge() > 0) ? 1 : 0;

      float tpc_n_sigma = GetTPCsigmas(track);
      float tof_n_sigma = iTof ? fPID->NumberOfSigmas(AliPIDResponse::kTOF, track, fParticle) : -999.f;

      for (int iR = iTof; iR >= 0; iR--) {
        /// TPC asymmetric cut to avoid contamination from protons in the DCA distributions. TOF sigma cut is set to 4
        /// to compensate for the shift in the sigma (to be rechecked in case of update of TOF PID response)
        if (tpc_n_sigma > -2. && tpc_n_sigma < 3. && (fabs(tof_n_sigma) < 4. || !iTof)) {
          fDCAxy[iR][iC]->Fill(centrality, pT, dca[0]);
          fDCAz[iR][iC]->Fill(centrality, pT, dca[1]);
        }
      }
      if (TMath::Abs(dca[0]) > fRequireMaxDCAxy) continue;
      fTPCcounts[iC]->Fill(centrality, pT, tpc_n_sigma);
      if (!pid_check || iTof == 0) continue;
      /// \f$ m = \frac{p}{\beta\gamma} \f$
      const float m2 = track->P() * track->P() * (1.f / (beta * beta) - 1.f);
      fTOFsignal[iC]->Fill(centrality, pT, m2 - fPDGMassOverZ * fPDGMassOverZ);

      if (fTOFfunctionPars.GetSize() == 4) {
        AliTOFPIDResponse& tofPid = fPID->GetTOFResponse();
        for (int iS = 0; iS < 5; ++iS) {
          const float expt0 = tofPid.GetExpectedSignal(track,kSpecies[iS]);
          const float sigma = tofPid.GetExpectedSigma(track->P(),expt0,kSpecies[iS]);
          const float smearing = fTOFfunction->GetRandom() + gRandom->Gaus(0., sqrt(sigma * sigma - tofPid.GetTimeResolution() * tofPid.GetTimeResolution()));
          const float expBeta = track->GetIntegratedLength() / ((expt0 + smearing) * LIGHT_SPEED);
          if (expBeta > EPS) {
            const float expM = track->P() * track->P() * (1.f / (expBeta * expBeta) - 1.f);
            fTOFtemplates[iS]->Fill(centrality,pT,expM - fPDGMassOverZ * fPDGMassOverZ);
          }
        }
      }
    }

  } // End AOD track loop

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
bool AliAnalysisTaskNucleiYield::AcceptTrack(AliAODTrack *track, Double_t dca[2]) {
  ULong_t status = track->GetStatus();
  fCutVec.SetPtEtaPhiM(track->Pt() * fCharge, track->Eta(), track->Phi(), fPDGMass);
  if (!(status & AliVTrack::kTPCrefit) && fRequireTPCrefit) return false;
  if (track->Eta() < fRequireEtaMin || track->Eta() > fRequireEtaMax) return false;
  if (fCutVec.Rapidity() < fRequireYmin || fCutVec.Rapidity() > fRequireYmax) return false;
  AliAODVertex *vtx1 = (AliAODVertex*)track->GetProdVertex();
  if(Int_t(vtx1->GetType()) == AliAODVertex::kKink && fRequireNoKinks) return false;
  if (track->Chi2perNDF() > fRequireMaxChi2) return false;
  if (track->GetTPCsignalN() < fRequireTPCsignal) return false;
  if (track->GetTPCsignal() < fRequireMinEnergyLoss) return false;
  if (fRequireMaxMomentum > 0 && track->P() > fRequireMaxMomentum) return false;

  /// ITS related cuts
  dca[0] = 0.;
  dca[1] = 0.;
  if (track->Pt() < fDisableITSatHighPt) {
    unsigned int nSPD = 0u, nSDD = 0u, nSSD = 0u;
    int nITS = GetNumberOfITSclustersPerLayer(track, nSPD, nSDD, nSSD);
    if (!(status & AliVTrack::kITSrefit) && fRequireITSrefit) return false;
    if (nITS < fRequireITSrecPoints) return false;
    if (nSPD < fRequireSPDrecPoints) return false;
    if (nSDD < fRequireSDDrecPoints) return false;
    if (fRequireVetoSPD && nSPD > 0) return false;
    Double_t cov[3];
    if (fPropagateTracks)
      if (!track->PropagateToDCA(fEventCut.GetPrimaryVertex(), fMagField, 100, dca, cov)) return false;
    if (TMath::Abs(dca[1]) > fRequireMaxDCAz) return false;
    //if (TMath::Abs(dca[0]) > fRequireMaxDCAxy) return false;
  }

  return true;
}

/// This function checks whether a track has or has not a prolongation in TOF.
///
/// \param track Track that has to be checked
/// \return \f$\beta\f$ of the particle, -1 means that there is no correct prolongation in TOF.
///
float AliAnalysisTaskNucleiYield::HasTOF(AliAODTrack *track, AliPIDResponse *pid) {
  bool hasTOFout  = track->GetStatus() & AliVTrack::kTOFout;
  bool hasTOFtime = track->GetStatus() & AliVTrack::kTIME;
  const float len = track->GetIntegratedLength();
  bool hasTOF = Bool_t(hasTOFout & hasTOFtime) && (len > 350.);

  if (!hasTOF) return -1.;
  const float p = track->GetTPCmomentum();
  const float tim = track->GetTOFsignal() - pid->GetTOFResponse().GetStartTime(p);
  if (tim < len / LIGHT_SPEED) return -1.;
  else {
    const float beta = len / (tim * LIGHT_SPEED);
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

/// This functions sets the \f$\mathrm{DCA}_{xy}\f$ bins used in the analysis
///
/// \param nbins Number of \f$\mathrm{DCA}_{xy}\f$ bins
/// \param bins Array with nbins + 1 elements contanining the edges of the bins
/// \return void
///
void AliAnalysisTaskNucleiYield::SetDCABins(Int_t nbins, Float_t *bins) {
  fDCABins.Set(nbins + 1, bins);
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

float AliAnalysisTaskNucleiYield::GetTPCsigmas(AliVTrack* t) {
  if (fCustomTPCpid.GetSize() < 6 || fIsMC) {
    AliTPCPIDResponse &tpcPidResp = fPID->GetTPCResponse();
    return tpcPidResp.GetNumberOfSigmas(t, fParticle);
  } else {
    const float p = t->GetTPCmomentum() / fPDGMassOverZ;
    const float r = fCharge * fCharge * AliExternalTrackParam::BetheBlochAleph(p, fCustomTPCpid[0], fCustomTPCpid[1],
        fCustomTPCpid[2], fCustomTPCpid[3],
        fCustomTPCpid[4]);
    return (t->GetTPCsignal() - r) / (fCustomTPCpid[5] * r);
  }


}

/// This function checks if the track passes the PID selection
///
/// \param t Track to be tested
/// \param sigmas Number of sigmas
/// \return Boolean value: true means that the track passes the PID selection
///
bool AliAnalysisTaskNucleiYield::PassesPIDSelection(AliAODTrack *t) {
  bool tofPID = true, itsPID = true, tpcPID = true;
  if (fRequireITSpidSigmas > 0 && t->Pt() < fDisableITSatHighPt) {
    AliITSPIDResponse &itsPidResp = fPID->GetITSResponse();
    itsPID = TMath::Abs(itsPidResp.GetNumberOfSigmas(t, fParticle)) < fRequireITSpidSigmas;
  }

  if (fRequireTOFpidSigmas > 0) {
    tofPID = TMath::Abs(fPID->NumberOfSigmasTOF(t, fParticle)) < fRequireTOFpidSigmas;
  }

  if (t->Pt() < fDisableTPCpidAtHighPt) {
    if (fCustomTPCpid.GetSize() < 6 || fIsMC) {
      AliTPCPIDResponse &tpcPidResp = fPID->GetTPCResponse();
      tpcPID = TMath::Abs(tpcPidResp.GetNumberOfSigmas(t, fParticle)) < fRequireTPCpidSigmas;
    } else {
      const float p = t->GetTPCmomentum() / fPDGMassOverZ;
      const float r = AliExternalTrackParam::BetheBlochAleph(p, fCustomTPCpid[0], fCustomTPCpid[1],
          fCustomTPCpid[2], fCustomTPCpid[3],
          fCustomTPCpid[4]);
      tpcPID = TMath::Abs(t->GetTPCsignal() - r) < fRequireTPCpidSigmas * fCustomTPCpid[5] * r;
    }
  }

  return itsPID && tpcPID && tofPID;
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
  fCharge = TMath::Abs(AliPID::ParticleCharge(fParticle));
}

/// This function provides the flattening of the centrality distribution.
/// Please check the hardcoded values! It is better to provide those number by yourself: the
/// probability is computed as \f[\mathrm{Probability}=\frac{C_{i}}{C_{ref}} \f] where \f$C_{i}\f$
/// is the centrality in the bin _i_ and \f$C_{ref}\f$ is the centrality of the reference bin
/// (namely the value around you want the centrality to fluctuate).
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
    if (cent >= fFlatteningProbs.GetSize()) return false;
    else return gRandom->Rndm() > fFlatteningProbs[int(cent)];
  } else {
    // This flattening is required since a strange peak in VOM distribution is observed in MC
    if (fFixForLHC14a6) {
      if (cent < 1.5e-3f) {
        return true;
      } else if (cent < 0.05 && gRandom->Rndm() < 0.5) {
        return true;
      }
    }
  }
  return false;
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

/// This function returns the number of clusters for each ITS subdetector
///
/// \param track
/// \param nSPD number of clusters in SPD
//  \param nSDD number of clusters in SDD
//  \param nSSD number of clusters in SSD
/// \return int number of clusters in ITS
///
int AliAnalysisTaskNucleiYield::GetNumberOfITSclustersPerLayer(AliVTrack *track, unsigned int &nSPD, unsigned int &nSDD, unsigned int &nSSD) {
  if (!track) return -1;
  nSPD = 0u;
  nSDD = 0u;
  nSSD = 0u;
  for (int i = 0; i < 6; ++i) {
    if (track->HasPointOnITSLayer(i)) {
      if (i < 2) nSPD++;
      else if (i < 4) nSDD++;
      else nSSD++;
    }
  }
  return nSPD + nSDD + nSSD;
}
