#include "AliAnalysisTaskNucleiYieldESD.h"

// ROOT includes
#include <TAxis.h>
#include <TChain.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TList.h>
#include <TMath.h>
#include <TParticle.h>
#include <TRandom3.h>

// ALIROOT includes
#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliESDtrack.h"
#include "AliInputEventHandler.h"
#include "AliMCEvent.h"
#include "AliMultSelection.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliMCEvent.h"
#include "AliTOFPIDResponse.h"
#include "AliTPCPIDResponse.h"
#include "AliVEvent.h"
#include "AliVEventHandler.h"
#include "AliVParticle.h"
#include "AliVTrack.h"
#include "AliVVertex.h"

#define LIGHT_SPEED 2.99792457999999984e-02 // in the units that TOF likes
#define EPS 1.e-15

const TString kLetters = "AM";
const TString kNames[5]= {"pion","kaon","proton","deuteron","triton"};
const AliPID::EParticleType kSpecies[5] = {AliPID::kPion, AliPID::kKaon, AliPID::kProton, AliPID::kDeuteron, AliPID::kTriton};

using TMath::TwoPi;

///\cond CLASSIMP
ClassImp(AliAnalysisTaskNucleiYieldESD);
///\endcond

double TOFsignal(double *x, double *par) {
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
  const int bins = axis->GetNbins();

  const Double_t from = axis->GetXmin();
  const Double_t to = axis->GetXmax();
  Double_t *newBins = new Double_t[bins + 1];

  newBins[0] = from;
  Double_t factor = pow(to / from, 1. / bins);

  for (int i = 1; i <= bins; i++) {
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
AliAnalysisTaskNucleiYieldESD::AliAnalysisTaskNucleiYieldESD(TString taskname) : AliAnalysisTaskSE(taskname.Data()),
  fTOFfunction(0x0),   
  fTOFtail(75.f),
  fList(0x0),  
  fYregion(0.5f),
  fTPCnSigmaCut(4.f),
  fParticle(AliPID::kDeuteron),   
  fPDG(AliPID::ParticleCode(AliPID::kDeuteron)),
  fPDGMass(AliPID::ParticleMass(AliPID::kDeuteron)),
  fPDGMassOverZ(AliPID::ParticleMassZ(AliPID::kDeuteron)),  
  fIsMC(false),
  fPID(0x0),
  fMagField(0.f),     
  fDCAzLimit(10.f),     
  fDCAzNbins(400),
  fPtCorrectionA(0),    
  fPtCorrectionM(0),   
  fTOFlowBoundary(-2.4),
  fTOFhighBoundary(3.6),
  fTOFnBins(75),
  fCentBins(),
  fDCABins(),       
  fPtBins(),        
  fFlatteningProbs(),   
  fCentralityClasses(0x0),
  fProduction(),
  fITS_TPC(),   
  fITS_TPC_TOF(),    
  fTotal(),  
  fPtCorrection(),
  fDCAPrimaryTPC(),
  fDCASecondaryTPC(),
  fDCAPrimaryTOF(),
  fDCASecondaryTOF(),
  fTOFsignal(),
  fTPCcounts(),
  fTPCdEdx(),
  fTPCdEdxTpcCut(),
  fTPCdEdxTofCut(),
  fDCAxyTPC(),      
  fDCAzTPC(),     
  fDCAxyTOF(),      
  fDCAzTOF(), 
  fTOFtemplates()
{
  AliESDtrackCuts* tmp = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(false, 0);
  fTrackCuts = *tmp;
  delete tmp;
  fTrackCuts.SetMaxChi2TPCConstrainedGlobal(36);

  fEventCuts.SetupLHC15o();
  
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}

/// Standard destructor
///
AliAnalysisTaskNucleiYieldESD::~AliAnalysisTaskNucleiYieldESD() {
  if (fList) delete fList;
}

/// This function creates all the histograms and all the objects in general used during the analysis
/// \return void
///
void AliAnalysisTaskNucleiYieldESD::UserCreateOutputObjects() {
  gRandom->SetSeed(0);

  fList = new TList();
  fList->SetOwner(kTRUE);

  const int nPtBins = fPtBins.GetSize() - 1;
  const int nCentBins = fCentBins.GetSize() - 1;
  const int nDCAbins = fDCABins.GetSize() - 1;
  const float *pTbins = fPtBins.GetArray();
  const float *centBins = fCentBins.GetArray();
  const float *dcaBins = fDCABins.GetArray();

  fCentralityClasses = new TH1F("fCentralityClasses",";Centrality classes(%%);Events / Class;",nCentBins,centBins);
  fList->Add(fCentralityClasses);

  if (fIsMC) {
    for (int iS = 0; iS < 2; ++iS) {
      fProduction[iS] = new TH1F(Form("f%cProduction",kLetters[iS]),";p (GeV/c);Entries",100,0.,10.);
      fTotal[iS] = new TH2F(Form("f%cTotal",kLetters[iS]),";Centrality (%%);p_{T} (GeV/c); Counts",
          nCentBins,centBins,nPtBins,pTbins);
      fITS_TPC[iS] = new TH2F(Form("f%cITS_TPC",kLetters[iS]),";Centrality (%%);p_{T} (GeV/c); Counts",
          nCentBins,centBins,nPtBins,pTbins);
      fITS_TPC_TOF[iS] = new TH2F(Form("f%cITS_TPC_TOF",kLetters[iS]),";Centrality (%%);p_{T} (GeV/c); Counts",
          nCentBins,centBins,nPtBins,pTbins);
      fDCAPrimaryTPC[iS] = new TH3F(Form("f%cDCAPrimaryTPC",kLetters[iS]),";Centrality (%%);p_{T} (GeV/c); DCA_{xy} (cm)",
          nCentBins,centBins,nPtBins,pTbins,nDCAbins,dcaBins);
      fDCASecondaryTPC[iS] = new TH3F(Form("f%cDCASecondaryTPC",kLetters[iS]),";Centrality (%%);p_{T} (GeV/c); DCA_{xy} (cm)",
          nCentBins,centBins,nPtBins,pTbins,nDCAbins,dcaBins);
      fDCAPrimaryTOF[iS] = new TH3F(Form("f%cDCAPrimaryTOF",kLetters[iS]),";Centrality (%%);p_{T} (GeV/c); DCA_{xy} (cm)",
          nCentBins,centBins,nPtBins,pTbins,nDCAbins,dcaBins);
      fDCASecondaryTOF[iS] = new TH3F(Form("f%cDCASecondaryTOF",kLetters[iS]),";Centrality (%%);p_{T} (GeV/c); DCA_{xy} (cm)",
          nCentBins,centBins,nPtBins,pTbins,nDCAbins,dcaBins);
      fPtCorrection[iS] = new TH2F(Form("f%cPtCorrection",kLetters[iS]),
          ";p_{T}^{rec} (GeV/c);p_{T}^{MC}-p_{T}^{rec} (GeV/c);Entries",
          160,0.4,6.,80,-1.,1.);

      fList->Add(fProduction[iS]);
      fList->Add(fTotal[iS]);
      fList->Add(fITS_TPC[iS]);
      fList->Add(fITS_TPC_TOF[iS]);
      fList->Add(fDCAPrimaryTPC[iS]);
      fList->Add(fDCASecondaryTPC[iS]);
      fList->Add(fDCAPrimaryTOF[iS]);
      fList->Add(fDCASecondaryTOF[iS]);
      fList->Add(fPtCorrection[iS]);
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

    const int nTPCptBins = 30;
    float tpcPtBins[nTPCptBins + 1];
    for (int i = 0; i <= nTPCptBins; ++i) {
      tpcPtBins[i] = 0.5 + 0.05 * i;
    }

    const int nTPCeLossBins = 800;
    float tpcElossBins[nTPCeLossBins + 1];
    for (int i = 0; i <= nTPCeLossBins; ++i) {
      tpcElossBins[i] = i * 2.f;
    }

    const int nSigmaBins = 240;
    float sigmaBins[nSigmaBins + 1];
    for (int i = 0; i <= nSigmaBins; ++i)
      sigmaBins[i] = -6.f + i * 0.05;

    for (int iS = 0; iS < 2; ++iS) {
      fTOFsignal[iS] = new TH3F(Form("f%cTOFsignal",kLetters[iS]),
          ";Centrality (%);p_{T} (GeV/c);m_{TOF}^{2}-m_{PDG}^{2} (GeV/c^{2})^{2}",
          nCentBins,centBins,nPtBins,pTbins,fTOFnBins,tofBins);
      fTPCcounts[iS] = new TH3F(Form("f%cTPCcounts",kLetters[iS]),";Centrality (%%);#it{p}_{T} (GeV/c); n_{#sigma} d",
          nCentBins,centBins,nPtBins,pTbins,nSigmaBins,sigmaBins);
      fDCAxyTPC[iS] = new TH3F(Form("f%cDCAxyTPC",kLetters[iS]),";Centrality (%%);#it{p}_{T} (GeV/c); DCA_{xy} (cm)",
          nCentBins,centBins,nPtBins,pTbins,nDCAbins,dcaBins);
      fDCAzTPC[iS] = new TH3F(Form("f%cDCAzTPC",kLetters[iS]),";Centrality (%%);#it{p}_{T} (GeV/c); DCA_{z} (cm)",
          nCentBins,centBins,nPtBins,pTbins,fDCAzNbins,dcazBins);
      fDCAxyTOF[iS] = new TH3F(Form("f%cDCAxyTOF",kLetters[iS]),";Centrality (%%);#it{p}_{T} (GeV/c); DCA_{xy} (cm)",
          nCentBins,centBins,nPtBins,pTbins,nDCAbins,dcaBins);
      fDCAzTOF[iS] = new TH3F(Form("f%cDCAzTOF",kLetters[iS]),";Centrality (%%);#it{p}_{T} (GeV/c); DCA_{z} (cm)",
          nCentBins,centBins,nPtBins,pTbins,fDCAzNbins,dcazBins);
      fTPCdEdx[iS] = new TH2F(Form("f%cTPCdEdx",kLetters[iS]),";#it{p} (GeV/c);TPC dE/dx (a.u.);Entries",
          196 * 2,0.2,10.,2400,0,2400);
      fTPCdEdxTpcCut[iS] = new TH2F(Form("f%cTPCdEdxTpcCut",kLetters[iS]),";#it{p} (GeV/c);TPC dE/dx (a.u.);Entries",
          196 * 2,0.2,10.,2400,0,2400);
      fTPCdEdxTofCut[iS] = new TH2F(Form("f%cTPCdEdxTofCut",kLetters[iS]),";#it{p} (GeV/c);TPC dE/dx (a.u.);Entries",
          196 * 2,0.2,10.,2400,0,2400);

      fList->Add(fTOFsignal[iS]);
      fList->Add(fTPCcounts[iS]);
      fList->Add(fDCAxyTPC[iS]);
      fList->Add(fDCAzTPC[iS]);
      fList->Add(fDCAxyTOF[iS]);
      fList->Add(fDCAzTOF[iS]);
      fList->Add(fTPCdEdx[iS]);
      fList->Add(fTPCdEdxTpcCut[iS]);
      fList->Add(fTPCdEdxTofCut[iS]);
    }

    for (int iS = 0; iS < 5; ++iS) {
      fTOFtemplates[iS] = new TH3F(Form("fTOFtemplates%i",iS),
          Form("%s;Centrality (%%);#it{p}_{T} (GeV/c);m_{TOF}^{2}-m_{PDG}^{2} (GeV/c^{2})^{2}",kNames[iS].Data()),
          nCentBins,centBins,nPtBins,pTbins,fTOFnBins,tofBins);
      fList->Add(fTOFtemplates[iS]);
    }
  }

  fEventCuts.AddQAplotsToList(fList);

  fTOFfunction = new TF1("fTOFfunction", TOFsignal, -2440., 2440., 4);

  PostData(1,fList);
}

/// This is the function that is evaluated for each event. The analysis code stays here.
///
/// \param options Deprecated parameter
/// \return void
///
void AliAnalysisTaskNucleiYieldESD::UserExec(Option_t *){
  /// The first check performed is on the particle type requested. If this type is AliPID::Unknown -
  /// the default one - the task does not process the information.
  if (fParticle == AliPID::kUnknown) {
    ::Error("AliAnalysisTaskNucleiYieldESD::UserExec", "No particle type set");
    PostData(1, fList);
    return;
  }

  /// If the setup is fine the current event is loaded in memory together with the analysis manager
  /// and the event handler.
  AliVEvent *ev = dynamic_cast<AliVEvent*> (InputEvent());
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* handl = (AliInputEventHandler*)mgr->GetInputEventHandler();

  /// Use the event cut class to apply the required selections
  if (!fEventCuts.AcceptEvent(ev)) {
    PostData(1, fList);
    return;
  }

  /// To perform the majority of the analysis - and also this one - the standard PID handler is
  /// required.
  fPID = handl->GetPIDResponse();
  
  /// Set TOF smearing parameters
  AliTOFPIDResponse& tofPid = fPID->GetTOFResponse();
  fTOFfunction->SetParameter(0, 1.);
  fTOFfunction->SetParameter(1, 0.);
  fTOFfunction->SetParameter(2, tofPid.GetTimeResolution());
  fTOFfunction->SetParameter(3, fTOFtail);

  float centrality = fEventCuts.GetCentrality();
  fCentralityClasses->Fill(centrality);

  AliMCEvent* mcEvent = nullptr;
  if (fIsMC) {
    // get branch "mcparticles"
    mcEvent = MCEvent();
    if (!mcEvent) {
      ::Fatal("AliAnalysisTaskNucleiYieldESD::UserExec", "Missing MC Event!");
      PostData(1, fList);
      return;
    }

    /// Making the list of the deuterons we want to measure
    for (int iMC = 0; iMC < mcEvent->GetNumberOfTracks(); ++iMC) {
      TParticle* part = mcEvent->Particle(iMC);
      const int pdg = part->GetPdgCode();
      int iC = pdg > 0 ? 1 : 0;
      if (TMath::Abs(part->Y()) > fYregion) continue;
      if (TMath::Abs(pdg) == fPDG) {
        fProduction[iC]->Fill(part->P());
        if (mcEvent->IsPhysicalPrimary(iMC)) fTotal[iC]->Fill(centrality,part->Pt());
      }
    }
  }

  /// Checking how many deuterons in acceptance are reconstructed well
  for (int iT = 0; iT < (int)ev->GetNumberOfTracks(); ++iT) {
    /// Create temporary track to avoid problems with other tasks on trains.
    AliESDtrack temp(*dynamic_cast<AliESDtrack*>(ev->GetTrack(iT)));
    AliESDtrack* track = &temp;

    if (track->GetID() <= 0) continue;
    int iC = track->Charge() > 0 ? 1 : 0;

    if (!fTrackCuts.AcceptTrack(track)) continue;

    const float energy = sqrt(track->P() * track->P() + fPDGMassOverZ * fPDGMassOverZ);
    const float pz = track->Pz();
    const float rapidity = (energy != TMath::Abs(pz)) ? 0.5* TMath::Log((energy+pz)/(energy-pz)) : -200.;
    if (fabs(rapidity) > fYregion) continue;
    const float beta = HasTOF(track);

    float dca[2], dcaCov[3];
    track->GetImpactParameters(dca,dcaCov);

    float pT = track->Pt();
    PtCorrection(pT,track->Charge() > 0);

    if (fIsMC) {
      const int id = TMath::Abs(track->GetLabel());
      TParticle* part = mcEvent->Particle(id);
      if (!part) continue;
      const int pdg = part->GetPdgCode();
      if (TMath::Abs(pdg) == fPDG) {
        if (mcEvent->IsPhysicalPrimary(id)) 
          fDCAPrimaryTPC[iC]->Fill(centrality,part->Pt(),dca[0]);
        else if (mcEvent->IsSecondaryFromMaterial(id))
          fDCASecondaryTPC[iC]->Fill(centrality,part->Pt(),dca[0]);
        fITS_TPC[iC]->Fill(centrality,part->Pt());
        fPtCorrection[iC]->Fill(pT,part->Pt() - pT);

        if (beta > 0) {
          fITS_TPC_TOF[iC]->Fill(centrality,part->Pt());
          if (mcEvent->IsPhysicalPrimary(id)) 
            fDCAPrimaryTOF[iC]->Fill(centrality,part->Pt(),dca[0]);
          else if (mcEvent->IsSecondaryFromMaterial(id))
            fDCASecondaryTOF[iC]->Fill(centrality,part->Pt(),dca[0]);
        }
      }
    } else {

      float tpc_n_sigma = fPID->GetTPCResponse().GetNumberOfSigmas(track, fParticle);
      fDCAxyTPC[iC]->Fill(centrality, pT, dca[0]);
      fDCAzTPC[iC]->Fill(centrality, pT, dca[1]);
      fTPCcounts[iC]->Fill(centrality, pT, tpc_n_sigma);
      if (beta < 1.e-10) continue;
      /// \f$ m = \frac{p}{\beta\gamma} \f$
      const float m = track->P() * track->P() * (1.f / (beta * beta) - 1.f);
      fTPCdEdx[iC]->Fill(track->GetTPCmomentum(),track->GetTPCsignal());
      if (fabs(fPID->NumberOfSigmasTOF(track,fParticle)) < 3)
        fTPCdEdxTofCut[iC]->Fill(track->GetTPCmomentum(),track->GetTPCsignal());
      if (fabs(tpc_n_sigma) > fTPCnSigmaCut) continue;
      fTPCdEdxTpcCut[iC]->Fill(track->GetTPCmomentum(),track->GetTPCsignal());
        
      fDCAxyTOF[iC]->Fill(centrality, pT, dca[0]);
      fDCAzTOF[iC]->Fill(centrality, pT, dca[1]);
      fTOFsignal[iC]->Fill(centrality, pT, m - fPDGMassOverZ * fPDGMassOverZ);
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

  } // End AOD track loop

  //  Post output data.
  PostData(1,fList);
}

/// Merge the output. Called once at the end of the query.
///
/// \return void
///
void AliAnalysisTaskNucleiYieldESD::Terminate(Option_t *) {
  return;
}

/// This function checks whether a track has or has not a prolongation in TOF.
///
/// \param track Track that has to be checked
/// \return \f$\beta\f$ of the particle, -1 means that there is no correct prolongation in TOF.
///
float AliAnalysisTaskNucleiYieldESD::HasTOF(AliVTrack *track) {
  bool hasTOFout  = track->GetStatus() & AliVTrack::kTOFout;
  bool hasTOFtime = track->GetStatus() & AliVTrack::kTIME;
  const float len = track->GetIntegratedLength();
  bool hasTOF = bool(hasTOFout & hasTOFtime) && (len > 350.f);

  if (!hasTOF) return -1.;
  const float p = track->P();
  const float tim = track->GetTOFsignal() - fPID->GetTOFResponse().GetStartTime(p);
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
void AliAnalysisTaskNucleiYieldESD::SetCentBins(int nbins, float *bins) {
  fCentBins.Set(nbins + 1, bins);
}

/// This functions sets the \f$\mathrm{DCA}_{xy}\f$ bins used in the analysis
///
/// \param nbins Number of \f$\mathrm{DCA}_{xy}\f$ bins
/// \param bins Array with nbins + 1 elements contanining the edges of the bins
/// \return void
///
void AliAnalysisTaskNucleiYieldESD::SetDCABins(int nbins, float min, float max) {
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
void AliAnalysisTaskNucleiYieldESD::SetDCABins(int nbins, float *bins) {
  fDCABins.Set(nbins + 1, bins);
}

/// This functions sets the \f$p_{\mathrm{T}}\f$ bins used in the analysis
///
/// \param nbins Number of \f$p_{\mathrm{T}}\f$ bins
/// \param bins Array with nbins + 1 elements contanining the edges of the bins
/// \return void
///
void AliAnalysisTaskNucleiYieldESD::SetPtBins(int nbins, float *bins) {
  fPtBins.Set(nbins + 1, bins);
}

/// This function sets the number of TOF bins and the boundaries of the histograms
///
/// \param nbins Number of bins
/// \param min Lower boundary of the histogram
/// \param max Higher boundary of the histogram
/// \return void
///
void AliAnalysisTaskNucleiYieldESD::SetTOFBins(int nbins, float min, float max) {
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
void AliAnalysisTaskNucleiYieldESD::SetDCAzBins(int nbins, float limit) {
  fDCAzNbins = nbins;
  fDCAzLimit = limit;
}

/// This function sets the particle type to be analysed
///
/// \param part Particle type
/// \return void
///
void AliAnalysisTaskNucleiYieldESD::SetParticleType(AliPID::EParticleType part) {
  fParticle = part;
  fPDGMass = AliPID::ParticleMass(part);
  fPDGMassOverZ = AliPID::ParticleMassZ(part);
  fPDG = AliPID::ParticleCode(part);
}

/// This function provides the correction for wrongly calculated \f$p_{\mathrm{T}}\f$.
///
/// \param pt \f$p_{\mathrm{T}}\f$ of the track
/// \param positiveCharge True if the track has positive sign.
/// \return void
///
void AliAnalysisTaskNucleiYieldESD::PtCorrection(float &pt, bool positiveCharge) {
  TArrayF& par = (positiveCharge) ? fPtCorrectionM : fPtCorrectionA;
  if (!par.GetSize()) return;
  const float correction = par[0] + par[1] * TMath::Exp(par[2] * pt);
  pt -= correction;
}

