#include "AliAnalysisTaskStrangenessRatios.h"

#include "AliAnalysisDataContainer.h"
#include "AliLog.h"

#include <algorithm>
#include <cmath>
#include <string>
using std::string;

// ROOT includes
#include <TAxis.h>
#include <TChain.h>
#include <TH2F.h>
#include <TList.h>
#include <TRandom3.h>
#include <TTree.h>

// ALIROOT includes
#include "AliAnalysisManager.h"
#include "AliAODTrack.h"
#include "AliAODcascade.h"
#include "AliAODMCParticle.h"
#include "AliMCEvent.h"
#include "AliInputEventHandler.h"
#include "AliPIDResponse.h"
#include "AliAODEvent.h"
#include "AliVEventHandler.h"
#include "AliVTrack.h"
#include "AliAODMCHeader.h"
#include "AliAnalysisUtils.h"

///\cond CLASSIMP
ClassImp(AliAnalysisTaskStrangenessRatios);
///\endcond

namespace
{

  double Sq(double x)
  {
    return x * x;
  }

  constexpr int kLambdaPdg{3122};
  constexpr double kLambdaMass{1.115683};
  constexpr int kXiPdg{3312};
  constexpr double kXiMass{1.32171};
  constexpr int kOmegaPdg{3334};
  constexpr double kOmegaMass{1.67245};
  constexpr double kcTauXi{4.91359839};
  constexpr double kcTauOmega{2.46129608};
  constexpr double kcTau[2]{kcTauXi, kcTauOmega};

}

/// Standard and default constructor of the class.
///
/// \param taskname Name of the task
/// \param partname Name of the analysed particle
///
AliAnalysisTaskStrangenessRatios::AliAnalysisTaskStrangenessRatios(bool isMC, TString taskname) : AliAnalysisTaskSE(taskname.Data()),
                                                                                                  fEventCut{false},
                                                                                                  fMC{isMC}
{
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
  DefineOutput(3, TTree::Class());
}

/// Standard destructor
///
AliAnalysisTaskStrangenessRatios::~AliAnalysisTaskStrangenessRatios()
{
  if (AliAnalysisManager::GetAnalysisManager()->IsProofMode())
    return;
  if (fList)
    delete fList;
  if (fTree)
    delete fTree;
  if (fTreeLambda)
    delete fTreeLambda;
  if (!fMC)
  {
    delete fRecCascade;
    delete fRecLambda;
  }
}

/// This function creates all the histograms and all the objects in general used during the analysis
/// \return void
///
void AliAnalysisTaskStrangenessRatios::UserCreateOutputObjects()
{

  fList = new TList();
  fList->SetOwner(kTRUE);
  fEventCut.AddQAplotsToList(fList);
  fRecCascade = fMC ? &fGenCascade : new MiniCascade;
  fRecLambda = fMC ? &fGenLambda : new MiniLambda;

  OpenFile(2);
  fTree = new TTree("XiOmegaTree", "Xi and Omega Tree");
  fTreeLambda = new TTree("LambdaTree", "Lambda");
  if (fMC)
  {
    fTree->Branch("MiniCascadeMC", &fGenCascade);
    fTreeLambda->Branch("MiniLambdaMC", &fGenLambda);
    fMCEvent = MCEvent();
  }
  else
  {
    fTree->Branch("MiniCascade", fRecCascade);
    fTreeLambda->Branch("MiniLambda", fRecLambda);
  }

  PostAllData();
}

/// This is the function that is evaluated for each event. The analysis code stays here.
///
/// \param options Deprecated parameter
/// \return void
///
void AliAnalysisTaskStrangenessRatios::UserExec(Option_t *)
{
  AliAODEvent *ev = (AliAODEvent *)InputEvent();
  if (!fEventCut.AcceptEvent(ev))
  {
    PostAllData();
    return;
  }

  double bField{ev->GetMagneticField()};
  auto pvObj = fEventCut.GetPrimaryVertex();
  double pv[3];
  pvObj->GetXYZ(pv);

  /// To perform the majority of the analysis - and also this one - the standard PID handler is
  /// required.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *handl = (AliInputEventHandler *)mgr->GetInputEventHandler();
  fPID = handl->GetPIDResponse();

  fRecCascade->centrality = fEventCut.GetCentrality();
  fRecLambda->centrality = fEventCut.GetCentrality();

  double rdmState{gRandom->Uniform()};

  std::vector<int> checkedLabel, checkedLambdaLabel;
  fGenCascade.isReconstructed = true;
  for (int iCasc = 0; iCasc < ev->GetNumberOfCascades(); iCasc++)
  {
    AliAODcascade *casc = ev->GetCascade(iCasc);
    if (!casc)
      continue;

    //cascade and V0 2D radii
    double vtxCasc[3]{casc->DecayVertexXiX(), casc->DecayVertexXiY(), casc->DecayVertexXiZ()};
    fRecCascade->radius = std::hypot(vtxCasc[0], vtxCasc[1]);
    fRecCascade->radiusV0 = casc->RadiusSecVtx();

    //get daughter tracks (positive, negative and bachelor)
    AliAODTrack *pTrackCasc = dynamic_cast<AliAODTrack *>(casc->GetDaughter(0));
    AliAODTrack *nTrackCasc = dynamic_cast<AliAODTrack *>(casc->GetDaughter(1));
    AliAODTrack *bTrackCasc = dynamic_cast<AliAODTrack *>(casc->GetDecayVertexXi()->GetDaughter(0));
    if (!pTrackCasc || !nTrackCasc || !bTrackCasc)
    {
      AliWarning("ERROR: Could not retrieve one of the 3 AOD daughter tracks of the cascade ...\n");
      continue;
    }

    if (!(pTrackCasc->GetStatus() & AliVTrack::kTPCrefit) || !(nTrackCasc->GetStatus() & AliVTrack::kTPCrefit) || !(bTrackCasc->GetStatus() & AliVTrack::kTPCrefit))
    {
      continue;
    }

    if (pTrackCasc->GetTPCsignalN() < 50 || nTrackCasc->GetTPCsignalN() < 50 || bTrackCasc->GetTPCsignalN() < 50)
    {
      continue;
    }

    if (std::abs(pTrackCasc->Eta()) > 0.8 || std::abs(nTrackCasc->Eta()) > 0.8 || std::abs(bTrackCasc->Eta()) > 0.8)
    {
      continue;
    }

    int labMothBac = -9999;
    int pdgCascade = -1;
    bool physPrim = true;
    if (fMC)
    {
      fGenCascade.pdg = 0;
      auto posPart = (AliAODMCParticle *)fMCEvent->GetTrack(std::abs(pTrackCasc->GetLabel()));
      auto negPart = (AliAODMCParticle *)fMCEvent->GetTrack(std::abs(nTrackCasc->GetLabel()));
      auto bacPart = (AliAODMCParticle *)fMCEvent->GetTrack(std::abs(bTrackCasc->GetLabel()));
      // Check lambda
      int labMothPos = posPart->GetMother();
      int labMothNeg = negPart->GetMother();
      labMothBac = bacPart->GetMother();
      auto lambda = (AliAODMCParticle *)fMCEvent->GetTrack(labMothNeg);
      auto cascade = (AliAODMCParticle *)fMCEvent->GetTrack(labMothBac);
      if (!cascade)
      {
        continue;
      }
      pdgCascade = std::abs(cascade->GetPdgCode());
      if (lambda && labMothNeg == labMothPos && std::abs(lambda->GetPdgCode()) == kLambdaPdg)
      {
        int labMothLam = lambda->GetMother();
        if (labMothLam == labMothBac && (pdgCascade == kXiPdg || pdgCascade == kOmegaPdg))
        {
          fGenCascade.pdg = cascade->GetPdgCode();
          fGenCascade.ptMC = cascade->Pt();
          fGenCascade.etaMC = cascade->Eta();
          fGenCascade.yMC = cascade->Y();
          physPrim = cascade->IsPhysicalPrimary();
          double pv[3], sv[3];
          cascade->XvYvZv(pv);
          bacPart->XvYvZv(sv);
          fGenCascade.ctMC = std::sqrt(Sq(pv[0] - sv[0]) + Sq(pv[1] - sv[1]) + Sq(pv[2] - sv[2])) * cascade->M() / cascade->P();
          checkedLabel.push_back(labMothBac);
        }
      }
      if (fOnlyTrueCandidates && fGenCascade.pdg == 0)
        continue;
    }

    fRecCascade->matter = casc->AlphaV0() > 0;
    fRecCascade->tpcNsigmaV0Pi = fPID->NumberOfSigmasTPC(fRecCascade->matter ? nTrackCasc : pTrackCasc, AliPID::kPion);
    fRecCascade->tpcNsigmaV0Pr = fPID->NumberOfSigmasTPC(fRecCascade->matter ? pTrackCasc : nTrackCasc, AliPID::kProton);

    fRecCascade->tpcClBach = bTrackCasc->GetTPCsignalN();
    fRecCascade->tpcClV0Pi = (fRecCascade->matter ? nTrackCasc : pTrackCasc)->GetTPCsignalN();
    fRecCascade->tpcClV0Pr = (fRecCascade->matter ? pTrackCasc : nTrackCasc)->GetTPCsignalN();

    //DCA info
    fRecCascade->dcaBachV0 = casc->DcaXiDaughters();
    fRecCascade->dcaBachPV = casc->DcaBachToPrimVertex();
    fRecCascade->dcaV0prPV = fRecCascade->matter ? casc->DcaPosToPrimVertex() : casc->DcaNegToPrimVertex();
    fRecCascade->dcaV0piPV = fRecCascade->matter ? casc->DcaNegToPrimVertex() : casc->DcaPosToPrimVertex();
    fRecCascade->dcaV0tracks = casc->DcaV0Daughters();
    fRecCascade->dcaV0PV = casc->DcaV0ToPrimVertex();

    //cascade and V0 cosine of pointing angle
    fRecCascade->cosPA = casc->CosPointingAngleXi(pv[0], pv[1], pv[2]);
    fRecCascade->cosPAV0 = casc->CosPointingAngle(pv);

    //TOF matching
    fRecCascade->hasTOFhit = !pTrackCasc->GetTOFBunchCrossing(bField) || !nTrackCasc->GetTOFBunchCrossing(bField) || !bTrackCasc->GetTOFBunchCrossing(bField);

    //track status: ( fCasc_NegTrackStatus & AliESDtrack::kITSrefit ) is the codition to check kITSrefit
    fRecCascade->hasITSrefit = (nTrackCasc->GetStatus() & AliVTrack::kITSrefit) || (pTrackCasc->GetStatus() & AliVTrack::kITSrefit) || (bTrackCasc->GetStatus() & AliVTrack::kITSrefit);

    fRecCascade->V0invMassDelta = ((fRecCascade->matter) ? casc->MassLambda() : casc->MassAntiLambda()) - kLambdaMass;

    //transverse momentum and eta
    fRecCascade->pt = std::sqrt(casc->Pt2Xi());

    //////////////////////////////
    //crossed raws
    double lCrosRawsPos = pTrackCasc->GetTPCClusterInfo(2, 1);
    double lCrosRawsNeg = nTrackCasc->GetTPCClusterInfo(2, 1);
    double lCrosRawsBac = bTrackCasc->GetTPCClusterInfo(2, 1);
    fCascLeastCRaws = (int)(lCrosRawsPos < lCrosRawsNeg ? std::min(lCrosRawsPos, lCrosRawsBac) : std::min(lCrosRawsNeg, lCrosRawsBac));
    //crossed raws / Findable clusters
    double lCrosRawsOvFPos = lCrosRawsPos / ((double)(pTrackCasc->GetTPCNclsF()));
    double lCrosRawsOvFNeg = lCrosRawsNeg / ((double)(nTrackCasc->GetTPCNclsF()));
    double lCrosRawsOvFBac = lCrosRawsBac / ((double)(bTrackCasc->GetTPCNclsF()));
    fCascLeastCRawsOvF = lCrosRawsOvFPos < lCrosRawsOvFNeg ? std::min(lCrosRawsOvFPos, lCrosRawsOvFBac) : std::min(lCrosRawsOvFNeg, lCrosRawsOvFBac);
    ///////////////////////////////

    //calculate DCA Bachelor-Baryon to remove "bump" structure in InvMass
    fRecCascade->bachBarCosPA = casc->BachBaryonCosPA();
    double p = std::sqrt(casc->Ptot2Xi());
    fRecCascade->eta = 0.5 * std::log((p + casc->MomXiZ()) / (p - casc->MomXiZ() + 1.e-16));

    //distance over total momentum
    double lOverP = std::sqrt((Sq(vtxCasc[0] - pv[0]) + Sq(vtxCasc[1] - pv[1]) + Sq(vtxCasc[2] - pv[2])) / (casc->Ptot2Xi() + 1e-10));
    double ctLambda = std::sqrt(Sq(vtxCasc[0] - casc->GetSecVtxX()) + Sq(vtxCasc[1] - casc->GetSecVtxY()) + Sq(vtxCasc[2] - casc->GetSecVtxZ())) / (casc->P() + 1e-10);

    if (std::abs(casc->MassOmega() - kOmegaMass) * 1000 < 30)
    {
      fRecCascade->mass = casc->MassOmega();
      fRecCascade->ct = lOverP * kOmegaMass;
      fRecCascade->tpcNsigmaBach = fPID->NumberOfSigmasTPC(bTrackCasc, AliPID::kKaon);
      fRecCascade->competingMass = std::abs(casc->MassXi() - kXiMass);
      if (physPrim && IsTopolSelected(true))
      {
        fRecCascade->isOmega = true;
        fTree->Fill();
      }
      else if (fMC && std::find(checkedLabel.begin(), checkedLabel.end(), labMothBac) != checkedLabel.end() && (pdgCascade == kOmegaPdg))
      {
        checkedLabel.erase(std::find(checkedLabel.begin(), checkedLabel.end(), labMothBac)); //checked particles that didn't pass the topological cut (have to be filled later)
      }
    }
    if (std::abs(casc->MassXi() - kXiMass) * 1000 < 30)
    {
      fRecCascade->mass = casc->MassXi();
      fRecCascade->ct = lOverP * kXiMass;
      fRecCascade->tpcNsigmaBach = fPID->NumberOfSigmasTPC(bTrackCasc, AliPID::kPion);
      fRecCascade->competingMass = std::abs(casc->MassOmega() - kOmegaMass);
      if (physPrim && IsTopolSelected(false))
      {
        fRecCascade->isOmega = false;
        fTree->Fill();
      }
      else if (fMC && std::find(checkedLabel.begin(), checkedLabel.end(), labMothBac) != checkedLabel.end() && (pdgCascade == kXiPdg))
      {
        checkedLabel.erase(std::find(checkedLabel.begin(), checkedLabel.end(), labMothBac)); //checked particles that didn't pass the topological cut (have to be filled later)
      }
    }
  }

  if (fFillLambdas && rdmState < fLambdaDownscaling)
  {
    fGenLambda.isReconstructed = true;
    for (int iV0{0}; iV0 < ev->GetNumberOfV0s(); ++iV0)
    {
      AliAODv0 *v0{ev->GetV0(iV0)};
      if (!v0)
        continue;

      fRecLambda->radius = v0->RadiusSecVtx();

      //get daughter tracks (positive, negative and bachelor)
      AliAODTrack *pTrack = dynamic_cast<AliAODTrack *>(v0->GetDaughter(0));
      AliAODTrack *nTrack = dynamic_cast<AliAODTrack *>(v0->GetDaughter(1));
      if (!pTrack || !nTrack)
      {
        AliWarning("ERROR: Could not retrieve one of the 2 AOD daughter tracks of the lambdas ...\n");
        continue;
      }

      if (!(pTrack->GetStatus() & AliVTrack::kTPCrefit) || !(nTrack->GetStatus() & AliVTrack::kTPCrefit) ||
          pTrack->GetTPCsignalN() < 50 || nTrack->GetTPCsignalN() < 50 ||
          std::abs(pTrack->Eta()) > 0.8 || std::abs(nTrack->Eta()) > 0.8)
      {
        continue;
      }

      int lambdaLabel{-1};
      if (fMC)
      {
        fGenLambda.pdg = 0;
        auto posPart = (AliAODMCParticle *)fMCEvent->GetTrack(std::abs(pTrack->GetLabel()));
        auto negPart = (AliAODMCParticle *)fMCEvent->GetTrack(std::abs(nTrack->GetLabel()));
        // Check lambda
        int labMothPos = posPart->GetMother();
        int labMothNeg = negPart->GetMother();
        auto lambda = (AliAODMCParticle *)fMCEvent->GetTrack(labMothNeg);
        if (lambda && labMothNeg == labMothPos && std::abs(lambda->GetPdgCode()) == kLambdaPdg)
        {
          lambdaLabel = labMothNeg;
          fGenLambda.pdg = lambda->GetPdgCode();
          fGenLambda.ptMC = lambda->Pt();
          fGenLambda.etaMC = lambda->Eta();
          fGenLambda.yMC = lambda->Y();
          fGenLambda.isPrimary = lambda->IsPhysicalPrimary();
          double pv[3], sv[3];
          lambda->XvYvZv(pv);
          posPart->XvYvZv(sv);
          fGenCascade.ctMC = std::sqrt(Sq(pv[0] - sv[0]) + Sq(pv[1] - sv[1]) + Sq(pv[2] - sv[2])) * lambda->M() / lambda->P();
        }
        if (fOnlyTrueLambdas && fGenLambda.pdg == 0)
          continue;
      }

      fRecLambda->matter = v0->AlphaV0() > 0;
      auto proton = fRecLambda->matter ? pTrack : nTrack;
      auto pion = fRecLambda->matter ? nTrack : pTrack;

      fRecLambda->pt = v0->Pt();
      fRecLambda->eta = v0->Eta();
      fRecLambda->mass = v0->MassLambda();
      fRecLambda->radius = v0->RadiusV0();
      fRecLambda->dcaV0PV = v0->DcaV0ToPrimVertex();
      fRecLambda->dcaPrPV = fRecLambda->matter ? v0->DcaPosToPrimVertex() : v0->DcaNegToPrimVertex();
      fRecLambda->dcaPiPV = fRecLambda->matter ? v0->DcaNegToPrimVertex() : v0->DcaPosToPrimVertex();
      fRecLambda->dcaV0tracks = v0->DcaV0Daughters();
      fRecLambda->cosPA = v0->CosPointingAngle(pv);
      fRecLambda->tpcNsigmaPi = fPID->NumberOfSigmasTPC(pion, AliPID::kPion);
      fRecLambda->tpcNsigmaPr = fPID->NumberOfSigmasTPC(proton, AliPID::kProton);
      fRecLambda->tpcClV0Pi = pion->GetTPCsignalN();
      fRecLambda->tpcClV0Pr = proton->GetTPCsignalN();
      fRecLambda->hasTOFhit = !pTrack->GetTOFBunchCrossing(bField) || !nTrack->GetTOFBunchCrossing(bField);
      fRecLambda->hasITSrefit = nTrack->GetStatus() & AliVTrack::kITSrefit || pTrack->GetStatus() & AliVTrack::kITSrefit;

      //crossed raws
      double lCrosRawsPos = pTrack->GetTPCClusterInfo(2, 1);
      double lCrosRawsNeg = nTrack->GetTPCClusterInfo(2, 1);
      fLambdaLeastCRaws = std::min(lCrosRawsPos, lCrosRawsNeg);
      //crossed raws / Findable clusters
      double lCrosRawsOvFPos = lCrosRawsPos / ((double)(pTrack->GetTPCNclsF()));
      double lCrosRawsOvFNeg = lCrosRawsNeg / ((double)(nTrack->GetTPCNclsF()));
      fLambdaLeastCRawsOvF = std::min(lCrosRawsOvFPos, lCrosRawsOvFNeg);

      if (IsTopolSelectedLambda())
      {
        if (lambdaLabel != -1)
        {
          checkedLambdaLabel.push_back(lambdaLabel);
        }
        fTreeLambda->Fill();
      }
    }
  }

  if (fMC)
  {
    fGenCascade.isReconstructed = false;
    //OOB pileup
    AliAODMCHeader *header = static_cast<AliAODMCHeader *>(ev->FindListObject(AliAODMCHeader::StdBranchName()));
    if (!header)
    {
      AliWarning("No header found.");
      PostAllData();
      return;
    }
    TClonesArray *MCTrackArray = dynamic_cast<TClonesArray *>(ev->FindListObject(AliAODMCParticle::StdBranchName()));
    if (MCTrackArray == NULL)
    {
      AliWarning("No MC track array found.");
      PostAllData();
      return;
    }
    //loop on generated
    for (int iT{0}; iT < fMCEvent->GetNumberOfTracks(); ++iT)
    {
      auto track = (AliAODMCParticle *)fMCEvent->GetTrack(iT);
      int pdg = std::abs(track->GetPdgCode());
      if (pdg != kXiPdg && pdg != kOmegaPdg)
      {
        continue;
      }
      if (std::find(checkedLabel.begin(), checkedLabel.end(), iT) != checkedLabel.end())
      {
        continue;
      }
      if (std::abs(track->Y()) > fCutY || !track->IsPhysicalPrimary() || AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(iT, header, MCTrackArray))
      { //removal of OOB pileup, cut on Y and PhysPrim
        continue;
      }
      fGenCascade.ptMC = track->Pt();
      fGenCascade.etaMC = track->Eta();
      fGenCascade.yMC = track->Y();
      fGenCascade.pdg = track->GetPdgCode();
      double pv[3], sv[3];
      track->XvYvZv(pv);
      for (int iD = track->GetDaughterFirst(); iD <= track->GetDaughterLast(); iD++)
      {
        auto daugh = (AliAODMCParticle *)fMCEvent->GetTrack(iD);
        if (!daugh)
        {
          continue;
        }
        if (std::abs(daugh->GetPdgCode()) == kLambdaPdg)
        {
          daugh->XvYvZv(sv);
          break;
        }
      }
      fGenCascade.ctMC = std::sqrt(Sq(pv[0] - sv[0]) + Sq(pv[1] - sv[1]) + Sq(pv[2] - sv[2])) * track->M() / track->P();
      fTree->Fill();
    }

    if (fFillLambdas && rdmState < fLambdaDownscaling)
    {
      fGenLambda.isReconstructed = false;
      //loop on generated
      for (int iT{0}; iT < fMCEvent->GetNumberOfTracks(); ++iT)
      {
        auto track = (AliAODMCParticle *)fMCEvent->GetTrack(iT);
        int pdg = std::abs(track->GetPdgCode());
        if (pdg != kLambdaPdg)
        {
          continue;
        }
        if (std::find(checkedLambdaLabel.begin(), checkedLambdaLabel.end(), iT) != checkedLambdaLabel.end())
        {
          continue;
        }

        if (std::abs(track->Y()) > fCutY || AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(iT, header, MCTrackArray))
        { //removal of OOB pileup, cut on Y and PhysPrim
          continue;
        }
        fGenLambda.ptMC = track->Pt();
        fGenLambda.etaMC = track->Eta();
        fGenLambda.yMC = track->Y();
        fGenLambda.pdg = track->GetPdgCode();
        fGenLambda.isPrimary = track->IsPhysicalPrimary();
        double pv[3], sv[3];
        track->XvYvZv(pv);
        for (int iD = track->GetDaughterFirst(); iD <= track->GetDaughterLast(); iD++)
        {
          auto daugh = (AliAODMCParticle *)fMCEvent->GetTrack(iD);
          if (!daugh)
          {
            continue;
          }
          if (std::abs(daugh->GetPdgCode()) == AliPID::ParticleCode(AliPID::kProton))
          {
            daugh->XvYvZv(sv);
            break;
          }
        }
        fGenLambda.ctMC = std::sqrt(Sq(pv[0] - sv[0]) + Sq(pv[1] - sv[1]) + Sq(pv[2] - sv[2])) * track->M() / track->P();
        fTreeLambda->Fill();
      }
    }
  }

  PostAllData();
}

AliAnalysisTaskStrangenessRatios *AliAnalysisTaskStrangenessRatios::AddTask(bool isMC, TString tskname, TString suffix)
{
  // Get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskStrangenessRatios", "No analysis manager found.");
    return nullptr;
  }

  // Check the analysis type using the event handlers connected to the analysis
  // manager.
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskStrangenessRatios", "This task requires an input event handler");
    return nullptr;
  }

  tskname.Append(suffix.Data());
  AliAnalysisTaskStrangenessRatios *task = new AliAnalysisTaskStrangenessRatios(isMC, tskname.Data());

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(
      Form("%s_summary", tskname.Data()), TList::Class(),
      AliAnalysisManager::kOutputContainer, "AnalysisResults.root");

  AliAnalysisDataContainer *coutput2 =
      mgr->CreateContainer(Form("%s_treeCascades", tskname.Data()), TTree::Class(),
                           AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
  coutput2->SetSpecialOutput();

  AliAnalysisDataContainer *coutput3 =
      mgr->CreateContainer(Form("%s_treeLambda", tskname.Data()), TTree::Class(),
                           AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
  coutput3->SetSpecialOutput();

  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, coutput1);
  mgr->ConnectOutput(task, 2, coutput2);
  mgr->ConnectOutput(task, 3, coutput3);
  return task;
}

//
//____________________________________________________________________________________________
bool AliAnalysisTaskStrangenessRatios::IsTopolSelected(bool isOmega)
{
  return fRecCascade->radius > fCutRadius[isOmega] &&
         fRecCascade->radiusV0 > fCutRadiusV0 &&
         fRecCascade->dcaBachPV > fCutDCABachToPV &&
         fRecCascade->dcaV0PV > fCutDCAV0toPV &&
         fRecCascade->dcaV0piPV > fCutDCAV0piToPV &&
         fRecCascade->dcaV0prPV > fCutDCAV0prToPV &&
         fRecCascade->dcaV0tracks < fCutDCAV0tracks &&
         fRecCascade->dcaBachV0 < fCutDCABachToV0[isOmega] &&
         fRecCascade->cosPA > fCutCosPA &&
         fRecCascade->cosPAV0 > fCutCosPAV0 &&
         std::abs(Eta2y(fRecCascade->pt, isOmega ? kOmegaMass : kXiMass, fRecCascade->eta)) < fCutY &&
         std::abs(fRecCascade->tpcNsigmaBach) < fCutNsigmaTPC &&
         std::abs(fRecCascade->tpcNsigmaV0Pr) < fCutNsigmaTPC &&
         std::abs(fRecCascade->tpcNsigmaV0Pi) < fCutNsigmaTPC &&
         fRecCascade->ct < fCutCt * kcTau[isOmega] &&
         fRecCascade->competingMass > fCutCompetingMass[isOmega] &&
         fRecCascade->tpcClBach > fCutTPCclu &&
         fRecCascade->tpcClV0Pi > fCutTPCclu &&
         fRecCascade->tpcClV0Pr > fCutTPCclu &&
         fCascLeastCRaws > fCutTPCrows &&
         fCascLeastCRawsOvF > fCutRowsOvF;
}

bool AliAnalysisTaskStrangenessRatios::IsTopolSelectedLambda()
{
  return fRecLambda->radius > fCutRadius[2] &&
         fRecLambda->cosPA > fCosPALambda &&
         fRecLambda->dcaPrPV > fCutDCALambdaPrToPV &&
         fRecLambda->dcaPiPV > fCutDCALambdaPiToPV &&
         fRecLambda->dcaV0tracks < fCutDCAV0tracks &&
         std::abs(Eta2y(fRecLambda->pt, kLambdaMass, fRecCascade->eta)) < fCutY &&
         fRecLambda->mass > fCutLambdaMass[0] && fRecLambda->mass < fCutLambdaMass[1] &&
         std::abs(fRecLambda->tpcNsigmaPi) < fCutNsigmaTPC &&
         std::abs(fRecLambda->tpcNsigmaPr) < fCutNsigmaTPC &&
         fLambdaLeastCRaws > fCutTPCrows &&
         fLambdaLeastCRawsOvF > fCutRowsOvF;
}

//
//____________________________________________________________________________________________
float AliAnalysisTaskStrangenessRatios::Eta2y(float pt, float m, float eta) const
{
  return std::asinh(pt / std::hypot(m, pt) * std::sinh(eta));
}

void AliAnalysisTaskStrangenessRatios::PostAllData()
{
  PostData(1, fList);
  PostData(2, fTree);
  PostData(3, fTreeLambda);
}

Bool_t AliAnalysisTaskStrangenessRatios::UserNotify() {
  TString cfn{CurrentFileName()};
  AliInfo(Form("Setting hash for file %s", cfn.Data()));

  gRandom->SetSeed(cfn.Hash());
  return true;
}
