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
#include "AliMCEvent.h"
#include "AliInputEventHandler.h"
#include "AliPIDResponse.h"
#include "AliAODEvent.h"
#include "AliVEventHandler.h"
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

  constexpr int kK0sPdg{310};
  constexpr double kK0sMass{0.497611};
  constexpr int kLambdaPdg{3122};
  constexpr double kLambdaMass{1.115683};
  constexpr int kXiPdg{3312};
  constexpr double kXiMass{1.32171};
  constexpr int kOmegaPdg{3334};
  constexpr double kOmegaMass{1.67245};
  constexpr double kcTauXi{4.91359839};
  constexpr double kcTauOmega{2.46129608};
  constexpr double kcTau[2]{kcTauXi, kcTauOmega};

  void getITScls(AliAODTrack* t, int &SPDcls, int &SDDSSDcls)
  {
    SPDcls = 0u;
    SDDSSDcls = 0u;
    for (int i = 0; i < 6; ++i)
    {
      if (t->HasPointOnITSLayer(i)){
        if (i < 2)
          SPDcls++;
        else
          SDDSSDcls++;
      }
    }
  }


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
  DefineOutput(4, TTree::Class());
  DefineOutput(5, TTree::Class());
  DefineOutput(6, TTree::Class());
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
  if (fTreeK0s)
    delete fTreeK0s;
  if (fTreeTrain)
    delete fTreeTrain;
  if (fTreeLambda)
    delete fTreeLambda;
  if (!fMC)
  {
    delete fRecCascade;
    delete fRecCascadeTrain;
    delete fRecLambda;
    delete fRecK0s;
    delete fRecLambdaBDTOut;
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
  fRecCascadeTrain = fMC ? &fGenCascadeTrain : new MiniCascadeTrain;
  fRecLambda = fMC ? &fGenLambda : new MiniLambda;
  fRecK0s = fMC ? &fGenK0s : new MiniK0s;
  fRecLambdaBDTOut = new MiniLambdaBDTOut;

  OpenFile(2);
  fTree = new TTree("XiOmegaTree", "Xi and Omega Tree");
  fTreeTrain = new TTree("XiOmegaTreeTrain", "Xi and Omega Tree (train)");
  fTreeLambda = new TTree("LambdaTree", "Lambda");
  fTreeLambdaBDTOut = new TTree("LambdaTreeBDTOut", "LambdaBDT");
  fTreeK0s = new TTree("K0sTree", "K0s");
  if (fMC)
  {
    fTree->Branch("MiniCascadeMC", &fGenCascade);
    fTreeTrain->Branch("MiniCascadeTrainMC", &fGenCascadeTrain);
    fTreeLambda->Branch("MiniLambdaMC", &fGenLambda);
    fTreeK0s->Branch("MiniK0sMC", &fGenK0s);
    fMCEvent = MCEvent();
  }
  else
  {
    fTree->Branch("MiniCascade", fRecCascade);
    fTreeTrain->Branch("MiniCascadeTrain", fRecCascadeTrain);
    fTreeLambda->Branch("MiniLambda", fRecLambda);
    fTreeK0s->Branch("MiniK0s", fRecK0s);
    fTreeLambdaBDTOut->Branch("MiniLambdaBDTOut", fRecLambdaBDTOut);
  }

  if (fFillLambdasBDTOut)
  {
    for (int iB=0; iB<(fCtBinsBDT.GetSize()-1); ++iB)
    {
      fBDT.push_back(new AliExternalBDT());
      if (!fBDT[iB]->LoadXGBoostModel(Form("%s%.0f_%.0f.model",fBDTPath.data(),fCtBinsBDT[iB],fCtBinsBDT[iB+1])))
      {
        fBDT[iB] = nullptr;
      }
    }
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

  double rdmState{gRandom->Uniform()};

  std::vector<int> checkedLabel, checkedLambdaLabel, checkedK0sLabel;

  if (fFillCascades || fFillCascadesTrain)
  {
    fRecCascade->centrality = fEventCut.GetCentrality();

    fGenCascade.isReconstructed = true;
    fGenCascadeTrain.isReconstructed = true;

    for (int iCasc = 0; iCasc < ev->GetNumberOfCascades(); iCasc++)
    {
      AliAODcascade *casc = ev->GetCascade(iCasc);
      if (!casc)
        continue;
      if (casc->GetOnFlyStatus() != fUseOnTheFly)
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
        fGenCascadeTrain.pdg = 0;
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
            fGenCascadeTrain.pdg = cascade->GetPdgCode();
            fGenCascade.ptMC = cascade->Pt();
            fGenCascade.etaMC = cascade->Eta();
            fGenCascade.yMC = cascade->Y();
            physPrim = cascade->IsPhysicalPrimary();
            fGenCascade.flag = 0u;
            if (physPrim)
              fGenCascade.flag |= kPrimary;
            else
              fGenCascade.flag |= cascade->IsSecondaryFromWeakDecay() ? kSecondaryFromWD : kSecondaryFromMaterial;
            double pv[3], sv[3];
            cascade->XvYvZv(pv);
            bacPart->XvYvZv(sv);
            fGenCascade.ctMC = std::sqrt(Sq(pv[0] - sv[0]) + Sq(pv[1] - sv[1]) + Sq(pv[2] - sv[2])) * cascade->M() / cascade->P();
            checkedLabel.push_back(labMothBac);

            //train cascade
            fGenCascadeTrain.ptMC = fGenCascade.ptMC;
            fGenCascadeTrain.etaMC = fGenCascade.etaMC;
            fGenCascadeTrain.ctMC = fGenCascade.ctMC;
            fGenCascadeTrain.yMC = fGenCascade.yMC;
            fGenCascadeTrain.flag = fGenCascade.flag;
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
      // double ctLambda = std::sqrt(Sq(vtxCasc[0] - casc->GetSecVtxX()) + Sq(vtxCasc[1] - casc->GetSecVtxY()) + Sq(vtxCasc[2] - casc->GetSecVtxZ())) / (casc->P() + 1e-10);

      //training cascade
      fRecCascadeTrain->dcaBachV0 = fRecCascade->dcaBachV0;
      fRecCascadeTrain->dcaBachPV = fRecCascade->dcaBachPV;
      fRecCascadeTrain->dcaV0prPV = fRecCascade->dcaV0prPV;
      fRecCascadeTrain->dcaV0piPV = fRecCascade->dcaV0piPV;
      fRecCascadeTrain->dcaV0tracks = fRecCascade->dcaV0tracks;
      fRecCascadeTrain->dcaV0PV = fRecCascade->dcaV0PV;
      fRecCascadeTrain->cosPA = fRecCascade->cosPA;
      fRecCascadeTrain->cosPAV0 = fRecCascade->cosPAV0;
      fRecCascadeTrain->tpcNsigmaV0Pr = fRecCascade->tpcNsigmaV0Pr;
      fRecCascadeTrain->centrality = fRecCascade->centrality;
      fRecCascadeTrain->pt = fRecCascade->matter ? fRecCascade->pt : -fRecCascade->pt;;

      if (std::abs(casc->MassOmega() - kOmegaMass) * 1000 < 30)
      {
        fRecCascade->mass = casc->MassOmega();
        fRecCascade->ct = lOverP * kOmegaMass;
        fRecCascade->tpcNsigmaBach = fPID->NumberOfSigmasTPC(bTrackCasc, AliPID::kKaon);
        fRecCascade->competingMass = std::abs(casc->MassXi() - kXiMass);

        fRecCascadeTrain->mass = fRecCascade->mass;
        if (physPrim && IsTopolSelected(true))
        {
          fRecCascade->isOmega = true;
          fRecCascadeTrain->mass = -fRecCascadeTrain->mass;
          if (fFillCascadesTrain)
            fTreeTrain->Fill();
          if (fFillCascades)
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

        fRecCascadeTrain->mass = fRecCascade->mass;
        if (physPrim && IsTopolSelected(false))
        {
          fRecCascade->isOmega = false;
          fRecCascadeTrain->mass = -fRecCascadeTrain->mass;
          if (fFillCascadesTrain)
            fTreeTrain->Fill();
          if (fFillCascades)
            fTree->Fill();
        }
        else if (fMC && std::find(checkedLabel.begin(), checkedLabel.end(), labMothBac) != checkedLabel.end() && (pdgCascade == kXiPdg))
        {
          checkedLabel.erase(std::find(checkedLabel.begin(), checkedLabel.end(), labMothBac)); //checked particles that didn't pass the topological cut (have to be filled later)
        }
      }
    }
  }

  if ((fFillLambdas||fFillLambdasBDTOut) && rdmState < fLambdaDownscaling)
  {
    fGenLambda.isReconstructed = true;
    for (int iV0{0}; iV0 < ev->GetNumberOfV0s(); ++iV0)
    {
      AliAODv0 *v0{ev->GetV0(iV0)};
      if (!v0)
        continue;
      if (v0->GetOnFlyStatus() != fUseOnTheFly)
        continue;

      fRecLambda->centrality = fEventCut.GetCentrality();
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
          std::abs(pTrack->Eta()) > 0.8 || std::abs(nTrack->Eta()) > 0.8 ||
          pTrack->Chi2perNDF() > 4 || nTrack->Chi2perNDF() > 4)
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
          fGenLambda.magFieldPolarity = bField > 0.;
          fGenLambda.pdg = lambda->GetPdgCode();
          fGenLambda.ptMC = lambda->Pt();
          fGenLambda.etaMC = lambda->Eta();
          fGenLambda.yMC = lambda->Y();
          fGenLambda.isPrimary = lambda->IsPhysicalPrimary();
          double ov[3], dv[3];
          lambda->XvYvZv(ov);
          posPart->XvYvZv(dv);
          fGenLambda.ctMC = std::sqrt(Sq(ov[0] - dv[0]) + Sq(ov[1] - dv[1]) + Sq(ov[2] - dv[2])) * lambda->M() / lambda->P();

          fGenLambda.flag = 0u;
          if (lambda->IsPrimary())
            fGenLambda.flag |= kPrimary;
          else if (lambda->IsSecondaryFromWeakDecay())
            FindWDLambdaMother(lambda);
          else fGenLambda.flag |= kSecondaryFromMaterial;
        }
        if (fOnlyTrueLambdas && fGenLambda.pdg == 0)
          continue;
      }

      fRecLambda->matter = v0->AlphaV0() > 0;
      auto proton = fRecLambda->matter ? pTrack : nTrack;
      auto pion = fRecLambda->matter ? nTrack : pTrack;

      int nSPDPr = 0u;
      int nSDDSSDPr = 0u;
      getITScls(proton, nSPDPr, nSDDSSDPr);
      if (nSDDSSDPr < fSDDSSDclsCut || nSPDPr < fSPDclsCut)
        continue;
      fRecLambda->pt = v0->Pt();
      fRecLambda->ptPr = proton->Pt();
      fRecLambda->eta = v0->Eta();
      fRecLambda->mass = fRecLambda->matter ? v0->MassLambda() : v0->MassAntiLambda();
      fRecLambda->ct = v0->Ct(kLambdaPdg, pv);
      fRecLambda->radius = v0->RadiusV0();
      fRecLambda->dcaV0PV = v0->DcaV0ToPrimVertex();
      fRecLambda->dcaPrPV = fRecLambda->matter ? v0->DcaPosToPrimVertex() : v0->DcaNegToPrimVertex();
      fRecLambda->dcaPiPV = fRecLambda->matter ? v0->DcaNegToPrimVertex() : v0->DcaPosToPrimVertex();
      fRecLambda->dcaV0tracks = v0->DcaV0Daughters();
      fRecLambda->cosPA = v0->CosPointingAngle(pv);
      fRecLambda->tpcNsigmaPi = fPID->NumberOfSigmasTPC(pion, AliPID::kPion);
      fRecLambda->tpcNsigmaPr = fPID->NumberOfSigmasTPC(proton, AliPID::kProton);
      fRecLambda->itsNsigmaPr = fPID->NumberOfSigmasITS(proton, AliPID::kProton);
      fRecLambda->tpcClV0Pi = pion->GetTPCsignalN();
      fRecLambda->tpcClV0Pr = proton->GetTPCsignalN();
      fRecLambda->hasTOFhit = !pTrack->GetTOFBunchCrossing(bField) || !nTrack->GetTOFBunchCrossing(bField);
      fRecLambda->hasITSrefit = nTrack->GetStatus() & AliVTrack::kITSrefit || pTrack->GetStatus() & AliVTrack::kITSrefit;

      //crossed raws
      double lCrosRawsPos = pTrack->GetTPCClusterInfo(2, 1);
      double lCrosRawsNeg = nTrack->GetTPCClusterInfo(2, 1);
      fV0LeastCRaws = std::min(lCrosRawsPos, lCrosRawsNeg);
      //crossed raws / Findable clusters
      double lCrosRawsOvFPos = lCrosRawsPos / ((double)(pTrack->GetTPCNclsF()));
      double lCrosRawsOvFNeg = lCrosRawsNeg / ((double)(nTrack->GetTPCNclsF()));
      fV0LeastCRawsOvF = std::min(lCrosRawsOvFPos, lCrosRawsOvFNeg);

      if (IsTopolSelectedLambda())
      {
        if (lambdaLabel != -1)
        {
          checkedLambdaLabel.push_back(lambdaLabel);
        }
        if (fFillLambdas)
          fTreeLambda->Fill();
        if (fFillLambdasBDTOut)
        {
          if (fRecLambda->ct < fCtPreselection || fRecLambda->ct > fMaxCt || fRecLambda->pt < fMinPt || fRecLambda->pt > fMaxPt || fRecLambda->radius < fRadiusPreselection || fRecLambda->radius > fRadiusOverflowCut || fRecLambda->dcaPiPV > fDCAV0piToPVOverflowCut || fRecLambda->dcaPrPV > fDCAV0prToPVOverflowCut || fRecLambda->dcaV0PV > fDCAV0toPVOverflowCut || fRecLambda->tpcClV0Pi < fTpcClV0PiPreselection || fRecLambda->tpcClV0Pr < fTpcClV0PrPreselection || fRecLambda->centrality < fMinCentrality || fRecLambda->centrality > fMaxCentrality) continue;
          int model_index = WhichBDT(fRecLambda->ct);
          if (model_index > (fCtBinsBDT.GetSize()-2)) {
            continue;
          }
          if (!fBDT[model_index])
          {
            AliError("ERROR: BDT not loaded, skip prediction ...\n");
            continue;
          }
          double features[]={fRecLambda->cosPA, fRecLambda->dcaV0tracks, fRecLambda->dcaPiPV, fRecLambda->dcaPrPV, fRecLambda->dcaV0PV, fRecLambda->tpcNsigmaPr, fRecLambda->radius};
          std::vector<double> bdt_out;
          fBDT[model_index]->Predict(features,7,bdt_out,false);
          if (bdt_out[0] > fBdtOutputBackgroundCut) {
            continue;
          }
          fRecLambdaBDTOut->bdtOutputBackground = bdt_out[0];
          fRecLambdaBDTOut->bdtOutputNonPrompt = bdt_out[1];
          fRecLambdaBDTOut->bdtOutputPrompt = bdt_out[2];
          fRecLambdaBDTOut->matter = fRecLambda->matter;
          fRecLambdaBDTOut->mass = fRecLambda->mass;
          fRecLambdaBDTOut->ct = fRecLambda->ct;
          fRecLambdaBDTOut->pt = fRecLambda->pt;
          fRecLambdaBDTOut->pileUpCheck = fRecLambda->hasTOFhit || fRecLambda->hasITSrefit;
          fTreeLambdaBDTOut->Fill();
        }
      }
    }
  }

  if (fFillK0s && rdmState < fK0sDownscaling)
  {
    fGenK0s.isReconstructed = true;
    for (int iV0{0}; iV0 < ev->GetNumberOfV0s(); ++iV0)
    {
      AliAODv0 *v0{ev->GetV0(iV0)};
      if (!v0)
        continue;
      if (v0->GetOnFlyStatus() != fUseOnTheFly)
        continue;

      fRecK0s->centrality = fEventCut.GetCentrality();
      fRecK0s->radius = v0->RadiusSecVtx();

      //get daughter tracks (positive, negative and bachelor)
      AliAODTrack *pTrack = dynamic_cast<AliAODTrack *>(v0->GetDaughter(0));
      AliAODTrack *nTrack = dynamic_cast<AliAODTrack *>(v0->GetDaughter(1));
      if (!pTrack || !nTrack)
      {
        AliWarning("ERROR: Could not retrieve one of the 2 AOD daughter tracks of the K0ss ...\n");
        continue;
      }

      if (!(pTrack->GetStatus() & AliVTrack::kTPCrefit) || !(nTrack->GetStatus() & AliVTrack::kTPCrefit) ||
          pTrack->GetTPCsignalN() < 50 || nTrack->GetTPCsignalN() < 50 ||
          std::abs(pTrack->Eta()) > 0.8 || std::abs(nTrack->Eta()) > 0.8 ||
          pTrack->Chi2perNDF() > 4 || nTrack->Chi2perNDF() > 4)
      {
        continue;
      }

      int K0sLabel{-1};
      if (fMC)
      {
        auto posPart = (AliAODMCParticle *)fMCEvent->GetTrack(std::abs(pTrack->GetLabel()));
        auto negPart = (AliAODMCParticle *)fMCEvent->GetTrack(std::abs(nTrack->GetLabel()));
        // Check K0s
        int labMothPos = posPart->GetMother();
        int labMothNeg = negPart->GetMother();
        auto K0s = (AliAODMCParticle *)fMCEvent->GetTrack(labMothNeg);
        if (K0s && labMothNeg == labMothPos && K0s->GetPdgCode() == kK0sPdg)
        {
          K0sLabel = labMothNeg;
          fGenK0s.ptMC = K0s->Pt();
          fGenK0s.etaMC = K0s->Eta();
          fGenK0s.yMC = K0s->Y();
          fGenK0s.isPrimary = K0s->IsPhysicalPrimary();
          double ov[3], dv[3];
          K0s->XvYvZv(ov);
          posPart->XvYvZv(dv);
          fGenK0s.ctMC = std::sqrt(Sq(ov[0] - dv[0]) + Sq(ov[1] - dv[1]) + Sq(ov[2] - dv[2])) * K0s->M() / K0s->P();

          fGenK0s.flag = 0u;
          if (K0s->IsPrimary())
            fGenK0s.flag |= kPrimary;
          else if (K0s->IsSecondaryFromWeakDecay())
            fGenK0s.flag |= kSecondaryFromWD;
          else fGenK0s.flag |= kSecondaryFromMaterial;
        }
        if (fOnlyTrueK0s && K0sLabel == -1)
          continue;
      }


      fRecK0s->pt = v0->Pt();
      fRecK0s->ptNeg = nTrack->Pt();
      fRecK0s->ptPos = pTrack->Pt();
      fRecK0s->eta = v0->Eta();
      fRecK0s->mass = v0->MassK0Short();
      fRecK0s->ct = v0->Ct(kK0sPdg, pv);
      fRecK0s->radius = v0->RadiusV0();
      fRecK0s->dcaV0PV = v0->DcaV0ToPrimVertex();
      fRecK0s->dcaPosPV = v0->DcaPosToPrimVertex();
      fRecK0s->dcaNegPV = v0->DcaNegToPrimVertex();
      fRecK0s->dcaV0tracks = v0->DcaV0Daughters();
      fRecK0s->cosPA = v0->CosPointingAngle(pv);
      fRecK0s->tpcNsigmaNeg = fPID->NumberOfSigmasTPC(nTrack, AliPID::kPion);
      fRecK0s->tpcNsigmaPos = fPID->NumberOfSigmasTPC(pTrack, AliPID::kPion);
      fRecK0s->tpcClV0Neg = nTrack->GetTPCsignalN();
      fRecK0s->tpcClV0Pos = pTrack->GetTPCsignalN();
      fRecK0s->hasTOFneg = HasTOF(nTrack);
      fRecK0s->hasTOFpos = HasTOF(pTrack);
      fRecK0s->hasTOFhit = !pTrack->GetTOFBunchCrossing(bField) || !nTrack->GetTOFBunchCrossing(bField);
      fRecK0s->hasITSrefit = nTrack->GetStatus() & AliVTrack::kITSrefit || pTrack->GetStatus() & AliVTrack::kITSrefit;

      //crossed raws
      double lCrosRawsPos = pTrack->GetTPCClusterInfo(2, 1);
      double lCrosRawsNeg = nTrack->GetTPCClusterInfo(2, 1);
      fV0LeastCRaws = std::min(lCrosRawsPos, lCrosRawsNeg);
      //crossed raws / Findable clusters
      double lCrosRawsOvFPos = lCrosRawsPos / ((double)(pTrack->GetTPCNclsF()));
      double lCrosRawsOvFNeg = lCrosRawsNeg / ((double)(nTrack->GetTPCNclsF()));
      fV0LeastCRawsOvF = std::min(lCrosRawsOvFPos, lCrosRawsOvFNeg);

      if (IsTopolSelectedK0s())
      {
        if (K0sLabel != -1)
        {
          checkedK0sLabel.push_back(K0sLabel);
        }
        if (fFillK0s)
          fTreeK0s->Fill();
      }
    }
  }


  if (fMC && (fFillCascades||fFillLambdas||fFillK0s||fFillCascadesTrain))
  {
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

    if (fFillCascades || fFillCascadesTrain)
    {
      fGenCascade.isReconstructed = false;
      fGenCascadeTrain.isReconstructed = false;
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
        fGenCascadeTrain.pdg = track->GetPdgCode();
        double pv[3], sv[3];
        track->XvYvZv(pv);
        bool goodDecay{false};
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
            goodDecay = true;
            break;
          }
        }
        if (!goodDecay)
          continue;
        fGenCascade.ctMC = std::sqrt(Sq(pv[0] - sv[0]) + Sq(pv[1] - sv[1]) + Sq(pv[2] - sv[2])) * track->M() / track->P();
        fGenCascade.flag = 0u;
        if (track->IsPrimary())
          fGenCascade.flag |= kPrimary;
        else
          fGenCascade.flag |= track->IsSecondaryFromWeakDecay() ? kSecondaryFromWD : kSecondaryFromMaterial;

        //train cascade
        fGenCascadeTrain.ptMC = fGenCascade.ptMC;
        fGenCascadeTrain.etaMC = fGenCascade.etaMC;
        fGenCascadeTrain.ctMC = fGenCascade.ctMC;
        fGenCascadeTrain.yMC = fGenCascade.yMC;
        fGenCascadeTrain.flag = fGenCascade.flag;

        if (fFillCascadesTrain)
          fTreeTrain->Fill();
        if (fFillCascades)
          fTree->Fill();
      }
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
        double ov[3], dv[3];
        track->XvYvZv(ov);
        bool neutralDecay{true};
        for (int iD = track->GetDaughterFirst(); iD <= track->GetDaughterLast(); iD++)
        {
          auto daugh = (AliAODMCParticle *)fMCEvent->GetTrack(iD);
          if (!daugh)
          {
            continue;
          }
          if (std::abs(daugh->GetPdgCode()) == AliPID::ParticleCode(AliPID::kProton))
          {
            neutralDecay = false;
            daugh->XvYvZv(dv);
            break;
          }
        }
        if (neutralDecay)
          continue;
        fGenLambda.flag = 0u;
        if (track->IsPrimary())
          fGenLambda.flag |= kPrimary;
        else if (track->IsSecondaryFromWeakDecay())
          FindWDLambdaMother(track);
        else fGenLambda.flag |= kSecondaryFromMaterial;

        fGenLambda.ctMC = std::sqrt(Sq(ov[0] - dv[0]) + Sq(ov[1] - dv[1]) + Sq(ov[2] - dv[2])) * track->M() / track->P();
        fTreeLambda->Fill();
      }
    }

    if (fFillK0s && rdmState < fK0sDownscaling)
    {
      fGenK0s.isReconstructed = false;
      //loop on generated
      for (int iT{0}; iT < fMCEvent->GetNumberOfTracks(); ++iT)
      {
        auto track = (AliAODMCParticle *)fMCEvent->GetTrack(iT);
        int pdg = std::abs(track->GetPdgCode());
        if (pdg != kK0sPdg)
        {
          continue;
        }
        if (std::find(checkedK0sLabel.begin(), checkedK0sLabel.end(), iT) != checkedK0sLabel.end())
        {
          continue;
        }

        if (std::abs(track->Y()) > fCutY || AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(iT, header, MCTrackArray))
        { //removal of OOB pileup, cut on Y and PhysPrim
          continue;
        }
        fGenK0s.ptMC = track->Pt();
        fGenK0s.etaMC = track->Eta();
        fGenK0s.yMC = track->Y();
        fGenK0s.isPrimary = track->IsPhysicalPrimary();
        double ov[3], dv[3];
        track->XvYvZv(ov);
        bool neutralDecay{true};
        for (int iD = track->GetDaughterFirst(); iD <= track->GetDaughterLast(); iD++)
        {
          auto daugh = (AliAODMCParticle *)fMCEvent->GetTrack(iD);
          if (!daugh)
          {
            continue;
          }
          if (std::abs(daugh->GetPdgCode()) == AliPID::ParticleCode(AliPID::kPion))
          {
            neutralDecay = false;
            daugh->XvYvZv(dv);
            break;
          }
        }
        if (neutralDecay)
          continue;
        fGenK0s.flag = 0u;
        if (track->IsPrimary())
          fGenK0s.flag |= kPrimary;
        else if (track->IsSecondaryFromWeakDecay())
          fGenK0s.flag |= kSecondaryFromWD;
        else fGenK0s.flag |= kSecondaryFromMaterial;

        fGenK0s.ctMC = std::sqrt(Sq(ov[0] - dv[0]) + Sq(ov[1] - dv[1]) + Sq(ov[2] - dv[2])) * track->M() / track->P();
        fTreeK0s->Fill();
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
      mgr->CreateContainer(Form("%s_treeCascadesTrain", tskname.Data()), TTree::Class(),
                           AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
  coutput3->SetSpecialOutput();

  AliAnalysisDataContainer *coutput4 =
      mgr->CreateContainer(Form("%s_treeLambda", tskname.Data()), TTree::Class(),
                           AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
  coutput4->SetSpecialOutput();

  AliAnalysisDataContainer *coutput5 =
      mgr->CreateContainer(Form("%s_treeLambdaBDTOut", tskname.Data()), TTree::Class(),
                           AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
  coutput5->SetSpecialOutput();

  AliAnalysisDataContainer *coutput6 =
      mgr->CreateContainer(Form("%s_treeK0s", tskname.Data()), TTree::Class(),
                           AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
  coutput6->SetSpecialOutput();

  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, coutput1);
  mgr->ConnectOutput(task, 2, coutput2);
  mgr->ConnectOutput(task, 3, coutput3);
  mgr->ConnectOutput(task, 4, coutput4);
  mgr->ConnectOutput(task, 5, coutput5);
  mgr->ConnectOutput(task, 6, coutput6);
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
         fCascLeastCRawsOvF > fCutRowsOvF &&
         (!fUseAdditionalCuts || (fUseAdditionalCuts && (fRecCascade->bachBarCosPA < fCutBachBarCosPA) && (fRecCascade->hasTOFhit || fRecCascade->hasITSrefit)));
}

bool AliAnalysisTaskStrangenessRatios::IsTopolSelectedLambda()
{
  return fRecLambda->radius > fCutRadius[2] &&
         fRecLambda->cosPA > fCosPALambda &&
         fRecLambda->dcaPrPV > fCutDCALambdaPrToPV &&
         fRecLambda->dcaPiPV > fCutDCALambdaPiToPV &&
         fRecLambda->dcaV0tracks < fCutDCAV0tracks &&
         std::abs(Eta2y(fRecLambda->pt, kLambdaMass, fRecLambda->eta)) < fCutY &&
         fRecLambda->mass > fCutLambdaMass[0] && fRecLambda->mass < fCutLambdaMass[1] &&
         std::abs(fRecLambda->tpcNsigmaPi) < fCutNsigmaTPC &&
         std::abs(fRecLambda->tpcNsigmaPr) < fCutNsigmaTPC &&
         fV0LeastCRaws > fCutTPCrows &&
         fV0LeastCRawsOvF > fCutRowsOvF;
}

bool AliAnalysisTaskStrangenessRatios::IsTopolSelectedK0s()
{
  return fRecK0s->radius > fCutRadius[2] &&
         fRecK0s->cosPA > fCosPAK0s &&
         fRecK0s->dcaPosPV > fCutDCAK0sProngToPV &&
         fRecK0s->dcaNegPV > fCutDCAK0sProngToPV &&
         fRecK0s->dcaV0tracks < fCutDCAV0tracks &&
         std::abs(Eta2y(fRecK0s->pt, kK0sMass, fRecK0s->eta)) < fCutY &&
         fRecK0s->mass > fCutK0sMass[0] && fRecK0s->mass < fCutK0sMass[1] &&
         std::abs(fRecK0s->tpcNsigmaPos) < fCutNsigmaTPC &&
         std::abs(fRecK0s->tpcNsigmaNeg) < fCutNsigmaTPC &&
         fV0LeastCRaws > fCutTPCrows &&
         fV0LeastCRawsOvF > fCutRowsOvF;
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
  PostData(3, fTreeTrain);
  PostData(4, fTreeLambda);
  PostData(5, fTreeLambdaBDTOut);
  PostData(6, fTreeK0s);
}

Bool_t AliAnalysisTaskStrangenessRatios::UserNotify()
{
  TString cfn{CurrentFileName()};
  AliInfo(Form("Setting hash for file %s", cfn.Data()));

  gRandom->SetSeed(cfn.Hash());
  return true;
}

int AliAnalysisTaskStrangenessRatios::WhichBDT(double ct)
{
  int iB=0;
  while ( (ct < fCtBinsBDT[iB] || ct > fCtBinsBDT[iB+1]) && iB < (fCtBinsBDT.GetSize()-1))
  {
    ++iB;
  }
  return iB;
}

void AliAnalysisTaskStrangenessRatios::FindWDLambdaMother(AliAODMCParticle *track)
{
  if (!track) return;
  int mothLambda = track->GetMother();
  auto mother = (AliAODMCParticle *)fMCEvent->GetTrack(mothLambda);
  switch (std::abs(mother->GetPdgCode())){
    case kXiPdg:
    {
      if (mother->IsPhysicalPrimary())
        fGenLambda.flag |= kSecondaryFromWDXi;
      else {
        int moth2Lambda = mother->GetMother();
        auto mother2 = (AliAODMCParticle *)fMCEvent->GetTrack(moth2Lambda);
        if (std::abs(mother2->GetPdgCode()) == kOmegaPdg && mother2->IsPhysicalPrimary()) {
          mother = (AliAODMCParticle *)fMCEvent->GetTrack(moth2Lambda);
          fGenLambda.flag |= kSecondaryFromWDOmega;
        }
        else
          fGenLambda.flag |= kSecondaryFromWD;
      }
    }
    break;
    case kOmegaPdg:
    {
      if (mother->IsPhysicalPrimary())
        fGenLambda.flag |= kSecondaryFromWDOmega;
      else
        fGenLambda.flag |= kSecondaryFromWD;
    }
    break;
    default:
      fGenLambda.flag |= kSecondaryFromWD;
    break;
  }
  fGenLambda.ptMotherMC = mother->Pt();
  double ovMother[3], dvMother[3];
  for (int iD = mother->GetDaughterFirst(); iD <= mother->GetDaughterLast(); iD++)
  {
    mother->XvYvZv(ovMother);
    auto daugh = (AliAODMCParticle *)fMCEvent->GetTrack(iD);
    if (!daugh)
    {
      continue;
    }
    if (std::abs(daugh->GetPdgCode()) == kLambdaPdg)
    {
      daugh->XvYvZv(dvMother);
      break;
    }
  }
  fGenLambda.ctMotherMC = std::sqrt(Sq(ovMother[0] - dvMother[0]) + Sq(ovMother[1] - dvMother[1]) + Sq(ovMother[2] - dvMother[2])) * mother->M() / mother->P();
}

bool AliAnalysisTaskStrangenessRatios::HasTOF(AliVTrack *track) {
  bool hasTOFout  = track->GetStatus() & AliVTrack::kTOFout;
  bool hasTOFtime = track->GetStatus() & AliVTrack::kTIME;
  const float len = track->GetIntegratedLength();
  return hasTOFout && hasTOFtime && (len > 350.);
}