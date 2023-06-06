#include "AliAnalysisTaskKaonXiCorrelation.h"

#include "AliAnalysisDataContainer.h"
#include "AliLog.h"

#include <algorithm>
#include <cmath>
#include <string>
using std::string;

// ROOT includes
#include <TAxis.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TList.h>
#include <TRandom3.h>
#include <TTree.h>

// ALIROOT includes
#include "AliAnalysisManager.h"
#include "AliAODTrack.h"
#include "AliMCEvent.h"
#include "AliInputEventHandler.h"
#include "AliPIDResponse.h"
#include "AliAODEvent.h"
#include "AliVEventHandler.h"
#include "AliAODMCHeader.h"
#include "AliAnalysisUtils.h"

#include "AliAODVZERO.h"

///\cond CLASSIMP
ClassImp(AliAnalysisTaskKaonXiCorrelation);
///\endcond

namespace
{

  double Sq(double x)
  {
    return x * x;
  }

  constexpr int kKaonPdg{321};
  constexpr int kXiPdg{3312};
  constexpr int kLambdaPdg(3122);
  constexpr double kOmegaMass{1.67245};
  constexpr double kXiMass{1.32171};
  constexpr double kLambdaMass{1.115683};
  constexpr double kKaonMass{0.493677};
  constexpr double kcTauXi{4.91359839};

}

/// Standard and default constructor of the class.
///
/// \param taskname Name of the task
/// \param partname Name of the analysed particle
///
AliAnalysisTaskKaonXiCorrelation::AliAnalysisTaskKaonXiCorrelation(bool isMC, TString taskname) : AliAnalysisTaskSE(taskname.Data()),
                                                                                                  fEventCuts{false},
                                                                                                  fMC{isMC}
{
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}

/// Standard destructor
///
AliAnalysisTaskKaonXiCorrelation::~AliAnalysisTaskKaonXiCorrelation()
{
  if (AliAnalysisManager::GetAnalysisManager()->IsProofMode())
    return;
  if (fList)
    delete fList;
  if (fTree)
    delete fTree;
  if (!fMC)
  {
    delete fXi;
    delete fKaon;
  }
}

/// This function creates all the histograms and all the objects in general used during the analysis
/// \return void
///
void AliAnalysisTaskKaonXiCorrelation::UserCreateOutputObjects()
{

  OpenFile(2);
  fXi = fMC ? &fGenXi : new MiniXi;
  fKaon = fMC ? &fGenKaon : new MiniKaon;
  fTree = new TTree("StrangenessTree", "Tree");
  fTree->Branch("MiniCollision", &fRecCollision);
  if (fMC)
  {
    fTree->Branch("MiniXiMC", &fGenCascades);
    fTree->Branch("MiniKaonMC", &fGenKaons);
  }
  else
  {
    fTree->Branch("MiniXi", &fRecCascades);
    fTree->Branch("MiniKaon", &fRecKaons);
  }

  for (int iB=0; iB<(fPtBinsBDT.GetSize()-1); ++iB)
  {
    fBDT.push_back(new AliExternalBDT());
    if (!fBDT[iB]->LoadXGBoostModel(Form("%s%.1f_%.1f.model", fBDTPath.data(), fPtBinsBDT[iB], fPtBinsBDT[iB+1])))
    {
      fBDT[iB] = nullptr;
    }
  }

  const char* det[] = {"ITS", "TPC", "TOF"};
  TFile calib_file(fCustomPidPath.data());
  for (int iD = 0; iD < 3; ++iD)
  {
    fCustomPidCalib[iD] = (TH3F*)calib_file.Get(Form("fCustom%spidCalib", det[iD]));
  }

  fList = new TList();
  fList->SetOwner(kTRUE);
  fEventCuts.AddQAplotsToList(fList);

  PostAllData();
}

/// This is the function that is evaluated for each event. The analysis code stays here.
///
/// \param options Deprecated parameter
/// \return void
///
void AliAnalysisTaskKaonXiCorrelation::UserExec(Option_t *)
{
  AliAODEvent *ev = (AliAODEvent *)InputEvent();
  if (!fEventCuts.AcceptEvent(ev))
  {
    PostAllData();
    return;
  }

  if (fMC)
  {
    fMCEvent = MCEvent();
  }

  double bField{ev->GetMagneticField()};
  fRecCollision.fCent = fEventCuts.GetCentrality(fEstimator);
  if (fRecCollision.fCent < fMinCentrality || fRecCollision.fCent > fMaxCentrality)
  {
    PostAllData();
    return;
  }

  double pv[3];
  fEventCuts.GetPrimaryVertex()->GetXYZ(pv);
  fRecCollision.fZ = pv[2];

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *handl = (AliInputEventHandler *)mgr->GetInputEventHandler();

  unsigned char tgr = 0x0;

  if (handl->IsEventSelected() & AliVEvent::kINT7)
    tgr |= kINT7;
  if (handl->IsEventSelected() & AliVEvent::kCentral)
    tgr |= kCentral;
  if (handl->IsEventSelected() & AliVEvent::kSemiCentral)
    tgr |= kSemiCentral;
  int magField = bField > 0 ? kPositiveB : 0;

  fRecCollision.fTrigger = tgr + magField;

  AliAODVZERO *vzero = ev->GetVZEROData();
  fRecCollision.fV0MAmp = 0;
  for (int i = 0; i < 64; ++i)
  {
    fRecCollision.fV0MAmp += vzero->GetMultiplicity(i);
  }
  fRecCollision.fNTrk = 0;

  fPID = handl->GetPIDResponse();

  std::vector<int> checkedLabelCasc, checkedLabelKaon;
  std::vector<XiDaughter> xiDaughters;

  fGenCascades.clear();
  fGenKaons.clear();
  fRecCascades.clear();
  fRecKaons.clear();

  if (fFillKaons)
  {
    for (int iT{0}; iT < ev->GetNumberOfTracks(); iT++)
    {
      AliAODTrack* aodTrack = dynamic_cast<AliAODTrack *>(ev->GetTrack(iT));
      if (!aodTrack)
      {
        AliWarning("ERROR: Could not retrieve track ...\n");
        continue;
      }

      int nSPD = 0;
      int nSDD = 0;
      int nSSD = 0;
      int nITS = GetITScls(aodTrack, nSPD, nSDD, nSSD);
      float dca[2]{0., 0.};
      aodTrack->GetImpactParameters(dca[0], dca[1]);
      double dcaMag = std::sqrt(dca[1] * dca[1] + dca[0] * dca[0]);
      bool tof = HasTOF(aodTrack);

      double itsNsigma = GetCustomNsigma(aodTrack, fRecCollision.fCent, 0);
      double tpcNsigma = GetCustomNsigma(aodTrack, fRecCollision.fCent, 1);
      double tofNsigma = tof ? GetCustomNsigma(aodTrack, fRecCollision.fCent, 2) : -999.f;

      if (!aodTrack->TestFilterBit(fFilterBit) && fFilterBit) continue;

      if (!(aodTrack->GetStatus() & AliVTrack::kTPCrefit) ||
          !(aodTrack->GetStatus() & AliVTrack::kITSrefit) ||
          std::abs(Eta2y(aodTrack->Pt(), kKaonMass, aodTrack->Eta())) > fCutY ||
          std::abs(aodTrack->Eta()) > 0.8 ||
          aodTrack->Chi2perNDF() > fCutChi2[2] ||
          aodTrack->GetTPCsignalN() < fCutTPCcls[2] ||
          aodTrack->GetITSchi2() > fCutMaxITSChi2 ||
          nITS < fCutITSrecPoints ||
          nSPD < fCutSPDrecPoints ||
          dcaMag > fCutDCA[2]
        )
      {
        continue;
      }
      
      fRecCollision.fNTrk += 1;

      if ( std::abs(tpcNsigma) > fCutKaonNsigmaTPC ||
          (aodTrack->Pt() > fPtTofCut && std::abs(tofNsigma) > fCutKaonNsigmaTOF) ||
          !( (!fUseITSpid || aodTrack->Pt() > fCutPtITSpid) || (fUseITSpid && std::abs(itsNsigma) < fCutKaonNsigmaITS && ( nSDD + nSSD > fCutSDDSSDrecPoints) && aodTrack->Pt() < fCutPtITSpid) )
        )
      {
        continue;
      }

      if (aodTrack->Pt() > fMaxPtKaon) continue;

      if (fMC)
      {
        auto part = (AliAODMCParticle *)fMCEvent->GetTrack(std::abs(aodTrack->GetLabel()));
        int pdg = part->GetPdgCode();
        if (std::abs(pdg) != kKaonPdg)
        {
          continue;
        }
        fGenKaon.fIsReconstructed = true;
        fGenKaon.fEtaMC = part->Eta();
        fGenKaon.fPtMC = pdg > 0 ? part->Pt() : -part->Pt();
        bool physPrim = part->IsPhysicalPrimary();
        fGenKaon.fFlag = 0u;
        if (physPrim)
          fGenKaon.fFlag |= kPrimary;
        else
          fGenKaon.fFlag |= part->IsSecondaryFromWeakDecay() ? kSecondaryFromWD : kSecondaryFromMaterial;
        checkedLabelKaon.push_back(std::abs(part->GetLabel()));
      }
      bool charge = aodTrack->Charge() > 0;
      fKaon->fPt = charge ? aodTrack->Pt() : -aodTrack->Pt();
      fKaon->fEta = aodTrack->Eta();
      fKaon->fNsigmaITS = itsNsigma;
      fKaon->fNsigmaTPC = tpcNsigma;
      fKaon->fNsigmaTOF = tofNsigma;
      fKaon->fCutBitMap = 0u;
      if (dcaMag < fCutDCA[0]) fKaon->fCutBitMap |= kDCAtightCut;
      else if (dcaMag < fCutDCA[1]) fKaon->fCutBitMap |= kDCAmidCut;
      if (aodTrack->GetTPCsignalN() > fCutTPCcls[0]) fKaon->fCutBitMap |= kTPCclsTightCut;
      else if (aodTrack->GetTPCsignalN() > fCutTPCcls[1]) fKaon->fCutBitMap |= kTPCclsMidCut;
      if (aodTrack->Chi2perNDF() < fCutChi2[0]) fKaon->fCutBitMap |= kChi2TightCut;
      else if (aodTrack->Chi2perNDF() < fCutChi2[1]) fKaon->fCutBitMap |= kChi2MidCut;

      if (!fMC)
        fRecKaons.push_back(*fKaon);
      else
        fGenKaons.push_back(fGenKaon);
      
      fRecCollision.fNTrk -= 1;
    }
  }

  if (fFillCascades)
  {
    for (int iCasc{0}; iCasc < ev->GetNumberOfCascades(); iCasc++)
    {
      AliAODcascade *casc = ev->GetCascade(iCasc);
      if (!casc)
        continue;
      if (casc->GetOnFlyStatus() != fUseOnTheFly)
        continue;

      // get daughter tracks (positive, negative and bachelor)
      AliAODTrack *pTrackCasc = dynamic_cast<AliAODTrack *>(casc->GetDaughter(0));
      AliAODTrack *nTrackCasc = dynamic_cast<AliAODTrack *>(casc->GetDaughter(1));
      AliAODTrack *bTrackCasc = dynamic_cast<AliAODTrack *>(casc->GetDecayVertexXi()->GetDaughter(0));
      if (!pTrackCasc || !nTrackCasc || !bTrackCasc)
      {
        AliWarning("ERROR: Could not retrieve one of the 3 AOD daughter tracks of the cascade ...\n");
        continue;
      }

      if (!(pTrackCasc->GetStatus() & AliVTrack::kTPCrefit) || !(nTrackCasc->GetStatus() & AliVTrack::kTPCrefit) || !(bTrackCasc->GetStatus() & AliVTrack::kTPCrefit) ||
          pTrackCasc->GetTPCsignalN() < 50 || nTrackCasc->GetTPCsignalN() < 50 || bTrackCasc->GetTPCsignalN() < 50 ||
          std::abs(pTrackCasc->Eta()) > 0.8 || std::abs(nTrackCasc->Eta()) > 0.8 || std::abs(bTrackCasc->Eta()) > 0.8 ||
          pTrackCasc->Chi2perNDF() > 4 || nTrackCasc->Chi2perNDF() > 4 || bTrackCasc->Chi2perNDF() > 4)
      {
        continue;
      }

      int labMothBac = -9999;
      int pdgCascade = -1;
      bool physPrim = true;
      if (fMC)
      {
        int pdg = 0;
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
          if (labMothLam == labMothBac && pdgCascade == kXiPdg)
          {
            pdg = cascade->GetPdgCode();
            fGenXi.fIsReconstructed = true;
            fGenXi.fPtMC = pdg > 0 ? cascade->Pt() : -cascade->Pt();
            fGenXi.fEtaMC = cascade->Eta();
            physPrim = cascade->IsPhysicalPrimary();
            fGenXi.fFlag = 0u;
            if (physPrim)
              fGenXi.fFlag |= kPrimary;
            else
              fGenXi.fFlag |= cascade->IsSecondaryFromWeakDecay() ? kSecondaryFromWD : kSecondaryFromMaterial;
            double pv[3], sv[3];
            cascade->XvYvZv(pv);
            bacPart->XvYvZv(sv);
            checkedLabelCasc.push_back(labMothBac);
          }
        }
        if (fOnlyTrueCandidates && pdg == 0)
          continue;
      }

      bool matter = casc->AlphaV0() > 0;
      double vtxCasc[3]{casc->DecayVertexXiX(), casc->DecayVertexXiY(), casc->DecayVertexXiZ()};
      double radius = std::hypot(vtxCasc[0], vtxCasc[1]);
      double radiusV0 = casc->RadiusSecVtx();

      double tpcNsigmaBach = fPID->NumberOfSigmasTPC(bTrackCasc, AliPID::kPion);
      double tpcNsigmaV0Pi = fPID->NumberOfSigmasTPC(matter ? nTrackCasc : pTrackCasc, AliPID::kPion);
      double tpcNsigmaV0Pr = fPID->NumberOfSigmasTPC(matter ? pTrackCasc : nTrackCasc, AliPID::kProton);

      double tpcClBach = bTrackCasc->GetTPCsignalN();
      double tpcClV0Pi = (matter ? nTrackCasc : pTrackCasc)->GetTPCsignalN();
      double tpcClV0Pr = (matter ? pTrackCasc : nTrackCasc)->GetTPCsignalN();

      // DCA info
      double dcaBachV0 = casc->DcaXiDaughters();
      double dcaBachPV = casc->DcaBachToPrimVertex();
      double dcaV0prPV = matter ? casc->DcaPosToPrimVertex() : casc->DcaNegToPrimVertex();
      double dcaV0piPV = matter ? casc->DcaNegToPrimVertex() : casc->DcaPosToPrimVertex();
      double dcaV0tracks = casc->DcaV0Daughters();
      double dcaV0PV = casc->DcaV0ToPrimVertex();

      // cascade and V0 cosine of pointing angle
      double cosPA = casc->CosPointingAngleXi(pv[0], pv[1], pv[2]);
      double cosPAV0 = casc->CosPointingAngle(pv);

      // TOF matching
      bool hasTOFhit = !pTrackCasc->GetTOFBunchCrossing(bField) || !nTrackCasc->GetTOFBunchCrossing(bField) || !bTrackCasc->GetTOFBunchCrossing(bField);

      // track status: ( fCasc_NegTrackStatus & AliESDtrack::kITSrefit ) is the codition to check kITSrefit
      bool hasITSrefit = (nTrackCasc->GetStatus() & AliVTrack::kITSrefit) || (pTrackCasc->GetStatus() & AliVTrack::kITSrefit) || (bTrackCasc->GetStatus() & AliVTrack::kITSrefit);

      double V0invMassDelta = ((matter) ? casc->MassLambda() : casc->MassAntiLambda()) - kLambdaMass;
      double competingMass = std::abs(casc->MassOmega() - kOmegaMass);

      // transverse momentum and eta
      double pt = std::sqrt(casc->Pt2Xi());

      //////////////////////////////
      // crossed rows
      double lCrosRowsPos = pTrackCasc->GetTPCClusterInfo(2, 1);
      double lCrosRowsNeg = nTrackCasc->GetTPCClusterInfo(2, 1);
      double lCrosRowsBac = bTrackCasc->GetTPCClusterInfo(2, 1);
      fCascLeastCRows = (int)(lCrosRowsPos < lCrosRowsNeg ? std::min(lCrosRowsPos, lCrosRowsBac) : std::min(lCrosRowsNeg, lCrosRowsBac));
      // crossed rows / Findable clusters
      double lCrosRowsOvFPos = lCrosRowsPos / ((double)(pTrackCasc->GetTPCNclsF()));
      double lCrosRowsOvFNeg = lCrosRowsNeg / ((double)(nTrackCasc->GetTPCNclsF()));
      double lCrosRowsOvFBac = lCrosRowsBac / ((double)(bTrackCasc->GetTPCNclsF()));
      fCascLeastCRowsOvF = lCrosRowsOvFPos < lCrosRowsOvFNeg ? std::min(lCrosRowsOvFPos, lCrosRowsOvFBac) : std::min(lCrosRowsOvFNeg, lCrosRowsOvFBac);
      ///////////////////////////////

      // calculate DCA Bachelor-Baryon to remove "bump" structure in InvMass
      double bachBarCosPA = casc->BachBaryonCosPA();
      double p = std::sqrt(casc->Ptot2Xi());
      double eta = 0.5 * std::log((p + casc->MomXiZ()) / (p - casc->MomXiZ() + 1.e-16));

      // distance over total momentum
      double lOverP = std::sqrt((Sq(vtxCasc[0] - pv[0]) + Sq(vtxCasc[1] - pv[1]) + Sq(vtxCasc[2] - pv[2])) / (casc->Ptot2Xi() + 1e-10));
      double ct = lOverP * kXiMass;

      bool isTopolSelected = radius > fCutRadius &&
            radiusV0 > fCutRadiusV0 &&
            dcaBachPV > fCutDCABachToPV &&
            dcaV0PV > fCutDCAV0toPV &&
            dcaV0piPV > fCutDCAV0piToPV &&
            dcaV0prPV > fCutDCAV0prToPV &&
            dcaV0tracks < fCutDCAV0tracks &&
            dcaBachV0 < fCutDCABachToV0 &&
            cosPA > fCutCosPA &&
            cosPAV0 > fCutCosPAV0 &&
            bachBarCosPA < fCutBachBarCosPA &&
            std::abs(Eta2y(pt, kXiMass, eta)) < fCutY &&
            std::abs(tpcNsigmaBach) < fCutNsigmaTPC &&
            std::abs(tpcNsigmaV0Pr) < fCutNsigmaTPC &&
            std::abs(tpcNsigmaV0Pi) < fCutNsigmaTPC &&
            ct < fCutCt * kcTauXi &&
            tpcClBach > fCutTPCclu &&
            tpcClV0Pi > fCutTPCclu &&
            tpcClV0Pr > fCutTPCclu &&
            fCascLeastCRows > fCutTPCrows &&
            fCascLeastCRowsOvF > fCutRowsOvF;

      if (std::abs(casc->MassXi() - kXiMass) * 1000 < 30)
      {

        if (pt < fMinPt || pt > fMaxPt ||
            radius > fRadiusOverflowCut || radiusV0 > fRadiusV0OverflowCut ||
            dcaV0piPV > fDCAV0piToPVOverflowCut || dcaV0prPV > fDCAV0prToPVOverflowCut || dcaBachPV > fDCABachToPVOverflowCut ||
            dcaV0PV > fDCAV0toPVOverflowCut)
        {
          isTopolSelected = false;
        }

        if (physPrim && isTopolSelected)
        {
          fXi->fPt = matter ? pt : -pt;
          fXi->fEta = eta;
          fXi->fMass = casc->MassXi();
          fXi->fBdtOut = -1.;
          fXi->fRecFlag = 0u;
          if (hasTOFhit || hasITSrefit)
            fXi->fRecFlag |= kHasTOFhitOrITSrefit;
          if (competingMass > fCutCompetingMass)
            fXi->fRecFlag |= kCompetingMassCut;

          if (fMC && !fApplyBdtToMC)
          {
            fGenCascades.push_back(fGenXi);
            continue;
          }
          else
          {
            int model_index = WhichBDT(pt);
            if (model_index > (fPtBinsBDT.GetSize()-2)) {
              continue;
            }
            if (!fBDT[model_index])
            {
              AliError("ERROR: BDT not loaded, skip prediction ...\n");
              continue;
            }
            double features[]={dcaV0tracks, dcaBachV0, cosPA, cosPAV0, tpcNsigmaV0Pr, dcaBachPV, dcaV0PV, dcaV0piPV, dcaV0prPV};
            std::vector<double> bdt_out;
            fBDT[model_index]->Predict(features, fNFeatures, bdt_out, false);
            if (bdt_out[0] < fBdtOutCut && !fMC)
            {
              continue;
            }
            fXi->fBdtOut = bdt_out[0];

            if (std::abs(fXi->fMass - kXiMass) < fCascMassWindow)
            {
              XiDaughter pDaughter(pTrackCasc->Px(), pTrackCasc->Py(), pTrackCasc->Pz(), 1);
              XiDaughter nDaughter(nTrackCasc->Px(), nTrackCasc->Py(), nTrackCasc->Pz(), -1);
              XiDaughter bDaughter(bTrackCasc->Px(), bTrackCasc->Py(), bTrackCasc->Pz(), matter ? -1 : 1);
              xiDaughters.push_back(pDaughter); xiDaughters.push_back(nDaughter); xiDaughters.push_back(bDaughter);
            }

            if (!fMC)
              fRecCascades.push_back(*fXi);
            else
              fGenCascades.push_back(fGenXi);
          }
        }
        else if (fMC && std::find(checkedLabelCasc.begin(), checkedLabelCasc.end(), labMothBac) != checkedLabelCasc.end() && (pdgCascade == kXiPdg))
        {
          checkedLabelCasc.erase(std::find(checkedLabelCasc.begin(), checkedLabelCasc.end(), labMothBac)); //checked particles that didn't pass the topological cut (have to be filled later)
        }
      }
    }
    if (HasTwoXiFromSameDaughters(xiDaughters)){
      fRecCollision.fTrigger |= kHasTwoXiFromSameDaughter;
    }
  }

  if (fMC && (fFillCascades || fFillKaons))
  {
    // OOB pileup
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

    if (fFillKaons)
    {
      fGenKaon.fIsReconstructed = false;
      // loop on generated
      for (int iT{0}; iT < fMCEvent->GetNumberOfTracks(); ++iT)
      {
        auto track = (AliAODMCParticle *)fMCEvent->GetTrack(iT);
        int pdg = track->GetPdgCode();
        if (std::abs(pdg) != kKaonPdg)
        {
          continue;
        }
        if (std::find(checkedLabelKaon.begin(), checkedLabelKaon.end(), iT) != checkedLabelKaon.end())
        {
          continue;
        }
        if (std::abs(track->Y()) > fCutY || AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(iT, header, MCTrackArray))
        { // removal of OOB pileup, cut on Y and PhysPrim
          continue;
        }
        fGenKaon.fPtMC = pdg > 0 ? track->Pt() : -track->Pt();
        fGenKaon.fEtaMC = track->Eta();
        fGenKaon.fFlag = 0u;
        if (track->IsPhysicalPrimary())
          fGenKaon.fFlag |= kPrimary;
        else
          fGenKaon.fFlag |= track->IsSecondaryFromWeakDecay() ? kSecondaryFromWD : kSecondaryFromMaterial;
        fGenKaons.push_back(fGenKaon);
      }
    }

    if (fFillCascades)
    {
      fGenXi.fIsReconstructed = false;
      // loop on generated
      for (int iT{0}; iT < fMCEvent->GetNumberOfTracks(); ++iT)
      {
        auto track = (AliAODMCParticle *)fMCEvent->GetTrack(iT);
        int pdg = track->GetPdgCode();
        if (std::abs(pdg) != kXiPdg)
        {
          continue;
        }
        if (std::find(checkedLabelCasc.begin(), checkedLabelCasc.end(), iT) != checkedLabelCasc.end())
        {
          continue;
        }
        if (std::abs(track->Y()) > fCutY || !track->IsPhysicalPrimary() || AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(iT, header, MCTrackArray))
        { // removal of OOB pileup, cut on Y and PhysPrim
          continue;
        }
        fGenXi.fPtMC = pdg > 0 ? track->Pt() : -track->Pt();
        fGenXi.fEtaMC = track->Eta();
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
        fGenXi.fFlag = 0u;
        if (track->IsPhysicalPrimary())
          fGenXi.fFlag |= kPrimary;
        else
          fGenXi.fFlag |= track->IsSecondaryFromWeakDecay() ? kSecondaryFromWD : kSecondaryFromMaterial;
        fGenCascades.push_back(fGenXi);
      }
    }
  }

  fTree->Fill();
  PostAllData();
}

AliAnalysisTaskKaonXiCorrelation *AliAnalysisTaskKaonXiCorrelation::AddTask(bool isMC, TString tskname, TString suffix)
{
  // Get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskKaonXiCorrelation", "No analysis manager found.");
    return nullptr;
  }

  // Check the analysis type using the event handlers connected to the analysis
  // manager.
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskKaonXiCorrelation", "This task requires an input event handler");
    return nullptr;
  }

  tskname.Append(suffix.Data());
  AliAnalysisTaskKaonXiCorrelation *task = new AliAnalysisTaskKaonXiCorrelation(isMC, tskname.Data());

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(
      Form("%s_summary", tskname.Data()), TList::Class(),
      AliAnalysisManager::kOutputContainer, "AnalysisResults.root");

  AliAnalysisDataContainer *coutput2 =
      mgr->CreateContainer(Form("%s_treeStrangeness", tskname.Data()), TTree::Class(),
                           AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
  coutput2->SetSpecialOutput();

  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, coutput1);
  mgr->ConnectOutput(task, 2, coutput2);
  return task;
}

//
//____________________________________________________________________________________________
float AliAnalysisTaskKaonXiCorrelation::Eta2y(float pt, float m, float eta) const
{
  return std::asinh(pt / std::hypot(m, pt) * std::sinh(eta));
}

void AliAnalysisTaskKaonXiCorrelation::PostAllData()
{
  PostData(1, fList);
  PostData(2, fTree);
}

Bool_t AliAnalysisTaskKaonXiCorrelation::UserNotify()
{
  TString cfn{CurrentFileName()};
  AliInfo(Form("Setting hash for file %s", cfn.Data()));

  gRandom->SetSeed(cfn.Hash());
  return true;
}

int AliAnalysisTaskKaonXiCorrelation::WhichBDT(double pt)
{
  int iB=0;
  if (fPtBinsBDT.GetSize() > 0)
  {
    while ( (pt < fPtBinsBDT[iB] || pt > fPtBinsBDT[iB+1]) && iB < (fPtBinsBDT.GetSize()-1))
    {
      ++iB;
    }
  }
  return iB;
}

int AliAnalysisTaskKaonXiCorrelation::GetITScls(AliAODTrack *track, int &nSPD, int &nSDD, int &nSSD)
{
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

bool AliAnalysisTaskKaonXiCorrelation::HasTOF(AliAODTrack *track)
{
  bool hasTOFout  = track->GetStatus() & AliVTrack::kTOFout;
  bool hasTOFtime = track->GetStatus() & AliVTrack::kTIME;
  const float len = track->GetIntegratedLength();
  bool hasTOF = hasTOFout && hasTOFtime && (len > 350.);
  return hasTOF;
}

bool AliAnalysisTaskKaonXiCorrelation::HasTwoXiFromSameDaughters(std::vector<XiDaughter> daughters)
{
  for (ULong_t iD = 0; iD < daughters.size(); ++iD)
  {
    for (ULong_t jD = iD + 1; jD < daughters.size(); ++jD)
    {
      if (daughters[iD].IsSame(daughters[jD]))
      {
        return true;
      }
    }
  }
  return false;
}

double AliAnalysisTaskKaonXiCorrelation::GetCustomNsigma(AliAODTrack *t, double cent, int det)
{
  double nsigma = 0.;
  switch (det) {
    case 0:
      nsigma = fPID->NumberOfSigmasITS(t, AliPID::kKaon);
      break;
    case 1:
      nsigma = fPID->NumberOfSigmasTPC(t, AliPID::kKaon);
      break;
    case 2:
      nsigma = fPID->NumberOfSigmasTOF(t, AliPID::kKaon);
      break;
  }
  if (fUseCustomPid && fCustomPidCalib[det]) {
    int binX = fCustomPidCalib[det]->GetXaxis()->FindBin(cent);
    int binY = fCustomPidCalib[det]->GetYaxis()->FindBin(t->Pt());
    int binZ = fCustomPidCalib[det]->GetZaxis()->FindBin(t->Eta());
    double muCalib = fCustomPidCalib[det]->GetBinContent(binX, binY, binZ);
    double sigCalib = fCustomPidCalib[det]->GetBinError(binX, binY, binZ);
    if (muCalib > 1.e-12 || sigCalib > 1.e-12){
      return ( nsigma - muCalib ) / sigCalib;
    }
  }
  return nsigma;
}