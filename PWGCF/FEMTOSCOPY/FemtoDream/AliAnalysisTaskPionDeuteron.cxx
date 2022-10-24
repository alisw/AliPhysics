#include "AliAnalysisTaskPionDeuteron.h"
#include "TChain.h"
#include "AliAODEvent.h"
#include "AliAODHeader.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "TMath.h"
#include "AliVTrack.h"
#include "TRandom3.h"

#include <utility>

ClassImp(AliAnalysisTaskPionDeuteron);

AliAnalysisTaskPionDeuteron::AliAnalysisTaskPionDeuteron(const char *name)
    : AliAnalysisTaskSE(name),
      fEventCuts{true},
      fPionCuts{BIT(7), 0.14, 4.0, 80, 70, 0.8, 0.3, 0.3, 3., 3., 3., 0.75},
      fProtonCuts{BIT(7), 0.5, 4.0, 80, 70, 0.8, 0.1, 0.1, 3., 3., 3., 0.75},
      fDeuteronCuts{BIT(8), 0.8, 4.0, 80, 70, 0.8, 0.2, 0.1, 3., 3., 3., 1.4},
      fOutputList{nullptr},
      fPID{nullptr},
      fEstimator{1},
      fP0{0.285},
      fMixingDepth{5},
      fSimpleCoalescence{true},
      fTwoGauss{false},
      fPrimaryPtBins{0},
      fKstarBins{0},
      fNormalisationHist{nullptr},
      hPionSpectrum{nullptr},
      hProtonSpectrum{nullptr},
      hDeuteronSpectrum{nullptr},
      hFakeDeuteronSpectrum{nullptr},
      hDeltaP{nullptr},
      hNparticles{nullptr},
      hNparticlesFake{nullptr},
      hSameEventPionProtonKstarLS{nullptr},
      hSameEventPionProtonKstarUS{nullptr},
      hMixedEventPionProtonKstarLS{nullptr},
      hMixedEventPionProtonKstarUS{nullptr},
      hSameEventPionDeuteronKstarLS{nullptr},
      hSameEventPionDeuteronKstarUS{nullptr},
      hMixedEventPionDeuteronKstarLS{nullptr},
      hMixedEventPionDeuteronKstarUS{nullptr},
      hSameEventPionFakeDeuteronKstarLS{nullptr},
      hSameEventPionFakeDeuteronKstarUS{nullptr},
      hMixedEventPionFakeDeuteronKstarLS{nullptr},
      hMixedEventPionFakeDeuteronKstarUS{nullptr},
      fMixingBufferProtons(),
      fMixingBufferAntiprotons(),
      fMixingBufferDeuterons(),
      fMixingBufferAntideuterons(),
      fMixingBufferFakeDeuterons(),
      fMixingBufferFakeAntideuterons()
{
  fZvtxArray = {-10., -8., -6., -4., -2., 0., 2., 4., 6., 8., 10.};
  fNZvtxBins = (int)fZvtxArray.size();
  fMultiplicityArray = {0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 80, 84, 88, 92, 96, 100};
  fNMultiplicityBins = (int)fMultiplicityArray.size();
  int BufferSize = fNZvtxBins * fNMultiplicityBins;
  fMixingBufferProtons.reserve(BufferSize);
  fMixingBufferAntiprotons.reserve(BufferSize);
  fMixingBufferDeuterons.reserve(BufferSize);
  fMixingBufferAntideuterons.reserve(BufferSize);
  fMixingBufferFakeDeuterons.reserve(BufferSize);
  fMixingBufferFakeAntideuterons.reserve(BufferSize);
  for (int i = 0; i < BufferSize; i++)
  {
    fMixingBufferProtons.push_back(CustomQueue<std::vector<TLorentzVector>>(fMixingDepth));
    fMixingBufferAntiprotons.push_back(CustomQueue<std::vector<TLorentzVector>>(fMixingDepth));
    fMixingBufferDeuterons.push_back(CustomQueue<std::vector<TLorentzVector>>(fMixingDepth));
    fMixingBufferAntideuterons.push_back(CustomQueue<std::vector<TLorentzVector>>(fMixingDepth));
    fMixingBufferFakeDeuterons.push_back(CustomQueue<std::vector<TLorentzVector>>(fMixingDepth));
    fMixingBufferFakeAntideuterons.push_back(CustomQueue<std::vector<TLorentzVector>>(fMixingDepth));
  }
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskPionDeuteron::~AliAnalysisTaskPionDeuteron()
{
  if (fOutputList)
  {
    delete fOutputList;
  }
}

void AliAnalysisTaskPionDeuteron::UserCreateOutputObjects()
{

  // Create Output List
  fOutputList = new TList();
  fOutputList->SetOwner(true);

  std::array<std::string, 4> norm_labels = {
      "No cuts",
      "Event selection",
      "Vertex reconstruction and quality",
      "Vertex position"};
  fNormalisationHist = new TH1F("fNormalisationHist", ";;Entries", norm_labels.size(), -.5, norm_labels.size() - 0.5);
  for (size_t iB = 1; iB <= norm_labels.size(); iB++)
  {
    fNormalisationHist->GetXaxis()->SetBinLabel(iB, norm_labels[iB - 1].data());
  }
  fOutputList->Add(fNormalisationHist);

  const int nPrimaryPtBins = fPrimaryPtBins.GetSize() - 1;
  const float *primaryPtBins = fPrimaryPtBins.GetArray();
  const int nKstarBins = fKstarBins.GetSize() - 1;
  const float *kStarBins = fKstarBins.GetArray();

  const char *charge_label[2] = {"Pos", "Neg"};

  for (int iCharge = 0; iCharge < 2; iCharge++)
  {
    hPionSpectrum[iCharge] = new TH1F(Form("h%sPionSpectrum", charge_label[iCharge]), ";#it{p} (GeV/#it{c});Entries", nPrimaryPtBins, primaryPtBins);
    fOutputList->Add(hPionSpectrum[iCharge]);
    hProtonSpectrum[iCharge] = new TH1F(Form("h%sProtonSpectrum", charge_label[iCharge]), ";#it{p} (GeV/#it{c});Entries", nPrimaryPtBins, primaryPtBins);
    fOutputList->Add(hProtonSpectrum[iCharge]);
    hDeuteronSpectrum[iCharge] = new TH1F(Form("h%sDeuteronSpectrum", charge_label[iCharge]), ";#it{p} (GeV/#it{c});Entries", nPrimaryPtBins, primaryPtBins);
    fOutputList->Add(hDeuteronSpectrum[iCharge]);
    hFakeDeuteronSpectrum[iCharge] = new TH1F(Form("h%sFakeDeuteronSpectrum", charge_label[iCharge]), ";#it{p} (GeV/#it{c});Entries", nPrimaryPtBins, primaryPtBins);
    fOutputList->Add(hFakeDeuteronSpectrum[iCharge]);
    hNparticles[iCharge] = new TH2F(Form("h%sNparticles", charge_label[iCharge]), ";N_{#pi};N_{d}", 51, -0.5, 50.5, 4, -0.5, 3.5);
    fOutputList->Add(hNparticles[iCharge]);
    hNparticlesFake[iCharge] = new TH2F(Form("h%sNparticlesFake", charge_label[iCharge]), ";N_{#pi};N_{d (fake)}", 51, -0.5, 50.5, 4, -0.5, 3.5);
    fOutputList->Add(hNparticlesFake[iCharge]);
    hSameEventPionProtonKstarLS[iCharge] = new TH1F(Form("h%sSameEventPionProtonKstarLS", charge_label[iCharge]), ";#it{k}* (GeV/#it{c});Entries", nKstarBins, kStarBins);
    fOutputList->Add(hSameEventPionProtonKstarLS[iCharge]);
    hSameEventPionProtonKstarUS[iCharge] = new TH1F(Form("h%sSameEventPionProtonKstarUS", charge_label[iCharge]), ";#it{k}* (GeV/#it{c});Entries", nKstarBins, kStarBins);
    fOutputList->Add(hSameEventPionProtonKstarUS[iCharge]);
    hMixedEventPionProtonKstarLS[iCharge] = new TH1F(Form("h%sMixedEventPionProtonKstarLS", charge_label[iCharge]), ";#it{k}* (GeV/#it{c});Entries", nKstarBins, kStarBins);
    fOutputList->Add(hMixedEventPionProtonKstarLS[iCharge]);
    hMixedEventPionProtonKstarUS[iCharge] = new TH1F(Form("h%sMixedEventPionProtonKstarUS", charge_label[iCharge]), ";#it{k}* (GeV/#it{c});Entries", nKstarBins, kStarBins);
    fOutputList->Add(hMixedEventPionProtonKstarUS[iCharge]);
    hSameEventPionDeuteronKstarLS[iCharge] = new TH1F(Form("h%sSameEventPionDeuteronKstarLS", charge_label[iCharge]), ";#it{k}* (GeV/#it{c});Entries", nKstarBins, kStarBins);
    fOutputList->Add(hSameEventPionDeuteronKstarLS[iCharge]);
    hSameEventPionDeuteronKstarUS[iCharge] = new TH1F(Form("h%sSameEventPionDeuteronKstarUS", charge_label[iCharge]), ";#it{k}* (GeV/#it{c});Entries", nKstarBins, kStarBins);
    fOutputList->Add(hSameEventPionDeuteronKstarUS[iCharge]);
    hMixedEventPionDeuteronKstarLS[iCharge] = new TH1F(Form("h%sMixedEventPionDeuteronKstarLS", charge_label[iCharge]), ";#it{k}* (GeV/#it{c});Entries", nKstarBins, kStarBins);
    fOutputList->Add(hMixedEventPionDeuteronKstarLS[iCharge]);
    hMixedEventPionDeuteronKstarUS[iCharge] = new TH1F(Form("h%sMixedEventPionDeuteronKstarUS", charge_label[iCharge]), ";#it{k}* (GeV/#it{c});Entries", nKstarBins, kStarBins);
    fOutputList->Add(hMixedEventPionFakeDeuteronKstarUS[iCharge]);
    hSameEventPionFakeDeuteronKstarLS[iCharge] = new TH1F(Form("h%sSameEventPionFakeDeuteronKstarLS", charge_label[iCharge]), ";#it{k}* (GeV/#it{c});Entries", nKstarBins, kStarBins);
    fOutputList->Add(hSameEventPionFakeDeuteronKstarLS[iCharge]);
    hSameEventPionFakeDeuteronKstarUS[iCharge] = new TH1F(Form("h%sSameEventPionFakeDeuteronKstarUS", charge_label[iCharge]), ";#it{k}* (GeV/#it{c});Entries", nKstarBins, kStarBins);
    fOutputList->Add(hSameEventPionFakeDeuteronKstarUS[iCharge]);
    hMixedEventPionFakeDeuteronKstarLS[iCharge] = new TH1F(Form("h%sMixedEventPionFakeDeuteronKstarLS", charge_label[iCharge]), ";#it{k}* (GeV/#it{c});Entries", nKstarBins, kStarBins);
    fOutputList->Add(hMixedEventPionFakeDeuteronKstarLS[iCharge]);
    hMixedEventPionFakeDeuteronKstarUS[iCharge] = new TH1F(Form("h%sMixedEventPionFakeDeuteronKstarUS", charge_label[iCharge]), ";#it{k}* (GeV/#it{c});Entries", nKstarBins, kStarBins);
    fOutputList->Add(hMixedEventPionFakeDeuteronKstarUS[iCharge]);
  }

  // DeltaP Distribution
  hDeltaP = new TH1F("hDeltaP", ";#Delta p (GeV/#it{c});Entries", 200, 0, 1.);
  hDeltaP->Sumw2();
  fOutputList->Add(hDeltaP);

  PostData(1, fOutputList);
}

void AliAnalysisTaskPionDeuteron::UserExec(Option_t *)
{

  AliAODEvent *ev = (AliAODEvent *)InputEvent();

  bool eventAccepted = fEventCuts.AcceptEvent(ev);

  std::array<AliEventCuts::NormMask, 4> norm_masks{
      AliEventCuts::kAnyEvent,
      AliEventCuts::kPassesNonVertexRelatedSelections,
      AliEventCuts::kHasReconstructedVertex,
      AliEventCuts::kPassesAllCuts};
  for (int iC = 0; iC < 4; ++iC)
  {
    if (fEventCuts.CheckNormalisationMask(norm_masks[iC]))
    {
      fNormalisationHist->Fill(iC);
    }
  }

  if (!eventAccepted)
  {
    PostData(1, fOutputList);
  }

  AliAODHeader *header = (AliAODHeader *)ev->GetHeader();
  if (!header)
  {
    ::Fatal("AliAnalysisTaskPionDeuteronMC::UserExec", "Header not found.");
  }

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *handl = (AliInputEventHandler *)mgr->GetInputEventHandler();
  fPID = handl->GetPIDResponse();

  std::vector<TLorentzVector> v_pospionTLV, v_negpionTLV, v_protonTLV, v_antiprotonTLV, v_deuteronTLV, v_antideuteronTLV, v_fakedeuteronTLV, v_fakeantideuteronTLV;
  std::vector<int> v_pospionID, v_negpionID, v_protonID, v_antiprotonID, v_deuteronID, v_antideuteronID;
  std::vector<std::pair<int, int>> v_fakedeuteronID, v_fakeantideuteronID;

  for (int iT = 0; iT < (int)ev->GetNumberOfTracks(); iT++)
  {
    AliAODTrack *track = dynamic_cast<AliAODTrack *>(ev->GetTrack(iT));

    TLorentzVector fourVector;
    float fDCAxy = 0.;

    if (applyTrackSelection(track, fPionCuts, fDCAxy) && applyStandardPID(track, fPionCuts, AliPID::kPion))
    {
      if (TMath::Abs(fDCAxy) < fPionCuts.dcaXY)
      {
        fourVector.SetPtEtaPhiM(track->Pt(), track->Eta(), track->Phi(), AliPID::ParticleMass(AliPID::kPion));
        if (track->Charge() > 0)
        {
          v_pospionID.push_back(iT);
          v_pospionTLV.push_back(fourVector);
          hPionSpectrum[0]->Fill(track->P());
        }
        else
        {
          v_negpionID.push_back(iT);
          v_negpionTLV.push_back(fourVector);
          hPionSpectrum[1]->Fill(track->P());
        }
      }
    }

    if (applyTrackSelection(track, fProtonCuts, fDCAxy) && applyStandardPID(track, fProtonCuts, AliPID::kProton))
    {
      if (TMath::Abs(fDCAxy) < fProtonCuts.dcaXY)
      {
        fourVector.SetPtEtaPhiM(track->Pt(), track->Eta(), track->Phi(), AliPID::ParticleMass(AliPID::kProton));
        if (track->Charge() > 0)
        {
          v_protonID.push_back(iT);
          v_protonTLV.push_back(fourVector);
          hProtonSpectrum[0]->Fill(track->P());
        }
        else
        {
          v_antiprotonID.push_back(iT);
          v_antiprotonTLV.push_back(fourVector);
          hProtonSpectrum[1]->Fill(track->P());
        }
      }
    }

    if (applyTrackSelection(track, fDeuteronCuts, fDCAxy) && applyStandardPID(track, fDeuteronCuts, AliPID::kDeuteron))
    {
      if (TMath::Abs(fDCAxy) < fDeuteronCuts.dcaXY)
      {
        fourVector.SetPtEtaPhiM(track->Pt(), track->Eta(), track->Phi(), AliPID::ParticleMass(AliPID::kDeuteron));
        if (track->Charge() > 0)
        {
          v_deuteronID.push_back(iT);
          v_deuteronTLV.push_back(fourVector);
          hDeuteronSpectrum[0]->Fill(track->P());
        }
        else
        {
          v_antideuteronID.push_back(iT);
          v_antideuteronTLV.push_back(fourVector);
          hDeuteronSpectrum[1]->Fill(track->P());
        }
      }
    }
  }

  DoFakeCoalescence(v_protonTLV, v_fakedeuteronTLV, v_fakedeuteronID, hFakeDeuteronSpectrum[0]);             // deuteron
  DoFakeCoalescence(v_antiprotonTLV, v_fakeantideuteronTLV, v_fakeantideuteronID, hFakeDeuteronSpectrum[1]); // antideuteron

  // count particles
  float nPosPions = (float)v_pospionTLV.size();
  float nNegPions = (float)v_negpionTLV.size();

  float nDeuterons = (float)v_deuteronTLV.size();
  float nAntideuterons = (float)v_antideuteronTLV.size();

  float nFakeDeuterons = (float)v_fakedeuteronTLV.size();
  float nFakeAntideuterons = (float)v_fakeantideuteronTLV.size();

  hNparticles[0]->Fill(nPosPions, nDeuterons);
  hNparticles[1]->Fill(nNegPions, nAntideuterons);

  hNparticlesFake[0]->Fill(nPosPions, nFakeDeuterons);
  hNparticlesFake[1]->Fill(nNegPions, nFakeAntideuterons);

  // same event
  for (int iPosPion = 0; iPosPion < (int)v_pospionTLV.size(); iPosPion++)
  {
    int pionID = v_pospionID[iPosPion];
    for (int iPosProton = 0; iPosProton < (int)v_protonTLV.size(); iPosProton++)
    {
      if (v_protonID[iPosProton] == pionID)
        continue;
      else
        hSameEventPionProtonKstarLS[0]->Fill(GetKstar(v_pospionTLV[iPosPion], v_protonTLV[iPosProton]));
    }
    for (int iNegProton = 0; iNegProton < (int)v_antiprotonTLV.size(); iNegProton++)
    {
      if (v_antiprotonID[iNegProton] == pionID)
        continue;
      else
        hSameEventPionProtonKstarUS[0]->Fill(GetKstar(v_pospionTLV[iPosPion], v_antiprotonTLV[iNegProton]));
    }
    for (int iPosDeuteron = 0; iPosDeuteron < (int)v_deuteronTLV.size(); iPosDeuteron++)
    {
      if (v_deuteronID[iPosDeuteron] == pionID)
        continue;
      else
        hSameEventPionDeuteronKstarLS[0]->Fill(GetKstar(v_pospionTLV[iPosPion], v_deuteronTLV[iPosDeuteron]));
    }
    for (int iNegDeuteron = 0; iNegDeuteron < (int)v_antideuteronTLV.size(); iNegDeuteron++)
    {
      if (v_deuteronID[iNegDeuteron] == pionID)
        continue;
      else
        hSameEventPionDeuteronKstarUS[0]->Fill(GetKstar(v_pospionTLV[iPosPion], v_antideuteronTLV[iNegDeuteron]));
    }
    for (int iPosFakeDeuteron = 0; iPosFakeDeuteron < (int)v_deuteronTLV.size(); iPosFakeDeuteron++)
    {
      if (v_fakedeuteronID[iPosFakeDeuteron].first == pionID || v_fakedeuteronID[iPosFakeDeuteron].second == pionID)
        continue;
      else
        hSameEventPionFakeDeuteronKstarLS[0]->Fill(GetKstar(v_pospionTLV[iPosPion], v_fakedeuteronTLV[iPosFakeDeuteron]));
    }
    for (int iNegFakeDeuteron = 0; iNegFakeDeuteron < (int)v_fakeantideuteronTLV.size(); iNegFakeDeuteron++)
    {
      if (v_fakeantideuteronID[iNegFakeDeuteron].first == pionID || v_fakeantideuteronID[iNegFakeDeuteron].second == pionID)
        continue;
      else
        hSameEventPionFakeDeuteronKstarUS[0]->Fill(GetKstar(v_pospionTLV[iPosPion], v_fakeantideuteronTLV[iNegFakeDeuteron]));
    }
  }

  for (int iNegPion = 0; iNegPion < (int)v_negpionTLV.size(); iNegPion++)
  {
    int pionID = v_negpionID[iNegPion];
    for (int iPosProton = 0; iPosProton < (int)v_protonTLV.size(); iPosProton++)
    {
      if (v_protonID[iPosProton] == pionID)
        continue;
      else
        hSameEventPionProtonKstarUS[1]->Fill(GetKstar(v_negpionTLV[iNegPion], v_protonTLV[iPosProton]));
    }
    for (int iNegProton = 0; iNegProton < (int)v_antiprotonTLV.size(); iNegProton++)
    {
      if (v_antiprotonID[iNegProton] == pionID)
        continue;
      else
        hSameEventPionProtonKstarLS[1]->Fill(GetKstar(v_negpionTLV[iNegPion], v_antiprotonTLV[iNegProton]));
    }
    for (int iPosDeuteron = 0; iPosDeuteron < (int)v_deuteronTLV.size(); iPosDeuteron++)
    {
      if (v_deuteronID[iPosDeuteron] == pionID)
        continue;
      else
        hSameEventPionDeuteronKstarUS[1]->Fill(GetKstar(v_negpionTLV[iNegPion], v_deuteronTLV[iPosDeuteron]));
    }
    for (int iNegDeuteron = 0; iNegDeuteron < (int)v_antideuteronTLV.size(); iNegDeuteron++)
    {
      if (v_antideuteronID[iNegDeuteron] == pionID)
        continue;
      else
        hSameEventPionDeuteronKstarLS[1]->Fill(GetKstar(v_negpionTLV[iNegPion], v_antideuteronTLV[iNegDeuteron]));
    }
    for (int iPosFakeDeuteron = 0; iPosFakeDeuteron < (int)v_deuteronTLV.size(); iPosFakeDeuteron++)
    {
      if (v_fakedeuteronID[iPosFakeDeuteron].first == pionID || v_fakedeuteronID[iPosFakeDeuteron].second == pionID)
        continue;
      else
        hSameEventPionFakeDeuteronKstarUS[1]->Fill(GetKstar(v_negpionTLV[iNegPion], v_fakedeuteronTLV[iPosFakeDeuteron]));
    }
    for (int iNegFakeDeuteron = 0; iNegFakeDeuteron < (int)v_fakeantideuteronTLV.size(); iNegFakeDeuteron++)
    {
      if (v_fakeantideuteronID[iNegFakeDeuteron].first == pionID || v_fakeantideuteronID[iNegFakeDeuteron].second == pionID)
        continue;
      else
        hSameEventPionDeuteronKstarLS[1]->Fill(GetKstar(v_negpionTLV[iNegPion], v_fakeantideuteronTLV[iNegFakeDeuteron]));
    }
  }

  // event mixing
  float zvtx = ev->GetPrimaryVertex()->GetZ();
  float mult = header->GetRefMultiplicityComb08();
  int index = FindBin(zvtx, mult);

  if (index >= 0)
  {
    // match pions with particles in the buffer
    if (v_pospionTLV.size() > 0)
    {
      FillMixedEvent(v_pospionTLV, fMixingBufferProtons[index], hMixedEventPionProtonKstarLS[0]);
      FillMixedEvent(v_pospionTLV, fMixingBufferAntiprotons[index], hMixedEventPionProtonKstarUS[0]);
      FillMixedEvent(v_pospionTLV, fMixingBufferDeuterons[index], hMixedEventPionDeuteronKstarLS[0]);
      FillMixedEvent(v_pospionTLV, fMixingBufferAntideuterons[index], hMixedEventPionDeuteronKstarUS[0]);
      FillMixedEvent(v_pospionTLV, fMixingBufferFakeDeuterons[index], hMixedEventPionFakeDeuteronKstarLS[0]);
      FillMixedEvent(v_pospionTLV, fMixingBufferFakeAntideuterons[index], hMixedEventPionFakeDeuteronKstarUS[0]);
    }
    if (v_negpionTLV.size() > 0)
    {
      FillMixedEvent(v_negpionTLV, fMixingBufferProtons[index], hMixedEventPionProtonKstarUS[1]);
      FillMixedEvent(v_negpionTLV, fMixingBufferAntiprotons[index], hMixedEventPionProtonKstarLS[1]);
      FillMixedEvent(v_negpionTLV, fMixingBufferDeuterons[index], hMixedEventPionDeuteronKstarUS[1]);
      FillMixedEvent(v_negpionTLV, fMixingBufferAntideuterons[index], hMixedEventPionDeuteronKstarLS[1]);
      FillMixedEvent(v_negpionTLV, fMixingBufferFakeDeuterons[index], hMixedEventPionFakeDeuteronKstarUS[1]);
      FillMixedEvent(v_negpionTLV, fMixingBufferFakeAntideuterons[index], hMixedEventPionFakeDeuteronKstarLS[1]);
    }
    // update particle buffer
    if (v_protonTLV.size() > 0)
    {
      fMixingBufferProtons[index].Fill(v_protonTLV);
    }
    if (v_antiprotonTLV.size() > 0)
    {
      fMixingBufferAntiprotons[index].Fill(v_antiprotonTLV);
    }
    if (v_deuteronTLV.size() > 0)
    {
      fMixingBufferDeuterons[index].Fill(v_deuteronTLV);
    }
    if (v_antideuteronTLV.size() > 0)
    {
      fMixingBufferAntideuterons[index].Fill(v_antideuteronTLV);
    }
    if (v_fakedeuteronTLV.size() > 0)
    {
      fMixingBufferFakeDeuterons[index].Fill(v_fakedeuteronTLV);
    }
    if (v_fakeantideuteronTLV.size() > 0)
    {
      fMixingBufferFakeAntideuterons[index].Fill(v_fakeantideuteronTLV);
    }
  }

  PostData(1, fOutputList);
}

bool AliAnalysisTaskPionDeuteron::applyTrackSelection(AliAODTrack *track, CutContainer &cuts, float &dcaXY)
{
  if (!track->TestFilterBit(cuts.filterBit))
    return false;
  if (TMath::Abs(track->Eta()) > cuts.eta)
    return false;
  if (track->Pt() < cuts.minPt)
    return false;
  if (track->Pt() > cuts.maxPt)
    return false;
  if (track->GetTPCNcls() < cuts.nTPCcls)
    return false;
  if (track->GetTPCnclsS() > 0)
    return false;
  if (track->GetTPCNCrossedRows() < (int)cuts.nCrossedRows)
    return false;
  if (track->GetTPCCrossedRows() / (float)track->GetTPCNclsF() < cuts.crossedRowsOverFindable)
    return false;
  if (track->GetTPCCrossedRows() / (float)track->GetTPCNclsF() < cuts.crossedRowsOverFindable)
    return false;
  float dcaZ = 0.;
  track->GetImpactParameters(dcaXY, dcaZ);
  if (TMath::Abs(dcaZ) > cuts.dcaZ)
    return false;
  return true;
}

bool AliAnalysisTaskPionDeuteron::applyStandardPID(AliAODTrack *track, CutContainer &cuts, AliPID::EParticleType kType)
{
  float nSigmaTPC = fPID->NumberOfSigmasTPC(track, kType);

  if (track->Pt() < cuts.ptPIDthreshold)
  {
    if (TMath::Abs(nSigmaTPC) > cuts.nSigmaTPC)
      return false;
    else
      return true;
  }
  else
  {
    float beta = hasTOF(track);
    if (beta > 0)
    {
      float nSigmaTOF = fPID->NumberOfSigmasTOF(track, kType);
      float nSigmaComb = TMath::Sqrt(nSigmaTPC * nSigmaTPC + nSigmaTOF * nSigmaTOF);
      if (TMath::Abs(nSigmaComb) > cuts.nSigmaComb)
        return false;
      else
        return true;
    }
    else
    {
      return false;
    }
  }
}

bool AliAnalysisTaskPionDeuteron::applyDeuteronPID(AliAODTrack *track, CutContainer &cuts, float &tofBeta)
{
  float nSigmaTPC = fPID->NumberOfSigmasTPC(track, AliPID::kDeuteron);
  if (track->Pt() < cuts.ptPIDthreshold)
  {
    if (TMath::Abs(nSigmaTPC) <= cuts.nSigmaTPC && fPID->NumberOfSigmasTPC(track, AliPID::kElectron) > 3.)
      return true;
    else
      return false;
  }
  else
  {
    float beta = hasTOF(track);
    if (beta > 0)
    {
      if (TMath::Abs(nSigmaTPC) > cuts.nSigmaTPC)
        return false;
      if (TMath::Abs(fPID->NumberOfSigmasTPC(track, AliPID::kPion)) < 3.)
        return false;
      float nSigmaTOF = fPID->NumberOfSigmasTOF(track, AliPID::kDeuteron);
      if (nSigmaTOF < -3.)
        return false;
      if (nSigmaTOF > 5.)
        return false;
      if (TMath::Abs(fPID->NumberOfSigmasTOF(track, AliPID::kElectron)) < 5.)
        return false;
      if (TMath::Abs(fPID->NumberOfSigmasTOF(track, AliPID::kPion)) < 5.)
        return false;
    }
    else
    {
      return false;
    }
  }
}

float AliAnalysisTaskPionDeuteron::hasTOF(AliAODTrack *track)
{
  bool hasTOFout = track->GetStatus() & AliVTrack::kTOFout;
  bool hasTOFtime = track->GetStatus() & AliVTrack::kTIME;
  const float len = track->GetIntegratedLength();
  bool hasTOF = hasTOFout && hasTOFtime && (len > 350.);

  if (!hasTOF)
    return -1.;
  const float time = track->GetTOFsignal() - fPID->GetTOFResponse().GetStartTime(track->GetTPCmomentum());
  const float beta = len / (time * LIGHT_SPEED);
  return beta;
}

void AliAnalysisTaskPionDeuteron::Terminate(Option_t *)
{
}

float AliAnalysisTaskPionDeuteron::GetKstar(TLorentzVector &p1, TLorentzVector &p2)
{
  TLorentzVector sum = p1 + p2;
  TVector3 boost_vector = sum.BoostVector();
  TLorentzVector p1_prf = p1;
  p1_prf.Boost(-boost_vector);
  TLorentzVector p2_prf = p2;
  p2_prf.Boost(-boost_vector);
  TLorentzVector kStar = p1_prf - p2_prf;
  return 0.5 * kStar.P();
}

void AliAnalysisTaskPionDeuteron::SetPrimaryPtBins(int nbins, float min, float max)
{
  const float delta = (max - min) / nbins;
  fPrimaryPtBins.Set(nbins + 1);
  for (int iB = 0; iB < nbins; ++iB)
  {
    fPrimaryPtBins[iB] = min + iB * delta;
  }
  fPrimaryPtBins[nbins] = max;
}

void AliAnalysisTaskPionDeuteron::SetPrimaryPtBins(int nbins, float *bins)
{
  fPrimaryPtBins.Set(nbins + 1, bins);
}

void AliAnalysisTaskPionDeuteron::SetKstarBins(int nbins, float min, float max)
{
  const float delta = (max - min) / nbins;
  fKstarBins.Set(nbins + 1);
  for (int iB = 0; iB < nbins; ++iB)
  {
    fKstarBins[iB] = min + iB * delta;
  }
  fKstarBins[nbins] = max;
}

void AliAnalysisTaskPionDeuteron::SetKstarBins(int nbins, float *bins)
{
  fKstarBins.Set(nbins + 1, bins);
}

void AliAnalysisTaskPionDeuteron::SetZvtxArray(std::vector<float> &vec)
{
  fZvtxArray.clear();
  for (auto p : vec)
  {
    fZvtxArray.push_back(p);
  }
  fNZvtxBins = (int)fZvtxArray.size();
}

void AliAnalysisTaskPionDeuteron::SetMultiplicityArray(std::vector<float> &vec)
{
  fMultiplicityArray.clear();
  for (auto p : vec)
  {
    fMultiplicityArray.push_back(p);
  }
  fNMultiplicityBins = (int)fMultiplicityArray.size();
}

void AliAnalysisTaskPionDeuteron::SetMixingDepth(unsigned int depth)
{
  for (auto p : fMixingBufferDeuterons)
  {
    p.SetDepth(depth);
  }
  for (auto p : fMixingBufferAntideuterons)
  {
    p.SetDepth(depth);
  }
  for (auto p : fMixingBufferFakeDeuterons)
  {
    p.SetDepth(depth);
  }
  for (auto p : fMixingBufferFakeAntideuterons)
  {
    p.SetDepth(depth);
  }
  for (auto p : fMixingBufferProtons)
  {
    p.SetDepth(depth);
  }
  for (auto p : fMixingBufferAntiprotons)
  {
    p.SetDepth(depth);
  }
}

int AliAnalysisTaskPionDeuteron::FindBin(float zvtx, float mult)
{
  if (zvtx < fZvtxArray[0] || zvtx > fZvtxArray[fNZvtxBins - 1] || zvtx < fMultiplicityArray[0] || zvtx > fMultiplicityArray[fNMultiplicityBins - 1])
    return -1;

  auto iter_zvtx = std::upper_bound(fZvtxArray.begin(), fZvtxArray.end(), zvtx,
                                    [](const float &comp1, const float &comp2)
                                    { return comp1 < comp2; });
  int index_zvtx = std::distance(fZvtxArray.begin(), iter_zvtx);
  auto iter_mult = std::upper_bound(fMultiplicityArray.begin(), fMultiplicityArray.end(), mult,
                                    [](const float &comp1, const float &comp2)
                                    { return comp1 < comp2; });
  int index_mult = std::distance(fMultiplicityArray.begin(), iter_mult);
  return index_zvtx * fNZvtxBins + index_mult;
}

void AliAnalysisTaskPionDeuteron::FillMixedEvent(std::vector<TLorentzVector> &vec, CustomQueue<std::vector<TLorentzVector>> &buffer, TH1F *histo)
{
  for (int i = 0; i < buffer.GetSize(); i++)
  {
    auto vectorBuffer = buffer.GetElement(i);
    for (auto &elementA : vec)
    {
      for (auto &elementB : vectorBuffer)
      {
        histo->Fill(GetKstar(elementA, elementB));
      }
    }
  }
}

void AliAnalysisTaskPionDeuteron::DoFakeCoalescence(std::vector<TLorentzVector> &v_proton, std::vector<TLorentzVector> &v_fakedeuteron, std::vector<std::pair<int, int>> &v_fakedeuteronID, TH1F *histo)
{

  int n_protons = (int)v_proton.size();
  std::vector<int> v_proton_status(n_protons, 0);

  for (int i = 0; i < (int)v_proton.size(); i++)
  {
    if (v_proton_status[i] != 0)
      continue;
    for (int j = i + 1; j < (int)v_proton.size(); j++)
    {
      if (v_proton_status[j] != 0)
        continue;
      TLorentzVector deuteron;
      deuteron.SetXYZM(v_proton[i].Px() + v_proton[j].Px(), v_proton[i].Py() + v_proton[j].Px(), v_proton[i].Pz() + v_proton[j].Pz(), AliPID::ParticleMass(AliPID::kDeuteron));
      TVector3 boost_vector = deuteron.BoostVector();

      // Lorentz Transformations (from Lab to Deuteron Frame)
      TLorentzVector proton_prime = v_proton[i];
      proton_prime.Boost(-boost_vector);
      TLorentzVector pseudoneutron_prime = v_proton[j];
      pseudoneutron_prime.Boost(-boost_vector);
      TLorentzVector deltaPvector = proton_prime - pseudoneutron_prime;
      double deltaP = deltaPvector.P();
      double mt = 0.5 * (v_proton[i].Mt() + v_proton[j].Mt());

      // Fill DeltaP Distribution
      hDeltaP->Fill(deltaP);

      if (fSimpleCoalescence)
      {
        // Simple Coalescence Condition
        if (deltaP < fP0)
        {
          v_fakedeuteron.push_back(deuteron);
          v_fakedeuteronID.push_back(std::make_pair(i, j));
          v_proton_status[j] = 1;
          histo->Fill(deuteron.P());
          break;
        }
      }
      else
      {
        float value = (float)gRandom->Uniform();
        float q = 0.5 * deltaP;
        float radius = ComputeRadius(mt);
        float prob = 1;
        if (fTwoGauss)
        {
          prob = 3 * (TMath::Power(3.979, 3) / TMath::Power(3.979 * 3.979 + 4 * radius * radius, 1.5) * TMath::Exp(-1 * q * q * 3.979 * 3.979 * 5.068 * 5.068) + (1. - 0.581) * TMath::Power(0.890, 3) / TMath::Power(0.890 * 0.890 + 4 * radius * radius, 1.5) * TMath::Exp(-1 * q * q * 0.890 * 0.890 * 5.068 * 5.068));
        }
        else
        {
          prob = 3 * TMath::Power(3.2, 3) / TMath::Power(3.2 * 3.2 + radius * radius, 1.5) * TMath::Exp(-1 * q * q * 3.2 * 3.2 * 5.068 * 5.068);
        }
        if (value < prob)
        {
          v_fakedeuteron.push_back(deuteron);
          v_fakedeuteronID.push_back(std::make_pair(i, j));
          v_proton_status[j] = 1;
          histo->Fill(deuteron.P());
          break;
        }
      }
    }
  }
}

float AliAnalysisTaskPionDeuteron::ComputeRadius(float mt)
{
  return 0.8 + TMath::Exp(1.4 - 1.9 * mt);
}
