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
      AliEasyFemto(),
      fEventCuts{true},
      fPionCuts{128, 0.14, 4.0, 80, 70, 0.8, 0.3, 0.3, 0.83, 3., 3., 3., 0.75},
      fProtonCuts{128, 0.5, 4.0, 80, 70, 0.8, 0.1, 0.1, 0.83, 3., 3., 3., 0.75},
      fDeuteronCuts{256, 0.8, 4.0, 80, 70, 0.8, 0.2, 0.1, 0.83, 3., 3., 3., 1.4},
      fOutputList{nullptr},
      fPID{nullptr},
      fEstimator{1},
      fNormalisationHist{nullptr},
      hPionTrackSelections{nullptr},
      hProtonTrackSelections{nullptr},
      hDeuteronTrackSelections{nullptr},
      hDeltaP{nullptr},
      hPionSpectrum{nullptr},
      hProtonSpectrum{nullptr},
      hDeuteronSpectrum{nullptr},
      hFakeDeuteronSpectrum{nullptr},
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

  hPionTrackSelections = new TH1F("hPionTrackSelections", ";;Entries", kTrackSelection::kNselections, -0.5, (int)kTrackSelection::kNselections - 0.5);
  hProtonTrackSelections = new TH1F("hProtonTrackSelections", ";;Entries", kTrackSelection::kNselections, -0.5, (int)kTrackSelection::kNselections - 0.5);
  hDeuteronTrackSelections = new TH1F("hDeuteronTrackSelections", ";;Entries", kTrackSelection::kNselections, -0.5, (int)kTrackSelection::kNselections - 0.5);
  for (int i = 0; i < kTrackSelection::kNselections; i++)
  {
    hPionTrackSelections->GetXaxis()->SetBinLabel(i + 1, vSelLabels[i]);
    hProtonTrackSelections->GetXaxis()->SetBinLabel(i + 1, vSelLabels[i]);
    hDeuteronTrackSelections->GetXaxis()->SetBinLabel(i + 1, vSelLabels[i]);
  }
  fOutputList->Add(hPionTrackSelections);
  fOutputList->Add(hProtonTrackSelections);
  fOutputList->Add(hDeuteronTrackSelections);

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

    if (applyTrackSelection(track, fPionCuts, fDCAxy, hPionTrackSelections))
    {
      if (TMath::Abs(fDCAxy) < fPionCuts.dcaXY)
      {
        hPionTrackSelections->Fill(AliEasyFemto::kDCAxy);
        if (applyStandardPID(track, fPID, fPionCuts, AliPID::kPion))
        {
          hPionTrackSelections->Fill(AliEasyFemto::kPID);
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
    }

    if (applyTrackSelection(track, fProtonCuts, fDCAxy, hProtonTrackSelections))
    {
      if (TMath::Abs(fDCAxy) < fProtonCuts.dcaXY)
      {
        hProtonTrackSelections->Fill(AliEasyFemto::kDCAxy);
        if (applyStandardPID(track, fPID, fProtonCuts, AliPID::kProton))
        {
          hProtonTrackSelections->Fill(AliEasyFemto::kPID);
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
    }

    if (applyTrackSelection(track, fDeuteronCuts, fDCAxy, hDeuteronTrackSelections))
    {
      if (TMath::Abs(fDCAxy) < fDeuteronCuts.dcaXY)
      {
        hDeuteronTrackSelections->Fill(AliEasyFemto::kDCAxy);
        if (applyStandardPID(track, fPID, fDeuteronCuts, AliPID::kDeuteron))
        {
          hDeuteronTrackSelections->Fill(AliEasyFemto::kPID);
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
  }

  DoFakeCoalescence(v_protonTLV, v_fakedeuteronTLV, v_fakedeuteronID, hFakeDeuteronSpectrum[0], hDeltaP);             // deuteron
  DoFakeCoalescence(v_antiprotonTLV, v_fakeantideuteronTLV, v_fakeantideuteronID, hFakeDeuteronSpectrum[1], hDeltaP); // antideuteron

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

bool AliAnalysisTaskPionDeuteron::applyDeuteronPID(AliAODTrack *track, AliPIDResponse *fPID, CutContainer &cuts, float &tofBeta)
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
    float beta = hasTOF(track, fPID);
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
      return true;
    }
    else
    {
      return false;
    }
  }
}

void AliAnalysisTaskPionDeuteron::Terminate(Option_t *)
{
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
