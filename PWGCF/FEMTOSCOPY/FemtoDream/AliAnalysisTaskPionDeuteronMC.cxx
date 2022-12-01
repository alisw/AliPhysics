#include "AliAnalysisTaskPionDeuteronMC.h"
#include "AliAODEvent.h"
#include "AliAODMCHeader.h"
#include "TClonesArray.h"
#include "AliAODMCParticle.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TDatabasePDG.h"
#include "AliPID.h"
#include "TChain.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TRandom3.h"

#include <vector>
#include <iostream>
#include <algorithm>
#include <array>

using namespace std;

///\cond CLASSIMP
ClassImp(AliAnalysisTaskPionDeuteronMC);
///\endcond

AliAnalysisTaskPionDeuteronMC::AliAnalysisTaskPionDeuteronMC(const char *name)
    : AliAnalysisTaskSE(name),
      AliEasyFemto(),
      fEventCuts{true},
      fOutputList{nullptr},
      fNormalisationHist{nullptr},
      hPionSpectrum{nullptr},
      hProtonSpectrum{nullptr},
      hNeutronSpectrum{nullptr},
      hDeuteronSpectrum{nullptr},
      hDeltaP{nullptr},
      hNparticles{nullptr},
      hSameEventPionProtonKstarLS{nullptr},
      hSameEventPionProtonKstarUS{nullptr},
      hMixedEventPionProtonKstarLS{nullptr},
      hMixedEventPionProtonKstarUS{nullptr},
      hSameEventPionDeuteronKstarLS{nullptr},
      hSameEventPionDeuteronKstarUS{nullptr},
      hMixedEventPionDeuteronKstarLS{nullptr},
      hMixedEventPionDeuteronKstarUS{nullptr},
      fMixingBufferProtons(),
      fMixingBufferAntiprotons(),
      fMixingBufferDeuterons(),
      fMixingBufferAntideuterons()
{
  int BufferSize = fNZvtxBins * fNMultiplicityBins;
  fMixingBufferProtons.reserve(BufferSize);
  fMixingBufferAntiprotons.reserve(BufferSize);
  fMixingBufferDeuterons.reserve(BufferSize);
  fMixingBufferAntideuterons.reserve(BufferSize);
  for (int i = 0; i < BufferSize; i++)
  {
    fMixingBufferProtons.push_back(CustomQueue<std::vector<TLorentzVector>>(fMixingDepth));
    fMixingBufferAntiprotons.push_back(CustomQueue<std::vector<TLorentzVector>>(fMixingDepth));
    fMixingBufferDeuterons.push_back(CustomQueue<std::vector<TLorentzVector>>(fMixingDepth));
    fMixingBufferAntideuterons.push_back(CustomQueue<std::vector<TLorentzVector>>(fMixingDepth));
  }
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskPionDeuteronMC::~AliAnalysisTaskPionDeuteronMC()
{
  if (fOutputList)
  {
    delete fOutputList;
  }
}

void AliAnalysisTaskPionDeuteronMC::UserCreateOutputObjects()
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
    hNeutronSpectrum[iCharge] = new TH1F(Form("h%sNeutronSpectrum", charge_label[iCharge]), ";#it{p} (GeV/#it{c});Entries", nPrimaryPtBins, primaryPtBins);
    fOutputList->Add(hNeutronSpectrum[iCharge]);
    hDeuteronSpectrum[iCharge] = new TH1F(Form("h%sDeuteronSpectrum", charge_label[iCharge]), ";#it{p} (GeV/#it{c});Entries", nPrimaryPtBins, primaryPtBins);
    fOutputList->Add(hDeuteronSpectrum[iCharge]);
    hNparticles[iCharge] = new TH2F(Form("h%sNparticles", charge_label[iCharge]), ";N_{#pi};N_{d}", 51, -0.5, 50.5, 4, -0.5, 3.5);
    fOutputList->Add(hNparticles[iCharge]);
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
    fOutputList->Add(hMixedEventPionDeuteronKstarUS[iCharge]);
  }

  // DeltaP Distribution
  hDeltaP = new TH1F("hDeltaP", ";#Delta p (GeV/#it{c});Entries", 200, 0, 1.);
  hDeltaP->Sumw2();
  fOutputList->Add(hDeltaP);

  PostData(1, fOutputList);
}

void AliAnalysisTaskPionDeuteronMC::UserExec(Option_t *)
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

  TClonesArray *stack = (TClonesArray *)ev->GetList()->FindObject(AliAODMCParticle::StdBranchName());
  if (!stack)
    ::Fatal("AliAnalysisTaskPionDeuteronMC::UserExec", "MC analysis requested on a sample without the MC particle array.");

  // Protons and Neutrons IDs
  vector<TLorentzVector> v_pospion, v_negpion, v_proton, v_antiproton, v_neutron, v_antineutron, v_deuteron, v_antideuteron;

  for (int iMC = 0; iMC < stack->GetEntriesFast(); iMC++)
  {
    AliAODMCParticle *part = (AliAODMCParticle *)stack->UncheckedAt(iMC);
    const int signed_pdg = part->GetPdgCode();
    float rapidity = TMath::Abs(part->Y());

    if (!part->IsPhysicalPrimary())
      continue;
    if (rapidity > 1.0)
      continue;

    float p = part->P();
    float px = part->Px();
    float py = part->Py();
    float pz = part->Pz();

    TLorentzVector fourvec;

    switch (signed_pdg)
    {
    case 2212:
      fourvec.SetXYZM(px, py, pz, TDatabasePDG::Instance()->GetParticle(2212)->Mass());
      v_proton.push_back(fourvec);
      if (rapidity <= 0.5)
      {
        hProtonSpectrum[0]->Fill(p);
      }
      break;
    case -2212:
      fourvec.SetXYZM(px, py, pz, TDatabasePDG::Instance()->GetParticle(2212)->Mass());
      v_antiproton.push_back(fourvec);
      if (rapidity <= 0.5)
      {
        hProtonSpectrum[1]->Fill(p);
      }
      break;
    case 2112:
      fourvec.SetXYZM(px, py, pz, TDatabasePDG::Instance()->GetParticle(2112)->Mass());
      v_neutron.push_back(fourvec);
      if (rapidity <= 0.5)
      {
        hNeutronSpectrum[0]->Fill(p);
      }
      break;
    case -2112:
      fourvec.SetXYZM(px, py, pz, TDatabasePDG::Instance()->GetParticle(2112)->Mass());
      v_antineutron.push_back(fourvec);
      if (rapidity <= 0.5)
      {
        hNeutronSpectrum[1]->Fill(p);
      }
      break;
    case 211:
      fourvec.SetXYZM(px, py, pz, TDatabasePDG::Instance()->GetParticle(211)->Mass());
      v_pospion.push_back(fourvec);
      if (rapidity <= 0.5)
      {
        hPionSpectrum[0]->Fill(p);
      }
      break;
    case -211:
      fourvec.SetXYZM(px, py, pz, TDatabasePDG::Instance()->GetParticle(211)->Mass());
      v_negpion.push_back(fourvec);
      if (rapidity <= 0.5)
      {
        hPionSpectrum[1]->Fill(p);
      }
      break;
    default:
      break;
    }
  }

  DoCoalescence(v_proton, v_neutron, v_deuteron, hDeuteronSpectrum[0], hDeltaP);             // deuteron
  DoCoalescence(v_antiproton, v_antineutron, v_antideuteron, hDeuteronSpectrum[1], hDeltaP); // antideuteron

  // count particles
  float nPosPions = (float)v_pospion.size();
  float nNegPions = (float)v_negpion.size();

  float nDeuterons = (float)v_deuteron.size();
  float nAntideuterons = (float)v_antideuteron.size();

  hNparticles[0]->Fill(nPosPions, nDeuterons);
  hNparticles[1]->Fill(nNegPions, nAntideuterons);

  // same event
  for (auto &pion : v_pospion)
  {
    for (auto &proton : v_proton)
    {
      hSameEventPionProtonKstarLS[0]->Fill(GetKstar(pion, proton));
    }
    for (auto &antiproton : v_antiproton)
    {
      hSameEventPionProtonKstarUS[0]->Fill(GetKstar(pion, antiproton));
    }
    for (auto &deuteron : v_deuteron)
    {
      hSameEventPionDeuteronKstarLS[0]->Fill(GetKstar(pion, deuteron));
    }
    for (auto &antideuteron : v_antideuteron)
    {
      hSameEventPionDeuteronKstarUS[0]->Fill(GetKstar(pion, antideuteron));
    }
  }

  for (auto &pion : v_negpion)
  {
    for (auto &proton : v_proton)
    {
      hSameEventPionProtonKstarUS[1]->Fill(GetKstar(pion, proton));
    }
    for (auto &antiproton : v_antiproton)
    {
      hSameEventPionProtonKstarLS[1]->Fill(GetKstar(pion, antiproton));
    }
    for (auto &deuteron : v_deuteron)
    {
      hSameEventPionDeuteronKstarUS[1]->Fill(GetKstar(pion, deuteron));
    }
    for (auto &antideuteron : v_antideuteron)
    {
      hSameEventPionDeuteronKstarLS[1]->Fill(GetKstar(pion, antideuteron));
    }
  }

  // event mixing
  float zvtx = ev->GetPrimaryVertex()->GetZ();
  float mult = header->GetRefMultiplicityComb08();
  int index = FindBin(zvtx, mult);

  if (index >= 0)
  {
    // match pions with particles in the buffer
    if (v_pospion.size() > 0)
    {
      FillMixedEvent(v_pospion, fMixingBufferProtons[index], hMixedEventPionProtonKstarLS[0]);
      FillMixedEvent(v_pospion, fMixingBufferAntiprotons[index], hMixedEventPionProtonKstarUS[0]);
      FillMixedEvent(v_pospion, fMixingBufferDeuterons[index], hMixedEventPionDeuteronKstarLS[0]);
      FillMixedEvent(v_pospion, fMixingBufferAntideuterons[index], hMixedEventPionDeuteronKstarUS[0]);
    }
    if (v_negpion.size() > 0)
    {
      FillMixedEvent(v_negpion, fMixingBufferProtons[index], hMixedEventPionProtonKstarUS[1]);
      FillMixedEvent(v_negpion, fMixingBufferAntiprotons[index], hMixedEventPionProtonKstarLS[1]);
      FillMixedEvent(v_negpion, fMixingBufferDeuterons[index], hMixedEventPionDeuteronKstarUS[1]);
      FillMixedEvent(v_negpion, fMixingBufferAntideuterons[index], hMixedEventPionDeuteronKstarLS[1]);
    }
    // update particle buffer
    if (v_proton.size() > 0)
    {
      fMixingBufferProtons[index].Fill(v_proton);
    }
    if (v_antiproton.size() > 0)
    {
      fMixingBufferAntiprotons[index].Fill(v_antiproton);
    }
    if (v_deuteron.size() > 0)
    {
      fMixingBufferDeuterons[index].Fill(v_deuteron);
    }
    if (v_antideuteron.size() > 0)
    {
      fMixingBufferAntideuterons[index].Fill(v_antideuteron);
    }
  }

  PostData(1, fOutputList);
}

void AliAnalysisTaskPionDeuteronMC::Terminate(Option_t *)
{
}

void AliAnalysisTaskPionDeuteronMC::SetMixingDepth(unsigned int depth)
{
  for (auto p : fMixingBufferDeuterons)
  {
    p.SetDepth(depth);
  }
  for (auto p : fMixingBufferAntideuterons)
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
