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
      fEventCuts{true},
      fOutputList{nullptr},
      fEstimator{1},
      fP0{0.285},
      fMixingDepth{5},
      fPrimaryPtBins{0},
      fKstarBins{0},
      fNormalisationHist{nullptr},
      hSourceSize{nullptr},
      hPionSpectrum{nullptr},
      hProtonSpectrum{nullptr},
      hNeutronSpectrum{nullptr},
      hDeuteronSpectrum{nullptr},
      hDeltaP{nullptr},
      hNparticles{nullptr},
      hSameEventKstarLS{nullptr},
      hSameEventKstarUS{nullptr},
      hMixedEventKstarLS{nullptr},
      hMixedEventKstarUS{nullptr},
      fMixingBufferDeuterons(),
      fMixingBufferAntideuterons()
{
  fZvtxArray = {-10., -8., -6., -4., -2., 0., 2., 4., 6., 8., 10.};
  fNZvtxBins = (int)fZvtxArray.size();
  fMultiplicityArray = {0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 80, 84, 88, 92, 96, 100};
  fNMultiplicityBins = (int)fMultiplicityArray.size();
  int BufferSize = fNZvtxBins * fNMultiplicityBins;
  fMixingBufferDeuterons.reserve(BufferSize);
  fMixingBufferAntideuterons.reserve(BufferSize);
  for (int i = 0; i < BufferSize; i++)
  {
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
    hSameEventKstarLS[iCharge] = new TH1F(Form("h%sSameEventKstarLS", charge_label[iCharge]), ";#it{k}* (GeV/#it{c});Entries", nKstarBins, kStarBins);
    fOutputList->Add(hSameEventKstarLS[iCharge]);
    hSameEventKstarUS[iCharge] = new TH1F(Form("h%sSameEventKstarUS", charge_label[iCharge]), ";#it{k}* (GeV/#it{c});Entries", nKstarBins, kStarBins);
    fOutputList->Add(hSameEventKstarUS[iCharge]);
    hMixedEventKstarLS[iCharge] = new TH1F(Form("h%sMixedEventKstarLS", charge_label[iCharge]), ";#it{k}* (GeV/#it{c});Entries", nKstarBins, kStarBins);
    fOutputList->Add(hMixedEventKstarLS[iCharge]);
    hMixedEventKstarUS[iCharge] = new TH1F(Form("h%sMixedEventKstarUS", charge_label[iCharge]), ";#it{k}* (GeV/#it{c});Entries", nKstarBins, kStarBins);
    fOutputList->Add(hMixedEventKstarUS[iCharge]);
  }

  // DeltaP Distribution
  hDeltaP = new TH1F("hDeltaP", ";#Delta p (GeV/#it{c});Entries", 200, 0, 1.);
  hDeltaP->Sumw2();
  fOutputList->Add(hDeltaP);

  PostData(1, fOutputList);
}

void AliAnalysisTaskPionDeuteronMC::UserExec(Option_t *)
{

  // Particle Masses in GeV
  double mpion = TDatabasePDG::Instance()->GetParticle(211)->Mass(); // Pion
  double mp = TDatabasePDG::Instance()->GetParticle(2212)->Mass();   // Proton
  double mn = TDatabasePDG::Instance()->GetParticle(2112)->Mass();   // Neutron
  double md = AliPID::ParticleMass(AliPID::kDeuteron);

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
  vector<int> proton_ID, antiproton_ID, neutron_ID, antineutron_ID;
  vector<int> neutron_status, antineutron_status;
  vector<TLorentzVector> v_pospion, v_negpion, v_deuteron, v_antideuteron;

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
      proton_ID.push_back(iMC);
      if (rapidity <= 0.5)
      {
        hProtonSpectrum[0]->Fill(p);
      }
      break;
    case -2212:
      antiproton_ID.push_back(iMC);
      if (rapidity <= 0.5)
      {
        hProtonSpectrum[1]->Fill(p);
      }
      break;
    case 2112:
      neutron_ID.push_back(iMC);
      neutron_status.push_back(0);
      if (rapidity <= 0.5)
      {
        hNeutronSpectrum[0]->Fill(p);
      }
      break;
    case -2112:
      antineutron_ID.push_back(iMC);
      antineutron_status.push_back(0);
      if (rapidity <= 0.5)
      {
        hNeutronSpectrum[1]->Fill(p);
      }
      break;
    case 211:
      fourvec.SetXYZM(px, py, pz, mpion);
      v_pospion.push_back(fourvec);
      if (rapidity <= 0.5)
      {
        hPionSpectrum[0]->Fill(p);
      }
      break;
    case -211:
      fourvec.SetXYZM(px, py, pz, mpion);
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

  //************************************************** SIMPLE COALESCENCE **************************************************//

  // deuteron

  for (int ip = 0; ip < (int)proton_ID.size(); ip++)
  {
    // Proton 4-Momentum
    AliAODMCParticle *proton = (AliAODMCParticle *)stack->UncheckedAt(proton_ID[ip]);
    TLorentzVector p_proton;
    p_proton.SetXYZM(proton->Px(), proton->Py(), proton->Pz(), mp);

    for (int in = 0; in < (int)neutron_ID.size(); in++)
    {

      // Neutron 4-Momentum
      AliAODMCParticle *neutron = (AliAODMCParticle *)stack->UncheckedAt(neutron_ID[in]);
      TLorentzVector p_neutron;
      p_neutron.SetXYZM(neutron->Px(), neutron->Py(), neutron->Pz(), mn);

      // Deuteron 4-Momentum
      TLorentzVector p_deuteron;
      p_deuteron.SetXYZM(p_proton.Px() + p_neutron.Px(), p_proton.Py() + p_neutron.Py(), p_proton.Pz() + p_neutron.Pz(), md);
      TVector3 boost_vector = p_deuteron.BoostVector();

      // Lorentz Transformations (from Lab to Deuteron Frame)
      TLorentzVector p_proton_prime = p_proton;
      p_proton_prime.Boost(-boost_vector);
      TLorentzVector p_neutron_prime = p_neutron;
      p_neutron_prime.Boost(-boost_vector);
      Double_t deltaP = (p_proton_prime - p_neutron_prime).P();

      // Fill DeltaP Distribution
      hDeltaP->Fill(deltaP);

      if (neutron_status[in] == 1)
        continue; // Skip already used neutrons

      // Simple Coalescence Condition
      if (deltaP < fP0)
      {
        v_deuteron.push_back(p_deuteron);
        neutron_status[in] = 1;
        double y_deuteron = p_deuteron.Rapidity();
        if (TMath::Abs(y_deuteron) < 0.5)
        {
          hDeuteronSpectrum[0]->Fill(p_deuteron.P());
        }
        break;
      }
    }
  }

  // antideuteron

  for (int ip = 0; ip < (int)antiproton_ID.size(); ip++)
  {
    // Antiproton 4-Momentum
    AliAODMCParticle *antiproton = (AliAODMCParticle *)stack->UncheckedAt(antiproton_ID[ip]);
    TLorentzVector p_antiproton;
    p_antiproton.SetXYZM(antiproton->Px(), antiproton->Py(), antiproton->Pz(), mp);

    for (int in = 0; in < (int)antineutron_ID.size(); in++)
    {

      // Antineutron 4-Momentum
      AliAODMCParticle *antineutron = (AliAODMCParticle *)stack->UncheckedAt(antineutron_ID[in]);
      TLorentzVector p_antineutron;
      p_antineutron.SetXYZM(antineutron->Px(), antineutron->Py(), antineutron->Pz(), mn);

      // Antiantideuteron 4-Momentum
      TLorentzVector p_antideuteron;
      p_antideuteron.SetXYZM(p_antiproton.Px() + p_antineutron.Px(), p_antiproton.Py() + p_antineutron.Py(), p_antiproton.Pz() + p_antineutron.Pz(), md);
      TVector3 boost_vector = p_antideuteron.BoostVector();

      // Lorentz Transformations (from Lab to Antideuteron Frame)
      TLorentzVector p_antiproton_prime = p_antiproton;
      p_antiproton_prime.Boost(-boost_vector);
      TLorentzVector p_antineutron_prime = p_antineutron;
      p_antineutron_prime.Boost(-boost_vector);
      Double_t deltaP = (p_antiproton_prime - p_antineutron_prime).P();

      // Fill DeltaP Distribution
      hDeltaP->Fill(deltaP);

      if (antineutron_status[in] == 1)
        continue; // Skip already used antineutrons

      // Simple Coalescence Condition
      if (deltaP < fP0)
      {
        v_antideuteron.push_back(p_antideuteron);
        antineutron_status[in] = 1;
        double y_antideuteron = p_antideuteron.Rapidity();
        if (TMath::Abs(y_antideuteron) < 0.5)
        {
          hDeuteronSpectrum[1]->Fill(p_antideuteron.P());
        }
        break;
      }
    }
  }

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
    for (auto &deuteron : v_deuteron)
    {
      hSameEventKstarLS[0]->Fill(GetKstar(pion, deuteron));
    }
    for (auto &antideuteron : v_antideuteron)
    {
      hSameEventKstarUS[0]->Fill(GetKstar(pion, antideuteron));
    }
  }

  for (auto &pion : v_negpion)
  {
    for (auto &antideuteron : v_antideuteron)
    {
      hSameEventKstarLS[1]->Fill(GetKstar(pion, antideuteron));
    }
    for (auto &deuteron : v_deuteron)
    {
      hSameEventKstarUS[1]->Fill(GetKstar(pion, deuteron));
    }
  }

  // event mixing

  float zvtx = ev->GetPrimaryVertex()->GetZ();
  float mult = header->GetRefMultiplicityComb08();
  int index = FindBin(zvtx, mult);

  if (index >= 0)
  {
    // match pions with deuterons in the buffer
    if (v_pospion.size() > 0)
    {
      FillMixedEvent(v_pospion, fMixingBufferDeuterons[index], hMixedEventKstarLS[0]);
      FillMixedEvent(v_pospion, fMixingBufferAntideuterons[index], hMixedEventKstarUS[0]);
    }
    if (v_negpion.size() > 0)
    {
      FillMixedEvent(v_negpion, fMixingBufferDeuterons[index], hMixedEventKstarUS[1]);
      FillMixedEvent(v_negpion, fMixingBufferAntideuterons[index], hMixedEventKstarLS[1]);
    }
    // update deuteron buffer
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

float AliAnalysisTaskPionDeuteronMC::GetKstar(TLorentzVector &p1, TLorentzVector &p2)
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

void AliAnalysisTaskPionDeuteronMC::SetPrimaryPtBins(int nbins, float min, float max)
{
  const float delta = (max - min) / nbins;
  fPrimaryPtBins.Set(nbins + 1);
  for (int iB = 0; iB < nbins; ++iB)
  {
    fPrimaryPtBins[iB] = min + iB * delta;
  }
  fPrimaryPtBins[nbins] = max;
}

void AliAnalysisTaskPionDeuteronMC::SetPrimaryPtBins(int nbins, float *bins)
{
  fPrimaryPtBins.Set(nbins + 1, bins);
}

void AliAnalysisTaskPionDeuteronMC::SetKstarBins(int nbins, float min, float max)
{
  const float delta = (max - min) / nbins;
  fKstarBins.Set(nbins + 1);
  for (int iB = 0; iB < nbins; ++iB)
  {
    fKstarBins[iB] = min + iB * delta;
  }
  fKstarBins[nbins] = max;
}

void AliAnalysisTaskPionDeuteronMC::SetKstarBins(int nbins, float *bins)
{
  fKstarBins.Set(nbins + 1, bins);
}

void AliAnalysisTaskPionDeuteronMC::SetZvtxArray(std::vector<int> &vec)
{
  fZvtxArray.clear();
  for (auto p : vec)
  {
    fZvtxArray.push_back(p);
  }
  fNZvtxBins = (int)fZvtxArray.size();
}

void AliAnalysisTaskPionDeuteronMC::SetMultiplicityArray(std::vector<int> &vec)
{
  fMultiplicityArray.clear();
  for (auto p : vec)
  {
    fMultiplicityArray.push_back(p);
  }
  fNMultiplicityBins = (int)fMultiplicityArray.size();
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
}

int AliAnalysisTaskPionDeuteronMC::FindBin(float zvtx, float mult)
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

void AliAnalysisTaskPionDeuteronMC::FillMixedEvent(std::vector<TLorentzVector> &vec, CustomQueue<std::vector<TLorentzVector>> &buffer, TH1F *histo)
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
