#include "AliEasyFemto.h"
#include "TMath.h"
#include "TChain.h"
#include "AliPIDResponse.h"
#include "TRandom3.h"

#include <utility>

ClassImp(AliEasyFemto);

AliEasyFemto::AliEasyFemto()
    : fP0{0.285},
      fMixingDepth{5},
      fSimpleCoalescence{true},
      fTwoGauss{false},
      fPrimaryPtBins{0},
      fKstarBins{0}
{
  fZvtxArray = {-10., -8., -6., -4., -2., 0., 2., 4., 6., 8., 10.};
  fNZvtxBins = (int)fZvtxArray.size();
  fMultiplicityArray = {0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 80, 84, 88, 92, 96, 100};
  fNMultiplicityBins = (int)fMultiplicityArray.size();
}

void AliEasyFemto::SetPrimaryPtBins(int nbins, float min, float max)
{
  const float delta = (max - min) / nbins;
  fPrimaryPtBins.Set(nbins + 1);
  for (int iB = 0; iB < nbins; ++iB)
  {
    fPrimaryPtBins[iB] = min + iB * delta;
  }
  fPrimaryPtBins[nbins] = max;
}

void AliEasyFemto::SetPrimaryPtBins(int nbins, float *bins)
{
  fPrimaryPtBins.Set(nbins + 1, bins);
}

void AliEasyFemto::SetKstarBins(int nbins, float min, float max)
{
  const float delta = (max - min) / nbins;
  fKstarBins.Set(nbins + 1);
  for (int iB = 0; iB < nbins; ++iB)
  {
    fKstarBins[iB] = min + iB * delta;
  }
  fKstarBins[nbins] = max;
}

void AliEasyFemto::SetKstarBins(int nbins, float *bins)
{
  fKstarBins.Set(nbins + 1, bins);
}

void AliEasyFemto::SetZvtxArray(std::vector<float> &vec)
{
  fZvtxArray.clear();
  for (auto p : vec)
  {
    fZvtxArray.push_back(p);
  }
  fNZvtxBins = (int)fZvtxArray.size();
}

void AliEasyFemto::SetMultiplicityArray(std::vector<float> &vec)
{
  fMultiplicityArray.clear();
  for (auto p : vec)
  {
    fMultiplicityArray.push_back(p);
  }
  fNMultiplicityBins = (int)fMultiplicityArray.size();
}

int AliEasyFemto::FindBin(float zvtx, float mult)
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

float AliEasyFemto::GetKstar(TLorentzVector &p1, TLorentzVector &p2)
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

void AliEasyFemto::FillMixedEvent(std::vector<TLorentzVector> &vec, CustomQueue<std::vector<TLorentzVector>> &buffer, TH1F *histo)
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

bool AliEasyFemto::applyTrackSelection(AliAODTrack *track, CutContainer &cuts, float &dcaXY, TH1F *hTrackSelections)
{
  if (hTrackSelections)
  {
    hTrackSelections->Fill(kTrackSelection::kNoSelection);
  }
  if (!track->TestFilterBit(cuts.filterBit))
  {
    return false;
  }
  if (hTrackSelections)
  {
    hTrackSelections->Fill(kTrackSelection::kFilterBit);
  }
  if (TMath::Abs(track->Eta()) > cuts.eta)
  {
    return false;
  }
  if (hTrackSelections)
  {
    hTrackSelections->Fill(kTrackSelection::kEta);
  }
  if (track->Pt() < cuts.minPt)
  {
    return false;
  }
  if (hTrackSelections)
  {
    hTrackSelections->Fill(kTrackSelection::kPtMin);
  }
  if (track->Pt() > cuts.maxPt)
  {
    return false;
  }
  if (hTrackSelections)
  {
    hTrackSelections->Fill(kTrackSelection::kPtMax);
  }
  if (track->GetTPCNcls() < cuts.nTPCcls)
  {
    return false;
  }
  if (hTrackSelections)
  {
    hTrackSelections->Fill(kTrackSelection::kTPCnCls);
  }
  if (track->GetTPCnclsS() > 0)
  {
    return false;
  }
  if (hTrackSelections)
  {
    hTrackSelections->Fill(kTrackSelection::kTPCnClsS);
  }
  if (track->GetTPCNCrossedRows() < (int)cuts.nCrossedRows)
  {
    return false;
  }
  if (hTrackSelections)
  {
    hTrackSelections->Fill(kTrackSelection::kTPCnCrossedRows);
  }
  if ((float)track->GetTPCNCrossedRows() / (float)track->GetTPCNclsF() < cuts.crossedRowsOverFindable)
  {
    return false;
  }
  if (hTrackSelections)
  {
    hTrackSelections->Fill(kTrackSelection::kTPCnCrossedOverFindable);
  }
  float dcaZ = 0.;
  track->GetImpactParameters(dcaXY, dcaZ);
  if (TMath::Abs(dcaZ) > cuts.dcaZ)
  {
    return false;
  }
  if (hTrackSelections)
  {
    hTrackSelections->Fill(kTrackSelection::kDCAz);
  }
  return true;
}

bool AliEasyFemto::applyStandardPID(AliAODTrack *track, AliPIDResponse *fPID, CutContainer &cuts, AliPID::EParticleType kType)
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
    float beta = hasTOF(track, fPID);
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

float AliEasyFemto::hasTOF(AliAODTrack *track, AliPIDResponse *fPID)
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

float AliEasyFemto::ComputeRadius(float mt)
{
  return 0.8 + TMath::Exp(1.4 - 1.9 * mt);
}

float AliEasyFemto::ComputeSimpleGaussProbability(float q, float radius, float d)
{
  return 3 * TMath::Power(d, 3) / TMath::Power(d * d + radius * radius, 1.5) * TMath::Exp(-1 * q * q * d * d * 5.068 * 5.068);
}

float AliEasyFemto::ComputeTwoGaussProbability(float q, float radius, float d1, float d2, float Delta)
{
  return 3 * (TMath::Power(d1, 3) / TMath::Power(d1 * d1 + 4 * radius * radius, 1.5) * TMath::Exp(-1 * q * q * d1 * d1 * 5.068 * 5.068) + (1. - Delta) * TMath::Power(d2, 3) / TMath::Power(d2 * d2 + 4 * radius * radius, 1.5) * TMath::Exp(-1 * q * q * d2 * d2 * 5.068 * 5.068));
}

TLorentzVector AliEasyFemto::BaseCoalescence(TLorentzVector &a, TLorentzVector &b, float &deltaP, float &mt)
{
  TLorentzVector deuteron;
  deuteron.SetXYZM(a.Px() + b.Px(), a.Py() + b.Px(), a.Pz() + b.Pz(), AliPID::ParticleMass(AliPID::kDeuteron));
  TVector3 boost_vector = deuteron.BoostVector();

  // Lorentz Transformations (from Lab to Deuteron Frame)
  TLorentzVector a_prime = a;
  a_prime.Boost(-boost_vector);
  TLorentzVector b_prime = b;
  b_prime.Boost(-boost_vector);
  TLorentzVector deltaPvector = a_prime - b_prime;
  deltaP = deltaPvector.P();
  mt = 0.5 * deltaPvector.Mt();

  return deuteron;
}

void AliEasyFemto::DoCoalescence(std::vector<TLorentzVector> &v_proton, std::vector<TLorentzVector> &v_neutron, std::vector<TLorentzVector> &v_deuteron, TH1F *hDeuteron, TH1F *hDeltaP)
{

  int n_neutrons = (int)v_neutron.size();
  std::vector<int> neutron_status(n_neutrons, 0);

  for (auto &prot : v_proton)
  {
    for (int in = 0; in < n_neutrons; in++)
    {
      float mt, deltaP;
      TLorentzVector deuteron = BaseCoalescence(prot, v_neutron[in], deltaP, mt);

      // Fill DeltaP Distribution
      hDeltaP->Fill(deltaP);

      if (neutron_status[in] == 1)
        continue; // Skip already used neutrons
      if (fSimpleCoalescence)
      {
        // Simple Coalescence Condition
        if (deltaP < fP0)
        {
          v_deuteron.push_back(deuteron);
          neutron_status[in] = 1;
          hDeuteron->Fill(deuteron.P());
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
          prob = ComputeTwoGaussProbability(q, radius, 3.979, 0.890, 0.581);
        }
        else
        {
          prob = ComputeSimpleGaussProbability(q, radius, 3.2);
        }
        if (value < prob)
        {
          v_deuteron.push_back(deuteron);
          neutron_status[in] = 1;
          hDeuteron->Fill(deuteron.P());
          break;
        }
      }
    }
  }
}

void AliEasyFemto::DoFakeCoalescence(std::vector<TLorentzVector> &v_proton, std::vector<TLorentzVector> &v_fakedeuteron, std::vector<std::pair<int, int>> &v_fakedeuteronID, TH1F *hDeuteron, TH1F *hDeltaP)
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
      float mt, deltaP;
      TLorentzVector deuteron = BaseCoalescence(v_proton[i], v_proton[j], deltaP, mt);

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
          hDeuteron->Fill(deuteron.P());
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
          prob = ComputeTwoGaussProbability(q, radius, 3.979, 0.890, 0.581);
        }
        else
        {
          prob = ComputeSimpleGaussProbability(q, radius, 3.2);
        }
        if (value < prob)
        {
          v_fakedeuteron.push_back(deuteron);
          v_fakedeuteronID.push_back(std::make_pair(i, j));
          v_proton_status[j] = 1;
          hDeuteron->Fill(deuteron.P());
          break;
        }
      }
    }
  }
}
