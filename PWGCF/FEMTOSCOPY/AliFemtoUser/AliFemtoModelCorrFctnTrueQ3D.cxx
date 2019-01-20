///
/// \file AliFemtoModelCorrFctnTrueQ3D.cxx
///

#include "AliFemtoModelCorrFctnTrueQ3D.h"
#include "AliFemtoPair.h"
#include "AliFemtoModelManager.h"
#include "AliFemtoModelHiddenInfo.h"

#include <TH3F.h>
#include <TList.h>
#include <TString.h>
#include <TRandom.h>

#include <tuple>
#include <iostream>

AliFemtoModelCorrFctnTrueQ3D::AliFemtoModelCorrFctnTrueQ3D():
  AliFemtoModelCorrFctnTrueQ3D("CF_TrueQ3D_")
{
}

AliFemtoModelCorrFctnTrueQ3D::AliFemtoModelCorrFctnTrueQ3D(const char *prefix):
  AliFemtoModelCorrFctnTrueQ3D(prefix, 56, 0.14)
{
}

AliFemtoModelCorrFctnTrueQ3D::AliFemtoModelCorrFctnTrueQ3D(const char *prefix,
                                                           UInt_t nbins,
                                                           Double_t qmax):
  AliFemtoModelCorrFctnTrueQ3D(prefix, nbins, -qmax, qmax)
{
}

AliFemtoModelCorrFctnTrueQ3D::AliFemtoModelCorrFctnTrueQ3D(const TString &prefix,
                                                           UInt_t nbins,
                                                           Double_t qmin,
                                                           Double_t qmax):
  AliFemtoCorrFctn()
  , fManager(nullptr)
  , fNumeratorGenerated(nullptr)
  , fNumeratorReconstructed(nullptr)
  , fNumeratorGenUnweighted(nullptr)
  , fNumeratorRecUnweighted(nullptr)
  , fDenominatorGenerated(nullptr)
  , fDenominatorReconstructed(nullptr)
  , fDenominatorGenWeighted(nullptr)
  , fDenominatorRecWeighted(nullptr)
  , fRng(new TRandom())
{
  fNumeratorGenerated = new TH3F(prefix + "NumGen",
                                 "Numerator (MC-Generated Momentum)",
                                 nbins, qmin, qmax,
                                 nbins, qmin, qmax,
                                 nbins, qmin, qmax);
  fNumeratorReconstructed = new TH3F(prefix + "NumRec",
                                     "Numerator (Reconstructed Momentum)",
                                     nbins, qmin, qmax,
                                     nbins, qmin, qmax,
                                     nbins, qmin, qmax);
  fNumeratorRecUnweighted = new TH3F(prefix + "NumRecUnweighted",
                                     "Numerator (Reconstructed Momentum - No femto weight)",
                                     nbins, qmin, qmax,
                                     nbins, qmin, qmax,
                                     nbins, qmin, qmax);
  fNumeratorGenUnweighted = new TH3F(prefix + "NumGenUnweighted",
                                     "Numerator (Generated Momentum - No femto weight)",
                                     nbins, qmin, qmax,
                                     nbins, qmin, qmax,
                                     nbins, qmin, qmax);
  fDenominatorGenerated = new TH3F(prefix + "DenGen",
                                   "Denominator (MC-Generated Momentum)",
                                   nbins, qmin, qmax,
                                   nbins, qmin, qmax,
                                   nbins, qmin, qmax);
  fDenominatorReconstructed = new TH3F(prefix + "DenRec",
                                       "Denominator (Reconstructed Momentum)",
                                       nbins, qmin, qmax,
                                       nbins, qmin, qmax,
                                       nbins, qmin, qmax);
  fDenominatorGenWeighted = new TH3F(prefix + "DenGenWeight",
                                     "Denominator (Generated Momentum - with femto weight)",
                                     nbins, qmin, qmax,
                                     nbins, qmin, qmax,
                                     nbins, qmin, qmax);
  fDenominatorRecWeighted = new TH3F(prefix + "DenRecWeight",
                                     "Numerator (Reconstructed Momentum - with femto weight)",
                                     nbins, qmin, qmax,
                                     nbins, qmin, qmax,
                                     nbins, qmin, qmax);
  fNumeratorGenerated->Sumw2();
  fNumeratorReconstructed->Sumw2();
  fDenominatorGenWeighted->Sumw2();
  fDenominatorRecWeighted->Sumw2();
}

AliFemtoModelCorrFctnTrueQ3D::AliFemtoModelCorrFctnTrueQ3D(UInt_t nbins, Double_t qmax):
  AliFemtoModelCorrFctnTrueQ3D(nbins, -abs(qmax), abs(qmax))
{
}

AliFemtoModelCorrFctnTrueQ3D::AliFemtoModelCorrFctnTrueQ3D(UInt_t nbins, Double_t qmin, Double_t qmax):
  AliFemtoModelCorrFctnTrueQ3D("", nbins, qmin, qmax)
{
}

AliFemtoModelCorrFctnTrueQ3D::AliFemtoModelCorrFctnTrueQ3D(const AliFemtoModelCorrFctnTrueQ3D::Parameters &params):
  AliFemtoModelCorrFctnTrueQ3D(params.prefix.Data(), params.bin_count, params.qmin, params.qmax)
{
  SetManager(params.mc_manager);
}


AliFemtoModelCorrFctnTrueQ3D::AliFemtoModelCorrFctnTrueQ3D(const AliFemtoModelCorrFctnTrueQ3D& orig):
  AliFemtoCorrFctn(orig)
  , fManager(orig.fManager)
  , fNumeratorGenerated(nullptr)
  , fNumeratorReconstructed(nullptr)
  , fNumeratorGenUnweighted(nullptr)
  , fNumeratorRecUnweighted(nullptr)
  , fDenominatorGenerated(nullptr)
  , fDenominatorReconstructed(nullptr)
  , fDenominatorGenWeighted(nullptr)
  , fDenominatorRecWeighted(nullptr)
  , fRng(new TRandom())
{
  fNumeratorGenerated = new TH3F(*orig.fNumeratorGenerated);
  fNumeratorReconstructed = new TH3F(*orig.fNumeratorReconstructed);
  fDenominatorGenerated = new TH3F(*orig.fDenominatorGenerated);
  fDenominatorReconstructed = new TH3F(*orig.fDenominatorReconstructed);

  fNumeratorGenUnweighted = new TH3F(*orig.fNumeratorGenUnweighted);
  fNumeratorRecUnweighted = new TH3F(*orig.fNumeratorRecUnweighted);
  fDenominatorGenWeighted = new TH3F(*orig.fDenominatorGenWeighted);
  fDenominatorRecWeighted = new TH3F(*orig.fDenominatorRecWeighted);
}


AliFemtoModelCorrFctnTrueQ3D&
AliFemtoModelCorrFctnTrueQ3D::operator=(const AliFemtoModelCorrFctnTrueQ3D&)
{
  return *this;
}

AliFemtoModelCorrFctnTrueQ3D::~AliFemtoModelCorrFctnTrueQ3D()
{
  delete fNumeratorGenerated;
  delete fNumeratorReconstructed;
  delete fDenominatorGenerated;
  delete fDenominatorReconstructed;
  delete fNumeratorGenUnweighted;
  delete fNumeratorRecUnweighted;
  delete fDenominatorGenWeighted;
  delete fDenominatorRecWeighted;
  delete fRng;
}


TList*
AliFemtoModelCorrFctnTrueQ3D::GetOutputList()
{
  TList *result = new TList();
  AppendOutputList(*result);
  return result;
}


TList*
AliFemtoModelCorrFctnTrueQ3D::AppendOutputList(TList &list)
{
  AddOutputObjectsTo(list);
  return &list;
}

void
AliFemtoModelCorrFctnTrueQ3D::AddOutputObjectsTo(TCollection &list)
{
  list.Add(fNumeratorGenerated);
  list.Add(fNumeratorReconstructed);
  list.Add(fDenominatorGenerated);
  list.Add(fDenominatorReconstructed);
  list.Add(fNumeratorGenUnweighted);
  list.Add(fNumeratorRecUnweighted);
  list.Add(fDenominatorGenWeighted);
  list.Add(fDenominatorRecWeighted);
}

/// Return q{Out-Side-Long} tuple, calculated from momentum vectors p1 & p2
static
std::tuple<Double_t, Double_t, Double_t>
Qcms(const AliFemtoLorentzVector &p1, const AliFemtoLorentzVector &p2)
{
  const AliFemtoLorentzVector p = p1 + p2,
                              d = p1 - p2;

  Double_t k1 = p.Perp(),
           k2 = d.x()*p.x() + d.y()*p.y();

  // relative momentum out component in lab frame
  Double_t qout = (k1 == 0) ? 0.0 : k2/k1;

  // relative momentum side component in lab frame
  Double_t qside = (k1 == 0) ? 0.0 : 2.0 * (p2.x()*p1.y() - p1.x()*p2.y())/k1;

  // relative momentum component in lab frame

  Double_t beta = p.z()/p.t(),
          gamma = 1.0 / TMath::Sqrt((1.0-beta)*(1.0+beta));

  Double_t qlong = gamma * (d.z() - beta*d.t());

  // double qlong = (p.t()*d.z() - p.z()*d.t()) / TMath::Sqrt(p.t()*p.t() - p.z()*p.z());

  return std::make_tuple(qout, qside, qlong);
}


static void
AddPair(const AliFemtoParticle &particle1,
        const AliFemtoParticle &particle2,
        TH3F &gen_hist,
        TH3F &rec_hist,
        TH3F &gen_hist_unw,
        TH3F &rec_hist_unw,
        Double_t weight)
{
  Double_t q_out, q_side, q_long;

  // Fill reconstructed histogram with "standard" particle momentum
  std::tie(q_out, q_side, q_long) = Qcms(particle1.FourMomentum(), particle2.FourMomentum());
  Int_t rec_bin = gen_hist.FindBin(q_out, q_long, q_side);
  if (!(gen_hist.IsBinOverflow(rec_bin) or gen_hist.IsBinUnderflow(rec_bin))) {
    rec_hist.Fill(q_out, q_side, q_long, weight);
    rec_hist_unw.Fill(q_out, q_side, q_long, 1.0);
  }

  // Get generated momentum from hidden info
  const AliFemtoModelHiddenInfo
    *info1 = dynamic_cast<const AliFemtoModelHiddenInfo*>(particle1.HiddenInfo()),
    *info2 = dynamic_cast<const AliFemtoModelHiddenInfo*>(particle2.HiddenInfo());

  if (info1 == nullptr || info2 == nullptr) {
    return;
  }

  const Float_t mass1 = info1->GetMass(),
                mass2 = info2->GetMass();

  // block all zero-mass particles from the correlation function
  if (mass1 == 0.0 || mass2 == 0.0) {
    return;
  }

  const AliFemtoThreeVector *true_momentum1 = info1->GetTrueMomentum(),
                            *true_momentum2 = info2->GetTrueMomentum();

  const Double_t e1 = sqrt(mass1 * mass1 + true_momentum1->Mag2()),
                 e2 = sqrt(mass2 * mass2 + true_momentum2->Mag2());

  const AliFemtoLorentzVector p1 = AliFemtoLorentzVector(e1, *true_momentum1),
                              p2 = AliFemtoLorentzVector(e2, *true_momentum2);

  // Fill generated-momentum histogram with "true" particle momentum
  std::tie(q_out, q_side, q_long) = Qcms(p1, p2);
  Int_t gen_bin = gen_hist.FindBin(q_out, q_long, q_side);
  if (!(gen_hist.IsBinOverflow(gen_bin) or gen_hist.IsBinUnderflow(gen_bin))) {
    gen_hist.Fill(q_out, q_side, q_long, weight);
    gen_hist_unw.Fill(q_out, q_side, q_long, weight);
  }
}

void
AliFemtoModelCorrFctnTrueQ3D::AddRealPair(AliFemtoPair *pair)
{
  const AliFemtoParticle *p1 = pair->Track1(),
                         *p2 = pair->Track2();

  // randomize to avoid ordering biases
  if (fRng->Uniform() >= 0.5) {
    std::swap(p1, p2);
  }
  AddPair(*p1, *p2,
          *fNumeratorGenUnweighted,
          *fNumeratorRecUnweighted,
          *fNumeratorGenerated,
          *fNumeratorReconstructed,
          fManager->GetWeight(pair));
}

void
AliFemtoModelCorrFctnTrueQ3D::AddMixedPair(AliFemtoPair *pair)
{
  const AliFemtoParticle *p1 = pair->Track1(),
                         *p2 = pair->Track2();

  // randomize to avoid ordering biases
  if (fRng->Uniform() >= 0.5) {
    std::swap(p1, p2);
  }

  AddPair(*p1, *p2,
          *fDenominatorGenWeighted,
          *fDenominatorRecWeighted,
          *fDenominatorGenerated,
          *fDenominatorReconstructed,
          fManager->GetWeight(pair));
}


AliFemtoString
AliFemtoModelCorrFctnTrueQ3D::Report()
{
  AliFemtoString report;
  return report;
}

void
AliFemtoModelCorrFctnTrueQ3D::Finish()
{
#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 10, 8)
  // remove {under,over}flow to shrink file sizes
  fNumeratorGenerated->ClearUnderflowAndOverflow();
  fNumeratorReconstructed->ClearUnderflowAndOverflow();
  fDenominatorGenerated->ClearUnderflowAndOverflow();
  fDenominatorReconstructed->ClearUnderflowAndOverflow();

  fNumeratorGenUnweighted->ClearUnderflowAndOverflow();
  fNumeratorRecUnweighted->ClearUnderflowAndOverflow();
  fDenominatorGenWeighted->ClearUnderflowAndOverflow();
  fDenominatorRecWeighted->ClearUnderflowAndOverflow();
#endif
}
