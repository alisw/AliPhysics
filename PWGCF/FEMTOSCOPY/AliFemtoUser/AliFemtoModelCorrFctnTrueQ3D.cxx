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

#define OPTIMIZE_EXPECT_FALSE(expr) __builtin_expect(expr, 0)


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
  AliFemtoModelCorrFctnTrueQ3D(prefix, nbins, -qmax, qmax, false, false)
{
}


AliFemtoModelCorrFctnTrueQ3D
  ::AliFemtoModelCorrFctnTrueQ3D(const TString &prefix,
                                 UInt_t nbins,
                                 Double_t qmin,
                                 Double_t qmax,
                                 Bool_t enable_extra_hists,
                                 Bool_t enable_extra_denominators):
  AliFemtoModelCorrFctnTrueQ3D(prefix,
                               nbins / 2 + 1, nbins, nbins,
                               qmax, qmax, qmax,
                               enable_extra_hists,
                               enable_extra_denominators)
{
}

AliFemtoModelCorrFctnTrueQ3D
  ::AliFemtoModelCorrFctnTrueQ3D(const TString &prefix,
                                 UInt_t nbo,
                                 UInt_t nbs,
                                 UInt_t nbl,
                                 Double_t qoutmax,
                                 Double_t qsidemax,
                                 Double_t qlongmax,
                                 Bool_t enable_extra_hists,
                                 Bool_t enable_extra_denominators):
  AliFemtoModelCorrFctnTrueQ3D(prefix,
                               nbo,
                               nbs,
                               nbl,
                               0.0,
                               qoutmax,
                               -qsidemax,
                               qsidemax,
                               -qlongmax,
                               qlongmax,
                               enable_extra_hists,
                               enable_extra_denominators)
{
}

AliFemtoModelCorrFctnTrueQ3D
  ::AliFemtoModelCorrFctnTrueQ3D(const TString &prefix,
                                 UInt_t nbo,
                                 UInt_t nbs,
                                 UInt_t nbl,
                                 Double_t qoutmin,
                                 Double_t qoutmax,
                                 Double_t qsidemin,
                                 Double_t qsidemax,
                                 Double_t qlongmin,
                                 Double_t qlongmax,
                                 Bool_t enable_extra_hists,
                                 Bool_t enable_extra_denominators):
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
{

  auto new_th3 = [&] (const TString &name, const TString &title)
    {
      return new TH3F(prefix + name,
                      title + "; q_{out} (GeV); q_{side} (GeV); q_{long} (Gev)",
                      nbo, 0.0, qoutmax,
                      nbs, qsidemin, qsidemax,
                      nbl, qlongmin, qlongmax);
    };

  fNumeratorGenerated = new_th3("NumGen", "Numerator (MC-Generated Momentum)");
  fNumeratorReconstructed = new_th3("NumRec", "Numerator (Reconstructed Momentum)");
  fDenominatorGenerated = new_th3("DenGen", "Denominator (MC-Generated Momentum)");
  fDenominatorReconstructed = new_th3("DenRec", "Denominator (Reconstructed Momentum)");

  fNumeratorGenerated->Sumw2();
  fNumeratorReconstructed->Sumw2();

  if (enable_extra_hists) {
    fNumeratorRecUnweighted = new_th3("NumRecUnweighted", "Numerator (Reconstructed Momentum - No femto weight)");
    fNumeratorGenUnweighted = new_th3("NumGenUnweighted", "Numerator (Generated Momentum - No femto weight)");
  }

  if (enable_extra_hists || enable_extra_denominators) {
    fDenominatorGenWeighted = new_th3("DenGenWeight", "Denominator (Generated Momentum - with femto weight)");
    fDenominatorRecWeighted = new_th3("DenRecWeight", "Denominator (Reconstructed Momentum - with femto weight)");

    fDenominatorGenWeighted->Sumw2();
    fDenominatorRecWeighted->Sumw2();
  }
}

AliFemtoModelCorrFctnTrueQ3D::AliFemtoModelCorrFctnTrueQ3D(UInt_t nbins, Double_t qmax):
  AliFemtoModelCorrFctnTrueQ3D(nbins, -abs(qmax), abs(qmax))
{
}

AliFemtoModelCorrFctnTrueQ3D::AliFemtoModelCorrFctnTrueQ3D(UInt_t nbins, Double_t qmin, Double_t qmax):
  AliFemtoModelCorrFctnTrueQ3D("", nbins, qmin, qmax)
{
}

AliFemtoModelCorrFctnTrueQ3D::AliFemtoModelCorrFctnTrueQ3D(const Parameters &params):
  AliFemtoModelCorrFctnTrueQ3D(params.prefix.Data(),
                               params.bin_count_out,
                               params.bin_count_side,
                               params.bin_count_long,
                               params.qomin,
                               params.qomax,
                               params.qsmin,
                               params.qsmax,
                               params.qlmin,
                               params.qlmax,
                               params.enable_extra_hists,
                               params.enable_extra_denoms)
{
  SetManager(params.mc_manager);
}


AliFemtoModelCorrFctnTrueQ3D::
  AliFemtoModelCorrFctnTrueQ3D(const TString &prefix,
                               const std::vector<double> &obins,
                               const std::vector<double> &sbins,
                               const std::vector<double> &lbins,
                               AliFemtoModelManager *mgr)
  : AliFemtoCorrFctn()
  , fManager(mgr)
  , fNumeratorGenerated(nullptr)
  , fNumeratorReconstructed(nullptr)
  , fNumeratorGenUnweighted(nullptr)
  , fNumeratorRecUnweighted(nullptr)
  , fDenominatorGenerated(nullptr)
  , fDenominatorReconstructed(nullptr)
  , fDenominatorGenWeighted(nullptr)
  , fDenominatorRecWeighted(nullptr)
{

  auto new_th3 = [&] (const TString &name, const TString &title)
    {
      return new TH3F(prefix + name,
                      title + "; q_{out} (GeV); q_{side} (GeV); q_{long} (Gev)",
                      obins.size()-1, obins.data(),
                      sbins.size()-1, sbins.data(),
                      lbins.size()-1, lbins.data());
    };

  fNumeratorGenerated = new_th3("NumGen", "Numerator (MC-Generated Momentum)");
  fNumeratorReconstructed = new_th3("NumRec", "Numerator (Reconstructed Momentum)");
  fDenominatorGenerated = new_th3("DenGen", "Denominator (MC-Generated Momentum)");
  fDenominatorReconstructed = new_th3("DenRec", "Denominator (Reconstructed Momentum)");

  fNumeratorGenerated->Sumw2();
  fNumeratorReconstructed->Sumw2();

  fNumeratorRecUnweighted = new_th3("NumRecUnweighted", "Numerator (Reconstructed Momentum - No femto weight)");
  fNumeratorGenUnweighted = new_th3("NumGenUnweighted", "Numerator (Generated Momentum - No femto weight)");

  fDenominatorGenWeighted = new_th3("DenGenWeight", "Denominator (Generated Momentum - with femto weight)");
  fDenominatorRecWeighted = new_th3("DenRecWeight", "Denominator (Reconstructed Momentum - with femto weight)");

  fDenominatorGenWeighted->Sumw2();
  fDenominatorRecWeighted->Sumw2();
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
{
  fNumeratorGenerated = new TH3F(*orig.fNumeratorGenerated);
  fNumeratorReconstructed = new TH3F(*orig.fNumeratorReconstructed);
  fDenominatorGenerated = new TH3F(*orig.fDenominatorGenerated);
  fDenominatorReconstructed = new TH3F(*orig.fDenominatorReconstructed);

  if (orig.fNumeratorGenUnweighted) {
    fNumeratorGenUnweighted = new TH3F(*orig.fNumeratorGenUnweighted);
    fNumeratorRecUnweighted = new TH3F(*orig.fNumeratorRecUnweighted);
  }
  if (orig.fDenominatorGenWeighted) {
    fDenominatorGenWeighted = new TH3F(*orig.fDenominatorGenWeighted);
    fDenominatorRecWeighted = new TH3F(*orig.fDenominatorRecWeighted);
  }
}


AliFemtoModelCorrFctnTrueQ3D&
AliFemtoModelCorrFctnTrueQ3D::operator=(const AliFemtoModelCorrFctnTrueQ3D &rhs)
{
  AliFemtoCorrFctn::operator=(rhs);

  fManager = rhs.fManager;

  *fNumeratorGenerated = *rhs.fNumeratorGenerated;
  *fNumeratorReconstructed = *rhs.fNumeratorReconstructed;
  *fDenominatorGenerated = *rhs.fDenominatorGenerated;
  *fDenominatorReconstructed = *rhs.fDenominatorReconstructed;

  auto copy_if_present = [] (const TH3F* src, TH3F *&dest)
    {
      if (src and dest) {
        *dest = *src;
      }
      else if (src) {
        dest = new TH3F(*src);
      }
      else if (dest) {
        delete dest;
        dest = nullptr;
      }
    };

  copy_if_present(rhs.fNumeratorGenUnweighted, fNumeratorGenUnweighted);
  copy_if_present(rhs.fNumeratorRecUnweighted, fNumeratorRecUnweighted);

  copy_if_present(rhs.fDenominatorGenWeighted, fDenominatorGenWeighted);
  copy_if_present(rhs.fDenominatorRecWeighted, fDenominatorRecWeighted);

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

  if (fNumeratorGenUnweighted) {
    list.Add(fNumeratorGenUnweighted);
    list.Add(fNumeratorRecUnweighted);
  }
  if (fDenominatorGenWeighted) {
    list.Add(fDenominatorGenWeighted);
    list.Add(fDenominatorRecWeighted);
  }
}

/// Return q{Out-Side-Long} tuple, calculated from momentum vectors p1 & p2
static
std::tuple<Double_t, Double_t, Double_t>
Qlcms(const AliFemtoLorentzVector &p1, const AliFemtoLorentzVector &p2)
{
  const AliFemtoLorentzVector p = p1 + p2,
                              d = p1 - p2;

  #define FAST_DIVIDE(num, den) __builtin_expect(den == 0.0, false) ? 0.0 : num / den

  const Double_t
    pt = p.Perp(),
    beta = p.z()/p.t(),
    gamma = 1.0 / TMath::Sqrt(1.0-beta*beta),

    qout = FAST_DIVIDE(d.x()*p.x() + d.y()*p.y(), pt),
    qside = FAST_DIVIDE(2.0 * (p2.x()*p1.y() - p1.x()*p2.y()), pt),
    qlong = gamma * (d.z() - beta*d.t());

  #undef FAST_DIVIDE

  // "flip" into positive qout region
  const Double_t f = std::copysign(1.0, qout);

  return std::make_tuple(f * qout, f * qside, f * qlong);
}


/// Calculates & fills historams with q-{Out,Side,Long} calculated
/// from the momentum given. First histogram gets weight, second gets
/// none
static void
fill_hists(TH3 *dest,
           TH3 *dest_unweighted,
           const AliFemtoLorentzVector &p1,
           const AliFemtoLorentzVector &p2,
           double weight)
{
  Double_t q_out, q_side, q_long;
  std::tie(q_out, q_side, q_long) = Qlcms(p1, p2);

  TH3 *hist = dest ? dest : dest_unweighted;

  const Int_t dest_bin = hist->FindBin(q_out, q_long, q_side);

  if (!hist->IsBinOverflow(dest_bin) and !hist->IsBinUnderflow(dest_bin)) {
    if (dest) {
      dest->Fill(q_out, q_side, q_long, weight);
    }
    if (dest_unweighted) {
      dest_unweighted->Fill(q_out, q_side, q_long, 1.0);
    }
  }
}


/// Adds pair to appropriate histograms
///
/// Used by AddRealPair & AddMixedPair
///
static void
AddPair(const AliFemtoParticle &particle1,
        const AliFemtoParticle &particle2,
        TH3F *gen_hist,
        TH3F *rec_hist,
        TH3F *gen_hist_unw,
        TH3F *rec_hist_unw,
        Double_t weight)
{
  // fill with reconstructed momentum
  fill_hists(rec_hist,
             rec_hist_unw,
             particle1.FourMomentum(),
             particle2.FourMomentum(),
             weight);

  // Get generated momentum from hidden info
  const AliFemtoModelHiddenInfo
    *info1 = dynamic_cast<const AliFemtoModelHiddenInfo*>(particle1.HiddenInfo()),
    *info2 = dynamic_cast<const AliFemtoModelHiddenInfo*>(particle2.HiddenInfo());

  if (OPTIMIZE_EXPECT_FALSE(info1 == nullptr || info2 == nullptr)) {
    return;
  }

  const Float_t mass1 = info1->GetMass(),
                mass2 = info2->GetMass();

  // block all zero-mass particles from the correlation function
  if (OPTIMIZE_EXPECT_FALSE(mass1 == 0.0 || mass2 == 0.0)) {
    return;
  }

  const AliFemtoThreeVector *true_momentum1 = info1->GetTrueMomentum(),
                            *true_momentum2 = info2->GetTrueMomentum();

  const Double_t e1 = sqrt(mass1 * mass1 + true_momentum1->Mag2()),
                 e2 = sqrt(mass2 * mass2 + true_momentum2->Mag2());

  // fill with monte-carlo generated momentum
  fill_hists(gen_hist,
             gen_hist_unw,
             AliFemtoLorentzVector(e1, *true_momentum1),
             AliFemtoLorentzVector(e2, *true_momentum2),
             weight);
}


void
AliFemtoModelCorrFctnTrueQ3D::AddRealPair(AliFemtoPair *pair)
{
  const AliFemtoParticle *p1 = pair->Track1(),
                         *p2 = pair->Track2();

  AddPair(*p1, *p2,
          fNumeratorGenerated,
          fNumeratorReconstructed,
          fNumeratorGenUnweighted,
          fNumeratorRecUnweighted,
          fManager->GetWeight(pair));
}


void
AliFemtoModelCorrFctnTrueQ3D::AddMixedPair(AliFemtoPair *pair)
{
  const AliFemtoParticle *p1 = pair->Track1(),
                         *p2 = pair->Track2();

  // skip weight calculation if we are not storing it
  const double weight = fDenominatorGenWeighted != nullptr
                      ? fManager->GetWeight(pair)
                      : 0.0;

  AddPair(*p1, *p2,
          fDenominatorGenWeighted,
          fDenominatorRecWeighted,
          fDenominatorGenerated,
          fDenominatorReconstructed,
          weight);
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
}
