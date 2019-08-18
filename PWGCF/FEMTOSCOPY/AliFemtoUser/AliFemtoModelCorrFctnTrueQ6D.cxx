///
/// \file AliFemtoModelCorrFctnTrueQ6D.cxx
///

#include "AliFemtoModelCorrFctnTrueQ6D.h"
#include "AliFemtoPair.h"
#include "AliFemtoModelHiddenInfo.h"

#include <TList.h>
#include <TString.h>
#include <TRandom.h>
#include <TAxis.h>

#include <array>
#include <tuple>
#include <iostream>
#include <exception>


void
AliFemtoModelCorrFctnTrueQ6D::UpdateQlimits()
{
  using Method = BinMethod;

  std::array<Int_t, 3> axes;

  switch (fBinMethod) {
  case Method::kGenRecLSO:
    axes = {2, 1, 0};
    break;
  case Method::kGenRecOSL:
    axes = {0, 1, 2};
    break;
  case Method::kGenLSORecOSL:
    axes = {2, 1, 0};
    break;

  case Method::kRecGenLSO:
    axes = {5, 4, 3};
    break;
  case Method::kRecGenOSL:
    axes = {3, 4, 5};
    break;
  case Method::kRecLSOGenOSL:
    axes = {3, 4, 5};
    break;

  case Method::kGroupedAxisLSO:
    axes = {5, 3, 1};
    break;
  case Method::kGroupedAxisOSL:
    axes = {0, 2, 4};
    break;

  default:
    throw std::invalid_argument(std::string(Form("Uknown binning method %d", fBinMethod)));
  }

  const auto aout = fHistogram->GetAxis(axes[0]),
             aside = fHistogram->GetAxis(axes[1]),
             along = fHistogram->GetAxis(axes[2]);

  fQlimits[0] = make_pair(aout->GetXmin(), aout->GetXmax());
  fQlimits[1] = make_pair(aside->GetXmin(), aside->GetXmax());
  fQlimits[2] = make_pair(along->GetXmin(), along->GetXmax());
}

AliFemtoModelCorrFctnTrueQ6D::Builder::Builder(AliFemtoConfigObject cfg)
  : Builder()
{
  cfg.pop_all()
    ("name", title)
    ("qmin", qmin)
    ("qmax", qmax)
    ("nbins", bin_count);
}

AliFemtoModelCorrFctnTrueQ6D::AliFemtoModelCorrFctnTrueQ6D()
  : AliFemtoModelCorrFctnTrueQ6D("CF_TrueQ6D")
{
}

AliFemtoModelCorrFctnTrueQ6D::AliFemtoModelCorrFctnTrueQ6D(const HistType &hist)
  : AliFemtoCorrFctn()
  , fHistogram(static_cast<HistType*>(hist.Clone()))
  , fBinMethod(kRecGenOSL)
  , fIgnoreZeroMassParticles(true)
{
  fBinMethod = GuessBinMethod(*fHistogram);
  UpdateQlimits();
}

AliFemtoModelCorrFctnTrueQ6D::AliFemtoModelCorrFctnTrueQ6D(HistType *&hist)
  : AliFemtoCorrFctn()
  , fHistogram(hist)
  , fBinMethod(kRecGenOSL)
  , fIgnoreZeroMassParticles(true)
{
  // take ownership of the pointer
  hist = nullptr;
  fBinMethod = GuessBinMethod(*fHistogram);
  UpdateQlimits();
}

template <typename T>
std::array<T, 6> sort_q(AliFemtoModelCorrFctnTrueQ6D::BinMethod bm,
   T q_out, T q_side, T q_long, T true_q_out, T true_q_side, T true_q_long)
{
  using Method = AliFemtoModelCorrFctnTrueQ6D::BinMethod;

  std::array<T, 6> result;

  switch (bm) {
  case Method::kGenRecLSO:
    result = {{true_q_long, true_q_side, true_q_out, q_long, q_side, q_out}};
    break;
  case Method::kGenRecOSL:
    result = {{true_q_out, true_q_side, true_q_long, q_out, q_side, q_long}};
    break;
  case Method::kGenLSORecOSL:
    result = {{true_q_long, true_q_side, true_q_out, q_out, q_side, q_long}};
    break;

  case Method::kRecGenLSO:
    result = {{q_long, q_side, q_out, true_q_long, true_q_side, true_q_out}};
    break;
  case Method::kRecGenOSL:
    result = {{q_out, q_side, q_long, true_q_out, true_q_side, true_q_long}};
    break;
  case Method::kRecLSOGenOSL:
    result = {{q_long, q_side, q_out, true_q_out, true_q_side, true_q_long}};
    break;

  case Method::kGroupedAxisLSO:
    result = {{q_long, true_q_long, q_side, true_q_side, q_out, true_q_out}};
    break;
  case Method::kGroupedAxisOSL:
    result = {{true_q_out, q_out, true_q_side, q_side, true_q_long, q_long}};
    break;

  default:
    throw std::invalid_argument(std::string(Form("Uknown binning method %d", bm)));
  }

  return result;
}

inline AliFemtoModelCorrFctnTrueQ6D::BinMethod
AliFemtoModelCorrFctnTrueQ6D::GuessBinMethod(const HistType &hist)
{
  using Method = AliFemtoModelCorrFctnTrueQ6D::BinMethod;

  std::array<TString, 6> titles {{ hist.GetAxis(0)->GetTitle(),
                                   hist.GetAxis(1)->GetTitle(),
                                   hist.GetAxis(2)->GetTitle(),
                                   hist.GetAxis(3)->GetTitle(),
                                   hist.GetAxis(4)->GetTitle(),
                                   hist.GetAxis(5)->GetTitle() }};

  const Method methods[] = {
    Method::kGenRecLSO,
    Method::kGenRecOSL,
    Method::kGenLSORecOSL,
    Method::kRecGenLSO,
    Method::kRecGenOSL,
    Method::kRecLSOGenOSL,
    Method::kGroupedAxisOSL,
    Method::kGroupedAxisLSO
  };

  for (auto method : methods) {
    if (titles == sort_q(method, TString("q_{o}"), TString("q_{s}"), TString("q_{l}"),
                                 TString("q_{t,o}"), TString("q_{t,s}"), TString("q_{t,l}"))) {
      return method;
    }
  }

 return Method::kGenLSORecOSL;
}


AliFemtoModelCorrFctnTrueQ6D::AliFemtoModelCorrFctnTrueQ6D(const TString &title)
  : AliFemtoModelCorrFctnTrueQ6D(title, 120, 0.3)
{
}


AliFemtoModelCorrFctnTrueQ6D::AliFemtoModelCorrFctnTrueQ6D(const TString &title,
                                                           UInt_t nbins,
                                                           Double_t qmax)
  : AliFemtoModelCorrFctnTrueQ6D(title, nbins, -qmax, qmax)
{
}


AliFemtoModelCorrFctnTrueQ6D
  ::AliFemtoModelCorrFctnTrueQ6D(const TString &prefix,
                                 UInt_t nbins,
                                 Double_t qmin,
                                 Double_t qmax,
                                 BinMethod binning)
  : AliFemtoModelCorrFctnTrueQ6D(prefix,
                                 nbins, qmin, qmax,
                                 nbins, qmin, qmax,
                                 nbins, qmin, qmax,
                                 binning)
{
}


AliFemtoModelCorrFctnTrueQ6D
  ::AliFemtoModelCorrFctnTrueQ6D(const TString &prefix,
                                 Int_t nbins_out, Double_t qout_lo, Double_t qout_hi,
                                 Int_t nbins_side, Double_t qside_lo, Double_t qside_hi,
                                 Int_t nbins_long, Double_t qlong_lo, Double_t qlong_hi,
                                 BinMethod binning)
  : AliFemtoModelCorrFctnTrueQ6D(prefix,
                                 nbins_out, qout_lo, qout_hi,
                                 nbins_side, qside_lo, qside_hi,
                                 nbins_long, qlong_lo, qlong_hi,
                                 nbins_out, qout_lo, qout_hi,
                                 nbins_side, qside_lo, qside_hi,
                                 nbins_long, qlong_lo, qlong_hi,
                                 binning)
{
}


AliFemtoModelCorrFctnTrueQ6D
  ::AliFemtoModelCorrFctnTrueQ6D(const TString &prefix,
                                 Int_t nbins_out, Double_t qout_lo, Double_t qout_hi,
                                 Int_t nbins_side, Double_t qside_lo, Double_t qside_hi,
                                 Int_t nbins_long, Double_t qlong_lo, Double_t qlong_hi,
                                 Int_t nbins_out_true, Double_t qout_lo_true, Double_t qout_hi_true,
                                 Int_t nbins_side_true, Double_t qside_lo_true, Double_t qside_hi_true,
                                 Int_t nbins_long_true, Double_t qlong_lo_true, Double_t qlong_hi_true,
                                 BinMethod binning)
  : AliFemtoCorrFctn()
  , fHistogram(nullptr)
  , fBinMethod(binning)
  , fIgnoreZeroMassParticles(true)
{
  auto nbins_v = sort_q(fBinMethod, nbins_out, nbins_side, nbins_long, nbins_out_true, nbins_side_true, nbins_long_true);
  auto hist_min = sort_q(fBinMethod, qout_lo, qside_lo, qlong_lo, qout_lo_true, qside_lo_true, qlong_lo_true);
  auto hist_max = sort_q(fBinMethod, qout_hi, qside_hi, qlong_hi, qout_hi_true, qside_hi_true, qlong_hi_true);

  auto axis_titles = sort_q(fBinMethod, "q_{o}", "q_{s}", "q_{l}", "q_{t,o}", "q_{t,s}", "q_{t,l}");

  TString title = "Momentum Correction Hypercube Histogram;";

  for (const TString &axis_title : axis_titles) {
    title += axis_title + ";";
  }

  fHistogram = new HistType(prefix + "HyperCube", title, 6, nbins_v.data(), hist_min.data(), hist_max.data());
  UpdateQlimits();
}


AliFemtoModelCorrFctnTrueQ6D::AliFemtoModelCorrFctnTrueQ6D(const Builder &params)
  : AliFemtoModelCorrFctnTrueQ6D(params.title.Data(), params.bin_count, params.qmin, params.qmax)
{
  fIgnoreZeroMassParticles = params.ignore_zeromass;
  UpdateQlimits();
  std::array<double, 2> a = {params.qout_range_min, params.qout_range_max},
                        b = {params.qside_range_min, params.qside_range_max},
                        c = {params.qlong_range_min, params.qlong_range_max};

  SetQrange(a.data(), b.data(), c.data());
}


AliFemtoModelCorrFctnTrueQ6D::
  AliFemtoModelCorrFctnTrueQ6D(const TString &prefix,
                               const std::vector<double> &obins,
                               const std::vector<double> &sbins,
                               const std::vector<double> &lbins,
                               BinMethod binning)
  : AliFemtoCorrFctn()
  , fHistogram(nullptr)
  , fBinMethod(binning)
  , fIgnoreZeroMassParticles(true)
{
  const Int_t
    nobins = obins.size()-1,
    nsbins = sbins.size()-1,
    nlbins = lbins.size()-1;

  auto nbins_v = sort_q(fBinMethod, nobins, nsbins, nlbins,
                                    nobins, nsbins, nlbins);

  auto hist_min = sort_q(fBinMethod, obins.front(), sbins.front(), lbins.front(),
                                     obins.front(), sbins.front(), lbins.front());

  auto hist_max = sort_q(fBinMethod, obins.front(), sbins.front(), lbins.front(),
                                     obins.front(), sbins.front(), lbins.front());

  auto axis_titles = sort_q(fBinMethod, "q_{o}", "q_{s}", "q_{l}", "q_{t,o}", "q_{t,s}", "q_{t,l}");

  TString title = "Momentum Correction Hypercube Histogram;";

  for (const TString &axis_title : axis_titles) {
    title += axis_title + ";";
  }

  fHistogram = new HistType(prefix + "HyperCube", title, 6, nbins_v.data(), hist_min.data(), hist_max.data());

  auto varbins = sort_q(fBinMethod, obins, sbins, lbins, obins, sbins, lbins);

  // change axis bins after histogram created... that's ok, right?
  for (int i=0; i<6; ++i) {
    fHistogram->GetAxis(i)->Set(varbins[i].size()-1, varbins[i].data());
  }

  UpdateQlimits();
}

AliFemtoModelCorrFctnTrueQ6D::AliFemtoModelCorrFctnTrueQ6D(const AliFemtoModelCorrFctnTrueQ6D &orig)
  : AliFemtoCorrFctn(orig)
  , fHistogram(reinterpret_cast<decltype(fHistogram)>(orig.fHistogram->Clone()))
  , fBinMethod(orig.fBinMethod)
  , fIgnoreZeroMassParticles(orig.fIgnoreZeroMassParticles)
{
  fQlimits[0] = orig.fQlimits[0];
  fQlimits[1] = orig.fQlimits[1];
  fQlimits[2] = orig.fQlimits[2];
}


AliFemtoModelCorrFctnTrueQ6D::~AliFemtoModelCorrFctnTrueQ6D()
{
  delete fHistogram;
}

AliFemtoModelCorrFctnTrueQ6D&
AliFemtoModelCorrFctnTrueQ6D::operator=(const AliFemtoModelCorrFctnTrueQ6D& rhs)
{
  if (this != &rhs) {
    delete fHistogram;
    fHistogram = static_cast<HistType*>(rhs.fHistogram->Clone());
    fBinMethod = rhs.fBinMethod;
    fIgnoreZeroMassParticles = rhs.fIgnoreZeroMassParticles;
    fQlimits[0] = rhs.fQlimits[0];
    fQlimits[1] = rhs.fQlimits[1];
    fQlimits[2] = rhs.fQlimits[2];
  }
  return *this;
}


TList*
AliFemtoModelCorrFctnTrueQ6D::GetOutputList()
{
  TList *result = new TList();
  AppendOutputList(*result);
  return result;
}


TList*
AliFemtoModelCorrFctnTrueQ6D::AppendOutputList(TList &list)
{
  list.Add(fHistogram);
  return &list;
}


static
std::tuple<Double_t, Double_t, Double_t>
Qcms(const AliFemtoLorentzVector &p1, const AliFemtoLorentzVector &p2)
{
  const AliFemtoLorentzVector p = p1 + p2,
                              d = p1 - p2;

  #define FAST_DIVIDE(num, den) __builtin_expect(den == 0.0, 0) ? 0.0 : num / den

  const Double_t
    pt = p.Perp(),
    E = p.e(),
    pz = p.z(),
    dE = d.e(),
    dz = d.z(),

    qout = FAST_DIVIDE(d.x()*p.x() + d.y()*p.y(), pt),

    qside = FAST_DIVIDE(p2.x()*p1.y() - p1.x()*p2.y(), pt),

    qlong = __builtin_expect(E == pz, 0) ? 0.0 : (E*dz - pz*dE) / std::sqrt(E*E - pz*pz);

  // Double_t beta = p.z()/p.t(),
  //         gamma = 1.0 / TMath::Sqrt((1.0-beta)*(1.0+beta)),
  //         qlong = gamma * (d.z() - beta*d.t());

  #undef FAST_DIVIDE

  const double factor = std::copysign(1.0, qout);

  return std::make_tuple(factor * qout, factor * qside, factor * qlong);
}

void
AliFemtoModelCorrFctnTrueQ6D::AddPair(const AliFemtoParticle &particle1,
                                      const AliFemtoParticle &particle2)
{
  auto out_of_bounds = [] (double q, const std::pair<double,double> limit) {
    return q < limit.first || limit.second < q;
  };

  // Get generated momentum from hidden info
  const AliFemtoModelHiddenInfo
    *info1 = static_cast<const AliFemtoModelHiddenInfo*>(particle1.HiddenInfo()),
    *info2 = static_cast<const AliFemtoModelHiddenInfo*>(particle2.HiddenInfo());

  if (info1 == nullptr || info2 == nullptr) {
    return;
  }

  const Float_t mass1 = info1->GetMass(),
                mass2 = info2->GetMass();

  // block all zero-mass particles from the correlation function
  if (fIgnoreZeroMassParticles && (mass1 == 0.0 || mass2 == 0.0)) {
    return;
  }

  const AliFemtoThreeVector &true_momentum1 = *info1->GetTrueMomentum(),
                            &true_momentum2 = *info2->GetTrueMomentum();

  const Double_t e1 = std::sqrt(mass1 * mass1 + true_momentum1.Mag2()),
                 e2 = std::sqrt(mass2 * mass2 + true_momentum2.Mag2());

  const AliFemtoLorentzVector p1(e1, true_momentum1),
                              p2(e2, true_momentum2);

  Double_t true_q_out, true_q_side, true_q_long;
  std::tie(true_q_out, true_q_side, true_q_long) = Qcms(p1, p2);

  if (out_of_bounds(true_q_out, fQlimits[0]) ||
      out_of_bounds(true_q_side, fQlimits[1]) ||
      out_of_bounds(true_q_long, fQlimits[2])) {
    return;
  }

  Double_t q_out, q_side, q_long;
  std::tie(q_out, q_side, q_long) = Qcms(particle1.FourMomentum(), particle2.FourMomentum());

  auto q = sort_q(fBinMethod,
                  q_out, q_side, q_long,
                  true_q_out, true_q_side, true_q_long);

  fHistogram->Fill(q.data());
}

void
AliFemtoModelCorrFctnTrueQ6D::AddRealPair(AliFemtoPair *pair)
{
}

void
AliFemtoModelCorrFctnTrueQ6D::AddMixedPair(AliFemtoPair *pair)
{
  const AliFemtoParticle *p1 = pair->Track1(),
                         *p2 = pair->Track2();

  AddPair(*p1, *p2);
}

AliFemtoString
AliFemtoModelCorrFctnTrueQ6D::Report()
{
  TString report;
  return AliFemtoString(report.Data());
}

void
AliFemtoModelCorrFctnTrueQ6D::SetQrange(const double o[2], const double s[2], const double l[2])
{
  fQlimits[0] = {o[0], o[1]};
  fQlimits[1] = {s[0], s[1]};
  fQlimits[2] = {l[0], l[1]};
}

AliFemtoModelCorrFctnTrueQ6D::Builder::Builder(AliFemtoModelCorrFctnTrueQ6D::Builder const &orig)
  : bin_count(orig.bin_count)
  , qmin(orig.qmin)
  , qmax(orig.qmax)
  , bin_method(orig.bin_method)
  , title(orig.title)
  , qout_range_min(orig.qout_range_min)
  , qout_range_max(orig.qout_range_max)
  , qside_range_min(orig.qside_range_min)
  , qside_range_max(orig.qside_range_max)
  , qlong_range_min(orig.qlong_range_min)
  , qlong_range_max(orig.qlong_range_max)
  , ignore_zeromass(orig.ignore_zeromass)
{
}


AliFemtoModelCorrFctnTrueQ6D::Builder&
AliFemtoModelCorrFctnTrueQ6D::Builder::operator=(AliFemtoModelCorrFctnTrueQ6D::Builder const &rhs)
{
  if (this != &rhs) {
    bin_count = rhs.bin_count;
    qmin = rhs.qmin;
    qmax = rhs.qmax;
    bin_method = rhs.bin_method;
    title = rhs.title;
    qout_range_min = rhs.qout_range_min;
    qout_range_max = rhs.qout_range_max;
    qside_range_min = rhs.qside_range_min;
    qside_range_max = rhs.qside_range_max;
    qlong_range_min = rhs.qlong_range_min;
    qlong_range_max = rhs.qlong_range_max;
    ignore_zeromass = rhs.ignore_zeromass;
  }
  return *this;
}
