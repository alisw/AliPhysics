///
/// \file AliFemtoModelCorrFctnTrueQ6D.cxx
///

#include "AliFemtoModelCorrFctnTrueQ6D.h"
#include "AliFemtoPair.h"
#include "AliFemtoModelManager.h"
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
  const auto aout = fHistogram->GetAxis(3),
             aside = fHistogram->GetAxis(4),
             along = fHistogram->GetAxis(5);

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

AliFemtoModelCorrFctnTrueQ6D::AliFemtoModelCorrFctnTrueQ6D(const HistType &hist, AliFemtoModelManager *mgr)
  : AliFemtoCorrFctn()
  , fManager(mgr)
  , fHistogram(reinterpret_cast<HistType*>(hist.Clone()))
  , fRng(new TRandom())
  , fBinMethod(kRecGenOSL)
  , fIgnoreZeroMassParticles(true)
  , fUseFemtoWeight(true)
{
  UpdateQlimits();
  fBinMethod = GuessBinMethod(*fHistogram);
}

AliFemtoModelCorrFctnTrueQ6D::AliFemtoModelCorrFctnTrueQ6D(HistType *&hist, AliFemtoModelManager *mgr)
  : AliFemtoCorrFctn()
  , fManager(mgr)
  , fHistogram(hist)
  , fRng(new TRandom())
  , fBinMethod(kRecGenOSL)
  , fIgnoreZeroMassParticles(true)
  , fUseFemtoWeight(true)
{
  // take ownership of the pointer
  hist = nullptr;
  UpdateQlimits();
  fBinMethod = GuessBinMethod(*fHistogram);
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

AliFemtoModelCorrFctnTrueQ6D::AliFemtoModelCorrFctnTrueQ6D(const TString &name,
                                                           UInt_t nbins,
                                                           Double_t qmin,
                                                           Double_t qmax)
  : AliFemtoCorrFctn()
  , fManager(nullptr)
  , fHistogram(nullptr)
  , fRng(new TRandom())
  , fBinMethod(kRecGenOSL)
  , fIgnoreZeroMassParticles(true)
  , fUseFemtoWeight(true)
{
  std::vector<Int_t> nbins_v(6, static_cast<Int_t>(nbins));
  std::vector<double> hist_min(6, qmin), hist_max(6, qmax);

  auto axis_titles = sort_q(fBinMethod,
      "q_{o}", "q_{s}", "q_{l}",
      "q_{t,o}", "q_{t,s}", "q_{t,l}");

  TString title = "Momentum Correction Hypercube Histogram ";

  for (const TString &axis_title : axis_titles) {
    title += "; " + axis_title;
  }

  fHistogram = new HistType(name +"_Histogram", title, 6, nbins_v.data(), hist_min.data(), hist_max.data());
  UpdateQlimits();
}


AliFemtoModelCorrFctnTrueQ6D::AliFemtoModelCorrFctnTrueQ6D(const Builder &params)
  : AliFemtoModelCorrFctnTrueQ6D(params.title.Data(), params.bin_count, params.qmin, params.qmax)
{
  fIgnoreZeroMassParticles = params.ignore_zeromass;
  SetManager(params.mc_manager);
  UpdateQlimits();
  std::array<double, 2> a = {params.qout_range_min, params.qout_range_max},
                        b = {params.qside_range_min, params.qside_range_max},
                        c = {params.qlong_range_min, params.qlong_range_max};

  SetQrange(a.data(), b.data(), c.data());
}


AliFemtoModelCorrFctnTrueQ6D::AliFemtoModelCorrFctnTrueQ6D(const AliFemtoModelCorrFctnTrueQ6D &orig)
  : AliFemtoCorrFctn(orig)
  , fManager(orig.fManager)
  , fHistogram(reinterpret_cast<decltype(fHistogram)>(orig.fHistogram->Clone()))
  , fRng(new TRandom())
  , fBinMethod(orig.fBinMethod)
  , fIgnoreZeroMassParticles(orig.fIgnoreZeroMassParticles)
  , fUseFemtoWeight(orig.fUseFemtoWeight)
{
  fQlimits[0] = orig.fQlimits[0];
  fQlimits[1] = orig.fQlimits[1];
  fQlimits[2] = orig.fQlimits[2];
}


AliFemtoModelCorrFctnTrueQ6D::~AliFemtoModelCorrFctnTrueQ6D()
{
  delete fHistogram;
  delete fRng;
}

AliFemtoModelCorrFctnTrueQ6D&
AliFemtoModelCorrFctnTrueQ6D::operator=(const AliFemtoModelCorrFctnTrueQ6D& rhs)
{
  if (this != &rhs) {
    delete fHistogram;
    fHistogram = static_cast<HistType*>(rhs.fHistogram->Clone());
    fBinMethod = rhs.fBinMethod;
    fIgnoreZeroMassParticles = rhs.fIgnoreZeroMassParticles;
    fUseFemtoWeight = rhs.fUseFemtoWeight;
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

  Double_t k1 = p.Perp(),
           k2 = d.x()*p.x() + d.y()*p.y();

  Double_t qout = (k1 == 0) ? 0.0 : k2/k1;

  Double_t qside = (k1 == 0) ? 0.0 : 2.0 * (p2.x()*p1.y() - p1.x()*p2.y())/k1;

  // Double_t beta = p.z()/p.t(),
  //         gamma = 1.0 / TMath::Sqrt((1.0-beta)*(1.0+beta)),
  //         qlong = gamma * (d.z() - beta*d.t());

  Double_t pt = p.t(),
           pz = p.z(),
           dt = d.t(),
           dz = d.z(),
           qlong = (pt*dz - pz*dt) / TMath::Sqrt(pt*pt - pz*pz);

  return std::make_tuple(qout, qside, qlong);
}

void
AliFemtoModelCorrFctnTrueQ6D::AddPair(const AliFemtoParticle &particle1,
                                      const AliFemtoParticle &particle2,
                                      double femto_weight /* = 1.0 */)
{
  // Fill reconstructed histogram with "standard" particle momentum
  Double_t q_out, q_side, q_long;
  std::tie(q_out, q_side, q_long) = Qcms(particle1.FourMomentum(), particle2.FourMomentum());

  auto out_of_bounds = [] (double q, const std::pair<double,double> limit) {
    return q < limit.first || limit.second < q;
  };

  if (out_of_bounds(q_out, fQlimits[0]) ||
      out_of_bounds(q_side, fQlimits[1]) ||
      out_of_bounds(q_long, fQlimits[2])) {
    return;
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
  if (fIgnoreZeroMassParticles && (mass1 == 0.0 || mass2 == 0.0)) {
    return;
  }

  const AliFemtoThreeVector &true_momentum1 = *info1->GetTrueMomentum(),
                            &true_momentum2 = *info2->GetTrueMomentum();

  const Double_t e1 = sqrt(mass1 * mass1 + true_momentum1.Mag2()),
                 e2 = sqrt(mass2 * mass2 + true_momentum2.Mag2());

  const AliFemtoLorentzVector p1(e1, true_momentum1),
                              p2(e2, true_momentum2);

  // Fill generated-momentum histogram with "true" particle momentum
  Double_t true_q_out, true_q_side, true_q_long;
  std::tie(true_q_out, true_q_side, true_q_long) = Qcms(p1, p2);

  auto q = sort_q(fBinMethod,
                  q_out, q_side, q_long,
                  true_q_out, true_q_side, true_q_long);

  fHistogram->Fill(q.data(), femto_weight);
}

void
AliFemtoModelCorrFctnTrueQ6D::AddRealPair(AliFemtoPair *pair)
{
  const AliFemtoParticle *p1 = pair->Track1(),
                         *p2 = pair->Track2();

  // randomize to avoid ordering biases
  if (fRng->Uniform() >= 0.5) {
    std::swap(p1, p2);
  }

  Double_t femto_weight = fUseFemtoWeight ? fManager->GetWeight(pair) : 1.0;
  AddPair(*p1, *p2, femto_weight);
}

void
AliFemtoModelCorrFctnTrueQ6D::AddMixedPair(AliFemtoPair *pair)
{
  const AliFemtoParticle *p1 = pair->Track1(),
                         *p2 = pair->Track2();

  // randomize to avoid ordering biases
  if (fRng->Uniform() >= 0.5) {
    std::swap(p1, p2);
  }

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
  , mc_manager(orig.mc_manager)
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
    mc_manager = rhs.mc_manager;
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
