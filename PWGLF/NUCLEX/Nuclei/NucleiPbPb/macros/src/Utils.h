#ifndef UTILS_H
#define UTILS_H

#include <cmath>
#include <limits.h>
#include <memory>
#include <vector>

#include <TColor.h>
#include <TF1.h>
#include <TGraph.h>
#include <TH1.h>
#include <TList.h>
#include <TObject.h>

namespace utils {

  template<class T> void Requires(T* obj, string msg = "") {
    if (!obj) {
      std::cout << "Missing object: " << msg.data() << "." << std::endl;
      abort();
    }
  }

  class TTList : public TList {
    public:
      TObject* Get(std::string name) {
        TObject* obj = this->FindObject(name.data());
        Requires(obj,name);
        return obj;
      }
  };

  /// Older version of ROOT 6 do not provide TMath::Sq
  template<typename T> constexpr T Sq(T x) {
    return ((x) * (x));
  }

  constexpr bool ValueInInterval(const double val, const double min, const double max) {
    return (val >= min && val < max);
  }

  double QuadratureSum(std::vector<double> &args) {
    double sum = 0.;
    for (double& val : args)
      sum += Sq(val);
    return std::sqrt(sum);
  }

  /// Workaround if C++14 is not available
  template<typename T, typename... Args>
  std::unique_ptr<T> make_unique(Args&&... args) {
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
  }

  double Ztest(const double mu0, const double sig0, const double mu1, const double sig1, const double corr = 0.) {
    const double sigma = std::sqrt(sig0 * sig0 + sig1 * sig1 + 2 * corr);
    if (sigma < FLT_MIN * 10.f) return FLT_MAX;
    else return (mu0 - mu1) / sigma;
  }

  void Divide(TH1* h, TGraph* gr) {
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
      if (h->GetBinContent(i) == 0 || std::abs(gr->Eval(h->GetBinCenter(i))) < FLT_MIN * 10.f) {
        continue;
      }
      h->SetBinContent(i, h->GetBinContent(i) / gr->Eval(h->GetBinCenter(i)));
      h->SetBinError(i, h->GetBinError(i) / gr->Eval(h->GetBinCenter(i)));
    }
  }

  TH1* ComputeEfficiency(TH1* tof, TH1* tot) {
    TH1* efftof = (TH1*)tof->Clone();
    for (int iBin = 1; iBin <= efftof->GetNbinsX(); ++iBin) {
      const double den = tot->GetBinContent(iBin);
      if (std::abs(den) < FLT_MIN * 10.f) continue;
      double eff = tof->GetBinContent(iBin) / den;
      efftof->SetBinContent(iBin, eff);
      efftof->SetBinError(iBin, std::sqrt(eff * (1. - eff) / den));
    }
    return efftof;
  }
}

#endif
