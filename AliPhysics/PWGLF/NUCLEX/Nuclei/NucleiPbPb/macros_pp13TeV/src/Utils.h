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

  bool replace(std::string& str, const std::string& from, const std::string& to) {
    size_t start_pos = str.find(from);
    if(start_pos == std::string::npos)
      return false;
    str.replace(start_pos, from.length(), to);
    return true;
  }

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

  void MeanRebin(TH1F*& hOriginal, int n_bins, const double* bin_limits, float Max_val = 100){
    int iBinCounter = 1;
    std::vector<float> content_tmp(n_bins,0.);
    std::vector<float> err2_tmp(n_bins,0.);
    std::vector<int> count_tmp(n_bins,0);
    TH1F* hTmp = (TH1F*)hOriginal->Clone("hTmp");
    const char* original_name = hOriginal->GetName();
    printf("original name: %s\n",original_name );
    const char* original_title = hOriginal->GetTitle();
    printf("original title: %s\n",original_title );
    const char* original_X = hOriginal->GetXaxis()->GetTitle();
    printf("original X: %s\n",original_X );
    const char* original_Y = hOriginal->GetYaxis()->GetTitle();
    printf("original Y: %s\n",original_Y );
    delete hOriginal;
    for(int iB=1; iB<=hTmp->GetNbinsX(); iB++){
      if(hTmp->GetBinCenter(iB) < bin_limits[0]) continue;
      if(hTmp->GetBinCenter(iB) > Max_val) break;
      if(hTmp->GetBinCenter(iB) < bin_limits[iBinCounter]){
        content_tmp[iBinCounter-1] += hTmp->GetBinContent(iB);
        err2_tmp[iBinCounter-1] += Sq(hTmp->GetBinError(iB));
        count_tmp[iBinCounter-1]++;
      }
      else{
        iBinCounter++;
        content_tmp[iBinCounter-1] += hTmp->GetBinContent(iB);
        err2_tmp[iBinCounter-1] += Sq(hTmp->GetBinError(iB));
        count_tmp[iBinCounter-1]++;
      }
    }
    hOriginal = new TH1F(original_name,Form("%s;%s;%s",original_title,original_X,original_Y),n_bins,bin_limits);
    hOriginal->SetDirectory(0);
    for(int i=0; i<n_bins; i++){
      if(hOriginal->GetBinCenter(i+1)>Max_val) continue;
      hOriginal->SetBinContent(i+1,content_tmp[i]/count_tmp[i]);
      hOriginal->SetBinError(i+1,TMath::Sqrt(err2_tmp[i])/count_tmp[i]);
    }
  }

  float zTest(const float mu0, const float sig0, const float mu1, const float sig1) {
    const float sigma = sqrt(sig0 * sig0 + sig1 * sig1);
    if (sigma < FLT_MIN * 10.f) return FLT_MAX;
    else return (mu0 - mu1) / std::abs(sig1-sig0);
  }

  void SmoothInRange(TH1F* histo, float low_limit, float up_limit, int n_times = 1){
    histo->GetXaxis()->SetRange(histo->FindBin(low_limit),histo->FindBin(up_limit));
    histo->Smooth(n_times,"R");
  }
}

#endif
