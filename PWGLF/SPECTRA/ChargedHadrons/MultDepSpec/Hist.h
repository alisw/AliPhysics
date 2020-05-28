#ifndef Hist_cxx
#define Hist_cxx

#include <iostream>
#include "TAxis.h"

namespace AnalysisHelpers
{

template<typename RootNdHist_t>
class Hist
{
public:
  Hist() : fAxes{}, fRawHist{nullptr}{}
  void AddAxis(std::string name, std::string title, const std::vector<double>& binEdges)
  {
    fAxes.push_back({name, title, binEdges});
  }
  void AddAxis(std::string name, std::string title, int nBins, double lowerEdge, double upperEdge)
  {
    std::vector<double> binEdges;
    for(int i = 0; i <= nBins; i++)
    {
      binEdges.push_back(lowerEdge + i*(upperEdge - lowerEdge)/nBins);
    }
    fAxes.push_back({name, title, binEdges});
  }
  RootNdHist_t* GenerateHist(std::string name, bool hasWeights = false)
  {
    const int MAX_DIM = 6;
    std::size_t nAxes = fAxes.size();
    if(nAxes == 0 || nAxes > MAX_DIM) return nullptr;

    int nBins[MAX_DIM] = {0};
    double lowerBounds[MAX_DIM] = {0.0};
    double upperBounds[MAX_DIM] = {0.0};

    std::string title = name + " [";
    // first figure out number of bins and dimensions
    for(int i = 0; i < nAxes; i++){
      nBins[i] = fAxes[i].binEdges.size()-1;
      lowerBounds[i] = fAxes[i].binEdges[0];
      upperBounds[i] = fAxes[i].binEdges[nBins[i]];
      title += fAxes[i].title;
      if(i < nAxes-1) title += " : "; else title += "]";
    }
    // create histogram
    if(fRawHist) delete fRawHist;
    fRawHist = new RootNdHist_t(name.c_str(), title.c_str(), nAxes, nBins, lowerBounds, upperBounds);
    if(!fRawHist) return nullptr;

    // set histogram axes
    for(int i = 0; i < nAxes; i++){
      fRawHist->SetBinEdges(i, fAxes[i].binEdges.data());
      fRawHist->GetAxis(i)->SetName((std::to_string(i) + "-" + fAxes[i].name).c_str());
      fRawHist->GetAxis(i)->SetTitle(fAxes[i].title.c_str());
    }
    if(hasWeights) fRawHist->Sumw2();

    fAxes.clear();
    return fRawHist;
  }

  template<std::size_t N>
  inline void Fill(const double (&values)[N]) {fRawHist->Fill(values);}
  
private:
  typedef struct{
    std::string name;
    std::string title;
    std::vector<double> binEdges;
  } Axis_t;
  std::vector<Axis_t> fAxes;
  RootNdHist_t* fRawHist;
};

} // end namespace AnalysisHelpers

#endif
