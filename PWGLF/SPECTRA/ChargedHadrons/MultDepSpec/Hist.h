#ifndef Hist_cxx
#define Hist_cxx

#include <iostream>

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THnBase.h"
#include "TAxis.h"

namespace AnalysisHelpers
{

template<typename RootHist_t>
class Hist
{
public:
  Hist() : fAxes{}, fRawHist{nullptr}{}
  Hist(const Hist&) = delete; // non construction-copyable
  Hist& operator=(const Hist&) = delete; // non copyable
  
  void AddAxis(const std::string& name, const std::string& title, const int& nBins, const double& lowerEdge, const double& upperEdge)
  {
    fAxes.push_back({name, title, {lowerEdge, upperEdge}, nBins});
  }
  void AddAxis(const std::string& name, const std::string& title, const std::vector<double>& binEdges)
  {
    fAxes.push_back({name, title, binEdges, 0});
  }

  // use SFINAE instead of 'if constexpr' since we are restricted to cpp11
  
  template<typename T = RootHist_t, typename std::enable_if<std::is_base_of<TH3, T>::value>::type* dummy = nullptr>
  RootHist_t* HistFactory (const std::string& name, const std::string& title, const int& nDim, const int nBins[], const double lowerBounds[], const double upperBounds[])
  {
    return (nDim != 3) ? nullptr : new RootHist_t(name.data(), title.data(), nBins[0], lowerBounds[0], upperBounds[0], nBins[1], lowerBounds[1], upperBounds[1], nBins[2], lowerBounds[2], upperBounds[2]);
  }
  template<typename T = RootHist_t, typename std::enable_if<std::is_base_of<TH2, T>::value>::type* dummy = nullptr>
  RootHist_t* HistFactory (const std::string& name, const std::string& title, const int& nDim, const int nBins[], const double lowerBounds[], const double upperBounds[])
  {
    return (nDim != 2) ? nullptr : new RootHist_t(name.data(), title.data(), nBins[0], lowerBounds[0], upperBounds[0], nBins[1], lowerBounds[1], upperBounds[1]);
  }
  template<typename T = RootHist_t, typename std::enable_if<!std::is_base_of<TH3, T>::value && !std::is_base_of<TH2, T>::value && std::is_base_of<TH1, T>::value>::type* dummy = nullptr>
  RootHist_t* HistFactory (const std::string& name, const std::string& title, const int& nDim, const int nBins[], const double lowerBounds[], const double upperBounds[])
  {
    return (nDim != 1) ? nullptr : new RootHist_t(name.data(), title.data(), nBins[0], lowerBounds[0], upperBounds[0]);
  }
  template<typename T = RootHist_t, typename std::enable_if<std::is_base_of<THnBase, T>::value>::type* dummy = nullptr>
  RootHist_t* HistFactory (const std::string& name, const std::string& title, const int& nDim, const int nBins[], const double lowerBounds[], const double upperBounds[])
  {
    return new RootHist_t(name.c_str(), title.c_str(), nDim, nBins, lowerBounds, upperBounds);
  }

  template<typename T = RootHist_t, typename std::enable_if<std::is_base_of<THnBase, T>::value>::type* dummy = nullptr>
  TAxis* GetAxis (const int& i)
  {
    return fRawHist->GetAxis(i);
  }
  template<typename T = RootHist_t, typename std::enable_if<!std::is_base_of<THnBase, T>::value>::type* dummy = nullptr>
  TAxis* GetAxis (const int& i)
  {
    return (i == 0) ? fRawHist->GetXaxis() : (i == 1) ? fRawHist->GetYaxis() : (i == 2) ? fRawHist->GetZaxis() : nullptr;
  }

  RootHist_t* GenerateHist(std::string name, bool hasWeights = false)
  {
    const std::size_t MAX_DIM = 6;
    const std::size_t nAxes = fAxes.size();
    if(nAxes == 0 || nAxes > MAX_DIM) return nullptr;

    int nBins[MAX_DIM] = {0};
    double lowerBounds[MAX_DIM] = {0.0};
    double upperBounds[MAX_DIM] = {0.0};

    // first figure out number of bins and dimensions
    std::string title = name + " [";
    for(std::size_t i = 0; i < nAxes; i++){
      nBins[i] = (fAxes[i].nBins) ? fAxes[i].nBins : fAxes[i].binEdges.size()-1;
      lowerBounds[i] = fAxes[i].binEdges[0];
      upperBounds[i] = fAxes[i].binEdges[fAxes[i].binEdges.size()-1];
      title += fAxes[i].name;
      if(i < nAxes-1) title += " : "; else title += "]";
    }
    
    // create histogram
    if(fRawHist) delete fRawHist;
    fRawHist = HistFactory(name, title, nAxes, nBins, lowerBounds, upperBounds);

    if(!fRawHist)
    {
      std::cout << "ERROR: The number of specified dimensions does not match the type." << std::endl;
      return nullptr;
    }

    // set axis properties
    for(std::size_t i = 0; i < nAxes; i++)
    {
      TAxis* axis = GetAxis(i);
      if(axis)
      {
        axis->SetName((std::to_string(i) + "-" + fAxes[i].name).c_str());
        axis->SetTitle(fAxes[i].title.c_str());
        // move bin edges in case variable binnining was requested
        if(!fAxes[i].nBins) axis->Set(nBins[i], fAxes[i].binEdges.data());
      }
    }
    if(hasWeights) fRawHist->Sumw2();
    fAxes.clear();
    return fRawHist;
  }
  
  
  template<typename T = RootHist_t, typename... Ts,
  typename std::enable_if<std::is_base_of<THnBase, T>::value>::type* dummy = nullptr>
  inline void Fill (const Ts&... position)
  {
    //if(fRawHist->GetNdimensions() != sizeof...(position)) return;
    double tempArray[] = {static_cast<double>(position)...};
    fRawHist->Fill(tempArray);
  }

  template<typename T = RootHist_t, typename... Ts, typename std::enable_if<
  (
  ((std::is_base_of<TH3, T>::value) && sizeof...(Ts) == 3) ||
  ((std::is_base_of<TH2, T>::value) && sizeof...(Ts) == 2) ||
  ((std::is_base_of<TH1, T>::value) && sizeof...(Ts) == 1)
  )
  >::type* dummy = nullptr>
  inline void Fill (const Ts&... position)
  {
    fRawHist->Fill(static_cast<double>(position)...);
  }
  
  template<typename T = RootHist_t, typename... Ts,
  typename std::enable_if<std::is_base_of<THnBase, T>::value>::type* dummy = nullptr>
  inline void FillWeight (const T& weight, const Ts&... position)
  {
    //if(fRawHist->GetNdimensions() != sizeof...(position)) return;
    double tempArray[] = {static_cast<double>(position)...};
    fRawHist->Fill(tempArray, static_cast<double>(weight));
  }

  template<typename T = RootHist_t, typename... Ts, \
  typename std::enable_if<
  (
  ((std::is_base_of<TH3, T>::value) && sizeof...(Ts) == 3) ||
  ((std::is_base_of<TH2, T>::value) && sizeof...(Ts) == 2) ||
  ((std::is_base_of<TH1, T>::value) && sizeof...(Ts) == 1)
  )
  >::type* dummy = nullptr>
  inline void FillWeight (const T& weight, const Ts&... position)
  {
    fRawHist->Fill(static_cast<double>(position)..., static_cast<double>(weight));
  }

  
private:
  typedef struct{
    std::string name;
    std::string title;
    std::vector<double> binEdges;
    int nBins; // 0 when bin edges are specified directly
  } Axis_t;

  std::vector<Axis_t> fAxes;
  RootHist_t* fRawHist;
};

} // end namespace AnalysisHelpers

#endif
