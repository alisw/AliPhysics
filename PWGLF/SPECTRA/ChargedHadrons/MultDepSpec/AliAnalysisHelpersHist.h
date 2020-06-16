#ifndef AliAnalysisHelpersHist_cxx
#define AliAnalysisHelpersHist_cxx

#include <iostream>

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THn.h"
#include "THnSparse.h"
#include "THnBase.h"
#include "TAxis.h"

namespace Hist
{

typedef struct{
  std::string name;
  std::string title;
  std::vector<double> binEdges;
  int nBins; // 0 when bin edges are specified directly
} Axis;


template<typename RootHist_t>
class Hist
{
public:
  Hist() : fAxes{}, fRawHist{nullptr}{}
  Hist(const Hist&) = delete; // non construction-copyable
  Hist& operator=(const Hist&) = delete; // non copyable

  void AddAxis(const Axis& axis)
  {
    fAxes.push_back(axis);
  }
  void AddAxes(const std::vector<Axis>& axes)
  {
    fAxes.insert(fAxes.end(), axes.begin(), axes.end());
  }
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
    return new RootHist_t(name.data(), title.data(), nDim, nBins, lowerBounds, upperBounds);
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
    const std::size_t MAX_DIM = 10;
    const std::size_t nAxes = fAxes.size();
    if(nAxes == 0 || nAxes > MAX_DIM) return nullptr;

    int nBins[MAX_DIM] = {0};
    double lowerBounds[MAX_DIM] = {0.0};
    double upperBounds[MAX_DIM] = {0.0};

    // first figure out number of bins and dimensions
    std::string title = "[ ";
    for(std::size_t i = 0; i < nAxes; i++){
      nBins[i] = (fAxes[i].nBins) ? fAxes[i].nBins : fAxes[i].binEdges.size()-1;
      lowerBounds[i] = fAxes[i].binEdges[0];
      upperBounds[i] = fAxes[i].binEdges[fAxes[i].binEdges.size()-1];
      title += fAxes[i].name;
      if(i < nAxes-1) title += " : "; else title += " ]";
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
        axis->SetTitle(fAxes[i].title.data());
        if(std::is_base_of<THnBase, RootHist_t>::value) axis->SetName((std::to_string(i) + "-" + fAxes[i].name).data());
        
        // move the bin edges in case a variable binnining was requested
        if(!fAxes[i].nBins)
        {
          if(!std::is_sorted(std::begin(fAxes[i].binEdges), std::end(fAxes[i].binEdges)))
          {
            std::cout << "ERROR: The bin edges specified for axis " << fAxes[i].name << " in histogram " << name << " are not in increasing order!" << std::endl;
            return nullptr;
          }
          axis->Set(nBins[i], fAxes[i].binEdges.data());
        }
      }
    }
    if(hasWeights) fRawHist->Sumw2();
    fAxes.clear();
    return fRawHist;
  }


  // fill functions

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

  template<typename T = RootHist_t, typename... Ts,
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


  // size functions

  template<typename T = RootHist_t, typename std::enable_if<std::is_base_of<TH1, T>::value>::type* dummy = nullptr>
  double GetSize()
  {
    return (!fRawHist) ? 0. : fRawHist->GetSize() * (sizeof(fRawHist->At(0)) + ((fRawHist->GetSumw2()->fN) ? sizeof(double) : 0.));
  }

  template<typename B>
  int GetBaseElementSize(THnT<B>* ptr)
  {
      return sizeof(B);
  };
  template<typename T = RootHist_t, typename std::enable_if<std::is_base_of<THn, T>::value>::type* dummy = nullptr>
  double GetSize()
  {
    return (!fRawHist) ? 0. : fRawHist->GetNbins() * (GetBaseElementSize(fRawHist) + ((fRawHist->GetSumw2() != -1.) ? sizeof(double) : 0.));
  }

  template<typename B>
  int GetBaseElementSize(THnSparseT<B>* ptr)
  {
      return sizeof(B);
  };
  template<typename T = RootHist_t, typename std::enable_if<std::is_base_of<THnSparse, T>::value>::type* dummy = nullptr>
  double GetSize(double fillFraction = 1.)
  {
    if(!fRawHist) return 0.;
    double nbinsTotal = 1.;
    for (Int_t d = 0; d < fRawHist->GetNdimensions(); ++d)
       nbinsTotal *= fRawHist->GetAxis(d)->GetNbins() + 2;

    Double_t overhead = 4.; // probably often less; unfortunatley cannot access fRawHist->GetCompactCoord()->GetBufferSize();

    return fillFraction * nbinsTotal * (GetBaseElementSize(fRawHist) + overhead + ((fRawHist->GetSumw2() != -1.) ? sizeof(double) : 0.));
  }

private:
  std::vector<Axis> fAxes;
  RootHist_t* fRawHist;
};

} // end namespace Hist

#endif
