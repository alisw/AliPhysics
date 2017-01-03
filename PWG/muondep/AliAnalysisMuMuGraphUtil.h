#ifndef ALIANALYSISMUMUGRAPHUTIL_H
#define ALIANALYSISMUMUGRAPHUTIL_H

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

#include <vector>
#include <string>
#include <set>

class TGraph;
class TGraphErrors;
class TCanvas;
class TObjArray;

#include "TAttFill.h"
#include "TAttMarker.h"
#include "TAttLine.h"
#include "TAttAxis.h"
#include "TString.h"

/**

  @ingroup pwg_muondep_mumu

  @class AliAnalysisMuMuGraphUtil

  @brief Utilities for graph drawing

  @author Laurent Aphecetche (Subatech)
*/

class AliAnalysisMuMuGraphUtil : public TObject
{
public:
  
  AliAnalysisMuMuGraphUtil(const char* ocdbpath="raw://");
  virtual ~AliAnalysisMuMuGraphUtil() {}

  static TGraphErrors* Combine(TObjArray& graph, Bool_t compact);
  
  static void Compact(TGraph& g);
  
  void DefaultStyle();

  TCanvas* DrawWith2Scales(TGraph& g1, TGraph& g2, const char* canvasName="c1");

  static Int_t GetRunNumber(const TGraph& g, Int_t i);

  void GetRuns(std::set<int>& runs, TGraph& graph) const;

  static Bool_t IsCompact(TGraph& g);
  
  void PlotSameWithLegend(TObjArray& a, Double_t ymin, Double_t ymax) const;

  void ShouldDrawPeriods(Bool_t value) { fShouldDrawPeriods = value; }

  void StyleGraph(TGraph& graph, UInt_t index) const;
  
  static void UnCompact(TGraph& g);
  
  static void GetYMinAndMax(TGraph& g, Double_t& ymin, Double_t& ymax);
  
  static TGraph* RelDif(TGraph& ga, TGraph& gb);

  
private:
  AliAnalysisMuMuGraphUtil(const AliAnalysisMuMuGraphUtil& rhs); // not implemented
  AliAnalysisMuMuGraphUtil& operator=(const AliAnalysisMuMuGraphUtil& rhs); // not implemented

  TString fOCDBPath; // OCDB path

  std::vector<TAttLine> fAttLine; // line attributes
  std::vector<TAttMarker> fAttMarker; // marker attributes
  std::vector<TAttFill> fAttFill; // fill attributes
  std::vector<TAttAxis> fAttXaxis; // x-axis attributes
  std::vector<TAttAxis> fAttYaxis; // y-axis attributes
  std::vector<std::string> fDrawOptions; // draw options
  
  Bool_t fShouldDrawPeriods; // draw period names on top of graphs
  
  ClassDef(AliAnalysisMuMuGraphUtil,0) // utility class to modify/plot graphs
};

#endif
