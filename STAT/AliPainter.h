#ifndef ALIPAINTER_H
#define ALIPAINTER_H

/*
 gSystem->AddIncludePath("-I$AliRoot_SRC/STAT");
.L $AliRoot_SRC/STAT/AliPainter.cxx+

TFile::SetCacheFileDir(".");
TFile *finput = TFile::Open("http://aliqatrkeos.web.cern.ch/aliqatrkeos/performance/AliPainterTest.root","CACHEREAD");
TTree *tree = (TTree *) finput->Get("hisPtAll");
hisArray = new TObjArray();
TList *keys = finput->GetListOfKeys();
for (Int_t iKey = 0; iKey<keys->GetEntries();iKey++) {
TObject *o = finput->Get(TString::Format("%s;%d", keys->At(iKey)->GetName(), ((TKey *) keys->At(iKey))->GetCycle()).Data());
hisArray->AddLast(o);
}
TCanvas *canvasQA1 = new TCanvas("canvasQA", "canvasQA", 1200, 800);
AliPainter::DivideTPad(canvasQA1, "<horizontal>[1]", "Canvas1");
canvasQA1->cd(1);
THnBase *hisN = hisArray->FindObject("hisK0DMassQPtTgl");
AliPainter::DrawHistogram(hisN, "hisName()(1)(name=gaus,option=W)(class=Raw,drawOpt=E)")

*/
/// \ingroup STAT
/// \class AliPainter
/*!
* \brief Class for generating QA reports
*  See the documentation in describing of functions.
* \author  <a href="marian.ivanov@cern.ch">Marian  Ivanov</a> , <a href="boris.rumyantsev@cern.ch">Boris Rumyantsev</a>
*/

#include <vector>
#include <map>
#include "TObjArray.h"
#include "TString.h"
#include "TPad.h"
#include "THnBase.h"
#include "TMultiGraph.h"
#include "TFormula.h"
#include "TObject.h"
#include "TLegend.h"

class AliPainter : public TObject {
  public:
    static TPad *DivideTPad(TPad *pad, const char *division, const char *classID="", const char *style="", Int_t verbose=0);
    static void SetMultiGraphTimeAxis(TMultiGraph *graph, TString option);
    static TPad *SetPadMargin(TPad *cPad, const char *position, const char *wMargin, const char *units, Double_t mValue, Int_t iCol, Int_t nCols);
    static TObjArray *PrepareHistogram(THnBase *hisN, const char *expression, TObjArray *&keepArray, TObjArray *metaData=nullptr, Int_t verbose=0);
    static void DrawHistogram(THnBase *hisN, const char *expression, TPad *pad=nullptr, TObjArray *keepArray=nullptr, TObjArray *metaData=nullptr, Int_t verbose=0);
    static void DrawHistogram(const TObjArray *histogramArray, const char *expression, TPad *pad=nullptr, TObjArray *keepArray=nullptr, TObjArray *metaData=nullptr, Int_t verbose=0);
    static TPad *GetNextPad(TPad *cPad, TPad *tempPad=nullptr, Int_t verbose=0);
    typedef std::map<Int_t, std::vector<TString> > axisRangesMap;
    static std::map<TString, TString> drawValues;
    static std::map<TString, TString> fitValues;
    static std::map<TString, TString> genValues;
    static std::vector<TString> rangesVec;
    static void FillAll(const char *, Int_t=0);
    static void ParseRanges(const TString, Int_t=0);
    static void SaveToKeepArray(TObject *, TObjArray *&, Int_t=0);
    static void SaveToKeepArray(TObjArray *, TObjArray *&, Int_t=0);
    static TObjArray *SetRanges(THnBase *, Int_t=0);
    static TObject *SetProjections(THnBase *, Int_t=0);
//    static TLegend BuildLegend(THnBase *, TObjArray *, Int_t=0);
    static Double_t *GetDataArray(TObjArray *, Long64_t &, Int_t=0);
    static void SetLimits(TObjArray *&, Int_t=0);
    static Double_t GetStatVal(Double_t *, Long64_t, const TString, Int_t=0);
    template <typename T>
      static void SetFitter(T *&, Int_t=0);
    template <typename T>
      static void SetDrawingOptions(T *&, Int_t=0);
    static std::vector<TString> ParseString(const char *iString,  const char *sep="()", Int_t verbose=0);
    static std::vector<TString> ParseOptionString(const char *optionsStr, Int_t defSize=0, const char sep=',', const char ignoreBrackets[2]="()", Int_t verbose=0);
    static void ParsePandasString(const TString optionsStr, std::map<TString, TString> &optMap, const char sep=',', const char ignoreBrackets[2]="[]", Int_t verbose=0);
    static void RegisterDefaultOptions();
    static void RangesToMap (TString range , Int_t axisNum, axisRangesMap &result);
    static void RangesMapToString(Int_t n, axisRangesMap result, std::vector<TString> arr, std::vector<TString> &res);
    ClassDef(AliPainter,1);
};
#endif
