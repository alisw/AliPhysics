#ifndef ALIPAINTER_H
#define ALIPAINTER_H

/// \ingroup STAT
/// \class AliPainter
/*!
 *
* \brief Class for generating QA reports
*  See the documentation in describing of functions.
* \author Marian  Ivanov marian.ivanov@cern.ch, Boris Rumyantsev boris.rumyantsev@cern.ch
*/

class TPad;
class TMultiGraph;
class TFormula;
#include "TObject.h"
#include <vector>
#include <map>
#include "TObjArray.h"
#include "TString.h"
#include "TPad.h"

class AliPainter : public TObject{
  public:
    static TPad *DivideTPad(const char *division, const char *classID="", const char *style="", TPad *pad=nullptr, Int_t verbose=0);
    static void SetMultiGraphTimeAxis(TMultiGraph *graph, TString option);
    static void DrawHistogram(char *expresion, const TObjArray *histogramArray, TPad *pad=nullptr, TObjArray *metaData=nullptr, TObjArray *keepArray=nullptr, Int_t verbose=4);
    static TPad *SetPadMargin(TPad *cPad, const char *position, const char *wMargin, const char *units, Double_t mValue, Int_t iCol, Int_t nCols);
    //static Double_t GetLimitValue(TString expression, TPad *cPad);
    //static void ApplyLimitValue();
    ClassDef(AliPainter,1);
  //TODO: it should be protected for inheritor test class @Boris
  //protected:
    typedef std::map<Int_t, std::vector<TString> > axisRangesMap;
    static std::map<TString, TString> optionValues;
    static std::map<TString, Double_t> statValues;
    static void ArgsParser(TString exprsn, TString &hisName, TString &projections, std::vector<TString> &fitOptions, std::vector<TString> &rangesStrings, Int_t verbose = 4);
    static std::vector<TString> OptionStringParser(const char *option, const char d[2], Int_t defSize);
    static void PandasOptionParser(const TString optionsStr);
    static void RegisterDefaultOptions();
    static std::vector<TString> RangesParser(TString inputRangesStr);
    static void RangesToMap (TString range , Int_t axisNum, axisRangesMap &result);
    static void RangesMapToString(Int_t n, axisRangesMap result, std::vector<TString> arr, std::vector<TString> &res);
};

#endif
