#ifndef ALIPAINTER_H
#define ALIPAINTER_H

/// \ingroup STAT
/// \class AliPainter
/*!
\brief AliPainter
\code
/// TODO add documentation
\endcode

  \author Marian  Ivanov marian.ivanov@cern.ch, Boris Rumyantsev boris.rumyantsev@cern.ch
*/

class TPad;
class TMultiGraph;
#include "TObject.h"
#include <vector>
#include <map>
#include "TObjArray.h"
#include "TString.h"
#include "TPad.h"

class AliPainter : public TObject{

public:
    static void       DivideTPad(TPad *pad, const char *division, const char *classID);
    static void       SetMultiGraphTimeAxis(TMultiGraph *graph, TString option);
    static void       DrawHistogram(char *expresion, const TObjArray *histogramArray, TPad *pad=NULL, TObjArray *metaData=NULL, TObjArray *keepArray=NULL, Int_t verbose=4);
//static TObject   *GetHistogram(TObject *hisN, TString range, std::vector<TString> fitOptions);
    static TPad      *SetPadMargin(TPad *cPad, const char *position, const char *wMargin, const char *units, Double_t mValue, Int_t iCol, Int_t nCols);
public:
    typedef std::map<int, std::vector<int> > axisRangesMap;
    static std::map<TString, TString> optionValues;
ClassDef(AliPainter,1);
    static void ArgsParser(TString exprsn, TString &hisName, TString &projections, std::vector<TString> &fitOptions, std::vector<TString> &rangesStrings, Int_t verbose = 4);
    static std::vector<TString> OptParser(const TString fitStr, const char d[2], Int_t defSize);
    static void OptionsParser(const TString optionsStr);
    static void  RegisterDefaultOptions();
    static std::vector<TString> RangesParser(TString inputRangesStr);
    static void RangesToMap (TString range , Int_t axisNum, axisRangesMap &result);
    static void RangesMapToString(Int_t n, axisRangesMap result, std::vector<TString> arr, std::vector<TString> &res);
};

#endif
