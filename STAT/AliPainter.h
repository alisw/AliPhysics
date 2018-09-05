#ifndef ALIPAINTER_H
#define ALIPAINTER_H

/// \ingroup STAT
/// \class AliPainter
/*!
* \brief Class for generating QA reports
*  See the documentation in describing of functions.
* \author  <a href="marian.ivanov@cern.ch">Marian  Ivanov</a> , <a href="boris.rumyantsev@cern.ch">Boris Rumyantsev</a>
*/


// TODO: add divFlag for inheritance Mother Class (ClassName) should be add to children (className)=>nameOfObject.class(ClassMother,ClassDaughter) @Boris
// TODO: perhaps, we should use some aliases? @Marian
// TODO: may be we should also use TVectorD or TMatrix instead enumeration? @Boris
// TODO: extend to TH2D, TH3D @Boris
// TODO: think how to combine  AliPainter::ParseString and AliPainter::ParseOptionString @Boris
// NOTE: we can provide list of brackets and remove ignoreBrackets from arguments. mb just add flag about do you want to ignore brackets or not.
// TODO - now we have 2 steps for parsing of ranges option to array of string with simple range values. 1. Initial string into map. 2. Map into array fo strings. First of all we can use vector of vector instead map, and the second we should think how we can avoid this 2 intermediate steps with map. @Boris
// TODO: add few classes with []
// TODO: change global maps to arguments of functions

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
#include "AliParser.h"

class AliPainter {
  public:
    static TPad *DivideTPad(TPad *pad, const char *division, const char *classID="", const char *style="", Int_t verbose=0);
    static void SetMultiGraphTimeAxis(TMultiGraph *graph, TString option);
    static void DrawHistogram(THnBase *hisN, const char *expression, TPad *pad=nullptr, TObjArray *keepArray=nullptr, TObjArray *metaData=nullptr, Int_t verbose=0);
    static void DrawHistogram(TObjArray *histogramArray, const char *expression, TPad *pad=nullptr, TObjArray *keepArray=nullptr, TObjArray *metaData=nullptr, Int_t verbose=0);
    static TPad *GetNextPad(TPad *cPad, TPad *tempPad=nullptr, Int_t verbose=0);
  private:
    static void SaveToKeepArray(TObject *, TObjArray *&, Int_t=0);
    static void SaveToKeepArray(TObjArray *, TObjArray *&, Int_t=0);
    static TPad *SetPadMargin(TPad *cPad, const char *position, const char *wMargin, const char *units, Double_t mValue, Int_t iCol, Int_t nCols);
    static TObjArray *PrepareHistogram(THnBase *hisN, const char *expression, TObjArray *&keepArray, TObjArray *metaData=nullptr, Int_t verbose=0);
    static TObjArray* SliceHistogram(THnBase *, TString, Int_t=0);
    static TObject *SetProjections(THnBase *, TString, Int_t=0);
//    static TLegend BuildLegend(THnBase *, TObjArray *, Int_t=0);
    static Double_t *GetDataArray(TObjArray *, Long64_t &, Int_t=0);
    static void SetLimits(TObjArray *&, TString, Int_t=0);
    static Double_t GetStatVal(Double_t *, Long64_t, const TString, Int_t=0);
    template <typename T>
      static void FitHistogram(T *&, std::map<TString, TString>, Int_t=0);
    template <typename T>
      static void SetDrawingOptions(T *&, std::map<TString, TString>, Int_t=0);
    ClassDef(AliPainter,1);
  private:
};
#endif
