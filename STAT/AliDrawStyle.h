#ifndef ALIDRAWSTYLE_H
#define ALIDRAWSTYLE_H
 /* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \ingroup STAT
/// \class AliDrawStyle
/// \brief AliDrawStyle
///  ### AliDrawStyle - Class to access drawing styles
///  * Several drawing styles can be registered and used in the same moment
///    * Styles are identified using strings as identifiers
///      * TStyle
///      * MarkerStyle[]  AliDrawStyle::GetMarkerStyle(const char *style, Int_t index);
///      * MarkerColors[] AliDrawStyle::GetMarkerColor(const char *style, Int_t index);
///      * FillColors[]   AliDrawStyle::GetFillColor(const char *style, Int_t index);
///  * Default styles are created  AliDrawStyle::SetDefaults()
///    * default style is based on the fig template -  https://twiki.cern.ch/twiki/pub/ALICE/ALICERecommendationsResultPresentationText/figTemplate.C
///    * users should be able to register their own styles (e.g in macros)
///  * Usage (work in progress)
///    * performance reports -  with styles as a parameter
///    * QA reports
///    * AliTreePlayer, and TStatToolkit
/// \author marian  Ivanov marian.ivanov@cern.ch
///
///  ## Example usage
/*!
\code
  AliDrawStyle::SetDefaults()
  // Style example
  //
  AliDrawStyle::PrintStyles(0,TPRegexp("."));
  AliDrawStyle::ApplyStyle("figTemplate");
  gPad->UseCurrentStyle();  // force current style for current data
  //
  // Standard ALICE latex symbols
  AliDrawStyle::PrintLatexSymbols(0,TPRegexp("."))
  AliDrawStyle::GetLatexAlice("qpt")
  AliDrawStyle::AddLatexSymbol("dphi", "#Delta#it#phi (unit)")
  AliDrawStyle::GetLatexAlice("dphi")
  //
  // Standard ALICE marker/colors arrays
  AliDrawStyle::GetMarkerStyle("figTemplate",0)
  AliDrawStyle::GetMarkerColor("figTemplate",0)
\endcode
*/



#include "TObject.h"
#include <map>
#include <vector>
#include <string>
#include "TString.h"
#include <iostream>
#include <fstream>
class TPRegexp;
class TStyle;
class TPad;
class TCanvas;
class TF1;
class TH1;
class TGraph;
class TLegend;
class TMultiGraph;
class TFrame;
class TPaveText;
class AliDrawStyle : public TObject{
public:
  static void ApplyStyle(const char * styleName);
  static const TStyle *GetStyle(const char * styleName) {return fStyleAlice[styleName];}
  static void SetDefaults();
  static void SetDefaultStyles(const char * styleName, const char* arrayName);
  static TString GetLatexAlice(const char * symbol);
  static void AddLatexSymbol(const char * symbolName, const char * symbolTitle);
  static const std::vector<int> &    GetMarkerStyles(const char *style){return AliDrawStyle::fMarkerStyles[style];};
  static const std::vector<float> &  GetMarkerSize(const char *style){return AliDrawStyle::fMarkerSize[style];};
  static const std::vector<int> &    GetMarkerColors(const char *style){return AliDrawStyle::fMarkerColors[style];};
  static const std::vector<float> &  GetLineWidth(const char *style){return AliDrawStyle::fLineWidth[style];};
  static const std::vector<int> &    GetFillColors(const char *style){return AliDrawStyle::fFillColors[style];};
  static Int_t   GetIntegerAt(const char * format, Int_t index, const char * separator=";");
  static Float_t GetFloatAt(const char * format, Int_t index, const char * separator=";");
  static Int_t   GetMarkerStyle(const char *style, Int_t index);
  static Float_t GetMarkerSize(const char *style, Int_t index);
  static Int_t   GetMarkerColor(const char *style, Int_t index);
  static Int_t   GetFillColor(const char *style, Int_t index);
  static Int_t   GetLineStyle(const char *style, Int_t index);
  static Int_t   GetLineColor(const char *style, Int_t index);
  static Float_t GetLineWidth(const char *style, Int_t index);
  static void PrintLatexSymbols(Option_t *option,TPRegexp& regExp);
  static void PrintStyles(Option_t *option, TPRegexp& regExp);
protected:
  static TString fDefaultTStyleID;                            ///< ID of the default TStyle
  static TString fDefaultArrayStyleID;                        ///< ID of the default array styles
  static std::map<TString, TString> fLatexAlice;              ///< map of predefined latex symbols - formatted according ALICE rules
  static std::map<TString, TStyle*>  fStyleAlice;             ///< map of Alice predefined styles (+user defined)
  static std::map<TString, TObjArray*>  fCssStyleAlice;       ///< map of Alice predefined styles corresponding to css notation
  static std::map<TString, std::vector<int> > fMarkerStyles;  ///< map of predefined marker styles arrays
  static std::map<TString, std::vector<int> > fMarkerColors;  ///< map of predefined colors  arrays
  static std::map<TString, std::vector<float> > fMarkerSize;  ///< map of predefined marker sizes ()
  static std::map<TString, std::vector<int> > fFillColors;    ///< map of predefined fill colors arrays
  static std::map<TString, std::vector<float> > fLineWidth;   ///< map of predefined line width
  static std::map<TString, std::vector<float> > fLineStyle;   ///< map of predefined line style
  static std::map<TString, std::vector<float> > fLineColor;   ///< map of predefined line color
  //
  static void  RegisterDefaultLatexSymbols();                 ///< initialize default LatexSymbols
  static void  RegisterDefaultStyle();                        ///< initialize default TStyles
  static void  RegisterDefaultMarkers();                      ///< initialize default Markers/Colors

public:
  //NEW code for AliDrawStyle:
  static const TObjArray *GetCssStyle(const char *styleName) {return fCssStyleAlice[styleName];}
  template <typename T>
    static void ObjectApplyStyle(const char* styleName, T *cObj, Int_t objNum=0, Int_t verbose=0);
  static void TGraphApplyStyle(const char* styleName, TGraph *cGraph, Int_t objNum=0, Int_t verbose=0);
  static void TH1ApplyStyle(const char* styleName, TH1 *cHis, Int_t objNum=0, Int_t verbose=0);
  static void TF1ApplyStyle(const char* styleName, TF1 *cFunc, Int_t objNum=0, Int_t verbose=0);
  static void TLegendApplyStyle(const char* styleName, TLegend *cLegend, Int_t objNum=0, Int_t verbose=0);
  static void TPadApplyStyle(const char* styleName, TPad *cPad, Int_t verbose=0);
  static void TCanvasApplyCssStyle(const char* styleName, TCanvas *cCanvas, Int_t verbose=0);
  //static void TFrameApplyCssStyle(const char* styleName, TFrame *cFrame);
  //static void TPaveTextApplyCssStyle(const char* styleName, TPaveText *cPaveText);
  static void ApplyCssStyle(TPad *pad, const char* styleName, Int_t verbose=0);
  static void RegisterCssStyle(const char *styleName, TObjArray*array ) { fCssStyleAlice[styleName]=array;}
    //parsing methods
  static Bool_t     IsSelected(TString selectors, TString elementID, TString classID, TString objectID);
  static TString    GetValue(const char * styleName, TString propertyName, TString elementID, TString classID, TString objectID, TString localStyle=TString(""), Int_t verbose=0);
  static TString    ParseDeclaration(const char *inputDec, const char *propertyName);
  static void       GetIds(TObject *cObj, TString &elementID, TString &classID, TString &objectID, TString &localStyle, Int_t verbose=0);
  template <typename T>
    static T GetNamedTypeAt(const char *inputStr, Bool_t &status, Int_t index=0, const char *propertyName="", Int_t verbose=0, const char sep=',', const char *ignoreBrackets="()");
  static Float_t    ConvertUnit(const char *inputValues, const char * option="", Int_t verbose=0);
  static Int_t      ConvertColor(const char *inputString, Int_t verbose=0);
  static Float_t    PixelsToFloat_t(const char *value, const char *option="", Int_t verbose=0);
  static Float_t    PercentToFloat_t(const char *value, Int_t verbose);
  static Int_t      RgbToColor_t(const char *inputString, Int_t verbose=0);
  static Int_t      HexToColor_t(const char *inputString, Int_t verbose=0);
  static TObjArray *ReadCSSFile(const char *  inputName, TObjArray * array=nullptr, Int_t verbose=0);
  static void       WriteCSSFile(TObjArray * cssArray, const char *  outputName, std::fstream *cssOut=nullptr);
  template <typename T>
    static T PrepareValue(const char* styleName, TString propertyName, TString elementID, TString classID, TString objectID, TString localStyle, Bool_t &status, Int_t objNum=0, Int_t verbose=0);
    // CSS like attribute fields parsing
  ClassDef(AliDrawStyle,1);
  private:
    static Int_t padNumber;
    static void SetPadNumber(Int_t num) {padNumber = num;};
    static Int_t GetPadNumber() {return padNumber;};
};

#endif