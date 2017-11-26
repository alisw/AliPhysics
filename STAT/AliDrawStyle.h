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
#include <iostream>     // std::cout
#include <fstream>
class TPRegexp;
class TStyle;
class TPad;
class TCanvas;
class TF1;
class TH1;
class TGraph;
class TMultiGraph;
class AliDrawStyle : public TObject{
public:
  static void ApplyStyle(const char * styleName);
  static const TStyle *GetStyle(const char * styleName) {return fStyleAlice[styleName];}
  /// \brief Get Css style by styleName.
  ///  Getter for css style.
  /// \param styleName  - name of style.
  /// \return           - TObjArray with pair selector and declaration from css file.
  static const TObjArray *GetCssStyle(const char *styleName){return fCssStyleAlice[styleName];}
  static void TGraphApplyStyle(const char* styleName, TGraph *cGraph);
  static void TH1ApplyStyle(const char* styleName, TH1 *cHis);
  static void TF1ApplyStyle(const char* styleName, TF1 *cFunc);
  static void TPadApplyStyle(const char* styleName, TPad *cPad);
  static void TCanvasApplyCssStyle(const char* styleName, TCanvas *cCanvas);
  static void ApplyCssStyle(TPad *pad, const char* styleName);
  /// \brief Set Css style by styleName.
  ///  Setter for css style.
  /// \param styleName  - name of style.
  /// \param array      - TObjArray with pair selector and declaration from css file.
  static void RegisterCssStyle(const char *styleName, TObjArray*array ){ fCssStyleAlice[styleName]=array;}
  static void SetDefaults();
  static void SetDefaultStyles(const char * styleName, const char* arrayName);
  static TString GetLatexAlice(const char * symbol);
  static void AddLatexSymbol(const char * symbolName, const char * symbolTitle);
  static const std::vector<int> &    GetMarkerStyles(const char *style){return AliDrawStyle::fMarkerStyles[style];};
  static const std::vector<float> &  GetMarkerSize(const char *style){return AliDrawStyle::fMarkerSize[style];};
  static const std::vector<int> &    GetMarkerColors(const char *style){return AliDrawStyle::fMarkerColors[style];};
  static const std::vector<float> &  GetLineWidth(const char *style){return AliDrawStyle::fLineWidth[style];};
  static const std::vector<int> &    GetFillColors(const char *style){return AliDrawStyle::fFillColors[style];};
  // CSS like attribute fields parsing
  static Bool_t  IsSelected(TString selectors, TString elementID, TString classID, TString objectID);
  static TString GetProperty(const char * styleName, TString propertyName, TString elementID, TString classID, TString objectID);
  static TString  GetPropertyValue(TString input, TString propertyName);
  //static Int_t    GetObjectIndex(TString &objName);
  static void     GetIds(TObject *cObj, TString &elementID, TString &classSet, TString &objectID, Int_t &objNum);
  static Int_t    GetNamedIntegerAt(TString input, TString propertyName, Int_t index, Bool_t &status);
  static Float_t  GetNamedFloatAt(TString input, TString propertyName, Int_t index, Bool_t &status);
  static Double_t UnitsConverter(TString value, Double_t k);
  static TObjArray * ReadCSSFile(const char *  inputName, TObjArray * array=NULL, Int_t verbose=0);
  static void    WriteCSSFile(TObjArray * cssArray, const char *  outputName, std::fstream *cssOut=NULL);
  //
  //
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

  ClassDef(AliDrawStyle,1);
};
