#ifndef ALIDRAWSTYLE_H
#define ALIDRAWSTYLE_H
 /* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \ingroup STAT
/// \class AliDrawStyle
/// \brief Class to access to drawing styles
///
///   _Class still experimental. It means the interface could be change._
///   Our goal is provide to users comfortable interface for changing style of drawing complicated reports,
///   which could be create via AliPainter.
///
///   We think that best tools for such tasks it's CSS. So this is the attempt to create analogue of CSS for aliroot.
///   All you need to know about CSS for usage AiDrawStyle is:
///   Selectors are used to target which HTML to style. Properties and values are used to set the style rules.
///   There are three kinds of selectors:
///       Name      | Description
///   ------------- | -------------
///       Tag       | An HTML tag such as h1 or p. (In our case this is a root class TH1, ...)
///      Class      | A class attribute of one or more elements, such as TH1F *fHist = new TH1F("exampleHisto.class(classA)", "anyTitle",100,-5,5). Referenced in CSS with a “.” before it.
///       ID        | An id attribute of a unique element, should only be used once, such as name of object("exampleHisto" from example above). Referenced in CSS with a “#” before it.
///
///   The sequence of work with AliDrawStyle:
///     1. Define the style (create file with description of style);
///     2. Register this style;
///     3. Apply it recursively for your objects on TPad or TCanvas;
///
///   Method of applying the style - AliDrawStyle::ApplyCssStyle() has three input parameters:
///     1. Name of style (it has a type - const char *, so for using TString use .Data());
///     2. TPad object in order to apply your style recursively according with you style;
///     3. Verbosity mode (4 for detailed output);
///
///   Also we support units and color formats. So, you can use rgb or hex notation for colors in your .css file and pixels or percents for sizes.
/// \author [Marian Ivanov](marian.ivanov@cern.ch), [Boris Rumyantsev](boris.rumyantsev@cern.ch)
///  ## Example usage:
///     First of all you should create and fill the file with extension -".css" with description
///     of your own styles according to your objects declaration. The general syntax of CSS you can see
///     <a href="https://www.w3schools.com/css/css_syntax.asp">here.</a>
///         + Selectors are patterns used to select the element(s) you want to style.
///         + Declaration consist from  properties and values.
///         We tried to use the same properties like in CSS
///         where possible and in the opposite case we use properties from root.
///
///     ### Let's see to examples:
///
///         First of all create the canvas with some histograms:
///
///                \code
///                 TCanvas *canvasT = new TCanvas("canvasT", "canvasT", 900, 500);
///                 canvasT->Divide(1,1);
///                 canvasT->cd(1);
///                 TH1F *his1 = new TH1F("his1.class(Error).style(line-color:#f30000;)", "his1", 100, -5,5);
///                 his1->FillRandom("gaus", 15000);
///                 his1->SetStats(kFALSE);
///                 TH1D *his2 = new TH1D("his2.class(Error)", "his2", 100, -5,5);
///                 his2->FillRandom("gaus", 10000);
///                 TH1I *his3 = new TH1I("his3.class(Error)", "his3", 100, -5,5);
///                 his3->FillRandom("gaus", 5000);
///                 his1->Draw();
///                 his2->Draw("same");
///                 his3->Draw("same");
///                 gStyle->SetOptTitle(0);
///                 TPaveText *pt = new TPaveText(-0.438183,694.575009,0.438183,740.053135);
///                 pt->AddText("Example of styling with using AliDrawStyle");
///                 pt->SetTextSize(0.04)
///                 pt->SetShadowColor(0)
///                 pt->SetBorderSize(0)
///                 pt->SetFillColor(0)
///                 pt->Draw()
///                 gPad->BuildLegend();
///                \endcode
///             ![Example for styling](AliDrawStyle_h_example1.png)
///
///             ##And now we will create the our own style:
///               1. Canvas properties:
///                 In case if you want to change some properties of canvas, you should specify one of selectors of this
///                 canvas. In our case we will use element it's "TCanvas", so first of all we add to our css file next
///                 strings:
///                 \code
///                   TCanvas   {
///                     height: 650;
///                     weight: 400;
///                   }
///                 \endcode
///               2. Pad properties:
///                 Let's also use element ("TPad") for pad selection, _but you should remember: in case you have few pads
///                 properties below selector "TPad" will be apply to each TPad from your canvas:_
///                 \code
///                   TPad {
///                     gridX: 1;
///                     gridY: 1;
///                     fill-color: rgb(254,254,254);
///                     tickX: 1;
///                     tickY: 1;
///                     margin-bottom: 0.15;
///                   }
///                 \endcode
///               3. Objects properties:
///                   We have 3 objects with such selectors:
///                   - elements - "TH1F","TH1D", "TH1I";
///                   - class   - "Error";
///                   - objects - "his1", "his2", "his3";
///
///                  It means, that we can be able to manage it via three different ways:
///
///               * <B>elements</B> - For elements you have two ways: you can use unique name of each elements
///                                   ("TH1F", "TH1D", "TH1I") and in this way you can apply your style
///                                   separately for each object or you can use regexp(<B>only star available</B>).
///                 Let's use this way with regexp:
///                 \code
///                     TH1* {
///                       line-width: 4;
///                       line-style: 1;
///                       line-color: #000000,#0f0fbb,#00d200;
///                       fill-color: #730000,#050545,#003900;
///                     }
///                 \endcode
///                  These properties will be common for all histograms.
///                  Pay your attention to this string: line-color: #f30000,#0f0fbb,#00d200;
///                  The style will be apply to each object recursively, and number of value
///                  will be equal to number of object.
///
///               * <B>classes</B> - In our case class also provide to access for styling all objects in one time:
///                 \code
///                     .Error  {
///                       title-size: 0.2;
///                       label-size: 0.05;
///                     }
///                 \endcode
///
///               * <B>objects</B> - Finally using names of objects ("his1", "his2", "his3")
///                                  you can manage to style of them separately.
///                 \code
///                     #his1 {
///                       fill-color: #730000;
///                     }
///
///                     #his2 {
///                       fill-color: #050545;
///                     }
///
///                     #his3 {
///                       fill-color: #003900;
///                     }
///                 \endcode
///
///               4. Other properties:
///                  Also we can manage other object such as TLegend, TF1, TGraph and classes which derived by them.
///                  \code
///                     TLegend  {
///                       line-width: 3;
///                       text-size: 0.04;
///                       text-color: rgb(50,50,50);
///                       x1: 2;
///                       x2: 4;
///                       y1: 500;
///                       y2: 600;
///                     }
///                  \endcode
///
///         One more important thing - we provide to users opportunity to specify local property
///         in the name of objects, such property has more priority than property from file.
///         If you defined the name of object like in example above,
///         you can specify priority attributes via .style("property:value;")
///         \code
///           TH1F *his1 = new TH1F("his1.class(Error).style(line-color:#f30000;)", "his1", 100, -5,5);
///         \endcode
///         As you can see here, in AliDrawStyleTutor.css for his1 we have black color for line, but color is red,
///         because local style in name has more priority.
///
///         Finally, content of our .css file will be so:
///         \code
///           TCanvas {
///             height: 1200;
///             width: 800;
///           }
///
///           TPad {
///             gridX: 1;
///             gridY: 1;
///             fill-color: rgb(254,254,254);
///             tickX: 1;
///             tickY: 1;
///             margin-bottom: 0.15;
///           }
///
///           TH1* {
///             line-width: 4;
///             line-style: 1;
///             line-color: #000000,#0f0fbb,#00d200;
///           }
///
///           .Error  {
///             title-size: 0.2;
///             label-size: 0.05;
///           }
///
///           TLegend  {
///             line-width: 3;
///             text-size: 0.04;
///             text-color: rgb(50,50,50);
///             x1: 2;
///             x2: 4;
///             y1: 500;
///             y2: 600;
///           }
///
///           #his1 {
///             fill-color: #730000;
///           }
///
///           #his2 {
///             fill-color: #050545;
///           }
///
///           #his3 {
///             fill-color: #003900;
///           }
///         \endcode
///         After applying
///         \code
///           AliDrawStyle::RegisterCssStyle("AliDrawStyleTutor", AliDrawStyle::ReadCSSFile("$AliRoot_SRC/STAT/test/AliDrawStyleTutor.css"));
///           AliDrawStyle::ApplyCssStyle(canvasT, "AliDrawStyleTutor");
///         \endcode
///         you will get such image:
///           ![Example after styling](AliDrawStyle_h_example2.png)
///         Example is available here: /../AliRoot/STAT/test/AliDrawStyleTutor.css

// TODO: Should we use some parser like bison, yacc or boost::spirit?
// TODO: create macro with handle style for test ApplyCssStyle
// TODO: change pad global indexes to local
// TODO: using TString - to be replaced by faster variant with rough pointers
// TODO: refactor godness function
// TODO: perhaps I should combine {TString elementID, TString classID, TString objectID, TString localStyle} to the array?
// TODO: add test for ReadCssString, ApplyCssStyle, WriteCssFile, private methods
// TODO: combine ElementSearch, ClassSearch, ObjectSearch into one function
// TODO: add parsing scheme

#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "TObject.h"
#include "TString.h"
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

class AliDrawStyle: public TObject {
  //OLD CODE
public:
  static void ApplyStyle(const char *styleName);
  static const TStyle *GetStyle(const char *styleName) { return fStyleAlice[styleName]; }
  static void SetDefaults();
  static void SetDefaultStyles(const char *styleName, const char *arrayName);
  static TString GetLatexAlice(const char *symbol);
  static void AddLatexSymbol(const char *symbolName, const char *symbolTitle);
  static const std::vector<int> &GetMarkerStyles(const char *style) { return AliDrawStyle::fMarkerStyles[style]; };
  static const std::vector<float> &GetMarkerSize(const char *style) { return AliDrawStyle::fMarkerSize[style]; };
  static const std::vector<int> &GetMarkerColors(const char *style) { return AliDrawStyle::fMarkerColors[style]; };
  static const std::vector<float> &GetLineWidth(const char *style) { return AliDrawStyle::fLineWidth[style]; };
  static const std::vector<int> &GetFillColors(const char *style) { return AliDrawStyle::fFillColors[style]; };
  static Int_t GetIntegerAt(const char *format, Int_t index, const char *separator = ";");
  static Float_t GetFloatAt(const char *format, Int_t index, const char *separator = ";");
  static Int_t GetMarkerStyle(const char *style, Int_t index);
  static Float_t GetMarkerSize(const char *style, Int_t index);
  static Int_t GetMarkerColor(const char *style, Int_t index);
  static Int_t GetFillColor(const char *style, Int_t index);
  static Int_t GetLineStyle(const char *style, Int_t index);
  static Int_t GetLineColor(const char *style, Int_t index);
  static Float_t GetLineWidth(const char *style, Int_t index);
  static void PrintLatexSymbols(Option_t *option, TPRegexp &regExp);
  static void PrintStyles(Option_t *option, TPRegexp &regExp);
protected:
  static TString fDefaultTStyleID;                            ///< ID of the default TStyle
  static TString fDefaultArrayStyleID;                        ///< ID of the default array styles
  static std::map <TString, TString> fLatexAlice;              ///< map of predefined latex symbols - formatted according ALICE rules
  static std::map<TString, TStyle *> fStyleAlice;             ///< map of Alice predefined styles (+user defined)
  static std::map<TString, TObjArray *> fCssStyleAlice;       ///< map of Alice predefined styles corresponding to css notation
  static std::map <TString, std::vector<int> > fMarkerStyles;  ///< map of predefined marker styles arrays
  static std::map <TString, std::vector<int> > fMarkerColors;  ///< map of predefined colors  arrays
  static std::map <TString, std::vector<float> > fMarkerSize;  ///< map of predefined marker sizes ()
  static std::map <TString, std::vector<int> > fFillColors;    ///< map of predefined fill colors arrays
  static std::map <TString, std::vector<float> > fLineWidth;   ///< map of predefined line width
  static std::map <TString, std::vector<float> > fLineStyle;   ///< map of predefined line style
  static std::map <TString, std::vector<float> > fLineColor;   ///< map of predefined line color
  static void RegisterDefaultLatexSymbols();                 ///< initialize default LatexSymbols
  static void RegisterDefaultStyle();                        ///< initialize default TStyles
  static void RegisterDefaultMarkers();                      ///< initialize default Markers/Colors


  // NEW CODE
public:
  static void       RegisterCssStyle(const char *styleName, TObjArray*array ) { fCssStyleAlice[styleName]=array;}
  static const TObjArray *GetCssStyle(const char *styleName) {return fCssStyleAlice[styleName];}
  static void       WriteCSSFile(TObjArray *cssArray, const char *  outputName, std::fstream *cssOut=nullptr);
  static TObjArray *ReadCSSFile(const char *inputName, TObjArray * array=nullptr, Int_t verbose=0);
  static TObjArray *ReadCssString(TString cssString, TObjArray *array=nullptr, Int_t verbose=0);
  static void       ApplyCssStyle(TPad *pad, const char *styleName, Int_t verbose=0);
//protected:
  template <typename T>
    static void     TObjectApplyStyle(const char* styleName, T *cObj, Int_t objNum=0, Int_t verbose=0);
  static void       TGraphApplyStyle(const char* styleName, TGraph *cGraph, Int_t objNum=0, Int_t verbose=0);
  static void       TH1ApplyStyle(const char* styleName, TH1 *cHis, Int_t objNum=0, Int_t verbose=0);
  static void       TF1ApplyStyle(const char* styleName, TF1 *cFunc, Int_t objNum=0, Int_t verbose=0);
  static void       TLegendApplyStyle(const char* styleName, TLegend *cLegend, Int_t objNum=0, Int_t verbose=0);
  static void       TPadApplyStyle(const char* styleName, TPad *cPad, Int_t verbose=0);
  static void       TCanvasApplyCssStyle(const char* styleName, TCanvas *cCanvas, Int_t verbose=0);

//public:
  static void       GetIds(TObject *cObj, TString &elementID, TString &classID, TString &objectID, TString &localStyle, Int_t verbose=0);
  static Float_t    PrepareValue(const char* styleName, TString propertyName, TString elementID, TString classID, TString objectID,                                        TString localStyle, Bool_t &status, Int_t objNum=0, Int_t verbose=0);
  static TString    ParseDeclaration(const char *inputDec, const char *propertyName);
  static TString    GetValue(const char * styleName, TString propertyName, TString elementID, TString classID, TString objectID, TString localStyle=TString(""), Int_t verbose=0);
  static Bool_t     IsSelected(TString selectors, const TString elementID, const TString classID, const TString objectID, Int_t verbose=0);
  static Bool_t     ElementSearch(const TString selectors, const TString elementID, Int_t verbose=0);
  static Bool_t     ClassSearch(const TString selectors, const TString classID, Int_t verbose=0);
  static Bool_t     ObjectSearch(const TString selectors, const TString objectID, Int_t verbose=0);
  static Float_t    GetNamedTypeAt(const char *inputStr, Bool_t &status, Int_t index=0, const char *propertyName="", Int_t verbose=0, const char sep=',', const char *ignoreBrackets="()");
  static Float_t    ConvertUnit(const char *inputValues, const char * option="", Int_t verbose=0);
  static Int_t      ConvertColor(const char *inputString, Int_t verbose=0);
  static Float_t    PixelsToFloat_t(const char *value, const char *option="", Int_t verbose=0);
  static Float_t    PercentToFloat_t(const char *value, Int_t verbose);
  static Int_t      RgbToColor_t(const char *inputString, Int_t verbose=0);
  static Int_t      HexToColor_t(const char *inputString, Int_t verbose=0);

  static Int_t      padNumber;
  static void       SetPadNumber(Int_t num) {padNumber = num;};
  static Int_t      GetPadNumber() {return padNumber;};
ClassDef(AliDrawStyle,1);
};

#endif
