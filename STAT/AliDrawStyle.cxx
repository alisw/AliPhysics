/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


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


#include "AliDrawStyle.h"
#include "TStyle.h"
#include "TError.h"
#include "TPRegexp.h"
#include "TColor.h"
#include "TMath.h"
#include <iostream>
//
TString AliDrawStyle::fDefaultTStyleID;                            ///< ID of the default TStyle
TString AliDrawStyle::fDefaultArrayStyleID;                        ///< ID of the default array styles
std::map<TString, TString>  AliDrawStyle::fLatexAlice;
std::map<TString, TStyle*>  AliDrawStyle::fStyleAlice;
std::map<TString, std::vector<int> > AliDrawStyle::fMarkerStyles;  // PLEASE LEAVE THE UNAESTHETIC SPACE
std::map<TString, std::vector<int> > AliDrawStyle::fMarkerColors;  // IN ORDER TO MAKE IT WORK WITH THE
std::map<TString, std::vector<float> > AliDrawStyle::fMarkerSize;  // NATIVE SLC6 COMPILER!!!
std::map<TString, std::vector<int> > AliDrawStyle::fFillColors;
std::map<TString, std::vector<float> > AliDrawStyle::fLineWidth;

void AliDrawStyle::SetDefaults(){
  AliDrawStyle::RegisterDefaultLatexSymbols();
  AliDrawStyle::RegisterDefaultStyle();
  AliDrawStyle::RegisterDefaultMarkers();
}

/// set AliDrawStyle::SetDefaultStyles
/// \param tstyleName  - default style to be used class function in case of empty style selection
/// \param arrayName   - default style to be used class function in case of empty array style selection
void AliDrawStyle::SetDefaultStyles(const char * tstyleName, const char* arrayName){
  fDefaultTStyleID=tstyleName;
  fDefaultArrayStyleID=arrayName;
}

TStyle* RegisterDefaultStyleFigTemplate(Bool_t grayScale);

/// Latex symbol section

/// \param symbol  -id name of latex symbol
/// \return latex symbol to be used in ROOT latex
TString AliDrawStyle::GetLatexAlice(const char * symbol){
  return  fLatexAlice[symbol];
}


/// Get integer from string at index
/// \param format    -  array string
/// \param index     -  element index
/// \param separator -  array separator
/// TODO: using TString - to be replaced by faster variant with rough pointers
/// Example usage:
/*!
\code
   AliDrawStyle::GetIntegerAt("1:4:8",1,":"); // return 4
   AliDrawStyle::GetIntegerAt("1;4;8",2,";"); // return  8
\endcode
 */
Int_t AliDrawStyle::GetIntegerAt(const char * format, Int_t index, const char * separator ){
  if (format==NULL) return -1;
  if (index<0) return -1;
  index++;
  TString sformat(format);
  TString token(format);
  Int_t position=0;
  Int_t counter=0;
  while (counter<index){
    if (sformat.Tokenize(token,position,separator)) {
      counter++;
    }else{
      break;
    }
  }
  if (counter==index) return token.Atoi();
  return -1;
}

/// Get integer from string
/// \param format    -  array string
/// \param index     -  element index
/// \param separator -  array separator
/// TODO: using TString - to be replaced by faster variant with rough pointers
/// Example usage:
/*!
\code
   AliDrawStyle::GetFloatAt("1.1:4.1:8.1",1,":"); // return 4.1
   AliDrawStyle::GetFloatAt("1.1;4.1;8.1",2,";"); // return  8.1
\endcode
*/
Float_t AliDrawStyle::GetFloatAt(const char * format, Int_t index, const char * separator ){
  if (format==NULL) return -1;
  if (index<0) return -1;
  index++;
  TString sformat(format);
  TString token(format);
  Int_t position=0;
  Int_t counter=0;
  while (counter<index){
    if (sformat.Tokenize(token,position,separator)) {
      counter++;
    }else{
      break;
    }
  }
  if (counter==index) return token.Atof();
  return -1;
}


// GetMarkerStyle associated to the style.
/// \param  style - name of style used
/// \param index  - marker index
/// \return marker style for given stylename, index
Int_t AliDrawStyle::GetMarkerStyle(const char *style, Int_t index){
  if (AliDrawStyle::fMarkerStyles[style].size() <= index) {
    return GetIntegerAt(style,index);
  }
  return  AliDrawStyle::fMarkerStyles[style][index];
}

/// GetMarkerColor associated to the style.
/// \param  style - name of style used
/// \param index  - marker index
/// \return marker color for given stylename, index
Int_t AliDrawStyle::GetMarkerColor(const char *style, Int_t index){
  if (AliDrawStyle::fMarkerColors[style].size() <= index) {
    return GetIntegerAt(style,index);
  }
  return  AliDrawStyle::fMarkerColors[style][index];
}
/// GetMarkerSize associated to the style.
/// \param  style - name of style used
/// \param index  - marker index
/// \return marker color for given stylename, index
Float_t AliDrawStyle::GetMarkerSize(const char *style, Int_t index){
  if (AliDrawStyle::fMarkerSize[style].size() <= index) {
    return GetIntegerAt(style,index);
  }
  return  AliDrawStyle::fMarkerSize[style][index];
}
/// GetFillColor associated to the style.
/// \param  style - name of style used
/// \param index  - marker index
/// \return fill color for given stylename, index
Int_t AliDrawStyle::GetFillColor(const char *style, Int_t index){
  if (AliDrawStyle::fFillColors[style].size() <= index) {
    return GetIntegerAt(style,index);
  }
  return  AliDrawStyle::fFillColors[style][index];
}
/// GetLineWidth associated to the style.
/// \param  style - name of style used
/// \param index  - marker index
/// \return fill color for given stylename, index
Float_t AliDrawStyle::GetLineWidth(const char *style, Int_t index){
  if (AliDrawStyle::fLineWidth[style].size() <= index) {
    return GetFloatAt(style,index);
  }
  return  AliDrawStyle::fLineWidth[style][index];
}



void AliDrawStyle::PrintLatexSymbols(Option_t */*option*/, TPRegexp& regExp){
  //print latex symbols
  typedef std::map<TString,TString>::const_iterator it_type;
  for(it_type iterator = fLatexAlice.begin(); iterator != fLatexAlice.end(); ++iterator) {
    if (regExp.Match(iterator->first.Data())==0) continue;
    std::cout<<iterator->first << " " << iterator->second << "\n";
  }
}

void AliDrawStyle::PrintStyles(Option_t *option, TPRegexp& regExp){
  //print latex symbols
  typedef std::map<TString,TStyle*>::const_iterator it_type;
  for(it_type iterator = fStyleAlice.begin(); iterator != fStyleAlice.end(); ++iterator) {
    if (regExp.Match(iterator->first.Data())==0) continue;
    if (option==NULL) std::cout << iterator->first << " " << iterator->second << "\n";
    if (option!=NULL) {
      iterator->second->Print(option);
      if (TString(option).Contains("dump")) iterator->second->Dump();
    }
  }
}


void AliDrawStyle::ApplyStyle(const char* styleName){
  //
  // Apply registered style
  //
  TStyle * style= fStyleAlice[styleName];
  if (style==NULL){
    ::Error("AliDrawStyle::ApplyStyle","Invalid style %s",styleName);
  }else{
    ::Info("AliDrawStyle::ApplyStyle","%s",styleName);
  }
  if (style) style->cd();
}


void  AliDrawStyle::AddLatexSymbol(const char * symbolName, const char * symbolTitle){
  fLatexAlice[symbolName]=symbolTitle;
}
void  AliDrawStyle::RegisterDefaultLatexSymbols(){
  //
  // Set default AliRoot/Latex/root shortcuts
  //
  fLatexAlice["qpt"]="#it{q}/#it{p}_{T} (GeV/#it{c})^{-1}";
  fLatexAlice["qpt0"]="#it{q}/#it{p}_{T}";
  fLatexAlice["pt"]="#it{p}_{T}  (GeV/#it{c}) ";
  fLatexAlice["pt0"]="#it{p}_{T} ";
  //  fLatexAlice["sqptmev"]="#sigma_{#it{q}/#it{p}_{T}}/#it{p}_{T}^{2} (MeV/#it{c})^{-1}";
  fLatexAlice["sqptmev"]="#sigma_{#it{q}/#it{p}_{T}} (MeV/#it{c})^{-1}";
  fLatexAlice["pbpb502"]="Pb#font[122]{-}Pb #sqrt{#it{s}_{NN}} =5.02 TeV";
  fLatexAlice["pp13"]="pp #sqrt{#it{s}} = 13 TeV ";
  fLatexAlice["drphi"]="#Delta_{#it{r#phi}} (cm)";
  fLatexAlice["srphi"]="#sigma_{#it{r#phi}} (cm)";
}

void   AliDrawStyle::RegisterDefaultStyle(){
  //
  fStyleAlice["figTemplate"]=RegisterDefaultStyleFigTemplate(kFALSE);
  fStyleAlice["figTemplateGrey"]=RegisterDefaultStyleFigTemplate(kFALSE);
  //
  TStyle *style=RegisterDefaultStyleFigTemplate(kFALSE);
  style->SetName("figTemplate2"); 
  style->SetTitleXSize(TMath::Power(2,0.5)*style->GetTitleXSize());
  style->SetTitleYSize(TMath::Power(2,0.5)*style->GetTitleYSize());
  style->SetLabelSize(TMath::Power(2,0.5)*style->GetLabelSize("X"),"X");
  style->SetLabelSize(TMath::Power(2,0.5)*style->GetLabelSize("Y"),"Y");
  style->SetLabelSize(TMath::Power(2,0.5)*style->GetLabelSize("Z"),"Z");  
  fStyleAlice["figTemplate2"]=style;
  //
  style=RegisterDefaultStyleFigTemplate(kFALSE);
  style->SetName("figTemplate3"); 
  style->SetTitleXSize(TMath::Power(3,0.5)*style->GetTitleXSize());
  style->SetTitleYSize(TMath::Power(3,0.5)*style->GetTitleYSize());
  style->SetLabelSize(TMath::Power(3,0.5)*style->GetLabelSize("X"),"X");
  style->SetLabelSize(TMath::Power(3,0.5)*style->GetLabelSize("Y"),"Y");
  style->SetLabelSize(TMath::Power(3,0.5)*style->GetLabelSize("Z"),"Z");  
  fStyleAlice["figTemplate3"]=style;



}

void  AliDrawStyle::RegisterDefaultMarkers(){
  //
  // Style source:
  // https://twiki.cern.ch/twiki/pub/ALICE/ALICERecommendationsResultPresentationText/figTemplate.C
  const Int_t fillColors[] = {kGray+1,  kRed-10, kBlue-9, kGreen-8, kMagenta-9, kOrange-9,kCyan-8,kYellow-7, kBlack, kRed+1 }; // for syst bands
  const Int_t colors[]     = {kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2,kGray+1,  kRed-10 };
  const Int_t markers[]    = {kFullCircle, kFullSquare,kOpenCircle,kOpenSquare,kOpenDiamond,kOpenCross,kFullCross,kFullDiamond,kFullStar,kOpenStar};
  //
  (fMarkerStyles["figTemplate"])=std::vector<int>(10);
  (fMarkerColors["figTemplate"])=std::vector<int>(10);
  (fMarkerSize["figTemplate"])=std::vector<float>(10);
  (fFillColors["figTemplate"])=std::vector<int>(10);
  (fLineWidth["figTemplate"])=std::vector<float>(10);
  for (Int_t i=0;i<10; i++){
    (fMarkerStyles["figTemplate"])[i]=markers[i];
    (fMarkerColors["figTemplate"])[i]=colors[i];
    (fMarkerSize["figTemplate"])[i]=1;
    (fFillColors["figTemplate"])[i]=fillColors[i];
    (fLineWidth["figTemplate"])[i]=0.5;
  }
  // style inspired by TRD performance paper
  Int_t colorsTRD[12]={0};
  const Int_t markersTRD[]    = {kOpenCircle,kFullCircle, kOpenSquare,kFullSquare, kOpenStar,kFullStar, kOpenDiamond,kFullDiamond, kOpenCross,kFullCross };
  const Float_t markerTRDSize[]    = {1,1, 0.9,0.9, 1.4,1.4, 1.1,1.1, 1.2,1.2 };
  colorsTRD[0]=TColor::GetColor("#0000DD");
  colorsTRD[1]=TColor::GetColor("#00EE00");
  colorsTRD[2]=TColor::GetColor("#FF0000");
  colorsTRD[3]=TColor::GetColor("#00EEDD");
  colorsTRD[4]=TColor::GetColor("#FFEE00");
  colorsTRD[5]=TColor::GetColor("#FF00DD");
  colorsTRD[6]=TColor::GetColor("#9999DD");
  colorsTRD[7]=TColor::GetColor("#99EE99");
  colorsTRD[8]=TColor::GetColor("#FF9999");
  colorsTRD[9]=TColor::GetColor("#66AADD");
  colorsTRD[10]=TColor::GetColor("#AAEE66");
  colorsTRD[11]=TColor::GetColor("#FF66AA");
  (fMarkerStyles["figTemplateTRD"])=std::vector<int>(10);
  (fMarkerColors["figTemplateTRD"])=std::vector<int>(10);
  (fMarkerSize["figTemplateTRD"])=std::vector<float>(10);
  (fFillColors["figTemplateTRD"])=std::vector<int>(10);
  (fLineWidth["figTemplateTRD"])=std::vector<float>(10);
  (fMarkerStyles["figTemplateTRDPair"])=std::vector<int>(10);
  (fMarkerColors["figTemplateTRDPair"])=std::vector<int>(10);
  (fMarkerSize["figTemplateTRDPair"])=std::vector<float>(10);
  (fFillColors["figTemplateTRDPair"])=std::vector<int>(10);
  (fLineWidth["figTemplateTRDPair"])=std::vector<float>(10);

  for (Int_t i=0; i<10; i++){
    (fMarkerStyles["figTemplateTRD"])[i]=markersTRD[i];
    (fMarkerColors["figTemplateTRD"])[i]=TColor::GetColorDark(colorsTRD[i]);
    (fMarkerSize["figTemplateTRD"])[i]=markerTRDSize[i];
    (fFillColors["figTemplateTRD"])[i]=fillColors[i];
    (fLineWidth["figTemplateTRD"])[i]=0.5;
    //
    (fMarkerStyles["figTemplateTRDPair"])[i]=markersTRD[i];
    (fMarkerColors["figTemplateTRDPair"])[i]=TColor::GetColorDark(colorsTRD[i/2]);
    (fMarkerSize["figTemplateTRDPair"])[i]=markerTRDSize[i];
    (fFillColors["figTemplateTRDPair"])[i]=fillColors[i/2];
    (fLineWidth["figTemplateTRDPair"])[i]=0.5;
  }

}

TStyle*  RegisterDefaultStyleFigTemplate(Bool_t graypalette) {
  // Style source:
  // https://twiki.cern.ch/twiki/pub/ALICE/ALICERecommendationsResultPresentationText/figTemplate.C
  //
  TStyle * figStyle = new TStyle;
  figStyle->Reset("Plain");
  figStyle->SetOptTitle(0);
  figStyle->SetOptStat(0);

  if(graypalette) figStyle->SetPalette(8,0);
  else figStyle->SetPalette(1);

  figStyle->SetCanvasColor(10);
  figStyle->SetCanvasBorderMode(0);
  figStyle->SetFrameLineWidth(1);
  figStyle->SetFrameFillColor(kWhite);
  figStyle->SetPadColor(10);
  figStyle->SetPadTickX(1);
  figStyle->SetPadTickY(1);
  figStyle->SetPadBottomMargin(0.15);
  figStyle->SetPadLeftMargin(0.15);
  figStyle->SetHistLineWidth(1);
  figStyle->SetHistLineColor(kRed);
  figStyle->SetFuncWidth(2);
  figStyle->SetFuncColor(kGreen);
  figStyle->SetLineWidth(2);
  figStyle->SetLabelSize(0.045,"xyz");
  figStyle->SetLabelOffset(0.01,"y");
  figStyle->SetLabelOffset(0.01,"x");
  figStyle->SetLabelColor(kBlack,"xyz");
  figStyle->SetTitleSize(0.05,"xyz");
  figStyle->SetTitleOffset(1.25,"y");
  figStyle->SetTitleOffset(1.2,"x");
  figStyle->SetTitleFillColor(kWhite);
  figStyle->SetTextSizePixels(26);
  figStyle->SetTextFont(42);
  figStyle->SetLegendBorderSize(0);
  figStyle->SetLegendFillColor(kWhite);
  figStyle->SetLegendFont(42);
  return figStyle;
}
