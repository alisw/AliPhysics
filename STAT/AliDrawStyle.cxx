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
///  * Several drawing styles can be regeistered and used in the same moment
///    * Styles are identified using strings as identifiets
///      * TStyle
///      * MarkerStyle[]  AliDrawStyle::GetMarkerStyle(const char *style, Int_t index);
///      * MarkerColors[] AliDrawStyle::GetMarkerColor(const char *style, Int_t index);
///      * FillColors[]   AliDrawStyle::GetFillColor(const char *style, Int_t index);
///  * Default styles are created  AliDrawStyle::SetDefaults()
///    * default style is based on the fig template -  https://twiki.cern.ch/twiki/pub/ALICE/ALICERecommendationsResultPresentationText/figTemplate.C
///    * users should be able to regiester their oun styles (e.g in macros)
///  * Usage (work in progress)
///    * performance reports -  with styles as a parameter
///    * QA reports
///    * AliTreePlayer, and TStatToolkit
/// \author marian  Ivanov marian.ivanov@cen.ch
///
///  ## Example usage
///
///  \code
///  AliDrawStyle::SetDefaults()
///  // Style example
///  //
///  AliDrawStyle::PrintStyles(0,TPRegexp("."));
///  AliDrawStyle::ApplyStyle("figTemplate");
///  //
///  // Standard ALICE latex symbols
///  AliDrawStyle::PrintLatexSymbols(0,TPRegexp("."))
///  AliDrawStyle::GetLatexAlice("qPt")
///  AliDrawStyle::AddLatexSymbol("dphi", "#Delta#it#phi (unit)")
///  AliDrawStyle::GetLatexAlice("dphi")
///  //
///  // Standard ALICE marker/colors arrays
///  AliDrawStyle::GetMarkerStyle("figTemplate",0)
///  AliDrawStyle::GetMarkerColor("figTemplate",0)
///  \endcode



#include "AliDrawStyle.h"
#include "TStyle.h"
#include "TError.h"
#include "TPRegexp.h"
#include "TColor.h"
#include <iostream>
//
std::map<TString, TString>  AliDrawStyle::fLatexAlice;
std::map<TString, TStyle*>  AliDrawStyle::fStyleAlice;
std::map<TString, std::vector<int> > AliDrawStyle::fMarkerStyles;
std::map<TString, std::vector<int> > AliDrawStyle::fMarkerColors;
std::map<TString, std::vector<int> > AliDrawStyle::fFillColors;

void AliDrawStyle::SetDefaults(){
  AliDrawStyle::RegisterDefaultLatexSymbols();
  AliDrawStyle::RegisterDefaultStyle();
  AliDrawStyle::RegisterDefaultMarkers();
}


TStyle* RegisterDefaultStyleFigTemplate(Bool_t grayScale);

/// Latex symbol section

/// \param symbol  -id name of latex symbol
/// \return latex symbol to be used in ROOT latex
TString AliDrawStyle::GetLatexAlice(const char * symbol){
  return  fLatexAlice[symbol];
}

/// \param  style - name of style used
/// \param index  - marker index
/// \return marker style for given stylename, index
Int_t AliDrawStyle::GetMarkerStyle(const char *style, Int_t index){
  return  AliDrawStyle::fMarkerStyles[style][index];
}

/// \param  style - name of style used
/// \param index  - marker index
/// \return marker color for given stylename, index
Int_t AliDrawStyle::GetMarkerColor(const char *style, Int_t index){
  return  AliDrawStyle::fMarkerColors[style][index];
}

/// \param  style - name of style used
/// \param index  - marker index
/// \return fill color for given stylename, index
Int_t AliDrawStyle::GetFillColor(const char *style, Int_t index){
  return  AliDrawStyle::fFillColors[style][index];
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
  fLatexAlice["qPt"]="#it{p}_{T} (GeV/#it{c})";
  fLatexAlice["Pt"]="#it{p}_{T}";
  fLatexAlice["sqPtMev"]="#sigma_{#it{q}/#it{p}_{T}}/#it{p}_{T}^{2} (MeV/#it{c})^{-1}";
  fLatexAlice["PbPb502"]="Pb#font[122]{-}Pb #sqrt{#it{s}_{NN}} =5.02 TeV";
  fLatexAlice["pp13"]="pp #sqrt{#it{s}} = 13 TeV ";
  fLatexAlice["drphi"]="#Delta_{#it{r#phi}} (cm)";
  fLatexAlice["srphi"]="#sigma_{#it{r#phi}} (cm)";
}

void   AliDrawStyle::RegisterDefaultStyle(){
  //
  fStyleAlice["figTemplate"]=RegisterDefaultStyleFigTemplate(kFALSE);
  fStyleAlice["figTemplateGrey"]=RegisterDefaultStyleFigTemplate(kFALSE);

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
  (fFillColors["figTemplate"])=std::vector<int>(10);
  (fMarkerStyles["figTemplatePair"])=std::vector<int>(10);
  (fMarkerColors["figTemplatePair"])=std::vector<int>(10);
  (fFillColors["figTemplatePair"])=std::vector<int>(10);
  for (Int_t i=0;i<10; i++){
    (fMarkerStyles["figTemplate"])[i]=markers[i];
    (fMarkerColors["figTemplate"])[i]=colors[i];
    (fFillColors["figTemplate"])[i]=fillColors[i];
    (fMarkerStyles["figTemplatePair"])[i]=markers[i];
    (fMarkerColors["figTemplatePair"])[i]=colors[i/2];
    (fFillColors["figTemplatePair"])[i]=fillColors[i/2];
    //
    (fMarkerStyles["figTemplateDark"])[i]=TColor::GetColorDark(markers[i]);
    (fMarkerColors["figTemplateDark"])[i]=TColor::GetColorDark(colors[i]);
    (fFillColors["figTemplateDark"])[i]=TColor::GetColorDark(fillColors[i/2]);
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
