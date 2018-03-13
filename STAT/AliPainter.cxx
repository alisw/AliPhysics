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
// TODO - ADD TESTS. check few properties | priority: most higher @Boris
// TODO: prepare presentation for offline week: | priority: high @Boris
//      -- AliDrawStyle
//      -- AliPainter
//      -- AliTMinuitToolkit
#include "AliPainter.h"
#include "AliTMinuitToolkit.h"
#include "TPad.h"
#include "TList.h"
#include "TAxis.h"
#include "TPRegexp.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include "TError.h"
#include "THn.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <vector>
#include <map>
#include <bitset>

//typedef std::map<int, std::vector<int>() > axisRangesMap;
std::map<TString, TString> AliPainter::optionValues;
std::map<TString, Double_t > AliPainter::statValues;

///
/// \brief Method allow to divide pad according to specify properties.
///
///
/// \param pad           - input pad to divide
/// \param division      - division string
///                         <NULL,vertical,horizontal>[div0,div1b, ...]
///                            divi - specify number of pads in row (resp. column)
///                                 - sharing parameter for axis
///                                 - btlrm  (bottom, left, top, right, middle) middle means rl for horizontal and tb for vertical
///                                 - set margin0 in case specified
///                                 - technically attribute can be added to the object name  axis-sharing=""
///                                 - 1lpx30 - means set for pad left margin equal 30pixels
/// \param classID        - adds classID to name of the pad
/*!
  #### Example use:

  1. Let's add some tree with Histogramm:
        For this you can use
        \code
          AliExternalInfo info("","",0);
          tree = info.GetTree("QA.TPC","LHC15o","pass1");
         \endcode
  2. Create canvas and divide it:
          \code
            TCanvas * canvasQA = new TCanvas("canvasQA","canvasQA", 1200,800);
            AliPainter::DivideTPad(canvasQA,"<horizontal>[1,1,1,1]", "Raw,Error");
            canvasQA->cd(1);
            TStatToolkit::DrawMultiGraph(TStatToolkit::MakeMultGraph(tree,"", "meanMIPele:run","1","25", "1",1,1,5,0),"ap");
            canvasQA->cd(2);
            TStatToolkit::DrawMultiGraph(TStatToolkit::MakeMultGraph(tree,"", "meanMIPele:run","1","25", "1",1,1,5,0),"ap");
            canvasQA->cd(3);
            TStatToolkit::DrawMultiGraph(TStatToolkit::MakeMultGraph(tree,"", "meanMIPele:meanMIP","1","25", "1",0,1,5,0),"ap");
            canvasQA->cd(4);
            TStatToolkit::DrawMultiGraph(TStatToolkit::MakeMultGraph(tree,"", "meanMIPele:meanMIP","1","25", "1",0,1,5,0),"ap");
         \endcode
         In this simple case we have four vertical pads:
         \image html $AliRoot_SRC/doxygen/canvasQA4v.jpg
         In case with horizontal position we will have next plot:
         \image html $AliRoot_SRC/doxygen/canvasQA4h.jpg
  3. You can sharing choosen axis (see meanings of flags in description.):
        \code
            TCanvas * canvasQA = new TCanvas("canvasQA","canvasQA", 1200,800);
            AliPainter::DivideTPad(canvasQA,"<vertical>[1r,1l,1r,1l]", "Raw,Error");
            ...
         \endcode
         \image html $AliRoot_SRC/doxygen/canvasQA4hbt.jpg

       and for vertical
       \code
            TCanvas * canvasQA = new TCanvas("canvasQA","canvasQA", 1200,800);
            AliPainter::DivideTPad(canvasQA,"<vertical>[1r,1l,1r,1l]", "Raw,Error");
            ...
       \endcode
        \image html $AliRoot_SRC/doxygen/canvasQA4vrl.jpg
    \code
  4. You can specify more than one pads in one raw:
         \code
            TCanvas * canvasQA = new TCanvas("canvasQA","canvasQA", 1200,800);
            AliPainter::DivideTPad(canvasQA,"<horizontal>[1,3,2,1]", "Raw,Error");
            canvasQA->cd(1);
            TStatToolkit::DrawMultiGraph(TStatToolkit::MakeMultGraph(tree,"", "meanMIPele:run","1","25", "1",1,1,5,0),"ap");
            canvasQA->cd(2);
            TStatToolkit::DrawMultiGraph(TStatToolkit::MakeMultGraph(tree,"", "meanMIPele:run","1","25", "1",1,1,5,0),"ap");
            canvasQA->cd(3);
            TStatToolkit::DrawMultiGraph(TStatToolkit::MakeMultGraph(tree,"", "meanMIPele:meanMIP","1","25", "1",0,1,5,0),"ap");
            canvasQA->cd(4);
            TStatToolkit::DrawMultiGraph(TStatToolkit::MakeMultGraph(tree,"", "meanMIPele:meanMIP","1","25", "1",0,1,5,0),"ap");
            canvasQA->cd(5);
            TStatToolkit::DrawMultiGraph(TStatToolkit::MakeMultGraph(tree,"", "meanMIPele:meanMIP","1","25", "1",0,1,5,0),"ap");
            canvasQA->cd(6);
            TStatToolkit::DrawMultiGraph(TStatToolkit::MakeMultGraph(tree,"", "meanMIPele:meanMIP","1","25", "1",0,1,5,0),"ap");
            canvasQA->cd(7);
            TStatToolkit::DrawMultiGraph(TStatToolkit::MakeMultGraph(tree,"", "meanMIPele:meanMIP","1","25", "1",0,1,5,0),"ap");
        \endcode
        \image html $AliRoot_SRC/doxygen/canvasQA13m2m1h.jpg
    or column:
        \code
            TCanvas * canvasQA = new TCanvas("canvasQA","canvasQA", 1200,800);
            AliPainter::DivideTPad(canvasQA,"<vertical>[1,3m,2m,1]", "Raw,Error");
             ...
        \endcode
        \image html $AliRoot_SRC/doxygen/canvasQA1321v.jpg

  5. You can specify special flag "m" if you want to join middle pads in column:
        \code
            TCanvas * canvasQA = new TCanvas("canvasQA","canvasQA", 1200,800);
            AliPainter::DivideTPad(canvasQA,"<vertical>[1,3m,2m,1]", "Raw,Error");
             ...
        \endcode
        \image html $AliRoot_SRC/doxygen/canvasQA13m2m1v.jpg

     or row:
             \code
            TCanvas * canvasQA = new TCanvas("canvasQA","canvasQA", 1200,800);
            AliPainter::DivideTPad(canvasQA,"<horizontal>[1,3m,2m,1]", "Raw,Error");
             ...
        \endcode
        \image html $AliRoot_SRC/doxygen/canvasQA13m2m1h.jpg

  6. All previous case set choosen margin equal 0, but you can set your own value in absolute units, now we are supporting
     only pixel("px" flag) and standart("st" flag) root values (percent from canvas size):
        \code
            TCanvas * canvasQA = new TCanvas("canvasQA","canvasQA", 1200,800);
            AliPainter::DivideTPad(canvasQA,"<horizontal>[1rpx200,1,1,1]", "Raw,Error");
             ...

            TCanvas * canvasQA = new TCanvas("canvasQA","canvasQA", 1200,800);
            AliPainter::DivideTPad(canvasQA,"<horizontal>[1rst0.2,1,1,1]", "Raw,Error");
             ...
        \endcode

//    TFile*f=new TFile("/Users/bdrum/Projects/alicesw/TPC_trending.root");
//    TTree*tree=(TTree*)f.Get("trending");



    canvasQA->cd(1);
    TStatToolkit::DrawMultiGraph(TStatToolkit::MakeMultGraph(tree,"", "meanMIPele:run","1","25", "1",1,1,5,0),"ap");
     canvasQA->cd(2);
    TStatToolkit::DrawMultiGraph(TStatToolkit::MakeMultGraph(tree,"", "meanMIPele:run","1","25", "1",1,1,5,0),"ap");
    canvasQA->cd(3);
    TStatToolkit::DrawMultiGraph(TStatToolkit::MakeMultGraph(tree,"", "meanMIPele:meanMIP","1","25", "1",0,1,5,0),"ap");
    canvasQA->cd(4);
    TStatToolkit::DrawMultiGraph(TStatToolkit::MakeMultGraph(tree,"", "meanMIPele:meanMIP","1","25", "1",0,1,5,0),"ap");
    canvasQA->cd(5);
    TStatToolkit::DrawMultiGraph(TStatToolkit::MakeMultGraph(tree,"", "meanMIPele:meanMIP","1","25", "1",0,1,5,0),"ap");
    canvasQA->cd(6);
    TStatToolkit::DrawMultiGraph(TStatToolkit::MakeMultGraph(tree,"", "meanMIPele:meanMIP","1","25", "1",0,1,5,0),"ap");
    canvasQA->cd(7);
    TStatToolkit::DrawMultiGraph(TStatToolkit::MakeMultGraph(tree,"", "meanMIPele:meanMIP","1","25", "1",0,1,5,0),"ap");
    canvasQA->cd(8);
    TStatToolkit::DrawMultiGraph(TStatToolkit::MakeMultGraph(tree,"", "meanMIPele:meanMIP","1","25", "1",0,1,5,0),"ap");
    canvasQA->cd(9);
    TStatToolkit::DrawMultiGraph(TStatToolkit::MakeMultGraph(tree,"", "meanMIPele:meanMIP","1","25", "1",0,1,5,0),"ap");
    canvasQA->cd(10);
    TStatToolkit::DrawMultiGraph(TStatToolkit::MakeMultGraph(tree,"", "meanMIPele:meanMIP","1","25", "1",0,1,5,0),"ap");
    canvasQA->cd(11);
    TStatToolkit::DrawMultiGraph(TStatToolkit::MakeMultGraph(tree,"", "meanMIPele:meanMIP","1","25", "1",0,1,5,0),"ap");
    canvasQA->cd(12);
    TStatToolkit::DrawMultiGraph(TStatToolkit::MakeMultGraph(tree,"", "meanMIPele:meanMIP","1","25", "1",0,1,5,0),"ap");
    \endcode
*/
//TODO: add divFlag for inheritance Mother Class (ClassName) should be add to children (className)=>nameOfObject.class(ClassMother,ClassDaughter) @Boris
//TODO: should we change the places of *pad and *classID? @Marian
//TODO: mb remove "<>" from division? @Marian
TPad *AliPainter::DivideTPad(const char *division, const char *classID, TPad *pad, Int_t verbose) {

  if (pad == nullptr) {
    TCanvas *canv = new TCanvas("AliPainterCanvas", "AliPainterCanvas",1200,800);
    ::Info("AliPainter::DivideTPad", "Input pad is null. We created new one.  pad = new TPad(\"AliPainter\", \"AliPainter\",0,0,1,1);");
    pad = canv;
  }

  Int_t    nPads    = 0, nRows = 0, maxNCols = 0;
  Double_t xl, yl, xu, yu, a, b, mValue;
  TString  position = "";
  TString  tempStr  = "";
  TString  wMargin  = ""; // margin zero
  TString  units    = ""; // set margin
  TString  padName  = "";

  TObjArray *padRows = TString(division).Tokenize("<>[],");
  position = padRows->At(0)->GetName();
  nRows = padRows->GetEntries() - 1;
  if (verbose == 4) ::Info("AliPainter::DivideTPad", "Number of rows is %d", nRows);

  for (Int_t iRow = 0; iRow < nRows; iRow++)
    if (maxNCols < TString(padRows->At(iRow + 1)->GetName()).Atoi())
      maxNCols = TString(padRows->At(iRow + 1)->GetName()).Atoi();

  for (Int_t iRow = 0; iRow < nRows; iRow++) {

    tempStr = TString(padRows->At(iRow + 1)->GetName());
    Int_t nCols = TString(tempStr(0, 1)).Atoi();
    if (verbose == 4) ::Info("AliPainter::DivideTPad", "Number of columns in %d row is %d", nRows, nCols);
    wMargin  = TString(tempStr(1, 1));
    units = TString(tempStr(2, 2));
    mValue = TString(tempStr(4, tempStr.Length())).Atof();

    for (Int_t iCol = 0; iCol < nCols; iCol++) {
      pad->cd();
      xl = iCol / Double_t(nCols);
      yl = (nRows - iRow - 1) / Double_t(nRows);
      xu = (iCol + 1) / Double_t(nCols);
      yu = (nRows - iRow) / Double_t(nRows);

      if (position == "vertical") {
        a = xl;
        b = yl;
        yl = a;       //yl -> xl
        xl = 1 - b;   //xl -> 1 - yl
        a = xu;
        b = yu;
        yu = a;       //yu -> xu
        xu = 1 - b;   //xu -> 1 - yu
      }
      if (classID != TString("")) padName = TString::Format("pad[%d].class(%s)", nPads, classID);
      else padName = TString::Format("pad[%d]", nPads);
      TPad *newPad = new TPad(padName.Data(), padName.Data(), xl, yl, xu, yu);
      if (verbose == 4) ::Info("AliPainter::DivideTPad", "New pad created: %s", padName.Data());
      if (position == "vertical") {
        newPad->SetTopMargin(nCols * newPad->GetLeftMargin() / maxNCols);
        newPad->SetBottomMargin(nCols * newPad->GetRightMargin() / maxNCols);
      }
      else {
        newPad->SetLeftMargin(nCols * newPad->GetLeftMargin() / maxNCols);
        newPad->SetRightMargin(nCols * newPad->GetRightMargin() / maxNCols);
      }
      newPad = AliPainter::SetPadMargin(newPad, position, wMargin, units, mValue, iCol, nCols);
      newPad->Draw();
      nPads++;
      newPad->SetNumber(nPads);
    }
  }
  if (classID != TString("")) padName = TString::Format("%s.class(%s)",pad->GetName(),classID);
  else padName = TString::Format("%s",pad->GetName());
  pad->SetName(padName);
  return pad;
}
///
/// \brief Function parses division string from AliPainter::DivideTPad and sets attributes.
/// \param cPad
/// \param position
/// \param units
/// \param value
TPad *AliPainter::SetPadMargin(TPad *cPad, const char *position, const char *wMargin, const char *units, Double_t mValue, Int_t iCol, Int_t nCols) {
  Float_t rlValue  = 0.0, btValue = 0.0;
  if (TString(units) == "px") {
    rlValue = mValue / cPad->GetWw();
    btValue = mValue / cPad->GetWh();
  }
  else if (TString(units) == "st") {
    rlValue = mValue;
    btValue = mValue;
  }
  else {
    rlValue = 0.0;
    btValue = 0.0;
  }

  if (TString(wMargin) == "r") cPad->SetRightMargin(rlValue);
  if (TString(wMargin) == "l") cPad->SetLeftMargin(rlValue);
  if (TString(wMargin) == "t") cPad->SetTopMargin(btValue);
  if (TString(wMargin) == "b") cPad->SetBottomMargin(btValue);

  if (TString(position) == "vertical") {
    if (nCols > 1 && TString(wMargin) == "m" && iCol == 0) cPad->SetTopMargin(btValue);  //top pad in the m
    if (nCols > 1 && TString(wMargin) == "m" && iCol > 0 && (iCol + 1) < nCols) { // middle pads in the m
      cPad->SetTopMargin(btValue);
      cPad->SetBottomMargin(btValue);
    }
    if (nCols > 1 && TString(wMargin) == "m" && (iCol + 1) == nCols) cPad->SetBottomMargin(btValue); //bottom in the m
    if (TString(wMargin) == "m" && nCols == 1) { // only 1 with m
      cPad->SetLeftMargin(rlValue);
      cPad->SetRightMargin(rlValue);
    }
  }
  else { //the same for horizontal
    if (nCols > 1 && TString(wMargin) == "m" && iCol == 0) cPad->SetRightMargin(rlValue);
    if (nCols > 1 && TString(wMargin) == "m" && iCol > 0 && (iCol + 1) < nCols) {
      cPad->SetRightMargin(0);
      cPad->SetLeftMargin(0);
    }
    if (nCols > 1 && TString(wMargin) == "m" && (iCol + 1) == nCols) cPad->SetLeftMargin(rlValue);
    if (TString(wMargin) == "m" && nCols == 1) {
      cPad->SetTopMargin(btValue);
      cPad->SetBottomMargin(btValue);
    }
  }
  return cPad;
}
///
/// \param graph
/// \param option
void AliPainter::SetMultiGraphTimeAxis(TMultiGraph *graph, TString option){
  TAxis *axis = NULL;
  for (Int_t i=0; i<graph->GetListOfGraphs()->GetEntries(); i++) {
    TGraph *cGraph=(TGraph *) graph->GetListOfGraphs()->At(i);
    if (option.Contains("X")) axis = cGraph->GetXaxis();
    if (option.Contains("Y")) axis = cGraph->GetYaxis();
    if (axis) {
      axis->SetNdivisions(510, kFALSE);
      axis->SetTimeDisplay(1);
      axis->SetTimeFormat("%d/%m");
    }
  }
}
// TODO: check perfomance and change symbols to character constant @Boris
/// \brief Private method for parsing arguments in AliPainter::DrawHistogram
/// \param exprsn - string with arguments
/// \return - vector of arguments
void AliPainter::ArgsParser(TString exprsn, TString &hisName, TString &projections, std::vector<TString> &fitOptions, std::vector<TString> &rangesStrings, Int_t verbose) {
  //check for match of brackets
  if (exprsn.CountChar('(') != exprsn.CountChar(')')) {
    ::Error("AliPainter::DrawHistogram", "check brackets in %s", exprsn.Data());
    return;
  }
  //save each argument from input expression into vector of string
  std::vector<TString> atts;
  TString verbStr = "";
  Int_t match = 0, startIndex = 0, finishIndex = 0;
  Bool_t isChange = kFALSE;
  for (Int_t i = 0; i < exprsn.Length(); i++) {
    if (exprsn(i) == '(' && match == 0) {
      match++;
      startIndex = i;
      isChange = kTRUE;
    } else if (exprsn(i) == '(' && match > 0) match++;
    else if (exprsn(i) == ')' && match == 1) {
      match--;
      finishIndex = i;
    } else if (exprsn(i) == ')' && match > 1) match--;
    if (match == 0 && isChange) atts.push_back(TString(exprsn(startIndex + 1, finishIndex - startIndex - 1)));
  }
  if (verbose == 4) {
    verbStr = "";
    for (Int_t i = 0; i < (Int_t) atts.size(); i++)
      verbStr += atts[i] + "|";
    ::Info("AliPainter::DrawHistogram", "Input expression - %s was transform into array %s", exprsn.Data(),
           verbStr.Data());
  }
  //get name of histogram
  hisName = exprsn(0, exprsn.Index("(", 0));
  //get string of projections
  projections = atts[1];
  //get options for fitting
  fitOptions = AliPainter::OptionStringParser(atts[2], "()", 6);
  if (verbose == 4 && atts[2] != "") {
    verbStr = "";
    for (Int_t i = 0; i < fitOptions.size(); i++)
      verbStr += fitOptions[i] + "|";
    ::Info("AliPainter::DrawHistogram", "Input fit option - %s was transform into array %s", atts[2].Data(),
           verbStr.Data());
  }
  //get ranges for hist looping
  if (atts[0].CountChar(',') > 5) {
    ::Error("AliPainter::DrawHistogram()", "AliPainter::ArgsParser error. rangesString has more parameters then expected.");
    return;
  }
  if (atts[0] != "") rangesStrings = AliPainter::RangesParser(atts[0]);
  if (verbose == 4 && atts[0] != "") {
    verbStr = "";
    for (Int_t i = 0; i < (Int_t) rangesStrings.size(); i++)
      verbStr += rangesStrings[i] + "|";
    ::Info("AliPainter::DrawHistogram", "Input ranges option - %s was transform into array %s", atts[0].Data(),
           verbStr.Data());
  }
  //get drawOption
  std::vector<TString> drawOptions = AliPainter::OptionStringParser(atts[3], "[]", 6);
  if (verbose == 4 && atts[3] != TString("")) {
    AliPainter::RegisterDefaultOptions();
    verbStr = "";
    for (Int_t i = 0; i < drawOptions.size(); i++) {
      if (drawOptions[i] != TString(""))
        AliPainter::PandasOptionParser(drawOptions[i]);
      verbStr += drawOptions[i] + "|";
    }
    ::Info("AliPainter::DrawHistogram", "Input draw option - %s was transform into array %s", atts[3].Data(),
           verbStr.Data());
  }
}
/// \brief Private method for parsing fit options in AliPainter::DrawHistogram
/// \param fitStr - string with fit options
/// \return array of values from inputOptions
std::vector<TString> AliPainter::OptionStringParser(const TString optStr, const char d[2], Int_t defSize) {
  Int_t arg = 0, startIndex = 0;
  std::vector<TString> vecOptions;
  for (Int_t i = 0; i < defSize+1; i++) vecOptions.push_back(TString(""));
  if (optStr == TString("")) {
    ::Error("AliPainter", "AliPainter::OptionStringParser(%s,%s,%d). Options string should not be empty.", optStr.Data(), d, defSize);
    return vecOptions;
  }
  for (Int_t i = 0; i <= optStr.Length(); i++) {
    if (arg > defSize) break;
    if (optStr(i) == ',' || i == optStr.Length()) {
      vecOptions[arg] = TString(optStr(startIndex, i - startIndex));
      arg++;
      startIndex = i + 1;
    } else if (optStr(i) == d[0]) {
      i = optStr.Index(d[1], i);
      continue;
    }
  }
  return vecOptions;
}

void AliPainter::RegisterDefaultOptions() {
  //TODO: perhaps, we should use some aliases? @Marian
  //basic option
  AliPainter::optionValues[TString("hisName")]       = TString("");
  AliPainter::optionValues[TString("x0Ranges")]      = TString("");
  AliPainter::optionValues[TString("x1Ranges")]      = TString("");
  AliPainter::optionValues[TString("x2Ranges")]      = TString("");
  AliPainter::optionValues[TString("proj")]          = TString("");
  //fit options
  AliPainter::optionValues[TString("fitName")]       = TString("");
  AliPainter::optionValues[TString("fitStrategy")]   = TString("");
  AliPainter::optionValues[TString("fitOption")]     = TString("");
  AliPainter::optionValues[TString("fitRanges")]     = TString("");
  AliPainter::optionValues[TString("fitInitParams")] = TString("");
  AliPainter::optionValues[TString("fitDOption")]    = TString("");
  //draw options
  AliPainter::optionValues[TString("div")]           = TString("");
  AliPainter::optionValues[TString("xlim")]          = TString("");
  AliPainter::optionValues[TString("ylim")]          = TString("");
  AliPainter::optionValues[TString("zlim")]          = TString("");
//  AliPainter::optionValues[TString("max")]           = TString("");
//  AliPainter::optionValues[TString("min")]           = TString("");
//  AliPainter::optionValues[TString("mean")]          = TString("");
//  AliPainter::optionValues[TString("rms")]           = TString("");
//  AliPainter::optionValues[TString("median")]        = TString("");
  AliPainter::optionValues[TString("class")]         = TString("");
  AliPainter::optionValues[TString("dOption")]       = TString("");

}
/// \brief Private method for parsing draw options in AliPainter::DrawHistogram
/// \param drawStr - string with draw options
/// \return array of values from inputOptions
void AliPainter::PandasOptionParser(const TString optionsStr) {
  TString key = "";
  TString value = "";
  key = TString(optionsStr(0, optionsStr.Index("="))).ReplaceAll(" ", "");
  value = TString(optionsStr(optionsStr.Index("=") + 1, optionsStr.Length())).ReplaceAll(" ", "");
  if (optionValues.find(key) == optionValues.end()) {
    ::Error("AliPainter::DrawHistogram", "key %s not found in default keys: hisName,x0Ranges,x1Ranges,x2Ranges,proj,fitName,fitStrategy,fitOption,fitRanges,fitInitParams,fitDOption,div,xlim,ylim,zlim,class,dOption", key.Data());
    return;
  }
  optionValues[key] = value;
}
// TODO - now we have 2 steps for parsing of ranges option to array of string with simple range values. 1. Initial string into map. 2. Map into array fo strings. First of all we can use vector of vector instead map, and the second we should think how we can avoid this 2 intermediate steps with map. @Boris
/// \brief Private method for parsing ranges options in AliPainter::DrawHistogram
/// \param ranges - String with ranges values
/// \return - map, where key is projection and values from ranges
std::vector<TString> AliPainter::RangesParser(TString ranges) {
  if (ranges == "") return std::vector<TString>();
  TString range = "";
  Ssiz_t fromStart0 = 0;
  Int_t axisNum = 0, numCnt = 1;
  axisRangesMap result;
  while (ranges.Tokenize(range, fromStart0, ",")) {
    if (range.IsFloat()) {
      if (numCnt % 2 == 0) axisNum--;
      numCnt++;
    }
    RangesToMap(range, axisNum, result);
    axisNum++;
  }
  std::vector<TString> arr;
  for (Int_t i = 0; i < (Int_t) result.size(); i++) {
    arr.push_back("");
  }
  std::vector<TString> res;
  AliPainter::RangesMapToString((Int_t) result.size(), result, arr, res);
  return res;
}
/// /brief Subsidiary method for RangesParser
/// \param n - number of dimensions
/// \param result - map with strings of axis ranges
/// \param arr - temprorary array
/// \param res - vector of string with values of ranges
void AliPainter::RangesMapToString(Int_t n, axisRangesMap result, std::vector<TString> arr, std::vector<TString> &res) {
  TString tempStr = "";
  if (n > 0) {
    for (Int_t i = 0; i < result[result.size() - n].size(); i += 2) {
      arr[result.size() - n] = TString::Format("%s,%s", result[result.size() - n][i].Data(), result[result.size() - n][i + 1].Data());
      RangesMapToString(n - 1, result, arr, res);
    }
  }
  else {
    for (Int_t j = 0; j < arr.size();j++) {
      tempStr += arr[j];
      if (j != arr.size() - 1) tempStr += ",";
    }
    res.push_back(tempStr);
  }
}
/// \brief Private method for filling result map (returned by  AliPainter::RangesParser())
/// \param range - concrete range value or python-like array
/// \param axisNum - num for parsing
/// \param result - output map
void AliPainter::RangesToMap(TString range, Int_t axisNum, axisRangesMap &result) {
  TString num = "";
  Ssiz_t fromStart1 = 0;
  Int_t i = 0;
  Int_t rangeArr[4] = {0,0,1,0}; // iLow, iHigh, step, delta
  std::vector<TString> vRanges;
  if (!range.IsFloat() && range.CountChar(':') > 0) {
    while (range.Tokenize(num, fromStart1, ":")) {
      rangeArr[i] = num.Atoi();
      i++;
    }
    for (Int_t j = rangeArr[0]; j <= rangeArr[1] - rangeArr[3]; j += rangeArr[2]) {
      vRanges.push_back(TString::Itoa(j,10));
      vRanges.push_back(TString::Itoa(j + rangeArr[3],10));
    }
    result[axisNum] = vRanges;
  }
  else if (range.IsFloat()) {
    result[axisNum].push_back(range);
  }
}
/// \brief Method finds histogram in inputArray and draw specified projection according to properties.
/// \param expresion        - histogram draw expression
///                         - syntax
///                           histogramName(<axisRanges>)(<projection string>)(<fitting string>)(<drawing string>)
///                             axisRange: @done
///                                if integer bin range   - TAxis::SetRange(min, max)
///                                if float   user range  - TAxis::SetRangeUser(min,max)
///                                if Expression is empty - do not specify anything
///                             projectionString: (i0,i1) @done
///                                new projection created THn his = hisInput->Projection(i0,i1....)
///                                at minimum one dimension should be specified, maximum 3D
///                             fitting string: (fitterName,fitOption,range,initialParam)
///                                fitterName - whatever fitter registered in the list of fitters
///                                             (defined in AliTMinuitToolkit , maybe also support
///                                              root fit functions)
///                                fitOption - see AliTMinuitToolkit fitOptions
///                                            we should put there checks of correctness of fit
///                                            options
///                                range - {x0min,x0max,x1min,xm1max,...}
///                                         in case not specified - range is not set
///                                intitialParam - {p0,p1;p2,...;ep0,ep1,...;minp0,minp1,...;
///                                                 maxp0,maxp1 ...} TODO: may be we should also use TVectorD or TMatrix instead enumeration? @Boris
///                                                 there are options- use only in case they are
///                                                 specified
///                                                 errors, min and max are optionals
///                             drawing string: (padDiv, lims, className, dOption)
///                                padDiv - allows to choose wich pad use for drawing:
///                                   divFlag = 0 - use the same pad for drawing; (default value)
///                                   divFlag = 1 - use differents pads for drawing;
///                                lims - allows to set limits for specified axis: (not specified by default)
///                                   xlim = [xmin, xmax] - set the minima and maxima for x-axes;
///                                   ylim = [ymin, ymax] - set the minima and maxima for y-axes;
///                                   zlim = [zmin, zmax] - set the minima and maxima for z-axes;
///                                className - adds name of class to each object of drawing. It need for applying css style - [AliDrawStyle:ApplyCssStyle();](link to docs):
///                                    class = [Raw,Error] - in the end of name of object will add ".class(Raw,Error)"
///                                    class = Raw - also as class=[Raw] will add .class(Raw)
///                                    class = [] - in this case nothing will add to the end of the name of object (default)
///                                dOption - root standard draw options [see docs](https://root.cern.ch/root/htmldoc/guides/users-guide/Histograms.html#draw-options)
/// \param inputArray       - array of input objects
///                         - Object to draw - histogramArray->FindObject(histogramName)
///                         - in case histogramArray is NULL or histogram not found gROOT->FindObject()
///                           will be used
/// \param metaData         - array with metadata describing histogram axis
///                         - for example in the trees we optionally keep metadata (array of TNamed ()tag,value) in the array "metaTable"
///                         - in case not specified -"metaTable" object from the histogramArray used
/// \param keepArray        - array to keep temporary objects
///
/// \return
/*!
   #### Example use:
    Data preparation:
    \code
    {
      TFile::SetCacheFileDir(".");
      TFile * finput = TFile::Open("http://aliqatrkeos.web.cern.ch/aliqatrkeos/performance/alice/data/2015/LHC15o/pass1/LHC15o.pass1.B1.Bin0/performanceHisto.root", "CACHEREAD");
      TTree* tree=(TTree*) finput.Get("hisPtAll");
      hisArray = new TObjArray();
      TList *keys = finput->GetListOfKeys();
      for (Int_t iKey = 0; iKey < keys->GetEntries(); iKey++) {
        TObject *o = finput->Get(TString::Format("%s;%d", keys->At(iKey)->GetName(), ((TKey *) keys->At(iKey))->GetCycle()).Data());
        hisArray->AddLast(o);
      }
    }
   \endcode
   Behaviour by default:
    \code
      AliPainter::DrawHistogram("hisK0DMassQPtTgl()(0)()()", hisArray);
    \endcode
    Using of ranges option:
    1. Simple values specified
      \code
        AliPainter::DrawHistogram("hisK0DMassQPtTgl(20,80)(0)()()", hisArray);
      \endcode
    2. Arrays of values specified
      2.1 Default behaviour
          By default draw use the same pad, but you can change it specify draw option.
        \code
          AliPainter::DrawHistogram("hisK0DMassQPtTgl([0,100], x1=[0:80:20:20],x2=[0,10])(0)()()", hisArray);
        \endcode
      2.2 Using AliPainter::DivideTPad();
          In case if you don't want to use the same pad you should specify it. Also we recommend to you use AliPainter::DivideTPad() for create a few pads on your canvas:
          Also here we add className(Raw) and you can see how it looks like via canvasQA->ls(). You can use this class name for applying your own styles in css file.
        \code
          TCanvas * canvasQA = new TCanvas("canvasQA","canvasQA", 1200,800);
          AliPainter::DivideTPad(canvasQA,"<horizontal>[1,1,1,1,1]", "Raw,Error");
          canvasQA->cd(1);
          AliPainter::DrawHistogram("hisK0DMassQPtTgl(0,100,0:80:20:20,0,10)(0)()(div=0,,class=[Raw,Error])", hisArray);
        \endcode
      2.3 Limitations to y-axis
          You can set the ranges to yaxis with using second parameter of drawOption:
        \code
          AliPainter::DrawHistogram("hisK0DMassQPtTgl(0,100,0:80:20:20,0,10)(0)()(div=0,ylim=[200,40000],class=Raw,dOption=E)", hisArray);
        \endcode
      2.4 Using of standard root draw options
        \code
          AliPainter::DrawHistogram("hisK0DMassQPtTgl(0,100,0:80:20:20,0,10)(0)()(div=0,ylim=[200,40000],class=Raw,dOption=E)", hisArray);
        \endcode
    3 Fitting
      3.1 Using of standard root fiting methods:
        \code
          AliPainter::DrawHistogram("hisK0DMassQPtTgl(0,100,0:80:20:20,0,10)(0)(gaus)(div=0,,class=Raw,dOption=E)", hisArray, canvasQA);
        \endcode
   \code
 // THn *hisn = hisArray->FindObject("hisK0DMassQPtTgl");
 // TH1 *his1 = hisn->Projection(1);
 // his1->Fit("gaus");
//  his1->Draw();
   //AliTMinuitToolkit::RegisterDefaultFitters()
 AliTMinuitToolkit *fitter = new AliTMinuitToolkit();
 TF1 *aFormExp = new TF1("formExp", "[0]*TMath::Exp(-[1]*x)");


 fitter->SetFitFunction(aFormExp, 0);
AliTMinuitToolkit::SetPredefinedFitter("ExpFit", fitter);
  TObjArray *keepArray = new TObjArray;
 // AliPainter::DivideTPad(canvasQA,"<horizontal>[1]", "Raw,Error");
 //AliPainter::DrawHistogram("hisK0DMassQPtTgl(0,100,0:60:20:20,0,10)(0)(pol8,,same,funOption(2,2,1),-0.01,0.01)()", hisArray, canvasQA, NULL, keepArray, 4);
 // AliPainter::DrawHistogram("hisK0DMassQPtTgl(0,100)(0)(gauss)()", hisArray, canvasQA, NULL, keepArray, 4);
   }
  \endcode
*/
//TODO: refactor the code @Boris
void AliPainter::DrawHistogram(char *expression, const TObjArray* histogramArray, TPad *pad, TObjArray *metaData, TObjArray *keepArray, Int_t verbose) {
  TString exprsn = expression;
  TString hisName = "";
  TString projections = "";
  TString uniqName = "";
  THn *hisN = nullptr;
  std::vector<TString> fitOptions;
  std::vector<TString> rangesString;

  //parsing arguments
  AliPainter::ArgsParser(exprsn, hisName, projections, fitOptions, rangesString);
  // check for existing of histogram
  hisN = (THn *) histogramArray->FindObject(hisName);
  if (hisN == nullptr)
    ::Info("AliPainter::DrawHistogram", "%s not found", (const char *) hisName);
  Int_t nDims = hisN->GetNdimensions();
  //checks  for number of dimensions
  if (nDims < projections.CountChar(',') + 1) {
    ::Error("AliPainter::DrawHistogram", "%s has only %d dimensions", (const char *) hisName, nDims);
    return;
  }
  else if (nDims > 3)
    ::Info("AliPainter::DrawHistogram", "You try to draw 4d histogram.");
  //check does the fitter exist if not register defaults:
  //TODO: how I can delete it if I use standart fitter. We should add clearing of fPredefinedFitters into AliTMinuitToolkit::ClearData() @Marian
  if (AliTMinuitToolkit::GetPredefinedFitter(fitOptions[0].Data()) == NULL) {
    AliTMinuitToolkit::RegisterDefaultFitters();
    if (verbose == 4)
      ::Info("AliPainter::DrawHistogram",
             "AliTMinuitToolkit::GetPredefinedFitter(%s) is NULL. AliTMinuitToolkit::RegisterDefaultFitters();",
             fitOptions[0].Data());
  }
  //pad->SetTitle(expression);
  Int_t histCnt = 1;
  Ssiz_t fromStart;
  Double_t xMin = fitOptions[4].Atof(), xMax = fitOptions[5].Atof();
  Double_t yMin, yMax;
  TString key;

  std::vector<TString> rangeVec;
  TString range = "";
  TString rangeString = "";
  TString fitStr = "";
  TString drawString = "";
  TH1D *hisArray[histCnt];
  TVirtualPad *cCanvas;
  TLegend *legend[histCnt];
  if (pad == nullptr && gPad == nullptr) {
    ::Error("AliPainter::DrawHistogram", "TPad object doesn't exist.");
    return;
  }
  if (pad == nullptr) cCanvas = gPad->GetMother();
  if (pad != nullptr && !pad->InheritsFrom("TCanvas")) cCanvas = pad->GetMother();
  TPad *nextPad;

  if (!rangesString.empty()) histCnt = (Int_t) rangesString.size();
  if (verbose == 4) ::Info("AliPainter::DrawHistogram", "Count of the histograms is %d", histCnt);
  //TODO: refactor it. We already have such strings. Let take it from ArgsParser(): @Boris see TODO in AliPainterTest() ArgsPareserTest()
  for (Int_t i = 0; i < fitOptions.size() && fitOptions[0] != ""; i++) {
    fitStr += fitOptions[i];
    if (i != fitOptions.size() - 1) fitStr += ",";
  }
  std::map<TString,TString>::iterator it;
  for (it=optionValues.begin(); it!=optionValues.end(); ++it) {
    if (it->second != TString(""))  drawString = it->first + TString("=") + it->second + TString(",");
    if (it->second != TString("") && it->first == TString("dOption")) drawString = it->first + TString("=") + it->second;
  }
  for (Int_t j = 0; j < histCnt; j++) {
    if (!rangesString.empty()) {
      rangeString = rangesString[j].Data();
      fromStart = 0;
      rangeVec.clear();
      while (rangesString[j].Tokenize(range, fromStart, ",")) rangeVec.push_back(range);
      for (Int_t i = 0; i < rangeVec.size(); i += 2) {
        if (rangeVec[i].CountChar('.') > 0 || rangeVec[i+1].CountChar('.') > 0) {
          hisN->GetAxis(i / 2)->SetRangeUser(rangeVec[i].Atof(), rangeVec[i + 1].Atof());
          if (verbose == 4) ::Info("AliPainter::DrawHistogram", "his->GetAxis(%d)->SetRangeUser(%s,%s);", i / 2, rangeVec[i].Data(), rangeVec[i + 1].Data());
        }
        else  {
          hisN->GetAxis(i / 2)->SetRange(rangeVec[i].Atoi(), rangeVec[i + 1].Atoi());
          if (verbose == 4) ::Info("AliPainter::DrawHistogram", "his->GetAxis(%d)->SetRange(%s,%s);", i / 2, rangeVec[i].Data(), rangeVec[i + 1].Data());
        }
      }
    }
    //fixme: such names don't work with AliDrawStyle::ApplyCssStyle() @Boris
    //uniqName = TString::Format("%s(%s)(%s)(%s)(%s)[%d]", hisName.Data(), rangeString.Data(), projections.Data(), fitStr.Data(), drawString.Data(), j).Data();
    uniqName = TString::Format("%s[%d]", hisName.Data(), j).Data();

    if (optionValues["class"] != TString("") && optionValues["class"].CountChar('[') == 0)
      uniqName = TString::Format("%s.class(%s)", uniqName.Data(), optionValues["class"].Data()).Data();
    else if (optionValues["class"] != TString("") && optionValues["class"].CountChar('[') > 0)
      uniqName = TString::Format("%s.class(%s)", uniqName.Data(), TString(optionValues["class"](1,optionValues["class"].Length()-2)).Data()).Data();
    // TH1
    if (projections.CountChar(',') + 1 == 1) {
      hisArray[j] = hisN->Projection(projections.Atoi());
      hisArray[j]->SetName(uniqName);
      hisArray[j]->SetTitle(uniqName);
      if (verbose == 4) ::Info("AliPainter::DrawHistogram", "hisArray[%d]->SetName(%s)", j, uniqName.Data());
      if (keepArray != NULL) {
        keepArray->AddLast((TObject *) hisArray[j]);
        if (verbose == 4) ::Info("AliPainter::DrawHistogram", "keepArray->AddLast((TObject*) %s)", uniqName.Data());
      }
        if (TString(optionValues["ylim"](1, optionValues["ylim"].Index(",") - 1)) != TString("")) {
          key = TString(optionValues["ylim"](1, optionValues["ylim"].Index(",") - 1));
          auto kit = AliPainter::statValues.find(key);
          if (kit == AliPainter::statValues.end()) yMin =  key.Atof();
          else yMin = AliPainter::statValues[key];
          hisArray[j]->SetMinimum(yMin);
          if (verbose == 4) ::Info("AliPainter::DrawHistogram", "hisArray[%d]->SetMinimum(%f)", j, yMin);
        }
        if (TString(optionValues["ylim"](optionValues["ylim"].Index(",") + 1,
                                         optionValues["ylim"].Index("]") -
                                         optionValues["ylim"].Index(",") - 1)) != TString("")) {
          key = TString(optionValues["ylim"](optionValues["ylim"].Index(",") + 1,
                                                     optionValues["ylim"].Index("]") -
                                                     optionValues["ylim"].Index(",") - 1));
          auto kit = AliPainter::statValues.find(key);
          if (kit == AliPainter::statValues.end()) yMax =  key.Atof();
          else yMax = AliPainter::statValues[key];
          hisArray[j]->SetMaximum(yMax);
          if (verbose == 4) ::Info("AliPainter::DrawHistogram", "hisArray[%d]->SetMaximum(%f)", j, yMax);
        }
      if (optionValues["div"] == TString("1")) {
        if (pad != NULL) {
          if (pad->InheritsFrom("TCanvas")) pad->cd(j + 1);
          else pad->cd();
        }
        hisArray[j]->SetStats(kFALSE);
        hisArray[j]->Draw(optionValues["dOption"].Data());
        if (pad != NULL) {
          if (!pad->InheritsFrom("TCanvas")) {
            if (cCanvas->GetListOfPrimitives()->After(pad) != NULL) {
              nextPad = (TPad *) cCanvas->GetListOfPrimitives()->After(pad);
              pad = nextPad;
            }
          }
        }
        if (verbose == 4) ::Info("AliPainter::DrawHistogram", "hisArray[%d]->Draw(%s)", j, optionValues["dOption"].Data());
        if (pad == NULL) {
          if (cCanvas->GetListOfPrimitives()->After(gPad) != NULL) {
            nextPad = (TPad *) cCanvas->GetListOfPrimitives()->After(gPad);
            nextPad->cd();
          }
        }
      }
      else {
        //fixme: here I should use more general way for check of empty pad @Boris
        if (gPad->GetListOfPrimitives()->GetEntries() > 1 || !gPad->GetListOfPrimitives()->At(0)->InheritsFrom("TFrame")) {
          hisArray[j]->SetStats(kFALSE);
          hisArray[j]->SetBit(TH1::kNoTitle);
          hisArray[j]->Draw(TString::Format("SAME%s", optionValues["dOption"].Data()).Data());
          if (verbose == 4) ::Info("AliPainter::DrawHistogram", "hisArray[%d]->Draw(\"SAME%s\")", j, optionValues["dOption"].Data());
        }
        else {
          hisArray[j]->SetStats(kFALSE);
          hisArray[j]->SetBit(TH1::kNoTitle);
          hisArray[j]->Draw(optionValues["dOption"].Data());
          if (verbose == 4)
            ::Info("AliPainter::DrawHistogram", "hisArray[%d]->Draw(%s)", j, optionValues["dOption"].Data());
        }
      }
      if (fitOptions[0] != TString("")) {
        //TODO: add check for standard names and change function to gaus and add standart root fitters @Boris
        if (AliTMinuitToolkit::GetPredefinedFitter(fitOptions[0].Data()) != NULL) {
          if (verbose == 4)
            ::Info("AliPainter::DrawHistogram", "AliTMinuitToolkit::Fit(%s)", fitStr.Data());
          AliTMinuitToolkit::Fit(hisArray[j], fitOptions[0], fitOptions[1], fitOptions[2], fitOptions[3],
                                xMin, xMax);
        }
        else {
          if (verbose == 4)
            ::Info("AliPainter::DrawHistogram",
                   "AliTMinuitToolkit::Fit(%s) - Non supported fitter %s. We will try to use standard root fitter.",
                   fitStr.Data(), fitOptions[0].Data());
          if (verbose == 4) ::Info("AliPainter::DrawHistogram", " hisArray[%d]->Fit(%s)", j, fitStr.Data());
          //BUG: fitOptions[n].Data() doesn't work for ubuntu 16.04, but works fine for macOS.
          hisArray[j]->Fit(fitOptions[0], fitOptions[1], fitOptions[2], xMin, xMax);
        }
      }
    }
    else if (projections.CountChar(',') + 1 > 3) {
      if (keepArray != NULL) keepArray->AddLast((TObject *) hisN->Projection(projections.CountChar(',') + 1));
      }
    else return;
  }
}
//TODO: add ApplyLimitsValue - call GetLimitsValue and apply it. priority: high @Boris
//Double_t AliPainter::GetLimitValue(TString expression, TPad *cPad) {
//  std::map<TString, TString> numNam;
//  //TODO: use regexp @Boris
//  numNam[TString("max")]    = TString("[0]");
//  numNam[TString("min")]    = TString("[1]");
//  numNam[TString("mean")]   = TString("[2]");
//  numNam[TString("rms")]    = TString("[3]");
//  numNam[TString("median")] = TString("[4]");
//
////TODO: combine objects into arrays and use TMath for calculating of parameters @Boris
//}
