/// \ingroup STAT/test
/// \brief  test methods of AliPainter
/// Example usage
/*
\code
.L $AliRoot_SRC/STAT/test/AliPainterTest.C+
AliPainterTest();
root.exe -b -q  $AliRoot_SRC/STAT/test/AliPainterTest.C+ | tee AliPainterTest.log
\endcode
*/
//TODO: if methods in AliPainter.h will be private we should create new class AliPainterTest inherits from Alipainter. @Boris
#include "AliPainter.h"
#include "TError.h"
#include "TCanvas.h"
#include "Riostream.h"
#include "Rtypes.h"
#include "AliDrawStyle.h"
#include "TSystem.h"
#include "TPad.h"
#include "TObjArray.h"
#include "TKey.h"
#include "TFile.h"
#include "TTree.h"
#include <cstring>
#include "TH1.h"
#include "THn.h"

void AliPainterTest_ParseRanges();
void AliPainterTest_ParseString();
void AliPainterTest_ParseOptionString();
void AliPainterTest_ParsePandasString();
void AliPainterTest_GetNextPad();
void AliPainterTest_DivideTPad();
//void AliPainterTest_SetMultiGraphTimeAxisTest();
void AliPainterTest_DrawHistogram();
void AliPainterTest_SetLimits();
void AliPainterTest_GenerateDoxyImages();


void AliPainterTest() {
  AliPainterTest_ParseRanges();
  AliPainterTest_ParseString();
  AliPainterTest_ParseOptionString();
  AliPainterTest_ParsePandasString();
  AliPainterTest_GetNextPad();
  AliPainterTest_DivideTPad();
//  AliPainterTest_SetMultiGraphTimeAxisTest();
  AliPainterTest_DrawHistogram();
  AliPainterTest_SetLimits();
  AliPainterTest_GenerateDoxyImages();
}

void AliPainterTest_ParseOptionString() {
  auto result = 0;
  TString input="gaus,W,fitFunction(1,2,3),E,10,200";
  std::vector<TString> optValuesHandle;
  std::vector<TString> optValues;
  optValuesHandle.push_back(TString("gaus"));
  optValuesHandle.push_back(TString("W"));
  optValuesHandle.push_back(TString("fitFunction(1,2,3)"));
  optValuesHandle.push_back(TString("E"));
  optValuesHandle.push_back(TString("10"));
  optValuesHandle.push_back(TString("200"));
  optValues = AliPainter::ParseOptionString(input, 6);
  for (Int_t i = 0; i < 6; i++) if (optValuesHandle[i] != optValues[i]) result++;
  if (result > 0) {
    ::Error("AliPainterTest","AliPainter::ParseOptionString(\"%s\",6)- FAILED", input.Data());
    return;
  }
  else{
    ::Info("AliPainterTest","AliPainter::ParseOptionString(\"%s\",6)- IsOK", input.Data());
  }
  input="gaus,,fitFunction(),,10.234";
  optValuesHandle.clear();
  optValues.clear();
  optValuesHandle.push_back(TString("gaus"));
  optValuesHandle.push_back(TString(""));
  optValuesHandle.push_back(TString("fitFunction()"));
  optValuesHandle.push_back(TString(""));
  optValuesHandle.push_back(TString("10.234"));
  optValuesHandle.push_back(TString(""));
  optValues = AliPainter::ParseOptionString(input, 6);
  for (Int_t i = 0; i < 6; i++) if (optValuesHandle[i] != optValues[i]) result++;
  if (result > 0) {
    ::Error("AliPainterTest","AliPainter::ParseOptionString(\"%s\",6)- FAILED", input.Data());
    return;
  }
  else{
    ::Info("AliPainterTest","AliPainter::ParseOptionString(\"%s\",6)- IsOK", input.Data());
  }
  input="div=0,class=Mass,dOption=E,xlim=[0,10000]";
  optValuesHandle.clear();
  optValues.clear();
  optValuesHandle.push_back(TString("div=0"));
  optValuesHandle.push_back(TString("class=Mass"));
  optValuesHandle.push_back(TString("dOption=E"));
  optValuesHandle.push_back(TString("xlim=[0,10000]"));
  optValuesHandle.push_back(TString(""));
  optValuesHandle.push_back(TString(""));
  optValues = AliPainter::ParseOptionString(input,6, ',', "[]");
  for (Int_t i = 0; i < 6; i++) if (optValuesHandle[i] != optValues[i]) result++;
  if (result > 0) {
    ::Error("AliPainterTest","AliPainter::ParseOptionString(\"%s\",6, \',\', \"[]\")- FAILED", input.Data());
    return;
  }
  else{
    ::Info("AliPainterTest","AliPainter::ParseOptionString(\"%s\",6, \',\', \"[]\")- IsOK", input.Data());
  }
  input="";
  optValuesHandle.clear();
  optValues.clear();
  optValues = AliPainter::ParseOptionString(input);
  if (optValues.size() != 1 && std::strncmp(optValues[0].Data(), "", 1) != 0) {
    ::Error("AliPainterTest","AliPainter::ParseOptionString(\"%s\",6, \',\', \"[]\")- FAILED", input.Data());
    return;
  }
  else{
    ::Info("AliPainterTest","AliPainter::ParseOptionString(\"%s\",6, \',\', \"[]\")- IsOK", input.Data());
  }
  input = "1;2;3;4";
  optValues.clear();
  optValues = AliPainter::ParseOptionString(input);
  if (optValues.size() != 1 && std::strncmp(optValues[0].Data(), "1;2;3;4", 7) != 0) {
    ::Error("AliPainterTest","AliPainter::ParseOptionString(\"%s\", 6)- FAILED", input.Data());
    return;
  }
  else{
    ::Info("AliPainterTest","AliPainter::ParseOptionString(\"%s\", 6)- IsOK", input.Data());
  }

  input="1,2,3,4,5";

  optValuesHandle.clear();
  optValues.clear();
  optValuesHandle.push_back(TString("1"));
  optValuesHandle.push_back(TString("2"));
  optValues = AliPainter::ParseOptionString(input,2);
  for (Int_t i = 0; i < 2; i++) if (optValuesHandle[i] != optValues[i]) result++;
  if (result > 0) {
    ::Error("AliPainterTest","AliPainter::ParseOptionString(\"%s\",6, \',\', \"[]\")- FAILED", input.Data());
    return;
  }
  else{
    ::Info("AliPainterTest","AliPainter::ParseOptionString(\"%s\",6, \',\', \"[]\")- IsOK", input.Data());
  }

}

void AliPainterTest_ParsePandasString() {
  AliPainter::RegisterDefaultOptions();
  TString input="div=0";
  AliPainter::ParsePandasString(input, AliPainter::drawValues);
  if (AliPainter::drawValues["div"] != TString("0")) {
    ::Error("AliPainterTest","AliPainter::PandasOptionParser(\"%s\")- FAILED", input.Data());
  }
  else{
    ::Info("AliPainterTest","AliPainter::PandasOptionParser(\"%s\")- IsOK", input.Data());
  }
  input="class=Mass";
  AliPainter::ParsePandasString(input, AliPainter::drawValues);
  if (AliPainter::drawValues["class"] != TString("Mass")) {
    ::Error("AliPainterTest","AliPainter::PandasOptionParser(\"%s\")- FAILED", input.Data());
  }
  else{
    ::Info("AliPainterTest","AliPainter::PandasOptionParser(\"%s\")- IsOK", input.Data());
  }
  input="drawOpt=E";
  AliPainter::ParsePandasString(input, AliPainter::drawValues);
  if (AliPainter::drawValues["drawOpt"] != TString("E")) {
    ::Error("AliPainterTest","AliPainter::PandasOptionParser(\"%s\")- FAILED", input.Data());
  }
  else{
    ::Info("AliPainterTest","AliPainter::PandasOptionParser(\"%s\")- IsOK", input.Data());
  }
  input="xlim=[10.123,20.435]";
  AliPainter::ParsePandasString(input, AliPainter::drawValues);
  if (AliPainter::drawValues["xlim"] != TString("[10.123,20.435]")) {
    ::Error("AliPainterTest","AliPainter::PandasOptionParser(\"%s\")- FAILED", input.Data());
  }
  else{
    ::Info("AliPainterTest","AliPainter::PandasOptionParser(\"%s\")- IsOK", input.Data());
  }
  input="zlim=[]";
  AliPainter::ParsePandasString(input, AliPainter::drawValues);
  if (AliPainter::drawValues["zlim"] != TString("[]")) {
    ::Error("AliPainterTest","AliPainter::PandasOptionParser(\"%s\")- FAILED", input.Data());
  }
  else{
    ::Info("AliPainterTest","AliPainter::PandasOptionParser(\"%s\")- IsOK", input.Data());
  }

  input="ylim=";
  AliPainter::ParsePandasString(input, AliPainter::drawValues);
  if (AliPainter::drawValues["ylim"] != TString("")) {
    ::Error("AliPainterTest","AliPainter::PandasOptionParser(\"%s\")- FAILED", input.Data());
  }
  else{
    ::Info("AliPainterTest","AliPainter::PandasOptionParser(\"%s\")- IsOK", input.Data());
  }

  input="zlim=[], xlim=[10.123,20.435], ylim=, class=Mass";
  AliPainter::ParsePandasString(input, AliPainter::drawValues);
  if (AliPainter::drawValues["zlim"] != TString("[]") && AliPainter::drawValues["xlim"] != TString("[10.123,20.435]") && \
      AliPainter::drawValues["ylim"] != TString("") && AliPainter::drawValues["class"] != TString("Mass")) {
    ::Error("AliPainterTest","AliPainter::PandasOptionParser(\"%s\")- FAILED", input.Data());
  }
  else{
    ::Info("AliPainterTest","AliPainter::PandasOptionParser(\"%s\")- IsOK", input.Data());
  }

}

void AliPainterTest_ParseRanges() {
  TString input = "10,20";
  AliPainter::ParseRanges(input);

  if (AliPainter::rangesVec[0] != input) {
    ::Error("AliPainterTest","AliPainter::RangesParser(\"%s\")- FAILED", input.Data());
  }
  else{
    ::Info("AliPainterTest","AliPainter::RangesParser(\"%s\")- IsOK", input.Data());
  }
  input = "10.23,20,54,32.45,43,65.34";
  AliPainter::ParseRanges(input);
  if (AliPainter::rangesVec[0] != input) {
    ::Error("AliPainterTest","AliPainter::RangesParser(\"%s\")- FAILED", input.Data());
  }
  else{
    ::Info("AliPainterTest","AliPainter::RangesParser(\"%s\")- IsOK", input.Data());
  }
  input = "0:20:10:10";
   AliPainter::ParseRanges(input);
  std::vector<TString> outPutHandle;
  outPutHandle.push_back("0,10");
  outPutHandle.push_back("10,20");
  auto result = 0;
  if (outPutHandle.size() != AliPainter::rangesVec.size()) {
    ::Error("AliPainterTest","AliPainter::RangesParser(\"%s\")- FAILED", input.Data());
    return;
  }
  for (Int_t i = 0; i < (Int_t) outPutHandle.size(); i++) if (outPutHandle[i] != AliPainter::rangesVec[i]) result++;
  if (result > 0) {
    ::Error("AliPainterTest","AliPainter::RangesParser(\"%s\")- FAILED", input.Data());
    return;
  }
  else{
    ::Info("AliPainterTest","AliPainter::RangesParser(\"%s\")- IsOK", input.Data());
  }
  input = "10:20:10,0:40:10:20,30,40";
  AliPainter::ParseRanges(input);
  outPutHandle.clear();
  outPutHandle.push_back("10,10,0,20,30,40");
  outPutHandle.push_back("10,10,10,30,30,40");
  outPutHandle.push_back("10,10,20,40,30,40");
  outPutHandle.push_back("20,20,0,20,30,40");
  outPutHandle.push_back("20,20,10,30,30,40");
  outPutHandle.push_back("20,20,20,40,30,40");
  result = 0;
  if (outPutHandle.size() != AliPainter::rangesVec.size()) {
    ::Error("AliPainterTest","AliPainter::RangesParser(\"%s\")- FAILED", input.Data());
    return;
  }
  for (Int_t i = 0; i < (Int_t) outPutHandle.size(); i++) if (outPutHandle[i] != AliPainter::rangesVec[i]) result++;
  if (result > 0) {
    ::Error("AliPainterTest","AliPainter::RangesParser(\"%s\")- FAILED", input.Data());
    return;
  }
  else{
    ::Info("AliPainterTest","AliPainter::RangesParser(\"%s\")- IsOK", input.Data());
  }
}

void AliPainterTest_ParseString() {
  TString input="hisK0DMassQPtTgl(20,80,0:80:20:20,0,10)(0)(name=gaus,drawOpt=W)(class=Mass,drawOpt=E)";
  std::vector<TString> args = AliPainter::ParseString(input);
  if (args[0] == TString("hisK0DMassQPtTgl") && args[2] == TString("0") && args[1] == TString("20,80,0:80:20:20,0,10") && \
      args[3] == TString("name=gaus,drawOpt=W") && args[4] == TString("class=Mass,drawOpt=E")) ::Info("AliPainterTest","AliPainter::ParseString(\"%s\")- IsOK", input.Data());
  else {
    ::Error("AliPainterTest","AliPainter::ParseString(\"%s\")- FAILED", input.Data());
    return;
  }

  input="hisK0DMassQPtTgl(10,20)()()()";
  args.clear();
  args = AliPainter::ParseString(input);
  if (args[0] == TString("hisK0DMassQPtTgl") && args[1] == TString("10,20") && args[2] == TString() && \
      args[3] == TString() && args[4] == TString()) ::Info("AliPainterTest","AliPainter::ParseString(\"%s\")- IsOK", input.Data());
  else {
    ::Error("AliPainterTest","AliPainter::ParseString(\"%s\")- FAILED", input.Data());
    return;
  }
  input = "<max>+3*<rms>";
  args.clear();
  if (args[0] == TString("max") && args[1] == TString("rms")) ::Info("AliPainterTest","AliPainter::ParseString(\"%s\")- IsOK", input.Data());
  else {
    ::Error("AliPainterTest","AliPainter::ParseString(\"%s\")- FAILED", input.Data());
    return;
  }
}

void AliPainterTest_DivideTPad() {
  TCanvas *canvasQA = new TCanvas("canvasQATest", "canvasQATest", 1200, 800);
  AliPainter::DivideTPad(canvasQA, "<horizontal>[1b,1m,1r300px,1lst0.3]");
  canvasQA->Print("canvasQADivideTPadTest.xml");
  canvasQA->Print("canvasQADivideTPadTestFixed.xml");

  auto nDiff = gSystem->GetFromPipe("diff canvasQADivideTPadTest.xml $AliRoot_SRC/STAT/test/canvasQADivideTPadTestFixed.xml  | wc -l").Atoi();
  if (nDiff - 6 <= 0) {
    ::Info("AliPainterTest","AliPainter::DivideTPad(\"canvasQATest\",\"<horizontal>[1b,1m,1m,1lst0.3]\",\"test\")- IsOK");
  }else{
    ::Error("AliPainterTest","AliDrawStyle::DivideTPad(\"canvasQATest\",\"<horizontal>[1b,1m,1m,1lst0.3]\",\"test\")- FAILED");
  }
}

void AliPainterTest_DrawHistogram() {
  TFile::SetCacheFileDir(".");
  TFile *finput = TFile::Open("http://aliqatrkeos.web.cern.ch/aliqatrkeos/performance/AliPainterTest.root","CACHEREAD");
  TTree *tree = (TTree *) finput->Get("hisPtAll");
  TObjArray *hisArray = new TObjArray();
  TList *keys = finput->GetListOfKeys();
  for (Int_t iKey = 0; iKey<keys->GetEntries();iKey++) {
    TObject *o = finput->Get(TString::Format("%s;%d", keys->At(iKey)->GetName(), ((TKey *) keys->At(iKey))->GetCycle()).Data());
    hisArray->AddLast(o);
  }
  TCanvas *canvasQA = new TCanvas("canvasQA", "canvasQA", 1200, 800);
  AliPainter::DivideTPad(canvasQA, "<horizontal>[1b,1t,1,1]", "Canvas41");
  canvasQA->cd(1);
  AliPainter::DrawHistogram(hisArray, "hisPtAll(0,10)(0)()(div=1,dOption=E,class=PtAll)");
  AliPainter::DrawHistogram(hisArray, "hisPtITS(0,10)(0)()(div=1,dOption=E,class=PtIts)");
  AliPainter::DrawHistogram(hisArray, "hisK0DMassQPtTgl(1,1)(2)()(div=1,dOption=E,class=Tgl)");
  AliPainter::DrawHistogram(hisArray, "hisK0DMassQPtTgl(20,80,40:80:20:20,0,10)(0)(gaus,W)(class=Mass,dOption=E)");
  canvasQA->Print("canvasQADrawHistogramTest.xml");

  auto nDiff = gSystem->GetFromPipe("diff canvasQADrawHistogramTest.xml $AliRoot_SRC/STAT/test/canvasQADrawHistogramTestFixed.xml  | wc -l").Atoi();
  if (nDiff - 6 <= 0) {
    ::Info("AliPainterTest",
           "AliPainter::DrawHistogram(\"hisPtAll(0,10)(0)()(div=1,dOption=E,class=PtAll)\",hisArray)- IsOK");
  } else {
    ::Error("AliPainterTest",
            "AliPainter::DrawHistogram(\"hisPtAll(0,10)(0)()(div=1,dOption=E,class=PtAll)\",hisArray)- FAILED");
  }
}

void AliPainterTest_GetNextPad() {
  TCanvas *canvasQA1 = new TCanvas("canvasQA1", "canvasQA1", 1200, 800);
  canvasQA1->Divide(4,1);
  TPad *pad1 = (TPad *) canvasQA1->cd(1);
  TPad *tPad = AliPainter::GetNextPad(pad1);
  if (std::strncmp(tPad->GetName(), "canvasQA1_2",11) == 0) {
    ::Info("AliPainterTest",
           "AliPainter::GetNextPad(\"%s\",nullptr,4)- IsOK", pad1->GetName());
  } else {
    ::Error("AliPainterTest",
            "AliPainter::GetNextPad(\"%s\",nullptr,4)- FAILED", pad1->GetName());
  }
  tPad = AliPainter::GetNextPad(tPad);
  if (std::strncmp(tPad->GetName(), "canvasQA1_3",11) == 0) {
    ::Info("AliPainterTest",
           "AliPainter::GetNextPad(\"%s\",nullptr,4)- IsOK", tPad->GetName());
  } else {
    ::Error("AliPainterTest",
            "AliPainter::GetNextPad(\"%s\",nullptr,4)- FAILED", tPad->GetName());
  }
  tPad = AliPainter::GetNextPad(tPad);
  if (std::strncmp(tPad->GetName(), "canvasQA1_4",11) == 0) {
    ::Info("AliPainterTest",
           "AliPainter::GetNextPad(\"%s\",nullptr,4)- IsOK", tPad->GetName());
  } else {
    ::Error("AliPainterTest",
            "AliPainter::GetNextPad(\"%s\",nullptr,4)- FAILED", tPad->GetName());
  }

  TPad *pad3 = (TPad *) canvasQA1->cd(2);
  pad3->Divide(2,1);
  tPad = (TPad *) pad3->cd(1);
  tPad = AliPainter::GetNextPad(tPad);
  if (std::strncmp(tPad->GetName(), "canvasQA1_2_2",11) == 0) {
    ::Info("AliPainterTest",
           "AliPainter::GetNextPad(\"%s\",nullptr,4)- IsOK", tPad->GetName());
  } else {
    ::Error("AliPainterTest",
            "AliPainter::GetNextPad(\"%s\",nullptr,4)- FAILED", tPad->GetName());
  }
  tPad = AliPainter::GetNextPad(tPad);
  if (std::strncmp(tPad->GetName(), "canvasQA1_3",11) == 0) {
    ::Info("AliPainterTest",
           "AliPainter::GetNextPad(\"%s\",nullptr,4)- IsOK", tPad->GetName());
  } else {
    ::Error("AliPainterTest",
            "AliPainter::GetNextPad(\"%s\",nullptr,4)- FAILED", tPad->GetName());
  }
}

void AliPainterTest_SetLimits() {
  TFile::SetCacheFileDir(".");
  TFile *finput = TFile::Open("http://aliqatrkeos.web.cern.ch/aliqatrkeos/performance/AliPainterTest.root","CACHEREAD");
  TTree *tree = (TTree *) finput->Get("hisPtAll");
  TObjArray *hisArray = new TObjArray();
  TList *keys = finput->GetListOfKeys();
  for (Int_t iKey = 0; iKey<keys->GetEntries();iKey++) {
    TObject *o = finput->Get(TString::Format("%s;%d", keys->At(iKey)->GetName(), ((TKey *) keys->At(iKey))->GetCycle()).Data());
    hisArray->AddLast(o);
  }
  THn *hisN = (THn *) hisArray->FindObject("hisK0DMassQPtTgl");
  TObjArray *keep = NULL;
  TObjArray *hisArr = AliPainter::PrepareHistogram(hisN, "hisK0DMassQPtTgl(20,80,40:80:20:20,0,10)(0)(name=gaus,option=W)(class=Mass,drawOpt=E,ylim=[<mean>,<max>])", keep);
  Long64_t comSize = 0;
  Double_t *combineArray = AliPainter::GetDataArray(hisArr,comSize);
  Double_t mean = TMath::Mean(comSize, combineArray);
  Double_t max = TMath::MaxElement(comSize, combineArray);
  AliPainter::SetLimits(hisArr, 4);

  if ((Int_t) mean == 2028 && (Int_t) max == 10730) {
    ::Info("AliPainterTest",
           "AliPainter::SetLimits()- IsOK");
  } else {
    ::Error("AliPainterTest",
            "AliPainter::SetLimits()- FAILED");
  }

}

void AliPainterTest_GenerateDoxyImages() {
  TFile::SetCacheFileDir(".");
  TFile *finput = TFile::Open("http://aliqatrkeos.web.cern.ch/aliqatrkeos/performance/AliPainterTest.root","CACHEREAD");
  TTree *tree = (TTree *) finput->Get("hisPtAll");
  TObjArray *hisArray = new TObjArray();
  TList *keys = finput->GetListOfKeys();
  for (Int_t iKey = 0; iKey<keys->GetEntries();iKey++) {
    TObject *o = finput->Get(TString::Format("%s;%d", keys->At(iKey)->GetName(), ((TKey *) keys->At(iKey))->GetCycle()).Data());
    hisArray->AddLast(o);
  }
  TCanvas *canvasQA = new TCanvas("canvasQA", "canvasQA", 1200, 800);

  AliPainter::DivideTPad(canvasQA, "<horizontal>[1,1,1,1]", "Canvas41");
  canvasQA->cd(1);
  AliPainter::DrawHistogram(hisArray, "hisPtAll(0,10)(0)()(div=1,drawOpt=E,class=PtAll)");
  AliPainter::DrawHistogram(hisArray, "hisPtITS(0,10)(0)()(div=1,drawOpt=E,class=PtIts)");
  AliPainter::DrawHistogram(hisArray, "hisK0DMassQPtTgl(1,1)(2)()(div=1,drawOpt=E,class=Tgl)");
  AliPainter::DrawHistogram(hisArray, "hisK0DMassQPtTgl(20,80,40:80:20:20,0,10)(0)(name=gaus,option=W)(class=Mass,drawOpt=E)");
  AliDrawStyle::ApplyCssStyle(canvasQA, "figTemplateHex");
  canvasQA->SaveAs("$AliRoot_SRC/STAT/imgdoc/AliPainter_cxx_example1.png");
  canvasQA->Clear();

  AliPainter::DivideTPad(canvasQA, "<horizontal>[1,1,1,1]", "Canvas41", "");
  canvasQA->cd(1);
  AliPainter::DrawHistogram(hisArray, "hisPtAll(0,10)(0)()(div=1,dOption=E,class=PtAll)");
  AliPainter::DrawHistogram(hisArray, "hisPtITS(0,10)(0)()(div=1,dOption=E,class=PtIts)");
  AliPainter::DrawHistogram(hisArray, "hisK0DMassQPtTgl(1,1)(2)()(div=1,dOption=E,class=Tgl)");
  AliPainter::DrawHistogram(hisArray, "hisK0DMassQPtTgl(20,80,40:80:20:20,0,10)(0)(gaus,W)(class=Mass,dOption=E)");
  AliDrawStyle::RegisterCssStyle("figTemplateHex", AliDrawStyle::ReadCSSFile("$AliRoot_SRC/STAT/test/figTemplateHex.css"));
  AliDrawStyle::ApplyCssStyle(canvasQA, "figTemplateHex");
  canvasQA->SaveAs("$AliRoot_SRC/STAT/imgdoc/AliPainter_cxx_example2.png");
  canvasQA->Clear();

  AliPainter::DivideTPad(canvasQA, "<vertical>[1r,1l,1r,1l]", "Raw,Error", "");
  canvasQA->cd(1);
  AliPainter::DrawHistogram(hisArray, "hisPtAll(0,10)(0)()(div=1,dOption=E,class=PtAll)");
  AliPainter::DrawHistogram(hisArray, "hisPtITS(0,10)(0)()(div=1,dOption=E,class=PtIts)");
  AliPainter::DrawHistogram(hisArray, "hisK0DMassQPtTgl(1,1)(2)()(div=1,dOption=E,class=Tgl)");
  AliPainter::DrawHistogram(hisArray, "hisK0DMassQPtTgl(20,80,40:80:20:20,0,10)(0)(gaus,W)(class=Mass,dOption=E)");
  AliDrawStyle::RegisterCssStyle("figTemplateHex", AliDrawStyle::ReadCSSFile("$AliRoot_SRC/STAT/test/figTemplateHex.css"));
  AliDrawStyle::ApplyCssStyle(canvasQA, "figTemplateHex");
  canvasQA->SaveAs("$AliRoot_SRC/STAT/imgdoc/AliPainter_cxx_example3.png");
  canvasQA->Clear();

  AliPainter::DivideTPad(canvasQA, "<vertical>[1,3m,2m,1]", "Raw,Error", "");
  canvasQA->cd(1);
  AliPainter::DrawHistogram(hisArray, "hisPtAll(0,10)(0)()(div=1,dOption=E,class=PtAll)");
  AliPainter::DrawHistogram(hisArray, "hisPtITS(0,10)(0)()(div=1,dOption=E,class=PtIts)");
  AliPainter::DrawHistogram(hisArray, "hisK0DMassQPtTgl(1,1)(2)()(div=1,dOption=E,class=Tgl)");
  AliPainter::DrawHistogram(hisArray, "hisPtAll(0,10)(0)()(div=1,dOption=E,class=PtAll)");
  AliPainter::DrawHistogram(hisArray, "hisK0DMassQPtTgl(20,80,40:80:20:20,0,10)(0)(gaus,W)(div=1,class=Mass,dOption=E)");
  AliPainter::DrawHistogram(hisArray, "hisPtAll(0,10)(0)()(div=1,dOption=E,class=PtAll)");
  AliDrawStyle::RegisterCssStyle("figTemplateHex", AliDrawStyle::ReadCSSFile("$AliRoot_SRC/STAT/test/figTemplateHex.css"));
  AliDrawStyle::ApplyCssStyle(canvasQA, "figTemplateHex");
  canvasQA->SaveAs("$AliRoot_SRC/STAT/imgdoc/AliPainter_cxx_example4.png");
  canvasQA->Clear();

  AliPainter::DivideTPad(canvasQA, "<horizontal>[1rpx200,1,1,1]", "Raw,Error", "");
  canvasQA->cd(1);
  AliPainter::DrawHistogram(hisArray, "hisPtAll(0,10)(0)()(div=1,dOption=E,class=PtAll)");
  AliPainter::DrawHistogram(hisArray, "hisPtITS(0,10)(0)()(div=1,dOption=E,class=PtIts)");
  AliPainter::DrawHistogram(hisArray, "hisK0DMassQPtTgl(1,1)(2)()(div=1,dOption=E,class=Tgl)");
  AliPainter::DrawHistogram(hisArray, "hisK0DMassQPtTgl(20,80,40:80:20:20,0,10)(0)(gaus,W)(class=Mass,dOption=E)");
  AliDrawStyle::RegisterCssStyle("figTemplateHex", AliDrawStyle::ReadCSSFile("$AliRoot_SRC/STAT/test/figTemplateHex.css"));
  AliDrawStyle::ApplyCssStyle(canvasQA, "figTemplateHex");
  canvasQA->SaveAs("$AliRoot_SRC/STAT/imgdoc/AliPainter_cxx_example5.png");
  canvasQA->Clear();

  AliPainter::DrawHistogram(hisArray, "hisK0DMassQPtTgl()(0)()()");
  canvasQA->SaveAs("$AliRoot_SRC/STAT/imgdoc/AliPainter_cxx_example6.png");
  canvasQA->Clear();

  AliPainter::DrawHistogram(hisArray, "hisK0DMassQPtTgl(20,80)(0)()()");
  canvasQA->SaveAs("$AliRoot_SRC/STAT/imgdoc/AliPainter_cxx_example7.png");
  canvasQA->Clear();

  AliPainter::DivideTPad(canvasQA, "<horizontal>[1]", "Canvas1");
  canvasQA->cd(1);
  AliPainter::DrawHistogram(hisArray, "hisK0DMassQPtTgl(20,80,40:80:20:20,0,10)(0)(name=gaus,option=W)(class=Mass,drawOpt=E,ylim=[0,])");
  AliDrawStyle::RegisterCssStyle("figTemplateHex", AliDrawStyle::ReadCSSFile("$AliRoot_SRC/STAT/test/figTemplateHex.css"));
  AliDrawStyle::ApplyCssStyle(canvasQA, "figTemplateHex");
  canvasQA->SaveAs("$AliRoot_SRC/STAT/imgdoc/AliPainter_cxx_example8.png");
  canvasQA->Clear();

  AliPainter::DivideTPad(canvasQA, "<horizontal>[1,1]", "Canvas1");
  canvasQA->cd(1);
  AliPainter::DrawHistogram(hisArray, "hisK0DMassQPtTgl(20,80,40:80:20:20,0,10)(0)(name=gaus,option=W)(class=Mass,drawOpt=E,ylim=[0,], div=1)");
  AliDrawStyle::ApplyCssStyle(canvasQA, "figTemplateHex");
  canvasQA->SaveAs("$AliRoot_SRC/STAT/imgdoc/AliPainter_cxx_example9.png");
  canvasQA->Clear();

  AliPainter::DivideTPad(canvasQA, "<horizontal>[1,1]", "Canvas1");
  canvasQA->cd(1);
  AliPainter::DrawHistogram(hisArray, "hisK0DMassQPtTgl(20,80,40:80:20:20,0,10)(0)(name=gaus,option=W)(class=Mass,drawOpt=E,ylim=[<min>+0.5*<mean>,<max>-0.5*<mean>], div=1)");
  AliDrawStyle::ApplyCssStyle(canvasQA, "figTemplateHex");
  canvasQA->SaveAs("$AliRoot_SRC/STAT/imgdoc/AliPainter_cxx_example10.png");
  canvasQA->Clear();
}