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

void AliPainterTest_OptionStringParser();
void AliPainterTest_PandasOptionParser();
void AliPainterTest_RangesParser();
void AliPainterTest_ArgsParser();
void AliPainterTest_DivideTPad();
//void AliPainterTest_SetMultiGraphTimeAxisTest();
void AliPainterTest_DrawHistogram();
//void AliPainterTest_GetLimitValueTest();
//void AliPainterTest_ApplyLimitValueTest();

void AliPainterTest(){
  AliPainterTest_OptionStringParser();
  AliPainterTest_PandasOptionParser();
  AliPainterTest_RangesParser();
  AliPainterTest_ArgsParser();
  AliPainterTest_DivideTPad();
//  AliPainterTest_SetMultiGraphTimeAxisTest();
  AliPainterTest_DrawHistogram();
//  AliPainterTest_GetLimitValueTest();
//  AliPainterTest_ApplyLimitValueTest();
}

void AliPainterTest_OptionStringParser(){
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
  optValues = AliPainter::OptionStringParser(input,"()",6);
  for (Int_t i = 0; i < 6; i++) if (optValuesHandle[i] != optValues[i]) result++;
  if (result > 0) {
    ::Error("AliPainterTest","AliPainter::OptionStringParser(\"%s\",\"()\",6)- FAILED", input.Data());
    return;
  }
  else{
    ::Info("AliPainterTest","AliPainter::OptionStringParser(\"%s\",\"()\",6)- IsOK", input.Data());
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
  optValues = AliPainter::OptionStringParser(input,"()",6);
  for (Int_t i = 0; i < 6; i++) if (optValuesHandle[i] != optValues[i]) result++;
  if (result > 0) {
    ::Error("AliPainterTest","AliPainter::OptionStringParser(\"%s\",\"()\",6)- FAILED", input.Data());
    return;
  }
  else{
    ::Info("AliPainterTest","AliPainter::OptionStringParser(\"%s\",\"()\",6)- IsOK", input.Data());
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
  optValues = AliPainter::OptionStringParser(input,"[]",6);
  for (Int_t i = 0; i < 6; i++) if (optValuesHandle[i] != optValues[i]) result++;
  if (result > 0) {
    ::Error("AliPainterTest","AliPainter::OptionStringParser(\"%s\",\"[]\",6)- FAILED", input.Data());
    return;
  }
  else{
    ::Info("AliPainterTest","AliPainter::OptionStringParser(\"%s\",\"[]\",6)- IsOK", input.Data());
  }
  input="";
  optValuesHandle.clear();
  optValues.clear();
  optValuesHandle.push_back(TString(""));
  optValuesHandle.push_back(TString(""));
  optValuesHandle.push_back(TString(""));
  optValuesHandle.push_back(TString(""));
  optValuesHandle.push_back(TString(""));
  optValuesHandle.push_back(TString(""));
  optValues = AliPainter::OptionStringParser(input,"[]",6);
  for (Int_t i = 0; i < 6; i++) if (optValuesHandle[i] != optValues[i]) result++;
  if (result > 0) {
    ::Error("AliPainterTest","AliPainter::OptionStringParser(\"%s\",\"[]\",6)- FAILED", input.Data());
    return;
  }
  else{
    ::Info("AliPainterTest","AliPainter::OptionStringParser(\"%s\",\"[]\",6)- IsOK", input.Data());
  }
  optValues.clear();
  optValues = AliPainter::OptionStringParser(input,"()",6);
  for (Int_t i = 0; i < 6; i++) if (optValuesHandle[i] != optValues[i]) result++;
  if (result > 0) {
    ::Error("AliPainterTest","AliPainter::OptionStringParser(\"%s\",\"()\",6)- FAILED", input.Data());
    return;
  }
  else{
    ::Info("AliPainterTest","AliPainter::OptionStringParser(\"%s\",\"()\",6)- IsOK", input.Data());
  }
}

void AliPainterTest_PandasOptionParser(){
  AliPainter::RegisterDefaultOptions();
  TString input="div=0";
  AliPainter::PandasOptionParser(input);
  if (AliPainter::optionValues["div"] != TString("0")) {
    ::Error("AliPainterTest","AliPainter::PandasOptionParser(\"%s\")- FAILED", input.Data());
  }
  else{
    ::Info("AliPainterTest","AliPainter::PandasOptionParser(\"%s\")- IsOK", input.Data());
  }
  input="class=Mass";
  AliPainter::PandasOptionParser(input);
  if (AliPainter::optionValues["class"] != TString("Mass")) {
    ::Error("AliPainterTest","AliPainter::PandasOptionParser(\"%s\")- FAILED", input.Data());
  }
  else{
    ::Info("AliPainterTest","AliPainter::PandasOptionParser(\"%s\")- IsOK", input.Data());
  }
  input="dOption=E";
  AliPainter::PandasOptionParser(input);
  if (AliPainter::optionValues["dOption"] != TString("E")) {
    ::Error("AliPainterTest","AliPainter::PandasOptionParser(\"%s\")- FAILED", input.Data());
  }
  else{
    ::Info("AliPainterTest","AliPainter::PandasOptionParser(\"%s\")- IsOK", input.Data());
  }
  input="xlim=[10.123,20.435]";
  AliPainter::PandasOptionParser(input);
  if (AliPainter::optionValues["xlim"] != TString("[10.123,20.435]")) {
    ::Error("AliPainterTest","AliPainter::PandasOptionParser(\"%s\")- FAILED", input.Data());
  }
  else{
    ::Info("AliPainterTest","AliPainter::PandasOptionParser(\"%s\")- IsOK", input.Data());
  }
  input="zlim=[]";
  AliPainter::PandasOptionParser(input);
  if (AliPainter::optionValues["zlim"] != TString("[]")) {
    ::Error("AliPainterTest","AliPainter::PandasOptionParser(\"%s\")- FAILED", input.Data());
  }
  else{
    ::Info("AliPainterTest","AliPainter::PandasOptionParser(\"%s\")- IsOK", input.Data());
  }

  input="ylim=";
  AliPainter::PandasOptionParser(input);
  if (AliPainter::optionValues["ylim"] != TString("")) {
    ::Error("AliPainterTest","AliPainter::PandasOptionParser(\"%s\")- FAILED", input.Data());
  }
  else{
    ::Info("AliPainterTest","AliPainter::PandasOptionParser(\"%s\")- IsOK", input.Data());
  }
}

void AliPainterTest_RangesParser() {
  TString input = "10,20";
  std::vector<TString> outPut = AliPainter::RangesParser(input);
  if (outPut[0] != input) {
    ::Error("AliPainterTest","AliPainter::RangesParser(\"%s\")- FAILED", input.Data());
  }
  else{
    ::Info("AliPainterTest","AliPainter::RangesParser(\"%s\")- IsOK", input.Data());
  }
  input = "10.23,20,54,32.45,43,65.34";
  outPut.clear();
  outPut = AliPainter::RangesParser(input);
  if (outPut[0] != input) {
    ::Error("AliPainterTest","AliPainter::RangesParser(\"%s\")- FAILED", input.Data());
  }
  else{
    ::Info("AliPainterTest","AliPainter::RangesParser(\"%s\")- IsOK", input.Data());
  }
  input = "0:20:10:10";
  outPut.clear();
  outPut = AliPainter::RangesParser(input);
  std::vector<TString> outPutHandle;
  outPutHandle.push_back("0,10");
  outPutHandle.push_back("10,20");
  auto result = 0;
  if (outPutHandle.size() != outPut.size()) {
    ::Error("AliPainterTest","AliPainter::RangesParser(\"%s\")- FAILED", input.Data());
    return;
  }
  for (Int_t i = 0; i < (Int_t) outPutHandle.size(); i++) if (outPutHandle[i] != outPut[i]) result++;
  if (result > 0) {
    ::Error("AliPainterTest","AliPainter::RangesParser(\"%s\")- FAILED", input.Data());
    return;
  }
  else{
    ::Info("AliPainterTest","AliPainter::RangesParser(\"%s\")- IsOK", input.Data());
  }
  input = "10:20:10,0:40:10:20,30,40";
  outPut.clear();
  outPut = AliPainter::RangesParser(input);
  outPutHandle.clear();
  outPutHandle.push_back("10,10,0,20,30,40");
  outPutHandle.push_back("10,10,10,30,30,40");
  outPutHandle.push_back("10,10,20,40,30,40");
  outPutHandle.push_back("20,20,0,20,30,40");
  outPutHandle.push_back("20,20,10,30,30,40");
  outPutHandle.push_back("20,20,20,40,30,40");
  result = 0;
  if (outPutHandle.size() != outPut.size()) {
    ::Error("AliPainterTest","AliPainter::RangesParser(\"%s\")- FAILED", input.Data());
    return;
  }
  for (Int_t i = 0; i < (Int_t) outPut.size(); i++) if (outPutHandle[i] != outPut[i]) result++;
  if (result > 0) {
    ::Error("AliPainterTest","AliPainter::RangesParser(\"%s\")- FAILED", input.Data());
    return;
  }
  else{
    ::Info("AliPainterTest","AliPainter::RangesParser(\"%s\")- IsOK", input.Data());
  }
}

void AliPainterTest_ArgsParser(){
  auto result=0;
  TString input="hisK0DMassQPtTgl(20,80,0:80:20:20,0,10)(0)(gaus,W)(class=Mass,dOption=E)";
  TString hisName = "";
  TString projections = "";
  std::vector<TString> fitOptions;
  std::vector<TString> rangesStrings;

  AliPainter::ArgsParser(input, hisName, projections, fitOptions, rangesStrings,4);
  if (hisName != TString("hisK0DMassQPtTgl")){
    ::Error("AliPainterTest","AliPainter::ArgsParser(\"%s\")- FAILED", input.Data());
    return;
  }
  else{
    ::Info("AliPainterTest","AliPainter::ArgsParser(\"%s\")- IsOK", input.Data());
  }
  if (projections != TString("0")){
    ::Error("AliPainterTest","AliPainter::ArgsParser(\"%s\")- FAILED", input.Data());
    return;
  }
  else{
    ::Info("AliPainterTest","AliPainter::ArgsParser(\"%s\")- IsOK", input.Data());
  }
}

void AliPainterTest_DivideTPad() {
  TCanvas *canvasQA = new TCanvas("canvasQATest", "canvasQATest", 1200, 800);
  AliPainter::DivideTPad("<horizontal>[1b,1m,1r300px,1lst0.3]", "", "", canvasQA);
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
  auto *finput = TFile::Open("$AliRoot_SRC/STAT/test/AliPainterTest.root");
  auto *tree = (TTree *) finput->Get("hisPtAll");
  auto *hisArray = new TObjArray();
  auto *keys = finput->GetListOfKeys();
  for (Int_t iKey = 0; iKey < keys->GetEntries(); iKey++) {
    auto *o = finput->Get(
            TString::Format("%s;%d", keys->At(iKey)->GetName(), ((TKey *) keys->At(iKey))->GetCycle()).Data());
    hisArray->AddLast(o);
  }
  auto *canvasQA = new TCanvas("canvasQA", "canvasQA", 1200, 800);
  AliPainter::DivideTPad("<horizontal>[1b,1t,1,1]", "Canvas41", "", canvasQA);
  canvasQA->cd(1);
  AliPainter::DrawHistogram((char *) "hisPtAll(0,10)(0)()(div=1,dOption=E,class=PtAll)", hisArray);
  AliPainter::DrawHistogram((char *) "hisPtITS(0,10)(0)()(div=1,dOption=E,class=PtIts)", hisArray);
  AliPainter::DrawHistogram((char *) "hisK0DMassQPtTgl(1,1)(2)()(div=1,dOption=E,class=Tgl)", hisArray);
  AliPainter::DrawHistogram((char *) "hisK0DMassQPtTgl(20,80,40:80:20:20,0,10)(0)(gaus,W)(class=Mass,dOption=E)",
                            hisArray);
  canvasQA->Print("canvasQADrawHistogramTest.xml");

  auto nDiff = gSystem->GetFromPipe("diff canvasQADrawHistogramTest.xml $AliRoot_SRC/STAT/test/canvasQADrawHistogramTestFixed.xml  | wc -l").Atoi();
  if (nDiff - 6 <= 0) {
    ::Info("AliPainterTest",
           "AliPainterTest::DrawHistogram(\"hisPtAll(0,10)(0)()(div=1,dOption=E,class=PtAll)\",hisArray)- IsOK");
  } else {
    ::Error("AliPainterTest",
            "AliPainterTest::DrawHistogram(\"hisPtAll(0,10)(0)()(div=1,dOption=E,class=PtAll)\",hisArray)- FAILED");
  }
}