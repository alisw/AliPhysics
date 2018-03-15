/// \ingroup STAT/test
/// \brief  test of AliDrawStyleTest macro
/// Example usage
/// Some benchmark of the performance doing parsing in loop
///      stopwatch.Start(); (for i; i<n; i++) doSomething(); formatted(dosomething, stopwatch.Print());
/*!
\code
.L $AliRoot_SRC/STAT/test/AliDrawStyleTest.C+
AliDrawStyleTest();
root.exe -b -q  $AliRoot_SRC/STAT/test/AliDrawStyleTest.C+ | tee AliDrawStyleTest.log
\endcode
*/

#include "Rtypes.h"
#include "TMath.h"
#include "AliDrawStyle.h"
#include "Riostream.h"
#include <iostream>
#include "TSystem.h"
#include "TStopwatch.h"
#include "TROOT.h"
#include "TH1.h"
#include "TGraph.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TPad.h"
#include "TLegend.h"
#include "AliPainter.h"

void AliDrawStyleTest_StyleArray();
void AliDrawStyleTest_Attributes();
void AliDrawStyleTest_ConvertColor();
void AliDrawStyleTest_GetIntValues();
void AliDrawStyleTest_GetFloatValues();
//  void AliDrawStyleTest_CSSReadWrite();
void AliDrawStyleTest_GetProperty();
void AliDrawStyleTest_ApplyCssStyle();
TCanvas *MakeTestPlot(Int_t nHis);

void AliDrawStyleTest(){
  AliDrawStyleTest_StyleArray();
  AliDrawStyleTest_Attributes();
  AliDrawStyleTest_ConvertColor();
  AliDrawStyleTest_GetIntValues();
  AliDrawStyleTest_GetFloatValues();
  //  AliDrawStyleTest_CSSReadWrite();
  AliDrawStyleTest_GetProperty();
  AliDrawStyleTest_ApplyCssStyle();
}

/// Test acces to the style indexed array
void AliDrawStyleTest_StyleArray(){
  //
  // Standard ALICE marker/colors arrays
  // TODO - extend it with the test of all line and marker atributes
  Int_t result=0;
  result = AliDrawStyle::GetMarkerStyle("1;2,3;4",0);
  if (result!=1){
    ::Error("AliDrawStyleTest","AliDrawStyle::GetMarkerStyle(\"1;2,3;4\",0)==%d should be 1-FAILED\n",result);
  }else{
    ::Info("AliDrawStyleTest","AliDrawStyle::GetMarkerStyle(\"1;2,3;4\",0)- IsOK");
  }
  result = AliDrawStyle::GetMarkerStyle("1;2,3;4",1);
  if (result!=2){
    ::Error("AliDrawStyleTest","AliDrawStyle::GetMarkerStyle(\"1;2,3;4\",1)==%d should be 2-FAILED\n",result);
  }else{
    ::Info("AliDrawStyleTest","AliDrawStyle::GetMarkerStyle(\"1;2,3;4\",1)- IsOK");
  }
  //
  result = AliDrawStyle::GetMarkerStyle("1;2,3;4",0);
  if (result!=1){
    ::Error("AliDrawStyleTest","AliDrawStyle::GetMarkerStyle(\"1;2,3;4\",0)==%d should be 1-FAILED\n",result);
  }else{
    ::Info("AliDrawStyleTest","AliDrawStyle::GetMarkerStyle(\"1;2,3;4\",0)- IsOK");
  }
  //
}

void AliDrawStyleTest_Attributes(){
  TString input="{marker_style:25,21,22,23; marker_color:1,2,4,5;}";
  if ( AliDrawStyle::GetPropertyValue(input,"marker_color").Contains("1,2,4,5")){
    ::Info("AliDrawStyleTest","AliDrawStyle::GetPropertyValue(%s, \"marker_color\")- IsOK", input.Data());
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetPropertyValue(%s, \"marker_color\")- FAILED", input.Data());
  }
  if ( AliDrawStyle::GetPropertyValue(input,"marker_style").Contains("25,21,22,23")){
    ::Info("AliDrawStyleTest","AliDrawStyle::GetPropertyValue(%s, \"marker_style\")- IsOK", input.Data());
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetPropertyValue(%s, \"marker_style\")- FAILED", input.Data());
  }
}

void AliDrawStyleTest_ConvertColor() {
  TString input="25,rgb(0,0,0),#0000FF,5,rgb(255,255,0),#00FFFF";
  Bool_t status;
  if ( AliDrawStyle::ConvertColor(input,0,status) == 25){
    ::Info("AliDrawStyleTest","AliDrawStyle::ConvertColor(\"%s\", 0, \"%d\")- IsOK", input.Data(), status);
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::ConvertColor(\"%s\", 0, \"%d\")- FAILED", input.Data(), status);
  }
  if ( AliDrawStyle::ConvertColor(input,1,status) == 1){
    ::Info("AliDrawStyleTest","AliDrawStyle::ConvertColor(\"%s\", 1, \"%d\")- IsOK", input.Data(), status);
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::ConvertColor(\"%s\", 1, \"%d\")- FAILED", input.Data(), status);
  }
  if ( AliDrawStyle::ConvertColor(input,2,status) == 4){
    ::Info("AliDrawStyleTest","AliDrawStyle::ConvertColor(\"%s\", 2, \"%d\")- IsOK", input.Data(), status);
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::ConvertColor(\"%s\", 2, \"%d\")- FAILED", input.Data(), status);
  }
  if ( AliDrawStyle::ConvertColor(input,3,status) == 5){
    ::Info("AliDrawStyleTest","AliDrawStyle::ConvertColor(\"%s\", 3, \"%d\")- IsOK", input.Data(), status);
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::ConvertColor(\"%s\", 3, \"%d\")- FAILED", input.Data(), status);
  }
  if ( AliDrawStyle::ConvertColor(input,4,status) == 5){
    ::Info("AliDrawStyleTest","AliDrawStyle::ConvertColor(\"%s\", 4, \"%d\")- IsOK", input.Data(), status);
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::ConvertColor(\"%s\", 4, \"%d\")- FAILED", input.Data(), status);
  }
  if ( AliDrawStyle::ConvertColor(input,5,status) == 7){
    ::Info("AliDrawStyleTest","AliDrawStyle::ConvertColor(\"%s\", 5, \"%d\")- IsOK", input.Data(), status);
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::ConvertColor(\"%s\", 5, \"%d\")- FAILED", input.Data(), status);
  }
}

/// test GetIntValues
void AliDrawStyleTest_GetIntValues(){
  TString input="{marker_style:25,21,22,23; marker_color:1,2,4,5;}";
  Bool_t status;
  if ( AliDrawStyle::GetNamedIntegerAt(input,"marker_color",0, status) == 1){
    ::Info("AliDrawStyleTest","AliDrawStyle::GetNamedIntegerAt(%s, \"marker_color\",0, \"%d\")- IsOK", input.Data(), status);
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetNamedIntegerAt(%s, \"marker_color\",0, \"%d\")- FAILED", input.Data(), status);
  }
  if ( AliDrawStyle::GetNamedIntegerAt(input,"marker_color",1, status) == 2){
    ::Info("AliDrawStyleTest","AliDrawStyle::GetNamedIntegerAt(%s, \"marker_color\",1, \"%d\")- IsOK", input.Data(), status);
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetNamedIntegerAt(%s, \"marker_color\",1, \"%d\")- FAILED", input.Data(), status);
  }
  if ( AliDrawStyle::GetNamedIntegerAt(input,"marker_color",2, status) == 4){
    ::Info("AliDrawStyleTest","AliDrawStyle::GetNamedIntegerAt(%s, \"marker_color\",2, \"%d\")- IsOK", input.Data(), status);
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetNamedIntegerAt(%s, \"marker_color\",2, \"%d\")- FAILED", input.Data(), status);
  }
  if ( AliDrawStyle::GetNamedIntegerAt(input,"marker_color",3, status) == 5){
    ::Info("AliDrawStyleTest","AliDrawStyle::GetNamedIntegerAt(%s, \"marker_color\",3, \"%d\")- IsOK", input.Data(), status);
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetNamedIntegerAt(%s, \"marker_color\",3, \"%d\")- FAILED", input.Data(), status);
  }
  if ( AliDrawStyle::GetNamedIntegerAt(input,"marker_style",0, status) == 25){
    ::Info("AliDrawStyleTest","AliDrawStyle::GetNamedIntegerAt(%s, \"marker_style\",0, \"%d\")- IsOK", input.Data(), status);
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetNamedIntegerAt(%s, \"marker_style\",0, \"%d\")- FAILED", input.Data(), status);
  }
  if ( AliDrawStyle::GetNamedIntegerAt(input,"marker_style",1, status) == 21){
    ::Info("AliDrawStyleTest","AliDrawStyle::GetNamedIntegerAt(%s, \"marker_style\",1, \"%d\")- IsOK", input.Data(), status);
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetNamedIntegerAt(%s, \"marker_style\",1, \"%d\")- FAILED", input.Data(), status);
  }
  if ( AliDrawStyle::GetNamedIntegerAt(input,"marker_style",2, status) == 22){
    ::Info("AliDrawStyleTest","AliDrawStyle::GetNamedIntegerAt(%s, \"marker_style\",2, \"%d\")- IsOK", input.Data(), status);
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetNamedIntegerAt(%s, \"marker_style\",2, \"%d\")- FAILED", input.Data(), status);
  }
  if ( AliDrawStyle::GetNamedIntegerAt(input,"marker_style",3, status) == 23){
    ::Info("AliDrawStyleTest","AliDrawStyle::GetNamedIntegerAt(%s, \"marker_style\",3, \"%d\")- IsOK", input.Data(), status);
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetNamedIntegerAt(%s, \"marker_style\",3, \"%d\")- FAILED", input.Data(), status);
  }
}

///
void AliDrawStyleTest_GetFloatValues(){
  TString input="{marker_style:25,21,22,23; marker_color:1,2,4,5;}";
  Bool_t status;
  if ( AliDrawStyle::GetNamedFloatAt(input,"marker_color",0, status) == 1){
    ::Info("AliDrawStyleTest","AliDrawStyle::GetNamedFloatAt(%s, \"marker_color\",0, \"%d\")- IsOK", input.Data(), status);
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetNamedFloatAt(%s, \"marker_color\",0, \"%d\")- FAILED", input.Data(), status);
  }
  if ( AliDrawStyle::GetNamedFloatAt(input,"marker_color",1, status) == 2){
    ::Info("AliDrawStyleTest","AliDrawStyle::GetNamedFloatAt(%s, \"marker_color\",1, \"%d\")- IsOK", input.Data(), status);
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetNamedFloatAt(%s, \"marker_color\",1, \"%d\")- FAILED", input.Data(), status);
  }
  if ( AliDrawStyle::GetNamedFloatAt(input,"marker_color",2, status) == 4){
    ::Info("AliDrawStyleTest","AliDrawStyle::GetNamedFloatAt(%s, \"marker_color\",2, \"%d\")- IsOK", input.Data(), status);
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetNamedFloatAt(%s, \"marker_color\",2, \"%d\")- FAILED", input.Data(), status);
  }
  if ( AliDrawStyle::GetNamedFloatAt(input,"marker_color",3, status) == 5){
    ::Info("AliDrawStyleTest","AliDrawStyle::GetNamedFloatAt(%s, \"marker_color\",3, \"%d\")- IsOK", input.Data(), status);
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetNamedFloatAt(%s, \"marker_color\",3, \"%d\")- FAILED", input.Data(), status);
  }
  if ( AliDrawStyle::GetNamedFloatAt(input,"marker_style",0, status) == 25){
    ::Info("AliDrawStyleTest","AliDrawStyle::GetNamedFloatAt(%s, \"marker_style\",0, \"%d\")- IsOK", input.Data(), status);
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetNamedFloatAt(%s, \"marker_style\",0, \"%d\")- FAILED", input.Data(), status);
  }
  if ( AliDrawStyle::GetNamedFloatAt(input,"marker_style",1, status) == 21){
    ::Info("AliDrawStyleTest","AliDrawStyle::GetNamedFloatAt(%s, \"marker_style\",1, \"%d\")- IsOK", input.Data(), status);
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetNamedFloatAt(%s, \"marker_style\",1, \"%d\")- FAILED", input.Data(), status);
  }
  if ( AliDrawStyle::GetNamedFloatAt(input,"marker_style",2, status) == 22){
    ::Info("AliDrawStyleTest","AliDrawStyle::GetNamedFloatAt(%s, \"marker_style\",2, \"%d\")- IsOK", input.Data(), status);
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetNamedFloatAt(%s, \"marker_style\",2, \"%d\")- FAILED", input.Data(), status);
  }
  if ( AliDrawStyle::GetNamedFloatAt(input,"marker_style",3, status) == 23){
    ::Info("AliDrawStyleTest","AliDrawStyle::GetNamedFloatAt(%s, \"marker_style\",3, \"%d\")- IsOK", input.Data(), status);
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetNamedFloatAt(%s, \"marker_style\",3, \"%d\")- FAILED", input.Data(), status);
  }
}


/// To test  - input CSS file to be read and than written
///          - diff between the files should be 0 except of the formatting
/// TODO test - ignoring commented fields in selector and in the declaration
// void AliDrawStyleTest_CSSReadWrite(){
//
//   if (gSystem->GetFromPipe(TString("[ -f ") + TString("$AliRoot_SRC/STAT/test/alirootTestStyle.css") +  TString(" ] && echo 1 || echo 0")) == "0") {
//     std::cout << "File doesn't exist1" << std::endl;
//     ::Info("AliDrawStyleTest","AliDrawStyleTest_CSSReadWrite()- ");
//     return;
//   }
//   TObjArray *cssArray = AliDrawStyle::ReadCSSFile("$AliRoot_SRC/STAT/test/alirootTestStyle.css",0);
//   if (cssArray == NULL) {
//     std::cout << "null-pointer error" << std::endl;
//     ::Info("AliDrawStyleTest","AliDrawStyleTest_CSSReadWrite()- ");
//     return;
//   }
//   AliDrawStyle::WriteCSSFilecssArray,"$AliRoot_SRC/STAT/test/test.css");
//   TObjArray *cssArrayFromTest = AliDrawStyle::ReadCSSFile("$AliRoot_SRC/STAT/test/test.css",0);
//   //AliDrawStyle::WriteCSSFile
//   //AliDrawStyle::WriteCSSFile
//
//   for (Int_t i = 0; i < cssArray->GetEntriesFast(); i++){
//     if ((TString(cssArray->At(i)->GetName()).ReplaceAll("\n", "") == TString(cssArrayFromTest->At(i)->GetName()).ReplaceAll("\n", "")) && TString(cssArray->At(i)->GetTitle()).ReplaceAll("\n", "") == TString(cssArrayFromTest->At(i)->GetTitle()).ReplaceAll("\n", "")) ::Info("AliDrawStyleTest","AliDrawStyleTest_CSSReadWrite()- ");
//     else ::Info("AliDrawStyleTest","AliDrawStyleTest_CSSReadWrite()- ");
//
//   }
//   gSystem->GetFromPipe("rm -f $AliRoot_SRC/STAT/test/test.css");
//
// }

void  AliDrawStyleTest_GetProperty(){
  AliDrawStyle::RegisterCssStyle("alirootTestStyle.css",AliDrawStyle::ReadCSSFile("$AliRoot_SRC/STAT/test/alirootTestStyle.css",0));
  if (AliDrawStyle::GetProperty("alirootTestStyle.css","marker-size", "TGraph", "Status", "TPC.QA.dcar_posA_1") == "1,2,3,4"){
    ::Info("AliDrawStyleTest","AliDrawStyle::GetProperty(\"alirootTestStyle.css\",\"marker-size\", \"TGraph\", \"Status\", \"TPC.QA.dcar_posA_1\")- IsOK");
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetProperty(\"alirootTestStyle.css\",\"marker-size\", \"TGraph\", \"Status\", \"TPC.QA.dcar_posA_1\")- FAILED");
  }
  if (AliDrawStyle::GetProperty("alirootTestStyle.css","marker-size", "TF1", "Status", "obj4") == "17,18,19,20"){
    ::Info("AliDrawStyleTest","AliDrawStyle::GetProperty(\"alirootTestStyle.css\",\"marker-size\", \"TF1\", \"Status\", \"obj4\")- IsOK");
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetProperty(\"alirootTestStyle.css\",\"marker-size\", \"TF1\", \"Status\", \"obj4\")- FAILED");
  }
  if (AliDrawStyle::GetProperty("alirootTestStyle.css","line-color", "TGraphErrors", "Warning", "asdasobj56") == "41,42,43,44"){
    ::Info("AliDrawStyleTest","AliDrawStyle::GetProperty(\"alirootTestStyle.css\",\"line-color\", \"TGraphErrors\", \"Warning\", \"asdasobj56\")- IsOK");
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetProperty(\"alirootTestStyle.css\",\"line-color\", \"TGraphErrors\", \"Warning\", \"asdasobj56\")- FAILED");
  }
  if (AliDrawStyle::GetProperty("alirootTestStyle.css","marker-color", "SomeNotExistingClass", "SomeNotExistingStatus", "obj3") == "37,38,39,40"){
    ::Info("AliDrawStyleTest","AliDrawStyle::GetProperty(\"alirootTestStyle.css\",\"marker-color\", \"SomeNotExistingClass\", \"SomeNotExistingStatus\", \"obj3\")- IsOK");
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetProperty(\"alirootTestStyle.css\",\"marker-color\", \"SomeNotExistingClass\", \"SomeNotExistingStatus\", \"obj3\")- FAILED");
  }
  if (AliDrawStyle::GetProperty("alirootTestStyle.css","line-color", "TH1", "Warning", "obj1") == "57,58,59,60"){
    ::Info("AliDrawStyleTest","AliDrawStyle::GetProperty(\"alirootTestStyle.css\",\"line-color\", \"TH1\", \"Warning\", \"obj1\")- IsOK");
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetProperty(\"alirootTestStyle.css\",\"line-color\", \"TH1\", \"Warning\", \"obj1\")- FAILED");
  }

}


/// Generate figut and apply 2 styles test1, -> test2 ->test1
/// resulting pad has to be the same
//SetStyle(styleName, pathToCssFile)
//specify into AliDrawStyle::ApplyCssStyle(pad or canvas, styleName) active pad or canvas, styleName
//for Apply a new style you should set new style and make AliDrawStyle::ApplyCssStyle() again with new styleName.

void AliDrawStyleTest_ApplyCssStyle(){
  TCanvas *canv = MakeTestPlot(3);
  AliDrawStyle::RegisterCssStyle("test1",AliDrawStyle::ReadCSSFile("$AliRoot_SRC/STAT/test/test1.css",0));
  AliDrawStyle::ApplyCssStyle(canv, "test1");
  canv->Print("test1.xml");
  AliDrawStyle::RegisterCssStyle("test2",AliDrawStyle::ReadCSSFile("$AliRoot_SRC/STAT/test/test2.css",0));
  AliDrawStyle::ApplyCssStyle(canv, "test2");
  AliDrawStyle::ApplyCssStyle(canv, "test1");
  canv->Print("test2-1.xml");
  Int_t nDiff = gSystem->GetFromPipe("diff  test1.xml test2-1.xml  | wc -l").Atoi();
  if (nDiff - 4  == 0) {
    ::Info("AliDrawStyleTest","AliDrawStyle::ApplyStyle(\"canv\",\"test1\")- IsOK");
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::ApplyStyle(\"canv\",\"test1\")- FAILED");
  }
}
///
TCanvas *MakeTestPlot(Int_t nHis) {
  TRandom r;
  TCanvas *exampleCanvas = new TCanvas("c1", "The AliDrawStyle::ApplyCssStyle example", 200, 10, 1200, 900);
  AliPainter::DivideTPad("<vertical>[1l,1r]", "Raw", exampleCanvas);
  exampleCanvas->cd(1);
  TH1F *hisArray[nHis];
  for (Int_t i = 0; i < nHis; i++) {
    hisArray[i] = new TH1F(TString::Format("his%d.class(Raw)", i).Data(),
                           TString::Format("his%d.class(Raw)", i).Data(), 100, -5, 5);
    hisArray[i]->SetStats(0);
    hisArray[i]->SetTitle(TString::Format("his%d", i).Data());
    hisArray[i]->SetMarkerStyle(1);
    hisArray[i]->FillRandom("gaus", 100000 / (i + 2));
    if (i == 0) hisArray[i]->Draw("err");
    else hisArray[i]->Draw("SAMEerr");
  }
  gPad->BuildLegend();

  exampleCanvas->cd(2);
  const Int_t n = 100;
  Double_t x[n], y[n];
  TGraph *grArray[nHis];
  for (Int_t j = 0; j < n; j++) {
    x[j] = j * 0.6 + 5;
    y[j] = (nHis) * log(x[j]);
  }

  for (Int_t i = 0; i < nHis; i++) {
    for (Int_t j = 0; j < n; j++) {
      x[j] = j * 0.6 + 5;
      y[j] = log(x[j]) / (i + 1);
    }
    grArray[i] = new TGraph(n, x, y);
    grArray[i]->SetName(TString::Format("graph%d.class(Raw)", i).Data());
    grArray[i]->SetTitle(TString::Format("gr%d", i).Data());
    if (i == 0) {
      grArray[i]->SetMinimum(0);
      grArray[i]->Draw("alp");
    }
    else grArray[i]->Draw("lp");
  }
  gPad->BuildLegend();
  return exampleCanvas;
}
