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

void AliDrawStyleTest_StyleArray();
void AliDrawStyleTest_Attributes();
void AliDrawStyleTest_GetIntValues();
void AliDrawStyleTest_GetFloatValues();
//void AliDrawStyleTest_CSSReadWrite();
Bool_t AliDrawStyleTest_IsSelected();

void AliDrawStyleTest(){
  AliDrawStyleTest_StyleArray();
  AliDrawStyleTest_Attributes();
  AliDrawStyleTest_GetIntValues();
  AliDrawStyleTest_GetFloatValues();
//  AliDrawStyleTest_CSSReadWrite();
  AliDrawStyleTest_IsSelected();
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
  result = AliDrawStyle::GetMarkerColor("1;2,3;4",0);
  if (result!=1){
    ::Error("AliDrawStyleTest","AliDrawStyle::GetMarkerStyle(\"1;2,3;4\",0)==%d should be 1-FAILED\n",result);
  }else{
     ::Info("AliDrawStyleTest","AliDrawStyle::GetMarkerStyle(\"1;2,3;4\",0)- IsOK");
  }
  //
}

/// test AliDrawStyle::GetPropertyValue
void AliDrawStyleTest_Attributes(){
  TString input="{\nmarker_style:25,21,22,23; \nmarker_color:1,2,4,5; \n}";
  if ( AliDrawStyle::GetPropertyValue(input,"marker_color").Contains("1,2,4,5")){
    ::Info("AliDrawStyleTest","AliDrawStyle::GetPropertyValue(input,\"marker_color\")- IsOK");
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetPropertyValue(input,\"marker_color\")- FAILED");
  }
  if ( AliDrawStyle::GetPropertyValue(input,"marker_style").Contains("25,21,22,23")){
    ::Info("AliDrawStyleTest","AliDrawStyle::GetPropertyValue(input,\"marker_style\")- IsOK");
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetPropertyValue(input,\"marker_style\")- FAILED");
  }
}

/// test GetIntValues
void AliDrawStyleTest_GetIntValues(){
  TString input="{\nmarker_style:25,21,22,23; \nmarker_color:1,2,4,5; \n}";
  if ( AliDrawStyle::GetNamedIntegerAt(input,"marker_color",0) == 1){
    ::Info("AliDrawStyleTest","AliDrawStyle::GetNamedIntegerAt(input,\"marker_color\",0)- IsOK");
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetAttributeValue(input,\"marker_color\",0)- FAILED");
  }
  if ( AliDrawStyle::GetNamedIntegerAt(input,"marker_color",1) == 2){
    ::Info("AliDrawStyleTest","AliDrawStyle::GetNamedIntegerAt(input,\"marker_color\",1)- IsOK");
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetAttributeValue(input,\"marker_color\",1)- FAILED");
  }
  if ( AliDrawStyle::GetNamedIntegerAt(input,"marker_color",2) == 4){
    ::Info("AliDrawStyleTest","AliDrawStyle::GetNamedIntegerAt(input,\"marker_color\",2)- IsOK");
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetAttributeValue(input,\"marker_color\",2)- FAILED");
  }
  if ( AliDrawStyle::GetNamedIntegerAt(input,"marker_color",3) == 5){
    ::Info("AliDrawStyleTest","AliDrawStyle::GetNamedIntegerAt(input,\"marker_color\",3)- IsOK");
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetAttributeValue(input,\"marker_color\",3)- FAILED");
  }
  if ( AliDrawStyle::GetNamedIntegerAt(input,"marker_style",0) == 25){
    ::Info("AliDrawStyleTest","AliDrawStyle::GetNamedIntegerAt(input,\"marker_style\",0)- IsOK");
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetNamedIntegerAt(input,\"marker_style\",0)- FAILED");
  }
  if ( AliDrawStyle::GetNamedIntegerAt(input,"marker_style",1) == 21){
    ::Info("AliDrawStyleTest","AliDrawStyle::GetNamedIntegerAt(input,\"marker_style\",1)- IsOK");
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetNamedIntegerAt(input,\"marker_style\",1)- FAILED");
  }
  if ( AliDrawStyle::GetNamedIntegerAt(input,"marker_style",2) == 22){
    ::Info("AliDrawStyleTest","AliDrawStyle::GetNamedIntegerAt(input,\"marker_style\",1)- IsOK");
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetNamedIntegerAt(input,\"marker_style\",1)- FAILED");
  }
  if ( AliDrawStyle::GetNamedIntegerAt(input,"marker_style",3) == 23){
    ::Info("AliDrawStyleTest","AliDrawStyle::GetNamedIntegerAt(input,\"marker_style\",1)- IsOK");
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetNamedIntegerAt(input,\"marker_style\",1)- FAILED");
  }
}

///
void AliDrawStyleTest_GetFloatValues(){
  TString input="{\nmarker_style:25,21,22,23; \nmarker_color:1,2,4,5; \n}";
  if ( AliDrawStyle::GetNamedFloatAt(input,"marker_color",0) == 1){
    ::Info("AliDrawStyleTest","AliDrawStyle::GetNamedFloatAt(input,\"marker_color\",0)- IsOK");
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetAttributeValue(input,\"marker_color\",0)- FAILED");
  }
  if ( AliDrawStyle::GetNamedFloatAt(input,"marker_color",1) == 2){
    ::Info("AliDrawStyleTest","AliDrawStyle::GetNamedFloatAt(input,\"marker_color\",1)- IsOK");
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetAttributeValue(input,\"marker_color\",1)- FAILED");
  }
  if ( AliDrawStyle::GetNamedFloatAt(input,"marker_color",2) == 4){
    ::Info("AliDrawStyleTest","AliDrawStyle::GetNamedFloatAt(input,\"marker_color\",2)- IsOK");
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetAttributeValue(input,\"marker_color\",2)- FAILED");
  }
  if ( AliDrawStyle::GetNamedFloatAt(input,"marker_color",3) == 5){
    ::Info("AliDrawStyleTest","AliDrawStyle::GetNamedFloatAt(input,\"marker_color\",3)- IsOK");
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetAttributeValue(input,\"marker_color\",3)- FAILED");
  }
  if ( AliDrawStyle::GetNamedFloatAt(input,"marker_style",0) == 25){
    ::Info("AliDrawStyleTest","AliDrawStyle::GetNamedFloatAt(input,\"marker_style\",0)- IsOK");
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetNamedFloatAt(input,\"marker_style\",0)- FAILED");
  }
  if ( AliDrawStyle::GetNamedFloatAt(input,"marker_style",1) == 21){
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetNamedFloatAt(input,\"marker_style\",1)- FAILED");
  }
  if ( AliDrawStyle::GetNamedFloatAt(input,"marker_style",2) == 22){
    ::Info("AliDrawStyleTest","AliDrawStyle::GetNamedFloatAt(input,\"marker_style\",1)- IsOK");
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetNamedFloatAt(input,\"marker_style\",1)- FAILED");
  }
  if ( AliDrawStyle::GetNamedFloatAt(input,"marker_style",3) == 23){
    ::Info("AliDrawStyleTest","AliDrawStyle::GetNamedFloatAt(input,\"marker_style\",1)- IsOK");
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetNamedFloatAt(input,\"marker_style\",1)- FAILED");
  }
}


/// To test  - input CSS file to be read and than written
///          - diff between the files should be 0 except of the formatting
/// TODO test - ignoring commented fields in selector and in the declaration
// void AliDrawStyleTest_CSSReadWrite(){
//   //TObjArray *cssArray = AliDrawStyle::ReadCSSFile("$AliRoot_SRC/STAT/test/alirootTestStyle.css");
// TObjArray *cssArray = AliDrawStyle::ReadCSSFile("/Users/bdrum/Projects/alicesw/AliRoot/STAT/test/alirootTestStyle.css");
//   AliDrawStyle::WriteCSSFile(cssArray,"test.css");
//   TString diff = gSystem->GetFromPipe("diff -w -B    test.css  $AliRoot_SRC/STAT/test/alirootTestStyle.css");
//   if (diff.Length()>0){
//     ::Error("AliDrawStyleTest","AliDrawStyleTestStyle_CSSReadWrite- FAILED");
//   }else{
//     ::Info("AliDrawStyleTest","AliDrawStyleTestStyle_CSSReadWrite- IsOK");
//   }
// }
/// TODO - add test for IsSelected @done
/// Add benchmark of is selected (timer, loop, print) @done if I get you right.
/// What time do we need?
/// Now we have
/// bench(100) = 0.001079.s || bench(1000) = 0.01162.s || bench(10000) == 0.06675.s bench(100000) = 0.6101.s
/// \param selector
/// \param className
/// \param attributeName
/// \return
Bool_t  AliDrawStyleTest_IsSelected(){//TString selector, TString className, TString attributeName){
  TString selectors = "TH1.Status#obj1, TH1.Warning#obj1, TH1.Warning#obj3 \tTGraph#obj1, TGraph.Status#TPC.QA.dcar_posA_1 \tTGraph.Warning#TPC.QA.dcar_posA_2 \tTF1.Status, .Status#obj1, #obj3";
  std::cout << selectors << std::endl;
  if (AliDrawStyle::IsSelected(selectors, "TF1", "Status", "anyObject")){
    ::Info("AliDrawStyleTest","AliDrawStyle::IsSelected(selectors, \"TF1\", \"Status\", \"anyObject\")- IsOK");
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::IsSelected(selectors, \"TF1\", \"Status\", \"anyObject\")- FAILED");
  }
  if (AliDrawStyle::IsSelected(selectors, "TGraphErrors", "AnyTag", "obj3")){
    ::Info("AliDrawStyleTest","AliDrawStyle::IsSelected(selectors, \"TGraphErrors\", \"AnyTag\", \"obj3\")- IsOK");
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::IsSelected(selectors, \"TGraphErrors\", \"AnyTag\", \"obj3\")- FAILED");
  }
  if (AliDrawStyle::IsSelected(selectors, "TGraph", "Warning", "TPC.QA.dcar_posA_2")){
    ::Info("AliDrawStyleTest","AliDrawStyle::IsSelected(selectors, \"TGraph\", \"Warning\", \"TPC.QA.dcar_posA_2\")- IsOK");
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::IsSelected(selectors, \"TGraph\", \"Warning\", \"TPC.QA.dcar_posA_2\")- FAILED");
  }
  if (!AliDrawStyle::IsSelected(selectors, "TH1", "Status", "obj2")){
    ::Info("AliDrawStyleTest","!AliDrawStyle::IsSelected(selectors, \"TH1\", \"Status\", \"obj2\")- IsOK");
  }else{
    ::Error("AliDrawStyleTest","!AliDrawStyle::IsSelected(selectors, \"TH1\", \"Status\", \"obj2\")- FAILED");
  }
  if (!AliDrawStyle::IsSelected(selectors, "Graph", "Warning", "obj1")){
    ::Info("AliDrawStyleTest","!AliDrawStyle::IsSelected(selectors, \"Graph\", \"Warning\", \"obj1\")- IsOK");
  }else{
    ::Error("AliDrawStyleTest","!AliDrawStyle::IsSelected(selectors, \"Graph\", \"Warning\", \"obj1\")- FAILED");
  }

  //should I make this like separate function?
  // Int_t n = 100000;
  // TStopwatch timer;
  // timer.Start();
  // for (Int_t i=0; i<n; i++){
  //   AliDrawStyle::IsSelected(selector, "tag123", "TGraphErrors","yyy");
  // }
  // timer.Stop();
  //
  // std::cout << "Benchmarking:" <<endl;
  // std::cout.precision(4);
  // std::cout << n << "calls of AliDrawStyle::IsSelected(tag*.TH*#*xxx \ttag*.*Errors#*yyy  , \"tag123\", \"TGraphErrors\",\"yyy\") tooks" << timer.RealTime() << ".s" <<endl;

}

Bool_t  AliDrawStyleTest_GetProperty(){
  AliDrawStyle::SetCssStyle("alirootTestStyle.css",AliDrawStyle::ReadCSSFile("$AliRoot_SRC/STAT/test/alirootTestStyle.css",0));
  if (AliDrawStyle::GetProperty("alirootTestStyle.css","marker_size", "TGraph", "Status", "TPC.QA.dcar_posA_1") == "   1,2,3,4"){
    ::Info("AliDrawStyleTest","AliDrawStyle::GetProperty(\"alirootTestStyle.css\",\"marker_size\", \"TGraph\", \"Status\", \"TPC.QA.dcar_posA_1\")- IsOK");
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetProperty(\"alirootTestStyle.css\",\"marker_size\", \"TGraph\", \"Status\", \"TPC.QA.dcar_posA_1\")- FAILED");
  }



}
