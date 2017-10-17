/// \ingroup STAT/test
/// \brief  test of AliDrawStyleTest macro
/// Example usage
/// Some benchmark of the performance doing parsing in loop
///      stopwatch.Start(); (for i; i<n; i++) doSomething(); formatted(dosomething, stopwatch.Print());
/*!
\code
.L $AliRoot_SRC/STAT/test/AliDrawStyleTest.C+
AliDrawStyleTestGetArray();
aliroot -b -q  $AliRoot_SRC/STAT/test/AliDrawStyleTest.C+ | tee AliDrawStyleTest.log
\endcode
*/

#include "Rtypes.h"
#include "TMath.h"
#include "AliDrawStyle.h"
#include "Riostream.h"
#include <iostream>
#include "TSystem.h"

void AliDrawStyleTest_StyleArray();
void AliDrawStyleTest_Attributes();
void AliDrawStyleTest_GetIntValues();
void AliDrawStyleTest_GetFloatValues();
void AliDrawStyleTest_CSSReadWrite();

void AliDrawStyleTest(){
  AliDrawStyleTest_StyleArray();
  AliDrawStyleTest_Attributes();
  AliDrawStyleTest_GetIntValues();
  AliDrawStyleTest_GetFloatValues();
  AliDrawStyleTest_CSSReadWrite();
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

/// test AliDrawStyle::GetAttributeValue
void AliDrawStyleTest_Attributes(){
  TString input="{\nmarker_style:25,21,22,23; \nmarker_color:1,2,4,5; \n}";
  if ( AliDrawStyle::GetAttributeValue(input,"marker_color").Contains("1,2,4,5")){
    ::Info("AliDrawStyleTest","AliDrawStyle::GetAttributeValue(input,\"marker_color\")- IsOK");
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetAttributeValue(input,\"marker_color\")- FAILED");
  }
  if ( AliDrawStyle::GetAttributeValue(input,"marker_style").Contains("25,21,22,23")){
    ::Info("AliDrawStyleTest","AliDrawStyle::GetAttributeValue(input,\"marker_style\")- IsOK");
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetAttributeValue(input,\"marker_style\")- FAILED");
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
void AliDrawStyleTest_CSSReadWrite(){
  TObjArray *cssArray = AliDrawStyle::ReadCSSFile("$AliRoot_SRC/STAT/test/alirootTestStyle.css");
  AliDrawStyle::WriteCSSFile(cssArray,"test.css");
  TString diff = gSystem->GetFromPipe("diff -w -B    test.css  $AliRoot_SRC/STAT/test/alirootTestStyle.css");
  if (diff.Length()>0){
    ::Error("AliDrawStyleTest","AliDrawStyleTestStyle_CSSReadWrite- FAILED");
  }else{
    ::Info("AliDrawStyleTest","AliDrawStyleTestStyle_CSSReadWrite- IsOK");
  }
}
/// TODO - add test for IsSelected
/// Add benchnchamrk of is selected (timer, loop, print)
/// \param selector
/// \param className
/// \param attributeName
/// \return
/// Comments:
///      * We will use root TRegexp for the pattern matching
///      * tab
Bool_t  AliDrawStyleTest_IsSelected(TString selector, TString className, TString attributeName){
  TString selector= ".TH*#*xxx .TH1*#yyy  ";
  // AliDrawStyle::IsSelected(selector, "TGraph","xxx");  //should be false
  // AliDrawStyle::IsSelected(selector, "TH2","xxx");     //should be true
  // AliDrawStyle::IsSelected(selector, "TH2","yyy");     //should be false
  // AliDrawStyle::IsSelected(selector, "TH1","xxx");     //should be true
  // AliDrawStyle::IsSelected(selector, "TH1","yyy");     //should be true
  // AliDrawStyle::IsSelected(selector, "TH1","yyy");     //should be true
  // TString className="TH2";
  // TString objectName="xxx";
  // Bool_t isSelected=0;
  // TObjArray * subArray = selector.Tokenize(" \t");
  // for (Int_t i=0; i<subArray->GetEntries(); i++){
  //   TString subExp=subArray->At(i)->GetName();
  // }
}