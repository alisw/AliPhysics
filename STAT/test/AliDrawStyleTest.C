/// \ingroup STAT/test
/// \brief  test of AliDrawStyleTest macro
/// Example usage
/*!
\code
.L $AliRoot_SRC/STAT/test/AliDrawStyleTest.C
AliDrawStyleTestGetArray();
aliroot -b -q  $AliRoot_SRC/STAT/test/AliDrawStyleTest.C+ | tee AliDrawStyleTest.log

\endcode


*/
#include "Rtypes.h"
#include "TMath.h"
#include "AliDrawStyle.h"
#include "Riostream.h"
#include <iostream>

void AliDrawStyleTestStyleArray();
void AliDrawStyleTestStyle_Attributes();
void AliDrawStyleTestStyle_GetIntValues();
void AliDrawStyleTestStyle_GetFloatValues();


void AliDrawStyleTest(){
  AliDrawStyleTestStyleArray();
  AliDrawStyleTestStyle_Attributes();
  AliDrawStyleTestStyle_GetIntValues();
  AliDrawStyleTestStyle_GetFloatValues();

}

void AliDrawStyleTestStyleArray(){
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

void AliDrawStyleTestStyle_PredefinedArray(){
  //
  // AliDrawStyle::GetMarkerStyle("figTemplate",1); //should be 20
  // AliDrawStyle::GetMarkerStyle("figTemplate",1); //should be 21
  // (AliDrawStyle::GetMarkerColor("figTemplate",0)==kBlack); // should be kBlack OK
  // (AliDrawStyle::GetMarkerColor("figTemplate",1)==kRed+1)
  // AliDrawStyle::GetMarkerColor("figTemplate",0); // should be
}

///
void AliDrawStyleTestStyle_Attributes(){
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


void AliDrawStyleTestStyle_GetIntValues(){
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



void AliDrawStyleTestStyle_GetFloatValues(){
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
    ::Info("AliDrawStyleTest","AliDrawStyle::GetNamedFloatAt(input,\"marker_style\",1)- IsOK");
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
