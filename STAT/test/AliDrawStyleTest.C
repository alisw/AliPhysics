/// \ingroup STAT/test
/// \brief  test of AliDrawStyleTest macro
/// Example usage
/*!
\code
.L $AliRoot_SRC/STAT/test/AliDrawStyleTest.C
AliDrawStyleTestGetArray();
aliroot -b -q  $AliRoot_SRC/STAT/test/AliDrawStyleTest.C | tee AliDrawStyleTest.log

\endcode


*/
#include "AliDrawStyle.h"
#include "Riostream.h"
#include <iostream>

void AliDrawStyleTestStyleArray();

void AliDrawStyleTest(){
  AliDrawStyleTestStyleArray();
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

