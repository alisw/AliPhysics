/// \ingroup STAT/test
/// \brief  test of AliDrawStyleTest macro
/// Example usage
/*!
 * For running:
\code
.L $AliRoot_SRC/STAT/test/AliDrawStyleTest.C+
AliDrawStyleTest();
\endcode
*/

#include "Rtypes.h"
#include "TMath.h"
#include "TStyle.h"
#include "TPaveText.h"
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
void AliDrawStyleTest_ParseDeclaration();
void AliDrawStyleTest_ConvertColor();
void AliDrawStyleTest_GetIntValues();
void AliDrawStyleTest_GetFloatValues();
void AliDrawStyleTest_IsSelected();
void AliDrawStyleTest_GetIds();
//  void AliDrawStyleTest_CSSReadWrite();
void AliDrawStyleTest_GetValue();
void AliDrawStyleTest_ApplyCssStyle();
TCanvas *MakeTestPlot(Int_t nHis);
void AliDrawStyleTest_GenerateDoxyImages();

void AliDrawStyleTest() {
  AliDrawStyleTest_StyleArray();
  AliDrawStyleTest_ParseDeclaration();
  AliDrawStyleTest_ConvertColor();
  AliDrawStyleTest_GetIntValues();
  AliDrawStyleTest_GetFloatValues();
  AliDrawStyleTest_IsSelected();
  AliDrawStyleTest_GetIds();
//  //  AliDrawStyleTest_CSSReadWrite();
  AliDrawStyleTest_GetValue();
  AliDrawStyleTest_ApplyCssStyle();
  AliDrawStyleTest_GenerateDoxyImages();
}

void AliDrawStyleTest_StyleArray() {
  Int_t result=0;
  result = AliDrawStyle::GetMarkerStyle("1;2,3;4",0);
  if (result!=1) {
    ::Error("AliDrawStyleTest","AliDrawStyle::GetMarkerStyle(\"1;2,3;4\",0)==%d should be 1-FAILED\n",result);
  }else{
    ::Info("AliDrawStyleTest","AliDrawStyle::GetMarkerStyle(\"1;2,3;4\",0)- IsOK");
  }
  result = AliDrawStyle::GetMarkerStyle("1;2,3;4",1);
  if (result!=2) {
    ::Error("AliDrawStyleTest","AliDrawStyle::GetMarkerStyle(\"1;2,3;4\",1)==%d should be 2-FAILED\n",result);
  }else{
    ::Info("AliDrawStyleTest","AliDrawStyle::GetMarkerStyle(\"1;2,3;4\",1)- IsOK");
  }
  result = AliDrawStyle::GetMarkerStyle("1;2,3;4",0);
  if (result!=1) {
    ::Error("AliDrawStyleTest","AliDrawStyle::GetMarkerStyle(\"1;2,3;4\",0)==%d should be 1-FAILED\n",result);
  }else{
    ::Info("AliDrawStyleTest","AliDrawStyle::GetMarkerStyle(\"1;2,3;4\",0)- IsOK");
  }
}

void AliDrawStyleTest_ParseDeclaration() {
  TString input="{marker-style:25,21,22,23; marker-color:1,2,4,5;}";
  if ( AliDrawStyle::ParseDeclaration(input,"marker-color").Contains("1,2,4,5")) {
    ::Info("AliDrawStyleTest","AliDrawStyle::ParseDeclaration(%s, \"marker-color\")- IsOK", input.Data());
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::ParseDeclaration(%s, \"marker-color\")- FAILED", input.Data());
  }
  if ( AliDrawStyle::ParseDeclaration(input,"marker-style").Contains("25,21,22,23")) {
    ::Info("AliDrawStyleTest","AliDrawStyle::ParseDeclaration(%s, \"marker-style\")- IsOK", input.Data());
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::ParseDeclaration(%s, \"marker-style\")- FAILED", input.Data());
  }
}

void AliDrawStyleTest_ConvertColor() {
  TString input="25";
  Bool_t status;
  if ( AliDrawStyle::ConvertColor(input.Data()) == 25) {
    ::Info("AliDrawStyleTest","AliDrawStyle::ConvertColor(\"%s\")- IsOK", input.Data());
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::ConvertColor(\"%s\")- FAILED", input.Data());
  }
  input = "rgb(0,0,0)";
  if ( AliDrawStyle::ConvertColor(input.Data()) == 1) {
    ::Info("AliDrawStyleTest","AliDrawStyle::ConvertColor(\"%s\")- IsOK", input.Data());
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::ConvertColor(\"%s\")- FAILED", input.Data());
  }
  input = "#0000FF";
  if ( AliDrawStyle::ConvertColor(input.Data()) == 4) {
    ::Info("AliDrawStyleTest","AliDrawStyle::ConvertColor(\"%s\")- IsOK", input.Data());
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::ConvertColor(\"%s\")- FAILED", input.Data());
  }
  input = "5";
  if ( AliDrawStyle::ConvertColor(input.Data()) == 5) {
    ::Info("AliDrawStyleTest","AliDrawStyle::ConvertColor(\"%s\")- IsOK", input.Data());
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::ConvertColor(\"%s\")- FAILED", input.Data());
  }
  input = "rgb(255,255,0)";
  if ( AliDrawStyle::ConvertColor(input.Data()) == 5) {
    ::Info("AliDrawStyleTest","AliDrawStyle::ConvertColor(\"%s\")- IsOK", input.Data());
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::ConvertColor(\"%s\")- FAILED", input.Data());
  }
  input = "#00FFFF";
  if ( AliDrawStyle::ConvertColor(input.Data()) == 7) {
    ::Info("AliDrawStyleTest","AliDrawStyle::ConvertColor(\"%s\")- IsOK", input.Data());
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::ConvertColor(\"%s\")- FAILED", input.Data());
  }
}

void AliDrawStyleTest_GetIntValues() {
  TString input="{marker_style:25,21,22,23; marker-color:1,2,4,5,rgb(123,123,123),#dfdfdf;}";
  Bool_t status;
  if ( AliDrawStyle::GetNamedTypeAt(input,status, 0, "marker-color") == 1) {
    ::Info("AliDrawStyleTest","AliDrawStyle::GetNamedTypeAt(%s, %d, 0, \"marker-color\")- IsOK", input.Data(), status);
  } else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetNamedTypeAt(%s, %d, 0, \"marker-color\")- FAILED", input.Data(), status);
  }
  if ( AliDrawStyle::GetNamedTypeAt(input,status, 1, "marker-color") == 2) {
    ::Info("AliDrawStyleTest","AliDrawStyle::GetNamedTypeAt(%s, %d, 1, \"marker-color\")- IsOK", input.Data(), status);
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetNamedTypeAt(%s, %d, 1, \"marker-color\")- FAILED", input.Data(), status);
  }
  if ( AliDrawStyle::GetNamedTypeAt(input,status, 2, "marker-color") == 4) {
    ::Info("AliDrawStyleTest","AliDrawStyle::GetNamedTypeAt(%s, %d, 2, \"marker-color\")- IsOK", input.Data(), status);
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetNamedTypeAt(%s, %d, 2, \"marker-color\")- FAILED", input.Data(), status);
  }
  if ( AliDrawStyle::GetNamedTypeAt(input,status, 3, "marker-color") == 5) {
    ::Info("AliDrawStyleTest","AliDrawStyle::GetNamedTypeAt(%s, %d, 3, \"marker-color\")- IsOK", input.Data(), status);
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetNamedTypeAt(%s, %d, 3, \"marker-color\")- FAILED", input.Data(), status);
  }
  if ( AliDrawStyle::GetNamedTypeAt(input,status, 0, "marker_style") == 25) {
    ::Info("AliDrawStyleTest","AliDrawStyle::GetNamedTypeAt(%s, %d, 0, \"marker_style\")- IsOK", input.Data(), status);
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetNamedTypeAt(%s, %d, 0, \"marker_style\")- FAILED", input.Data(), status);
  }
  if ( AliDrawStyle::GetNamedTypeAt(input,status, 1, "marker_style") == 21) {
    ::Info("AliDrawStyleTest","AliDrawStyle::GetNamedTypeAt(%s, %d, 1, \"marker_style\")- IsOK", input.Data(), status);
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetNamedTypeAt(%s, %d, 1, \"marker_style\")- FAILED", input.Data(), status);
  }
  if ( AliDrawStyle::GetNamedTypeAt(input,status, 2, "marker_style") == 22) {
    ::Info("AliDrawStyleTest","AliDrawStyle::GetNamedTypeAt(%s, %d, 2, \"marker_style\")- IsOK", input.Data(), status);
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetNamedTypeAt(%s, %d, 2, \"marker_style\")- FAILED", input.Data(), status);
  }
  if ( AliDrawStyle::GetNamedTypeAt(input,status, 3, "marker_style") == 23) {
    ::Info("AliDrawStyleTest","AliDrawStyle::GetNamedTypeAt(%s, %d, 3, \"marker_style\")- IsOK", input.Data(), status);
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetNamedTypeAt(%s, %d, 3, \"marker_style\")- FAILED", input.Data(), status);
  }
  if ( AliDrawStyle::GetNamedTypeAt(input,status, 4, "marker-color") == 202 || AliDrawStyle::GetNamedTypeAt(input,status, 4, "marker-color") == 14) {
    ::Info("AliDrawStyleTest","AliDrawStyle::GetNamedTypeAt(%s, %d, 4, \"marker-color\")- IsOK", input.Data(), status);
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetNamedTypeAt(%s, %d, 4, \"marker-color\")- FAILED", input.Data(), status);
  }
  if ( AliDrawStyle::GetNamedTypeAt(input,status, 5, "marker-color") == 1179 || AliDrawStyle::GetNamedTypeAt(input,status, 5, "marker-color") == 924 || AliDrawStyle::GetNamedTypeAt(input,status, 5, "marker-color") == 18) {
    ::Info("AliDrawStyleTest","AliDrawStyle::GetNamedTypeAt(%s, %d, 5, \"marker-color\")- IsOK", input.Data(), status);
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetNamedTypeAt(%s, %d, 5, \"marker-color\")- FAILED", input.Data(), status);
  }
}

void AliDrawStyleTest_GetFloatValues() {
  TString input="{margin-top:300%,474px,5; marker-size:8px,2,4,500%;}";
  Bool_t status;

  if ( AliDrawStyle::GetNamedTypeAt(input,status,0, "marker-size") - 1 <= 0.01) {
    ::Info("AliDrawStyleTest","AliDrawStyle::GetNamedTypeAt(\"%s\", %d, 0, \"marker-size\")- IsOK", input.Data(), status);
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetNamedTypeAt(\"%s\", %d, 0, \"marker-size\")- FAILED", input.Data(), status);
  }
  if ( AliDrawStyle::GetNamedTypeAt(input,status,1, "marker-size") - 2 <= 0.01) {
    ::Info("AliDrawStyleTest","AliDrawStyle::GetNamedTypeAt(\"%s\", %d, 1, \"marker-size\")- IsOK", input.Data(), status);
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetNamedTypeAt(\"%s\", %d, 1, \"marker-size\")- FAILED", input.Data(), status);
  }
  if (AliDrawStyle::GetNamedTypeAt(input,status,2, "marker-size") - 4 <= 0.01) {
    ::Info("AliDrawStyleTest","AliDrawStyle::GetNamedTypeAt(\"%s\", %d, 2, \"marker-size\")- IsOK", input.Data(), status);
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetNamedTypeAt(\"%s\", %d, 2, \"marker-size\")- FAILED", input.Data(), status);
  }
  if ( AliDrawStyle::GetNamedTypeAt(input,status,3, "marker-size") - 5 <= 0.01) {
    ::Info("AliDrawStyleTest","AliDrawStyle::GetNamedTypeAt(\"%s\", %d, 3, \"marker-size\")- IsOK", input.Data(), status);
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetNamedTypeAt(\"%s\", %d, 3, \"marker-size\")- FAILED", input.Data(), status);
  }
  TCanvas *exampleCanvas = new TCanvas("c1", "The AliDrawStyle::ApplyCssStyle example", 200, 10, 1200, 900);
  TPad *pad = new TPad("testPad", "testPad",0,0,1,1);
  pad->Draw();
  if ( AliDrawStyle::GetNamedTypeAt(input,status,0, "margin-top") - 3 <= 0.01) {
    ::Info("AliDrawStyleTest","AliDrawStyle::GetNamedTypeAt(\"%s\", %d, 0, \"margin-top\")- IsOK", input.Data(), status);
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetNamedTypeAt(\"%s\", %d, 0, \"margin-top\")- FAILED", input.Data(), status);
  }
  if ( AliDrawStyle::GetNamedTypeAt(input,status,1, "margin-top") - 1 <= 0.01) {
    ::Info("AliDrawStyleTest","AliDrawStyle::GetNamedTypeAt(\"%s\", %d, 1, \"margin-top\")- IsOK", input.Data(), status);
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetNamedTypeAt(\"%s\", %d, 1, \"margin-top\")- FAILED", input.Data(), status);
  }
  if ( AliDrawStyle::GetNamedTypeAt(input,status,2, "margin-top") - 5 <= 0.01) {
    ::Info("AliDrawStyleTest","AliDrawStyle::GetNamedTypeAt(\"%s\", %d, 2, \"margin-top\")- IsOK", input.Data(), status);
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetNamedTypeAt(\"%s\", %d, 2, \"margin-top\")- FAILED", input.Data(), status);
  }
}

void AliDrawStyleTest_IsSelected() {
  TString selectors = "TH1.Status#obj1, TH1.Warning#obj1, TH1.Warning#obj3 \tTGraph#obj1, TGraph.Status#TPC.QA.dcar_posA_1 \tTGraph.Warning#TPC.QA.dcar_posA_2 \tTF1.Status, .Status#obj1, #obj3\t .deadBand";

  if (!AliDrawStyle::IsSelected(selectors, "TH1", "deadBand", "")) {
    ::Info("AliDrawStyleTest","!AliDrawStyle::IsSelected(\"%s\",\"TH1\", \"deadBand\", \"\")- IsOK", selectors.Data());
  } else{
    ::Error("AliDrawStyleTest","!AliDrawStyle::IsSelected(\"%s\",\"TH1\", \"deadBand\", \"\")- FAILED", selectors.Data());
  }
  if(AliDrawStyle::IsSelected(selectors, "TGraph", "deadBand", "")) {
    ::Info("AliDrawStyleTest","AliDrawStyle::IsSelected(\"%s\",\"TGraph\", \"deadBand\", \"\")- IsOK", selectors.Data());
  } else{
    ::Error("AliDrawStyleTest","AliDrawStyle::IsSelected(\"%s\",\"TGraph\", \"deadBand\", \"\")- FAILED", selectors.Data());
  }
  if(AliDrawStyle::IsSelected(selectors, "", "deadBand", "")) {
    ::Info("AliDrawStyleTest","AliDrawStyle::IsSelected(\"%s\",\"\", \"deadBand\", \"\")- IsOK", selectors.Data());
  } else{
    ::Error("AliDrawStyleTest","AliDrawStyle::IsSelected(\"%s\",\"\", \"deadBand\", \"\")- FAILED", selectors.Data());
  }
  if(AliDrawStyle::IsSelected(selectors, "", "deadBand2", "")) {
    ::Info("AliDrawStyleTest","AliDrawStyle::IsSelected(\"%s\",\"\", \"deadBand2\", \"\")- IsOK", selectors.Data());
  } else{
    ::Error("AliDrawStyleTest","AliDrawStyle::IsSelected(\"%s\",\"\", \"deadBand2\", \"\")- FAILED", selectors.Data());
  }
  if(!AliDrawStyle::IsSelected(selectors, "TH1", "deadBand2", "")) {
    ::Info("AliDrawStyleTest","!AliDrawStyle::IsSelected(\"%s\",\"TH1\", \"deadBand2\", \"\")- IsOK", selectors.Data());
  } else{
    ::Error("AliDrawStyleTest","!AliDrawStyle::IsSelected(\"%s\",\"TH1\", \"deadBand2\", \"\")- FAILED", selectors.Data());
  }
  if(AliDrawStyle::IsSelected(selectors, "TGraph", "Status", "")) {
    ::Info("AliDrawStyleTest","AliDrawStyle::IsSelected(\"%s\",\"TGraph\", \"Status\", \"\")- IsOK", selectors.Data());
  } else{
    ::Error("AliDrawStyleTest","AliDrawStyle::IsSelected(\"%s\",\"TGraph\", \"Status\", \"\")- FAILED", selectors.Data());
  }
  if(!AliDrawStyle::IsSelected(selectors, "TGraph", "Warning", "obj3")) {
    ::Info("AliDrawStyleTest","!AliDrawStyle::IsSelected(\"%s\",\"TGraph\", \"Warning\", \"obj3\")- IsOK", selectors.Data());
  } else{
    ::Error("AliDrawStyleTest","!AliDrawStyle::IsSelected(\"%s\",\"TGraph\", \"Warning\", \"obj3\")- FAILED", selectors.Data());
  }
}

void AliDrawStyleTest_GetIds() {
  TString elementID = "";
  TString classID   = "";
  TString objectID  = "";
  TString styleID   = "";

  TH1F *his = new TH1F("his1.style(345)", "title", 100,-5,5);
  TObject *obj = (TObject *) his;
  AliDrawStyle::GetIds(obj, elementID, classID, objectID, styleID);
  if(elementID == TString("TH1F") && classID == TString("") && objectID == TString("his1") && styleID == TString("345")) {
    ::Info("AliDrawStyleTest","AliDrawStyle::GetIds(\"%s\",\"%s\", \"%s\", \"%s\", \"%s\")- IsOK", obj->GetName(), elementID.Data(), classID.Data(), objectID.Data(), styleID.Data());
  } else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetIds(\"%s\",\"%s\", \"%s\", \"%s\", \"%s\")- FAILED", obj->GetName(), elementID.Data(), classID.Data(), objectID.Data(), styleID.Data());
  }

  TH1F *his8 = new TH1F("his1.class(123)", "title", 100,-5,5);
  obj = (TObject *) his8;
  AliDrawStyle::GetIds(obj, elementID, classID, objectID, styleID);
  if(elementID == TString("TH1F") && classID == TString("123") && objectID == TString("his1") && styleID == TString("")) {
    ::Info("AliDrawStyleTest","AliDrawStyle::GetIds(\"%s\",\"%s\", \"%s\", \"%s\", \"%s\")- IsOK", obj->GetName(), elementID.Data(), classID.Data(), objectID.Data(), styleID.Data());
  } else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetIds(\"%s\",\"%s\", \"%s\", \"%s\", \"%s\")- FAILED", obj->GetName(), elementID.Data(), classID.Data(), objectID.Data(), styleID.Data());
  }

  elementID = "";
  classID   = "";
  objectID  = "";
  styleID   = "";
  TH1F *his2 = new TH1F("his1.class(123).style(345)", "title", 100,-5,5);
  obj = (TObject *) his2;
  AliDrawStyle::GetIds(obj, elementID, classID, objectID, styleID);
  if(elementID == TString("TH1F") && classID == TString("123") && objectID == TString("his1") && styleID == TString("345")) {
    ::Info("AliDrawStyleTest","AliDrawStyle::GetIds(\"%s\",\"%s\", \"%s\", \"%s\", \"%s\")- IsOK", obj->GetName(), elementID.Data(), classID.Data(), objectID.Data(), styleID.Data());
  } else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetIds(\"%s\",\"%s\", \"%s\", \"%s\", \"%s\")- FAILED", obj->GetName(), elementID.Data(), classID.Data(), objectID.Data(), styleID.Data());
  }

  elementID = "";
  classID   = "";
  objectID  = "";
  styleID   = "";
  TH1F *his3 = new TH1F("his1.TPC.RAW.class(123).style(345)", "title", 100,-5,5);
  obj = (TObject *) his3;
  AliDrawStyle::GetIds(obj, elementID, classID, objectID, styleID);
  if(elementID == TString("TH1F") && classID == TString("123") && objectID == TString("his1.TPC.RAW") && styleID == TString("345")) {
    ::Info("AliDrawStyleTest","AliDrawStyle::GetIds(\"%s\",\"%s\", \"%s\", \"%s\", \"%s\")- IsOK", obj->GetName(), elementID.Data(), classID.Data(), objectID.Data(), styleID.Data());
  } else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetIds(\"%s\",\"%s\", \"%s\", \"%s\", \"%s\")- FAILED", obj->GetName(), elementID.Data(), classID.Data(), objectID.Data(), styleID.Data());
  }

  elementID = "";
  classID   = "";
  objectID  = "";
  styleID   = "";
  TH1F *his4 = new TH1F("his1.TPC.RAW.style(345).class(123)", "title", 100,-5,5);
  obj = (TObject *) his4;
  AliDrawStyle::GetIds(obj, elementID, classID, objectID, styleID);
  if(elementID == TString("TH1F") && classID == TString("123") && objectID == TString("his1.TPC.RAW") && styleID == TString("345")) {
    ::Info("AliDrawStyleTest","AliDrawStyle::GetIds(\"%s\",\"%s\", \"%s\", \"%s\", \"%s\")- IsOK", obj->GetName(), elementID.Data(), classID.Data(), objectID.Data(), styleID.Data());
  } else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetIds(\"%s\",\"%s\", \"%s\", \"%s\", \"%s\")- FAILED", obj->GetName(), elementID.Data(), classID.Data(), objectID.Data(), styleID.Data());
  }

  elementID = "";
  classID   = "";
  objectID  = "";
  styleID   = "";
  TH1F *his5 = new TH1F("his.TPC.RAW[1].style(345).class(123)", "title", 100,-5,5);
  obj = (TObject *) his5;
  AliDrawStyle::GetIds(obj, elementID, classID, objectID, styleID);
  if(elementID == TString("TH1F") && classID == TString("123") && objectID == TString("his.TPC.RAW") && styleID == TString("345")) {
    ::Info("AliDrawStyleTest","AliDrawStyle::GetIds(\"%s\",\"%s\", \"%s\", \"%s\", \"%s\")- IsOK", obj->GetName(), elementID.Data(), classID.Data(), objectID.Data(), styleID.Data());
  } else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetIds(\"%s\",\"%s\", \"%s\", \"%s\", \"%s\")- FAILED", obj->GetName(), elementID.Data(), classID.Data(), objectID.Data(), styleID.Data());
  }

  elementID = "";
  classID   = "";
  objectID  = "";
  styleID   = "";
  TH1F *his6 = new TH1F("his[1]", "title", 100,-5,5);
  obj = (TObject *) his6;
  AliDrawStyle::GetIds(obj, elementID, classID, objectID, styleID);
  if(elementID == TString("TH1F") && classID == TString("") && objectID == TString("his") && styleID == TString("")) {
    ::Info("AliDrawStyleTest","AliDrawStyle::GetIds(\"%s\",\"%s\", \"%s\", \"%s\", \"%s\")- IsOK", obj->GetName(), elementID.Data(), classID.Data(), objectID.Data(), styleID.Data());
  } else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetIds(\"%s\",\"%s\", \"%s\", \"%s\", \"%s\")- FAILED", obj->GetName(), elementID.Data(), classID.Data(), objectID.Data(), styleID.Data());
  }

  elementID = "";
  classID   = "";
  objectID  = "";
  styleID   = "";
  TH1F *his7 = new TH1F("his7", "title", 100,-5,5);
  obj = (TObject *) his7;
  AliDrawStyle::GetIds(obj, elementID, classID, objectID, styleID);
  if(elementID == TString("TH1F") && classID == TString("") && objectID == TString("his7") && styleID == TString("")) {
    ::Info("AliDrawStyleTest","AliDrawStyle::GetIds(\"%s\",\"%s\", \"%s\", \"%s\", \"%s\")- IsOK", obj->GetName(), elementID.Data(), classID.Data(), objectID.Data(), styleID.Data());
  } else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetIds(\"%s\",\"%s\", \"%s\", \"%s\", \"%s\")- FAILED", obj->GetName(), elementID.Data(), classID.Data(), objectID.Data(), styleID.Data());
  }

}

/// To test  - input CSS file to be read and than written
///          - diff between the files should be 0 except of the formatting
/// TODO test - ignoring commented fields in selector and in the declaration
// void AliDrawStyleTest_CSSReadWrite() {
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
//   for (Int_t i = 0; i < cssArray->GetEntriesFast(); i++) {
//     if ((TString(cssArray->At(i)->GetName()).ReplaceAll("\n", "") == TString(cssArrayFromTest->At(i)->GetName()).ReplaceAll("\n", "")) && TString(cssArray->At(i)->GetTitle()).ReplaceAll("\n", "") == TString(cssArrayFromTest->At(i)->GetTitle()).ReplaceAll("\n", "")) ::Info("AliDrawStyleTest","AliDrawStyleTest_CSSReadWrite()- ");
//     else ::Info("AliDrawStyleTest","AliDrawStyleTest_CSSReadWrite()- ");
//
//   }
//   gSystem->GetFromPipe("rm -f $AliRoot_SRC/STAT/test/test.css");
//
// }

void  AliDrawStyleTest_GetValue() {
  AliDrawStyle::RegisterCssStyle("alirootTestStyle.css",AliDrawStyle::ReadCSSFile("$AliRoot_SRC/STAT/test/alirootTestStyle.css",0));
  if (AliDrawStyle::GetValue("alirootTestStyle.css","marker-size", "TGraph", "Status", "TPC.QA.dcar_posA_1", "") == "1,2,3,4") {
    ::Info("AliDrawStyleTest","AliDrawStyle::GetValue(\"alirootTestStyle.css\",\"marker-size\", \"TGraph\", \"Status\", \"TPC.QA.dcar_posA_1\", \"\")- IsOK");
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetValue(\"alirootTestStyle.css\",\"marker-size\", \"TGraph\", \"Status\", \"TPC.QA.dcar_posA_1\", \"\")- FAILED");
  }
  if (AliDrawStyle::GetValue("alirootTestStyle.css","marker-size", "TF1", "Status", "obj4", "") == "17,18,19,20") {
    ::Info("AliDrawStyleTest","AliDrawStyle::GetValue(\"alirootTestStyle.css\",\"marker-size\", \"TF1\", \"Status\", \"obj4\", \"\")- IsOK");
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetValue(\"alirootTestStyle.css\",\"marker-size\", \"TF1\", \"Status\", \"obj4\", \"\")- FAILED");
  }
  if (AliDrawStyle::GetValue("alirootTestStyle.css","line-color", "TGraphErrors", "Warning", "asdasobj56", "") == "41,42,43,44") {
    ::Info("AliDrawStyleTest","AliDrawStyle::GetValue(\"alirootTestStyle.css\",\"line-color\", \"TGraphErrors\", \"Warning\", \"asdasobj56\", \"\")- IsOK");
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetValue(\"alirootTestStyle.css\",\"line-color\", \"TGraphErrors\", \"Warning\", \"asdasobj56\", \"\")- FAILED");
  }
  if (AliDrawStyle::GetValue("alirootTestStyle.css","marker-color", "TObject", "SomeNotExistingStatus", "obj3", "") == "37,38,39,40") {
    ::Info("AliDrawStyleTest","AliDrawStyle::GetValue(\"alirootTestStyle.css\",\"marker-color\", \"SomeNotExistingClass\", \"SomeNotExistingStatus\", \"obj3\", \"\")- IsOK");
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetValue(\"alirootTestStyle.css\",\"marker-color\", \"SomeNotExistingClass\", \"SomeNotExistingStatus\", \"obj3\", \"\")- FAILED");
  }
  if (AliDrawStyle::GetValue("alirootTestStyle.css","line-color", "TH1", "Warning", "obj1", "") == "57,58,59,60") {
    ::Info("AliDrawStyleTest","AliDrawStyle::GetValue(\"alirootTestStyle.css\",\"line-color\", \"TH1\", \"Warning\", \"obj1\", \"\")- IsOK");
  }else{
    ::Error("AliDrawStyleTest","AliDrawStyle::GetValue(\"alirootTestStyle.css\",\"line-color\", \"TH1\", \"Warning\", \"obj1\", \"\")- FAILED");
  }
}

void AliDrawStyleTest_ApplyCssStyle() {
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

TCanvas *MakeTestPlot(Int_t nHis) {
  TRandom r;
  TCanvas *exampleCanvas = new TCanvas("c1", "The AliDrawStyle::ApplyCssStyle example", 200, 10, 1200, 900);
  exampleCanvas->Divide(1,2);
  //AliPainter::DivideTPad("<vertical>[1l,1r]", "Raw", "gridX:0;gridY:0;tickX:0;tickY:0;", exampleCanvas);
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
 // gPad->BuildLegend();

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
    grArray[i]->SetName(TString::Format("graph%d.class(Raw).style(marker-size:%d;marker-color:%d;)", i,i/3 + 1, i + 1).Data());
    grArray[i]->SetTitle(TString::Format("gr%d", i).Data());
    if (i == 0) {
      grArray[i]->SetMinimum(0);
      grArray[i]->Draw("alp");
    }
    else grArray[i]->Draw("lp");
  }
 // gPad->BuildLegend();
  return exampleCanvas;
}

void AliDrawStyleTest_GenerateDoxyImages() {
  TCanvas *canvasT = new TCanvas("canvasT", "canvasT", 900, 500);
  canvasT->Divide(1,1);
  canvasT->cd(1);
  TH1F *his1 = new TH1F("his1.class(Error).style(line-color:#f30000;)", "his1", 100, -5, 5);
  his1->FillRandom("gaus", 15000);
  his1->SetStats(kFALSE);
  TH1D *his2 = new TH1D("his2.class(Error)", "his2", 100, -5, 5);
  his2->FillRandom("gaus", 10000);
  TH1I *his3 = new TH1I("his3.class(Error)", "his3", 100, -5, 5);
  his3->FillRandom("gaus", 5000);
  his1->Draw();
  his2->Draw("same");
  his3->Draw("same");
  gStyle->SetOptTitle(0);
  TPaveText *pt = new TPaveText(-0.438183, 694.575009, 0.438183, 740.053135);
  pt->AddText("Example of styling with using AliDrawStyle");
  pt->SetTextSize(0.04);
  pt->SetShadowColor(0);
  pt->SetBorderSize(0);
  pt->SetFillColor(0);
  pt->Draw();
  gPad->BuildLegend();
  canvasT->SaveAs("$AliRoot_SRC/STAT/imgdoc/AliDrawStyle_h_example1.png");
  canvasT->Clear();

  canvasT->Divide(1,1);
  canvasT->cd(1);
  his1->Draw();
  his2->Draw("same");
  his3->Draw("same");
  gStyle->SetOptTitle(0);
  pt->Draw();
  gPad->BuildLegend();

  AliDrawStyle::RegisterCssStyle("AliDrawStyleTutor", AliDrawStyle::ReadCSSFile("$AliRoot_SRC/STAT/test/AliDrawStyleTutor.css"));
  AliDrawStyle::ApplyCssStyle(canvasT, "AliDrawStyleTutor");
  canvasT->SaveAs("$AliRoot_SRC/STAT/imgdoc/AliDrawStyle_h_example2.png");
  canvasT->Clear();
 }