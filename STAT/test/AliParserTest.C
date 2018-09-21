/// \ingroup STAT/test
/// \brief  test methods of AliParser
/// Example usage
/*
\code
.L $AliRoot_SRC/STAT/test/AliParserTest.C+
AliParserTest();
root.exe -b -q  $AliRoot_SRC/STAT/test/AliParserTest.C+ 2>&1 | tee AliParserTest.log
\endcode
*/

#include "AliParser.h"
#include "TString.h"
#include "TMatrixD.h"
#include <vector>
#include <map>
#include "TError.h"
#include "Riostream.h"
#include "Rtypes.h"
#include "TSystem.h"

void AliParserTest_Split();
void AliParserTest_ExtractSurroundingBy();
void AliParserTest_ExtractBetween();
void AliParserTest_Parse();
void AliParserTest_Slice2IArray();
void AliParserTest_Slice2Matrix();
TString Vec2TString(std::vector<TString>);
TString Map2TString(std::map<TString, TString>);
TString VecD2TString(std::vector<Double_t> parsedVector);
TString VecI2TString(std::vector<Int_t> parsedVector);
TString TMatrixD2TString(TMatrixD matrix);


void AliParserTest() {
  AliParserTest_Split();
    AliParserTest_ExtractSurroundingBy();
  AliParserTest_Parse();
  AliParserTest_ExtractBetween();
  AliParserTest_Slice2IArray();
  AliParserTest_Slice2Matrix();
}

TString Vec2TString(std::vector<TString> parsedVector) {
  TString infoString = "";
  for (std::vector<TString>::iterator it = parsedVector.begin(); it != parsedVector.end(); ++it)
    infoString += *it + "|";
  return infoString;
}

TString VecD2TString(std::vector<Double_t> parsedVector) {
  TString infoString = "";
  for (std::vector<Double_t >::iterator it = parsedVector.begin(); it != parsedVector.end(); ++it)
    infoString += TString::Format("%.2g|",*it);
  return infoString;
}

TString VecI2TString(std::vector<Int_t> parsedVector) {
  TString infoString = "";
  for (std::vector<Int_t >::iterator it = parsedVector.begin(); it != parsedVector.end(); ++it)
    infoString += TString::Format("%d|",*it);
  return infoString;
}

TString TMatrixD2TString (TMatrixD matrix) {
  TString infoString = "\n";
  for(auto i=0;i<matrix.GetNrows();++i) {
    for(auto j=0;j<matrix.GetNcols();++j) {
      infoString += TString::Format("%d|", (Int_t) TMatrixDRow(matrix,i)[j]);
    }
    infoString += "\n";
  }
  return infoString;
}

TString Map2TString(std::map<TString, TString> parsedMap) {
  TString infoString = "";
  for (std::map<TString, TString>::iterator it = parsedMap.begin(); it != parsedMap.end(); it++)
    infoString += it->first + "=" + it->second + "|";
  return infoString;
}

void AliParserTest_Split() {
  TString input = "gaus,W,fitFunction(1,2,3),E,10,200";
  if (Vec2TString(AliParser::Split(input.Data())) == TString("gaus|W|fitFunction(1,2,3)|E|10|200|"))
    ::Info("AliParserTest", "AliParser::Split(\"%s\")- IsOK", input.Data());
  else
    ::Error("AliParserTest", "AliParser::Split(\"%s\")- FAILED", input.Data());

  input = "1,2,3,4";
  if (Vec2TString(AliParser::Split(input.Data())) == TString("1|2|3|4|"))
    ::Info("AliParserTest", "AliParser::Split(\"%s\")- IsOK", input.Data());
  else
    ::Error("AliParserTest", "AliParser::Split(\"%s\")- FAILED", input.Data());

  input = "a,b,c(1,2,3),d";
  if (Vec2TString(AliParser::Split(input.Data())) == TString("a|b|c(1,2,3)|d|"))
    ::Info("AliParserTest", "AliParser::Split(\"%s\")- IsOK", input.Data());
  else
    ::Error("AliParserTest", "AliParser::Split(\"%s\")- FAILED", input.Data());

  input = "a,b,c[1,2,3],d";
  if (Vec2TString(AliParser::Split(input.Data())) == TString("a|b|c[1,2,3]|d|"))
    ::Info("AliParserTest", "AliParser::Split(\"%s\")- IsOK", input.Data());
  else
    ::Error("AliParserTest", "AliParser::Split(\"%s\")- FAILED", input.Data());

  input = "(1,2,3),4,5,6";
  if (Vec2TString(AliParser::Split(input.Data())) == TString("(1,2,3)|4|5|6|"))
    ::Info("AliParserTest", "AliParser::Split(\"%s\")- IsOK", input.Data());
  else
    ::Error("AliParserTest", "AliParser::Split(\"%s\")- FAILED", input.Data());

}

void AliParserTest_ExtractSurroundingBy() {
  TString input = "[1][2][3][4]";
  if (Vec2TString(AliParser::ExtractSurroundingBy(input.Data(), '[', ']')) == TString("|1|2|3|4|"))
    ::Info("AliParserTest", "AliParser::ExtractSurroundingBy(\"%s\")- IsOK", input.Data());
  else
    ::Error("AliParserTest", "AliParser::ExtractSurroundingBy(\"%s\")- FAILED", input.Data());

  input = "(1)(2)(3)(4)";
  if (Vec2TString(AliParser::ExtractSurroundingBy(input.Data())) == TString("|1|2|3|4|"))
    ::Info("AliParserTest", "AliParser::ExtractSurroundingBy(\"%s\")- IsOK", input.Data());
  else
    ::Error("AliParserTest", "AliParser::ExtractSurroundingBy(\"%s\")- FAILED", input.Data());

  input = "(a(1,2,3))(b)(c)(d(4,5,6))";
  if (Vec2TString(AliParser::ExtractSurroundingBy(input.Data())) == TString("|a(1,2,3)|b|c|d(4,5,6)|"))
    ::Info("AliParserTest", "AliParser::ExtractSurroundingBy(\"%s\")- IsOK", input.Data());
  else
    ::Error("AliParserTest", "AliParser::ExtractSurroundingBy(\"%s\")- FAILED", input.Data());


  input = "((1,2,3))(2)(3)(4)";
  if (Vec2TString(AliParser::ExtractSurroundingBy(input.Data())) == TString("|(1,2,3)|2|3|4|"))
    ::Info("AliParserTest", "AliParser::ExtractSurroundingBy(\"%s\")- IsOK", input.Data());
  else
    ::Error("AliParserTest", "AliParser::ExtractSurroundingBy(\"%s\")- FAILED", input.Data());

  input = "abc(a)(b[10,20])(c[1,2,3])(d)";
  if (Vec2TString(AliParser::ExtractSurroundingBy(input.Data())) == TString("abc|a|b[10,20]|c[1,2,3]|d|"))
    ::Info("AliParserTest", "AliParser::ExtractSurroundingBy(\"%s\")- IsOK", input.Data());
  else
    ::Error("AliParserTest", "AliParser::ExtractSurroundingBy(\"%s\")- FAILED", input.Data());

}

void AliParserTest_ExtractBetween() {

  TString input = "pi*<max>+345*<min>";
  if (Vec2TString(AliParser::ExtractBetween(input.Data(), "<", ">", 0)) == TString("max|min|"))
    ::Info("AliParserTest", "AliParser::Split(\"%s\")- IsOK", input.Data());
  else
    ::Error("AliParserTest", "AliParser::Split(\"%s\")- FAILED", input.Data());

  input = "{marker-style:25,21,22,23; marker-color:1,2,4,5; }";
  if (Vec2TString(AliParser::ExtractBetween(input.Data(), "marker-style:", ";", 0)) == TString("25,21,22,23|"))
    ::Info("AliParserTest", "AliParser::ExtractBetween(\"%s\", \"marker-style:\")- IsOK", input.Data());
  else
    ::Error("AliParserTest", "AliParser::ExtractBetween(\"%s\", \"marker-style:\")- FAILED", input.Data());

  if (Vec2TString(AliParser::ExtractBetween(input.Data(), "marker-color:", ";", 0)) == TString("1,2,4,5|"))
    ::Info("AliParserTest", "AliParser::ExtractBetween(\"%s\", \"marker-color:\")- IsOK", input.Data());
  else
    ::Error("AliParserTest", "AliParser::ExtractBetween(\"%s\", \"marker-color:\")- FAILED", input.Data());

  input = "{marker-style:25,21,22,23; marker-color:rgb(1,2,3),#f0f0f0,1; }";
  if (Vec2TString(AliParser::ExtractBetween(input.Data(), "marker-color:", ";", 0)) == TString("rgb(1,2,3),#f0f0f0,1|"))
    ::Info("AliParserTest", "AliParser::ExtractBetween(\"%s\")- IsOK", input.Data());
  else
    ::Error("AliParserTest", "AliParser::ExtractBetween(\"%s\")- FAILED", input.Data());
}

void AliParserTest_Parse() {
  TString input = "div=0";
  if (Map2TString(AliParser::Parse(input.Data())) == TString("div=0|"))
    ::Info("AliParserTest", "AliParser::Parse(\"%s\")- IsOK", input.Data());
  else
    ::Error("AliParserTest", "AliParser::Parse(\"%s\")- FAILED", input.Data());

  input = "class=[Mass,PtAll]";
  if (Map2TString(AliParser::Parse(input.Data())) == TString("class=[Mass,PtAll]|"))
    ::Info("AliParserTest", "AliParser::Parse(\"%s\")- IsOK", input.Data());
  else
    ::Error("AliParserTest", "AliParser::Parse(\"%s\")- FAILED", input.Data());

  input = "div=0, class=[Mass,PtAll], xlim=[10.123,20.435]";
  if (Map2TString(AliParser::Parse(input.Data())) == TString("class=[Mass,PtAll]|div=0|xlim=[10.123,20.435]|"))
    ::Info("AliParserTest", "AliParser::Parse(\"%s\")- IsOK", input.Data());
  else
    ::Error("AliParserTest", "AliParser::Parse(\"%s\")- FAILED", input.Data());

  input = "ylim=";
  if (Map2TString(AliParser::Parse(input.Data())) == TString("ylim=|"))
    ::Info("AliParserTest", "AliParser::Parse(\"%s\")- IsOK", input.Data());
  else
    ::Error("AliParserTest", "AliParser::Parse(\"%s\")- FAILED", input.Data());

  input = "zlim=[],xlim=[10.123,20.435],ylim=,class=Mass";
  if (Map2TString(AliParser::Parse(input.Data())) == TString("class=Mass|xlim=[10.123,20.435]|ylim=|zlim=[]|"))
    ::Info("AliParserTest", "AliParser::Parse(\"%s\")- IsOK", input.Data());
  else
    ::Error("AliParserTest", "AliParser::Parse(\"%s\")- FAILED", input.Data());
  input = "zlim=(),xlim=(10.123,20.435),ylim=,class=Mass";
  if (Map2TString(AliParser::Parse(input.Data())) == TString("class=Mass|xlim=(10.123,20.435)|ylim=|zlim=()|"))
    ::Info("AliParserTest", "AliParser::Parse(\"%s\")- IsOK", input.Data());
  else
    ::Error("AliParserTest", "AliParser::Parse(\"%s\")- FAILED", input.Data());
}

void AliParserTest_Slice2IArray() {
  TString input = "10:40:10:10";
  if (VecI2TString(AliParser::Slice2IArray(input.Data())) == TString("10|20|20|30|30|40|"))
    ::Info("AliParserTest", "AliParser::Slice2IArray(\"%s\")- IsOK", input.Data());
  else
    ::Error("AliParserTest", "AliParser::Slice2IArray(\"%s\")- FAILED", input.Data());
  input = "1:4";
  if (VecI2TString(AliParser::Slice2IArray(input.Data())) == TString("1|2|2|3|3|4|"))
    ::Info("AliParserTest", "AliParser::Slice2IArray(\"%s\")- IsOK", input.Data());
  else
    ::Error("AliParserTest", "AliParser::Slice2IArray(\"%s\")- FAILED", input.Data());
  input = "1";
  if (VecI2TString(AliParser::Slice2IArray(input.Data())) == TString())
    ::Info("AliParserTest", "AliParser::Slice2IArray(\"%s\")- IsOK", input.Data());
  else
    ::Error("AliParserTest", "AliParser::Slice2IArray(\"%s\")- FAILED", input.Data());
  input = "";
  if (VecI2TString(AliParser::Slice2IArray(input.Data())) == TString())
    ::Info("AliParserTest", "AliParser::Slice2IArray(\"%s\")- IsOK", input.Data());
  else
    ::Error("AliParserTest", "AliParser::Slice2IArray(\"%s\")- FAILED", input.Data());
  input = "1:4:1:1:1:1:1";
  if (VecI2TString(AliParser::Slice2IArray(input.Data())) == TString("1|2|2|3|3|4|"))
    ::Info("AliParserTest", "AliParser::Slice2IArray(\"%s\")- IsOK", input.Data());
  else
    ::Error("AliParserTest", "AliParser::Slice2IArray(\"%s\")- FAILED", input.Data());
  input = "1:4:1";
  if (VecI2TString(AliParser::Slice2IArray(input.Data())) == TString("1|2|2|3|3|4|"))
    ::Info("AliParserTest", "AliParser::Slice2IArray(\"%s\")- IsOK", input.Data());
  else
    ::Error("AliParserTest", "AliParser::Slice2IArray(\"%s\")- FAILED", input.Data());
  input = "1:4:1:0";
  if (VecI2TString(AliParser::Slice2IArray(input.Data())) == TString("1|1|2|2|3|3|4|4|"))
    ::Info("AliParserTest", "AliParser::Slice2IArray(\"%s\")- IsOK", input.Data());
  else
    ::Error("AliParserTest", "AliParser::Slice2IArray(\"%s\")- FAILED", input.Data());
  input = "1.2:4.3";
  if (VecI2TString(AliParser::Slice2IArray(input.Data())) == TString("1|2|2|3|3|4|"))
    ::Info("AliParserTest", "AliParser::Slice2IArray(\"%s\")- IsOK", input.Data());
  else
    ::Error("AliParserTest", "AliParser::Slice2IArray(\"%s\")- FAILED", input.Data());

}

void AliParserTest_Slice2Matrix() {
  TString input = "10,20,30:50:10:10,60,70";
  if (TMatrixD2TString(AliParser::Slice2Matrix(input.Data())) == TString("\n10|20|30|40|60|70|\n10|20|40|50|60|70|\n0|0|0|0|0|0|\n"))
    ::Info("AliParserTest", "AliParser::Slice2Matrix(\"%s\")- IsOK", input.Data());
  else
    ::Error("AliParserTest", "AliParser::Slice2Matrix(\"%s\")- FAILED", input.Data());
  input = "10,20,30,50,60,70";
  if (TMatrixD2TString(AliParser::Slice2Matrix(input.Data())) == TString("\n10|20|30|50|60|70|\n0|0|0|0|0|0|\n"))
    ::Info("AliParserTest", "AliParser::Slice2Matrix(\"%s\")- IsOK", input.Data());
  else
    ::Error("AliParserTest", "AliParser::Slice2Matrix(\"%s\")- FAILED", input.Data());
//  input = "";
//  if (TMatrixD2TString(AliParser::Slice2Matrix(input.Data())) == TString("\n"))
//    ::Info("AliParserTest", "AliParser::Slice2Matrix(\"%s\")- IsOK", input.Data());
//  else
//    ::Error("AliParserTest", "AliParser::Slice2Matrix(\"%s\")- FAILED", input.Data());
}
