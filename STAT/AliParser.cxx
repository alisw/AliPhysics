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

#include "AliParser.h"
#include <vector>
#include <map>
#include <iostream>
#include "TString.h"
#include "TError.h"
#include "TMatrixD.h"

/// \brief  Extracts content between specified patterns (startStr, endStr).
///
/// The analogue of this method [exists in matlab](https://www.mathworks.com/help/matlab/ref/extractbetween.html)
/// \param str - input string
/// \param startStr - start pattern
/// \param endStr - end pattern
/// \param verbose
/// \return - vector of TString with content between startStr and endStr
///
/// \b Example \b usage:
/// \code
/// root[] AliParser::ExtractBetween("{\nmarker-style:25,21,22,23; \nmarker-color:1,2,4,5; \n}", "marker-style:", ";")
/// (std::vector<TString>) { "25,21,22,23" }
/// root [] AliParser::ExtractBetween("3*<mean>-<max>/2,<rms>*10+<min>", "<", ">")
/// (std::vector<TString>) { "mean", "max", "rms", "mean" }
/// \endcode
std::vector<TString> AliParser::ExtractBetween(const char *str, const char *startStr, const char *endStr, Int_t verbose) {
  TString iTStr = TString(str);
  Int_t startPos = iTStr.Index(startStr);
  Int_t pos=0;
  Int_t endPos = iTStr.Index(endStr, startPos);
  std::vector<TString> res;

  while (startPos != -1) {
    pos = startPos+TString(startStr).Length();
    res.push_back(iTStr(pos, endPos - pos));
    startPos = iTStr.Index(startStr, endPos);
    endPos = iTStr.Index(endStr, startPos+1);
  }

  if (verbose == 4) {
    TString infoString = "";
    for (std::vector<TString>::iterator it = res.begin(); it != res.end(); ++it)
      infoString += *it + "|";
    ::Info("AliParser::ParseString", "Input string \"%s\" was parsed to %s", iTStr.Data(), infoString.Data());
  }
  return res;
}


/// \brief  Splits input string to array according to specified char delimiter.
///
///  Widespread function in popular languages like [python](https://docs.python.org/2/library/stdtypes.html#str.split), or [javascript](https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/String/split).
///  The main difference our method from aforementioned  is ignoring specified delimiter inside of such - "(){}[]<>" brackets.
///  Template of input string: <value1|delimiter|value2(value2.1|delimiter|value2.2)> ==> {value1,value2(value2.1|delimiter|value2.2)}
/// \param inputExpr - input string
/// \param del - delimiter. by default - ','
/// \param verbose
/// \return - array with parsed content
///
/// \b Example \b usage:
/// \code
///   root [] AliParser::Split("a(1,2,3),b[1,2,3],c{1,2,3}")
///   (std::vector<TString>) { "a(1,2,3)", "b[1,2,3]", "c{1,2,3}"}
///   root [] AliParser::Split("1,2,3,4,5")
///   (std::vector<TString>) { "1", "2", "3", "4", "5"}
/// \endcode
std::vector<TString> AliParser::Split(const char *inputExpr, const char del, Int_t verbose) {
  std::vector<TString> res;
  TString inputTStr(inputExpr);
  Int_t startIndex = 0;
  std::map<char, char> ignoredBrackets;
  ignoredBrackets.insert(std::make_pair('(', ')'));
  ignoredBrackets.insert(std::make_pair('[', ']'));
  ignoredBrackets.insert(std::make_pair('{', '}'));
  ignoredBrackets.insert(std::make_pair('<', '>'));

  for (UShort_t i = 0; i <= inputTStr.Length(); ++i) {
    if (inputTStr(i) == TString(del) || i == inputTStr.Length()) {
      res.push_back(TString(inputTStr(startIndex, i - startIndex)));
      startIndex = i + 1;
    } else if (ignoredBrackets.find(inputTStr(i)) != ignoredBrackets.end()) {
      i = inputTStr.Index(ignoredBrackets[inputTStr(i)], i);
      continue;
    }
  }

  if (verbose == 4) {
    TString infoString = "";
    for (std::vector<TString>::iterator it = res.begin(); it != res.end(); ++it)
      infoString += *it + "|";
    ::Info("AliPainter::ParseString", "Input string \"%s\" was parsed to %s", inputExpr, infoString.Data());
  }
  return res;
}

/// \brief Extracts content from specified parentheses.
///
///  The similar method exist in java library from google - [guava](https://github.com/google/guava/issues/1615)
///  Template of input string: <(content1)(content2(content3))> ==> {content1,content2(content3)}
///  This function very similar to AliParser::ExtractBetween. The differences:
///  1. In AliParser::ExtractBetween separators could be any length;
///  2. In case such string "abc<1>+..." will be input string of both methods, that all elements of returned arrays
///     from AliParser::ExtractBetween and AliParser::ExtractSurroundingBy will be the same except first element.
///     AliParser::ExtractBetween(...)[0] will be "1"(!), but AliParser::ExtractSurroundingBy(...)[0] will be "abc"(!).
/// \param inpString
/// \param begin - by default - '('
/// \param end - by default - ')'
/// \param verbose
/// \return - arrays with contents of parentheses
///
/// \b Example \b usage:
/// \code
///  root [] AliParser::ExtractSurroundingBy("(a(1,2,3))(b[1,2,3])(c{1,2,3})")
///  (std::vector<TString>) {  "a(1,2,3)", "b[1,2,3]", "c{1,2,3}" }
///  AliParser::ExtractSurroundingBy("nameOfSomething(param1,param2)(param3,param4)(functionName(arg1,arg2))")
/// (std::vector<TString>) { "nameOfSomething", "param1,param2", "param3,param4", "functionName(arg1,arg2)" }
/// \endcode
std::vector<TString> AliParser::ExtractSurroundingBy(const char *inpString, const char begin, const char end, Int_t verbose) {
  std::vector<TString> parsedArr;
  TString exprsn(inpString);
  if (exprsn.CountChar(begin) != exprsn.CountChar(end)) {
    ::Error("AliPainter::DrawHistogram", "check brackets in %s", exprsn.Data());
    return parsedArr;
  }

  if (exprsn.Index(begin) != 0)
    parsedArr.push_back(TString(exprsn(0, exprsn.Index(begin))));
  else
    parsedArr.push_back(TString());

  TString verbStr = "";
  Int_t match = 0, startIndex = 0, finishIndex = 0;
  Bool_t isChange = kFALSE;

  for (Int_t i = 0; i < exprsn.Length(); ++i) {
    if (exprsn(i) == begin && match == 0) {
      match++;
      startIndex = i;
      isChange = kTRUE;
    } else if (exprsn(i) == begin && match > 0) match++;
    else if (exprsn(i) == end && match == 1) {
      match--;
      finishIndex = i;
    } else if (exprsn(i) == end && match > 1) match--;

    if (match == 0 && isChange) {
      parsedArr.push_back(TString(exprsn(startIndex + 1, finishIndex - startIndex - 1)));
      isChange = kFALSE;
    }
  }

  if (verbose == 4) {
    TString infoString = "";
    for (std::vector<TString>::iterator it = parsedArr.begin(); it != parsedArr.end(); ++it)
      infoString += *it + "|";
    ::Info("AliPainter::ParseString", "Input string \"%s\" was parsed to %s", inpString, infoString.Data());
  }
  return parsedArr;
}

//TODO: add also support to another options like -b, -n10, -n 10, --tools "...", ...
/// \brief Parses string with named arguments.
///
///        Template of input string: <nameOfArg1=value1,nameOfArg2=value2,nameOfArg3=[value3.1,value3.2]...>
///        In some cases could be useful usage of named parameters. You can study more information in [python argparse module.](https://docs.python.org/3/library/argparse.html)
/// \param optionsStr - input string
/// \param verbose
/// \param defKeys - the list of predefined keys, if parsed key will not found in this list - warning'll be generate.
/// \return
///
/// \b Example \b usage:
/// \code
///  root[] AliParser::Parse("name=gaus, strategy=misac(10,20),xlim=[10,20],class=Abc",4)
///  (std::map<TString, TString>) { "name" => "gaus", "strategy" => "misac(10,20)", "xlim" => "[10,20]", "class" => "Abc" }
/// \endcode
std::map<TString, TString> AliParser::Parse(const char *iStr, Int_t verbose, std::vector<TString> defKeys) {
  std::map<TString, TString> optMap;
  TString str(iStr);
  std::vector<TString> options = AliParser::Split(iStr, ',', verbose);
  for (UShort_t i = 0; i < options.size(); i++) {
    TString optionStr = options[i];
    TString key = "";
    TString value = "";
    key = TString(optionStr(0, optionStr.Index("="))).ReplaceAll(" ", "");
    value = TString(optionStr(optionStr.Index("=") + 1, optionStr.Length())).ReplaceAll(" ", "");
//    if (std::find(defKeys.begin(), defKeys.end(), key) == defKeys.end() && key != TString() && defKeys.size() > 0) {
//      TString defaultKeys = "";
//      for (std::vector<TString>::iterator it = defKeys.begin(); it != defKeys.end(); ++it)
//        defaultKeys += *it + ",";
//      ::Warning("AliPainter::DrawHistogram", "key \"%s\" not found in the list of default keys: \"%s\"", key.Data(),
//                defaultKeys.Data());
//    }
    optMap[key] = value;
  }
  return optMap;
}

/// \brief Returns array according with python-like interface.
///
///        In python for getting values from array, you can use ":" [see python docs](https://docs.python.org/3.7/library/functions.html?highlight=slice#slice)
///        Here we provide the same functionality, but output array will generate from string.
/// \param range
/// \return
///
/// \b Example \b usage:
/// \code
///  root []  AliParser::Slice2IArray("10:40:10:10")
///    (std::vector<Int_t>) { 10, 20, 20, 30, 30, 40 }
///  root []  AliParser::Slice2IArray("1:4")
///    (std::vector<Int_t>) { 1, 2, 2, 3, 3, 4 }
///  root [] AliParser::Slice2IArray("1:4:1:0")
///    (std::vector<Int_t>) { 1, 1, 2, 2, 3, 3, 4, 4 }
/// \endcode
std::vector<Int_t> AliParser::Slice2IArray(const char *inputString) {
  Int_t initArr[4] = {0, 0, 1, 1}; // start, stop, step, delta
  std::vector<TString> initValues = AliParser::Split(TString(inputString), ':');
  for (Int_t i=0; i < 4 && i < (Int_t) initValues.size(); ++i)
    initArr[i] = initValues[i].Atoi();
  std::vector<Int_t> vRanges;
  for (Int_t j = initArr[0]; j <= initArr[1] - initArr[3]; j += initArr[2]) {
    vRanges.push_back(j);
    vRanges.push_back(j + initArr[3]);
  }
  return  vRanges;
}

// TODO: code will crash if iRanges.size() is not even. Crash when (it+1) row 229
//::Error("AliParser", "SliceRanges: count of values should be even.");
//return TMatrixD();
/// \brief Returns TMatrixD from input string.
///
///   The idea is providing method for generating values of ranges with slicer:
///   Such string "10,20,30:60:10:10,70:110:10:10" should be transform to the next map:
///  |axisNum  |   0    |   1    |   2    |
///  |:-------:|:-----:|:-----:|:-----:|
///  |values   |   10   |   30   |   70   |
///  |          |   20   |   40   |   80   |
///  |          |        |   40   |   80   |
///  |          |        |   50   |   90   |
///  |          |        |   50   |   90   |
///  |          |        |   60   |   100  |
///  |          |        |        |   100  |
///  |          |        |        |   110  |
/// And then such map will be transform to the matrix:
///
/// \b nCols - (count of axes) * 2
///
/// \b nRows - 1 + multiplication of sizes of each array from map divided by 2
///
/// \b last \b row - floatFlag of each column. In case float flag is 1 SetRangeUser will be applied else SetRange
///
///  |x0Min|x0Max|x1Min|x1Max|x2Min|x2Max|
///  |:---:|:---:|:---:|:---:|:---:|:---:|
///  | 10  |  20 |  30 |  40 |  70 |  80 |
///  |  10 |  20 |  30 |  40 |  80 |  90 |
///  |  10 |  20 |  30 |  40 |  90 |  100|
///  |  10 |  20 |  30 |  40 |  100|  110|
///  |  10 |  20 |  40 |  50 |  70 |  80 |
///  |  10 |  20 |  40 |  50 |  80 |  90 |
///  |  10 |  20 |  40 |  50 |  90 |  100|
///  |  10 |  20 |  40 |  50 |  100|  110|
///  |  10 |  20 |  50 |  60 |  70 |  80 |
///  |  10 |  20 |  50 |  60 |  80 |  90 |
///  |  10 |  20 |  50 |  60 |  90 |  100|
///  |  10 |  20 |  50 |  60 |  100|  110|
///  |  0  |  0  |  0  |  0  |  0  |  0  |
///
/// \param inputString - string for transforming
/// \param verbose
/// \return matrix of parsed values
///
/// \b Example \b usage:
/// \code
/*
auto m = AliParser::Slice2Matrix("10,20,30:60:10:10,70:110:10:10")
{
  for(auto i=0;i<m.GetNrows();++i) {
    for(auto j=0;j<m.GetNcols();++j) {
      std::cout << " " << TMatrixDRow(m,i)[j];
    }
    std::cout << std::endl;
  }
}
 */
/// \endcode
TMatrixD AliParser::Slice2Matrix(const char *iStr, Int_t verbose) {
  TString inputString(iStr);
  if (inputString == TString()) return TMatrixD();
  std::vector<TString> initRanges = AliParser::Split(iStr, ',', verbose);
  std::map<Int_t, std::vector<Double_t> > iRanges;
  std::vector<Double_t> tempVector;
  std::vector<Int_t> sliceVec;
  std::vector<Double_t> floatFlagArray;

  Int_t axisNum = 0;
  std::vector<TString>::iterator it = initRanges.begin();
  try {
    while (it < initRanges.end()) {
      tempVector.clear();
      sliceVec.clear();
      if (it->Contains(':')) {
        sliceVec = AliParser::Slice2IArray(it->Data());
        tempVector.insert(tempVector.end(), sliceVec.begin(), sliceVec.end());
        AliParser::FillFloatFlagArray(floatFlagArray, *it);
        AliParser::FillFloatFlagArray(floatFlagArray, *it);
        it++;
      } else {
        tempVector.push_back(it->Atof());
        tempVector.push_back((it + 1)->Atof());
        AliParser::FillFloatFlagArray(floatFlagArray, *it);
        AliParser::FillFloatFlagArray(floatFlagArray, *(it + 1));
        it += 2;
      }
      iRanges[axisNum] = tempVector;
      axisNum++;
    }
  }
  catch (std::exception &e) {
    std::cerr << "Exception catched : " << e.what() << std::endl;
    std::cout << "Most probably you have wrong number of range values. It should be even number." << std::endl;
    return TMatrixD();
  }
  Int_t *ind = new Int_t[iRanges.size()]();
  std::vector<Double_t> darr;
  AliParser::Map2Array(darr, iRanges, ind, 1);
  delete [] ind;
  darr.insert(darr.end(), floatFlagArray.begin(), floatFlagArray.end());
  Int_t rowCnt = 1;
  for (Int_t r = 0; r < (Int_t) iRanges.size(); ++r)
    rowCnt *= iRanges[r].size() / 2;
  TMatrixD matrix(rowCnt + 1, iRanges.size() * 2, &darr[0]);
  return matrix;
}

///
/// \param array
/// \param iRanges
/// \param indexes
/// \param cnt
void AliParser::Map2Array(std::vector<Double_t> &array, std::map<Int_t, std::vector<Double_t> > iRanges, Int_t *indexes, Int_t cnt) {
  std::vector<Double_t> tempArray;
  Int_t d = 0;
  for (d = 0; d < (Int_t) iRanges.size(); ++d) {
    tempArray.push_back(iRanges[d][indexes[d]]);
    tempArray.push_back(iRanges[d][indexes[d] + 1]);
  }
  array.insert(array.end(), tempArray.begin(), tempArray.end());
  for (Int_t i = d - 1; i >= 0; --i) {
    if (indexes[i] + 2 >= (Int_t) iRanges[i].size())
      indexes[i] = 0;
    else {
      indexes[i] = indexes[i] + 2;
      break;
    }
  }
  Int_t rowCnt = 1;
  for (Int_t r = 0; r < (Int_t) iRanges.size(); ++r)
    rowCnt *= iRanges[r].size()/2;
  if (cnt == rowCnt) return;
  cnt++;
  AliParser::Map2Array(array, iRanges, indexes, cnt);
}

///
/// \param floatFlagArray
/// \param str
void AliParser::FillFloatFlagArray(std::vector<Double_t> &floatFlagArray, TString str) {
  if (str.Contains('.')) {
    floatFlagArray.push_back(1.);
  }
  else {
    floatFlagArray.push_back(0.);
  }
}
