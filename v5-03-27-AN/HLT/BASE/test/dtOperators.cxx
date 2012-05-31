// $Id: AliHLTComponent.cxx 22571 2007-11-28 09:54:41Z richterm $

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   dtOperators.cxx
    @author Matthias Richter
    @date   
    @brief  Test program for data type operators
 */

#include <iostream>
#include "AliHLTDataTypes.h"
#include "AliHLTComponent.h"

using namespace std;

typedef struct test_t {
  AliHLTComponentDataType dt1;
  AliHLTComponentDataType dt2;
  int type; // 1 -> ==, 0 -> !=
  bool result;
} test_t;

int main(int /*argc*/, const char** /*argv*/)
{
  AliHLTComponentDataType testdt1={
    sizeof(AliHLTComponentDataType),
    {'D','U','M','M','Y','D','A','T'},
    {'T','E','S','T'}};

  AliHLTComponentDataType testdt2={
    sizeof(AliHLTComponentDataType),
    {'_','S','O','_','E','V','E','R'},
    {'W','H','A','T'}};

  AliHLTComponentDataType testdt3={
    sizeof(AliHLTComponentDataType),
    {'D','U','M','M','Y','D','A','T'},
    kAliHLTDataOriginAny};

  AliHLTComponentDataType testdt4={
    sizeof(AliHLTComponentDataType),
    kAliHLTAnyDataTypeID,
    {'T','E','S','T'}};

  AliHLTComponentDataType testdt5={
    sizeof(AliHLTComponentDataType),
    {'D','D','L','_','R','A','W',' '},
    {'T','R','D',' '}};

  test_t tests[] = {
    {kAliHLTDataTypeDDLRaw|"TRD " , testdt5            , 1 , true },
    {kAliHLTDataTypeDDLRaw|"TRD " , testdt5            , 0 , false},

    {kAliHLTDataTypeDDLRaw , testdt5            , 1 , true },
    {kAliHLTDataTypeDDLRaw , testdt5            , 0 , false},

//     {testdt1             , testdt1              , 1 , true },
//     {testdt1             , testdt1              , 0 , false},

    {testdt1             , testdt2              , 1 , false},
    {testdt1             , testdt2              , 0 , true },

    {testdt1             , testdt3              , 1 , true },
    {testdt1             , testdt3              , 0 , false},

    {testdt1             , testdt4              , 1 , true },
    {testdt1             , testdt4              , 0 , false},

    {testdt3             , testdt4              , 1 , true },
    {testdt3             , testdt4              , 0 , false},

    {testdt2             , testdt4              , 1 , false},
    {testdt2             , testdt4              , 0 , true },

    {testdt1             , kAliHLTAnyDataType   , 1 , true },
    {testdt1             , kAliHLTAnyDataType   , 0 , false},

    //
    {kAliHLTAnyDataType  , kAliHLTAnyDataType   , 1 , true },
    {kAliHLTAnyDataType  , kAliHLTAnyDataType   , 0 , false},

    {kAliHLTVoidDataType , kAliHLTAnyDataType   , 1 , false},
    {kAliHLTVoidDataType , kAliHLTAnyDataType   , 0 , true },

    {kAliHLTVoidDataType , kAliHLTVoidDataType  , 1 , true },
    {kAliHLTVoidDataType , kAliHLTVoidDataType  , 0 , false},

    {kAliHLTVoidDataType , testdt5  , 1 , false },
    {kAliHLTVoidDataType , testdt5  , 0 , true},

    {kAliHLTVoidDataType|"TRD " , kAliHLTVoidDataType  , 1 , false },
    {kAliHLTVoidDataType|"TRD " , kAliHLTVoidDataType  , 0 , true  },

    {kAliHLTVoidDataType|"TRD " , testdt5  , 1 , false },
    {kAliHLTVoidDataType|"TRD " , testdt5  , 0 , true },
  };

  int result=0;
  cout << "checking data type operators" << endl;
  for (unsigned int i=0; i<(sizeof(tests)/sizeof(test_t)); i++) {
    //if (!tests[i].type) continue;
    const char* op=(tests[i].type)?" == ":" != ";
    cout << "checking: " << AliHLTComponent::DataType2Text(tests[i].dt1).c_str() << op << AliHLTComponent::DataType2Text(tests[i].dt2).c_str() << " -> ";
    if (tests[i].type) result=tests[i].dt1==tests[i].dt2;
    else  result=tests[i].dt1!=tests[i].dt2;
    if (result) {
      if (tests[i].type) result=tests[i].dt2==tests[i].dt1;
      else  result=tests[i].dt2!=tests[i].dt1;
      if (!result) cout << " interchanged sequence ";
    }
    cout << result << endl;
    if (result^tests[i].result) return 1;
  }

  return 0;
}
