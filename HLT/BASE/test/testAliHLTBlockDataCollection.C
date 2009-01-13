// $Id$

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

/** @file   testAliHLTBlockDataCollection.C
    @author Matthias Richter
    @date   
    @brief  Test macro/program for the AliHLTBlockDataCollection
 */

#ifndef __CINT__
#include <ostream>
#include <vector>
#include "AliHLTBlockDataCollection.h"
#include "AliHLTComponent.h"
#else
#define EINVAL 22
#define EPROTO 71
#endif //__CINT__

typedef struct init_t {
  const char* arguments;
  int result;
} init_t;

typedef struct select_t {
  const char* id;
  const char* origin;
  UInt_t specification;
  int result;
} select_t;

init_t gInits[]={
  {"-datatype 'ESD_TREE' 'TPC ' ", 3},
  {"-datatype 'ESD_TREE' 'TPC ' "
   "-origin PHOS "
   "-verbosity "
   "-origin 'TRD ' "
   "-dataspec 0xdeadbeef ", 9},
  {"-datatype 'ESD_TREE'", -EPROTO},
  {NULL, 0}
};

select_t gSelections[]={
  {"ESD_TREE", "TPC ", 0, 1},
  {"DDL_RAW ", "TPC ", 0xaffe, 0},
  {"ESD_TREE", "TRD ", 0xaffe, 0},
  {"ESD_TREE", "TRD ", 0xdeadbeef, 1},
  {NULL, NULL, 0, 0}
};

int gVerbosity=0;

int InitCollection(AliHLTBlockDataCollection& collection, const char* init)
{
  // build argument vector
  int iResult=0;
  char* arguments=new char[strlen(init)+1];
  strcpy(arguments, init);
  vector<int> positions;
  bool bQuote=false;
  int i=0;
  int iStart=0;
  for (;arguments[i]!=0; i++) {
    if (arguments[i]=='\'') {
      if (bQuote=!bQuote) {
	// opening quote, set start
      } else {
	// closing quote, add argument
	arguments[i]=0;
	if (i-iStart>0) positions.push_back(iStart);
      }
      iStart=i+1;
    } else if ((arguments[i]==' ' && !bQuote) ||
	       arguments[i]==0) {
      arguments[i]=0;
      if (i-iStart>0) positions.push_back(iStart);
      iStart=i+1;
    }
  }
  if (i-iStart>0) positions.push_back(iStart);

  int argc=positions.size();
  const char** argv=new const char*[argc];
  for (int j=0; j<argc ; j++) {
    argv[j]=&arguments[positions[j]];
  }

  // test argument scan
  for (i=0; i<argc; i++) {
    if (gVerbosity>1) cout << "scanning " << argc-i << " arguments: " << argv[i] << endl;
    int result=collection.ScanArgument(argc-i, &argv[i]);
    if (result>0) {
      iResult+=result;
      i+=result-1;
    } else if (result<0) {
      iResult=result;
      break;
    }
  }

  // cleanup
  delete [] argv;
  delete [] arguments;

  return iResult;
}

int testAliHLTBlockDataCollection(int verbosity=0)
{
  gVerbosity=verbosity;
  AliHLTBlockDataCollection collection;
  for (int initNo=0; gInits[initNo].arguments!=NULL; initNo++) {
    if (gVerbosity>0) cout << "checking: " << gInits[initNo].arguments << endl;
    int result=InitCollection(collection, gInits[initNo].arguments);
    if (result!=gInits[initNo].result) {
      cerr << "failed: " << initNo << " (" << gInits[initNo].arguments << ") " << ": result " << result << ", expected " << gInits[initNo].result << endl;
      return -1;
    } else {
      if (gVerbosity>1) cout << "init " << initNo << " (" << gInits[initNo].arguments << ") " << "succeeded: result " << result << endl;
    }
  }

  AliHLTComponentBlockData bd;
  AliHLTComponent::FillBlockData(bd);

  for (int selectNo=0; gSelections[selectNo].id!=NULL; selectNo++) {
    AliHLTComponent::SetDataType(bd.fDataType, gSelections[selectNo].id, gSelections[selectNo].origin );
    bd.fSpecification=gSelections[selectNo].specification;
    int result=collection.IsSelected(bd);
    if (result!=gSelections[selectNo].result) {
      cerr << "failed: block " << gSelections[selectNo].id << ":" << gSelections[selectNo].origin << " 0x" << hex << gSelections[selectNo].specification << ": result " << result << ", expected " << gSelections[selectNo].result << endl;
      return -1;
    } else {
      if (gVerbosity>0) cout << "checking: block " << gSelections[selectNo].id << ":" << gSelections[selectNo].origin << " 0x" << hex << gSelections[selectNo].specification << (result==1?" selected":" not selected") << " (correctly)" << endl;
    }
  }

  return 0;
}

int main(int /*argc*/, const char** /*argv*/)
{
  int iResult=0;
  if ((iResult=testAliHLTBlockDataCollection(1))<0) {
    cerr << "<<<<< re-run with higher verbosity >>>>>>>" << endl;
    testAliHLTBlockDataCollection(2);
  }
  return iResult;
}
