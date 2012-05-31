// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/** 
 * @file   dumpTPCDigits.C
 * @author Matthias.Richter@ift.uib.no
 * @date   
 * @brief  Convert TPC digit info to ascii and print to stdout.
 *
 * The macro reads the simulated TPC digits from the RunLoader and
 * prints the digit info to stout.
 *
 * aliroot -b -q dumpTPCDigits.C | tee digits.log
 *
 * The macro asumes the data to be already simulated.
 */

#ifndef __CINT__
#include "TSystem.h"
#include "TTree.h"
#include "AliRunLoader.h"
#include "AliSimDigits.h"
#endif //__CINT__

int dumpDigits(AliRunLoader* runloader, int event);

int dumpTPCDigits()
{
  const char* galice_file="galice.root";
  TString param=galice_file;
  param+="?filetype=raw";
  TFile file(param);
  if (file.IsZombie()) {
    cerr << "can not open file " << galice_name << ", skipping test" << endl;
    return 0;
  }

  AliRunLoader* rl=AliRunLoader::Open(galice_file);
  if (!rl) {
    cerr << "can not open RunLoader" << endl;
    return -1;
  }

  gSystem->Load("libAliHLTTPC.so");

  dumpDigits(rl, 0);

  return 0;
}

int dumpDigits(AliRunLoader* runloader, int event)
{
  int iResult=0;
  AliLoader* tpcLoader = runloader->GetLoader("TPCLoader");
  if(!tpcLoader){
    cerr << "error: can not get TPC loader" << endl;
    return -1;
  }

  runloader->GetEvent(event);
  tpcLoader->LoadDigits();

  TTree *digitsTree = tpcLoader->TreeD();
  if(!digitsTree) {
    cerr << "error: can not get digits tree for event no " << event << endl;
    return -1;
  }

  AliSimDigits *digits;
  digitsTree->GetBranch("Segment")->SetAddress(&digits);

  AliTPCParam *tpcParam;
  runloader->CdGAFile();
  tpcParam = (AliTPCParam*)gFile->Get("75x40_100x60_150x60");
  if(!tpcParam){
    cerr << "error: can not retrieve parameters" << endl;
  }

  Int_t iPrintedSlice=-1;
  Int_t iPrintedPart=-1;

  for(Int_t n=0; n<digitsTree->GetEntries(); n++) {
    Int_t sector, row; // coordinates of the simualted data
    Int_t lslice,lrow; // local coordinates
    digitsTree->GetEvent(n);
    tpcParam->AdjustSectorRow(digits->GetID(),sector,row);
    if(!AliHLTTPCTransform::Sector2Slice(lslice,lrow,sector,row)){
      cerr << "error: conversion of coordinates sector/row " << sector << "/" << row << "failed" << endl;
      return -1;
    }

    Int_t part=AliHLTTPCTransform::GetPatch(lrow);

    if (iPrintedSlice!=lslice || iPrintedPart!=part) {
      iPrintedSlice=lslice;
      iPrintedPart=part;
      cout << "====================================================================" << endl;
      cout << "    Slice: " << iPrintedSlice << "   Partition: " << iPrintedPart << "   digit tree entry " << n << endl;
    }

    cout << "--------------------------------------------------------------------" << endl;
    cout << "Row: " << lrow << endl;

    digits->First();
    Int_t lasttime=-1;
    Int_t lastpad=-1;
    do {
      Int_t time=digits->CurrentRow();
      Int_t pad=digits->CurrentColumn();
      Int_t dig = digits->GetDigit(time,pad);
      if (lastpad!=pad) {
	if (lastpad!=-1) cout << "        -> " << lasttime << endl;
	lastpad=pad;
	cout << "Row: " << lrow << "  Pad: " << lastpad << endl;
	lasttime=-1;
      }
      if (lasttime!=time+1 && lasttime!=time-1 ) {
	if (lasttime!=-1) cout << "        -> " << lasttime << endl;
	cout << "                     Time " << time << ":  ";
      }
      lasttime=time;
      cout << "  " << dig;
    } while (digits->Next());
    if (lasttime) cout << "        -> " << lasttime << endl;
  }
}
