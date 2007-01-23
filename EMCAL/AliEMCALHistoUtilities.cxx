/**************************************************************************
 * Copyright(c) 1998-2002, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

//_________________________________________________________________________
// This is a set of histogram
// utilities for the EMCAL
// to make some common
// functions easier
//
//*-- Authors: J.L. Klay (LLNL) & Aleksei Pavlinov (WSU) 

#include "AliEMCALHistoUtilities.h"

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TList.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TROOT.h>
#include <TRegexp.h>
#include <TString.h>

ClassImp(AliEMCALHistoUtilities)

AliEMCALHistoUtilities::AliEMCALHistoUtilities(const char *name, const char *tit) : TNamed(name,tit)
{
  // constructor
}

AliEMCALHistoUtilities::~AliEMCALHistoUtilities()
{
	//destructor
}  

TList* AliEMCALHistoUtilities::MoveHistsToList(const char* name, Bool_t putToBrowser)
{
  // Move HIST to list
  gROOT->cd();
  TIter nextHist(gDirectory->GetList());
  TList *list = new TList;
  list->SetName(name);
  TObject *objHist;
  while((objHist=nextHist())){
    if (!objHist->InheritsFrom("TH1")) continue;
    ((TH1*)objHist)->SetDirectory(0); // Remove from gROOT
    list->Add(objHist);
  }
  if(putToBrowser) gROOT->GetListOfBrowsables()->Add((TObject*)list);
  return list;
}

void AliEMCALHistoUtilities::FillH1(TList *l, Int_t ind, Double_t x, Double_t w)
{
  //fill 1d histogram
  static TH1* hid=0;
  if(l == 0) return;
  if(ind < l->GetSize()){
    hid = (TH1*)l->At(ind);
    hid->Fill(x,w);
  }
}

void AliEMCALHistoUtilities::FillH2(TList *l, Int_t ind, Double_t x, Double_t y, Double_t w)
{
  //fill 2d histogram
  static TH2* hid=0;
  if(l == 0) return;
  if(ind < l->GetSize()){
    hid = (TH2*)l->At(ind);
    hid->Fill(x,y,w);
  }
}

int AliEMCALHistoUtilities::SaveListOfHists(TList *mylist,const char* name,Bool_t kSingleKey,const char* opt)
{
  //write histograms to file
  printf(" Name of out file |%s|\n", name); 
  int save = 0;
  if(mylist && mylist->GetSize() && strlen(name)){
    TString nf(name); 
    if(nf.Contains(".root") == kFALSE) nf += ".root";
    TFile file(nf.Data(),opt);
    TIter nextHist(mylist);
    TObject* objHist=0;
    int nh=0;
    if(kSingleKey) {
       file.cd();
       mylist->Write(mylist->GetName(),TObject::kSingleKey);
       mylist->ls();
       save = 1;
    } else {
      while((objHist=nextHist())) { // loop over list 
        if(objHist->InheritsFrom("TH1")) {
          TH1* hid = (TH1*)objHist;
          file.cd();
          hid->Write();
          nh++;
          printf("Save hist. %s \n",hid ->GetName());
        }
      }
      printf("%i hists. save to file -> %s\n", nh,file.GetName());
      if(nh>0) save = 1;
    }
    file.Close();
  } else {
    printf("AliEMCALHistoUtilities::SaveListOfHists : N O  S A V I N G \n");
  }
  return save;
}

int AliEMCALHistoUtilities::ParseString(const TString &topt, TObjArray &Opt)
{ 
  // Moved from AliEMCALGeometry
  // Feb 06, 2006
  Ssiz_t begin, index, end, end2;
  begin = index = end = end2 = 0;
  TRegexp separator("[^ ;,\\t\\s/]+");
  while ( (begin < topt.Length()) && (index != kNPOS) ) {
    // loop over given options
    index = topt.Index(separator,&end,begin);
    if (index >= 0 && end >= 1) {
      TString substring(topt(index,end));
      Opt.Add(new TObjString(substring.Data()));
    }
    begin += end+1;
  }
  return Opt.GetEntries();
}
