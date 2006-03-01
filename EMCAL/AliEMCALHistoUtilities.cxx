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

/*
$Log$
Revision 1.1  2006/02/28 21:55:11  jklay
add histogram utilities class, correct package definitions

*/

//*-- Authors: J.L. Klay (LLNL) & Aleksei Pavlinov (WSU) 

//*

#include <TBrowser.h>
#include <TFile.h>
#include <TList.h>
#include <TH1.h>
#include <TH2.h>
#include <TROOT.h>
#include <TString.h>

#include "AliEMCALHistoUtilities.h"

ClassImp(AliEMCALHistoUtilities)

AliEMCALHistoUtilities::AliEMCALHistoUtilities(const char *name, const char *tit) : TNamed(name,tit)
{
	//constructor
  fDebug = 0;
  gROOT->cd();
  fListHist = MoveHistsToList("Hist For AliEMCALHistoUtilities", kFALSE); 
}

AliEMCALHistoUtilities::~AliEMCALHistoUtilities()
{
	//destructor
}  

void AliEMCALHistoUtilities::Browse(TBrowser* b)
{
  // Browse
   if(fListHist)  b->Add((TObject*)fListHist);
   //   TObject::Browse(b);
}

Bool_t  AliEMCALHistoUtilities::IsFolder() const
{
  // Is folder
  if(fListHist) return kTRUE;
  else                   return kFALSE;
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
  static TH1* hid=0;
  if(l == 0) return;
  if(ind < l->GetSize()){
    hid = (TH1*)l->At(ind);
    hid->Fill(x,w);
  }
}

void AliEMCALHistoUtilities::FillH2(TList *l, Int_t ind, Double_t x, Double_t y, Double_t w)
{
  static TH2* hid=0;
  if(l == 0) return;
  if(ind < l->GetSize()){
    hid = (TH2*)l->At(ind);
    hid->Fill(x,y,w);
  }
}

int AliEMCALHistoUtilities::SaveListOfHists(TList *mylist,const char* name,Bool_t kSingleKey,const char* opt)
{
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
