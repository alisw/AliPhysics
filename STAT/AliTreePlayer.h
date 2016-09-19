#ifndef ALITREEPLAYER_H
#define ALITREEPLAYER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Alice extension of TTreePlayer            //
////////////////////////////////////////////////

 
#include "AliTreePlayer.h"

class AliTreePlayer : public TNamed {

public:
  AliTreePlayer() {}
  AliTreePlayer(const char *name, const char *title);
  static TObjArray  * selectMetadata(TTree * tree, TString query, Int_t verbose);
  static TObjArray  * selectTreeInfo(TTree* tree, TString query,Int_t verbose); 
  static void selectWhatWhereOrderBy(TTree * tree, TString what, TString where, TString orderBy,  Int_t firstentry, Int_t nentries, TString outputFormat, TString outputName);
  static TString  printSelectedTreeInfo(TTree*tree, TString infoType,  TString regExpFriend, TString regExpTag, Int_t verbose);
  // to add
  // Metadata query from the TStatToolkit (GetMetadata, AddMetadata)
  // TH1* DrawSorted(const char * expression, const char weights, Bool_t down, Int_t cut);  // draw sorted version of histogram 
  // Friend tree queries for time series using (time in entry or in array)
  //       - calculate nearest entry
  //       -- cache adding  branch   
  //       - calculate (index, mean, median, rms ...) entry in interval
  //       -- cache adding branch
  // To solve/consider: 
  //       - index branches have to have the same names - indeces to be created in all trees
  //       - can we change indeces? is it needed ?
  //       - implemntations for indeces er entry and internal indices within array will be different
  //       -- are both implementation needed - in which order ... 
  
  //static TH1* DrawSorted(TTree * tree, const char* varexp, const TCut& selection, Bool_t nbins=20, Long64_t nentries = 1000000000, Long64_t firstentry = 0);
protected:
  
  ClassDef(AliTreePlayer,1)  // Extension of the TTreePlayer
};

#endif

