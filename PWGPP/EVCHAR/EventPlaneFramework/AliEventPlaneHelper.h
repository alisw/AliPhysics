/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/*
***********************************************************

 Event plane framework with helper functions wrapped in a namespace
 Contact: Jaap Onderwaater, j.onderwaater@gsi.de, jacobus.onderwaater@cern.ch
 
   Based on:
   AliDielectronHelper
   Dielectron helper functions wrapped in a namespace
   Authors: 
   Jens Wiechula <Jens.Wiechula@cern.ch> 
   Frederick Kramer <Frederick.Kramer@cern.ch> 
   Julian Book <Julian.Book@cern.ch>
***********************************************************
*/
#ifndef ALIEVENTPLANEHELPER_H
#define ALIEVENTPLANEHELPER_H


#include <TVectorDfwd.h>
#include <TDirectory.h>
#include <THashList.h>
#include <TObjArray.h>
#include <TChain.h>
#include <TArray.h>
#include <TFile.h>
#include <TKey.h>
#include <TAxis.h>
#include <THn.h>


class AliEventPlaneHelper : public TNamed{
public:

      
  AliEventPlaneHelper();
  virtual ~AliEventPlaneHelper();


  //TObjArray* gHistLists=0x0;   // main histogram list for the current running process

  static TVectorD* MakeLogBinning(Int_t nbinsX, Double_t xmin, Double_t xmax);
  static TVectorD* MakeLinBinning(Int_t nbinsX, Double_t xmin, Double_t xmax);
  static TVectorD* MakeArbitraryBinning(const char* bins);
  static Double_t* MakeBins(Int_t nbins, Double_t min, Double_t max);
  
 
  static void CloseFile() ;
  static void InitFile(const Char_t* filename) ;
  static TChain* GetChain(const Char_t* filename, Int_t howMany, Int_t offset, Long64_t& entries,
                                     const Char_t* friendpath, TChain* friendChain/*=0x0*/, const Char_t* friendChainFile/*=0x0*/) ;
  static TChain* GetChain(const Char_t* filename, Int_t howMany, Int_t offset, Long64_t& entries,
                                     TChain* friendChain=0x0, const Char_t* friendChainFile=0x0) ;
  static TArrayD * ArrayConversion(TClonesArray* array, Int_t nvectors);
  static TArrayD * AppendArray(TArrayD* array, TArrayD arr, Int_t length);
  static Double_t * GetElements(TClonesArray* array, Int_t vector);
  static TList* GetHistogramList(const Char_t* listname);    // get a histogram list
  static TObject* GetHistogram(const Char_t* listname, const Char_t* hname);  // get a histogram from an old output
  static THnF* AddHistogram( const Char_t* name, const Char_t* title,Int_t nDimensions,TArrayD* binLimits);
  static THnF* AddHistogram( const Char_t* name, const Char_t* title,Int_t nDimensions,TAxis* binLimits);

  static TDirectoryFile* fgHistCali;  // main directory for a standard tree analysis output (used for calibration, plotting etc.)
  static TFile* fgHistCaliFile;      // pointer to a TFile opened for reading

private:

  ClassDef(AliEventPlaneHelper,1)
};

#endif
