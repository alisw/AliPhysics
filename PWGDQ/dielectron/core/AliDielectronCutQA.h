#ifndef ALIDIELECTRONCUTQA_H
#define ALIDIELECTRONCUTQA_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//#################################################################
//#                                                               #
//#             Class AliDielectronCutQA                          #
//#              Dielectron Group of cuts                         #
//#                                                               #
//#  Authors:                                                     #
//#   Julian    Book,     Uni Ffm / Julian.Book@cern.ch           #
//#                                                               #
//#################################################################

#include <TNamed.h>
#include <AliAnalysisFilter.h>
#include <TH1F.h>
#include <TList.h>
#include <TObjArray.h>

class TCollection;

class AliDielectronCutQA : public TNamed {
  
public:
  enum { kEvent=0, kTrack, kPair, kNtypes };

  AliDielectronCutQA();
  AliDielectronCutQA(const char*name, const char* title);
  
  virtual ~AliDielectronCutQA();

  //Analysis cuts interface
  //
  void Init();
  void AddTrackFilter(     AliAnalysisFilter *trkFilter);
  /* void AddPrePairFilter(   AliAnalysisFilter *prePairFilter); */
  /* void AddPrePairLegFilter(AliAnalysisFilter *prePairLegFilter); */
  void AddPairFilter(      AliAnalysisFilter *pairFilter);
  void AddEventFilter(     AliAnalysisFilter *eventFilter);

  //  void Fill(AliAnalysisCuts *cut);
  void Fill(UInt_t mask, TObject *obj);
  void FillAll(TObject *obj);// { fCutQA->Fill(0); }

  const TObjArray * GetQAHistArray() const { return &fQAHistArray; }


private:

  TObjArray fQAHistArray;              //-> array of QA histograms
  TH1F *fCutQA[kNtypes];               // qa histogram
  Int_t fNCuts[kNtypes];               // number of cuts
  const char* fCutNames[20][kNtypes];  // cut names
  const char* fTypeKeys[kNtypes];      // type names


  UInt_t GetObjIndex(TObject *obj);    // return object index

  AliDielectronCutQA(const AliDielectronCutQA &);
  AliDielectronCutQA &operator=(const AliDielectronCutQA &);
  
  ClassDef(AliDielectronCutQA,1) //Group of cuts
};

#endif
