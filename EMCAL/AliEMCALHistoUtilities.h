#ifndef AliEMCALHistoUtilities_H
#define AliEMCALHistoUtilities_H
/* Copyright(c) 1998-2002, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

/* $Id$ */

//_________________________________________________________________________
// This is just set of static methods for common using
//                  
//*-- Authors: J.L. Klay (LLNL) & Aleksei Pavlinov (WSU)

#include <TNamed.h>
class TList;

class AliEMCALHistoUtilities: public TNamed {
  public:
  AliEMCALHistoUtilities(const char *name="emcalUtilitiesRoutines",
  const char *tit="EMCAL utility routines");
  AliEMCALHistoUtilities(const  AliEMCALHistoUtilities &) : TNamed("", ""){
    Fatal("cpy ctor", "not implemented") ; }
  virtual ~AliEMCALHistoUtilities();

  // service routine
  static TList *MoveHistsToList(const char* name="ListOfHists", Bool_t putToBrowser=kTRUE);
  static void FillH1(TList *l=0, Int_t ind=0, Double_t x=-99999., Double_t w=1.);
  static void FillH2(TList *l=0, Int_t ind=0, Double_t x=-99999., Double_t y=-99999., Double_t w=1.);
  static int  SaveListOfHists(TList *mylist=0, const char* name="test", Bool_t kSingleKey=kFALSE,
  const char* opt="RECREATE");
  // 
  static int ParseString(const TString &topt, TObjArray &Opt); 

  AliEMCALHistoUtilities & operator = (const AliEMCALHistoUtilities &) {
    Fatal("operator =", "not implemented") ; return *this ; }
  
  ClassDef(AliEMCALHistoUtilities,1) // EMCAL utility routines
};

#endif // AliEMCALHistoUtilities_H
