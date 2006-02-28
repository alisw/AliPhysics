#ifndef AliEMCALHistoUtilities_H
#define AliEMCALHistoUtilities_H
/* Copyright(c) 1998-2002, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

/* $Id$ */

//_________________________________________________________________________
//                  
//*-- Authors: J.L. Klay (LLNL) & Aleksei Pavlinov (WSU)

#include <TNamed.h>
class TList;
class TBrowser;

class AliEMCALHistoUtilities: public TNamed {
  public:
  AliEMCALHistoUtilities(const char *name="emcalHistoUtilities",
  const char *tit="Histogram Utilities methods for EMCAL");
  AliEMCALHistoUtilities(const  AliEMCALHistoUtilities &) : TNamed("", ""){
    Fatal("cpy ctor", "not implemented") ; }
  virtual ~AliEMCALHistoUtilities();

  void    SetDebug(Int_t flag) {fDebug = flag;}
  Float_t GetDebug() const  {return fDebug;}
  virtual Bool_t  IsFolder() const;
  virtual void Browse(TBrowser* b) const ;

  // service routine
  static TList *MoveHistsToList(const char* name="ListOfHists", Bool_t putToBrowser=kTRUE);
  static void FillH1(TList *l=0, Int_t ind=0, Double_t x=-99999., Double_t w=1.);
  static void FillH2(TList *l=0, Int_t ind=0, Double_t x=-99999., Double_t y=-99999., Double_t w=1.);
  static int  SaveListOfHists(TList *list=0, const char* name="test", Bool_t kSingleKey=kFALSE,
  const char* opt="RECREATE");

  AliEMCALHistoUtilities & operator = (const AliEMCALHistoUtilities &) {
    Fatal("operator =", "not implemented") ; return *this ; }
  
  private:
  Int_t   fDebug;	// debug flag
  TList*  fListHist;    //!

  ClassDef(AliEMCALHistoUtilities,1) // EMCAL Histogram service routines
};

#endif // AliEMCALHistoUtilities_H
