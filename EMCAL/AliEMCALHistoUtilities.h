#ifndef AliEMCALHistoUtilities_H
#define AliEMCALHistoUtilities_H
/* Copyright(c) 1998-2002, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

/* $Id$ */

//_________________________________________________________________________
// This is a set of histogram
// utilities for the EMCAL
// to make some common
// functions easier
//                  
//*-- Authors: J.L. Klay (LLNL) & Aleksei Pavlinov (WSU)

#include <TNamed.h>

class TList;
class TH1;
class TF1;
class TLatex;
class TChain;
class TLorentzVector;

class AliESDCaloCluster;
class AliEMCALRecPoint;

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
  static void AddToNameAndTitle(TH1   *h=0, const char *name=0, const char *title=0);
  static void AddToNameAndTitleToList(TList *l=0, const char *name=0, const char *title=0);
  static TLatex *lat(const char *text="", Float_t x=0.0,Float_t y=0.0, Int_t align=12, Float_t tsize=0.05, short tcolor = 1); 
  // TChain
  static void InitChain(TChain *chain=0, const char* nameListOfFiles=0, Int_t nFileMax=0); 
  // 
  static int ParseString(const TString &topt, TObjArray &Opt); 
  // Analysis utilites
  static Bool_t GetLorentzVectorFromESDCluster(TLorentzVector &v, const AliESDCaloCluster *cl);
  static Bool_t GetLorentzVectorFromRecPoint(TLorentzVector &v, const AliEMCALRecPoint  *rp);
  // Drawing 
  static void DrawHist(TH1* hid=0,int lineWidth=1,int lineColor=1,const char* opt="",int lineStyle=1);
  // Fitting:
  static TF1* Gausi(const char *addName, double xmi, double xma, double N, double mean, double sig, 
  double width);
  static TF1* Gausi(const char *addName, double xmi, double xma, TH1 *h);

  static TF1* GausiPol2(const char *addName, double xmi, double xma, TF1 *g, TF1* bg);
  //
  static Double_t Gi(Double_t *x, Double_t *par);
  static Double_t GiPol2(Double_t *x, Double_t *par);

  AliEMCALHistoUtilities & operator = (const AliEMCALHistoUtilities &) {
    Fatal("operator =", "not implemented") ; return *this ; }
  
  ClassDef(AliEMCALHistoUtilities,1) // EMCAL utility routines
};

#endif // AliEMCALHistoUtilities_H
