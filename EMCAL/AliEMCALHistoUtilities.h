#ifndef ALIEMCALHISTOUTILITIES_H
#define ALIEMCALHISTOUTILITIES_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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
class TString;
class TH1;
class TGraph;
class TGraphErrors;
class TF1;
class TLatex;
class TChain;
//class TLorentzVector;
class TArrayF;

//class AliESDCaloCluster;
//class AliEMCALRecPoint;
//class AliRunLoader;

class AliEMCALHistoUtilities: public TNamed {
  public:  AliEMCALHistoUtilities(const char *name="emcalUtilitiesRoutines",
  const char *tit="EMCAL utility routines");
  AliEMCALHistoUtilities(const  AliEMCALHistoUtilities &) : TNamed("", ""){
    Fatal("cpy ctor", "not implemented") ; }
  virtual ~AliEMCALHistoUtilities();

  // service routine
  static TList *MoveHistsToList(const char* name="ListOfHists", Bool_t putToBrowser=kTRUE);
  static void FillH1(TList *l=0, Int_t ind=0, Double_t x=-99999., Double_t w=1., Double_t error=0.);
  static void FillH2(TList *l=0, Int_t ind=0, Double_t x=-99999., Double_t y=-99999., Double_t w=1.);
  static int  SaveListOfHists(TList *mylist=0, const char* name="test", Bool_t kSingleKey=kFALSE,
  const char* opt="RECREATE");
  static void AddToNameAndTitle(TH1   *h=0, const char *name=0, const char *title=0);
  static void AddToNameAndTitleToList(TList *l=0, const char *name=0, const char *title=0);
  static void ResetListOfHists(TList *l);

  static TLatex *Lat(const char *text="", Float_t x=0.0,Float_t y=0.0, Int_t align=12, Float_t tsize=0.05, short tcolor = 1); 
  static TGraph *DrawGraph(Int_t n=4, Double_t *x=0, Double_t *y=0, Int_t markerColor=4,  
  Int_t markerStyle=4, const char* opt="", const char* tit="", const char* xTit="  jet E_{t}  ",
  const char* yTit="", Int_t ifun=0, const char *optFit="W+", const char *fun="");
  static TGraphErrors *DrawGraphErrors(const Int_t n=4,Double_t *x=0,Double_t *y=0,Double_t *ex=0, 
  Double_t *ey=0, Int_t markerColor=4,Int_t markerStyle=4, const char* opt="", 
  const char* tit="", const char* xTit="  jet E_{t}  ",
  const char* yTit="", Int_t ifun=0, const char *optFit="W+", const char *fun="");
  // TChain
  static void InitChain(TChain *chain=0, const char* nameListOfFiles=0, Int_t nFileMax=0); 
  //static AliRunLoader* InitKinematics(const Int_t nev=0, const char* galiceName="galice.root");
  //static AliRunLoader* GetRunLoader(const Int_t nev, const Char_t* galiceName,
	//			 const Char_t* eventFolderName, AliRunLoader* rlOld);
  //
  static Double_t GetMomentum(const char* nameListOfFiles); 
  static int ParseString(const TString &topt, TObjArray &Opt); 
  // Analysis utilites
  //static Bool_t GetLorentzVectorFromESDCluster(TLorentzVector &v, const AliESDCaloCluster *cl);
  //static Bool_t GetLorentzVectorFromRecPoint(TLorentzVector &v, const AliEMCALRecPoint  *rp);
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
  // Calibration stuff
  static Double_t GetCorrectionCoefficientForGamma1(const Double_t eRec);
  static Double_t GetCorrectedEnergyForGamma1(const Double_t eRec);
  static TF1* GetResolutionFunction(const char *opt, TString &latexName);
  //
  // Analysis
  //
  // Trigger 
  static TList* GetTriggersListOfHists(const Int_t scale=0, const Int_t nTrig=5, const Bool_t toBrowser=kFALSE);
  static void   FillTriggersListOfHists(TList *l=0, TArrayF *triggerPosition=0, TArrayF *triggerAmplitudes=0);
  // Jet(s) kinematics
  static TList* GetJetsListOfHists(Int_t njet=2, Bool_t toBrowser=kFALSE);
  //static void   FillJetKineListOfHists(TList *l, AliRunLoader* rl, TLorentzVector &goodJet);

  AliEMCALHistoUtilities & operator = (const AliEMCALHistoUtilities &) {
    Fatal("operator =", "not implemented") ; return *this ; }
  
  ClassDef(AliEMCALHistoUtilities,1) // EMCAL utility routines
};

#endif // ALIEMCALHISTOUTILITIES_H
