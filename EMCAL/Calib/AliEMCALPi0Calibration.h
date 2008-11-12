#ifndef ALIEMCALPI0CALIBRATION_H
#define ALIEMCALPI0CALIBRATION_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

/* $Id: AliEMCALPi0Calibration.h 21465 2007-10-16 14:38:13Z pavlinov $ */
 
//*--  Authors: Aleksei Pavlinov (WSU)
//  Pi0 calibration
//  Tuning parameters of coordinate calculations 
//  May 2008 - move to AliAnalysisTaskSE

#include "AliAnalysisTaskSE.h"

#include <TObjArray.h>

class AliEMCALGeometry;
class AliEMCALFolder;
class AliRunLoader;
class AliEMCALRecPoint;
class AliEMCALCellInfo;
class AliEMCALCellIndexes;

class TList;
class TCanvas;
class TH1F;
class TH2F;
class TBrowser;
class TChain;
class TFolder;
class TArrayI;

class AliEMCALPi0Calibration :  public AliAnalysisTaskSE {
  public:
    AliEMCALPi0Calibration();
    AliEMCALPi0Calibration(const char* name);
    //AliEMCALPi0Calibration(const AliAnalysisTaskSE& obj);
    //AliEMCALPi0Calibration& operator=(const AliAnalysisTaskSE& other);

    virtual ~AliEMCALPi0Calibration();

    virtual void UserExec(const Option_t* opt);
  //
  void InitStructure(Int_t it);
  static TList *DefineHistsOfRP(const char *name="RP",  Double_t p=110.0, Int_t keyOpt=0);
  static TList *DefineHistsOfKineVsRP(const char *name="KineVsRP",  Double_t p=110.0, Int_t keyOpt=0);
  static TList *DefineHistsForShowerProfile(const char *name="ProfY", Double_t p=1.);
  static void   FillHistsOfKineVsRP(TList *l, AliRunLoader* RL, TClonesArray &lvM);
  static void   FillHistsForShowerProfile(TList *l, AliEMCALRecPoint *rp, AliEMCALCellInfo* t);
  static void   EvalLocalPhiPosition(const Double_t wlog, const AliEMCALRecPoint *rp, const AliEMCALCellInfo* t, Double_t &xcog, Int_t &phiSize, AliEMCALCellIndexes &rMax); 
  //
  TList *GetListKineVsRP()      {return fLKineVsRP;}
  TList *GetListShowerProfile() {return fLShowerProfile;}
  //  static TList *DefineHistsOfTowers(const char *name="towers");
  // 
  void FitEffMassHist(); // *MENU*  
  void PrintInfo();      // *MENU*
  //
  void    SetChain(TChain *chain)  {fChain = chain;}
  TChain* GetChain()               {return fChain;}
  void    SetMomentum(Double_t p);
  Double_t GetMomentum() const {return fPmom;}

  AliEMCALFolder* CreateEmcalFolder(const Int_t it);
  void SetEmcalFolder(AliEMCALFolder* folder); 
  void SetEmcalOldFolder(AliEMCALFolder* folder); 
  AliEMCALFolder* GetEmcalOldFolder(const Int_t nsm);
  AliEMCALCellInfo* GetCellsInfo() {return fCellsInfo;}
  //
  void      SetStringOfRunOpts(const char *st) {fRunOpts = st;}
  TObjArray GetOptsArray() const {return fArrOpts;}
  Int_t     GetKeyOptsValue(Int_t key);  // *MENU*
  void      CheckRunOpts();
  //
  virtual void Browse(TBrowser* b);
  virtual Bool_t  IsFolder() const;
  //
  void    Save(Int_t ver=0, const char *optIO="NEW");   // *MENU*
  static  AliEMCALPi0Calibration* ReadSelector(const char* nf = "/home/pavlinov/ALICE/SHISHKEBAB/RF/CALIB/PROF/PROFILE_0.root");
  

  //
  //// Pictures staf - Jun 26, 2007
  //
  void ReadAllEmcalFolders();
  void PictVsIterNumber(const Int_t ind=0, const Int_t nsm=0);
  // Gamma
  TH1F* FitHistOfRecPointEnergy(const char *opt="CLONE");
  static TCanvas *Linearity(TList *l, Int_t ifun=3);
  static TCanvas *DrawKineVsRP(TList *l);
  // Profile
  static TCanvas *DrawMeanVsLog(TH2F *h2);
  // Geometry staff
  TCanvas *DrawPhiEtaAnglesDistribution(const char *gn="SHISH_TRD1_CURRENT_2X2"); // *MENU*
  // Geometry constants 
  TCanvas *DrawDeffVsEnergy();  // *MENU*
  TCanvas *DrawDeffVsEnergy2(const char *opt="fit1"); // *MENU*
  void     ReadParsDeffAndW0(const char *dirName="/data/r22b/ALICE/CALIB/FIT/",
	   double *deff=0, double *edeff=0, double *w0=0, double *ew0=0, const Int_t pri=0);
  TCanvas *DrawSpaceResolution();
  // 
  static AliEMCALFolder* GetEmcalFolder() {return fgEMCAL;}
  static AliEMCALFolder* GetEmcalOldFolder() {return fgEMCALOld;}
  static void SetFitParameters(Double_t deff, Double_t w0, Double_t slope) 
  {
    fgDistEff = deff; fgW0 = w0; fgSlopePhiShift = slope;
  }
  static void GetFitParameters(Double_t &deff, Double_t &w0, Double_t &slope)
  {
    deff = fgDistEff; w0 = fgW0; slope = fgSlopePhiShift;
  } 
  void ResetAllListOfHists();
  void ReloadChain(Long64_t entry=0);
  void GetInitialParsForFit(const Int_t var, Double_t &deff, Double_t &w0, Double_t &phislope, const int phiCase=0);

 protected:
  static AliEMCALFolder*  fgEMCAL;      // current  EMCAL object
  static AliEMCALFolder*  fgEMCALOld;   // previous EMCAL object
  //
  static Double_t fgDistEff;  // effective depth of electromagnetic shower
  static Double_t fgW0;       // parameter of log. methods 
  static Double_t fgSlopePhiShift; // phi shift of cluster = fSlopePhiShift * phi

  Double_t fPmom; // positive if defined
  //
  TChain* fChain; //! chain if ESD files
  TList* fLofHistsPC; // list of histograms of pseudo clusters 
  TList* fLofHistsRP; // list of histograms of rec.points 
  TList* fLKineVsRP;  // list of histograms kinematics vs rec.points 
  TList* fLShowerProfile;  // list of histograms for shower profile business
  //
  AliEMCALCellInfo *fCellsInfo; // pointer to current cell
  TFolder*         fEmcalPool;  // folder of EMCAL objects
  //
  // Options - Jul 10, 2007
  //
  TString   fRunOpts;        // String of running options
  TObjArray fArrOpts;        // Array of options 
  // Options keys
  TArrayI  *fKeyOpts;        // optins key; 0-disable, 1-enable
  // Static parameters
 private:
  AliEMCALPi0Calibration(const AliEMCALPi0Calibration&);
  AliEMCALPi0Calibration& operator=(const AliEMCALPi0Calibration&);
  //
  static AliEMCALGeometry* fgEmcalGeo; //! pointer to EMCAL geometry
  static Int_t fgNmaxCell;  //! max number of cells
  static const Char_t **fgAnaOpt; //! aray of options
  static Int_t fgNanaOpt;   //! number of options

  ClassDef(AliEMCALPi0Calibration, 1);
};
#endif
