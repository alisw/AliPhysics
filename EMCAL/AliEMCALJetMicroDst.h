#ifndef ALIEMCALJETMICRODST_H
#define ALIEMCALJETMICRODST_H
/* Copyright(c) 1998-2002, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

/* $Id$ */

//_________________________________________________________________________
//                  
//*-- Author: Aleksei Pavlinov (WSU)
#include <TNamed.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>

class AliGenHijingEventHeader;
class AliRun;
class AliEMCALJetFinder;
class TVector3;
class TBrowser;

class AliEMCALJetMicroDst: public TNamed {

  private:
  Int_t   fDebug;

  public:
  TFile*  fFile;
  TTree*  fTree;
  TString fName;
  TList*  fListHist;    //!
  TString fFileName;    // for convenience

  Float_t decone;   //! for EMCAL
  Float_t ptcone;   //! for ch.particles 
  // For partons after hard scattering
  Int_t   npart;
  Float_t xpt[4];  //[npart]
  Float_t xeta[4]; //[npart]
  Float_t xphi[4]; //[npart]
  // Jet 
  Int_t   njet;
  Float_t jet[10];   //[njet]
  Float_t jetal[10]; //[njet]
  Float_t jphil[10]; //[njet]
  Float_t jetaw[10]; //[njet]
  Float_t jphiw[10]; //[njet]
  // Charge particle in jet ??
  // eT in EMCAL itself - 24-jan-2003
  Int_t   ncell;         // 96*144 =13824 
  Int_t   idcell[13824]; //[ncell]
  Float_t etcell[13824]; //[ncell] : de = det*sf
  // eT in EMCAL grid for jet finder 
  Int_t   ngrid;         // 96*144 =13824 
  Int_t   idgrid[13824]; //[ngrid]
  Float_t etgrid[13824]; //[ngrid]
  // charge particle which hit to EMCAL - 28-jan-2003
  Int_t   nchp;
  Int_t   pid[20000];  //[nchp]
  Float_t ppt[20000];  //[nchp]
  Float_t peta[20000]; //[nchp]
  Float_t pphi[20000]; //[nchp]

  public:
  AliEMCALJetMicroDst(char *name="jetMicroDst",
  char *tit="jet Micro Dst for preparation of proposal");
  virtual ~AliEMCALJetMicroDst();
  Bool_t  Create(TFile *file);
  Bool_t  Create(const char  *fname);
  void    Fill(AliRun *run=0, AliEMCALJetFinder* jetFinder=0, Int_t modeFilling=0);
  void    FillPartons(AliGenHijingEventHeader *header);
  void    FillPartons();
  void    FillJets(AliEMCALJetFinder* jetFinder);
  void    FillEtForEMCAL(AliEMCALJetFinder* jetFinder);
  void    FillEtForGrid(AliEMCALJetFinder* jetFinder);
  void    FillArrays(TH2* hid, Int_t &n, Int_t *id, Float_t *et);
  void    FillChargeParticles(AliEMCALJetFinder* jetFinder);  
  
  void    FillJetsControl(); // 18-jan-2003

  Bool_t  Open(const Int_t mode=1) {return Open(DefineName(mode));}  // *MENU* 
  Bool_t  Open(const char  *fname);                                  // *MENU* 
  const Char_t* DefineName(const Int_t mode=1);                      // *MENU*
  Bool_t  Initialize(TFile *file);
  void    Print(Option_t* option="") const;                          // *MENU* 
  Int_t   GetEntry(Int_t entry);
  void    Test();
  Int_t   GetNpart() {return npart;}
  Bool_t  GetParton(Int_t i, Float_t& pt, Float_t& eta, Float_t& phi);
  Bool_t  GetParton(Int_t i, TVector3& vec);
  Int_t   GetNjet() {return njet;} 
  Bool_t  GetJet(Int_t i,Int_t mode, Float_t& pt,Float_t& eta,Float_t& phi);
  Bool_t  GetJet(Int_t i,Int_t mode, TVector3& vec);
  static  void FillVector(Float_t pt, Float_t eta, Float_t phi, TVector3& vec);
  void    GetEtaPhi(Int_t id, Double_t &eta, Double_t &phi);
  TVector3& GetCellVector(Int_t i);
  TVector3& GetGridVector(Int_t i);
  // 13-apr-2003
  Double_t GetSumInCone(TVector3 &jet, Int_t nc, Float_t *et,Float_t *eta,Float_t *phi, Double_t cellEtCut, Double_t rJet);
  Double_t GetEmcalEtInCone(TVector3 &jet, Double_t cellEtCut=0.0, Double_t rJet=0.5);
  Double_t GetTpcPtInCone(TVector3 &jet, Double_t cellEtCut=0.0, Double_t rJet=0.5);
  Double_t GetSum(Int_t n, Float_t *ar, Double_t cut=0.0);
  Double_t GetSumEmcal(Double_t cut=0.0) {return GetSum(ncell, etcell, cut);}
  Double_t GetSumTpc(Double_t cut=0.0) {return GetSum(nchp, ppt, cut);}

  void    SetDebug(Int_t flag) {fDebug = flag;}
  Float_t GetDebug() const  {return fDebug;}

  TTree* GetTree() {return fTree;}
  TFile* GetFile() {return fFile;}
  void   Close();

  Bool_t  IsPythiaDst();
  virtual Bool_t  IsFolder() const;
  virtual void Browse(TBrowser* b);

  static TList *MoveHistsToList(char* name="List of Hist", Bool_t putToBrowser=kTRUE);

  ClassDef(AliEMCALJetMicroDst,1) // Micro Dst for jet analysis
};

#endif // AliEMCALJETMICRODST_H
/*
What to do
1. Common info about event
 */
