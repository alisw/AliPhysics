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
  Int_t  fDebug; 

  protected:
  TFile*  fFile;
  TTree*  fTree;
  TString fName;
  TList*  fListHist;    //!

  // For partons after hard scattering
  Int_t   fNpart;
  Float_t fXpt[4]; 
  Float_t fXeta[4]; 
  Float_t fXphi[4];
  // Jet 
  Int_t   fNjet;
  Float_t fJet[10]; 
  Float_t fJetaL[10]; 
  Float_t fJphiL[10];
  Float_t fJetaW[10]; 
  Float_t fJphiW[10];
  // Charge particle in jet ??

  public:
  AliEMCALJetMicroDst(char *name="jetMicroDst",
  char *tit="jet Micro Dst for preparation of proposal");
  virtual ~AliEMCALJetMicroDst();
  virtual Bool_t  Create(TFile *file);
  virtual Bool_t  Create(const char  *fname);
  virtual void    Fill(AliRun *run, AliEMCALJetFinder* jetFinder);
  virtual void    FillPartons(AliGenHijingEventHeader *header);
  virtual void    FillPartons();
  virtual void    FillJets(AliEMCALJetFinder* jetFinder);

  virtual Bool_t  Open(const char  *fname);    // *MENU* 
  virtual Bool_t  Initialize(TFile *file);
  virtual void    Print(Option_t* option="") const;    // *MENU* 
  virtual Int_t   GetEntry(Int_t entry);
  virtual void    Test();
  Int_t   GetNpart() {return fNpart;}
  Bool_t  GetParton(Int_t i, Float_t& pt, Float_t& eta, Float_t& phi);
  Bool_t  GetParton(Int_t i, TVector3& vec);
  Int_t   GetNjet() {return fNjet;} 
  Bool_t  GetJet(Int_t i,Int_t mode, Float_t& pt,Float_t& eta,Float_t& phi);
  Bool_t  GetJet(Int_t i,Int_t mode, TVector3& vec);
  static  void FillVector(Float_t pt, Float_t eta, Float_t phi, TVector3& vec);

  void    SetDebug(Int_t flag) {fDebug = flag;}
  Float_t GetDebug() const  {return fDebug;}

  TTree* GetTree() {return fTree;}
  TFile* GetFile() {return fFile;}
  void   Close();

  virtual Bool_t  IsFolder() const;
  virtual void Browse(TBrowser* b);

  static TList *MoveHistsToList(char* name="List of Hist", Bool_t putToBrowser=kTRUE);

  ClassDef(AliEMCALJetMicroDst,1) // Micro Dst for jet analysis
};

#endif // AliEMCALJETMICRODST_H
