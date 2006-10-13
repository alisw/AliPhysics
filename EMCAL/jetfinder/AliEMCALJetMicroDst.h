#ifndef ALIEMCALJETMICRODST_H
#define ALIEMCALJETMICRODST_H
/* Copyright(c) 1998-2002, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

/* $Id$ */

//_________________________________________________________________________
//                  
//*-- Author: Aleksei Pavlinov (WSU)
#include <TNamed.h>
class TFile;
class TTree;
class TH1;
class TH1F;
class TH2F;
class TH2;

class AliGenHijingEventHeader;
class AliRun;
class AliEMCALJetFinder;
class TVector3;
class TBrowser;

class AliEMCALJetMicroDst: public TNamed {
  public:
  AliEMCALJetMicroDst(const char *name="jetMicroDst",
  const char *tit="jet Micro Dst for preparation of proposal");
  AliEMCALJetMicroDst(const  AliEMCALJetMicroDst &);

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

  Bool_t  Open(Int_t mode=1) {return Open(DefineName(mode));}  // *MENU* 
  Bool_t  Open(const char  *fname);                                  // *MENU* 
  const Char_t* DefineName(Int_t mode=1);                      // *MENU*
  Bool_t  Initialize(TFile *file);
  void    Print(Option_t* option="") const;                          // *MENU* 
  Int_t   GetEntry(Int_t entry);
  void    Test();
  Int_t   GetNpart() const {return fnpart;}
  Bool_t  GetParton(Int_t i, Float_t& pt, Float_t& eta, Float_t& phi) const ;
  Bool_t  GetParton(Int_t i, TVector3& vec) const ;
  Int_t   GetNjet() const {return fnjet;} 
  Bool_t  GetJet(Int_t i,Int_t mode, Float_t& pt,Float_t& eta,Float_t& phi) const ;
  Bool_t  GetJet(Int_t i,Int_t mode, TVector3& vec) const ;
  static  void FillVector(Float_t pt, Float_t eta, Float_t phi, TVector3& vec);
  void    GetEtaPhi(Int_t id, Double_t &eta, Double_t &phi) const ;
  TVector3& GetCellVector(Int_t i) const ;
  TVector3& GetGridVector(Int_t i) const ;
  // 13-apr-2003
  Double_t GetSumInCone(TVector3 &jet, Int_t nc, Float_t *et,Float_t *eta,Float_t *phi, Double_t cellEtCut, Double_t rJet) const ;
  Double_t GetEmcalEtInCone(TVector3 &jet, Double_t cellEtCut=0.0, Double_t rJet=0.5) ;
  Double_t GetTpcPtInCone(TVector3 &jet, Double_t cellEtCut=0.0, Double_t rJet=0.5) ;
  Double_t GetSum(Int_t n, Float_t *ar, Double_t cut=0.0) const ;
  Double_t GetSumEmcal(Double_t cut=0.0) {return GetSum(fncell, fetcell, cut);}
  Double_t GetSumTpc(Double_t cut=0.0) {return GetSum(fnchp, fppt, cut);}

  void    SetDebug(Int_t flag) {fDebug = flag;}
  Float_t GetDebug() const  {return fDebug;}

  TTree* GetTree() const {return fTree;}
  TFile* GetFile() const {return fFile;}
  void   Close();

  Bool_t  IsPythiaDst() const ;
  virtual Bool_t  IsFolder() const;
  virtual void Browse(TBrowser* b) ;

  // service routine
  static TList *MoveHistsToList(const char* name="ListOfHists", Bool_t putToBrowser=kTRUE);
  static void FillH1(TList *l=0, Int_t ind=0, Double_t x=-99999., Double_t w=1.);
  static void FillH2(TList *l=0, Int_t ind=0, Double_t x=-99999., Double_t y=-99999., Double_t w=1.);
  static int  SaveListOfHists(TList *list=0, const char* name="test", Bool_t kSingleKey=kFALSE,
  const char* opt="RECREATE");
  static void Sgpdge(Int_t pdgId, Int_t &gId); 

  AliEMCALJetMicroDst & operator = (const AliEMCALJetMicroDst &) {
    Fatal("operator =", "not implemented") ; return *this ; }
  
  private:
  
  Float_t fpphi[20000]; //[nchp]
  Int_t   fDebug;	// debug flag
  TFile*  fFile;	// filename
  TTree*  fTree;	// Tree pointer
  TList*  fListHist;    //!
  TString fFileName;    // for convenience

  Float_t fdecone;   //! for EMCAL
  Float_t fptcone;   //! for ch.particles 
  // For partons after hard scattering
  Int_t   fnpart;	//npartons
  Float_t fxpt[4];  //[npart]
  Float_t fxeta[4]; //[npart]
  Float_t fxphi[4]; //[npart]
  // Jet 
  Int_t   fnjet;	// number of jets
  Float_t fjet[10];   //[njet]
  Float_t fjetal[10]; //[njet]
  Float_t fjphil[10]; //[njet]
  Float_t fjetaw[10]; //[njet]
  Float_t fjphiw[10]; //[njet]
  // Charge particle in jet ??
  // eT in EMCAL itself - 24-jan-2003
  Int_t   fncell;         // 96*144 =13824 
  Int_t   fidcell[13824]; //[ncell]
  Float_t fetcell[13824]; //[ncell] : de = det*sf
  // eT in EMCAL grid for jet finder 
  Int_t   fngrid;         // 96*144 =13824 
  Int_t   fidgrid[13824]; //[ngrid]
  Float_t fetgrid[13824]; //[ngrid]
  // charge particle which hit to EMCAL - 28-jan-2003
  Int_t   fnchp;	//number of charged particles
  Int_t   fpid[20000];  //[nchp]
  Float_t fppt[20000];  //[nchp]
  Float_t fpeta[20000]; //[nchp]
  TH1F*  fhPtPart; //hist
  TH1F*	 fhNJet;	//hist
  TH1F*  fhPtJet;	//hist
  TH2F*  fhEtaPhiPart;	//hist
  TH2F*  fhEtaPhiJet;	//hist
  TH1F*  fhNcell;	//hist
  TH1F*	 fhCellId;	//hist
  TH1F*	 fhCellEt;	//hist
  	
  TH1F*	 fhSumEt;	//hist
  TH1F*  fhNgrid;	//hist
  TH1F*	 fhGridId;	//hist
  TH1F*  fhGridEt;	//hist
  TH1F*  fhSumEtGrForJF;	//hist


  ClassDef(AliEMCALJetMicroDst,3) // Micro Dst for jet analysis
};

#endif // AliEMCALJETMICRODST_H

typedef AliEMCALJetMicroDst sv; // for convinience

/*
What to do
1. Common info about event
 */
