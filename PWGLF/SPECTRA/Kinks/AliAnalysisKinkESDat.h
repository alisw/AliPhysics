#ifndef ALIANALYSISKINKESDat_H
#define ALIANALYSISKINKESDat_H

/*  See cxx source for full Copyright notice */

//-----------------------------------------------------------------
//                 AliAnalysisKinkESDat class
//         This task is an example of an analysis task
//                  for kink topology Study
//          Authors: Martha Spyropoulou-Stassinaki
//           and members of the Greek group at the
//          Physics Department of Athens University
//                    mspyrop@phys.uoa.gr
//-----------------------------------------------------------------
class AliPIDResponse;
class AliESDVertex;
class AliESDtrack;
class TF1;
class TH1F;
class TH2F;
class TH3F;
class TH1D;
class TH2D;
class TList;
class AliESDtrackCuts;

#include "AliAnalysisTaskSE.h"

class AliAnalysisKinkESDat : public AliAnalysisTaskSE {
 public:
  AliAnalysisKinkESDat(const char *name = "AliAnalysisKinkESDat");
  virtual ~AliAnalysisKinkESDat() {}

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
  const AliESDVertex *GetEventVertex(const AliESDEvent* esd) const;
//           Marek multiplicity bins
  void SetMulCut(Int_t low, Int_t up){fLowMulcut=low;fUpMulcut=up;}	  
  void SetKinkRadius(Float_t lRadiusKLow, Float_t lRadiusKUp)  { fKinkRadLow=lRadiusKLow; fKinkRadUp=lRadiusKUp;}	  

  void SetNClusterCut(Int_t lowCluster){fLowCluster=lowCluster;}	  
  void SetQtCut(Float_t   lowQt){fLowQt=lowQt;}	  
 private:
   TH1F        *fHistPtESD; //Pt spectrum of all ESD inside eta, Pt cuts
   TH1F        *fHistPt; //Pt spectrum of all ESD tracks
   TH1F        *fHistQtAll; //Qt spectrum of all kinks
   TH1F        *fHistQt1; //Qt spectrum of Kaon selected sample
   TH1F        *fHistQt2; //Qt spectrum in Qt region of kaons
   TH1F        *fHistPtKaon; //Pt Kaon spectrum of clean sample
   TH1F        *fHistPtKPDG; //Pt Kaon spectrum , confirmed by  PDG,inside kaon Qt region
   TH1F        *fHistEta; //Eta spectrum of all kinks
   TH1F        *fHistEtaK; //Eta spectrum of kaons selected by kink topology
   TH1F        *fptKMC; //Pt Kaon spectrum MC, inside eta and pt cuts 
   TH1F        *fMultiplMC; //charge multipl MC 
   TH1F        *fESDMult; //ESD charged mult
   TH1F        *fgenpt; //Pt Kaon-Kink->mu  spectrum , MC, inside eta, Pt, radius cuts
   TH1F        *frad; //radius of kinks,  MC , inside the eta nad Pt cuts 
   TH1F        *fKinkKaon; //Pt of PDG Kaons inside the selcted ones by the KInk topology 
   TH1F        *fKinKRbn; //Pt of PDG Kaons inside the selcted ones by the KInk topology 
   TH1F        *fKinkKaonBg; //Pt of the BG inside the kink-Kaon identified spectrum
   TH1F        *fM1kaon; //inv mass of kink-tracks taken as kaons decaying to  mu + neutrino
   TH1F        *fPtKink; //Pt  spectrum   of all kinks  from track bank
   TH1F        *fptKink; //Pt  spectrum of all kinks from kink bank
   TH2F        *fAngMomK; // Decay angle vrs Mother Mom for pdg kaons
   TH2F        *fAngMomPi; // Decay angle vrs Mother Mom for pdg pions
   TH2F        *fAngMomKC; //Decay angle vrs Mother Mom for pdg kaons, inside the selected sample
   TH1F        *fMultESDK; //ESD charged mult
   TH1F        *fMultMCK; //MC K charged mult
   TH2F        *fSignPtNcl;//signPt vrs number of clusters in TPC for kaons from kink sele sample
   TH2F        *fSignPtEta;//signPt vrs Eta  in TPC for kaons from kink sele sample
   TH2F        *fEtaNcl;//Eta    vrs Nclu in TPC for kaons from kink sele sample
   TH1F        *fSignPt;//signPt  in TPC for kaons from kink sele sample
   TH2F        *fChi2NclTPC;//chi2 vrs TPC Nclusters for kaons from kink sele sample
   TH1F        *fRatChi2Ncl;// Ratio chi2/ Ncl  TPC  for kaons from kink sele sample
   TH2F        *fRadiusNcl;//kink  Radius      Ncl  TPC  for kaons from kink sele sample
   TH2F        *fTPCSgnlP;//kink  Radius      Ncl  TPC  for kaons from kink sele sample
   TH2F        *fTPCSgnlPa;//kink  Radius      Ncl  TPC  for  kink sele sample
   TH1D        *fRpr; // Radius of VTX at Y , X plane              
   TH1D        *fZpr; //Z distrio of main vertex                  
   TH1D        *fdcatoVxXY; //dca to Vertex XY  distrio                   
   TH1D        *fnSigmToVx; //nSigma to Vertex  distrio of main vertex                  
   TH2F        *fKinkMothDau; //Mother vrs Daughter                                       
   TH2F        *fZvXv; //two dime of Z vrs X of vtx main           
   TH2F        *fZvYv; // two dime of Z vrs Y of vtx main           
   TH2F        *fXvYv; // two dime of X vrs Y of main tracks vtx main           
   TH1F        *fHistPtKaoP; //Pt Kaon spectrum of clean sample pos
   TH1F        *fHistPtKaoN; //Pt Kaon spectrum of clean sample neg
   TH1F        *frapiKESD;// rapidi K      
   TH1F        *flifetime;//radius of kinks,  MC , inside the eta nad Pt cuts 
   TH1F        *fradLK;//Length  of kinks,  MC , inside the eta nad Pt cuts 
   TH3F        *fradPtRpDt;//radius of kinks,  MC , inside the eta nad Pt cuts 
   TH1F        *fInvMuNuAll;//radius of kinks,  MC , inside the eta nad Pt cuts 
   TH2F        *fQtInvM;// 
   TH1F        *fDCAkink;//!MC dcs kink
   TH2F        *fPosiKink;//!MC position  kink
   TH2F        *fPosiKinkK;//!MC position  kink
   TH2F        *fPosiKinKXZ;//!MC position  kink
   TH2F        *fPosiKinKYZ;//!MC position  kink
   TH2F        *fPosiKinKBg;//!MC position  kink
   TH2F        *fQtMothP;//!qt vrs p mother  
   TH2F        *fTPCSgnlPtpc;//Kink mother moment vrs TPC signal                 
   TH2F        *fTPCMomNSgnl;//kink  mother TPC momentum vrs nsigmas of dEdx                    
   TH2F        *fMothKinkMomSgnl;//kink  mother TPC momentum vrs nsigmas of dEdx                    
   TH1F        *fNSigmTPC;//kink  mother TPC momentum vrs nsigmas of dEdx                    
   TH2F        *fTPCSgnlKinkDau;//Kink mother moment vrs TPC signal                 
   TH1F        *fPtKinkPos; //Pos K Pt  spectrum   of all kinks  from track bank, K0 bins
   TH1F        *fPtKinkNeg; //Neg K Pt  spectrum   of all kinks  from track bank, K0 bins
   TH2F        *fRadNclCln;//kink  Radius      Ncl  TPC  for kaons from kink clean sample
   TH1F        *fRatioCrossedRows; //ratio  crossed rows                                           
   TH1F        *fRatioCrossedRowsKink; //ratio  crossed rows  for kinks                                         
   TH2F        *fRadiusPt;//kinks,  Radius      vs Pt                                        
   TH2F        *fRadiusPtcln;//kinks,  Radius      vs Pt    for clean kaons                                     
   TH2F        *fInvMassMuNuPt;//kinks,Invariant Mass MuNu     vs Pt                                         
   TH1F        *fPtCut1; //K Pt  spectrum   of all kinks  from track bank, K0 bins
   TH1F        *fPtCut2; //K Pt  spectrum   of all kinks  from track bank, K0 bins
   TH1F        *fPtCut3; //K Pt  spectrum   of all kinks  from track bank, K0 bins
   TH2F        *fAngMomKKinks;//kinks,  Angle vs Momentum  for K-kinks                                      



   TF1         *f1;
   TF1         *f2;
  TList        *fListOfHistos; // list of histos
//          Marek Multiplicity
Int_t fLowMulcut;  // 
Int_t fUpMulcut;
Int_t fKinkRadUp;
Int_t fKinkRadLow;
Int_t fLowCluster;
Float_t  fLowQt;       
AliESDtrackCuts* fCutsMul;
AliESDtrackCuts* fMaxDCAtoVtxCut;  
AliPIDResponse *fPIDResponse;     //! PID response object

  AliAnalysisKinkESDat(const AliAnalysisKinkESDat&); // not implemented
  AliAnalysisKinkESDat& operator=(const AliAnalysisKinkESDat&); // not implemented

  ClassDef(AliAnalysisKinkESDat, 1); // example of analysis
};

#endif
