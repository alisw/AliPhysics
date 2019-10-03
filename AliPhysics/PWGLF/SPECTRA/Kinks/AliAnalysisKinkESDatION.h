#ifndef ALIANALYSISKINKESDatION_H
#define ALIANALYSISKINKESDatION_H

/*  See cxx source for full Copyright notice */

//-----------------------------------------------------------------
//                 AliAnalysisKinkESDatION class
//         This task is an example of an analysis task
//          for kink topology study in Pb-Pb collisions
//          Authors: Martha Spyropoulou-Stassinaki
//           and members of the Greek group at the
//          Physics Department of Athens University
//                    mspyrop@phys.uoa.gr
//-----------------------------------------------------------------
//class AliTPCPIDResponse;
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
class AliPhysicsSelection;


#include "AliAnalysisTaskSE.h"

class AliAnalysisKinkESDatION : public AliAnalysisTaskSE {
 public:
 // AliAnalysisKinkESDatION();
  AliAnalysisKinkESDatION(const char *name = "AliAnalysisKinkESDatION");
  virtual ~AliAnalysisKinkESDatION() {}

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
//  3/6   Float_t  GetSigmaToVertex(AliESDtrack* esdTrack) const;
  const AliESDVertex *GetEventVertex(const AliESDEvent* esd) const;
  
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
   TH2F        *fKinkMothDau; //Mother vrs Daughter                                       
   TH2F        *fZvXv; //two dime of Z vrs X of vtx main           
   TH2F        *fZvYv; // two dime of Z vrs Y of vtx main           
   TH2F        *fXvYv; // two dime of X vrs Y of main tracks vtx main           
   TH1F        *fPtPrKink; // pt of Primary PDG kaons inside the selected ones by the kink topology              
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
   TH2F        *fQtMothP;//!qt vrs p mother  
   TH1F        *fKinkKPt05; //Pt of  Kaons  selcted   centrality 0-5               
   TH1F        *fKinkKPt510; //Pt of  Kaons Selcted   centrality 5-10 
   TH1F        *fKinkKPt1020; //Pt of aons  selcted  centrality 10-20
   TH1F        *fKinkKPt2030; //Pt of Kaons  selcted with centrality 20-30
   TH1F        *fKinkKPt3040; //Pt of Kaons   selcted  inside centrality  30-40 
   TH1F        *fKinkKPt4050; //Pt of Kaons   selcted  inside centrality  40-50 
   TH1F        *fKinkKPt5060; //Pt of Kaons   selcted  inside centrality  50-60 
   TH1F        *fKinkKPt6070; //Pt of Kaons   selcted  inside centrality  60-70 
   TH1F        *fKinkKPt7080; //Pt of Kaons   selcted  inside centrality  70-80 
   TH1F        *fKinkKPt8090; //Pt of Kaons   selcted  inside centrality  80-90 
   TH1F        *fKinkMul05; //Pt of  Kaons  selcted   centrality 0-5               
   TH1F        *fKinkMul510; //Pt of  Kaons Selcted   centrality 5-10 
   TH1F        *fKinkMul1020; //Pt of aons  selcted  centrality 10-20
   TH1F        *fKinkMul2030; //Pt of Kaons  selcted with centrality 20-30
   TH1F        *fKinkMul3040; //Pt of Kaons   selcted  inside centrality  30-40 
   TH1F        *fKinkMul4050; //Pt of Kaons   selcted  inside centrality  40-50 
   TH1F        *fKinkMul5060; //Pt of Kaons   selcted  inside centrality  50-60 
   TH1F        *fKinkMul6070; //Pt of Kaons   selcted  inside centrality  60-70 
   TH1F        *fKinkMul7080; //Pt of Kaons   selcted  inside centrality  70-80 
   TH1F        *fKinkMul8090; //Pt of Kaons   selcted  inside centrality  80-90 
   TH2F        *fRadiusNclAll;//kink  Radius      Ncl  TPC  for kaons from kink sele sample
   TH2F        *fRadiusNclK;//kink  Radius      Ncl  TPC  for kaons from kink sele sample
   TH2F        *fRadiusNclClean;//kink  Radius      Ncl  TPC  for kaons from kink sele sample
   TH1F        *fRatioCrossedRows;//kink      
   TH1F        *fRatioCrossedRowsKink;//kink      

   TF1         *f1;
   TF1         *f2;
  TList        *fListOfHistos; // list of histos
  AliESDtrackCuts* fMaxDCAtoVtxCut;
 AliPIDResponse *fPIDResponse;     //! PID response object

  AliAnalysisKinkESDatION(const AliAnalysisKinkESDatION&); // not implemented
  AliAnalysisKinkESDatION& operator=(const AliAnalysisKinkESDatION&); // not implemented

  ClassDef(AliAnalysisKinkESDatION, 1); // example of analysis
};

#endif
