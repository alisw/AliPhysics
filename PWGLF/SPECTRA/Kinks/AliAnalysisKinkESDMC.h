#ifndef ALIANALYSISKINKESDMC_H
#define ALIANALYSISKINKESDMC_H

/*  See cxx source for full Copyright notice */

//-----------------------------------------------------------------
//                 AliAnalysisKinkESDMC class
//         This task is an example of an analysis task
//                  for kink topology Study
//          Authors: Martha Spyropoulou-Stassinaki
//           and members of the Greek group at the
//          Physics Department of Athens University
//                    mspyrop@phys.uoa.gr
//-----------------------------------------------------------------

class AliPIDResponse;
class AliESDEvent;
class AliESDVertex;
class AliESDtrack;
class TF1;
class TH1F;
class TH2F;
class TH1D;
class TH2D;
class TList;
class AliESDtrackCuts;
class AliPhysicsSelection;

#include "AliAnalysisTaskSE.h"

class AliAnalysisKinkESDMC : public AliAnalysisTaskSE {
 public:
  AliAnalysisKinkESDMC(const char *name = "AliAnalysisKinkESDMC");
  virtual ~AliAnalysisKinkESDMC() {}

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
  const AliESDVertex *GetEventVertex(const AliESDEvent* esd) const;
//  const AliESDVertex *GetEventVertex(AliESDEvent* esd) ;
  void SetMulCut(Int_t low, Int_t up){fLowMulcut=low;fUpMulcut=up;}

  void SetKinkRadius(Float_t lRadiusKLow, Float_t lRadiusKUp)  { fKinkRadLow=lRadiusKLow; fKinkRadUp=lRadiusKUp;}

  
 private:
   TH1F        *fHistPtESD; //!Pt spectrum of all ESD inside eta, Pt cuts
   TH1F        *fHistPt; //!Pt spectrum of all ESD tracks
   TH1F        *fHistQtAll; //!Qt spectrum of all kinks
   TH1F        *fHistQt1; //!Qt spectrum of Kaon selected sample
   TH1F        *fHistQt2; //!Qt spectrum in Qt region of kaons
   TH1F        *fHistPtKaon; //!Pt Kaon spectrum of clean sample
   TH1F        *fHistPtKPDG; //!Pt Kaon spectrum , confirmed by  PDG,inside kaon Qt region
   TH1F        *fHistEta; //!Eta spectrum of all kinks
   TH1F        *fHistEtaK; //!Eta spectrum of kaons selected by kink topology
   TH1F        *fptKMC; //!Pt Kaon spectrum MC, inside eta and pt cuts 
   TH1F        *fMultiplMC; //!charge multipl MC 
   TH1F        *fESDMult; //!ESD charged mult
   TH1F        *frad; //!radius of kinks,  MC , inside the eta nad Pt cuts 
   TH1F        *fradMC; //!radius of kinks,  MC , inside the eta nad Pt cuts 
   TH1F        *fKinkKaon; //!Pt of PDG Kaons inside the selcted ones by the KInk topology 
   TH1F        *fKinkKaonBg; //!Pt of the BG inside the kink-Kaon identified spectrum
   TH1F        *fM1kaon; //!inv mass of kink-tracks taken as kaons decaying to  mu + neutrino
   TH1F        *fgenPtEtR; //!MC Pt spectrum of kaons decaying to muon+neutrino and pi +pi, inside eta,Pt,Rad cuts
   TH1F        *fPtKink; //!Pt  spectrum   of all kinks  from track bank
   TH2F        *fcodeH; //!PDG code(mother)  vrs PDG dcode(daughter) of kinks with Qt <0.12 (fake)
   TH2F        *fdcodeH; //!inks, code  vrs dcode of BG,if mother code is 321 and daughter code > 
   TH2F        *fAngMomK; //! Decay angle vrs Mother Mom for pdg kaons
   TH2F        *fAngMomPi; //! Decay angle vrs Mother Mom for pdg pions
   TH2F        *fAngMomKC; //!Decay angle vrs Mother Mom for pdg kaons, inside the selected sample
   TH1F        *fMultESDK; //!ESD charged mult
   TH1F        *fMultMCK; //!MC K charged mult
   TH2F        *fSignPtNcl;//!signPt vrs number of clusters in TPC for kaons from kink sele sample
   TH2F        *fSignPtEta;//!signPt vrs Eta  in TPC for kaons from kink sele sample
   TH2F        *fSignPtEtaMC;//!signPt vrs Eta  in TPC for kaons from kink sele sample
   TH1F        *fSignPtMC;//!signPt   in TPC for kaons 
   TH2F        *fEtaNcl;//!Eta    vrs Nclu in TPC for kaons from kink sele sample
   TH1F        *fSignPt;//!signPt  in TPC for kaons from kink sele sample
   TH2F        *fChi2NclTPC;//!chi2 vrs TPC Nclusters for kaons from kink sele sample
   TH1F        *fRatChi2Ncl;//! Ratio chi2/ Ncl  TPC  for kaons from kink sele sample
   TH2F        *fRadiusNcl;//! Radis  f kink vetex   for kaons from kink sele sample
   TH2F        *fTPCSgnlP;//! TPC de/dx signal      for kaons from kink sele sample
   TH2F        *fTPCSgnlPa;//! TPC de/dx signal      for  kink  sample
   TH1F        *fSignPtGen;//!signPt   in TPC for kaonsgenerated 
   TH1D        *fRpr;//! Radius of VTX at Y , X plane
   TH1D        *fZpr;//!Z distrio of main vertex                  
   TH1D        *fdcatoVxXY;//! dca to Vertex XY distr          
   TH1F        *fMCEtaKaon;//!MC eta for kaons                                     
   TH2F        *fZvXv;//! two dime of Z vrs X of vtx main           
   TH2F        *fZvYv;//! two dime of Z vrs Y of vtx main           
   TH2F        *fXvYv;//! two dime of X vrs Y of main tracks vtx main           
   TH1F        *fPtPrKink;//! pt of Primary PDG kaons inside the selected ones by the kink topology              
   TH1F        *fgenPtEtRP;//!MC Pt spectrum of kaons decaying to muon+neutrino and pi +pi, inside eta,Pt,Rad cuts
   TH1F        *fgenPtEtRN;//!MC Pt spectrum of kaons decaying to muon+neutrino and pi +pi, inside eta,Pt,Rad cuts
   TH1F        *fkinkKaonP;//!MC Pt spectrum of kaons decaying to muon+neutrino and pi +pi, inside eta,Pt,Rad cuts
   TH1F        *fkinkKaonN;//!MC Pt spectrum of kaons decaying to muon+neutrino and pi +pi, inside eta,Pt,Rad cuts
   TH1F        *frapidESDK;//!MC ESD K  rapidity distr.                                                                   
   TH1F        *frapidKMC;//!MC MC K rapidity dis
   TH1F        *fPtKPlMC;//!MC MC K rapidity dis
   TH1F        *fPtKMnMC;//!MC MC K rapidity dis
   TH1F        *fHistPtKaoP;//!MC MC K rapidity dis
   TH1F        *fHistPtKaoN;//!MC MC K rapidity dis
   TH1F        *fHiPtKPDGP;//!MC MC K rapidity dis
   TH1F        *fHiPtKPDGN;//!MC MC K rapidity dis
   TH1F        *fKinKBGP;//!MC MC K rapidity dis
   TH1F        *fKinKBGN;//!MC MC K rapidity dis
   TH1F        *fQtKMu;//!MC MC K Qt K to mu
   TH1F        *fQtKPi;//!MC MC K Qt K to mu
   TH1F        *fQtKEl;//!MC MC K Qt K to mu
   TH1F        *fFakepipi;//!MC Fake pipi
   TH1F        *fFakeKPi;//!MC Fake Kpi
   TH1F        *fDCAkink;//!MC dcs kink
   TH1F        *fDCAkinkBG;//!MC dcs kink
   TH2F        *fPosiKink;//!MC position  kink
   TH2F        *fPosiKinkK;//!MC position  kink
   TH2F        *fPosiKinKXZ;//!MC position  kink
   TH2F        *fPosiKinKYZ;//!MC position  kink
   TH2F        *fPosiKinKBgZY;//!MC position  kink
   TH2F        *fcode2;//!PDG code(mother)  vrs PDG dcode(daughter) of kinks with Qt <0.12 (fake)
   TH2F        *fcode4;//!PDG code(mother)  vrs PDG dcode(daughter) of kinks with Qt <0.12 (fake)
   TH2F        *fZkinkZDau;//!   z-position of kink z position of daughter bg                       )
   TH1F        *fQtKMuMC;//!MC MC K Qt K to mu
   TH1F        *fQtKElMC;//!MC MC K Qt K to mu
   TH1F        *fQtKPiMC;//!MC MC K Qt K to mu
   TH1F        *fQtK3PiP;//!EDS   K Qt K to 3Pi
   TH1F        *fQtK3PiM;//!EDS   K Qt K to 3Pi
   TH2F        *fmaxAngMomKmu; //!Decay angle vrs Mother Mom for pdg kaons, inside the selected sample
   TH2F        *fPosiKinKBgZX;//!MC position  kink
   TH2F        *fPosiKinKBgXY;//!MC position  kink
   TH1F        *fMinvPi;//!MC R life time     
   TH1F        *fMinvKa;//!MC R life time     
   TH1F        *fMinvPr;//!MC R life time     
   TH2F        *fTPCSgnlPtpc;//Kink mother moment vrs TPC signal                 
   TH2F        *fTPCMomNSgnl;//kink  mother TPC momentum vrs nsigmas of dEdx                    
   TH2F        *fMothKinkMomSgnl;//kink  mother TPC momentum vrs nsigmas of dEdx                    
   TH1F        *fNSigmTPC;//kink  mother TPC momentum vrs nsigmas of dEdx                    
   TH2F        *fTPCSgnlKinkDau;//Kink mother moment vrs TPC signal 
   TH2F        *fcodeDau1;//!PDG code(mother)  vrs PDG dcode(daughter) of kinks   daughters          
   TH2F        *fcodeDau2;//!PDG code(mother)  vrs PDG dcode(daughter) of kinks   daughters , Bg study
   TH2F        *fMothKinkMomSgnlD;//kink  mother TPC momentum vrs nsigmas of dEdx                    
   TH1F        *fInvMassMuNuAll;//kinks,  Inv Mass all kinks  MuNu                                       
   TH2F        *fInvMassMuNuPt;//kinks,Invariant Mass MuNu     vs Pt                                         
   TH1F        *fRatioCrossedRows; //ratio  crossed rows                                           
   TH1F        *fRatioCrossedRowsKink; //ratio  crossed rows  for kinks                                      
   TH2F        *fRadiusPt;//kinks,  Radius      vs Pt                                        
   TH2F        *fRadiusPtcln;//kinks,  Radius      vs Pt    for clean kaons                                  
   TH1F        *fPtCut1; //K Pt  spectrum   of all kinks  from track bank, K0 bins
   TH1F        *fPtCut2; //K Pt  spectrum   of all kinks  from track bank, K0 bins
   TH1F        *fPtCut3; //K Pt  spectrum   of all kinks  from track bank, K0 bins
   TH2F        *fAngMomKKinks;//kinks,  Angle vs Momentum  for K-kinks                                      
    TH1F        *flengthMCK;//!MC R life time
   TH1F        *flifetiMCK;//!MC R life time
   TH1F        *flifetim2;//!MC R life time
   TH1F        *fLHelESDK;//!MC R life time
   TH1F        *flifeInt;//!MC R life time
   TH1F        *flifeYuri;//!MC R life time
   TH1F        *flenYuri;//!MC R life time
   TH1F        *flenTrRef;//!MC R life time
   TH1F        *flifeSmall;//!MC R life time
   TH1F        *flifetime;//!MC R life time
   TH1F        *flifTiESDK;//!MC R life time
   TH1F        *flifeKink;//!MC R life time
   TH1F        *flenHelx;//!MC R life time
   TH3F        *fradPtRapMC;//!MC R life time
   TH3F        *fradPtRapDC;//!MC R life time
   TH3F        *fradPtRapESD;//!MC R life time
   TH2F        *fRadNclcln;//!MC R life time


    
   TF1         *f1;
   TF1         *f2;
  TList        *fListOfHistos; //! list of histos

 Int_t fLowMulcut; //
Int_t fUpMulcut;
Int_t fKinkRadUp;
Int_t fKinkRadLow;
AliESDtrackCuts*  fCutsMul;

     AliESDtrackCuts* fMaxDCAtoVtxCut;
 AliPIDResponse *fPIDResponse;     //! PID response object

  AliAnalysisKinkESDMC(const AliAnalysisKinkESDMC&); // not implemented
  AliAnalysisKinkESDMC& operator=(const AliAnalysisKinkESDMC&); // not implemented

  ClassDef(AliAnalysisKinkESDMC, 1); // example of analysis
};

#endif
