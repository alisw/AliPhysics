#ifndef ALIANALYSISKINKTASKPbPb2018_H
#define ALIANALYSISKINKTASKPbPb2018_H

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
#include "AliEventCuts.h"

class AliPIDResponse;
class AliESDVertex;
class AliESDEvent;
class AliESDtrack;
class TF1;
class TH1F;
class TH2F;
class TH3F;
class TH1D;
class TH2D;
class TList;
//class AliAnalysisUtils;
class AliESDtrackCuts;

#include "AliAnalysisTaskSE.h"

class AliAnalysisKinkTaskPbPb2018 : public AliAnalysisTaskSE {
 public:
  AliAnalysisKinkTaskPbPb2018(const char *name = "AliAnalysisKinkTaskPbPb2018");
  virtual ~AliAnalysisKinkTaskPbPb2018() {}

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
  const AliESDVertex *GetEventVertex(const AliESDEvent* esd) const;
   //
       Bool_t IsGoodSPDvertexRes(const AliESDVertex * spdVertex = NULL);
          Bool_t selectVertex2015pp(AliESDEvent *esd,Bool_t checkSPDres, Bool_t 
          requireSPDandTrk, Bool_t checkProximity);
   


//*
//           Marek multiplicity bins
 // void SetMulCut(Int_t low, Int_t up){fLowMulcut=low;fUpMulcut=up;}	  
  void SetKinkRadius(Float_t lRadiusKLow, Float_t lRadiusKUp)  { fKinkRadLow=lRadiusKLow; fKinkRadUp=lRadiusKUp;}	  

  void SetNClusterCut(Int_t lowCluster){fLowCluster=lowCluster;}	  
  void SetQtCut(Float_t   lowQt){fLowQt=lowQt;}	  
  void SetYKRange(Float_t  RapidityK){fRapiK=RapidityK;}

  AliEventCuts fEventCuts; 
//*/
 private:
   TH1F        *fMultiplicityBeforeCuts; //charge multiplicity ESD
   TH1F        *fIncompletEv; //charge multiplicity ESD ater Incompl. cut
   TH1F        *fMultiplicityAfterTriggerBit; // charge  multiplicity after trigger
   TH1F        *fMultiplicityAfterPileup; // c
   TH1D        *fZMainVx;//ZMain vertex distibution                                
   TH1F        *fMultiplicityAfterVertexCut;//ESD Multiplicity inside ZMain vertex                              
//
   TH1F        *fTrackPtAfterTrackCuts; //Pt spectrum of all ESD inside eta, Pt cuts
   TH1F        *fTrackPtAll; //Pt spectrum of all ESD tracks
   TH1F        *fHistQtAll; //Qt spectrum of all kinks
   TH1F        *fQtAfterAcceptance; //Qt spectrum in Qt region of kaons
   TH1F        *fPtSelectedKaons; //Pt Kaon spectrum of clean sample
   TH1F        *fKaonKinkPtAfterKinkNclCut; //Pt Kaon spectrum , confirmed by  PDG,inside kaon Qt region
   TH1F        *fEtaAfterAcceptance; //Eta spectrum of all kinks
   TH1F        *fTrackEtaSelectedKaons; //Eta spectrum of kaons selected by kink topology
   TH1F        *fPtAllKinks; //Pt Kaon spectrum MC, inside eta and pt cuts 
   TH1F        *fgenpt; //Pt Kaon-Kink->mu  spectrum , MC, inside eta, Pt, radius cuts
   TH1F        *fRadiusSelectedKinks; //radius of kinks,  MC , inside the eta nad Pt cuts 
   TH1F        *fKinKRbn; //Pt of PDG Kaons inside the selcted ones by the KInk topology 
   TH1F        *fKinkMomFromMother; //Pt of the BG inside the kink-Kaon identified spectrum
   TH1F        *fInvMassKaonInR; //inv mass of kink-tracks taken as kaons decaying to  mu + neutrino
   TH1F        *fPtFromMotherAllKinks; //inv
   TH2F        *fAngMomK; // Decay angle vrs Mother Mom for pdg kaons
   TH2F        *fAngleVsMomentumKaonsInR; //Decay angle
   TH2F        *fAngleMomKaonsinTPC; // Decay angle vrs Mother Mom for pdg pions
   TH2F        *fNclinTPCVsSignedPtSelectedKaons;//signPt vrs number of clusters in TPC for kaons from kink sele sample
   TH2F        *fRapidityVsSignedPtSelectedKons;//signPt vrs Eta  in TPC for kaons from kink sele sample
   TH2F        *fNclVsRapiditySelectedKAons;//Eta    vrs Nclu in TPC for kaons from kink sele sample
   TH1F        *fSignedPtSelectedKaons;//signPt  in TPC for kaons from kink sele sample
   TH2F        *fNclVsChi2SelectedKaons;//chi2 vrs TPC Nclusters for kaons from kink sele sample
   TH1F        *fChi2perTPCclusterSelectedKaons;// Ratio chi2/ Ncl  TPC  for kaons from kink sele sample
   TH2F        *fRadiusVsNclInR;//kink  Radius      Ncl  TPC  for kaons from kink sele sample
   TH2F        *fTPCSignalVsMomSelectedKaons;//kink  Radius      Ncl  TPC  for kaons from kink sele sample
   TH2F        *fTPCSgnlPa;//kink  Radius      Ncl  TPC  for  kink sele sample
   TH1D        *fDCAz; // Radius of VTX at Y , X plane              
   TH1D        *fDCAxy; //dca to Vertex XY  distrio                   
   TH1D        *fnSigmToVx; //nSigma to Vertex  distrio of main vertex                  
   TH2F        *fKinkMothDau; //Mother vrs Daughter                                       
   TH2F        *fZvXv; //two dime of Z vrs X of vtx main           
   TH2F        *fZvYv; // two dime of Z vrs Y of vtx main           
   TH2F        *fXvYv; // two dime of X vrs Y of main tracks vtx main           
   TH1F        *fRapiditySelectedKaons;// rapidi K      
   TH1F        *fLifetimeSelectedKaons;//radius of kinks,  MC , inside the eta nad Pt cuts 
   TH1F        *fKaonLifetimeSelectedKaons;//Length  of kinks,  MC , inside the eta nad Pt cuts 
   //TH3F        *fradPtRpDt;//radius of kinks,  MC , inside the eta nad Pt cuts 
   TH1F        *fInvMassMuNuKaonTPC;//radius of kinks,  MC , inside the eta nad Pt cuts 
   TH2F        *fQtInvMassKaonTPC;// 
   TH1F        *fDCAkinkSelectedKaons;//!MC dcs kink
   TH2F        *fxyKinkPosition;//!MC position  kink
   TH2F        *fPosiKinkK;//!MC position  kink
   TH2F        *fPosiKinKXZ;//!MC position  kink
   TH2F        *fPosiKinKYZ;//!MC position  kink
   TH2F        *fZKinkProductionVsKinkRadSelectedKaons;//!MC position  kink
   TH2F        *fQtVsKinkMomAfterAcceptance;//!qt vrs p mother  
   TH2F        *fdedxMthVsTPCMomSelectedKaons;//Kink mother moment vrs TPC signal                 
   TH2F        *fNsigmaVsTPCmomSelectedKaons;//kink  mother TPC momentum vrs nsigmas of dEdx                    
   TH2F        *fSignalMthVsSignalDaughterSelectedKaons;//kink  mother TPC momentum vrs nsigmas of dEdx                    
   TH1F        *fNsigmaSelectedKaons;//kink  mother TPC momentum vrs nsigmas of dEdx                    
   TH2F        *fSignalDaughterVsDaughterMomSelectedKaons;//Kink mother moment vrs TPC signal                 
   TH1F        *fPtPositiveSelectedKaons; //Pos K Pt  spectrum   of all kinks  from track bank, K0 bins
   TH1F        *fPtNegativeSelectedKaons; //Neg K Pt  spectrum   of all kinks  from track bank, K0 bins
   TH2F        *fNclustersVsRadiusSelectedKaons;//kink  Radius      Ncl  TPC  for kaons from kink clean sample
   TH1F        *fRatioCrossedRows; //ratio  crossed rows                                           
   TH1F        *fRatioCrossedRowsKink; //ratio  crossed rows  for kinks                                         
   TH2F        *fRadiusVsPtKaonTPC;//kinks,  Radius      vs Pt                                        
   TH2F        *fPtVsRadiusSelectedKaons;//kinks
   TH2F        *fPtVsInvMassSelectedKaons;//kinks
   TH2F        *fInvMassPtKaonTPC;//kinks,Invariant Mass MuNu     vs Pt                                         
   TH1F        *fPtKaonInTPC; //K Pt  spectrum   of all kinks  from track bank, K0 bins
   TH1F        *fPtPreSelectedkinks; //K Pt  spectrum   of all kinks  from track bank, K0 bins
   TH1F        *fPtKinkBeforedEdx; //K
   TH2F        *fAngleVsMomPreSelectedKinks; //K
   TH1F        *fPtKinkK0; //Sum K Pt  spectrum   of all kinks  from track bank, K0 bins
   TH1F        *fPtKinkK0P; //Pos K Pt  spectrum   of all kinks  from track bank, K0 bins
   TH1F        *fPtKinkK0N; //Neg K Pt  spectrum   of all kinks  from track bank, K0 bins
   TH1F        *fPtKinkGyu; //Pt  spectrum   of all kinks with High pt binning 
   TH1F        *fPtKinkGyuP;//Pt  spectrum   of Pos kinks with High pt binning 
   TH1F        *fPtKinkGyuN;//Pt  spectrum   of Neg kinks with High pt binning 

   TF1         *f1;
   TF1         *f2;
   TList       *fListOfHistos; // list of histos
   
   Float_t fKinkRadUp;
   Float_t fKinkRadLow;
   Int_t fLowCluster;
   Float_t  fLowQt;       
   Float_t  fRapiK;
   //AliESDtrackCuts* fCutsMul;
   //AliESDtrackCuts* fMaxDCAtoVtxCut;  
   //AliPIDResponse *fPIDResponse;     //! PID response object
   //*/
   AliPIDResponse *fPIDResponse;     //! PID response object
   TH1F        *fMultiBin1ChargedMulti; //  charged multiplicity in multiplicity bins
   TH1F        *fMultiBin2ChargedMulti; //  charged multiplicity in multiplicity bins
   TH1F        *fMultiBin3ChargedMulti; //  charged multiplicity in multiplicity bins

   TH1F        *fMultiBin1KaonKinksPos; //  pos kaons in multiplicity bins
   TH1F        *fMultiBin2KaonKinksPos; //  pos kaons in multiplicity bins
   TH1F        *fMultiBin3KaonKinksPos; //  pos kaons in multiplicity bins

   TH1F        *fMultiBin1KaonKinksNeg; //  neg kaons in multiplicity bins
   TH1F        *fMultiBin2KaonKinksNeg; //  neg kaons in multiplicity bins
   TH1F        *fMultiBin3KaonKinksNeg; //  neg kaons in multiplicity bins

//   AliAnalysisUtils   *fUtils; // Analysis Utils
   TH1D        *fMultiBin1Vertex;       // multi vtx
   TH1D        *fMultiBin2Vertex;       // multi vtx
   TH1D        *fMultiBin3Vertex;       // multi vtx

   TH1D        *fVertexNet;       // multi vtx
   TH1F        *fhPS;       // multi vtx
   TH1F        *fhvtxP;       // multi vtx
   AliESDtrackCuts * fESDtrackCuts;
   TH1F        *fpercentile; // percentile distribution   
   
   AliAnalysisKinkTaskPbPb2018(const AliAnalysisKinkTaskPbPb2018&); // not implemented
   AliAnalysisKinkTaskPbPb2018& operator=(const AliAnalysisKinkTaskPbPb2018&); // not implemented
   
   ClassDef(AliAnalysisKinkTaskPbPb2018, 1); // example of analysis
};

#endif
