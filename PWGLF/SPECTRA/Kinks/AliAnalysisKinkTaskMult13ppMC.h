#ifndef ALIANALYSISKINKTASKMULT13ppMC_H
#define ALIANALYSISKINKTASKMULT13ppMC_H

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
class AliESDEvent;
class AliESDtrack;
class TF1;
class TH1F;
class TH2F;
class TH3F;
class TH1D;
class TH2D;
class TList;
class AliPIDResponse;
//class AliESDtrackCuts;

#include "AliAnalysisTaskSE.h"

class AliAnalysisKinkTaskMult13ppMC : public AliAnalysisTaskSE {
 public:
  AliAnalysisKinkTaskMult13ppMC(const char *name = "AliAnalysisKinkTaskMult13ppMC");
  virtual ~AliAnalysisKinkTaskMult13ppMC() {}

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
//*/
 private:
   TH1F        *fMultiplicityBeforeCuts; //!charge multiplicity ESD   okkkk
   TH1F        *fIncompletEv; //!charge multiplicity ESD ater Incompl. cut  okkkkk
   TH1F        *fMultiplicityAfterTriggerBit; //! charge  multiplicity after trigger   okkkk
   TH1F        *fMultiplicityAfterPileup; //! c   okkkk
   TH1D        *fZMainVx;//!ZMain vertex distibution    okkkk                             
   TH1F        *fMultiplicityAfterVertexCut;//!ESD Multiplicity inside ZMain vertex      okkkkk                        
//

   TF1         *f1;//!
   TF1         *f2;//!
   TList       *fListOfHistos; //! list of histos
  
   Float_t fKinkRadUp;
   Float_t fKinkRadLow;
   Int_t fLowCluster;
   Float_t  fLowQt;       
   Float_t  fRapiK;
   //AliESDtrackCuts* fCutsMul;
   //AliESDtrackCuts* fMaxDCAtoVtxCut;  
   //AliPIDResponse *fPIDResponse;     //! PID response object
   //*/
   // AliPIDResponse *fPIDResponse;     //! PID response object


   TH1F        *fptKMC; //!Pt Kaon spectrum MC, inside eta and pt cuts   okkkkk
   TH1F        *fPtKPlMC;//!MC MC K rapidity dis   okkk
   TH1F        *fPtKMnMC;//!MC MC K rapidity dis    okkkk
   TH1F        *frapidKMC;//!MC MC K rapidity dis   okkkkk
   TH1F        *flenTrRef;//!MC R life time    okkkk
   TH1F        *flifeSmall;//!MC R life time   okkkk
   TH1F        *fLifeMCProcess13;//!MC R life time   okkkk
   TH1F        *fLifeP;//!MC R life time   okkk
   TH1F        *fLengthRZ;//!MC R life time   okkkk
   TH1F        *fLengthP;//!MC R life time   okkkk
   TH1F        *flengthKMC;//!MC R life time    okkkk
   TH1F        *fLifeMCProcess4;//!MC R life time    okkkk
   //   TH3F        *fradPtRapMC;//!MC R life time
   TH1F        *fLifeDECAYmuonORelectron;//!MC R life time   okkkk
   TH2F        *fmaxAngleVsMomKmuKel; //!Decay angle vrs Mother Mom    okkkk

   //   TH3F        *fradiusPtRapidityKplusMu;//!MC R life time
   TH1F        *fradiusKplusMu;//!MC   okkkkk
   TH1F        *fQtKplusMu;//!MC     okkkkk
   TH2F        *fPtRapidityKplusMu;//!MC    okkkk
   TH1F        *fLifetimeKplusMu;//!MC   okkkk

   //  TH3F        *fradiusPtRapidityKminusMu;//!MC R life time
   TH1F        *fradiusKminusMu;//!MC  okkkk
   TH1F        *fQtKminusMu;//!MC   okkkk
   TH2F        *fPtRapidityKminusMu;//!MC    okkkk
   TH1F        *fLifetimeKminusMu;//!MC    okkkk

   TH1F        *fMultiBin1KaonMCPos; //!  pos kaons in multiplicity bins
   TH1F        *fMultiBin2KaonMCPos; //!  pos kaons in multiplicity bins
   TH1F        *fMultiBin3KaonMCPos; //!  pos kaons in multiplicity bins
   TH1F        *fMultiBin4KaonMCPos; //!  pos kaons in multiplicity bins
   TH1F        *fMultiBin5KaonMCPos; //!  pos kaons in multiplicity bins
   TH1F        *fMultiBin6KaonMCPos; //!  pos kaons in multiplicity bins
   TH1F        *fMultiBin7KaonMCPos; //!  pos kaons in multiplicity bins
   TH1F        *fMultiBin8KaonMCPos; //!  pos kaons in multiplicity bins
   TH1F        *fMultiBin9KaonMCPos; //!  pos kaons in multiplicity bins
   TH1F        *fMultiBin10KaonMCPos; //!  pos kaons in multiplicity bins

   TH1F        *fMultiBin1KaonMCNeg; //!  pos kaons in multiplicity bins
   TH1F        *fMultiBin2KaonMCNeg; //!  pos kaons in multiplicity bins
   TH1F        *fMultiBin3KaonMCNeg; //!  pos kaons in multiplicity bins   
   TH1F        *fMultiBin4KaonMCNeg; //!  pos kaons in multiplicity bins
   TH1F        *fMultiBin5KaonMCNeg; //!  pos kaons in multiplicity bins
   TH1F        *fMultiBin6KaonMCNeg; //!  pos kaons in multiplicity bins
   TH1F        *fMultiBin7KaonMCNeg; //!  pos kaons in multiplicity bins
   TH1F        *fMultiBin8KaonMCNeg; //!  pos kaons in multiplicity bins
   TH1F        *fMultiBin9KaonMCNeg; //!  pos kaons in multiplicity bins
   TH1F        *fMultiBin10KaonMCNeg; //!  pos kaons in multiplicity bins
   TH1F        *flifetime;//!!MC R life time
   TH1F        *fMultiBin1fPtKplusMu; //! !multi vtx
   TH1F        *fMultiBin1fPtKminusMu; //! !multi vtx
   TH1F        *fMultiBin1fPtKplusEl; //! !multi vtx
   TH1F        *fMultiBin1fPtKminusEl; //! !multi vtx
   TH1F        *fMultiBin1fPtKPiPlus; //! !multi vtx
   TH1F        *fMultiBin1fPtKPiMinus; //! !multi vtx

   TH1F        *fMultiBin2fPtKplusMu; //! !multi vtx
   TH1F        *fMultiBin2fPtKminusMu; //! !multi vtx
   TH1F        *fMultiBin2fPtKplusEl; //! !multi vtx
   TH1F        *fMultiBin2fPtKminusEl; //! !multi vtx
   TH1F        *fMultiBin2fPtKPiPlus; //! !multi vtx
   TH1F        *fMultiBin2fPtKPiMinus; //! !multi vtx

   TH1F        *fMultiBin3fPtKplusMu; //! !multi vtx
   TH1F        *fMultiBin3fPtKminusMu; //! !multi vtx
   TH1F        *fMultiBin3fPtKplusEl; //! !multi vtx
   TH1F        *fMultiBin3fPtKminusEl; //! !multi vtx
   TH1F        *fMultiBin3fPtKPiPlus; //! !multi vtx
   TH1F        *fMultiBin3fPtKPiMinus; //! !multi vtx

   TH1F        *fMultiBin4fPtKplusMu; //! !multi vtx
   TH1F        *fMultiBin4fPtKminusMu; //! !multi vtx
   TH1F        *fMultiBin4fPtKplusEl; //! !multi vtx
   TH1F        *fMultiBin4fPtKminusEl; //! !multi vtx
   TH1F        *fMultiBin4fPtKPiPlus; //! !multi vtx
   TH1F        *fMultiBin4fPtKPiMinus; //! !multi vtx

   TH1F        *fMultiBin5fPtKplusMu; //! !multi vtx
   TH1F        *fMultiBin5fPtKminusMu; //! !multi vtx
   TH1F        *fMultiBin5fPtKplusEl; //! !multi vtx
   TH1F        *fMultiBin5fPtKminusEl; //! !multi vtx
   TH1F        *fMultiBin5fPtKPiPlus; //! !multi vtx
   TH1F        *fMultiBin5fPtKPiMinus; //! !multi vtx

   TH1F        *fMultiBin6fPtKplusMu; //! !multi vtx
   TH1F        *fMultiBin6fPtKminusMu; //! !multi vtx
   TH1F        *fMultiBin6fPtKplusEl; //! !multi vtx
   TH1F        *fMultiBin6fPtKminusEl; //! !multi vtx
   TH1F        *fMultiBin6fPtKPiPlus; //! !multi vtx
   TH1F        *fMultiBin6fPtKPiMinus; //! !multi vtx

   TH1F        *fMultiBin7fPtKplusMu; //! !multi vtx
   TH1F        *fMultiBin7fPtKminusMu; //! !multi vtx
   TH1F        *fMultiBin7fPtKplusEl; //! !multi vtx
   TH1F        *fMultiBin7fPtKminusEl; //! !multi vtx
   TH1F        *fMultiBin7fPtKPiPlus; //! !multi vtx
   TH1F        *fMultiBin7fPtKPiMinus; //! !multi vtx

   TH1F        *fMultiBin8fPtKplusMu; //! !multi vtx
   TH1F        *fMultiBin8fPtKminusMu; //! !multi vtx
   TH1F        *fMultiBin8fPtKplusEl; //! !multi vtx
   TH1F        *fMultiBin8fPtKminusEl; //! !multi vtx
   TH1F        *fMultiBin8fPtKPiPlus; //! !multi vtx
   TH1F        *fMultiBin8fPtKPiMinus; //! !multi vtx

   TH1F        *fMultiBin9fPtKplusMu; //! !multi vtx
   TH1F        *fMultiBin9fPtKminusMu; //! !multi vtx
   TH1F        *fMultiBin9fPtKplusEl; //! !multi vtx
   TH1F        *fMultiBin9fPtKminusEl; //! !multi vtx
   TH1F        *fMultiBin9fPtKPiPlus; //! !multi vtx
   TH1F        *fMultiBin9fPtKPiMinus; //! !multi vtx

   TH1F        *fMultiBin10fPtKplusMu; //! !multi vtx
   TH1F        *fMultiBin10fPtKminusMu; //! !multi vtx
   TH1F        *fMultiBin10fPtKplusEl; //! !multi vtx
   TH1F        *fMultiBin10fPtKminusEl; //! !multi vtx
   TH1F        *fMultiBin10fPtKPiPlus; //! !multi vtx
   TH1F        *fMultiBin10fPtKPiMinus; //! !multi vtx

   TH1F        *fradiusKplusEl;//!MC
   TH1F        *fQtKplusEl;//!MC
   TH2F        *fPtRapidityKplusEl;//!MC
   TH1F        *fLifetimeKplusEl;//!MC

//   TH3F        *fradiusPtRapidityKminusEl;//!MC R life time
   TH1F        *fradiusKminusEl;//!MC
   TH1F        *fQtKminusEl;//!MC
   TH2F        *fPtRapidityKminusEl;//!MC
   TH1F        *fLifetimeKminusEl;//!MC

   TH1F        *fradiusKPiPlus;//!MC
   TH1F        *fQtKPiPlus;//!MC
   TH2F        *fPtRapidityKPiPlus;//!MC
   TH1F        *fLifetimeKPiPlus;//!MC

   //   TH3F        *fradiusPtRapidityKPiMinus;//!MC R life time
   TH1F        *fradiusKPiMinus;//!MC
   TH1F        *fQtKPiMinus;//!MC
   TH2F        *fPtRapidityKPiMinus;//!MC
   TH1F        *fLifetimeKPiMinus;//!MC
   AliPIDResponse *fPIDResponse;//!
   TH1F        *fTrackPtAll; //!Pt spectrum of all ESD tracks
   TH1F        *fRatioCrossedRows; //!ratio  crossed rows
   TH2F        *fZvXv; //!two dime of Z vrs X of vtx main           
   TH2F        *fZvYv; //! two dime of Z vrs Y of vtx main           
   TH2F        *fXvYv; //! two dime of X vrs Y of main tracks vtx main  
   TH1D        *fDCAz; //! Radius of VTX at Y , X plane              
   TH1D        *fDCAxy; //!dca to Vertex XY  distrio
   TH1F        *fMultiBin1ChargedMulti; //!  charged multiplicity in multiplicity bins
   TH1F        *fMultiBin2ChargedMulti; //!  charged multiplicity in multiplicity bins
   TH1F        *fMultiBin3ChargedMulti; //!  charged multiplicity in multiplicity bins
   TH1F        *fMultiBin4ChargedMulti; //!  charged multiplicity in multiplicity bins
   TH1F        *fMultiBin5ChargedMulti; //!  charged multiplicity in multiplicity bins
   TH1F        *fMultiBin6ChargedMulti; //!  charged multiplicity in multiplicity bins
   TH1F        *fMultiBin7ChargedMulti; //!  charged multiplicity in multiplicity bins
   TH1F        *fMultiBin8ChargedMulti; //!  charged multiplicity in multiplicity bins
   TH1F        *fMultiBin9ChargedMulti; //!  charged multiplicity in multiplicity bins
   TH1F        *fMultiBin10ChargedMulti; //! charged multiplicity in multiplicity bins

   TH1D        *fMultiBin1Vertex;       //! multi vtx
   TH1D        *fMultiBin2Vertex;       //! multi vtx
   TH1D        *fMultiBin3Vertex;       //! multi vtx
   TH1D        *fMultiBin4Vertex;       //! multi vtx
   TH1D        *fMultiBin5Vertex;       //! multi vtx
   TH1D        *fMultiBin6Vertex;       //! multi vtx
   TH1D        *fMultiBin7Vertex;       //! multi vtx
   TH1D        *fMultiBin8Vertex;       //! multi vtx
   TH1D        *fMultiBin9Vertex;       //! multi vtx
   TH1D        *fMultiBin10Vertex;       //! multi vtx
   TH1F        *fTrackPtAfterTrackCuts; //!Pt spectrum of all ESD inside eta, Pt cuts
   TH1F        *fPtAllKinks; //!Pt Kaon spectrum MC, inside eta and pt cuts
   TH1F        *fRatioCrossedRowsKink; //!ratio  crossed rows  for kinks  
   TH2F        *fxyKinkPosition;//!!MC position  kink
   TH1F        *fHistQtAll; //!Qt spectrum of all kinks
   TH1F        *fPtFromMotherAllKinks; //!inv
   TH2F        *fQtVsKinkMomAfterAcceptance;//!!qt vrs p mother
   TH1F        *fQtAfterAcceptance; //!Qt spectrum in Qt region of kaons
   TH1F        *fEtaAfterAcceptance; //!Eta spectrum of all kinks
   TH1F        *fKinkMomFromMother; //!Pt of the BG inside the kink-Kaon identified spectrum
   TH1F        *fQtBeforeAngleCut; //!!Qt spectrum of Kaon selected sample
   TH1F        *fFakepipi;//!!MC Fake pipi
   TH1F        *fFakeKPi;//!!MC Fake Kpi
   TH1F        *fPtKaonInTPC; //!K Pt  spectrum   of all kinks  from track bank, K0 bins
   TH2F        *fAngleMomKaonsinTPC; //! Decay angle vrs Mother Mom for pdg pions
   TH2F        *fRadiusPtPion;//!kinks,  Radius      vs Pt    for clean kaons
   TH2F        *fRadiusPtKaon;//!kinks,  Radius      vs Pt    for clean kaons 
   TH1F        *fQtKMu;//!MC MC K Qt K to mu
   TH1F        *fQtKPi;//!MC MC K Qt K to mu
   TH1F        *fQtKEl;//!MC MC K Qt K to mu
   TH1F        *fQtK3PiP;//!EDS   K Qt K to 3Pi
   TH1F        *fQtK3PiM;//!EDS   K Qt K to 3Pi
   TH1F        *fPtKaonPDG; //!Pt Kaon
   TH1F        *fPtKaonPDGPos; //!Pt Kaon
   TH1F        *fPtKaonPDGNeg; //!Pt Kaon
   TH1F        *fKaonPDGEta; //!Eta spectrum
   TH1F        *fKaonPDGrapidity;//!MC
   TH2F        *fKaonPDGpTvsRadius;//kinks
   TH1F        *fKaonPDGqT; //!Qt spectrum
   TH2F        *fAngMomKaonPDG; // Decay angle vrs Mother Mom for pdg kaons
   TH2F        *fAngMomPi; //! Decay angle vrs Mother Mom for pdg pions
   TH2F        *fQtInvMassKaonTPC;//
   TH2F        *fInvMassPtKaonTPC;//kinks,Invariant Mass MuNu     vs Pt
   TH2F        *fRadiusVsPtKaonTPC;//kinks,  Radius      vs Pt 
   TH2F        *fSignPtNcl;//!signPt vrs number of clusters in TPC for kaons from kink sele sample
   TH2F        *fAngleVsMomentumKaonsInR; //Decay angle
   TH1F        *fInvMassKaonInR; //inv mass of kink-tracks taken as kaons decaying to  mu + neutrino
   TH2F        *fRadiusVsNclInR;//kink  Radius      Ncl  TPC  for kaons from kink sele sample
   TH1F        *fkaonToMu; //!inv mass of kink-tracks taken as kaons decaying to  mu + neutrino
   TH1F        *fkaonToPi; //!inv mass of kink-tracks taken as kaons decaying to  mu + neutrino
   TH1F        *fkaonToKa; //!inv mass of kink-tracks
   TH2F        *fRadiusNcl;//! Radis  f kink vetex   for kaons from kink selected sample
   TH2F        *fcodeMotherVsDaughter; //!PDG code(mother)  vrs PDG dcode(daughter)
   TH2F        *fZMotherVsDaughter;//!   z-position of kink z position of daughter bg
   TH2F        *fPtVsRadiusFake;//!kinks, Radius vs Pt for fake kaons
   TH1F        *fKaonKinkPtAfterKinkNclCut; //!Pt Kaon spectrum , confirmed by  PDG,inside kaon Qt region
   TH2F        *fAngleVsMomPreSelectedKinks; //!
   TH1F        *fPtPreSelectedkinks; //!K Pt  spectrum   of all kinks  from track bank, K0 bins
   TH1F        *fPtKinkBeforedEdx; //!K
   TH2F        *fTPCSgnlPa;//!kink  Radius      Ncl  TPC  for  kink sele sample
   TH2F        *fTPCSignalVsMomSelectedKaons;//!kink  Radius      Ncl  TPC  for kaons from kink sele sample
   TH2F        *fNclustersVsRadiusSelectedKaons;//!kink  Radius      Ncl  TPC  for kaons from kink clean sample
   TH2F        *fPtVsRadiusSelectedKaons;//!kinks
   TH2F        *fPtVsInvMassSelectedKaons;//!kinks
   TH1F        *fPtKaonRECpos;//!MC MC K rapidity dis
   TH1F        *fPtKaonRECneg;//!MC MC K rapidity dis
   TH1F        *fPtKaonREC;//!MC MC
   TH2F        *fdedxMthVsTPCMomSelectedKaons;//!Kink mother moment vrs TPC signal
   TH2F        *fSignalMthVsSignalDaughterSelectedKaons;//!kink  mother TPC momentum vrs nsigmas of dEdx
   TH2F        *fSignalDaughterVsDaughterMomSelectedKaons;//!Kink mother moment vrs TPC signal 
   TH2F        *fNsigmaVsTPCmomSelectedKaons;//!kink  mother TPC momentum vrs nsigmas of dEdx
   TH1F        *fRadiusSelectedKinks; //!radius of kinks,  MC , inside the eta nad Pt cuts
   TH1F        *fKaonLifetimeSelectedKaons;//!Length  of kinks,  MC , inside the eta nad Pt cuts
   TH1F        *fTrackEtaSelectedKaons; //!Eta spectrum of kaons selected by kink topology
   TH1F        *fRapiditySelectedKaons;//! rapidi K
   TH2F        *fZKinkProductionVsKinkRadSelectedKaons;//!MC position  kink
   TH2F        *fNclinTPCVsSignedPtSelectedKaons;//!signPt vrs number of clusters in TPC for kaons from kink sele sample
   TH2F        *fRapidityVsSignedPtSelectedKons;//!signPt vrs Eta  in TPC for kaons from kink sele sample
   TH2F        *fNclVsRapiditySelectedKAons;//!Eta    vrs Nclu in TPC for kaons from kink sele sample
   TH2F        *fNclVsChi2SelectedKaons;//!chi2 vrs TPC Nclusters for kaons from kink sele sample
   TH1F        *fChi2perTPCclusterSelectedKaons;//! Ratio chi2/ Ncl  TPC  for kaons from kink sele sample
   TH1F        *fLifetimeSelectedKaons;//!radius of kinks,  MC , inside the eta nad Pt cuts
   TH1F        *fPtSelectedKaons; //!Pt Kaon spectrum of clean sample
   TH1F        *fDCAkinkSelectedKaons;//!MC dcs kink
   TH1F        *fPtPositiveSelectedKaons; //!Pos K Pt  spectrum   of all kinks  from track bank, K0 bins
   TH1F        *fPtNegativeSelectedKaons; //!Neg K Pt  spectrum   of all kinks  from track bank, K0 bins
   TH1F        *fkinkKaonPDG;//!MC MC
   TH1F        *fkinkKaonPDGpos;//!MC MC
   TH1F        *fkinkKaonPDGneg;//!MC MC
   TH1F        *fkinkKaonPDGposMulti1;//!MC MC
   TH1F        *fkinkKaonPDGnegMulti1;//!MC MC
   TH1F        *fkinkKaonPDGposMulti2;//!MC MC
   TH1F        *fkinkKaonPDGnegMulti2;//!MC MC
   TH1F        *fkinkKaonPDGposMulti3;//!MC MC
   TH1F        *fkinkKaonPDGnegMulti3;//!MC MC
   TH1F        *fkinkKaonPDGposMulti4;//!MC MC
   TH1F        *fkinkKaonPDGnegMulti4;//!MC MC
   TH1F        *fkinkKaonPDGposMulti5;//!MC MC
   TH1F        *fkinkKaonPDGnegMulti5;//!MC MC
   TH1F        *fkinkKaonPDGposMulti6;//!MC MC
   TH1F        *fkinkKaonPDGnegMulti6;//!MC MC
   TH1F        *fkinkKaonPDGposMulti7;//!MC MC
   TH1F        *fkinkKaonPDGnegMulti7;//!MC MC
   TH1F        *fkinkKaonPDGposMulti8;//!MC MC
   TH1F        *fkinkKaonPDGnegMulti8;//!MC MC
   TH1F        *fkinkKaonPDGposMulti9;//!MC MC
   TH1F        *fkinkKaonPDGnegMulti9;//!MC MC
   TH1F        *fkinkKaonPDGposMulti10;//!MC MC
   TH1F        *fkinkKaonPDGnegMulti10;//!MC MC
   TH1F        *fkinkKaonPDGBkg;//!MC MC
   TH1F        *fkinkKaonPDGBkgMulti1;//!MC MC
   TH1F        *fkinkKaonPDGBkgMulti2;//!MC MC
   TH1F        *fkinkKaonPDGBkgMulti3;//!MC MC
   TH1F        *fkinkKaonPDGBkgMulti4;//!MC MC
   TH1F        *fkinkKaonPDGBkgMulti5;//!MC MC
   TH1F        *fkinkKaonPDGBkgMulti6;//!MC MC
   TH1F        *fkinkKaonPDGBkgMulti7;//!MC MC
   TH1F        *fkinkKaonPDGBkgMulti8;//!MC MC
   TH1F        *fkinkKaonPDGBkgMulti9;//!MC MC
   TH1F        *fkinkKaonPDGBkgMulti10;//!MC M

   TH1F        *fMultiBin1KaonKinksPos; //!  pos kaons in multiplicity bins
   TH1F        *fMultiBin2KaonKinksPos; //!  pos kaons in multiplicity bins
   TH1F        *fMultiBin3KaonKinksPos; //!  pos kaons in multiplicity bins
   TH1F        *fMultiBin4KaonKinksPos; //!  pos kaons in multiplicity bins
   TH1F        *fMultiBin5KaonKinksPos; //!  pos kaons in multiplicity bins
   TH1F        *fMultiBin6KaonKinksPos; //!  pos kaons in multiplicity bins
   TH1F        *fMultiBin7KaonKinksPos; //!  pos kaons in multiplicity bins
   TH1F        *fMultiBin8KaonKinksPos; //!  pos kaons in multiplicity bins
   TH1F        *fMultiBin9KaonKinksPos; //!  pos kaons in multiplicity bins
   TH1F        *fMultiBin10KaonKinksPos; //! pos kaons in multiplicity bins

   TH1F        *fMultiBin1KaonKinksNeg; //!  neg kaons in multiplicity bins
   TH1F        *fMultiBin2KaonKinksNeg; //!  neg kaons in multiplicity bins
   TH1F        *fMultiBin3KaonKinksNeg; //!  neg kaons in multiplicity bins
   TH1F        *fMultiBin4KaonKinksNeg; //!  neg kaons in multiplicity bins
   TH1F        *fMultiBin5KaonKinksNeg; //!  neg kaons in multiplicity bins
   TH1F        *fMultiBin6KaonKinksNeg; //!  neg kaons in multiplicity bins
   TH1F        *fMultiBin7KaonKinksNeg; //!  neg kaons in multiplicity bins
   TH1F        *fMultiBin8KaonKinksNeg; //!  neg kaons in multiplicity bins
   TH1F        *fMultiBin9KaonKinksNeg; //!  neg kaons in multiplicity bins
   TH1F        *fMultiBin10KaonKinksNeg; //! neg kaons in multiplicity bins

   TH2F        *fAngMomK; //! Decay angle vrs Mother Mom for pdg kaons
   TH2F        *fPosiKinkK;//!MC position  kink
   TH2F        *fPosiKinKXZ;//!MC position  kink
   TH2F        *fPosiKinKYZ;//!MC position  kink

   
   AliAnalysisKinkTaskMult13ppMC(const AliAnalysisKinkTaskMult13ppMC&); // not implemented
   AliAnalysisKinkTaskMult13ppMC& operator=(const AliAnalysisKinkTaskMult13ppMC&); // not implemented
   
   ClassDef(AliAnalysisKinkTaskMult13ppMC, 1); // example of analysis
};

#endif
