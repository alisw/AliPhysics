#ifndef ALIANACHARGEDPARTICLES_H
#define ALIANACHARGEDPARTICLES_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
//
// Class for track selection and identification (not done now)
// Tracks from the CTS are kept in the AOD.
// Few histograms produced.
//
//-- Author: Gustavo Conesa (INFN-LNF)

// Root system
class TH2F; 

// Analysis system
#include "AliAnaCaloTrackCorrBaseClass.h"
 
class AliAnaChargedParticles : public AliAnaCaloTrackCorrBaseClass {
  
 public: 
  AliAnaChargedParticles() ; // default ctor
  virtual ~AliAnaChargedParticles() { ; } //virtual dtor

  TList * GetCreateOutputObjects();
  
  void    Init();
  
  void    InitParameters();
  
  void    Print(const Option_t * opt) const;
  
  void    MakeAnalysisFillAOD()  ;
  
  void    MakeAnalysisFillHistograms() ; 
  
  Int_t   GetPdgOfSelectedCharged()  const  { return fPdg ; }
  void    SelectChargedWithPdg( Int_t pdg ) { fPdg = pdg  ; }
  
  void    SwitchOnFillPileUpHistograms()    { fFillPileUpHistograms = kTRUE  ; }
  void    SwitchOffFillPileUpHistograms()   { fFillPileUpHistograms = kFALSE ; }
  
 private:
  
  Int_t  fPdg ;                  // identified particle id
  Bool_t fFillPileUpHistograms;  // Fill pile-up related histograms

  //Histograms
  TH1F * fhNtracks;     //! track multiplicity distribution
  TH1F * fhPt;          //! pT distribution
  TH1F * fhPtPileUp[7]; //! pT distribution, pile-up defined events
  TH2F * fhPhiNeg;      //! phi distribution vs pT, negative
  TH2F * fhEtaNeg;      //! eta distribution vs pT, negative
  TH2F * fhPhiPos;      //! phi distribution vs pT, positive
  TH2F * fhEtaPos;      //! eta distribution vs pT, positive
  TH2F * fhEtaPhiPos;   //! eta vs phi distribution of positive charge  
  TH2F * fhEtaPhiNeg;   //! eta vs phi distribution of negative charge
  
  //MC
  TH1F * fhPtPion;      //! pT distribution
  TH2F * fhPhiPion;     //! phi distribution vs pT
  TH2F * fhEtaPion;     //! eta distribution vs pT
  
  TH1F * fhPtProton;    //! pT distribution
  TH2F * fhPhiProton;   //! phi distribution vs pT
  TH2F * fhEtaProton;   //! eta distribution vs pT
  
  TH1F * fhPtElectron;  //! pT distribution
  TH2F * fhPhiElectron; //! phi distribution vs pT
  TH2F * fhEtaElectron; //! eta distribution vs pT
  
  TH1F * fhPtKaon;      //! pT distribution
  TH2F * fhPhiKaon;     //! phi distribution vs pT
  TH2F * fhEtaKaon;     //! eta distribution vs pT
  
  TH1F * fhPtUnknown;   //! pT distribution
  TH2F * fhPhiUnknown;  //! phi distribution vs pT
  TH2F * fhEtaUnknown;  //! eta distribution vs pT
  
  // TOF
  TH1F * fhTOFSignal;                    //! TOF signal
  TH1F * fhTOFSignalPtCut;               //! TOF signal pt and acceptance cut
  TH1F * fhTOFSignalBCOK;                //! TOF signal pt and acceptance cut
  TH2F * fhPtTOFSignal;                  //! TOF signal vs track pT, good status
  TH2F * fhPtTOFSignalPileUp[7];         //! TOF signal vs track pT, good status, pile-up
  TH1F * fhPtTOFStatus0;                 //! pT of tracks not passing TOF status selection
  TH2F * fhEtaPhiTOFStatus0;             //! eta/phi of tracks not passing TOF status selection
  TH2F * fhEtaPhiTOFBC0;                 //! eta/phi of tracks passing TOF status selection, tracks in BC=0
  TH2F * fhEtaPhiTOFBCPlus;              //! eta/phi of tracks passing TOF status selection, tracks in BC>0
  TH2F * fhEtaPhiTOFBCMinus;             //! eta/phi of tracks passing TOF status selection, tracks in BC<0
  TH2F * fhEtaPhiTOFBC0PileUpSPD;        //! eta/phi of tracks passing TOF status selection, tracks in BC=0, pile-up spd
  TH2F * fhEtaPhiTOFBCPlusPileUpSPD;     //! eta/phi of tracks passing TOF status selection, tracks in BC>0, pile-up spd
  TH2F * fhEtaPhiTOFBCMinusPileUpSPD;    //! eta/phi of tracks passing TOF status selection, tracks in BC<0, pile-up spd
//  TH1F * fhProductionVertexBC;           //!  check BC of production vertex

  TH2F * fhPtDCA[3];                     //! DCA (xy,z,constrained) of all tracks
  //TH2F * fhPtDCAVtxOutBC0[3];            //! DCA (xy,z,constrained) of all tracks, vertex BC!=0
  TH2F * fhPtDCAPileUp[3];               //! DCA (xy,z,constrained) of all tracks, SPD pile-up
  //TH2F * fhPtDCAVtxOutBC0PileUp[3];      //! DCA (xy,z,constrained) of all tracks, vertex BC!=0, SPD pile-up

  TH2F * fhPtDCATOFBC0[3];               //! DCA (xy,z,constrained) of all tracks, hit in TOF and BC=0
  TH2F * fhPtDCAPileUpTOFBC0[3];         //! DCA (xy,z,constrained) of all tracks, hit in TOF and BC=0

  TH2F * fhPtDCANoTOFHit[3];                //! DCA (xy,z,constrained) of all tracks, no hit in TOF
  //TH2F * fhPtDCAVtxOutBC0NoTOFHit[3];       //! DCA (xy,z,constrained) of all tracks, vertex BC!=0, no hit in TOF
  TH2F * fhPtDCAPileUpNoTOFHit[3];          //! DCA (xy,z,constrained) of all tracks, SPD pile-up, no hit in TOF
  //TH2F * fhPtDCAVtxOutBC0PileUpNoTOFHit[3]; //! DCA (xy,z,constrained) of all tracks, vertex BC!=0, SPD pile-up, no hit in TOF
  
  
  AliAnaChargedParticles(              const AliAnaChargedParticles & ch) ; // cpy ctor
  AliAnaChargedParticles & operator = (const AliAnaChargedParticles & ch) ; // cpy assignment
  
  ClassDef(AliAnaChargedParticles,5)

} ;


#endif //ALIANACHARGEDPARTICLES_H



