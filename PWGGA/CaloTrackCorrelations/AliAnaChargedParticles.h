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

  Bool_t  AcceptDCA(const Float_t pt, const Float_t dca) ;
  
  TList * GetCreateOutputObjects();
  
  Int_t   GetVertexBC(const AliVVertex * vtx);
  
  void    Init();
  
  void    InitParameters();
  
  void    Print(const Option_t * opt) const;
  
  void    MakeAnalysisFillAOD()  ;
  
  void    MakeAnalysisFillHistograms() ; 
  
  void    SwitchOnFillPileUpHistograms()     { fFillPileUpHistograms    = kTRUE  ; }
  void    SwitchOffFillPileUpHistograms()    { fFillPileUpHistograms    = kFALSE ; }
  
  void    SwitchOnFillVertexBC0Histograms()  { fFillVertexBC0Histograms = kTRUE  ; }
  void    SwitchOffFillVertexBC0Histograms() { fFillVertexBC0Histograms = kFALSE ; }

  void    SwitchOnRecalculateVertexBC()      { fRecalculateVertexBC     = kTRUE  ; }
  void    SwitchOffRecalculateVertexBC()     { fRecalculateVertexBC     = kFALSE ; }

  void    SetDCACutParameters(Int_t i, Float_t par) { if(i >= 0 && i < 3) fDCACutParam[i] = par ; }
  
 private:
  
  Bool_t  fFillPileUpHistograms;    // Fill pile-up related histograms
  Bool_t  fFillVertexBC0Histograms; // Fill histograms for tracks with vertex BC=0 or not related histograms
  Bool_t  fRecalculateVertexBC;     // Recalculate vertex BC for older AODs
  Float_t fDCACutParam[3];          // DCA cut function parameters
  
  //Histograms
  TH1F * fhNtracks;     //! track multiplicity distribution
  TH1F * fhPt;          //! pT distribution
  TH1F * fhPtNoCut;     //! pT distribution, no cut
  TH1F * fhPtCutDCA;    //! pT distribution, Apply DCA cut
  TH1F * fhPtCutDCABCOK;//! pT distribution, Apply DCA cut, BC=0 or -100

  TH1F * fhPtPileUp[7]; //! pT distribution, pile-up defined events
  TH2F * fhPhiNeg;      //! phi distribution vs pT, negative
  TH2F * fhEtaNeg;      //! eta distribution vs pT, negative
  TH2F * fhPhiPos;      //! phi distribution vs pT, positive
  TH2F * fhEtaPos;      //! eta distribution vs pT, positive
  TH2F * fhEtaPhiPos;   //! eta vs phi distribution of positive charge  
  TH2F * fhEtaPhiNeg;   //! eta vs phi distribution of negative charge
  
  TH1F * fhPtVtxOutBC0;    //! pT distribution of tracks from a vertex with BC!=0
  TH2F * fhEtaPhiVtxOutBC0;//! eta/phi distribution of tracks from a vertex with BC!=0
  TH1F * fhPtVtxInBC0;     //! pT distribution of tracks from a vertex with BC=0
  TH2F * fhEtaPhiVtxInBC0; //! eta/phi distribution of tracks from a vertex with BC=0

  TH1F * fhPtSPDRefit;     //! pT distribution of tracks with SPD and ITS refit
  TH1F * fhPtNoSPDRefit;   //! pT distribution of constrained tracks no SPD and with ITSRefit
  TH1F * fhPtNoSPDNoRefit; //! pT distribution of constrained tracks with no SPD requierement and without ITSRefit

  TH2F * fhEtaPhiSPDRefitPt02;     //! eta-phi distribution of tracks with SPD and ITS refit, 0 < pT < 2 GeV
  TH2F * fhEtaPhiNoSPDRefitPt02;   //! eta-phi distribution of constrained tracks no SPD and with ITSRefit,  0 < pT < 2 GeV
  TH2F * fhEtaPhiNoSPDNoRefitPt02; //! eta-phi distribution of constrained tracks with no SPD requierement and without ITSRefit,  0 < pT < 2 GeV

  TH2F * fhEtaPhiSPDRefitPt3;     //! eta-phi distribution of tracks with SPD and ITS refit, pT > 3 GeV
  TH2F * fhEtaPhiNoSPDRefitPt3;   //! eta-phi distribution of constrained tracks no SPD and with ITSRefit,  pT > 3 GeV
  TH2F * fhEtaPhiNoSPDNoRefitPt3; //! eta-phi distribution of constrained tracks with no SPD requierement and without ITSRefit,  pT > 3 GeV

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
  TH2F * fhPtTOFSignalVtxOutBC0;         //! TOF signal vs track pT, good status
  TH2F * fhPtTOFSignalVtxOutBC0PileUp[7];//! TOF signal vs track pT, good status, pile-up
  TH2F * fhPtTOFSignalVtxInBC0;          //! TOF signal vs track pT, good status
  TH2F * fhPtTOFSignalVtxInBC0PileUp[7]; //! TOF signal vs track pT, good status, pile-up
  TH1F * fhPtTOFStatus0;                 //! pT of tracks not passing TOF status selection
  TH2F * fhEtaPhiTOFStatus0;             //! eta/phi of tracks not passing TOF status selection
  TH2F * fhEtaPhiTOFBC0;                 //! eta/phi of tracks passing TOF status selection, tracks in BC=0
  TH2F * fhEtaPhiTOFBCPlus;              //! eta/phi of tracks passing TOF status selection, tracks in BC>0
  TH2F * fhEtaPhiTOFBCMinus;             //! eta/phi of tracks passing TOF status selection, tracks in BC<0
  TH2F * fhEtaPhiTOFBC0PileUpSPD;        //! eta/phi of tracks passing TOF status selection, tracks in BC=0, pile-up spd
  TH2F * fhEtaPhiTOFBCPlusPileUpSPD;     //! eta/phi of tracks passing TOF status selection, tracks in BC>0, pile-up spd
  TH2F * fhEtaPhiTOFBCMinusPileUpSPD;    //! eta/phi of tracks passing TOF status selection, tracks in BC<0, pile-up spd
  
  TH1F * fhProductionVertexBC;           //!  check BC of production vertex
  TH1F * fhProductionVertexBCPileUp[7];  //!  check BC of production vertex, pile-up

  TH2F * fhPtDCA[3];                     //! DCA (xy,z,constrained) of all tracks
  
  TH2F * fhPtDCASPDRefit[3];             //! DCA (xy,z,constrained) of tracks with SPD and ITS refit
  TH2F * fhPtDCANoSPDRefit[3];           //! DCA (xy,z,constrained) of constrained tracks no SPD and with ITSRefit
  TH2F * fhPtDCANoSPDNoRefit[3];         //! DCA (xy,z,constrained) of constrained tracks with no SPD requierement and without ITSRefit

  TH2F * fhPtDCAVtxOutBC0[3];            //! DCA (xy,z,constrained) of all tracks, vertex BC!=0
  TH2F * fhPtDCAVtxInBC0[3];             //! DCA (xy,z,constrained) of all tracks, vertex BC==0
  TH2F * fhPtDCAPileUp[3];               //! DCA (xy,z,constrained) of all tracks, SPD pile-up
  TH2F * fhPtDCAVtxOutBC0PileUp[3];      //! DCA (xy,z,constrained) of all tracks, vertex BC!=0, SPD pile-up
  TH2F * fhPtDCAVtxInBC0PileUp[3];       //! DCA (xy,z,constrained) of all tracks, vertex BC==0, SPD pile-up

  TH2F * fhPtDCATOFBC0[3];               //! DCA (xy,z,constrained) of all tracks, hit in TOF and BC=0
  TH2F * fhPtDCAPileUpTOFBC0[3];         //! DCA (xy,z,constrained) of all tracks, hit in TOF and BC=0
  TH2F * fhPtDCATOFBCOut[3];             //! DCA (xy,z,constrained) of all tracks, hit in TOF and BC!=0

  TH2F * fhPtDCANoTOFHit[3];                //! DCA (xy,z,constrained) of all tracks, no hit in TOF
  TH2F * fhPtDCAVtxOutBC0NoTOFHit[3];       //! DCA (xy,z,constrained) of all tracks, vertex BC!=0, no hit in TOF
  TH2F * fhPtDCAVtxInBC0NoTOFHit[3];        //! DCA (xy,z,constrained) of all tracks, vertex BC=0, no hit in TOF
  TH2F * fhPtDCAPileUpNoTOFHit[3];          //! DCA (xy,z,constrained) of all tracks, SPD pile-up, no hit in TOF
  TH2F * fhPtDCAVtxOutBC0PileUpNoTOFHit[3]; //! DCA (xy,z,constrained) of all tracks, vertex BC!=0, SPD pile-up, no hit in TOF
  TH2F * fhPtDCAVtxInBC0PileUpNoTOFHit[3];  //! DCA (xy,z,constrained) of all tracks, vertex BC=0, SPD pile-up, no hit in TOF
  
  
  AliAnaChargedParticles(              const AliAnaChargedParticles & ch) ; // cpy ctor
  AliAnaChargedParticles & operator = (const AliAnaChargedParticles & ch) ; // cpy assignment
  
  ClassDef(AliAnaChargedParticles,7)

} ;


#endif //ALIANACHARGEDPARTICLES_H



