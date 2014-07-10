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
  
  void    SwitchOnFillPileUpHistograms()     { fFillPileUpHistograms    = kTRUE  ; }
  void    SwitchOffFillPileUpHistograms()    { fFillPileUpHistograms    = kFALSE ; }
  
  void    SwitchOnFillTrackBCHistograms()    { fFillVertexBC0Histograms = kTRUE  ; }
  void    SwitchOffFillTrackBCHistograms()   { fFillVertexBC0Histograms = kFALSE ; }

  void    SwitchOnFillVertexBC0Histograms()  { fFillVertexBC0Histograms = kTRUE  ; }
  void    SwitchOffFillVertexBC0Histograms() { fFillVertexBC0Histograms = kFALSE ; }

 private:
  
  Bool_t  fFillPileUpHistograms;    // Fill pile-up related histograms
  Bool_t  fFillTrackBCHistograms;   // Fill histograms for tracks with TOF BC=0 or not related histograms
  Bool_t  fFillVertexBC0Histograms; // Fill histograms for tracks with vertex BC=0 or not related histograms
  
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
  
  enum mvType{kmcPion = 0, kmcProton = 1, kmcKaon = 2, kmcMuon = 3, kmcElectron = 4, kmcUnknown = 4 };

  TH1F * fhPtMCPart [6];     //! pT distribution, 6 hadron ID
  TH2F * fhPhiMCPart[6];     //! phi distribution vs pT, 6 hadron ID
  TH2F * fhEtaMCPart[6];     //! eta distribution vs pT, 6 hadron ID
  
  TH1F * fhPtMCPrimPart [6]; //! Number of generated charged hadrons vs pT coming from MC particle, 6 hadron ID
  TH2F * fhPhiMCPrimPart[6]; //! Number of generated charged hadrons vs phi coming from MC particle, 6 hadron ID
  TH2F * fhEtaMCPrimPart[6]; //! Number of generated charged hadrons vs eta coming from MC particle, 6 hadron ID

  // TOF
  TH1F * fhTOFSignal;                    //! TOF signal
  TH1F * fhTOFSignalPtCut;               //! TOF signal pt and acceptance cut
  TH1F * fhTOFSignalBCOK;                //! TOF signal pt and acceptance cut
  TH2F * fhPtTOFSignal;                  //! TOF signal vs track pT, good status
  TH2F * fhPtTOFSignalDCACut;            //! TOF signal vs track pT, good status
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
  
  TH2F * fhPtNPileUpSPDVtx;	              //! cluster pt vs number of spd pile-up vertices
  TH2F * fhPtNPileUpTrkVtx;               //! cluster pt vs number of track pile-up vertices
  TH2F * fhPtNPileUpSPDVtxBC0;	          //! cluster pt vs number of spd pile-up vertices, track in BC=0
  TH2F * fhPtNPileUpTrkVtxBC0;            //! cluster pt vs number of track pile-up vertices, track in BC=0

  AliAnaChargedParticles(              const AliAnaChargedParticles & ch) ; // cpy ctor
  AliAnaChargedParticles & operator = (const AliAnaChargedParticles & ch) ; // cpy assignment
  
  ClassDef(AliAnaChargedParticles,10)

} ;


#endif //ALIANACHARGEDPARTICLES_H



