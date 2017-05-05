#ifndef ALIANACHARGEDPARTICLES_H
#define ALIANACHARGEDPARTICLES_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
/// \class AliAnaChargedParticles
/// \ingroup CaloTrackCorrelationsAnalysis 
/// \brief Track selection for correlation analysis.
///
/// Class for track selection and identification (not done now)
/// Tracks from the CTS are kept in the AOD.
/// Few histograms produced.
///
/// More information can be found in this [twiki](https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PhotonHadronCorrelations)
/// and particularly in this [section](https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PhotonHadronCorrelations#AliAnaChargedParticles).
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-IN2P3-CNRS
//_________________________________________________________________________

// Root system
class TH2F; 

// Analysis system
#include "AliAnaCaloTrackCorrBaseClass.h"
 
class AliAnaChargedParticles : public AliAnaCaloTrackCorrBaseClass {
  
 public:
    
  AliAnaChargedParticles() ;
    
  /// Virtual destructor.
  virtual ~AliAnaChargedParticles() { ; }
  
  void    FillPrimaryHistograms();
  
  TList * GetCreateOutputObjects();
    
  Bool_t  GetTrackSide  (Float_t trackEta) const { if ( trackEta >=0 ) return kTRUE ; else return kFALSE ; }
  Int_t   GetTrackSector(Float_t trackPhi) const { return Int_t (trackPhi / 0.349065850398865896) ; } // 20*pi/180
  
  void    Init();
  
  void    InitParameters();
  
  void    Print(const Option_t * opt) const;
  
  void    MakeAnalysisFillAOD()  ;
  
  void    MakeAnalysisFillHistograms() ;
  
  void    SwitchOnFillTrackBCHistograms()        { fFillTrackBCHistograms   = kTRUE  ; }
  void    SwitchOffFillTrackBCHistograms()       { fFillTrackBCHistograms   = kFALSE ; }

  void    SwitchOnFillVertexBC0Histograms()      { fFillVertexBC0Histograms = kTRUE  ; }
  void    SwitchOffFillVertexBC0Histograms()     { fFillVertexBC0Histograms = kFALSE ; }

  void    SwitchOnFillEtaPhiRegionsHistograms()  { fFillEtaPhiRegionHistograms = kTRUE  ; }
  void    SwitchOffFillEtaPhiRegionsHistograms() { fFillEtaPhiRegionHistograms = kFALSE ; }

  void    SwitchOnFillTrackMultiplicityHistograms()  { fFillTrackMultHistograms = kTRUE  ; }
  void    SwitchOffFillTrackMultiplicityHistograms() { fFillTrackMultHistograms = kFALSE ; }
  
  void    SwitchOnFillTrackDCAHistograms()       { fFillTrackDCAHistograms = kTRUE  ; }
  void    SwitchOffFillTrackDCAHistograms()      { fFillTrackDCAHistograms = kFALSE ; }
  
 private:
  
  Bool_t  fFillTrackBCHistograms;           ///<  Fill histograms for tracks with TOF BC=0 or not related histograms
  Bool_t  fFillVertexBC0Histograms;         ///<  Fill histograms for tracks with vertex BC=0 or not related histograms
  Bool_t  fFillEtaPhiRegionHistograms;      ///<  Fill track pT spectrum histograms in different eta-phi windows
  Bool_t  fFillTrackMultHistograms;         ///<  Fill track pT spectrum histograms vs track multiplicity or track sum pt
  Bool_t  fFillTrackDCAHistograms;          ///<  Fill track DCA histograms 
  
  TLorentzVector fMomentum;                 //!<! Temporary momentum container
  
  //Histograms
  TH2F * fhNTracks;                         //!<! Track multiplicity distribution per event, different pT cuts
  TH2F * fhSumPtTracks;                     //!<! Track sum pT distribution per event, different pT cuts
  
  TH2F * fhPtTrackNTracks    [10];          //!<! Track multiplicity distribution per event vs track pT, different pT cuts
  TH2F * fhPtTrackSumPtTracks[10];          //!<! Track sum pT distribution per event vs track pT, different pT cuts  
  
  TH1F * fhPt;                              //!<! pT distribution
  TH1F * fhPtNoCut;                         //!<! pT distribution, no cut
  TH1F * fhPtCutDCA;                        //!<! pT distribution, Apply DCA cut
  TH1F * fhPtCutDCABCOK;                    //!<! pT distribution, Apply DCA cut, BC=0 or -100

  TH1F * fhPtPerRegion   [18][2];           //!<! pT distribution in TPC regions
  TH1F * fhSumPtPerRegion[18][2];           //!<! pT distribution in TPC regions
  
  TH1F * fhPtNotPrimary;                    //!<! pT spectra of tracks not declared as primary (AOD)
  TH1F * fhPtNotSharedClusterCut;           //!<! pT spectra of tracks not passing the shared clusters cut (AOD)
  
  TH1F * fhPtPileUp[7];                     //!<! pT distribution, pile-up defined events
  TH2F * fhPhiNeg;                          //!<! phi distribution vs pT, negative
  TH2F * fhEtaNeg;                          //!<! eta distribution vs pT, negative
  TH2F * fhPhiPos;                          //!<! phi distribution vs pT, positive
  TH2F * fhEtaPos;                          //!<! eta distribution vs pT, positive
  TH2F * fhEtaPhiPos;                       //!<! eta vs phi distribution of positive charge
  TH2F * fhEtaPhiNeg;                       //!<! eta vs phi distribution of negative charge
  
  TH2F * fhTrackResolution;                 //!<! track resolution sigma pT vs pT, ESDs
  
  TH1F * fhPtVtxOutBC0;                     //!<! pT distribution of tracks from a vertex with BC!=0
  TH2F * fhEtaPhiVtxOutBC0;                 //!<! eta/phi distribution of tracks from a vertex with BC!=0
  TH1F * fhPtVtxInBC0;                      //!<! pT distribution of tracks from a vertex with BC=0
  TH2F * fhEtaPhiVtxInBC0;                  //!<! eta/phi distribution of tracks from a vertex with BC=0

  TH1F * fhPtSPDRefit;                      //!<! pT distribution of tracks with SPD and ITS refit
  TH1F * fhPtNoSPDRefit;                    //!<! pT distribution of constrained tracks no SPD and with ITSRefit
  TH1F * fhPtNoSPDNoRefit;                  //!<! pT distribution of constrained tracks with no SPD requierement and without ITSRefit

  TH2F * fhEtaPhiSPDRefitPt02;              //!<! eta-phi distribution of tracks with SPD and ITS refit, 0 < pT < 2 GeV
  TH2F * fhEtaPhiNoSPDRefitPt02;            //!<! eta-phi distribution of constrained tracks no SPD and with ITSRefit,  0 < pT < 2 GeV
  TH2F * fhEtaPhiNoSPDNoRefitPt02;          //!<! eta-phi distribution of constrained tracks with no SPD requierement and without ITSRefit,  0 < pT < 2 GeV

  TH2F * fhEtaPhiSPDRefitPt3;               //!<! eta-phi distribution of tracks with SPD and ITS refit, pT > 3 GeV
  TH2F * fhEtaPhiNoSPDRefitPt3;             //!<! eta-phi distribution of constrained tracks no SPD and with ITSRefit,  pT > 3 GeV
  TH2F * fhEtaPhiNoSPDNoRefitPt3;           //!<! eta-phi distribution of constrained tracks with no SPD requierement and without ITSRefit,  pT > 3 GeV

  // MC
  /// Indeces for histograms depending on the origin of the track at MC level
  enum mvType { kmcPion = 0, kmcProton = 1, kmcKaon = 2, kmcMuon = 3, kmcElectron = 4, kmcUnknown = 5 };

  TH1F * fhPtMCPart [6];                    //!<! pT distribution, 6 hadron ID
  TH2F * fhPhiMCPart[6];                    //!<! phi distribution vs pT, 6 hadron ID
  TH2F * fhEtaMCPart[6];                    //!<! eta distribution vs pT, 6 hadron ID
  
  TH1F * fhPtMCPrimPart [6];                //!<! Number of generated charged hadrons vs pT coming from MC particle, 6 hadron ID
  TH2F * fhPhiMCPrimPart[6];                //!<! Number of generated charged hadrons vs phi coming from MC particle, 6 hadron ID
  TH2F * fhEtaMCPrimPart[6];                //!<! Number of generated charged hadrons vs eta coming from MC particle, 6 hadron ID

  // TOF and BC
  TH1F * fhTOFSignal;                       //!<! TOF signal
  TH1F * fhTOFSignalPtCut;                  //!<! TOF signal pt and acceptance cut
  TH1F * fhTOFSignalBCOK;                   //!<! TOF signal pt and acceptance cut
  TH2F * fhPtTOFSignal;                     //!<! TOF signal vs track pT, good status
  TH2F * fhPtTOFSignalDCACut;               //!<! TOF signal vs track pT, good status
  TH2F * fhPtTOFSignalPileUp[7];            //!<! TOF signal vs track pT, good status, pile-up
  TH2F * fhPtTOFSignalVtxOutBC0;            //!<! TOF signal vs track pT, good status
  TH2F * fhPtTOFSignalVtxOutBC0PileUp[7];   //!<! TOF signal vs track pT, good status, pile-up
  TH2F * fhPtTOFSignalVtxInBC0;             //!<! TOF signal vs track pT, good status
  TH2F * fhPtTOFSignalVtxInBC0PileUp[7];    //!<! TOF signal vs track pT, good status, pile-up
  TH1F * fhPtTOFStatus0;                    //!<! pT of tracks not passing TOF status selection
  TH2F * fhEtaPhiTOFStatus0;                //!<! eta/phi of tracks not passing TOF status selection
  TH2F * fhEtaPhiTOFBC0;                    //!<! eta/phi of tracks passing TOF status selection, tracks in BC=0
  TH2F * fhEtaPhiTOFBCPlus;                 //!<! eta/phi of tracks passing TOF status selection, tracks in BC>0
  TH2F * fhEtaPhiTOFBCMinus;                //!<! eta/phi of tracks passing TOF status selection, tracks in BC<0
  TH2F * fhEtaPhiTOFBC0PileUpSPD;           //!<! eta/phi of tracks passing TOF status selection, tracks in BC=0, pile-up spd
  TH2F * fhEtaPhiTOFBCPlusPileUpSPD;        //!<! eta/phi of tracks passing TOF status selection, tracks in BC>0, pile-up spd
  TH2F * fhEtaPhiTOFBCMinusPileUpSPD;       //!<! eta/phi of tracks passing TOF status selection, tracks in BC<0, pile-up spd
  
  TH1F * fhProductionVertexBC;              //!<! Check BC of production vertex
  TH1F * fhProductionVertexBCPileUp[7];     //!<! Check BC of production vertex, pile-up

  // DCA
  TH2F * fhPtDCA[3];                        //!<! DCA (xy,z,constrained) of all tracks
  
  TH2F * fhPtDCASPDRefit[3];                //!<! DCA (xy,z,constrained) of tracks with SPD and ITS refit
  TH2F * fhPtDCANoSPDRefit[3];              //!<! DCA (xy,z,constrained) of constrained tracks no SPD and with ITSRefit
  TH2F * fhPtDCANoSPDNoRefit[3];            //!<! DCA (xy,z,constrained) of constrained tracks with no SPD requierement and without ITSRefit

  TH2F * fhPtDCAVtxOutBC0[3];               //!<! DCA (xy,z,constrained) of all tracks, vertex BC!=0
  TH2F * fhPtDCAVtxInBC0[3];                //!<! DCA (xy,z,constrained) of all tracks, vertex BC==0
  TH2F * fhPtDCAPileUp[3];                  //!<! DCA (xy,z,constrained) of all tracks, SPD pile-up
  TH2F * fhPtDCAVtxOutBC0PileUp[3];         //!<! DCA (xy,z,constrained) of all tracks, vertex BC!=0, SPD pile-up
  TH2F * fhPtDCAVtxInBC0PileUp[3];          //!<! DCA (xy,z,constrained) of all tracks, vertex BC==0, SPD pile-up

  TH2F * fhPtDCATOFBC0[3];                  //!<! DCA (xy,z,constrained) of all tracks, hit in TOF and BC=0
  TH2F * fhPtDCAPileUpTOFBC0[3];            //!<! DCA (xy,z,constrained) of all tracks, hit in TOF and BC=0
  TH2F * fhPtDCATOFBCOut[3];                //!<! DCA (xy,z,constrained) of all tracks, hit in TOF and BC!=0

  TH2F * fhPtDCANoTOFHit[3];                //!<! DCA (xy,z,constrained) of all tracks, no hit in TOF
  TH2F * fhPtDCAVtxOutBC0NoTOFHit[3];       //!<! DCA (xy,z,constrained) of all tracks, vertex BC!=0, no hit in TOF
  TH2F * fhPtDCAVtxInBC0NoTOFHit[3];        //!<! DCA (xy,z,constrained) of all tracks, vertex BC=0, no hit in TOF
  TH2F * fhPtDCAPileUpNoTOFHit[3];          //!<! DCA (xy,z,constrained) of all tracks, SPD pile-up, no hit in TOF
  TH2F * fhPtDCAVtxOutBC0PileUpNoTOFHit[3]; //!<! DCA (xy,z,constrained) of all tracks, vertex BC!=0, SPD pile-up, no hit in TOF
  TH2F * fhPtDCAVtxInBC0PileUpNoTOFHit[3];  //!<! DCA (xy,z,constrained) of all tracks, vertex BC=0, SPD pile-up, no hit in TOF
  
  // pile-up
  TH2F * fhPtNPileUpSPDVtx;	                //!<! cluster pt vs number of spd pile-up vertices
  TH2F * fhPtNPileUpTrkVtx;                 //!<! cluster pt vs number of track pile-up vertices
  TH2F * fhPtNPileUpSPDVtxBC0;	            //!<! cluster pt vs number of spd pile-up vertices, track in BC=0
  TH2F * fhPtNPileUpTrkVtxBC0;              //!<! cluster pt vs number of track pile-up vertices, track in BC=0

  /// Copy constructor not implemented.
  AliAnaChargedParticles(              const AliAnaChargedParticles & ch) ;

  /// Assignment operator not implemented.
  AliAnaChargedParticles & operator = (const AliAnaChargedParticles & ch) ;
  
  /// \cond CLASSIMP
  ClassDef(AliAnaChargedParticles,12) ;
  /// \endcond

} ;


#endif //ALIANACHARGEDPARTICLES_H



