#ifndef ALIANACLUSTERPILEUP_H
#define ALIANACLUSTERPILEUP_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
/// \class AliAnaClusterPileUp
/// \ingroup CaloTrackCorrelationsAnalysis 
/// \brief Fill histograms for cluster spectra dependence on pile-up.
///
/// Class for the study of Pile-up effect on
/// Calorimeter clusters.
/// Open time cuts in reader.
///
/// More information can be found in this [twiki](https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PhotonHadronCorrelations)
/// and particularly in this [section](https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PhotonHadronCorrelations#AliAnaClusterPileUp).
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-IN2P3-CNRS
//_________________________________________________________________________


// --- ROOT system ---
class TH2F ;
class TH1F;
class TObjString;
class TList ;

// --- ANALYSIS system ---
#include "AliAnaCaloTrackCorrBaseClass.h"

class AliAnaClusterPileUp : public AliAnaCaloTrackCorrBaseClass {

 public: 
           AliAnaClusterPileUp() ;
    
  /// Virtual destructor.
  virtual ~AliAnaClusterPileUp() { ; }
	
  //---------------------------------------
  // General analysis frame methods
  //---------------------------------------
  
  TObjString * GetAnalysisCuts();
  
  TList      * GetCreateOutputObjects();
  
  void         Init();

  void         InitParameters();

  void         MakeAnalysisFillHistograms() ; 
  
  void         Print(const Option_t * opt)const;

  // Analysis parameters setters getters
  
  void         SetNCellCut(Int_t n)                   { fNCellsCut = n             ; }
  Double_t     GetNCellCut()                    const { return fNCellsCut          ; }
  
private:
 
  Int_t   fNCellsCut ;                              ///<  Accept for the analysis clusters with more than fNCellsCut cells

  TLorentzVector fMomentum;                         //!<! Cluster momentum
  
  //Histograms
  
  TH1F * fhPtPileUp[7];                             //!<! pT distribution of clusters before any selection
  TH1F * fhPtNeutralPileUp[7];                      //!<! pT distribution of track matched clusters
  TH2F * fhLambda0PileUp[7];                        //!<! E vs M02 distribution of clusters, before any selection
  TH2F * fhLambda0NeutralPileUp[7];                 //!<! E vs M02 distribution of clusters, track matched clusters
  TH2F * fhClusterCellTimePileUp[7];                //!<! E vs Time inside cluster, before any selection, not max cell
  TH2F * fhClusterTimeDiffPileUp[7];                //!<! E vs Time difference inside cluster, before any selection
  TH2F * fhClusterTimeDiffNeutralPileUp[7];         //!<! E vs Time difference inside cluster for track matched clusters
  TH2F * fhClusterEFracLongTimePileUp[7];           //!<! E vs fraction of cluster energy from cells with large time
  TH2F * fhTimePtNoCut;                             //!<! Time of cluster vs Pt, no cut
  TH2F * fhTimePtSPD;                               //!<! Time of cluster vs Pt, IsSPDPileUp
  TH2F * fhTimeNPileUpVertSPD;                      //!<! Time of cluster vs n pile-up vertices from SPD
  TH2F * fhTimeNPileUpVertTrack;                    //!<! Time of cluster vs n pile-up vertices from Tracks
  TH2F * fhTimeNPileUpVertContributors;             //!<! Time of cluster vs n pile-up vertex from SPD contributors
  TH2F * fhTimePileUpMainVertexZDistance;           //!<! Time of cluster vs difference of z main vertex and pile-up vertex
  TH2F * fhTimePileUpMainVertexZDiamond;            //!<! Time of cluster vs difference of z diamond and pile-up vertex
  TH2F * fhClusterMultSPDPileUp[4];                 //!<! E max cluster vs event cluster multiplicity, for tmax-tdiff cuts, pile up event
  TH2F * fhClusterMultNoPileUp[4];                  //!<! E max cluster vs event cluster multiplicity, for tmax-tdiff cuts, not pile up event
  TH2F * fhEtaPhiBC0;                               //!<! eta/phi of clusters in BC=0
  TH2F * fhEtaPhiBCPlus;                            //!<! eta/phi of clusters in BC>0
  TH2F * fhEtaPhiBCMinus;                           //!<! eta/phi of clusters in BC<0
  TH2F * fhEtaPhiBC0PileUpSPD;                      //!<! eta/phi of clusters in BC=0, SPD pile-up
  TH2F * fhEtaPhiBCPlusPileUpSPD;                   //!<! eta/phi of clusters in BC>0, SPD pile-up
  TH2F * fhEtaPhiBCMinusPileUpSPD;                  //!<! eta/phi of clusters in BC<0, SPD pile-up

  TH2F * fhPtNPileUpSPDVtx;	                        //!<! Cluster pt vs number of spd pile-up vertices
  TH2F * fhPtNPileUpTrkVtx;                         //!<! Cluster pt vs number of track pile-up vertices
  TH2F * fhPtNPileUpSPDVtxTimeCut;	                //!<! Cluster pt vs number of spd pile-up vertices, time cut +-25 ns
  TH2F * fhPtNPileUpTrkVtxTimeCut;                  //!<! Cluster pt vs number of track pile-up vertices, time cut +- 25 ns
  TH2F * fhPtNPileUpSPDVtxTimeCut2;	                //!<! Cluster pt vs number of spd pile-up vertices, time cut +-75 ns
  TH2F * fhPtNPileUpTrkVtxTimeCut2;                 //!<! Cluster pt vs number of track pile-up vertices, time cut +- 75 ns
	
  /// Copy constructor not implemented.
  AliAnaClusterPileUp(              const AliAnaClusterPileUp & pu) ;
    
  /// Assignment operator not implemented.
  AliAnaClusterPileUp & operator = (const AliAnaClusterPileUp & pu) ;
  
  /// \cond CLASSIMP
  ClassDef(AliAnaClusterPileUp,2) ;
  /// \endcond

} ;
 
#endif//ALIANACLUSTERPILEUP_H



