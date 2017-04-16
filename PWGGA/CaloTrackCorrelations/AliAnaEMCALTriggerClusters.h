#ifndef ALIANAEMCALTRIGGERCLUSTERS_H
#define ALIANAEMCALTRIGGERCLUSTERS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
/// \class AliAnaEMCALTriggerClusters
/// \ingroup CaloTrackCorrelationsAnalysis 
/// \brief Class for study of EMCAL trigger behaviour.
///
/// Class for study of EMCAL trigger behaviour.
///
/// More information can be found in this [twiki](https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PhotonHadronCorrelations)
/// and particularly in this [section](https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PhotonHadronCorrelations#AliAnaEMCALTriggerClusters).
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-IN2P3-CNRS
//_________________________________________________________________________

// --- ROOT system ---
class TH2F ;
class TH1F;
class TString ;
class TObjString;
class TList ;

// --- ANALYSIS system ---
#include "AliAnaCaloTrackCorrBaseClass.h"

class AliAnaEMCALTriggerClusters : public AliAnaCaloTrackCorrBaseClass {

 public: 
           AliAnaEMCALTriggerClusters() ;
    
  /// Virtual destructor.
  virtual ~AliAnaEMCALTriggerClusters() { ; }
	
  //---------------------------------------
  // General analysis frame methods
  //---------------------------------------
  
  TObjString * GetAnalysisCuts();
  
  TList      * GetCreateOutputObjects();
  
  void         Init();

  void         InitParameters();

  void         MakeAnalysisFillHistograms() ; 
  
  void         Print(const Option_t * opt)const;
  
  void         FillBadTriggerEventHistogram();
  
  void         FillRawClusterTriggerBCHistograms(Int_t idcalo,       Float_t ecluster,   Float_t tofcluster,
                                                 Float_t etacluster, Float_t phicluster);
  
  // Analysis parameters setters getters

  void         SetNCellCut(Int_t n)                { fNCellsCut = n               ; }
  Double_t     GetNCellCut()                 const { return fNCellsCut            ; }
  
  void         SetM02(Float_t min, Float_t max)    { fMinM02 = min; fMaxM02 = max ; }
  Float_t      GetM02Min()                   const { return fMinM02               ; }
  Float_t      GetM02Max()                   const { return fMaxM02               ; }
  
  Bool_t       IsTrackMatchRejectionOn()     const { return fRejectTrackMatch     ; }
  void         SwitchOnTrackMatchRejection()       { fRejectTrackMatch = kTRUE    ; }
  void         SwitchOffTrackMatchRejection()      { fRejectTrackMatch = kFALSE   ; }
  
  private:
 
  Bool_t  fRejectTrackMatch ;                         ///<  Reject clusters which have an associated TPC track
  Int_t   fNCellsCut ;                                ///<  Accept for the analysis clusters with more than fNCellsCut cells
  Float_t fMinM02  ;                                  ///<  Remove clusters with large M02
  Float_t fMaxM02  ;                                  ///<  Remove clusters with small M02
  
  TLorentzVector fMomentum ;                          //!<! Cluster momentum
  
  // Histograms
  
  TH1F * fhE               ;                          //!<! Raw clusters E
  TH1F * fhESelected       ;                          //!<! Selected custers E
  TH2F * fhEtaPhi          ;                          //!<! Raw Pseudorapidity vs Phi of clusters for E > 0.5
  TH2F * fhEtaPhiSelected  ;                          //!<! Pseudorapidity vs Phi of clusters for E > 0.5
  TH2F * fhEtaPhiEMCALBC0  ;                          //!<! Pseudorapidity vs Phi of clusters for E > 0.5
  TH2F * fhEtaPhiEMCALBC1  ;                          //!<! Pseudorapidity vs Phi of clusters for E > 0.5
  TH2F * fhEtaPhiEMCALBCN  ;                          //!<! Pseudorapidity vs Phi of clusters for E > 0.5

  TH2F * fhEtaPhiTriggerEMCALBC[11] ;                 //!<! Pseudorapidity vs Phi of clusters for E > 2
  TH2F * fhTimeTriggerEMCALBC  [11] ;                 //!<! Time distribution of clusters, when trigger is in a given BC
  TH2F * fhTimeTriggerEMCALBCPileUpSPD[11];           //!<! Time distribution of clusters, when trigger is in a given BC, tagged as pile-up SPD

  TH2F * fhEtaPhiTriggerEMCALBCUM[11] ;               //!<! Pseudorapidity vs Phi of clusters for E > 2, not matched to trigger
  TH2F * fhTimeTriggerEMCALBCUM  [11] ;               //!<! Time distribution of clusters, when trigger is in a given BC, not matched to trigger
  
  TH2F * fhEtaPhiTriggerEMCALBCCluster  [11] ;        //!<! Pseudorapidity vs Phi of trigger clusters
  TH2F * fhTimeTriggerEMCALBCCluster ;                //!<! Time distribution of clusters, when trigger cluster is in a given BC
  TH2F * fhEtaPhiTriggerEMCALBCUMCluster[11] ;        //!<! Pseudorapidity vs Phi of highest E cluster  in event, not matched to trigger
  TH2F * fhTimeTriggerEMCALBCUMCluster ;              //!<! Time distribution of highest energy cluster in event, when trigger is in a given BC, not
  
  TH2F * fhEtaPhiTriggerEMCALBCClusterOverTh     ;    //!<! Pseudorapidity vs Phi of trigger clusters, over nominal threshold
  TH2F * fhEtaPhiTriggerEMCALBCUMClusterOverTh   ;    //!<! Pseudorapidity vs Phi of highest E cluster  in event, not matched to trigger, over nominal threshold
  TH2F * fhEtaPhiTriggerEMCALBCClusterBelowTh1   ;    //!<! Pseudorapidity vs Phi of trigger clusters, 1 GeV below nominal threshold
  TH2F * fhEtaPhiTriggerEMCALBCUMClusterBelowTh1 ;    //!<! Pseudorapidity vs Phi of highest E cluster  in event, not matched to trigger, 2 GeV below nominal threshold
  TH2F * fhEtaPhiTriggerEMCALBCClusterBelowTh2   ;    //!<! Pseudorapidity vs Phi of trigger clusters, 1 GeV below nominal threshold
  TH2F * fhEtaPhiTriggerEMCALBCUMClusterBelowTh2 ;    //!<! Pseudorapidity vs Phi of highest E cluster  in event, not matched to trigger, 2 GeV below nominal threshold

  TH2F * fhEtaPhiTriggerEMCALBCExotic            ;    //!<! Pseudorapidity vs Phi of trigger exotic clusters
  TH2F * fhTimeTriggerEMCALBCExotic              ;    //!<! Time distribution of clusters, when trigger exotic cluster
  TH2F * fhEtaPhiTriggerEMCALBCUMExotic          ;    //!<! Pseudorapidity vs Phi of highest E exotic cluster  in event, not matched to trigger
  TH2F * fhTimeTriggerEMCALBCUMExotic            ;    //!<! Time distribution of highest energy exotic cluster in event, not matched to trigger

  TH2F * fhEtaPhiTriggerEMCALBCBad               ;    //!<! Pseudorapidity vs Phi of trigger exotic clusters
  TH2F * fhTimeTriggerEMCALBCBad                 ;    //!<! Time distribution of clusters, when trigger exotic
  TH2F * fhEtaPhiTriggerEMCALBCUMBad             ;    //!<! Pseudorapidity vs Phi of highest E exotic cluster  in event, not matched to trigger
  TH2F * fhTimeTriggerEMCALBCUMBad               ;    //!<! Time distribution of highest energy exotic cluster in event, not matched to trigger
  
  TH2F * fhEtaPhiTriggerEMCALBCBadExotic         ;    //!<! Pseudorapidity vs Phi of trigger exotic and bad clusters
  TH2F * fhTimeTriggerEMCALBCBadExotic           ;    //!<! Time distribution of clusters, when trigger exotic and bad cluster
  TH2F * fhEtaPhiTriggerEMCALBCUMBadExotic       ;    //!<! Pseudorapidity vs Phi of highest E exotic cluster  in event, not matched to trigger
  TH2F * fhTimeTriggerEMCALBCUMBadExotic         ;    //!<! Time distribution of highest energy exotic cluster in event, not matched to trigger
  
  TH2F * fhEtaPhiTriggerEMCALBCExoticCluster     ;    //!<! Pseudorapidity vs Phi of trigger exotic clusters
  TH2F * fhTimeTriggerEMCALBCExoticCluster       ;    //!<! Time distribution of clusters, when trigger exotic cluster
  TH2F * fhEtaPhiTriggerEMCALBCUMExoticCluster   ;    //!<! Pseudorapidity vs Phi of highest E exotic cluster  in event, not matched to trigger
  TH2F * fhTimeTriggerEMCALBCUMExoticCluster     ;    //!<! Time distribution of highest energy exotic cluster in event, not matched to trigger
  
  TH2F * fhEtaPhiTriggerEMCALBCBadCluster        ;    //!<! Pseudorapidity vs Phi of trigger bad clusters
  TH2F * fhTimeTriggerEMCALBCBadCluster          ;    //!<! Time distribution of clusters, when trigger bad cluster is in a given BC
  TH2F * fhEtaPhiTriggerEMCALBCUMBadCluster      ;    //!<! Pseudorapidity vs Phi of highest E bad cluster  in event, not matched to trigger
  TH2F * fhTimeTriggerEMCALBCUMBadCluster        ;    //!<! Time distribution of highest energy bad cluster in event, when trigger is in a given BC, not

  TH2F * fhEtaPhiTriggerEMCALBCBadExoticCluster  ;    //!<! Pseudorapidity vs Phi of trigger exotic and bad clusters
  TH2F * fhTimeTriggerEMCALBCBadExoticCluster    ;    //!<! Time distribution of clusters, when trigger exotic and bad cluster
  TH2F * fhEtaPhiTriggerEMCALBCUMBadExoticCluster;    //!<! Pseudorapidity vs Phi of highest E exotic and bad cluster in event, not matched to trigger
  TH2F * fhTimeTriggerEMCALBCUMBadExoticCluster  ;    //!<! Time distribution of highest energy exotic and bad cluster in event, not matched to trigger
  
  TH2F * fhTimeTriggerEMCALBCBadMaxCell          ;    //!<! Time distribution of trigger clusters, when trigger bad max cell
  TH2F * fhTimeTriggerEMCALBCUMBadMaxCell        ;    //!<! Time distribution of highest energy bad max cell cluster in event, when trigger is not found
  TH2F * fhTimeTriggerEMCALBCBadMaxCellExotic    ;    //!<! Time distribution of trigger clusters, when trigger exotic cluster with bad max cell
  TH2F * fhTimeTriggerEMCALBCUMBadMaxCellExotic  ;    //!<! Time distribution of highest energy exotic with bad max cell cluster in event, when trigger is not found
  
  TH2F * fhEtaPhiTriggerEMCALBCUMReMatchOpenTimeCluster ;  //!<! Pseudorapidity vs Phi of highest E bad cluster  in event, not matched to trigger, rematched open time trigger
  TH2F * fhTimeTriggerEMCALBCUMReMatchOpenTimeCluster   ;  //!<! Time distribution of highest energy bad max cell cluster in event, when trigger is not found, rematched open time trigger
  TH2F * fhEtaPhiTriggerEMCALBCUMReMatchCheckNeighCluster; //!<! Pseudorapidity vs Phi of highest E bad cluster  in event, not matched to trigger, rematched with neigbour patchs
  TH2F * fhTimeTriggerEMCALBCUMReMatchCheckNeighCluster ;  //!<! Time distribution of highest energy bad max cell cluster in event, when trigger is not found, rematched with neigbour patchs
  TH2F * fhEtaPhiTriggerEMCALBCUMReMatchBothCluster;  //!<! Pseudorapidity vs Phi of highest E bad cluster  in event, not matched to trigger, rematched open both
  TH2F * fhTimeTriggerEMCALBCUMReMatchBothCluster ;   //!<! Time distribution of highest energy bad max cell cluster in event, when trigger is not found, rematched open both
  
  TH2F * fhTimeTriggerEMCALBC0UMReMatchOpenTime   ;   //!<! Time distribution of clusters, not matched to trigger, rematched open time trigger
  TH2F * fhTimeTriggerEMCALBC0UMReMatchCheckNeigh ;   //!<! Time distribution of clusters, not matched to trigger, rematched with neighbour patchs
  TH2F * fhTimeTriggerEMCALBC0UMReMatchBoth       ;   //!<! Time distribution of clusters, not matched to trigger, rematched open both
  
  TH2F * fhEtaPhiNoTrigger ;                          //!<! Pseudorapidity vs Phi of highest E exotic cluster  in event, no trigger at all
  TH2F * fhTimeNoTrigger   ;                          //!<! Time distribution of highest energy exotic cluster in event, no trigger at all
  
  TH2F * fhEtaPhiSelectedEMCALBC0  ;                  //!<! Pseudorapidity vs Phi of identified  photon for E > 0.5
  TH2F * fhEtaPhiSelectedEMCALBC1  ;                  //!<! Pseudorapidity vs Phi of identified  photon for E > 0.5
  TH2F * fhEtaPhiSelectedEMCALBCN  ;                  //!<! Pseudorapidity vs Phi of identified  photon for E > 0.5
  TH2F * fhEtaPhiSelectedTriggerEMCALBC[11];          //!<! Pseudorapidity vs Phi of photons for E > 0.5
  TH2F * fhTimeSelectedTriggerEMCALBC  [11];          //!<! Time distribution of photons, when trigger is in a given BC
  TH2F * fhTimeSelectedTriggerEMCALBCPileUpSPD[11] ;  //!<! Time distribution of photons, when trigger is in a given BC, tagged as pile-up SPD
  TH2F * fhEtaPhiSelectedTriggerEMCALBCUM[11];        //!<! Pseudorapidity vs Phi of photons for E > 2, not matched to trigger
  TH2F * fhTimeSelectedTriggerEMCALBCUM  [11];        //!<! Time distribution of photons, when trigger is in a given BC, not matched to trigger

  TH2F * fhTimeSelectedTriggerEMCALBC0UMReMatchOpenTime   ;  //!<! Time distribution of photons in event, when trigger is not found, rematched open time trigger
  TH2F * fhTimeSelectedTriggerEMCALBC0UMReMatchCheckNeigh ;  //!<! Time distribution of photons in event, when trigger is not found, rematched with neigbour patchs
  TH2F * fhTimeSelectedTriggerEMCALBC0UMReMatchBoth       ;  //!<! Time distribution of photons in event, when trigger is not found, rematched open both
  
  /// Copy constructor not implemented.
  AliAnaEMCALTriggerClusters(              const AliAnaEMCALTriggerClusters & g) ;
    
  /// Assignment operator not implemented.
  AliAnaEMCALTriggerClusters & operator = (const AliAnaEMCALTriggerClusters & g) ;
  
  /// \cond CLASSIMP
  ClassDef(AliAnaEMCALTriggerClusters,2) ;
  /// \endcond

} ;
 
#endif//ALIANAEMCALTRIGGERCLUSTERS_H



