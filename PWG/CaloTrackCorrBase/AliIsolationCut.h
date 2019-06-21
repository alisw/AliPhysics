#ifndef ALIISOLATIONCUT_H
#define ALIISOLATIONCUT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
/// \class AliIsolationCut
/// \ingroup CaloTrackCorrelationsBase
/// \brief Class with utils to perform Isolation Cuts.
///
/// Class containing methods for the isolation cut.
/// An AOD candidate (AliCaloTrackParticleCorrelation type)
/// is passed. Look in a cone around the candidate and study
/// the hadronic activity inside to decide if the candidate is isolated
///
/// More information can be found in this [twiki](https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PhotonHadronCorrelations).
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-IN2P3-CNRS
//_________________________________________________________________________

// --- ROOT system ---
#include <TObject.h>
class TObjArray ;
class TList   ;
#include <TLorentzVector.h>

// --- ANALYSIS system ---
class AliCaloTrackParticleCorrelation ;
class AliCaloTrackReader ;
class AliCaloPID;
class AliHistogramRanges;

class AliIsolationCut : public TObject {

 public:

  AliIsolationCut() ;  // default ctor

  /// Virtual destructor.
  virtual ~AliIsolationCut() { ; }

  // Enums

  enum type       { kPtThresIC, kSumPtIC, kPtFracIC, kSumPtFracIC, kSumDensityIC, kSumBkgSubIC } ;

  enum partInCone { kNeutralAndCharged=0, kOnlyNeutral=1, kOnlyCharged=2 } ;

  // Main Methods
  
  void       InitParameters() ;

  TString    GetICParametersList() ;

  Float_t    GetCellDensity(  AliCaloTrackParticleCorrelation * pCandidate,
                              AliCaloTrackReader * reader) const ;

  TList *    GetCreateOutputObjects();  
  
  void       MakeIsolationCut(AliCaloTrackParticleCorrelation  * pCandidate, AliCaloTrackReader * reader,
                              Bool_t bFillAOD, Bool_t useRefs, TString aodObjArrayName,
                              TObjArray * bgTrk, TObjArray * bgCls,
                              Int_t calorimeter, AliCaloPID * pid,
                              Int_t &n, Int_t & nfrac, Float_t &ptSum, Float_t &ptLead, Bool_t & isolated, 
                              Double_t histoWeight = 1) ;

  void       Print(const Option_t * opt) const ;

  Float_t    Radius(Float_t etaCandidate, Float_t phiCandidate, Float_t eta, Float_t phi) const ;

  // Cone content calculation
  
  void       CalculateCaloSignalInCone (AliCaloTrackParticleCorrelation * aodParticle, AliCaloTrackReader * reader, 
                                        Bool_t    bFillAOD    , Bool_t    useRefs, 
                                        TString   refArrayName, TObjArray *bgTrk,
                                        Int_t     calorimeter , AliCaloPID * pid,
                                        Int_t   & nPart       , Int_t   & nfrac,
                                        Float_t & coneptsum   , Float_t & coneptLead,
                                        Float_t & etaBandPtSum, Float_t & phiBandPtSum,  
                                        Double_t histoWeight = 1) ;
  
  void       CalculateTrackSignalInCone(AliCaloTrackParticleCorrelation * aodParticle, AliCaloTrackReader * reader,
                                        Bool_t    bFillAOD    , Bool_t    useRefs, 
                                        TString   refArrayName, TObjArray *bgCls,
                                        Int_t   & nPart       , Int_t   & nfrac,
                                        Float_t & coneptsum   , Float_t & coneptLead,  
                                        Float_t & etaBandPtSum, Float_t & phiBandPtSum, 
                                        Float_t & perpBandPtSum,Double_t  histoWeight = 1) ;
  
  // Cone background studies medthods

  Float_t    CalculateExcessAreaFraction(Float_t excess) const ;

  void       CalculateUEBandClusterNormalization(AliCaloTrackReader * reader,     Float_t   etaC, Float_t phiC,
                                                 Float_t   phiUEptsumCluster,     Float_t   etaUEptsumCluster,
                                                 Float_t & phiUEptsumClusterNorm, Float_t & etaUEptsumClusterNorm,
                                                 Float_t & excessFracPhi,         Float_t & excessFracEta         ) const ;

  void       CalculateUEBandTrackNormalization  (AliCaloTrackReader * reader,     Float_t   etaC,  Float_t phiC,
                                                 Float_t   phiUEptsumTrack,       Float_t   etaUEptsumTrack    , 
                                                 Float_t & phiUEptsumTrackNorm,   Float_t & etaUEptsumTrackNorm,
                                                 Float_t & excessFracPhi,         Float_t & excessFracEta        ) const ;

  void 	     GetCoeffNormBadCell(AliCaloTrackParticleCorrelation * pCandidate,
                                 AliCaloTrackReader * reader,
                                 Float_t & coneBadCellsCoeff,
                                 Float_t & etaBandBadCellsCoeff  , Float_t & phiBandBadCellsCoeff) ;


  // Parameter setters and getters

  Float_t    GetConeSize()            const { return fConeSize       ; }
  Float_t    GetPtThreshold()         const { return fPtThreshold    ; }
  Float_t    GetPtThresholdMax()      const { return fPtThresholdMax ; }
  Float_t    GetSumPtThreshold()      const { return fSumPtThreshold ; }
  Float_t    GetSumPtThresholdMax()   const { return fSumPtThresholdMax ; }
  Float_t    GetPtFraction()          const { return fPtFraction     ; }
  Int_t      GetICMethod()            const { return fICMethod       ; }
  Int_t      GetParticleTypeInCone()  const { return fPartInCone     ; }
  Int_t      GetDebug()               const { return fDebug          ; }
  Bool_t     GetFracIsThresh()        const { return fFracIsThresh   ; }
  Float_t    GetMinDistToTrigger()    const { return fDistMinToTrigger ; }

  void       SetConeSize(Float_t r)                            { fConeSize          = r    ; }
  void       SetPtThreshold(Float_t pt)                        { fPtThreshold       = pt   ; }
  void       SetPtThresholdMax(Float_t pt)                     { fPtThresholdMax    = pt   ; }
  void       SetSumPtThreshold(Float_t s)                      { fSumPtThreshold    = s    ; }
  void       SetSumPtThresholdMax(Float_t s)                   { fSumPtThresholdMax = s    ; }
  void       SetPtFraction(Float_t pt)                         { fPtFraction        = pt   ; }
  void       SetICMethod(Int_t i )                             { fICMethod          = i    ; }
  void       SetParticleTypeInCone(Int_t i)                    { fPartInCone        = i    ; }
  void       SetDebug(Int_t d)                                 { fDebug             = d    ; }
  void       SetFracIsThresh(Bool_t f )                        { fFracIsThresh      = f    ; }
  void       SetTrackMatchedClusterRejectionInCone(Bool_t tm)  { fIsTMClusterInConeRejected = tm ; }
  void       SetMinDistToTrigger(Float_t md)                   { fDistMinToTrigger  = md   ; }
  void       SetHistogramRanges(AliHistogramRanges * range)    { fHistoRanges       = range; } 
  
  void       SwitchOnFillEtaPhiHistograms ()                   { fFillEtaPhiHistograms = kTRUE  ;}
  void       SwitchOffFillEtaPhiHistograms()                   { fFillEtaPhiHistograms = kFALSE ;}
  
 private:

  Bool_t     fFillHistograms;    ///< Fill histograms if GetCreateOuputObjects() was called. 
  Bool_t     fFillEtaPhiHistograms; ///< Fill histograms if GetCreateOuputObjects() was called with eta/phi or band related histograms 

  Float_t    fConeSize ;         ///< Size of the isolation cone

  Float_t    fPtThreshold ;      ///< Minimum pt of the particles in the cone or sum in cone (UE pt mean in the forward region cone)

  Float_t    fPtThresholdMax ;   ///< Maximum pt of the particles outside the cone (needed to fit shower distribution isolated/non-isolated particles)

  Float_t    fSumPtThreshold ;   ///< Minimum of sum pt of the particles in the cone (UE sum in the forward region cone)

  Float_t    fSumPtThresholdMax ;///< Maximum of sum pt of the particles in the cone (UE sum in the forward region cone)

  Float_t    fPtFraction ;       ///< Fraction of the momentum of particles in cone or sum in cone.

  Int_t      fICMethod ;         ///< Isolation cut method to be used: kPtIC, kSumPtIC, kPtFracIC, kSumPtFracIC.

  Int_t      fPartInCone;        ///< Type of particles inside cone: kNeutralAndCharged, kOnlyNeutral, kOnlyCharged.

  Bool_t     fFracIsThresh;      ///< Use threshold instead of fraction when pt leading is small.

  Bool_t     fIsTMClusterInConeRejected; ///< Enable to remove the Track matching removal of clusters in cone sum pt calculation in case of kNeutralAndCharged analysis
  
  Float_t    fDistMinToTrigger;  ///<  Minimal distance between isolation candidate particle and particles in cone to count them for this isolation. Do not count in cone particles close to the trigger.
  
  Int_t      fDebug;             ///< Debug level.

  TLorentzVector fMomentum;      //!<! Momentum of cluster, temporal object.

  TVector3   fTrackVector;       //!<! Track moment, temporal object.
  
  // Histograms
  
  AliHistogramRanges * fHistoRanges;                   ///!  Histogram bins and ranges  data-base
  
  TH2F *   fhPtInCone ;                                //!<! Cluster/track Pt in the cone.
  TH2F *   fhPtClusterInCone ;                         //!<! Cluster Pt in the cone.
  TH2F *   fhPtTrackInCone ;                           //!<! Track Pt in the cone.
  TH2F *   fhPtInPerpCone ;                            //!<! Particle Pt  in cone at the perpendicular phi region to trigger axis  (phi +90).
  
  TH2F *   fhConeSumPt ;                               //!<! Cluster and tracks Sum Pt in the cone.
  TH2F *   fhConeSumPtCluster ;                        //!<! Clusters Sum Pt in the cone.
  TH2F *   fhConeSumPtTrack ;                          //!<! Tracks Sum Pt in the cone.
  TH2F *   fhPerpConeSumPt ;                           //!<! Sum Pt in cone at the perpendicular phi region to trigger axis  (phi +90).

  TH2F *   fhConePtLead ;                              //!<! Cluster and tracks leading pt in the cone.
  TH2F *   fhConePtLeadCluster ;                       //!<! Clusters leading pt in the cone.
  TH2F *   fhConePtLeadTrack ;                         //!<! Tracks leading pt in the cone.
  
  TH2F *   fhConeSumPtClustervsTrack ;                 //!<! Cluster vs tracks Sum Pt Sum Pt in the cone.
  
  TH2F *   fhConeSumPtClusterTrackFrac ;               //!<! Cluster / tracks Sum Pt Sum Pt in the cone.
  
  TH2F *   fhConePtLeadClustervsTrack;                 //!<! Tracks vs Clusters leading pt.
  
  TH2F *   fhConePtLeadClusterTrackFrac;               //!<! Trigger pt vs cluster/track leading pt.
  
  TH2F *   fhConeSumPtTrigEtaPhi ;                     //!<! Cluster and tracks Sum Pt Sum Pt in the cone, per eta-phi bin of trigger.
  
  TH2F *   fhEtaPhiCluster ;                           //!<! Eta vs. phi of all clusters.
  TH2F *   fhEtaPhiTrack ;                             //!<! Eta vs. phi of all tracks.
  TH2F *   fhEtaPhiInConeCluster ;                     //!<! Eta vs. phi of clusters in cone.
  TH2F *   fhEtaPhiInConeTrack ;                       //!<! Eta vs. phi of tracks in cone.
  TH2F *   fhEtaPhiInPerpCone ;                        //!<! Eta vs. phi of tracks in perpendicular cone

  // UE bands
  
  TH2F *   fhEtaBandClusterPt ;                        //!<! pT in Eta band to estimate UE in cone, only clusters.
  TH2F *   fhPhiBandClusterPt ;                        //!<! pT in Phi band to estimate UE in cone, only clusters.
  TH2F *   fhEtaBandTrackPt   ;                        //!<! pT in Eta band to estimate UE in cone, only tracks.
  TH2F *   fhPhiBandTrackPt   ;                        //!<! pT in Phi band to estimate UE in cone, only tracks.
  
  TH2F *   fhConeSumPtEtaBandUECluster;                //!<! Cluster Sum Pt in the eta band for clusters, before normalization.
  TH2F *   fhConeSumPtPhiBandUECluster;                //!<! Cluster Sum Pt in the phi band for clusters, before normalization.
  TH2F *   fhConeSumPtEtaBandUETrack;                  //!<! Track Sum Pt in the eta band for tracks, before normalization.
  TH2F *   fhConeSumPtPhiBandUETrack;                  //!<! Track Sum Pt in the phi badn for tracks, before normalization.
  
  TH2F *   fhEtaBandClusterEtaPhi ;                    //!<! Eta vs Phi in Eta band to estimate UE in cone, only clusters. 
  TH2F *   fhPhiBandClusterEtaPhi ;                    //!<! Eta vs Phi in Phi band to estimate UE in cone, only clusters.
  TH2F *   fhEtaBandTrackEtaPhi   ;                    //!<! Eta vs Phi in Eta band to estimate UE in cone, only tracks.
  TH2F *   fhPhiBandTrackEtaPhi   ;                    //!<! Eta vs Phi in Phi band to estimate UE in cone, only tracks. 
  
  TH2F *   fhConeSumPtEtaBandUEClusterTrigEtaPhi;      //!<! Cluster Sum Pt in the eta band for clusters, per eta-phi bin of trigger,before normalization.
  TH2F *   fhConeSumPtPhiBandUEClusterTrigEtaPhi;      //!<! Cluster Sum Pt in the phi band for clusters, per eta-phi bin of trigger, before normalization.
  TH2F *   fhConeSumPtEtaBandUETrackTrigEtaPhi;        //!<! Track Sum Pt in the eta band for tracks, per eta-phi bin of trigger, before normalization.
  TH2F *   fhConeSumPtPhiBandUETrackTrigEtaPhi;        //!<! Track Sum Pt in the phi badn for tracks, per eta-phi bin of trigger, before normalization.
  
  TH2F *   fhConeSumPtVSUETracksEtaBand;               //!<! Tracks, eta band: sum pT in cone vs bkg to subtract.
  TH2F *   fhConeSumPtVSUETracksPhiBand;               //!<! Tracks, phi band:  sum pT in cone vs bkg to subtract.
  TH2F *   fhConeSumPtVSUEClusterEtaBand;              //!<! Clusters, eta band: sum pT in cone vs bkg to subtract.
  TH2F *   fhConeSumPtVSUEClusterPhiBand;              //!<! Clusters, phi band:  sum pT in cone vs bkg to subtract.
  
  TH2F *   fhConeSumPtEtaUESub;                        //!<! Cluster and tracks Sum Pt in the cone after bkg subtraction, vs pT trigger.
  TH2F *   fhConeSumPtPhiUESub;                        //!<! Cluster and tracks Sum Pt in the cone after bkg subtraction, vs pT trigger.
  TH2F *   fhConeSumPtEtaUESubTrigEtaPhi;              //!<! Cluster and tracks Sum Pt in the cone after bkg subtraction, vs eta-phi trigger.
  TH2F *   fhConeSumPtPhiUESubTrigEtaPhi;              //!<! Cluster and tracks Sum Pt in the cone after bkg subtraction, vs eta-phi trigger.
  
  TH2F *   fhConeSumPtEtaUENormCluster;                //!<! Cluster Sum Pt in the normalized eta UE cone vs pT trigger.
  TH2F *   fhConeSumPtPhiUENormCluster;                //!<! Cluster Sum Pt in the normalized phi UE cone vs pT trigger.
  TH2F *   fhConeSumPtEtaUESubCluster;                 //!<! Cluster Sum Pt in the cone after bkg subtraction, vs pT trigger.
  TH2F *   fhConeSumPtPhiUESubCluster;                 //!<! Cluster Sum Pt in the cone after bkg subtraction, vs pT trigger.
  TH2F *   fhConeSumPtEtaUESubClusterTrigEtaPhi;       //!<! Cluster Sum Pt in the cone after bkg subtraction, vs eta-phi trigger.
  TH2F *   fhConeSumPtPhiUESubClusterTrigEtaPhi;       //!<! Cluster Sum Pt in the cone after bkg subtraction, vs eta-phi trigger.

  TH2F *   fhConeSumPtEtaUENormTrack;                  //!<! Track Sum Pt in the normalized eta UE cone vs pT trigger.
  TH2F *   fhConeSumPtPhiUENormTrack;                  //!<! Track Sum Pt in the normalized phi UE cone vs pT trigger.
  TH2F *   fhConeSumPtEtaUESubTrack;                   //!<! Track Sum Pt in the cone after bkg subtraction, vs pT trigger.
  TH2F *   fhConeSumPtPhiUESubTrack;                   //!<! Track Sum Pt in the cone after bkg subtraction, vs pT trigger.
  TH2F *   fhConeSumPtEtaUESubTrackTrigEtaPhi;         //!<! Track Sum Pt in the cone after bkg subtraction, vs eta-phi trigger.
  TH2F *   fhConeSumPtPhiUESubTrackTrigEtaPhi;         //!<! Track Sum Pt in the cone after bkg subtraction, vs eta-phi trigger.
  
  TH2F *   fhFractionTrackOutConeEta;                  //!<! Fraction of cone out of tracks acceptance in eta.
  TH2F *   fhFractionTrackOutConeEtaTrigEtaPhi;        //!<! Fraction of cone out of tracks acceptance in eta, vs trigger eta-phi.
  TH2F *   fhFractionClusterOutConeEta;                //!<! Fraction of cone out of clusters acceptance in eta.
  TH2F *   fhFractionClusterOutConeEtaTrigEtaPhi;      //!<! Fraction of cone out of clusters acceptance in eta, vs trigger eta-phi.
  TH2F *   fhFractionClusterOutConePhi;                //!<! Fraction of cone out of clusters acceptance in phi.
  TH2F *   fhFractionClusterOutConePhiTrigEtaPhi;      //!<! Fraction of cone out of clusters acceptance in phi, vs trigger eta-phi.
  
  TH2F *   fhConeSumPtEtaUESubClustervsTrack ;         //!<! Cluster vs tracks Sum Pt Sum Pt in the cone, after subtraction in eta band.
  TH2F *   fhConeSumPtPhiUESubClustervsTrack ;         //!<! Cluster vs tracks Sum Pt Sum Pt in the cone, after subtraction in phi band.
  
  TH2F *   fhEtaBandClustervsTrack ;                   //!<! Accumulated pT in Eta band to estimate UE in cone, clusters vs tracks.
  TH2F *   fhPhiBandClustervsTrack ;                   //!<! Accumulated pT in Phi band to estimate UE in cone, clusters vs tracks.
  TH2F *   fhEtaBandNormClustervsTrack ;               //!<! Accumulated pT in Eta band to estimate UE in cone, normalized to cone size, clusters vs tracks.
  TH2F *   fhPhiBandNormClustervsTrack ;               //!<! Accumulated pT in Phi band to estimate UE in cone, normalized to cone size, clusters vs tracks.
  
  TH2F *   fhConeSumPtSubvsConeSumPtTotPhiTrack;       //!<! Tracks, phi band: sum pT in cone after bkg sub vs sum pT in cone before bkg sub.
  TH2F *   fhConeSumPtSubNormvsConeSumPtTotPhiTrack;   //!<! Tracks, phi band: sum pT in cone after bkg sub normalized by sum pT in cone before bkg sub vs sum pT in cone before bkg sub.
  TH2F *   fhConeSumPtSubvsConeSumPtTotEtaTrack;       //!<! Tracks, eta band: sum pT in cone after bkg sub vs sum pT in cone before bkg sub
  TH2F *   fhConeSumPtSubNormvsConeSumPtTotEtaTrack;   //!<! Tracks, eta band: sum pT in cone after bkg sub normalized by sum pT in cone before bkg sub vs sum pT in cone before bkg sub.
  TH2F *   fhConeSumPtSubvsConeSumPtTotPhiCluster;     //!<! Clusters, phi band: sum pT in cone after bkg sub vs sum pT in cone before bkg sub.
  TH2F *   fhConeSumPtSubNormvsConeSumPtTotPhiCluster; //!<! Clusters, phi band: sum pT in cone after bkg sub normalized by sum pT in cone before bkg sub vs sum pT in cone before bkg sub.
  TH2F *   fhConeSumPtSubvsConeSumPtTotEtaCluster;     //!<! Clusters, eta band: sum pT in cone after bkg sub vs sum pT in cone before bkg sub
  TH2F *   fhConeSumPtSubNormvsConeSumPtTotEtaCluster; //!<! Clusters, eta band: sum pT in cone after bkg sub normalized by sum pT in cone before bkg sub vs sum pT in cone before bkg sub.
  
  /// Copy constructor not implemented.
  AliIsolationCut(              const AliIsolationCut & g) ;

  /// Assignment operator not implemented.
  AliIsolationCut & operator = (const AliIsolationCut & g) ; 

  /// \cond CLASSIMP
  ClassDef(AliIsolationCut,12) ;
  /// \endcond

} ;

#endif //ALIISOLATIONCUT_H



