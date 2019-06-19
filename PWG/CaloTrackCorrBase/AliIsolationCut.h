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
                                        Double_t histoWeight = 1) ;
  
  // Cone background studies medthods

  Float_t    CalculateExcessAreaFraction(Float_t excess) const ;

  void       CalculateUEBandClusterNormalization(AliCaloTrackReader * reader,     Float_t   etaC, Float_t phiC,
                                                 Float_t   phiUEptsumCluster,     Float_t   etaUEptsumCluster,
                                                 Float_t & phiUEptsumClusterNorm, Float_t & etaUEptsumClusterNorm,
                                                 Float_t & excessFracEta,         Float_t & excessFracPhi              ) const ;

  void       CalculateUEBandTrackNormalization  (AliCaloTrackReader * reader,     Float_t   etaC,  Float_t phiC,
                                                 Float_t   phiUEptsumTrack,       Float_t   etaUEptsumTrack,
                                                 Float_t & phiUEptsumTrackNorm,   Float_t & etaUEptsumTrackNorm,
                                                 Float_t & excessFracEta,         Float_t & excessFracPhi              ) const ;

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
  
 private:

  Bool_t     fFillHistograms;    ///< Fill histograms if GetCreateOuputObjects() was called. 

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
  
  Float_t    fDistMinToTrigger;  ///<  Minimal distance between isolation candidate particle and particles in cone to count them for this isolation.
  
  Int_t      fDebug;             ///< Debug level.

  TLorentzVector fMomentum;      //!<! Momentum of cluster, temporal object.

  TVector3   fTrackVector;       //!<! Track moment, temporal object.
  
  // Histograms
  
  AliHistogramRanges * fHistoRanges;                   ///!  Histogram bins and ranges  data-base
  
  TH2F *   fhPtInCone ;                                //!<! Cluster/track Pt in the cone.
  TH2F *   fhPtClusterInCone ;                         //!<! Cluster Pt in the cone.
  TH2F *   fhPtTrackInCone ;                           //!<! Track Pt in the cone.

  TH2F *   fhConeSumPt ;                               //!<! Cluster and tracks Sum Pt in the cone.
  TH2F *   fhConeSumPtCluster ;                        //!<! Clusters Sum Pt in the cone.
  TH2F *   fhConeSumPtTrack ;                          //!<! Tracks Sum Pt in the cone.
  
  TH2F *   fhConePtLead ;                              //!<! Cluster and tracks leading pt in the cone.
  TH2F *   fhConePtLeadCluster ;                       //!<! Clusters leading pt in the cone.
  TH2F *   fhConePtLeadTrack ;                         //!<! Tracks leading pt in the cone.
  
  TH2F *   fhConeSumPtClustervsTrack ;                 //!<! Cluster vs tracks Sum Pt Sum Pt in the cone.
  
  TH2F *   fhConeSumPtClusterTrackFrac ;               //!<! Cluster / tracks Sum Pt Sum Pt in the cone.
  
  TH2F *   fhConePtLeadClustervsTrack;                 //!<! Tracks vs Clusters leading pt.
  
  TH2F *   fhConePtLeadClusterTrackFrac;               //!<! Trigger pt vs cluster/track leading pt.
  
  TH2F *   fhConeSumPtTrigEtaPhi ;                     //!<! Cluster and tracks Sum Pt Sum Pt in the cone, per eta-phi bin of trigger.

  
  
  /// Copy constructor not implemented.
  AliIsolationCut(              const AliIsolationCut & g) ;

  /// Assignment operator not implemented.
  AliIsolationCut & operator = (const AliIsolationCut & g) ; 

  /// \cond CLASSIMP
  ClassDef(AliIsolationCut,12) ;
  /// \endcond

} ;

#endif //ALIISOLATIONCUT_H



