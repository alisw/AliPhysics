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
/// An AOD candidate (AliAODPWG4ParticleCorrelation type)
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
#include <TLorentzVector.h>

// --- ANALYSIS system ---
class AliAODPWG4ParticleCorrelation ;
class AliCaloTrackReader ;
class AliCaloPID;

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

  Float_t    GetCellDensity(  AliAODPWG4ParticleCorrelation * pCandidate,
                              AliCaloTrackReader * reader) const ;

  void       MakeIsolationCut(TObjArray * plCTS, TObjArray * plNe,
                              AliCaloTrackReader * reader,
                              AliCaloPID * pid,
                              Bool_t bFillAOD,
                              AliAODPWG4ParticleCorrelation  * pCandidate, TString aodObjArrayName,
                              Int_t &n, Int_t & nfrac, Float_t &ptSum, Float_t &ptLead, Bool_t & isolated) ;

  void       Print(const Option_t * opt) const ;

  Float_t    Radius(Float_t etaCandidate, Float_t phiCandidate, Float_t eta, Float_t phi) const ;

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

  void 	     GetCoeffNormBadCell(AliAODPWG4ParticleCorrelation * pCandidate,
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
    
 private:

  Float_t    fConeSize ;         ///< Size of the isolation cone

  Float_t    fPtThreshold ;      ///< Minimum pt of the particles in the cone or sum in cone (UE pt mean in the forward region cone)

  Float_t    fPtThresholdMax ;   ///< Maximum pt of the particles outside the cone (needed to fit shower distribution isolated/non-isolated particles)

  Float_t    fSumPtThreshold ;   ///< Minimum of sum pt of the particles in the cone (UE sum in the forward region cone)

  Float_t    fSumPtThresholdMax ;///< Maximum of sum pt of the particles in the cone (UE sum in the forward region cone)

  Float_t    fPtFraction ;       ///< Fraction of the momentum of particles in cone or sum in cone.

  Int_t      fICMethod ;         ///< Isolation cut method to be used: kPtIC, kSumPtIC, kPtFracIC, kSumPtFracIC.

  Int_t      fPartInCone;        ///< Type of particles inside cone: kNeutralAndCharged, kOnlyNeutral, kOnlyCharged.

  Int_t      fDebug;             ///< Debug level.

  Bool_t     fFracIsThresh;      ///< Use threshold instead of fraction when pt leading is small.

  Bool_t     fIsTMClusterInConeRejected; ///< Enable to remove the Track matching removal of clusters in cone sum pt calculation in case of kNeutralAndCharged analysis
  
  Float_t    fDistMinToTrigger;  ///<  Minimal distance between isolation candidate particle and particles in cone to count them for this isolation.
  
  TLorentzVector fMomentum;      //!<! Momentum of cluster, temporal object.

  TVector3   fTrackVector;       //!<! Track moment, temporal object.

  /// Copy constructor not implemented.
  AliIsolationCut(              const AliIsolationCut & g) ;

  /// Assignment operator not implemented.
  AliIsolationCut & operator = (const AliIsolationCut & g) ; 

  /// \cond CLASSIMP
  ClassDef(AliIsolationCut,11) ;
  /// \endcond

} ;

#endif //ALIISOLATIONCUT_H



