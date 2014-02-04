#ifndef ALIISOLATIONCUT_H
#define ALIISOLATIONCUT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
// Class containing methods for the isolation cut. 
// An AOD candidate (AliAODPWG4ParticleCorrelation type)
// is passed. Look in a cone around the candidate and study
// the hadronic activity inside to decide if the candidate is isolated
//
// -- Author: Gustavo Conesa (INFN-LNF)
//
// -- Yaxian Mao (add the possibility for different IC method with different pt range, 01/10/2010)

// --- ROOT system --- 
#include <TObject.h>
class TObjArray ;

// --- ANALYSIS system ---
class AliAODPWG4ParticleCorrelation ;
class AliCaloTrackReader ;
class AliCaloPID; 

class AliIsolationCut : public TObject {
  
 public: 
  
  AliIsolationCut() ;            // default ctor
  virtual ~AliIsolationCut() {;} // virtual dtor
 
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
                              Int_t &n, Int_t & nfrac, Float_t &ptsum, Bool_t & isolated) const ;  
  
  void       Print(const Option_t * opt) const ;
  
  Float_t    Radius(Float_t etaCandidate, Float_t phiCandidate, Float_t eta, Float_t phi) const ;
  
  // Cone background studies medthods
  
  Float_t  CalculateExcessAreaFraction(Float_t excess) const ;

  void     CalculateUEBandClusterNormalization(AliCaloTrackReader * reader,     Float_t   etaC, Float_t phiC,
                                               Float_t   phiUEptsumCluster,     Float_t   etaUEptsumCluster,
                                               Float_t & phiUEptsumClusterNorm, Float_t & etaUEptsumClusterNorm,
                                               Float_t & excessFracEta,         Float_t & excessFracPhi              ) const ;
  
  void     CalculateUEBandTrackNormalization  (AliCaloTrackReader * reader,     Float_t   etaC,  Float_t phiC,
                                               Float_t   phiUEptsumTrack,       Float_t   etaUEptsumTrack,
                                               Float_t & phiUEptsumTrackNorm,   Float_t & etaUEptsumTrackNorm,
                                               Float_t & excessFracEta,         Float_t & excessFracPhi              )   const ;

  void 	   GetCoeffNormBadCell(AliAODPWG4ParticleCorrelation * pCandidate,
                               AliCaloTrackReader * reader, 
                               Float_t & coneBadCellsCoeff,
                               Float_t & etaBandBadCellsCoeff  , Float_t & phiBandBadCellsCoeff) ;

  
  // Parameter setters and getters
  
  Float_t    GetConeSize()            const { return fConeSize       ; }
  Float_t    GetPtThreshold()         const { return fPtThreshold    ; }
  Float_t    GetPtThresholdMax()      const { return fPtThresholdMax    ; }
  Float_t    GetSumPtThreshold()      const { return fSumPtThreshold ; }
  Float_t    GetPtFraction()          const { return fPtFraction     ; }
  Int_t      GetICMethod()            const { return fICMethod       ; }
  Int_t      GetParticleTypeInCone()  const { return fPartInCone     ; }
  Bool_t     GetFracIsThresh()        const { return fFracIsThresh   ; }
	
  void       SetConeSize(Float_t r)         { fConeSize       = r    ; }
  void       SetPtThreshold(Float_t pt)     { fPtThreshold    = pt   ; }
  void       SetPtThresholdMax(Float_t pt)  { fPtThresholdMax    = pt   ; }
  void       SetSumPtThreshold(Float_t s)   { fSumPtThreshold = s    ; }
  void       SetPtFraction(Float_t pt)      { fPtFraction     = pt   ; }
  void       SetICMethod(Int_t i )          { fICMethod       = i    ; }
  void       SetParticleTypeInCone(Int_t i) { fPartInCone     = i    ; }
  void       SetDebug(Int_t d)              { fDebug          = d    ; }
  void       SetFracIsThresh(Bool_t f )     { fFracIsThresh   = f    ; }
 private:
  
  Float_t    fConeSize ;       // Size of the isolation cone 
  Float_t    fPtThreshold ;    // Mimium pt of the particles in the cone or sum in cone (UE pt mean in the forward region cone)
  Float_t    fPtThresholdMax ; // Maximum pt of the particles outside the cone (needed to fit shower distribution isolated/non-isolated particles)
  Float_t    fSumPtThreshold ; // Minium of sum pt of the particles in the cone (UE sum in the forward region cone)
  Float_t    fPtFraction ;     // Fraction of the momentum of particles in cone or sum in cone
  Int_t      fICMethod ;       // Isolation cut method to be used
                               // kPtIC: Pt threshold method
                               // kSumPtIC: Cone pt sum method
                               // kPtFracIC:   Pt threshold, fraction of candidate pt, method
                               // kSumPtFracIC:   Cone pt sum , fraction of cone sum, method
  Int_t      fPartInCone;      // Type of particles inside cone:
                               // kNeutralAndCharged, kOnlyNeutral, kOnlyCharged

  Int_t      fDebug;           // Debug level
  Bool_t     fFracIsThresh;    // Use threshold instead of fraction when pt leading is small
  
  AliIsolationCut(              const AliIsolationCut & g) ; // cpy ctor
  AliIsolationCut & operator = (const AliIsolationCut & g) ; // cpy assignment
  
  ClassDef(AliIsolationCut,6)
} ;


#endif //ALIISOLATIONCUT_H



