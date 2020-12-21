// Class for pair information
// Author: Ionut-Cristian Arsene (iarsene@cern.ch)
// 

#ifndef ALIREDUCEDPAIRINFO_H
#define ALIREDUCEDPAIRINFO_H

#include <TMath.h>
#include "AliReducedBaseTrack.h"


//_____________________________________________________________________
class AliReducedPairInfo : public AliReducedBaseTrack {

  friend class AliAnalysisTaskReducedTreeMaker;  // friend analysis task which fills the object

  public:
  enum CandidateType {
    kGammaConv=0,            // gamma        -> e+ e-  (conversion)
    kK0sToPiPi,              // K0s          -> pi+ pi-
    kLambda0ToPPi,           // Lambda0      -> p pi-
    kALambda0ToPPi,          // anti-Lambda0 -> anti-p pi+
    kPhiToKK,                // phi          -> K+ K-
    kJpsiToEE,               // J/psi        -> e+ e-
    kPsi2SToEE,              // psi(2S)      -> e+ e- 
    kUpsilon,                // Upsilon      -> e+ e-
    kDplusToK0sPiplus,       // D+           -> K0s pi+
    kDplusToK0sKplus,        // D+           -> K0s K+
    kDplusToPhiPiplus,       // D+           -> phi pi+
    kDminusToK0sPiminus,     // D-           -> K0s pi-
    kDminusToK0sKminus,      // D-           -> K0s K-
    kDminusToPhiPiminus,     // D-           -> phi pi-
    kDzeroToKminusPiplus,    // D0           -> K- pi+
    kADzeroToKplusPiminus,   // anti-D0      -> K+ pi-
    kDsplusToK0sKplus,       // Ds+          -> K0s K+
    kDsminusToK0sKminus,     // Ds-          -> K0s K-
    kNMaxCandidateTypes
  };
  static const Char_t* fgkDecayChannelNames[kNMaxCandidateTypes][4];
  
  AliReducedPairInfo();
  AliReducedPairInfo(const AliReducedPairInfo &c);
  virtual ~AliReducedPairInfo();

  // setters
  void CandidateId(Char_t id) {fCandidateId=id;}
  void PairType(Char_t type) {fPairType=type;}
  void PairTypeSPD(Char_t type) {fPairTypeSPD=type;}
  void SetMass(Float_t m) {fMass[0]=m;}
  void SetLegIds(UShort_t leg1, UShort_t leg2) {fLegIds[0]=leg1;fLegIds[1]=leg2;}
  void SetLxy(Float_t lxy) {fLxy = lxy;}
  void SetPseudoProper(Float_t lpsproper) {fPsProper = lpsproper;}
  void SetPointingAngle(Float_t pa) {fPointingAngle = pa;}
  void SetChisquare(Float_t chi2) {fChisquare = chi2;}
  
  // getters
  Char_t   CandidateId()         const {return fCandidateId;}
  Char_t   PairType()            const {return fPairType;}
  Char_t   PairTypeSPD()            const {return fPairTypeSPD;}
  Int_t    LegId(Int_t leg)      const {return (leg==0 || leg==1 ? fLegIds[leg] : -1);}
  Float_t  Mass(Int_t idx=0)     const {return (idx>=0 && idx<4 ? fMass[idx] : -999.);}
  Float_t  Energy()              const;
  Float_t  Rapidity()            const;
  Float_t  Lxy()                 const {return fLxy;}
  Float_t  PsProper()                 const {return fPsProper;}
  Float_t  DecayRadius()         const {return fLxy;}
  Float_t  PointingAngle()       const {return fPointingAngle;}
  Float_t  Chi2()                const {return fChisquare;}
  Bool_t   IsOnTheFly()          const {return fPairType;}
  Bool_t   IsPureV0K0s()         const {return (fQualityFlags&(UInt_t(1)<<1));}
  Bool_t   IsPureV0Lambda()      const {return (fQualityFlags&(UInt_t(1)<<2));}
  Bool_t   IsPureV0ALambda()     const {return (fQualityFlags&(UInt_t(1)<<3));}
  Bool_t   IsPureV0Gamma()       const {return (fQualityFlags&(UInt_t(1)<<4));}
  UInt_t   QualityFlags()        const {return fQualityFlags;}
  
 protected:
  Char_t  fCandidateId;         // candidate type (K0s, Lambda, J/psi, phi, etc)
  Char_t  fPairType;            // 0- offline, 1- on the fly for V0 candidates; 0 ++; 1 +-; 2 -- for other pairs; 
  Char_t  fPairTypeSPD;         // 2 both / 1 one / 0 none of the two legs has an hit in the first layers of the SPD; 
  UShort_t fLegIds[2];          // leg ids 
  Float_t fMass[4];             // invariant mass for pairs (3 extra mass values for other V0 pid assumptions)
                                // idx=0 -> K0s assumption; idx=1 -> Lambda; idx=2 -> anti-Lambda; idx=3 -> gamma conversion
  Float_t fLxy;                 // pseudo-proper decay length (pair candidates) or radius of the secondary vertex for V0s 
  Float_t fPsProper;                 // pseudo-proper decay length (pair candidates) or radius of the secondary vertex for V0s 
  Float_t fPointingAngle;       // angle between the pair momentum vector and the secondary vertex position vector
  Float_t fChisquare;           // chi2 for the legs matching
  
  AliReducedPairInfo& operator= (const AliReducedPairInfo &c);

  ClassDef(AliReducedPairInfo, 2);
};

//_______________________________________________________________________________
inline Float_t AliReducedPairInfo::Energy() const 
{
  //
  // Return the energy
  //
  Float_t mass=fMass[0];
  switch (fCandidateId) {
    case kK0sToPiPi:
      mass = fMass[0];
      break;
    case kLambda0ToPPi:
      mass = fMass[1];
      break;
    case kALambda0ToPPi:
      mass = fMass[2];
      break;
    case kGammaConv:
      mass = fMass[3];
      break;
    default:
      mass = fMass[0];
      break;    
  }
  Float_t p = P();
  return TMath::Sqrt(mass*mass+p*p);
}


//_______________________________________________________________________________
inline Float_t AliReducedPairInfo::Rapidity() const
{
  //
  // return rapidity
  //
  Float_t e = Energy();
  Float_t pz = Pz();
  if(e-TMath::Abs(pz)>1.0e-10)
    return 0.5*TMath::Log((e+pz)/(e-pz));
  else 
    return -999.;
}

#endif
