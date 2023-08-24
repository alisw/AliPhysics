// Base class for track information
// Author: Ionut-Cristian Arsene (iarsene@cern.ch)
// 

#ifndef ALIREDUCEDBASETRACK_H
#define ALIREDUCEDBASETRACK_H

#include <TMath.h>
#include <TObject.h>

//_____________________________________________________________________
class AliReducedBaseTrack : public TObject {
  
  friend class AliAnalysisTaskReducedTreeMaker;  // friend analysis task which fills the object
  
  public:
    AliReducedBaseTrack();
    AliReducedBaseTrack(const AliReducedBaseTrack &c);
    virtual ~AliReducedBaseTrack();
  
    // getters
    UShort_t TrackId()                     const {return fTrackId;}
    Float_t Px()  const {return (fIsCartesian ? fP[0] : TMath::Abs(fP[0])*TMath::Cos(fP[1]));}
    Float_t Py()  const {return (fIsCartesian ? fP[1] : TMath::Abs(fP[0])*TMath::Sin(fP[1]));}
    Float_t Pz()  const {return (fIsCartesian ? fP[2] : TMath::Abs(fP[0])*TMath::SinH(fP[2]));}
    Float_t P()   const {return (fIsCartesian ? TMath::Sqrt(fP[0]*fP[0]+fP[1]*fP[1]+fP[2]*fP[2]) : TMath::Abs(fP[0])*TMath::CosH(fP[2]));}
    Float_t Phi() const;
    Float_t Pt()  const {return (fIsCartesian ? TMath::Sqrt(fP[0]*fP[0]+fP[1]*fP[1]) : fP[0]);}
    Float_t Eta() const;
    Float_t Rapidity(Float_t massAssumption) const;
    Float_t Theta() const;
    Float_t Energy(Float_t massAssumption) const {return TMath::Sqrt(massAssumption*massAssumption+P()*P());}
    Bool_t  TestFlag(UShort_t iflag)       const {return ((iflag<(8*sizeof(ULong_t))) ? fFlags&(ULong_t(1)<<iflag) : kFALSE);} 
    ULong_t GetFlags()                     const {return fFlags;}
    Int_t  Charge()                        const {return fCharge;} 
    Bool_t IsCartesian() const {return fIsCartesian;}           
    
    ULong_t GetQualityFlags()             const {return fQualityFlags;}
    UInt_t  GetFirstHalfOfQualityFlags() const;
    UInt_t  GetSecondHalfOfQualityFlags() const;
    Bool_t UsedForQvector()               const {return fQualityFlags&(UShort_t(1)<<0);}
    Bool_t TestQualityFlag(UShort_t iflag)  const {return ((iflag<(8*sizeof(ULong_t))) ? fQualityFlags&(ULong_t(1)<<iflag) : kFALSE);}
    Bool_t IsMCTruth()                        const {return (fIsMCTruth ? kTRUE : kFALSE);}
    //Bool_t HasMCTruthInfo()               const {return (fMCFlags ? kTRUE : kFALSE);}
    Bool_t IsTRDmatch()                     const {return fQualityFlags&(ULong_t(1)<<26);}
    Bool_t IsGammaLeg()                   const {return fQualityFlags&(ULong_t(1)<<1);}
    Bool_t IsPureGammaLeg()            const {return fQualityFlags&(ULong_t(1)<<8);}
    Bool_t IsK0sLeg()                           const {return fQualityFlags&(ULong_t(1)<<2);}
    Bool_t IsPureK0sLeg()                   const {return fQualityFlags&(ULong_t(1)<<9);}
    Bool_t IsLambdaLeg()                   const {return fQualityFlags&(ULong_t(1)<<3);}
    Bool_t IsPureLambdaLeg()            const {return fQualityFlags&(ULong_t(1)<<10);}
    Bool_t IsALambdaLeg()                 const {return fQualityFlags&(ULong_t(1)<<4);}
    Bool_t IsPureALambdaLeg()          const {return fQualityFlags&(ULong_t(1)<<11);}
    Bool_t IsKink(Int_t i=0)               const {return (i>=0 && i<3 ? (fQualityFlags&(ULong_t(1)<<(5+i))) : kFALSE);}
    Bool_t IsKinkNegativeLabel(Int_t i=0)     const {return (i>=0 && i<3 ? (fQualityFlags&(ULong_t(1)<<(12+i))) : kFALSE);}
    Float_t GetBayesProb(Int_t specie)  const { return (fQualityFlags&(ULong_t(1)<<(15+specie)) ? (fQualityFlags&(ULong_t(1)<<21) ? 0.9 : (fQualityFlags&(ULong_t(1)<<20) ? 0.8 : (fQualityFlags&(ULong_t(1)<<19)           ? 0.7 : 0.5)))   : 0.0);}
    
    UInt_t GetMCFlags()                      const {return fMCFlags;}
    Bool_t TestMCFlag(UShort_t iflag) const {return ((iflag<(8*sizeof(UInt_t))) ? fMCFlags & (UInt_t(1)<<iflag) : kFALSE);}
    Bool_t IsMCKineParticle()              const {return (fIsMCTruth ? kTRUE : kFALSE);}
    
    // setters
    void   Px(Float_t px) {fP[0] = px; fIsCartesian=kTRUE;}
    void   Py(Float_t py) {fP[1] = py; fIsCartesian=kTRUE;}
    void   Pz(Float_t pz) {fP[2] = pz; fIsCartesian=kTRUE;}
    void   PxPyPz(Float_t px, Float_t py, Float_t pz) {fP[0]=px;fP[1]=py;fP[2]=pz;fIsCartesian=kTRUE;}
    void   Pt(Float_t pt) {fP[0] = pt; fIsCartesian=kFALSE;}
    void   Phi(Float_t phi) {fP[1] = phi; fIsCartesian=kFALSE;}
    void   Eta(Float_t eta) {fP[2] = eta; fIsCartesian=kFALSE;}
    void   PtPhiEta(Float_t pt, Float_t phi, Float_t eta) {fP[0]=pt;fP[1]=phi;fP[2]=eta;fIsCartesian=kFALSE;}
    void   Charge(Int_t ch) {fCharge=ch;}
    void   ResetFlags() {fFlags=0;}
    void   SetFlags(ULong_t flags) {fFlags=flags;}
    Bool_t SetFlag(UShort_t iflag)  {if(iflag>=8*sizeof(ULong_t)) return kFALSE; fFlags|=(ULong_t(1)<<iflag); return kTRUE;}
    Bool_t UnsetFlag(UShort_t iflag) {if(iflag>=8*sizeof(ULong_t)) return kFALSE; if(TestFlag(iflag)) fFlags^=(ULong_t(1)<<iflag); return kTRUE;}
    void   ResetQualityFlags() {fQualityFlags=0;}
    void   SetQualityFlags(ULong_t flags) {fQualityFlags=flags;}
    Bool_t SetQualityFlag(UShort_t iflag)      {if (iflag>=8*sizeof(ULong_t)) return kFALSE; fQualityFlags|=(ULong_t(1)<<iflag); return kTRUE;}
    Bool_t UnsetQualityFlag(UShort_t iflag)  {if (iflag>=8*sizeof(ULong_t)) return kFALSE; if(TestQualityFlag(iflag)) fQualityFlags^=(ULong_t(1)<<iflag); return kTRUE;}
    void    SetMCFlags(UInt_t flags) {fMCFlags = flags;}
    Bool_t SetMCFlag(UShort_t iflag) {if(iflag>=8*sizeof(UInt_t)) return kFALSE; fMCFlags |= (UInt_t(1)<<iflag); return kTRUE;}
    Bool_t UnsetMCFlag(UShort_t iflag) {if(iflag>=8*sizeof(UInt_t)) return kFALSE; if(TestMCFlag(iflag)) fMCFlags ^= (UInt_t(1)<<iflag); return kTRUE;}
    void SetIsMCTruth(Bool_t flag=kTRUE) {fIsMCTruth = flag;}
   
  protected:
    UShort_t fTrackId;            // track id 
    Float_t fP[3];         // 3-momentum vector
    Bool_t  fIsCartesian;  // if false then the 3-momentum vector is in spherical coordinates (pt,phi,eta)
    Char_t  fCharge;       // electrical charge
    ULong_t fFlags;        // flags reserved for various operations
    ULong_t fQualityFlags;          // BIT0 toggled if track used for TPC event plane
                                                   // BIT1 toggled if track belongs to a gamma conversion
                                                   // BIT2 toggled if track belongs to a K0s
                                                   // BIT3 toggled if track belongs to a Lambda
                                                   // BIT4 toggled if track belongs to an Anti-Lambda
                                                   // BIT5 toggled if the track has kink0 index > 0
                                                   // BIT6 toggled if the track has kink1 index > 0
                                                   // BIT7 toggled if the track has kink2 index > 0
                                                   // BIT8 toggled if track belongs to a pure gamma conversion
                                                   // BIT9 toggled if track belongs to a pure K0s
                                                   // BIT10 toggled if track belongs to a pure Lambda
                                                   // BIT11 toggled if track belongs to a pure ALambda
                                                   // BIT12 toggled if the track has kink0 index < 0
                                                   // BIT13 toggled if the track has kink1 index < 0
                                                   // BIT14 toggled if the track has kink2 index < 0
                                                   // AOD
                                                   // BIT(15+i) toggled if track has filter bit 0+i , 0 <= i <= 10
                                                   // BIT26 toggled if this is a track matched in TRD
                                                   // BIT27 toggled if the track was used in the vertex determination
                                                   
                                                   // For AliReducedPairInfo objects
                                                   // BIT1 toggled for pure V0 K0s candidates
                                                   // BIT2 toggled for pure V0 Lambda candidates
                                                   // BIT3 toggled for pure V0 anti-Lambda candidates
                                                   // BIT4 toggled for pure V0 photon candidates
    UInt_t fMCFlags;                     // flags used to encode Monte-Carlo information
    Bool_t fIsMCTruth;                  // TRUE: if the particle is pure MC track; FALSE: if the particle is a reconstructed track
        
    AliReducedBaseTrack& operator= (const AliReducedBaseTrack &c);
    
    ClassDef(AliReducedBaseTrack, 6)
};

//_______________________________________________________________________________
inline Float_t AliReducedBaseTrack::Phi() const {
  //
  // Return the azimuthal angle of this particle
  //
  if(!fIsCartesian) return fP[1];
  Float_t phi=TMath::ATan2(fP[1],fP[0]); 
  if(phi>=0.0) 
    return phi;
  else 
    return (TMath::TwoPi()+phi);
}

//_______________________________________________________________________________
inline Float_t AliReducedBaseTrack::Theta() const {
  //
  // Return the polar angle for this particle
  //
  Float_t p=P(); 
  if(p>=1.0e-6) 
    return TMath::ACos(Pz()/p);
  else 
    return 0.0;
}

//_______________________________________________________________________________
inline Float_t AliReducedBaseTrack::Eta() const {
  //
  // Return the pseudorapidity of this particle
  //
  if(!fIsCartesian) return fP[2];
  Float_t eta = TMath::Tan(0.5*Theta());
  if(eta>1.0e-6) 
    return -1.0*TMath::Log(eta);
  else 
    return 0.0;
}

//_______________________________________________________________________________
inline Float_t AliReducedBaseTrack::Rapidity(Float_t massAssumption) const {
  //
  // Return the rapidity of this particle using a massAssumption
  //
  Float_t e = Energy(massAssumption);
  Float_t factor = e-Pz();
  if(TMath::Abs(factor)<1.0e-6) return 0.0;
  factor = (e+Pz())/factor;
  if(factor<1.0e-6) return 0.0;
  return 0.5*TMath::Log(factor);
}

#endif
