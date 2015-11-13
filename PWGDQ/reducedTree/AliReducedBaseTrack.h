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
    virtual ~AliReducedBaseTrack();
  
    // getters
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
    Bool_t  TestFlag(UShort_t iflag)       const {return (iflag<8*sizeof(ULong_t) ? fFlags&(ULong_t(1)<<iflag) : kFALSE);} 
    ULong_t GetFlags()                     const {return fFlags;}
    Int_t  Charge()                        const {return fCharge;} 
    Bool_t IsCartesian() const {return fIsCartesian;}           
    
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
        
  protected:
    Float_t fP[3];         // 3-momentum vector
    Bool_t  fIsCartesian;  // if false then the 3-momentum vector is in spherical coordinates (pt,phi,eta)
    Char_t  fCharge;       // electrical charge
    ULong_t fFlags;        // flags reserved for various operations
        
    AliReducedBaseTrack(const AliReducedBaseTrack &c);      
    AliReducedBaseTrack& operator= (const AliReducedBaseTrack &c);
    
    ClassDef(AliReducedBaseTrack, 1)
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
