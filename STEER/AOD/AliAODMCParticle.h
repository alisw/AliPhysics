#ifndef AliAODMCPARTICLE_H
#define AliAODMCPARTICLE_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/// \class AliAODMCParticle
/// \brief AliVParticle realisation for MC Particles in the AOD
///
/// Stripped dow AliMCParticle
///
/// \author Christian Klein Bösing, CERN

#include <Rtypes.h>
#include <TParticlePDG.h>
#include <TExMap.h>
#include <TString.h>


#include "AliTrackReference.h"
#include "AliVParticle.h"
#include "AliMCParticle.h"

class AliAODEvent;
class TParticle;
class TClonesArray;

class AliAODMCParticle: public AliVParticle {
 public:
  AliAODMCParticle();
  AliAODMCParticle(const AliMCParticle* part, Int_t label=-1,Int_t flag = 0);
  virtual ~AliAODMCParticle(){};
  AliAODMCParticle(const AliAODMCParticle& mcPart); 
  AliAODMCParticle& operator=(const AliAODMCParticle& mcPart);
  
  // Kinematics
    virtual Double_t Px()        const;
    virtual Double_t Py()        const;
    virtual Double_t Pz()        const;
    virtual Double_t Pt()        const;
    virtual Double_t P()         const;
    virtual Bool_t   PxPyPz(Double_t p[3]) const;
  
    virtual void     Momentum(TLorentzVector & lv)  { lv.SetPxPyPzE(Px(),Py(),Pz(),E()) ; }
  
    virtual Double_t OneOverPt() const;
    virtual Double_t Phi()       const;
    virtual Double_t Theta()     const;
    
    virtual Double_t Xv() const;
    virtual Double_t Yv() const;
    virtual Double_t Zv() const;
    virtual Bool_t   XvYvZv(Double_t x[3]) const;  
    virtual Double_t T()  const;
    virtual Double_t Tv() const;

    virtual Double_t E()          const;
    virtual Double_t M()          const;
    
    virtual Double_t Eta()        const;
    virtual Double_t Y()          const;
    
    virtual Short_t Charge()      const;

    virtual Int_t   Label()       const;
    virtual Int_t   GetLabel()    const {return Label();} 
	    
    // PID
    virtual const Double_t *PID() const {return 0;} // return PID object (to be defined, still)

    // 
    virtual Double_t GetCalcMass() const;
    virtual void  SetDaughter(Int_t i,Int_t id){if(i<2)fDaughter[i] = id;}
    virtual Int_t GetDaughterLabel(Int_t i) const { return i ? fDaughter[1] : fDaughter[0];}
    virtual Int_t GetNDaughters  () const { return fDaughter[1]>0 ? fDaughter[1]-fDaughter[0]+1 : (fDaughter[0]>0 ? 1:0 ) ;}
    virtual void  SetMother(Int_t im){fMother = im;}
    virtual Int_t GetMother() const {return fMother;}

    virtual Int_t   GetDaughterFirst()   const {return fDaughter[0];}
    virtual Int_t   GetDaughterLast()    const {return fDaughter[1];}

    virtual void Print(const Option_t *opt = "") const;
    virtual Int_t GetPdgCode() const { return fPdgCode;}
    virtual Int_t PdgCode()    const { return GetPdgCode();}
    virtual void    SetGeneratorIndex(Short_t i) {fGeneratorIndex = i;}
    virtual Short_t GetGeneratorIndex() const {return fGeneratorIndex;}
    enum { kPrimary = 1<<0, kPhysicalPrim = 1<<1, kSecondaryFromWeakDecay = 1<<2, kSecondaryFromMaterial = 1 <<3,
	   kFromSubsidiary = 1<<4 }; // use only the first 8bits!
    virtual void SetFlag(UInt_t flag){fFlag = flag;} // carefull flag encodes three different types of information 
    virtual UInt_t GetFlag() const {return fFlag;}


    // for the status we use the upper 16 bits/2 bytes of the flag word
    void SetStatus(Int_t status){
      if(status<0)return; // a TParticle can have a negative status, catch this here and do nothing
      fFlag &= 0xffff;    // reset the upper bins keep the lower bins
      fFlag |= (((UInt_t)status)<<16); // bit shift by 16
    }

    virtual ULong64_t GetStatus() const {
      /// bit shift by 16

      return ULong64_t(fFlag>>16);
    }

    void        SetMCStatusCode(Int_t status) { SetStatus(status)  ; }
    virtual UInt_t MCStatusCode()       const { return GetStatus() ; }

    // Bitwise operations
    void SetPrimary(Bool_t b = kTRUE){
      if(b)fFlag |= kPrimary;
      else fFlag &= ~kPrimary;
    }
    virtual Bool_t IsPrimary() const {return ((fFlag&kPrimary)==kPrimary);} 

    void SetPhysicalPrimary(Bool_t b = kTRUE){
     if(b)fFlag |= kPhysicalPrim;
     else fFlag &= ~kPhysicalPrim; 
    } 
    Bool_t IsPhysicalPrimary() const {return ((fFlag&kPhysicalPrim)==kPhysicalPrim);} 

    void SetSecondaryFromWeakDecay(Bool_t b = kTRUE){
     if(b)fFlag |= kSecondaryFromWeakDecay;
     else fFlag &= ~kSecondaryFromWeakDecay; 
    } 
    Bool_t IsSecondaryFromWeakDecay() const {return ((fFlag&kSecondaryFromWeakDecay)==kSecondaryFromWeakDecay);} 

    void SetSecondaryFromMaterial(Bool_t b = kTRUE){
     if(b)fFlag |= kSecondaryFromMaterial;
     else fFlag &= ~kSecondaryFromMaterial; 
    } 
    Bool_t IsSecondaryFromMaterial() const {return ((fFlag&kSecondaryFromMaterial)==kSecondaryFromMaterial);} 

    virtual void   SetFromSubsidiaryEvent(Bool_t b = kTRUE){
      if(b)fFlag |= kFromSubsidiary;
      else fFlag &= ~kFromSubsidiary;
    }
    virtual Bool_t IsFromSubsidiaryEvent() const {return ((fFlag&kFromSubsidiary)==kFromSubsidiary);} 

    void SetMCProcessCode(UInt_t mcProcess){
      if(mcProcess>1<<7)return; // should not be larger than 46 (see TMCProcess) allow up to 128
      fFlag &= 0xffff00ff; // keep the upper bins and the lower bins just reset 9-16
      fFlag |= (mcProcess<<8); // bit shift by 8
    }

    UInt_t GetMCProcessCode(){
      return ((fFlag&0xff00)>>8); // just return bit shifted bits 9-16 
    }
    

    
    

    static const char* StdBranchName(){return fgkStdBranchName.Data();}

 private:

    static TString fgkStdBranchName;      ///< Standard branch name


  Int_t            fPdgCode;              ///< PDG code of the particle
  UInt_t           fFlag;                 ///< Flag for indication of primary etc, Status code in the upper 16 bits 17-32, MC process id (AKA UniqueID) in bins 16-9)
  Int_t            fLabel;                ///< Label of the original MCParticle
  Int_t            fMother;               ///< Index of the mother particles
  Int_t            fDaughter[2];          ///< Indices of the daughter particles
  Double32_t       fPx;                   ///< x component of momentum
  Double32_t       fPy;                   ///< y component of momentum
  Double32_t       fPz;                   ///< z component of momentum
  Double32_t       fE;                    ///< Energy

  /// x of production vertex
  Double32_t       fVx;                   //[0.,0.,12]
  /// y of production vertex
  Double32_t       fVy;                   //[0.,0.,12]
  /// z of production vertex
  Double32_t       fVz;                   //[0.,0.,12]
  /// t of production vertex
  Double32_t       fVt;                   //[0.,0.,12]
  Short_t          fGeneratorIndex;       //!<! Index of generator in cocktail
  // Copy the uniquID to another data member? unique ID is correctly handled 
  // via TOBject Copy construct but not by AliVParticle ctor (no passing of 
  // TParticles
  // Need a flag for primaries?

  /*
    const TMCProcess kMCprocesses[kMaxMCProcess] =
    {
     kPNoProcess, kPMultipleScattering, kPEnergyLoss, kPMagneticFieldL,
     kPDecay, kPPair, kPCompton, kPPhotoelectric, kPBrem, kPDeltaRay,
     kPAnnihilation, kPHadronic, kPNoProcess, kPEvaporation, kPNuclearFission,
     kPNuclearAbsorption, kPPbarAnnihilation, kPNCapture, kPHElastic,
     kPHInhelastic, kPMuonNuclear, kPTOFlimit,kPPhotoFission, kPNoProcess,
     kPRayleigh, kPNoProcess, kPNoProcess, kPNoProcess, kPNull, kPStop
    };
  */

  ClassDef(AliAODMCParticle,9)  // AliVParticle realisation for AODMCParticles

};

inline Double_t AliAODMCParticle::Px()        const {return fPx;}
inline Double_t AliAODMCParticle::Py()        const {return fPy;}
inline Double_t AliAODMCParticle::Pz()        const {return fPz;}
inline Double_t AliAODMCParticle::Pt()        const {return TMath::Sqrt(fPx*fPx+fPy*fPy);}
inline Double_t AliAODMCParticle::P()         const {return TMath::Sqrt(fPx*fPx+fPy*fPy+fPz*fPz); }
inline Double_t AliAODMCParticle::OneOverPt() const {return 1. / Pt();}
inline Bool_t   AliAODMCParticle::PxPyPz(Double_t p[3]) const { p[0] = fPx; p[1] = fPy; p[2] = fPz; return kTRUE; }
inline Double_t AliAODMCParticle::Phi()       const {return TMath::Pi()+TMath::ATan2(-fPy,-fPx); }  // note that Phi() returns an angle between 0 and 2pi
inline Double_t AliAODMCParticle::Theta()     const {return (fPz==0)?TMath::PiOver2():TMath::ACos(fPz/P()); }
inline Double_t AliAODMCParticle::Xv()        const {return fVx;}
inline Double_t AliAODMCParticle::Yv()        const {return fVy;}
inline Double_t AliAODMCParticle::Zv()        const {return fVz;}
inline Bool_t   AliAODMCParticle::XvYvZv(Double_t x[3]) const { x[0] = fVx; x[1] = fVy; x[2] = fVz; return kTRUE; }
inline Double_t AliAODMCParticle::T()         const {return fVt;}
inline Double_t AliAODMCParticle::Tv()        const {return fVt;}
inline Double_t AliAODMCParticle::E()         const {return fE;}
inline Double_t AliAODMCParticle::Eta()       const {  
  Double_t pmom = P();
  if (pmom != TMath::Abs(fPz)) return 0.5*TMath::Log((pmom+fPz)/(pmom-fPz));
  else                         return 1.e30;
}


inline Double_t AliAODMCParticle::Y()         const 
{
    Double_t e  = E();
    Double_t pz = Pz();
    
    if (e > TMath::Abs(pz)) { 
	return 0.5*TMath::Log((e+pz)/(e-pz));
    } else { 
	return -999.;
    }
}

inline Int_t AliAODMCParticle::Label()       const {return fLabel;}

inline Double_t AliAODMCParticle::GetCalcMass() const {

  Double_t m2 = E()*E()-Px()*Px()-Py()*Py()-Pz()*Pz();
  if(m2<0)return 0;
  return TMath::Sqrt(m2);
}


#endif
