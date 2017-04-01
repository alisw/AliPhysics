#ifndef AliPicoTrack_H
#define AliPicoTrack_H

#include <TMath.h>

#include "AliVTrack.h"
#include "AliAODMCParticle.h"

class AliVCluster;
class AliAODTrack;

class AliPicoTrack: public AliVTrack {
 public:
  AliPicoTrack();
  AliPicoTrack(Double_t pt, Double_t eta, Double_t phi, Byte_t q, Int_t label, Byte_t type,
               Double_t etaemc=0, Double_t phiemc=0, Double_t ptemc=0, Bool_t ise=0, Double_t mass=0.13957);
  AliPicoTrack(const AliPicoTrack &pc); 
  AliPicoTrack &operator=(const AliPicoTrack &pc);

  Double_t        Px()                        const { return fPt*TMath::Cos(fPhi);  }
  Double_t        Py()                        const { return fPt*TMath::Sin(fPhi);  }
  Double_t        Pz()                        const { return fPt*TMath::SinH(fEta); }
  Double_t        Pt()                        const { return fPt;                   }
  Double_t        P()                         const { return fPt*TMath::CosH(fEta); }
  
  Bool_t          PxPyPz(Double_t p[3])       const { p[0]=Px();p[1]=Py();p[2]=Pz(); return 1;         }
  Bool_t          GetPxPyPz(Double_t p[3])    const { p[0]=Px();p[1]=Py();p[2]=Pz(); return 1;         }
  Double_t        Xv()                        const { return 0.;      }
  Double_t        Yv()                        const { return 0.;      }
  Double_t        Zv()                        const { return 0.;      }
  Bool_t          XvYvZv(Double_t x[3])       const { x[0]=0;x[1]=0;x[2]=0; return 1;                  }
  Double_t        OneOverPt()                 const { return 1./fPt;  }
  Double_t        Phi()                       const { return fPhi;    }
  Double_t        Theta()                     const { return 2*TMath::ATan(TMath::Exp(-fEta));         }
  Double_t        E()                         const { Double_t p=P(); return TMath::Sqrt(M()*M()+p*p); }
  Double_t        M()                         const { return fM; }
  Double_t        Eta()                       const { return fEta;    }
  Double_t        Y()                         const { return 0.5*TMath::Log((E()+Pz())/(E()-Pz()));  }
  Short_t         Charge()                    const { return (char)fQ;}
  Double_t        GetSign()                   const { return fQ;}
  Int_t           GetLabel()                  const { return fLabel;  }
  void            SetLabel(Int_t label)             { fLabel = label; }
  Byte_t          GetTrackType()              const { return fTrackType;}
  void            SetTrackType(Byte_t type)         { fTrackType = type;}
  Int_t           PdgCode()                   const { return 0;       }
  const Double_t *PID()                const { return 0;       }
  Int_t           GetID()                     const { return 0;       }
  UChar_t         GetITSClusterMap()          const { return 0;       }
  Int_t           GetEMCALcluster()           const { return fClusId; }
  void            SetEMCALcluster(Int_t id)         { fClusId = id;   }
  Double_t        GetEtaEmc()                 const { return GetTrackEtaOnEMCal(); }
  Double_t        GetPhiEmc()                 const { return GetTrackPhiOnEMCal(); }
  Bool_t          IsEMCAL()                   const { return fEmcal;  }
  ULong_t         GetStatus()                 const { return 0;       }
  Bool_t          GetXYZ(Double_t *v)         const { v[0]=0; v[1]=0; v[2]=0; return 0; }
  using AliVTrack::GetXYZ;
  Double_t        GetBz()                     const { return 0;       }
  void            GetBxByBz(Double_t b[3])    const { b[0]=0;b[1]=0;b[2]=0; }
  Bool_t          GetCovarianceXYZPxPyPz(Double_t /*cv*/[21]) const { return 0; }
  AliVTrack      *GetTrack()                  const { return fOrig;   }
  void            SetTrack(AliVTrack *t)            { fOrig = t;      }
  Double_t        GetTrackPhiOnEMCal()        const { return fPhiEmc ; }
  Double_t        GetTrackEtaOnEMCal()        const { return fEtaEmc ; }
  Double_t        GetTrackPtOnEMCal()         const { return fPtEmc  ; }
  Double_t        GetTrackPOnEMCal()          const { return fPtEmc*TMath::CosH(fEtaEmc);}
  void            SetTrackPhiEtaPtOnEMCal(Double_t eta, Double_t phi, Double_t pt) { fEtaEmc = eta; fPhiEmc = phi; fPtEmc = pt;}
  void            Clear(Option_t * /*option*/ ="")  { fClusId = -1; fOrig = 0; } 
  Int_t           Compare(const TObject* obj) const;
  Bool_t          PropagateToDCA(const AliVVertex *, Double_t, Double_t, Double_t *, Double_t *) { return 0; }
  void            SetFlag(UInt_t flag) {fFlag = flag;}
  UInt_t          GetFlag() const {return fFlag;}
  void SetPrimary(Bool_t b = kTRUE){
    if(b)fFlag |= AliAODMCParticle::kPrimary;
    else fFlag &= ~AliAODMCParticle::kPrimary;
  }
  Bool_t IsPrimary() const {return ((fFlag&AliAODMCParticle::kPrimary)==AliAODMCParticle::kPrimary);} 
  
  void SetPhysicalPrimary(Bool_t b = kTRUE){
    if(b)fFlag |= AliAODMCParticle::kPhysicalPrim;
    else fFlag &= ~AliAODMCParticle::kPhysicalPrim; 
  } 
  Bool_t IsPhysicalPrimary() const {return ((fFlag&AliAODMCParticle::kPhysicalPrim)==AliAODMCParticle::kPhysicalPrim);} 
  
  void SetSecondaryFromWeakDecay(Bool_t b = kTRUE){
    if(b)fFlag |= AliAODMCParticle::kSecondaryFromWeakDecay;
    else fFlag &= ~AliAODMCParticle::kSecondaryFromWeakDecay; 
  } 
  Bool_t IsSecondaryFromWeakDecay() const {return ((fFlag&AliAODMCParticle::kSecondaryFromWeakDecay)==AliAODMCParticle::kSecondaryFromWeakDecay);} 
  
  void SetSecondaryFromMaterial(Bool_t b = kTRUE){
    if(b)fFlag |= AliAODMCParticle::kSecondaryFromMaterial;
    else fFlag &= ~AliAODMCParticle::kSecondaryFromMaterial; 
  } 
  Bool_t IsSecondaryFromMaterial() const {return ((fFlag&AliAODMCParticle::kSecondaryFromMaterial)==AliAODMCParticle::kSecondaryFromMaterial);} 

  Bool_t  IsFromHIJING() const { return (fGeneratorIndex == 0); }
  void    SetGeneratorIndex(Short_t i) {fGeneratorIndex = i;}
  Short_t GetGeneratorIndex() const {return fGeneratorIndex;}

  static void     GetEtaPhiDiff(const AliVTrack *t, const AliVCluster *v, Double_t &phidiff, Double_t &etadiff);
  static Byte_t   GetTrackType(const AliVTrack *t);
  static Byte_t   GetTrackType(const AliAODTrack *aodTrack, UInt_t filterBit1, UInt_t filterBit2);

 protected:
  Double32_t      fPt;             //[0,0,12]   pt at vertex
  Double32_t      fEta;            //[-1,1,12]  eta at vertex
  Double32_t      fPhi;            //[0,6.3,12] phi at vertex
  Double32_t      fM;              //           mass
  Byte_t          fQ;              //           charge
  Int_t           fLabel;          //           label  
  Byte_t          fTrackType;      //           0=global track; 1=w/o SPD, w/ ITS refit; 2=w/o SPD, w/o ITS refit
  Double32_t      fEtaEmc;         //[-1,1,12]  eta at emcal surface
  Double32_t      fPhiEmc;         //[0,6.3,12] phi at emcal surface
  Double32_t      fPtEmc;          //[0,0,12]   pt at emcal surface
  Bool_t          fEmcal;          //           is true if track propagated to emcal
  UInt_t          fFlag;           //!          Flag for indication of primary etc (MC)
  Short_t         fGeneratorIndex; //!          Index of generator in cocktail
  Short_t         fClusId;         //!          cluster id of matched cluster; -1 if not set
  AliVTrack      *fOrig;           //!          ptr to original track

  ClassDef(AliPicoTrack, 8) // Pico track class
};
#endif
