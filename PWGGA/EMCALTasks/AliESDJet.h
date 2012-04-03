#ifndef AliESDJet_H
#define AliESDJet_H

// $Id$

#include "AliVParticle.h"
#include <TLorentzVector.h>
#include <TMath.h>

class AliESDJet : public AliVParticle
{
 public:
  AliESDJet() : AliVParticle(), fPt(0), fEta(0), fPhi(0), fM(0), fNEF(0), fArea(0), 
                fNch(0), fNn(0), fMaxCPt(0), fMaxNPt(0) {;}
  AliESDJet(Double_t pt, Double_t eta, Double_t phi, Double_t m);
  AliESDJet(Double_t px, Double_t py, Double_t pz);
  AliESDJet(const AliESDJet &jet); 
  AliESDJet& operator=(const AliESDJet &jet);

  Double_t    Px()                         const { return fPt*TMath::Cos(fPhi);  }
  Double_t    Py()                         const { return fPt*TMath::Cos(fPhi);  }
  Double_t    Pz()                         const { return fPt*TMath::SinH(fEta); }
  Double_t    Pt()                         const { return fPt;                   }
  Double_t    P()                          const { return fPt*TMath::CosH(fEta); }
  Bool_t      PxPyPz(Double_t p[3])        const { p[0]=Px();p[1]=Py();p[2]=Pz(); return 1;         }
  Double_t    Xv()                         const { return 0.;      }
  Double_t    Yv()                         const { return 0.;      }
  Double_t    Zv()                         const { return 0.;      }
  Bool_t      XvYvZv(Double_t x[3])        const { x[0]=0;x[1]=0;x[2]=0; return 1;                  }
  Double_t    OneOverPt()                  const { return 1./fPt;  }
  Double_t    Phi()                        const { return fPhi;    }
  Double_t    Theta()                      const { return 2*TMath::ATan(TMath::Exp(-fEta));         }
  Double_t    E()                          const { Double_t p=P(); return TMath::Sqrt(M()*M()+p*p); }
  Double_t    M()                          const { return 0.13957; }
  Double_t    Eta()                        const { return fEta;    }
  Double_t    Y()                          const { return 0.5*TMath::Log((E()+Pz())/(E()-Pz()));    }
  Short_t     Charge()                     const { return 0;       }
  Int_t       GetLabel()                   const { return -1;      }
  Int_t       PdgCode()                    const { return 0;       }
  const Double_t *PID()                    const { return 0;       }
  void        GetMom(TLorentzVector &vec)  const;
  void        Print(Option_t* option = "") const;

  Double_t    Area()                       const { return fArea;   }
  Double_t    NEF()                        const { return fNEF;    }
  UShort_t    N()                          const { return fNch+fNn;}
  Double32_t  MaxTrackPt()                 const { return fMaxCPt; }
  Double32_t  MaxClusterPt()               const { return fMaxNPt; }

  void        SetArea(Double_t a)                { fArea   = a;    }
  void        SetNEF(Double_t nef)               { fNEF    = nef;  }
  void        SetMaxTrackPt(Double32_t t)        { fMaxCPt = t;    }
  void        SetMaxClusterPt(Double32_t t)      { fMaxNPt = t;    }

 protected:
  Double32_t  fPt;           //[0,0,12]   pt 
  Double32_t  fEta;          //[-1,1,12]  eta
  Double32_t  fPhi;          //[0,6.3,12] phi
  Double32_t  fM;            //[0,0,8]    mass
  Double32_t  fNEF;          //[0,1,8]    neutral energy fraction
  Double32_t  fArea;         //[0,0,12]   area
  UShort_t    fNch;          //           number of charged constituents
  UShort_t    fNn;           //           number of neutral constituents
  Double32_t  fMaxCPt;       //[0,0,12]   pt of maximum track
  Double32_t  fMaxNPt;       //[0,0,12]   pt of maximum cluster

  ClassDef(AliESDJet,1) // ESD jet class in cylindrical coordinates
};
#endif
