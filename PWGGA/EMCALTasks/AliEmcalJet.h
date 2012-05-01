#ifndef AliEmcalJet_H
#define AliEmcalJet_H

// $Id$

#include <TArrayS.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include "AliVParticle.h"

class AliEmcalJet : public AliVParticle
{
 public:
  AliEmcalJet() : AliVParticle(), fPt(0), fEta(0), fPhi(0), fM(0), fNEF(0), fArea(0), 
                  fMaxCPt(0), fMaxNPt(0), fClusterIDs(), fTrackIDs() {;}
  AliEmcalJet(Double_t pt, Double_t eta, Double_t phi, Double_t m);
  AliEmcalJet(Double_t px, Double_t py, Double_t pz);
  AliEmcalJet(const AliEmcalJet &jet); 
  AliEmcalJet& operator=(const AliEmcalJet &jet);

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

  Double_t    Area()                       const { return fArea;                 }
  Short_t     ClusterAt(Int_t idx)         const { return fClusterIDs.At(idx);   }
  UShort_t    GetNumberOfClusters()        const { return Nn();                  }
  UShort_t    GetNumberOfTracks()          const { return Nch();                 }
  Double_t    MaxClusterPt()               const { return fMaxNPt;               }
  Double_t    MaxTrackPt()                 const { return fMaxCPt;               }
  Double_t    NEF()                        const { return fNEF;                  }
  UShort_t    Nn()                         const { return fClusterIDs.GetSize(); }
  UShort_t    Nch()                        const { return fTrackIDs.GetSize();   }
  UShort_t    N()                          const { return Nch()+Nn();            }
  Short_t     TrackAt(Int_t idx)           const { return fTrackIDs.At(idx);     }
  void        AddClusterAt(Int_t clus, Int_t idx){ fClusterIDs.AddAt(clus, idx); }
  void        AddTrackAt(Int_t track, Int_t idx) { fTrackIDs.AddAt(track, idx);  }
  void        SetArea(Double_t a)                { fArea   = a;                  }
  void        SetMaxClusterPt(Double32_t t)      { fMaxNPt = t;                  }
  void        SetMaxTrackPt(Double32_t t)        { fMaxCPt = t;                  }
  void        SetNEF(Double_t nef)               { fNEF    = nef;                }
  void        SetNumberOfClusters(Int_t n)       { fClusterIDs.Set(n);           }
  void        SetNumberOfTracks(Int_t n)         { fTrackIDs.Set(n);             }

  void        SortConstituents();

 protected:
  Double32_t  fPt;           //[0,0,12]   pt 
  Double32_t  fEta;          //[-1,1,12]  eta
  Double32_t  fPhi;          //[0,6.3,12] phi
  Double32_t  fM;            //[0,0,8]    mass
  Double32_t  fNEF;          //[0,1,8]    neutral energy fraction
  Double32_t  fArea;         //[0,0,12]   area
  Double32_t  fMaxCPt;       //[0,0,12]   pt of maximum track
  Double32_t  fMaxNPt;       //[0,0,12]   pt of maximum cluster
  TArrayS     fClusterIDs;   //           array of cluster constituents  
  TArrayS     fTrackIDs;     //           array of track constituents    

  ClassDef(AliEmcalJet,2) // ESD jet class in cylindrical coordinates
};
#endif
