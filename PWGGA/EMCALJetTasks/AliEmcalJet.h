#ifndef AliEmcalJet_H
#define AliEmcalJet_H

// $Id$

#include <TArrayS.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TClonesArray.h>

#include "AliVParticle.h"
#include "AliVCluster.h"

class AliEmcalJet : public AliVParticle
{
 public:
  AliEmcalJet();
  AliEmcalJet(Double_t px, Double_t py, Double_t pz);
  AliEmcalJet(Double_t pt, Double_t eta, Double_t phi, Double_t m);
  AliEmcalJet(const AliEmcalJet &jet); 
  AliEmcalJet& operator=(const AliEmcalJet &jet);

  Double_t          Px()                         const { return fPt*TMath::Cos(fPhi);  }
  Double_t          Py()                         const { return fPt*TMath::Cos(fPhi);  }
  Double_t          Pz()                         const { return fPt*TMath::SinH(fEta); }
  Double_t          Pt()                         const { return fPt;                   }
  Double_t          P()                          const { return fPt*TMath::CosH(fEta); }
  Bool_t            PxPyPz(Double_t p[3])        const { p[0]=Px();p[1]=Py();p[2]=Pz(); return 1;         }
  Double_t          Xv()                         const { return 0.;      }
  Double_t          Yv()                         const { return 0.;      }
  Double_t          Zv()                         const { return 0.;      }
  Bool_t            XvYvZv(Double_t x[3])        const { x[0]=0;x[1]=0;x[2]=0; return 1;                  }
  Double_t          OneOverPt()                  const { return 1./fPt;  }
  Double_t          Phi()                        const { return fPhi;    }
  Double_t          Theta()                      const { return 2*TMath::ATan(TMath::Exp(-fEta));         }
  Double_t          E()                          const { Double_t p=P(); return TMath::Sqrt(M()*M()+p*p); }
  Double_t          M()                          const { return 0.13957; }
  Double_t          Eta()                        const { return fEta;    }
  Double_t          Y()                          const { return 0.5*TMath::Log((E()+Pz())/(E()-Pz()));    }
  Short_t           Charge()                     const { return 0;       }
  Int_t             GetLabel()                   const { return -1;      }
  Int_t             PdgCode()                    const { return 0;       }
  const Double_t   *PID()                        const { return 0;       }
  void              GetMom(TLorentzVector &vec)  const;
  void              Print(Option_t* option = "") const;

  Double_t          Area()                       const { return fArea;                     }
  Double_t          AreaEmc()                    const { return fAreaEmc;                  }
  Bool_t            AxisInEmcal()                const { return fAxisInEmcal;              }
  UShort_t          GetNumberOfClusters()        const { return fClusterIDs.GetSize();     }
  Short_t           ClusterAt(Int_t idx)         const { return fClusterIDs.At(idx);       }
  AliVCluster      *ClusterAt(Int_t idx, TClonesArray *clusarray)  const { return dynamic_cast<AliVCluster*>(clusarray->At(ClusterAt(idx))); }
  UShort_t          GetNumberOfTracks()          const { return fTrackIDs.GetSize();       }
  Short_t           TrackAt(Int_t idx)           const { return fTrackIDs.At(idx);         }
  AliVParticle     *TrackAt(Int_t idx, TClonesArray *trackarray)   const { return dynamic_cast<AliVParticle*>(trackarray->At(TrackAt(idx))); } 
  Double_t          FracEmcalArea()              const { return fAreaEmc/fArea;            }
  Bool_t            IsInsideEmcal()              const { return (fAreaEmc/fArea>0.999);    }
  Bool_t            IsInEmcal()                  const { return (fAreaEmc>0);              }
  Double_t          MaxNeutralPt()               const { return fMaxNPt;                   }
  Double_t          MaxChargedPt()               const { return fMaxCPt;                   }
  Double_t          NEF()                        const { return fNEF;                      }
  UShort_t          Nn()                         const { return fNn;                       }
  UShort_t          Nch()                        const { return fNch;                      }
  UShort_t          N()                          const { return Nch()+Nn();                }
  Double_t          MCPt()                       const { return fMCPt;                     }
  Bool_t            IsMC()                       const { return (Bool_t)(MCPt() > E() - MCPt()); }
  AliEmcalJet*      ClosestJet()                 const { return fClosestJets[0];           }
  Double_t          ClosestJetDistance()         const { return fClosestJetsDist[0];       }
  AliEmcalJet*      SecondClosestJet()           const { return fClosestJets[1];           }
  Double_t          SecondClosestJetDistance()   const { return fClosestJetsDist[1];       }
  AliEmcalJet*      MatchedJet()                 const { return fMatched < 2 ? fClosestJets[fMatched] : 0; }
  // Only for backward compatibility
  Double_t          MaxClusterPt()               const { return MaxNeutralPt();            }
  Double_t          MaxTrackPt()                 const { return MaxChargedPt();            }

  void              AddClusterAt(Int_t clus, Int_t idx){ fClusterIDs.AddAt(clus, idx);     }
  void              AddTrackAt(Int_t track, Int_t idx) { fTrackIDs.AddAt(track, idx);      }
  void              Clear(Option_t */*option*/="")     { fClusterIDs.Set(0); 
                                                         fTrackIDs.Set(0);                 }
  void              SetArea(Double_t a)                { fArea    = a;                     }
  void              SetAreaEmc(Double_t a)             { fAreaEmc = a;                     }
  void              SetAxisInEmcal(Bool_t b)           { fAxisInEmcal = b;                 }
  void              SetMaxNeutralPt(Double32_t t)      { fMaxNPt  = t;                     }
  void              SetMaxChargedPt(Double32_t t)      { fMaxCPt  = t;                     }
  void              SetNEF(Double_t nef)               { fNEF     = nef;                   }
  void              SetNumberOfClusters(Int_t n)       { fClusterIDs.Set(n);               }
  void              SetNumberOfTracks(Int_t n)         { fTrackIDs.Set(n);                 }
  void              SetNumberOfCharged(Int_t n)        { fNch = n;                         }
  void              SetNumberOfNeutrals(Int_t n)       { fNn = n;                          }
  void              SetMCPt(Double_t p)                { fMCPt = p;                        }
  void              SortConstituents();
  void              SetClosestJet(AliEmcalJet *j, Double_t d)       { fClosestJets[0] = j; fClosestJetsDist[0] = d; }
  void              SetSecondClosestJet(AliEmcalJet *j, Double_t d) { fClosestJets[1] = j; fClosestJetsDist[1] = d; }
  void              SetMatchedToClosest()                           { fMatched = 0;        }
  void              SetMatchedToSecondClosest()                     { fMatched = 1;        }

 protected:
  Double32_t   fPt;                  //[0,0,12]   pt 
  Double32_t   fEta;                 //[-1,1,12]  eta
  Double32_t   fPhi;                 //[0,6.3,12] phi
  Double32_t   fM;                   //[0,0,8]    mass
  Double32_t   fNEF;                 //[0,1,8]    neutral energy fraction
  Double32_t   fArea;                //[0,0,12]   area
  Double32_t   fAreaEmc;             //[0,0,12]   area on EMCAL surface (determined from ghosts)
  Bool_t       fAxisInEmcal;         //           =true if jet axis inside EMCAL acceptance
  Double32_t   fMaxCPt;              //[0,0,12]   pt of maximum charged constituent
  Double32_t   fMaxNPt;              //[0,0,12]   pt of maximum neutral constituent
  Double32_t   fMCPt;                //           pt from MC particles contributing to the jet
  Int_t        fNn;                  //           number of neutral constituents
  Int_t        fNch;                 //           number of charged constituents
  TArrayS      fClusterIDs;          //           array of cluster constituents  
  TArrayS      fTrackIDs;            //           array of track constituents   
  AliEmcalJet *fClosestJets[2];      //!          if this is MC it contains the two closest detector level jets in order of distance and viceversa
  Double32_t   fClosestJetsDist[2];  //!          distance to closest jets (see above)
  UShort_t     fMatched;             //!          0,1 if it is matched with one of the closest jets; 2 if it is not matched

  ClassDef(AliEmcalJet,6) // Emcal jet class in cylindrical coordinates
};
#endif
