#ifndef AliEmcalJet_H
#define AliEmcalJet_H

#include <Riosfwd.h>
#include <vector>
#include <algorithm>
#include <utility>

#include <TArrayS.h>
#include <TMath.h>
#include <TClonesArray.h>
#include <TVector2.h>
#include <TLorentzVector.h>
#include <TString.h>

#include "AliVParticle.h"
#include "AliVCluster.h"
#include "AliVEvent.h"

class AliEmcalJet : public AliVParticle
{
 public:
     enum EFlavourTag{
       kDStar = 1<<0,
       kD0 = 1<<1,
       kSig1 = 1<<2,
       kSig2 = 1<<3,
       kBckgrd1 = 1<<4,
       kBckgrd2 = 1<<5,
       kBckgrd3 = 1<<6
       //.....
    };

  AliEmcalJet();
  AliEmcalJet(Double_t px, Double_t py, Double_t pz);
  AliEmcalJet(Double_t pt, Double_t eta, Double_t phi, Double_t m);
  AliEmcalJet(const AliEmcalJet &jet);
  AliEmcalJet& operator=(const AliEmcalJet &jet);
  friend std::ostream &operator<<(std::ostream &in, const AliEmcalJet &jet);

  std::ostream &Print(std::ostream &in) const;
  TString toString() const;

  Double_t          Px()                         const { return fPt*TMath::Cos(fPhi);  }
  Double_t          Py()                         const { return fPt*TMath::Sin(fPhi);  }
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
  Double_t          Phi_0_2pi()                  const { return TVector2::Phi_0_2pi(fPhi);    }
  Double_t          Theta()                      const { return 2*TMath::ATan(TMath::Exp(-fEta));         }
  Double_t          E()                          const { Double_t p=P(); return TMath::Sqrt(fM*fM+p*p); }
  Double_t          M()                          const { return fM; }
  Double_t          Eta()                        const { return fEta;    }
  Double_t          Y()                          const { Double_t e = E(); Double_t pz = Pz(); return 0.5*TMath::Log((e+pz)/(e-pz));    }
  Short_t           Charge()                     const { return 0;       }
  Int_t             GetLabel()                   const { return fLabel;  }
  Int_t             PdgCode()                    const { return 0;       }
  const Double_t   *PID()                        const { return 0;       }
  void              GetMom(TLorentzVector &vec)  const {GetMomentum(vec);}
  void              GetMomentum(TLorentzVector &vec) const;

  Double_t          Area()                       const { return fArea;                     }
  Double_t          AreaPt()                     const { return fArea;                     }
  Double_t          AreaEta()                    const { return fAreaEta;                  }
  Double_t          AreaPhi()                    const { return fAreaPhi;                  }
  Double_t          AreaE()                      const { return fAreaE;                    }
  Double_t          AreaEmc()                    const { return fAreaEmc;                  }
  Bool_t            AxisInEmcal()                const { return fAxisInEmcal;              }
  Int_t             Compare(const TObject* obj)  const;
  Short_t           ClusterAt(Int_t idx)         const { return fClusterIDs.At(idx);       }
  AliVCluster      *ClusterAt(Int_t idx, TClonesArray *ca)  const { if (!ca) return 0; return dynamic_cast<AliVCluster*>(ca->At(ClusterAt(idx))); }
  Int_t             ContainsCluster(AliVCluster* cluster, TClonesArray* clusters) const { return clusters==NULL||cluster==NULL ? 0 : ContainsCluster(clusters->IndexOf(cluster)); }
  Int_t             ContainsCluster(Int_t ic)    const;
  AliVCluster      *GetLeadingCluster(TClonesArray *clusters) const;
  UShort_t          GetNumberOfClusters()        const { return fClusterIDs.GetSize();     }
  UShort_t          GetNumberOfTracks()          const { return fTrackIDs.GetSize();       }
  UShort_t          GetNumberOfConstituents()    const { return GetNumberOfClusters()+GetNumberOfTracks();       }
  Double_t          FracEmcalArea()              const { return fAreaEmc/fArea;            }
  Bool_t            IsInsideEmcal()              const { return (fAreaEmc/fArea>0.999);    }
  Bool_t            IsInEmcal()                  const { return (fAreaEmc>0);              }
  Bool_t            IsMC()                       const { return (Bool_t)(MCPt() > 0);      }
  Bool_t            IsSortable()                 const { return kTRUE;                     }
  Double_t          MaxNeutralPt()               const { return fMaxNPt;                   }
  Double_t          MaxChargedPt()               const { return fMaxCPt;                   }
  Double_t          NEF()                        const { return fNEF;                      }
  UShort_t          Nn()                         const { return fNn;                       }
  UShort_t          Nch()                        const { return fNch;                      }
  UShort_t          N()                          const { return Nch()+Nn();                }
  Int_t             NEmc()                       const { return fNEmc;                     }
  Double_t          MCPt()                       const { return fMCPt;                     }
  Double_t          MaxClusterPt()               const { return MaxNeutralPt();            }
  Double_t          MaxTrackPt()                 const { return MaxChargedPt();            }
  Double_t          MaxPartPt()                  const { return fMaxCPt < fMaxNPt ? fMaxNPt : fMaxCPt;     }
  Double_t          PtEmc()                      const { return fPtEmc;                    }
  Double_t          PtSub()                      const { return fPtSub;                    }
  Double_t          PtSubVect()                  const { return fPtSubVect;                }
  Double_t          PtSub(Double_t rho, Bool_t save = kFALSE);
  Double_t          PtSubVect(Double_t rho, Bool_t save = kFALSE);
  TLorentzVector    SubtractRhoVect(Double_t rho, Bool_t save = kFALSE);
  Short_t           TrackAt(Int_t idx)           const { return fTrackIDs.At(idx);         }
  AliVParticle     *TrackAt(Int_t idx, TClonesArray *ta)  const { if (!ta) return 0; return dynamic_cast<AliVParticle*>(ta->At(TrackAt(idx))); }
  Int_t             ContainsTrack(AliVParticle* track, TClonesArray* tracks) const { return tracks==NULL||track==NULL ? 0 : ContainsTrack(tracks->IndexOf(track)); }
  Int_t             ContainsTrack(Int_t it)      const;
  AliVParticle     *GetLeadingTrack(TClonesArray *tracks) const;


  void              AddClusterAt(Int_t clus, Int_t idx){ fClusterIDs.AddAt(clus, idx);     }
  void              AddTrackAt(Int_t track, Int_t idx) { fTrackIDs.AddAt(track, idx);      }
  void              Clear(Option_t */*option*/="")     { fClusterIDs.Set(0); fTrackIDs.Set(0); fClosestJets[0] = 0; fClosestJets[1] = 0;
                                                         fClosestJetsDist[0] = 0; fClosestJetsDist[1] = 0; fMatched = 0; fPtSub = 0;
                                                         fGhosts.clear(); fHasGhost = kFALSE; }
  Double_t          DeltaR(const AliVParticle* part) const;
  Double_t          GetZ ( const Double_t trkPx, const Double_t trkPy, const Double_t trkPz ) const; // Get Z of constituent trk
  Double_t          GetZ ( const AliVParticle* trk )       const; // Get Z of constituent trk
  Double_t          GetXi ( const AliVParticle* trk )      const { return TMath::Log ( 1/GetZ (trk) ); } // Get Xi of constituent trk
  Double_t          GetXi ( const Double_t trkPx, const Double_t trkPy, const Double_t trkPz ) const { return TMath::Log ( 1/GetZ (trkPx, trkPy, trkPz ) ); } // Get Xi of constituent trk

  void              SetLabel(Int_t l)                  { fLabel = l;                       }
  void              SetArea(Double_t a)                { fArea    = a;                     }
  void              SetAreaEta(Double_t a)             { fAreaEta = a;                     }
  void              SetAreaPhi(Double_t a)             { fAreaPhi = TVector2::Phi_0_2pi(a); }
  void              SetAreaE(Double_t a)               { fAreaE = a;                       }
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
  std::vector<int>  SortConstituentsPt(TClonesArray *tracks) const;
  void              SetNEmc(Int_t n)                   { fNEmc           = n;              }
  void              SetPtEmc(Double_t pt)              { fPtEmc          = pt;             }
  void              SetPtSub(Double_t ps)              { fPtSub          = ps;             }
  void              SetPtSubVect(Double_t ps)          { fPtSubVect      = ps;             }

  // Trigger
  Bool_t            IsTriggerJet(UInt_t trigger=AliVEvent::kEMCEJE) const   { return (Bool_t)((fTriggers & trigger) != 0); }
  void              SetTrigger(UInt_t trigger)                              { fTriggers  = trigger;                        }
  void              AddTrigger(UInt_t trigger)                              { fTriggers |= trigger;                        }

  // Matching
  void              SetClosestJet(AliEmcalJet *j, Double_t d)       { fClosestJets[0] = j; fClosestJetsDist[0] = d    ; }
  void              SetSecondClosestJet(AliEmcalJet *j, Double_t d) { fClosestJets[1] = j; fClosestJetsDist[1] = d    ; }
  void              SetMatchedToClosest(UShort_t m)                 { fMatched        = 0; fMatchingType       = m    ; }
  void              SetMatchedToSecondClosest(UShort_t m)           { fMatched        = 1; fMatchingType       = m    ; }
  void              ResetMatching();
  AliEmcalJet*      ClosestJet()                              const { return fClosestJets[0]                          ; }
  Double_t          ClosestJetDistance()                      const { return fClosestJetsDist[0]                      ; }
  AliEmcalJet*      SecondClosestJet()                        const { return fClosestJets[1]                          ; }
  Double_t          SecondClosestJetDistance()                const { return fClosestJetsDist[1]                      ; }
  AliEmcalJet*      MatchedJet()                              const { return fMatched < 2 ? fClosestJets[fMatched] : 0; }
  UShort_t          GetMatchingType()                         const { return fMatchingType                            ; }

  void              SetTaggedJet(AliEmcalJet *j)                    { fTaggedJet = j                                  ; }
  void              SetTagStatus(Int_t i)                           { fTagStatus = i                                  ; }
  AliEmcalJet*      GetTaggedJet()                            const { return fTaggedJet                               ; }
  Int_t             GetTagStatus()                            const { return fTagStatus                               ; }

  //jet shape derivatives
  //jet mass
  void              SetFirstDerivative(Double_t d)                  { fJetShapeMassFirstDer = d                       ; }
  void              SetSecondDerivative(Double_t d)                 { fJetShapeMassSecondDer = d                      ; }
  void              SetFirstOrderSubtracted(Double_t d)             { fJetShapeMassFirstSub = d                       ; }
  void              SetSecondOrderSubtracted(Double_t d)            { fJetShapeMassSecondSub = d                      ; }
  Double_t          GetFirstDerivative()                      const { return fJetShapeMassFirstDer                    ; }
  Double_t          GetSecondDerivative()                     const { return fJetShapeMassSecondDer                   ; }
  Double_t          GetFirstOrderSubtracted()                 const { return fJetShapeMassFirstSub                    ; }
  Double_t          GetSecondOrderSubtracted()                const { return fJetShapeMassSecondSub                   ; }

  //jet structure function
  TArrayF           GetGRNumerator()                          const { return fGRNumerator                             ; }
  TArrayF           GetGRDenominator()                        const { return fGRDenominator                           ; }
  TArrayF           GetGRNumeratorSub()                       const { return fGRNumeratorSub                          ; }
  TArrayF           GetGRDenominatorSub()                     const { return fGRDenominatorSub                        ; }
  void              AddGRNumAt(Float_t num, Int_t idx)              { fGRNumerator.AddAt(num, idx)                    ; }
  void              AddGRDenAt(Float_t den, Int_t idx)              { fGRDenominator.AddAt(den, idx)                  ; }
  void              SetGRNumSize(UInt_t s)                          { fGRNumerator.Set(s)                             ; }
  void              SetGRDenSize(UInt_t s)                          { fGRDenominator.Set(s)                           ; }

  void              AddGRNumSubAt(Float_t num, Int_t idx)           { fGRNumeratorSub.AddAt(num, idx)                 ; }
  void              AddGRDenSubAt(Float_t den, Int_t idx)           { fGRDenominatorSub.AddAt(den, idx)               ; }
  void              SetGRNumSubSize(UInt_t s)                       { fGRNumeratorSub.Set(s)                          ; }
  void              SetGRDenSubSize(UInt_t s)                       { fGRDenominatorSub.Set(s)                        ; }
  void              PrintGR();

  //Angularity
  void              SetFirstDerivativeAngularity(Double_t d)        { fJetShapeAngularityFirstDer = d                 ; }
  void              SetSecondDerivativeAngularity(Double_t d)       { fJetShapeAngularitySecondDer = d                ; }
  void              SetFirstOrderSubtractedAngularity(Double_t d)   { fJetShapeAngularityFirstSub = d                 ; }
  void              SetSecondOrderSubtractedAngularity(Double_t d)  { fJetShapeAngularitySecondSub = d                ; }
  Double_t          GetFirstDerivativeAngularity()            const { return fJetShapeAngularityFirstDer              ; }
  Double_t          GetSecondDerivativeAngularity()           const { return fJetShapeAngularitySecondDer             ; }
  Double_t          GetFirstOrderSubtractedAngularity()       const { return fJetShapeAngularityFirstSub              ; }
  Double_t          GetSecondOrderSubtractedAngularity()      const { return fJetShapeAngularitySecondSub             ; }

  //pTD
  void              SetFirstDerivativepTD(Double_t d)               { fJetShapepTDFirstDer = d                        ; }
  void              SetSecondDerivativepTD(Double_t d)              { fJetShapepTDSecondDer = d                       ; }
  void              SetFirstOrderSubtractedpTD(Double_t d)          { fJetShapepTDFirstSub = d                        ; }
  void              SetSecondOrderSubtractedpTD(Double_t d)         { fJetShapepTDSecondSub = d                       ; }
  Double_t          GetFirstDerivativepTD()                   const { return fJetShapepTDFirstDer                     ; }
  Double_t          GetSecondDerivativepTD()                  const { return fJetShapepTDSecondDer                    ; }
  Double_t          GetFirstOrderSubtractedpTD()              const { return fJetShapepTDFirstSub                     ; }
  Double_t          GetSecondOrderSubtractedpTD()             const { return fJetShapepTDSecondSub                    ; }

  //Circularity
  void              SetFirstDerivativeCircularity(Double_t d)       { fJetShapeCircularityFirstDer = d                ; }
  void              SetSecondDerivativeCircularity(Double_t d)      { fJetShapeCircularitySecondDer = d               ; }
  void              SetFirstOrderSubtractedCircularity(Double_t d)  { fJetShapeCircularityFirstSub = d                ; }
  void              SetSecondOrderSubtractedCircularity(Double_t d) { fJetShapeCircularitySecondSub = d               ; }
  Double_t          GetFirstDerivativeCircularity()           const { return fJetShapeCircularityFirstDer             ; }
  Double_t          GetSecondDerivativeCircularity()          const { return fJetShapeCircularitySecondDer            ; }
  Double_t          GetFirstOrderSubtractedCircularity()      const { return fJetShapeCircularityFirstSub             ; }
  Double_t          GetSecondOrderSubtractedCircularity()     const { return fJetShapeCircularitySecondSub            ; }

  //Sigma2
  void              SetFirstDerivativeSigma2(Double_t d)            { fJetShapeSigma2FirstDer = d                     ; }
  void              SetSecondDerivativeSigma2(Double_t d)           { fJetShapeSigma2SecondDer = d                    ; }
  void              SetFirstOrderSubtractedSigma2(Double_t d)       { fJetShapeSigma2FirstSub = d                     ; }
  void              SetSecondOrderSubtractedSigma2(Double_t d)      { fJetShapeSigma2SecondSub = d                    ; }
  Double_t          GetFirstDerivativeSigma2()           const      { return fJetShapeSigma2FirstDer                  ; }
  Double_t          GetSecondDerivativeSigma2()          const      { return fJetShapeSigma2SecondDer                 ; }
  Double_t          GetFirstOrderSubtractedSigma2()      const      { return fJetShapeSigma2FirstSub                  ; }
  Double_t          GetSecondOrderSubtractedSigma2()     const      { return fJetShapeSigma2SecondSub                 ; }


  //number of contituents
  void              SetFirstDerivativeConstituent(Double_t d)       { fJetShapeConstituentFirstDer = d                ; }
  void              SetSecondDerivativeConstituent(Double_t d)      { fJetShapeConstituentSecondDer = d               ; }
  void              SetFirstOrderSubtractedConstituent(Double_t d)  { fJetShapeConstituentFirstSub = d                ; }
  void              SetSecondOrderSubtractedConstituent(Double_t d) { fJetShapeConstituentSecondSub = d               ; }
  Double_t          GetFirstDerivativeConstituent()           const { return fJetShapeConstituentFirstDer             ; }
  Double_t          GetSecondDerivativeConstituent()          const { return fJetShapeConstituentSecondDer            ; }
  Double_t          GetFirstOrderSubtractedConstituent()      const { return fJetShapeConstituentFirstSub             ; }
  Double_t          GetSecondOrderSubtractedConstituent()     const { return fJetShapeConstituentSecondSub            ; }

  //leading minus subleading constituent
  void              SetFirstDerivativeLeSub(Double_t d)             { fJetShapeLeSubFirstDer = d                      ; }
  void              SetSecondDerivativeLeSub(Double_t d)            { fJetShapeLeSubSecondDer = d                     ; }
  void              SetFirstOrderSubtractedLeSub(Double_t d)        { fJetShapeLeSubFirstSub = d                      ; }
  void              SetSecondOrderSubtractedLeSub(Double_t d)       { fJetShapeLeSubSecondSub = d                     ; }
  Double_t          GetFirstDerivativeLeSub()                 const { return fJetShapeLeSubFirstDer                   ; }
  Double_t          GetSecondDerivativeLeSub()                const { return fJetShapeLeSubSecondDer                  ; }
  Double_t          GetFirstOrderSubtractedLeSub()            const { return fJetShapeLeSubFirstSub                   ; }
  Double_t          GetSecondOrderSubtractedLeSub()           const { return fJetShapeLeSubSecondSub                  ; }

  //heavy-flavor jets
  Int_t             GetFlavour()                 const { return fFlavourTagging;                       }
  void              AddFlavourTag(Int_t tag)           { fFlavourTagging |= tag;                       }
  void              SetFlavour(Int_t flavour)          { fFlavourTagging = flavour;                    }
  Bool_t            TestFlavourTag(Int_t tag)    const { return (Bool_t)((tag & fFlavourTagging) !=0); }
  void              AddFlavourTrack(AliVParticle* hftrack){ if (!fFlavourTracks) fFlavourTracks = new TObjArray(); fFlavourTracks->Add(hftrack);               }
  AliVParticle     *GetFlavourTrack(Int_t i=0)   const { return fFlavourTracks && i >= 0 && i < fFlavourTracks->GetEntriesFast() ? static_cast<AliVParticle*>(fFlavourTracks->At(i)) : 0; }
  Double_t          GetFlavourTrackZ(Int_t i=0)  const;
  AliVParticle     *RemoveFlavourTrack(Int_t i=0)      { return fFlavourTracks && i >= 0 && i < fFlavourTracks->GetEntriesFast() ? static_cast<AliVParticle*>(fFlavourTracks->RemoveAt(i)) : 0; }
  void              ClearFlavourTracks()               { if (fFlavourTracks) fFlavourTracks->Clear(); }
  
  void AddGhost(const Double_t dPx, const Double_t dPy, const Double_t dPz, const Double_t dE) {
    TLorentzVector ghost(dPx, dPy, dPz, dE);
    fGhosts.push_back(ghost);
    if (!fHasGhost) fHasGhost = kTRUE;
    return;
  }

  Bool_t HasGhost() const { return fHasGhost; }
  const std::vector<TLorentzVector> GetGhosts() const { return fGhosts; }

  void Print(Option_t* /*opt*/ = "") const;
  void PrintConstituents(TClonesArray* tracks, TClonesArray* clusters) const;

 protected:
  Double32_t        fPt;                  //[0,0,12]   pt
  Double32_t        fEta;                 //[-1,1,12]  eta
  Double32_t        fPhi;                 //[0,6.3,12] phi
  Double32_t        fM;                   //[0,0,8]    mass
  Double32_t        fNEF;                 //[0,1,8]    neutral energy fraction
  Double32_t        fArea;                //[0,0,12]   area (transverse)
  Double32_t        fAreaEta;             //[0,0,12]   area eta
  Double32_t        fAreaPhi;             //[0,0,12]   area phi
  Double32_t        fAreaE;               //[0,0,12]   temporal area component
  Double32_t        fAreaEmc;             //[0,0,12]   area on EMCAL surface (determined from ghosts)
  Bool_t            fAxisInEmcal;         //           =true if jet axis inside EMCAL acceptance
  Int_t             fFlavourTagging;      //           tag jet with a flavour, bit 0 = no tag; bit 1= Dstar; bit 2 = D0
  TObjArray        *fFlavourTracks;       //!          heavy flavour candidate tracks matched to the jet
  Double32_t        fMaxCPt;              //[0,0,12]   pt of maximum charged constituent
  Double32_t        fMaxNPt;              //[0,0,12]   pt of maximum neutral constituent
  Double32_t        fMCPt;                //           pt from MC particles contributing to the jet
  Int_t             fNn;                  //           number of neutral constituents
  Int_t             fNch;                 //           number of charged constituents
  Double32_t        fPtEmc;               //[0,0,12]   pt in EMCAL acceptance
  Int_t             fNEmc;                //           number of constituents in EMCAL acceptance
  TArrayI           fClusterIDs;          //           array containing ids of cluster constituents
  TArrayI           fTrackIDs;            //           array containing ids of track constituents
  AliEmcalJet      *fClosestJets[2];      //!          if this is MC it contains the two closest detector level jets in order of distance and viceversa
  Double32_t        fClosestJetsDist[2];  //!          distance to closest jets (see above)
  UShort_t          fMatched;             //!          0,1 if it is matched with one of the closest jets; 2 if it is not matched
  UShort_t          fMatchingType;        //!          matching type
  AliEmcalJet      *fTaggedJet;           //!          jet tagged to this jet
  Int_t             fTagStatus;           //!          status of tagging -1: NA 0: not tagged 1: tagged
  Double_t          fPtSub;               //!          background subtracted pt (not stored set from outside)
  Double_t          fPtSubVect;           //!          background vector subtracted pt (not stored set from outside)
  UInt_t            fTriggers;            //!          triggers that the jet might have fired (AliVEvent::EOfflineTriggerTypes)

  Double_t          fJetShapeMassFirstDer;         //!   result from shape derivatives for jet mass: 1st derivative
  Double_t          fJetShapeMassSecondDer;        //!   result from shape derivatives for jet mass: 2nd derivative
  Double_t          fJetShapeMassFirstSub;         //!   result from shape derivatives for jet mass: 1st order subtracted
  Double_t          fJetShapeMassSecondSub;        //!   result from shape derivatives for jet mass: 2nd order subtracted
  Int_t             fLabel;                        //    label to inclusive jet for constituent subtracted jet

  TArrayF           fGRNumerator;                  //!   array with angular structure function numerator
  TArrayF           fGRDenominator;                //!   array with angular structure function denominator
  TArrayF           fGRNumeratorSub;               //!   array with angular structure function numerator
  TArrayF           fGRDenominatorSub;             //!   array with angular structure function denominator

  Double_t          fJetShapeAngularityFirstDer;   //!   result from shape derivatives for jet Angularity: 1st derivative
  Double_t          fJetShapeAngularitySecondDer;  //!   result from shape derivatives for jet Angularity: 2nd derivative
  Double_t          fJetShapeAngularityFirstSub;   //!   result from shape derivatives for jet Angularity: 1st order subtracted
  Double_t          fJetShapeAngularitySecondSub;  //!   result from shape derivatives for jet Angularity: 2nd order subtracted

  Double_t          fJetShapepTDFirstDer;          //!   result from shape derivatives for jet pTD: 1st derivative
  Double_t          fJetShapepTDSecondDer;         //!   result from shape derivatives for jet pTD: 2nd derivative
  Double_t          fJetShapepTDFirstSub;          //!   result from shape derivatives for jet pTD: 1st order subtracted
  Double_t          fJetShapepTDSecondSub;         //!   result from shape derivatives for jet pTD: 2nd order subtracted

  Double_t          fJetShapeCircularityFirstDer;  //!   result from shape derivatives for jet circularity: 1st derivative
  Double_t          fJetShapeCircularitySecondDer; //!   result from shape derivatives for jet circularity: 2nd derivative
  Double_t          fJetShapeCircularityFirstSub;  //!   result from shape derivatives for jet circularity: 1st order subtracted
  Double_t          fJetShapeCircularitySecondSub; //!   result from shape derivatives for jetcircularity: 2nd order subtracted

  Double_t          fJetShapeSigma2FirstDer;       //!   result from shape derivatives for jet sigma2: 1st derivative
  Double_t          fJetShapeSigma2SecondDer;      //!   result from shape derivatives for jet sigma2: 2nd derivative
  Double_t          fJetShapeSigma2FirstSub;       //!   result from shape derivatives for jet sigma2: 1st order subtracted
  Double_t          fJetShapeSigma2SecondSub;      //!   result from shape derivatives for jetsigma2: 2nd order subtracted

  Double_t          fJetShapeConstituentFirstDer;  //!   result from shape derivatives for jet const: 1st derivative
  Double_t          fJetShapeConstituentSecondDer; //!   result from shape derivatives for jet const: 2nd derivative
  Double_t          fJetShapeConstituentFirstSub;  //!   result from shape derivatives for jet const: 1st order subtracted
  Double_t          fJetShapeConstituentSecondSub; //!   result from shape derivatives for jet const: 2nd order subtracted

  Double_t          fJetShapeLeSubFirstDer;        //!   result from shape derivatives for jet LeSub: 1st derivative
  Double_t          fJetShapeLeSubSecondDer;       //!   result from shape derivatives for jet LeSub: 2nd derivative
  Double_t          fJetShapeLeSubFirstSub;        //!   result from shape derivatives for jet LeSub: 1st order subtracted
  Double_t          fJetShapeLeSubSecondSub;       //!   result from shape derivatives for jet LeSub: 2nd order subtracted

  Bool_t fHasGhost;
  std::vector<TLorentzVector> fGhosts;

  private:
    struct sort_descend
        { // sort in decreasing order
          // first value of the pair is Pt and the second is entry index
        bool operator () (const std::pair<Double_t, Int_t>& p1, const std::pair<Double_t, Int_t>& p2)  { return p1.first > p2.first ; }
        };

  ClassDef(AliEmcalJet,17) // Emcal jet class in cylindrical coordinates
};

std::ostream &operator<<(std::ostream &in, const AliEmcalJet &jet);
#endif
