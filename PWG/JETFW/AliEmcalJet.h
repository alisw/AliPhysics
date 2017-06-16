#ifndef ALIEMCALJET_H
#define ALIEMCALJET_H
/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <vector>
#include <algorithm>
#include <utility>

#include <iosfwd>
#include <TArrayI.h>
#include <TMath.h>
#include <TClonesArray.h>
#include <TVector2.h>
#include <TLorentzVector.h>
#include <TString.h>

#include <AliVParticle.h>
#include <AliVCluster.h>
#include <AliVEvent.h>

#include "AliEmcalJetShapeProperties.h"

/**
 * @class AliEmcalJet
 * @brief Represent a jet reconstructed using the EMCal jet framework
 * @ingroup EMCALCOREFW
 * @author Salvatore Aiola <salvatore.aiola@yale.edu>, Yale University
 * @author Constantin Loizides <cloizides@lbl.gov>, Lawrence Berkeley National Laboratory
 *
 * This class encapsulates a jet reconstructed using the EMCal jet framework.
 * It can represent charged (tracks), neutral (EMCal clusters) or full (tracks+clusters) jet,
 * or a particle-level jet reconstructed from a Monte Carlo simulation.
 * Information contained in the class includes:
 * - reconstructed jet axis in cylindrical coordinates (eta, phi, pT);
 * - jet area (used for background subtraction);
 * - jet constituents (cluster, tracks, particles);
 * - flavor tagging;
 * - jet shape properties;
 * - matching with other reconstructed jets (e.g. detector level with particole level).
 * The class implements also a number of service function to calculate other observable,
 * such as fragmentation functions, subtracted jet momentum, etc.
 */
class AliEmcalJet : public AliVParticle
{
 public:
  
  /**
   * @enum JetAcceptanceType
   * @brief Bit definition for jet geometry acceptance. Cut implemented in AliJetContainer
   * by comparing jet's bits (set in jet finder) to container's bits (set by user).
   * If user doesn't set jet acceptance cut value, no cut is performed (equivalent to kUser).
   * The boundaries defined for each bit should be taken as approximate (within a couple
   * cells) -- the user should verify the definitions if precision is crucial.
   * If you create jets outside of the standard jet finder, you may have to manually set these
   * acceptance bits if you want to use the acceptance selection cut in the jet container
   * e.g. "jet->SetJetAcceptanceType(fJetTask->FindJetAcceptanceType(eta,phi,r));".
   */
  enum JetAcceptanceType {
    kTPC              = 1<<0,     ///< TPC acceptance
    kTPCfid           = 1<<1,     ///< TPC fiducial acceptance (each eta edge narrowed by jet R)
    kEMCAL            = 1<<2,     ///< EMCal acceptance
    kEMCALfid         = 1<<3,     ///< EMCal fiducial acceptance (each eta, phi edge narrowed by jet R)
    kDCAL             = 1<<4,     ///< DCal acceptance -- spans entire rectangular region in eta-phi (including most of PHOS)
    kDCALfid          = 1<<5,     ///< DCal fiducial acceptance (each eta, phi edge narrowed by jet R)
    kDCALonly         = 1<<6,     ///< DCal acceptance -- spans ONLY DCal (no PHOS or gap)
    kDCALonlyfid      = 1<<7,     ///< DCal fiducial acceptance (each eta, phi edge narrowed by jet R)
    kPHOS             = 1<<8,     ///< PHOS acceptance
    kPHOSfid          = 1<<9,     ///< PHOS fiducial acceptance (each eta, phi edge narrowed by jet R)
    kUser             = 1<<10     ///< Full acceptance, i.e. no acceptance cut applied -- left to user
  };
  
  /**
   * @enum EFlavourTag
   * @brief Bit definition for the flavor tagging
   */
  enum EFlavourTag {
    kDStar   = 1<<0,    ///< Jet is tagged to contain a D* meson
    kD0      = 1<<1,    ///< Jet is tagged to contain a D0 meson
    kSig1    = 1<<2,    ///< Generic signal 1
    kSig2    = 1<<3,    ///< Generic signal 2
    kBckgrd1 = 1<<4,    ///< Generic background 1
    kBckgrd2 = 1<<5,    ///< Generic background 2
    kBckgrd3 = 1<<6     ///< Generic background 3
  };

  AliEmcalJet();
  AliEmcalJet(Double_t px, Double_t py, Double_t pz);
  AliEmcalJet(Double_t pt, Double_t eta, Double_t phi, Double_t m);
  AliEmcalJet(const AliEmcalJet &jet);
  AliEmcalJet& operator=(const AliEmcalJet &jet);
  virtual ~AliEmcalJet();
  friend std::ostream &operator<<(std::ostream &in, const AliEmcalJet &jet);
  Int_t Compare(const TObject* obj)  const;
  std::ostream &Print(std::ostream &in) const;
  TString toString() const;

  // Implementation of AliVParticle interface
  Double_t          Px()                         const { return fPt*TMath::Cos(fPhi) ; }
  Double_t          Py()                         const { return fPt*TMath::Sin(fPhi) ; }
  Double_t          Pz()                         const { return fPt*TMath::SinH(fEta); }
  Double_t          Pt()                         const { return fPt                  ; }
  Double_t          P()                          const { return fPt*TMath::CosH(fEta); }
  Bool_t            PxPyPz(Double_t p[3])        const { p[0]=Px();p[1]=Py();p[2]=Pz(); return kTRUE; }
  Double_t          Xv()                         const { return 0.;      }
  Double_t          Yv()                         const { return 0.;      }
  Double_t          Zv()                         const { return 0.;      }
  Bool_t            XvYvZv(Double_t x[3])        const { x[0]=0;x[1]=0;x[2]=0         ; return kTRUE; }
  Double_t          OneOverPt()                  const { return 1./fPt ; }
  Double_t          Phi()                        const { return fPhi   ; }
  Double_t          Theta()                      const { return 2*TMath::ATan(TMath::Exp(-fEta))      ; }
  Double_t          E()                          const { Double_t p=P(); return TMath::Sqrt(fM*fM+p*p); }
  Double_t          M()                          const { return fM     ; }
  Double_t          Eta()                        const { return fEta   ; }
  Double_t          Y()                          const { Double_t e = E(); Double_t pz = Pz(); return 0.5*TMath::Log((e+pz)/(e-pz));    }
  Short_t           Charge()                     const { return 0      ; }
  Int_t             GetLabel()                   const { return fLabel ; }
  Int_t             PdgCode()                    const { return 0;       }
  const Double_t   *PID()                        const { return 0;       }

  // Other kinematic and jet properties
  Double_t          Phi_0_2pi()                  const { return TVector2::Phi_0_2pi(fPhi); }
  Double_t          Area()                       const { return fArea                    ; }
  Double_t          AreaPt()                     const { return fArea                    ; }
  Double_t          AreaEta()                    const { return fAreaEta                 ; }
  Double_t          AreaPhi()                    const { return fAreaPhi                 ; }
  Double_t          AreaE()                      const { return fAreaE                   ; }
  Double_t          AreaEmc()                    const { return fAreaEmc                 ; }
  Bool_t            AxisInEmcal()                const { return fAxisInEmcal             ; }
  Int_t             ClusterAt(Int_t idx)         const { return fClusterIDs.At(idx)      ; }
  UShort_t          GetNumberOfClusters()        const { return fClusterIDs.GetSize()    ; }
  UShort_t          GetNumberOfTracks()          const { return fTrackIDs.GetSize()      ; }
  UShort_t          GetNumberOfConstituents()    const { return GetNumberOfClusters()+GetNumberOfTracks(); }
  Double_t          FracEmcalArea()              const { return fAreaEmc/fArea           ; }
  Bool_t            IsInsideEmcal()              const { return (fAreaEmc/fArea>0.999)   ; }
  Bool_t            IsInEmcal()                  const { return (Bool_t)(fAreaEmc > 0)   ; }
  Bool_t            IsMC()                       const { return (Bool_t)(MCPt() > 0)     ; }
  Bool_t            IsSortable()                 const { return kTRUE                    ; }
  Double_t          MaxNeutralPt()               const { return fMaxNPt                  ; }
  Double_t          MaxChargedPt()               const { return fMaxCPt                  ; }
  Double_t          NEF()                        const { return fNEF                     ; }
  UShort_t          Nn()                         const { return fNn                      ; }
  UShort_t          Nch()                        const { return fNch                     ; }
  UShort_t          N()                          const { return Nch()+Nn()               ; }
  Int_t             NEmc()                       const { return fNEmc                    ; }
  Double_t          MCPt()                       const { return fMCPt                    ; }
  Double_t          MaxClusterPt()               const { return MaxNeutralPt()           ; }
  Double_t          MaxTrackPt()                 const { return MaxChargedPt()           ; }
  Double_t          MaxPartPt()                  const { return fMaxCPt < fMaxNPt ? fMaxNPt : fMaxCPt; }
  Double_t          PtEmc()                      const { return fPtEmc                   ; }
  Double_t          PtSub()                      const { return fPtSub                   ; }
  Double_t          PtSubVect()                  const { return fPtSubVect               ; }
  Int_t             TrackAt(Int_t idx)           const { return fTrackIDs.At(idx)        ; }

  // Background subtraction
  Double_t          PtSub(Double_t rho, Bool_t save = kFALSE)          ;
  Double_t          PtSubVect(Double_t rho, Bool_t save = kFALSE)      ;
  TLorentzVector    SubtractRhoVect(Double_t rho, Bool_t save = kFALSE);

  // Jet constituents
  AliVCluster      *Cluster(Int_t idx)                                             const;
  AliVCluster      *ClusterAt(Int_t idx, TClonesArray *ca)                         const;
  Int_t             ContainsCluster(AliVCluster* cluster, TClonesArray* clusters)  const;
  Int_t             ContainsCluster(Int_t ic)                                      const;
  AliVCluster      *GetLeadingCluster(TClonesArray *clusters)                      const;
  AliVParticle     *Track(Int_t idx)                                               const;
  AliVParticle     *TrackAt(Int_t idx, TClonesArray *ta)                           const;
  Int_t             ContainsTrack(AliVParticle* track, TClonesArray* tracks)       const;
  Int_t             ContainsTrack(Int_t it)                                        const;
  AliVParticle     *GetLeadingTrack(TClonesArray *tracks)                          const;
  Bool_t            IsGhost()                                                      const { return fNn + fNch == 0; }

  // Fragmentation function
  Double_t          GetZ(const Double_t trkPx, const Double_t trkPy, const Double_t trkPz)  const;
  Double_t          GetZ(const AliVParticle* trk )                                          const;
  Double_t          GetXi(const AliVParticle* trk )                                         const;
  Double_t          GetXi(const Double_t trkPx, const Double_t trkPy, const Double_t trkPz) const;

  // Other service methods
  void              GetMomentum(TLorentzVector &vec)                                        const;
  Double_t          DeltaR(const AliVParticle* part)                                        const;

  // Setters
  void              SetLabel(Int_t l)                  { fLabel   = l;                     }
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
  void              SetNEmc(Int_t n)                   { fNEmc           = n;              }
  void              SetPtEmc(Double_t pt)              { fPtEmc          = pt;             }
  void              SetPtSub(Double_t ps)              { fPtSub          = ps;             }
  void              SetPtSubVect(Double_t ps)          { fPtSubVect      = ps;             }
  void              AddClusterAt(Int_t clus, Int_t idx){ fClusterIDs.AddAt(clus, idx);     }
  void              AddTrackAt(Int_t track, Int_t idx) { fTrackIDs.AddAt(track, idx);      }
  void              Clear(Option_t */*option*/="");

  // Sorting methods
  void              SortConstituents();
  std::vector<int>  GetPtSortedTrackConstituentIndexes(TClonesArray *tracks) const;

  // Trigger
  Bool_t            IsTriggerJet(UInt_t trigger=AliVEvent::kEMCEJE) const   { return (Bool_t)((fTriggers & trigger) != 0); }
  void              SetTrigger(UInt_t trigger)                              { fTriggers  = trigger;                        }
  void              AddTrigger(UInt_t trigger)                              { fTriggers |= trigger;                        }

  // Matching
  void              ResetMatching();
  void              SetClosestJet(AliEmcalJet *j, Double_t d)       { fClosestJets[0] = j; fClosestJetsDist[0] = d    ; }
  void              SetSecondClosestJet(AliEmcalJet *j, Double_t d) { fClosestJets[1] = j; fClosestJetsDist[1] = d    ; }
  void              SetMatchedToClosest(UShort_t m)                 { fMatched        = 0; fMatchingType       = m    ; }
  void              SetMatchedToSecondClosest(UShort_t m)           { fMatched        = 1; fMatchingType       = m    ; }
  AliEmcalJet*      ClosestJet()                              const { return fClosestJets[0]                          ; }
  Double_t          ClosestJetDistance()                      const { return fClosestJetsDist[0]                      ; }
  AliEmcalJet*      SecondClosestJet()                        const { return fClosestJets[1]                          ; }
  Double_t          SecondClosestJetDistance()                const { return fClosestJetsDist[1]                      ; }
  AliEmcalJet*      MatchedJet()                              const { return fMatched < 2 ? fClosestJets[fMatched] : 0; }
  UShort_t          GetMatchingType()                         const { return fMatchingType                            ; }

  // Jet tagging
  void              SetTaggedJet(AliEmcalJet *j)                    { fTaggedJet = j                                  ; }
  void              SetTagStatus(Int_t i)                           { fTagStatus = i                                  ; }
  AliEmcalJet*      GetTaggedJet()                            const { return fTaggedJet                               ; }
  Int_t             GetTagStatus()                            const { return fTagStatus                               ; }

  // Ghosts
  void AddGhost(const Double_t dPx, const Double_t dPy, const Double_t dPz, const Double_t dE);
  Bool_t HasGhost() const                               { return fHasGhost; }
  const std::vector<TLorentzVector> GetGhosts()   const { return fGhosts  ; }

  // Debug printouts
  void Print(Option_t* /*opt*/ = "") const;
  void PrintConstituents(TClonesArray* tracks, TClonesArray* clusters) const;

  //heavy-flavor jets
  Int_t             GetFlavour()                 const { return fFlavourTagging;                       }
  void              AddFlavourTag(Int_t tag)           { fFlavourTagging |= tag;                       }
  void              SetFlavour(Int_t flavour)          { fFlavourTagging = flavour;                    }
  Bool_t            TestFlavourTag(Int_t tag)    const { return (Bool_t)((tag & fFlavourTagging) !=0); }
  void              ClearFlavourTracks()               { if (fFlavourTracks) fFlavourTracks->Clear(); }
  void              AddFlavourTrack(AliVParticle* hftrack);
  AliVParticle     *GetFlavourTrack(Int_t i=0)       const;
  Double_t          GetFlavourTrackZ(Int_t i=0)      const;
  AliVParticle     *RemoveFlavourTrack(Int_t i=0)         ;
  
  // Jet shape
  AliEmcalJetShapeProperties* GetShapeProperties() const{ return fJetShapeProperties; }
  AliEmcalJetShapeProperties* GetShapeProperties() { if (!fJetShapeProperties) CreateShapeProperties(); return fJetShapeProperties; }
  void CreateShapeProperties() { if (fJetShapeProperties) delete fJetShapeProperties; fJetShapeProperties = new AliEmcalJetShapeProperties(); }
  
  // Jet geometrical acceptance
  void     SetJetAcceptanceType(UInt_t type)         { fJetAcceptanceType = type;}
  UInt_t   GetJetAcceptanceType()              const { return fJetAcceptanceType; }

 protected:
  /// Jet transverse momentum
  Double32_t        fPt;                  //[0,0,12]
  /// Jet pseudo-rapidity
  Double32_t        fEta;                 //[-1,1,12]
  /// Jet axis azimuthal angle
  Double32_t        fPhi;                 //[0,6.3,12]
  /// Jet mass
  Double32_t        fM;                   //[0,0,8]
  /// Jet Neutral Energy Fraction
  Double32_t        fNEF;                 //[0,1,8]
  /// Jet transverse area
  Double32_t        fArea;                //[0,0,12]
  /// Jet eta area
  Double32_t        fAreaEta;             //[0,0,12]
  /// Jet phi area
  Double32_t        fAreaPhi;             //[0,0,12]
  /// Jet temporal area component
  Double32_t        fAreaE;               //[0,0,12]
  /// Area on EMCAL surface (determined by ghosts in EMCal acceptance)
  Double32_t        fAreaEmc;             //[0,0,12]
  Bool_t            fAxisInEmcal;         ///<  Whether the jet axis is inside the EMCAL acceptance
  Int_t             fFlavourTagging;      ///<  Tag jet with a flavor (use bits defined in enum EFlavourTag)
  TObjArray        *fFlavourTracks;       //!<! Heavy flavor candidate tracks matched to the jet
  /// Pt of maximum charged constituent
  Double32_t        fMaxCPt;              //[0,0,12]
  /// Pt of maximum neutral constituent
  Double32_t        fMaxNPt;              //[0,0,12]
  Double32_t        fMCPt;                ///<  Pt from MC particles contributing to the jet
  Int_t             fNn;                  ///<  Number of neutral constituents
  Int_t             fNch;                 ///<  Number of charged constituents
  /// Pt in EMCAL acceptance
  Double32_t        fPtEmc;               //[0,0,12]
  Int_t             fNEmc;                ///<  Number of constituents in EMCAL acceptance
  TArrayI           fClusterIDs;          ///<  Array containing ids of cluster constituents
  TArrayI           fTrackIDs;            ///<  Array containing ids of track constituents
  AliEmcalJet      *fClosestJets[2];      //!<! If this is MC it contains the two closest detector level jets in order of distance and viceversa
  Double32_t        fClosestJetsDist[2];  //!<! Distance from the two closest jets
  UShort_t          fMatched;             //!<! 0 or 1 if it is matched with one of the closest jets; 2 if it is not matched
  UShort_t          fMatchingType;        //!<! Matching type
  AliEmcalJet      *fTaggedJet;           //!<! Jet tagged to this jet
  Int_t             fTagStatus;           //!<! Status of tagging -1: NA 0: not tagged 1: tagged
  Double_t          fPtSub;               //!<! Background subtracted pt (not stored set from outside)
  Double_t          fPtSubVect;           //!<! Background vector subtracted pt (not stored set from outside)
  UInt_t            fTriggers;            //!<! Triggers that the jet might have fired (AliVEvent::EOfflineTriggerTypes)
  Int_t             fLabel;               //!<! Label to inclusive jet for constituent subtracted jet

  Bool_t            fHasGhost;            //!<! Whether ghost particle are included within the constituents
  std::vector<TLorentzVector> fGhosts;    //!<! Vector containing the ghost particles

  AliEmcalJetShapeProperties *fJetShapeProperties; //!<! Pointer to the jet shape properties
  UInt_t fJetAcceptanceType;    //!<!  Jet acceptance type (stored bitwise)

 private:
  /**
   * @struct sort_descend
   * @brief Simple C structure to allow sorting in descending order
   */
  struct sort_descend {
    // first value of the pair is Pt and the second is entry index
    bool operator () (const std::pair<Double_t, Int_t>& p1, const std::pair<Double_t, Int_t>& p2)  { return p1.first > p2.first ; }
  };

  /// \cond CLASSIMP
  ClassDef(AliEmcalJet,19);
  /// \endcond
};

std::ostream &operator<<(std::ostream &in, const AliEmcalJet &jet);
#endif
