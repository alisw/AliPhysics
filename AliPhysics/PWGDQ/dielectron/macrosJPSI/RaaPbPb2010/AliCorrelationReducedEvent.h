// Classes used for creating a reduced information tree
// Author: Ionut-Cristian Arsene (i.c.arsene@gsi.de)
// 
//  Basic structure:
//  1. Event wise information
//  2. List of tracks in the event
//  3. List of resonance candidates

#ifndef ALICORRELATIONREDUCEDEVENT_H
#define ALICORRELATIONREDUCEDEVENT_H

#include <TClonesArray.h>
#include <TBits.h>
#include <TMath.h>


const Int_t fgkNMaxHarmonics = 10;
/*class AliCorrelationReducedTrack;
class AliCorrelationReducedPair;
class AliCorrelationReducedEventFriend;
class AliCorrelationReducedEvent;
class AliCorrelationReducedCaloCluster;*/

//_____________________________________________________________________
class AliCorrelationReducedTrack : public TObject {

  friend class AliAnalysisTaskCorrelationTree;

 public:
  AliCorrelationReducedTrack();
  ~AliCorrelationReducedTrack();

  // getters
  UShort_t TrackId()                     const {return fTrackId;}
  ULong_t  Status()                      const {return fStatus;}
  Bool_t   CheckTrackStatus(UInt_t flag) const {return (flag<8*sizeof(ULong_t) ? (fStatus&(1<<flag)) : kFALSE);}
  Int_t    Charge()                      const {return (fGlobalPt>0.0 ? +1 : -1);}
  Float_t  Px()                          const {return TMath::Abs(fGlobalPt)*TMath::Cos(fGlobalPhi);}
  Float_t  Py()                          const {return TMath::Abs(fGlobalPt)*TMath::Sin(fGlobalPhi);}
  Float_t  Pz()                          const {return TMath::Abs(fGlobalPt)*TMath::SinH(fGlobalEta);}
  Float_t  P()                           const {return TMath::Abs(fGlobalPt)*TMath::CosH(fGlobalEta);};
  Float_t  Phi()                         const {return fGlobalPhi;}
  Float_t  Pt()                          const {return TMath::Abs(fGlobalPt);}
  Float_t  Eta()                         const {return fGlobalEta;}
  Float_t  Theta()                       const {return TMath::ACos(TMath::TanH(fGlobalEta));}
  Float_t  PxTPC()                       const {return fTPCPt*TMath::Cos(fTPCPhi);}
  Float_t  PyTPC()                       const {return fTPCPt*TMath::Sin(fTPCPhi);}
  Float_t  PzTPC()                       const {return fTPCPt*TMath::SinH(fTPCEta);}
  Float_t  PTPC()                        const {return fTPCPt*TMath::CosH(fTPCEta);};
  Float_t  PhiTPC()                      const {return fTPCPhi;}
  Float_t  PtTPC()                       const {return fTPCPt;}
  Float_t  EtaTPC()                      const {return fTPCEta;}
  Float_t  ThetaTPC()                    const {return TMath::ACos(TMath::TanH(fTPCEta));}
  Float_t  Pin()                         const {return fMomentumInner;}
  Float_t  DCAxy()                       const {return fDCA[0];}
  Float_t  DCAz()                        const {return fDCA[1];}
  
  UShort_t ITSncls()                const;
  UChar_t  ITSclusterMap()          const {return fITSclusterMap;}
  Bool_t   ITSLayerHit(Int_t layer) const {return (layer>=0 && layer<6 ? (fITSclusterMap&(1<<layer)) : kFALSE);};
  Float_t  ITSsignal()              const {return fITSsignal;}
  
  UChar_t TPCncls()                        const {return fTPCNcls;}
  UChar_t TPCFindableNcls()                const {return fTPCNclsF;}
  UChar_t TPCCrossedRows()                 const {return fTPCCrossedRows;}
  UChar_t TPCnclsIter1()                   const {return fTPCNclsIter1;}
  UChar_t TPCClusterMap()                  const {return fTPCClusterMap;}
  Int_t   TPCClusterMapBitsFired()         const;
  Bool_t  TPCClusterMapBitFired(Int_t bit) const {return (bit>=0 && bit<8 ? (fTPCClusterMap&(1<<bit)) : kFALSE);};
  Float_t TPCsignal()                      const {return fTPCsignal;}
  Float_t TPCnSig(Int_t specie)            const {return (specie>=0 && specie<=3 ? fTPCnSig[specie] : -999.);}
  
  Float_t  TOFbeta()             const {return fTOFbeta;}    
  Float_t  TOFnSig(Int_t specie) const {return (specie>=0 && specie<=3 ? fTOFnSig[specie] : -999.);}
  
  Int_t    TRDntracklets(Int_t type)  const {return (type==0 || type==1 ? fTRDntracklets[type] : -1);}
  Float_t  TRDpid(Int_t specie)       const {return (specie>=0 && specie<=1 ? fTRDpid[specie] : -999.);}
  
  Int_t    CaloClusterId() const {return fCaloClusterId;}
  //Float_t  CaloEnergy(AliCorrelationReducedEvent* event) const {if(fCaloClusterId>0) return event->GetCaloCluster(fCaloClusterId)->Energy();}
  //Float_t  CaloDx(AliCorrelationReducedEvent* event) const {if(fCaloClusterId>0) return event->GetCaloCluster(fCaloClusterId)->Dx();}
  //Float_t  CaloDz(AliCorrelationReducedEvent* event) const {if(fCaloClusterId>0) return event->GetCaloCluster(fCaloClusterId)->Dz();}
  
  Float_t  BayesPID(Int_t specie) const {return (specie>=0 && specie<=2 ? fBayesPID[specie] : -999.);}
  Bool_t   UsedForQvector()       const {return fFlags&(1<<0);}
  
 private:
  UShort_t fTrackId;            // track id 
  ULong_t fStatus;              // tracking status
  Float_t fGlobalPhi;           // phi at the vertex from global track, in the [0,2pi) interval
  Float_t fGlobalPt;            // pt*charge at the vertex from global track
  Float_t fGlobalEta;           // eta at the vertex from global track
  Float_t fTPCPhi;              // phi at the vertex from TPC alone tracking , in the [0,2pi) interval
  Float_t fTPCPt;               // pt at the vertex from TPC alone tracking  
  Float_t fTPCEta;              // eta at the vertex from TPC alone tracking 
  Float_t fMomentumInner;       // inner param momentum (only the magnitude)
  Float_t fDCA[2];              // DCA xy,z
  
  // ITS
  UChar_t  fITSclusterMap;      // ITS cluster map
  Float_t  fITSsignal;          // ITS signal
  
  // TPC
  UChar_t fTPCNcls;            // TPC ncls                          
  UChar_t fTPCCrossedRows;     // TPC crossed rows                  
  UChar_t fTPCNclsF;           // TPC findable ncls                 
  UChar_t fTPCNclsIter1;       // TPC no clusters after first iteration
  UChar_t fTPCClusterMap;      // TPC cluster distribution map
  Float_t fTPCsignal;          // TPC de/dx
  Float_t fTPCnSig[4];         // 0-electron; 1-pion; 2-kaon; 3-proton
    
  // TOF
  Float_t fTOFbeta;             // TOF pid info
  Float_t fTOFnSig[4];          // TOF n-sigma deviation from expected signal
  
  // TRD
  UChar_t fTRDntracklets[2];       // 0 - AliESDtrack::GetTRDntracklets(); 1 - AliESDtrack::GetTRDntrackletsPID()   TODO: use only 1 char
  Float_t fTRDpid[2];              // TRD pid probabilities, [0]-electron, [1]-pion
  
  // EMCAL/PHOS
  Int_t  fCaloClusterId;          // ID for the calorimeter cluster (if any)
  
  // Bayesian PID
  Float_t fBayesPID[3];            // Combined Bayesian PID   pi/K/p
  
  UShort_t fFlags;                // BIT0 toggled if track used for TPC event plane   TODO combine with other posible flags, use for MC pid?
  // TODO flag for which TPC part used for pid  --> Char_t  used in 2011 data
  
  AliCorrelationReducedTrack(const AliCorrelationReducedTrack &c);
  AliCorrelationReducedTrack& operator= (const AliCorrelationReducedTrack &c);

  ClassDef(AliCorrelationReducedTrack, 2);
};


//_____________________________________________________________________
class AliCorrelationReducedPair : public TObject {

  friend class AliAnalysisTaskCorrelationTree;

 public:
  enum CandidateType {
    kK0sToPiPi=0,
    kPhiToKK,
    kLambda0ToPPi,
    kALambda0ToPPi,
    kJpsiToEE,
    kNMaxCandidateTypes
  };
  AliCorrelationReducedPair();
  AliCorrelationReducedPair(const AliCorrelationReducedPair &c);
  ~AliCorrelationReducedPair();

  // getters
  Char_t   CandidateId()         const {return fCandidateId;}
  Char_t   PairType()            const {return fPairType;}
  Int_t    LegId(Int_t leg)      const {return (leg==0 || leg==1 ? fLegIds[leg] : -1);}
  Float_t  Mass(Int_t idx=0)     const {return (idx>=0 && idx<3 ? fMass[idx] : -999.);}
  Float_t  Px()                  const {return fPt*TMath::Cos(fPhi);}
  Float_t  Py()                  const {return fPt*TMath::Sin(fPhi);}
  Float_t  Pz()                  const {return fPt*TMath::SinH(fEta);}
  Float_t  P()                   const {return fPt*TMath::CosH(fEta);}
  Float_t  Phi()                 const {return fPhi;}
  Float_t  Pt()                  const {return fPt;}
  Float_t  Eta()                 const {return fEta;}
  Float_t  Energy()              const;
  Float_t  Rapidity()            const;
  Float_t  Theta()               const {return TMath::ACos(TMath::TanH(fEta));}
  Float_t  Lxy()                 const {return fLxy;}
  Float_t  OpeningAngle()        const {return fOpeningAngle;}
  Bool_t   IsOnTheFly()          const {return fPairType;}
  UInt_t   MCid()                const {return fMCid;}
  Bool_t   CheckMC(const Int_t flag) const {return (flag<32 ? (fMCid&(1<<flag)) : kFALSE);}
  
 private:
  Char_t  fCandidateId;         // candidate type (K0s, Lambda, J/psi, phi, etc)
  Char_t  fPairType;            // 0 ++; 1 +-; 2 -- for dielectron pairs; 0- offline, 1- on the fly for V0 candidates
  UShort_t fLegIds[2];          // leg ids 
  Float_t fMass[3];             // invariant mass for pairs (2 extra mass values for other V0 pid assumptions)
                                // idx=0 -> K0s assumption; idx=1 -> Lambda; idx=2 -> anti-Lambda
  Float_t fPhi;                 // pair phi in the [0,2*pi) interval
  Float_t fPt;                  // pair pt
  Float_t fEta;                 // pair eta 
  Float_t fLxy;                 // pseudo-proper decay length
  Float_t fOpeningAngle;        // opening angle                TODO remove   ???
  UInt_t  fMCid;                // Bit map with Monte Carlo info about the pair

  AliCorrelationReducedPair& operator= (const AliCorrelationReducedPair &c);

  ClassDef(AliCorrelationReducedPair, 1);
};


//_________________________________________________________________________
class AliCorrelationReducedEventFriend : public TObject {
  
  friend class AliAnalysisTaskCorrelationTree;
  
 public: 
  enum EventPlaneStatus {
    kRaw=0,
    kCalibrated,
    kRecentered,
    kShifted,
    kNMaxFlowFlags
  };
  enum EventPlaneDetector {
    kTPC=0,       
    kTPCptWeights,
    kTPCpos,
    kTPCneg,
    kVZEROA,
    kVZEROC,
    kFMD,
    kZDCA,
    kZDCC,
    kNdetectors
  };
  
  AliCorrelationReducedEventFriend();
  ~AliCorrelationReducedEventFriend();
  
  Double_t Qx(Int_t det, Int_t harmonic)  const {return (det>=0 && det<kNdetectors && harmonic>0 && harmonic<=fgkNMaxHarmonics ? fQvector[det][harmonic-1][0] : -999.);}
  Double_t Qy(Int_t det, Int_t harmonic)  const {return (det>=0 && det<kNdetectors && harmonic>0 && harmonic<=fgkNMaxHarmonics ? fQvector[det][harmonic-1][1] : -999.);}
  Double_t EventPlane(Int_t det, Int_t h) const;
  UChar_t GetEventPlaneStatus(Int_t det, Int_t h) const {return (det>=0 && det<kNdetectors && h>0 && h<=fgkNMaxHarmonics ? fEventPlaneStatus[det][h] : -999.);} 
  Bool_t  CheckEventPlaneStatus(Int_t det, Int_t h, EventPlaneStatus flag) const;
  void    CopyEvent(AliCorrelationReducedEventFriend* event);

  void SetQx(Int_t det, Int_t harmonic, Float_t qx) { if(det>=0 && det<kNdetectors && harmonic>0 && harmonic<=fgkNMaxHarmonics) fQvector[det][harmonic-1][0]=qx;}
  void SetQy(Int_t det, Int_t harmonic, Float_t qy) { if(det>=0 && det<kNdetectors && harmonic>0 && harmonic<=fgkNMaxHarmonics) fQvector[det][harmonic-1][1]=qy;}
  void SetEventPlaneStatus(Int_t det, Int_t harmonic, EventPlaneStatus status) { 
    if(det>=0 && det<kNdetectors && harmonic>0 && harmonic<=fgkNMaxHarmonics) 
      fEventPlaneStatus[det][harmonic-1] |= (1<<status);
  }
  
 private:
  // Q-vectors for the first 10 harmonics from TPC, VZERO, FMD and ZDC detectors
  Double_t fQvector[kNdetectors][fgkNMaxHarmonics][2];     // Q vector components for all detectors and 6 harmonics
  UChar_t fEventPlaneStatus[kNdetectors][fgkNMaxHarmonics];  // Bit maps for the event plane status (1 char per detector and per harmonic)
   
  void ClearEvent();
  AliCorrelationReducedEventFriend(const AliCorrelationReducedEventFriend &c);
  AliCorrelationReducedEventFriend& operator= (const AliCorrelationReducedEventFriend &c);

  ClassDef(AliCorrelationReducedEventFriend, 1);
};


//_________________________________________________________________________
class AliCorrelationReducedCaloCluster : public TObject {
  
  friend class AliAnalysisTaskCorrelationTree;
  
 public:
  enum ClusterType {
    kUndefined=0, kEMCAL, kPHOS  
  };
   
  AliCorrelationReducedCaloCluster();
  ~AliCorrelationReducedCaloCluster();
  
  Bool_t  IsEMCAL() const {return (fType==kEMCAL ? kTRUE : kFALSE);}
  Bool_t  IsPHOS()  const {return (fType==kPHOS ? kTRUE : kFALSE);}
  Float_t Energy()  const {return fEnergy;}
  Float_t Dx()      const {return fTrackDx;}
  Float_t Dz()      const {return fTrackDz;}
  
 private:
  Char_t  fType;         // cluster type (EMCAL/PHOS)
  Float_t fEnergy;       // cluster energy
  Float_t fTrackDx;      // distance to closest track in phi
  Float_t fTrackDz;      // distance to closest track in z
  
  AliCorrelationReducedCaloCluster(const AliCorrelationReducedCaloCluster &c);
  AliCorrelationReducedCaloCluster& operator= (const AliCorrelationReducedCaloCluster &c);

  ClassDef(AliCorrelationReducedCaloCluster, 1);
};


//_________________________________________________________________________
class AliCorrelationReducedEvent : public TObject {

  friend class AliAnalysisTaskCorrelationTree;

 public:
  AliCorrelationReducedEvent();
  ~AliCorrelationReducedEvent();

  // getters
  Int_t     RunNo()                           const {return fRunNo;}
  UShort_t  BC()                              const {return fBC;}
  ULong64_t TriggerMask()                     const {return fTriggerMask;}
  Bool_t    IsPhysicsSelection()              const {return fIsPhysicsSelection;}
  Float_t   Vertex(Int_t axis)                const {return (axis>=0 && axis<=2 ? fVtx[axis] : 0);}
  Int_t     VertexNContributors()             const {return fNVtxContributors;}
  Float_t   VertexTPC(Int_t axis)             const {return (axis>=0 && axis<=2 ? fVtxTPC[axis] : 0);}
  Int_t     VertexTPCContributors()           const {return fNVtxTPCContributors;}
  Float_t   CentralityVZERO()                 const {return fCentrality[0];}
  Float_t   CentralitySPD()                   const {return fCentrality[1];}
  Float_t   CentralityTPC()                   const {return fCentrality[2];}
  Float_t   CentralityZEMvsZDC()              const {return fCentrality[3];}
  Int_t     CentralityQuality()               const {return fCentQuality;}
  Int_t     NV0CandidatesTotal()              const {return fNV0candidates[0];}
  Int_t     NV0Candidates()                   const {return fNV0candidates[1];}
  Int_t     NDielectrons()                    const {return fNDielectronCandidates;}
  Int_t     NTracksTotal()                    const {return fNtracks[0];}
  Int_t     NTracks()                         const {return fNtracks[1];}
  Int_t     SPDntracklets()                   const {return fSPDntracklets;}
  
  Float_t   MultChannelVZERO(Int_t channel)   const {return (channel>=0 && channel<=63 ? fVZEROMult[channel] : -999.);}
  Float_t   MultVZEROA()                      const;
  Float_t   MultVZEROC()                      const;
  Float_t   MultVZERO()                       const;
  Float_t   MultRingVZEROA(Int_t ring)        const;
  Float_t   MultRingVZEROC(Int_t ring)        const;
  
  Float_t   EnergyZDC(Int_t channel)   const {return (channel>=0 && channel<8 ? fZDCnEnergy[channel] : -999.);}
  Float_t   EnergyZDCnA(Int_t channel) const {return (channel>=0 && channel<4 ? fZDCnEnergy[channel+4] : -999.);}
  Float_t   EnergyZDCnC(Int_t channel) const {return (channel>=0 && channel<4 ? fZDCnEnergy[channel] : -999.);}
  
  AliCorrelationReducedTrack* GetTrack(Int_t i)         const {return (i<fNtracks[1] ? (AliCorrelationReducedTrack*)fTracks->At(i) : 0x0);}
  AliCorrelationReducedPair* GetV0Pair(Int_t i)         const {return (i>=0 && i<fNV0candidates[1] ? (AliCorrelationReducedPair*)fCandidates->At(i) : 0x0);}
  AliCorrelationReducedPair* GetDielectronPair(Int_t i) const {return (i>=0 && i<fNDielectronCandidates ? (AliCorrelationReducedPair*)fCandidates->At(i+fNV0candidates[1]) : 0x0);}
  TClonesArray* GetPairs()                              const {return fCandidates;}
  TClonesArray* GetTracks()                             const {return fTracks;}

  Int_t GetNCaloClusters() const {return fNCaloClusters;}
  AliCorrelationReducedCaloCluster* GetCaloCluster(Int_t i) const {return (i>=0 && i<fNCaloClusters ? (AliCorrelationReducedCaloCluster*)fCaloClusters->At(i) : 0x0);}
  
  void GetQvector(Double_t Qvec[][2], Int_t det);
  void GetTPCQvector(Double_t Qvec[][2], Int_t det);
  void GetVZEROQvector(Double_t Qvec[][2], Int_t det);
  void GetVZEROQvector(Double_t Qvec[][2], Int_t det, Float_t* vzeroMult);
  void GetZDCQvector(Double_t Qvec[][2], Int_t det);

 private:
  Int_t     fRunNo;                 // run number
  UShort_t  fBC;                    // bunch crossing
  ULong64_t fTriggerMask;           // trigger mask
  Bool_t    fIsPhysicsSelection;    // PhysicsSelection passed event
  Float_t   fVtx[3];                // global event vertex vector in cm
  Int_t     fNVtxContributors;      // global event vertex contributors
  Float_t   fVtxTPC[3];             // TPC only event vertex           
  Int_t     fNVtxTPCContributors;   // TPC only event vertex contributors
  Float_t   fCentrality[4];         // centrality; 0-VZERO, 1-SPD, 2-TPC, 3-ZEMvsZDC 
  Int_t     fCentQuality;           // quality flag for the centrality 
  Int_t     fNV0candidates[2];      // number of V0 candidates, [0]-total, [1]-selected for the tree
  Int_t     fNDielectronCandidates; // number of pairs selected as dielectrons
  Int_t     fNtracks[2];            // number of tracks, [0]-total, [1]-selected for the tree
  Int_t     fSPDntracklets;         // number of SPD tracklets in |eta|<1.0 

  Float_t   fVZEROMult[64];         // VZERO multiplicity in all 64 channels
  Float_t   fZDCnEnergy[8];         // neutron ZDC energy in all 8 channels
    
  TClonesArray* fTracks;            //->   array containing global tracks
  static TClonesArray* fgTracks;
  
  TClonesArray* fCandidates;        //->   array containing pair candidates
  static TClonesArray* fgCandidates;
  
  Int_t     fNCaloClusters;         // number of calorimeter clusters  
  TClonesArray* fCaloClusters;        //->   array containing calorimeter clusters
  static TClonesArray* fgCaloClusters;
  
  void ClearEvent();
  AliCorrelationReducedEvent(const AliCorrelationReducedEvent &c);
  AliCorrelationReducedEvent& operator= (const AliCorrelationReducedEvent &c);

  ClassDef(AliCorrelationReducedEvent, 2);
};

//_______________________________________________________________________________
inline UShort_t AliCorrelationReducedTrack::ITSncls() const
{
  //
  // ITS number of clusters from the cluster map
  //
  UShort_t ncls=0;
  for(Int_t i=0; i<6; ++i) ncls += (ITSLayerHit(i) ? 1 : 0);
  return ncls;
}


//_______________________________________________________________________________
inline Int_t AliCorrelationReducedTrack::TPCClusterMapBitsFired()  const
{
  //
  // Count the number of bits fired in the TPC cluster map
  //
  Int_t nbits=0;
  for(Int_t i=0; i<8; ++i) nbits += (TPCClusterMapBitFired(i) ? 1 : 0);
  return nbits;
}


//_______________________________________________________________________________
inline Float_t AliCorrelationReducedPair::Energy() const 
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
    default:
      mass = fMass[0];
      break;    
  }
  Float_t p = P();
  return TMath::Sqrt(mass*mass+p*p);
}


//_______________________________________________________________________________
inline Float_t AliCorrelationReducedPair::Rapidity() const
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


//_______________________________________________________________________________
inline Double_t AliCorrelationReducedEventFriend::EventPlane(Int_t det, Int_t harmonic) const
{
  //
  // Event plane from detector "det" and harmonic "harmonic"
  //
  if(det<0 || det>=kNdetectors || harmonic<1 || harmonic>fgkNMaxHarmonics) return -999.;
  return TMath::ATan2(fQvector[det][harmonic-1][1], fQvector[det][harmonic-1][0])/Double_t(harmonic);
}

//_______________________________________________________________________________
inline Bool_t AliCorrelationReducedEventFriend::CheckEventPlaneStatus(Int_t det, Int_t h, EventPlaneStatus flag) const {
  //
  // Check the status of the event plane for a given detector and harmonic
  //
  if(det<0 || det>=kNdetectors || h<1 || h>fgkNMaxHarmonics) return kFALSE;
  return (flag<kNMaxFlowFlags ? (fEventPlaneStatus[det][h]&(1<<flag)) : kFALSE);
}

#endif
