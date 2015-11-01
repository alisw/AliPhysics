// Classes used for creating a reduced information tree
// Author: Ionut-Cristian Arsene (i.c.arsene@gsi.de)
// 
//  Basic structure:
//  1. Event wise information
//  2. List of tracks in the event
//  3. List of resonance candidates

#ifndef ALIREDUCEDEVENT_H
#define ALIREDUCEDEVENT_H

#include <TClonesArray.h>
#include <TBits.h>
#include <TMath.h>


const Int_t fgkNMaxHarmonics = 10;
const Float_t gZdcNalpha= 1.0;

//_____________________________________________________________________
class AliReducedTrack : public TObject {

  friend class AliAnalysisTaskReducedTree;  // friend analysis task which fills the object

 public:
  AliReducedTrack();
  ~AliReducedTrack();

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
  Float_t  TrackLength()                 const {return fTrackLength;}
  
  UShort_t ITSncls()                const;
  UChar_t  ITSclusterMap()          const {return fITSclusterMap;}
  Bool_t   ITSLayerHit(Int_t layer) const {return (layer>=0 && layer<6 ? (fITSclusterMap&(1<<layer)) : kFALSE);};
  Float_t  ITSsignal()              const {return fITSsignal;}
  Float_t  ITSnSig(Int_t specie)    const {return (specie>=0 && specie<=3 ? fITSnSig[specie] : -999.);}
  Float_t  ITSchi2()                const {return fITSchi2;}
  
  UChar_t TPCncls()                        const {return fTPCNcls;}
  UChar_t TPCFindableNcls()                const {return fTPCNclsF;}
  UChar_t TPCCrossedRows()                 const {return fTPCCrossedRows;}
  UChar_t TPCnclsIter1()                   const {return fTPCNclsIter1;}
  UChar_t TPCClusterMap()                  const {return fTPCClusterMap;}
  Int_t   TPCClusterMapBitsFired()         const;
  Bool_t  TPCClusterMapBitFired(Int_t bit) const {return (bit>=0 && bit<8 ? (fTPCClusterMap&(1<<bit)) : kFALSE);};
  Float_t TPCsignal()                      const {return fTPCsignal;}
  UChar_t TPCsignalN()                     const {return fTPCsignalN;}
  Float_t TPCnSig(Int_t specie)            const {return (specie>=0 && specie<=3 ? fTPCnSig[specie] : -999.);}
  Float_t TPCchi2()                        const {return fTPCchi2;}
  
  Float_t  TOFbeta()             const {return fTOFbeta;}    
  Float_t  TOFnSig(Int_t specie) const {return (specie>=0 && specie<=3 ? fTOFnSig[specie] : -999.);}
  Short_t  TOFdeltaBC()          const {return fTOFdeltaBC;}
  
  Int_t    TRDntracklets(Int_t type)  const {return (type==0 || type==1 ? fTRDntracklets[type] : -1);}
  Float_t  TRDpid(Int_t specie)       const {return (specie>=0 && specie<=1 ? fTRDpid[specie] : -999.);}
  Float_t  TRDpidLQ1D(Int_t specie)   const {return (specie>=0 && specie<=1 ? fTRDpid[specie] : -999.);}
  Float_t  TRDpidLQ2D(Int_t specie)   const {return (specie>=0 && specie<=1 ? fTRDpidLQ2D[specie] : -999.);}
  
  Int_t    CaloClusterId() const {return fCaloClusterId;}
    
  Float_t  BayesPID(Int_t specie) const {return (specie>=0 && specie<=2 ? fBayesPID[specie] : -999.);}
  
  Bool_t UsedForQvector()             const {return fFlags&(UShort_t(1)<<0);}
  Bool_t TestFlag(UShort_t iflag)     const {return (iflag<8*sizeof(UShort_t) ? fFlags&(UShort_t(1)<<iflag) : kFALSE);}
  Bool_t SetFlag(UShort_t iflag)            {if (iflag>=8*sizeof(UShort_t)) return kFALSE; fFlags|=(UShort_t(1)<<iflag); return kTRUE;}
  Bool_t IsGammaLeg()                 const {return fFlags&(UShort_t(1)<<1);}
  Bool_t IsPureGammaLeg()             const {return fFlags&(UShort_t(1)<<8);}
  Bool_t IsK0sLeg()                   const {return fFlags&(UShort_t(1)<<2);}
  Bool_t IsPureK0sLeg()               const {return fFlags&(UShort_t(1)<<9);}
  Bool_t IsLambdaLeg()                const {return fFlags&(UShort_t(1)<<3);}
  Bool_t IsPureLambdaLeg()            const {return fFlags&(UShort_t(1)<<10);}
  Bool_t IsALambdaLeg()               const {return fFlags&(UShort_t(1)<<4);}
  Bool_t IsPureALambdaLeg()           const {return fFlags&(UShort_t(1)<<11);}
  Bool_t IsKink(Int_t i=0)            const {return (i>=0 && i<3 ? fFlags&(UShort_t(1)<<(5+i)) : kFALSE);}
  Bool_t TestFlagMore(UShort_t iflag) const {return (iflag<8*sizeof(ULong_t) ? fMoreFlags&(ULong_t(1)<<iflag) : kFALSE);}
  Bool_t SetFlagMore(UShort_t iflag)        {if(iflag>=8*sizeof(ULong_t)) return kFALSE; fMoreFlags|=(ULong_t(1)<<iflag); return kTRUE;}
  Bool_t UnsetFlagMore(UShort_t iflag)      {if(iflag>=8*sizeof(ULong_t)) return kFALSE; fMoreFlags^=(ULong_t(1)<<iflag); return kTRUE;}
  ULong_t GetFlagsMore()               const {return fMoreFlags;}
  
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
  Float_t fTrackLength;         // track length
  
  // ITS
  UChar_t  fITSclusterMap;      // ITS cluster map
  Float_t  fITSsignal;          // ITS signal
  Float_t  fITSnSig[4];         // 0-electron; 1-pion; 2-kaon; 3-proton
  Float_t  fITSchi2;            // ITS chi2 / cls
  
  // TPC
  UChar_t fTPCNcls;            // TPC ncls                          
  UChar_t fTPCCrossedRows;     // TPC crossed rows                  
  UChar_t fTPCNclsF;           // TPC findable ncls                 
  UChar_t fTPCNclsIter1;       // TPC no clusters after first iteration
  UChar_t fTPCClusterMap;      // TPC cluster distribution map
  Float_t fTPCsignal;          // TPC de/dx
  UChar_t fTPCsignalN;         // TPC no clusters de/dx
  Float_t fTPCnSig[4];         // 0-electron; 1-pion; 2-kaon; 3-proton
  Float_t fTPCchi2;            // TPC chi2 / cls
    
  // TOF
  Float_t fTOFbeta;             // TOF pid info
  Float_t fTOFnSig[4];          // TOF n-sigma deviation from expected signal
  Short_t fTOFdeltaBC;          // BC(event) - BC(track) estimated by TOF
  
  // TRD
  UChar_t fTRDntracklets[2];       // 0 - AliESDtrack::GetTRDntracklets(); 1 - AliESDtrack::GetTRDntrackletsPID()   TODO: use only 1 char
  Float_t fTRDpid[2];              // TRD pid 1D likelihoods, [0]-electron , [1]- pion
  Float_t fTRDpidLQ2D[2];          // TRD pid 2D likelihoods, [0]-electron , [1]- pion
  
  // EMCAL/PHOS
  Int_t  fCaloClusterId;          // ID for the calorimeter cluster (if any)
  
  // Bayesian PID
  Float_t fBayesPID[3];           // Combined Bayesian PID   pi/K/p
  
  UShort_t fFlags;                // BIT0 toggled if track used for TPC event plane
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
  ULong_t  fMoreFlags;            // Space reserved for more information which might be needed later for analysis
    
  AliReducedTrack(const AliReducedTrack &c);
  AliReducedTrack& operator= (const AliReducedTrack &c);

  ClassDef(AliReducedTrack, 4);
};


//_____________________________________________________________________
class AliReducedPair : public TObject {

  friend class AliAnalysisTaskReducedTree;  // friend analysis task which fills the object

 public:
  enum CandidateType {
    kK0sToPiPi=0,
    kPhiToKK,
    kLambda0ToPPi,
    kALambda0ToPPi,
    kJpsiToEE,
    kUpsilon,
    kGammaConv,
    kNMaxCandidateTypes
  };
  AliReducedPair();
  AliReducedPair(const AliReducedPair &c);
  ~AliReducedPair();

  // getters
  Char_t   CandidateId()         const {return fCandidateId;}
  Char_t   PairType()            const {return fPairType;}
  Int_t    LegId(Int_t leg)      const {return (leg==0 || leg==1 ? fLegIds[leg] : -1);}
  Float_t  Mass(Int_t idx=0)     const {return (idx>=0 && idx<4 ? fMass[idx] : -999.);}
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
  Float_t  DecayRadius()         const {return fLxy;}
  Float_t  LxyErr()              const {return fLxyErr;}
  Float_t  PointingAngle()       const {return fPointingAngle;}
  Float_t  Chi2()                const {return fChisquare;}
  Bool_t   IsOnTheFly()          const {return fPairType;}
  Bool_t   IsPureV0K0s()         const {return (fMCid&(UInt_t(1)<<1));}
  Bool_t   IsPureV0Lambda()      const {return (fMCid&(UInt_t(1)<<2));}
  Bool_t   IsPureV0ALambda()     const {return (fMCid&(UInt_t(1)<<3));}
  Bool_t   IsPureV0Gamma()       const {return (fMCid&(UInt_t(1)<<4));}
  UInt_t   MCid()                const {return fMCid;}
  Bool_t   CheckMC(const Int_t flag) const {return (flag<32 ? (fMCid&(1<<flag)) : kFALSE);}
  
 private:
  Char_t  fCandidateId;         // candidate type (K0s, Lambda, J/psi, phi, etc)
  Char_t  fPairType;            // 0 ++; 1 +-; 2 -- for dielectron pairs; 0- offline, 1- on the fly for V0 candidates
  UShort_t fLegIds[2];          // leg ids 
  Float_t fMass[4];             // invariant mass for pairs (3 extra mass values for other V0 pid assumptions)
                                // idx=0 -> K0s assumption; idx=1 -> Lambda; idx=2 -> anti-Lambda; idx=3 -> gamma conversion
  Float_t fPhi;                 // pair phi in the [0,2*pi) interval
  Float_t fPt;                  // pair pt
  Float_t fEta;                 // pair eta 
  Float_t fLxy;                 // pseudo-proper decay length (pair candidates) or radius of the secondary vertex for V0s 
  Float_t fLxyErr;              // error on Lxy
  Float_t fPointingAngle;       // angle between the pair momentum vector and the secondary vertex position vector
  Float_t fChisquare;           // chi2 for the legs matching
  UInt_t  fMCid;                // Bit map with Monte Carlo info about the pair

  AliReducedPair& operator= (const AliReducedPair &c);

  ClassDef(AliReducedPair, 3);
};


//_____________________________________________________________________
class AliReducedFMD : public TObject {

  friend class AliAnalysisTaskReducedTree;  // friend analysis task which fills the object

 public:
  AliReducedFMD();
  ~AliReducedFMD();

  // getters
  UShort_t Multiplicity()	    const {return fMultiplicity;}
  //Float_t Eta()	 		          const {return fEta;}
  UShort_t Id()	    	        const {return fId;}

  //void SetIgnoreSteamer();
  
 private:

  UShort_t fMultiplicity;        
  //Float_t fEta;            
  UShort_t fId;           

    
  AliReducedFMD(const AliReducedFMD &c);
  AliReducedFMD& operator= (const AliReducedFMD &c);

  ClassDef(AliReducedFMD, 1);
};


//_________________________________________________________________________
class AliReducedEventFriend : public TObject {
  
  friend class AliAnalysisTaskReducedTree;    // friend analysis task which fills the object
  
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
  
  AliReducedEventFriend();
  ~AliReducedEventFriend();
  
  Double_t Qx(Int_t det, Int_t harmonic)  const {return (det>=0 && det<kNdetectors && harmonic>0 && harmonic<=fgkNMaxHarmonics ? fQvector[det][harmonic-1][0] : -999.);}
  Double_t Qy(Int_t det, Int_t harmonic)  const {return (det>=0 && det<kNdetectors && harmonic>0 && harmonic<=fgkNMaxHarmonics ? fQvector[det][harmonic-1][1] : -999.);}
  Double_t EventPlane(Int_t det, Int_t h) const;
  UChar_t GetEventPlaneStatus(Int_t det, Int_t h) const {return (det>=0 && det<kNdetectors && h>0 && h<=fgkNMaxHarmonics ? fEventPlaneStatus[det][h-1] : UChar_t(255));} 
  Bool_t  CheckEventPlaneStatus(Int_t det, Int_t h, EventPlaneStatus flag) const;
  void    CopyEvent(const AliReducedEventFriend* event);

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
  AliReducedEventFriend(const AliReducedEventFriend &c);
  AliReducedEventFriend& operator= (const AliReducedEventFriend &c);

  ClassDef(AliReducedEventFriend, 1);
};


//_________________________________________________________________________
class AliReducedCaloCluster : public TObject {
  
  friend class AliAnalysisTaskReducedTree;         // friend analysis task which fills the object
  
 public:
  enum ClusterType {
    kUndefined=0, kEMCAL, kPHOS  
  };
   
  AliReducedCaloCluster();
  ~AliReducedCaloCluster();
  
  Bool_t  IsEMCAL()    const {return (fType==kEMCAL ? kTRUE : kFALSE);}
  Bool_t  IsPHOS()     const {return (fType==kPHOS ? kTRUE : kFALSE);}
  Float_t Energy()     const {return fEnergy;}
  Float_t Dx()         const {return fTrackDx;}
  Float_t Dz()         const {return fTrackDz;}
  Float_t M20()        const {return fM20;}
  Float_t M02()        const {return fM02;}
  Float_t Dispersion() const {return fDispersion;}
  
 private:
  Char_t  fType;         // cluster type (EMCAL/PHOS)
  Float_t fEnergy;       // cluster energy
  Float_t fTrackDx;      // distance to closest track in phi
  Float_t fTrackDz;      // distance to closest track in z
  Float_t fM20;          // short axis
  Float_t fM02;          // long axis
  Float_t fDispersion;   // dispersion
  
  AliReducedCaloCluster(const AliReducedCaloCluster &c);
  AliReducedCaloCluster& operator= (const AliReducedCaloCluster &c);

  ClassDef(AliReducedCaloCluster, 2);
};


//_________________________________________________________________________
class AliReducedEvent : public TObject {

  friend class AliAnalysisTaskReducedTree;     // friend analysis task which fills the object

 public:
  AliReducedEvent();
  AliReducedEvent(const Char_t* name);
  ~AliReducedEvent();

  // getters
  ULong64_t EventTag()                        const {return fEventTag;}
  Bool_t    EventTag(UShort_t bit)            const {return (bit<8*sizeof(ULong64_t) ? (fEventTag&(ULong64_t(1)<<bit)) : kFALSE);}
  Int_t     EventNumberInFile()               const {return fEventNumberInFile;}
  UInt_t    L0TriggerInputs()                 const {return fL0TriggerInputs;}
  Bool_t    L0TriggerInput(UShort_t bit)         const {return (bit<8*sizeof(UInt_t) ? (fL0TriggerInputs&(UInt_t(1)<<bit)) : kFALSE);}
  UInt_t    L1TriggerInputs()                 const {return fL1TriggerInputs;}
  Bool_t    L1TriggerInput(UShort_t bit)         const {return (bit<8*sizeof(UInt_t) ? (fL1TriggerInputs&(UInt_t(1)<<bit)) : kFALSE);}
  UShort_t  L2TriggerInputs()                 const {return fL2TriggerInputs;}
  Bool_t    L2TriggerInput(UShort_t bit)         const {return (bit<8*sizeof(UShort_t) ? (fL2TriggerInputs&(UShort_t(1)<<bit)) : kFALSE);}
  Int_t     RunNo()                           const {return fRunNo;}
  UShort_t  BC()                              const {return fBC;}
  UInt_t    TimeStamp()                       const {return fTimeStamp;}
  UInt_t    EventType()                       const {return fEventType;}
  ULong64_t TriggerMask()                     const {return fTriggerMask;}
  Bool_t    IsPhysicsSelection()              const {return fIsPhysicsSelection;}
  Bool_t    IsSPDPileup()                     const {return fIsSPDPileup;}
  Bool_t    IsSPDPileupMultBins()             const {return fIsSPDPileupMultBins;}
  Int_t     IRIntClosestIntMap(Int_t id)      const {return (id>=0 && id<2 ? fIRIntClosestIntMap[id] : -999);}
  Bool_t    IsFMDReduced()                    const {return fIsFMDReduced;}
  Float_t   Vertex(Int_t axis)                const {return (axis>=0 && axis<=2 ? fVtx[axis] : 0);}
  Int_t     VertexNContributors()             const {return fNVtxContributors;}
  Float_t   VertexTPC(Int_t axis)             const {return (axis>=0 && axis<=2 ? fVtxTPC[axis] : 0);}
  Int_t     VertexTPCContributors()           const {return fNVtxTPCContributors;}
  Float_t   VertexTZERO()                     const {return fT0zVertex;}
  Int_t     NpileupSPD()                      const {return fNpileupSPD;}
  Int_t     NpileupTracks()                   const {return fNpileupTracks;}
  Int_t     NPMDtracks()                      const {return fNPMDtracks;}
  Int_t     NTRDtracks()                      const {return fNTRDtracks;}
  Int_t     NTRDtracklets()                   const {return fNTRDtracklets;}
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
  Int_t     SPDntracklets(Int_t bin)          const {return (bin>=0 && bin<32 ? fSPDntrackletsEta[bin] : -999);}
  Int_t     TracksPerTrackingFlag(Int_t flag) const {return (flag>=0 && flag<32 ? fNtracksPerTrackingFlag[flag] : -999);}
  UShort_t  NFMDchannels(Int_t det)	      const {return (det>=0 && det<5 ? fNFMDchannels[det] : 0);}
  UShort_t  FMDtotalMult(Int_t det)           const {return (det>=0 && det<5 ? fFMDtotalMult[det] : 0);}
  
  Float_t   MultChannelVZERO(Int_t channel)   const {return (channel>=0 && channel<=63 ? fVZEROMult[channel] : -999.);}
  Float_t   MultVZEROA()                      const;
  Float_t   MultVZEROC()                      const;
  Float_t   MultVZERO()                       const;
  Float_t   MultRingVZEROA(Int_t ring)        const;
  Float_t   MultRingVZEROC(Int_t ring)        const;
  
  Float_t   AmplitudeTZEROA()                       const;
  Float_t   AmplitudeTZEROC()                       const;
  Float_t   AmplitudeTZERO()                        const;
  Float_t   AmplitudeTZEROch(Int_t ch)              const {return (ch>=0 && ch<=25 ? fT0amplitude[ch] : -999.);}
  Float_t   EventTZEROStartTime()                   const {return fT0start;}
  Float_t   EventTZEROStartTimeTOFfirst(Int_t side) const {return (side>=0 && side<3 ? fT0TOF[side] : -999.);}
  Float_t   EventTZEROStartTimeTOFbest(Int_t side)  const {return (side>=0 && side<3 ? fT0TOFbest[side] : -999.);}
  Bool_t    IsPileupTZERO()                         const {return fT0pileup;}
  Bool_t    IsSatteliteCollisionTZERO()             const {return fT0sattelite;}
  
  Float_t   EnergyZDCnTree(UShort_t channel)  const {return (channel<10 ? fZDCnEnergy[channel] : -999.);};
  Float_t   EnergyZDCpTree(UShort_t channel)  const {return (channel<10 ? fZDCpEnergy[channel] : -999.);};
  Float_t   EnergyZDCn(Int_t channel)  const;
  Float_t   EnergyZDCA() const;
  Float_t   EnergyZDCC() const;
  Float_t   EnergyZDC() const;
  
  Bool_t    TestEventTag(UShort_t iflag)  const {return (iflag<8*sizeof(ULong64_t) ? fEventTag&(ULong64_t(1)<<iflag) : kFALSE);}
  Bool_t    SetEventTag(UShort_t iflag)         {if (iflag>=8*sizeof(ULong64_t)) return kFALSE; fEventTag|=(ULong64_t(1)<<iflag); return kTRUE;}
  
  AliReducedTrack* GetTrack(Int_t i)         const 
    {return (i<fNtracks[1] ? (AliReducedTrack*)fTracks->At(i) : 0x0);}
  AliReducedPair* GetV0Pair(Int_t i)         const 
    {return (i>=0 && i<fNV0candidates[1] ? (AliReducedPair*)fCandidates->At(i) : 0x0);}
  AliReducedPair* GetDielectronPair(Int_t i) const 
    {return (i>=0 && i<fNDielectronCandidates ? (AliReducedPair*)fCandidates->At(i+fNV0candidates[1]) : 0x0);}
  TClonesArray* GetPairs()                              const {return fCandidates;}
  TClonesArray* GetTracks()                             const {return fTracks;}
  TClonesArray* GetFMD(Int_t det)                       const {return (det==0 ? fFMD1 :
                                                                 (det==1 ? fFMD2I :
                                                                  (det==2 ? fFMD2O :
                                                                   (det==3 ? fFMD3I :
                                                                    (det==4 ? fFMD3O :
                                                                     0)))));}
  TClonesArray* GetFMD1()                                const {return fFMD1;}
  TClonesArray* GetFMD2I()                                const {return fFMD2I;}
  TClonesArray* GetFMD2O()                                const {return fFMD2O;}
  TClonesArray* GetFMD3I()                                const {return fFMD3I;}
  TClonesArray* GetFMD3O()                                const {return fFMD3O;}
  AliReducedFMD* GetFMD1Channel(UShort_t ch) const {return (ch<fNFMDchannels[0] ? (AliReducedFMD*)fFMD1->At(ch) : 0x0);}
  AliReducedFMD* GetFMD2IChannel(UShort_t ch) const {return (ch<fNFMDchannels[1] ? (AliReducedFMD*)fFMD2I->At(ch) : 0x0);}
  AliReducedFMD* GetFMD2OChannel(UShort_t ch) const {return (ch<fNFMDchannels[2] ? (AliReducedFMD*)fFMD2O->At(ch) : 0x0);}
  AliReducedFMD* GetFMD3IChannel(UShort_t ch) const {return (ch<fNFMDchannels[3] ? (AliReducedFMD*)fFMD3I->At(ch) : 0x0);}
  AliReducedFMD* GetFMD3OChannel(UShort_t ch) const {return (ch<fNFMDchannels[4] ? (AliReducedFMD*)fFMD3O->At(ch) : 0x0);}
  
  Int_t GetNCaloClusters() const {return fNCaloClusters;}
  AliReducedCaloCluster* GetCaloCluster(Int_t i) const 
    {return (i>=0 && i<fNCaloClusters ? (AliReducedCaloCluster*)fCaloClusters->At(i) : 0x0);}
  
  void  GetQvector(Double_t Qvec[][2], Int_t det, Float_t etaMin=-0.8, Float_t etaMax=+0.8, Bool_t (*IsTrackSelected)(AliReducedTrack*)=NULL);
  Int_t GetTPCQvector(Double_t Qvec[][2], Int_t det, Float_t etaMin=-0.8, Float_t etaMax=+0.8, Bool_t (*IsTrackSelected)(AliReducedTrack*)=NULL);
  void  GetVZEROQvector(Double_t Qvec[][2], Int_t det) ;
  void  GetVZEROQvector(Double_t Qvec[][2], Int_t det, Float_t* vzeroMult);
  void  GetZDCQvector(Double_t Qvec[][2], Int_t det) const;
  void  GetZDCQvector(Double_t Qvec[][2], Int_t det, const Float_t* zdcEnergy) const;
  void  SubtractParticleFromQvector(AliReducedTrack* particle, Double_t Qvec[][2], Int_t det,
                                    Float_t etaMin=-0.8, Float_t etaMax=+0.8,
	 			    Bool_t (*IsTrackSelected)(AliReducedTrack*)=NULL);

 private:
  ULong64_t fEventTag;              // Event tags to be used either during analysis or to filter events
  Int_t     fEventNumberInFile;     // Event number in ESD file
  UInt_t    fL0TriggerInputs;       // L0 trigger inputs
  UInt_t    fL1TriggerInputs;       // L1 trigger inputs
  UShort_t  fL2TriggerInputs;       // L2 trigger inputs
  Int_t     fRunNo;                 // run number
  UShort_t  fBC;                    // bunch crossing
  UInt_t    fTimeStamp;             // time stamp of the event                
  UInt_t    fEventType;             // event type                             
  ULong64_t fTriggerMask;           // trigger mask
  Bool_t    fIsPhysicsSelection;    // PhysicsSelection passed event
  Bool_t    fIsSPDPileup;           // identified as pileup event by SPD
  Bool_t    fIsSPDPileupMultBins;   // identified as pileup event by SPD in multiplicity bins
  Int_t     fIRIntClosestIntMap[2]; // out of bunch interactions, [0]-Int1, [1]-Int2 
  Bool_t    fIsFMDReduced;          // FMD info, if present, is reduced       (NEW)
  Float_t   fVtx[3];                // global event vertex vector in cm
  Int_t     fNVtxContributors;      // global event vertex contributors
  Float_t   fVtxTPC[3];             // TPC only event vertex           
  Int_t     fNVtxTPCContributors;   // TPC only event vertex contributors
  Int_t     fNpileupSPD;            // number of pileup vertices from SPD     
  Int_t     fNpileupTracks;         // number of pileup vertices from tracks  
  Int_t     fNPMDtracks;            // number of PMD tracks                   
  Int_t     fNTRDtracks;            // number of TRD tracks                   
  Int_t     fNTRDtracklets;         // number of TRD tracklets                
  Float_t   fCentrality[4];         // centrality; 0-VZERO, 1-SPD, 2-TPC, 3-ZEMvsZDC 
  Int_t     fCentQuality;           // quality flag for the centrality 
  Int_t     fNV0candidates[2];      // number of V0 candidates, [0]-total, [1]-selected for the tree
  Int_t     fNDielectronCandidates; // number of pairs selected as dielectrons
  Int_t     fNtracks[2];            // number of tracks, [0]-total, [1]-selected for the tree
  Int_t     fSPDntracklets;         // number of SPD tracklets in |eta|<1.0 
  Int_t     fSPDntrackletsEta[32];  // number of SPD tracklets in equal eta bins between -1.6 --> +1.6    
  Int_t     fNtracksPerTrackingFlag[32];  // number of tracks for each tracking status bit                
  
  Float_t   fVZEROMult[64];         // VZERO multiplicity in all 64 channels
  Float_t   fZDCnEnergy[10];         // neutron ZDC energy in all 8 channels
  Float_t   fZDCpEnergy[10];         // neutron ZDC energy in all 8 channels
  Float_t   fT0amplitude[26];         // T0 amplitude in all 24 channels
  Float_t   fT0TOF[3];               // T0 timing for A&C, A, and C (first time)
  Float_t   fT0TOFbest[3];           // T0 timing for A&C, A, and C (best time)
  Float_t   fT0zVertex;                // T0 z vertex estimation
  Float_t   fT0start;                // T0 timing
  Bool_t    fT0pileup;               // TZERO pileup flag
  Bool_t    fT0sattelite;            // TZERO flag for collisions from sattelite bunches
    
  TClonesArray* fTracks;            //->   array containing global tracks
  static TClonesArray* fgTracks;    //       global tracks
  
  TClonesArray* fCandidates;        //->   array containing pair candidates
  static TClonesArray* fgCandidates;  // pair candidates
  
  UShort_t   fNFMDchannels[5];           // number of FMD channels read out 	      (NEW)
  UShort_t   fFMDtotalMult[5];           // number of FMD channels read out 	      (NEW)
  TClonesArray* fFMD1;            //->   array containing fmd readout          (NEW)
  static TClonesArray* fgFMD1;    //       fmd readout			      (NEW)
  TClonesArray* fFMD2I;            //->   array containing fmd readout          (NEW)
  static TClonesArray* fgFMD2I;    //       fmd readout			      (NEW)
  TClonesArray* fFMD2O;            //->   array containing fmd readout          (NEW)
  static TClonesArray* fgFMD2O;    //       fmd readout			      (NEW)
  TClonesArray* fFMD3I;            //->   array containing fmd readout          (NEW)
  static TClonesArray* fgFMD3I;    //       fmd readout			      (NEW)
  TClonesArray* fFMD3O;            //->   array containing fmd readout          (NEW)
  static TClonesArray* fgFMD3O;    //       fmd readout			      (NEW)
  
  Int_t     fNCaloClusters;         // number of calorimeter clusters  
  TClonesArray* fCaloClusters;        //->   array containing calorimeter clusters
  static TClonesArray* fgCaloClusters;     // calorimeter clusters
  
  void ClearEvent();
  AliReducedEvent(const AliReducedEvent &c);
  AliReducedEvent& operator= (const AliReducedEvent &c);

  ClassDef(AliReducedEvent, 5);
};

//_______________________________________________________________________________
inline UShort_t AliReducedTrack::ITSncls() const
{
  //
  // ITS number of clusters from the cluster map
  //
  UShort_t ncls=0;
  for(Int_t i=0; i<6; ++i) ncls += (ITSLayerHit(i) ? 1 : 0);
  return ncls;
}


//_______________________________________________________________________________
inline Int_t AliReducedTrack::TPCClusterMapBitsFired()  const
{
  //
  // Count the number of bits fired in the TPC cluster map
  //
  Int_t nbits=0;
  for(Int_t i=0; i<8; ++i) nbits += (TPCClusterMapBitFired(i) ? 1 : 0);
  return nbits;
}


//_______________________________________________________________________________
inline Float_t AliReducedPair::Energy() const 
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
    case kGammaConv:
      mass = fMass[3];
      break;
    default:
      mass = fMass[0];
      break;    
  }
  Float_t p = P();
  return TMath::Sqrt(mass*mass+p*p);
}


//_______________________________________________________________________________
inline Float_t AliReducedPair::Rapidity() const
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
inline Double_t AliReducedEventFriend::EventPlane(Int_t det, Int_t harmonic) const
{
  //
  // Event plane from detector "det" and harmonic "harmonic"
  //
  if(det<0 || det>=kNdetectors || harmonic<1 || harmonic>fgkNMaxHarmonics) return -999.;
  return TMath::ATan2(fQvector[det][harmonic-1][1], fQvector[det][harmonic-1][0])/Double_t(harmonic);
}

//_______________________________________________________________________________
inline Bool_t AliReducedEventFriend::CheckEventPlaneStatus(Int_t det, Int_t h, EventPlaneStatus flag) const {
  //
  // Check the status of the event plane for a given detector and harmonic
  //
  if(det<0 || det>=kNdetectors || h<1 || h>fgkNMaxHarmonics) return kFALSE;
  return (flag<kNMaxFlowFlags ? (fEventPlaneStatus[det][h-1]&(1<<flag)) : kFALSE);
}

#endif
