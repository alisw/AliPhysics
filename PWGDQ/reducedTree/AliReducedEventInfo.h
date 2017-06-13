// Class for reduced event information
// Author: Ionut-Cristian Arsene (iarsene@cern.ch)
// 

#ifndef ALIREDUCEDEVENTINFO_H
#define ALIREDUCEDEVENTINFO_H

#include "AliReducedBaseEvent.h"
#include "AliReducedEventPlaneInfo.h"

class AliReducedCaloClusterInfo;
class AliReducedPairInfo;
class AliReducedTrackInfo;

//_________________________________________________________________________
class AliReducedEventInfo : public AliReducedBaseEvent {

  friend class AliAnalysisTaskReducedTreeMaker;     // friend analysis task which fills the object

 public:
  AliReducedEventInfo();
  AliReducedEventInfo(const Char_t* name, Int_t trackOption = AliReducedBaseEvent::kNoInit);
  virtual ~AliReducedEventInfo();

  // getters
  Int_t     EventNumberInFile()               const {return fEventNumberInFile;}
  UInt_t    L0TriggerInputs()                 const {return fL0TriggerInputs;}
  Bool_t    L0TriggerInput(UShort_t bit)      const {return (bit<8*sizeof(UInt_t) ? (fL0TriggerInputs&(UInt_t(1)<<bit)) : kFALSE);}
  UInt_t    L1TriggerInputs()                 const {return fL1TriggerInputs;}
  Bool_t    L1TriggerInput(UShort_t bit)      const {return (bit<8*sizeof(UInt_t) ? (fL1TriggerInputs&(UInt_t(1)<<bit)) : kFALSE);}
  UShort_t  L2TriggerInputs()                 const {return fL2TriggerInputs;}
  Bool_t    L2TriggerInput(UShort_t bit)      const {return (bit<8*sizeof(UShort_t) ? (fL2TriggerInputs&(UShort_t(1)<<bit)) : kFALSE);}
  UShort_t  BC()                              const {return fBC;}
  UInt_t    TimeStamp()                       const {return fTimeStamp;}
  UInt_t    EventType()                       const {return fEventType;}
  ULong64_t TriggerMask()                     const {return fTriggerMask;}
  Bool_t    IsPhysicsSelection()              const {return fIsPhysicsSelection;}
  Bool_t    IsSPDPileup()                     const {return fIsSPDPileup;}
  Bool_t    IsSPDPileupMultBins()             const {return fIsSPDPileupMultBins;}
  Int_t     IRIntClosestIntMap(Int_t id)      const {return (id>=0 && id<2 ? fIRIntClosestIntMap[id] : -999);}
  Float_t   VertexCovMatrix(Int_t iCov = 0)   const {return (iCov>=0 && iCov<6 ? fVtxCovMatrix[iCov] : 0.0);}
  Float_t   VertexTPC(Int_t axis)             const {return (axis>=0 && axis<=2 ? fVtxTPC[axis] : 0);}
  Int_t     VertexTPCContributors()           const {return fNVtxTPCContributors;}
  Float_t   VertexSPD(Int_t axis)             const {return (axis>=0 && axis<=2 ? fVtxSPD[axis] : 0);}
  Int_t     VertexSPDContributors()           const {return fNVtxSPDContributors;}
  Int_t     NTPCClusters()                    const {return fNTPCclusters;}
  Float_t   VertexTZERO()                     const {return fT0zVertex;}
  Int_t     NpileupSPD()                      const {return fNpileupSPD;}
  Int_t     NpileupTracks()                   const {return fNpileupTracks;}
  Int_t     NPMDtracks()                      const {return fNPMDtracks;}
  Int_t     NTRDtracks()                      const {return fNTRDtracks;}
  Int_t     NTRDtracklets()                   const {return fNTRDtracklets;}
  Int_t     SPDntracklets()                   const {return fSPDntracklets;}
  Int_t     SPDntracklets(Int_t bin)          const {return (bin>=0 && bin<32 ? fSPDntrackletsEta[bin] : -999);}
  Short_t   SPDFiredChips(Int_t layer)        const {return (layer==1 || layer==2 ? fSPDFiredChips[layer-1] : -999);}
  UInt_t    ITSClusters(Int_t layer)          const {return (layer>=1 && layer<=6 ? fITSClusters[layer-1] : 0);}
  Int_t     SPDnSingleClusters()              const {return fSPDnSingle;}
  Int_t     TracksPerTrackingFlag(Int_t flag) const {return (flag>=0 && flag<32 ? fNtracksPerTrackingFlag[flag] : -999);}
  
  Float_t   MultEstimatorOnlineV0M()   const {return fMultiplicityEstimators[0];}
  Float_t   MultEstimatorOnlineV0A()   const {return fMultiplicityEstimators[1];}
  Float_t   MultEstimatorOnlineV0C()   const {return fMultiplicityEstimators[2];}
  Float_t   MultEstimatorADM()   const {return fMultiplicityEstimators[3];}
  Float_t   MultEstimatorADA()   const {return fMultiplicityEstimators[4];}
  Float_t   MultEstimatorADC()   const {return fMultiplicityEstimators[5];}
  Float_t   MultEstimatorSPDClusters()   const {return fMultiplicityEstimators[6];}
  Float_t   MultEstimatorSPDTracklets()   const {return fMultiplicityEstimators[7];}
  Float_t   MultEstimatorRefMult05()   const {return fMultiplicityEstimators[8];}
  Float_t   MultEstimatorRefMult08()   const {return fMultiplicityEstimators[9];}
  Float_t   MultEstimatorPercentileOnlineV0M()   const {return fMultiplicityEstimatorPercentiles[0];}
  Float_t   MultEstimatorPercentileOnlineV0A()   const {return fMultiplicityEstimatorPercentiles[1];}
  Float_t   MultEstimatorPercentileOnlineV0C()   const {return fMultiplicityEstimatorPercentiles[2];}
  Float_t   MultEstimatorPercentileADM()   const {return fMultiplicityEstimatorPercentiles[3];}
  Float_t   MultEstimatorPercentileADA()   const {return fMultiplicityEstimatorPercentiles[4];}
  Float_t   MultEstimatorPercentileADC()   const {return fMultiplicityEstimatorPercentiles[5];}
  Float_t   MultEstimatorPercentileSPDClusters()   const {return fMultiplicityEstimatorPercentiles[6];}
  Float_t   MultEstimatorPercentileSPDTracklets()   const {return fMultiplicityEstimatorPercentiles[7];}
  Float_t   MultEstimatorPercentileRefMult05()   const {return fMultiplicityEstimatorPercentiles[8];}
  Float_t   MultEstimatorPercentileRefMult08()   const {return fMultiplicityEstimatorPercentiles[9];}
  
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
  TClonesArray* GetFMD() const {return fFMD;}
  
  Bool_t    TestEventTag(UShort_t iflag)  const {return (iflag<8*sizeof(ULong64_t) ? fEventTag&(ULong64_t(1)<<iflag) : kFALSE);}
  Bool_t    SetEventTag(UShort_t iflag)         {if (iflag>=8*sizeof(ULong64_t)) return kFALSE; fEventTag|=(ULong64_t(1)<<iflag); return kTRUE;}

  Bool_t    IsTriggerFired(UShort_t iflag)      {if (iflag>=8*sizeof(ULong64_t)) return kFALSE; return (fTriggerMask&(ULong64_t(1)<<iflag) ? kTRUE : kFALSE);}

  Double_t GetQvectorFMD(Int_t c, Double_t etamin, Double_t etamax);

  Int_t GetNCaloClusters() const {return fNCaloClusters;}
  AliReducedCaloClusterInfo* GetCaloCluster(Int_t i) const 
    {return (i>=0 && i<fNCaloClusters ? (AliReducedCaloClusterInfo*)fCaloClusters->At(i) : 0x0);}
  
  void  GetQvector(Double_t Qvec[][2], Int_t det, Float_t etaMin=-0.8, Float_t etaMax=+0.8, Bool_t (*IsTrackSelected)(AliReducedTrackInfo*)=NULL);
  Int_t GetTPCQvector(Double_t Qvec[][2], Int_t det, Float_t etaMin=-0.8, Float_t etaMax=+0.8, Bool_t (*IsTrackSelected)(AliReducedTrackInfo*)=NULL);
  void  GetVZEROQvector(Double_t Qvec[][2], Int_t det) ;
  void  GetVZEROQvector(Double_t Qvec[][2], Int_t det, Float_t* vzeroMult);
  void  GetZDCQvector(Double_t Qvec[][2], Int_t det) const;
  void  GetZDCQvector(Double_t Qvec[][2], Int_t det, const Float_t* zdcEnergy) const;
  void  SubtractParticleFromQvector(AliReducedTrackInfo* particle, Double_t Qvec[][2], Int_t det,
                                    Float_t etaMin=-0.8, Float_t etaMax=+0.8,
	 			    Bool_t (*IsTrackSelected)(AliReducedTrackInfo*)=NULL);
  
  // Event plane information handling for the case when event plane information is written directly in the trees
  //void SetEventPlane(const AliReducedEventPlaneInfo* ep) 
  //    {if(ep) {fEventPlane=new AliReducedEventPlaneInfo(); fEventPlane->CopyEvent(ep);}};
  //AliReducedEventPlaneInfo* GetEventPlane() const {return fEventPlane;};
  void SetEventPlane(const AliReducedEventPlaneInfo* ep) {if(ep) fEventPlane.CopyEvent(ep);}
  
  //Double_t GetEventPlane(Int_t detector, Int_t harmonic) const 
   //   {if(fEventPlane) return fEventPlane->EventPlane(detector, harmonic); return 0.0;};
  Double_t GetEventPlane(Int_t detector, Int_t harmonic) const 
      {return fEventPlane.EventPlane(detector, harmonic); return 0.0;};    
      

  virtual void ClearEvent();
  
  static const Float_t fgkZdcNalpha;
  
 protected:
  Int_t     fEventNumberInFile;     // Event number in ESD file
  UInt_t    fL0TriggerInputs;       // L0 trigger inputs
  UInt_t    fL1TriggerInputs;       // L1 trigger inputs
  UShort_t  fL2TriggerInputs;       // L2 trigger inputs
  UShort_t  fBC;                    // bunch crossing
  UInt_t    fTimeStamp;             // time stamp of the event                
  UInt_t    fEventType;             // event type                             
  ULong64_t fTriggerMask;           // trigger mask
  Float_t   fMultiplicityEstimators[10];   // multiplicity estimators: "OnlineV0M", "OnlineV0A", "OnlineV0C", "ADM", "ADA", "ADC", "SPDClusters", "SPDTracklets", "RefMult05", "RefMult08"
  Float_t   fMultiplicityEstimatorPercentiles[10];   // multiplicity estimators: "OnlineV0M", "OnlineV0A", "OnlineV0C", "ADM", "ADA", "ADC", "SPDClusters", "SPDTracklets", "RefMult05", "RefMult08"
  Bool_t    fIsPhysicsSelection;    // PhysicsSelection passed event
  Bool_t    fIsSPDPileup;           // identified as pileup event by SPD
  Bool_t    fIsSPDPileupMultBins;   // identified as pileup event by SPD in multiplicity bins
  Int_t     fIRIntClosestIntMap[2]; // out of bunch interactions, [0]-Int1, [1]-Int2 
  Float_t   fVtxCovMatrix[6];       // Covariance matrix of the event vertex
  Float_t   fVtxTPC[3];             // TPC only event vertex       
  Int_t     fNVtxTPCContributors;   // TPC only event vertex contributors
  Float_t   fVtxSPD[3];             // SPD only event vertex
  Int_t     fNVtxSPDContributors;  // SPD only event vertex contributors
  Int_t     fNpileupSPD;            // number of pileup vertices from SPD     
  Int_t     fNpileupTracks;         // number of pileup vertices from tracks  
  Int_t     fNTPCclusters;          // number of TPC clusters
  Int_t     fNPMDtracks;            // number of PMD tracks                   
  Int_t     fNTRDtracks;            // number of TRD tracks                   
  Int_t     fNTRDtracklets;         // number of TRD tracklets                
  Int_t     fSPDntracklets;         // number of SPD tracklets in |eta|<1.0 
  Int_t     fSPDntrackletsEta[32];  // number of SPD tracklets in equal eta bins between -1.6 --> +1.6    
  Short_t   fSPDFiredChips[2];      // number of fired chips in the two layers
  UInt_t    fITSClusters[6];        // number of ITS clusters per layer
  Int_t     fSPDnSingle;            // number of clusters in SPD layer 1, not associated to a tracklet on SPD layer 2
  Int_t     fNtracksPerTrackingFlag[32];  // number of tracks for each tracking status bit                
  Float_t   fVZEROMult[64];         // VZERO multiplicity in all 64 channels
  Float_t   fVZEROTotalMult[2];    // Total VZERO multiplicity
  Float_t   fZDCnEnergy[10];         // neutron ZDC energy in all 8 channels
  Float_t   fZDCpEnergy[10];         // proton ZDC energy in all 8 channels
  Float_t   fZDCnTotalEnergy[2];   // total neutron ZDC energy
  Float_t   fZDCpTotalEnergy[2];  // total proton ZDC energy
  Float_t   fT0amplitude[26];         // T0 amplitude in all 24 channels
  Float_t   fT0TOF[3];               // T0 timing for A&C, A, and C (first time)
  Float_t   fT0TOFbest[3];           // T0 timing for A&C, A, and C (best time)
  Float_t   fT0zVertex;                // T0 z vertex estimation
  Float_t   fT0start;                // T0 timing
  Bool_t    fT0pileup;               // TZERO pileup flag
  Bool_t    fT0sattelite;            // TZERO flag for collisions from sattelite bunches
    
  Int_t     fNCaloClusters;         // number of calorimeter clusters  
  TClonesArray* fCaloClusters;        //->   array containing calorimeter clusters
  static TClonesArray* fgCaloClusters;     // calorimeter clusters
  TClonesArray* fFMD;            //->   array containing fmd readout          (NEW)
  static TClonesArray* fgFMD;    //       fmd readout			      (NEW)

  //AliReducedEventPlaneInfo* fEventPlane;     //-> container for event plane information
  AliReducedEventPlaneInfo fEventPlane;     // container for event plane information
  
  AliReducedEventInfo(const AliReducedEventInfo &c);
  AliReducedEventInfo& operator= (const AliReducedEventInfo &c);

  ClassDef(AliReducedEventInfo, 6);
};

#endif
