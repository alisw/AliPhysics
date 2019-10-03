#ifndef ALITRACKLETALG_H
#define ALITRACKLETALG_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//_________________________________________________________________________
// 
//        Implementation of the ITS-SPD trackleter class
//   Clone version of the AliITSMultReconstructor class (October 2010) 
//   that can be used in an AliAnalysisTask
//
//_________________________________________________________________________
#include "AliTrackleter.h"

class TBits;
class TTree;
class TH1F;
class TH2F; 
class AliITSDetTypeRec;
class AliITSgeom;
class AliESDEvent;
class AliESDtrack;
class AliVertex;
class AliESDVertex;
class AliMultiplicity;

class AliTrackletAlg : public AliTrackleter
{
public:
  //
  enum {kClTh,kClPh,kClZ,kClMC0,kClMC1,kClMC2,kClNPar};
  enum {kTrTheta,kTrPhi,kTrDPhi,kTrDTheta,kTrLab1,kTrLab2,kClID1,kClID2,kTrNPar};
  enum {kSCTh,kSCPh,kSCLab,kSCID,kSCNPar};   
  enum {kITSTPC,kITSSAP,kITSTPCBit=BIT(kITSTPC),kITSSAPBit=BIT(kITSSAP)}; // RS
  AliTrackletAlg();
  virtual ~AliTrackletAlg();

  void Reconstruct(AliESDEvent* esd, TTree* treeRP);
  void Reconstruct(TTree* tree, Float_t* vtx, Float_t* vtxRes);   // old reconstructor invocation
  void FindTracklets(const Float_t* vtx); 
  void LoadClusterFiredChips(TTree* tree);
  void FlagClustersInOverlapRegions(Int_t ic1,Int_t ic2);
  void FlagTrackClusters(Int_t id);
  void FlagIfSecondary(AliESDtrack* track, const AliVertex* vtx);
  void FlagV0s(const AliESDVertex *vtx);
  void ProcessESDTracks();
  Bool_t  CanBeElectron(const AliESDtrack* trc) const;
  
  void CreateMultiplicityObject();
  //
  void SetPhiWindow(Float_t w=0.08) {fPhiWindow=w;}
  void SetThetaWindow(Float_t w=0.025) {fThetaWindow=w;}
  void SetPhiShift(Float_t w=0.0045) {fPhiShift=w;}
  void SetRemoveClustersFromOverlaps(Bool_t b = kFALSE) {fRemoveClustersFromOverlaps = b;}
  void SetPhiOverlapCut(Float_t w=0.005) {fPhiOverlapCut=w;}
  void SetZetaOverlapCut(Float_t w=0.05) {fZetaOverlapCut=w;}
  void SetPhiRotationAngle(Float_t w=0.0) {fPhiRotationAngle=w;}

  Int_t GetNClustersLayer1() const {return fNClustersLay1;}
  Int_t GetNClustersLayer2() const {return fNClustersLay2;}
  Int_t GetNTracklets() const {return fNTracklets;}
  Int_t GetNSingleClusters() const {return fNSingleCluster;}
  Short_t GetNFiredChips(Int_t layer) const {return fNFiredChips[layer];}

  Float_t* GetClusterLayer1(Int_t n) {return &fClustersLay1[n*kClNPar];}
  Float_t* GetClusterLayer2(Int_t n) {return &fClustersLay2[n*kClNPar];}

  Float_t* GetTracklet(Int_t n) {return fTracklets[n];}
  Float_t* GetCluster(Int_t n) {return fSClusters[n];}

  void SetHistOn(Bool_t b=kFALSE) {fHistOn=b;}
  void SaveHists();

  AliITSDetTypeRec *GetDetTypeRec() const {return fDetTypeRec;}
  void SetDetTypeRec(AliITSDetTypeRec *ptr){fDetTypeRec = ptr;}
  //
  void    SetCutPxDrSPDin(Float_t v=0.1)             { fCutPxDrSPDin = v;}
  void    SetCutPxDrSPDout(Float_t v=0.15)           { fCutPxDrSPDout = v;}
  void    SetCutPxDz(Float_t v=0.2)                  { fCutPxDz = v;}
  void    SetCutDCArz(Float_t v=0.5)                 { fCutDCArz = v;}
  void    SetCutMinElectronProbTPC(Float_t v=0.5)    { fCutMinElectronProbTPC = v;}
  void    SetCutMinElectronProbESD(Float_t v=0.1)    { fCutMinElectronProbESD = v;}
  void    SetCutMinP(Float_t v=0.05)                 { fCutMinP = v;}
  void    SetCutMinRGamma(Float_t v=2.)              { fCutMinRGamma = v;}
  void    SetCutMinRK0(Float_t v=1.)                 { fCutMinRK0 = v;}
  void    SetCutMinPointAngle(Float_t v=0.98)        { fCutMinPointAngle = v;}
  void    SetCutMaxDCADauther(Float_t v=0.5)         { fCutMaxDCADauther = v;}
  void    SetCutMassGamma(Float_t v=0.03)            { fCutMassGamma = v;}
  void    SetCutMassGammaNSigma(Float_t v=5.)        { fCutMassGammaNSigma = v;}
  void    SetCutMassK0(Float_t v=0.03)               { fCutMassK0 = v;}
  void    SetCutMassK0NSigma(Float_t v=5.)           { fCutMassK0NSigma = v;}
  void    SetCutChi2cGamma(Float_t v=2.)             { fCutChi2cGamma = v;}
  void    SetCutChi2cK0(Float_t v=2.)                { fCutChi2cK0 = v;}
  void    SetCutGammaSFromDecay(Float_t v=-10.)      { fCutGammaSFromDecay = v;}
  void    SetCutK0SFromDecay(Float_t v=-10.)         { fCutK0SFromDecay = v;}
  void    SetCutMaxDCA(Float_t v=1.)                 { fCutMaxDCA = v;}
  //
  Float_t GetCutPxDrSPDin()                    const {return fCutPxDrSPDin;}
  Float_t GetCutPxDrSPDout()                   const {return fCutPxDrSPDout;}
  Float_t GetCutPxDz()                         const {return fCutPxDz;}
  Float_t GetCutDCArz()                        const {return fCutDCArz;}
  Float_t GetCutMinElectronProbTPC()           const {return fCutMinElectronProbTPC;}
  Float_t GetCutMinElectronProbESD()           const {return fCutMinElectronProbESD;}
  Float_t GetCutMinP()                         const {return fCutMinP;}
  Float_t GetCutMinRGamma()                    const {return fCutMinRGamma;}
  Float_t GetCutMinRK0()                       const {return fCutMinRK0;}
  Float_t GetCutMinPointAngle()                const {return fCutMinPointAngle;}
  Float_t GetCutMaxDCADauther()                const {return fCutMaxDCADauther;}
  Float_t GetCutMassGamma()                    const {return fCutMassGamma;}
  Float_t GetCutMassGammaNSigma()              const {return fCutMassGammaNSigma;}
  Float_t GetCutMassK0()                       const {return fCutMassK0;}
  Float_t GetCutMassK0NSigma()                 const {return fCutMassK0NSigma;}
  Float_t GetCutChi2cGamma()                   const {return fCutChi2cGamma;}
  Float_t GetCutChi2cK0()                      const {return fCutChi2cK0;}
  Float_t GetCutGammaSFromDecay()              const {return fCutGammaSFromDecay;}
  Float_t GetCutK0SFromDecay()                 const {return fCutK0SFromDecay;}
  Float_t GetCutMaxDCA()                       const {return fCutMaxDCA;}

  //
protected:
  AliTrackletAlg(const AliTrackletAlg& mr);
  AliTrackletAlg& operator=(const AliTrackletAlg& mr);
  AliITSDetTypeRec* fDetTypeRec;            //! pointer to DetTypeRec
  AliESDEvent*      fESDEvent;              //! pointer to ESD event
  TTree*            fTreeRP;                //! ITS recpoints

  UInt_t*       fUsedClusLay1;               // RS: flag of clusters usage in ESD tracks: 0=unused, else ID+1 in word0=TPC/ITS+ITSSA, word1=ITSSA_Pure
  UInt_t*       fUsedClusLay2;               // RS: flag of clusters usage in ESD tracks: 0=unused, else ID+1 word0=TPC/ITS+ITSSA, word1=ITSSA_Pure

  Float_t*      fClustersLay1;               // clusters in the 1st layer of ITS 
  Float_t*      fClustersLay2;               // clusters in the 2nd layer of ITS 
  Int_t*        fDetectorIndexClustersLay1;  // module index for clusters 1st ITS layer
  Int_t*        fDetectorIndexClustersLay2;  // module index for clusters 2nd ITS layer
  Bool_t*       fOverlapFlagClustersLay1;    // flag for clusters in the overlap regions 1st ITS layer
  Bool_t*       fOverlapFlagClustersLay2;    // flag for clusters in the overlap regions 2nd ITS layer 

  Float_t**     fTracklets;            // tracklets 
  Float_t**     fSClusters;            // single clusters (unassociated)
  
  Int_t         fNClustersLay1;        // Number of clusters (Layer1)
  Int_t         fNClustersLay2;        // Number of clusters (Layer2)
  Int_t         fNTracklets;           // Number of tracklets
  Int_t         fNSingleCluster;       // Number of unassociated clusters
  Short_t       fNFiredChips[2];       // Number of fired chips in the two SPD layers
  //
  // Following members are set via AliITSRecoParam
  //
  Float_t       fPhiWindow;                    // Search window in phi
  Float_t       fThetaWindow;                  // Search window in theta
  Float_t       fPhiShift;                     // Phi shift reference value (at 0.5 T) 
  Bool_t        fRemoveClustersFromOverlaps;   // Option to skip clusters in the overlaps
  Float_t       fPhiOverlapCut;                // Fiducial window in phi for overlap cut
  Float_t       fZetaOverlapCut;               // Fiducial window in eta for overlap cut
  Float_t       fPhiRotationAngle;             // Angle to rotate the inner layer cluster for combinatorial reco only 

  // cuts for secondaries identification
  Float_t       fCutPxDrSPDin;                 // max P*DR for primaries involving at least 1 SPD
  Float_t       fCutPxDrSPDout;                // max P*DR for primaries not involving any SPD
  Float_t       fCutPxDz;                      // max P*DZ for primaries
  Float_t       fCutDCArz;                     // max DR or DZ for primares
  //
  // cuts for flagging tracks in V0s
  Float_t       fCutMinElectronProbTPC;     // min probability for e+/e- PID involving TPC
  Float_t       fCutMinElectronProbESD;     // min probability for e+/e- PID not involving TPC
  //
  Float_t       fCutMinP;                   // min P of V0
  Float_t       fCutMinRGamma;              // min transv. distance from ESDVertex to V0 for gammas
  Float_t       fCutMinRK0;                 // min transv. distance from ESDVertex to V0 for K0s
  Float_t       fCutMinPointAngle;          // min pointing angle cosine
  Float_t       fCutMaxDCADauther;          // max DCA of daughters at V0
  Float_t       fCutMassGamma;              // max gamma mass
  Float_t       fCutMassGammaNSigma;        // max standard deviations from 0 for gamma
  Float_t       fCutMassK0;                 // max K0 mass difference from PGD value
  Float_t       fCutMassK0NSigma;           // max standard deviations for K0 mass from PDG value
  Float_t       fCutChi2cGamma;             // max constrained chi2 cut for gammas
  Float_t       fCutChi2cK0;                // max constrained chi2 cut for K0s
  Float_t       fCutGammaSFromDecay;        // min path*P for gammas
  Float_t       fCutK0SFromDecay;           // min path*P for K0s
  Float_t       fCutMaxDCA;                 // max DCA for V0 at ESD vertex  

  Bool_t        fHistOn;               // Option to define and fill the histograms 

  TH1F*         fhClustersDPhiAcc;     // Phi2 - Phi1 for tracklets 
  TH1F*         fhClustersDThetaAcc;   // Theta2 - Theta1 for tracklets 
  TH1F*         fhClustersDPhiAll;     // Phi2 - Phi1 all the combinations 
  TH1F*         fhClustersDThetaAll;   // Theta2 - Theta1 all the combinations
 
  TH2F*         fhDPhiVsDThetaAll;     // 2D plot for all the combinations  
  TH2F*         fhDPhiVsDThetaAcc;     // same plot for tracklets 

  TH1F*         fhetaTracklets;        // Pseudorapidity distr. for tracklets 
  TH1F*         fhphiTracklets;        // Azimuthal (Phi) distr. for tracklets  
  TH1F*         fhetaClustersLay1;     // Pseudorapidity distr. for Clusters L. 1
  TH1F*         fhphiClustersLay1;     // Azimuthal (Phi) distr. for Clusters L. 1 


  void LoadClusterArrays(TTree* tree);

  ClassDef(AliTrackletAlg,1) 
};

#endif
