#ifndef ALIITSMULTRECONSTRUCTOR_H
#define ALIITSMULTRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//_________________________________________________________________________
// 
//        Implementation of the ITS-SPD trackleter class
//
// It retrieves clusters in the pixels (theta and phi) and finds tracklets.
// These can be used to extract charged particle multiplicity from the ITS.
//
// A tracklet consists of two ITS clusters, one in the first pixel layer and 
// one in the second. The clusters are associated if the differences in 
// Phi (azimuth) and Theta (polar angle) are within fiducial windows.
// In case of multiple candidates the candidate with minimum
// distance is selected. 
//_________________________________________________________________________
#include "AliTrackleter.h"
#include "AliITSsegmentationSPD.h"
#include "TMath.h"

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
class AliRefArray;

class AliITSMultReconstructor : public AliTrackleter
{
public:
  //
  enum {kClTh,kClPh,kClZ,kClMC0,kClMC1,kClMC2,kClNPar};
  enum {kTrTheta,kTrPhi,kTrDPhi,kTrDTheta,kTrLab1,kTrLab2,kClID1,kClID2,kTrNPar};
  enum {kSCTh,kSCPh,kSCLab,kSCID,kSCNPar};   
  enum {kITSTPC,kITSSAP,kITSTPCBit=BIT(kITSTPC),kITSSAPBit=BIT(kITSSAP)}; // RS
  AliITSMultReconstructor();
  virtual ~AliITSMultReconstructor();

  void Reconstruct(AliESDEvent* esd, TTree* treeRP);
  void Reconstruct(TTree* tree, Float_t* vtx, Float_t* vtxRes=0);   // old reconstructor invocation
  void ReconstructMix(TTree* clusterTree, TTree* clusterTreeMix, const Float_t* vtx, Float_t* vtrRes=0);
  void FindTracklets(const Float_t* vtx); 
  void LoadClusterFiredChips(TTree* tree);
  void FlagClustersInOverlapRegions(Int_t ic1,Int_t ic2);
  void FlagTrackClusters(Int_t id);
  void FlagIfSecondary(AliESDtrack* track, const AliVertex* vtx);
  void FlagV0s(const AliESDVertex *vtx);
  void ProcessESDTracks();
  Bool_t  CanBeElectron(const AliESDtrack* trc) const;
  
  virtual void CreateMultiplicityObject();
  //
  // Following members are set via AliITSRecoParam
  void SetPhiWindow(Float_t w=0.08)    {fDPhiWindow=w;   fDPhiWindow2 = w*w;}
  void SetThetaWindow(Float_t w=0.025) {fDThetaWindow=w; fDThetaWindow2=w*w;}
  void SetPhiShift(Float_t w=0.0045) {fPhiShift=w;}
  void SetRemoveClustersFromOverlaps(Bool_t b = kFALSE) {fRemoveClustersFromOverlaps = b;}
  void SetPhiOverlapCut(Float_t w=0.005) {fPhiOverlapCut=w;}
  void SetZetaOverlapCut(Float_t w=0.05) {fZetaOverlapCut=w;}
  void SetPhiRotationAngle(Float_t w=0.0) {fPhiRotationAngle=w;}

  Int_t GetNClustersLayer1() const {return fNClustersLay[0];}
  Int_t GetNClustersLayer2() const {return fNClustersLay[1];}
  Int_t GetNClustersLayer(Int_t i) const {return fNClustersLay[i];}
  Int_t GetNTracklets() const {return fNTracklets;}
  Int_t GetNSingleClusters() const {return fNSingleCluster;}
  Short_t GetNFiredChips(Int_t layer) const {return fNFiredChips[layer];}

  Float_t* GetClusterLayer1(Int_t n) {return &fClustersLay[0][n*kClNPar];}
  Float_t* GetClusterLayer2(Int_t n) {return &fClustersLay[1][n*kClNPar];}
  Float_t* GetClusterOfLayer(Int_t lr,Int_t n) {return &fClustersLay[lr][n*kClNPar];}
  
  Float_t* GetTracklet(Int_t n) {return fTracklets[n];}
  Float_t* GetCluster(Int_t n) {return fSClusters[n];}

  void     SetScaleDThetaBySin2T(Bool_t v=kTRUE)         {fScaleDTBySin2T = v;}
  Bool_t   GetScaleDThetaBySin2T()               const   {return fScaleDTBySin2T;}
  //
  void     SetNStdDev(Float_t f=1.)                      {fNStdDev = f<0.01 ? 0.01 : f; fNStdDevSq=TMath::Sqrt(fNStdDev);}
  Float_t  GetNStdDev()                          const   {return fNStdDev;}
  //
  void SetHistOn(Bool_t b=kFALSE) {fHistOn=b;}
  void SaveHists();
  //
  void   SetBuildRefs(Bool_t v=kTRUE)                    {fBuildRefs = v;}
  Bool_t GetBuildRefs()                          const   {return fBuildRefs;}

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
  void  InitAux();
  void  ClusterPos2Angles(const Float_t *vtx);
  void  ClusterPos2Angles(Float_t *clPar, const Float_t *vtx) const;
  Int_t AssociateClusterOfL1(Int_t iC1);
  Int_t StoreTrackletForL2Cluster(Int_t iC2);
  void  StoreL1Singles();
  TClonesArray* GetClustersOfLayer(Int_t il)   const {return fClArr[il];}
  void  LoadClusters()                               {LoadClusterArrays(fTreeRP);}
  void  SetTreeRP(TTree* rp)                         {fTreeRP    = rp;}
  void  SetTreeRPMix(TTree* rp=0)                    {fTreeRPMix = rp;}
  Bool_t AreClustersLoaded()                   const {return fClustersLoaded;}
  Bool_t GetCreateClustersCopy()               const {return fCreateClustersCopy;}
  Bool_t IsRecoDone()                          const {return fRecoDone;}
  void  SetCreateClustersCopy(Bool_t v=kTRUE)        {fCreateClustersCopy=v;}
  //
  //  Float_t* GetClustersArray(Int_t lr)          const {return (Float_t*) (lr==0) ? fClustersLay[0]:fClustersLay[1];}
  Float_t* GetClustersArray(Int_t lr)          const {if(lr==0){return fClustersLay[0];} 
                                                      else {return fClustersLay[1];}}
  Int_t*   GetPartnersOfL2()                   const {return (Int_t*)fPartners;}
  Float_t* GetMinDistsOfL2()                   const {return (Float_t*)fMinDists;}
  Double_t GetDPhiShift()                      const {return fDPhiShift;}
  Double_t GetDPhiWindow2()                    const {return fDPhiWindow2;}
  Double_t GetDThetaWindow2()                  const {return fDThetaWindow2;}
  Double_t CalcDist(Double_t dphi, Double_t dtheta, Double_t theta) const; 
  //
 protected:
  void   SetClustersLoaded(Bool_t v=kTRUE)           {fClustersLoaded = v;}
  AliITSMultReconstructor(const AliITSMultReconstructor& mr);
  AliITSMultReconstructor& operator=(const AliITSMultReconstructor& mr);
  void              CalcThetaPhi(float dx,float dy,float dz,float &theta,float &phi) const;
  AliITSDetTypeRec* fDetTypeRec;            //! pointer to DetTypeRec
  AliESDEvent*      fESDEvent;              //! pointer to ESD event
  TTree*            fTreeRP;                //! ITS recpoints
  TTree*            fTreeRPMix;             //! ITS recpoints for mixing
  AliRefArray*      fUsedClusLay[2][2];     //! RS: clusters usage in ESD tracks
  //
  Float_t*      fClustersLay[2];            //! clusters in the SPD layers of ITS 
  Int_t*        fDetectorIndexClustersLay[2];  //! module index for clusters in ITS layers
  Bool_t*       fOverlapFlagClustersLay[2];  //! flag for clusters in the overlap regions in ITS layers

  Float_t**     fTracklets;            //! tracklets 
  Float_t**     fSClusters;            //! single clusters (unassociated)
  
  Int_t         fNClustersLay[2];      // Number of clusters on each layer
  Int_t         fNTracklets;           // Number of tracklets
  Int_t         fNSingleCluster;       // Number of unassociated clusters
  Short_t       fNFiredChips[2];       // Number of fired chips in the two SPD layers
  //
  // Following members are set via AliITSRecoParam
  //
  Float_t       fDPhiWindow;                   // Search window in phi
  Float_t       fDThetaWindow;                 // Search window in theta
  Float_t       fPhiShift;                     // Phi shift reference value (at 0.5 T) 
  Bool_t        fRemoveClustersFromOverlaps;   // Option to skip clusters in the overlaps
  Float_t       fPhiOverlapCut;                // Fiducial window in phi for overlap cut
  Float_t       fZetaOverlapCut;               // Fiducial window in eta for overlap cut
  Float_t       fPhiRotationAngle;             // Angle to rotate the inner layer cluster for combinatorial reco only 
  //
  Bool_t        fScaleDTBySin2T;               // use in distance definition
  Float_t       fNStdDev;                      // number of standard deviations to keep
  Float_t       fNStdDevSq;                    // sqrt of number of standard deviations to keep
  //
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

  // temporary stuff for single event trackleting
  Double_t      fDPhiShift;            // shift in dphi due to the curvature
  Double_t      fDPhiWindow2;          // phi window^2
  Double_t      fDThetaWindow2;        // theta window^2
  Int_t*        fPartners;             //! L2 partners of L1
  Int_t*        fAssociatedLay1;       //! association flag
  Float_t*      fMinDists;             //! smallest distances for L2->L1
  AliRefArray*  fBlackList;            //! blacklisted cluster references
  Bool_t        fStoreRefs[2][2];      //! which cluster to track refs to store
  //
  // this is for the analysis mode only
  TClonesArray  *fClArr[2];            //! original clusters
  Bool_t        fCreateClustersCopy;   //  read and clone clusters directly from the tree
  Bool_t        fClustersLoaded;       // flag of clusters loaded
  Bool_t        fRecoDone;             // flag that reconstruction is done
  Bool_t        fBuildRefs;            // build cluster to tracks references
  //
  AliITSsegmentationSPD fSPDSeg;       // SPD segmentation model
  //
  void LoadClusterArrays(TTree* tree, TTree* treeMix=0);
  void LoadClusterArrays(TTree* tree,int il);

  ClassDef(AliITSMultReconstructor,11)
};

//____________________________________________________________________
inline void AliITSMultReconstructor::ClusterPos2Angles(Float_t *clPar, const Float_t *vtx) const
{
  // convert cluster coordinates to angles wrt vertex
  Float_t x = clPar[kClTh] - vtx[0];
  Float_t y = clPar[kClPh] - vtx[1];
  Float_t z = clPar[kClZ]  - vtx[2];
  Float_t r    = TMath::Sqrt(x*x + y*y + z*z);
  clPar[kClTh] = TMath::ACos(z/r);                   // Store Theta
  clPar[kClPh] = TMath::Pi() + TMath::ATan2(-y,-x);  // Store Phi 
  //
}

//____________________________________________________________________
inline Double_t AliITSMultReconstructor::CalcDist(Double_t dphi, Double_t dtheta, Double_t theta) const
{
  // calculate eliptical distance. theta is the angle of cl1, dtheta = tht(cl1)-tht(cl2)
  dphi = TMath::Abs(dphi) - fDPhiShift;
  if (fScaleDTBySin2T) {
    double sinTI = TMath::Sin(theta-dtheta/2);
    sinTI *= sinTI;
    dtheta /= sinTI>1.e-6 ? sinTI : 1.e-6;
  }
  return dphi*dphi/fDPhiWindow2 + dtheta*dtheta/fDThetaWindow2;
}

//____________________________________________________________________
inline void AliITSMultReconstructor::CalcThetaPhi(float x, float y,float z,float &theta,float &phi) const
{
  // get theta and phi in tracklet convention
  theta = TMath::ACos(z/TMath::Sqrt(x*x + y*y + z*z));
  phi   = TMath::Pi() + TMath::ATan2(-y,-x);
}


#endif
