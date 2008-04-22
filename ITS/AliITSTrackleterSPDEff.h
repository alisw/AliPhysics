#ifndef ALIITSTRACKLETERSPDEFF_H
#define ALIITSTRACKLETERSPDEFF_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//____________________________________________________________________
// 
// AliITSTrackleterSPDEff - find SPD chips efficiencies by using tracklets.
// 
// This class has been derived from AliITSMultReconstructor (see
// it for more details). It is the class for the Trackleter used to estimate
// SPD plane efficiency.
// The trackleter prediction is built using the vertex and 1 cluster.

//
// 
//  Author :  Giuseppe Eugenio Bruno, based on the skeleton of Reconstruct method  provided by Tiziano Virgili
//  email:    giuseppe.bruno@ba.infn.it
//  
//____________________________________________________________________

#include "AliITSMultReconstructor.h"
#include "AliITSPlaneEffSPD.h"

class AliStack;

class AliITSTrackleterSPDEff : public AliITSMultReconstructor 
{
public:
  AliITSTrackleterSPDEff();
  virtual ~AliITSTrackleterSPDEff();

  void Reconstruct(TTree* tree, Float_t* vtx, Float_t* vtxRes, AliStack* pStack=0x0);

  void SetPhiWindowL1(Float_t w=0.08) {fPhiWindowL1=w;}
  void SetZetaWindowL1(Float_t w=1.) {fZetaWindowL1=w;}
  void SetOnlyOneTrackletPerC1(Bool_t b = kTRUE) {fOnlyOneTrackletPerC1 = b;}
  
  AliITSPlaneEffSPD* GetPlaneEff() const {return fPlaneEffSPD;}
  
  void SetMC(Bool_t mc=kTRUE) {fMC=mc; InitPredictionMC(); return;}
  Bool_t GetMC() const {return fMC;}
  void SetUseOnlyPrimaryForPred(Bool_t flag=kTRUE) {CallWarningMC(); fUseOnlyPrimaryForPred = flag; }
  void SetUseOnlySecondaryForPred(Bool_t flag=kTRUE) {CallWarningMC(); fUseOnlySecondaryForPred = flag;}
  void SetUseOnlySameParticle(Bool_t flag=kTRUE) {CallWarningMC(); fUseOnlySameParticle = flag;}
  void SetUseOnlyDifferentParticle(Bool_t flag=kTRUE) {CallWarningMC(); fUseOnlyDifferentParticle = flag;}
  void SetUseOnlyStableParticle(Bool_t flag=kTRUE) {CallWarningMC(); fUseOnlyStableParticle = flag;}
  Bool_t GetUseOnlyPrimaryForPred() const {CallWarningMC(); return fUseOnlyPrimaryForPred; }
  Bool_t GetUseOnlySecondaryForPred() const {CallWarningMC(); return fUseOnlySecondaryForPred;}
  Bool_t GetUseOnlySameParticle() const {CallWarningMC(); return fUseOnlySameParticle;}
  Bool_t GetUseOnlyDifferentParticle() const {CallWarningMC(); return fUseOnlyDifferentParticle;}
  Bool_t GetUseOnlyStableParticle() const {CallWarningMC(); return fUseOnlyStableParticle;}
  Int_t GetPredictionPrimary(const UInt_t key) const;
  Int_t GetPredictionSecondary(const UInt_t key) const;
  Int_t GetClusterPrimary(const UInt_t key) const;
  Int_t GetClusterSecondary(const UInt_t key) const;
  Int_t GetPredictionPrimary(const UInt_t mod, const UInt_t chip) const
        {return GetPredictionPrimary(fPlaneEffSPD->GetKey(mod,chip));};
  Int_t GetPredictionSecondary(const UInt_t mod, const UInt_t chip) const
        {return GetPredictionSecondary(fPlaneEffSPD->GetKey(mod,chip));};
  Int_t GetClusterPrimary(const UInt_t mod, const UInt_t chip) const
        {return GetClusterPrimary(fPlaneEffSPD->GetKey(mod,chip));};
  Int_t GetClusterSecondary(const UInt_t mod, const UInt_t chip) const
        {return GetClusterSecondary(fPlaneEffSPD->GetKey(mod,chip));};
  void SavePredictionMC(TString filename="TrackletsMCpred.txt") const;
  void ReadPredictionMC(TString filename="TrackletsMCpred.txt");
  // Print some class info in ascii form to stream (cut values and MC statistics)
  virtual void PrintAscii(ostream *os)const;
  // Read some class info in ascii form from stream (cut values and MC statistics)
  virtual void ReadAscii(istream *is);
  Bool_t GetHistOn() const {return fHistOn;}; // return status of histograms
  Bool_t WriteHistosToFile(TString filename="TrackleterSPDHistos.root",Option_t* option = "RECREATE");
  void SetHistOn(Bool_t his=kTRUE) {AliITSMultReconstructor::SetHistOn(his); 
         if(GetHistOn()) {DeleteHistos(); BookHistos();} else DeleteHistos(); return;}

protected:
  AliITSTrackleterSPDEff(const AliITSTrackleterSPDEff& mr);
  AliITSTrackleterSPDEff& operator=(const AliITSTrackleterSPDEff& mr);

  Bool_t*       fAssociationFlag1;    // flag for the associations (Layer 1)
  UInt_t*       fChipPredOnLay2;      // prediction for the chip traversed by the tracklet 
                                      // based on vtx and ClusterLay1 (to be used in extrapolation)
  UInt_t*       fChipPredOnLay1;      // prediction for the chip traversed by the tracklet 
                                      // based on vtx and ClusterLay2 (to be used in interpolation)
  Int_t         fNTracklets1;   // Number of tracklets layer 1
  Float_t       fPhiWindowL1;     // Search window in phi (Layer 1)
  Float_t       fZetaWindowL1;    // SEarch window in zeta (Layer 1)
  Bool_t        fOnlyOneTrackletPerC1; // only one tracklet per cluster in L. 1
  AliITSPlaneEffSPD* fPlaneEffSPD; // pointer to SPD plane efficiency class
  Bool_t   fMC; // Boolean to access Kinematics (only for MC events )
  Bool_t   fUseOnlyPrimaryForPred; // Only for MC: if this is true, build tracklet prediction using only primary particles
  Bool_t   fUseOnlySecondaryForPred; // Only for MC: if this is true build tracklet prediction using only secondary particles
  Bool_t   fUseOnlySameParticle; // Only for MC: if this is true, assign a success only if clusters from same particles 
                                 // (i.e. PP or SS) otherwise ignore the combination
  Bool_t   fUseOnlyDifferentParticle; // Only for MC: if this is true, assign a success only if clusters from different particles 
                                      // (i.e. PP' or PS or SS') otherwise ignore the combination
  Bool_t   fUseOnlyStableParticle; // Only for MC: if this is kTRUE then method PrimaryTrackChecker return kTRUE only 
                                //              for particles decaying (eventually) after pixel layers
  Int_t *fPredictionPrimary;  // those for correction of bias from secondaries
  Int_t *fPredictionSecondary; // chip_by_chip: number of Prediction built with primaries/secondaries
  Int_t *fClusterPrimary;  //   number of clusters on a given chip fired by (at least) a primary
  Int_t *fClusterSecondary; //  number of clusters on a given chip fired by (only) secondaries
 // extra histograms with respect to the base class AliITSMultReconstructor
  TH1F*         fhClustersDPhiInterpAcc;   // Phi2 - Phi1 for tracklets (interpolation phase)
  TH1F*         fhClustersDThetaInterpAcc; // Theta2 - Theta1 for tracklets (interpolation phase)
  TH1F*         fhClustersDZetaInterpAcc;  // z2 - z1projected for tracklets (interpolation phase)
  TH1F*         fhClustersDPhiInterpAll;   // Phi2 - Phi1 all the combinations (interpolation phase)
  TH1F*         fhClustersDThetaInterpAll; // Theta2 - Theta1 all the combinations (interpolation phase)
  TH1F*         fhClustersDZetaInterpAll;  // z2 - z1projected all the combinations (interpolation phase)
  TH2F*         fhDPhiVsDThetaInterpAll; // 2D plot for all the combinations
  TH2F*         fhDPhiVsDThetaInterpAcc; // same plot for tracklets
  TH2F*         fhDPhiVsDZetaInterpAll;  // 2d plot for all the combination
  TH2F*         fhDPhiVsDZetaInterpAcc;  // same plot for tracklets
  TH1F*         fhetaClustersLay2; // Pseudorapidity distr. for Clusters L. 2
  TH1F*         fhphiClustersLay2; // Azimuthal (Phi) distr. for Clusters L. 2
//
  Double_t GetRLayer(Int_t layer); // return average radius of layer (0,1) from Geometry
  Bool_t PrimaryTrackChecker(Int_t ipart,AliStack* stack=0x0);  // check if a MC particle is primary (need AliStack)
  Int_t DecayingTrackChecker(Int_t ipart,AliStack* stack=0x0);  // For a primary particle, check if it is stable (see cxx)
  void InitPredictionMC();
  // method to locate a chip using current vtx and polar coordinate od tracklet w.r.t. to vtx (zVtx may not be given)
  Bool_t FindChip(UInt_t &key, Int_t layer,  Float_t* vtx, Float_t thetaVtx, Float_t phiVtx, Float_t zVtx=999.); 
  // method to transform from Global Cilindrical coordinate to local (module) Cartesian coordinate
  Bool_t FromGloCilToLocCart(Int_t ilayer,Int_t idet, Double_t r, Double_t phi, Double_t z,
                           Float_t &xloc, Float_t &zloc);
  // method to obtain the module (detector) index using global coordinates
  Int_t FindDetectorIndex(Int_t layer, Double_t phi, Double_t z);
  // this method gives you the intersections between a line and a circle (centred in the origin) 
  // using polar coordinates
  Bool_t FindIntersectionPolar(Double_t vtx[2],Double_t phiVtx, Double_t R,Double_t &phi);
  Bool_t SetAngleRange02Pi(Double_t &angle); // set the range of angle in [0,2pi[ 
  Bool_t SetAngleRange02Pi(Float_t  &angle) 
  {Double_t tmp=(Double_t)angle; Bool_t ret=SetAngleRange02Pi(tmp);angle=(Float_t)tmp;return ret;};  
  void CallWarningMC() const {if(!fMC) AliWarning("You can use this method only for MC! Call SetMC() first");}
  Bool_t SaveHists();
  void BookHistos(); // booking of extra histograms w.r.t. base class
  void DeleteHistos(); //delete histos from memory

  ClassDef(AliITSTrackleterSPDEff,1)
};
// Input and output function for standard C++ input/output (for the cut values and MC statistics).
ostream &operator<<(ostream &os,const AliITSTrackleterSPDEff &s);
istream &operator>>(istream &is, AliITSTrackleterSPDEff &s);
#endif
