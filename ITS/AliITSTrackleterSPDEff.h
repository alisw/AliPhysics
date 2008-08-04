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

class AliStack;
#include "AliITSMultReconstructor.h"
#include "AliITSPlaneEffSPD.h"

class AliITSTrackleterSPDEff : public AliITSMultReconstructor 
{
public:
  AliITSTrackleterSPDEff();
  virtual ~AliITSTrackleterSPDEff();
  // Main method to perform the trackleter and the SPD efficiency evaluation
  void Reconstruct(TTree* tree, Float_t* vtx, Float_t* vtxRes, AliStack* pStack=0x0, TTree* tRef=0x0);

  void SetReflectClusterAroundZAxisForLayer(Int_t ilayer,Bool_t b=kTRUE){  // method to study residual background:
    if(b) AliInfo(Form("All clusters on layer %d will be rotated by 180 deg around z",ilayer)); 
    if(ilayer==0) fReflectClusterAroundZAxisForLayer0=b;                   // a rotation by 180degree around the Z axis  
    else if(ilayer==1) fReflectClusterAroundZAxisForLayer1=b;              // (x->-x; y->-y) to all RecPoints on a 
    else AliInfo("Nothing done: input argument (ilayer) either 0 or 1");   // given layer is applied. In such a way 
  }                                                                        // you remove all the true tracklets.

  void SetPhiWindowL1(Float_t w=0.08) {fPhiWindowL1=w;}  // method to set the cuts in the interpolation
  void SetZetaWindowL1(Float_t w=1.) {fZetaWindowL1=w;}  // phase; use method of the base class for extrap.
  void SetOnlyOneTrackletPerC1(Bool_t b = kTRUE) {fOnlyOneTrackletPerC1 = b;} // as in the base class but 
									      // for the inner layer
  void SetUpdateOncePerEventPlaneEff(Bool_t b = kTRUE) {fUpdateOncePerEventPlaneEff = b;}
  
  AliITSPlaneEffSPD* GetPlaneEff() const {return fPlaneEffSPD;}  // return a pointer to the AliITSPlaneEffSPD
  
  void SetMC(Bool_t mc=kTRUE) {fMC=mc; fMC? InitPredictionMC() : DeletePredictionMC(); return;}  // switch on access to MC true 
  Bool_t GetMC() const {return fMC;}  // check the access to MC true
  // Only for MC: use only "primary" particles (according to PrimaryTrackChecker) for the tracklet prediction
  void SetUseOnlyPrimaryForPred(Bool_t flag=kTRUE) {CallWarningMC(); fUseOnlyPrimaryForPred = flag; } 
  // Only for MC: use only "secondary" particles (according to PrimaryTrackChecker) for the tracklet prediction
  void SetUseOnlySecondaryForPred(Bool_t flag=kTRUE) {CallWarningMC(); fUseOnlySecondaryForPred = flag;}
  // Only for MC: associate a cluster to the tracklet prediction if  from the same particle
  void SetUseOnlySameParticle(Bool_t flag=kTRUE) {CallWarningMC(); fUseOnlySameParticle = flag;}
  // Only for MC: associate a cluster to the tracklet prediction if  from different particles
  void SetUseOnlyDifferentParticle(Bool_t flag=kTRUE) {CallWarningMC(); fUseOnlyDifferentParticle = flag;}
  //  Only for MC: re-define "primary" a particle if it is also "stable" (according to definition in method DecayingTrackChecker)
  void SetUseOnlyStableParticle(Bool_t flag=kTRUE) {CallWarningMC(); fUseOnlyStableParticle = flag;}
  // only for MC: Getters relative to the above setters
  Bool_t GetUseOnlyPrimaryForPred() const {CallWarningMC(); return fUseOnlyPrimaryForPred; }
  Bool_t GetUseOnlySecondaryForPred() const {CallWarningMC(); return fUseOnlySecondaryForPred;}
  Bool_t GetUseOnlySameParticle() const {CallWarningMC(); return fUseOnlySameParticle;}
  Bool_t GetUseOnlyDifferentParticle() const {CallWarningMC(); return fUseOnlyDifferentParticle;}
  Bool_t GetUseOnlyStableParticle() const {CallWarningMC(); return fUseOnlyStableParticle;}
  // Getters for the data members related to MC true statisitcs (see below)
  Int_t GetPredictionPrimary(const UInt_t key) const;
  Int_t GetPredictionSecondary(const UInt_t key) const;
  Int_t GetClusterPrimary(const UInt_t key) const;
  Int_t GetClusterSecondary(const UInt_t key) const;
  Int_t GetSuccessPP(const UInt_t key) const;
  Int_t GetSuccessTT(const UInt_t key) const;
  Int_t GetSuccessS(const UInt_t key) const;
  Int_t GetSuccessP(const UInt_t key) const;
  Int_t GetFailureS(const UInt_t key) const;
  Int_t GetFailureP(const UInt_t key) const;
  Int_t GetRecons(const UInt_t key) const;
  Int_t GetNonRecons(const UInt_t key) const;
  Int_t GetPredictionPrimary(const UInt_t mod, const UInt_t chip) const
        {return GetPredictionPrimary(fPlaneEffSPD->GetKey(mod,chip));};
  Int_t GetPredictionSecondary(const UInt_t mod, const UInt_t chip) const
        {return GetPredictionSecondary(fPlaneEffSPD->GetKey(mod,chip));};
  Int_t GetClusterPrimary(const UInt_t mod, const UInt_t chip) const
        {return GetClusterPrimary(fPlaneEffSPD->GetKey(mod,chip));};
  Int_t GetClusterSecondary(const UInt_t mod, const UInt_t chip) const
        {return GetClusterSecondary(fPlaneEffSPD->GetKey(mod,chip));};
  Int_t GetSuccessPP(const UInt_t mod, const UInt_t chip) const
        {return GetSuccessPP(fPlaneEffSPD->GetKey(mod,chip));};
  Int_t GetSuccessTT(const UInt_t mod, const UInt_t chip) const
       {return GetSuccessTT(fPlaneEffSPD->GetKey(mod,chip));};
  Int_t GetSuccessS(const UInt_t mod, const UInt_t chip) const
       {return GetSuccessS(fPlaneEffSPD->GetKey(mod,chip));};
  Int_t GetSuccessP(const UInt_t mod, const UInt_t chip) const
       {return GetSuccessP(fPlaneEffSPD->GetKey(mod,chip));};
  Int_t GetFailureS(const UInt_t mod, const UInt_t chip) const
       {return GetFailureS(fPlaneEffSPD->GetKey(mod,chip));};
  Int_t GetFailureP(const UInt_t mod, const UInt_t chip) const
       {return GetFailureP(fPlaneEffSPD->GetKey(mod,chip));};
  Int_t GetRecons(const UInt_t mod, const UInt_t chip) const
       {return GetRecons(fPlaneEffSPD->GetKey(mod,chip));};
  Int_t GetNonRecons(const UInt_t mod, const UInt_t chip) const
       {return GetNonRecons(fPlaneEffSPD->GetKey(mod,chip));};
  // methods to write/reas cuts and MC statistics into/from file 
  // if filename contains  ".root", then data are stored into histograms (->root file). 
  void SavePredictionMC(TString filename="TrackletsMCpred.txt") const;
  void ReadPredictionMC(TString filename="TrackletsMCpred.txt");
  // Print some class info in ascii form to stream (cut values and MC statistics)
  virtual void PrintAscii(ostream *os)const;
  // Read some class info in ascii form from stream (cut values and MC statistics)
  virtual void ReadAscii(istream *is);
  Bool_t GetHistOn() const {return fHistOn;}; // return status of histograms
  // write histograms into a root file on disk
  Bool_t WriteHistosToFile(TString filename="TrackleterSPDHistos.root",Option_t* option = "RECREATE");
  // switch on/off the extra histograms
  void SetHistOn(Bool_t his=kTRUE) {AliITSMultReconstructor::SetHistOn(his); 
         if(GetHistOn()) {DeleteHistos(); BookHistos();} else DeleteHistos(); return;}

protected:
  AliITSTrackleterSPDEff(const AliITSTrackleterSPDEff& mr); // protected method: no copy allowed from outside
  AliITSTrackleterSPDEff& operator=(const AliITSTrackleterSPDEff& mr);

  Bool_t*       fAssociationFlag1;    // flag for the associations (Layer 1)
  UInt_t*       fChipPredOnLay2;      // prediction for the chip traversed by the tracklet 
                                      // based on vtx and ClusterLay1 (to be used in extrapolation)
  UInt_t*       fChipPredOnLay1;      // prediction for the chip traversed by the tracklet 
                                      // based on vtx and ClusterLay2 (to be used in interpolation)
  Int_t         fNTracklets1;   // Number of tracklets layer 1
  // possible cuts :
  Float_t       fPhiWindowL1;     // Search window in phi (Layer 1)
  Float_t       fZetaWindowL1;    // SEarch window in zeta (Layer 1)
  Bool_t        fOnlyOneTrackletPerC1; // only one tracklet per cluster in L. 1
  Bool_t        fUpdateOncePerEventPlaneEff;  //  If this is kTRUE, then you can update the chip efficiency only once
                                              //  per event in that chip. This to avoid double counting from the
                                              //  same tracklets which has two rec-points on one layer.
  Bool_t*       fChipUpdatedInEvent;          //  boolean (chip by chip) to flag which chip has been updated its efficiency
                                              //  in that event
  AliITSPlaneEffSPD* fPlaneEffSPD; // pointer to SPD plane efficiency class
  Bool_t   fReflectClusterAroundZAxisForLayer0;  // if kTRUE, then a 180degree rotation around Z is applied to all 
  Bool_t   fReflectClusterAroundZAxisForLayer1;  // clusters on that layer (x->-x; y->-y)
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
  Int_t *fSuccessPP;     // number of successes by using the same primary track (vs. chip of the success)
  Int_t *fSuccessTT;     // number of successes by using the same track (either a primary or a secondary) (vs. chip of the success)
  Int_t *fSuccessS;      // number of successes by using a secondary for the prediction (vs. chip of the success)
  Int_t *fSuccessP;      // number of successes by using a primary for the prediction (vs. chip of the success)
  Int_t *fFailureS;      // number of failures by using a secondary for the prediction (vs. chip of the failure)
  Int_t *fFailureP;      // number of failures by using a primary for the prediction (vs. chip of the failure)
  Int_t *fRecons;        // number of particle which can be reconstructed (only for MC from TrackRef)
  Int_t *fNonRecons;     // unmber of particle which cannot be reconstructed (only for MC from TrackRef)
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
// check if a MC particle is reconstructable
  Bool_t IsReconstructableAt(Int_t layer,Int_t iC,Int_t ipart,Float_t* vtx,AliStack* stack=0x0,TTree* ref=0x0);
  void InitPredictionMC(); // allocate memory for cuts and MC data memebers
  void DeletePredictionMC(); // deallocate memory
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
  // Method to apply a rotation by 180degree to all RecPoints (x->-x; y->-y) on a given layer
  void ReflectClusterAroundZAxisForLayer(Int_t ilayer); // to be used for backgnd estimation on real data 

  ClassDef(AliITSTrackleterSPDEff,3)
};
// Input and output function for standard C++ input/output (for the cut values and MC statistics).
ostream &operator<<(ostream &os,const AliITSTrackleterSPDEff &s);
istream &operator>>(istream &is, AliITSTrackleterSPDEff &s);
#endif
