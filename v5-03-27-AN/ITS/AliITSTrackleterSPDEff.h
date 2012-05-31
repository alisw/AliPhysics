#ifndef ALIITSTRACKLETERSPDEFF_H
#define ALIITSTRACKLETERSPDEFF_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//____________________________________________________________________
// 
// AliITSTrackleterSPDEff - find SPD chips efficiencies by using tracklets.
// 
// This class was originally derived from AliITSMultReconstructor (see
// it for more details). Later on, the inheritance was changed to AliTracker
// It is the class for the Trackleter used to estimate
// SPD plane efficiency.
// The trackleter prediction is built using the vertex and 1 cluster.

//
// 
//  Author :  Giuseppe Eugenio Bruno, based on the skeleton of Reconstruct method  provided by Tiziano Virgili
//  email:    giuseppe.bruno@ba.infn.it
//  
//____________________________________________________________________

class AliStack;
class TTree;
class TH1F;
class TH2F;
class AliPlaneEff;

#include "AliTracker.h"
#include "AliITSPlaneEffSPD.h"

class AliITSTrackleterSPDEff : public  AliTracker
{
public:
  AliITSTrackleterSPDEff();
  virtual ~AliITSTrackleterSPDEff();
  Int_t Clusters2Tracks(AliESDEvent *esd);
  Int_t PostProcess(AliESDEvent *);

  virtual Int_t PropagateBack(AliESDEvent*) {return 0;}
  virtual Int_t RefitInward(AliESDEvent*) {return 0;}
  Int_t LoadClusters(TTree* cl) {LoadClusterArrays(cl); return 0;} // see implementation in AliITSMultReconstructor
  virtual void UnloadClusters() {return;}
  virtual AliCluster *GetCluster(Int_t) const {return NULL;}

  // Main method to perform the trackleter and the SPD efficiency evaluation
  void Reconstruct(AliStack* pStack=0x0, TTree* tRef=0x0, Bool_t lbkg=kFALSE);

  void SetReflectClusterAroundZAxisForLayer(Int_t ilayer,Bool_t b=kTRUE);  // method to study residual background:
                                                                           // a rotation by 180degree around the Z axis  
                                                                           // (x->-x; y->-y) to all RecPoints on a 
                                                                           // given layer is applied. In such a way 
                                                                           // you remove all the true tracklets.
  void SetLightBkgStudyInParallel(Bool_t b = kTRUE); // if you set this on, then the estimation of the 
						     // SPD efficiency is done as usual for data, but in 
						     // parallel a light (i.e. without control histograms, etc.) 
						     // evaluation of combinatorial background is performed
						     // with the usual ReflectClusterAroundZAxisForLayer method.
  Bool_t GetLightBkgStudyInParallel() const {return fLightBkgStudyInParallel;}
  void SetOnlyOneTrackletPerC2(Bool_t b = kTRUE) {fOnlyOneTrackletPerC2 = b;}
  void SetPhiWindowL2(Float_t w=0.08) {fPhiWindowL2=w;}
  void SetZetaWindowL2(Float_t w=1.) {fZetaWindowL2=w;}

  void SetPhiWindowL1(Float_t w=0.08) {fPhiWindowL1=w;}  // method to set the cuts in the interpolation
  void SetZetaWindowL1(Float_t w=1.) {fZetaWindowL1=w;}  // phase; use method of the base class for extrap.
  void SetOnlyOneTrackletPerC1(Bool_t b = kTRUE) {fOnlyOneTrackletPerC1 = b;} // as in the base class but 
  void SetMinContVtx(Int_t min=3) {fMinContVtx=min;} // set minimum n. of contributors to vertex

  Int_t GetNClustersLayer1() const {return fNClustersLay1;}
  Int_t GetNClustersLayer2() const {return fNClustersLay2;}
  Int_t GetNTracklets() const {return fNTracklets;}

  Float_t* GetClusterLayer1(Int_t n) {return fClustersLay1[n];}
  Float_t* GetClusterLayer2(Int_t n) {return fClustersLay2[n];}
  Float_t* GetTracklet(Int_t n) {return fTracklets[n];}
									      // for the inner layer
  void SetUpdateOncePerEventPlaneEff(Bool_t b = kTRUE) {fUpdateOncePerEventPlaneEff = b;}
  
  AliITSPlaneEffSPD* GetPlaneEffSPD() const {return fPlaneEffSPD;}  // return a pointer to the AliITSPlaneEffSPD
  AliPlaneEff *GetPlaneEff() {return (AliPlaneEff*)fPlaneEffSPD;}   // return the pointer to AliPlaneEff
  
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
  void SavePredictionMC(TString filename="TrackletsMCpred.root") const;
  void ReadPredictionMC(TString filename="TrackletsMCpred.root");
  // Print some class info in ascii form to stream (cut values and MC statistics)
  virtual void PrintAscii(ostream *os)const;
  // Read some class info in ascii form from stream (cut values and MC statistics)
  virtual void ReadAscii(istream *is);
  Bool_t GetHistOn() const {return fHistOn;}; // return status of histograms
  // write histograms into a root file on disk
  Bool_t WriteHistosToFile(TString filename="TrackleterSPDHistos.root",Option_t* option = "RECREATE");
  // switch on/off the extra histograms
  void SetHistOn(Bool_t his=kTRUE) {fHistOn=his; 
         if(GetHistOn()) {DeleteHistos(); BookHistos();} else DeleteHistos(); return;}

protected:
  AliITSTrackleterSPDEff(const AliITSTrackleterSPDEff& mr); // protected method: no copy allowed from outside
  AliITSTrackleterSPDEff& operator=(const AliITSTrackleterSPDEff& mr);
//
//// From AliITSMultReconstructor
//
  Float_t**     fClustersLay1;               //! clusters in the 1st layer of ITS
  Float_t**     fClustersLay2;               //! clusters in the 2nd layer of ITS

  Float_t**     fTracklets;            //! tracklets
  Bool_t*       fAssociationFlag;      //! flag for the associations

  Int_t         fNClustersLay1;        // Number of clusters (Layer1)
  Int_t         fNClustersLay2;        // Number of clusters (Layer2)
  Int_t         fNTracklets;           // Number of tracklets

  // Following members are set via AliITSRecoParam
  Bool_t        fOnlyOneTrackletPerC2;         // Allow only one tracklet per cluster in the outer layer
  Float_t       fPhiWindowL2;                    // Search window in phi
  Float_t       fZetaWindowL2;                   // Search window in eta
  Float_t       fPhiOverlapCut;                // Fiducial window in phi for overlap cut
  Float_t       fZetaOverlapCut;               // Fiducial window in eta for overlap cut

  Bool_t        fHistOn;               // Option to define and fill the histograms

  TH1F*         fhClustersDPhiAcc;     //! Phi2 - Phi1 for tracklets
  TH1F*         fhClustersDThetaAcc;   //! Theta2 - Theta1 for tracklets
  TH1F*         fhClustersDZetaAcc;    //! z2 - z1projected for tracklets
  TH1F*         fhClustersDPhiAll;     //! Phi2 - Phi1 all the combinations
  TH1F*         fhClustersDThetaAll;   //! Theta2 - Theta1 all the combinations
  TH1F*         fhClustersDZetaAll;    //! z2 - z1projected all the combinations

  TH2F*         fhDPhiVsDThetaAll;     //! 2D plot for all the combinations
  TH2F*         fhDPhiVsDThetaAcc;     //! same plot for tracklets
  TH2F*         fhDPhiVsDZetaAll;      //! 2d plot for all the combination
  TH2F*         fhDPhiVsDZetaAcc;      //! same plot for tracklets

  TH1F*         fhetaTracklets;        //! Pseudorapidity distr. for tracklets
  TH1F*         fhphiTracklets;        //! Azimuthal (Phi) distr. for tracklets
  TH1F*         fhetaClustersLay1;     //! Pseudorapidity distr. for Clusters L. 1
  TH1F*         fhphiClustersLay1;     //! Azimuthal (Phi) distr. for Clusters L. 1
//
// 
  Bool_t*       fAssociationFlag1;    //! flag for the associations (Layer 1)
  UInt_t*       fChipPredOnLay2;      //! prediction for the chip traversed by the tracklet 
                                      //  based on vtx and ClusterLay1 (to be used in extrapolation)
  UInt_t*       fChipPredOnLay1;      //! prediction for the chip traversed by the tracklet 
                                      // based on vtx and ClusterLay2 (to be used in interpolation)
  Int_t         fNTracklets1;   // Number of tracklets layer 1
  // possible cuts :
  Float_t       fPhiWindowL1;     // Search window in phi (Layer 1)
  Float_t       fZetaWindowL1;    // SEarch window in zeta (Layer 1)
  Bool_t        fOnlyOneTrackletPerC1; // only one tracklet per cluster in L. 1
  Bool_t        fUpdateOncePerEventPlaneEff;  //  If this is kTRUE, then you can update the chip efficiency only once
  Int_t		fMinContVtx;  // minimum number of contributors (tracklets) to the vertex for the event to be used 
                                              //  per event in that chip. This to avoid double counting from the
                                              //  same tracklets which has two rec-points on one layer.
  Bool_t*       fChipUpdatedInEvent;          //!  boolean (chip by chip) to flag which chip has been updated its efficiency
                                              //  in that event
  AliITSPlaneEffSPD* fPlaneEffSPD; //! pointer to SPD plane efficiency class
  AliITSPlaneEffSPD* fPlaneEffBkg; //! pointer to SPD plane efficiency class for background evaluation
  Bool_t   fReflectClusterAroundZAxisForLayer0;  // if kTRUE, then a 180degree rotation around Z is applied to all 
  Bool_t   fReflectClusterAroundZAxisForLayer1;  // clusters on that layer (x->-x; y->-y)
  Bool_t   fLightBkgStudyInParallel; // if this is kTRUE, the basic and correct evaluation of background is performed
                                     // in paralell to standard SPD efficiency evaluation
  Bool_t   fMC; // Boolean to access Kinematics (only for MC events )
  Bool_t   fUseOnlyPrimaryForPred; // Only for MC: if this is true, build tracklet prediction using only primary particles
  Bool_t   fUseOnlySecondaryForPred; // Only for MC: if this is true build tracklet prediction using only secondary particles
  Bool_t   fUseOnlySameParticle; // Only for MC: if this is true, assign a success only if clusters from same particles 
                                 // (i.e. PP or SS) otherwise ignore the combination
  Bool_t   fUseOnlyDifferentParticle; // Only for MC: if this is true, assign a success only if clusters from different particles 
                                      // (i.e. PP' or PS or SS') otherwise ignore the combination
  Bool_t   fUseOnlyStableParticle; // Only for MC: if this is kTRUE then method PrimaryTrackChecker return kTRUE only 
                                //              for particles decaying (eventually) after pixel layers
  Int_t *fPredictionPrimary;  //! those for correction of bias from secondaries
  Int_t *fPredictionSecondary; //! chip_by_chip: number of Prediction built with primaries/secondaries
  Int_t *fClusterPrimary;  //!   number of clusters on a given chip fired by (at least) a primary
  Int_t *fClusterSecondary; //!  number of clusters on a given chip fired by (only) secondaries
  Int_t *fSuccessPP;     //! number of successes by using the same primary track (vs. chip of the success)
  Int_t *fSuccessTT;     //! number of successes by using the same track (either a primary or a secondary) (vs. chip of the success)
  Int_t *fSuccessS;      //! number of successes by using a secondary for the prediction (vs. chip of the success)
  Int_t *fSuccessP;      //! number of successes by using a primary for the prediction (vs. chip of the success)
  Int_t *fFailureS;      //! number of failures by using a secondary for the prediction (vs. chip of the failure)
  Int_t *fFailureP;      //! number of failures by using a primary for the prediction (vs. chip of the failure)
  Int_t *fRecons;        //! number of particle which can be reconstructed (only for MC from TrackRef)
  Int_t *fNonRecons;     //! unmber of particle which cannot be reconstructed (only for MC from TrackRef)
 // extra histograms with respect to the base class AliITSMultReconstructor
  TH1F*         fhClustersDPhiInterpAcc;   //! Phi2 - Phi1 for tracklets (interpolation phase)
  TH1F*         fhClustersDThetaInterpAcc; //! Theta2 - Theta1 for tracklets (interpolation phase)
  TH1F*         fhClustersDZetaInterpAcc;  //! z2 - z1projected for tracklets (interpolation phase)
  TH1F*         fhClustersDPhiInterpAll;   //! Phi2 - Phi1 all the combinations (interpolation phase)
  TH1F*         fhClustersDThetaInterpAll; //! Theta2 - Theta1 all the combinations (interpolation phase)
  TH1F*         fhClustersDZetaInterpAll;  //! z2 - z1projected all the combinations (interpolation phase)
  TH2F*         fhDPhiVsDThetaInterpAll; //! 2D plot for all the combinations
  TH2F*         fhDPhiVsDThetaInterpAcc; //! same plot for tracklets
  TH2F*         fhDPhiVsDZetaInterpAll;  //! 2d plot for all the combination
  TH2F*         fhDPhiVsDZetaInterpAcc;  //! same plot for tracklets
  TH1F*         fhetaClustersLay2; //! Pseudorapidity distr. for Clusters L. 2
  TH1F*         fhphiClustersLay2; //! Azimuthal (Phi) distr. for Clusters L. 2
  TH1F*         fhClustersInChip; //! number of fired clusters versus chip number [0,1199]
  TH2F**        fhClustersInModuleLay1; //! distribution of cluster in the module Lay 1 (sub-chip scale)
  TH2F**        fhClustersInModuleLay2; //! distribution of cluster in the module Lay 2 (sub-chip scale)
//
  void Init(); // initialize pointers and allocate memory 
  Double_t GetRLayer(Int_t layer); // return average radius of layer (0,1) from Geometry
  Bool_t PrimaryTrackChecker(Int_t ipart,AliStack* stack=0x0);  // check if a MC particle is primary (need AliStack)
  Int_t DecayingTrackChecker(Int_t ipart,AliStack* stack=0x0);  // For a primary particle, check if it is stable (see cxx)
// check if a MC particle is reconstructable
  Bool_t IsReconstructableAt(Int_t layer,Int_t iC,Int_t ipart,const Float_t* vtx,const AliStack* stack=0x0,TTree* ref=0x0);
  void InitPredictionMC(); // allocate memory for cuts and MC data memebers
  void DeletePredictionMC(); // deallocate memory
  // method to locate a chip using current vtx and polar coordinate od tracklet w.r.t. to vtx (zVtx may not be given)
  Bool_t FindChip(UInt_t &key, Int_t layer,const Float_t* vtx, Float_t thetaVtx, Float_t phiVtx, Float_t zVtx=999.); 
  // method to transform from Global Cilindrical coordinate to local (module) Cartesian coordinate
  Bool_t FromGloCilToLocCart(Int_t ilayer,Int_t idet, Double_t r, Double_t phi, Double_t z,
                           Float_t &xloc, Float_t &zloc);
  // method to obtain the module (detector) index using global coordinates
  Int_t FindDetectorIndex(Int_t layer, Double_t phi, Double_t z);
  // this method gives you the intersections between a line and a circle (centred in the origin) 
  // using polar coordinates
  Bool_t FindIntersectionPolar(Double_t vtx[2],Double_t phiVtx, Double_t R,Double_t &phi);
  Bool_t SetAngleRange02Pi(Double_t &angle) const; // set the range of angle in [0,2pi[ 
  Bool_t SetAngleRange02Pi(Float_t  &angle) const 
  {Double_t tmp=(Double_t)angle; Bool_t ret=SetAngleRange02Pi(tmp);angle=(Float_t)tmp;return ret;};  
  void CallWarningMC() const {if(!fMC) AliWarning("You can use this method only for MC! Call SetMC() first");}
  Bool_t SaveHists();
  void BookHistos(); // booking of extra histograms w.r.t. base class
  void DeleteHistos(); //delete histos from memory
  // Method to apply a rotation by 180degree to all RecPoints (x->-x; y->-y) on a given layer
  void ReflectClusterAroundZAxisForLayer(Int_t ilayer); // to be used for backgnd estimation on real data 

  void LoadClusterArrays(TTree* tree);

  ClassDef(AliITSTrackleterSPDEff,6)
};
// Input and output function for standard C++ input/output (for the cut values and MC statistics).
ostream &operator<<(ostream &os,const AliITSTrackleterSPDEff &s);
istream &operator>>(istream &is, AliITSTrackleterSPDEff &s);
#endif

