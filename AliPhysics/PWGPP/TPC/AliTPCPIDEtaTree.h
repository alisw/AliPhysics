#ifndef ALITPCPIDETATREE_H
#define ALITPCPIDETATREE_H

/*
This task determines the eta dependence of the TPC signal.
For this purpose, only tracks fulfilling some standard quality cuts are taken into account.
The obtained data can be used to derive the functional behaviour of the eta dependence.
Such a function can be plugged into this task to correct for the eta dependence and to see
if there is then really no eta dependence left.

Class written by Benjamin Hess.
Contact: bhess@cern.ch
*/

class TTree;
class TObjArray;
class THnSparse;
class TH2I;

#include "AliTPCPIDBase.h"

class AliTPCPIDEtaTree : public AliTPCPIDBase {
 public:
  enum PIDtype { kMCid = 0, kTPCid = 1, kV0idNoTOF = 2, kTPCandTOFid = 3, kV0idPlusTOFaccepted = 4, kV0idPlusTOFrejected = 5 };
  
  AliTPCPIDEtaTree();
  AliTPCPIDEtaTree(const char *name);
  virtual ~AliTPCPIDEtaTree();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(const Option_t*);
  
  Bool_t GetCorrectdEdxEtaDependence() const { return fCorrectdEdxEtaDependence; };
  void SetCorrectdEdxEtaDependence(Bool_t flag) { fCorrectdEdxEtaDependence = flag; };
  
  Bool_t GetCorrectdEdxMultiplicityDependence() const { return fCorrectdEdxMultiplicityDependence; };
  void SetCorrectdEdxMultiplicityDependence(Bool_t flag) { fCorrectdEdxMultiplicityDependence = flag; };
  
  Bool_t GetDoAdditionalQA() const { return fDoAdditionalQA; };
  void SetDoAdditionalQA(Bool_t doAdditionalQA = kTRUE) { fDoAdditionalQA = doAdditionalQA; };
  
  Bool_t GetStoreMultiplicity() const  { return fStoreMultiplicity; };
  void SetStoreMultiplicity(Bool_t storeMultiplicity = kTRUE) { fStoreMultiplicity = storeMultiplicity; };
  
  Bool_t GetStoreNumOfSubthresholdclusters() const  { return fStoreNumOfSubthresholdclusters; };
  void SetStoreNumOfSubthresholdclusters(Bool_t storeNumOfSubthresholdclusters = kTRUE)
    { fStoreNumOfSubthresholdclusters = storeNumOfSubthresholdclusters; };
    
  Bool_t GetStoreNumClustersInActiveVolume() const  { return fStoreNumClustersInActiveVolume; };
  void SetStoreNumClustersInActiveVolume(Bool_t storeNumClustersInActiveVolume = kTRUE)
    { fStoreNumClustersInActiveVolume = storeNumClustersInActiveVolume; };
  
  Double_t GetPtpcPionCut() const { return fPtpcPionCut; };
  void SetPtpcPionCut(Double_t pTPCpionCut) { fPtpcPionCut = pTPCpionCut; };
  
 private:
  Short_t fNumEtaCorrReqErrorsIssued;  // Number of times the error about eta correction issues have been displayed
  Short_t fNumMultCorrReqErrorsIssued; // Number of times the error about multiplicity correction issues have been displayed
  
  Bool_t fStoreMultiplicity; // Store multiplicity in tree?
  Bool_t fStoreNumOfSubthresholdclusters; // Store number of subthreshold clusters in tree?
  Bool_t fStoreNumClustersInActiveVolume; // Store number of clusters in active volume in tree?
  Bool_t fDoAdditionalQA; // Save output for additional QA, like TOF QA?
  Double_t fPtpcPionCut; // Cut on pions with lower tpc momentum
  
  Double_t fPtpc; // TPC momentum
  Double_t fPt; // Transverse momentum
  Double_t fDeDx; // Measured dE/dx
  Double_t fDeDxExpected; // Expected dE/dx according to parametrisation
  Double_t fTanTheta; // Tangens of (local) theta at TPC inner wall
  //Double_t fSinAlpha; // Sine of (local) phi at TPC inner wall
  //Double_t fY; // Local Y at TPC inner wall
  Double_t fPhiPrime; // Phi prime
  UShort_t fTPCsignalN; // Number of TPC clusters for PID
  UShort_t fTPCsignalNsubthreshold; // Number of TPC subthreshold clusters for PID
  Double_t fNumTPCClustersInActiveVolume; // Number of TPC clusters in active volume
  UChar_t  fPIDtype; // Type of identification (TPC dEdx, V0, ...)
  
  // In case of PbpB
  Int_t fMultiplicity; // Multiplicity in case of PbPb
  
  Bool_t fCorrectdEdxEtaDependence;    // Correct eta dependence for dEdxExpected
  Bool_t fCorrectdEdxMultiplicityDependence;    // Correct multiplicity dependence for dEdxExpected
  
  TTree* fTree; //! data tree
  TTree* fTreePions; //! data tree pions
  
  TObjArray* fOutputContainer; //! Output data container for TOF qa
  THnSparseI* fhTOFqa; //! THnSparse with TOF qa data
  TH2I* fhMultiplicityQA; //! QA histo for multiplicity
  
  AliTPCPIDEtaTree(const AliTPCPIDEtaTree&); // not implemented
  AliTPCPIDEtaTree& operator=(const AliTPCPIDEtaTree&); // not implemented
  
  ClassDef(AliTPCPIDEtaTree, 3); 
};

#endif
