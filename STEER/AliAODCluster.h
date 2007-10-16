#ifndef AliAODCluster_H
#define AliAODCluster_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//     AOD cluster base class
//     Author: Markus Oldenburg, CERN
//-------------------------------------------------------------------------

#include <TObject.h>

class AliAODCluster : public TObject {

 public:
  
  enum AODClu_t {kUndef = -1, 
		 kPHOSNeutral, 
		 kPHOSCharged,
		 kEMCALClusterv1,
		 kPMDNeutral, 
		 kPMDCharged};

  enum AODCluPID_t {
    kUnknown = 0, 
    kPhoton  = 1, 
    kPi0     = 2, 
    kNeutron = 3, 
    kKaon0   = 4,
    kEleCon  = 5, 
    kCharged = 6, 
    kNeutral = 7 , 
    kOther   = 8};

  AliAODCluster();
  AliAODCluster(Int_t id,
		UInt_t nLabel,
		Int_t *label,
		Double_t energy,
		Double_t x[3],
		Double_t pid[9],
		Char_t ttype=kUndef,
		UInt_t selectInfo=0);

   AliAODCluster(Int_t id,
		 UInt_t nLabel,
		 Int_t *label,
		 Float_t energy,
		 Float_t x[3],
		 Float_t pid[9],
		 Char_t ttype=kUndef,
		 UInt_t selectInfo=0);

  virtual ~AliAODCluster();
  AliAODCluster(const AliAODCluster& clus); 
  AliAODCluster& operator=(const AliAODCluster& clus);

  Double_t Chi2() const { return fChi2; }

  virtual Double_t E() const { return fEnergy; }

  // PID
  virtual const Double_t *PID() const { return fPID; }
  AODCluPID_t GetMostProbablePID() const;
 
  template <class T> void GetPID(T *pid) const {
    for(Int_t i=0; i<9; ++i) pid[i]=fPID[i];}
 
  template <class T> void SetPID(const T *pid) {
    if(pid) for(Int_t i=0; i<9; ++i) fPID[i]=pid[i];
    else {for(Int_t i=0; i<9; fPID[i++]=0);} fPID[AliAODCluster::kUnknown]=1.;}

  Int_t  GetID() const { return fID; }
  Int_t  GetLabel(UInt_t i) const;
  UInt_t GetNLabel() const { return (UInt_t)fNLabel; }
  Bool_t TestFilterBit(UInt_t filterBit) const { return (Bool_t) ((filterBit & fFilterMap) != 0); }
  Char_t GetType() const { return fType; }

  template <class T> Bool_t GetPosition(T *x) const {
    x[0]=fPosition[0]; x[1]=fPosition[1]; x[2]=fPosition[2];
    return kTRUE;}

  Bool_t IsEMCALCluster() {if(fType == kEMCALClusterv1) return kTRUE;
    else return kFALSE;}
  Bool_t IsPHOSCluster() {if(fType == kPHOSCharged || fType == kPHOSNeutral) return kTRUE;
    else return kFALSE;}

  // print
  void  Print(const Option_t *opt = "") const;

  // setters
  void SetID(Int_t id) { fID = id; }
  void SetType(AODClu_t ttype) { fType=ttype; }
  void SetLabel(Int_t *label, UInt_t size);  
  void RemoveLabel();
 
  template <class T> void SetPosition(const T *x);

  void SetChi2(Double_t chi2) { fChi2 = chi2; }

 private :

  // Energy & position
  Double32_t    fEnergy;         // energy
  Double32_t    fPosition[3];    // position of the cluster

  Double32_t    fChi2;           // chi2 (probably not necessary for PMD)
  Double32_t    fPID[9];         // [0.,1.,8] pointer to PID object

  Int_t         fID;             // unique cluster ID, points back to the ESD cluster
  Int_t         fNLabel;         // number of original track for this cluster      
  Int_t        *fLabel;          // [fNLabel] particle label, points back to MC tracks
  UInt_t        fFilterMap;      // filter information, one bit per set of cuts
  
  Char_t        fType;           // cluster type

  ClassDef(AliAODCluster,4);
};

#endif
