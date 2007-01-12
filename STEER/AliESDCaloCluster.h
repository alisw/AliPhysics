#ifndef ALIESDCALOCLUSTER_H
#define ALIESDCALOCLUSTER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
/* $Log $ */

//-------------------------------------------------------------------------
//                          Class AliESDCaloCluster
//   This is the class to deal with during the physics analysis of data
//
//   New container for calorimeter clusters, which are the effective 
//   "tracks" for calorimeter detectors.  Can be used by PHOS and EMCAL
//
//     J.L. Klay (LLNL)
//-------------------------------------------------------------------------

#include <TObject.h>
#include "AliPID.h"

class TLorentzVector;

class AliESDCaloCluster : public TObject {

public:

  AliESDCaloCluster();
  AliESDCaloCluster(const AliESDCaloCluster& clus);
  AliESDCaloCluster & operator=(const AliESDCaloCluster& source);
  virtual ~AliESDCaloCluster();

  void SetID(Int_t id) {fID = id;}
  Int_t GetID() const {return fID;}

  enum ClusterType {kPseudoCluster, kClusterv1};//Two types of clusters stored
                                                //in EMCAL.
  void SetClusterType(Int_t type) { fClusterType = type; }
  Int_t GetClusterType() const {return fClusterType; }

  void SetEMCAL(Bool_t emc) { fEMCALCluster = emc;}
  Bool_t IsEMCAL() const {return fEMCALCluster;}

  void SetPHOS(Bool_t phos) { fPHOSCluster = phos;}
  Bool_t IsPHOS() const {return fPHOSCluster;}

  void SetGlobalPosition(const Float_t *pos) {
    fGlobalPos[0] = pos[0]; fGlobalPos[1] = pos[1]; fGlobalPos[2] = pos[2];
  }
  void GetGlobalPosition(Float_t *pos) const {
    pos[0] = fGlobalPos[0]; pos[1] = fGlobalPos[1]; pos[2] = fGlobalPos[2];
  }

  void SetClusterEnergy(Float_t ene) { fEnergy = ene;}
  Float_t GetClusterEnergy() const   { return fEnergy;}

  void SetClusterDisp(Float_t disp)  { fDispersion = disp; }
  Float_t GetClusterDisp() const     { return fDispersion; }

  void SetClusterChi2(Float_t chi2)  { fChi2 = chi2; }
  Float_t GetClusterChi2() const     { return fChi2; }

  void SetPid(const Float_t *p);
  Float_t *GetPid() {return fPID;}

  void SetPrimaryIndex(Int_t primary)     { fPrimaryIndex = primary; }
  Int_t GetPrimaryIndex() const           { return fPrimaryIndex; }

  void SetM20(Float_t m20)                { fM20 = m20; }
  Float_t GetM20() const                  { return fM20; }

  void SetM02(Float_t m02)                { fM02 = m02; }
  Float_t GetM02() const                  { return fM02; }

  void SetM11(Float_t m11)                { fM11 = m11; }
  Float_t GetM11() const                  { return fM11; }

  void SetNExMax(UShort_t nExMax)         { fNExMax = nExMax; }
  UShort_t GetNExMax() const              { return fNExMax; }

  void SetEmcCpvDistance(Float_t dEmcCpv) { fEmcCpvDistance = dEmcCpv; }
  Float_t GetEmcCpvDistance() const       { return fEmcCpvDistance; }

  void SetNumberOfDigits(Int_t ndig)      { fNumberOfDigits = ndig; }
  Int_t GetNumberOfDigits() const         { return fNumberOfDigits; }
  
  void SetDigitAmplitude(UShort_t *adc)   { fDigitAmplitude = adc;}
  UShort_t *GetDigitAmplitude() const     { return fDigitAmplitude;}

  void SetDigitTime(UShort_t *time)       { fDigitTime = time;}
  UShort_t *GetDigitTime() const          { return fDigitTime;}

  void SetDigitIndex(UShort_t *digit)     { fDigitIndex = digit;}
  UShort_t *GetDigitIndex() const         { return fDigitIndex; }

  void GetMomentum(TLorentzVector& p);

protected:

  Int_t     fID;               // Unique Id of the cluster
  Int_t     fClusterType;      // Flag for different clustering versions
  Bool_t    fEMCALCluster;     // Is this is an EMCAL cluster?
  Bool_t    fPHOSCluster;      // Is this is a PHOS cluster?
  Float_t   fGlobalPos[3];     // position in global coordinate system
  Float_t   fEnergy;           // energy measured by calorimeter
  Float_t   fDispersion;       // cluster dispersion, for shape analysis
  Float_t   fChi2;             // chi2 of cluster fit
  Float_t   fPID[AliPID::kSPECIESN]; //"detector response probabilities" (for the PID)
  Int_t     fPrimaryIndex;     // primary track number associated with this cluster
  Float_t   fM20;              // 2-nd moment along the main eigen axis
  Float_t   fM02;              // 2-nd moment along the second eigen axis
  Float_t   fM11;              // 2-nd mixed moment Mxy
  UShort_t  fNExMax ;          // number of (Ex-)maxima before unfolding
  Float_t   fEmcCpvDistance;   // the distance from PHOS EMC rec.point to the closest CPV rec.point



  Int_t     fNumberOfDigits;   // number of calorimeter digits in cluster
                               // Very important! The streamer needs to
                               // know how big these arrays are for
                               // each event that is written out: 
  UShort_t* fDigitAmplitude;   //[fNumberOfDigits] digit energy (integer units)
  UShort_t* fDigitTime;        //[fNumberOfDigits] time of this digit (integer units)
  UShort_t* fDigitIndex;       //[fNumberOfDigits] calorimeter digit index

  ClassDef(AliESDCaloCluster,1)  //ESDCaloCluster 
};

#endif 

