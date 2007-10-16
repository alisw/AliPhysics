#ifndef AliAODCaloCluster_H
#define AliAODCaloCluster_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//     AOD calorimeter cluster class (for PHOS and EMCAL)
//     Author: Markus Oldenburg, CERN, 
//             Gustavo Conesa, INFN
//-------------------------------------------------------------------------

#include "AliAODCluster.h"

#include <TRefArray.h>
#include <TArrayS.h>

class AliAODCaloCluster : public AliAODCluster {

 public:
  
  AliAODCaloCluster();
  AliAODCaloCluster(Int_t id,
		    UInt_t nLabel,
		    Int_t *label,
		    Double_t energy,
		    Double_t x[3],
		    Double_t pid[9],
		    Char_t ttype=kUndef,
		    UInt_t selectInfo=0);
  
  AliAODCaloCluster(Int_t id,
		    UInt_t nLabel,
		    Int_t *label,
		    Float_t energy,
		    Float_t x[3],
		    Float_t pid[9],
		    Char_t ttype=kUndef,
		    UInt_t selectInfo=0);
  
  virtual ~AliAODCaloCluster();
  AliAODCaloCluster(const AliAODCaloCluster& clus); 
  AliAODCaloCluster& operator=(const AliAODCaloCluster& clus);

  
  // getters
  Double_t GetDistToBadChannel() const { return fDistToBadChannel; }
  Double_t GetDispersion() const { return fDispersion; }
  Double_t GetM20() const { return fM20; }
  Double_t GetM01() const { return fM01; }
  Double_t GetM02() const { return fM02; }
  Double_t GetM11() const { return fM11; }
  Double_t GetEmcCpvDistance() const { return fEmcCpvDistance; }
  UShort_t GetNExMax() const { return fNExMax; }

  Int_t    GetNTracksMatched() const { return fTracksMatched.GetEntriesFast(); }
  TObject *GetTrackMatched(Int_t i) const { return fTracksMatched.At(i); }
  Int_t    GetNCellNumbers() const { return fCellNumber.GetSize(); }
  UShort_t GetCellNumber(Int_t i) const { return fCellNumber.At(i); }

  // setters
  void SetDistToBadChannel(Double_t dist) { fDistToBadChannel = dist; }
  void SetDispersion(Double_t disp) { fDispersion = disp; }
  void SetM20(Double_t m20) { fM20 = m20; }
  void SetM01(Double_t m01) { fM01 = m01; }
  void SetM02(Double_t m02) { fM02 = m02; }
  void SetM11(Double_t m11) { fM11 = m11; }
  void SetEmcCpvDistance(Double_t emcCpvDist) { fEmcCpvDistance = emcCpvDist; }
  void SetNExMax(UShort_t nExMax) { fNExMax = nExMax; }

  void SetCaloCluster(Double_t dist = -999., 
		      Double_t disp = -1., 
		      Double_t m20 = 0., 
		      Double_t m01 = 0., 
		      Double_t m02 = 0., 
		      Double_t m11 = 0., 
		      Double_t emcCpvDist = -999., 
		      UShort_t nExMax = 0) 
  {
    fDistToBadChannel = dist;
    fDispersion = disp;
    fM20 = m20;
    fM01 = m01;
    fM02 = m02;
    fM11 = m11;
    fEmcCpvDistance = emcCpvDist;
    fNExMax = nExMax;
  }

  void AddTrackMatched(TObject *trk) { fTracksMatched.Add(trk); }
  void RemoveTrackMatched(TObject *trk) { fTracksMatched.Remove(trk); }
  Bool_t HasTrackMatched(TObject *trk) const;

 private :

  Double32_t   fDistToBadChannel; // Distance to nearest bad channel
  Double32_t   fDispersion;       // cluster dispersion, for shape analysis
  Double32_t   fM20;              // 2-nd moment along the main eigen axis
  Double32_t   fM01;              // 
  Double32_t   fM02;              // 2-nd moment along the second eigen axis
  Double32_t   fM11;              // 2-nd mixed moment Mxy
  Double32_t   fEmcCpvDistance;   // the distance from PHOS EMC rec.point to the closest CPV rec.point
  UShort_t     fNExMax;           // number of (Ex-)maxima before unfolding

  TRefArray    fTracksMatched;    // references to tracks close to cluster. First entry is the most likely match.
  TArrayS      fCellNumber;       // fired calorimeter cell numbers

  ClassDef(AliAODCaloCluster,1);
};

#endif
