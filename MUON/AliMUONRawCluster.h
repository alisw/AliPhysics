#ifndef ALIMUONRAWCLUSTER_H
#define ALIMUONRAWCLUSTER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

/// \ingroup base
/// \class AliMUONRawCluster
/// \brief MUON raw cluster
///
/// Class for the MUON RecPoint
/// It contains the properties of the physics cluters found in the tracking chambers
/// RawCluster contains also the information from the both cathode of the chambers.


#include "AliMUONVCluster.h"
#include <TMath.h> // because of inline funtion GetRadius
#include <TArrayF.h>

class AliMUONRawCluster : public AliMUONVCluster {

public:
   AliMUONRawCluster();
   virtual ~AliMUONRawCluster();
   
           /// Clear method (used by TClonesArray)
   virtual void Clear(Option_t* = "") {}
   
           /// Create a copy of the current cluster
   virtual AliMUONRawCluster* Clone(const char* = "") const {return new AliMUONRawCluster(*this);}
   
           /// Set coordinates (cm)
  virtual void     SetXYZ(Double_t x, Double_t y, Double_t z) {fX[0] = x; fY[0] = y; fZ[0] = z;}
           /// Return coordinate X (cm)
  virtual Double_t GetX() const {return fX[0];}
           /// Return coordinate Y (cm)
  virtual Double_t GetY() const {return fY[0];}
           /// Return coordinate Z (cm)
  virtual Double_t GetZ() const {return fZ[0];}
  
	   /// Set resolution (cm) on coordinates (X,Y)
  virtual void     SetErrXY(Double_t errX, Double_t errY) {fErrXY[0] = errX; fErrXY[1] = errY;}
           /// Return resolution (cm) on coordinate X
  virtual Double_t GetErrX() const {return fErrXY[0];}
           /// Return resolution**2 (cm**2) on coordinate X
  virtual Double_t GetErrX2() const {return fErrXY[0] * fErrXY[0];}
           /// Return resolution (cm) on coordinate Y
  virtual Double_t GetErrY() const {return fErrXY[1];}
           /// Return resolution**2 (cm**2) on coordinate Y
  virtual Double_t GetErrY2() const {return fErrXY[1] * fErrXY[1];}
  
           /// Set the cluster charge
  virtual void     SetCharge(Double_t q) {fQ[0] = q;}
           /// Set the cluster charge
  virtual Double_t GetCharge() const {return fQ[0];}
  
           /// Return chamber Id
  virtual Int_t    GetChamberId() const {return fDetElemId/100 - 1;}
           /// Set detection element Id
          void     SetDetElemId(Int_t id) {fDetElemId = id;}
           /// Return detection element Id
  virtual Int_t    GetDetElemId() const {return fDetElemId;}
  
  virtual void     SetDigitsId(Int_t nDigits, const UInt_t *digitsId);
           /// Add a digit Id to the array of associated digits
  virtual void     AddDigitId(UInt_t id) {fIndexMap[fMultiplicity[0]++][0] = id;}

           /// Return number of associated digits
  virtual Int_t    GetNDigits() const {return fMultiplicity[0];}
           /// Return Id of digits i
  virtual UInt_t   GetDigitId(Int_t i) const {return (i < fMultiplicity[0] && i < 50) ? (UInt_t)fIndexMap[i][0] : 0;}
  
           /// Set chi2 of cluster
  virtual void     SetChi2( Double_t chi2) {fChi2[0] = chi2;}
           /// Return chi2 of cluster
  virtual Double_t GetChi2() const {return fChi2[0];}
   
           /// Set the corresponding MC track number
  virtual void     SetMCLabel(Int_t label) {SetTrack(0, label);}
           /// Return the corresponding MC track number
  virtual Int_t    GetMCLabel() const {return GetTrack(0);}
  
  /// Return radius
   Float_t      GetRadius(Int_t i) const {return TMath::Sqrt(fX[i]*fX[i]+fY[i]*fY[i]);}
   /// Return true as the function Compare() is implemented
   Bool_t       IsSortable() const {return kTRUE;}
   Int_t        Compare(const TObject *obj) const;
   Int_t        PhysicsContribution() const;
   virtual void Print(Option_t* opt="") const;
   static Int_t BinarySearch(Float_t r, TArrayF ccord, Int_t from, Int_t upto);
   static void  SortMin(Int_t *idx,Float_t *xdarray, Float_t *xarray, Float_t *yarray, Float_t *qarray,Int_t ntr);
   void         DumpIndex();

   Int_t        AddCharge(Int_t i, Float_t Q);
   Int_t        AddX(Int_t i, Float_t X);
   Int_t        AddY(Int_t i, Float_t Y);
   Int_t        AddZ(Int_t i, Float_t Z);

   Float_t      GetCharge(Int_t i) const;
   Float_t      GetX(Int_t i) const;
   Float_t      GetY(Int_t i) const;
   Float_t      GetZ(Int_t i) const;
   Int_t        GetTrack(Int_t i=0) const;
   Float_t      GetPeakSignal(Int_t i=0) const;
   Int_t        GetMultiplicity(Int_t i=0) const;
   Int_t        GetClusterType() const;
   Int_t        GetGhost() const;
   Int_t        GetNcluster(Int_t i=0) const;
   Float_t      GetChi2(Int_t i) const;
   Int_t        GetIndex(Int_t i, Int_t j) const;
   Int_t        GetOffset(Int_t i, Int_t j) const;
   Float_t      GetContrib(Int_t i, Int_t j) const;
   Int_t        GetPhysics(Int_t i) const;

   Int_t        SetCharge(Int_t i, Float_t Q);
   Int_t        SetX(Int_t i, Float_t X);
   Int_t        SetY(Int_t i, Float_t Y);
   Int_t        SetZ(Int_t i, Float_t Z);
   Int_t        SetTrack(Int_t i, Int_t track);
   Int_t        SetPeakSignal(Int_t i, Float_t peaksignal);
   Int_t        SetMultiplicity(Int_t i, Int_t mul);
   Int_t        SetClusterType(Int_t type);
   Int_t        SetGhost(Int_t ghost);
   Int_t        SetNcluster(Int_t i, Int_t ncluster);
   Int_t        SetChi2(Int_t i, Float_t chi2);
   void         SetIndex(Int_t i, Int_t j, Int_t index);
   void         SetOffset(Int_t i, Int_t j, Int_t offset);
   void         SetContrib(Int_t i, Int_t j, Float_t contrib);
   void         SetPhysics(Int_t i, Int_t physics);

private:
   Int_t       fIndexMap[50][2];  ///< Indices of digits
   Int_t       fOffsetMap[50][2]; ///< Emmanuel special
   Float_t     fContMap[50][2];   ///< Contribution from digit
   Int_t       fPhysicsMap[50];   ///< Distinguish signal and background contr.
  
   Float_t     fQ[2]  ;           ///< Q of cluster (in ADC counts)     
   Float_t     fX[2]  ;           ///< X of cluster
   Float_t     fY[2]  ;           ///< Y of cluster
   Float_t     fZ[2]  ;           ///< Z of cluster
   Int_t       fTracks[3];        ///< Labels of overlapped tracks
   Float_t     fPeakSignal[2];    ///< Peak signal 
   Int_t       fMultiplicity[2];  ///< Cluster multiplicity
   Int_t       fClusterType;      ///< Cluster type
   Int_t       fGhost;            ///< Ghost info
                                  // 0 if not a ghost or ghost problem solved
                                  // >0 if ghost problem remains because
                                  // 1 both (true and ghost) satify 
                                  //   charge chi2 compatibility
                                  // 2 none give satisfactory chi2
   Int_t       fNcluster[2];      ///< Number of clusters
   Float_t     fChi2[2];          ///< Chi**2 of fit
   Int_t       fDetElemId;        ///< ID number of the detection element (slat) on which the cluster is found. 
   Float_t     fErrXY[2];         ///< coordinate errors
   
   ClassDef(AliMUONRawCluster,3)  //Cluster class for MUON
};

// inline functions

/// Return Indices of digits
inline  Int_t  AliMUONRawCluster::GetIndex(Int_t i, Int_t j) const
{ return fIndexMap[i][j]; }

/// Return Emmanuel special offset map
inline  Int_t  AliMUONRawCluster::GetOffset(Int_t i, Int_t j) const
{ return fOffsetMap[i][j]; }

/// Return Contribution from digit
inline  Float_t  AliMUONRawCluster::GetContrib(Int_t i, Int_t j) const
{ return fContMap[i][j]; }

/// Return Distinguish signal and background contr.
inline  Int_t  AliMUONRawCluster::GetPhysics(Int_t i) const
{ return fPhysicsMap[i]; }

/// Set Indices of digits
inline  void  AliMUONRawCluster::SetIndex(Int_t i, Int_t j, Int_t index)
{ fIndexMap[i][j] = index; }

/// Set Emmanuel special offset map
inline  void  AliMUONRawCluster::SetOffset(Int_t i, Int_t j, Int_t offset)
{ fOffsetMap[i][j] = offset; }

/// Set Contribution from digit
inline  void  AliMUONRawCluster::SetContrib(Int_t i, Int_t j, Float_t contrib)
{ fContMap[i][j] = contrib; }

/// Set Distinguish signal and background contr.
inline  void  AliMUONRawCluster::SetPhysics(Int_t i, Int_t physics)
{ fPhysicsMap[i] = physics; }


#endif

