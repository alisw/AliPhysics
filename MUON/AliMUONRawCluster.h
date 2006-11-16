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


#include <TObject.h>
#include <TMath.h> // because of inline funtion GetRadius
#include <TArrayF.h>
class AliMUONRawCluster : public TObject {

public:
   AliMUONRawCluster();
   virtual ~AliMUONRawCluster() { }
   Float_t      GetRadius(Int_t i) {return TMath::Sqrt(fX[i]*fX[i]+fY[i]*fY[i]);}
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

   Float_t        GetCharge(Int_t i=0) const;
   Float_t      GetX(Int_t i=0) const;
   Float_t      GetY(Int_t i=0) const;
   Float_t      GetZ(Int_t i=0) const;
   Int_t        GetTrack(Int_t i=0) const;
   Float_t        GetPeakSignal(Int_t i=0) const;
   Int_t        GetMultiplicity(Int_t i=0) const;
   Int_t        GetClusterType() const;
   Int_t        GetGhost() const;
   Int_t        GetNcluster(Int_t i=0) const;
   Float_t      GetChi2(Int_t i=0) const;
   Int_t        GetIndex(Int_t i, Int_t j) const;
   Int_t        GetOffset(Int_t i, Int_t j) const;
   Float_t      GetContrib(Int_t i, Int_t j) const;
   Int_t        GetPhysics(Int_t i) const;
   Int_t        GetDetElemId() const ; 
   Float_t      GetError(Int_t iXY) const;
   Float_t      GetErrX() const;
   Float_t      GetErrY() const;

   Int_t        SetCharge(Int_t i, Float_t Q);
   Int_t        SetX(Int_t i, Float_t X);
   Int_t        SetY(Int_t i, Float_t Y);
   Int_t        SetZ(Int_t i, Float_t Z);
   void         SetDetElemId(Int_t Id); 
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
   void         SetError(Int_t iXY, Float_t err);
   void         SetErrX(Float_t err);
   void         SetErrY(Float_t err);

private:
   Int_t       fIndexMap[50][2];  ///< Indices of digits
   Int_t       fOffsetMap[50][2]; ///< Emmanuel special
   Float_t     fContMap[50][2];   ///< Contribution from digit
   Int_t     fPhysicsMap[50];   ///< Distinguish signal and background contr.
  
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
   ClassDef(AliMUONRawCluster,2)  //Cluster class for MUON
};

// inline functions

inline  Int_t  AliMUONRawCluster::GetIndex(Int_t i, Int_t j) const
{ return fIndexMap[i][j]; }

inline  Int_t  AliMUONRawCluster::GetOffset(Int_t i, Int_t j) const
{ return fOffsetMap[i][j]; }

inline  Float_t  AliMUONRawCluster::GetContrib(Int_t i, Int_t j) const
{ return fContMap[i][j]; }

inline  Int_t  AliMUONRawCluster::GetPhysics(Int_t i) const
{ return fPhysicsMap[i]; }

inline  void  AliMUONRawCluster::SetIndex(Int_t i, Int_t j, Int_t index)
{ fIndexMap[i][j] = index; }

inline  void  AliMUONRawCluster::SetOffset(Int_t i, Int_t j, Int_t offset)
{ fOffsetMap[i][j] = offset; }

inline  void  AliMUONRawCluster::SetContrib(Int_t i, Int_t j, Float_t contrib)
{ fContMap[i][j] = contrib; }

inline  void  AliMUONRawCluster::SetPhysics(Int_t i, Int_t physics)
{ fPhysicsMap[i] = physics; }

inline void AliMUONRawCluster::SetDetElemId(Int_t Id)
{ fDetElemId = Id; }

inline Int_t AliMUONRawCluster::GetDetElemId() const
{ return fDetElemId;}

inline void AliMUONRawCluster::SetError(Int_t iXY, Float_t err)
{ fErrXY[iXY] = err; }

inline void AliMUONRawCluster::SetErrX(Float_t err)
{ SetError(0, err); }

inline void AliMUONRawCluster::SetErrY(Float_t err)
{ SetError(1, err); }

inline Float_t AliMUONRawCluster::GetError(Int_t iXY) const
{ return fErrXY[iXY]; }

inline Float_t AliMUONRawCluster::GetErrX() const
{ return GetError(0); }

inline Float_t AliMUONRawCluster::GetErrY() const
{ return GetError(1); }

#endif






