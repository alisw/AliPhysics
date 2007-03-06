#ifndef ALICLUSTERTGEO_H
#define ALICLUSTERTGEO_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//                         Class AliClusterTGeo
// This is the future base for managing the clusters in barrel detectors.
// It is fully interfaced with the ROOT geometrical modeller TGeo.
// Each cluster contains XYZ coordinates in the local tracking c.s. and
// the unique ID of the sensitive detector element which continas the
// cluster. The coordinates in global c.s. are computed using the interface
// to TGeo and will be not overwritten by the derived sub-detector cluster
// classes.
//
// cvetan.cheshkov@cern.ch  & jouri.belikov@cern.ch     5/3/2007
//-------------------------------------------------------------------------

#include <TObject.h>

class TGeoHMatrix;
class TGeoPNEntry;

class AliClusterTGeo : public TObject {
 public:
  AliClusterTGeo();
  AliClusterTGeo(UShort_t volId, const Float_t *hit, Float_t x = 0, Float_t sigyz = 0, const Int_t *lab = NULL);
  virtual ~AliClusterTGeo() {;}

  Int_t    GetLabel(Int_t i) const {return fTracks[i];}
  Float_t  GetX()            const {return fX;}
  Float_t  GetY()            const {return fY;}
  Float_t  GetZ()            const {return fZ;}
  Float_t  GetSigmaY2()      const {return fSigmaY2;}
  Float_t  GetSigmaZ2()      const {return fSigmaZ2;}
  UShort_t GetVolumeId()     const {return fVolumeId;}

  virtual void Use(Int_t = 0) {;}

  Bool_t   GetGlobalXYZ(Float_t xyz[3]) const;
  Bool_t   GetGlobalCov(Float_t cov[6]) const;
  Bool_t   GetXRefPlane(Float_t &xref) const;

  Bool_t   Misalign();

 protected:

  AliClusterTGeo(const AliClusterTGeo& cluster);
  AliClusterTGeo &operator=(const AliClusterTGeo& cluster);

 private:

  TGeoHMatrix*         GetMatrix(Bool_t original = kFALSE) const;
  const TGeoHMatrix*   GetTracking2LocalMatrix() const;
  TGeoPNEntry*         GetPNEntry() const;

  Int_t    fTracks[3];//MC labels
  Float_t  fX;        // X of the cluster in the tracking c.s.
  Float_t  fY;        // Y of the cluster in the tracking c.s.
  Float_t  fZ;        // Z of the cluster in the tracking c.s.
  Float_t  fSigmaY2;  // Sigma Y square of cluster
  Float_t  fSigmaZ2;  // Sigma Z square of cluster
  Float_t  fSigmaYZ;  // Non-diagonal element of cov.matrix
  UShort_t fVolumeId; // Volume ID of the detector element
  Bool_t   fIsMisaligned; // Cluster was misagned or not?

  ClassDef(AliClusterTGeo,1) // Barrel detectors cluster
};

#endif
