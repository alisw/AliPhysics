#ifndef ALIITSGEOMETRYITSV_H
#define ALIITSGEOMETRYITSV_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
  $Id$
 */

/*
  ITS Mother volume for ITS Geometry v11.
*/
#include "AliITSBaseGeometry.h"
class TVector3;
class AliITS;
 
class  AliITSGeometryITSV : public AliITSBaseGeometry {
 public:
    AliITSGeometryITSV();
    AliITSGeometryITSV(AliITS *its,const char *moth);
    virtual ~AliITSGeometryITSV(){};
    void CreateG3Geometry();
    void PositionGeometry(const char *moth,Int_t copy,TVector3 &t,
			  Int_t irot=0);
    void CreateG3Materials();
    void BuildDisplayGeometry();
    AliITSBaseVolParams& GetParams(){return fA;}// Returns parameters of
    // this logical volume.
    void PolyCone(AliITSPConeData &d,Int_t med); // Special version for this
 private:
    AliITSPConeData fA; // Poly-cone Volume A.
    Int_t fAir; // ITS mother volume material number.

    ClassDef(AliITSGeometryITSV,1)// ITS Mother Volume for Geometry v11.
};
 
#endif

