#ifndef ALIITSGEOM_H
#define ALIITSGEOM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////////////////////////
//  ITS geometry manipulation routines.
//  Created April 15 1999.
//  version: 0.0.0
//  By: Bjorn S. Nilsen
//
//     A package of geometry routines to do transformations between
// local, detector active area, and ALICE global coordinate system in such
// a way as to allow for detector alignment studies and the like. All of
// the information needed to do the coordinate transformation are kept in
// a specialized structure for ease of implementation.
/////////////////////////////////////////////////////////////////////////
#include <fstream.h>
#include <TObjArray.h>
#include <TVector.h>

#include "AliITSgeomSPD.h"
#include "AliITSgeomSDD.h"
#include "AliITSgeomSSD.h"


struct AliITSgeomS {
    Int_t   fShapeIndex; // Shape index for this volume
    Float_t fx0,fy0,fz0; // Translation vector
    Float_t frx,fry,frz; // Rotation about axis, angle radians
    Float_t fr[9];       // the rotation matrix
    Float_t angles[6];   // module center, theta and phi
    Double_t rottrack[3][3]; // the tracking rotation matrix
};

//_______________________________________________________________________

class AliITSgeom : public TObject {

 public:
    AliITSgeom();                            // Default constructor
    AliITSgeom(const char *filename);        // Constructor
    AliITSgeom(const AliITSgeom &source);    // Copy constructor
    //    void operator=(const AliITSgeom &source);// = operator
    AliITSgeom& operator=(const AliITSgeom &source);// = operator

    virtual ~AliITSgeom();                   // Default destructor

    Int_t GetNdetectors(Int_t layer) const {return fNdet[layer-1];}// return number of detector a ladder has
    Int_t GetNladders(Int_t layer)   const {return fNlad[layer-1];}// return number of laders a layer has
    Int_t GetNlayers()               const {return fNlayers;} // return number of layer the ITS has
    void GetAngles(Int_t lay,Int_t lad,Int_t det,
			  Float_t &rx,Float_t &ry,Float_t &rz)const {
                          rx = fGm[lay-1][fNdet[lay-1]*(lad-1)+det-1].frx;
                          ry = fGm[lay-1][fNdet[lay-1]*(lad-1)+det-1].fry;
                          rz = fGm[lay-1][fNdet[lay-1]*(lad-1)+det-1].frz;} // Get agnles of roations for a give module
    void GetTrans(Int_t lay,Int_t lad,Int_t det,
			 Float_t &x,Float_t &y,Float_t &z)const {
                         x = fGm[lay-1][fNdet[lay-1]*(lad-1)+det-1].fx0;
                         y = fGm[lay-1][fNdet[lay-1]*(lad-1)+det-1].fy0;
                         z = fGm[lay-1][fNdet[lay-1]*(lad-1)+det-1].fz0;} // get translation for a given module
    void SetByAngles(Int_t lay,Int_t lad,Int_t det,
		     Float_t rx,Float_t ry,Float_t rz);
    void SetByAngles(Int_t index,Double_t angl[]);
    void SetTrans(Int_t lay,Int_t lad,Int_t det,
			 Float_t x,Float_t y,Float_t z){
                         fGm[lay-1][fNdet[lay-1]*(lad-1)+det-1].fx0 = x;
                         fGm[lay-1][fNdet[lay-1]*(lad-1)+det-1].fy0 = y;
                         fGm[lay-1][fNdet[lay-1]*(lad-1)+det-1].fz0 = z;}// Set translation vector for a give module
    void SetTrans(Int_t index,Double_t x[]);
    void GetRotMatrix(Int_t lay,Int_t lad,Int_t det,Float_t *mat);
    void GetRotMatrix(Int_t lay,Int_t lad,Int_t det,Double_t *mat);
    void GetRotMatrix(Int_t index,Float_t *mat);
    void GetRotMatrix(Int_t index,Double_t *mat);
    Int_t GetStartDet(Int_t dtype );
    Int_t GetLastDet(Int_t dtype);
    Int_t GetStartSPD() {return GetModuleIndex(1,1,1);} // return starting index for SPD
    Int_t GetLastSPD() {return GetModuleIndex(2,fNlad[1],fNdet[1]);}// return Ending index for SPD
    Int_t GetStartSDD() {return GetModuleIndex(3,1,1);} // return starting index for SDD
    Int_t GetLastSDD() {return GetModuleIndex(4,fNlad[3],fNdet[3]);}// return Ending index for SDD
    Int_t GetStartSSD() {return GetModuleIndex(5,1,1);} // return starting index for SSD
    Int_t GetLastSSD() {return GetModuleIndex(6,fNlad[5],fNdet[5]);}// return Ending index for SSD
    Int_t GetIndexMax() {return GetModuleIndex(fNlayers,fNlad[fNlayers-1],
					       fNdet[fNlayers-1])+1;}// return Ending index for all ITS
    void GtoL(Int_t lay,Int_t lad,Int_t det,const Float_t *g,Float_t *l);
    void GtoL(const Int_t *id,const Float_t *g,Float_t *l);
    void GtoL(const Int_t index,const Float_t *g,Float_t *l);
    void GtoL(Int_t lay,Int_t lad,Int_t det,const Double_t *g,Double_t *l);
    void GtoL(const Int_t *id,const Double_t *g,Double_t *l);
    void GtoL(const Int_t index,const Double_t *g,Double_t *l);
    void GtoLMomentum(Int_t lay,Int_t lad,Int_t det,
		      const Float_t *g,Float_t *l);
    void GtoLMomentum(Int_t lay,Int_t lad,Int_t det,
		      const Double_t *g,Double_t *l);
    void LtoG(Int_t lay,Int_t lad,Int_t det,const Float_t *l,Float_t *g);
    void LtoG(const Int_t *id,const Float_t *l,Float_t *g);
    void LtoG(const Int_t index,const Float_t *l,Float_t *g);
    void LtoG(Int_t lay,Int_t lad,Int_t det,const Double_t *l,Double_t *g);
    void LtoG(const Int_t *id,const Double_t *l,Double_t *g);
    void LtoG(const Int_t index,const Double_t *l,Double_t *g);
    void LtoGMomentum(Int_t lay,Int_t lad,Int_t det,
		      const Float_t *l,Float_t *g);
    void LtoGMomentum(Int_t lay,Int_t lad,Int_t det,
		      const Double_t *l,Double_t *g);
    void LtoL(const Int_t *id1,const Int_t *id2,Double_t *l1,Double_t *l2);
    void LtoL(const Int_t index1,const Int_t index2,Double_t *l1,Double_t *l2);
    void LtoLMomentum(const Int_t *id1,const Int_t *id2,
		      const Double_t *l1,Double_t *l2);
    void GtoLErrorMatrix(const Int_t index,Double_t **g,Double_t **l);
    void LtoGErrorMatrix(const Int_t index,Double_t **l,Double_t **g);
    void LtoLErrorMatrix(const Int_t index1,const Int_t index2,
			 Double_t **l1,Double_t **l2);
    Int_t GetModuleIndex(Int_t lay,Int_t lad,Int_t det);
    void GetModuleId(Int_t index,Int_t &lay,Int_t &lad,Int_t &det);
    void GlobalChange(Float_t  *tran,Float_t  *rot);
    void GlobalCylindericalChange(Float_t *tran,Float_t *rot);
    void RandomChange(Float_t *stran,Float_t *srot);
    void RandomCylindericalChange(Float_t *stran,Float_t *srot);
    void PrintComparison(FILE *fp,AliITSgeom *other);
    void PrintData(FILE *fp,Int_t lay,Int_t lad,Int_t det);
    ofstream &PrintGeom(ofstream &out);
    ifstream &ReadGeom(ifstream &in);
    virtual Int_t IsVersion() const {return 1;} // return version number
    void AddShape(TObject *shp){fShape->AddLast(shp);} // Add Shapes
    virtual TObject *GetShape(Int_t lay,Int_t lad,Int_t det)
	const {return fShape->At(fGm[lay-1][fNdet[lay-1]*(lad-1)+det-1].
				 fShapeIndex);} // return specific TShape

    TObjArray *Shape() {return fShape;} // return Shapes array

    void GeantToTracking(AliITSgeom &source); // This converts the geometry
    // transformations from that used by the ITS and it's Monte Carlo to that
    // used by the track finding code.
    // Usage:
    // AliITSgeom *gm,*gt;
    // gm = ((AliITS *) ITS)->GetITSgeom();
    // gt->GeantToTracking(*gm);
    // This allocates and fills gt with the geometry transforms between the
    // global coordinate system to the local coordinate system used to do
    // tracking.

    void GtoLtracking(Int_t lay,Int_t lad,Int_t det,const Double_t *g,Double_t *l);
    void LtoGtracking(Int_t lay,Int_t lad,Int_t det,const Double_t *l,Double_t *g);
    void GetCenterThetaPhi(Int_t lay,Int_t lad,Int_t det, TVector &x) const {
                       x(0) = fGm[lay-1][fNdet[lay-1]*(lad-1)+det-1].fx0;
		       x(1) = fGm[lay-1][fNdet[lay-1]*(lad-1)+det-1].fy0;
		       x(2) = fGm[lay-1][fNdet[lay-1]*(lad-1)+det-1].fz0;	                         		 
		       x(3) = fGm[lay-1][fNdet[lay-1]*(lad-1)+det-1].angles[0];
		       x(4) = fGm[lay-1][fNdet[lay-1]*(lad-1)+det-1].angles[1];
		       x(5) = fGm[lay-1][fNdet[lay-1]*(lad-1)+det-1].angles[2];
		       x(6) = fGm[lay-1][fNdet[lay-1]*(lad-1)+det-1].angles[3];
		       x(7) = fGm[lay-1][fNdet[lay-1]*(lad-1)+det-1].angles[4];
		       x(8) = fGm[lay-1][fNdet[lay-1]*(lad-1)+det-1].angles[5];

}

 private:
    Int_t        fNlayers; // The number of layers.
    Int_t        *fNlad;   // Array of the number of ladders/layer(layer)
    Int_t        *fNdet;   // Array of the number of detectors/ladder(layer)
    AliITSgeomS  **fGm;     // Structure of translation and rotation.
    TObjArray    *fShape;  // Array of shapes and detector information.
    
    ClassDef(AliITSgeom,1) // ITS geometry class
};

#endif
