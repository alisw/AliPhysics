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
#include <Riostream.h>
#include <TObject.h>
#include <TObjArray.h>
#include <TVector.h>

#include "AliITSgeomMatrix.h"

typedef enum {kSPD=0, kSDD=1, kSSD=2, kSSDp=3,kSDDp=4} AliITSDetector;

//_______________________________________________________________________

class AliITSgeom : public TObject {

 public:
    AliITSgeom();                      // Default constructor
    AliITSgeom(const char *filename);  // Constructor
    AliITSgeom(Int_t itype,Int_t nlayers,Int_t *nlads,Int_t *ndets,
	       Int_t nmods); // Constructor
    AliITSgeom(const AliITSgeom &source);    // Copy constructor
    AliITSgeom& operator=(const AliITSgeom &source);// = operator
    virtual ~AliITSgeom();             // Default destructor
    // this function allocates a AliITSgeomMatrix for a particular module.
    void CreatMatrix(Int_t mod,Int_t lay,Int_t lad,Int_t det,
		     AliITSDetector idet,const Double_t tran[3],
		     const Double_t rot[10]);
    void ReadNewFile(const char *filename);  // Constructor for new format.
    void WriteNewFile(const char *filename); // Output for new format.
// Getters
    Int_t GetTransformationType() const {return fTrans;}
//
    // returns kTRUE if the transformation defined by this class is
    // for Global GEANT coordinate system to the local GEANT coordinate system
    // of the detector. These are the transformation used by GEANT.
    Bool_t IsGeantToGeant()     const {return (fTrans == 0);}
    // returns kTRUE if the transformation defined by this class is
    // for Global GEANT coordinate system to the local "Tracking" coordinate
    // system of the detector. These are the transformation used by the
    // Tracking code.
    Bool_t IsGeantToTracking()  const {return ((fTrans&&0xfffe)!= 0);}
    // returns kTRUE if the transformation defined by this class is
    // for Global GEANT coordinate system to the local GEANT coordinate system
    // of the detector but may have been displaced by some typically small
    // amount. These are modified transformation similar to that used by GEANT.
    Bool_t IsGeantToDisplaced() const {return ((fTrans&&0xfffd)!= 0);}
    // returns kTRUE if the shape defined by ishape has been defined in this
    // set of transformations. Typical values of ishape are kSPD, kSDD, kSSD,
    // SSD2.
    Bool_t IsShapeDefined(Int_t ishape)const {
	if(fShape!=0){return ((fShape->At(ishape))!=0);}else return kFALSE;}
//
    //     This function returns a pointer to the particular AliITSgeomMatrix
    // class for a specific module index.
    AliITSgeomMatrix *GetGeomMatrix(Int_t index){
	return (AliITSgeomMatrix*)(fGm->At(index));}
    //     This function returns the number of detectors/ladder for a give 
    // layer. In particular it returns fNdet[layer-1].
    Int_t GetNdetectors(Int_t lay) const {return fNdet[lay-1];}
    //     This function returns the number of ladders for a give layer. In
    // particular it returns fNlad[layer-1].
    Int_t GetNladders(Int_t lay)   const {return fNlad[lay-1];};
    //     This function returns the number of layers defined in the ITS
    // geometry. In particular it returns fNlayers.
    Int_t GetNlayers()                   const {return fNlayers;}
    Int_t GetModuleIndex(Int_t lay,Int_t lad,Int_t det);
    //     This function returns the module index number given the layer,
    // ladder and detector numbers put into the array id[3].
    Int_t GetModuleIndex(const Int_t *id){
	return GetModuleIndex(id[0],id[1],id[2]);}
    void  GetModuleId(Int_t index,Int_t &lay,Int_t &lad,Int_t &det);
    // Returns the detector type
    Int_t GetModuleType(Int_t index){
	                   return GetGeomMatrix(index)->GetDetectorIndex();}
    // Returns the detector type as a string
    Char_t * GetModuleTypeName(Int_t index){switch(GetModuleType(index)) {
    case kSPD: return "kSPD";case kSDD: return "kSDD";case kSSD: return "kSSD";
    case kSSDp: return"kSSDp";case kSDDp: return "kSDDp"; default: return "";}}
//
    Int_t GetStartDet(Int_t dtype );
    Int_t GetLastDet(Int_t dtype);
    //     Returns the starting module index number for SPD detector,
    // assuming the modules are placed in the "standard" cylindrical
    // ITS structure.
    Int_t GetStartSPD() {return GetModuleIndex(1,1,1);}
    //     Returns the ending module index number for SPD detector,
    // assuming the modules are placed in the "standard" cylindrical
    // ITS structure.
    Int_t GetLastSPD()  {return GetModuleIndex(2,fNlad[1],fNdet[1]);}
    //     Returns the starting module index number for SDD detector,
    // assuming the modules are placed in the "standard" cylindrical
    // ITS structure.
    Int_t GetStartSDD() {return GetModuleIndex(3,1,1);}
    //     Returns the ending module index number for SDD detector,
    // assuming the modules are placed in the "standard" cylindrical
    // ITS structure.
    Int_t GetLastSDD()  {return GetModuleIndex(4,fNlad[3],fNdet[3]);}
    //     Returns the starting module index number for SSD detector,
    // assuming the modules are placed in the "standard" cylindrical
    // ITS structure.
    Int_t GetStartSSD() {return GetModuleIndex(5,1,1);}
    //     Returns the ending module index number for SSD detector,
    // assuming the modules are placed in the "standard" cylindrical
    // ITS structure.
    Int_t GetLastSSD()  {return GetModuleIndex(6,fNlad[5],fNdet[5]);}
    //     Returns the last module index number.
    Int_t GetIndexMax() const {return fNmodules;}
//
    //     This function returns the rotation angles for a give module 
    // in the Double point array ang[3]. The angles are in radians
    void  GetAngles(Int_t index,Double_t *ang) {
                    GetGeomMatrix(index)->GetAngles(ang);}
    //     This function returns the rotation angles for a give module
    // in the three floating point variables provided. rx = frx,
    // fy = fry, rz = frz. The angles are in radians
    void  GetAngles(Int_t index,Float_t &rx,Float_t &ry,Float_t &rz) {
                    Double_t a[3];GetAngles(index,a);
                    rx = a[0];ry = a[1];rz = a[2];}
    //     This function returns the rotation angles for a give detector on
    // a give ladder in a give layer in the three floating point variables
    // provided. rx = frx, fy = fry, rz = frz. The angles are in radians
    void  GetAngles(Int_t lay,Int_t lad,Int_t det,
                    Float_t &rx,Float_t &ry,Float_t &rz) {
                    GetAngles(GetModuleIndex(lay,lad,det),rx,ry,rz);}
//
    //     This function returns the 6 GEANT rotation angles for a give 
    // module in the double point array ang[3]. The angles are in degrees
    void  GetGeantAngles(Int_t index,Double_t *ang){
	GetGeomMatrix(index)->SixAnglesFromMatrix(ang);}
//
    //     This function returns the Cartesian translation for a give
    // module in the Double array t[3]. The units are
    // those of the Monte Carlo, generally cm.
    void  GetTrans(Int_t index,Double_t *t) {
                   GetGeomMatrix(index)->GetTranslation(t);}
    //     This function returns the Cartesian translation for a give
    // module index in the three floating point variables provided.
    // x = fx0, y = fy0, z = fz0. The units are those of the Mont
    // Carlo, generally cm.
    void  GetTrans(Int_t index,Float_t &x,Float_t &y,Float_t &z) {
                   Double_t t[3];GetTrans(index,t);
                   x = t[0];y = t[1];z = t[2];}
    //     This function returns the Cartesian translation for a give
    // detector on a give ladder in a give layer in the three floating
    // point variables provided. x = fx0, y = fy0, z = fz0. The units are
    // those of the Monte Carlo, generally cm.
    void  GetTrans(Int_t lay,Int_t lad,Int_t det,
                   Float_t &x,Float_t &y,Float_t &z) {
                   GetTrans(GetModuleIndex(lay,lad,det),x,y,z);}
//
    //     This function returns the Cartesian translation for a give
    // module in the Double array t[3]. The units are
    // those of the Monte Carlo, generally cm.
    void  GetTransCyln(Int_t index,Double_t *t) {
                   GetGeomMatrix(index)->GetTranslationCylinderical(t);}
    //     This function returns the Cartesian translation for a give
    // module index in the three floating point variables provided.
    // x = fx0, y = fy0, z = fz0. The units are those of the Mont
    // Carlo, generally cm.
    void  GetTransCyln(Int_t index,Float_t &x,Float_t &y,Float_t &z) {
                   Double_t t[3];GetTransCyln(index,t);
                   x = t[0];y = t[1];z = t[2];}
    //     This function returns the Cartesian translation for a give
    // detector on a give ladder in a give layer in the three floating
    // point variables provided. x = fx0, y = fy0, z = fz0. The units are
    // those of the Monte Carlo, generally cm.
    void  GetTransCyln(Int_t lay,Int_t lad,Int_t det,
                   Float_t &x,Float_t &y,Float_t &z) {
                   GetTransCyln(GetModuleIndex(lay,lad,det),x,y,z);}
//
    //      This function returns the Cartesian translation [cm] and the
    // 6 GEANT rotation angles [degrees]for a given layer ladder and
    // detector number, in the TVector x (at least 9 elements large).
    // This function is required to be inlined for speed.
    void  GetCenterThetaPhi(Int_t lay,Int_t lad,Int_t det,TVector &x){
    Double_t t[3],a[6];Int_t i=GetModuleIndex(lay,lad,det);GetTrans(i,t);
    GetGeantAngles(i,a);x(0)=t[0];x(1)=t[1];x(2)=t[2];x(3)=a[0];x(4)=a[1];
    x(5)=a[2];x(6)=a[3];x(7)=a[4];x(8)=a[5];}
//
    //     This function returns the rotation matrix in Double
    // precision for a given module.
    void  GetRotMatrix(Int_t index,Double_t mat[3][3]){
          GetGeomMatrix(index)->GetMatrix(mat);}
    //     This function returns the rotation matrix in a Double
    // precision pointer for a given module. mat[i][j] => mat[3*i+j].
    void  GetRotMatrix(Int_t index,Double_t *mat){
          Double_t rot[3][3];GetRotMatrix(index,rot);
          for(Int_t i=0;i<3;i++)for(Int_t j=0;j<3;j++) mat[3*i+j] = rot[i][j];}
    //     This function returns the rotation matrix in a floating 
    // precision pointer for a given layer ladder and detector module.
    // mat[i][j] => mat[3*i+j].
    void  GetRotMatrix(Int_t lay,Int_t lad,Int_t det,
                 Float_t *mat){GetRotMatrix(GetModuleIndex(lay,lad,det),mat);}
    //     This function returns the rotation matrix in a Double
    // precision pointer for a given layer ladder and detector module.
    // mat[i][j] => mat[3*i+j].
    void  GetRotMatrix(Int_t lay,Int_t lad,Int_t det,
                Double_t *mat){GetRotMatrix(GetModuleIndex(lay,lad,det),mat);}
    //     This function returns the rotation matrix in a floating
    // precision pointer for a given module. mat[i][j] => mat[3*i+j].
    void  GetRotMatrix(Int_t index,Float_t *mat){
          Double_t rot[3][3];
	  GetGeomMatrix(index)->GetMatrix(rot);
          for(Int_t i=0;i<3;i++)for(Int_t j=0;j<3;j++) mat[3*i+j] = rot[i][j];}
//
    //     Will define fShape if it isn't already defined.
    void DefineShapes(Int_t size=4)
	{if(fShape==0) fShape = new TObjArray(size);else fShape->Expand(size);}
    //     this function returns a pointer to the class describing a particular
    // detector type based on AliITSDetector value. This will return a pointer
    // to one of the classes AliITSgeomSPD, AliITSgeomSDD, or AliITSgeomSSD,
    // for example.
    virtual TObject *GetShape(AliITSDetector idet)
	{return fShape->At((Int_t)idet);};
    //     This function returns a pointer to the class describing the
    // detector for a particular module index. This will return a pointer
    // to one of the classes AliITSgeomSPD, AliITSgeomSDD, or AliITSgeomSSD,
    // for example.
    virtual TObject *GetShape(Int_t index){
	return fShape->At(GetGeomMatrix(index)->
			  GetDetectorIndex());}
    //     This function returns a pointer to the class describing the
    // detector for a particular layer ladder and detector numbers. This
    // will return a pointer to one of the classes AliITSgeomSPD,
    // AliITSgeomSDD, or AliITSgeomSSD, for example.
    virtual TObject *GetShape(Int_t lay,Int_t lad,Int_t det)
	                     {return GetShape(GetModuleIndex(lay,lad,det));}
//
//  Setters
    //     Sets the rotation angles and matrix for a give module index
    // via the double precision array a[3] [radians].
    void SetByAngles(Int_t index,const Double_t a[]){
	GetGeomMatrix(index)->SetAngles(a);}
    //     Sets the rotation angles and matrix for a give module index
    // via the 3 floating precision variables rx, ry, and rz [radians].
    void SetByAngles(Int_t index,
		     Float_t rx, Float_t ry, Float_t rz) {
                     Double_t a[3];a[0] = rx;a[1] = ry;a[2] = rz;
                     GetGeomMatrix(index)->SetAngles(a);}
    //     Sets the rotation angles and matrix for a give layer, ladder,
    // and detector numbers via the 3 floating precision variables rx,
    // ry, and rz [radians].
    void SetByAngles(Int_t lay,Int_t lad,Int_t det,
                     Float_t rx, Float_t ry, Float_t rz) {
                     SetByAngles(GetModuleIndex(lay,lad,det),rx,ry,rz);}
//
    //     Sets the rotation angles and matrix for a give module index
    // via the Double precision array a[6] [degree]. The angles are those
    // defined by GEANT 3.12.
    void SetByGeantAngles(Int_t index,const Double_t *ang){
	GetGeomMatrix(index)->MatrixFromSixAngles(ang);}
    //     Sets the rotation angles and matrix for a give layer, ladder
    // and detector, in the array id[3] via the Double precision array
    // a[6] [degree]. The angles are those defined by GEANT 3.12.
    void SetByGeantAngles(const Int_t *id,const Double_t *ang){
	SetByGeantAngles(GetModuleIndex(id),ang);}
    //     Sets the rotation angles and matrix for a give layer, ladder
    // and detector, via the Double precision array a[6] [degree]. The
    // angles are those defined by GEANT 3.12.
    void SetByGeantAngles(Int_t lay,Int_t lad,Int_t det,
			  const Double_t *ang){
	SetByGeantAngles(GetModuleIndex(lay,lad,det),ang);}
//
    //     This function sets a new translation vector, given by the
    // array x[3], for the Cartesian coordinate transformation
    // for a give module index.
    void SetTrans(Int_t index,Double_t x[]){
	GetGeomMatrix(index)->SetTranslation(x);}
    //     This function sets a new translation vector, given by the three
    // variables x, y, and z, for the Cartesian coordinate transformation
    // for the detector defined by layer, ladder and detector.
    void SetTrans(Int_t lay,Int_t lad,Int_t det,
                  Float_t x,Float_t y,Float_t z){Double_t t[3];
                  t[0] = x;t[1] = y;t[2] = z;
                  SetTrans(GetModuleIndex(lay,lad,det),t);}
//
    //     This function adds one more shape element to the TObjArray
    // fShape. It is primarily used in the constructor functions of the
    // AliITSgeom class. The pointer *shape can be the pointer to any
    // class that is derived from TObject (this is true for nearly every
    // ROOT class). This does not appear to be working properly at this time.
    void AddShape(TObject *shp){fShape->AddLast(shp);}
    //     This function deletes an existing shape element, of type TObject,
    // and replaces it with the one specified. This is primarily used to
    // changes the parameters to the segmentation class for a particular
    // type of detector.
    void ReSetShape(Int_t dtype,TObject *shp){
         fShape->RemoveAt(dtype);fShape->AddAt(shp,dtype);}
//
//  transformations
    //     Transforms from the ALICE Global coordinate system
    // to the detector local coordinate system for the detector
    // defined by the layer, ladder, and detector numbers. The
    // global and local coordinate are given in two floating point
    // arrays g[3], and l[3].
    void GtoL(Int_t lay,Int_t lad,Int_t det,
	      const Float_t *g,Float_t *l){
         GtoL(GetModuleIndex(lay,lad,det),g,l);}
    //     Transforms from the ALICE Global coordinate system
    // to the detector local coordinate system for the detector
    // defined by the id[0], id[1], and id[2] numbers. The
    // global and local coordinate are given in two floating point
    // arrays g[3], and l[3].
    void GtoL(const Int_t *id,const Float_t *g,Float_t *l){
         GtoL(GetModuleIndex(id),g,l);}
    //     Transforms from the ALICE Global coordinate system
    // to the detector local coordinate system for the detector
    // module index number. The global and local coordinate are
    // given in two floating point arrays g[3], and l[3].
    void GtoL(Int_t index,const Float_t *g,Float_t *l){
         Double_t dg[3],dl[3];Int_t i;for(i=0;i<3;i++) dg[i] = g[i];
         GetGeomMatrix(index)->GtoLPosition(dg,dl);
         for(i=0;i<3;i++) l[i] =dl[i];}
    //     Transforms from the ALICE Global coordinate system
    // to the detector local coordinate system for the detector
    // defined by the layer, ladder, and detector numbers. The
    // global and local coordinate are given in two Double point
    // arrays g[3], and l[3].
    void GtoL(Int_t lay,Int_t lad,Int_t det,
	      const Double_t *g,Double_t *l){
         GtoL(GetModuleIndex(lay,lad,det),g,l);}
    //     Transforms from the ALICE Global coordinate system
    // to the detector local coordinate system for the detector
    // defined by the id[0], id[1], and id[2] numbers. The
    // global and local coordinate are given in two Double point
    // arrays g[3], and l[3].
    void GtoL(const Int_t *id,const Double_t *g,Double_t *l){
         GtoL(GetModuleIndex(id),g,l);}
    //     Transforms from the ALICE Global coordinate system
    // to the detector local coordinate system for the detector
    // module index number. The global and local coordinate are
    // given in two Double point arrays g[3], and l[3].
    void GtoL(Int_t index,const Double_t *g,Double_t *l){
         Double_t dg[3],dl[3];Int_t i;for(i=0;i<3;i++) dg[i] = g[i];
         GetGeomMatrix(index)->GtoLPosition(dg,dl);
         for(i=0;i<3;i++) l[i] =dl[i];}
//
    //     Transforms from the ALICE Global coordinate system
    // to the detector local coordinate system (used for ITS tracking)
    // for the detector module index number. The global and local
    // coordinate are given in two Double point arrays g[3], and l[3].
    void GtoLtracking(Int_t index,const Double_t *g,Double_t *l){
	 if(IsGeantToTracking()) GtoL(index,g,l);
	 else GetGeomMatrix(index)->GtoLPositionTracking(g,l);}
    //     Transforms from the ALICE Global coordinate system
    // to the detector local coordinate system (used for ITS tracking)
    // for the detector id[3]. The global and local
    // coordinate are given in two Double point arrays g[3], and l[3].
    void GtoLtracking(const Int_t *id,const Double_t *g,Double_t *l){
	 GtoLtracking(GetModuleIndex(id),g,l);}
    //     Transforms from the ALICE Global coordinate system
    // to the detector local coordinate system (used for ITS tracking)
    // for the detector layer ladder and detector numbers. The global
    // and local coordinate are given in two Double point arrays g[3],
    // and l[3].
    void GtoLtracking(Int_t lay,Int_t lad,Int_t det,
		      const Double_t *g,Double_t *l){
	 GtoLtracking(GetModuleIndex(lay,lad,det),g,l);}
//
    //     Transforms of momentum types of quantities from the ALICE
    // Global coordinate system to the detector local coordinate system
    // for the detector layer ladder and detector numbers. The global
    // and local coordinate are given in two float point arrays g[3],
    // and l[3].
    void GtoLMomentum(Int_t lay,Int_t lad,Int_t det,
		      const Float_t *g,Float_t *l){
                         GtoLMomentum(GetModuleIndex(lay,lad,det),g,l);}
    //     Transforms of momentum types of quantities from the ALICE
    // Global coordinate system to the detector local coordinate system
    // for the detector module index number. The global and local
    // coordinate are given in two float point arrays g[3], and l[3].
    void GtoLMomentum(Int_t index,const Float_t *g,Float_t *l){
         Double_t dg[3],dl[3];Int_t i;for(i=0;i<3;i++) dg[i] = g[i];
         GetGeomMatrix(index)->GtoLMomentum(dg,dl);
         for(i=0;i<3;i++) l[i] =dl[i];}
    //     Transforms of momentum types of quantities from the ALICE
    // Global coordinate system to the detector local coordinate system
    // for the detector layer ladder and detector numbers. The global
    // and local coordinate are given in two Double point arrays g[3],
    // and l[3].
    void GtoLMomentum(Int_t lay,Int_t lad,Int_t det,
		      const Double_t *g,Double_t *l){
         GtoLMomentum(GetModuleIndex(lay,lad,det),g,l);}
    //     Transforms of momentum types of quantities from the ALICE
    // Global coordinate system to the detector local coordinate system
    // for the detector module index number. The global and local
    // coordinate are given in two Double point arrays g[3], and l[3].
    void GtoLMomentum(Int_t index,const Double_t *g,Double_t *l){
         Double_t dg[3],dl[3];Int_t i;for(i=0;i<3;i++) dg[i] = g[i];
         GetGeomMatrix(index)->GtoLMomentum(dg,dl);
         for(i=0;i<3;i++) l[i] =dl[i];}
//
    //     Transforms of momentum types of quantities from the ALICE
    // Global coordinate system to the detector local coordinate system
    // (used for ITS tracking) for the detector module index number.
    // The global and local coordinate are given in two Double point
    // arrays g[3], and l[3].
    void GtoLMomentumTracking(Int_t index,const Double_t *g,Double_t *l){
         if(IsGeantToTracking()) GtoLMomentum(index,g,l);
         else GetGeomMatrix(index)->GtoLMomentumTracking(g,l);}
    //     Transforms of momentum types of quantities from the ALICE
    // Global coordinate system to the detector local coordinate system
    // (used for ITS tracking) for the detector id[3].
    // The global and local coordinate are given in two Double point
    // arrays g[3], and l[3].
    void GtoLMomentumTracking(const Int_t *id,const Double_t *g,Double_t *l){
                 GtoLMomentumTracking(GetModuleIndex(id),g,l);}
    //     Transforms of momentum types of quantities from the ALICE
    // Global coordinate system to the detector local coordinate system
    // (used for ITS tracking) for the detector layer ladder and detector
    // numbers. The global and local coordinate are given in two Double point
    // arrays g[3], and l[3].
    void GtoLMomentumTracking(Int_t lay,Int_t lad,Int_t det,
			      const Double_t *g,Double_t *l){
                        GtoLMomentumTracking(GetModuleIndex(lay,lad,det),g,l);}
//
    //     Transforms from the detector local coordinate system
    // to the ALICE Global coordinate  system for the detector
    // defined by the layer, ladder, and detector numbers. The
    // global and local coordinate are given in two floating point
    // arrays g[3], and l[3].
    void LtoG(Int_t lay,Int_t lad,Int_t det,
	      const Float_t *l,Float_t *g){
                     LtoG(GetModuleIndex(lay,lad,det),l,g);}
    //     Transforms from the detector local coordinate system
    // to the ALICE Global coordinate system for the detector
    // defined by the id[0], id[1], and id[2] numbers. The
    // global and local coordinate are given in two floating point
    // arrays g[3], and l[3].
    void LtoG(const Int_t *id,const Float_t *l,Float_t *g){
                     LtoG(GetModuleIndex(id),l,g);}
    //     Transforms from the detector local coordinate system
    // to the ALICE Global coordinate system for the detector
    // module index number. The global and local coordinate are
    // given in two floating point arrays g[3], and l[3].
    void LtoG(Int_t index,const Float_t *l,Float_t *g){
         Double_t dg[3],dl[3];Int_t i;for(i=0;i<3;i++) dl[i] = l[i];
         GetGeomMatrix(index)->LtoGPosition(dl,dg);
         for(i=0;i<3;i++) g[i] =dg[i];}
    //     Transforms from the detector local coordinate system
    // to the ALICE Global coordinate system for the detector
    // defined by the layer, ladder, and detector numbers. The
    // global and local coordinate are given in two Double point
    // arrays g[3], and l[3].
    void LtoG(Int_t lay,Int_t lad,Int_t det,
	      const Double_t *l,Double_t *g){
                      LtoG(GetModuleIndex(lay,lad,det),l,g);}
    //     Transforms from the detector local coordinate system
    // to the ALICE Global coordinate system for the detector
    // defined by the id[0], id[1], and id[2] numbers. The
    // global and local coordinate are given in two Double point
    // arrays g[3], and l[3].
    void LtoG(const Int_t *id,const Double_t *l,Double_t *g){
                       LtoG(GetModuleIndex(id),l,g);}
    //     Transforms from the detector local coordinate system
    // to the ALICE Global coordinate system for the detector
    // module index number. The global and local coordinate are
    // given in two Double point arrays g[3], and l[3].
    void LtoG(Int_t index,const Double_t *l,Double_t *g){
         Double_t dg[3],dl[3];Int_t i;for(i=0;i<3;i++) dl[i] = l[i];
         GetGeomMatrix(index)->LtoGPosition(dl,dg);
         for(i=0;i<3;i++) g[i] =dg[i];}
//
    //     Transforms from the detector local coordinate system (used
    // for ITS tracking) to the ALICE Global coordinate system 
    // for the detector module index number. The global and local
    // coordinate are given in two Double point arrays g[3], and l[3].
    void LtoGtracking(Int_t index,const Double_t *l,Double_t *g){
	 if(IsGeantToTracking()) LtoG(index,l,g);
	 else GetGeomMatrix(index)->LtoGPositionTracking(l,g);}
    //     Transforms from the detector local coordinate system (used
    // for ITS tracking) to the ALICE Global coordinate system 
    // for the detector id[3]. The global and local
    // coordinate are given in two Double point arrays g[3], and l[3].
    void LtoGtracking(const Int_t *id,const Double_t *l,Double_t *g){
	 LtoGtracking(GetModuleIndex(id),l,g);}
    //     Transforms from the detector local coordinate system (used
    // for ITS tracking) to the detector local coordinate system
    // for the detector layer ladder and detector numbers. The global
    // and local coordinate are given in two Double point arrays g[3],
    // and l[3].
    void LtoGtracking(Int_t lay,Int_t lad,Int_t det,
		      const Double_t *l,Double_t *g){
	 LtoGtracking(GetModuleIndex(lay,lad,det),l,g);}
//
    //     Transforms of momentum types of quantities from the detector
    // local coordinate system to the ALICE Global coordinate system
    // for the detector layer ladder and detector numbers. The global
    // and local coordinate are given in two float point arrays g[3],
    // and l[3].
    void LtoGMomentum(Int_t lay,Int_t lad,Int_t det,
		      const Float_t *l,Float_t *g){
         LtoGMomentum(GetModuleIndex(lay,lad,det),l,g);}
    //     Transforms of momentum types of quantities from the detector
    // local coordinate system to the ALICE Global coordinate system
    // for the detector module index number. The global and local
    // coordinate are given in two float point arrays g[3], and l[3].
    void LtoGMomentum(Int_t index,const Float_t *l,Float_t *g){
         Double_t dg[3],dl[3];Int_t i;for(i=0;i<3;i++) dl[i] = l[i];
         GetGeomMatrix(index)->LtoGMomentum(dl,dg);
         for(i=0;i<3;i++) g[i] =dg[i];}
    //     Transforms of momentum types of quantities from the detector
    // local coordinate system to the ALICE Global coordinate system
    // for the detector layer ladder and detector numbers. The global
    // and local coordinate are given in two Double point arrays g[3],
    // and l[3].
    void LtoGMomentum(Int_t lay,Int_t lad,Int_t det,
			   const Double_t *l,Double_t *g){
                        LtoGMomentum(GetModuleIndex(lay,lad,det),l,g);}
    //     Transforms of momentum types of quantities from the detector
    // local coordinate system to the ALICE Global coordinate system
    // for the detector module index number. The global and local
    // coordinate are given in two Double point arrays g[3], and l[3].
    void LtoGMomentum(Int_t index,const Double_t *l,Double_t *g){
         GetGeomMatrix(index)->LtoGMomentum(l,g);}
//
    //     Transforms of momentum types of quantities from the detector 
    // local coordinate system (used for ITS tracking) to the detector
    // system ALICE Global for the detector module index number.
    // The global and local coordinate are given in two Double point
    // arrays g[3], and l[3].
    void LtoGMomentumTracking(Int_t index,const Double_t *l,Double_t *g){
         if(IsGeantToTracking()) LtoGMomentum(index,l,g);
         else GetGeomMatrix(index)->LtoGMomentumTracking(l,g);}
    //     Transforms of momentum types of quantities from the detector
    // local coordinate system (used for ITS tracking) to the ALICE
    // Global coordinate system for the detector id[3].
    // The global and local coordinate are given in two Double point
    // arrays g[3], and l[3].
    void LtoGMomentumTracking(const Int_t *id,const Double_t *l,Double_t *g){
                 LtoGMomentumTracking(GetModuleIndex(id),l,g);}
    //     Transforms of momentum types of quantities from the detector
    // local coordinate system (used for ITS tracking) to the ALICE
    // Global coordinate system for the detector layer ladder and detector
    // numbers. The global and local coordinate are given in two Double point
    // arrays g[3], and l[3].
    void LtoGMomentumTracking(Int_t lay,Int_t lad,Int_t det,
			      const Double_t *l,Double_t *g){
                        LtoGMomentumTracking(GetModuleIndex(lay,lad,det),l,g);}
//
    //     Transforms from one detector local coordinate system
    // to another detector local coordinate system for the detector
    // module index1 number to the detector module index2 number. The
    //  local coordinates are given in two Double point arrays l1[3],
    // and l2[3].
    void LtoL(Int_t index1,Int_t index2,Double_t *l1,Double_t *l2){
         Double_t g[3]; LtoG(index1,l1,g);GtoL(index2,g,l2);}
    //     Transforms from one detector local coordinate system
    // to another detector local coordinate system for the detector
    // id1[3] to the detector id2[3]. The local coordinates are given
    // in two Double point arrays l1[3], and l2[3].
    void LtoL(const Int_t *id1,const Int_t *id2,Double_t *l1,Double_t *l2){
         LtoL(GetModuleIndex(id1[0],id1[1],id1[2]),
              GetModuleIndex(id2[0],id2[1],id2[2]),l1,l2);}
//
    //     Transforms from one detector local coordinate system (used for
    // ITS tracking) to another detector local coordinate system (used
    // for ITS tracking) for the detector module index1 number to the
    // detector module index2 number. The local coordinates are given
    // in two Double point arrays l1[3], and l2[3].
    void LtoLtracking(Int_t index1,Int_t index2,
			   Double_t *l1,Double_t *l2){
         Double_t g[3]; LtoGtracking(index1,l1,g);GtoLtracking(index2,g,l2);}
    //     Transforms from one detector local coordinate system (used for
    // ITS tracking) to another detector local coordinate system (used
    // for ITS tracking) for the detector id1[3] to the detector id2[3].
    // The local coordinates are given in two Double point arrays l1[3],
    // and l2[3].
    void LtoLtracking(const Int_t *id1,const Int_t *id2,
			   Double_t *l1,Double_t *l2){
         LtoLtracking(GetModuleIndex(id1[0],id1[1],id1[2]),
              GetModuleIndex(id2[0],id2[1],id2[2]),l1,l2);}
//
    //     Transforms of momentum types of quantities from one detector
    // local coordinate system to another detector local coordinate
    // system for the detector module index1 number to the detector
    // module index2 number. The local coordinates are given in two
    // Double point arrays l1[3], and l2[3].
    void LtoLMomentum(Int_t index1,Int_t index2,
		      const Double_t *l1,Double_t *l2){
         Double_t g[3]; LtoGMomentum(index1,l1,g);GtoLMomentum(index2,g,l2);}
    //     Transforms of momentum types of quantities from one detector
    // local coordinate system to another detector local coordinate
    // system for the detector id1[3] to the detector id2[3]. The local
    // coordinates are given in two Double point arrays l1[3], and l2[3].
    void LtoLMomentum(const Int_t *id1,const Int_t *id2,
		      const Double_t *l1,Double_t *l2){
         LtoLMomentum(GetModuleIndex(id1[0],id1[1],id1[2]),
                      GetModuleIndex(id2[0],id2[1],id2[2]),l1,l2);}
//
    //     Transforms of momentum types of quantities from one detector
    // local coordinate system (used by ITS tracking) to another detector
    // local coordinate system (used by ITS tracking) for the detector
    // module index1 number to the detector module index2 number. The
    // local coordinates are given in two Double point arrays l1[3],
    // and l2[3].
    void LtoLMomentumTracking(Int_t index1,Int_t index2,
			   Double_t *l1,Double_t *l2){
         Double_t g[3]; LtoGMomentumTracking(index1,l1,g);
                        GtoLMomentumTracking(index2,g,l2);}
    //     Transforms of momentum types of quantities from one detector
    // local coordinate system (used by ITS tracking) to another detector
    // local coordinate system (used by ITS tracking) for the detector
    // id1[3] to the detector id2[3]. The local coordinates are given in
    // two Double point arrays l1[3], and l2[3].
    void LtoLMomentumTracking(const Int_t *id1,const Int_t *id2,
			   Double_t *l1,Double_t *l2){
         LtoLMomentumTracking(GetModuleIndex(id1[0],id1[1],id1[2]),
                              GetModuleIndex(id2[0],id2[1],id2[2]),l1,l2);}
//
    //     Transforms a matrix, like an Uncertainty or Error matrix from
    // the ALICE Global coordinate system to a detector local coordinate
    // system. The specific detector is determined by the module index
    // number.
    void GtoLErrorMatrix(Int_t index,const Double_t **g,Double_t **l){
         GetGeomMatrix(index)->GtoLPositionError((Double_t (*)[3])g,(Double_t (*)[3])l);}
//
    //     Transforms a matrix, like an Uncertainty or Error matrix from
    // the ALICE Global coordinate system to a detector local coordinate
    // system (used by ITS tracking). The specific detector is determined
    // by the module index number.
    void GtoLErrorMatrixTracking(Int_t index,const Double_t **g,
				 Double_t **l){
	if(IsGeantToTracking()) GetGeomMatrix(index)->GtoLPositionError((
	    Double_t (*)[3])g,(Double_t (*)[3])l);
	else GetGeomMatrix(index)->GtoLPositionErrorTracking(
	     (Double_t (*)[3])g,(Double_t (*)[3])l);}
//
    //     Transforms a matrix, like an Uncertainty or Error matrix from
    // the detector local coordinate system to a ALICE Global coordinate
    // system. The specific detector is determined by the module index
    // number.
    void LtoGErrorMatrix(Int_t index,const Double_t **l,Double_t **g){
         GetGeomMatrix(index)->LtoGPositionError((Double_t (*)[3])l,(Double_t (*)[3])g);}
//
    //     Transforms a matrix, like an Uncertainty or Error matrix from
    // the detector local coordinate system (used by ITS tracking) to a
    // ALICE Global coordinate system. The specific detector is determined
    // by the module index number.
    void LtoGErrorMatrixTracking(Int_t index,const Double_t **l,
				 Double_t **g){
         if(IsGeantToTracking()) GetGeomMatrix(index)->LtoGPositionError((
	    Double_t (*)[3])g,(Double_t (*)[3])l);
	else GetGeomMatrix(index)->LtoGPositionErrorTracking((Double_t (*)[3])l,
					       (Double_t (*)[3])g);}
//
    //     Transforms a matrix, like an Uncertainty or Error matrix from
    // one detector local coordinate system to another detector local
    // coordinate system. The specific detector is determined by the
    // two module index number index1 and index2.
    void LtoLErrorMatrix(Int_t index1,Int_t index2,
			 const Double_t **l1,Double_t **l2){
	Double_t g[3][3];
                  LtoGErrorMatrix(index1,l1,(Double_t **)g);
                  GtoLErrorMatrix(index2,(const Double_t **)g,l2);}
//
    //     Transforms a matrix, like an Uncertainty or Error matrix from
    // one detector local coordinate system (used by ITS tracking) to
    // another detector local coordinate system (used by ITS tracking).
    // The specific detector is determined by the two module index number
    // index1 and index2.
    void LtoLErrorMatrixTraking(Int_t index1,Int_t index2,
			 const Double_t **l1,Double_t **l2){Double_t g[3][3];
                  LtoGErrorMatrixTracking(index1,l1,(Double_t **)g);
                  GtoLErrorMatrixTracking(index2,(const Double_t **)g,l2);}
//  Find Specific Modules
    // Locate the nearest module to the point g, in ALICE global Cartesian
    // coordinates [cm] in a give layer. If layer = 0 then it search in
    // all layers.
    Int_t    GetNearest(const Double_t g[3],Int_t lay=0);
    // Locates the nearest 27 modules, in nearest order, to the point g, in
    // ALICE global Cartesian coordinates [cm] in a give layer. If layer = 0
    // then it searches in all layers. (there are 27 elements in a 3x3x3
    // cube.
    void     GetNearest27(const Double_t g[3],Int_t n[27],Int_t lay=0);
    // Returns the distance [cm] between the point g[3] and the center of
    // the detector/module specified by the the module index number.
    Double_t Distance(Int_t index,const Double_t g[3]){
         return  TMath::Sqrt(GetGeomMatrix(index)->Distance2(g));}
//  Geometry manipulation
    // This function performs a Cartesian translation and rotation of
    // the full ITS from its default position by an amount determined by
    // the three element arrays tran and rot.
    void GlobalChange(const Float_t  *tran,const Float_t  *rot);
    // This function performs a Cylindrical translation and rotation of
    // the full ITS from its default position by an amount determined by
    // the three element arrays tran and rot.
    void GlobalCylindericalChange(const Float_t *tran,const Float_t *rot);
    // This function performs a Gaussian random displacement and/or
    // rotation about the present global position of each active
    // volume/detector of the ITS with variances given by stran and srot.
    void RandomChange(const Float_t *stran,const Float_t *srot);
    // This function performs a Gaussian random displacement and/or
    // rotation about the present global position of each active
    // volume/detector of the ITS with variances given by stran and srot.
    // But in Cylindrical coordinates.
    void RandomCylindericalChange(const Float_t *stran,const Float_t *srot);
    // This function converts these transformations from Alice global and
    // local to Tracking global and local.
    void GeantToTracking(AliITSgeom &source); // This converts the geometry
//  Other routines.
    // This routine prints, to a file, the difference between this class
    // and "other".
    void PrintComparison(FILE *fp,AliITSgeom *other);
    // This routine prints, to a file, the contents of this class.
    void PrintData(FILE *fp,Int_t lay,Int_t lad,Int_t det);
    // This function prints out this class in a single stream. This steam
    // can be read by ReadGeom.
    ofstream &PrintGeom(ofstream &out);
    // This function reads in that single steam printed out by PrintGeom.
    ifstream &ReadGeom(ifstream &in);

 private:
    char       fVersion[20];// Transformation version.
    Int_t      fTrans;   // Flag to keep track of which transformation 
    Int_t      fNmodules;// The total number of modules
    Int_t      fNlayers; // The number of layers.
    Int_t     *fNlad;  //[fNlayers] Array of the number of ladders/layer(layer)
    Int_t     *fNdet;//[fNlayers] Array of the number of detector/ladder(layer)
    TObjArray *fGm;      // Structure of translation. and rotation.
    TObjArray *fShape;   // Array of shapes and detector information.

    ClassDef(AliITSgeom,2) // ITS geometry class
};

#endif
