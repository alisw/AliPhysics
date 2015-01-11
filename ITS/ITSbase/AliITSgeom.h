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
#include <TString.h>
#include <TArrayI.h>
#include <TMath.h>
//
#include "AliITSgeomMatrix.h"
#include "AliLog.h"

typedef enum {kND=-1,kSPD=0, kSDD=1, kSSD=2, kSSDp=3,kSDDp=4, kUPG=5} AliITSDetector;

//_______________________________________________________________________

class AliITSgeom : public TObject {

 public:
    AliITSgeom();                      // Default constructor
    AliITSgeom(Int_t itype,Int_t nlayers,const Int_t *nlads,const Int_t *ndets,
               Int_t nmods); // Constructor
    AliITSgeom(const AliITSgeom &source);    // Copy constructor
    AliITSgeom& operator=(const AliITSgeom &source);// = operator
    virtual ~AliITSgeom();             // Default destructor
    // Zero and reinitilizes this class.
    void Init(Int_t itype,Int_t nlayers,const Int_t *nlads,
              const Int_t *ndets,Int_t mods);
    // this function allocates a AliITSgeomMatrix for a particular module.
    void CreateMatrix(Int_t mod,Int_t lay,Int_t lad,Int_t det,
                      AliITSDetector idet,const Double_t tran[3],
                      const Double_t rot[10]);
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
    Bool_t IsGeantToTracking()  const {return ((fTrans&0xfffe)!= 0);}
    // returns kTRUE if the transformation defined by this class is
    // for Global GEANT coordinate system to the local GEANT coordinate system
    // of the detector but may have been displaced by some typically small
    // amount. These are modified transformation similar to that used by GEANT.
    Bool_t IsGeantToDisplaced() const {return ((fTrans&0xfffd)!= 0);}
    //
    //     This function returns a pointer to the particular AliITSgeomMatrix
    // class for a specific module index.
    AliITSgeomMatrix *GetGeomMatrix(Int_t index){if(index<fGm.GetSize()&&index>=0)
        return (AliITSgeomMatrix*)(fGm.At(index));else 
            Error("GetGeomMatrix","index=%d<0||>=GetSize()=%d",index,fGm.GetSize());return 0;}
    AliITSgeomMatrix *GetGeomMatrix(Int_t index)const{if(index<fGm.GetSize()&&index>=0)
        return (AliITSgeomMatrix*)(fGm.At(index));else 
            Error("GetGeomMatrix","index=%d<0||>=GetSize()=%d",index,fGm.GetSize());return 0;}
    // This function find and return the number of detector types only.
    Int_t GetNDetTypes()const{Int_t max;return GetNDetTypes(max);};
    // This function find and return the number of detector types and the
    // maximum det type value.
    Int_t GetNDetTypes(Int_t &max)const;
    // This function finds and return the number of detector types and the
    // and the number of each type in the TArrayI and their types.
    Int_t GetNDetTypes(TArrayI &maxs,AliITSDetector *types)const;
    //     This function returns the number of detectors/ladder for a give 
    // layer. In particular it returns fNdet[layer-1].
    Int_t GetNdetectors(Int_t lay) const {return fNdet[lay-1];}
    //     This function returns the number of ladders for a give layer. In
    // particular it returns fNlad[layer-1].
    Int_t GetNladders(Int_t lay)   const {return fNlad[lay-1];};
    //     This function returns the number of layers defined in the ITS
    // geometry. In particular it returns fNlayers.
    Int_t GetNlayers()                   const {return fNlayers;}
    Int_t GetModuleIndex(Int_t lay,Int_t lad,Int_t det)const;
    //     This function returns the module index number given the layer,
    // ladder and detector numbers put into the array id[3].
    Int_t GetModuleIndex(const Int_t *id)const{
        return GetModuleIndex(id[0],id[1],id[2]);}
    void  GetModuleId(Int_t index,Int_t &lay,Int_t &lad,Int_t &det)const;
    // Returns the detector type
    //Int_t GetModuleType(Int_t index)const{
    //    return GetGeomMatrix(index)->GetDetectorIndex();}
    AliITSDetector GetModuleType(Int_t index)const{
        return (AliITSDetector)(GetGeomMatrix(index)->GetDetectorIndex());}
    // Returns the detector type as a string
    const char * GetModuleTypeName(Int_t index)const{
        return GetDetectorTypeName(GetModuleType(index));} 
    // Returns the detector type as a string
    const char * GetDetectorTypeName(Int_t index)const{switch(index) {
    case kSPD : return "kSPD" ; case kSDD : return "kSDD" ;
    case kSSD : return "kSSD" ; case kSSDp: return "kSSDp";
    case kSDDp: return "kSDDp"; default   : return "Undefined";};}
    //
    Int_t GetStartDet(Int_t dtype )const;
    Int_t GetLastDet(Int_t dtype)const;
    //     Returns the starting module index number for SPD detector,
    // assuming the modules are placed in the "standard" cylindrical
    // ITS structure.
    Int_t GetStartSPD()const{return GetStartDet(kSPD);}
    //     Returns the ending module index number for SPD detector,
    // assuming the modules are placed in the "standard" cylindrical
    // ITS structure.
    Int_t GetLastSPD()const{return GetLastDet(kSPD);}
    //     Returns the starting module index number for SDD detector,
    // assuming the modules are placed in the "standard" cylindrical
    // ITS structure.
    Int_t GetStartSDD()const{return GetStartDet(kSDD);}
    //     Returns the ending module index number for SDD detector,
    // assuming the modules are placed in the "standard" cylindrical
    // ITS structure.
    Int_t GetLastSDD()const{return GetLastDet(kSDD);}
    //     Returns the starting module index number for SSD detector,
    // assuming the modules are placed in the "standard" cylindrical
    // ITS structure.
    Int_t GetStartSSD()const{return GetStartDet(kSSD);}
    //     Returns the ending module index number for SSD detector,
    // assuming the modules are placed in the "standard" cylindrical
    // ITS structure.
    Int_t GetLastSSD()const{return GetLastDet(kSSD);}
    //     Returns the last module index number.
    Int_t GetIndexMax() const {return fNmodules;}
    //
    //     This function returns the rotation angles for a give module 
    // in the Double point array ang[3]. The angles are in radians
    void  GetAngles(Int_t index,Double_t *ang)const{
        GetGeomMatrix(index)->GetAngles(ang);}
    //     This function returns the rotation angles for a give module
    // in the three floating point variables provided. rx = frx,
    // fy = fry, rz = frz. The angles are in radians
    void  GetAngles(Int_t index,Float_t &rx,Float_t &ry,Float_t &rz)const{
        Double_t a[3];GetAngles(index,a);rx = a[0];ry = a[1];rz = a[2];}
    //     This function returns the rotation angles for a give detector on
    // a give ladder in a give layer in the three floating point variables
    // provided. rx = frx, fy = fry, rz = frz. The angles are in radians
    void  GetAngles(Int_t lay,Int_t lad,Int_t det,
                    Float_t &rx,Float_t &ry,Float_t &rz)const{
        GetAngles(GetModuleIndex(lay,lad,det),rx,ry,rz);}
    //
    //     This function returns the 6 GEANT rotation angles for a give 
    // module in the double point array ang[3]. The angles are in degrees
    void  GetGeantAngles(Int_t index,Double_t *ang)const{
	GetGeomMatrix(index)->SixAnglesFromMatrix(ang);}
    //
    //     This function returns the Cartesian translation for a give
    // module in the Double array t[3]. The units are
    // those of the Monte Carlo, generally cm.
    void  GetTrans(Int_t index,Double_t *t)const{
        GetGeomMatrix(index)->GetTranslation(t);}
    //     This function returns the Cartesian translation for a give
    // module index in the three floating point variables provided.
    // x = fx0, y = fy0, z = fz0. The units are those of the Mont
    // Carlo, generally cm.
    void  GetTrans(Int_t index,Float_t &x,Float_t &y,Float_t &z)const{
        Double_t t[3];GetTrans(index,t);x = t[0];y = t[1];z = t[2];}
    //     This function returns the Cartesian translation for a give
    // detector on a give ladder in a give layer in the three floating
    // point variables provided. x = fx0, y = fy0, z = fz0. The units are
    // those of the Monte Carlo, generally cm.
    void  GetTrans(Int_t lay,Int_t lad,Int_t det,
                   Float_t &x,Float_t &y,Float_t &z)const{
        GetTrans(GetModuleIndex(lay,lad,det),x,y,z);}
    //
    //     This function returns the Cartesian translation for a give
    // module in the Double array t[3]. The units are
    // those of the Monte Carlo, generally cm.
    void  GetTransCyln(Int_t index,Double_t *t)const{
        GetGeomMatrix(index)->GetTranslationCylinderical(t);}
    //     This function returns the Cartesian translation for a give
    // module index in the three floating point variables provided.
    // x = fx0, y = fy0, z = fz0. The units are those of the Mont
    // Carlo, generally cm.
    void  GetTransCyln(Int_t index,Float_t &x,Float_t &y,Float_t &z)const{
        Double_t t[3];GetTransCyln(index,t);x = t[0];y = t[1];z = t[2];}
    //     This function returns the Cartesian translation for a give
    // detector on a give ladder in a give layer in the three floating
    // point variables provided. x = fx0, y = fy0, z = fz0. The units are
    // those of the Monte Carlo, generally cm.
    void  GetTransCyln(Int_t lay,Int_t lad,Int_t det,
                       Float_t &x,Float_t &y,Float_t &z)const{
        GetTransCyln(GetModuleIndex(lay,lad,det),x,y,z);}
    //
    //      This function returns the Cartesian translation [cm] and the
    // 6 GEANT rotation angles [degrees]for a given layer ladder and
    // detector number, in the TVector x (at least 9 elements large).
    // This function is required to be in-lined for speed.
    void  GetCenterThetaPhi(Int_t lay,Int_t lad,Int_t det,TVector &x)const{
        Double_t t[3],a[6];Int_t i=GetModuleIndex(lay,lad,det);GetTrans(i,t);
        GetGeantAngles(i,a);x(0)=t[0];x(1)=t[1];x(2)=t[2];x(3)=a[0];x(4)=a[1];
        x(5)=a[2];x(6)=a[3];x(7)=a[4];x(8)=a[5];}
    //
    //     This function returns the rotation matrix in Double
    // precision for a given module.
    void  GetRotMatrix(Int_t index,Double_t mat[3][3])const{
        GetGeomMatrix(index)->GetMatrix(mat);}
    //     This function returns the rotation matrix in a Double
    // precision pointer for a given module. mat[i][j] => mat[3*i+j].
    void  GetRotMatrix(Int_t index,Double_t *mat)const{
        Double_t rot[3][3];GetRotMatrix(index,rot);
        for(Int_t i=0;i<3;i++)for(Int_t j=0;j<3;j++) mat[3*i+j] = rot[i][j];}
    //     This function returns the rotation matrix in a floating 
    // precision pointer for a given layer ladder and detector module.
    // mat[i][j] => mat[3*i+j].
    void  GetRotMatrix(Int_t lay,Int_t lad,Int_t det,Float_t *mat)const{
        GetRotMatrix(GetModuleIndex(lay,lad,det),mat);}
    //     This function returns the rotation matrix in a Double
    // precision pointer for a given layer ladder and detector module.
    // mat[i][j] => mat[3*i+j].
    void  GetRotMatrix(Int_t lay,Int_t lad,Int_t det,Double_t *mat)const{
        GetRotMatrix(GetModuleIndex(lay,lad,det),mat);}
    //     This function returns the rotation matrix in a floating
    // precision pointer for a given module. mat[i][j] => mat[3*i+j].
    void  GetRotMatrix(Int_t index,Float_t *mat)const{
          Double_t rot[3][3];
          GetGeomMatrix(index)->GetMatrix(rot);
          for(Int_t i=0;i<3;i++)for(Int_t j=0;j<3;j++) mat[3*i+j] = rot[i][j];}
    //     This function sets the rotation matrix in a Double
    // precision pointer for a given module. mat[i][j] => mat[3*i+j].
    void  SetRotMatrix(Int_t index,const Double_t *mat){Double_t rot[3][3];
          for(Int_t i=0;i<3;i++)for(Int_t j=0;j<3;j++) rot[i][j]=mat[3*i+j];
          GetGeomMatrix(index)->SetMatrix(rot);}
    // Return the normal for a specific module
    void GetGlobalNormal(Int_t index,Double_t n[3]){
        GetGeomMatrix(index)->GetGlobalNormal(n[0],n[1],n[2]);}
    //
    //
    //  Setters
    //     Sets the rotation angles and matrix for a give module index
    // via the double precision array a[3] [radians].
    void SetByAngles(Int_t index,const Double_t a[]){
        GetGeomMatrix(index)->SetAngles(a);}
    //     Sets the rotation angles and matrix for a give module index
    // via the 3 floating precision variables rx, ry, and rz [radians].
    void SetByAngles(Int_t index,Float_t rx, Float_t ry, Float_t rz) {
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
    //  transformations
    //     Transforms from the ALICE Global coordinate system
    // to the detector local coordinate system for the detector
    // defined by the layer, ladder, and detector numbers. The
    // global and local coordinate are given in two floating point
    // arrays g[3], and l[3].
    void GtoL(Int_t lay,Int_t lad,Int_t det,
              const Float_t *g,Float_t *l)const{
        GtoL(GetModuleIndex(lay,lad,det),g,l);}
    //     Transforms from the ALICE Global coordinate system
    // to the detector local coordinate system for the detector
    // defined by the id[0], id[1], and id[2] numbers. The
    // global and local coordinate are given in two floating point
    // arrays g[3], and l[3].
    void GtoL(const Int_t *id,const Float_t *g,Float_t *l)const{
        GtoL(GetModuleIndex(id),g,l);}
    //     Transforms from the ALICE Global coordinate system
    // to the detector local coordinate system for the detector
    // module index number. The global and local coordinate are
    // given in two floating point arrays g[3], and l[3].
    void GtoL(Int_t index,const Float_t *g,Float_t *l)const{
        Double_t dg[3],dl[3];Int_t i;for(i=0;i<3;i++) dg[i] = g[i];
        GetGeomMatrix(index)->GtoLPosition(dg,dl);
        for(i=0;i<3;i++) l[i] =dl[i];}
    //     Transforms from the ALICE Global coordinate system
    // to the detector local coordinate system for the detector
    // defined by the layer, ladder, and detector numbers. The
    // global and local coordinate are given in two Double point
    // arrays g[3], and l[3].
    void GtoL(Int_t lay,Int_t lad,Int_t det,
              const Double_t *g,Double_t *l)const{
        GtoL(GetModuleIndex(lay,lad,det),g,l);}
    //     Transforms from the ALICE Global coordinate system
    // to the detector local coordinate system for the detector
    // defined by the id[0], id[1], and id[2] numbers. The
    // global and local coordinate are given in two Double point
    // arrays g[3], and l[3].
    void GtoL(const Int_t *id,const Double_t *g,Double_t *l)const{
        GtoL(GetModuleIndex(id),g,l);}
    //     Transforms from the ALICE Global coordinate system
    // to the detector local coordinate system for the detector
    // module index number. The global and local coordinate are
    // given in two Double point arrays g[3], and l[3].
    void GtoL(Int_t index,const Double_t g[3],Double_t l[3])const{
        GetGeomMatrix(index)->GtoLPosition(g,l);}
 
    //  Find Specific Modules
    // Locate the nearest module to the point g, in ALICE global Cartesian
    // coordinates [cm] in a give layer. If layer = 0 then it search in
    // all layers.
    Int_t    GetNearest(const Double_t g[3],Int_t lay=0)const;
    // Locates the nearest 27 modules, in nearest order, to the point g, in
    // ALICE global Cartesian coordinates [cm] in a give layer. If layer = 0
    // then it searches in all layers. (there are 27 elements in a 3x3x3
    // cube.
    void     GetNearest27(const Double_t g[3],Int_t n[27],Int_t lay=0)const;
    // Returns the distance [cm] between the point g[3] and the center of
    // the detector/module specified by the the module index number.
    Double_t Distance(Int_t index,const Double_t g[3])const{
        return  TMath::Sqrt(GetGeomMatrix(index)->Distance2(g));}
    // This routine prints, to a file, the contents of this class.
    void PrintData(FILE *fp,Int_t lay,Int_t lad,Int_t det)const;
    // This function prints out this class in a single stream. This steam
    // can be read by ReadGeom.
    // (Coverity warnings)    void PrintGeom(ostream *out)const;

    //Conversion from det. local coordinates to local ("V2") coordinates
    //used for tracking

    void DetLToTrackingV2(Int_t md,Float_t xin,Float_t zin,
                          Float_t &yout, Float_t &zout) const ; 

    void TrackingV2ToDetL(Int_t md,Float_t yin,Float_t zin,
                          Float_t &xout,Float_t &zout) const ;

 private:
    TString    fVersion; // Transformation version.
    Int_t      fTrans;   // Flag to keep track of which transformation 
    Int_t      fNmodules;// The total number of modules
    Int_t      fNlayers; // The number of layers.
    TArrayI    fNlad;    // Array of the number of ladders/layer(layer)
    TArrayI    fNdet;    // Array of the number of detector/ladder(layer)
    TObjArray  fGm;      // Structure of translation. and rotation.

    ClassDef(AliITSgeom,4) // ITS geometry class
}; 
// Input and output function for standard C++ input/output.

#endif
