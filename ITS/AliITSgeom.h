#ifndef ITSGEOM_H
#define ITSGEOM_H
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
#include "TObjArray.h"
#include "AliITSgeomSPD.h"
#include "AliITSgeomSDD.h"
#include "AliITSgeomSSD.h"

////////////////////////////////////////////////////////////////////////
// The structure ITS_geom:
//     The structure ITS_geom has been defined to hold all of the
// information necessary to do the coordinate transformations for one
// detector between the ALICE Cartesian global and the detector local
// coordinate systems. The rotations are implemented in the following
// order, Rz*Ry*Rx*(Vglobal-Vtrans)=Vlocal (in matrix notation). 
// In addition it contains an index to the TObjArray containing all of
// the information about the shape of the active detector volume, and
// any other useful detector parameters. See the definition of *fShape
// below and the classes AliITSgeomSPD, AliITSgeomSDD, and AliITSgeomSSD
// for a full description. This structure is not available outside of 
// these routines.
//
// Int_t fShapeIndex
//     The index to the array of detector shape information. In this way
// only an index is needed to be stored and not all of the shape
// information. This saves much space since most, if not all, of the
// detectors of a give type have the same shape information and are only
// placed in a different spot in the ALICE/ITS detector.
//
// Float_t fx0,fy0,fz0
//     The Cartesian translation vector used to define part of the
// coordinate transformation. The units of the translation are kept
// in the Monte Carlo distance units, usually cm.
//
// Float_t frx,fry,frz
//     The three rotation angles that define the rotation matrix. The
// angles are, frx the rotation about the x axis. fry the rotation about
// the "new" or "rotated" y axis. frz the rotation about the "new" or
// "rotated" z axis. These angles, although redundant with the rotation
// matrix fr, are kept for speed. This allows for their retrieval without
// having to compute them each and every time. The angles are kept in
// radians
//
// Float_t fr[9]
//     The 3x3 rotation matrix defined by the angles frx, fry, and frz,
// for the Global to Local transformation is
//    |fr[0] fr[1] fr[2]| | cos(frz)  sin(frz) 0| | cos(fry) 0  sin(fry)|
// fr=|fr[3] fr[4] fr[4]|=|-sin(frz)  cos(frz) 0|*|   0      1    0     |
//    |fr[6] fr[7] fr[8]| |   0         0      1| |-sin(fry) 0  cos(fry)|
//
//    |1    0        0     |
//   *|0  cos(frx) sin(frx)|
//    |0 -sin(frx) cos(frx)|
//
// Even though this information is redundant with the three rotation
// angles, because this transformation matrix can be used so much it is
// kept to speed things up a lot. The coordinate system used is Cartesian.
////////////////////////////////////////////////////////////////////////

struct ITS_geom {
    Int_t   fShapeIndex; // Shape index for this volume
    Float_t fx0,fy0,fz0; // Translation vector
    Float_t frx,fry,frz; // Rotation about axis, angle radians
    Float_t fr[9];       // the rotation matrix
};

//_______________________________________________________________________

class AliITSgeom : public TObject {
////////////////////////////////////////////////////////////////////////
//
// version: 0
// Written by Bjorn S. Nilsen
//
// Data Members:
//
// Int_t fNlayers
//     The number of ITS layers for this geometry. By default this
//  is 6, but can be modified by the creator function if there are
// more layers defined.
//
// Int_t *fNlad
//     A pointer to an array fNlayers long containing the number of 
// ladders for each layer. This array is typically created and filled 
// by the AliITSgeom creator function.
//
// Int_t *fNdet
//     A pointer to an array fNlayers long containing the number of
// active detector volumes for each ladder. This array is typically
// created and filled by the AliITSgeom creator function.
//
// ITS_geom **fg
//     A pointer to an array of pointers pointing to the ITS_geom
// structure containing the coordinate transformation information.
// The ITS_geom structure corresponding to layer=lay, ladder=lad,
// and detector=det is gotten by fg[lay-1][(fNlad[lay-1]*(lad-1)+det-1)].
// In this way a lot of space is saved over trying to keep a three
// dimensional array fNlayersXmax(fNlad)Xmax(fNdet), since the number
// of detectors typically increases with layer number.
//
// TObjArray *fShape
//     A pointer to an array of TObjects containing the detailed shape
// information for each type of detector used in the ITS. For example
// I have created AliITSgeomSPD, AliITSgeomSDD, and AliITSgeomSSD as
// example structures, derived from TObjects, to hold the detector
// information. I would recommend that one element in each of these
// structures, that which describes the shape of the active volume,
// be one of the ROOT classes derived from TShape. In this way it would
// be easy to have the display program display the correct active
// ITS volumes. See the example classes AliITSgeomSPD, AliITSgeomSDD,
// and AliITSgeomSSD for a more detailed example.
//
// Member Functions:
//
// AliITSgeom()
//     The default constructor for the AliITSgeom class. It, by default,
// sets fNlayers to zero and zeros all pointers.
//
// AliITSgeom(const char *filename)
//     The constructor for the AliITSgeom class. All of the data to fill
// this structure is read in from the file given my the input filename.
//
// AliITSgeom(AliITSgeom &source)
//     The copy constructor for the AliITSgeom class. It calls the
// = operator function. See the = operator function for more details.
//
// void operator=(AliITSgeom &source)
//     The = operator function for the AliITSgeom class. It makes an
// independent copy of the class in such a way that any changes made
// to the copied class will not affect the source class in any way.
// This is required for many ITS alignment studies where the copied
// class is then modified by introducing some misalignment.
//
// ~AliITSgeom()
//     The destructor for the AliITSgeom class. If the arrays fNlad,
// fNdet, or fg have had memory allocated to them, there pointer values
// are non zero, then this memory space is freed and they are set
// to zero. In addition, fNlayers is set to zero. The destruction of
// TObjArray fShape is, by default, handled by the TObjArray destructor.
//
// Int_t GetNdetectors(Int_t layer)
//     This function returns the number of detectors/ladder for a give 
// layer. In particular it returns fNdet[layer-1].
//
// Int_t GetNladders(Int_t layer)
//     This function returns the number of ladders for a give layer. In
// particular it returns fNlad[layer-1].
//
// Int_t GetNlayers()
//     This function returns the number of layers defined in the ITS
// geometry. In particular it returns fNlayers.
//
// GetAngles(Int_t layer,Int_t ladder,Int_t detector,
//           Float_t &rx, Float_t &ry, Float_t &rz)
//     This function returns the rotation angles for a give detector on
// a give ladder in a give layer in the three floating point variables
// provided. rx = frx, fy = fry, rz = frz. The angles are in radians
//
// GetTrans(Int_t layer,Int_t ladder,Int_t detector,
//          Float_t &x, Float_t &y, Float_t &z)
//     This function returns the Cartesian translation for a give
// detector on a give ladder in a give layer in the three floating
// point variables provided. x = fx0, y = fy0, z = fz0. The units are
// those of the Monte Carlo, generally cm.
//
// SetByAngles(Int_t layer,Int_t ladder,Int_t detector,
//            Float_t &rx, Float_t &ry, Float_t &rz)
//     This function computes a new rotation matrix based on the angles
// rx, ry, and rz (in radians) for a give detector on the give ladder
// in the give layer. A new
// fg[layer-1][(fNlad[layer-1]*(ladder-1)+detector-1)].fr[] array is
// computed.
//
// SetTrans(Int_t layer,Int_t ladder,Int_t detector,
//          Float_t x, Float_t y, Float_t z)
//     This function sets a new translation vector, given by the three
// variables x, y, and z, for the Cartesian coordinate transformation
// for the detector defined by layer, ladder and detector.
//
// GetRotMatrix(Int_t layer, Int_t ladder, Int_t detector, Float_t *mat)
//     Returns, in the Float_t array pointed to by mat, the full rotation
// matrix for the give detector defined by layer, ladder, and detector.
// It returns all nine elements of fr in the ITS_geom structure. See the
// description of the ITS_geom structure for further details of this
// rotation matrix.
//
// GtoL(Int_t layer, Int_t ladder, Int_t detector,
//       const Float_t *g, Float_t *l)
//     The function that does the global ALICE Cartesian coordinate
// to local active volume detector Cartesian coordinate transformation.
// The local detector coordinate system is determined by the layer, 
// ladder, and detector numbers. The global coordinates are entered by
// the three element Float_t array g and the local coordinate values
// are returned by the three element Float_t array l. The order of the 
// three elements are g[0]=x, g[1]=y, and g[2]=z, similarly for l.
//
// GtoL(const Int_t *Id, const Float_t *g, Float_t *l)
//     The function that does the global ALICE Cartesian coordinate
// to local active volume detector Cartesian coordinate transformation.
// The local detector coordinate system is determined by the three
// element array Id containing as it's three elements Id[0]=layer, 
// Id[1]=ladder, and Id[2]=detector numbers. The global coordinates
// are entered by the three element Float_t array g and the local
// coordinate values are returned by the three element Float_t array l.
// The order of the three elements are g[0]=x, g[1]=y, and g[2]=z,
// similarly for l.
//
//  LtoG(Int_t layer, Int_t ladder, Int_t detector,
//       const Float_t *l, Float_t *g)
//     The function that does the local active volume detector Cartesian
// coordinate to global ALICE Cartesian coordinate transformation.
// The local detector coordinate system is determined by the layer, 
// ladder, and detector numbers. The local coordinates are entered by
// the three element Float_t array l and the global coordinate values
// are returned by the three element Float_t array g. The order of the 
// three elements are l[0]=x, l[1]=y, and l[2]=z, similarly for g.
//
// LtoG(const Int_t *Id, const Float_t *l, Float_t *g)
//     The function that does the local active volume detector Cartesian
// coordinate to global ALICE Cartesian coordinate transformation.
// The local detector coordinate system is determined by the three
// element array Id containing as it's three elements Id[0]=layer, 
// Id[1]=ladder, and Id[2]=detector numbers. The local coordinates
// are entered by the three element Float_t array l and the global
// coordinate values are returned by the three element Float_t array g.
// The order of the three elements are l[0]=x, l[1]=y, and l[2]=z,
// similarly for g.
//
// Int_t IsVersion()
//     This function returns the version number of this AliITSgeom
// class.
//
// AddShape(TObject *shape)
//     This function adds one more shape element to the TObjArray
// fShape. It is primarily used in the constructor functions of the
// AliITSgeom class. The pointer *shape can be the pointer to any
// class that is derived from TObject (this is true for nearly every
// ROOT class). This does not appear to be working properly at this time.
//
// PrintComparison(FILE *fp, AliITSgeom *other)
//     This function was primarily created for diagnostic reasons. It
// print to a file pointed to by the file pointer fp the difference
// between two AliITSgeom classes. The format of the file is basicly,
// define d? to be the difference between the same element of the two
// classes. For example dfrx = this->fg[i][j].frx - other->fg[i][j].frx.
// if(at least one of dfx0, dfy0, dfz0,dfrx,dfry,dfrz are non zero) then print
// layer ladder detector dfx0 dfy0 dfz0 dfrx dfry dfrz
// if(at least one of the 9 elements of dfr[] are non zero) then print
// layer ladder detector dfr[0] dfr[1] dfr[2]
//                       dfr[3] dfr[4] dfr[5]
//                       dfr[6] dfr[7] dfr[8]
// Only non zero values are printed to save space. The differences are
// typical written to a file because there are usually a lot of numbers
// printed out and it is usually easier to read them in some nice editor
// rather than zooming quickly past you on a screen. fprintf is used to
// do the printing. The fShapeIndex difference is not printed at this time.
//
// PrintData(FILE *fp, Int_t layer, Int_t ladder, Int_t detector)
//     This function prints out the coordinate transformations for
// the particular detector defined by layer, ladder, and detector
// to the file pointed to by the File pointer fp. fprinf statements
// are used to print out the numbers. The format is
// layer ladder detector Trans= fx0 fy0 fz0 rot= frx fry frz Shape=fShapeIndex
//                         dfr= fr[0] fr[1] fr[2]
//                         dfr= fr[3] fr[4] fr[5]
//                         dfr= fr[6] fr[7] fr[8]
// By indicating which detector, some control over the information 
// is given to the user. The output it written to the file pointed
// to by the file pointer fp. This can be set to stdout if you want.
//
// Streamer(TBuffer &R__b)
//     The default Streamer function "written by ROOT" doesn't write out
// the arrays referenced by pointers. Therefore, a specific Streamer function
// has to be written. This function should not be modified but instead added
// on to so that older versions can still be read. The proper handling of
// the version dependent streamer function hasn't been written do to the lack
// of finding an example at the time of writting.
//
//----------------------------------------------------------------------
//
//     The following member functions are defined to modify an existing
// AliITSgeom data structure. They were developed for the use in doing
// alignment studies of the ITS.
//
// GlobalChange(Float_t *dtranslation, Float_t *drotation)
//     This function performs a Cartesian translation and rotation of
// the full ITS from its default position by an amount determined by
// the three element arrays dtranslation and drotation. If every element
// of dtranslation and drotation are zero then there is no change made
// the geometry. The change is global in that the exact same translation
// and rotation is done to every detector element in the exact same way.
// The units of the translation are those of the Monte Carlo, usually cm,
// and those of the rotation are in radians. The elements of dtranslation
// are dtranslation[0] = x, dtranslation[1] = y, and dtranslation[2] = z.
// The elements of drotation are drotation[0] = rx, drotation[1] = ry, and
// drotation[2] = rz. A change in x will move the hole ITS in the ALICE
// global x direction, the same for a change in y. A change in z will
// result in a translation of the ITS as a hole up or down the beam line.
// A change in the angles will result in the inclination of the ITS with
// respect to the beam line, except for an effective rotation about the
// beam axis which will just rotate the ITS as a hole about the beam axis.
//
// GlobalCylindericalChange(Float_t *dtranslation, Float_t *drotation)
//     This function performs a cylindrical translation and rotation of
// each ITS element by a fixed about in radius, rphi, and z from its
// default position by an amount determined by the three element arrays
// dtranslation and drotation. If every element of dtranslation and
// drotation are zero then there is no change made the geometry. The
// change is global in that the exact same distance change in translation
// and rotation is done to every detector element in the exact same way.
// The units of the translation are those of the Monte Carlo, usually cm,
// and those of the rotation are in radians. The elements of dtranslation
// are dtranslation[0] = r, dtranslation[1] = rphi, and dtranslation[2] = z.
// The elements of drotation are drotation[0] = rx, drotation[1] = ry, and
// drotation[2] = rz. A change in r will results in the increase of the
// radius of each layer by the same about. A change in rphi will results in
// the rotation of each layer by a different angle but by the same
// circumferential distance. A change in z will result in a translation
// of the ITS as a hole up or down the beam line. A change in the angles
// will result in the inclination of the ITS with respect to the beam
// line, except for an effective rotation about the beam axis which will
// just rotate the ITS as a hole about the beam axis.
//
// RandomChange(Float_t *stranslation, Float_t *srotation)
//     This function performs a Gaussian random displacement and/or
// rotation about the present global position of each active
// volume/detector of the ITS. The sigma of the random displacement
// is determined by the three element array stranslation, for the
// x y and z translations, and the three element array srotation,
// for the three rotation about the axis x y and z.
//
// RandomCylindericalChange(Float_t *stranslation, Float_t *srotation)
//     This function performs a Gaussian random displacement and/or
// rotation about the present global position of each active
// volume/detector of the ITS. The sigma of the random displacement
// is determined by the three element array stranslation, for the
// r rphi and z translations, and the three element array srotation,
// for the three rotation about the axis x y and z. This random change
// in detector position allow for the simulation of a random uncertainty
// in the detector positions of the ITS.
////////////////////////////////////////////////////////////////////////
 private:
    Int_t     fNlayers; // The number of layers.
    Int_t     *fNlad;   // Array of the number of ladders/layer(layer)
    Int_t     *fNdet;   // Array of the number of detectors/ladder(layer)
    ITS_geom  **fg;     // Structure of translation and rotation.
    TObjArray *fShape;  // Array of shapes and detector information.

 public:
    AliITSgeom();                      // Default constructor
    AliITSgeom(const char *filename);  // Constructor
    AliITSgeom(AliITSgeom &source);    // Copy constructor
    void operator=(AliITSgeom &source);// = operator
    virtual ~AliITSgeom();             // Default destructor
    // this is a dummy routine for now.
    inline Int_t GetNdetectors(Int_t layer) {return fNdet[layer-1];}
    inline Int_t GetNladders(Int_t layer)   {return fNlad[layer-1];}
    inline Int_t GetNlayers()               {return fNlayers;}
    inline void GetAngles(Int_t lay,Int_t lad,Int_t det,
			  Float_t &rx,Float_t &ry,Float_t &rz){
                          rx = fg[lay-1][fNdet[lay-1]*(lad-1)+det-1].frx;
                          ry = fg[lay-1][fNdet[lay-1]*(lad-1)+det-1].fry;
                          rz = fg[lay-1][fNdet[lay-1]*(lad-1)+det-1].frz;}
    inline void GetTrans(Int_t lay,Int_t lad,Int_t det,
			 Float_t &x,Float_t &y,Float_t &z){
                         x = fg[lay-1][fNdet[lay-1]*(lad-1)+det-1].fx0;
                         y = fg[lay-1][fNdet[lay-1]*(lad-1)+det-1].fy0;
                         z = fg[lay-1][fNdet[lay-1]*(lad-1)+det-1].fz0;}
    void SetByAngles(Int_t lay,Int_t lad,Int_t det,
		     Float_t rx,Float_t ry,Float_t rz);
    inline void SetTrans(Int_t lay,Int_t lad,Int_t det,
			 Float_t x,Float_t y,Float_t z){
                         fg[lay-1][fNdet[lay-1]*(lad-1)+det-1].fx0 = x;
                         fg[lay-1][fNdet[lay-1]*(lad-1)+det-1].fy0 = y;
                         fg[lay-1][fNdet[lay-1]*(lad-1)+det-1].fz0 = z;}
    void GetRotMatrix(Int_t lay,Int_t lad,Int_t det,Float_t *mat);
    void GtoL(Int_t lay,Int_t lad,Int_t det,const Float_t *g,Float_t *l);
    void GtoL(const Int_t *id,const Float_t *g,Float_t *l);
    void GtoL(const Int_t index,const Float_t *g,Float_t *l);
    void GtoLMomentum(Int_t lay,Int_t lad,Int_t det,const Float_t *g,Float_t *l);
    void LtoG(Int_t lay,Int_t lad,Int_t det,const Float_t *l,Float_t *g);
    void LtoG(const Int_t *id,const Float_t *l,Float_t *g);
    void LtoG(const Int_t index,const Float_t *l,Float_t *g);
    void LtoGMomentum(Int_t lay,Int_t lad,Int_t det,const Float_t *l,Float_t *g);
    Int_t GetModuleIndex(Int_t lay,Int_t lad,Int_t det);
    void  GetModuleId(Int_t index,Int_t &lay,Int_t &lad,Int_t &det);
    void GlobalChange(Float_t  *tran,Float_t  *rot);
    void GlobalCylindericalChange(Float_t *tran,Float_t *rot);
    void RandomChange(Float_t *stran,Float_t *srot);
    void RandomCylindericalChange(Float_t *stran,Float_t *srot);
    void PrintComparison(FILE *fp,AliITSgeom *other);
    void PrintData(FILE *fp,Int_t lay,Int_t lad,Int_t det);
    ofstream &PrintGeom(ofstream &out);
    ifstream &ReadGeom(ifstream &in);
    virtual Int_t IsVersion() const {return 0;}
    inline void AddShape(TObject *shp){fShape->AddLast(shp);}

  ClassDef(AliITSgeom,1)
};

#endif
