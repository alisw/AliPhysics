/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/*
$Log$
Revision 1.4.4.5  2000/03/04 23:42:39  nilsen
Updated the comments/documentations and improved the maintainability of the
code.

Revision 1.4.4.4  2000/03/02 21:27:07  nilsen
Added two functions, SetByAngles and SetTrans.

Revision 1.4.4.3  2000/01/23 03:09:10  nilsen
// fixed compiler warnings for new function LtLErrorMatrix(...)

Revision 1.4.4.2  2000/01/19 23:18:20  nilsen
Added transformations of Error matrix to AliITSgeom and fixed some typos
in AliITS.h and AliITShitIndex.h

Revision 1.4.4.1  2000/01/12 19:03:32  nilsen
This is the version of the files after the merging done in December 1999.
See the ReadMe110100.txt file for details

Revision 1.4  1999/10/15 07:03:20  fca
Fixed bug in GetModuleId(Int_t index,Int_t &lay,Int_t &lad, Int_t &det) and
a typo in the creator. aliroot need to be rerun to get a fixed geometry.

Revision 1.3  1999/10/04 15:20:12  fca
Correct syntax accepted by g++ but not standard for static members, remove minor warnings

Revision 1.2  1999/09/29 09:24:20  fca
Introduction of the Copyright and cvs Log

*/

///////////////////////////////////////////////////////////////////////
// ITS geometry manipulation routines.                               //
// Created April 15 1999.                                            //
// version: 0.0.0                                                    //
// By: Bjorn S. Nilsen                                               //
// version: 0.0.1                                                    //
// Updated May 27 1999.                                              //
// Added Cylindrical random and global based changes.               //
// Added  function PrintComparison.                                  //
///////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////
// The structure AliITSgeomS:
//     The structure AliITSgeomS has been defined to hold all of the
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
//
//     The local coordinate system by, default, is show in the following
// figures. Also shown are the ladder numbering scheme.
//Begin_Html
/*
<img src="picts/ITS/its1+2_convention_front_5.gif">
</pre>
<br clear=left>
<font size=+2 color=blue>
<p>This shows the front view of the SPDs and the orientation of the local
pixel coordinate system. Note that the inner pixel layer has its y coordinate
in the opposite direction from all of the other layers.
</font>
<pre>

<pre>
<img src="picts/ITS/its3+4_convention_front_5.gif">
</pre>
<br clear=left>
<font size=+2 color=blue>
<p>This shows the front view of the SDDs and the orientation of the local
pixel coordinate system.
</font>
<pre>

<pre>
<img src="picts/ITS/its5+6_convention_front_5.gif">
</pre>
<br clear=left>
<font size=+2 color=blue>
<p>This shows the front view of the SSDs and the orientation of the local
pixel coordinate system.
</font>
<pre>
*/
//End_Html

////////////////////////////////////////////////////////////////////////

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
// AliITSgeomS **fGm
//     A pointer to an array of pointers pointing to the AliITSgeomS
// structure containing the coordinate transformation information.
// The AliITSgeomS structure corresponding to layer=lay, ladder=lad,
// and detector=det is gotten by fGm[lay-1][(fNlad[lay-1]*(lad-1)+det-1)].
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
// Inlined Member Functions:
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
// SetTrans(Int_t layer,Int_t ladder,Int_t detector,
//          Float_t x, Float_t y, Float_t z)
//     This function sets a new translation vector, given by the three
// variables x, y, and z, for the Cartesian coordinate transformation
// for the detector defined by layer, ladder and detector.
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
// Int_t GetStartSPD()
//     This functions returns the starting module index number for the
// silicon pixels detectors (SPD). Typically this is zero. To loop over all
// of the pixel detectors do: for(Int_t i=GetStartSPD();i<=GetLastSPD();i++)
//
// Int_t GetLastSPD()
//     This functions returns the last module index number for the
// silicon pixels detectors (SPD). To loop over all of the pixel detectors 
// do: for(Int_t i=GetStartSPD();i<=GetLastSPD();i++)
//
// Int_t GetStartSDD()
//     This functions returns the starting module index number for the
// silicon drift detectors (SDD). To loop over all of the drift detectors 
// do: for(Int_t i=GetStartSDD();i<=GetLastSDD();i++)
//
// Int_t GetLastSDD()
//     This functions returns the last module index number for the
// silicon drift detectors (SDD). To loop over all of the drift detectors 
// do: for(Int_t i=GetStartSDD();i<=GetLastSDD();i++)
//
// Int_t GetStartSSD()
//     This functions returns the starting module index number for the
// silicon strip detectors (SSD). To loop over all of the strip detectors 
// do: for(Int_t i=GetStartSSD();i<=GetLastSSD();i++)
//
// Int_t GetStartSSD()
//     This functions returns the last module index number for the
// silicon strip detectors (SSD). To loop over all of the strip detectors 
// do: for(Int_t i=GetStartSSD();i<=GetLastSSD();i++)
//
// TObject *GetShape(Int_t lay,Int_t lad,Int_t det)
//     This functions returns the shape object AliITSgeomSPD, AliITSgeomSDD,
// or AliITSgeomSSD for that particular module designated by lay, lad, and
// detector. In principle there can be additional shape objects. In this
// way a minimum of shape objects are created since one AliITSgeomS?D shape
// object is used for all modules of that type.
////////////////////////////////////////////////////////////////////////

#include <iostream.h>
#include <fstream.h>
#include <iomanip.h>
#include <stdio.h>
#include "AliITSgeom.h"
#include "AliITSgeomSPD300.h"
#include "AliITSgeomSPD425.h"
#include "AliITSgeomSDD.h"
#include "AliITSgeomSSD.h"
#include "TRandom.h"

ClassImp(AliITSgeom)

//_____________________________________________________________________
AliITSgeom::AliITSgeom(){
////////////////////////////////////////////////////////////////////////
//     The default constructor for the AliITSgeom class. It, by default,
// sets fNlayers to zero and zeros all pointers.
////////////////////////////////////////////////////////////////////////
  // Default constructor.
  // Do not allocate anything zero everything
   fNlayers = 0;
   fNlad    = 0;
   fNdet    = 0;
   fGm       = 0;
   fShape   = 0;
   return;
}

//_____________________________________________________________________
AliITSgeom::~AliITSgeom(){
////////////////////////////////////////////////////////////////////////
//     The destructor for the AliITSgeom class. If the arrays fNlad,
// fNdet, or fGm have had memory allocated to them, there pointer values
// are non zero, then this memory space is freed and they are set
// to zero. In addition, fNlayers is set to zero. The destruction of
// TObjArray fShape is, by default, handled by the TObjArray destructor.
////////////////////////////////////////////////////////////////////////
  // Default destructor.
  // if arrays exist delete them. Then set everything to zero.
   if(fGm!=0){
      for(Int_t i=0;i<fNlayers;i++) delete[] fGm[i];
      delete[] fGm;
   } // end if fGm!=0
   if(fNlad!=0) delete[] fNlad;
   if(fNdet!=0) delete[] fNdet;
   fNlayers = 0;
   fNlad    = 0;
   fNdet    = 0;
   fGm       = 0;
   return;
}

//_____________________________________________________________________
AliITSgeom::AliITSgeom(const char *filename){
////////////////////////////////////////////////////////////////////////
//     The constructor for the AliITSgeom class. All of the data to fill
// this structure is read in from the file given my the input filename.
////////////////////////////////////////////////////////////////////////
   FILE     *pf;
   Int_t    i;
   AliITSgeomS *g;
   Int_t    l,a,d;
   Float_t  x,y,z,o,p,q,r,s,t;
   Double_t oor,pr,qr,rr,sr,tr; // Radians
   Double_t lr[9];
   Double_t si; // sin(angle)
   Double_t pi = TMath::Pi(), byPI = pi/180.;

   pf = fopen(filename,"r");

   fNlayers = 6; // set default number of ladders
   fNlad    = new Int_t[fNlayers];
   fNdet    = new Int_t[fNlayers];
   // find the number of ladders and detectors in this geometry.
   for(i=0;i<fNlayers;i++){fNlad[i]=fNdet[i]=0;} // zero out arrays
   for(;;){ // for ever loop
      i = fscanf(pf,"%d %d %d %f %f %f %f %f %f %f %f %f",
                     &l,&a,&d,&x,&y,&z,&o,&p,&q,&r,&s,&t);
      if(i==EOF) break;
      if(l<1 || l>fNlayers) {
         printf("error in file %s layer=%d min is 1 max is %d/n",
                 filename,l,fNlayers);
         continue;
      }// end if l
      if(fNlad[l-1]<a) fNlad[l-1] = a;
      if(fNdet[l-1]<d) fNdet[l-1] = d;
   } // end for ever loop
   // counted the number of ladders and detectors now allocate space.
   fGm = new AliITSgeomS* [fNlayers];
   for(i=0;i<fNlayers;i++){
      fGm[i] = 0;
      l = fNlad[i]*fNdet[i];
      fGm[i] = new AliITSgeomS[l]; // allocate space for transforms
   } // end for i

   // Set up Shapes for a default configuration of 6 layers.
   fShape = new TObjArray(3);
   AddShape((TObject *) new AliITSgeomSPD300());  // shape 0
   AddShape((TObject *) new AliITSgeomSDD());  // shape 1
   AddShape((TObject *) new AliITSgeomSSD());  // shape 2

   // prepare to read in transforms
   rewind(pf); // start over reading file
   for(;;){ // for ever loop
      i = fscanf(pf,"%d %d %d %f %f %f %f %f %f %f %f %f",
                     &l,&a,&d,&x,&y,&z,&o,&p,&q,&r,&s,&t);
      if(i==EOF) break;
      if(l<1 || l>fNlayers) {
         printf("error in file %s layer=%d min is 1 max is %d/n",
                 filename,l,fNlayers);
         continue;
      }// end if l
      l--; a--; d--; // shift layer, ladder, and detector counters to zero base
      i = d + a*fNdet[l]; // position of this detector
      g = &(fGm[l][i]);

      oor = byPI*o;
      pr = byPI*p;
      qr = byPI*q;
      rr = byPI*r;
      sr = byPI*s;
      tr = byPI*t;

      g->fx0   = x;
      g->fy0   = y;
      g->fz0   = z;
//
      si    = sin(oor);if(o== 90.0) si = +1.0;
                      if(o==270.0) si = -1.0;
                      if(o==  0.0||o==180.) si = 0.0;
      lr[0] = si * cos(pr);
      lr[1] = si * sin(pr);
      lr[2] = cos(oor);if(o== 90.0||o==270.) lr[2] = 0.0;
                      if(o== 0.0)           lr[2] = +1.0;
                      if(o==180.0)          lr[2] = -1.0;
//
      si    =  sin(qr);if(q== 90.0) si = +1.0; 
                       if(q==270.0) si = -1.0;
                       if(q==  0.0||q==180.) si = 0.0;
      lr[3] = si * cos(rr);
      lr[4] = si * sin(rr);
      lr[5] = cos(qr);if(q== 90.0||q==270.) lr[5] = 0.0;
                      if(q==  0.0)          lr[5] = +1.0;
                      if(q==180.0)          lr[5] = -1.0;
//
      si    = sin(sr);if(s== 90.0) si = +1.0;
                      if(s==270.0) si = -1.0;
                      if(s==  0.0||s==180.) si = 0.0;
      lr[6] = si * cos(tr);
      lr[7] = si * sin(tr);
      lr[8] = cos(sr);if(s== 90.0||s==270.0) lr[8] =  0.0;
                      if(s==  0.0)           lr[8] = +1.0;
                      if(s==180.0)           lr[8] = -1.0;
      // Normalize these elements
      for(a=0;a<3;a++){// reuse float Si and integers a and d.
         si = 0.0;
         for(d=0;d<3;d++) si += lr[3*a+d]*lr[3*a+d];
         si = TMath::Sqrt(1./si);
         for(d=0;d<3;d++) g->fr[3*a+d] = lr[3*a+d] = si*lr[3*a+d];
      } // end for a
      // get angles from matrix up to a phase of 180 degrees.
      oor     = atan2(lr[7],lr[8]);if(oor<0.0) oor += 2.0*pi;
      pr     = asin(lr[2]);       if(pr<0.0) pr += 2.0*pi;
      qr     = atan2(lr[3],lr[0]);if(qr<0.0) qr += 2.0*pi;
      g->frx = oor;
      g->fry = pr;
      g->frz = qr;
      // l = layer-1 at this point.
           if(l==0||l==1) g->fShapeIndex = 0; // SPD's
      else if(l==2||l==3) g->fShapeIndex = 1; // SDD's
      else if(l==4||l==5) g->fShapeIndex = 2; // SSD's
   } // end for ever loop
   fclose(pf);
}

//________________________________________________________________________
AliITSgeom::AliITSgeom(const AliITSgeom &source){
////////////////////////////////////////////////////////////////////////
//     The copy constructor for the AliITSgeom class. It calls the
// = operator function. See the = operator function for more details.
////////////////////////////////////////////////////////////////////////

    *this = source;  // Just use the = operator for now.

    return;
}

//________________________________________________________________________
/*void AliITSgeom::operator=(const AliITSgeom &source){
////////////////////////////////////////////////////////////////////////
//     The = operator function for the AliITSgeom class. It makes an
// independent copy of the class in such a way that any changes made
// to the copied class will not affect the source class in any way.
// This is required for many ITS alignment studies where the copied
// class is then modified by introducing some misalignment.
////////////////////////////////////////////////////////////////////////
   Int_t i,j,k;

   if(this == &source) return; // don't assign to ones self.

   // if there is an old structure allocated delete it first.
   if(fGm != 0){
      for(i=0;i<fNlayers;i++) delete[] fGm[i];
      delete[] fGm;
   } // end if fGm != 0 
   if(fNlad != 0) delete[] fNlad;
   if(fNdet != 0) delete[] fNdet;

   fNlayers = source.fNlayers;
   fNlad = new Int_t[fNlayers];
   for(i=0;i<fNlayers;i++) fNlad[i] = source.fNlad[i];
   fNdet = new Int_t[fNlayers];
   for(i=0;i<fNlayers;i++) fNdet[i] = source.fNdet[i];
   fShape = new TObjArray(*(source.fShape));//This does not make a proper copy.
   fGm = new AliITSgeomS* [fNlayers];
   for(i=0;i<fNlayers;i++){
      fGm[i] = new AliITSgeomS[fNlad[i]*fNdet[i]];
      for(j=0;j<(fNlad[i]*fNdet[i]);j++){
	  fGm[i][j].fShapeIndex = source.fGm[i][j].fShapeIndex;
	  fGm[i][j].fx0 = source.fGm[i][j].fx0;
	  fGm[i][j].fy0 = source.fGm[i][j].fy0;
	  fGm[i][j].fz0 = source.fGm[i][j].fz0;
	  fGm[i][j].frx = source.fGm[i][j].frx;
	  fGm[i][j].fry = source.fGm[i][j].fry;
	  fGm[i][j].frz = source.fGm[i][j].frz;
	  for(k=0;k<9;k++) fGm[i][j].fr[k] = source.fGm[i][j].fr[k];
      } // end for j
   } // end for i
   return;
   }*/
//________________________________________________________________________
AliITSgeom& AliITSgeom::operator=(const AliITSgeom &source){
////////////////////////////////////////////////////////////////////////
//     The = operator function for the AliITSgeom class. It makes an
// independent copy of the class in such a way that any changes made
// to the copied class will not affect the source class in any way.
// This is required for many ITS alignment studies where the copied
// class is then modified by introducing some misalignment.
////////////////////////////////////////////////////////////////////////
   Int_t i,j,k;

   if(this == &source) return *this; // don't assign to ones self.

   // if there is an old structure allocated delete it first.
   if(fGm != 0){
      for(i=0;i<fNlayers;i++) delete[] fGm[i];
      delete[] fGm;
   } // end if fGm != 0 
   if(fNlad != 0) delete[] fNlad;
   if(fNdet != 0) delete[] fNdet;

   fNlayers = source.fNlayers;
   fNlad = new Int_t[fNlayers];
   for(i=0;i<fNlayers;i++) fNlad[i] = source.fNlad[i];
   fNdet = new Int_t[fNlayers];
   for(i=0;i<fNlayers;i++) fNdet[i] = source.fNdet[i];
   fShape = new TObjArray(*(source.fShape));//This does not make a proper copy.
   fGm = new AliITSgeomS* [fNlayers];
   for(i=0;i<fNlayers;i++){
      fGm[i] = new AliITSgeomS[fNlad[i]*fNdet[i]];
      for(j=0;j<(fNlad[i]*fNdet[i]);j++){
	  fGm[i][j].fShapeIndex = source.fGm[i][j].fShapeIndex;
	  fGm[i][j].fx0 = source.fGm[i][j].fx0;
	  fGm[i][j].fy0 = source.fGm[i][j].fy0;
	  fGm[i][j].fz0 = source.fGm[i][j].fz0;
	  fGm[i][j].frx = source.fGm[i][j].frx;
	  fGm[i][j].fry = source.fGm[i][j].fry;
	  fGm[i][j].frz = source.fGm[i][j].frz;
	  for(k=0;k<9;k++) fGm[i][j].fr[k] = source.fGm[i][j].fr[k];
      } // end for j
   } // end for i
   return *this;
}
//________________________________________________________________________
void AliITSgeom::GtoL(Int_t lay,Int_t lad,Int_t det,
                       const Double_t *g,Double_t *l){
////////////////////////////////////////////////////////////////////////
//     The function that does the global ALICE Cartesian coordinate
// to local active volume detector Cartesian coordinate transformation.
// The local detector coordinate system is determined by the layer, 
// ladder, and detector numbers. The global coordinates are entered by
// the three element Double_t array g and the local coordinate values
// are returned by the three element Double_t array l. The order of the 
// three elements are g[0]=x, g[1]=y, and g[2]=z, similarly for l.
////////////////////////////////////////////////////////////////////////
   Double_t x,y,z;
   AliITSgeomS *gl;

   lay--; lad--; det--;
   gl = &(fGm[lay][fNdet[lay]*lad+det]);

   x    = g[0] - gl->fx0;
   y    = g[1] - gl->fy0;
   z    = g[2] - gl->fz0;
   l[0] = gl->fr[0]*x + gl->fr[1]*y + gl->fr[2]*z;
   l[1] = gl->fr[3]*x + gl->fr[4]*y + gl->fr[5]*z;
   l[2] = gl->fr[6]*x + gl->fr[7]*y + gl->fr[8]*z;
   return;
}
//________________________________________________________________________
void AliITSgeom::GtoL(const Int_t *id,const Double_t *g,Double_t *l){
////////////////////////////////////////////////////////////////////////
//     The function that does the local active volume detector Cartesian
// coordinate to global ALICE Cartesian coordinate transformation.
// The local detector coordinate system is determined by the id[0]=layer, 
// id[1]=ladder, and id[2]=detector numbers. The local coordinates are
// entered by the three element Double_t array l and the global coordinate
// values are returned by the three element Double_t array g. The order of the 
// three elements are l[0]=x, l[1]=y, and l[2]=z, similarly for g.
////////////////////////////////////////////////////////////////////////
    GtoL(id[0],id[1],id[2],g,l);
    return;
}
//________________________________________________________________________
void AliITSgeom::GtoL(const Int_t index,const Double_t *g,Double_t *l){
////////////////////////////////////////////////////////////////////////
//     The function that does the local active volume detector Cartesian
// coordinate to global ALICE Cartesian coordinate transformation.
// The local detector coordinate system is determined by the detector
// index numbers (see GetModuleIndex and GetModuleID). The local 
// coordinates are entered by the three element Double_t array l and the 
// global coordinate values are returned by the three element Double_t array g.
// The order of the three elements are l[0]=x, l[1]=y, and l[2]=z, similarly 
// for g.
////////////////////////////////////////////////////////////////////////
    Int_t    lay,lad,det;

    this->GetModuleId(index,lay,lad,det);

    GtoL(lay,lad,det,g,l);
    return;
}
//________________________________________________________________________
void AliITSgeom::GtoL(Int_t lay,Int_t lad,Int_t det,
                       const Float_t *g,Float_t *l){
////////////////////////////////////////////////////////////////////////
//     The function that does the global ALICE Cartesian coordinate
// to local active volume detector Cartesian coordinate transformation.
// The local detector coordinate system is determined by the layer, 
// ladder, and detector numbers. The global coordinates are entered by
// the three element Float_t array g and the local coordinate values
// are returned by the three element Float_t array l. The order of the 
// three elements are g[0]=x, g[1]=y, and g[2]=z, similarly for l.
////////////////////////////////////////////////////////////////////////
    Int_t    i;
    Double_t gd[3],ld[3];

    for(i=0;i<3;i++) gd[i] = (Double_t) g[i];
    GtoL(lay,lad,det,(Double_t *)gd,(Double_t *)ld);
    for(i=0;i<3;i++) l[i] = (Float_t) ld[i];
    return;
}
//________________________________________________________________________
void AliITSgeom::GtoL(const Int_t *id,const Float_t *g,Float_t *l){
////////////////////////////////////////////////////////////////////////
//     The function that does the local active volume detector Cartesian
// coordinate to global ALICE Cartesian coordinate transformation.
// The local detector coordinate system is determined by the Int_t array id,
// id[0]=layer, id[1]=ladder, and id[2]=detector numbers. The local 
// coordinates are entered by the three element Float_t array l and the
// global coordinate values are returned by the three element Float_t array g.
// The order of the three elements are l[0]=x, l[1]=y, and l[2]=z, similarly
// for g. The order of the three elements are g[0]=x, g[1]=y, and g[2]=z,
// similarly for l.
////////////////////////////////////////////////////////////////////////
    Int_t    i;
    Double_t gd[3],ld[3];

    for(i=0;i<3;i++) gd[i] = (Double_t) g[i];
    GtoL(id[0],id[1],id[2],(Double_t *)gd,(Double_t *)ld);
    for(i=0;i<3;i++) l[i] = (Float_t) ld[i];
    return;
}
//________________________________________________________________________
void AliITSgeom::GtoL(const Int_t index,const Float_t *g,Float_t *l){
////////////////////////////////////////////////////////////////////////
//     The function that does the local active volume detector Cartesian
// coordinate to global ALICE Cartesian coordinate transformation.
// The local detector coordinate system is determined by the detector
// index numbers (see GetModuleIndex and GetModuleID). The local 
// coordinates are entered by the three element Float_t array l and the 
// global coordinate values are returned by the three element Float_t array g.
// The order of the three elements are l[0]=x, l[1]=y, and l[2]=z, similarly 
// for g.
////////////////////////////////////////////////////////////////////////
    Int_t    lay,lad,det;
    Int_t    i;
    Double_t gd[3],ld[3];

    this->GetModuleId(index,lay,lad,det);

    for(i=0;i<3;i++) gd[i] = (Double_t) g[i];
    GtoL(lay,lad,det,(Double_t *)gd,(Double_t *)ld);
    for(i=0;i<3;i++) l[i] = (Float_t) ld[i];
    return;
}
//________________________________________________________________________
void AliITSgeom::LtoG(Int_t lay,Int_t lad,Int_t det,
		      const Double_t *l,Double_t *g){
////////////////////////////////////////////////////////////////////////
//     The function that does the local active volume detector Cartesian
// coordinate to global ALICE Cartesian coordinate transformation.
// The local detector coordinate system is determined by the layer, 
// ladder, and detector numbers. The local coordinates are entered by
// the three element Float_t array l and the global coordinate values
// are returned by the three element Float_t array g. The order of the 
// three elements are l[0]=x, l[1]=y, and l[2]=z, similarly for g.
////////////////////////////////////////////////////////////////////////
   Double_t x,y,z;
   AliITSgeomS *gl;

   lay--; lad--; det--;
   gl   = &(fGm[lay][fNdet[lay]*lad+det]);

   x    = gl->fr[0]*l[0] + gl->fr[3]*l[1] + gl->fr[6]*l[2];
   y    = gl->fr[1]*l[0] + gl->fr[4]*l[1] + gl->fr[7]*l[2];
   z    = gl->fr[2]*l[0] + gl->fr[5]*l[1] + gl->fr[8]*l[2];
   g[0] = x + gl->fx0;
   g[1] = y + gl->fy0;
   g[2] = z + gl->fz0;
   return;
}
//________________________________________________________________________
void AliITSgeom::LtoG(const Int_t *id,const Double_t *l,Double_t *g){
////////////////////////////////////////////////////////////////////////
//     The function that does the local active volume detector Cartesian
// coordinate to global ALICE Cartesian coordinate transformation.
// The local detector coordinate system is determined by the three
// element array Id containing as it's three elements Id[0]=layer, 
// Id[1]=ladder, and Id[2]=detector numbers. The local coordinates
// are entered by the three element Double_t array l and the global
// coordinate values are returned by the three element Double_t array g.
// The order of the three elements are l[0]=x, l[1]=y, and l[2]=z,
// similarly for g.
////////////////////////////////////////////////////////////////////////
    LtoG(id[0],id[1],id[2],l,g);
    return;
}
//________________________________________________________________________
void AliITSgeom::LtoG(const Int_t index,const Double_t *l,Double_t *g){
////////////////////////////////////////////////////////////////////////
//     The function that does the local active volume detector Cartesian
// coordinate to global ALICE Cartesian coordinate transformation.
// The local detector coordinate system is determined by the detector  
// index number (see GetModuleIndex and GetModuleId). The local coordinates
// are entered by the three element Double_t array l and the global
// coordinate values are returned by the three element Double_t array g.
// The order of the three elements are l[0]=x, l[1]=y, and l[2]=z,
// similarly for g.
////////////////////////////////////////////////////////////////////////
    Int_t    lay,lad,det;

    this->GetModuleId(index,lay,lad,det);

    LtoG(lay,lad,det,l,g);
    return;
}
//________________________________________________________________________
void AliITSgeom::LtoG(Int_t lay,Int_t lad,Int_t det,
		      const Float_t *l,Float_t *g){
////////////////////////////////////////////////////////////////////////
//     The function that does the local active volume detector Cartesian
// coordinate to global ALICE Cartesian coordinate transformation.
// The local detector coordinate system is determined by the layer, 
// ladder, and detector numbers. The local coordinates are entered by
// the three element Float_t array l and the global coordinate values
// are returned by the three element Float_t array g. The order of the 
// three elements are l[0]=x, l[1]=y, and l[2]=z, similarly for g.
////////////////////////////////////////////////////////////////////////
    Int_t    i;
    Double_t gd[3],ld[3];

    for(i=0;i<3;i++) ld[i] = (Double_t) l[i];
    LtoG(lay,lad,det,(Double_t *)ld,(Double_t *)gd);
    for(i=0;i<3;i++) g[i] = (Float_t) gd[i];
    return;
}
//________________________________________________________________________
void AliITSgeom::LtoG(const Int_t *id,const Float_t *l,Float_t *g){
////////////////////////////////////////////////////////////////////////
//     The function that does the local active volume detector Cartesian
// coordinate to global ALICE Cartesian coordinate transformation.
// The local detector coordinate system is determined by the three
// element array Id containing as it's three elements Id[0]=layer, 
// Id[1]=ladder, and Id[2]=detector numbers. The local coordinates
// are entered by the three element Float_t array l and the global
// coordinate values are returned by the three element Float_t array g.
// The order of the three elements are l[0]=x, l[1]=y, and l[2]=z,
// similarly for g.
////////////////////////////////////////////////////////////////////////
    Int_t    i;
    Double_t gd[3],ld[3];

    for(i=0;i<3;i++) ld[i] = (Double_t) l[i];
    LtoG(id[0],id[1],id[2],(Double_t *)ld,(Double_t *)gd);
    for(i=0;i<3;i++) g[i] = (Float_t) gd[i];
    return;
}
//________________________________________________________________________
void AliITSgeom::LtoG(const Int_t index,const Float_t *l,Float_t *g){
////////////////////////////////////////////////////////////////////////
//     The function that does the local active volume detector Cartesian
// coordinate to global ALICE Cartesian coordinate transformation.
// The local detector coordinate system is determined by the detector  
// index number (see GetModuleIndex and GetModuleId). The local coordinates
// are entered by the three element Float_t array l and the global
// coordinate values are returned by the three element Float_t array g.
// The order of the three elements are l[0]=x, l[1]=y, and l[2]=z,
// similarly for g.
////////////////////////////////////////////////////////////////////////
    Int_t    i,lay,lad,det;
    Double_t gd[3],ld[3];

    this->GetModuleId(index,lay,lad,det);

    for(i=0;i<3;i++) ld[i] = (Double_t) l[i];
    LtoG(lay,lad,det,(Double_t *)ld,(Double_t *)gd);
    for(i=0;i<3;i++) g[i] = (Float_t) gd[i];
    return;
}
//______________________________________________________________________
void AliITSgeom::LtoL(const Int_t *id1,const Int_t *id2,
		      Double_t *l1,Double_t *l2){
////////////////////////////////////////////////////////////////////////
//     The function that does the local active volume detector Cartesian
// coordinate to a different local active volume detector Cartesian coordinate
// transformation. The original local detector coordinate system is determined
// by the detector array id1, id1[0]=layer, id1[1]=ladder, and id1[2]=detector
// and the new coordinate system is determined by the detector array id2,
// id2[0]=layer, id2[1]=ladder, and id2[2]=detector. The original local
// coordinates are entered by the three element Double_t array l1 and the
// other new local coordinate values are returned by the three element
// Double_t array l2. The order of the three elements are l1[0]=x, l1[1]=y,
// and l1[2]=z, similarly for l2.
////////////////////////////////////////////////////////////////////////
    Double_t g[3];

    LtoG(id1,l1,g);
    GtoL(id2,g,l2);
    return;
}
//______________________________________________________________________
void AliITSgeom::LtoL(const Int_t index1,const Int_t index2,
		      Double_t *l1,Double_t *l2){
////////////////////////////////////////////////////////////////////////
//     The function that does the local active volume detector Cartesian
// coordinate to a different local active volume detector Cartesian coordinate
// transformation. The original local detector coordinate system is determined
// by the detector index number index1, and the new coordinate system is
// determined by the detector index number index2, (see GetModuleIndex and
// GetModuleId). The original local coordinates are entered by the three
// element Double_t array l1 and the other new local coordinate values are
// returned by the three element Double_t array l2. The order of the three
// elements are l1[0]=x, l1[1]=y, and l1[2]=z, similarly for l2.
////////////////////////////////////////////////////////////////////////
    Double_t g[3];

    LtoG(index1,l1,g);
    GtoL(index2,g,l2);
    return;
}
//________________________________________________________________________
void AliITSgeom::GtoLMomentum(Int_t lay,Int_t lad,Int_t det,
			      const Double_t *g,Double_t *l){
////////////////////////////////////////////////////////////////////////
//     The function that does the global ALICE Cartesian momentum
// to local active volume detector Cartesian momentum transformation.
// The local detector coordinate system is determined by the layer, 
// ladder, and detector numbers. The global momentums are entered by
// the three element Double_t array g and the local momentums values
// are returned by the three element Double_t array l. The order of the 
// three elements are g[0]=x, g[1]=y, and g[2]=z, similarly for l.
////////////////////////////////////////////////////////////////////////
   Double_t px,py,pz;
   AliITSgeomS *gl;

   lay--; lad--; det--;
   gl = &(fGm[lay][fNdet[lay]*lad+det]);

   px   = g[0];
   py   = g[1];
   pz   = g[2];
   l[0] = gl->fr[0]*px + gl->fr[1]*py + gl->fr[2]*pz;
   l[1] = gl->fr[3]*px + gl->fr[4]*py + gl->fr[5]*pz;
   l[2] = gl->fr[6]*px + gl->fr[7]*py + gl->fr[8]*pz;
   return;
}
//________________________________________________________________________
void AliITSgeom::GtoLMomentum(Int_t lay,Int_t lad,Int_t det,
			      const Float_t *g,Float_t *l){
////////////////////////////////////////////////////////////////////////
//     The function that does the global ALICE Cartesian momentum
// to local active volume detector Cartesian momentum transformation.
// The local detector coordinate system is determined by the layer, 
// ladder, and detector numbers. The global momentums are entered by
// the three element Float_t array g and the local momentums values
// are returned by the three element Float_t array l. The order of the 
// three elements are g[0]=x, g[1]=y, and g[2]=z, similarly for l.
////////////////////////////////////////////////////////////////////////
    Int_t i;
    Double_t gd[3],ld[3];

    for(i=0;i<3;i++) gd[i] = (Double_t) g[i];
    GtoLMomentum(lay,lad,det,(Double_t *)gd,(Double_t *)ld);
    for(i=0;i<3;i++) l[i] = (Float_t) ld[i];
    return;
}
//________________________________________________________________________
void AliITSgeom::LtoGMomentum(Int_t lay,Int_t lad,Int_t det,
			      const Double_t *l,Double_t *g){
////////////////////////////////////////////////////////////////////////
//     The function that does the local active volume detector Cartesian
// momentum to global ALICE Cartesian momentum transformation.
// The local detector momentum system is determined by the layer, 
// ladder, and detector numbers. The local momentums are entered by
// the three element Double_t array l and the global momentum values
// are returned by the three element Double_t array g. The order of the 
// three elements are l[0]=x, l[1]=y, and l[2]=z, similarly for g.
////////////////////////////////////////////////////////////////////////
   Double_t px,py,pz;
   AliITSgeomS *gl;

   lay--; lad--; det--;
   gl   = &(fGm[lay][fNdet[lay]*lad+det]);

   px   = gl->fr[0]*l[0] + gl->fr[3]*l[1] + gl->fr[6]*l[2];
   py   = gl->fr[1]*l[0] + gl->fr[4]*l[1] + gl->fr[7]*l[2];
   pz   = gl->fr[2]*l[0] + gl->fr[5]*l[1] + gl->fr[8]*l[2];
   g[0] = px;
   g[1] = py;
   g[2] = pz;
   return;
}
//________________________________________________________________________
void AliITSgeom::LtoGMomentum(Int_t lay,Int_t lad,Int_t det,
			      const Float_t *l,Float_t *g){
////////////////////////////////////////////////////////////////////////
//     The function that does the local active volume detector Cartesian
// momentum to global ALICE Cartesian momentum transformation.
// The local detector momentum system is determined by the layer, 
// ladder, and detector numbers. The local momentums are entered by
// the three element Float_t array l and the global momentum values
// are returned by the three element Float_t array g. The order of the 
// three elements are l[0]=x, l[1]=y, and l[2]=z, similarly for g.
////////////////////////////////////////////////////////////////////////
    Int_t i;
    Double_t gd[3],ld[3];

    for(i=0;i<3;i++) ld[i] = (Double_t) l[i];
    LtoGMomentum(lay,lad,det,(Double_t *)ld,(Double_t *)gd);
    for(i=0;i<3;i++) g[i] = (Float_t) gd[i];
    return;
}
//______________________________________________________________________
void AliITSgeom::LtoLMomentum(const Int_t *id1,const Int_t *id2,
			      const Double_t *l1,Double_t *l2){
////////////////////////////////////////////////////////////////////////
//     The function that does the local active volume detector Cartesian
// momentum to a different local active volume detector Cartesian momentum
// transformation. The original local detector momentum system is determined
// by the Int_t array id1 (id1[0]=lay, id1[1]=lad, id1[2]=det). The new local
// coordinate system id determined by the Int_t array id2. The local
// momentums are entered by the three element Double_t array l1 and the other
// local momentum values are returned by the three element Double_t array l2.
// The order of the three elements are l1[0]=x, l1[1]=y, and l1[2]=z,
// similarly for l2.
////////////////////////////////////////////////////////////////////////
    Double_t g[3];

    LtoGMomentum(id1[0],id1[1],id1[2],l1,g);
    GtoLMomentum(id2[0],id2[1],id2[2],g,l2);
    return;
}
//______________________________________________________________________
void AliITSgeom::GtoLErrorMatrix(const Int_t index,Double_t **g,Double_t **l){
////////////////////////////////////////////////////////////////////////
//      This converts an error matrix, expressed in global coordinates
// into an error matrix expressed in local coordinates. Since the 
// translations do not change the error matrix they are not included.
// Definition: if GtoL is l[i] = T[i][j]*g[j], then from the definition
// of the transformation matrix above T[i][j] = fr[3*i+j]. Then for a 
// matrix l[i][l] = T[i][j]*g[j][k]*T[l][k] (sum over repeated indexes). 
// Where T[l][k] is the transpose of T[k][l].
////////////////////////////////////////////////////////////////////////
    Double_t lR[3][3],lRt[3][3];
    Int_t    lay,lad,det,i,j,k,n;
    AliITSgeomS *gl;

    GetModuleId(index,lay,lad,det);
    lay--;lad--;det--;
    gl = &(fGm[lay][fNdet[lay]*lad+det]);

    for(i=0;i<3;i++)for(j=0;j<3;j++){
	lR[i][j] = lRt[j][i] = gl->fr[3*i+j];
    } // end for i,j

    for(i=0;i<3;i++)for(j=0;j<3;j++)for(k=0;k<3;k++)for(n=0;n<3;n++){
	l[i][n] = lR[i][j]*g[j][k]*lRt[k][n];
    } // end for i,j,k,l
    return;
}
//______________________________________________________________________
void AliITSgeom::LtoGErrorMatrix(const Int_t index,Double_t **l,Double_t **g){
////////////////////////////////////////////////////////////////////////
//      This converts an error matrix, expressed in local coordinates
// into an error matrix expressed in global coordinates. Since the 
// translations do not change the error matrix they are not included.
// Definition: if GtoL is l[i] = T[i][j]*g[j], then from the definition
// of the transformation matrix above T[i][j] = fr[3*i+j]. Then for a 
// matrix g[i][l] = T[j][i]*l[j][k]*T[k][l] (sum over repeated indexes). 
// Where T[j][i] is the transpose of T[i][j].
////////////////////////////////////////////////////////////////////////
    Double_t lR[3][3],lRt[3][3];
    Int_t    lay,lad,det,i,j,k,n;
    AliITSgeomS *gl;

    GetModuleId(index,lay,lad,det);
    lay--;lad--;det--;
    gl = &(fGm[lay][fNdet[lay]*lad+det]);

    for(i=0;i<3;i++)for(j=0;j<3;j++){
	lR[i][j] = lRt[j][i] = gl->fr[3*i+j];
    } // end for i,j

    for(i=0;i<3;i++)for(j=0;j<3;j++)for(k=0;k<3;k++)for(n=0;n<3;n++){
	g[i][n] = lRt[i][j]*l[j][k]*lR[k][n];
    } // end for i,j,k,l
    return;
}
//______________________________________________________________________
void AliITSgeom::LtoLErrorMatrix(const Int_t index1,const Int_t index2,
				 Double_t **l1,Double_t **l2){
////////////////////////////////////////////////////////////////////////
//      This converts an error matrix, expressed in one local coordinates
// into an error matrix expressed in different local coordinates. Since  
// the translations do not change the error matrix they are not included.
// This is done by going through the global coordinate system for 
// simplicity and constancy.
////////////////////////////////////////////////////////////////////////
    Double_t g[3][3];

    this->LtoGErrorMatrix(index1,l1,(Double_t **)g);
    this->GtoLErrorMatrix(index2,(Double_t **)g,l2);
    return;
}
//______________________________________________________________________
Int_t AliITSgeom::GetModuleIndex(Int_t lay,Int_t lad,Int_t det){
////////////////////////////////////////////////////////////////////////
//      This routine computes the module index number from the layer,
// ladder, and detector numbers. The number of ladders and detectors
// per layer is determined when this geometry package is constructed,
// see AliITSgeom(const char *filename) for specifics.
////////////////////////////////////////////////////////////////////////
    Int_t i,j,k;

    i = fNdet[lay-1] * (lad-1) + det - 1;
    j = 0;
    for(k=0;k<lay-1;k++) j += fNdet[k]*fNlad[k];
    return (i+j);
}
//___________________________________________________________________________
void AliITSgeom::GetModuleId(Int_t index,Int_t &lay,Int_t &lad,Int_t &det){
////////////////////////////////////////////////////////////////////////
//      This routine computes the layer, ladder and detector number 
// given the module index number. The number of ladders and detectors
// per layer is determined when this geometry package is constructed,
// see AliITSgeom(const char *filename) for specifics.
////////////////////////////////////////////////////////////////////////
    Int_t i,j,k;

    j = 0;
    for(k=0;k<fNlayers;k++){
	j += fNdet[k]*fNlad[k];
	if(j>index)break;
    } // end for k
    lay = k+1;
    i = index -j + fNdet[k]*fNlad[k];
    j = 0;
    for(k=0;k<fNlad[lay-1];k++){
	j += fNdet[lay-1];
	if(j>i)break;
    } // end for k
    lad = k+1;
    det = 1+i-fNdet[lay-1]*k;
    return;
}
//___________________________________________________________________________
void AliITSgeom::GetRotMatrix(Int_t lay,Int_t lad,Int_t det,Double_t *mat){
////////////////////////////////////////////////////////////////////////
//     Returns, in the Double_t array pointed to by mat, the full rotation
// matrix for the give detector defined by layer, ladder, and detector.
// It returns all nine elements of fr in the AliITSgeomS structure. See the
// description of the AliITSgeomS structure for further details of this
// rotation matrix.
////////////////////////////////////////////////////////////////////////
   Int_t    i;
   AliITSgeomS *g;

   lay--; lad--; det--; // shift to base 0
   g = &(fGm[lay][fNdet[lay]*lad+det]);
   for(i=0;i<9;i++) mat[i] = g->fr[i];
   return;
}
//___________________________________________________________________________
void AliITSgeom::GetRotMatrix(Int_t index,Double_t *mat){
////////////////////////////////////////////////////////////////////////
//     Returns, in the Double_t array pointed to by mat, the full rotation
// matrix for the give detector defined by the module index number.
// It returns all nine elements of fr in the AliITSgeomS structure. See the
// description of the AliITSgeomS structure for further details of this
// rotation matrix.
////////////////////////////////////////////////////////////////////////
   Int_t    lay,lad,det;

   this->GetModuleId(index,lay,lad,det);
   GetRotMatrix(lay,lad,det,mat);
   return;
}
//___________________________________________________________________________
void AliITSgeom::GetRotMatrix(Int_t lay,Int_t lad,Int_t det,Float_t *mat){
////////////////////////////////////////////////////////////////////////
//     Returns, in the Float_t array pointed to by mat, the full rotation
// matrix for the give detector defined by layer, ladder, and detector.
// It returns all nine elements of fr in the AliITSgeomS structure. See the
// description of the AliITSgeomS structure for further details of this
// rotation matrix.
////////////////////////////////////////////////////////////////////////
   Int_t    i;
   Double_t matd[9];

   GetRotMatrix(lay,lad,det,(Double_t *)matd);
   for(i=0;i<9;i++) mat[i] = (Float_t) matd[i];
   return;
}

//___________________________________________________________________________
void AliITSgeom::GetRotMatrix(Int_t index,Float_t *mat){
////////////////////////////////////////////////////////////////////////
//     Returns, in the Float_t array pointed to by mat, the full rotation
// matrix for the give detector defined by module index number.
// It returns all nine elements of fr in the AliITSgeomS structure. See the
// description of the AliITSgeomS structure for further details of this
// rotation matrix.
////////////////////////////////////////////////////////////////////////
   Int_t    i,lay,lad,det;
   Double_t matd[9];

   this->GetModuleId(index,lay,lad,det);
   GetRotMatrix(lay,lad,det,(Double_t *)matd);
   for(i=0;i<9;i++) mat[i] = (Float_t) matd[i];
   return;
}

//___________________________________________________________________________
Int_t AliITSgeom::GetStartDet(Int_t id){
  /////////////////////////////////////////////////////////////////////////
  // returns the starting module index value for a give type of detector id
  /////////////////////////////////////////////////////////////////////////
  Int_t first;
  switch(id)
  {
  case 0:
     first = GetModuleIndex(1,1,1);
     break;
  case 1:
     first = GetModuleIndex(3,1,1);
     break;
  case 2:
     first = GetModuleIndex(5,1,1);
     break;
  default:
     printf("<AliITSgeom::GetFirstDet> undefined detector type\n");
     first = 0;

  }
  return first;
}

//___________________________________________________________________________
Int_t AliITSgeom::GetLastDet(Int_t id){
  /////////////////////////////////////////////////////////////////////////
  // returns the last module index value for a give type of detector id
  /////////////////////////////////////////////////////////////////////////
  Int_t last;
  switch(id)
  {
  case 0:
     last = GetLastSPD();
     break;
   case 1:
     last = GetLastSDD();
     break;
   case 2:
     last = GetLastSSD();
     break;
   default:
     printf("<AliITSgeom::GetLastDet> undefined detector type\n");
     last = 0;
  }
  return last;
}

//___________________________________________________________________________
void AliITSgeom::PrintComparison(FILE *fp,AliITSgeom *other){
////////////////////////////////////////////////////////////////////////
//     This function was primarily created for diagnostic reasons. It
// print to a file pointed to by the file pointer fp the difference
// between two AliITSgeom classes. The format of the file is basicly,
// define d? to be the difference between the same element of the two
// classes. For example dfrx = this->fGm[i][j].frx - other->fGm[i][j].frx.
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
////////////////////////////////////////////////////////////////////////
   Int_t    i,j,k,l;
   Double_t xt,yt,zt,xo,yo,zo;
   Double_t rxt,ryt,rzt,rxo,ryo,rzo;  // phi in radians
   AliITSgeomS *gt,*go;
   Bool_t   t;

   for(i=0;i<this->fNlayers;i++){
      for(j=0;j<this->fNlad[i];j++) for(k=0;k<this->fNdet[i];k++){
	 l   = this->fNdet[i]*j+k; // resolved index
         gt  = &(this->fGm[i][l]);
	 go  = &(other->fGm[i][l]);
         xt  = gt->fx0; yt  = gt->fy0; zt  = gt->fz0;
         xo  = go->fx0; yo  = go->fy0; zo  = go->fz0;
         rxt = gt->frx; ryt = gt->fry; rzt = gt->frz;
         rxo = go->frx; ryo = go->fry; rzo = go->frz;
	 if(!(xt==xo&&yt==yo&&zt==zo&&rxt==rxo&&ryt==ryo&&rzt==rzo))
	 fprintf(fp,"%1.1d %2.2d %2.2d dTrans=%f %f %f drot=%f %f %f\n",
		 i+1,j+1,k+1,xt-xo,yt-yo,zt-zo,rxt-rxo,ryt-ryo,rzt-rzo);
	 t = kFALSE;
	 for(i=0;i<9;i++) t = gt->fr[i] != go->fr[i];
	 if(t){
	     fprintf(fp,"%1.1d %2.2d %2.2d dfr= %e %e %e\n",i+1,j+1,k+1,
                 gt->fr[0]-go->fr[0],gt->fr[1]-go->fr[1],gt->fr[2]-go->fr[2]);
	     fprintf(fp,"        dfr= %e %e %e\n",
                 gt->fr[3]-go->fr[3],gt->fr[4]-go->fr[4],gt->fr[5]-go->fr[5]);
	     fprintf(fp,"        dfr= %e %e %e\n",
                 gt->fr[6]-go->fr[6],gt->fr[7]-go->fr[7],gt->fr[8]-go->fr[8]);
	 }
      } // end for j,k
   } // end for i
   return;
}

//___________________________________________________________________________
void AliITSgeom::PrintData(FILE *fp,Int_t lay,Int_t lad,Int_t det){
////////////////////////////////////////////////////////////////////////
//     This function prints out the coordinate transformations for
// the particular detector defined by layer, ladder, and detector
// to the file pointed to by the File pointer fp. fprintf statements
// are used to print out the numbers. The format is
// layer ladder detector Trans= fx0 fy0 fz0 rot= frx fry frz Shape=fShapeIndex
//                         dfr= fr[0] fr[1] fr[2]
//                         dfr= fr[3] fr[4] fr[5]
//                         dfr= fr[6] fr[7] fr[8]
// By indicating which detector, some control over the information 
// is given to the user. The output it written to the file pointed
// to by the file pointer fp. This can be set to stdout if you want.
////////////////////////////////////////////////////////////////////////
   Int_t    i,j,k,l;
   AliITSgeomS *gt;

   i  = lay-1;
   j  = lad-1;
   k  = det-1;
   l  = this->fNdet[i]*j+k; // resolved index
   gt = &(this->fGm[i][l]);
   fprintf(fp,"%1.1d %2.2d %2.2d Trans=%f %f %f rot=%f %f %f Shape=%d\n",
	   i+1,j+1,k+1,gt->fx0,gt->fy0,gt->fz0,gt->frx,gt->fry,gt->frz,
           gt->fShapeIndex);
   fprintf(fp,"        dfr= %e %e %e\n",gt->fr[0],gt->fr[1],gt->fr[2]);
   fprintf(fp,"        dfr= %e %e %e\n",gt->fr[3],gt->fr[4],gt->fr[5]);
   fprintf(fp,"        dfr= %e %e %e\n",gt->fr[6],gt->fr[7],gt->fr[8]);
   return;
}
//___________________________________________________________________________
ofstream & AliITSgeom::PrintGeom(ofstream &lRb){
////////////////////////////////////////////////////////////////////////
//     The default Streamer function "written by ROOT" doesn't write out
// the arrays referenced by pointers. Therefore, a specific Streamer function
// has to be written. This function should not be modified but instead added
// on to so that older versions can still be read. The proper handling of
// the version dependent streamer function hasn't been written do to the lack
// of finding an example at the time of writing.
////////////////////////////////////////////////////////////////////////
   // Stream an object of class AliITSgeom.
    Int_t i,j,k;

    lRb.setf(ios::scientific);
    lRb << fNlayers << " ";
    for(i=0;i<fNlayers;i++) lRb << fNlad[i] << " ";
    for(i=0;i<fNlayers;i++) lRb << fNdet[i] << "\n";
    for(i=0;i<fNlayers;i++) for(j=0;j<fNlad[i]*fNdet[i];j++){
	lRb <<setprecision(16) << fGm[i][j].fShapeIndex << " ";
	lRb <<setprecision(16) << fGm[i][j].fx0 << " ";
	lRb <<setprecision(16) << fGm[i][j].fy0 << " ";
	lRb <<setprecision(16) << fGm[i][j].fz0 << " ";
	lRb <<setprecision(16) << fGm[i][j].frx << " ";
	lRb <<setprecision(16) << fGm[i][j].fry << " ";
	lRb <<setprecision(16) << fGm[i][j].frz << "\n";
	for(k=0;k<9;k++) lRb <<setprecision(16) << fGm[i][j].fr[k] << " ";
	lRb << "\n";
      } // end for i,j
//      lRb << fShape;
      return lRb;
}
//___________________________________________________________________________
ifstream & AliITSgeom::ReadGeom(ifstream &lRb){
////////////////////////////////////////////////////////////////////////
//     The default Streamer function "written by ROOT" doesn't write out
// the arrays referenced by pointers. Therefore, a specific Streamer function
// has to be written. This function should not be modified but instead added
// on to so that older versions can still be read. The proper handling of
// the version dependent streamer function hasn't been written do to the lack
// of finding an example at the time of writing.
////////////////////////////////////////////////////////////////////////
   // Stream an object of class AliITSgeom.
    Int_t i,j,k;

      lRb >> fNlayers;
      if(fNlad!=0) delete[] fNlad;
      if(fNdet!=0) delete[] fNdet;
      fNlad = new Int_t[fNlayers];
      fNdet = new Int_t[fNlayers];
      for(i=0;i<fNlayers;i++) lRb >> fNlad[i];
      for(i=0;i<fNlayers;i++) lRb >> fNdet[i];
      if(fGm!=0){
	  for(i=0;i<fNlayers;i++) delete[] fGm[i];
	  delete[] fGm;
      } // end if fGm!=0
      fGm = new AliITSgeomS*[fNlayers];
      for(i=0;i<fNlayers;i++){
	  fGm[i] = new AliITSgeomS[fNlad[i]*fNdet[i]];
	  for(j=0;j<fNlad[i]*fNdet[i];j++){
	      lRb >> fGm[i][j].fShapeIndex;
	      lRb >> fGm[i][j].fx0;
	      lRb >> fGm[i][j].fy0;
	      lRb >> fGm[i][j].fz0;
	      lRb >> fGm[i][j].frx;
	      lRb >> fGm[i][j].fry;
	      lRb >> fGm[i][j].frz;
	      for(k=0;k<9;k++) lRb >> fGm[i][j].fr[k];
	  } // end for j
      } // end for i
//      lRb >> fShape;
      return lRb;
}
//______________________________________________________________________
//     The following routines modify the transformation of "this"
// geometry transformations in a number of different ways.
//______________________________________________________________________
void AliITSgeom::SetByAngles(Int_t lay,Int_t lad,Int_t det,
			     Float_t rx,Float_t ry,Float_t rz){
////////////////////////////////////////////////////////////////////////
//     This function computes a new rotation matrix based on the angles
// rx, ry, and rz (in radians) for a give detector on the give ladder
// in the give layer. A new
// fGm[layer-1][(fNlad[layer-1]*(ladder-1)+detector-1)].fr[] array is
// computed.
////////////////////////////////////////////////////////////////////////
   AliITSgeomS *g;
   Double_t  sx,cx,sy,cy,sz,cz;

   lay--; lad--; det--; // set to zero base now.
   g = &(fGm[lay][fNdet[lay]*lad+det]);

   sx = sin(rx); cx = cos(rx);
   sy = sin(ry); cy = cos(ry);
   sz = sin(rz); cz = cos(rz);
   g->frx   = rx;
   g->fry   = ry;
   g->frz   = rz;
   g->fr[0] =  cz*cy;
   g->fr[1] = -cz*sy*sx - sz*cx;
   g->fr[2] = -cz*sy*cx + sz*sx;
   g->fr[3] =  sz*cy;
   g->fr[4] = -sz*sy*sx + cz*cx;
   g->fr[5] = -sz*sy*cx - cz*sx;
   g->fr[6] =  sy;
   g->fr[7] =  cy*sx;
   g->fr[8] =  cy*cx;
   return;
}
//______________________________________________________________________
void AliITSgeom::SetByAngles(Int_t index,Double_t angl[]){
////////////////////////////////////////////////////////////////////////
//     Sets the coordinate rotation transformation for a given module
// as determined by the module index number.
////////////////////////////////////////////////////////////////////////
    Int_t lay,lad,det;
    Float_t x,y,z;

    GetModuleId(index,lay,lad,det);
    x = (Float_t) angl[0];
    y = (Float_t) angl[1];
    z = (Float_t) angl[2];
    SetByAngles(lay,lad,det,x,y,z);
    return;
}
//______________________________________________________________________
void AliITSgeom::SetTrans(Int_t index,Double_t v[]){
////////////////////////////////////////////////////////////////////////
//     Sets the coordinate translation for a given module as determined
// by the module index number.
////////////////////////////////////////////////////////////////////////
    Int_t lay,lad,det;
    Float_t x,y,z;

    GetModuleId(index,lay,lad,det);
    x = (Float_t) v[0];
    y = (Float_t) v[1];
    z = (Float_t) v[2];
    SetTrans(lay,lad,det,x,y,z);
    return;
}
//___________________________________________________________________________
void AliITSgeom::GlobalChange(Float_t *tran,Float_t *rot){
////////////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////////
   Int_t    i,j,k,l;
   Double_t rx,ry,rz;
   Double_t sx,cx,sy,cy,sz,cz;
   AliITSgeomS *gl;

   for(i=0;i<fNlayers;i++){
      for(j=0;j<fNlad[i];j++) for(k=0;k<fNdet[i];k++){
	 l = fNdet[i]*j+k; // resolved index
         gl = &(fGm[i][l]);
         gl->fx0 += tran[0];
         gl->fy0 += tran[1];
         gl->fz0 += tran[2];
         gl->frx +=  rot[0];
         gl->fry +=  rot[1];
         gl->frz +=  rot[2];
         rx = gl->frx; ry = gl->fry; rz = gl->frz;
         sx = sin(rx); cx = cos(rx);
         sy = sin(ry); cy = cos(ry);
         sz = sin(rz); cz = cos(rz);
         gl->fr[0] =  cz*cy;
         gl->fr[1] = -cz*sy*sx - sz*cx;
         gl->fr[2] = -cz*sy*cx + sz*sx;
         gl->fr[3] =  sz*cy;
         gl->fr[4] = -sz*sy*sx + cz*cx;
         gl->fr[5] = -sz*sy*cx - cz*sx;
         gl->fr[6] =  sy;
         gl->fr[7] =  cy*sx;
         gl->fr[8] =  cy*cx;
      } // end for j,k
   } // end for i
   return;
}

//___________________________________________________________________________
void AliITSgeom::GlobalCylindericalChange(Float_t *tran,Float_t *rot){
////////////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////////
   Int_t    i,j,k,l;
   Double_t rx,ry,rz,r,phi,rphi; // phi in radians
   Double_t sx,cx,sy,cy,sz,cz,r0;
   AliITSgeomS *gl;

   for(i=0;i<fNlayers;i++){
      for(j=0;j<fNlad[i];j++) for(k=0;k<fNdet[i];k++){
	 l     = fNdet[i]*j+k; // resolved index
         gl    = &(fGm[i][l]);
	 r = r0= TMath::Hypot(gl->fy0,gl->fx0);
	 phi   = atan2(gl->fy0,gl->fx0);
	 rphi  = r0*phi;
	 r    += tran[0];
	 rphi += tran[1];
	 phi   = rphi/r0;
         gl->fx0  = r*TMath::Cos(phi);
         gl->fy0  = r*TMath::Sin(phi);
         gl->fz0 += tran[2];
         gl->frx +=  rot[0];
         gl->fry +=  rot[1];
         gl->frz +=  rot[2];
         rx = gl->frx; ry = gl->fry; rz = gl->frz;
         sx = sin(rx); cx = cos(rx);
         sy = sin(ry); cy = cos(ry);
         sz = sin(rz); cz = cos(rz);
         gl->fr[0] =  cz*cy;
         gl->fr[1] = -cz*sy*sx - sz*cx;
         gl->fr[2] = -cz*sy*cx + sz*sx;
         gl->fr[3] =  sz*cy;
         gl->fr[4] = -sz*sy*sx + cz*cx;
         gl->fr[5] = -sz*sy*cx - cz*sx;
         gl->fr[6] =  sy;
         gl->fr[7] =  cy*sx;
         gl->fr[8] =  cy*cx;
      } // end for j,k
   } // end for i
   return;
}

//___________________________________________________________________________
void AliITSgeom::RandomChange(Float_t *stran,Float_t *srot){
////////////////////////////////////////////////////////////////////////
//     This function performs a Gaussian random displacement and/or
// rotation about the present global position of each active
// volume/detector of the ITS. The sigma of the random displacement
// is determined by the three element array stran, for the
// x y and z translations, and the three element array srot,
// for the three rotation about the axis x y and z.
////////////////////////////////////////////////////////////////////////
   Int_t    i,j,k,l;
   Double_t rx,ry,rz;
   Double_t sx,cx,sy,cy,sz,cz;
   TRandom  ran;
   AliITSgeomS *gl;

   for(i=0;i<fNlayers;i++){
      for(j=0;j<fNlad[i];j++) for(k=0;k<fNdet[i];k++){
	 l = fNdet[i]*j+k; // resolved index
         gl = &(fGm[i][l]);
         gl->fx0 += ran.Gaus(0.0,stran[0]);
         gl->fy0 += ran.Gaus(0.0,stran[1]);
         gl->fz0 += ran.Gaus(0.0,stran[2]);
         gl->frx += ran.Gaus(0.0, srot[0]);
         gl->fry += ran.Gaus(0.0, srot[1]);
         gl->frz += ran.Gaus(0.0, srot[2]);
         rx = gl->frx; ry = gl->fry; rz = gl->frz;
         sx = sin(rx); cx = cos(rx);
         sy = sin(ry); cy = cos(ry);
         sz = sin(rz); cz = cos(rz);
         gl->fr[0] =  cz*cy;
         gl->fr[1] = -cz*sy*sx - sz*cx;
         gl->fr[2] = -cz*sy*cx + sz*sx;
         gl->fr[3] =  sz*cy;
         gl->fr[4] = -sz*sy*sx + cz*cx;
         gl->fr[5] = -sz*sy*cx - cz*sx;
         gl->fr[6] =  sy;
         gl->fr[7] =  cy*sx;
         gl->fr[8] =  cy*cx;
      } // end for j,k
   } // end for i
   return;
}

//___________________________________________________________________________
void AliITSgeom::RandomCylindericalChange(Float_t *stran,Float_t *srot){
////////////////////////////////////////////////////////////////////////
//     This function performs a Gaussian random displacement and/or
// rotation about the present global position of each active
// volume/detector of the ITS. The sigma of the random displacement
// is determined by the three element array stran, for the
// r rphi and z translations, and the three element array srot,
// for the three rotation about the axis x y and z. This random change
// in detector position allow for the simulation of a random uncertainty
// in the detector positions of the ITS.
////////////////////////////////////////////////////////////////////////
   Int_t     i,j,k,l;
   Double_t  rx,ry,rz,r,phi,x,y;  // phi in radians
   Double_t  sx,cx,sy,cy,sz,cz,r0;
   TRandom   ran;
   AliITSgeomS  *gl;

   for(i=0;i<fNlayers;i++){
      for(j=0;j<fNlad[i];j++) for(k=0;k<fNdet[i];k++){
	 l     = fNdet[i]*j+k; // resolved index
         gl    = &(fGm[i][l]);
	 x     = gl->fx0;
	 y     = gl->fy0;
	 r = r0= TMath::Hypot(y,x);
	 phi   = TMath::ATan2(y,x);
	 r    += ran.Gaus(0.0,stran[0]);
	 phi  += ran.Gaus(0.0,stran[1])/r0;
         gl->fx0  = r*TMath::Cos(phi);
         gl->fy0  = r*TMath::Sin(phi);
         gl->fz0 += ran.Gaus(0.0,stran[2]);
         gl->frx += ran.Gaus(0.0, srot[0]);
         gl->fry += ran.Gaus(0.0, srot[1]);
         gl->frz += ran.Gaus(0.0, srot[2]);
         rx = gl->frx; ry = gl->fry; rz = gl->frz;
         sx = sin(rx); cx = cos(rx);
         sy = sin(ry); cy = cos(ry);
         sz = sin(rz); cz = cos(rz);
         gl->fr[0] =  cz*cy;
         gl->fr[1] = -cz*sy*sx - sz*cx;
         gl->fr[2] = -cz*sy*cx + sz*sx;
         gl->fr[3] =  sz*cy;
         gl->fr[4] = -sz*sy*sx + cz*cx;
         gl->fr[5] = -sz*sy*cx - cz*sx;
         gl->fr[6] =  sy;
         gl->fr[7] =  cy*sx;
         gl->fr[8] =  cy*cx;
      } // end for j,k
   } // end for i
   return;
}
//______________________________________________________________________
void AliITSgeom::GeantToTracking(AliITSgeom &source){
/////////////////////////////////////////////////////////////////////////
//     Copy the geometry data but change it to make coordinate systems
// changes between the Global to the Local coordinate system used for 
// ITS tracking. Basicly the difference is that the direction of the
// y coordinate system for layer 1 is rotated about the z axis 180 degrees
// so that it points in the same direction as it does in all of the other
// layers.
// Fixed for bug and new calulation of tracking coordiantes. BSN June 8 2000.
////////////////////////////////////////////////////////////////////////////
   Double_t oor,pr,qr;
   Int_t    i,j,k;
   Double_t pi = TMath::Pi();

   if(this == &source) return; // don't assign to ones self.

   // if there is an old structure allocated delete it first.
   if(fGm != 0){
      for(i=0;i<fNlayers;i++) delete[] fGm[i];
      delete[] fGm;
   } // end if fGm != 0 
   if(fNlad != 0) delete[] fNlad;
   if(fNdet != 0) delete[] fNdet;

   fNlayers = source.fNlayers;
   fNlad = new Int_t[fNlayers];
   for(i=0;i<fNlayers;i++) fNlad[i] = source.fNlad[i];
   fNdet = new Int_t[fNlayers];
   for(i=0;i<fNlayers;i++) fNdet[i] = source.fNdet[i];
   fShape = new TObjArray(*(source.fShape));//This does not make a proper copy.
   fGm = new AliITSgeomS* [fNlayers];
   for(i=0;i<fNlayers;i++){
      fGm[i] = new AliITSgeomS[fNlad[i]*fNdet[i]];
      for(j=0;j<(fNlad[i]*fNdet[i]);j++){
	  fGm[i][j].fShapeIndex = source.fGm[i][j].fShapeIndex;
	  fGm[i][j].fx0 = source.fGm[i][j].fx0;
	  fGm[i][j].fy0 = source.fGm[i][j].fy0;
	  fGm[i][j].fz0 = source.fGm[i][j].fz0;
	  fGm[i][j].frx = source.fGm[i][j].frx;
	  fGm[i][j].fry = source.fGm[i][j].fry;
	  fGm[i][j].frz = source.fGm[i][j].frz;
	  for(k=0;k<9;k++) fGm[i][j].fr[k] = source.fGm[i][j].fr[k];
	  if(i==0) { // layer=1 is placed up side down
	      // mupliply by -1  0 0
	      //              0 -1 0
	      //              0  0 1.
	      fGm[i][j].fr[0] = -source.fGm[i][j].fr[0];
	      fGm[i][j].fr[1] = -source.fGm[i][j].fr[1];
	      fGm[i][j].fr[2] = -source.fGm[i][j].fr[2];
	      fGm[i][j].fr[3] = -source.fGm[i][j].fr[3];
              fGm[i][j].fr[4] = -source.fGm[i][j].fr[4];
              fGm[i][j].fr[5] = -source.fGm[i][j].fr[5];
	  } // end if i=1
	  // get angles from matrix up to a phase of 180 degrees.
	  oor     = atan2(fGm[i][j].fr[7],fGm[i][j].fr[8]);
	  if(oor<0.0) oor += 2.0*pi;
	  pr     = asin(fGm[i][j].fr[2]);
	  if(pr<0.0) pr += 2.0*pi;
	  qr     = atan2(fGm[i][j].fr[3],fGm[i][j].fr[0]);
	  if(qr<0.0) qr += 2.0*pi;
	  fGm[i][j].frx = oor;
	  fGm[i][j].fry = pr;
	  fGm[i][j].frz = qr;
      } // end for j
   } // end for i
   return;
}
//___________________________________________________________________________
void AliITSgeom::Streamer(TBuffer &lRb){
////////////////////////////////////////////////////////////////////////
//     The default Streamer function "written by ROOT" doesn't write out
// the arrays referenced by pointers. Therefore, a specific Streamer function
// has to be written. This function should not be modified but instead added
// on to so that older versions can still be read. The proper handling of
// the version dependent streamer function hasn't been written do to the lack
// of finding an example at the time of writing.
////////////////////////////////////////////////////////////////////////
   // Stream an object of class AliITSgeom.
    Int_t i,j,k,n;


   printf("AliITSgeomStreamer starting\n");
   if (lRb.IsReading()) {
      Version_t lRv = lRb.ReadVersion(); if (lRv) { }
      TObject::Streamer(lRb);
      printf("AliITSgeomStreamer reading fNlayers\n");
      lRb >> fNlayers;
      if(fNlad!=0) delete[] fNlad;
      if(fNdet!=0) delete[] fNdet;
      fNlad = new Int_t[fNlayers];
      fNdet = new Int_t[fNlayers];
      printf("AliITSgeomStreamer fNlad\n");
      for(i=0;i<fNlayers;i++) lRb >> fNlad[i];
      printf("AliITSgeomStreamer fNdet\n");
      for(i=0;i<fNlayers;i++) lRb >> fNdet[i];
      if(fGm!=0){
	  for(i=0;i<fNlayers;i++) delete[] fGm[i];
	  delete[] fGm;
      } // end if fGm!=0
      fGm = new AliITSgeomS*[fNlayers];
      printf("AliITSgeomStreamer AliITSgeomS\n");
      for(i=0;i<fNlayers;i++){
	  n     = fNlad[i]*fNdet[i];
	  fGm[i] = new AliITSgeomS[n];
	  for(j=0;j<n;j++){
	      lRb >> fGm[i][j].fShapeIndex;
	      lRb >> fGm[i][j].fx0;
	      lRb >> fGm[i][j].fy0;
	      lRb >> fGm[i][j].fz0;
	      lRb >> fGm[i][j].frx;
	      lRb >> fGm[i][j].fry;
	      lRb >> fGm[i][j].frz;
	      for(k=0;k<9;k++) lRb >> fGm[i][j].fr[k];
	  } // end for j
      } // end for i
      /*
      if(fShape!=0){
	  delete fShape;
      } // end if
      printf("AliITSgeomStreamer reading fShape\n");
      lRb >> fShape;
      */
      //if (fShape) fShape->Streamer(lRb);
   } else {
      lRb.WriteVersion(AliITSgeom::IsA());
      TObject::Streamer(lRb);
      lRb << fNlayers;
      for(i=0;i<fNlayers;i++) lRb << fNlad[i];
      for(i=0;i<fNlayers;i++) lRb << fNdet[i];
      for(i=0;i<fNlayers;i++) for(j=0;j<fNlad[i]*fNdet[i];j++){
	  lRb << fGm[i][j].fShapeIndex;
	  lRb << fGm[i][j].fx0;
	  lRb << fGm[i][j].fy0;
	  lRb << fGm[i][j].fz0;
	  lRb << fGm[i][j].frx;
	  lRb << fGm[i][j].fry;
	  lRb << fGm[i][j].frz;
	  for(k=0;k<9;k++) lRb << fGm[i][j].fr[k];
      } // end for i,j
      // lRb << fShape;
      //if (fShape) fShape->Streamer(lRb);
   } // end if reading
   printf("AliITSgeomStreamer Finished\n");
}
