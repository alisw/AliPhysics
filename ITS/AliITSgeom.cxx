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
Revision 1.3  1999/10/04 15:20:12  fca
Correct syntax accepted by g++ but not standard for static members, remove minor warnings

Revision 1.2  1999/09/29 09:24:20  fca
Introduction of the Copyright and cvs Log

*/

///////////////////////////////////////////////////////////////////////
// ITS geometry manimulaiton routines.                               //
// Created April 15 1999.                                            //
// version: 0.0.0                                                    //
// By: Bjorn S. Nilsen                                               //
// version: 0.0.1                                                    //
// Updated May 27 1999.                                              //
// Added Cylinderical random and global based changes.               //
// Added  function PrintComparison.                                  //
///////////////////////////////////////////////////////////////////////
#include <iostream.h>
#include <fstream.h>
#include <iomanip.h>
#include <stdio.h>
#include "AliITSgeom.h"
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
   fg       = 0;
   fShape   = 0;
   return;
}

//_____________________________________________________________________
AliITSgeom::~AliITSgeom(){
////////////////////////////////////////////////////////////////////////
//     The destructor for the AliITSgeom class. If the arrays fNlad,
// fNdet, or fg have had memory allocated to them, there pointer values
// are non zero, then this memory space is freed and they are set
// to zero. In addition, fNlayers is set to zero. The destruction of
// TObjArray fShape is, by default, handled by the TObjArray destructor.
////////////////////////////////////////////////////////////////////////
  // Default destructor.
  // if arrays exist delet them. Then set everything to zero.
   if(fg!=0){
      for(Int_t i=0;i<fNlayers;i++) delete[] fg[i];
      delete[] fg;
   } // end if fg!=0
   if(fNlad!=0) delete[] fNlad;
   if(fNdet!=0) delete[] fNdet;
   fNlayers = 0;
   fNlad    = 0;
   fNdet    = 0;
   fg       = 0;
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
   ITS_geom *g;
   Int_t    l,a,d;
   Float_t  x,y,z,o,p,q,r,s,t;
   Double_t oor,pr,qr,rr,sr,tr; // Radians
   Double_t lr[9];
   Double_t si; // sin(angle)
   Double_t PI = TMath::Pi(), byPI = PI/180.;

   pf = fopen(filename,"r");

   fNlayers = 6; // set default number of ladders
   fNlad    = new Int_t[fNlayers];
   fNdet    = new Int_t[fNlayers];
   // find the number of laders and detectors in this geometry.
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
   // counted the number of laders and detectors now allocate space.
   fg = new ITS_geom* [fNlayers];
   for(i=0;i<fNlayers;i++){
      fg[i] = 0;
      l = fNlad[i]*fNdet[i];
      fg[i] = new ITS_geom[l]; // allocate space for transforms
   } // end for i

   // Set up Shapes for a default configuration of 6 layers.
   fShape = new TObjArray;
   AddShape((TObject *) new AliITSgeomSPD());  // shape 0
   AddShape((TObject *) new AliITSgeomSDD());  // shape 1
   AddShape((TObject *) new AliITSgeomSPD());  // shape 2

   // prepair to read in transforms
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
      l--; a--; d--; // shift layer, lader, and detector counters to zero base
      i = d + a*fNdet[l]; // position of this detector
      g = &(fg[l][i]);

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
      for(a=0;a<3;a++){// reuse float si and integers a and d.
         si = 0.0;
         for(d=0;d<3;d++) si += lr[3*a+d]*lr[3*a+d];
         si = TMath::Sqrt(1./si);
         for(d=0;d<3;d++) g->fr[3*a+d] = lr[3*a+d] = si*lr[3*a+d];
      } // end for a
      // get angles from matrix up to a phase of 180 degrees.
      oor     = atan2(lr[7],lr[8]);if(oor<0.0) oor += 2.0*PI;
      pr     = asin(lr[2]);       if(pr<0.0) pr += 2.0*PI;
      qr     = atan2(lr[3],lr[0]);if(qr<0.0) qr += 2.0*PI;
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
AliITSgeom::AliITSgeom(AliITSgeom &source){
////////////////////////////////////////////////////////////////////////
//     The copy constructor for the AliITSgeom class. It calls the
// = operator function. See the = operator function for more details.
////////////////////////////////////////////////////////////////////////
   source = *this;  // Just use the = operator for now.
   return;
}

//________________________________________________________________________
void AliITSgeom::operator=(AliITSgeom &source){
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
   if(fg != 0){
      for(i=0;i<fNlayers;i++) delete[] fg[i];
      delete[] fg;
   } // end if fg != 0 
   if(fNlad != 0) delete[] fNlad;
   if(fNdet != 0) delete[] fNdet;

   fNlayers = source.fNlayers;
   fNlad = new Int_t[fNlayers];
   for(i=0;i<fNlayers;i++) fNlad[i] = source.fNlad[i];
   fNdet = new Int_t[fNlayers];
   for(i=0;i<fNlayers;i++) fNdet[i] = source.fNdet[i];
   fShape = new TObjArray(*(source.fShape));//This does not make a proper copy.
   fg = new ITS_geom* [fNlayers];
   for(i=0;i<fNlayers;i++){
      fg[i] = new ITS_geom[fNlad[i]*fNdet[i]];
      for(j=0;j<(fNlad[i]*fNdet[i]);j++){
	  fg[i][j].fShapeIndex = source.fg[i][j].fShapeIndex;
	  fg[i][j].fx0 = source.fg[i][j].fx0;
	  fg[i][j].fy0 = source.fg[i][j].fy0;
	  fg[i][j].fz0 = source.fg[i][j].fz0;
	  fg[i][j].frx = source.fg[i][j].frx;
	  fg[i][j].fry = source.fg[i][j].fry;
	  fg[i][j].frz = source.fg[i][j].frz;
	  for(k=0;k<9;k++) fg[i][j].fr[k] = source.fg[i][j].fr[k];
      } // end for j
   } // end for i
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
   Double_t x,y,z;
   ITS_geom *gl;

   lay--; lad--; det--;
   gl = &(fg[lay][fNdet[lay]*lad+det]);

   x    = g[0] - gl->fx0;
   y    = g[1] - gl->fy0;
   z    = g[2] - gl->fz0;
   l[0] = gl->fr[0]*x + gl->fr[1]*y + gl->fr[2]*z;
   l[1] = gl->fr[3]*x + gl->fr[4]*y + gl->fr[5]*z;
   l[2] = gl->fr[6]*x + gl->fr[7]*y + gl->fr[8]*z;
   return;
}

//________________________________________________________________________
void AliITSgeom::GtoL(const Int_t *id,const Float_t *g,Float_t *l){
////////////////////////////////////////////////////////////////////////
//     The function that does the local active volume detector Cartesian
// coordinate to global ALICE Cartesian coordinate transformation.
// The local detector coordinate system is determined by the layer, 
// ladder, and detector numbers. The local coordinates are entered by
// the three element Float_t array l and the global coordinate values
// are returned by the three element Float_t array g. The order of the 
// three elements are l[0]=x, l[1]=y, and l[2]=z, similarly for g.
////////////////////////////////////////////////////////////////////////
   Int_t    lay,lad,det;
   Double_t x,y,z;
   ITS_geom *gl;

   lay = id[0]; lad = id[1]; det = id[2];
   lay--; lad--; det--;
   gl = &(fg[lay][fNdet[lay]*lad+det]);

   x    = g[0] - gl->fx0;
   y    = g[1] - gl->fy0;
   z    = g[2] - gl->fz0;
   l[0] = gl->fr[0]*x + gl->fr[1]*y + gl->fr[2]*z;
   l[1] = gl->fr[3]*x + gl->fr[4]*y + gl->fr[5]*z;
   l[2] = gl->fr[6]*x + gl->fr[7]*y + gl->fr[8]*z;
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
   Double_t x,y,z;
   ITS_geom *gl;

   this->GetModuleId(index,lay,lad,det);
   lay--; lad--; det--;
   gl = &(fg[lay][fNdet[lay]*lad+det]);

   x    = g[0] - gl->fx0;
   y    = g[1] - gl->fy0;
   z    = g[2] - gl->fz0;
   l[0] = gl->fr[0]*x + gl->fr[1]*y + gl->fr[2]*z;
   l[1] = gl->fr[3]*x + gl->fr[4]*y + gl->fr[5]*z;
   l[2] = gl->fr[6]*x + gl->fr[7]*y + gl->fr[8]*z;
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
   Double_t x,y,z;
   ITS_geom *gl;

   lay--; lad--; det--;
   gl   = &(fg[lay][fNdet[lay]*lad+det]);

   x    = gl->fr[0]*l[0] + gl->fr[3]*l[1] + gl->fr[6]*l[2];
   y    = gl->fr[1]*l[0] + gl->fr[4]*l[1] + gl->fr[7]*l[2];
   z    = gl->fr[2]*l[0] + gl->fr[5]*l[1] + gl->fr[8]*l[2];
   g[0] = x + gl->fx0;
   g[1] = y + gl->fy0;
   g[2] = z + gl->fz0;
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
   Int_t    lay,lad,det;
   Double_t x,y,z;
   ITS_geom *gl;

   lay = id[0]; lad = id[1]; det = id[2];
   lay--; lad--; det--;
   gl   = &(fg[lay][fNdet[lay]*lad+det]);

   x    = gl->fr[0]*l[0] + gl->fr[3]*l[1] + gl->fr[6]*l[2];
   y    = gl->fr[1]*l[0] + gl->fr[4]*l[1] + gl->fr[7]*l[2];
   z    = gl->fr[2]*l[0] + gl->fr[5]*l[1] + gl->fr[8]*l[2];
   g[0] = x + gl->fx0;
   g[1] = y + gl->fy0;
   g[2] = z + gl->fz0;
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
   Int_t    lay,lad,det;
   Double_t x,y,z;
   ITS_geom *gl;

   this->GetModuleId(index,lay,lad,det);
   lay--; lad--; det--;
   gl   = &(fg[lay][fNdet[lay]*lad+det]);

   x    = gl->fr[0]*l[0] + gl->fr[3]*l[1] + gl->fr[6]*l[2];
   y    = gl->fr[1]*l[0] + gl->fr[4]*l[1] + gl->fr[7]*l[2];
   z    = gl->fr[2]*l[0] + gl->fr[5]*l[1] + gl->fr[8]*l[2];
   g[0] = x + gl->fx0;
   g[1] = y + gl->fy0;
   g[2] = z + gl->fz0;
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
   Double_t px,py,pz;
   ITS_geom *gl;

   lay--; lad--; det--;
   gl = &(fg[lay][fNdet[lay]*lad+det]);

   px   = g[0];
   py   = g[1];
   pz   = g[2];
   l[0] = gl->fr[0]*px + gl->fr[1]*py + gl->fr[2]*pz;
   l[1] = gl->fr[3]*px + gl->fr[4]*py + gl->fr[5]*pz;
   l[2] = gl->fr[6]*px + gl->fr[7]*py + gl->fr[8]*pz;
   return;
}
//________________________________________________________________________
void AliITSgeom::LtoGMomentum(Int_t lay,Int_t lad,Int_t det,
			      const Float_t *l,Float_t *g){
////////////////////////////////////////////////////////////////////////
//     The function that does the local active volume detector Cartesian
// momentum to global ALICE Cartesian momentum transformation.
// The local detector momentum system is determined by the layer, 
// ladder, and detector numbers. The locall momentums are entered by
// the three element Float_t array l and the global momentum values
// are returned by the three element Float_t array g. The order of the 
// three elements are l[0]=x, l[1]=y, and l[2]=z, similarly for g.
////////////////////////////////////////////////////////////////////////
   Double_t px,py,pz;
   ITS_geom *gl;

   lay--; lad--; det--;
   gl   = &(fg[lay][fNdet[lay]*lad+det]);

   px   = gl->fr[0]*l[0] + gl->fr[3]*l[1] + gl->fr[6]*l[2];
   py   = gl->fr[1]*l[0] + gl->fr[4]*l[1] + gl->fr[7]*l[2];
   pz   = gl->fr[2]*l[0] + gl->fr[5]*l[1] + gl->fr[8]*l[2];
   g[0] = px;
   g[1] = py;
   g[2] = pz;
   return;
}
//___________________________________________________________________________
Int_t AliITSgeom::GetModuleIndex(Int_t lay,Int_t lad,Int_t det){
    Int_t i,j,k;

    i = fNdet[lay-1] * (lad-1) + det - 1;
    j = 0;
    for(k=0;k<lay-1;k++) j += fNdet[k]*fNlad[k];
    return (i+j);
}
//___________________________________________________________________________
void AliITSgeom::GetModuleId(Int_t index,Int_t &lay,Int_t &lad,Int_t &det){
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
   ITS_geom *gl;

   for(i=0;i<fNlayers;i++){
      for(j=0;j<fNlad[i];j++) for(k=0;k<fNdet[i];k++){
	 l = fNdet[i]*j+k; // resolved index
         gl = &(fg[i][l]);
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
   ITS_geom *gl;

//   printf("trans=%f %f %f rot=%f %f %f\n",tran[0],tran[1],tran[2],
//	  rot[0],rot[1],rot[2]);
   for(i=0;i<fNlayers;i++){
      for(j=0;j<fNlad[i];j++) for(k=0;k<fNdet[i];k++){
	 l     = fNdet[i]*j+k; // resolved index
         gl    = &(fg[i][l]);
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
// is determined by the three element array stranslation, for the
// x y and z translations, and the three element array srotation,
// for the three rotation about the axis x y and z.
////////////////////////////////////////////////////////////////////////
   Int_t    i,j,k,l;
   Double_t rx,ry,rz;
   Double_t sx,cx,sy,cy,sz,cz;
   TRandom  ran;
   ITS_geom *gl;

   for(i=0;i<fNlayers;i++){
      for(j=0;j<fNlad[i];j++) for(k=0;k<fNdet[i];k++){
	 l = fNdet[i]*j+k; // resolved index
         gl = &(fg[i][l]);
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
// is determined by the three element array stranslation, for the
// r rphi and z translations, and the three element array srotation,
// for the three rotation about the axis x y and z. This random change
// in detector position allow for the simulation of a random uncertainty
// in the detector positions of the ITS.
////////////////////////////////////////////////////////////////////////
   Int_t     i,j,k,l;
   Double_t  rx,ry,rz,r,phi,x,y;  // phi in radians
   Double_t  sx,cx,sy,cy,sz,cz,r0;
   TRandom   ran;
   ITS_geom  *gl;

//   printf("trans=%f %f %f rot=%f %f %f\n",stran[0],stran[1],stran[2],
//	  srot[0],srot[1],srot[2]);
   for(i=0;i<fNlayers;i++){
      for(j=0;j<fNlad[i];j++) for(k=0;k<fNdet[i];k++){
	 l     = fNdet[i]*j+k; // resolved index
         gl    = &(fg[i][l]);
	 x     = gl->fx0;
	 y     = gl->fy0;
	 r = r0= TMath::Hypot(y,x);
	 phi   = TMath::ATan2(y,x);
//	 if(phi<0.0) phi += 2.0*TMath::Pi();
	 r    += ran.Gaus(0.0,stran[0]);
	 phi  += ran.Gaus(0.0,stran[1])/r0;
//	 printf("fx0=%f fy0=%f rcos(phi)=%f rsin(phi)=%f\n",gl->fx0,gl->fy0,
//		r*TMath::Cos(phi),r*TMath::Sin(phi));
         gl->fx0  = r*TMath::Cos(phi);
         gl->fy0  = r*TMath::Sin(phi);
//	 printf("r0=%f r=%f hypot=%f phi0=%f phi=%f ATan2=%f\n",
//		r0,r,TMath::Hypot(gl->fy0,gl->fx0),
//		phi0,phi,TMath::ATan2(gl->fy0,gl->fx0));
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
void AliITSgeom::SetByAngles(Int_t lay,Int_t lad,Int_t det,
			     Float_t rx,Float_t ry,Float_t rz){
////////////////////////////////////////////////////////////////////////
//     This function computes a new rotation matrix based on the angles
// rx, ry, and rz (in radians) for a give detector on the give ladder
// in the give layer. A new
// fg[layer-1][(fNlad[layer-1]*(ladder-1)+detector-1)].fr[] array is
// computed.
////////////////////////////////////////////////////////////////////////
   ITS_geom *g;
   Double_t  sx,cx,sy,cy,sz,cz;

   lay--; lad--; det--; // set to zero base now.
   g = &(fg[lay][fNdet[lay]*lad+det]);

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

//___________________________________________________________________________
void AliITSgeom::GetRotMatrix(Int_t lay,Int_t lad,Int_t det,Float_t *mat){
////////////////////////////////////////////////////////////////////////
//     Returns, in the Float_t array pointed to by mat, the full rotation
// matrix for the give detector defined by layer, ladder, and detector.
// It returns all nine elements of fr in the ITS_geom structure. See the
// description of the ITS_geom structure for further details of this
// rotation matrix.
////////////////////////////////////////////////////////////////////////
   Int_t    i;
   ITS_geom *g;

   lay--; lad--; det--; // shift to base 0
   g = &(fg[lay][fNdet[lay]*lad+det]);
   for(i=0;i<9;i++) mat[i] = g->fr[i];
   return;
}

//___________________________________________________________________________
void AliITSgeom::PrintComparison(FILE *fp,AliITSgeom *other){
////////////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////////
   Int_t    i,j,k,l;
   Double_t xt,yt,zt,xo,yo,zo;
   Double_t rxt,ryt,rzt,rxo,ryo,rzo;  // phi in radians
   ITS_geom *gt,*go;
   Bool_t   t;

   for(i=0;i<this->fNlayers;i++){
      for(j=0;j<this->fNlad[i];j++) for(k=0;k<this->fNdet[i];k++){
	 l   = this->fNdet[i]*j+k; // resolved index
         gt  = &(this->fg[i][l]);
	 go  = &(other->fg[i][l]);
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
// to the file pointed to by the File pointer fp. fprinf statements
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
   ITS_geom *gt;

   i  = lay-1;
   j  = lad-1;
   k  = det-1;
   l  = this->fNdet[i]*j+k; // resolved index
   gt = &(this->fg[i][l]);
   fprintf(fp,"%1.1d %2.2d %2.2d Trans=%f %f %f rot=%f %f %f Shape=%d\n",
	   i+1,j+1,k+1,gt->fx0,gt->fy0,gt->fz0,gt->frx,gt->fry,gt->frz,
           gt->fShapeIndex);
   fprintf(fp,"        dfr= %e %e %e\n",gt->fr[0],gt->fr[1],gt->fr[2]);
   fprintf(fp,"        dfr= %e %e %e\n",gt->fr[3],gt->fr[4],gt->fr[5]);
   fprintf(fp,"        dfr= %e %e %e\n",gt->fr[6],gt->fr[7],gt->fr[8]);
   return;
}
//___________________________________________________________________________
void AliITSgeom::Streamer(TBuffer &R__b){
////////////////////////////////////////////////////////////////////////
//     The default Streamer function "written by ROOT" doesn't write out
// the arrays referenced by pointers. Therefore, a specific Streamer function
// has to be written. This function should not be modified but instead added
// on to so that older versions can still be read. The proper handling of
// the version dependent streamer function hasn't been written do to the lack
// of finding an example at the time of writting.
////////////////////////////////////////////////////////////////////////
   // Stream an object of class AliITSgeom.
    Int_t i,j,k;

   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(); if (R__v) { }
      TObject::Streamer(R__b);
      R__b >> fNlayers;
      if(fNlad!=0) delete[] fNlad;
      if(fNdet!=0) delete[] fNdet;
      fNlad = new Int_t[fNlayers];
      fNdet = new Int_t[fNlayers];
      for(i=0;i<fNlayers;i++) R__b >> fNlad[i];
      for(i=0;i<fNlayers;i++) R__b >> fNdet[i];
      if(fg!=0){
	  for(i=0;i<fNlayers;i++) delete[] fg[i];
	  delete[] fg;
      } // end if fg!=0
      fg = new ITS_geom*[fNlayers];
      for(i=0;i<fNlayers;i++){
	  fg[i] = new ITS_geom[fNlad[i]*fNdet[i]];
	  for(j=0;j<fNlad[i]*fNdet[i];j++){
	      R__b >> fg[i][j].fShapeIndex;
	      R__b >> fg[i][j].fx0;
	      R__b >> fg[i][j].fy0;
	      R__b >> fg[i][j].fz0;
	      R__b >> fg[i][j].frx;
	      R__b >> fg[i][j].fry;
	      R__b >> fg[i][j].frz;
	      for(k=0;k<9;k++) R__b >> fg[i][j].fr[k];
	  } // end for j
      } // end for i
      R__b >> fShape;
   } else {
      R__b.WriteVersion(AliITSgeom::IsA());
      TObject::Streamer(R__b);
      R__b << fNlayers;
      for(i=0;i<fNlayers;i++) R__b << fNlad[i];
      for(i=0;i<fNlayers;i++) R__b << fNdet[i];
      for(i=0;i<fNlayers;i++) for(j=0;j<fNlad[i]*fNdet[i];j++){
	  R__b << fg[i][j].fShapeIndex;
	  R__b << fg[i][j].fx0;
	  R__b << fg[i][j].fy0;
	  R__b << fg[i][j].fz0;
	  R__b << fg[i][j].frx;
	  R__b << fg[i][j].fry;
	  R__b << fg[i][j].frz;
	  for(k=0;k<9;k++) R__b << fg[i][j].fr[k];
      } // end for i,j
      R__b << fShape;
   }
}

//___________________________________________________________________________
ofstream & AliITSgeom::PrintGeom(ofstream &R__b){
////////////////////////////////////////////////////////////////////////
//     The default Streamer function "written by ROOT" doesn't write out
// the arrays referenced by pointers. Therefore, a specific Streamer function
// has to be written. This function should not be modified but instead added
// on to so that older versions can still be read. The proper handling of
// the version dependent streamer function hasn't been written do to the lack
// of finding an example at the time of writting.
////////////////////////////////////////////////////////////////////////
   // Stream an object of class AliITSgeom.
    Int_t i,j,k;

    R__b.setf(ios::scientific);
    R__b << fNlayers << " ";
    for(i=0;i<fNlayers;i++) R__b << fNlad[i] << " ";
    for(i=0;i<fNlayers;i++) R__b << fNdet[i] << "\n";
    for(i=0;i<fNlayers;i++) for(j=0;j<fNlad[i]*fNdet[i];j++){
	R__b <<setprecision(16) << fg[i][j].fShapeIndex << " ";
	R__b <<setprecision(16) << fg[i][j].fx0 << " ";
	R__b <<setprecision(16) << fg[i][j].fy0 << " ";
	R__b <<setprecision(16) << fg[i][j].fz0 << " ";
	R__b <<setprecision(16) << fg[i][j].frx << " ";
	R__b <<setprecision(16) << fg[i][j].fry << " ";
	R__b <<setprecision(16) << fg[i][j].frz << "\n";
	for(k=0;k<9;k++) R__b <<setprecision(16) << fg[i][j].fr[k] << " ";
	R__b << "\n";
      } // end for i,j
//      R__b << fShape;
      return R__b;
}

//___________________________________________________________________________
ifstream & AliITSgeom::ReadGeom(ifstream &R__b){
////////////////////////////////////////////////////////////////////////
//     The default Streamer function "written by ROOT" doesn't write out
// the arrays referenced by pointers. Therefore, a specific Streamer function
// has to be written. This function should not be modified but instead added
// on to so that older versions can still be read. The proper handling of
// the version dependent streamer function hasn't been written do to the lack
// of finding an example at the time of writting.
////////////////////////////////////////////////////////////////////////
   // Stream an object of class AliITSgeom.
    Int_t i,j,k;

      R__b >> fNlayers;
      if(fNlad!=0) delete[] fNlad;
      if(fNdet!=0) delete[] fNdet;
      fNlad = new Int_t[fNlayers];
      fNdet = new Int_t[fNlayers];
      for(i=0;i<fNlayers;i++) R__b >> fNlad[i];
      for(i=0;i<fNlayers;i++) R__b >> fNdet[i];
      if(fg!=0){
	  for(i=0;i<fNlayers;i++) delete[] fg[i];
	  delete[] fg;
      } // end if fg!=0
      fg = new ITS_geom*[fNlayers];
      for(i=0;i<fNlayers;i++){
	  fg[i] = new ITS_geom[fNlad[i]*fNdet[i]];
	  for(j=0;j<fNlad[i]*fNdet[i];j++){
	      R__b >> fg[i][j].fShapeIndex;
	      R__b >> fg[i][j].fx0;
	      R__b >> fg[i][j].fy0;
	      R__b >> fg[i][j].fz0;
	      R__b >> fg[i][j].frx;
	      R__b >> fg[i][j].fry;
	      R__b >> fg[i][j].frz;
	      for(k=0;k<9;k++) R__b >> fg[i][j].fr[k];
	  } // end for j
      } // end for i
//      R__b >> fShape;
      return R__b;
}
