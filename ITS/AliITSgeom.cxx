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
Revision 1.12  2001/01/17 15:41:27  barbera
Some corrections suggested by Peter Hristov added

Revision 1.11  2000/10/02 16:32:35  barbera
Forward declaration added

Revision 1.4.4.15  2000/10/02 15:52:05  barbera
Forward declaration added

Revision 1.10  2000/09/05 14:25:50  nilsen
Made fixes for HP compiler. All function parameter default values placed
in .h file. Fixed the usual problem with HP comilers and the "for(Int_t i..."
business. Replaced casting (Double_t [3][3]) to (Double_t (*)[3]) for HP.
Lastly removed all "const" before function parameters which were 2 dim. arrays,
because on HP root generates some strange code (?). Thanks Peter for the
changes.

Revision 1.9  2000/08/29 20:19:03  nilsen
Removed dependancy on structure AliITSeomS and replaced it with class
AliITSgeomMatrix. Added many new functions with many new arguments. Most
in the form of in line functions for speed.

Revision 1.4.4.6  2000/06/04 16:33:32  Nilsen
A restructured AliITSgeom class. Now used AliITSgeomMatrix.

Revision 1.4.4.5  2000/03/04 23:42:39  Nilsen
Updated the comments/documentations and improved the maintainability of the
code.

Revision 1.4.4.4  2000/03/02 21:27:07  Nilsen
Added two functions, SetByAngles and SetTrans.

Revision 1.4.4.3  2000/01/23 03:09:10  Nilsen
// fixed compiler warnings for new function LtLErrorMatrix(...)

Revision 1.4.4.2  2000/01/19 23:18:20  Nilsen
Added transformations of Error matrix to AliITSgeom and fixed some typos
in AliITS.h and AliITShitIndex.h

Revision 1.4.4.1  2000/01/12 19:03:32  Nilsen
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
//     The local coordinate system by, default, is show in the following
// figures. Also shown are the ladder numbering scheme.
//Begin_Html
/*
<img src="picts/ITS/AliITSgeomMatrix_L1.gif">
</pre>
<br clear=left>
<font size=+2 color=blue>
<p>This shows the relative geometry differences between the ALICE Global
coordinate system and the local detector coordinate system.
</font>
<pre>

<pre>
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
//
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
// AliITSgeomMatrix *fGm
//     A pointer to an array of AliITSgeomMatrix classes. One element 
// per module (detector) in the ITS. AliITSgeomMatrix basicly contains
// all of the necessary information about the detector and it's coordinate
// transformations.
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
////////////////////////////////////////////////////////////////////////
#include <iostream.h>
#include <fstream.h>
#include <iomanip.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <TObject.h>
#include <TRandom.h>

#include "AliITSgeom.h"
#include "AliITSgeomMatrix.h"
#include "AliITSgeomSPD.h"
#include "AliITSgeomSDD.h"
#include "AliITSgeomSSD.h"

ClassImp(AliITSgeom)

//_____________________________________________________________________
AliITSgeom::AliITSgeom(){
////////////////////////////////////////////////////////////////////////
//     The default constructor for the AliITSgeom class. It, by default,
// sets fNlayers to zero and zeros all pointers.
////////////////////////////////////////////////////////////////////////
  // Default constructor.
  // Do not allocate anything zero everything
   fTrans   = 0; // standard GEANT global/local coordinate system.
   fNlayers = 0;
   fNlad    = 0;
   fNdet    = 0;
   fGm      = 0;
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
      for(Int_t i=0;i<fNmodules;i++) delete fGm[i];
      delete [] fGm;
   } // end if fGm!=0
   if(fNlad!=0) delete[] fNlad;
   if(fNdet!=0) delete[] fNdet;
   fNlayers = 0;
   fNlad    = 0;
   fNdet    = 0;
   fGm      = 0;
   return;
}
//______________________________________________________________________
void AliITSgeom::ReadNewFile(const char *filename){
    printf("New file format not defined yet\n");
    return;
}
//_____________________________________________________________________
AliITSgeom::AliITSgeom(const char *filename){
////////////////////////////////////////////////////////////////////////
//     The constructor for the AliITSgeom class. All of the data to fill
// this structure is read in from the file given my the input filename.
////////////////////////////////////////////////////////////////////////
   FILE     *pf;
   Int_t    i,lm=0,id[3];
   Int_t    l,a,d;
   Float_t  x,y,z,o,p,q,r,s,t;
   Double_t rot6[6],tran[3];
   char     buf[200],*buff=0; // input character buffer;

   pf = fopen(filename,"r");

   fNlayers = 6; // set default number of ladders
TryAgain:
   fNlad    = new Int_t[fNlayers];
   fNdet    = new Int_t[fNlayers];
   fNmodules = 0;
   // find the number of ladders and detectors in this geometry.
   for(i=0;i<fNlayers;i++){fNlad[i]=fNdet[i]=0;} // zero out arrays
   while(fgets(buf,200,pf)!=NULL){ // for ever loop
      for(i=0;i<200;i++)if(buf[i]!=' '){ // remove blank spaces.
           buff = &(buf[i]);
           break;
      } // end for i
      // remove blank lines and comments.
      if(buff[0]=='\n'||buff[0]=='#'||buff[0]=='!'||
         (buff[0]=='/'&&buff[1]=='/')) continue;
      if(isalpha(buff[0])) { // must be the new file formated file.
            fclose(pf);
            delete[] fNlad;delete[] fNdet;
            ReadNewFile(filename);
            return;
      } // end if isalpha(buff[0])
      sscanf(buff,"%d %d %d %f %f %f %f %f %f %f %f %f",
                  &l,&a,&d,&x,&y,&z,&o,&p,&q,&r,&s,&t);
      if(l>lm) lm = l;
      if(l<1 || l>fNlayers) {
         printf("error in file %s layer=%d min. is 1 max is %d/n",
                 filename,l,fNlayers);
         continue;
      }// end if l
      fNmodules++;
      if(l<=fNlayers&&fNlad[l-1]<a) fNlad[l-1] = a;
      if(l<=fNlayers&&fNdet[l-1]<d) fNdet[l-1] = d;
   } // end while ever loop
   if(lm>fNlayers){
	delete[] fNlad;
	delete[] fNdet;
	fNlayers = lm;
	goto TryAgain;
   } // end if lm>fNlayers
   // counted the number of ladders and detectors now allocate space.
   fGm = new AliITSgeomMatrix*[fNmodules];

   // Set up Shapes for a default configuration of 6 layers.
   fTrans   = 0; // standard GEANT global/local coordinate system.
   fShape = new TObjArray(3);
   AddShape((TObject *) new AliITSgeomSPD());  // shape 0
   AddShape((TObject *) new AliITSgeomSDD());  // shape 1
   AddShape((TObject *) new AliITSgeomSPD());  // shape 2

   // prepare to read in transforms
   lm = 0; // reuse lm as counter of modules.
   rewind(pf); // start over reading file
   while(fgets(buf,200,pf)!=NULL){ // for ever loop
      for(i=0;i<200;i++)if(buf[i]!=' '){ // remove blank spaces.
           buff = &(buf[i]);
           break;
      } // end for i
      // remove blank lines and comments.
      if(buff[0]=='\n'||buff[0]=='#'||buff[0]=='!'||
        (buff[0]=='/'&&buff[1]=='/')) continue;
      x = y = z = o = p = q = r = s = t = 0.0;
      sscanf(buff,"%d %d %d %f %f %f %f %f %f %f %f %f",
                  &l,&a,&d,&x,&y,&z,&o,&p,&q,&r,&s,&t);
      if(l<1 || l>fNlayers) {
         printf("error in file %s layer=%d min. is 1 max is %d/n",
                 filename,l,fNlayers);
         continue;
      }// end if l
      id[0] = l;id[1] = a;id[2] = d;
      tran[0] = tran[1] = tran[2]  = 0.0;
      tran[0] = (Double_t)x;tran[1] = (Double_t)y;tran[2] = (Double_t)z;
      rot6[0] = rot6[1] = rot6[2] = rot6[3] = rot6[4] = rot6[5] =0.0;
      rot6[0] = (Double_t)o;rot6[1] = (Double_t)p;rot6[2] = (Double_t)q;
      rot6[3] = (Double_t)r;rot6[4] = (Double_t)s;rot6[5] = (Double_t)t;
      switch (l){
      case 1: case 2: // layer 1 or2 SPD
          fGm[lm++] = new AliITSgeomMatrix(rot6,0,id,tran);
	  break;
      case 3: case 4: // layer 3 or 4 SDD
          fGm[lm++] = new AliITSgeomMatrix(rot6,1,id,tran);
	  break;
      case 5: case 6: // layer 5 or 6 SSD
          fGm[lm++] = new AliITSgeomMatrix(rot6,2,id,tran);
	  break;
      } // end switch
   } // end while ever loop
   fclose(pf);
}

//________________________________________________________________________
AliITSgeom::AliITSgeom(AliITSgeom &source){
////////////////////////////////////////////////////////////////////////
//     The copy constructor for the AliITSgeom class. It calls the
// = operator function. See the = operator function for more details.
////////////////////////////////////////////////////////////////////////

    *this = source;  // Just use the = operator for now.

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
   Int_t i;

   if(this == &source) return; // don't assign to ones self.

   // if there is an old structure allocated delete it first.
   if(this->fGm != 0){
      for(i=0;i<this->fNmodules;i++) delete this->fGm[i];
      delete this->fGm;
   } // end if fGm != 0 
   if(fNlad != 0) delete[] fNlad;
   if(fNdet != 0) delete[] fNdet;

   this->fTrans    = source.fTrans;
   this->fNmodules = source.fNmodules;
   this->fNlayers = source.fNlayers;
   this->fNlad = new Int_t[fNlayers];
   for(i=0;i<this->fNlayers;i++) this->fNlad[i] = source.fNlad[i];
   this->fNdet = new Int_t[fNlayers];
   for(i=0;i<this->fNlayers;i++) this->fNdet[i] = source.fNdet[i];
   this->fShape = new TObjArray(*(source.fShape));//This does not make a proper copy.
   this->fGm = new AliITSgeomMatrix*[this->fNmodules];
   for(i=0;i<this->fNmodules;i++){
       this->fGm[i] = new AliITSgeomMatrix(*(source.fGm[i]));
   } // end for i
   return;
}//_____________________________________________________________________
Int_t AliITSgeom::GetModuleIndex(const Int_t lay,const Int_t lad,
				 const Int_t det){
////////////////////////////////////////////////////////////////////////
//      This routine computes the module index number from the layer,
// ladder, and detector numbers. The number of ladders and detectors
// per layer is determined when this geometry package is constructed,
// see AliITSgeom(const char *filename) for specifics.
////////////////////////////////////////////////////////////////////////
    Int_t i,j,k,id[3];

    i = fNdet[lay-1] * (lad-1) + det - 1;
    j = 0;
    for(k=0;k<lay-1;k++) j += fNdet[k]*fNlad[k];
    i = i+j;
    fGm[i]->GetIndex(id);
    if(id[0]==lay&&id[1]==lad&&id[2]==det) return i;
    // Array of modules fGm is not in expected order. Search for this index
    for(i=0;i<fNmodules;i++){
	fGm[i]->GetIndex(id);
	if(id[0]==lay&&id[1]==lad&&id[2]==det) return i;
    } // end for i
    // This layer ladder and detector combination does not exist return -1.
    return -1;
}
//______________________________________________________________________
void AliITSgeom::GetModuleId(const Int_t index,
			     Int_t &lay,Int_t &lad,Int_t &det){
////////////////////////////////////////////////////////////////////////
//      This routine computes the layer, ladder and detector number 
// given the module index number. The number of ladders and detectors
// per layer is determined when this geometry package is constructed,
// see AliITSgeom(const char *filename) for specifics.
////////////////////////////////////////////////////////////////////////
    Int_t id[3];

    fGm[index]->GetIndex(id);
    lay = id[0]; lad = id[1]; det = id[2];
    return;

    // The old way kept for posterity.
/*
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
*/
}
//___________________________________________________________________________
Int_t AliITSgeom::GetStartDet(const Int_t dtype){
  /////////////////////////////////////////////////////////////////////////
  // returns the starting module index value for a give type of detector id
  /////////////////////////////////////////////////////////////////////////

  switch(dtype){
  case 0:
     return GetModuleIndex(1,1,1);
     break;
  case 1:
     return GetModuleIndex(3,1,1);
     break;
  case 2:
     return GetModuleIndex(5,1,1);
     break;
  default:
     printf("<AliITSgeom::GetFirstDet> undefined detector type\n");
     return 0;
  } // end switch

  printf("<AliITSgeom::GetFirstDet> undefined detector type\n");
  return 0;
}

//___________________________________________________________________________
Int_t AliITSgeom::GetLastDet(const Int_t dtype){
  /////////////////////////////////////////////////////////////////////////
  // returns the last module index value for a give type of detector id
  /////////////////////////////////////////////////////////////////////////

  switch(dtype){
  case 0:
     return GetLastSPD();
     break;
   case 1:
     return GetLastSDD();
     break;
   case 2:
     return GetLastSSD();
     break;
   default:
     printf("<AliITSgeom::GetLastDet> undefined detector type\n");
     return 0;
  } // end switch

  printf("<AliITSgeom::GetLastDet> undefined detector type\n");
  return 0;
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
   Int_t    i,j,idt[3],ido[3];
   Double_t tt[3],to[3];  // translation
   Double_t rt[3],ro[3];  // phi in radians
   Double_t mt[3][3],mo[3][3]; // matrixes
   AliITSgeomMatrix *gt,*go;
   Bool_t   t;

   for(i=0;i<this->fNmodules;i++){
         gt  =  this->GetGeomMatrix(i);
	 go  = other->GetGeomMatrix(i);
         gt->GetIndex(idt);
         go->GetIndex(ido);
         t = kFALSE;
         for(i=0;i<3;i++) t = t&&idt[i]!=ido[i];
         if(t) fprintf(fp,"%4.4d %1.1d %2.2d %2.2d %1.1d %2.2d %2.2d\n",i,
                       idt[0],idt[1],idt[2],ido[0],ido[1],ido[2]);
         gt->GetTranslation(tt);
         go->GetTranslation(to);
         gt->GetAngles(rt);
         go->GetAngles(ro);
         t = kFALSE;
         for(i=0;i<3;i++) t = t&&tt[i]!=to[i];
	 if(t) fprintf(fp,"%1.1d %2.2d %2.2d dTrans=%f %f %f drot=%f %f %f\n",
		       idt[0],idt[1],idt[2],
                       tt[0]-to[0],tt[1]-to[1],tt[2]-to[2],
                       rt[0]-ro[0],rt[1]-ro[1],rt[2]-ro[2]);
	 t = kFALSE;
         gt->GetMatrix(mt);
         go->GetMatrix(mo);
	 for(i=0;i<3;i++)for(j=0;j<3;j++)  t = mt[i][j] != mo[i][j];
	 if(t){
	     fprintf(fp,"%1.1d %2.2d %2.2d dfr= %e %e %e\n",
                     idt[0],idt[1],idt[2],
                 mt[0][0]-mo[0][0],mt[0][1]-mo[0][1],mt[0][2]-mo[0][2]);
	     fprintf(fp,"        dfr= %e %e %e\n",
                 mt[1][0]-mo[1][0],mt[1][1]-mo[1][1],mt[1][2]-mo[1][2]);
	     fprintf(fp,"        dfr= %e %e %e\n",
                 mt[2][0]-mo[2][0],mt[2][1]-mo[2][1],mt[2][2]-mo[2][2]);
	 } // end if t
   } // end for i
   return;
}

//___________________________________________________________________________
void AliITSgeom::PrintData(FILE *fp,
			   const Int_t lay,const Int_t lad,const Int_t det){
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
   AliITSgeomMatrix *gt;
   Double_t t[3],r[3],m[3][3];

   gt = this->GetGeomMatrix(GetModuleIndex(lay,lad,det));
   gt->GetTranslation(t);
   gt->GetAngles(r);
   fprintf(fp,"%1.1d %2.2d %2.2d Trans=%f %f %f rot=%f %f %f Shape=%d\n",
	   lay,lad,det,t[0],t[1],t[2],r[0],r[1],r[2],
           gt->GetDetectorIndex());
   gt->GetMatrix(m);
   fprintf(fp,"        dfr= %e %e %e\n",m[0][0],m[0][1],m[0][2]);
   fprintf(fp,"        dfr= %e %e %e\n",m[1][0],m[1][1],m[1][2]);
   fprintf(fp,"        dfr= %e %e %e\n",m[2][0],m[2][1],m[2][2]);
   return;
}
//___________________________________________________________________________
ofstream & AliITSgeom::PrintGeom(ofstream &R__b){
////////////////////////////////////////////////////////////////////////
//     The default Streamer function "written by ROOT" doesn't write out
// the arrays referenced by pointers. Therefore, a specific Streamer function
// has to be written. This function should not be modified but instead added
// on to so that older versions can still be read. The proper handling of
// the version dependent streamer function hasn't been written do to the lack
// of finding an example at the time of writing.
////////////////////////////////////////////////////////////////////////
   // Stream an object of class AliITSgeom.
    Int_t i;

    R__b.setf(ios::scientific);
    R__b << fTrans << " ";
    R__b << fNmodules << " ";
    R__b << fNlayers << " ";
    for(i=0;i<fNlayers;i++) R__b << fNlad[i] << " ";
    for(i=0;i<fNlayers;i++) R__b << fNdet[i] << "\n";
    for(i=0;i<fNmodules;i++) {
	R__b <<setprecision(16) << *(fGm[i]) << "\n";
    } // end for i
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
// of finding an example at the time of writing.
////////////////////////////////////////////////////////////////////////
   // Stream an object of class AliITSgeom.
      Int_t i;

      fNlad = new Int_t[fNlayers];
      fNdet = new Int_t[fNlayers];
      if(fGm!=0){
	  for(i=0;i<fNmodules;i++) delete fGm[i];
	  delete fGm;
      } // end if fGm!=0

      R__b >> fTrans >> fNmodules >> fNlayers;
      fNlad = new Int_t[fNlayers];
      fNdet = new Int_t[fNlayers];
      for(i=0;i<fNlayers;i++) R__b >> fNlad[i];
      for(i=0;i<fNlayers;i++) R__b >> fNdet[i];
      fGm = new AliITSgeomMatrix*[fNmodules];
      for(i=0;i<fNmodules;i++){
	  fGm[i] = new AliITSgeomMatrix;
	  R__b >> *(fGm[i]);
      } // end for i
      return R__b;
}

//______________________________________________________________________
//     The following routines modify the transformation of "this"
// geometry transformations in a number of different ways.
//______________________________________________________________________
void AliITSgeom::GlobalChange(const Float_t *tran,const Float_t *rot){
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
   Int_t    i,j;
   Double_t t[3],r[3];
   AliITSgeomMatrix *g;

   fTrans = (fTrans && 0xfffd) + 2;  // set bit 1 true.
   for(i=0;i<fNmodules;i++){
         g = this->GetGeomMatrix(i);
         g->GetTranslation(t);
         g->GetAngles(r);
         for(j=0;j<3;j++){
              t[j] += tran[j];
              r[j] += rot[j];
         } // end for j
         g->SetTranslation(t);
         g->SetAngles(r);
   } // end for i
   return;
}
//___________________________________________________________________________
void AliITSgeom::GlobalCylindericalChange(const Float_t *tran,const Float_t *rot){
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
   Int_t    i,j;
   Double_t t[3],ro[3],r,r0,phi,rphi;
   AliITSgeomMatrix *g;

   fTrans = (fTrans && 0xfffd) + 2;  // set bit 1 true.
   for(i=0;i<fNmodules;i++){
         g = this->GetGeomMatrix(i);
         g->GetTranslation(t);
         g->GetAngles(ro);
	 r = r0= TMath::Hypot(t[1],t[0]);
	 phi   = TMath::ATan2(t[1],t[0]);
	 rphi  = r0*phi;
	 r    += tran[0];
	 rphi += tran[1];
	 phi   = rphi/r0;
         t[0]  = r*TMath::Cos(phi);
         t[1]  = r*TMath::Sin(phi);
         t[2] += tran[2];
         for(j=0;j<3;j++){
              ro[j] += rot[j];
         } // end for j
         g->SetTranslation(t);
         g->SetAngles(ro);
   } // end for i
   return;
}
//___________________________________________________________________________
void AliITSgeom::RandomChange(const Float_t *stran,const Float_t *srot){
////////////////////////////////////////////////////////////////////////
//     This function performs a Gaussian random displacement and/or
// rotation about the present global position of each active
// volume/detector of the ITS. The sigma of the random displacement
// is determined by the three element array stran, for the
// x y and z translations, and the three element array srot,
// for the three rotation about the axis x y and z.
////////////////////////////////////////////////////////////////////////
   Int_t    i,j;
   Double_t t[3],r[3];
   TRandom ran;
   AliITSgeomMatrix *g;

   fTrans = (fTrans && 0xfffd) + 2;  // set bit 1 true.
   for(i=0;i<fNmodules;i++){
         g = this->GetGeomMatrix(i);
         g->GetTranslation(t);
         g->GetAngles(r);
         for(j=0;j<3;j++){
              t[j] += ran.Gaus(0.0,stran[j]);
              r[j] += ran.Gaus(0.0, srot[j]);
         } // end for j
         g->SetTranslation(t);
         g->SetAngles(r);
   } // end for i
   return;
}
//___________________________________________________________________________
void AliITSgeom::RandomCylindericalChange(const Float_t *stran,
					  const Float_t *srot){
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
   Int_t    i,j;
   Double_t t[3],ro[3],r,r0,phi,rphi;
   TRandom ran;
   AliITSgeomMatrix *g;

   fTrans = (fTrans && 0xfffd) + 2;  // set bit 1 true.
   for(i=0;i<fNmodules;i++){
         g = this->GetGeomMatrix(i);
         g->GetTranslation(t);
         g->GetAngles(ro);
	 r = r0= TMath::Hypot(t[1],t[0]);
	 phi   = TMath::ATan2(t[1],t[0]);
	 rphi  = r0*phi;
	 r    += ran.Gaus(0.0,stran[0]);
	 rphi += ran.Gaus(0.0,stran[1]);
	 phi   = rphi/r0;
         t[0]  = r*TMath::Cos(phi);
         t[1]  = r*TMath::Sin(phi);
         t[2] += ran.Gaus(0.0,stran[2]);
         for(j=0;j<3;j++){
              ro[j] += ran.Gaus(0.0, srot[j]);
         } // end for j
         g->SetTranslation(t);
         g->SetAngles(ro);
   } // end for i
   return;
}
//______________________________________________________________________
void AliITSgeom::GeantToTracking(AliITSgeom &source){
/////////////////////////////////////////////////////////////////////////
//     Copy the geometry data but change it to go between the ALICE
// Global coordinate system to that used by the ITS tracking. A slightly
// different coordinate system is used when tracking. This coordinate 
// system is only relevant when the geometry represents the cylindrical
// ALICE ITS geometry. For tracking the Z axis is left alone but X-> -Y
// and Y-> X such that X always points out of the ITS cylinder for every
// layer including layer 1 (where the detectors are mounted upside down).
//Begin_Html
/*
<img src="picts/ITS/AliITSgeomMatrix_T1.gif">
*/
//End_Html
////////////////////////////////////////////////////////////////////////
   Int_t    i,j,k,l,id[3];
   Double_t R0[3][3],R1[3][3];
   Double_t A0[3][3] = {{0.,+1.,0.},{-1.,0.,0.},{0.,0.,+1.}};
   Double_t A1[3][3] = {{0.,-1.,0.},{+1.,0.,0.},{0.,0.,+1.}};

   *this = source;  // copy everything
   for(i=0;i<GetIndexMax();i++){
       fGm[i]->GetIndex(id);
       fGm[i]->GetMatrix(R0);
       if(id[0]==1){ // Layer 1 is treated different from the others.
	   for(j=0;j<3;j++) for(k=0;k<3;k++){
	       R1[j][k] = 0.;
	       for(l=0;l<3;l++) R1[j][k] += A0[j][l]*R0[l][k];
	   } // end for j,k
       }else{
	   for(j=0;j<3;j++) for(k=0;k<3;k++){
	       R1[j][k] = 0.;
	       for(l=0;l<3;l++) R1[j][k] += A1[j][l]*R0[l][k];
	   } // end for j,k
       } // end if
       fGm[i]->SetMatrix(R1);
   } // end for i
   this->fTrans = (this->fTrans && 0xfffe) + 1;  // set bit 0 true.
   return;
}
//______________________________________________________________________
Int_t AliITSgeom::GetNearest(const Double_t g[3],const Int_t lay){
////////////////////////////////////////////////////////////////////////
//      Finds the Detector (Module) that is nearest the point g [cm] in
// ALICE Global coordinates. If layer !=0 then the search is restricted
// to Detectors (Modules) in that particular layer.
////////////////////////////////////////////////////////////////////////
     Int_t    i,l,a,e,in=0;
     Double_t d,dn=1.0e10;
     Bool_t   t=lay!=0; // skip if lay = 0 default value check all layers.

     for(i=0;i<fNmodules;i++){
          if(t){GetModuleId(i,l,a,e);if(l!=lay) continue;}
          if((d=fGm[i]->Distance2(g))<dn){
               dn = d;
               in = i;
          } // end if
     } // end for i
     return in;
}
//______________________________________________________________________
void AliITSgeom::GetNearest27(const Double_t g[3],Int_t n[27],const Int_t lay){
////////////////////////////////////////////////////////////////////////
//      Finds 27 Detectors (Modules) that are nearest the point g [cm] in
// ALICE Global coordinates. If layer !=0 then the search is restricted
// to Detectors (Modules) in that particular layer. The number 27 comes 
// from including the nearest detector and all those around it (up, down,
// left, right, forwards, backwards, and the corners).
////////////////////////////////////////////////////////////////////////
     Int_t    i,l,a,e,in[27]={0,0,0,0,0,0,0,0,0,
                              0,0,0,0,0,0,0,0,0,
                              0,0,0,0,0,0,0,0,0,};
     Double_t d,dn[27]={1.0e10,1.0e10,1.0e10,1.0e10,1.0e10,1.0e10,
                        1.0e10,1.0e10,1.0e10,1.0e10,1.0e10,1.0e10,
                        1.0e10,1.0e10,1.0e10,1.0e10,1.0e10,1.0e10,
                        1.0e10,1.0e10,1.0e10,1.0e10,1.0e10,1.0e10,
                        1.0e10,1.0e10,1.0e10};
     Bool_t   t=(lay!=0); // skip if lay = 0 default value check all layers.

     for(i=0;i<fNmodules;i++){
          if(t){GetModuleId(i,l,a,e);if(l!=lay) continue;}
          for(a=0;a<27;a++){
               d = fGm[i]->Distance2(g);
               if(d<dn[a]){
                   for(e=26;e>a;e--){dn[e] = dn[e-1];in[e] = in[e-1];}
                   dn[a] = d; in[a] = i;
               } // end if d<dn[i]
          } // end for a
     } // end for i
     for(i=0;i<27;i++) n[i] = in[i];
}
//----------------------------------------------------------------------
