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

/* $Id$ */

///////////////////////////////////////////////////////////////////////
// ITS geometry manipulation routines.                               //
// Created April 15 1999.                                            //
// version: 0.0.0                                                    //
// By: Bjorn S. Nilsen                                               //
// version: 0.0.1                                                    //
// Updated May 27 1999.                                              //
// Added Cylindrical random and global based changes.                //
//                                                                   //
// Modified and added functions Feb. 7 2006                          //
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
// TString    fVersion 
//     Transformation version.
// Int_t      fTrans
//     Flag to keep track of which transformation 
// Int_t      fNmodules
//      The total number of modules
// Int_t fNlayers
//     The number of ITS layers for this geometry. By default this
//  is 6, but can be modified by the creator function if there are
// more layers defined.
//
// TArrayI fNlad
//     A pointer to an array fNlayers long containing the number of 
// ladders for each layer. This array is typically created and filled 
// by the AliITSgeom creator function.
//
// TArrayI fNdet
//     A pointer to an array fNlayers long containing the number of
// active detector volumes for each ladder. This array is typically
// created and filled by the AliITSgeom creator function.
//
// TObjArray fGm containing objects of type AliITSgeomMatrix
//     A pointer to an array of AliITSgeomMatrix classes. One element 
// per module (detector) in the ITS. AliITSgeomMatrix basicly contains
// all of the necessary information about the detector and it's coordinate
// transformations.
//
////////////////////////////////////////////////////////////////////////
#include <Riostream.h>
#include <ctype.h>

#include <TRandom.h>
#include <TSystem.h>
#include <TArrayI.h>

#include "AliITSgeom.h"
#include "AliLog.h"

ClassImp(AliITSgeom)

//______________________________________________________________________
AliITSgeom::AliITSgeom():
TObject(),
fVersion("GEANT"),// Transformation version.
fTrans(0),       // Flag to keep track of which transformation 
fNmodules(0),    // The total number of modules
fNlayers(0),     // The number of layers.
fNlad(),         //[] Array of the number of ladders/layer(layer)
fNdet(),         //[] Array of the number of detector/ladder(layer)
fGm(0,0)        // Structure of translation. and rotation.
{
    //     The default constructor for the AliITSgeom class. It, by default,
    // sets fNlayers to zero and zeros all pointers.
    // Do not allocate anything zero everything.
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    a zeroed AliITSgeom object.

    fGm.SetOwner(kTRUE);
    return;
}

//______________________________________________________________________
AliITSgeom::AliITSgeom(Int_t itype,Int_t nlayers,const Int_t *nlads,
                       const Int_t *ndets,Int_t mods):
TObject(),
fVersion("GEANT"),    // Transformation version.
fTrans(itype),       // Flag to keep track of which transformation 
fNmodules(mods),     // The total number of modules
fNlayers(nlayers),   // The number of layers.
fNlad(nlayers,nlads),//[] Array of the number of ladders/layer(layer)
fNdet(nlayers,ndets),//[] Array of the number of detector/ladder(layer)
fGm(mods,0)         // Structure of translation. and rotation.
{
    //     A simple constructor to set basic geometry class variables
    // Inputs:
    //      Int_t itype   the type of transformation kept.
    //                    bit 0 => Standard GEANT
    //                    bit 1 => ITS tracking
    //                    bit 2 => A change in the coordinate system 
    //                    has been made. others are still to be defined 
    //                    as needed.
    //      Int_t nlayers The number of ITS layers also set the size of 
    //                    the arrays
    //      Int_t *nlads  an array of the number of ladders for each 
    //                    layer. This array must be nlayers long.
    //      Int_t *ndets  an array of the number of detectors per ladder
    //                    for each layer. This array must be nlayers long.
    //      Int_t mods    The number of modules. Typically the sum of all the 
    //                    detectors on every layer and ladder.
    // Outputs:
    //     none
    // Return:
    //     A properly inilized AliITSgeom object.

    fGm.SetOwner(kTRUE);
    return;
}
//______________________________________________________________________
void AliITSgeom::Init(Int_t itype,Int_t nlayers,const Int_t *nlads,
                      const Int_t *ndets,Int_t mods){
    //     A simple Inilizer to set basic geometry class variables
    // Inputs:
    //      Int_t itype   the type of transformation kept.
    //                    bit 0 => Standard GEANT
    //                    bit 1 => ITS tracking
    //                    bit 2 => A change in the coordinate system 
    //                    has been made. others are still to be defined 
    //                    as needed.
    //      Int_t nlayers The number of ITS layers also set the size of 
    //                    the arrays
    //      Int_t *nlads  an array of the number of ladders for each 
    //                    layer. This array must be nlayers long.
    //      Int_t *ndets  an array of the number of detectors per ladder 
    //                    for each layer. This array must be nlayers long.
    //      Int_t mods    The number of modules. Typically the sum of all the 
    //                    detectors on every layer and ladder.
    // Outputs:
    //     none
    // Return:
    //     A properly inilized AliITSgeom object.

    fVersion  = "GEANT";     // Transformation version.
    fTrans    = itype;       // Flag to keep track of which transformation 
    fNmodules = mods;        // The total number of modules
    fNlayers  = nlayers;     // The number of layers.
    fNlad.Set(nlayers,nlads);//[] Array of the number of ladders/layer(layer)
    fNdet.Set(nlayers,ndets);//[] Array of the number of detector/ladder(layer)
    fGm.Clear();
    fGm.Expand(mods);        // Structure of translation. and rotation.
    fGm.SetOwner(kTRUE);
    return;
}
//______________________________________________________________________
void AliITSgeom::CreateMatrix(Int_t mod,Int_t lay,Int_t lad,Int_t det,
                              AliITSDetector idet,const Double_t tran[3],
                              const Double_t rot[10]){
    // Given the translation vector tran[3] and the rotation matrix rot[1],
    // this function creates and adds to the TObject Array fGm the
    // AliITSgeomMatrix object.
    // The rot[10] matrix is set up like:
    /*   / rot[0]  rot[1]  rot[2] \
    //  |  rot[3]  rot[4]  rot[5]  |
    //   \ rot[6]  rot[7]  rot[8] /  if(rot[9]!=0) then the Identity matrix
    // is used regardless of the values in rot[0]-rot[8].
    */
    // Inputs:
    //    Int_t           mod     The module number. The location in TObjArray
    //    Int_t           lay     The layer where this module is
    //    Int_t           lad     On which ladder this module is
    //    Int_t           det     Which detector on this ladder this module is
    //    AliITSDetector idet     The type of detector see AliITSgeom.h
    //    Double_t       tran[3]  The translation vector
    //    Double_t       rot[10]  The rotation matrix.
    // Outputs:
    //    none
    // Return:
    //    none.
    Int_t id[3];
    Double_t r[3][3] = {{1.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};

    if(mod<0||mod>=fGm.GetSize()){ 
	Error("CreateMatrix","mod=%d is out of bounds max value=%d",mod,
	      fGm.GetSize());
	return;
    } // end if
    delete fGm.At(mod);
    id[0] = lay; id[1] = lad; id[2] = det;
    if(rot[9]!=0.0) { // null rotation
        r[0][0] = rot[0]; r[0][1] = rot[1]; r[0][2] = rot[2];
        r[1][0] = rot[3]; r[1][1] = rot[4]; r[1][2] = rot[5];
        r[2][0] = rot[6]; r[2][1] = rot[7]; r[2][2] = rot[8];
    } // end if
    fGm.AddAt(new AliITSgeomMatrix(idet,id,r,tran),mod);
}
//______________________________________________________________________
AliITSgeom::~AliITSgeom(){
    //     The destructor for the AliITSgeom class. If the arrays fNlad,
    // fNdet, or fGm have had memory allocated to them, there pointer values
    // are non zero, then this memory space is freed and they are set
    // to zero. In addition, fNlayers is set to zero. The destruction of
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    none.

    return;
}
//______________________________________________________________________
AliITSgeom::AliITSgeom(const AliITSgeom &source) : 
TObject(source),
fVersion(source.fVersion), // Transformation version.
fTrans(source.fTrans),   // Flag to keep track of which transformation
fNmodules(source.fNmodules),// The total number of modules
fNlayers(source.fNlayers), // The number of layers.
fNlad(source.fNlad),    // Array of the number of ladders/layer(layer)
fNdet(source.fNdet),    // Array of the number of detector/ladder(layer)
fGm(source.fGm.GetSize(),source.fGm.LowerBound())// Structure of 
                                                  // translation and rotation.
{
    //     The copy constructor for the AliITSgeom class. It calls the
    // = operator function. See the = operator function for more details.
    // Inputs:
    //     AliITSgeom &source  The AliITSgeom class with which to make this
    //                         a copy of.
    // Outputs:
    //     none.
    // Return:
    //     none.
    Int_t i,n;

    n = source.fGm.GetLast()+1;
    for(i=source.fGm.LowerBound();i<n;i++){
        fGm.AddAt(new AliITSgeomMatrix(*((AliITSgeomMatrix*)(
                                             source.fGm.At(i)))),i);
    } // end for i
    fGm.SetOwner(kTRUE);
    return;
}
//______________________________________________________________________
AliITSgeom& AliITSgeom::operator=(const AliITSgeom &source){
    //     The = operator function for the AliITSgeom class. It makes an
    // independent copy of the class in such a way that any changes made
    // to the copied class will not affect the source class in any way.
    // This is required for many ITS alignment studies where the copied
    // class is then modified by introducing some misalignment.
    // Inputs:
    //     AliITSgeom &source  The AliITSgeom class with which to make this
    //                         a copy of.
    // Outputs:
    //     none.
    // Return:
    //     *this The a new copy of source.
    Int_t i;

    if(this == &source) return *this; // don't assign to ones self.

    // if there is an old structure allocated delete it first.
    this->fGm.Clear();

    this->fVersion  = source.fVersion;
    this->fTrans    = source.fTrans;
    this->fNmodules = source.fNmodules;
    this->fNlayers  = source.fNlayers;
    this->fNlad     = source.fNlad;
    this->fNdet     = source.fNdet;
    this->fGm.Expand(this->fNmodules);
    for(i=source.fGm.LowerBound();i<source.fGm.GetLast();i++){
        fGm.AddAt(new AliITSgeomMatrix(*((AliITSgeomMatrix*)(
                                             source.fGm.At(i)))),i);
    } // end for i
    fGm.SetOwner(kTRUE);
    return *this;
}
//______________________________________________________________________
Int_t AliITSgeom::GetModuleIndex(Int_t lay,Int_t lad,Int_t det)const{
    //      This routine computes the module index number from the layer,
    // ladder, and detector numbers. The number of ladders and detectors
    // per layer is determined when this geometry package is constructed,
    // see AliITSgeom(const char *filename) for specifics.
    // Inputs:
    //    Int_t lay  The layer number. Starting from 1.
    //    Int_t lad  The ladder number. Starting from 1.
    //    Int_t det  The detector number. Starting from 1.
    // Outputs:
    //    none.
    // Return:
    //    the module index number, starting from zero.
    Int_t i,j,k,id[3];

    i = fNdet[lay-1] * (lad-1) + det - 1;
    j = 0;
    for(k=0;k<lay-1;k++) j += fNdet[k]*fNlad[k];
    i = i+j;
    if(i>=fNmodules) return -1;
    GetGeomMatrix(i)->GetIndex(id);
    if(id[0]==lay&&id[1]==lad&&id[2]==det) return i;
    // Array of modules fGm is not in expected order. Search for this index
    for(i=0;i<fNmodules;i++){
        GetGeomMatrix(i)->GetIndex(id);
        if(id[0]==lay&&id[1]==lad&&id[2]==det) return i;
    } // end for i
    // This layer ladder and detector combination does not exist return -1.
    return -1;
}
//______________________________________________________________________
void AliITSgeom::GetModuleId(Int_t index,Int_t &lay,Int_t &lad,Int_t &det)
const{
    //      This routine computes the layer, ladder and detector number 
    // given the module index number. The number of ladders and detectors
    // per layer is determined when this geometry package is constructed,
    // see AliITSgeom(const char *filename) for specifics.
    // Inputs:
    //     Int_t index  The module index number, starting from zero.
    // Outputs:
    //     Int_t lay    The layer number. Starting from 1.
    //     Int_t lad    The ladder number. Starting from 1.
    //     Int_t det    The detector number. Starting from 1.
    // Return:
    //     none.
    Int_t id[3];
    AliITSgeomMatrix *g = GetGeomMatrix(index);

    if (g == 0x0){
        Error("GetModuleId","Can not get GeoMatrix for index = %d",index);
        lay = -1; lad = -1; det = -1;
    }else{
        g->GetIndex(id);
        lay = id[0]; lad = id[1]; det = id[2];
    }// End if
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
//______________________________________________________________________
Int_t AliITSgeom::GetNDetTypes(Int_t &max)const{
    // Finds and returns the number of detector types used and the
    // maximum detector type value. Only counts id >=0 (no undefined
    // values. See AliITSgeom.h for list of AliITSDetecor enumerated types.
    // Inputs:
    //    none.
    // Outputs:
    //    The maximum detector type used
    // Return:
    //    The number of detector types used
    Int_t i,*n,id;

    max = -1;
    for(i=0;i<GetIndexMax();i++){
        id = GetModuleType(i);
        if(id>max) max=id;
    } // end for i
    n = new Int_t[max+1];
    for(i=0;i<max;i++) n[i] = 0;
    for(i=0;i<GetIndexMax();i++){
        id = GetModuleType(i);
        if(id>-1)n[id]++; // note id=-1 => undefined.
    } // end for i
    id = 0;
    for(i=0;i<max;i++) if(n[i]!=0) id++;
    delete[] n;
    return id+1;
}
//______________________________________________________________________
Int_t AliITSgeom::GetNDetTypes(TArrayI &maxs,AliITSDetector *types)const{
    // Finds and returns the number of detector types used and the
    // number of each detector type. Only counts id >=0 (no undefined
    // values. See AliITSgeom.h for list of AliITSDetecor enumerated types.
    // Inputs:
    //    none.
    // Outputs:
    //    The maximum detector type used
    // Return:
    //    The number of detector types used
    Int_t i,j,*n,id,max;

    max = -1;
    for(i=0;i<GetIndexMax();i++){
        id = GetModuleType(i);
        if(id>max) max=id;
    } // end for i
    n = new Int_t[max+1];
    for(i=0;i<max;i++) n[i] = 0;
    for(i=0;i<GetIndexMax();i++){
        id = GetModuleType(i);
        if(id>-1)n[id]++; // note id=-1 => undefined.
    } // end for i
    id = 0;
    for(i=0;i<=max;i++) if(n[i]!=0) id++;
    maxs.Set(id);
    j = 0;
    for(i=0;i<=max;i++) if(n[i]!=0){
        maxs[j] = n[i];
        types[j++] = (AliITSDetector) i;
    } // end for i/end if
    delete[] n;
    return id;
}
//______________________________________________________________________
Int_t AliITSgeom::GetStartDet(Int_t dtype)const{
    // returns the starting module index value for a give type of detector id.
    // This assumes that the detector types are different on different layers
    // and that they are not mixed up.
    // Inputs:
    //    Int_t dtype A detector type number. 0 for SPD, 1 for SDD, 
    //                and 2 for SSD.
    // Outputs:
    //    none.
    // Return:
    //    the module index for the first occurrence of that detector type.

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
        Warning("GetStartDet","undefined detector type %d",dtype);
        return 0;
    } // end switch

    Warning("GetStartDet","undefined detector type %d",dtype);
    return 0;
}
//______________________________________________________________________
Int_t AliITSgeom::GetLastDet(Int_t dtype)const{
    // returns the last module index value for a give type of detector id.
    // This assumes that the detector types are different on different layers
    // and that they are not mixed up.
    // Inputs:
    //     Int_t dtype A detector type number. 0 for SPD, 1 for SDD, 
    //                 and 2 for SSD.
    // Outputs:
    // Return:
    //     the module index for the last occurrence of that detector type.

    switch((AliITSDetector)dtype){
    case kSPD:
        return GetModuleIndex(3,1,1)-1;
        break;
    case kSDD:
        return GetModuleIndex(5,1,1)-1;
        break;
    case kSSD:
        return GetIndexMax()-1;
        break;
    case kSSDp: case kSDDp: case kND:
    default:
        Warning("GetLastDet","undefined detector type %d",dtype);
        return 0;
    } // end switch

    Warning("GetLastDet","undefined detector type %d",dtype);
    return 0;
}

//______________________________________________________________________
void AliITSgeom::PrintData(FILE *fp,Int_t lay,Int_t lad,Int_t det)const{
    //     This function prints out the coordinate transformations for
    // the particular detector defined by layer, ladder, and detector
    // to the file pointed to by the File pointer fp. fprintf statements
    // are used to print out the numbers. The format is
    // layer ladder detector Trans= fx0 fy0 fz0 rot= frx fry frz 
    // Shape=fShapeIndex
    //                         dfr= fr[0] fr[1] fr[2]
    //                         dfr= fr[3] fr[4] fr[5]
    //                         dfr= fr[6] fr[7] fr[8]
    // By indicating which detector, some control over the information 
    // is given to the user. The output it written to the file pointed
    // to by the file pointer fp. This can be set to stdout if you want.
    // Inputs:
    //     FILE *fp           A file pointer to an opened file for 
    //                        writing in which the results of the 
    //                        comparison will be written.
    //     Int_t lay          The layer number. Starting from 1.
    //     Int_t lad          The ladder number. Starting from 1.
    //     Int_t det          The detector number. Starting from 1.
    // Outputs:
    //     none
    // Return:
    //     none.
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

//______________________________________________________________________
Int_t AliITSgeom::GetNearest(const Double_t g[3],Int_t lay)const{
    //      Finds the Detector (Module) that is nearest the point g [cm] in
    // ALICE Global coordinates. If layer !=0 then the search is restricted
    // to Detectors (Modules) in that particular layer.
    // Inputs:
    //     Double_t g[3]  The ALICE Cartesian global coordinate from which the
    //                    distance is to be calculated with.
    //     Int_t lay      The layer to restrict the search to. If layer=0 then
    //                    all layers are searched. Default is lay=0.
    // Output:
    //     none.
    // Return:
    //     The module number representing the nearest module.
    Int_t    i,l,a,e,in=0;
    Double_t d,dn=1.0e10;
    Bool_t   t=lay!=0; // skip if lay = 0 default value check all layers.

    for(i=0;i<fNmodules;i++){
        if(t){GetModuleId(i,l,a,e);if(l!=lay) continue;}
        if((d=GetGeomMatrix(i)->Distance2(g))<dn){
            dn = d;
            in = i;
        } // end if
    } // end for i
    return in;
}
//______________________________________________________________________
void AliITSgeom::GetNearest27(const Double_t g[3],Int_t n[27],Int_t lay)const{
    //      Finds 27 Detectors (Modules) that are nearest the point g [cm] in
    // ALICE Global coordinates. If layer !=0 then the search is restricted
    // to Detectors (Modules) in that particular layer. The number 27 comes 
    // from including the nearest detector and all those around it (up, down,
    // left, right, forwards, backwards, and the corners).
    // Input:
    //     Double_t g[3]  The ALICE Cartesian global coordinate from which the
    //                    distance is to be calculated with.
    //     Int_t lay      The layer to restrict the search to. If layer=0 then
    //                    all layers are searched. Default is lay=0.
    // Output:
    //     Int_t n[27]    The module number representing the nearest 27 modules
    //                    in order.
    // Return:
    //     none.
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
            d = GetGeomMatrix(i)->Distance2(g);
            if(d<dn[a]){
                for(e=26;e>a;e--){dn[e] = dn[e-1];in[e] = in[e-1];}
                dn[a] = d; in[a] = i;
            } // end if d<dn[i]
        } // end for a
    } // end for i
    for(i=0;i<27;i++) n[i] = in[i];
}
//_______________________________________________________________________
void AliITSgeom::DetLToTrackingV2(Int_t md,Float_t xin,Float_t zin,
                                  Float_t &yout,Float_t &zout) const {

    //Conversion from local coordinates on detectors to local
    //coordinates used for tracking ("v2")
    // Inputs:
    //   Int_t   md      Module number
    //   Float_t xin     Standard local coordinate x
    //   Float_t zin     Standard local coordinate z
    // Output:
    //   Float_t yout    Tracking local coordinate y
    //   Float_t zout    Tracking local coordinate z
    // Return:
    //   none.
    Float_t x,y,z;
    Double_t rt[9],al;

    GetTrans(md,x,y,z);
    GetRotMatrix(md,rt);
    al = TMath::ATan2(rt[1],rt[0])+TMath::Pi();
    yout = -(-xin+(x*((Float_t)TMath::Cos(al))+y*((Float_t)TMath::Sin(al))));
    if(md<(GetModuleIndex(2,1,1))) yout *= -1; 
    zout = -zin+z; 
}
//_______________________________________________________________________
void AliITSgeom::TrackingV2ToDetL(Int_t md,Float_t yin,Float_t zin,
                                  Float_t &xout,Float_t &zout) const {
    //Conversion from local coordinates used for tracking ("v2") to
    //local detector coordinates  
    // Inputs:
    //   Int_t   md      Module number
    //   Float_t yin     Tracking local coordinate y
    //   Float_t zin     Tracking local coordinate z
    // Output:
    //   Float_t xout    Standard local coordinate x
    //   Float_t zout    Standard local coordinate z
    // Return:
    //   none.
    Float_t x,y,z;
    Double_t rt[9],al;

    GetTrans(md,x,y,z);
    GetRotMatrix(md,rt);
    al = TMath::ATan2(rt[1],rt[0])+TMath::Pi();
    xout = yin;
    if(md<(GetModuleIndex(2,1,1))) xout = -xout;
    xout += (x*((Float_t)TMath::Cos(al))+y*((Float_t)TMath::Sin(al)));
    zout  = -zin+z; 
}
//----------------------------------------------------------------------

