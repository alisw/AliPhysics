/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                          *
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
Revision 1.8  2000/10/02 16:32:57  barbera
Forward declarations added and formatting

Revision 1.3.4.8  2000/10/02 15:55:26  barbera
Forward declarations added and formatting

Revision 1.7  2000/09/22 12:36:38  nilsen
Minor changes to improve compilation and create less NOISE.

Revision 1.6  2000/07/10 16:07:18  fca
Release version of ITS code

Revision 1.3.4.2  2000/03/02 21:42:29  nilsen 
Linked AliDetector::fDigit to AliITSmodule::fDigitsM and AliITS::fITSRecPoints
to AliITSmodule::fRecPointsM. Renamed AliITSmodule::fPointsM to fRecPointsM.
Removed the deletion of fDigitsM from the distructor since it is only a copy
of what is in AliDetector. Fixed a bug in the functions LineSegmentL and
LineSegmentG. Added two new versions of LineSegmentL and LineSegmentG to
additionaly return track number from the hit. Removed FastPoint function,
haven't found anywere it was used, also it had very many problems, should
just call FastPointSPD... .

Revision 1.3.4.1  2000/01/12 19:03:32  nilsen
This is the version of the files after the merging done in December 1999.
See the ReadMe110100.txt file for details

Revision 1.3  1999/10/04 15:20:12  fca
Correct syntax accepted by g++ but not standard for static members, remove minor warnings

Revision 1.2  1999/09/29 09:24:20  fca
Introduction of the Copyright and cvs Log

*/

#include <TArrayI.h>

#include <stdlib.h>

#include "AliRun.h"
#include "AliITS.h"
#include "AliITShit.h"
#include "AliITSmodule.h"
#include "AliITSgeom.h"

ClassImp(AliITSmodule)

//_______________________________________________________________________
//
// Impementation of class AliITSmodule
//
// created by: A. Bouchm, W. Peryt, S. Radomski, P. Skowronski 
//             R. Barbers, B. Batyunia, B. S. Nilsen
// ver 1.0     CERN 16.09.1999
//_______________________________________________________________________
//________________________________________________________________________
// 
// Constructors and deconstructor
//________________________________________________________________________
//
AliITSmodule::AliITSmodule() {
  // constructor
    fIndex       = 0;
    fHitsM       = 0;
    fTrackIndex  = 0;
    fHitIndex    = 0;
    fITS         = 0;

}
//_________________________________________________________________________
AliITSmodule::AliITSmodule(Int_t index) {
  // constructor

    fIndex      = index;
    fHitsM      = new TObjArray();
    fTrackIndex = new TArrayI(16);
    fHitIndex   = new TArrayI(16);
    fITS        = (AliITS*)(gAlice->GetDetector("ITS"));
}
//__________________________________________________________________________
AliITSmodule::~AliITSmodule() {
    // The destructor for AliITSmodule. Before destoring AliITSmodule
    // we must first destroy all of it's members.

    fIndex   = 0;
    if(fHitsM){
	for(Int_t i=0;i<fHitsM->GetEntriesFast();i++) 
	    delete ((AliITShit *)(fHitsM->At(i)));
	// must delete each object in the TObjArray.
	delete fHitsM;
    } // end if
    delete fTrackIndex;
    delete fHitIndex;
    fITS = 0; // We don't delete this pointer since it is just a copy.
}
//____________________________________________________________________________
AliITSmodule::AliITSmodule(const AliITSmodule &source){
////////////////////////////////////////////////////////////////////////
//     Copy Constructor 
////////////////////////////////////////////////////////////////////////
  printf("AliITSmodule error: AliITSmodule class has not to be copied! Exit.\n");
  exit(1);
}

//_____________________________________________________________________________
AliITSmodule& AliITSmodule::operator=(const AliITSmodule &source){
////////////////////////////////////////////////////////////////////////
//    Assignment operator 
////////////////////////////////////////////////////////////////////////
  printf("AliITSmodule error: AliITSmodule class has not to be copied! Exit.\n");
  exit(1);
  return *this; // fake return neded on Sun
} 

//_________________________________________________________________________
// 
// Hits management
//__________________________________________________________________________
Int_t AliITSmodule::AddHit(AliITShit* hit,Int_t t,Int_t h) {
// Hits management

  //printf("AddHit: beginning hit %p t h %d %d\n",hit,t,h);
    fHitsM->AddLast(new AliITShit(*hit));
    Int_t fNhitsM = fHitsM->GetEntriesFast();
    if(fNhitsM-1>=fTrackIndex->GetSize()){ // need to expand the TArrayI
	TArrayI *p = new TArrayI(fNhitsM+64);
	for(Int_t i=0;i<fTrackIndex->GetSize();i++) 
	    (*p)[i] = fTrackIndex->At(i);
	delete fTrackIndex;
	fTrackIndex = p;
    } // end if
    if(fNhitsM-1>=fHitIndex->GetSize()){ // need to expand the TArrayI
	TArrayI *p = new TArrayI(fNhitsM+64);
	for(Int_t i=0;i<fHitIndex->GetSize();i++) 
	    (*p)[i] = fHitIndex->At(i);
	delete fHitIndex;
	fHitIndex = p;
    } // end if
    (*fTrackIndex)[fNhitsM-1] = t;
    (*fHitIndex)[fNhitsM-1]   = h;
    return fNhitsM;
}

//___________________________________________________________________________
void AliITSmodule::MedianHitG(Int_t index,
			      Float_t hitx1,Float_t hity1,Float_t hitz1,
			      Float_t hitx2,Float_t hity2,Float_t hitz2,
			      Float_t &xMg, Float_t &yMg, Float_t &zMg){
  // median hit
   AliITSgeom *gm = fITS->GetITSgeom();
   Float_t x1l,y1l,z1l;
   Float_t x2l,y2l,z2l;
   Float_t xMl,yMl=0,zMl;
   Float_t l[3], g[3];

   g[0] = hitx1;
   g[1] = hity1;
   g[2] = hitz1;
   gm->GtoL(index,g,l);
   x1l = l[0];
   y1l = l[1];
   z1l = l[2];

   g[0] = hitx2;
   g[1] = hity2;
   g[2] = hitz2;
   gm->GtoL(index,g,l);
   x2l = l[0];
   y2l = l[1];
   z2l = l[2];

   xMl = (-y1l / (y2l-y1l))*(x2l-x1l) + x1l;
   zMl = (-y1l / (y2l-y1l))*(z2l-z1l) + z1l;

   l[0] = xMl;
   l[1] = yMl;
   l[2] = zMl;
   gm->LtoG(index,l,g);
   xMg = g[0];
   yMg = g[1];
   zMg = g[2];
}
//___________________________________________________________________________
Double_t AliITSmodule::PathLength(Int_t index,AliITShit *itsHit1,
				  AliITShit *itsHit2){
  // path lenght
   Float_t  x1g,y1g,z1g;   
   Float_t  x2g,y2g,z2g;
   Double_t s;

   itsHit1->GetPositionG(x1g,y1g,z1g);
   itsHit2->GetPositionG(x2g,y2g,z2g);

   s = TMath::Sqrt( ((Double_t)(x2g-x1g)*(Double_t)(x2g-x1g)) +
		    ((Double_t)(y2g-y1g)*(Double_t)(y2g-y1g)) +
		    ((Double_t)(z2g-z1g)*(Double_t)(z2g-z1g))  );
   return s;
}
//___________________________________________________________________________
void AliITSmodule::PathLength(Int_t index,
			      Float_t x,Float_t y,Float_t z,
			      Int_t status,Int_t &nseg,
			      Float_t &x1,Float_t &y1,Float_t &z1,
			      Float_t &dx1,Float_t &dy1,Float_t &dz1,
			      Int_t &flag){
  // path length
    static Float_t x0,y0,z0;

    if (status == 66){ // entering
	x0 = x;
	y0 = y;
	z0 = z;
	nseg = 0;
	flag = 1;
    }else{
	x1 = x0;
	y1 = y0;
	z1 = z0;
	dx1 = x-x1;
	dy1 = y-y1;
	dz1 = z-z1;
	nseg++;
	if (status == 68) flag = 0; //exiting
	else flag = 2; //inside
	x0 = x;
	y0 = y;
	z0 = z;
    } // end if
}
//___________________________________________________________________________
Bool_t AliITSmodule::LineSegmentL(Int_t hitindex,Double_t &a,Double_t &b,
				  Double_t &c,Double_t &d,
				  Double_t &e,Double_t &f,Double_t &de){
  // line segment
    static Int_t hitindex0;
    AliITShit *h0,*h1;

    if(hitindex>= fHitsM->GetEntriesFast()) return kFALSE;

    h1 = (AliITShit *) (fHitsM->At(hitindex));
    if(h1->StatusEntering()){ // if track entering volume, get index for next
	                      // step
	hitindex0 = hitindex;
	return kFALSE;
    } // end if StatusEntering()
    // else stepping
    h0 = (AliITShit *) (fHitsM->At(hitindex0));
    de = h1->GetIonization();
    h0->GetPositionL(a,c,e);
    h1->GetPositionL(b,d,f);
    b = b - a;
    d = d - c;
    f = f - e;
    hitindex0 = hitindex;
    return kTRUE;
}
//___________________________________________________________________________
Bool_t AliITSmodule::LineSegmentG(Int_t hitindex,Double_t &a,Double_t &b,
				  Double_t &c,Double_t &d,
				  Double_t &e,Double_t &f,Double_t &de){
  // line segment
    static Int_t hitindex0;
    AliITShit *h0,*h1;

    if(hitindex>= fHitsM->GetEntriesFast()) return kFALSE;

    h1 = (AliITShit *) (fHitsM->At(hitindex));
    if(h1->StatusEntering()){ // if track entering volume, get index for next
	                      // step
	hitindex0 = hitindex;
	return kFALSE;
    } // end if StatusEntering()
    // else stepping
    h0 = (AliITShit *) (fHitsM->At(hitindex0));
    de = h1->GetIonization();
    h0->GetPositionG(a,c,e);
    h1->GetPositionG(b,d,f);
    b = b - a;
    d = d - c;
    f = f - e;
    hitindex0 = hitindex;
    return kTRUE;
}
//___________________________________________________________________________
Bool_t AliITSmodule::LineSegmentL(Int_t hitindex,Double_t &a,Double_t &b,
				  Double_t &c,Double_t &d,
				  Double_t &e,Double_t &f,
				  Double_t &de,Int_t &track){
  // line segmente
    static Int_t hitindex0;
    AliITShit *h0,*h1;

    if(hitindex>= fHitsM->GetEntriesFast()) return kFALSE;

    h1 = (AliITShit *) (fHitsM->At(hitindex));
    if(h1->StatusEntering()){ // if track entering volume, get index for next
	                      // step
	hitindex0 = hitindex;
	track = h1->GetTrack();
	return kFALSE;
    } // end if StatusEntering()
    // else stepping
    h0 = (AliITShit *) (fHitsM->At(hitindex0));
    de = h1->GetIonization();
    h0->GetPositionL(a,c,e);
    h1->GetPositionL(b,d,f);
    b = b - a;
    d = d - c;
    f = f - e;
    hitindex0 = hitindex;
    track = h1->GetTrack();
    return kTRUE;
}
//___________________________________________________________________________
Bool_t AliITSmodule::LineSegmentG(Int_t hitindex,Double_t &a,Double_t &b,
				  Double_t &c,Double_t &d,
				  Double_t &e,Double_t &f,
				  Double_t &de,Int_t &track){
  // line segment
    static Int_t hitindex0;
    AliITShit *h0,*h1;

    if(hitindex>= fHitsM->GetEntriesFast()) return kFALSE;

    h1 = (AliITShit *) (fHitsM->At(hitindex));
    if(h1->StatusEntering()){ // if track entering volume, get index for next
	                      // step
	hitindex0 = hitindex;
	track = h1->GetTrack();
	return kFALSE;
    } // end if StatusEntering()
    // else stepping
    h0 = (AliITShit *) (fHitsM->At(hitindex0));
    de = h1->GetIonization();
    h0->GetPositionG(a,c,e);
    h1->GetPositionG(b,d,f);
    b = b - a;
    d = d - c;
    f = f - e;
    hitindex0 = hitindex;
    track = h1->GetTrack();
    return kTRUE;
}
//___________________________________________________________________________
void AliITSmodule::MedianHitL(Int_t index, 
               		     AliITShit *itsHit1, 
	       	    	     AliITShit *itsHit2, 
	       		     Float_t &xMl, Float_t &yMl, Float_t &zMl){
  // median hit
   Float_t x1l,y1l,z1l;
   Float_t x2l,y2l,z2l;

   itsHit1->GetPositionL(x1l,y1l,z1l);
   itsHit2->GetPositionL(x2l,y2l,z2l);

   xMl = (-y1l / (y2l-y1l))*(x2l-x1l) + x1l;
   yMl = 0.0;
   zMl = (-y1l / (y2l-y1l))*(z2l-z1l) + z1l;	     
}
//___________________________________________________________________________
void AliITSmodule::MedianHit(Int_t index,
			     Float_t xg,Float_t yg,Float_t zg,
			     Int_t status,
			     Float_t &xMg,Float_t &yMg,Float_t &zMg,
			     Int_t &flag){
  // median hit
   static Float_t x1,y1,z1;

   if (status == 66){ // entering
       x1 = xg;
       y1 = yg;
       z1 = zg;
       flag = 1;
   } // end if
   if (status == 68){ // exiting
       MedianHitG(index,x1,y1,z1,xg,yg,zg,xMg,yMg,zMg);
       flag = 0;
   } // end if
   if ((status != 66) && (status != 68)) flag = 1;
}
//___________________________________________________________________________
void AliITSmodule::GetID(Int_t &lay,Int_t &lad,Int_t &det){
  // get ID
	fITS->GetITSgeom()->GetModuleId(fIndex,lay,lad,det);
	return ;
}

