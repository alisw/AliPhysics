/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/* $Id$ */
/* $Author$ */
/* $Date$ */
/* $Name$ */
/* $Header$ */
/*
   $Log$
   Revision 1.1.2.2  2000/06/04 16:35:37  nilsen
   One more try to fix log comments.

   Revision 1.1.2.1  2000/03/02 20:13:52  nilsen
   A new class useful for ITS detector alignment studdies.
 */
/* $Revision$ */

// Standard Root Libraries
#include <TMath.h>

// ITS libraries
#include "AliITSgeom.h"
#include "AliITSAlignmentTrack.h"
#include "AliITSstatistics.h"
#include "AliITSstatistics2.h"

ClassImp(AliITSAlignmentTrack)

//______________________________________________________________________
AliITSAlignmentTrack::AliITSAlignmentTrack(){

    ftrack=fnclust=0,fnclustMax=0;
    ffunc=-1;
    fclust=0;
    for(Int_t i=0;i<10;i++) fpar[i]=0.0;
    fpx=fpy=fpz=fp=fpt=0.0;
    fChi2=-1.0;
}
//______________________________________________________________________
AliITSAlignmentTrack::~AliITSAlignmentTrack(){

    ftrack=fnclust=0,fnclustMax=0;
    ffunc=-1;
    delete[] fclust;
    fclust=0;
    for(Int_t i=0;i<10;i++) fpar[i]=0.0;
    fpx=fpy=fpz=fp=fpt=0.0;
    fChi2=-1.0;
}
//______________________________________________________________________
void AliITSAlignmentTrack::func0(Double_t *go,Double_t *gi){
    Double_t x,y,z;

    x = gi[0];
    y = gi[1];
    z = gi[2];
    x = fpar[0]+fpar[1]*z;
    y = fpar[2]+fpar[3]*z;
    go[0] = x;
    go[1] = y;
    go[2] = z;
    return;
}
//______________________________________________________________________
void AliITSAlignmentTrack::func1(Double_t *go,Double_t *gi){
    Double_t x,y,z;

    x = gi[0];
    y = gi[1];
    z = gi[2];
    x = fpar[0]+fpar[1]*y;
    z = fpar[2]+fpar[3]*y;
    go[0] = x;
    go[1] = y;
    go[2] = z;
    return;
}
//______________________________________________________________________
void AliITSAlignmentTrack::func2(Double_t *go,Double_t *gi){
    Double_t x,y,z,r,th;

    x = gi[0];
    y = gi[1];
    z = gi[2];
    th = TMath::ATan2(y-fpar[1],x-fpar[0]);
    r  = TMath::Hypot(x-fpar[0],y-fpar[1]);
    if(th<0.0) th += 2.0*TMath::Pi();
    x = fpar[0]+fpar[2]*TMath::Cos(th);
    y = fpar[1]+fpar[2]*TMath::Sin(th);
    z = fpar[3]+fpar[4]*r;
    go[0] = x;
    go[1] = y;
    go[2] = z;
    return;
}
//______________________________________________________________________
void AliITSAlignmentTrack::func(Double_t *go,Double_t *gi){

    switch (ffunc){
    case 0:
	func0(go,gi);
	return;
    case 1:
	func1(go,gi);
	return;
    case 2:
	func2(go,gi);
	return;
    } // end switch
}
//______________________________________________________________________
Double_t AliITSAlignmentTrack::ComputeChi2(){
    Int_t    i,j,k,l;
    Double_t chi2=0.0,go[3],gi[3];

    switch (ffunc) {
    case -1:
	return -1.0;
	break;
    case 0:
	for(i=0;i<fnclust;i++){
	    for(j=0;j<3;j++)gi[j] = fclust[i].fxg[j];
	    func0(go,gi);
	    for(k=0;k<2;k++) for(l=0;l<2;l++) {
		//if(k==2 || l==2) continue;
		chi2 += (go[k] - gi[k])*fclust[i].fExg[k][l]*(go[l] - gi[l]);
	    } // end for k,l
	} // end for i
	break;
    case 1:
	for(i=0;i<fnclust;i++){
	    for(j=0;j<3;j++)gi[j] = fclust[i].fxg[j];
	    func1(go,gi);
	    for(k=0;k<3;k++) for(l=0;l<3;l++){
		if(k==1 || l==1) continue;
		chi2 += (go[k] - gi[k])*fclust[i].fExg[k][l]*(go[l] - gi[l]);
	    } // end for k,l
	} // end for i
	break;
    case 2:
	for(i=0;i<fnclust;i++){
	    for(j=0;j<3;j++)gi[j] = fclust[i].fxg[j];
	    func0(go,gi);
	    for(k=0;k<3;k++) for(l=0;l<3;l++){
		chi2 += (go[k] - gi[k])*fclust[i].fExg[k][l]*(go[l] - gi[l]);
	    } // end for k,l
	} // end for i
	break;
    } // end switch

    fChi2 = (Float_t) chi2;

    return  chi2;
}
//______________________________________________________________________
Int_t AliITSAlignmentTrack::FindCircleCenter(Double_t *xc,Double_t *x1,
		                            Double_t  *x2,Double_t *x3){
////////////////////////////////////////////////////////////////////////
//     This was derived as folows. Given three non-linear points find
// the circle that is therefor defined by those three non-linar points.
// Assume that the circle is centers at xc,yc and has a radous R. Then
// (1) R^2 = (x1-xc)^2 + (y1-yc)^2
// (2) R^2 = (x2-xc)^2 + (y2-yc)^2
// (3) R^2 = (x3-xc)^2 + (y3-yc)^2.
// Now consider the two equations derived from the above
// (1) - (2) = x1^2 - x2^2 -2xc(x1-x2) + y1^2 - y2y2 -2yc(y1-y2) = 0
// (1) - (3) = x1^2 - x3^2 -2xc(x1-x3) + y1^2 - y3y2 -2yc(y1-y3) = 0
// solving these two equations for x0 and y0 gives
// xc = +{(y1-y2)*(y1-y3)*(y1-y3)+x1*x1*(y1-y3)+x2*x2*(y3-y1)+x3*x3*(y1-y2)}/2d
// yc = -{(x1-x2)*(x1-x3)*(x1-x3)+y1*y1*(x1-x3)+y2*y2*(x3-x1)+y3*y3*(x1-x2)}/2d
// with d = (x1-x2)*(y1-y3) - (x1-x3)*(y1-y2)
////////////////////////////////////////////////////////////////////////
    Double_t d;

    d = (x1[0]-x2[0])*(x1[1]-x3[1]) - (x1[0]-x3[0])*(x1[1]-x2[1]);
    if(d==0.0) return 0;  // fits to a line!

    xc[0] = (x1[1]-x2[1])*(x1[1]-x3[1])*(x1[1]-x3[1])+
	     x1[0]*x1[0]*(x1[1]-x3[1])+
             x2[0]*x2[0]*(x3[1]-x1[1])+
             x3[0]*x3[0]*(x1[1]-x2[1]);
    xc[1] = (x1[0]-x2[0])*(x1[0]-x3[0])*(x1[0]-x3[0])+
             x1[1]*x1[1]*(x1[0]-x3[0])+
	     x2[1]*x2[1]*(x3[0]-x1[0])+
             x3[1]*x3[1]*(x1[0]-x2[0]);
    xc[0] = +0.5*xc[0]/d;
    xc[1] = -0.5*xc[1]/d;

    return 1;
}
//______________________________________________________________________
Int_t AliITSAlignmentTrack::FitTrackToLineG(){
// Xg = fpar[0] + fpar[1] * Zg
// Yg = fpar[2] + fpar[3] * Zg;
// Local Variables
    Int_t   i;
    Double_t x,y,z,wx,wy;
    Double_t a,b,c,d;
    AliITSstatistics2 *sx = new AliITSstatistics2(4);
    AliITSstatistics2 *sy = new AliITSstatistics2(4);
    Double_t b0,d0;

    fChi2 = -1.0;
    ffunc = 0;
    if(fnclust<3) return -1;

    b  = (fclust[0].fxg[0]-fclust[fnclust].fxg[0])/
         (fclust[0].fxg[2]-fclust[fnclust].fxg[2]);
    d  = (fclust[0].fxg[1]-fclust[fnclust].fxg[1])/
         (fclust[0].fxg[2]-fclust[fnclust].fxg[2]);
    do{
	b0 = b;
	d0 = d;
	for(i=0;i<fnclust;i++){
	    sx->Reset();
	    sy->Reset();
	    x    = fclust[i].fxg[0];
	    y    = fclust[i].fxg[1];
	    z    = fclust[i].fxg[2];
	    wx   = 1./(1./fclust[i].fExg[0][0] +
		       (1./fclust[i].fExg[2][2])*b0*b0);// 1.0/rms^2
	    wy   = 1./(1./fclust[i].fExg[1][1] +
		       (1./fclust[i].fExg[2][2])*d0*d0);// 1.0/rms^2
	    sx->AddValue(x,z,wx);
	    sy->AddValue(y,z,wy);
       } // end for i
	fChi2  = sx->FitToLine(a,b);
	fChi2 += sy->FitToLine(c,d);
	//} while(fabs(b0-b)<1.E-5 && fabs(d0-d)<1.E-5);
    } while(TMath::Abs(b0-b)<1.E-5 && TMath::Abs(d0-d)<1.E-5);
    fpar[0] = a;
    fpar[1] = b;
    fpar[2] = c;
    fpar[3] = d;
    return 0;
}
//______________________________________________________________________
Int_t AliITSAlignmentTrack::FitTrackToLineL(AliITSgeom *gm){
// X = fpar[0] + fpar[1] * y;
// Z = fpar[2] + fpar[3] * y;
// in the local coordinate system of the detector fclust[0].findex.
   // Local Variables
   Int_t   i,j,k;
   Double_t wx/*,wy*/,wz,Exll[3][3];
   Double_t a,b,c,d;
   Double_t xg[3],xl[3],x2g[3],x2l[3];
   AliITSstatistics2 *Fx  = new AliITSstatistics2(2);
   AliITSstatistics2 *Fz  = new AliITSstatistics2(2);

   fChi2 = -1.0;
   ffunc = 1;
   if(fnclust<3) return -1;

   Int_t Npts = fnclust;
   for(i=0;i<Npts;i++){
       for(j=0;j<3;j++)for(k=0;k<3;k++) Exll[j][k] = fclust[i].fExl[j][k];
       for(j=0;j<3;j++){x2l[j] = fclust[i].fxl[j];}
       gm->LtoL(fclust[i].findex,fclust[0].findex,x2l,xl);
       gm->LtoLErrorMatrix(fclust[i].findex,fclust[0].findex,
			   (Double_t **) fclust[i].fExl,(Double_t **) Exll);
       wx = Exll[0][0];
       wz = Exll[2][2];
       Fx->AddValue(xl[0],xl[1],wx);
       Fz->AddValue(xl[2],xl[1],wz);
   } // end for i
   fChi2  = Fx->FitToLine(a,b);
   fChi2 += Fz->FitToLine(c,d);
   fpar[0] = a;
   fpar[1] = b;
   fpar[2] = c;
   fpar[3] = d;
   // convert to global if posible.
   xl[0]  = a;
   xl[1]  = 0.0;
   xl[2]  = c;
   x2l[0] = a+b;
   x2l[1] = 1.0;
   x2l[2] = c+d;
   gm->LtoG(fclust[0].findex,xl,xg);
   gm->LtoG(fclust[0].findex,x2l,x2g);
   c = xg[2] - x2g[2];
   if(c!=0.0){
       b = (xg[0] - x2g[0])/c;
       d = (xg[1] - x2g[1])/c;
       a = xg[0] - b*xg[2];
       c = xg[1] - d*xg[2];
       fpar[4] = a;
       fpar[5] = b;
       fpar[6] = c;
       fpar[7] = d;
   }else{
       fpar[4] = 0.0;
       fpar[5] = 0.0;
       fpar[6] = 0.0;
       fpar[7] = 0.0;
       return -1;
   }// end if c!=0.0
   return 0;
}
//______________________________________________________________________
void AliITSAlignmentTrack::FitToFunction(Int_t n,AliITSgeom *gm){
    Int_t i,j,k,l;
    Double_t r,w,xc[3],x1[3],x2[3],x3[3];
    AliITSstatistics2 *sa = new AliITSstatistics2(4);
    AliITSstatistics2 *sb = new AliITSstatistics2(4);

    ffunc = -1;
    if(fnclust<3) return;

    switch (n){
    case -1:
	return;
	break;
    case 0:
	ffunc = 0;
	FitTrackToLineG();
	return;
	break;
    case 1:
	ffunc = 1;
	FitTrackToLineL(gm);
	return;
	break;
    case 2:
	ffunc = 2;
	sa->Reset();
	sb->Reset();
	for(i=0;i<fnclust;i++){
	    for(l=0;l<3;l++) x1[l] = (Double_t) fclust[i].fxg[l];
	    r = TMath::Hypot((Double_t)x1[0],(Double_t)x1[1]);
	    w = (Double_t) (fclust[i].fExg[2][2]);
	    sb->AddValue((Double_t)fclust[i].fxg[2],r,w);
	    if(i<fnclust-2)for(j=i+1;j<fnclust-1;j++)for(k=j+1;k<fnclust;k++){
		for(l=0;l<3;l++) x2[l] = (Double_t) fclust[j].fxg[l];
		for(l=0;l<3;l++) x3[l] = (Double_t) fclust[k].fxg[l];
		w = 0.0;
		w += TMath::Hypot((Double_t)fclust[i].fExg[0][0],
				  (Double_t)fclust[i].fExg[1][1]);
		w += TMath::Hypot(w,
                           TMath::Hypot((Double_t)fclust[j].fExg[0][0],
                                        (Double_t)fclust[j].fExg[1][1]));
		w += TMath::Hypot(w,
                           TMath::Hypot((Double_t)fclust[k].fExg[0][0],
                                        (Double_t)fclust[k].fExg[1][1]));
		if(FindCircleCenter(xc,x1,x2,x3)==0){ // Can't find circle
		    FitToFunction(1,gm); // fit to line.
		    return;
		} // end if
		sa->AddValue(xc[1],xc[0],w);
	    } // end for j,k
	} // end for i
	fpar[0] = sa->GetMeanX();
	fpar[1] = sa->GetMeanY();
	fpar[2] = sb->GetMeanX();
	fChi2   = sb->FitToLine(fpar[3],fpar[4]);
	fChi2   = (Float_t) ComputeChi2();
	return;
	break;
    default:
	return;
	break;
    } // end switch
    return;
}
//______________________________________________________________________
void AliITSAlignmentTrack::Streamer(TBuffer &R__b){
   // Stream an object of class AliITSAlignmentTrack.

   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(); if (R__v) { }
      TObject::Streamer(R__b);
      R__b >> ftrack;
      R__b >> fnclust;
      R__b >> fnclustMax;
//      R__b >> fclust;
      R__b >> ffunc;
      R__b.ReadStaticArray(fpar);
      R__b >> fChi2;
      R__b >> fpx;
      R__b >> fpy;
      R__b >> fpz;
      R__b >> fp;
      R__b >> fpt;
   } else {
      R__b.WriteVersion(AliITSAlignmentTrack::IsA());
      TObject::Streamer(R__b);
      R__b << ftrack;
      R__b << fnclust;
      R__b << fnclustMax;
//      R__b << fclust;
      R__b << ffunc;
      R__b.WriteArray(fpar, 10);
      R__b << fChi2;
      R__b << fpx;
      R__b << fpy;
      R__b << fpz;
      R__b << fp;
      R__b << fpt;
   }
}
