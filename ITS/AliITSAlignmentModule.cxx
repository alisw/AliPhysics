/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/* $Id$ */
/* $Author$ */
/* $Date$ */
/* $Name$ */
/* $Header$ */
/*
   $Log$
   Revision 1.1.2.3  2000/06/04 16:35:09  nilsen
   One more try to fix log comments.

   Revision 1.1.2.2  2000/03/04 23:39:36  nilsen
   Fixed the logs???
 */
/* Revision 1.1.2.1  2000/03/02 20:12:23  nilsen */
/* A new class usefull for ITS detector Alignment Studdies */
/* */
/* $Revision$ */

// Standard C & C++ libraries
#include <TObject.h>

// Standard Root Libraries
#include <TParticle.h>

// ITS libraries
#include "AliITSgeom.h"
#include "AliITSAlignmentTrack.h"
#include "AliITSAlignmentModule.h"

ClassImp(AliITSAlignmentModule)

//______________________________________________________________________
AliITSAlignmentModule::AliITSAlignmentModule(){

    findex = -1;
    flay   = -1;
    flad   = -1;
    fdet   = -1;
    fx0[0] = 0.0;
    fx0[1] = 0.0;
    fx0[2] = 0.0;
    fangles[0] = 0.0;
    fangles[1] = 0.0;
    fangles[2] = 0.0;
    fM[0][0]   = 0.0;
    fM[0][1]   = 0.0;
    fM[0][2]   = 0.0;
    fM[1][0]   = 0.0;
    fM[1][1]   = 0.0;
    fM[1][2]   = 0.0;
    fM[2][0]   = 0.0;
    fM[2][1]   = 0.0;
    fM[2][2]   = 0.0;
    ftrksM =0;
    fChi2 = 0.0;
}
//______________________________________________________________________
AliITSAlignmentModule::AliITSAlignmentModule(Int_t index,AliITSgeom *gm,
                                         Int_t ntrk,AliITSAlignmentTrack *trk){
    Float_t x,y,z,n;
    Int_t   i,j;

    findex = index;
    gm->GetModuleId(index,flay,flad,fdet);
    gm->GetRotMatrix(index,(Double_t *) &(fM[0][0]));
    gm->GetTrans(flay,flad,fdet,x,y,z);
    fx0[0] = (Double_t )x;
    fx0[1] = (Double_t )y;
    fx0[2] = (Double_t )z;
    gm->GetAngles(flay,flad,fdet,x,y,z);
    fangles[0] = (Double_t )x;
    fangles[1] = (Double_t )y;
    fangles[2] = (Double_t )z;
    ftrksM = new TObjArray();
    fChi2 = 0.0;

    for(i=0;i<ntrk;i++){
	n = 0;
	for(j=0;j<trk[i].GetNumberOfClustersSl();j++){
	    if(trk[i].GetIndex(j,index)==findex){
		ftrksM->AddAtFree((TObject *) &trk[i]);
		break; // break out of the j loop
	    } // end if
	} // end for j
    } // end for i
}
//______________________________________________________________________
Double_t AliITSAlignmentModule::ComputeChi2(){
    Float_t n;
    Int_t   i,j,ntrk,ntr;

    ntrk = ftrksM->GetEntriesFast();
    fChi2 = 0.0;
    for(i=0;i<ntrk;i++){
	n = 0;
	ntr =((AliITSAlignmentTrack *)(ftrksM->At(i)))->GetNumberOfClustersSl();
	for(j=0;j<ntr;j++){
		fChi2 += (Double_t) ((AliITSAlignmentTrack *)(ftrksM->At(i)))->
		                     GetChi2();
		n++;
	} // end for j
	if(n==0){
	    fChi2 = -1.0;
	}else{
	    fChi2 /= (Double_t) n;
	} // end if n==0
    } // end for i
    return fChi2;
}
//______________________________________________________________________
void AliITSAlignmentModule::lnsrch(Int_t npar,Double_t *xold,Double_t fold,
				   Double_t *g,Double_t *p,Double_t *x,
				   Double_t &f, Double_t stpmax,Int_t &check){
    Double_t ALF = 1.0e-4, TOLX = 1.0E-7;

    Int_t    i;
    Double_t a,alam,alam2=0.0,alamin,b,disc,f2=0.0,rhs1,rhs2,slope,sum,temp,
	     test,tmplam;

    check = 0;
    for(sum=0.0,i=0;i<npar;i++) sum += p[i]*p[i];
    sum = TMath::Sqrt(sum);
    if(sum>stpmax) for(i=0;i<npar;i++) p[i] *= stpmax/sum;
    for(slope=0.0,i=0;i<npar;i++) slope += g[i]*p[i];
    if(slope >=0.0) printf("Error: round off problem in lnsrch.\n");
    test = 0.0;
    for(i=0;i<npar;i++){
	temp = TMath::Abs(p[i])/TMath::Max(TMath::Abs(xold[i]),1.0);
	if(temp > test) test = temp;
    } // end for i
    alamin = TOLX/test;
    alam = 1.0;
    for(;;){
	for(i=0;i<npar;i++) x[i] = xold[i] + alam*p[i];
	f = Chi2(x);
	if(alam < alamin){
	    for(i=0;i<npar;i++) x[i] = xold[i];
	    check = 1;
	    return;
	}else if(f <= fold+ALF*alam*slope) return;
	else{
	    if(alam == 1.0) tmplam = -slope/(2.0*(f-fold-slope));
	    else{
		rhs1 = f-fold-alam*slope;
		rhs2 = f2-fold-alam2*slope;
		a = (rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
		b = (-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam*alam);
		if(a==0.0) tmplam = -slope/(2.0*b);
		else{
		    disc = b*b - 3.0*a*slope;
		    if(disc<0.0) tmplam = 0.5*alam;
		    else if(b<=0.0) tmplam = (-b+TMath::Sqrt(disc))/(3.0*a);
		    else tmplam = -slope/(b+TMath::Sqrt(disc));
		} // end if a == 0.0
		if(tmplam > 0.5*alam) tmplam = 0.5*alam;
	    } // end if alam == 1.0
	} // end if alam < alamin
	alam2 = alam;
	f2 = f;
	alam = TMath::Max(tmplam,0.1*alam);
    } // end for ever loop

}
//______________________________________________________________________
void AliITSAlignmentModule::MRVMminimization(Int_t npar,Double_t *p,
					     Double_t &fret,Double_t gtol,
					     Int_t &iter){
    Int_t ITMAX = 200;
    Double_t EPS = 3.0e-8, TOLX = 4.0*EPS, STPMX = 100.0;

    Int_t    check,i,its,j;
    Double_t den,fac,fad,fae,fp,stpmax,sum=0.0,sumdg,sumxi,temp,test;
    Double_t *dg,*g,*hdg,**hessin,*pnew,*xi;

    // allocate space.
    dg     = new Double_t[npar];
    g      = new Double_t[npar];
    hdg    = new Double_t[npar];
    hessin = new Double_t * [npar];
    for(i=0;i<npar;i++) hessin[i] = new Double_t[npar];
    pnew   = new Double_t[npar];
    xi     = new Double_t[npar];

    // init function values
    fp = Chi2(p);
    dChi2(p,g);

    for(i=0;i<npar;i++){
	for(j=0;j<npar;j++) hessin[i][j] = 0.0;
	hessin[i][i] = 1.0;
	xi[i] = -g[i];
	sum += p[i]*p[i];
    } // end for i
    stpmax = STPMX*TMath::Max(TMath::Sqrt(sum),(Double_t)(npar+1));

    for(its=0;its<ITMAX;its++){
	iter = its;
	lnsrch(npar,p,fp,g,xi,pnew,fret,stpmax,check);
	fp = fret;
	for(i=0;i<npar;i++){
	    xi[i] = pnew[i] - p[i];
	    p[i] = pnew[i];
	} // end for i
	test = 0.0;
	for(i=0;i<npar;i++){
	    temp = TMath::Abs(xi[i])/TMath::Max(TMath::Abs(p[i]),1.0);
	    if(temp > test) test = temp;
	} // end for i
	if(test<TOLX) break;
	for(i=0;i<npar;i++) dg[i] = g[i];
	dChi2(p,g);
	test = 0.0;
	den = TMath::Max(fret,1.0);
	for(i=0;i<npar;i++){
	    temp = TMath::Abs(g[i])*TMath::Max(TMath::Abs(p[i]),1.0)/den;
	    if(temp>test) test=temp;
	} // end for i
	if(test<gtol) break;
	for(i=0;i<npar;i++) dg[i] = g[i] - dg[i];
	for(i=0;i<npar;i++){
	    hdg[i] = 0.0;
	    for(j=0;j<npar;j++) hdg[i] += hessin[i][j]*dg[j];
	} // end for i
	fac = fae = sumdg = sumxi = 0.0;
	for(i=0;i<npar;i++){
	    fac += dg[i]*xi[i];
	    fae += dg[i]*hdg[i];
	    sumdg += TMath::Sqrt(dg[i]);
	    sumxi += TMath::Sqrt(xi[i]);
	} // end for i
	if(fac>TMath::Sqrt(EPS*sumdg*sumxi)){
	    fac = 1.0/fac;
	    fad = 1.0/fae;
	    for(i=0;i<npar;i++) dg[i] = fac*xi[i]-fad*hdg[i];
	    for(i=0;i<npar;i++) for(j=i;j<npar;j++){
		hessin[i][j] += fac*xi[i]*xi[j] - fad*hdg[i]*hdg[j] + 
		                                  fae*dg[i]*dg[j];
		hessin[j][i] = hessin[i][j];
	    } // end for i,j
	} // end if fac>...
	for(i=0;i<npar;i++){
	    xi[i] = 0.0;
	    for(j=0;j<npar;j++) xi[i] -= hessin[i][j]*g[j];
	} // end for i
    } // end for its
    if(its==ITMAX) printf("Error: too many iterations\n");
    delete[] dg;
    delete[] g;
    delete[] hdg;
    delete[] pnew;
    delete[] xi;
    for(i=0;i<npar;i++) delete[] hessin[i];
    delete[] hessin;
}
//______________________________________________________________________
void AliITSAlignmentModule::SetByAngles(Double_t *th){
   Double_t  sx,cx,sy,cy,sz,cz;

   sx = TMath::Sin(th[0]); cx = TMath::Cos(th[0]);
   sy = TMath::Sin(th[1]); cy = TMath::Cos(th[1]);
   sz = TMath::Sin(th[2]); cz = TMath::Cos(th[2]);
   for(Int_t i=0;i<3;i++) fangles[i]   = th[i];
   fM[0][0] =  cz*cy;
   fM[0][1] = -cz*sy*sx - sz*cx;
   fM[0][2] = -cz*sy*cx + sz*sx;
   fM[1][0] =  sz*cy;
   fM[1][1] = -sz*sy*sx + cz*cx;
   fM[1][2] = -sz*sy*cx - cz*sx;
   fM[2][0] =  sy;
   fM[2][1] =  cy*sx;
   fM[2][2] =  cy*cx;
}
//______________________________________________________________________
void AliITSAlignmentModule::dfMdthx(Double_t dfMx[3][3]){
   Double_t  sx,cx,sy,cy,sz,cz;

   sx = TMath::Sin(fangles[0]); cx = TMath::Cos(fangles[0]);
   sy = TMath::Sin(fangles[1]); cy = TMath::Cos(fangles[1]);
   sz = TMath::Sin(fangles[2]); cz = TMath::Cos(fangles[2]);
   dfMx[0][0] =  0.0;
   dfMx[0][1] = -cz*sy*cx + sz*sx;
   dfMx[0][2] =  cz*sy*sx + sz*cx;
   dfMx[1][0] =  0.0;
   dfMx[1][1] = -sz*sy*cx - cz*sx;
   dfMx[1][2] =  sz*sy*sx - cz*cx;
   dfMx[2][0] =  0.0;
   dfMx[2][1] =  cy*cx;
   dfMx[2][2] = -cy*sx;
}
//______________________________________________________________________
void AliITSAlignmentModule::dfMdthy(Double_t dfMx[3][3]){
   Double_t  sx,cx,sy,cy,sz,cz;

   sx = TMath::Sin(fangles[0]); cx = TMath::Cos(fangles[0]);
   sy = TMath::Sin(fangles[1]); cy = TMath::Cos(fangles[1]);
   sz = TMath::Sin(fangles[2]); cz = TMath::Cos(fangles[2]);
   dfMx[0][0] = -cz*sy;
   dfMx[0][1] = -cz*cy*sx;
   dfMx[0][2] = -cz*cy*cx;
   dfMx[1][0] = -sz*sy;
   dfMx[1][1] = -sz*cy*sx;
   dfMx[1][2] = -sz*cy*cx;
   dfMx[2][0] =  cy;
   dfMx[2][1] = -sy*sx;
   dfMx[2][2] = -sy*cx;
}
//______________________________________________________________________
void AliITSAlignmentModule::dfMdthz(Double_t dfMx[3][3]){
   Double_t  sx,cx,sy,cy,sz,cz;

   sx = TMath::Sin(fangles[0]); cx = TMath::Cos(fangles[0]);
   sy = TMath::Sin(fangles[1]); cy = TMath::Cos(fangles[1]);
   sz = TMath::Sin(fangles[2]); cz = TMath::Cos(fangles[2]);
   dfMx[0][0] = -sz*cy;
   dfMx[0][1] =  sz*sy*sx - cz*cx;
   dfMx[0][2] =  sz*sy*cx + cz*sx;
   dfMx[1][0] =  cz*cy;
   dfMx[1][1] = -cz*sy*sx - sz*cx;
   dfMx[1][2] = -cz*sy*cx + sz*sx;
   dfMx[2][0] =  0.0;
   dfMx[2][1] =  0.0;
   dfMx[2][2] =  0.0;
}
//______________________________________________________________________
void AliITSAlignmentModule::LtoG(Double_t xl[],Double_t xg[]){

    Int_t i;
    for(i=0;i<3;i++) xg[i] = fx0[i];
    for(i=0;i<3;i++)for(Int_t j=0;j<3;j++) xg[i] += fM[j][i]*xl[j];
}
//______________________________________________________________________
void AliITSAlignmentModule::GtoL(Double_t xg[],Double_t xl[]){
    Int_t i;
    for(i=0;i<3;i++) xl[i] = 0.0;
    for(i=0;i<3;i++)for(Int_t j=0;j<3;j++)xl[i]+=fM[i][j]*(xg[j]-fx0[j]);
}
//______________________________________________________________________
Double_t AliITSAlignmentModule::Chi2(Double_t p[]){
    Int_t    i,j,k,l,n=0,indx;
    Double_t chi=0.0,xo[3],xi[3],xg[3],Ex[3][3];
    AliITSAlignmentTrack *tr;

    for(i=0;i<3;i++) fx0[i] = p[i];
    SetByAngles(&(p[3]));

    for(i=0;i<ftrksM->GetEntriesFast();i++) {
	tr = (AliITSAlignmentTrack *)(ftrksM->At(i));
	for(j=0;j<tr->GetNumberOfClustersSl();j++){
	    tr->GetIndex(j,indx);
	    if(indx==findex){
		n++;
		tr->GetPointG(j,(Double_t *)xi);
		tr->func(xi,xo);
		tr->GetPointL(j,(Double_t *)xi);
		tr->GetErrorG(j,(Double_t **)Ex);
		LtoG(xi,xg);
		for(k=0;k<3;k++)for(l=0;l<3;l++)
		    chi += (xg[k] - xo[k])*Ex[k][l]*(xg[l]-xo[l]);
	    } // end if indx==findex
	} // end for j
    } // end for i
    if(n<7) return chi;
    return chi/(Double_t)(n-6);
}
//______________________________________________________________________
void AliITSAlignmentModule::dChi2(Double_t p[],Double_t dChi2[]){
    Int_t    i,j,k,l,m,n=0,indx;
    Double_t chi[6]={0.0,0.0,0.0,0.0,0.0,0.0},xo[3],xi[3],xg[3],Ex[3][3];
    Double_t dxdp[3][6],fMx[3][3],fMy[3][3],fMz[3][3];
    AliITSAlignmentTrack *tr;

    for(i=0;i<3;i++) fx0[i] = p[i];
    SetByAngles(&(p[3]));
    dfMdthx(fMx);
    dfMdthy(fMy);
    dfMdthz(fMz);
    for(i=0;i<3;i++)for(j=0;j<6;j++)  dxdp[i][j] = 0.0;
    dxdp[0][0] = 1.0; // dx/dx
    dxdp[1][1] = 1.0; // dy/dy
    dxdp[2][2] = 1.0; // dz/dz

    for(i=0;i<ftrksM->GetEntriesFast();i++) {
	tr = (AliITSAlignmentTrack *)(ftrksM->At(i));
	for(j=0;j<tr->GetNumberOfClustersSl();j++){
	    tr->GetIndex(j,indx);
	    if(indx==findex){
		n++;
		tr->GetPointG(j,(Double_t *)xi);
		tr->func(xi,xo);
		tr->GetPointL(j,(Double_t *)xi);
		tr->GetErrorG(j,(Double_t **)Ex);
		LtoG(xi,xg);
		for(m=0;m<3;m++) for(k=0;k<3;k++){
		    dxdp[m][3] += fMx[m][k]*xi[k];
		    dxdp[m][4] += fMy[m][k]*xi[k];
		    dxdp[m][5] += fMz[m][k]*xi[k];
		} // end for m
		for(m=0;m<6;m++){
		    for(k=0;k<3;k++)for(l=0;l<3;l++)
			chi[m] += (xg[k] - xo[k])*Ex[k][l]*dxdp[l][m];
		} // end for m
	    } // end if indx==findex
	} // end for j
    } // end for i
    if(n<7) return ;
    for(m=0;m<6;m++) chi[m] /= (Double_t)(n-6);
    return;
}
//______________________________________________________________________
void AliITSAlignmentModule::AlignModule(){
    static Int_t npar=6;
    Int_t    iter,i;
    //Double_t p[npar],fret,gtol=1.0E-5;
    Double_t p[6],fret,gtol=1.0E-5;

    for(i=0;i<3;i++) {p[i] = fx0[i]; p[i+3] = fangles[i];}
    MRVMminimization(npar,(Double_t *)p,fret,gtol,iter);
    for(i=0;i<3;i++) {fx0[i] = p[i]; fangles[i] = p[i+3];}
    printf("AlignModule #%d:Xt=(%e,%e,%e) cm angles=(%e,%e,%e)rad,"
	   " Chi2=%e loops=%d\n",findex,
	   fx0[0],fx0[1],fx0[2],fangles[0],fangles[1],fangles[2],fret,iter);
}
//______________________________________________________________________
void AliITSAlignmentModule::Streamer(TBuffer &R__b){
   // Stream an object of class AliITSAlignmentModule.

   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(); if (R__v) { }
      TObject::Streamer(R__b);
      R__b >> findex;
      R__b >> flay;
      R__b >> flad;
      R__b >> fdet;
      R__b >> ftrksM;
      R__b >> fChi2;
//      R__b.ReadStaticArray(fx0);
//      R__b.ReadStaticArray((double*)fM);
//      R__b.ReadStaticArray(fangles);
   } else {
      R__b.WriteVersion(AliITSAlignmentModule::IsA());
      TObject::Streamer(R__b);
      R__b << findex;
      R__b << flay;
      R__b << flad;
      R__b << fdet;
      R__b << ftrksM;
      R__b << fChi2;
//      R__b.WriteArray(fx0, 3);
//      R__b.WriteArray((double*)fM, 9);
//      R__b.WriteArray(fangles, 3);
   }
}
