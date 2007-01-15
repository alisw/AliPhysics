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

//-------------------------------------------------------------------------
//   Field with Magnetic Field map
//   Used by AliRun class
//   Author:
//-------------------------------------------------------------------------

#include <TMath.h>
#include <TSystem.h>

#include "AliLog.h"
#include "AliMagFDM.h"

ClassImp(AliMagFDM)

//_______________________________________________________________________
AliMagFDM::AliMagFDM():
  fSolenoid(0),
  fInd(0),
  fZmin(0),
  fZmax(0),
  fYmax(0),
  fYmin(0),
  fZpmx(0),
  fZpmn(0),
  fRmax(0),
  fRmin(0),
  fXdel(0),
  fYdel(0),
  fZdel(0),
  fRdel(0),
  fPhid(0),
  fZpdl(0),
  fCx1(0),
  fCx2(0),
  fAx1(0),
  fAx2(0),
  fXl(0),
  fYl(0),
  fZl(0),    
  fRn(0),
  fPhin(0),
  fZpl(0)  
{
  //
  // Default constructor for the Dipole field
  //
}

//_______________________________________________________________________
AliMagFDM::AliMagFDM(const char *name, const char *title, Int_t integ,
                     Float_t factor, Float_t fmax):
  AliMagFC(name,title,integ,factor,fmax),
  fSolenoid(0),
  fInd(0),
  fZmin(0),
  fZmax(0),
  fYmax(0),
  fYmin(0),
  fZpmx(0),
  fZpmn(0),
  fRmax(0),
  fRmin(0),
  fXdel(0),
  fYdel(0),
  fZdel(0),
  fRdel(0),
  fPhid(0),
  fZpdl(0),
  fCx1(0),
  fCx2(0),
  fAx1(0),
  fAx2(0),
  fXl(0),
  fYl(0),
  fZl(0),    
  fRn(0),
  fPhin(0),
  fZpl(0)  
{
  //
  // Standard constructor for the Dipole field
  //
  fType = kDipoMap;
  fMap  = 3;
  SetSolenoidField();
  
  AliDebug(1, Form(
       "Field Map for Muon Arm from IP till muon filter %s created: map= %d, integ= %d, factor= %f, file=%s",
       fName.Data(), fMap ,integ,factor,fTitle.Data()));
 
}

//_______________________________________________________________________
void AliMagFDM::Field(Float_t *xfi, Float_t *b) const
{
  //
  // Main routine to compute the field in a point
  //
  const Double_t keps=0.1E-06;
  const Double_t kPI2=2.*TMath::Pi();
  const Double_t kone=1;
  
  const    Int_t  kiip=33; 
  const    Int_t  kmiip=0;    
  const    Int_t  kliip=0;
  
  const    Int_t  kiic=0;
  const    Int_t  kmiic=0;    
  const    Int_t  kliic=0;       
  
  const Double_t    kfZbg=502.92;  // Start of Map using in z
  const Double_t    kfZL3=600;  // Beginning of L3 door in z

  Double_t   x[3];   
  Double_t   xL3[3]; 
  Double_t   bint[3]; 
  
  Double_t r0;
  Int_t iKvar,jb;

  Double_t zp1, zp2,xp1,xp2,yp1,yp2; 
  Double_t zz1, zz2,yy1,yy2,x2,x1; 
     
// --- start the map fiel from z = 502.92 cm ---
//
// This map has been calculated in a coordinate system in which the muon spectrometer sits at z > 0
// Transfor correspondingly.

  x[0] = - xfi[0];
  x[1] =   xfi[1];
  x[2] = - xfi[2];
  b[0]=b[1]=b[2]=0;
//
  //       printf("x[0]  %f,x[1] %f,x[2]  %f\n",x[0],x[1],x[2]); 

  Double_t rr=TMath::Sqrt(x[0]*x[0]+x[1]*x[1]);
           r0=rr/100;
  Double_t rPmax;
  rPmax=fRmax; 
  if ( (-700<x[2] && x[2]<=kfZbg && 
        (x[0]*x[0]+(x[1]+30)*(x[1]+30))< 560*560)
       || (kfZbg<x[2] && x[2]<=kfZL3 && (rr>rPmax*100 && rr< 560)) )
       {
        b[0]=b[1]=0;
        b[2]=fSolenoid;
       }

  xL3[0]=x[0]/100;
  xL3[1]=(x[1]+30)/100;
  xL3[2]=x[2]/100; 

 
  if (TMath::Abs(xL3[0])<=0.1E-06) xL3[0]=TMath::Sign(0.1E-05,xL3[0]);
  if (TMath::Abs(xL3[1])<=0.1E-06) xL3[1]=TMath::Sign(0.1E-05,xL3[1]);

  Double_t xminn=xL3[2]*fAx1+fCx1;
  Double_t xmaxx=xL3[2]*fAx2+fCx2;
  Double_t zCmin,zCmax,yCmin,yCmax;

  zCmin=fZmin;
  zCmax=fZmax;
  yCmin=fYmin;
  yCmax=fYmax;
if ((kfZbg/100<xL3[2] && xL3[2]<=zCmin && r0<=rPmax) || ((zCmin<xL3[2] && xL3[2] <= zCmax ) && (yCmin<xL3[1] && xL3[1]<= yCmax) &&  (xminn<xL3[0] && xL3[0] <= xmaxx)))
    {
     if(fMap==3) 
      { 
        if (kfZbg/100<xL3[2] && xL3[2]<=zCmin && r0<=rPmax)
       {
       //---------------------  Polar part ----------------------
       
       Double_t yyp,ph0;
       Int_t  kp0, lp0, mp0;
       Int_t kpi,lpi,mpi;
       Double_t alp1,alp2,alp3;
       Double_t zpz,rp,fip,cphi; 

       kpi=kiip; 
       lpi=kliip;
       mpi=kmiip;
   
       zpz=xL3[2];
       kp0=FZ(zpz, fZp ,fZpdl,kpi,fZpl) ;
       zp2=(zpz-fZp[kp0])/(fZp[kp0+1]-fZp[kp0]);
       zp1= kone-zp2;
       yyp=xL3[1]- 0.3;
       cphi=TMath::Abs(yyp/r0);
       Int_t kcphi=0;
       if (cphi > kone) {
        AliDebug(2,Form("xL3[0] %e, xL3[1] %e, xL3[2] %e, yyp %e, r0 %e, cphi %e",xL3[0],xL3[1],xL3[2],yyp,r0,cphi));
        cphi =kone;
        kcphi=777;
       } 
       ph0=TMath::ACos(cphi);
       if (xL3[0] < 0 && yyp > 0 ) {ph0=kPI2/2 - ph0;}  
       if (xL3[0] < 0 && yyp < 0 ) {ph0=kPI2/2 + ph0;} 
       if (xL3[0] > 0 && yyp < 0 ) {ph0=kPI2 - ph0;}  
       if (ph0 > kPI2) {       ph0=ph0 - kPI2;}
       if (kcphi==777) {
        AliDebug(2,Form("xL3[0] %e, xL3[1] %e, xL3[2] %e, yyp %e, r0 %e, ph0 %e",xL3[0],xL3[1],xL3[2],yyp,r0,ph0));
       }  
       fip=ph0; 
       mp0=FZ(fip,fPhi,fPhid ,mpi,fPhin);
       xp2=(fip-fPhi[mp0])/(fPhi[mp0+1]-fPhi[mp0]);
       xp1=kone-xp2;

       Double_t rDel;
       rDel=fRdel;

       if (r0<= fRdel)
        {

         if(r0< keps) 
         {

          bint[0]=(zp1*fBpx[kp0][0][0] + zp2*fBpx[kp0+1][0][0])*10;     
          bint[1]=(zp1*fBpy[kp0][0][0] + zp2*fBpy[kp0+1][0][0])*10;      
          bint[2]=(zp1*fBpz[kp0][0][0] + zp2*fBpz[kp0+1][0][0])*10; 

         } else { 
      
        alp2= fB[0][0][mp0]*yyp + fB[0][1][mp0]*xL3[0]; 
        alp3= fB[1][0][mp0]*yyp + fB[1][1][mp0]*xL3[0];      
        alp1= kone - alp2 - alp3;
     
        for (jb=0; jb<3 ; jb++) 
         {
          iKvar=jb;     
          bint[jb] = Ba(iKvar,zp1,zp2,alp1,alp2,alp3, kp0,mp0)*10 ;
         } 
        } // end   of keps <=r0 
       }
        else
       {
        rp=r0;

        lp0=FZ(rp,fR ,fRdel,lpi,fRn);
        yp2=(rp-fR[lp0])/(fR[lp0+1]-fR[lp0]);
	yp1=kone-yp2;

        for (jb=0; jb<3 ; jb++) 
        {
         iKvar=jb;
         bint[jb] = Bb(zp1,zp2,yp1,yp2,xp1,xp2,iKvar,kp0,lp0,mp0)*10 ;
        }
       }

       b[0]=-bint[0];
       b[1]=bint[1];
       b[2]=-bint[2];

   }  
   else 
   {
   //-------------- Cartensian part ------------------
          
   Double_t zzc,yyc; 
   Int_t k0, l0,m0;
   Double_t  xx1, xx2,dx, xxx ,xXl; 
   Int_t kci,mci,lci;

   kci=kiic;
   lci=kliic;
   mci=kmiic;

   xx1 = fAx1*xL3[2] + fCx1;
   xx2 = fAx2*xL3[2] + fCx2;

   zzc=xL3[2];
   k0=FZ(zzc, fZc ,fZdel, kci, fZl); 
   zz2=(zzc-fZc[k0])/(fZc[k0+1]-fZc[k0]);
   zz1=kone - zz2;;    

   yyc=xL3[1];
   l0=FZ(yyc, fY , fYdel,lci,fYl);  
   yy2=(yyc-fY[l0])/(fY[l0+1]-fY[l0]);;
   yy1=kone - yy2;    
   xXl = fXl-kone; 
   dx = (xx2-xx1)/xXl;
   xxx= xL3[0]-xx1;
   m0 = int(xxx/dx);   

   if(xL3[0]<(xx1+m0*dx) || xL3[0] >(xx1+(m0+1)*dx)) 
    {
      m0=m0+1;   
      AliDebug(2,Form(" m0 %d, m0+1 %d\n",m0,m0+1));  
    }

   x2=(xL3[0]-( xx1+m0*dx))/dx;
   x1=kone-x2;
   m0=m0-1;
   for (jb=3; jb<6; jb++) 
     {
       iKvar=jb;     
       bint[jb-3] = Bb(zz1,zz2,yy1,yy2,x1,x2,iKvar,k0, l0, m0)*10 ; 
     }    
 
   b[0]=-bint[0];
   b[1]=bint[1];
   b[2]=-bint[2]; 

   } 


  } else {
        AliError(Form("Unknown map of Dipole region %d",fMap));
 }
           
} else {
    ZDCField(xfi,b);

    }

  if(fFactor!=1) {
    b[0]*=fFactor;
    b[1]*=fFactor;
    b[2]*=fFactor;
  }
}

//_______________________________________________________________________
Int_t AliMagFDM::FZ(Double_t temp, const Float_t *Ar, 
                    Float_t delu, Int_t ik, Int_t nk) const
{
  //
  // Quest of a point position at x,y,z (Cartensian) and R,Phi,z (Polar) axises
  //
  Int_t l,ikj;
  Double_t ddu,ar;

  Int_t ku,kf;
  kf=0;
  ar=Ar[ik];
  ddu=temp-ar;

  ku=int(ddu/delu);
  ikj=ik+ku+1;
  if (ddu<=0) ikj=0;

   for(l=ikj; l<nk; l++)
   {

    if(temp <= Ar[l])
     {
       
      kf=l-1;
      break;
     }
    }
    return kf;
  }

//_______________________________________________________________________
Double_t AliMagFDM::Ba(Int_t kaai,Double_t zaa1, Double_t zaa2, 
                       Double_t alf1, Double_t alf2, Double_t alf3, 
                       Int_t kaa, Int_t maa) const
{
  //
  // Calculation of field componet for case (keps <r0<= fRdel) at a given axis
  //
  Double_t fa11=0,fa12=0,fa13=0;
  Double_t fa21=0,fa22=0,fa23=0;
  Double_t faY1=0,faY2=0;
  Double_t bba=0;

  switch (kaai) {
  case 0:
    fa11 = fBpx[kaa][0][0];
    fa12 = fBpx[kaa][0][maa];
    fa13 = fBpx[kaa][0][maa+1];
    fa21 = fBpx[kaa+1][0][0];
    fa22 = fBpx[kaa+1][0][maa];               
    fa23 = fBpx[kaa+1][0][maa+1]; 
    break;
  case 1:
    fa11 = fBpy[kaa][0][0];
    fa12 = fBpy[kaa][0][maa];
    fa13 = fBpy[kaa][0][maa+1];
    fa21 = fBpy[kaa+1][0][0];
    fa22 = fBpy[kaa+1][0][maa];               
    fa23 = fBpy[kaa+1][0][maa+1]; 
    break;
  case 2:
    fa11 = fBpz[kaa][0][0];
    fa12 = fBpz[kaa][0][maa];
    fa13 = fBpz[kaa][0][maa+1];
    fa21 = fBpz[kaa+1][0][0];
    fa22 = fBpz[kaa+1][0][maa];               
    fa23 = fBpz[kaa+1][0][maa+1]; 
    break;
  default:
    AliFatal(Form("Invalid value of kaai %d",kaai));
  }                            
  faY1=alf1*fa11+alf2*fa12+alf3*fa13;
  faY2=alf1*fa21+alf2*fa22+alf3*fa23;
  bba =  zaa1*faY1+zaa2*faY2;    
  return bba;  
}

  
//_______________________________________________________________________
Double_t AliMagFDM::Bb(Double_t z1,Double_t z2, Double_t y1,Double_t y2, 
                       Double_t x1,Double_t x2, Int_t kv, Int_t k, Int_t l, Int_t m) const
{  
  //
  // Calculation of field componet at a given axis (general case)
  //
  Double_t fy1=0,fy2=0,ffy=0;
  Double_t gy1=0,gy2=0,ggy=0;
  Double_t bbi=0;
  Double_t bf11=0,bf12=0,bf21=0,bf22=0;
  Double_t bg11=0,bg12=0,bg21=0,bg22=0;

  
  /*-----------------Polar part ------------------*/
  
  switch (kv) {
  case 0:
    bf11=fBpx[k][l][m];
    bf12=fBpx[k+1][l][m];
    bf21=fBpx[k+1][l+1][m];
    bf22=fBpx[k][l+1][m];
    
    bg11=fBpx[k][l][m+1];
    bg12=fBpx[k+1][l][m+1];
    bg21=fBpx[k+1][l+1][m+1];
    bg22=fBpx[k][l+1][m+1];
    break;

  case 1:
    bf11=fBpy[k][l][m];
    bf12=fBpy[k+1][l][m];
    bf21=fBpy[k+1][l+1][m];
    bf22=fBpy[k][l+1][m];
    
    bg11=fBpy[k][l][m+1];
    bg12=fBpy[k+1][l][m+1];
    bg21=fBpy[k+1][l+1][m+1];
    bg22=fBpy[k][l+1][m+1];
    break;

  case 2:
    bf11=fBpz[k][l][m];
    bf12=fBpz[k+1][l][m];
    bf21=fBpz[k+1][l+1][m];
    bf22=fBpz[k][l+1][m];
    
    bg11=fBpz[k][l][m+1];
    bg12=fBpz[k+1][l][m+1];
    bg21=fBpz[k+1][l+1][m+1];
    bg22=fBpz[k][l+1][m+1]; 
    break;
  /*-----------------Cartensian part ---------------*/ 
  
  case 3:
    bf11=fBcx[k][l][m];
    bf12=fBcx[k+1][l][m];
    bf21=fBcx[k+1][l+1][m];
    bf22=fBcx[k][l+1][m];
    
    bg11=fBcx[k][l][m+1];
    bg12=fBcx[k+1][l][m+1];
    bg21=fBcx[k+1][l+1][m+1];
    bg22=fBcx[k][l+1][m+1];
    break;

  case 4:
    bf11=fBcy[k][l][m];
    bf12=fBcy[k+1][l][m];
    bf21=fBcy[k+1][l+1][m];
    bf22=fBcy[k][l+1][m];
    
    bg11=fBcy[k][l][m+1];
    bg12=fBcy[k+1][l][m+1];
    bg21=fBcy[k+1][l+1][m+1];
    bg22=fBcy[k][l+1][m+1];
    break;

  case 5:
    bf11=fBcz[k][l][m];
    bf12=fBcz[k+1][l][m];
    bf21=fBcz[k+1][l+1][m];
    bf22=fBcz[k][l+1][m];
    
    bg11=fBcz[k][l][m+1];
    bg12=fBcz[k+1][l][m+1];
    bg21=fBcz[k+1][l+1][m+1];
    bg22=fBcz[k][l+1][m+1];
    break;

  default:
    AliFatal(Form("Invalid value of kv %d",kv));
  }  
  
  
  
  fy1=z1*bf11+z2*bf12;
  fy2=z2*bf21+z1* bf22;
  ffy=y1*fy1+ y2*fy2; 
  
  
  gy1 = z1*bg11+z2*bg12;
  gy2 = z2*bg21+z1*bg22;
  ggy= y1*gy1 +  y2*gy2;
  
  bbi =  x1*ffy+x2*ggy;

  return bbi;
  
}     

//_______________________________________________________________________
void AliMagFDM::ReadField()
{
  //
  // Method to read the magnetic field from file
  //
  FILE *magfile;

  Int_t ik, il, im;
  Float_t zzp, rr,phii;
  Float_t zz, yy, bx,by,bz,bb;

  char *fname;
  AliDebug(1,Form("Reading Magnetic Field %s from file %s",fName.Data(),fTitle.Data()));
  fname = gSystem->ExpandPathName(fTitle.Data());
  magfile=fopen(fname,"r");
  delete [] fname;

  AliDebug(2,"Cartensian part");
 
  if (magfile) {
  
//  Cartensian part 
 
    fscanf(magfile,"%d %d %d ",&fYl, &fXl, &fZl); 
    
    AliDebug(3,Form("fYl %d, fXl %d, fZl %d",fYl, fXl, fZl));     
    
    for (ik=0; ik<fZl; ik++)
    { 
   
      fscanf(magfile, " %e ", &zz);
      fZc[ik]=zz; 

    } 
   
    for (ik=0; ik<fYl; ik++)
    {    
       fscanf(magfile, " %e ", &yy); 
       fY[ik]=yy;
 
    } 
    for (ik=0; ik<81; ik++)
    {    
           AliDebug(4,Form("fZc %e,fY %e", fZc[ik],fY[ik])); 
    }   
             
    fscanf(magfile," %e %e %e %e %e %e %e %e %e %e %e ", &fYdel,&fXdel,&fZdel,&fZmax,&fZmin,&fYmax,&fYmin,&fAx1,&fCx1,&fAx2,&fCx2); 

AliDebug(3,Form("fYdel %e, fXdel %e, fZdel %e",fYdel,fXdel,fZdel));
AliDebug(3,Form("fZmax %e, fZmin %e, fYmax %e,fYmin %e",fZmax,fZmin,fYmax,fYmin));
AliDebug(3,Form("fAx1 %e, fCx1 %e, fAx2 %e, fCx %e",fAx1,fCx1,fAx2,fCx2));

    for (il=0; il<44; il++)  { 
     for (im=0; im<81; im++)  {      
      for (ik=0; ik<81; ik++)  {      
      
      fscanf(magfile, " %e ", &by); 
      fBcy[ik][im][il]=by;        
      }
     }
    } 

    for (il=0; il<44; il++)  { 
     for (im=0; im<81; im++)  {      
      for (ik=0; ik<81; ik++)  {      
      
      fscanf(magfile, " %e ", &bx); 
      fBcx[ik][im][il]=bx;        
      }    
     }     
    }
    
   for (il=0; il<44; il++)  { 
     for (im=0; im<81; im++)  {      
      for (ik=0; ik<81; ik++)  {      
      
      fscanf(magfile, " %e ", &bz); 
      fBcz[ik][im][il]=bz;          
      }              
     }     
    } 
//----------------------   Polar part ---------------------------------

    AliDebug(2,"Polar part");
    fscanf(magfile,"%d %d %d ", &fZpl, &fRn, &fPhin); 
    AliDebug(3,Form("fZpl %d, fRn %d, fPhin %d",fZpl,fRn,fPhin));   

    AliDebug(4," fZp array"); 
     
    for (ik=0; ik<51; ik++) 
    {    
     fscanf(magfile, " %e ", &zzp);
     fZp[ik]=zzp; 
     AliDebug(4,Form(" %e",fZp[ik]));      
    } 
  
    AliDebug(4," fR array"); 
         
    for (ik=0; ik<10; ik++) 
    {    
     fscanf(magfile, " %e ", &rr); 
     fR[ik]=rr;
     AliDebug(4,Form(" %e",fR[ik]));
    } 
    
//    AliDebug(4,"fPhi array"); 
     
     for (il=0; il<33; il++)  
     {
       fscanf(magfile, " %e ", &phii); 
       fPhi[il]=phii; 
//        AliDebug(4,Form(" %e",fPhi[il]));          
     }

    fscanf(magfile," %e %e %e %e %e %e %e ",&fZpdl,&fPhid,&fRdel,&fZpmx,&fZpmn,&fRmax, &fRmin); 

AliDebug(3,Form("fZpdl %e, fPhid %e, fRdel %e, fZpmx %e, fZpmn %e,fRmax %e,fRmin %e", fZpdl,fPhid, fRdel,fZpmx, fZpmn,fRmax, fRmin));

                      
    for (il=0; il<33; il++)  { 
     for (im=0; im<10; im++)  {      
      for (ik=0; ik<51; ik++)  {
      fscanf(magfile, " %e ", &by); 
        fBpy[ik][im][il]=by;        
      }
     }
    } 
    
    for (il=0; il<33; il++)  { 
     for (im=0; im<10; im++)  {      
      for (ik=0; ik<51; ik++)  {
      fscanf(magfile, " %e ", &bx); 
        fBpx[ik][im][il]=bx;                     
      }    
     }     
    }      


    for (il=0; il<33; il++)  { 
     for (im=0; im<10; im++)  {      
      for (ik=0; ik<51; ik++)  {
      fscanf(magfile, " %e ", &bz); 
        fBpz[ik][im][il]=bz;                      
      }              
     }     
    } 
    
    
    for (il=0; il<32; il++) { 
     for (im=0; im<2; im++)  {      
      for (ik=0; ik<2; ik++)  {
      fscanf(magfile, " %e ", &bb);    
        fB[ik][im][il]=bb;  
      }              
     } 
    }
//
  } else { 
    AliFatal(Form("File %s not found !",fTitle.Data()));
  }
}
