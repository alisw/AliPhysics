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
Revision 1.6  2000/11/10 18:09:55  fca
New field map for the ZDC

Revision 1.5  2000/10/27 14:17:04  morsch
- Bug causing segmentation violation during muon reconstruction corrected
- Coding rule violations corrected.
(Galina Chabratova)

Revision 1.4  2000/10/02 21:28:14  fca
Removal of useless dependecies via forward declarations

Revision 1.3  2000/07/13 16:19:09  fca
Mainly coding conventions + some small bug fixes

Revision 1.2  2000/07/12 08:56:25  fca
Coding convention correction and warning removal

Revision 1.1  2000/07/11 18:24:59  fca
Coding convention corrections + few minor bug fixes

*/

#include <stdlib.h>

#include "AliMagFDM.h"
#include "TSystem.h"
    

ClassImp(AliMagFDM)

//________________________________________
AliMagFDM::AliMagFDM(const char *name, const char *title, const Int_t integ,
const Int_t map, const Float_t factor, const Float_t fmax)
  : AliMagF(name,title,integ,map,factor,fmax)
  
{
  //
  // Standard constructor for the Dipole field
  //
  fType = kDipoMap;

  printf("Field Map for Muon Arm from IP till muon filter %s created: map= %d, factor= %f, file=%s\n",fName.Data(),map,factor,fTitle.Data());
  
}

//________________________________________

void AliMagFDM::Field(Float_t *xfi, Float_t *b)
{
  //
  // Main routine to compute the field in a point
  //
  static  const Double_t keps=0.1E-06;
  static  const Double_t kpi2=.6283185E+01;
  static  const Double_t kone=1;

  static  const    Int_t  kiip=33; 
  static  const    Int_t  kmiip=0;    
  static  const    Int_t  kliip=0;

  static  const    Int_t  kiic=0;
  static  const    Int_t  kmiic=0;    
  static  const    Int_t  kliic=0;       

  static  const Double_t    kfZbg=502.92;  // Start of Map using in z
  static  const Double_t    kfZL3=600;  // Beginning of L3 door in z

  Double_t   x[3];   
  Double_t   xL3[3]; 
  Double_t   bint[3]; 
  
  Double_t r0;
  
  Double_t bbj;
  Int_t iKvar,jb;

  Double_t zp1, zp2,xp1,xp2,yp1,yp2; 
  Double_t zz1, zz2,yy1,yy2,x2,x1; 
     
// --- start the map fiel from z = 502.92 cm ---

  x[0] = xfi[0];
  x[1] = xfi[1];
  x[2] = xfi[2];
    
  /*  
  if (x[0]==0) x[0]=0.1E-06;
  if (x[1]==0) x[1]=0.1E-06;
  */ 
 b[0]=b[1]=b[2]=0;
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
        b[2]=2;
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
  /*          
if ((kfZbg/100<xL3[2] && xL3[2]<zCmin && r0<rPmax) || ((zCmin<=xL3[2] && xL3[2] <= zCmax ) && (yCmin<=xL3[1] && xL3[1]<= yCmax) &&  (xminn <= xL3[0] && xL3[0] <= xmaxx)))
  */
if ((kfZbg/100<xL3[2] && xL3[2]<=zCmin && r0<=rPmax) || ((zCmin<xL3[2] && xL3[2] <= zCmax ) && (yCmin<xL3[1] && xL3[1]<= yCmax) &&  (xminn<xL3[0] && xL3[0] <= xmaxx)))
    {
     if(fMap==3) 
      { 
        if (kfZbg/100<xL3[2] && xL3[2]<=zCmin && r0<=rPmax)
	//       if (xL3[2]<zCmin && r0<rPmax)   
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

       FZ(&zpz, fZp ,&fZpdl,&kpi,&kp0,&zp1 ,&zp2,&fZpl) ;

       yyp=xL3[1]- 0.3;
       cphi=yyp/r0;
       ph0=TMath::ACos(cphi);
       if (xL3[0]< 0) {ph0=kpi2 - ph0;}  
                                  
       fip=ph0;
       FZ(&fip,fPhi,&fPhid ,&mpi,&mp0, &xp1,&xp2,&fPhin); 

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
          FRfuncBi(&iKvar,&zp1,&zp2,&alp1,&alp2,&alp3, &kp0,&mp0, &bbj);
          bint[jb] = bbj*10 ;
         } 
        } // end   of keps <=r0 
       }
        else
       {
        rp=r0;

        FZ(&rp,fR ,&fRdel,&lpi,&lp0,&yp1,&yp2,&fRn);

        for (jb=0; jb<3 ; jb++) 
        {
         iKvar=jb;
         FGfuncBi(&zp1,&zp2,&yp1,&yp2,&xp1,&xp2,&iKvar,&kp0,&lp0,&mp0,&bbj); 

         bint[jb] = bbj*10 ;
        }
       }

       b[0]=bint[0];
       b[1]=bint[1];
       b[2]=bint[2];

//    fprintf(fitest,"-------------   Freg2 run -------------\n");       
  
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
   FZ(&zzc, fZc ,&fZdel, &kci,&k0, &zz1, &zz2, &fZl);     

   yyc=xL3[1];
   FZ(&yyc, fY , &fYdel,&lci, &l0, &yy1, &yy2,&fYl);  

   xXl = fXl-kone; 
   dx = (xx2-xx1)/xXl;
   xxx= xL3[0]-xx1;
   //     xm = xxx/dx; 
   m0 = int(xxx/dx);   

   if(xL3[0]<(xx1+m0*dx) || xL3[0] >(xx1+(m0+1)*dx)) 
    {
      m0=m0+1;   
      printf(" m0 %d, m0+1 %d\n",m0,m0+1);  
    }

   x2=(xL3[0]-( xx1+m0*dx))/dx;
   x1=kone-x2;
   m0=m0-1;
   for (jb=3; jb<6; jb++) 
     {
       iKvar=jb;     
       FGfuncBi(&zz1,&zz2,&yy1,&yy2,&x1,&x2,&iKvar,&k0, &l0, &m0, &bbj); 
       bint[jb-3] = bbj*10 ; 
     }    
 
   b[0]=bint[0];
   b[1]=bint[1];
   b[2]=bint[2]; 

//   fprintf(fitest,"------------   Freg1 run -----------------\n");            
   } 

  } else {
        printf("Unknown map of Dipole region %d\n",fMap);
 }
           
} else {

//This is the ZDC part
    Float_t rad2=x[0]*x[0]+x[1]*x[1];
    if(x[2]>kCORBEG2 && x[2]<kCOREND2){
      if(rad2<kCOR2RA2){
        b[0] = kFCORN2;
      }
    }
    else if(x[2]>kZ1BEG && x[2]<kZ1END){  
      if(rad2<kZ1RA2){
        b[0] = -kG1*x[1];
        b[1] = -kG1*x[0];
      }
    }
    else if(x[2]>kZ2BEG && x[2]<kZ2END){  
      if(rad2<kZ2RA2){
        b[0] = kG1*x[1];
        b[1] = kG1*x[0];
      }
    }
    else if(x[2]>kZ3BEG && x[2]<kZ3END){  
      if(rad2<kZ3RA2){
        b[0] = kG1*x[1];
        b[1] = kG1*x[0];
      }
    }
    else if(x[2]>kZ4BEG && x[2]<kZ4END){  
      if(rad2<kZ4RA2){
        b[0] = -kG1*x[1];
        b[1] = -kG1*x[0];
      }
    }
    else if(x[2]>kD1BEG && x[2]<kD1END){ 
      if(rad2<kD1RA2){
        b[1] = -kFDIP;
      }
    }
    else if(x[2]>kD2BEG && x[2]<kD2END){
      if(((x[0]-kXCEN1D2)*(x[0]-kXCEN1D2)+(x[1]-kYCEN1D2)*(x[1]-kYCEN1D2))<kD2RA2
        || ((x[0]-kXCEN2D2)*(x[0]-kXCEN2D2)+(x[1]-kYCEN2D2)*(x[1]-kYCEN2D2))<kD2RA2){
	b[1] = kFDIP;
      }
    }
      }
    

  if(fFactor!=1) {
    b[0]*=fFactor;
    b[1]*=fFactor;
    b[2]*=fFactor;
  }
}

//_________________________________________

void AliMagFDM::FZ(Double_t *u, Float_t  *Ar, Float_t *du,Int_t *ki,Int_t *kf,Double_t *a1,Double_t *a2 ,Int_t *nu)
{
  //
  // Z component of the field
  //
  static  const Double_t kone=1;
  Int_t l,ik,ikj;
  Double_t temp;
  Double_t ddu,delu,ar;

  Int_t nk,ku;
  temp=*u;
  nk=*nu;
  ik=*ki;
  delu=*du; 

  ar=Ar[ik];
  ddu=temp-ar;

  ku=int(ddu/delu);
  //  ikj=ik+ku; //falsh
    ikj=ik+ku+1; //correct
  if (ddu<=0) ikj=0;

   for(l=ikj; l<nk; l++)
   {
     //    if(temp < Ar[l])
    if(temp <= Ar[l])
     {
       
      *kf=l-1;
      *a2=(temp-Ar[l-1])/(Ar[l]-Ar[l-1]);
      
       /*
      *kf=l;
      *a2=(temp-Ar[l])/(Ar[l+1]-Ar[l]);
      */
       break;
     }
    }
   //    *a2=(temp-Ar[*kf])/(Ar[*kf+1]-Ar[*kf]);
    *a1= kone - *a2;
  }

/*-------------FRfuncBi----------------*/

void AliMagFDM::FRfuncBi(Int_t *kai,Double_t *za1, Double_t *za2, Double_t *al1, Double_t *al2, Double_t *al3, Int_t *ka, Int_t *ma, Double_t *ba)
{
  //
  // This method needs to be commented
  //
  Double_t fa11,fa12,fa13;
  Double_t fa21,fa22,fa23;
  Double_t faY1,faY2;
  Double_t bba;
  
  Double_t zaa1,zaa2,alf1,alf2,alf3;
  Int_t kaai,kaa,maa;
  kaai=*kai;
  kaa=*ka;
  maa=*ma;
  zaa1=*za1;
  zaa2=*za2;
  alf1=*al1;  
  alf2=*al2;
  alf3=*al3; 
  
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
    Fatal("FRfuncBi","Invalid value of kaai %d\n",kaai);
    exit(1);
  }                            
  faY1=alf1*fa11+alf2*fa12+alf3*fa13;
  faY2=alf1*fa21+alf2*fa22+alf3*fa23;
  bba =  zaa1*faY1+zaa2*faY2;    
  *ba=bba;
  
}

  
/*----------- FGfuncBi------------*/

void AliMagFDM::FGfuncBi(Double_t *zz1,Double_t *zz2, Double_t *yy1,Double_t *yy2, Double_t *xx1,Double_t *xx2, Int_t *kvr, Int_t *kk, Int_t *ll, Int_t *mm, Double_t *bb)
{  
  //
  // This method needs to be commented
  //
  Double_t fy1, fy2, ffy;
  Double_t gy1,gy2,ggy;
  Double_t z1,z2,y1,y2,x1,x2;
  
  Int_t k,l,m,kv;
  Double_t bbi;
  
  Double_t bf11,bf12,bf21,bf22;
  Double_t bg11,bg12,bg21,bg22;
  k=*kk;
  l=*ll;
  m=*mm;
  
  kv=*kvr;
  
  z1=*zz1;
  z2=*zz2;
  y1=*yy1;
  y2=*yy2;
  x1=*xx1;
  x2=*xx2; 
  
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
    Fatal("FGfuncBi","Invalid value of kv %d\n",kv);
    exit(1);
  }  
  
  
  
  fy1=z1*bf11+z2*bf12;
  fy2=z2*bf21+z1* bf22;
  ffy=y1*fy1+ y2*fy2; 
  
  
  gy1 = z1*bg11+z2*bg12;
  gy2 = z2*bg21+z1*bg22;
  ggy= y1*gy1 +  y2*gy2;
  
  bbi =  x1*ffy+x2*ggy;
  
  *bb=bbi; 
  
}     
//____________________________________________

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
  printf("Reading Magnetic Field %s from file %s\n",fName.Data(),fTitle.Data());
  fname = gSystem->ExpandPathName(fTitle.Data());
  magfile=fopen(fname,"r");
  delete [] fname;

  printf("Cartensian part\n");
 
  if (magfile) {
  
//  Cartensian part 
 
    fscanf(magfile,"%d %d %d ",&fYl, &fXl, &fZl); 
    
    printf("fYl %d, fXl %d, fZl %d\n",fYl, fXl, fZl);     
    
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
           printf("fZc %e,fY %e\n", fZc[ik],fY[ik]); 
    }   
             
    fscanf(magfile," %e %e %e %e %e %e %e %e %e %e %e ", &fYdel,&fXdel,&fZdel,&fZmax,&fZmin,&fYmax,&fYmin,&fAx1,&fCx1,&fAx2,&fCx2); 

printf("fYdel %e, fXdel %e, fZdel %e\n",fYdel,fXdel,fZdel);
printf("fZmax %e, fZmin %e, fYmax %e,fYmin %e\n",fZmax,fZmin,fYmax,fYmin);
printf("fAx1 %e, fCx1 %e, fAx2 %e, fCx %e\n",fAx1,fCx1,fAx2,fCx2);

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

    printf("Polar part\n");
    fscanf(magfile,"%d %d %d ", &fZpl, &fRn, &fPhin); 
    printf("fZpl %d, fRn %d, fPhin %d\n",fZpl,fRn,fPhin);   

    printf(" fZp array\n"); 
     
    for (ik=0; ik<51; ik++) 
    {    
     fscanf(magfile, " %e ", &zzp);
     fZp[ik]=zzp; 
     printf(" %e\n",fZp[ik]);      
    } 
  
    printf(" fR array\n"); 
         
    for (ik=0; ik<10; ik++) 
    {    
     fscanf(magfile, " %e ", &rr); 
     fR[ik]=rr;
     printf(" %e\n",fR[ik]);
    } 
    
//    printf("fPhi array\n"); 
     
     for (il=0; il<33; il++)  
     {
       fscanf(magfile, " %e ", &phii); 
       fPhi[il]=phii; 
//        printf(" %e\n",fPhi[il]);          
     }

    fscanf(magfile," %e %e %e %e %e %e %e ",&fZpdl,&fPhid,&fRdel,&fZpmx,&fZpmn,&fRmax, &fRmin); 

printf("fZpdl %e, fPhid %e, fRdel %e, fZpmx %e, fZpmn %e,fRmax %e,fRmin %e \n", fZpdl,fPhid, fRdel,fZpmx, fZpmn,fRmax, fRmin);

                      
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
    printf("File %s not found !\n",fTitle.Data());
    exit(1);
  }
}
//________________________________  
