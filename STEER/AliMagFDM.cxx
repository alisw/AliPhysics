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
*/

#include "AliMagFDM.h"
#include "TSystem.h"
    

ClassImp(AliMagFDM)

//________________________________________
AliMagFDM::AliMagFDM(const char *name, const char *title, const Int_t integ,
const Int_t map, const Float_t factor, const Float_t fmax)
  : AliMagF(name,title,integ,map,factor,fmax)
  
{
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

  static  const Double_t    kfdZbg=502.92;  // Start of Map using in z
  static  const Double_t    kfdZL3=600;  // Beginning of L3 door in z

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
  b[0]=b[1]=b[2]=0;
  //       printf("x[0]  %f,x[1] %f,x[2]  %f\n",x[0],x[1],x[2]); 

  Double_t rr=TMath::Sqrt(x[0]*x[0]+x[1]*x[1]);
           r0=rr/100;
  Double_t rPmax;
  rPmax=fdRmax; 
  if ( (-700<x[2] && x[2]<=kfdZbg && 
        (x[0]*x[0]+(x[1]+30)*(x[1]+30))< 560*560)
       || (kfdZbg<x[2] && x[2]<=kfdZL3 && rr>=rPmax*100) )
       {
        b[2]=2;
       }

  xL3[0]=x[0]/100;
  xL3[1]=(x[1]+30)/100;
  xL3[2]=x[2]/100; 
 
  Double_t xminn=xL3[2]*fdAx1+fdCx1;
  Double_t xmaxx=xL3[2]*fdAx2+fdCx2;
  Double_t zCmin,zCmax,yCmin,yCmax;

  zCmin=fdZmin;
  zCmax=fdZmax;
  yCmin=fdYmin;
  yCmax=fdYmax;
          
if ((kfdZbg/100<xL3[2] && xL3[2]<zCmin && r0<rPmax) || ((zCmin<=xL3[2] && xL3[2] <= zCmax ) && (yCmin<=xL3[1] && xL3[1]<= yCmax) &&  (xminn <= xL3[0] && xL3[0] <= xmaxx)))
    {
     if(fMap==3) 
      { 
       if (xL3[2]<zCmin && r0<rPmax)   
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

       FZ(&zpz, fdZp ,&fdZpdl,&kpi,&kp0,&zp1 ,&zp2,&fdZpl) ;

       yyp=xL3[1]- 0.3;
       cphi=yyp/r0;
       ph0=TMath::ACos(cphi);
       if (xL3[0]< 0) {ph0=kpi2 - ph0;}  
                                  
       fip=ph0;
       FZ(&fip,fdPhi,&fdPhid ,&mpi,&mp0, &xp1,&xp2,&fdPhin); 

       Double_t rDel;
       rDel=fdRdel;

       if (r0<= fdRdel)
        {

         if(r0< keps) 
         {

          bint[0]=(zp1*fdBpx[kp0][0][0] + zp2*fdBpx[kp0+1][0][0])*10;     
          bint[1]=(zp1*fdBpy[kp0][0][0] + zp2*fdBpy[kp0+1][0][0])*10;      
          bint[2]=(zp1*fdBpz[kp0][0][0] + zp2*fdBpz[kp0+1][0][0])*10; 

         }  
      
        alp2= fdB[0][0][mp0]*yyp + fdB[0][1][mp0]*xL3[0]; 
        alp3= fdB[1][0][mp0]*yyp + fdB[1][1][mp0]*xL3[0];      
        alp1= kone - alp2 - alp3;
     
        for (jb=0; jb<3 ; jb++) 
        {
         iKvar=jb;     
         FRfuncBi(&iKvar,&zp1,&zp2,&alp1,&alp2,&alp3, &kp0,&mp0, &bbj);
         bint[jb] = bbj*10 ;
        }  
       }
        else
       {
        rp=r0;

        FZ(&rp,fdR ,&fdRdel,&lpi,&lp0,&yp1,&yp2,&fdRn);

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

   xx1 = fdAx1*xL3[2] + fdCx1;
   xx2 = fdAx2*xL3[2] + fdCx2;

   zzc=xL3[2];
   FZ(&zzc, fdZc ,&fdZdel, &kci,&k0, &zz1, &zz2, &fdZl);     

   yyc=xL3[1];
   FZ(&yyc, fdY , &fdYdel,&lci, &l0, &yy1, &yy2,&fdYl);  

   xXl = fdXl-kone; 
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
      if(rad2<kD2RA2) {
        if(x[2]>kD2BEG) {
          
//    Separator Dipole D2
          if(x[2]<kD2END) b[1]=kFDIP;
        } else if(x[2]>kD1BEG) {

//    Separator Dipole D1
          if(x[2]<kD1END) b[1]=-kFDIP;
        }
        if(rad2<kCORRA2) {

//    First quadrupole of inner triplet de-focussing in x-direction
//    Inner triplet
          if(x[2]>kZ4BEG) {
            if(x[2]<kZ4END) {
              
//    2430 <-> 3060
              b[0]=-kG1*x[1];
              b[1]=-kG1*x[0];
            }
          } else if(x[2]>kZ3BEG) {
            if(x[2]<kZ3END) {

//    1530 <-> 2080
              b[0]=kG1*x[1];
              b[1]=kG1*x[0];
            }
          } else if(x[2]>kZ2BEG) {
            if(x[2]<kZ2END) {

//    890 <-> 1430
              b[0]=kG1*x[1];
              b[1]=kG1*x[0];
            }
          } else if(x[2]>kZ1BEG) {
            if(x[2]<kZ1END) {

//    0 <->  630
              b[0]=-kG1*x[1];
              b[1]=-kG1*x[0];
            }
          } else if(x[2]>kCORBEG) {
            if(x[2]<kCOREND) {
//    Corrector dipole (because of dimuon arm)
//            b[0]=kFCORN;
              b[0]=-kFCORN;
            }
          }
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
  ikj=ik+ku;
  if (ddu<=0) ikj=0;

   for(l=ikj; l<nk; l++)
   {

    if(temp < Ar[l])
     {
      *kf=l;
      *a2=(temp-Ar[l])/(Ar[l+1]-Ar[l]);
      *a1= kone - *a2;
      break;     
     }
    }
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
  
  if (kaai==0 ) {
    fa11 = fdBpx[kaa][0][0];
    fa12 = fdBpx[kaa][0][maa];
    fa13 = fdBpx[kaa][0][maa+1];
    fa21 = fdBpx[kaa+1][0][0];
    fa22 = fdBpx[kaa+1][0][maa];               
    fa23 = fdBpx[kaa+1][0][maa+1]; 
  }
  if (kaai==1 ) {
    fa11 = fdBpy[kaa][0][0];
    fa12 = fdBpy[kaa][0][maa];
    fa13 = fdBpy[kaa][0][maa+1];
    fa21 = fdBpy[kaa+1][0][0];
    fa22 = fdBpy[kaa+1][0][maa];               
    fa23 = fdBpy[kaa+1][0][maa+1]; 
  }
  if (kaai==2 ) {
    fa11 = fdBpz[kaa][0][0];
    fa12 = fdBpz[kaa][0][maa];
    fa13 = fdBpz[kaa][0][maa+1];
    fa21 = fdBpz[kaa+1][0][0];
    fa22 = fdBpz[kaa+1][0][maa];               
    fa23 = fdBpz[kaa+1][0][maa+1]; 
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
  
  if(kv==0) {
    bf11=fdBpx[k][l][m];
    bf12=fdBpx[k+1][l][m];
    bf21=fdBpx[k+1][l+1][m];
    bf22=fdBpx[k][l+1][m];
    
    bg11=fdBpx[k][l][m+1];
    bg12=fdBpx[k+1][l][m+1];
    bg21=fdBpx[k+1][l+1][m+1];
    bg22=fdBpx[k][l+1][m+1];
  }
  if(kv==1) {
    bf11=fdBpy[k][l][m];
    bf12=fdBpy[k+1][l][m];
    bf21=fdBpy[k+1][l+1][m];
    bf22=fdBpy[k][l+1][m];
    
    bg11=fdBpy[k][l][m+1];
    bg12=fdBpy[k+1][l][m+1];
    bg21=fdBpy[k+1][l+1][m+1];
    bg22=fdBpy[k][l+1][m+1];
  }  
  
  if(kv==2) {
    bf11=fdBpz[k][l][m];
    bf12=fdBpz[k+1][l][m];
    bf21=fdBpz[k+1][l+1][m];
    bf22=fdBpz[k][l+1][m];
    
    bg11=fdBpz[k][l][m+1];
    bg12=fdBpz[k+1][l][m+1];
    bg21=fdBpz[k+1][l+1][m+1];
    bg22=fdBpz[k][l+1][m+1];           
  } 
  /*-----------------Cartensian part ---------------*/ 
  
  if(kv==3) {
    bf11=fdBcx[k][l][m];
    bf12=fdBcx[k+1][l][m];
    bf21=fdBcx[k+1][l+1][m];
    bf22=fdBcx[k][l+1][m];
    
    bg11=fdBcx[k][l][m+1];
    bg12=fdBcx[k+1][l][m+1];
    bg21=fdBcx[k+1][l+1][m+1];
    bg22=fdBcx[k][l+1][m+1];
  }          
  
  if(kv==4) {
    bf11=fdBcy[k][l][m];
    bf12=fdBcy[k+1][l][m];
    bf21=fdBcy[k+1][l+1][m];
    bf22=fdBcy[k][l+1][m];
    
    bg11=fdBcy[k][l][m+1];
    bg12=fdBcy[k+1][l][m+1];
    bg21=fdBcy[k+1][l+1][m+1];
    bg22=fdBcy[k][l+1][m+1];
  }          
  if(kv==5) {
    bf11=fdBcz[k][l][m];
    bf12=fdBcz[k+1][l][m];
    bf21=fdBcz[k+1][l+1][m];
    bf22=fdBcz[k][l+1][m];
    
    bg11=fdBcz[k][l][m+1];
    bg12=fdBcz[k+1][l][m+1];
    bg21=fdBcz[k+1][l+1][m+1];
    bg22=fdBcz[k][l+1][m+1];
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
 
    fscanf(magfile,"%d %d %d ",&fdYl, &fdXl, &fdZl); 
    
    printf("fdYl %d, fdXl %d, fdZl %d\n",fdYl, fdXl, fdZl);     
    
    for (ik=0; ik<fdZl; ik++)
    { 
   
      fscanf(magfile, " %e ", &zz);
      fdZc[ik]=zz; 

    } 
   
    for (ik=0; ik<fdYl; ik++)
    {    
       fscanf(magfile, " %e ", &yy); 
       fdY[ik]=yy;
 
    } 
    for (ik=0; ik<81; ik++)
    {    
           printf("fdZc %e,fdY %e\n", fdZc[ik],fdY[ik]); 
    }   
             
    fscanf(magfile," %e %e %e %e %e %e %e %e %e %e %e ", &fdYdel,&fdXdel,&fdZdel,&fdZmax,&fdZmin,&fdYmax,&fdYmin,&fdAx1,&fdCx1,&fdAx2,&fdCx2); 

printf("fdYdel %e, fdXdel %e, fdZdel %e\n",fdYdel,fdXdel,fdZdel);
printf("fdZmax %e, fdZmin %e, fdYmax %e,fdYmin %e\n",fdZmax,fdZmin,fdYmax,fdYmin);
printf("fdAx1 %e, fdCx1 %e, fdAx2 %e, fdCx %e\n",fdAx1,fdCx1,fdAx2,fdCx2);

    for (il=0; il<44; il++)  { 
     for (im=0; im<81; im++)  {      
      for (ik=0; ik<81; ik++)  {      
      
      fscanf(magfile, " %e ", &by); 
      fdBcy[ik][im][il]=by;        
      }
     }
    } 

    for (il=0; il<44; il++)  { 
     for (im=0; im<81; im++)  {      
      for (ik=0; ik<81; ik++)  {      
      
      fscanf(magfile, " %e ", &bx); 
      fdBcx[ik][im][il]=bx;        
      }    
     }     
    }
    
   for (il=0; il<44; il++)  { 
     for (im=0; im<81; im++)  {      
      for (ik=0; ik<81; ik++)  {      
      
      fscanf(magfile, " %e ", &bz); 
      fdBcz[ik][im][il]=bz;          
      }              
     }     
    } 
//----------------------   Polar part ---------------------------------

    printf("Polar part\n");
    fscanf(magfile,"%d %d %d ", &fdZpl, &fdRn, &fdPhin); 
    printf("fdZpl %d, fdRn %d, fdPhin %d\n",fdZpl,fdRn,fdPhin);   

    printf(" fdZp array\n"); 
     
    for (ik=0; ik<51; ik++) 
    {    
     fscanf(magfile, " %e ", &zzp);
     fdZp[ik]=zzp; 
     printf(" %e\n",fdZp[ik]);      
    } 
  
    printf(" fdR array\n"); 
         
    for (ik=0; ik<10; ik++) 
    {    
     fscanf(magfile, " %e ", &rr); 
     fdR[ik]=rr;
     printf(" %e\n",fdR[ik]);
    } 
    
//    printf("fdPhi array\n"); 
     
     for (il=0; il<33; il++)  
     {
       fscanf(magfile, " %e ", &phii); 
       fdPhi[il]=phii; 
//        printf(" %e\n",fdPhi[il]);          
     }

    fscanf(magfile," %e %e %e %e %e %e %e ",&fdZpdl,&fdPhid,&fdRdel,&fdZpmx,&fdZpmn,&fdRmax, &fdRmin); 

printf("fdZpdl %e, fdPhid %e, fdRdel %e, fdZpmx %e, fdZpmn %e,fdRmax %e,fdRmin %e \n", fdZpdl,fdPhid, fdRdel,fdZpmx, fdZpmn,fdRmax, fdRmin);

                      
    for (il=0; il<33; il++)  { 
     for (im=0; im<10; im++)  {      
      for (ik=0; ik<51; ik++)  {
      fscanf(magfile, " %e ", &by); 
        fdBpy[ik][im][il]=by;        
      }
     }
    } 
    
    for (il=0; il<33; il++)  { 
     for (im=0; im<10; im++)  {      
      for (ik=0; ik<51; ik++)  {
      fscanf(magfile, " %e ", &bx); 
        fdBpx[ik][im][il]=bx;                     
      }    
     }     
    }      


    for (il=0; il<33; il++)  { 
     for (im=0; im<10; im++)  {      
      for (ik=0; ik<51; ik++)  {
      fscanf(magfile, " %e ", &bz); 
        fdBpz[ik][im][il]=bz;                      
      }              
     }     
    } 
    
    
    for (il=0; il<32; il++) { 
     for (im=0; im<2; im++)  {      
      for (ik=0; ik<2; ik++)  {
      fscanf(magfile, " %e ", &bb);    
        fdB[ik][im][il]=bb;  
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
