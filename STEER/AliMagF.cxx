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

Revision 1.4  2000/03/28 12:40:24  fca
Introduce factor for magnetic field


Revision 1.3  1999/09/29 09:24:29  fca
Introduction of the Copyright and cvs Log

*/


#include "AliMagF.h"
#include "TSystem.h"

#include <stdlib.h>
#include <stdio.h>

//ZDC part -------------------------------------------------------------------

  static const Float_t G1=20.03;
  static const Float_t FDIP=-37.34;
  static const Float_t FDIMU=6.;
  static const Float_t FCORN=11.72;
//
// ZBEG       Beginning of the inner triplet
// D1BEG      Beginning of separator dipole 1
// D2BEG      Beginning of separator dipole 2
// CORBEG     Corrector dipole beginning (because of dimuon arm)
//
  static const Float_t CORBEG=1920,COREND=CORBEG+190, CORRA2=4.5*4.5;
//
  static const Float_t ZBEG=2300;
  static const Float_t Z1BEG=ZBEG+   0,Z1END=Z1BEG+630,Z1RA2=3.5*3.5;
  static const Float_t Z2BEG=ZBEG+ 880,Z2END=Z2BEG+550,Z2RA2=3.5*3.5;
  static const Float_t Z3BEG=ZBEG+1530,Z3END=Z3BEG+550,Z3RA2=3.5*3.5;
  static const Float_t Z4BEG=ZBEG+2430,Z4END=Z4BEG+630,Z4RA2=3.5*3.5;
  static const Float_t D1BEG=5843.5   ,D1END=D1BEG+945,D1RA2=4.5*4.5;
  static const Float_t D2BEG=12113.2  ,D2END=D2BEG+945,D2RA2=4.5*.5;

//ZDC part -------------------------------------------------------------------

ClassImp(AliMagF)

//________________________________________
AliMagF::AliMagF(const char *name, const char *title, const Int_t integ, const Int_t map, 
		 const Float_t factor, const Float_t fmax)
  : TNamed(name,title)
{
  fMap = map;
  fType = Undef;
  fInteg = integ;
  fFactor = factor;
  fMax = fmax;
}

//________________________________________
void AliMagF::Field(Float_t*, Float_t *b)
{
  printf("Undefined MagF Field called, returning 0\n");
  b[0]=b[1]=b[2]=0;
}
      
ClassImp(AliMagFC)

//________________________________________
AliMagFC::AliMagFC(const char *name, const char *title, const Int_t integ, const Int_t map, 
		 const Float_t factor, const Float_t fmax)
  : AliMagF(name,title,integ,map,factor,fmax)
{
  printf("Constant Field %s created: map= %d, factor= %f\n",fName.Data(),map,factor);
  fType = Const;
}

//________________________________________
void AliMagFC::Field(Float_t *x, Float_t *b)
{
  b[0]=b[1]=b[2]=0;
  if(fMap==1) {
    if(TMath::Abs(x[2])<700 && x[0]*x[0]+(x[1]+30)*(x[1]+30) < 560*560) {
      b[2]=2;
    } else {
      if ( 725 <= x[2] && x[2] <= 1225 ) {
	Float_t dz = TMath::Abs(975-x[2])*0.01;
	b[0]=(1-0.1*dz*dz)*7;
      }
      else {
//This is the ZDC part
	Float_t rad2=x[0]*x[0]+x[1]*x[1];
	if(rad2<D2RA2) {
	  if(x[2]>D2BEG) {
	    
//    Separator Dipole D2
	    if(x[2]<D2END) b[1]=FDIP;
	  } else if(x[2]>D1BEG) {
	    
//    Separator Dipole D1
	    if(x[2]<D1END) b[1]=-FDIP;
	  }
	  if(rad2<CORRA2) {

//    First quadrupole of inner triplet de-focussing in x-direction
//    Inner triplet
	    if(x[2]>Z4BEG) {
	      if(x[2]<Z4END) {
	      
//    2430 <-> 3060
		b[0]=-G1*x[1];
		b[1]=-G1*x[0];
	      }
	    } else if(x[2]>Z3BEG) {
	      if(x[2]<Z3END) {

//    1530 <-> 2080
		b[0]=G1*x[1];
		b[1]=G1*x[0];
	      }
	    } else if(x[2]>Z2BEG) {
	      if(x[2]<Z2END) {
	      
//    890 <-> 1430
		b[0]=G1*x[1];
		b[1]=G1*x[0];
	      }
	    } else if(x[2]>Z1BEG) {
	      if(x[2]<Z1END) {

//    0 <->  630
		b[0]=-G1*x[1];
		b[1]=-G1*x[0];
	      }
	    } else if(x[2]>CORBEG) {
	      if(x[2]<COREND) {
//    Corrector dipole (because of dimuon arm)
		b[0]=FCORN;
	      }
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
  } else {
    printf("Invalid field map for constant field %d\n",fMap);
    exit(1);
  }
}
    
ClassImp(AliMagFCM)

//________________________________________
AliMagFCM::AliMagFCM(const char *name, const char *title, const Int_t integ, const Int_t map, 
		 const Float_t factor, const Float_t fmax)
  : AliMagF(name,title,integ,map,factor,fmax)
{
  fType = ConMesh;
  printf("Constant Mesh Field %s created: map= %d, factor= %f, file= %s\n",fName.Data(),map,factor,fTitle.Data());
}

//________________________________________
void AliMagFCM::Field(Float_t *x, Float_t *b)
{
  Double_t ratx, raty, ratz, hix, hiy, hiz, ratx1, raty1, ratz1, 
    bhyhz, bhylz, blyhz, blylz, bhz, blz, xl[3];
  const Double_t one=1;
  Int_t ix, iy, iz;
    
  // --- find the position in the grid ---

  b[0]=b[1]=b[2]=0;
  if(-700<x[2] && x[2]<fZbeg && x[0]*x[0]+(x[1]+30)*(x[1]+30) < 560*560) {
    b[2]=2;
  } else  {
    Bool_t infield=(fZbeg<=x[2] && x[2]<fZbeg+fZdel*(fZn-1)
       &&  ( fXbeg <= TMath::Abs(x[0]) && TMath::Abs(x[0]) < fXbeg+fXdel*(fXn-1) )
       &&  ( fYbeg <= TMath::Abs(x[1]) && TMath::Abs(x[1]) < fYbeg+fYdel*(fYn-1) ));
      if(infield) {
      xl[0]=TMath::Abs(x[0])-fXbeg;
      xl[1]=TMath::Abs(x[1])-fYbeg;
      xl[2]=x[2]-fZbeg;
      
    // --- start with x
    
      hix=xl[0]*fXdeli;
      ratx=hix-int(hix);
      ix=int(hix);
      
      hiy=xl[1]*fYdeli;
      raty=hiy-int(hiy);
      iy=int(hiy);
      
      hiz=xl[2]*fZdeli;
      ratz=hiz-int(hiz);
      iz=int(hiz);
    
      if(fMap==2) {
      // ... simple interpolation
	ratx1=one-ratx;
	raty1=one-raty;
	ratz1=one-ratz;
	bhyhz = Bx(ix  ,iy+1,iz+1)*ratx1+Bx(ix+1,iy+1,iz+1)*ratx;
	bhylz = Bx(ix  ,iy+1,iz  )*ratx1+Bx(ix+1,iy+1,iz  )*ratx;
	blyhz = Bx(ix  ,iy  ,iz+1)*ratx1+Bx(ix+1,iy  ,iz+1)*ratx;
	blylz = Bx(ix  ,iy  ,iz  )*ratx1+Bx(ix+1,iy  ,iz  )*ratx;
	bhz   = blyhz             *raty1+bhyhz             *raty;
	blz   = blylz             *raty1+bhylz             *raty;
	b[0]  = blz               *ratz1+bhz               *ratz;
	//
	bhyhz = By(ix  ,iy+1,iz+1)*ratx1+By(ix+1,iy+1,iz+1)*ratx;
	bhylz = By(ix  ,iy+1,iz  )*ratx1+By(ix+1,iy+1,iz  )*ratx;
	blyhz = By(ix  ,iy  ,iz+1)*ratx1+By(ix+1,iy  ,iz+1)*ratx;
	blylz = By(ix  ,iy  ,iz  )*ratx1+By(ix+1,iy  ,iz  )*ratx;
	bhz   = blyhz             *raty1+bhyhz             *raty;
	blz   = blylz             *raty1+bhylz             *raty;
	b[1]  = blz               *ratz1+bhz               *ratz;
	//
	bhyhz = Bz(ix  ,iy+1,iz+1)*ratx1+Bz(ix+1,iy+1,iz+1)*ratx;
	bhylz = Bz(ix  ,iy+1,iz  )*ratx1+Bz(ix+1,iy+1,iz  )*ratx;
	blyhz = Bz(ix  ,iy  ,iz+1)*ratx1+Bz(ix+1,iy  ,iz+1)*ratx;
	blylz = Bz(ix  ,iy  ,iz  )*ratx1+Bz(ix+1,iy  ,iz  )*ratx;
	bhz   = blyhz             *raty1+bhyhz             *raty;
	blz   = blylz             *raty1+bhylz             *raty;
	b[2]  = blz               *ratz1+bhz               *ratz;
      //printf("ratx,raty,ratz,b[0],b[1],b[2] %f %f %f %f %f %f\n",
      //ratx,raty,ratz,b[0],b[1],b[2]);
      //
    // ... use the dipole symmetry
	if (x[0]*x[1] < 0) b[1]=-b[1];
	if (x[0]<0) b[2]=-b[2];
      } else {
	printf("Invalid field map for constant mesh %d\n",fMap);
      }
    } else {
//This is the ZDC part
      Float_t rad2=x[0]*x[0]+x[1]*x[1];
      if(rad2<D2RA2) {
	if(x[2]>D2BEG) {
	  
//    Separator Dipole D2
	  if(x[2]<D2END) b[1]=FDIP;
	} else if(x[2]>D1BEG) {

//    Separator Dipole D1
	  if(x[2]<D1END) b[1]=-FDIP;
	}
	if(rad2<CORRA2) {

//    First quadrupole of inner triplet de-focussing in x-direction
//    Inner triplet
	  if(x[2]>Z4BEG) {
	    if(x[2]<Z4END) {
	      
//    2430 <-> 3060
	      b[0]=-G1*x[1];
	      b[1]=-G1*x[0];
	    }
	  } else if(x[2]>Z3BEG) {
	    if(x[2]<Z3END) {

//    1530 <-> 2080
	      b[0]=G1*x[1];
	      b[1]=G1*x[0];
	    }
	  } else if(x[2]>Z2BEG) {
	    if(x[2]<Z2END) {

//    890 <-> 1430
	      b[0]=G1*x[1];
	      b[1]=G1*x[0];
	    }
	  } else if(x[2]>Z1BEG) {
	    if(x[2]<Z1END) {

//    0 <->  630
	      b[0]=-G1*x[1];
	      b[1]=-G1*x[0];
	    }
	  } else if(x[2]>CORBEG) {
	    if(x[2]<COREND) {
//    Corrector dipole (because of dimuon arm)
	      b[0]=FCORN;
	    }
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

//________________________________________
void AliMagFCM::ReadField()
{
  FILE *magfile;
  Int_t ix, iy, iz, ipx, ipy, ipz;
  Float_t bx, by, bz;
  char *fname;
  printf("Reading Magnetic Field %s from file %s\n",fName.Data(),fTitle.Data());
  fname = gSystem->ExpandPathName(fTitle.Data());
  magfile=fopen(fname,"r");
  delete [] fname;
  if (magfile) {
    fscanf(magfile,"%d %d %d %f %f %f %f %f %f",
    	   &fXn, &fYn, &fZn, &fXdel, &fYdel, &fZdel, &fXbeg, &fYbeg, &fZbeg);
    printf("fXn %d, fYn %d, fZn %d, fXdel %f, fYdel %f, fZdel %f, fXbeg %f, fYbeg %f, fZbeg %f\n",
    	   fXn, fYn, fZn, fXdel, fYdel, fZdel, fXbeg, fYbeg, fZbeg);
    fXdeli=1./fXdel;
    fYdeli=1./fYdel;
    fZdeli=1./fZdel;
    fB = new TVector(3*fXn*fYn*fZn);
    for (iz=0; iz<fZn; iz++) {
      ipz=iz*3*(fXn*fYn);
      for (iy=0; iy<fYn; iy++) {
	ipy=ipz+iy*3*fXn;
	for (ix=0; ix<fXn; ix++) {
	  ipx=ipy+ix*3;
	  fscanf(magfile,"%f %f %f",&bz,&by,&bx);
	  (*fB)(ipx+2)=bz;
	  (*fB)(ipx+1)=by;
	  (*fB)(ipx  )=bx;
	}
      }
    }
  } else { 
    printf("File %s not found !\n",fTitle.Data());
    exit(1);
  }
}
// ------------------------------------------------------- 

ClassImp(AliMagFDM)

//________________________________________
AliMagFDM::AliMagFDM(const char *name, const char *title, const Int_t integ,
const Int_t map, const Float_t factor, const Float_t fmax)
  : AliMagF(name,title,integ,map,factor,fmax)
  
{
  fType = DipoMap;

  printf("Field Map for Muon Arm from IP till muon filter %s created: map= %d, factor= %f, file=%s\n",fName.Data(),map,factor,fTitle.Data());
  
}

//________________________________________

void AliMagFDM::Field(Float_t *xfi, Float_t *b)
{
  static  const Double_t eps=0.1E-06;
  static  const Double_t pi2=.6283185E+01;
  static  const Double_t one=1;
  static  const Double_t fdYaxi = 0.3;

  static  const    Int_t  kiip=33; 
  static  const    Int_t  miip=0;    
  static  const    Int_t  liip=0;

  static  const    Int_t  kiic=0;
  static  const    Int_t  miic=0;    
  static  const    Int_t  liic=0;       

  static  const Double_t    fdZbg=502.92;  // Start of Map using in z
  static  const Double_t    fdZL3=600;  // Beginning of L3 door in z

  Double_t   x[3];   
  Double_t   xL3[3]; 
  Double_t   bint[3]; 
  
  Double_t r0;
  
  Double_t bbj;
  Int_t Kvar,jb;

  Double_t Zp1, Zp2,Xp1,Xp2,Yp1,Yp2; 
  Double_t Zz1, Zz2,Yy1,Yy2,X2,X1; 
     
// --- start the map fiel from z = 502.92 cm ---

  x[0] = xfi[0];
  x[1] = xfi[1];
  x[2] = xfi[2];
  b[0]=b[1]=b[2]=0;
  //       printf("x[0]  %f,x[1] %f,x[2]  %f\n",x[0],x[1],x[2]); 

  Double_t rr=TMath::Sqrt(x[0]*x[0]+x[1]*x[1]);
           r0=rr/100;
  Double_t Rpmax;
  Rpmax=fdRmax; 
  if ( (-700<x[2] && x[2]<=fdZbg && 
        (x[0]*x[0]+(x[1]+30)*(x[1]+30))< 560*560)
       || (fdZbg<x[2] && x[2]<=fdZL3 && rr>=Rpmax*100) )
       {
        b[2]=2;
       }

  xL3[0]=x[0]/100;
  xL3[1]=(x[1]+30)/100;
  xL3[2]=x[2]/100; 
 
  Double_t xminn=xL3[2]*fdAx1+fdCx1;
  Double_t xmaxx=xL3[2]*fdAx2+fdCx2;
  Double_t Zcmin,Zcmax,Ycmin,Ycmax;

  Zcmin=fdZmin;
  Zcmax=fdZmax;
  Ycmin=fdYmin;
  Ycmax=fdYmax;
          
if ((fdZbg/100<xL3[2] && xL3[2]<Zcmin && r0<Rpmax) || ((Zcmin<=xL3[2] && xL3[2] <= Zcmax ) && (Ycmin<=xL3[1] && xL3[1]<= Ycmax) &&  (xminn <= xL3[0] && xL3[0] <= xmaxx)))
    {
     if(fMap==3) 
      { 
       if (xL3[2]<Zcmin && r0<Rpmax)   
       {
       //---------------------  Polar part ----------------------
       
       Double_t yyp,ph0;
       Int_t  kp0, lp0, mp0;
       Int_t kpi,lpi,mpi;
       Double_t alp1,alp2,alp3;
       Double_t zpz,rp,fip,cphi; 

       kpi=kiip; 
       lpi=liip;
       mpi=miip;
   
       zpz=xL3[2];

       FZ(&zpz, fdZp ,&fdZpdl,&kpi,&kp0,&Zp1 ,&Zp2,&fdZpl) ;

       yyp=xL3[1]- 0.3;
       cphi=yyp/r0;
       ph0=TMath::ACos(cphi);
       if (xL3[0]< 0) {ph0=pi2 - ph0;}  
                                  
       fip=ph0;
       FZ(&fip,fdPhi,&fdPhid ,&mpi,&mp0, &Xp1,&Xp2,&fdPhin); 

       Double_t Rdel;
       Rdel=fdRdel;

       if (r0<= fdRdel)
        {

         if(r0< eps) 
         {

          bint[0]=(Zp1*fdBpx[kp0][0][0] + Zp2*fdBpx[kp0+1][0][0])*10;     
          bint[1]=(Zp1*fdBpy[kp0][0][0] + Zp2*fdBpy[kp0+1][0][0])*10;      
          bint[2]=(Zp1*fdBpz[kp0][0][0] + Zp2*fdBpz[kp0+1][0][0])*10; 

         }  
      
        alp2= fdB[0][0][mp0]*yyp + fdB[0][1][mp0]*xL3[0]; 
        alp3= fdB[1][0][mp0]*yyp + fdB[1][1][mp0]*xL3[0];      
        alp1= one - alp2 - alp3;
     
        for (jb=0; jb<3 ; jb++) 
        {
         Kvar=jb;     
         FRfuncBi(&Kvar,&Zp1,&Zp2,&alp1,&alp2,&alp3, &kp0,&mp0, &bbj);
         bint[jb] = bbj*10 ;
        }  
       }
        else
       {
        rp=r0;

        FZ(&rp,fdR ,&fdRdel,&lpi,&lp0,&Yp1,&Yp2,&fdRn);

        for (jb=0; jb<3 ; jb++) 
        {
         Kvar=jb;
         FGfuncBi(&Zp1,&Zp2,&Yp1,&Yp2,&Xp1,&Xp2,&Kvar,&kp0,&lp0,&mp0,&bbj); 

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
   lci=liic;
   mci=miic;

   xx1 = fdAx1*xL3[2] + fdCx1;
   xx2 = fdAx2*xL3[2] + fdCx2;

   zzc=xL3[2];
   FZ(&zzc, fdZc ,&fdZdel, &kci,&k0, &Zz1, &Zz2, &fdZl);     

   yyc=xL3[1];
   FZ(&yyc, fdY , &fdYdel,&lci, &l0, &Yy1, &Yy2,&fdYl);  

   xXl = fdXl-one; 
   dx = (xx2-xx1)/xXl;
   xxx= xL3[0]-xx1;
   //     xm = xxx/dx; 
   m0 = int(xxx/dx);   

   if(xL3[0]<(xx1+m0*dx) || xL3[0] >(xx1+(m0+1)*dx)) 
    {
      m0=m0+1;   
      printf(" m0 %d, m0+1 %d\n",m0,m0+1);  
    }

   X2=(xL3[0]-( xx1+m0*dx))/dx;
   X1=one-X2;
   m0=m0-1;
   for (jb=3; jb<6; jb++) 
     {
       Kvar=jb;     
       FGfuncBi(&Zz1,&Zz2,&Yy1,&Yy2,&X1,&X2,&Kvar,&k0, &l0, &m0, &bbj); 
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
      if(rad2<D2RA2) {
        if(x[2]>D2BEG) {
          
//    Separator Dipole D2
          if(x[2]<D2END) b[1]=FDIP;
        } else if(x[2]>D1BEG) {

//    Separator Dipole D1
          if(x[2]<D1END) b[1]=-FDIP;
        }
        if(rad2<CORRA2) {

//    First quadrupole of inner triplet de-focussing in x-direction
//    Inner triplet
          if(x[2]>Z4BEG) {
            if(x[2]<Z4END) {
              
//    2430 <-> 3060
              b[0]=-G1*x[1];
              b[1]=-G1*x[0];
            }
          } else if(x[2]>Z3BEG) {
            if(x[2]<Z3END) {

//    1530 <-> 2080
              b[0]=G1*x[1];
              b[1]=G1*x[0];
            }
          } else if(x[2]>Z2BEG) {
            if(x[2]<Z2END) {

//    890 <-> 1430
              b[0]=G1*x[1];
              b[1]=G1*x[0];
            }
          } else if(x[2]>Z1BEG) {
            if(x[2]<Z1END) {

//    0 <->  630
              b[0]=-G1*x[1];
              b[1]=-G1*x[0];
            }
          } else if(x[2]>CORBEG) {
            if(x[2]<COREND) {
//    Corrector dipole (because of dimuon arm)
//            b[0]=FCORN;
              b[0]=-FCORN;
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
  static  const Double_t one=1;
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
      *a1= one - *a2;
      break;     
     }
    }
  }

/*-------------FRfuncBi----------------*/

void AliMagFDM::FRfuncBi(Int_t *kai,Double_t *za1, Double_t *za2, Double_t *al1, Double_t *al2, Double_t *al3, Int_t *ka, Int_t *ma, Double_t *ba)
  
{
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


  




