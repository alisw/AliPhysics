
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

  

