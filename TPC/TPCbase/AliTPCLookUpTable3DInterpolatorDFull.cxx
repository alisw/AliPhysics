#include "AliTPCLookUpTable3DInterpolatorDFull.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TDecompLU.h"


/// \cond CLASSIMP3
ClassImp(AliTPCLookUpTable3DInterpolatorDFull)
/// \endcond



AliTPCLookUpTable3DInterpolatorDFull::AliTPCLookUpTable3DInterpolatorDFull()
{	
	fOrder = 1;
	fIsAllocatingLookUp = kFALSE;
}



AliTPCLookUpTable3DInterpolatorDFull::AliTPCLookUpTable3DInterpolatorDFull
(
	Int_t nr, 
	Double_t rmin, 
	Double_t rmax, 
	Int_t nphi, 
	Double_t phimin, 
	Double_t phimax, 
	Int_t nz , 
	Double_t zmin, 
	Double_t zmax
)
{
	fOrder = 1;
	fIsAllocatingLookUp = kTRUE;
	
	fNR = nr;
	fNPhi = nphi;
	fNZ = nz;
	
	fLookUpR = new TMatrixD*[fNPhi];
	fLookUpPhi = new TMatrixD*[fNPhi];
	fLookUpZ = new TMatrixD*[fNPhi];
	
	for (Int_t m=0;m < fNPhi;m++) {
		fLookUpR[m] = new TMatrixD(fNR,fNZ);   ///< Array to store distterortion followind the drift
		fLookUpPhi[m] = new TMatrixD(fNR,fNZ); ;   ///< Array to store distortion followind the drift
		fLookUpZ[m] = new TMatrixD(fNR,fNZ); ;   ///< Array to store distortion followind the drift
	}	

	// list now contains different position
	fRlist = new TMatrixD *[fNPhi];
	fPhilist = new TMatrixD *[fNPhi];
	fZlist = new TMatrixD *[fNPhi];


	/**
	Double_t dr = (rmax - rmin) / fNR;
	Double_t dphi = (phimax - phimin) / fNPhi;
	Double_t dz = (zmax - zmin) / fNPhi;
	
	for (Int_t m=0;m<fNPhi;m++)  fPhilist[m] = phimin + dphi * m;
	for (Int_t m=0;m<fNR;m++)  fRlist[m] = rmin + dr * m;
	for (Int_t m=0;m<fNZ;m++)  fZlist[m] = zmin + dz * m;
	**/
	
	
	
}


AliTPCLookUpTable3DInterpolatorDFull::AliTPCLookUpTable3DInterpolatorDFull
(
	Int_t nr, 
	TMatrixD**lookupr, 
	TMatrixD**rlist, 
	Double_t * rlistn,
	Int_t nphi, 
	TMatrixD**lookupphi, 
	TMatrixD**philist,
	Double_t * philistn,
	Int_t nz, 
	TMatrixD**lookupz, 
	TMatrixD**zlist,
	Double_t * zlistn,
	Int_t order,
	Int_t stepR,
	Int_t stepZ,
	Int_t stepPhi,
	Int_t type
)
{
	fIsAllocatingLookUp = kFALSE;
	
	
	SetNR(nr);
	SetLookUpR(lookupr);	
	SetRList(rlist);
	
	SetNPhi(nphi);
	SetLookUpPhi(lookupphi);
	SetPhiList(philist);
	
	SetNZ(nz);
	SetLookUpZ(lookupz);
	SetZList(zlist);	
	SetOrder(order);


	fInterpolatorR = new AliTPC3DCylindricalInterpolatorFull(nr,nz,nphi,stepR,stepZ,stepPhi,type);
	fInterpolatorZ = new AliTPC3DCylindricalInterpolatorFull(nr,nz,nphi,stepR,stepZ,stepPhi,type);
	fInterpolatorPhi = new AliTPC3DCylindricalInterpolatorFull(nr,nz,nphi,stepR,stepZ,stepPhi,type);
	
	fInterpolatorR->SetNR(nr);
	fInterpolatorR->SetNZ(nz);
	fInterpolatorR->SetNPhi(nphi);

	
	fInterpolatorR->SetRListNormalized(rlistn);
	fInterpolatorR->SetZListNormalized(zlistn);
	fInterpolatorR->SetPhiListNormalized(philistn);
       
	
	fInterpolatorR->SetOrder(order);
	
	//interpolatorR->SetValue(lookupr);
	
	fInterpolatorZ->SetNR(nr);
	fInterpolatorZ->SetNZ(nz);
	fInterpolatorZ->SetNPhi(nphi);

	
	fInterpolatorZ->SetRListNormalized(rlistn);
	fInterpolatorZ->SetZListNormalized(zlistn);
	fInterpolatorZ->SetPhiListNormalized(philistn);
       

	
	fInterpolatorZ->SetOrder(order);
	
	//interpolatorZ->SetValue(lookupz);
	
	
	fInterpolatorPhi->SetNR(nr);
	fInterpolatorPhi->SetNZ(nz);
	fInterpolatorPhi->SetNPhi(nphi);

	
	fInterpolatorPhi->SetRListNormalized(rlistn);
	fInterpolatorPhi->SetZListNormalized(zlistn);
	fInterpolatorPhi->SetPhiListNormalized(philistn);
       
	fInterpolatorPhi->SetOrder(order);
	//interpolatorZ->SetValue(lookupphi);
	
}

AliTPCLookUpTable3DInterpolatorDFull::~AliTPCLookUpTable3DInterpolatorDFull()
{
	
	if  (fIsAllocatingLookUp) {
		for (Int_t m=0;m < fNPhi;m++) {
			delete fLookUpR[m];
			delete fLookUpPhi[m];
			delete fLookUpZ[m];
			delete fRlist[m];
			delete fPhilist[m];
			delete fZlist[m];
		}
		delete fLookUpR;
		delete fLookUpPhi;
		delete fLookUpZ;
		delete fRlist;
		delete fPhilist;
		delete fZlist;
	
	}

	delete fInterpolatorR;
	delete fInterpolatorZ;
	delete fInterpolatorPhi;
	
}


Double_t AliTPCLookUpTable3DInterpolatorDFull::Interpolate3DTableCyl
( 
	Int_t order, 
	Double_t r,   
	Double_t z,
	Double_t phi,   	
	Int_t  nr,    
	Int_t  nz,
	Int_t  nphi,    	
	const Double_t rlist[], 
	const Double_t zlist[],
	const Double_t philist[], 
	
	TMatrixD **arrayofArrays 
) 
{
  /// Interpolate table (TMatrix format) - 3D interpolation
  /// Float version (in order to decrease the OCDB size)

  static  Int_t ilow = 0, jlow = 0, klow = 0, m=0;
  Float_t saveArray[6]= {0.,0.,0.,0.,0.,0.};
  Float_t savedArray[6]= {0.,0.,0.,0.,0.,0.} ;

	while (phi < 0.0) phi = TMath::TwoPi() + phi;
	while (phi > TMath::TwoPi()) phi = phi - TMath::TwoPi() ;

  Search( nr, rlist, r, ilow   ) ;
  Search( nz, zlist, z, jlow   ) ;  
  Search( nphi, philist, phi, klow   ) ;

	
	
  if ( ilow < 0 ) ilow = 0 ;   // check if out of range
  if ( jlow < 0 ) jlow = 0 ;  
	if ( klow < 0 ) klow = nphi + klow ;
  
  // cubic spline
	if (order > 2) {
		if (ilow >= (order-2))
		   ilow = ilow - (order - 2);
		if (jlow >= (order-2))
		   jlow = jlow - (order - 2);
		klow = klow - (order - 2);		   			
		if ( klow < 0 ) klow = nphi + klow ;	   
	}

	//if ( klow < 0 ) {
	//	klow = nphi + klow ;
	//}

	

  if ( ilow + order  >=    nr - 1 ) ilow =   nr- 1 - order ;
  if ( jlow + order  >=    nz - 1 ) jlow =   nz - 1 - order ;
  
  
  
  

  for ( Int_t k = 0 ; k <  order + 1 ; k++ )
    {
			m = (klow + k) % nphi;
			
      TMatrixD &table = *arrayofArrays[m] ;
      
      for ( Int_t i = ilow ; i < ilow + order + 1 ; i++ )
			{
				saveArray[i-ilow] = Interpolate( &zlist[jlow], &table(i,jlow), order, z )   ;
			}
      savedArray[k] = Interpolate( &rlist[ilow], saveArray, order, r )  ;
      //table.Print();
    }
  return( InterpolatePhi( &philist[0], klow, nphi,  savedArray, order, phi ) )   ;
}


Double_t AliTPCLookUpTable3DInterpolatorDFull::InterpolatePhi
( 
	const Double_t xArray[], 
	const Int_t ilow,
	const Int_t nx,
	const Float_t yArray[],
	Int_t order, 
	Double_t x 
)
{
  /// Interpolate function Y(x) using linear (order=1) or quadratic (order=2) interpolation.

	Int_t i0 = ilow;
	Double_t xi0 = xArray[ilow];
	
	Int_t i1 = (ilow + 1) % nx;
	Double_t xi1 = xArray[i1];
	Int_t i2 = (ilow + 2) % nx;
	Double_t xi2 = xArray[i2];
	if ((ilow + 1) >= nx) {
		xi1 += TMath::TwoPi();
	}
	
	
	if ((ilow + 2) >= nx) {
		xi2 += TMath::TwoPi();
		
	}
  Double_t y ;
  
  if (order > 2) {
		Double_t dphi = xArray[1] - xArray[0];
		
		Double_t * philist = new Double_t[order + 1];
		
		for (Int_t i=0;i < order +1;i++)
			philist[i] = xArray[ilow] + i*dphi;
		
		Double_t * y2Array = new Double_t[order + 1];
		
		Spline3(philist,yArray,order+1,y2Array);
		if (x < philist[0]) x = TMath::TwoPi() + x;
		y = SplineInt3(philist,yArray,y2Array,order+1,x);		
//		printf("phi:(%f,%f) (%f,%f)\n",philist[0],philist[order],x,y);
		delete[] philist;
		delete[] y2Array;
	}
	else if ( order == 2 ) {                // Quadratic Interpolation = 2		
    y  = (x-xi1) * (x-xi2) * yArray[0] / ( (xi0-xi1) * (xi0-xi2) ) ;
    y += (x-xi2) * (x-xi0) * yArray[1] / ( (xi1-xi2) * (xi1-xi0) ) ;
    y += (x-xi0) * (x-xi1) * yArray[2] / ( (xi2-xi0) * (xi2-xi1) ) ;
  } else {                           // Linear Interpolation = 1
    y  = yArray[0] + ( yArray[1]-yArray[0] ) * ( x-xArray[i0] ) / (xi1 - xArray[i0] ) ;
  }

  return (y);

}



void AliTPCLookUpTable3DInterpolatorDFull::Search
( 
	Int_t n, 
	const Double_t xArray[], 
	Double_t x, 
	Int_t &low 
) 
{
  /// Search an ordered table by starting at the most recently used point

  Long_t middle, high ;
  Int_t  ascend = 0, increment = 1 ;

  if ( xArray[n-1] >= xArray[0] ) ascend = 1 ;  // Ascending ordered table if true

  if ( low < 0 || low > n-1 ) {
    low = -1 ; high = n ;
  } else {                                            // Ordered Search phase
    if ( (Int_t)( x >= xArray[low] ) == ascend )  {
      if ( low == n-1 ) return ;
      high = low + 1 ;
      while ( (Int_t)( x >= xArray[high] ) == ascend ) {
	low = high ;
	increment *= 2 ;
	high = low + increment ;
	if ( high > n-1 )  {  high = n ; break ;  }
      }
    } else {
      if ( low == 0 )  {  low = -1 ;  return ;  }
      high = low - 1 ;
      while ( (Int_t)( x < xArray[low] ) == ascend ) {
	high = low ;
	increment *= 2 ;
	if ( increment >= high )  {  low = -1 ;  break ;  }
	else  low = high - increment ;
      }
    }
  }

  while ( (high-low) != 1 ) {                     // Binary Search Phase
    middle = ( high + low ) / 2 ;
    if ( (Int_t)( x >= xArray[middle] ) == ascend )
      low = middle ;
    else
      high = middle ;
  }

  if ( x == xArray[n-1] ) low = n-2 ;
  if ( x == xArray[0]   ) low = 0 ;

}


Float_t AliTPCLookUpTable3DInterpolatorDFull::Interpolate
( 
	const Double_t xArray[], 
	const Float_t yArray[],
	Int_t order, 
	Double_t x 
) 
{
  /// Interpolate function Y(x) using linear (order=1) or quadratic (order=2) interpolation.
  /// Float version (in order to decrease the OCDB size)
  

  Float_t y ;
  
  if ( order > 2) {
		
		Double_t * y2Array = new Double_t[order + 1];
		Spline3(xArray,yArray,order+1,y2Array);
		y = SplineInt3(xArray,yArray,y2Array,order+1,x);		
		//printf("(%f,%f)(%f,%f)\n",xArray[0],xArray[order],x,y);
		delete[] y2Array;
		
	}else  if ( order == 2 ) {                // Quadratic Interpolation = 2
    y  = (x-xArray[1]) * (x-xArray[2]) * yArray[0] / ( (xArray[0]-xArray[1]) * (xArray[0]-xArray[2]) ) ;
    y += (x-xArray[2]) * (x-xArray[0]) * yArray[1] / ( (xArray[1]-xArray[2]) * (xArray[1]-xArray[0]) ) ;
    y += (x-xArray[0]) * (x-xArray[1]) * yArray[2] / ( (xArray[2]-xArray[0]) * (xArray[2]-xArray[1]) ) ;
  } else {                           // Linear Interpolation = 1
    y  = yArray[0] + ( yArray[1]-yArray[0] ) * ( x-xArray[0] ) / ( xArray[1] - xArray[0] ) ;
  }

  return (y);

}


Double_t AliTPCLookUpTable3DInterpolatorDFull::Interpolate
( 
	const Double_t xArray[], 
	const Double_t yArray[],
	Int_t order, 
	Double_t x 
) 
{
  /// Interpolate function Y(x) using linear (order=1) or quadratic (order=2) interpolation.
  
 
  Double_t y ;
  
  if ( order > 2) {
		Double_t * y2Array = new Double_t[order + 1];
		Spline3(xArray,yArray,order+1,y2Array);
		y = SplineInt3(xArray,yArray,y2Array,order+1,x);		
		//printf("(%f,%f)(%f,%f)\n",xArray[0],xArray[order],x,y);
		delete[] y2Array;
		
	} else  if ( order == 2 ) {                // Quadratic Interpolation = 2
    y  = (x-xArray[1]) * (x-xArray[2]) * yArray[0] / ( (xArray[0]-xArray[1]) * (xArray[0]-xArray[2]) ) ;
    y += (x-xArray[2]) * (x-xArray[0]) * yArray[1] / ( (xArray[1]-xArray[2]) * (xArray[1]-xArray[0]) ) ;
    y += (x-xArray[0]) * (x-xArray[1]) * yArray[2] / ( (xArray[2]-xArray[0]) * (xArray[2]-xArray[1]) ) ;
  } else {                           // Linear Interpolation = 1
    y  = yArray[0] + ( yArray[1]-yArray[0] ) * ( x-xArray[0] ) / ( xArray[1] - xArray[0] ) ;
  }

  return (y);

}

void AliTPCLookUpTable3DInterpolatorDFull::CopyVals()  {
	fInterpolatorR->SetVals(fLookUpR,fRlist,fPhilist,fZlist);	
	fInterpolatorZ->SetVals(fLookUpZ,fRlist,fPhilist,fZlist);
	fInterpolatorPhi->SetVals(fLookUpPhi,fRlist,fPhilist,fZlist);
}


void AliTPCLookUpTable3DInterpolatorDFull::CopyValsCartesian()  {
	fInterpolatorR->SetValsCartesian(fLookUpR,fRlist,fPhilist,fZlist);	
	fInterpolatorZ->SetValsCartesian(fLookUpZ,fRlist,fPhilist,fZlist);
	fInterpolatorPhi->SetValsCartesian(fLookUpPhi,fRlist,fPhilist,fZlist);
}

// copy vals irregular grid at z == j
// 
void AliTPCLookUpTable3DInterpolatorDFull::CopyVals(Int_t j)  {
  fInterpolatorR->SetVals(fLookUpR,fRlist,fPhilist,fZlist,j);	
  fInterpolatorZ->SetVals(fLookUpZ,fRlist,fPhilist,fZlist,j);
  fInterpolatorPhi->SetVals(fLookUpPhi,fRlist,fPhilist,fZlist,j);

  /**
  if (fOrder > 2) 
    {
      fInterpolatorR->InitCubicSpline(j);
      fInterpolatorZ->InitCubicSpline(j);
      fInterpolatorPhi->InitCubicSpline(j);
    }
  **/
	
}


void AliTPCLookUpTable3DInterpolatorDFull::GetValue
(
	Double_t r, 
	Double_t phi, 
	Double_t z,
	Double_t &vr,
	Double_t &vphi,
	Double_t &vz,
	Int_t rindex,
	Int_t phiindex,
	Int_t zindex,
	Int_t stepR,
	Int_t stepPhi,
	Int_t stepZ
) 
{
	
	//vr = InterpolateCylindrical(fOrder, r,z,phi, fNR, fNZ, fNPhi, fListR, fListZ, fListPhi, fLookUpR);
  vr =  fInterpolatorR->GetValue(r,phi,z,rindex,phiindex,zindex,stepR,stepPhi,stepZ);
	//vphi = InterpolateCylindrical(fOrder, r,z,phi, fNR, fNZ, fNPhi, fListR, fListZ, fListPhi, fLookUpPhi);
  vphi =  fInterpolatorPhi->GetValue(r,phi,z,rindex,phiindex,zindex,stepR,stepPhi,stepZ);
	//InterpolateCylindrical(fOrder, r,z,phi, fNR, fNZ, fNPhi, fListR, fListZ, fListPhi, fLookUpPhi);
	//vz = InterpolateCylindrical(fOrder, r,z,phi, fNR, fNZ, fNPhi, fListR, fListZ, fListPhi, fLookUpZ);
  vz =  fInterpolatorZ->GetValue(r,phi,z,rindex,phiindex,zindex,stepR,stepPhi,stepZ);
	//InterpolateCylindrical(fOrder, r,z,phi, fNR, fNZ, fNPhi, fListR, fListZ, fListPhi, fLookUpZ);
	
	
}



void AliTPCLookUpTable3DInterpolatorDFull::GetValueCartesian
(
	Double_t r, 
	Double_t phi, 
	Double_t z,
	Double_t &vr,
	Double_t &vphi,
	Double_t &vz,
	Int_t rindex,
	Int_t phiindex,
	Int_t zindex,
	Int_t stepR,
	Int_t stepPhi,
	Int_t stepZ
) 
{
	
	//vr = InterpolateCylindrical(fOrder, r,z,phi, fNR, fNZ, fNPhi, fListR, fListZ, fListPhi, fLookUpR);
  vr =  fInterpolatorR->GetValueCartesian(r,phi,z,rindex,phiindex,zindex,stepR,stepPhi,stepZ);
	//vphi = InterpolateCylindrical(fOrder, r,z,phi, fNR, fNZ, fNPhi, fListR, fListZ, fListPhi, fLookUpPhi);
  vphi =  fInterpolatorPhi->GetValueCartesian(r,phi,z,rindex,phiindex,zindex,stepR,stepPhi,stepZ);
	//InterpolateCylindrical(fOrder, r,z,phi, fNR, fNZ, fNPhi, fListR, fListZ, fListPhi, fLookUpPhi);
	//vz = InterpolateCylindrical(fOrder, r,z,phi, fNR, fNZ, fNPhi, fListR, fListZ, fListPhi, fLookUpZ);
  vz =  fInterpolatorZ->GetValueCartesian(r,phi,z,rindex,phiindex,zindex,stepR,stepPhi,stepZ);
	//InterpolateCylindrical(fOrder, r,z,phi, fNR, fNZ, fNPhi, fListR, fListZ, fListPhi, fLookUpZ);
	
	
}

// add max zindex
void AliTPCLookUpTable3DInterpolatorDFull::GetValue
(
	Double_t r, 
	Double_t phi, 
	Double_t z,
	Double_t &vr,
	Double_t &vphi,
	Double_t &vz,
	Int_t rindex,
	Int_t phiindex,
	Int_t zindex,
	Int_t stepR,
	Int_t stepPhi,
	Int_t stepZ,
	Int_t minzindex
) 
{
	
	//vr = InterpolateCylindrical(fOrder, r,z,phi, fNR, fNZ, fNPhi, fListR, fListZ, fListPhi, fLookUpR);
  vr =  fInterpolatorR->GetValue(r,phi,z,rindex,phiindex,zindex,stepR,stepPhi,stepZ,minzindex);
	//vphi = InterpolateCylindrical(fOrder, r,z,phi, fNR, fNZ, fNPhi, fListR, fListZ, fListPhi, fLookUpPhi);
  vphi =  fInterpolatorPhi->GetValue(r,phi,z,rindex,phiindex,zindex,stepR,stepPhi,stepZ,minzindex);
	//InterpolateCylindrical(fOrder, r,z,phi, fNR, fNZ, fNPhi, fListR, fListZ, fListPhi, fLookUpPhi);
	//vz = InterpolateCylindrical(fOrder, r,z,phi, fNR, fNZ, fNPhi, fListR, fListZ, fListPhi, fLookUpZ);
  vz =  fInterpolatorZ->GetValue(r,phi,z,rindex,phiindex,zindex,stepR,stepPhi,stepZ,minzindex);
	//InterpolateCylindrical(fOrder, r,z,phi, fNR, fNZ, fNPhi, fListR, fListZ, fListPhi, fLookUpZ);
	
	
}


// add max zindex
void AliTPCLookUpTable3DInterpolatorDFull::GetValueCartesian
(
	Double_t r, 
	Double_t phi, 
	Double_t z,
	Double_t &vr,
	Double_t &vphi,
	Double_t &vz,
	Int_t rindex,
	Int_t phiindex,
	Int_t zindex,
	Int_t stepR,
	Int_t stepPhi,
	Int_t stepZ,
	Int_t minzindex
) 
{
	
	//vr = InterpolateCylindrical(fOrder, r,z,phi, fNR, fNZ, fNPhi, fListR, fListZ, fListPhi, fLookUpR);
  vr =  fInterpolatorR->GetValueCartesian(r,phi,z,rindex,phiindex,zindex,stepR,stepPhi,stepZ,minzindex);
	//vphi = InterpolateCylindrical(fOrder, r,z,phi, fNR, fNZ, fNPhi, fListR, fListZ, fListPhi, fLookUpPhi);
  vphi =  fInterpolatorPhi->GetValueCartesian(r,phi,z,rindex,phiindex,zindex,stepR,stepPhi,stepZ,minzindex);
	//InterpolateCylindrical(fOrder, r,z,phi, fNR, fNZ, fNPhi, fListR, fListZ, fListPhi, fLookUpPhi);
	//vz = InterpolateCylindrical(fOrder, r,z,phi, fNR, fNZ, fNPhi, fListR, fListZ, fListPhi, fLookUpZ);
  vz =  fInterpolatorZ->GetValueCartesian(r,phi,z,rindex,phiindex,zindex,stepR,stepPhi,stepZ,minzindex);
	//InterpolateCylindrical(fOrder, r,z,phi, fNR, fNZ, fNPhi, fListR, fListZ, fListPhi, fLookUpZ);
	
	
}

void AliTPCLookUpTable3DInterpolatorDFull::GetValue
(
	Double_t r, 
	Double_t phi, 
	Double_t z,
	Double_t &vr,
	Double_t &vphi,
	Double_t &vz
) 
{
	
	//vr = InterpolateCylindrical(fOrder, r,z,phi, fNR, fNZ, fNPhi, fListR, fListZ, fListPhi, fLookUpR);
  vr =  fInterpolatorR->GetValue(r,phi,z);
	//vphi = InterpolateCylindrical(fOrder, r,z,phi, fNR, fNZ, fNPhi, fListR, fListZ, fListPhi, fLookUpPhi);
  vphi =  fInterpolatorPhi->GetValue(r,phi,z);
	//InterpolateCylindrical(fOrder, r,z,phi, fNR, fNZ, fNPhi, fListR, fListZ, fListPhi, fLookUpPhi);
	//vz = InterpolateCylindrical(fOrder, r,z,phi, fNR, fNZ, fNPhi, fListR, fListZ, fListPhi, fLookUpZ);
  vz =  fInterpolatorZ->GetValue(r,phi,z);
	//InterpolateCylindrical(fOrder, r,z,phi, fNR, fNZ, fNPhi, fListR, fListZ, fListPhi, fLookUpZ);
	
	
}




void AliTPCLookUpTable3DInterpolatorDFull::GetValue
(
	Double_t r, 
	Double_t phi, 
	Double_t z,
	Float_t &vr,
	Float_t &vphi,
	Float_t &vz,
	Int_t rindex,
	Int_t phiindex,
	Int_t zindex,
	Int_t startR,
	Int_t startPhi,
	Int_t startZ
) 
{
	//fLookUpR[0]->Print();
	//printf("before ddr=%f,ddrphi=%f,ddz=%f,r=%f,phi=%f,z=%f\n",vr,vphi,vz,r,phi,z);
	//printf("(%d,%d,%d)\n",fNR,fNZ,fNPhi);
	
	
  vr =  fInterpolatorR->GetValue(r,phi,z,rindex,phiindex,zindex,startR,startPhi,startZ);
  vphi =  fInterpolatorPhi->GetValue(r,phi,z,rindex,phiindex,zindex,startR,startPhi,startZ);
  vz =  fInterpolatorZ->GetValue(r,phi,z,rindex,phiindex,zindex,startR,startPhi,startZ);
	//vr   = InterpolateCylindrical(fOrder, r,z,phi, fNR, fNZ, fNPhi, fListR, fListZ, fListPhi, fLookUpR);
	//vphi = InterpolateCylindrical(fOrder, r,z,phi, fNR, fNZ, fNPhi, fListR, fListZ, fListPhi, fLookUpPhi);
	//vz 	 = InterpolateCylindrical(fOrder, r,z,phi, fNR, fNZ, fNPhi, fListR, fListZ, fListPhi, fLookUpZ);
	
	//fInterpolatorR->;
	//printf("after ddr=%f,ddrphi=%f,ddz=%f\n",vr,vphi,vz);
	
	
}

// spline from recepi
void AliTPCLookUpTable3DInterpolatorDFull::Spline3
(	
	const Double_t xArray[], 
	const Double_t yArray[],
	const Int_t n,
	Double_t y2Array[]
) 
{
	Double_t u[n];
	Double_t sig,p,qn,un;
	y2Array[n - 1] = 0.0;
	u[n-1] = 0.0;
	y2Array[0] = 0.0;
	u[0] = 0.0;
	for (Int_t i=1;i<=n-2;i++) {
		sig = (xArray[i] - xArray[i-1])/(xArray[i+1] - xArray[i-1]);
		p = sig * y2Array[i-1]+ 2.0;
		y2Array[i] = (sig - 1.0)/p;
		u[i] = (yArray[i+1] - yArray[i])/(xArray[i+1] - xArray[i]) - (yArray[i] - yArray[i-1])/(xArray[i] - xArray[i-1]);
		u[i] = (6.0 * u[i]/(xArray[i+1] - xArray[i-1]) - sig * u[i-1])/p;	
	}
	qn = un = 0.0;
	y2Array[n-1] = (un - qn * u[n-2])/(qn*y2Array[n-2]+1.0);
	for (Int_t k=n-2;k>=1;k--)
		y2Array[k] = y2Array[k] * y2Array[k+1] + u[k];	
}



// spline from recepi
Double_t AliTPCLookUpTable3DInterpolatorDFull::SplineInt3
(	
	const Double_t xArray[], 
	const Double_t yArray[],
	Double_t y2Array[],	
	const Int_t n,
	Double_t x
) 
{
	Int_t klo,khi,k;
	Float_t h,b,a;
	
	
	klo = 0;
	khi = n-1;
	while (khi-klo > 1) {
		k = (khi + klo) >> 1;
		if (xArray[k] > x) khi=k;
		else klo = k;
	}
	
	h = xArray[khi] - xArray[klo];
	
	
	if (h < 1e-20)  {
		printf("not found\n");
		return 0.0;
	}
	
	a = (xArray[khi] - x) / h;
	b = (x - xArray[klo]) / h;
	
	
	Double_t y = a *yArray[klo] + b * yArray[khi] + ((a*a*a -a) * y2Array[klo] + (b*b*b -b) * y2Array[khi]) * (h*h)/6.0;
	return y;
	
}


// spline from recepi
void AliTPCLookUpTable3DInterpolatorDFull::Spline3
(	
	const Double_t xArray[], 
	const Float_t yArray[],
	const Int_t n,
	Double_t y2Array[]
) 
{
	Double_t u[n];
	Double_t sig,p,qn,un;
	y2Array[n - 1] = 0.0;
	u[n-1] = 0.0;
	y2Array[0] = 0.0;
	u[0] = 0.0;
	for (Int_t i=1;i<=n-2;i++) {
		sig = (xArray[i] - xArray[i-1])/(xArray[i+1] - xArray[i-1]);
		p = sig * y2Array[i-1]+ 2.0;
		y2Array[i] = (sig - 1.0)/p;
		u[i] = (yArray[i+1] - yArray[i])/(xArray[i+1] - xArray[i]) - (yArray[i] - yArray[i-1])/(xArray[i] - xArray[i-1]);
		u[i] = (6.0 * u[i]/(xArray[i+1] - xArray[i-1]) - sig * u[i-1])/p;	
	}
	qn = un = 0.0;
	y2Array[n-1] = (un - qn * u[n-2])/(qn*y2Array[n-2]+1.0);
	for (Int_t k=n-2;k>=1;k--)
		y2Array[k] = y2Array[k] * y2Array[k+1] + u[k];	
}



// spline from recepi
Double_t AliTPCLookUpTable3DInterpolatorDFull::SplineInt3
(	
	const Double_t xArray[], 
	const Float_t yArray[],
	Double_t y2Array[],	
	const Int_t n,
	Double_t x
) 
{
	Int_t klo,khi,k;
	Float_t h,b,a;
	
	
	klo = 0;
	khi = n-1;
	while (khi-klo > 1) {
		k = (khi + klo) >> 1;
		if (xArray[k] > x) khi=k;
		else klo = k;
	}
	
	h = xArray[khi] - xArray[klo];
	
	
	if (h < 1e-20)  {
		printf("not found\n");
		return 0.0;
	}
	
	a = (xArray[khi] - x) / h;
	b = (x - xArray[klo]) / h;
	
	
	Double_t y = a *yArray[klo] + b * yArray[khi] + ((a*a*a -a) * y2Array[klo] + (b*b*b -b) * y2Array[khi]) * (h*h)/6.0;
	return y;
	
}
