//
//
// This class is based on Jim Thomas's implementation of the ExB effects.
// It was shortened and made to comply with AliROOT coding conventions,
// but not changed in functionallity.
//
//


/***********************************************************************
 *
 * Author: Jim Thomas   11/4/2009
 *
 ***********************************************************************
 *
 * Description: 
 *    Utilities for handling SpaceCharge and a few other distortions 
 *    that may occur in the ALICE TPC
 *
 ***********************************************************************/

#include "AliTPCDistortions.h"
#include <AliMagF.h>
#include <TMath.h>

//
//  Standard maps for E and B Field Distortions 
//
//  Note the careful steps in Z around the Central Electrode due to possible discontinuities at CE.  
//  Needed for interpolation tools on the grid.  Also note that whenever we interpolate this grid, 
//  we explicitly do not allow Z to get closer to CE than 0.2 cm.  This gives three points for 
//  quadratic interpolation in all cases (if you need it).

const Double_t AliTPCDistortions::fgkRList[AliTPCDistortions::kNR] =  {   84.0,   84.5,   85.0,   85.5,   86.0,  87.0,    88.0,
			      	   90.0,   92.0,   94.0,   96.0,   98.0,  100.0,  102.0,  104.0,  106.0,  108.0, 
			      	  110.0,  112.0,  114.0,  116.0,  118.0,  120.0,  122.0,  124.0,  126.0,  128.0, 
				  130.0,  132.0,  134.0,  136.0,  138.0,  140.0,  142.0,  144.0,  146.0,  148.0, 
				  150.0,  152.0,  154.0,  156.0,  158.0,  160.0,  162.0,  164.0,  166.0,  168.0, 
			       	  170.0,  172.0,  174.0,  176.0,  178.0,  180.0,  182.0,  184.0,  186.0,  188.0,
				  190.0,  192.0,  194.0,  196.0,  198.0,  200.0,  202.0,  204.0,  206.0,  208.0,
				  210.0,  212.0,  214.0,  216.0,  218.0,  220.0,  222.0,  224.0,  226.0,  228.0,
				  230.0,  232.0,  234.0,  236.0,  238.0,  240.0,  242.0,  244.0,  246.0,  248.0,
		       		  249.0,  249.5,  250.0,  251.5,  252.0  } ;
  
const Double_t AliTPCDistortions::fgkPhiList[AliTPCDistortions::kNPhi] = {  0.0000,           
				  2.0*TMath::Pi()* 1.0/18, 2.0*TMath::Pi()* 2.0/18, 2.0*TMath::Pi()* 3.0/18,   
				  2.0*TMath::Pi()* 4.0/18, 2.0*TMath::Pi()* 5.0/18, 2.0*TMath::Pi()* 6.0/18,  
				  2.0*TMath::Pi()* 7.0/18, 2.0*TMath::Pi()* 8.0/18, 2.0*TMath::Pi()* 9.0/18, 
				  2.0*TMath::Pi()*10.0/18, 2.0*TMath::Pi()*11.0/18, 2.0*TMath::Pi()*12.0/18,  
				  2.0*TMath::Pi()*13.0/18, 2.0*TMath::Pi()*14.0/18, 2.0*TMath::Pi()*15.0/18,  
				  2.0*TMath::Pi()*16.0/18, 2.0*TMath::Pi()*17.0/18, 2.0*TMath::Pi() } ;
                                  // 19 planes of phi - 18+1 so can wrap around

const Double_t AliTPCDistortions::fgkZList[AliTPCDistortions::kNZ]     = 
                               { -249.5, -249.0, -248.5, -248.0, -247.0, -246.0, -245.0, -243.0, -242.0, -241.0,
				 -240.0, -238.0, -236.0, -234.0, -232.0, -230.0, -228.0, -226.0, -224.0, -222.0,
				 -220.0, -218.0, -216.0, -214.0, -212.0, -210.0, -208.0, -206.0, -204.0, -202.0,
				 -200.0, -198.0, -196.0, -194.0, -192.0, -190.0, -188.0, -186.0, -184.0, -182.0,
				 -180.0, -178.0, -176.0, -174.0, -172.0, -170.0, -168.0, -166.0, -164.0, -162.0,
				 -160.0, -158.0, -156.0, -154.0, -152.0, -150.0, -148.0, -146.0, -144.0, -142.0,
				 -140.0, -138.0, -136.0, -134.0, -132.0, -130.0, -128.0, -126.0, -124.0, -122.0,
				 -120.0, -118.0, -116.0, -114.0, -112.0, -110.0, -108.0, -106.0, -104.0, -102.0,
				 -100.0,  -98.0,  -96.0,  -94.0,  -92.0,  -90.0,  -88.0,  -86.0,  -84.0,  -82.0,
				  -80.0,  -78.0,  -76.0,  -74.0,  -72.0,  -70.0,  -68.0,  -66.0,  -64.0,  -62.0,
				  -60.0,  -58.0,  -56.0,  -54.0,  -52.0,  -50.0,  -48.0,  -46.0,  -44.0,  -42.0,
				  -40.0,  -38.0,  -36.0,  -34.0,  -32.0,  -30.0,  -28.0,  -26.0,  -24.0,  -22.0,
			          -20.0,  -18.0,  -16.0,  -14.0,  -12.0,  -10.0,   -8.0,   -6.0,   -4.0,   -2.0,
 			           -1.0,   -0.5,   -0.2,   -0.1,  -0.05,   0.05,    0.1,    0.2,    0.5,    1.0, 
				    2.0,    4.0,    6.0,    8.0,   10.0,   12.0,   14.0,   16.0,   18.0,   20.0, 
			           22.0,   24.0,   26.0,   28.0,   30.0,   32.0,   34.0,   36.0,   38.0,   40.0, 
			           42.0,   44.0,   46.0,   48.0,   50.0,   52.0,   54.0,   56.0,   58.0,   60.0, 
			           62.0,   64.0,   66.0,   68.0,   70.0,   72.0,   74.0,   76.0,   78.0,   80.0, 
			           82.0,   84.0,   86.0,   88.0,   90.0,   92.0,   94.0,   96.0,   98.0,  100.0, 
			          102.0,  104.0,  106.0,  108.0,  110.0,  112.0,  114.0,  116.0,  118.0,  120.0, 
			          122.0,  124.0,  126.0,  128.0,  130.0,  132.0,  134.0,  136.0,  138.0,  140.0, 
			          142.0,  144.0,  146.0,  148.0,  150.0,  152.0,  154.0,  156.0,  158.0,  160.0, 
			          162.0,  164.0,  166.0,  168.0,  170.0,  172.0,  174.0,  176.0,  178.0,  180.0, 
			          182.0,  184.0,  186.0,  188.0,  190.0,  192.0,  194.0,  196.0,  198.0,  200.0,
				  202.0,  204.0,  206.0,  208.0,  210.0,  212.0,  214.0,  216.0,  218.0,  220.0,
				  222.0,  224.0,  226.0,  228.0,  230.0,  232.0,  234.0,  236.0,  238.0,  240.0,
			          242.0,  243.0,  244.0,  245.0,  246.0,  247.0,  248.0,  248.5,  249.0,  249.5   } ;

const Double_t AliTPCDistortions::fgkIFCRadius= 83.06; // Mean Radius of the Inner Field Cage ( 82.43 min,  83.70 max) (cm)
const Double_t AliTPCDistortions::fgkOFCRadius=254.5;  // Mean Radius of the Outer Field Cage (252.55 min, 256.45 max) (cm)
const Double_t AliTPCDistortions::fgkTPC_Z0=249.7;     // Z location of TPC Gated Grid (cm)
const Double_t AliTPCDistortions::fgkZOffSet= 0.2;     // Offset from CE: calculate all distortions closer to CE as if at this point
const Double_t AliTPCDistortions::fgkCathodeV   =-100000.0; // Cathode Voltage (volts)
const Double_t AliTPCDistortions::fgkGG         =-70.0;     // Gating Grid voltage (volts)
const Double_t AliTPCDistortions::fgkAliceDriftV= 2.73;     // Drift Velocity (cm/microSec) Magnitude




AliTPCDistortions* AliTPCDistortions::fgInstance = 0;




//_ singleton implementation __________________________________________________
AliTPCDistortions* AliTPCDistortions::Instance()
{
  //
  // Singleton implementation
  // Returns an instance of this class, it is created if neccessary
  //
  if (fgInstance == 0){
    fgInstance = new AliTPCDistortions();
    //fgInstance->Init();
  }
  return fgInstance;
}





AliTPCDistortions::AliTPCDistortions()
  :fOmegaTau(0.),fT1(1.),fT2(1.),fC0(0.),fC1(1.),fC2(1.),
   fBField(0),fXTwist(0.),fYTwist(0.),fIFCShift(0.),fDeltaVGGA(0.),fDeltaVGGC(0.),
   fJLow(0),fKLow(0)
{ 
  //
  // default constructor (all corrections are set to have zero effect)
  //

  InitIFCShiftDistortion();
  InitGGVoltErrorDistortion();
}

void AliTPCDistortions::RecalculateCs() {
  //
  // helper function to recalculate the coefficients C0-C2
  // after a change of omega-tau or the tensor terms
  //
  fC0=1.                          / ( 1. + fT2*fT2*fOmegaTau*fOmegaTau ) ;
  fC1=fT1*fOmegaTau               / ( 1. + fT1*fT1*fOmegaTau*fOmegaTau ) ;
  fC2=fT2*fT2*fOmegaTau*fOmegaTau / ( 1. + fT2*fT2*fOmegaTau*fOmegaTau ) ;
}

void AliTPCDistortions::UndoTwistDistortion(const Double_t x[],Double_t xprime[],Int_t roc)
{
  //
  // Twist distortion
  //
  // Remove the effects of a simple "twist" of the TPC in the magnet.  If there is
  // an angle between the E and B fields, there will be a distortion in the recorded
  // tracks.  This routine takes out that distortion.
  //
  // The parameters are:
  // twist[0]: twist angle in x-z
  // twist[1]: twist angle in y-z
  //

  Double_t        zdrift ;
  Int_t           sign ;
  
  Double_t z = x[2] ;                                 // Creat temporary copy of x[2]

  if ( roc <  36 ) sign =  1 ;                       // (TPC End A)
  else             sign = -1 ;                       // (TPC End C)
  if ( roc <= 35 && z <  fgkZOffSet ) z =  fgkZOffSet ;    // Protect against discontinuity at CE
  if ( roc >= 36 && z > -fgkZOffSet ) z = -fgkZOffSet ;    // Protect against discontinuity at CE

  zdrift    = sign * ( fgkTPC_Z0 - TMath::Abs(z) ) ;
  
  xprime[0] = x[0] - (     fC1 * fYTwist - fC2 * fXTwist ) * zdrift/1000 ;
  xprime[1] = x[1] - ( -1* fC1 * fXTwist - fC2 * fYTwist ) * zdrift/1000 ;
  xprime[2] = x[2] ;                                       // Subtract to undo the distortion 
}


void AliTPCDistortions::InitIFCShiftDistortion() {
  //      cout << "AliTPCDistortions::IFCShift   Please wait for the tables to fill ...  ~5 seconds" << endl ;
  Int_t nterms = 100 ;
  Double_t r,z;
  for ( Int_t i = 0 ; i < kNZ ; ++i ) {
    z = TMath::Abs( fgkZList[i] ) ;
    for ( Int_t j = 0 ; j < kNR ; ++j ) {
      r = fgkRList[j] ;
      fShiftER[i][j] = 0.0 ; 	    
      Double_t intz = 0.0 ;
      for ( Int_t n = 1 ; n < nterms ; ++n ) {
	Double_t k  = (2*n-1) * TMath::Pi() / fgkTPC_Z0 ;
	Double_t Cn = -4.0 / ( k * fgkTPC_Z0 ) ;
	Double_t numerator =
	  TMath::BesselK0( k*fgkOFCRadius ) * TMath::BesselI1( k*r ) +
	  TMath::BesselK1( k*r )         * TMath::BesselI0( k*fgkOFCRadius ) ;
	Double_t denominator =
	  TMath::BesselK0( k*fgkOFCRadius ) * TMath::BesselI0( k*fgkIFCRadius ) -
	  TMath::BesselK0( k*fgkIFCRadius ) * TMath::BesselI0( k*fgkOFCRadius ) ;
	Double_t zterm = 1 + TMath::Cos( k*z ) ;
	Double_t qwe = numerator / denominator ;
	intz += Cn * zterm * qwe ;
	if ( n>10 && fabs(intz)*1.e-10 > fabs(qwe) ) break;
      }
      if  ( fgkZList[i] < 0 )  intz = -1 * intz ;  // Force AntiSymmetry of solutions in Z
      fShiftER[i][j] = intz ; 	    
    }
  }
}

void AliTPCDistortions::UndoIFCShiftDistortion(const Double_t x[],Double_t xprime[],Int_t roc)
{ 
  //
  // IFC Shift Distortion
  // 
  // The Inner field cage of the TPC may not be perfectly aligned with the outer field cage 
  // of the TPC.  They can be shifted along the Z axis by up to a 1 mm.  This causes a tilting 
  // of the equi-potential lines inside the TPC and therefore a DCA error at the vertex.  
  // The distortion is anti-symmetric in Z. 
  //

  Int_t   order     = 1     ;               // Linear interpolation = 1, Quadratic = 2         

  Double_t  intEr, intEphi ;
  Double_t r, phi, z ;
  Int_t    sign ;

  r      =  TMath::Sqrt( x[0]*x[0] + x[1]*x[1] ) ;
  phi    =  TMath::ATan2(x[1],x[0]) ;
  if ( phi < 0 ) phi += TMath::TwoPi() ;                   // Table uses phi from 0 to 2*Pi
  z      =  x[2] ;                                         // Create temporary copy of x[2]

  if ( roc <  36 ) sign =  1 ;                       // (TPC End A)
  else                   sign = -1 ;                       // (TPC End C)
  if ( roc <= 35 && z <  fgkZOffSet ) z =  fgkZOffSet ;    // Protect against discontinuity at CE
  if ( roc >= 36 && z > -fgkZOffSet ) z = -fgkZOffSet ;    // Protect against discontinuity at CE

  Interpolate2DEdistortion( order, r, z, fShiftER, intEr ) ;
  intEphi = 0.0 ;  // Efield is symmetric in phi

  // Subtract to Undo the distortions
  if ( r > 0.0 ) 
    {
      phi =  phi - fIFCShift*( fC0*intEphi - fC1*intEr ) / r ;      
      r   =  r   - fIFCShift*( fC0*intEr   + fC1*intEphi ) ;  
    }

  xprime[0] = r * TMath::Cos(phi) ;
  xprime[1] = r * TMath::Sin(phi) ;
  xprime[2] = x[2] ;

}

void AliTPCDistortions::InitGGVoltErrorDistortion() {
  //      cout << "AliTPCDistortions::UndoGG VE  Please wait for the tables to fill ...  ~5 seconds" << endl ;
  Double_t r,z;
  Int_t nterms = 100 ;
  for ( Int_t i = 0 ; i < kNZ ; ++i ) {
    z = fgkZList[i] ;
    for ( Int_t j = 0 ; j < kNR ; ++j ) {
      r = fgkRList[j] ;
      fGGVoltErrorER[i][j] = 0.0 ; 	    
      Double_t intz = 0.0 ;
      for ( Int_t n = 1 ; n < nterms ; ++n ) {
	Double_t k    =  n * TMath::Pi() / fgkTPC_Z0 ;
	Double_t ein  =  0 ;                    // Error potential on the IFC
	Double_t eout =  0 ;                    // Error potential on the OFC
	if ( z < 0 ) {
	  ein   =  2.0 * ((z-fgkTPC_Z0)/fgkTPC_Z0) / ( k * (fgkCathodeV - fgkGG) ) ;       
	  eout  =  2.0 * ((z-fgkTPC_Z0)/fgkTPC_Z0) / ( k * (fgkCathodeV - fgkGG) ) ;       
	}
	if ( z == 0 ) continue ;
	if ( z > 0 ) {
	  ein   =  2.0 * ((fgkTPC_Z0-z)/fgkTPC_Z0) / ( k * (fgkCathodeV - fgkGG) ) ;       
	  eout  =  2.0 * ((fgkTPC_Z0-z)/fgkTPC_Z0) / ( k * (fgkCathodeV - fgkGG) ) ;       
	}
	Double_t an   =  ein  * TMath::BesselK0( k*fgkOFCRadius ) - eout * TMath::BesselK0( k*fgkIFCRadius ) ;
	Double_t bn   =  eout * TMath::BesselI0( k*fgkIFCRadius ) - ein  * TMath::BesselI0( k*fgkOFCRadius ) ;
	Double_t numerator =
	  an * TMath::BesselI1( k*r ) - bn * TMath::BesselK1( k*r ) ;
	Double_t denominator =
	  TMath::BesselK0( k*fgkOFCRadius ) * TMath::BesselI0( k*fgkIFCRadius ) -
	  TMath::BesselK0( k*fgkIFCRadius ) * TMath::BesselI0( k*fgkOFCRadius ) ;
	Double_t zterm = TMath::Cos( k*(fgkTPC_Z0-TMath::Abs(z)) ) - 1 ;
	intz += zterm * numerator / denominator ;
	// Assume series converges, break if small terms
	if ( n>10 && fabs(intz)*1.e-10 > fabs(numerator/denominator) ) break;   
      }
      fGGVoltErrorER[i][j] = intz ;
    }
  }
}



void AliTPCDistortions::UndoGGVoltErrorDistortion(const Double_t x[], Double_t xprime[],Int_t roc)
{ 
  //
  // Gated Grid Voltage Error
  //
  // Calculate the effect of having an incorrect voltage on the A or C end plate Gated Grids.
  //
  //  Electrostatic Equations from StarNote SN0253 by Howard Wieman.
  //  Note that we use a funny coordinate system where Z==0 at the GG.
  //
  
  Int_t   order     = 1     ;               // Linear interpolation = 1, Quadratic = 2         
 
  Double_t  intEr, intEphi ;
  Double_t r, phi, z ;
  Int_t    sign ;

  Double_t deltaVGG;
  
  r      =  TMath::Sqrt( x[0]*x[0] + x[1]*x[1] ) ;
  phi    =  TMath::ATan2(x[1],x[0]) ;
  if ( phi < 0 ) phi += TMath::TwoPi() ;                   // Table uses phi from 0 to 2*Pi
  z      =  x[2] ;

  if ( roc <  36 ) sign =  1 ;                       // (TPC End A)
  else             sign = -1 ;                       // (TPC End C)
  if ( roc <= 35 && z <  fgkZOffSet ) z =  fgkZOffSet ;    // Protect against discontinuity at CE
  if ( roc >= 36 && z > -fgkZOffSet ) z = -fgkZOffSet ;    // Protect against discontinuity at CE

  if ( roc <  36 ) deltaVGG = fDeltaVGGA ;            // (TPC End A)
  else             deltaVGG = fDeltaVGGC ;            // (TPC End C)

  Interpolate2DEdistortion( order, r, z, fGGVoltErrorER, intEr ) ;
  intEphi = 0.0 ;  // Efield is symmetric in phi

  // Subtract to Undo the distortions
  if ( r > 0.0 ) 
    {
      phi =  phi - deltaVGG*( fC0*intEphi - fC1*intEr ) / r ;      
      r   =  r   - deltaVGG*( fC0*intEr   + fC1*intEphi ) ;  
    }

  xprime[0] = r * TMath::Cos(phi) ;
  xprime[1] = r * TMath::Sin(phi) ;
  xprime[2] = x[2] ;

}

void AliTPCDistortions::UndoExBShapeDistortion( const Double_t x[], Double_t xprime[], Int_t roc)
{
  //
  // Example of using Ruben's integrals to do ExB Shape distortions
  //
  //

  xprime[0] = x[0] ;  xprime[1] = x[1] ;  xprime[2] = x[2]  ;

  if (!fBField) return;

  Double_t intBStart[3], intBEnd[3] , xStart[3], xEnd[3] ;
  Double_t intBxOverBz, intByOverBz;
  //Double_t denominator ;
  Int_t    sign;
  if ( roc <  36 ) sign =  1 ;                    // (TPC End A)
  else             sign = -1 ;                    // (TPC End C)

  xStart[0] = x[0] ;  xStart[1] = x[1] ;  xStart[2] = x[2]          ;
  xEnd  [0] = x[0] ;  xEnd  [1] = x[1] ;  xEnd  [2] = sign * fgkTPC_Z0 ;
  
  fBField -> GetTPCInt(xStart, intBStart) ;
  fBField -> GetTPCInt(xEnd  , intBEnd  ) ;
  
  if ( TMath::Abs(intBStart[2] - intBEnd[2]) < 0.1 ) return ;     // Protect against divide by zero, below
  
  intBxOverBz = (intBStart[0] - intBEnd[0]) ;  
  intByOverBz = (intBStart[1] - intBEnd[1]) ;  
  
  xprime[0] +=  ( fC2*intBxOverBz - fC1*intByOverBz )  ;
  xprime[1] +=  ( fC2*intByOverBz + fC1*intBxOverBz )  ;
  xprime[2] +=  0.0 ;

}


void AliTPCDistortions::Interpolate2DEdistortion( const Int_t order, const Double_t r, const Double_t z, 
						  const Double_t er[kNZ][kNR], Double_t &er_value )
{
  //
  // Interpolate table - 2D interpolation
  //
  Double_t save_er[10] ;

  Search( kNZ,   fgkZList,  z,   fJLow   ) ;
  Search( kNR,   fgkRList,  r,   fKLow   ) ;
  if ( fJLow < 0 ) fJLow = 0 ;   // check if out of range
  if ( fKLow < 0 ) fKLow = 0 ;
  if ( fJLow + order  >=    kNZ - 1 ) fJLow =   kNZ - 1 - order ;
  if ( fKLow + order  >=    kNR - 1 ) fKLow =   kNR - 1 - order ;

  for ( Int_t j = fJLow ; j < fJLow + order + 1 ; j++ )
    {
      save_er[j-fJLow]     = Interpolate( &fgkRList[fKLow], &er[j][fKLow], order, r )   ;
    }
  er_value = Interpolate( &fgkZList[fJLow], save_er, order, z )   ;

}


Double_t AliTPCDistortions::Interpolate( const Double_t xArray[], const Double_t yArray[], 
				       const Int_t order, const Double_t x )
{
  //
  // Interpolate function Y(x) using linear (order=1) or quadratic (order=2) interpolation.
  //

  Double_t y ;

  if ( order == 2 )                // Quadratic Interpolation = 2 

    {
      y  = (x-xArray[1]) * (x-xArray[2]) * yArray[0] / ( (xArray[0]-xArray[1]) * (xArray[0]-xArray[2]) ) ; 
      y += (x-xArray[2]) * (x-xArray[0]) * yArray[1] / ( (xArray[1]-xArray[2]) * (xArray[1]-xArray[0]) ) ; 
      y += (x-xArray[0]) * (x-xArray[1]) * yArray[2] / ( (xArray[2]-xArray[0]) * (xArray[2]-xArray[1]) ) ; 
      
    }

  else                             // Linear Interpolation = 1

    {
      y  = yArray[0] + ( yArray[1]-yArray[0] ) * ( x-xArray[0] ) / ( xArray[1] - xArray[0] ) ;
    }

  return (y) ;

}


void AliTPCDistortions::Search( const Int_t n, const Double_t xArray[], const Double_t x, Int_t &low )
{
  //
  // Search an ordered table by starting at the most recently used point
  //

  Long_t middle, high ;
  Int_t  ascend = 0, increment = 1 ;

  if ( xArray[n-1] >= xArray[0] ) ascend = 1 ;  // Ascending ordered table if true
  
  if ( low < 0 || low > n-1 ) { low = -1 ; high = n ; }

  else                                            // Ordered Search phase
    {
      if ( (Int_t)( x >= xArray[low] ) == ascend ) 
	{
	  if ( low == n-1 ) return ;          
	  high = low + 1 ;
	  while ( (Int_t)( x >= xArray[high] ) == ascend )  
	    {
	      low = high ;
	      increment *= 2 ;
	      high = low + increment ;
	      if ( high > n-1 )  {  high = n ; break ;  }
	    }
	}
      else
	{
	  if ( low == 0 )  {  low = -1 ;  return ;  }
	  high = low - 1 ;
	  while ( (Int_t)( x < xArray[low] ) == ascend )
	    {
	      high = low ;
	      increment *= 2 ;
	      if ( increment >= high )  {  low = -1 ;  break ;  }
	      else  low = high - increment ;
	    }
	}
    }

  while ( (high-low) != 1 )                      // Binary Search Phase
    {
      middle = ( high + low ) / 2 ;
      if ( (Int_t)( x >= xArray[middle] ) == ascend )
	low = middle ;
      else
	high = middle ;
    }

  if ( x == xArray[n-1] ) low = n-2 ;
  if ( x == xArray[0]   ) low = 0 ;

}

ClassImp(AliTPCDistortions);
