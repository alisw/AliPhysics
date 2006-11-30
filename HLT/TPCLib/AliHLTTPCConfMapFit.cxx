// @(#) $Id$
// Original: AliHLTConfMapFit.cxx,v 1.14 2005/06/14 10:55:21 cvetan 

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "AliHLTTPCRootTypes.h"
#include "AliHLTTPCLogging.h"
#include "AliHLTTPCVertex.h"
#include "AliHLTTPCConfMapTrack.h"
#include "AliHLTTPCConfMapPoint.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCConfMapFit.h"

/** \class AliHLTTPCConfMapFit
<pre>
//_____________________________________________________________
// AliHLTTPCConfMapFit
//
// Fit class for conformal mapping tracking
</pre>
*/

#if __GNUC__ >= 3
using namespace std;
#endif

ClassImp(AliHLTTPCConfMapFit)


AliHLTTPCConfMapFit::AliHLTTPCConfMapFit(AliHLTTPCConfMapTrack *track,AliHLTTPCVertex *vertex)
{
  //constructor
  fTrack = track;
  fVertex = vertex;
}

Int_t AliHLTTPCConfMapFit::FitHelix()
{
  //fit the helix
  if(FitCircle())
    {
      LOG(AliHLTTPCLog::kError,"AliHLTTPCConfMapFit::FitHelix","TrackFit")<<AliHLTTPCLog::kDec<<
	"Problems during circle fit"<<ENDLOG;
      return 1;
    }
  if(FitLine())
    {
      LOG(AliHLTTPCLog::kError,"AliHLTTPCConfMapFit::FitHelix","TrackFit")<<AliHLTTPCLog::kDec<<
	"Problems during line fit"<<ENDLOG;
      return 1;
    }
  return 0;
}

// #### -B0-CHANGE-START == JMT

Int_t AliHLTTPCConfMapFit::FitStraightLine() {
    //fit the straight line 
    if(FitLineXY()) {
	LOG(AliHLTTPCLog::kError,"AliHLTTPCConfMapFit::FitStraightLine","TrackFit")<<AliHLTTPCLog::kDec<<
	    "Problems during stright line fit in XY plane"<<ENDLOG;
	return 1;
    }
    if(FitLineSZ()){
	LOG(AliHLTTPCLog::kError,"AliHLTTPCConfMapFit::FitStraightLine","TrackFit")<<AliHLTTPCLog::kDec<<
	    "Problems during stright line fit in SZ plane"<<ENDLOG;
	return 1;
    }
    return 0;
}

// #### -B0-CHANGE-END == JMT

Int_t AliHLTTPCConfMapFit::FitCircle()
{
  //-----------------------------------------------------------------
  //Fits circle parameters using algorithm
  //described by ChErnov and Oskov in Computer Physics
  //Communications.
  // 
  //Written in FORTRAN by Jawluen Tang, Physics department , UT-Austin 
  //Moved to C by Pablo Yepes
  //Moved to AliROOT by ASV.
  //------------------------------------------------------------------
  
  Double_t wsum  = 0.0 ;
  Double_t xav   = 0.0 ;
  Double_t yav   = 0.0 ;
  
  Int_t num_of_hits = fTrack->GetNumberOfPoints();
  //
  //     Loop over hits calculating average
  Int_t co=0;
  
  for(fTrack->StartLoop(); fTrack->LoopDone(); fTrack->GetNextHit())
    {
      co++;
      AliHLTTPCConfMapPoint *cHit = (AliHLTTPCConfMapPoint*)fTrack->GetCurrentHit();
      cHit->SetXYWeight( 1./ (Double_t)(cHit->GetXerr()*cHit->GetXerr() + cHit->GetYerr()*cHit->GetYerr()) );
      wsum      += cHit->GetXYWeight() ;
      xav       += cHit->GetXYWeight() * cHit->GetX() ;
      yav       += cHit->GetXYWeight() * cHit->GetY() ;
    }
  if(co!=num_of_hits) 
    LOG(AliHLTTPCLog::kError,"AliHLTTPCConfMapFit::FitCircle","TrackFit")<<AliHLTTPCLog::kDec<<
      "Mismatch of hits. Counter: "<<co<<" nHits: "<<num_of_hits<<ENDLOG;
  if (fTrack->ComesFromMainVertex() == true)
    {    
      wsum += fVertex->GetXYWeight() ;
      xav  += fVertex->GetX() ;
      yav  += fVertex->GetY() ;
    }
  
  xav = xav / wsum ;
  yav = yav / wsum ;
//
//  CALCULATE <X**2>, <XY>, AND <Y**2> WITH <X> = 0, & <Y> = 0
//
  Double_t xxav  = 0.0 ;
  Double_t xyav  = 0.0 ; 
  Double_t yyav  = 0.0 ;
  Double_t xi, yi ;

  for(fTrack->StartLoop(); fTrack->LoopDone(); fTrack->GetNextHit())
    { 
      //AliHLTTPCConfMapPoint *cHit = (AliHLTTPCConfMapPoint *)hits->At(hit_counter);
      AliHLTTPCConfMapPoint *cHit = (AliHLTTPCConfMapPoint*)fTrack->GetCurrentHit();
      xi        = cHit->GetX() - xav ;
      yi        = cHit->GetY() - yav ;
      xxav     += xi * xi * cHit->GetXYWeight() ;
      xyav     += xi * yi * cHit->GetXYWeight() ;
      yyav     += yi * yi * cHit->GetXYWeight() ;
    }
  
  if (fTrack->ComesFromMainVertex() == true)
    {
      xi        = fVertex->GetX() - xav ;
      yi        = fVertex->GetY() - yav ;
      xxav     += xi * xi * fVertex->GetXYWeight() ;
      xyav     += xi * yi * fVertex->GetXYWeight() ;
      yyav     += yi * yi * fVertex->GetXYWeight() ; 
    }
  xxav = xxav / wsum ;
  xyav = xyav / wsum ;
  yyav = yyav / wsum ;
//
//-->  ROTATE COORDINATES SO THAT <XY> = 0
//
//-->  SIGN(C**2 - S**2) = SIGN(XXAV - YYAV) >
//-->  &                                     > ==> NEW : (XXAV-YYAV) > 0
//-->  SIGN(S) = SIGN(XYAV)                  >

  Double_t a = fabs( xxav - yyav ) ;
  Double_t b = 4.0 * xyav * xyav ;

  Double_t asqpb  = a * a + b  ;
  Double_t rasqpb = sqrt ( asqpb) ;

  Double_t splus  = 1.0 + a / rasqpb ;
  Double_t sminus = b / (asqpb * splus) ;

  splus  = sqrt (0.5 * splus ) ;
  sminus = sqrt (0.5 * sminus) ;
//
//->  FIRST REQUIRE : SIGN(C**2 - S**2) = SIGN(XXAV - YYAV)
//
  Double_t sinrot, cosrot ;
  if ( xxav <= yyav ) {
	 cosrot = sminus ;
	 sinrot = splus  ;
  }
  else {
	  cosrot = splus ;
	  sinrot = sminus ;
  }
//
//->  REQUIRE : SIGN(S) = SIGN(XYAV) * SIGN(C) (ASSUMING SIGN(C) > 0)
//
  if ( xyav < 0.0 ) sinrot = - sinrot ;
//
//-->  WE NOW HAVE THE SMALLEST ANGLE THAT GUARANTEES <X**2> > <Y**2>
//-->  TO GET THE SIGN OF THE CHARGE RIGHT, THE NEW X-AXIS MUST POINT
//-->  OUTWARD FROM THE ORGIN.  WE ARE FREE TO CHANGE SIGNS OF BOTH
//-->  COSROT AND SINROT SIMULTANEOUSLY TO ACCOMPLISH THIS.
//
//-->  CHOOSE SIGN OF C WISELY TO BE ABLE TO GET THE SIGN OF THE CHARGE
//
  if ( cosrot*xav+sinrot*yav < 0.0 ) {
	  cosrot = -cosrot ;
	  sinrot = -sinrot ;
  }
//
//->  NOW GET <R**2> AND RSCALE= SQRT(<R**2>)
//
  Double_t rrav   = xxav + yyav ;
  Double_t rscale = sqrt(rrav) ;

  xxav   = 0.0 ;
  yyav   = 0.0 ;
  xyav   = 0.0 ;
  Double_t xrrav = 0.0 ;
  Double_t yrrav = 0.0 ;
  Double_t rrrrav  = 0.0 ;

  Double_t xixi, yiyi, riri, wiriri, xold, yold ;
  
  //for (hit_counter=0; hit_counter<num_of_hits; hit_counter++) 
  for(fTrack->StartLoop(); fTrack->LoopDone(); fTrack->GetNextHit())  
    { 
      //AliHLTTPCConfMapPoint *cHit = (AliHLTTPCConfMapPoint*)hits->At(hit_counter);  
      AliHLTTPCConfMapPoint* cHit = (AliHLTTPCConfMapPoint*)fTrack->GetCurrentHit();

      xold = cHit->GetX() - xav ;
      yold = cHit->GetY() - yav ;
      //
      //-->  ROTATE SO THAT <XY> = 0 & DIVIDE BY RSCALE SO THAT <R**2> = 1
      //
      xi = (  cosrot * xold + sinrot * yold ) / rscale ;
      yi = ( -sinrot * xold + cosrot * yold ) / rscale ;
      
      xixi   = xi * xi ;
      yiyi   = yi * yi ;
      riri   = xixi + yiyi ;
      wiriri = cHit->GetXYWeight() * riri ;
      
      xyav   += cHit->GetXYWeight() * xi * yi ;
      xxav   += cHit->GetXYWeight() * xixi ;
      yyav   += cHit->GetXYWeight() * yiyi ;
      
      xrrav  += wiriri * xi ;
      yrrav  += wiriri * yi ;
      rrrrav += wiriri * riri ;
    }
  //
//   Include vertex if required
//
  if (fTrack->ComesFromMainVertex() == true)
    {
	xold = fVertex->GetX() - xav ;
	yold = fVertex->GetY() - yav ;
	//
	//-->  ROTATE SO THAT <XY> = 0 & DIVIDE BY RSCALE SO THAT <R**2> = 1
	//
	xi = (  cosrot * xold + sinrot * yold ) / rscale ;
	yi = ( -sinrot * xold + cosrot * yold ) / rscale ;
	
	xixi   = xi * xi ;
	yiyi   = yi * yi ;
	riri   = xixi + yiyi ;
	wiriri = fVertex->GetXYWeight() * riri ;

	xyav   += fVertex->GetXYWeight() * xi * yi ;
	xxav   += fVertex->GetXYWeight() * xixi ;
	yyav   += fVertex->GetXYWeight() * yiyi ;

	xrrav  += wiriri * xi ;
	yrrav  += wiriri * yi ;
	rrrrav += wiriri * riri ;
  }
  //
  //    
  //
  //-->  DIVIDE BY WSUM TO MAKE AVERAGES
  //
  xxav    = xxav   / wsum ;
  yyav    = yyav   / wsum ;
  xrrav   = xrrav  / wsum ;
  yrrav   = yrrav  / wsum ;
  rrrrav  = rrrrav / wsum ;
  xyav    = xyav   / wsum ;

  Int_t const ntry = 5 ;
//
//-->  USE THESE TO GET THE COEFFICIENTS OF THE 4-TH ORDER POLYNIMIAL
//-->  DON'T PANIC - THE THIRD ORDER TERM IS ZERO !
//
  Double_t xrrxrr = xrrav * xrrav ;
  Double_t yrryrr = yrrav * yrrav ;
  Double_t rrrrm1 = rrrrav - 1.0  ;
  Double_t xxyy   = xxav  * yyav  ;        

  Double_t c0  =          rrrrm1*xxyy - xrrxrr*yyav - yrryrr*xxav ;
  Double_t c1  =        - rrrrm1      + xrrxrr      + yrryrr   - 4.0*xxyy ;        
  Double_t c2  =   4.0  + rrrrm1                               - 4.0*xxyy ;           
  Double_t c4  = - 4.0  ;                
//
//-->  COEFFICIENTS OF THE DERIVATIVE - USED IN NEWTON-RAPHSON ITERATIONS
//
  Double_t c2d =   2.0 * c2 ;
  Double_t c4d =   4.0 * c4 ;
//
//-->  0'TH VALUE OF LAMDA - LINEAR INTERPOLATION BETWEEN P(0) & P(YYAV)
//
//   LAMDA = YYAV * C0 / (C0 + YRRSQ * (XXAV-YYAV))
  Double_t lamda  = 0.0 ;
  Double_t dlamda = 0.0 ;
//
  Double_t chiscl = wsum * rscale * rscale ;
  Double_t dlamax = 0.001 / chiscl ;   
   
  Double_t p, pd ;
  for ( int itry = 1 ; itry <= ntry ; itry++ ) {
     p      = c0 + lamda * (c1 + lamda * (c2 + lamda * lamda * c4 )) ;
     pd     = (c1 + lamda * (c2d + lamda * lamda * c4d)) ;
     dlamda = -p / pd ;
     lamda  = lamda + dlamda ;
     if (fabs(dlamda)<   dlamax) break ;
  }

  Double_t chi2 = (Double_t)(chiscl * lamda) ;
 
  fTrack->SetChiSq1(chi2);
  // Double_t dchisq = chiscl * dlamda ;	     
//
//-->  NOW CALCULATE THE MATRIX ELEMENTS FOR ALPHA, BETA & KAPPA
//
  Double_t h11   = xxav  -     lamda ;
  Double_t h14   = xrrav ;
  Double_t h22   = yyav  -     lamda ; 
  Double_t h24   = yrrav ;
  Double_t h34   = 1.0   + 2.0*lamda ;
  if ( h11 == 0.0 || h22 == 0.0 ){
    LOG(AliHLTTPCLog::kError,"AliHLTTPCConfMapFit::FitCircle","TrackFit")<<AliHLTTPCLog::kDec<<
      "Problems fitting circle"<<ENDLOG;
    return 1 ;
  }
  Double_t rootsq = (h14*h14)/(h11*h11) + 4.0*h34 ;

  Double_t ratio, kappa, beta ;
  if ( fabs(h22) > fabs(h24) ) {
     ratio  = h24 / h22 ;
     rootsq = ratio * ratio + rootsq ;
     kappa = 1.0 / sqrt(rootsq) ;
     beta  = - ratio * kappa ;
  }
  else {
     ratio  = h22 / h24 ;
     rootsq = 1.0 + ratio * ratio * rootsq ;
     beta  = 1.0 / sqrt(rootsq) ;
     if ( h24 > 0 ) beta = - beta ;
     kappa = -ratio * beta ;
  }            
  Double_t alpha = - (h14/h11) * kappa ;
//
//-->  transform these into the lab coordinate system
//-->  first get kappa and back to real dimensions
//
  Double_t kappa1 = kappa / rscale ;
  Double_t dbro   = 0.5   / kappa1 ;
//
//-->  next rotate alpha and beta and scale
//
  Double_t alphar = (cosrot * alpha - sinrot * beta)* dbro ;
  Double_t betar  = (sinrot * alpha + cosrot * beta)* dbro ;
//
//-->  then translate by (xav,yav)
//
  Double_t acent  = (double)(xav - alphar) ;
  Double_t bcent  = (double)(yav - betar ) ;
  Double_t radius = (double)dbro ;
//
//   Get charge
//
  Int_t q = ( ( yrrav < 0 ) ? 1 : -1 ) ;

  fTrack->SetCharge(q);
  
  
  //Set the first point on the track to the space point coordinates of the innermost track
  //This will be updated to lie on the fit later on (AliHLTTPCTrack::UpdateToFirstPoint).
  Double_t x0,y0,psi,pt ;
  AliHLTTPCConfMapPoint *lHit = (AliHLTTPCConfMapPoint*)fTrack->GetLastHit();
  x0 = lHit->GetX();
  y0 = lHit->GetY();
  fTrack->SetFirstPoint(x0,y0,0); //Z-value is set in FitLine

  psi  = (Double_t)atan2(bcent-y0,acent-x0) ;
  psi  = psi + q * AliHLTTPCTransform::PiHalf();
  if ( psi < 0 ) psi = psi + AliHLTTPCTransform::TwoPi();
  pt   = (Double_t)(AliHLTTPCTransform::GetBFieldValue() * radius ) ;
  
  //Update the track parameters with the parameters from this fit:
  fTrack->SetPsi(psi);
  fTrack->SetPt(pt);
  fTrack->SetRadius(radius);
  fTrack->SetCenterX(acent);
  fTrack->SetCenterY(bcent);

  //
//    Get errors from fast fit
//
  //if ( getPara()->getErrors ) getErrorsCircleFit ( acent, bcent, radius ) ;
//
  return 0 ;
  
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//    Fit Line in s-z plane
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t AliHLTTPCConfMapFit::FitLine ( )
{
  //
  //Initialization 
  //
  Double_t sum = 0.F ;
  Double_t ss  = 0.F ;
  Double_t sz  = 0.F ;
  Double_t sss = 0.F ;
  Double_t ssz = 0.F ;
  //
  //find sum , sums ,sumz, sumss 
  // 
  Double_t dx, dy ;
  Double_t radius = (Double_t)(fTrack->GetPt() / AliHLTTPCTransform::GetBFieldValue() ) ;

  //TObjArray *hits = fTrack->GetHits();
  //Int_t num_of_hits = fTrack->GetNumberOfPoints();

  if (0)// fTrack->ComesFromMainVertex() == true ) 
    {
      dx = ((AliHLTTPCConfMapPoint*)fTrack->GetFirstHit())->GetX() - fVertex->GetX();
      dy = ((AliHLTTPCConfMapPoint*)fTrack->GetFirstHit())->GetY() - fVertex->GetY() ;
    }
  else 
    {
      dx = ((AliHLTTPCConfMapPoint *)fTrack->GetFirstHit())->GetX() - ((AliHLTTPCConfMapPoint *)fTrack->GetLastHit())->GetX() ;
      dy = ((AliHLTTPCConfMapPoint *)fTrack->GetFirstHit())->GetY() - ((AliHLTTPCConfMapPoint *)fTrack->GetLastHit())->GetY() ;
      //dx = ((AliHLTTPCConfMapPoint *)hits->First())->GetX() - ((AliHLTTPCConfMapPoint *)hits->Last())->GetX() ;
      //dy = ((AliHLTTPCConfMapPoint *)hits->First())->GetY() - ((AliHLTTPCConfMapPoint *)hits->Last())->GetY() ;
    }
  
  Double_t localPsi = 0.5F * sqrt ( dx*dx + dy*dy ) / radius ;
  Double_t total_s ;
  
  if ( fabs(localPsi) < 1. ) 
    {
      total_s = 2.0 * radius * asin ( localPsi ) ;
    } 
  else 
    { 
      total_s = 2.0 * radius * AliHLTTPCTransform::Pi() ;
    } 
  
  AliHLTTPCConfMapPoint *previousHit = NULL;

  // FtfBaseHit *previousHit = 0  ;
  
  //for ( startLoop() ; done() ; nextHit() ) {
  Double_t dpsi,s;

  //  for(hit_counter=0; hit_counter<num_of_hits; hit_counter++)
  for(fTrack->StartLoop(); fTrack->LoopDone(); fTrack->GetNextHit())  
    {
      // AliHLTTPCConfMapPoint *cHit = (AliHLTTPCConfMapPoint*)hits->At(hit_counter);
      AliHLTTPCConfMapPoint *cHit = (AliHLTTPCConfMapPoint*)fTrack->GetCurrentHit();
      // if ( GetCurrentHit() != GetFirstHit() ) 
      if(cHit != fTrack->GetFirstHit())//  hits->First())
	{
	  dx   = cHit->GetX() - previousHit->GetX() ;
	  dy   = cHit->GetY() - previousHit->GetY() ;
	  dpsi = 0.5 * (Double_t)sqrt ( dx*dx + dy*dy ) / radius ;
	  fTrack->SetPsierr(dpsi);
	  s = previousHit->GetS() - 2.0 * radius * (Double_t)asin ( dpsi ) ;
	  cHit->SetS(s);
	}
      else
	cHit->SetS(total_s);
      //	cHit->s = total_s ;
    
      sum += cHit->GetZWeight() ;
      ss  += cHit->GetZWeight() * cHit->GetS() ;
      sz  += cHit->GetZWeight() * cHit->GetZ() ;
      sss += cHit->GetZWeight() * cHit->GetS() * cHit->GetS() ;
      ssz += cHit->GetZWeight() * cHit->GetS() * cHit->GetZ() ;
      previousHit = cHit ;
    }
  
  Double_t chi2,det = sum * sss - ss * ss;
  if ( fabs(det) < 1e-20)
    { 
      chi2 = 99999.F ;
      fTrack->SetChiSq2(chi2);
      return 0 ;
    }
  
  //Compute the best fitted parameters A,B
  Double_t tanl,z0,dtanl,dz0;

  tanl = (Double_t)((sum * ssz - ss * sz ) / det );
  z0   = (Double_t)((sz * sss - ssz * ss ) / det );

  fTrack->SetTgl(tanl);
  fTrack->SetZ0(z0);
  
  //     calculate chi-square 
  
  chi2 = 0.;
  Double_t r1 ;
  
  //for(hit_counter=0; hit_counter<num_of_hits; hit_counter++)
  for(fTrack->StartLoop(); fTrack->LoopDone(); fTrack->GetNextHit())  
    {
      //AliHLTTPCConfMapPoint *cHit = (AliHLTTPCConfMapPoint*)hits->At(hit_counter);
      AliHLTTPCConfMapPoint *cHit = (AliHLTTPCConfMapPoint*)fTrack->GetCurrentHit();
      r1   = cHit->GetZ() - tanl * cHit->GetS() - z0 ;
      chi2 += (Double_t) ( (Double_t)cHit->GetZWeight() * (r1 * r1) );
    }
  fTrack->SetChiSq2(chi2);
  //
  //     calculate estimated variance
  //      varsq=chi/(double(n)-2.) 
  //     calculate covariance matrix 
  //      siga=sqrt(varsq*sxx/det) 
  //      sigb=sqrt(varsq*sum/det) 
  //
  dtanl = (Double_t) ( sum / det );
  dz0   = (Double_t) ( sss / det );
  
  fTrack->SetTglerr(dtanl);
  fTrack->SetZ0err(dz0);
  
  return 0 ;
} 


// #### -B0-CHANGE-START == JMT

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//    Straight Line Fit  in x-y plane
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t AliHLTTPCConfMapFit::FitLineXY ( ){
    // -----------------------------------------------------------------------------
    // Implementation after Numerical Recipes in C, 2nd Edtion, chapter 15.2, p. 661
    // with y = b*x + a
    // and Data Analysis for Physical Science Students, Luis Lyons, chapter 2.4 p.51 
    // with y = a' + bx' , x' = x - <x>
    // -----------------------------------------------------------------------------

    Double_t S = 0.;
    Double_t Sx = 0.;

    Double_t SPrime = 0.;
    Double_t SxPrime = 0.;
    Double_t SxxPrime = 0.;
    Double_t SyPrime = 0.;
    Double_t SxyPrime = 0.;

    Double_t chi2 = 0.;

    Int_t num_of_hits = fTrack->GetNumberOfPoints();

    Int_t co=0;
    
    // - Loop over hits calculating average : xav
    for(fTrack->StartLoop(); fTrack->LoopDone(); fTrack->GetNextHit()) {
	co++;
	AliHLTTPCConfMapPoint *currentHit = (AliHLTTPCConfMapPoint*)fTrack->GetCurrentHit();
	// ** maybe not necessary, already done in ConfMapPoint
	currentHit->SetXYWeight( 1./ (Double_t)(currentHit->GetXerr()*currentHit->GetXerr() + currentHit->GetYerr()*currentHit->GetYerr()) );
	// **
	S   += currentHit->GetXYWeight();
	Sx  += currentHit->GetXYWeight() * currentHit->GetX();
    }   
    
    if(co!=num_of_hits) 
	LOG(AliHLTTPCLog::kError,"AliHLTTPCConfMapFit::FitLineXY","TrackFit") << "Mismatch of hits. Counter: "<<co<<" nHits: "<<num_of_hits<<ENDLOG;

    Double_t xav = (Double_t)Sx / S;
    
    // Calculate weighted means
    for(fTrack->StartLoop(); fTrack->LoopDone(); fTrack->GetNextHit()) {
	AliHLTTPCConfMapPoint *currentHit = (AliHLTTPCConfMapPoint*)fTrack->GetCurrentHit();

	Double_t xPrime =  currentHit->GetX() - xav;
	SPrime   += currentHit->GetXYWeight();
	SxPrime  += currentHit->GetXYWeight() * xPrime;
	SxxPrime += currentHit->GetXYWeight() * xPrime * xPrime;
	SyPrime  += currentHit->GetXYWeight() * currentHit->GetY();
	SxyPrime += currentHit->GetXYWeight() * xPrime * currentHit->GetY();
    }

    Double_t det = SPrime*SxxPrime + SxPrime*SxPrime;

    if (fabs(det) < 1e-20) { 
	LOG(AliHLTTPCLog::kDebug,"AliHLTTPCConfMapFit::FitLineXY","TrackFit") << "Determinant == 0" << ENDLOG;	
	chi2 = 99999.F ;
	fTrack->SetChiSq1(chi2);
	return -1 ;
    }

    Double_t b   = (Double_t)(SPrime*SxyPrime - SxPrime*SyPrime) / det;        // line parameter b
    Double_t aPrime   = (Double_t)(SxxPrime*SyPrime - SxPrime*SxyPrime) / det; // line parameter a

    Double_t sigma2b = (Double_t)1. / SxxPrime;
    //-- Double_t sigma2aprime = (Double_t)1. /SPrime;

    // Get gradient angle psi of line in xy plane
    Double_t psi  = (Double_t) atan(b) ; 

    // Calculate chi2
    for(fTrack->StartLoop(); fTrack->LoopDone(); fTrack->GetNextHit()) {
	AliHLTTPCConfMapPoint *currentHit = (AliHLTTPCConfMapPoint*)fTrack->GetCurrentHit();
	Double_t tempchi = currentHit->GetY() - aPrime - b*(currentHit->GetX() - xav);
	chi2 += tempchi*tempchi*currentHit->GetXYWeight() ;
    }

    Double_t a = aPrime - b*xav;


    // Set TrackParameter
    fTrack->SetChiSq1(chi2);
    fTrack->SetPsi(psi);
    fTrack->SetPsierr(sigma2b);   
    fTrack->SetCenterX(0.);    // Set to point on the track (for UpdateToFirstPoint)
    fTrack->SetCenterY(a);     // Set to point on the track (for UpdateToFirstPoint)

    //Set the first point on the track to the space point coordinates of the innermost track
    //This will be updated to lie on the fit later on (AliHLTTPCTrack::UpdateToFirstPoint).
    AliHLTTPCConfMapPoint *lastHit = (AliHLTTPCConfMapPoint*)fTrack->GetLastHit();
    Double_t x0 = lastHit->GetX();
    Double_t y0 = lastHit->GetY();
    fTrack->SetFirstPoint(x0,y0,0); //Z-value is set in FitLineSZ
    

    //Set Defaults
    fTrack->SetRadius(-1.);
    fTrack->SetCharge(1);
    fTrack->SetPt(-1.);
  

    return 0;
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//    Straight Line Fit  in s-z plane
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t AliHLTTPCConfMapFit::FitLineSZ ( ){
    // -----------------------------------------------------------------------------
    // Implementation after Numerical Recipes in C, 2nd Edtion, chapter 15.2, p. 661
    // with z = b*s + a
    // and Data Analysis for Physical Science Students, Luis Lyons, chapter 2.4 p.51 
    // with z = a' + bs' , s' = s - <s>
    // -----------------------------------------------------------------------------

    Double_t S = 0.;
    Double_t Ss = 0.;

    Double_t SPrime = 0.;
    Double_t SsPrime = 0.;
    Double_t SssPrime = 0.;
    Double_t SzPrime = 0.;
    Double_t SszPrime = 0.;

    Double_t chi2 = 0.;
    Double_t s = 0.;

    AliHLTTPCConfMapPoint *previousHit = NULL;
  
    // - Loop over hits calculating length in xy-plane: s
    // - Loop over hits calculating average : sav
    for(fTrack->StartLoop(); fTrack->LoopDone(); fTrack->GetNextHit()) {
	AliHLTTPCConfMapPoint *currentHit = (AliHLTTPCConfMapPoint*)fTrack->GetCurrentHit();
	if(currentHit != fTrack->GetFirstHit()) {
	    Double_t dx = currentHit->GetX() - previousHit->GetX() ;
	    Double_t dy = currentHit->GetY() - previousHit->GetY() ;
	    s = previousHit->GetS() - (Double_t)sqrt ( dx*dx + dy*dy );
	}
	else{
	    Double_t dx = ((AliHLTTPCConfMapPoint *)fTrack->GetFirstHit())->GetX() - ((AliHLTTPCConfMapPoint *)fTrack->GetLastHit())->GetX(); 
	    Double_t dy = ((AliHLTTPCConfMapPoint *)fTrack->GetFirstHit())->GetY() - ((AliHLTTPCConfMapPoint *)fTrack->GetLastHit())->GetY();
	    s = (Double_t)sqrt ( dx*dx + dy*dy );
	}

	currentHit->SetS(s);

	S   += currentHit->GetZWeight();
	Ss  += currentHit->GetZWeight() * currentHit->GetS();
    }

    Double_t sav = (Double_t)Ss / S;

    // Calculate weighted means
    for(fTrack->StartLoop(); fTrack->LoopDone(); fTrack->GetNextHit()) {
	AliHLTTPCConfMapPoint *currentHit = (AliHLTTPCConfMapPoint*)fTrack->GetCurrentHit();

	Double_t sPrime =  currentHit->GetS() - sav;
	SPrime   += currentHit->GetZWeight();
	SsPrime  += currentHit->GetZWeight() * sPrime;
	SssPrime += currentHit->GetZWeight() * sPrime * sPrime;
	SzPrime  += currentHit->GetZWeight() * currentHit->GetZ();
	SszPrime += currentHit->GetZWeight() * sPrime * currentHit->GetZ();
    }

    Double_t det = SPrime*SssPrime + SsPrime*SsPrime;

    if (fabs(det) < 1e-20) { 
	LOG(AliHLTTPCLog::kDebug,"AliHLTTPCConfMapFit::FitLineSZ","TrackFit") << "Determinant == 0" << ENDLOG;	
	chi2 = 99999.F ;
	fTrack->SetChiSq2(chi2);
	return -1 ;
    }

    Double_t b   = (Double_t)(SPrime*SszPrime - SsPrime*SzPrime) / det;        // line parameter b
    Double_t aPrime   = (Double_t)(SssPrime*SzPrime - SsPrime*SszPrime) / det; // line parameter a

    Double_t a = aPrime - b*sav;

    Double_t sigma2b = (Double_t) 1. / SssPrime;
    Double_t sigma2aprime = (Double_t) 1. /SPrime;

    Double_t sigma2a = sigma2aprime + sav*sav * sigma2b*sigma2b;

    // Calculate chi2
    for(fTrack->StartLoop(); fTrack->LoopDone(); fTrack->GetNextHit()) {
	AliHLTTPCConfMapPoint *currentHit = (AliHLTTPCConfMapPoint*)fTrack->GetCurrentHit();
	Double_t tempchi = currentHit->GetZ() - aPrime - b*(currentHit->GetS() - sav);
	chi2 += tempchi*tempchi*currentHit->GetZWeight() ;
    }

    // Set TrackParameter
    fTrack->SetChiSq2(chi2);
    fTrack->SetTgl(b);
    fTrack->SetZ0(a);
    fTrack->SetTglerr(sigma2b);
//  fTrack->SetZ0err(sigma2aprime);   // maybe subject to check
    fTrack->SetZ0err(sigma2a);        // maybe subject to check
    return 0;
}

// #### -B0-CHANGE-END == JMT
