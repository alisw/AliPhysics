// @(#) $Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "AliHLTStandardIncludes.h"
#include "AliHLTRootTypes.h"
#include "AliHLTLogging.h"
#include "AliHLTVertex.h"
#include "AliHLTConfMapTrack.h"
#include "AliHLTConfMapPoint.h"
#include "AliHLTTransform.h"
#include "AliHLTConfMapFit.h"

/** \class AliHLTConfMapFit
<pre>
//_____________________________________________________________
// AliHLTConfMapFit
//
// Fit class for conformal mapping tracking
</pre>
*/

#if __GNUC__ >= 3
using namespace std;
#endif

ClassImp(AliHLTConfMapFit)


AliHLTConfMapFit::AliHLTConfMapFit(AliHLTConfMapTrack *track,AliHLTVertex *vertex)
{
  //constructor
  fTrack = track;
  fVertex = vertex;
}

Int_t AliHLTConfMapFit::FitHelix()
{
  //fit the helix
  if(FitCircle())
    {
      LOG(AliHLTLog::kError,"AliHLTConfMapFit::FitHelix","TrackFit")<<AliHLTLog::kDec<<
	"Problems during circle fit"<<ENDLOG;
      return 1;
    }
  if(FitLine())
    {
      LOG(AliHLTLog::kError,"AliHLTConfMapFit::FitHelix","TrackFit")<<AliHLTLog::kDec<<
	"Problems during line fit"<<ENDLOG;
      return 1;
    }
  return 0;
}

Int_t AliHLTConfMapFit::FitCircle()
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
      AliHLTConfMapPoint *cHit = (AliHLTConfMapPoint*)fTrack->GetCurrentHit();
      cHit->SetXYWeight( 1./ (Double_t)(cHit->GetXerr()*cHit->GetXerr() + cHit->GetYerr()*cHit->GetYerr()) );
      wsum      += cHit->GetXYWeight() ;
      xav       += cHit->GetXYWeight() * cHit->GetX() ;
      yav       += cHit->GetXYWeight() * cHit->GetY() ;
    }
  if(co!=num_of_hits) 
    LOG(AliHLTLog::kError,"AliHLTConfMapFit::FitCircle","TrackFit")<<AliHLTLog::kDec<<
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
      //AliHLTConfMapPoint *cHit = (AliHLTConfMapPoint *)hits->At(hit_counter);
      AliHLTConfMapPoint *cHit = (AliHLTConfMapPoint*)fTrack->GetCurrentHit();
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
      //AliHLTConfMapPoint *cHit = (AliHLTConfMapPoint*)hits->At(hit_counter);  
      AliHLTConfMapPoint* cHit = (AliHLTConfMapPoint*)fTrack->GetCurrentHit();

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
    LOG(AliHLTLog::kError,"AliHLTConfMapFit::FitCircle","TrackFit")<<AliHLTLog::kDec<<
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
  //This will be updated to lie on the fit later on (AliHLTTrack::UpdateToFirstPoint).
  Double_t x0,y0,psi,pt ;
  AliHLTConfMapPoint *lHit = (AliHLTConfMapPoint*)fTrack->GetLastHit();
  x0 = lHit->GetX();
  y0 = lHit->GetY();
  fTrack->SetFirstPoint(x0,y0,0); //Z-value is set in FitLine

  psi  = (Double_t)atan2(bcent-y0,acent-x0) ;
  psi  = psi + q * AliHLTTransform::PiHalf();
  if ( psi < 0 ) psi = psi + AliHLTTransform::TwoPi();
  pt   = (Double_t)(AliHLTTransform::GetBFieldValue() * radius ) ;
  
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
Int_t AliHLTConfMapFit::FitLine ( )
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
  Double_t radius = (Double_t)(fTrack->GetPt() / AliHLTTransform::GetBFieldValue() ) ;

  //TObjArray *hits = fTrack->GetHits();
  //Int_t num_of_hits = fTrack->GetNumberOfPoints();

  if (0)// fTrack->ComesFromMainVertex() == true ) 
    {
      dx = ((AliHLTConfMapPoint*)fTrack->GetFirstHit())->GetX() - fVertex->GetX();
      dy = ((AliHLTConfMapPoint*)fTrack->GetFirstHit())->GetY() - fVertex->GetY() ;
    }
  else 
    {
      dx = ((AliHLTConfMapPoint *)fTrack->GetFirstHit())->GetX() - ((AliHLTConfMapPoint *)fTrack->GetLastHit())->GetX() ;
      dy = ((AliHLTConfMapPoint *)fTrack->GetFirstHit())->GetY() - ((AliHLTConfMapPoint *)fTrack->GetLastHit())->GetY() ;
      //dx = ((AliHLTConfMapPoint *)hits->First())->GetX() - ((AliHLTConfMapPoint *)hits->Last())->GetX() ;
      //dy = ((AliHLTConfMapPoint *)hits->First())->GetY() - ((AliHLTConfMapPoint *)hits->Last())->GetY() ;
    }
  
  Double_t localPsi = 0.5F * sqrt ( dx*dx + dy*dy ) / radius ;
  Double_t total_s ;
  
  if ( fabs(localPsi) < 1. ) 
    {
      total_s = 2.0 * radius * asin ( localPsi ) ;
    } 
  else 
    { 
      total_s = 2.0 * radius * AliHLTTransform::Pi() ;
    } 
  
  AliHLTConfMapPoint *previousHit = NULL;

  // FtfBaseHit *previousHit = 0  ;
  
  //for ( startLoop() ; done() ; nextHit() ) {
  Double_t dpsi,s;

  //  for(hit_counter=0; hit_counter<num_of_hits; hit_counter++)
  for(fTrack->StartLoop(); fTrack->LoopDone(); fTrack->GetNextHit())  
    {
      // AliHLTConfMapPoint *cHit = (AliHLTConfMapPoint*)hits->At(hit_counter);
      AliHLTConfMapPoint *cHit = (AliHLTConfMapPoint*)fTrack->GetCurrentHit();
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
      //AliHLTConfMapPoint *cHit = (AliHLTConfMapPoint*)hits->At(hit_counter);
      AliHLTConfMapPoint *cHit = (AliHLTConfMapPoint*)fTrack->GetCurrentHit();
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
