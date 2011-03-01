// @(#) $Id$
// Original: AliHLTFitter.cxx,v 1.14 2005/06/14 10:55:21 cvetan 

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Anders Vestbo, maintained by                          *
 *                  Matthias Richter <Matthias.Richter@ift.uib.no>        *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   AliHLTTPCFitter.cxx
    @author Anders Vestbo, maintained by Matthias Richter
    @date   
    @brief  Fit class HLT for helix
*/

#include <math.h>
#include "AliHLTTPCLogging.h"
#include "AliHLTTPCFitter.h"
#include "AliHLTTPCVertex.h"
#include "AliHLTTPCTrack.h"
#include "AliHLTTPCSpacePointData.h"
#include "AliHLTTPCTransform.h"
//#include "AliHLTTPC.h"

#if __GNUC__ >= 3
using namespace std;
#endif

ClassImp(AliHLTTPCFitter)


AliHLTTPCFitter::AliHLTTPCFitter()
  :
  fTrack(NULL),
  fVertex(NULL),
  fVertexConstraint(0)
{
  //constructor
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  memset(fClusters,0,36*6*sizeof(AliHLTTPCSpacePointData*));
  memset(fNcl,0,36*6*sizeof(UInt_t));
}

AliHLTTPCFitter::AliHLTTPCFitter(AliHLTTPCVertex *vertex,Bool_t vertexconstraint)
  :
  fTrack(NULL),
  fVertex(vertex),
  fVertexConstraint(vertexconstraint)
{
  //constructor
  fTrack=0;
  memset(fClusters,0,36*6*sizeof(AliHLTTPCSpacePointData*));
}

AliHLTTPCFitter::~AliHLTTPCFitter()
{
  //destructor
  for(Int_t i=0; i<36; i++)
    {
      for(Int_t j=0; j<6; j++)
	{
	  if(fClusters[i][j])
	    delete [] fClusters[i][j];
	}
    }
}

void AliHLTTPCFitter::SortTrackClusters(AliHLTTPCTrack *track) const
{
  //Sort the internal cluster list in each track with respect to row numbering.
  //This may be necessary when no conventional track follower has been
  //applied, in which the cluster list has been maintained in a more
  //arbitrary fashion.

  Int_t nhits = track->GetNHits();
  Int_t *ids = (Int_t*)track->GetHitNumbers();
  Int_t *origids = new Int_t[nhits];
  Int_t *mk = new Int_t[nhits];
  Int_t k;

  for(k=0; k<nhits; k++) {origids[k] = ids[k]; mk[k] = -1;}
  
  Int_t slice,patch,id,padrow,maxrow,maxk;
  UInt_t pos;
  for(Int_t j=0; j<nhits; j++)
    {
      maxrow=-1;
      maxk=200;
      for(k=0; k<nhits; k++)
	{
	  id=ids[k];
	  if(id < 0) continue;
	  slice = AliHLTTPCSpacePointData::GetSlice(id);
	  patch = AliHLTTPCSpacePointData::GetPatch(id);
	  pos = AliHLTTPCSpacePointData::GetNumber(id);
	  AliHLTTPCSpacePointData *points = fClusters[slice][patch];
	  padrow = points[pos].fPadRow;
	  if(padrow > maxrow)
	    {
	      maxrow = padrow;
	      maxk=k;
	    }
	}
      mk[j]=maxk;
      ids[maxk]=-1;
    }
    
  for(k=0; k<nhits; k++)
    ids[k] = origids[mk[k]];
  delete [] origids;
  delete [] mk;
}

Int_t AliHLTTPCFitter::FitHelix(AliHLTTPCTrack *track)
{
  //fit helix parameters
  fTrack = track;
  if(FitCircle())
    {
      LOG(AliHLTTPCLog::kError,"AliHLTTPCFitter::FitHelix","TrackFit")<<AliHLTTPCLog::kDec<<
	"Problems during circle fit"<<ENDLOG;
      return 1;
    }
  if(FitLine())
    {
      LOG(AliHLTTPCLog::kError,"AliHLTTPCFitter::FitHelix","TrackFit")<<AliHLTTPCLog::kDec<<
	"Problems during line fit"<<ENDLOG;
      return 1;
    }
  return 0;
}

Int_t AliHLTTPCFitter::FitCircle()
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

  if (!fTrack) return -1;
  
  Double_t wsum  = 0.0 ;
  Double_t xav   = 0.0 ;
  Double_t yav   = 0.0 ;
  
  //
  //     Loop over hits calculating average
  Double_t * fXYWeight = new Double_t[(fTrack->GetNHits())];
  UInt_t *hitnum = fTrack->GetHitNumbers();
  if (!fXYWeight || !hitnum) {
    if (fXYWeight) delete [] fXYWeight;
    return -1;
  }
  memset(fXYWeight, 0, fTrack->GetNHits()*sizeof(Double_t));
  for(Int_t i=0; i<fTrack->GetNHits(); i++)
    {
      UInt_t id = hitnum[i];
      Int_t slice = AliHLTTPCSpacePointData::GetSlice(id);
      Int_t patch = AliHLTTPCSpacePointData::GetPatch(id);
      UInt_t pos = AliHLTTPCSpacePointData::GetNumber(id);
      AliHLTTPCSpacePointData *points = fClusters[slice][patch];
      fXYWeight[i] = 1./ (Double_t)(points[pos].fSigmaY2 + points[pos].fSigmaY2);
      wsum += fXYWeight[i];
      xav += fXYWeight[i]*points[pos].fX;
      yav += fXYWeight[i]*points[pos].fY;
    }
  if (fVertexConstraint == kTRUE)
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
  
  for(Int_t i=0; i<fTrack->GetNHits(); i++)
    { 
      UInt_t id = hitnum[i];
      Int_t slice = AliHLTTPCSpacePointData::GetSlice(id);
      Int_t patch = AliHLTTPCSpacePointData::GetPatch(id);
      UInt_t pos = AliHLTTPCSpacePointData::GetNumber(id);
      AliHLTTPCSpacePointData *points = fClusters[slice][patch];

      xi = points[pos].fX -xav;
      yi        = points[pos].fY - yav ;
      xxav     += xi * xi * fXYWeight[i];
      xyav     += xi * yi * fXYWeight[i];
      yyav     += yi * yi * fXYWeight[i];
    }
  
  if (fVertexConstraint == kTRUE)
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
  
  for(Int_t i=0; i<fTrack->GetNHits(); i++)
    { 
      UInt_t id = hitnum[i];
      Int_t slice = AliHLTTPCSpacePointData::GetSlice(id);
      Int_t patch = AliHLTTPCSpacePointData::GetPatch(id);
      UInt_t pos = AliHLTTPCSpacePointData::GetNumber(id);
      AliHLTTPCSpacePointData *points = fClusters[slice][patch];
      
      xold = points[pos].fX - xav ;
      yold = points[pos].fY - yav ;
      //
      //-->  ROTATE SO THAT <XY> = 0 & DIVIDE BY RSCALE SO THAT <R**2> = 1
      //
      xi = (  cosrot * xold + sinrot * yold ) / rscale ;
      yi = ( -sinrot * xold + cosrot * yold ) / rscale ;
      
      xixi   = xi * xi ;
      yiyi   = yi * yi ;
      riri   = xixi + yiyi ;
      wiriri = fXYWeight[i] * riri ;
      
      xyav   += fXYWeight[i] * xi * yi ;
      xxav   += fXYWeight[i] * xixi ;
      yyav   += fXYWeight[i] * yiyi ;
      
      xrrav  += wiriri * xi ;
      yrrav  += wiriri * yi ;
      rrrrav += wiriri * riri ;
    }
//
//   Include vertex if required
//
  if (fVertexConstraint == kTRUE)
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

  Int_t const kntry = 5 ;
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
  for ( int itry = 1 ; itry <= kntry ; itry++ ) {
     p      = c0 + lamda * (c1 + lamda * (c2 + lamda * lamda * c4 )) ;
     pd     = (c1 + lamda * (c2d + lamda * lamda * c4d)) ;
     dlamda = -p / pd ;
     lamda  = lamda + dlamda ;
     if (fabs(dlamda)<   dlamax) break ;
  }

  //Double_t chi2 = (Double_t)(chiscl * lamda) ;
  //fTrack->SetChiSq1(chi2);
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
    LOG(AliHLTTPCLog::kError,"AliHLTTPCFitter::FitCircle","TrackFit")<<AliHLTTPCLog::kDec<<
      "Problems fitting circle"<<ENDLOG;
    delete [] fXYWeight;
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
  Int_t lastid=fTrack->GetNHits()-1;
  UInt_t id = hitnum[lastid];
  Int_t slice = AliHLTTPCSpacePointData::GetSlice(id);
  Int_t patch = AliHLTTPCSpacePointData::GetPatch(id);
  UInt_t pos = AliHLTTPCSpacePointData::GetNumber(id);
  AliHLTTPCSpacePointData *points = fClusters[slice][patch];
  x0   =  points[pos].fX;
  y0   =  points[pos].fY;
  fTrack->SetFirstPoint(x0,y0,0); //Z-value is set in FitLine
  
  //Set the remaining fit parameters
  psi  = (Double_t)atan2(bcent-y0,acent-x0) ;
  psi  = psi + q * 0.5F * AliHLTTPCTransform::Pi() ;
  if ( psi < 0 ) psi = psi + 2*AliHLTTPCTransform::Pi();
  
  pt   = (Double_t)(AliHLTTPCTransform::GetBFact() * AliHLTTPCTransform::GetBField() * radius ) ;
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
  delete [] fXYWeight;
  return 0 ;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//    Fit Line in s-z plane
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t AliHLTTPCFitter::FitLine ( )
{
  // Fit Line in s-z plane
  if (!fTrack) return -1;

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
  Double_t radius = (Double_t)(fTrack->GetPt() / ( AliHLTTPCTransform::GetBFact() * AliHLTTPCTransform::GetBField() ) ) ;

  Double_t * fS = new Double_t[(fTrack->GetNHits())];
  Double_t *fZWeight = new Double_t[fTrack->GetNHits()];
  UInt_t *hitnum = fTrack->GetHitNumbers();
  if (!fS || !fZWeight || !hitnum) {
    if (fS) delete [] fS;
    if (fZWeight) delete [] fZWeight;
    return -1;
  }
  memset(fS, 0, fTrack->GetNHits()*sizeof(Double_t));
  memset(fZWeight, 0, fTrack->GetNHits()*sizeof(Double_t));
  if (0)//fVertexConstraint==kTRUE)
    {
      UInt_t id = hitnum[0];
      Int_t slice = AliHLTTPCSpacePointData::GetSlice(id);
      Int_t patch = AliHLTTPCSpacePointData::GetPatch(id);
      UInt_t pos = AliHLTTPCSpacePointData::GetNumber(id);
      AliHLTTPCSpacePointData *points = fClusters[slice][patch];
      
      dx = points[pos].fX - fVertex->GetX();
      dy = points[pos].fY - fVertex->GetY();
    }
  else 
    {
      UInt_t id = hitnum[0];
      Int_t slice = AliHLTTPCSpacePointData::GetSlice(id);
      Int_t patch = AliHLTTPCSpacePointData::GetPatch(id);
      UInt_t posf = AliHLTTPCSpacePointData::GetNumber(id);
      AliHLTTPCSpacePointData *pointsf = fClusters[slice][patch];
      id = hitnum[(fTrack->GetNHits()-1)];
      slice = AliHLTTPCSpacePointData::GetSlice(id);
      patch = AliHLTTPCSpacePointData::GetPatch(id);
      UInt_t posl = AliHLTTPCSpacePointData::GetNumber(id);
      AliHLTTPCSpacePointData *pointsl = fClusters[slice][patch];
      dx = pointsf[posf].fX - pointsl[posl].fX;
      dy = pointsf[posf].fY - pointsl[posl].fY;
    }

  Double_t localPsi = 0.5F * sqrt ( dx*dx + dy*dy ) / radius ;
  Double_t totals ;
  
  if ( fabs(localPsi) < 1. ) 
    {
      totals = 2.0 * radius * asin ( localPsi ) ;
    } 
  else 
    { 
      totals = 2.0 * radius * AliHLTTPCTransform::Pi() ;
    } 
  
  Double_t dpsi,s;
  
  for(Int_t i=0; i<fTrack->GetNHits(); i++)
    { 
      UInt_t id = hitnum[i];
      Int_t slice = AliHLTTPCSpacePointData::GetSlice(id);
      Int_t patch = AliHLTTPCSpacePointData::GetPatch(id);
      UInt_t pos = AliHLTTPCSpacePointData::GetNumber(id);
      AliHLTTPCSpacePointData *points = fClusters[slice][patch];
      
      fZWeight[i] = 1./(Double_t)(points[pos].fSigmaZ2);
      if(i>0)
	{
	  id = hitnum[i-1];
	  slice = AliHLTTPCSpacePointData::GetSlice(id);
	  patch = AliHLTTPCSpacePointData::GetPatch(id);
	  UInt_t lastpos = AliHLTTPCSpacePointData::GetNumber(id);
	  AliHLTTPCSpacePointData *lastpoints = fClusters[slice][patch];
	  dx = points[pos].fX -lastpoints[lastpos].fX;
	  dy = points[pos].fY -lastpoints[lastpos].fY;
	  dpsi = 0.5 * (Double_t)sqrt ( dx*dx + dy*dy ) / radius ;
	  if(fabs(dpsi) > 1) {
	    delete [] fS;
	    delete [] fZWeight;
	    return 1;
	  }
	  fTrack->SetPsierr(dpsi);
	  s = fS[i-1] - 2.0 * radius * (Double_t)asin ( dpsi ) ;
	  fS[i]=s;
	}
      else
	fS[i]=totals;
      
      sum += fZWeight[i];
      ss  += fZWeight[i] * fS[i];
      sz  += fZWeight[i] * points[pos].fZ;
      sss += fZWeight[i] * fS[i] * fS[i];
      ssz += fZWeight[i] * fS[i] * points[pos].fZ;
      
    }
  
  
  Double_t chi2,det = sum * sss - ss * ss;
  if ( fabs(det) < 1e-20)
    { 
      chi2 = 99999.F ;
      //fTrack->SetChiSq2(chi2);
      delete [] fS;
      delete [] fZWeight;
      return 0 ;
    }
  
  //Compute the best fitted parameters A,B
  Double_t tanl,z0,dtanl,dz0;

  tanl = (Double_t)((sum * ssz - ss * sz ) / det );
  z0   = (Double_t)((sz * sss - ssz * ss ) / det );

  fTrack->SetTgl(tanl);
  fTrack->SetZ0(z0);
  
  //calculate chi-square 
  chi2 = 0.;
  Double_t r1 ;
  
  for(Int_t i=0; i<fTrack->GetNHits(); i++)
    { 
      UInt_t id = hitnum[i];
      Int_t slice = AliHLTTPCSpacePointData::GetSlice(id);
      Int_t patch = AliHLTTPCSpacePointData::GetPatch(id);
      UInt_t pos = AliHLTTPCSpacePointData::GetNumber(id);
      AliHLTTPCSpacePointData *points = fClusters[slice][patch];
      r1   = points[pos].fZ - tanl * fS[i] - z0 ;
      chi2 += (Double_t) ( (Double_t)(fZWeight[i]) * (r1 * r1) );
    }
  
  //fTrack->SetChiSq2(chi2);
  //
  //calculate estimated variance
  //varsq=chi/(double(n)-2.) 
  //calculate covariance matrix 
  //siga=sqrt(varsq*sxx/det) 
  //sigb=sqrt(varsq*sum/det) 
  //
  dtanl = (Double_t) ( sum / det );
  dz0   = (Double_t) ( sss / det );
  
  fTrack->SetTglerr(dtanl);
  fTrack->SetZ0err(dz0);
  delete [] fZWeight;
  delete [] fS;
  return 0 ;
} 
