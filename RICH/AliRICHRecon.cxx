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

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliRICHRecon                                                         //
//                                                                      //
// RICH class to perfom pattern recognition based on Hough transfrom    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TMath.h>
#include <TRotation.h>
#include <TVector3.h>

#include "AliRICH.h"
#include "AliRICHParam.h"
#include "AliRICHRecon.h"
#include "AliRICHHelix.h"
#include <AliLog.h>

#define NPointsOfRing 201

//__________________________________________________________________________________________________
AliRICHRecon::AliRICHRecon(AliRICHHelix *pHelix,TClonesArray *pClusters,Int_t iMipId)
             :TTask("RichRec","RichPat")
{
// main ctor
  SetFreonScaleFactor(1);
  fIsWEIGHT = kFALSE;
  fThetaBin=750; fThetaMin = 0.0; fThetaMax = 0.75; 
  fDTheta       = 0.001;   fWindowWidth  = 0.045;
  fRadiatorWidth = AliRICHParam::Zfreon();
  fQuartzWidth   = AliRICHParam::Zwin();
  fGapWidth      = AliRICHParam::Freon2Pc() - fRadiatorWidth - fQuartzWidth;
  fXmin         = -AliRICHParam::PcSizeX()/2.;
  fXmax         =  AliRICHParam::PcSizeX()/2.;
  fYmin         = -AliRICHParam::PcSizeY()/2.;
  fYmax         =  AliRICHParam::PcSizeY()/2.; 
  SetTrackTheta(pHelix->Ploc().Theta());
  SetTrackPhi(pHelix->Ploc().Phi());
  SetMipIndex(iMipId);
  SetShiftX(pHelix->PosRad().X());
  SetShiftY(pHelix->PosRad().Y());
  fpClusters = pClusters;
}
//__________________________________________________________________________________________________
Double_t AliRICHRecon::ThetaCerenkov()
{
// Pattern recognition method based on Hough transform
// Return theta Cerenkov for a given track and list of clusters which are set in ctor  

  if(fpClusters->GetEntries()==0) return -1;//no clusters at all for a given track
  Bool_t kPatRec = kFALSE;  
    
  AliDebug(1,Form("---Track Parameters--- Theta: %f , Phi: %f ",GetTrackTheta()*TMath::RadToDeg(),GetTrackPhi()*TMath::RadToDeg()));

  Int_t candidatePhotons = 0;

  SetThetaCerenkov(999.);
  SetHoughPhotons(0);
  SetHoughPhotonsNorm(0);

  for (Int_t j=0; j < fpClusters->GetEntries(); j++){//clusters loop
    SetPhotonIndex(j);
    SetPhotonFlag(0);
    SetPhotonEta(-999.);
    SetPhotonWeight(0.);
    if (j == GetMipIndex()) continue; // do not consider MIP cluster as a candidate photon
    Float_t xtoentr = ((AliRICHCluster*)fpClusters->UncheckedAt(j))->X() - GetShiftX();
    Float_t ytoentr = ((AliRICHCluster*)fpClusters->UncheckedAt(j))->Y() - GetShiftY();
    SetEntranceX(xtoentr);
    SetEntranceY(ytoentr);
    FindPhiPoint();
//      Int_t photonStatus = PhotonInBand();
//      if(photonStatus == 0) continue;
    SetPhotonFlag(1);
    FindThetaPhotonCerenkov();
    Float_t thetaPhotonCerenkov = GetThetaPhotonCerenkov();
    AliDebug(1,Form("THETA CERENKOV ---> %f",thetaPhotonCerenkov));
    SetPhotonEta(thetaPhotonCerenkov);
    candidatePhotons++;
  }//clusters loop

  if(candidatePhotons >= 1) kPatRec = kTRUE;

  if(!kPatRec) return -1;

  SetPhotonsNumber(fpClusters->GetEntries());

  HoughResponse();
  
  fNrings++;

  FlagPhotons();
  Int_t nPhotonHough = GetHoughPhotons();
 
  if(nPhotonHough < 1) 
    {
      SetThetaCerenkov(999.);
      SetHoughPhotonsNorm(0.);
      return -1;
    }

  FindThetaCerenkov();

  AliDebug(1,Form("Number of clusters accepted --->  %i",nPhotonHough));
  
//  Float_t thetaCerenkov = GetThetaCerenkov();  
//  SetThetaOfRing(thetaCerenkov);
//  FindAreaAndPortionOfRing();

//  Float_t nPhotonHoughNorm = ((Float_t)nPhotonHough)/GetPortionOfRing();
//  SetHoughPhotonsNorm(nPhotonHoughNorm);

  // Calculate the area where the photon are accepted...
/*
  Float_t thetaInternal = thetaCerenkov - 0.5*fWindowWidth; 
  SetThetaOfRing(thetaInternal);
  FindAreaAndPortionOfRing();
  Float_t internalArea = GetAreaOfRing();

  Float_t thetaExternal = thetaCerenkov + 0.5*fWindowWidth; 
  SetThetaOfRing(thetaExternal);
  FindAreaAndPortionOfRing();
  Float_t externalArea = GetAreaOfRing();

  Float_t houghArea = externalArea - internalArea;

  SetHoughArea(houghArea);
*/
  return GetThetaCerenkov();

}//ThetaCerenkov()
//__________________________________________________________________________________________________
void AliRICHRecon::FindEmissionPoint()
{
  //estimate the emission point in radiator

// Find emission point

  Float_t absorbtionLenght=7.83*fRadiatorWidth; //absorption length in the freon (cm)
  // 7.83 = -1/ln(T0) where 
  // T0->Trasmission freon at 180nm = 0.88 (Eph=6.85eV)
  Float_t photonLenght, photonLenghtMin, photonLenghtMax;

  photonLenght=exp(-fRadiatorWidth/(absorbtionLenght*cos(fCerenkovAnglePad)));
  photonLenghtMin=fRadiatorWidth*photonLenght/(1.-photonLenght);
  photonLenghtMax=absorbtionLenght*cos(fCerenkovAnglePad);
  Float_t emissionPoint = fRadiatorWidth + photonLenghtMin - photonLenghtMax;

  SetEmissionPoint(emissionPoint);
  SetEmissionPoint(fRadiatorWidth/2); // tune the emission point
}
//__________________________________________________________________________________________________
Int_t AliRICHRecon::PhotonInBand()
{
  //search band fro photon candidates

  //  Float_t massOfParticle;
  Float_t nfreon;

  Float_t thetacer;

  Float_t xtoentr = GetEntranceX();
  Float_t ytoentr = GetEntranceY();

  Float_t innerRadius;
  Float_t outerRadius;

  Float_t phpad = GetPhiPoint();


  // inner radius //
  SetPhotonEnergy(5.6);
  SetEmissionPoint(fRadiatorWidth -0.0001);
  SetFreonRefractiveIndex();

  nfreon = GetFreonRefractiveIndex();
  thetacer = 0.;

  AliDebug(1,Form("thetacer in photoninband min %f",thetacer));

  FindThetaAtQuartz(thetacer);

  if(thetacer == 999. || GetThetaAtQuartz() == 999.)
    {
      innerRadius = -999.;
      SetXInnerRing(-999.);
      SetYInnerRing(-999.);
      SetRadiusInnerRing(-999.);
    }
  else
    {
      SetThetaPhotonInDRS(GetThetaAtQuartz());
      SetPhiPhotonInDRS(phpad);

      innerRadius = FromEmissionToCathode();
       if(innerRadius == 999.) innerRadius = -999.;
      
      SetXInnerRing(GetXPointOnCathode());
      SetYInnerRing(GetYPointOnCathode());
      SetRadiusInnerRing(innerRadius);
    }
  
  // outer radius //
  SetPhotonEnergy(7.7);
  SetEmissionPoint(0.);
//  SetMassHypotesis(0.139567);
  SetFreonRefractiveIndex();

  nfreon = GetFreonRefractiveIndex();

  thetacer = Cerenkovangle(nfreon,1);

  //  thetacer = 0.75;

  AliDebug(1,Form("thetacer in photoninband max %f",thetacer));

  FindThetaAtQuartz(thetacer);

  if(thetacer == 999. || GetThetaAtQuartz() == 999.)
    {
      outerRadius = 999.;
      SetXOuterRing(999.);
      SetYOuterRing(999.);
      SetRadiusOuterRing(999.);
    }
  else
    {
      SetThetaPhotonInDRS(GetThetaAtQuartz());
      SetPhiPhotonInDRS(phpad);

      outerRadius = FromEmissionToCathode();
//      cout << " outerRadius " << outerRadius << endl;
      SetXOuterRing(GetXPointOnCathode());
      SetYOuterRing(GetYPointOnCathode());
      SetRadiusOuterRing(outerRadius);
    }

  Float_t padradius = sqrt(TMath::Power(xtoentr,2)+TMath::Power(ytoentr,2));
  
  AliDebug(1,Form("rmin %f r %f rmax %f",innerRadius,padradius,outerRadius));

  if(padradius>=innerRadius && padradius<=outerRadius) return 1;
  return 0;
}
//__________________________________________________________________________________________________
void AliRICHRecon::FindThetaAtQuartz(Float_t thetaCerenkov)
{
  //find the theta at the quartz plate

  if(thetaCerenkov == 999.) 
    {
      SetThetaAtQuartz(999.);
      return;
    }

  Float_t thetaAtQuartz = 999.;

  Float_t trackTheta = GetTrackTheta();

  if(trackTheta == 0) {

    thetaAtQuartz = thetaCerenkov;
    SetThetaAtQuartz(thetaAtQuartz);
    return;
  }

  Float_t trackPhi   = GetTrackPhi();
  Float_t phiPoint = GetPhiPoint();

  Double_t den = TMath::Sin((Double_t)trackTheta)
    *TMath::Cos((Double_t)trackPhi)
    *TMath::Cos((Double_t)phiPoint) +
    TMath::Sin((Double_t)trackTheta)
    *TMath::Sin((Double_t)trackPhi)
    *TMath::Sin((Double_t)phiPoint); 
  Double_t b = TMath::Cos((Double_t)trackTheta)/den;
  Double_t c = -TMath::Cos((Double_t)thetaCerenkov)/den;

  Double_t underSqrt = 1 + b*b - c*c;

  if(underSqrt < 0) {
    SetThetaAtQuartz(999.);
    return;
  }

  Double_t sol1 = (1+TMath::Sqrt(underSqrt))/(b-c);
  Double_t sol2 = (1-TMath::Sqrt(underSqrt))/(b-c);

  Double_t thetaSol1 = 2*TMath::ATan(sol1);
  Double_t thetaSol2 = 2*TMath::ATan(sol2);

  if(thetaSol1>0 && thetaSol1 < TMath::Pi()) thetaAtQuartz = (Float_t)thetaSol1;
  if(thetaSol2>0 && thetaSol2 < TMath::Pi()) thetaAtQuartz = (Float_t)thetaSol2;

//  AliDebug(1,Form(" Theta @ quartz window %f ",thetaAtQuartz));

  SetThetaAtQuartz(thetaAtQuartz);
}
//__________________________________________________________________________________________________
void AliRICHRecon::FindThetaPhotonCerenkov()
{
  //find theta cerenkov of ring

  Float_t thetaCerMin = 0.;
  Float_t thetaCerMax = 0.75;
  Float_t thetaCerMean;

  Float_t radiusMin, radiusMax, radiusMean;
  Int_t nIteration = 0;

  const Float_t kTollerance = 0.05;


  Float_t phiPoint = GetPhiPoint();

  SetPhotonEnergy(AliRICHParam::MeanCkovEnergy());
  SetEmissionPoint(fRadiatorWidth/2);

  Float_t xPoint = GetEntranceX();
  Float_t yPoint = GetEntranceY();
  Float_t distPoint = TMath::Sqrt(xPoint*xPoint + yPoint*yPoint);

//  AliDebug(1,Form(" DistPoint %f ",distPoint));

  // Star minimization...

  // First value...

  FindThetaAtQuartz(thetaCerMin);
  
  if(GetThetaAtQuartz() == 999.)
    {
      radiusMin = -999.;
    }
  else
    {
      SetThetaPhotonInDRS(GetThetaAtQuartz());
      SetPhiPhotonInDRS(phiPoint);
      
      radiusMin = FromEmissionToCathode();
    }

  // Second value...

  FindThetaAtQuartz(thetaCerMax);
  if(GetThetaAtQuartz() == 999.)
    {
      radiusMax = 999.;
    }
  else
    {
      SetThetaPhotonInDRS(GetThetaAtQuartz());
      SetPhiPhotonInDRS(phiPoint);
      
      radiusMax = FromEmissionToCathode();
    }
  // Mean value...

  thetaCerMean = (thetaCerMax + thetaCerMin)/2;

  FindThetaAtQuartz(thetaCerMean);
  if(GetThetaAtQuartz() == 999.)
    {
      radiusMean = 999.;
    }
  else
    {
      SetThetaPhotonInDRS(GetThetaAtQuartz());
      SetPhiPhotonInDRS(phiPoint);
      
      radiusMean = FromEmissionToCathode();
    }

//  AliDebug(1,Form(" r1 %f rmean %f r2 %f",radiusMin,radiusMean,radiusMax));

  while (TMath::Abs(radiusMean-distPoint) > kTollerance)
    {

      if((radiusMin-distPoint)*(radiusMean-distPoint) < 0) thetaCerMax = thetaCerMean;
      if((radiusMin-distPoint)*(radiusMean-distPoint) > 0) {

	thetaCerMin = thetaCerMean;

	FindThetaAtQuartz(thetaCerMin);
	SetThetaPhotonInDRS(GetThetaAtQuartz());
	SetPhiPhotonInDRS(phiPoint);

	radiusMin =FromEmissionToCathode();
      }

      thetaCerMean = (thetaCerMax + thetaCerMin)/2;

      FindThetaAtQuartz(thetaCerMean);
      SetThetaPhotonInDRS(GetThetaAtQuartz());
      SetPhiPhotonInDRS(phiPoint);

      radiusMean = FromEmissionToCathode();

      nIteration++;
      if(nIteration>=50) {
//	AliDebug(1,Form(" max iterations in FindPhotonCerenkov ",nIteration));
	SetThetaPhotonCerenkov(999.);
	return;
      }
    }

//  AliDebug(1,Form(" distpoint %f radius %f ",distPoint,radiusMean));
  SetThetaPhotonCerenkov(thetaCerMean);

}
//__________________________________________________________________________________________________
void AliRICHRecon::FindAreaAndPortionOfRing()
{
  //find fraction of the ring accepted by the RICH

  Float_t xPoint[NPointsOfRing], yPoint[NPointsOfRing];

  //  Float_t xtoentr = GetEntranceX();
  //  Float_t ytoentr = GetEntranceY();
  Float_t shiftX = GetShiftX();
  Float_t shiftY = GetShiftY();

  Float_t xemiss = GetXCoordOfEmission(); 
  Float_t yemiss = GetYCoordOfEmission(); 

  Float_t x0 = xemiss + shiftX;
  Float_t y0 = yemiss + shiftY;


  SetPhotonEnergy(AliRICHParam::MeanCkovEnergy());
  SetFreonRefractiveIndex();

  SetEmissionPoint(fRadiatorWidth/2.);

  Float_t theta = GetThetaOfRing();
  
  Int_t nPoints = 0;
  Int_t nPsiAccepted = 0;
  Int_t nPsiTotal = 0;

  for(Int_t i=0;i<NPointsOfRing-1;i++)
    {

      Float_t psi = 2*TMath::Pi()*i/NPointsOfRing;
      
      SetThetaPhotonInTRS(theta);
      SetPhiPhotonInTRS(psi);
      FindPhotonAnglesInDRS();
      
      Float_t radius = FromEmissionToCathode();
      if (radius == 999.) continue;
      
      nPsiTotal++;

      Float_t xPointRing = GetXPointOnCathode() + shiftX;
      Float_t yPointRing = GetYPointOnCathode() + shiftY;
      
      SetDetectorWhereX(xPointRing);
      SetDetectorWhereY(yPointRing);
      
      Int_t zone = CheckDetectorAcceptance();


      if (zone != 0) 
	{
	  FindIntersectionWithDetector();
	  xPoint[nPoints] = GetIntersectionX();
	  yPoint[nPoints] = GetIntersectionY();
	}
      else
	{
	  xPoint[nPoints] = xPointRing;
	  yPoint[nPoints] = yPointRing;
	  nPsiAccepted++;
	}

      nPoints++;

    }

  xPoint[nPoints] = xPoint[0];
  yPoint[nPoints] = yPoint[0];
  
  // find area...

  Float_t area = 0;

  for (Int_t i = 0; i < nPoints; i++)
    {
      area += TMath::Abs((xPoint[i]-x0)*(yPoint[i+1]-y0) - (xPoint[i+1]-x0)*(yPoint[i]-y0));
    }
  
  area *= 0.5;
  
  Float_t portionOfRing = ((Float_t)nPsiAccepted)/((Float_t)(nPsiTotal));


  SetAreaOfRing(area);
  SetPortionOfRing(portionOfRing);
}
//__________________________________________________________________________________________________
void AliRICHRecon::FindIntersectionWithDetector()
{
  // find ring intersection with CsI edges

  Float_t xIntersect, yIntersect;
  Float_t x1, x2, y1, y2;

  Float_t shiftX = GetShiftX();
  Float_t shiftY = GetShiftY();

  Float_t xPoint = GetXPointOnCathode() + shiftX;
  Float_t yPoint = GetYPointOnCathode() + shiftY;

  Float_t xemiss = GetXCoordOfEmission(); 
  Float_t yemiss = GetYCoordOfEmission(); 

  Float_t phi = GetPhiPhotonInDRS();
  Float_t m = tan(phi);

  Float_t x0 = xemiss + shiftX;
  Float_t y0 = yemiss + shiftY;

  if(xPoint > x0)
    {
      x1 = x0;
      x2 = xPoint;
    }
  else
    {
      x2 = x0;
      x1 = xPoint;
    }
  if(yPoint > y0)
    {
      y1 = y0;
      y2 = yPoint;
    }
  else
    {
      y2 = y0;
      y1 = yPoint;
    }
  //
  xIntersect = fXmax;
  yIntersect = m*(xIntersect - x0) + y0;
  if (yIntersect >= fYmin && yIntersect <= fYmax && xIntersect >= x1 && xIntersect <= x2)
    {
      SetIntersectionX(xIntersect);
      SetIntersectionY(yIntersect);
      return;
    }
  //
  xIntersect = fXmin;
  yIntersect = m*(xIntersect - x0) + y0;
  if (yIntersect >= fYmin && yIntersect <= fYmax && xIntersect >= x1 && xIntersect <= x2)
    {
      SetIntersectionX(xIntersect);
      SetIntersectionY(yIntersect);
      return;
    }
  //
  yIntersect = fYmax;
  xIntersect = (yIntersect - y0)/m + x0;
  if (xIntersect >= fXmin && xIntersect <= fXmax && yIntersect >= y1 && yIntersect <= y2)
    {
      SetIntersectionX(xIntersect);
      SetIntersectionY(yIntersect);
      return;
    }
  //
  yIntersect = fYmin;
  xIntersect = (yIntersect - y0)/m + x0;
  if (xIntersect >= fXmin && xIntersect <= fXmax && yIntersect >= y1 && yIntersect <= y2)
    {
      SetIntersectionX(xIntersect);
      SetIntersectionY(yIntersect);
      return;
    }
  
  cout << " sono fuori!!!!!!" << endl;
  
}
//__________________________________________________________________________________________________
Int_t AliRICHRecon::CheckDetectorAcceptance() const
{
  // check for the acceptance

  // crosses X -2.6 2.6 cm
  // crosses Y -1 1 cm

  Float_t xcoord = GetDetectorWhereX();
  Float_t ycoord = GetDetectorWhereY();

  if(xcoord > fXmax)
    {
      if(ycoord > fYmax) return 2;
      if(ycoord > fYmin && ycoord < fYmax) return 3;
      if(ycoord < fYmin) return 4;
    }
  if(xcoord < fXmin)
    {
      if(ycoord > fYmax) return 8;
      if(ycoord > fYmin && ycoord < fYmax) return 7;
      if(ycoord < fYmin) return 6;
    }
  if(xcoord > fXmin && xcoord < fXmax)
    {
      if(ycoord > fYmax) return 1;
      if(ycoord > fYmin && ycoord < fYmax) return 0;
      if(ycoord < fYmin) return 5;
    }
  return 999;
}
//__________________________________________________________________________________________________
void AliRICHRecon::FindPhotonAnglesInDRS()
{
  // Setup the rotation matrix of the track...

  TRotation mtheta;
  TRotation mphi;
  TRotation minv;
  TRotation mrot;
  
  Float_t trackTheta = GetTrackTheta();
  Float_t trackPhi = GetTrackPhi();

  mtheta.RotateY(trackTheta);
  mphi.RotateZ(trackPhi);
  
  mrot = mphi * mtheta;
  //  minv = mrot.Inverse();

  TVector3 photonInRadiator(1,1,1);

  Float_t thetaCerenkov = GetThetaPhotonInTRS();
  Float_t phiCerenkov   = GetPhiPhotonInTRS();

  photonInRadiator.SetTheta(thetaCerenkov);
  photonInRadiator.SetPhi(phiCerenkov);
  photonInRadiator = mrot * photonInRadiator;
  Float_t theta = photonInRadiator.Theta();
  Float_t phi = photonInRadiator.Phi();
  SetThetaPhotonInDRS(theta);
  SetPhiPhotonInDRS(phi);

}
//__________________________________________________________________________________________________
Float_t AliRICHRecon::FromEmissionToCathode()
{
  // trace from emission point to cathode

  Float_t nfreon, nquartz, ngas; 

  SetFreonRefractiveIndex();
  SetQuartzRefractiveIndex();
  SetGasRefractiveIndex();

  nfreon  = GetFreonRefractiveIndex();
  nquartz = GetQuartzRefractiveIndex();
  ngas    = GetGasRefractiveIndex();

  Float_t trackTheta = GetTrackTheta();
  Float_t trackPhi = GetTrackPhi();
  Float_t lengthOfEmissionPoint = GetEmissionPoint();

  Float_t theta = GetThetaPhotonInDRS();
  Float_t phi   = GetPhiPhotonInDRS();

//   cout << " Theta " << Theta << " Phi " << Phi << endl;

  Float_t xemiss = lengthOfEmissionPoint*tan(trackTheta)*cos(trackPhi);
  Float_t yemiss = lengthOfEmissionPoint*tan(trackTheta)*sin(trackPhi);

  SetXCoordOfEmission(xemiss);
  SetYCoordOfEmission(yemiss);
  
  Float_t thetaquar = SnellAngle(nfreon, nquartz, theta);

  if(thetaquar == 999.) 
    {
      SetXPointOnCathode(999.);
      SetYPointOnCathode(999.);
      return thetaquar;
    }

  Float_t thetagap  = SnellAngle( nquartz, ngas, thetaquar);

  if(thetagap == 999.) 
    {
      SetXPointOnCathode(999.);
      SetYPointOnCathode(999.);
      return thetagap;
    }

  Float_t xw = (fRadiatorWidth - lengthOfEmissionPoint)*cos(phi)*tan(theta);
  Float_t xq = fQuartzWidth*cos(phi)*tan(thetaquar);
  Float_t xg = fGapWidth*cos(phi)*tan(thetagap);
  Float_t yw = (fRadiatorWidth - lengthOfEmissionPoint)*sin(phi)*tan(theta);
  Float_t yq = fQuartzWidth*sin(phi)*tan(thetaquar);
  Float_t yg = fGapWidth*sin(phi)*tan(thetagap);


  Float_t xtot = xemiss + xw + xq + xg;
  Float_t ytot = yemiss + yw + yq + yg;

  SetXPointOnCathode(xtot);
  SetYPointOnCathode(ytot);


  Float_t distanceFromEntrance = sqrt(TMath::Power(fPhotonLimitX,2)
				    +TMath::Power(fPhotonLimitY,2)); 

  return distanceFromEntrance;

}
//__________________________________________________________________________________________________
void AliRICHRecon::FindPhiPoint()
{
  //find phi of generated point 

  Float_t xtoentr = GetEntranceX();
  Float_t ytoentr = GetEntranceY();

  Float_t trackTheta = GetTrackTheta();
  Float_t trackPhi = GetTrackPhi();

  Float_t emissionPoint = GetEmissionPoint();

  Float_t argY = ytoentr - emissionPoint*tan(trackTheta)*sin(trackPhi);
  Float_t argX = xtoentr - emissionPoint*tan(trackTheta)*cos(trackPhi);
  Float_t phi = atan2(argY,argX);

  SetPhiPoint(phi);

}
//__________________________________________________________________________________________________
Float_t AliRICHRecon::Cerenkovangle(Float_t n, Float_t beta)
{
  // cerenkov angle from n and beta

// Compute the cerenkov angle

  Float_t thetacer;

  if((n*beta)<1.) {
    thetacer = 999.;
    //    cout << " warning in Cerenkoangle !!!!!! " << endl;
    return thetacer;
  }

  thetacer = acos (1./(n*beta));
  return thetacer;
}
//__________________________________________________________________________________________________
Float_t AliRICHRecon::SnellAngle(Float_t n1, Float_t n2, Float_t theta1)
{
  // Snell law

// Compute the Snell angle

  Float_t sinrefractangle;
  Float_t refractangle;

  sinrefractangle = (n1/n2)*sin(theta1);

  if(sinrefractangle>1.) {
    //    cout << " PROBLEMS IN SNELL ANGLE !!!!! " << endl;
    refractangle = 999.;
    return refractangle;
  }
  
  refractangle = asin(sinrefractangle);  
  return refractangle;
}
//__________________________________________________________________________________________________
void AliRICHRecon::HoughResponse()
{
  //Hough response

// Implement Hough response pat. rec. method

  Float_t *hCSspace;

  int 		bin=0;
  int           bin1=0;
  int           bin2=0;
  int           i, j, k, nCorrBand;
  float         hcs[750],hcsw[750];
  float         angle, weight;
  float         lowerlimit,upperlimit;

  float         etaPeak[100];

  int           nBin;

  float etaPeakPos  = -1;

  Int_t   etaPeakCount = -1;
  
  Float_t thetaCerenkov = 0.;
    
  nBin = (int)(0.5+fThetaMax/(fDTheta));
  nCorrBand = (int)(0.5+ fWindowWidth/(2 * fDTheta)); 

  memset ((void *)hcs, 0, fThetaBin*sizeof(float));
  memset ((void *)hcsw, 0, fThetaBin*sizeof(float));

  Int_t nPhotons = GetPhotonsNumber();

  Int_t weightFlag = 0;

  for (k=0; k< nPhotons; k++) {

    SetPhotonIndex(k);

    angle = GetPhotonEta();

    if(angle == -999.) continue;

    if (angle>=fThetaMin && angle<= fThetaMax) 

      {

	bin = (int)(0.5+angle/(fDTheta));

	bin1= bin-nCorrBand;
	bin2= bin+nCorrBand;

	// calculate weights

	if(fIsWEIGHT)
	  {
	    lowerlimit = ((Float_t)bin1)*fDTheta + 0.5*fDTheta;
	    SetThetaOfRing(lowerlimit);
	    FindAreaAndPortionOfRing();
	    Float_t area1 = GetAreaOfRing();
	    
	    upperlimit = ((Float_t)bin2)*fDTheta + 0.5*fDTheta;
	    SetThetaOfRing(upperlimit);
	    FindAreaAndPortionOfRing();
	    Float_t area2 = GetAreaOfRing();
	    
	    //	    cout << "lowerlimit" << lowerlimit << "upperlimit " << upperlimit << endl;
            Float_t diffarea = area2 - area1;

            if(diffarea>0)
              {
	        weight = 1./(area2-area1);
              }
            else
              {
                weightFlag = 1;
		weight = 1.;
              }

	    //	    cout <<" low "<< lowerlimit << " up " << upperlimit << 
	    //	      " area1 " << area1 << " area2 " << area2 << " weight " << weight << endl;
	    
	  }
	else
	  {
	    weight = 1.;
	  }

	SetPhotonWeight(weight);
	
	//	cout << "weight..." << weight << endl;


	if (bin1<0)    bin1=0;
	if (bin2>nBin) bin2=nBin;
      
	for (j=bin1; j<bin2; j++) 
	  {
	    hcs[j] += 1; 
	    hcsw[j] += weight;
	  }
      }
  }
  

  if(weightFlag == 0) 
    {
      hCSspace = hcsw;
    }
  else
    {
      hCSspace = hcs;
      //      cout << " probems with weight...normal procedure adopted " << endl;
    }

  HoughFiltering(hCSspace);

  for (bin=0; bin <nBin; bin++) {
    angle = (bin+0.5) * (fDTheta);
    if (hCSspace[bin] && hCSspace[bin] > etaPeakPos) {
      etaPeakCount = 0;
      etaPeakPos = hCSspace[bin];
      etaPeak[0]=angle;
    }
    else { 
      if (hCSspace[bin] == etaPeakPos) {
	etaPeak[++etaPeakCount] = angle;
      }
    }
  } 

  for (i=0; i<etaPeakCount+1; i++) {
    thetaCerenkov += etaPeak[i];
  }
  if (etaPeakCount>=0) {
    thetaCerenkov /= etaPeakCount+1;
    fThetaPeakPos = etaPeakPos;
  }

  SetThetaCerenkov(thetaCerenkov);
}
//__________________________________________________________________________________________________
void AliRICHRecon::HoughFiltering(float hcs[])
{
  // filter for Hough

// hough filtering

   float hcsFilt[750];
   float k[5] = {0.05, 0.25, 0.4, 0.25, 0.05};
   int nx, i, nxDx;
   int sizeHCS;
   int nBin;

   nBin =  (int)(1+fThetaMax/fDTheta); 
   sizeHCS = fThetaBin*sizeof(float);

   memset ((void *)hcsFilt, 0, sizeHCS); 

   for (nx = 0; nx < nBin; nx++) {
      for (i = 0; i < 5; i++)	{
        nxDx = nx + (i-2);
	if (nxDx> -1 && nxDx<nBin)
             hcsFilt[nx] +=  hcs[nxDx] * k[i];
      }      
   }
     
   for (nx = 0; nx < nBin; nx++) {
     hcs[nx] = hcsFilt[nx];
   }
}
//__________________________________________________________________________________________________
void AliRICHRecon::FindThetaCerenkov()
{
  // manage with weight for photons

  Float_t wei = 0.;
  Float_t weightThetaCerenkov = 0.;

  Int_t nPhotons = GetPhotonsNumber();
  for(Int_t i=0;i<nPhotons;i++)
    {
      SetPhotonIndex(i);

      if(GetPhotonFlag() == 2)
	{
	  Float_t photonEta = GetPhotonEta();
	  Float_t photonWeight = GetPhotonWeight();
	  weightThetaCerenkov += photonEta*photonWeight;
	  wei += photonWeight;
	}
    }

  if(wei != 0.) 
    {
      weightThetaCerenkov /= wei;
    }
  else
    {
      weightThetaCerenkov = 0.;
    }
  
  SetThetaCerenkov(weightThetaCerenkov);

  AliDebug(1,Form(" thetac weighted -> %f",weightThetaCerenkov));
}
//__________________________________________________________________________________________________
void AliRICHRecon::FlagPhotons()
{
  // flag photons

  Int_t nPhotonHough = 0;

  Float_t thetaCerenkov = GetThetaCerenkov();
  AliDebug(1,Form(" fThetaCerenkov %f ",thetaCerenkov));

  Float_t thetaDist= thetaCerenkov - fThetaMin;
  Int_t steps = (Int_t)(thetaDist / fDTheta);

  Float_t tmin = fThetaMin + (Float_t)(steps - 1)*fDTheta;
  Float_t tmax = fThetaMin + (Float_t)(steps)*fDTheta;
  Float_t tavg = 0.5*(tmin+tmax);

  tmin = tavg - 0.5*fWindowWidth;
  tmax = tavg + 0.5*fWindowWidth;

  //  Int_t candidatePhotonsNumber = GetCandidatePhotonsNumber();

  Int_t nPhotons = GetPhotonsNumber();

  //  for(Int_t i=0;i<candidatePhotonsNumber;i++)

  for(Int_t i=0;i<nPhotons;i++)
    {
      SetPhotonIndex(i);

      Float_t photonEta = GetPhotonEta();

      if(photonEta == -999.) continue;

      if(photonEta >= tmin && photonEta <= tmax)
	{
	  SetPhotonFlag(2);
	  nPhotonHough++;
	}
    }
  SetHoughPhotons(nPhotonHough);
}

