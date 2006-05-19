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
// for single chamber                                                   //
//////////////////////////////////////////////////////////////////////////

#include "AliRICHRecon.h"  //class header
#include <TMath.h>
#include <TRotation.h>
#include <TVector3.h>
#include <TH1F.h>

#include "AliRICHCluster.h" //ThetaCerenkov()
#include "AliRICHParam.h"
#include "AliRICHHelix.h"   //ThetaCerenkov()
#include <AliLog.h>

#define NPointsOfRing 201

//__________________________________________________________________________________________________
AliRICHRecon::AliRICHRecon()
             :TTask       ("RichRec","RichPat")
{
// main ctor
  fThetaMin = 0.0; fThetaMax = 0.75; 
  fDTheta       = 0.001;   fWindowWidth  = 0.045;
  fMinNumPhots = 3;
  fParam=AliRICHParam::Instance(); //get the pointer to AliRICHParam
}
//__________________________________________________________________________________________________
Double_t AliRICHRecon::ThetaCerenkov(AliRICHHelix *pHelix,TClonesArray *pClusters,Int_t &iMipId)
{
// Pattern recognition method based on Hough transform
// Return theta Cerenkov for a given track and list of clusters which are set in ctor  
// Remeber that list of clusters must contain more then 1 cluster. This considiration implies that normally we have 1 mip cluster and few photon clusters per track.  
// Argume
//   Returns: Track theta ckov in rad, nPhot contains number of photon candidates accepted for reconstruction track theta ckov   
  SetTrackTheta(pHelix->Ploc().Theta());  SetTrackPhi(pHelix->Ploc().Phi());
  SetShiftX(pHelix->PosRad().X());        SetShiftY(pHelix->PosRad().Y());
  fClusters = pClusters;
  if(pClusters->GetEntries()>200) fIsWEIGHT = kTRUE; // offset to take into account bkg in reconstruction
  else                            fIsWEIGHT = kFALSE;



  SetThetaCerenkov(-1);   

  //
  // Photon Flag:  Flag = 0 initial set; Flag = 1 good candidate (charge compatible with photon); Flag = 2 photon used for the ring;
  //
  
  for (Int_t iClu=0; iClu<fClusters->GetEntriesFast();iClu++){//clusters loop
    if(iClu == iMipId) continue; // do not consider MIP cluster as a photon candidate
    SetPhotonIndex(iClu);
    SetPhotonFlag(0);
    SetPhotonEta(-999.);
    SetPhotonWeight(0.);
    AliRICHCluster *pClu=(AliRICHCluster*)fClusters->UncheckedAt(iClu);                      //get pointer to current cluster
    if(pClu->Q()>AliRICHParam::QthMIP()) continue;                                           //avoid MIP clusters from bkg
    SetEntranceX(pClu->X() - GetShiftX());    SetEntranceY(pClu->Y() - GetShiftY());         //cluster position with respect to track intersection
    FindPhiPoint();
    SetPhotonFlag(1); 
    FindThetaPhotonCerenkov();
    Float_t thetaPhotonCerenkov = GetThetaPhotonCerenkov();
    AliDebug(1,Form("Track Theta=%5.2f deg, Phi=%5.2f deg Photon clus=%2i ThetaCkov=%5.2f rad",GetTrackTheta()*TMath::RadToDeg(),GetTrackPhi()*TMath::RadToDeg()
                                                                                              ,iClu,thetaPhotonCerenkov ));
    SetPhotonEta(thetaPhotonCerenkov);
  }//clusters loop

  SetPhotonsNumber(fClusters->GetEntries());

  if((iMipId=FlagPhotons(HoughResponse()))<1) return -11; //flag photons according to individual theta ckov with respect to most probable track theta ckov


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
  FindThetaCerenkov();
  return GetThetaCerenkov();
}//ThetaCerenkov()
//__________________________________________________________________________________________________
Int_t AliRICHRecon::PhotonInBand()
{
// Define valid band for photon candidates. For that photons with ThetaMin and ThetaMax are traced up to photcathode

  Float_t nfreon;

  Float_t thetacer;

  Float_t xtoentr = GetEntranceX();
  Float_t ytoentr = GetEntranceY();

  Float_t innerRadius;
  Float_t outerRadius;

  Float_t phpad = GetPhiPoint();


  // inner radius //
  SetEmissionPoint(AliRICHParam::RadThick() -0.0001);

  nfreon = fParam->IdxC6F14(fParam->EckovMin());
  thetacer = 0.;

  AliDebug(1,Form("thetacer in photoninband min %f",thetacer));

  FindThetaAtQuartz(thetacer);

  if(thetacer == 999. || GetThetaAtQuartz() == 999.) {
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
  SetEmissionPoint(0.);
//  SetMassHypotesis(0.139567);

  nfreon = fParam->IdxC6F14(fParam->EckovMax());

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
}//PhotonInBand()
//__________________________________________________________________________________________________
void AliRICHRecon::FindThetaAtQuartz(Float_t thetaCerenkov)
{
// find the theta at the quartz plate

  if(thetaCerenkov == 999.) { SetThetaAtQuartz(999.); return; }

  Float_t thetaAtQuartz = 999.;

  Float_t trackTheta = GetTrackTheta();

  if(trackTheta == 0) {
    thetaAtQuartz = thetaCerenkov;
    SetThetaAtQuartz(thetaAtQuartz);
    return;
  }

  Float_t trackPhi   = GetTrackPhi();
  Float_t phiPoint = GetPhiPoint();

  Double_t den = TMath::Sin((Double_t)trackTheta)*TMath::Cos((Double_t)trackPhi)*TMath::Cos((Double_t)phiPoint) +
    TMath::Sin((Double_t)trackTheta)*TMath::Sin((Double_t)trackPhi)*TMath::Sin((Double_t)phiPoint); 
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

  SetEmissionPoint(AliRICHParam::RadThick()/2);

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



  SetEmissionPoint(AliRICHParam::RadThick()/2.);

  Float_t theta = GetThetaOfRing();
  
  Int_t nPoints = 0;
  Int_t nPsiAccepted = 0;
  Int_t nPsiTotal = 0;

  for(Int_t i=0;i<NPointsOfRing-1;i++){
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

       AliDebug(1,Form("acceptance to detector zone -> %d",zone));	    

      if (zone != 0){
	      FindIntersectionWithDetector();
	      xPoint[nPoints] = GetIntersectionX();	  yPoint[nPoints] = GetIntersectionY();
    	}else{
	      xPoint[nPoints] = xPointRing;	      yPoint[nPoints] = yPointRing;
	      nPsiAccepted++;
	    }
      nPoints++;
  }

  xPoint[nPoints] = xPoint[0];  yPoint[nPoints] = yPoint[0];
  
  // find area...

  Float_t area = 0;

  for (Int_t i = 0; i < nPoints; i++)
    {
      area += TMath::Abs((xPoint[i]-x0)*(yPoint[i+1]-y0) - (xPoint[i+1]-x0)*(yPoint[i]-y0));
    }
  
  area *= 0.5;
  
  Float_t portionOfRing = 0;
  if (nPsiTotal>0) 
    portionOfRing = ((Float_t)nPsiAccepted)/((Float_t)(nPsiTotal));


  SetAreaOfRing(area);
  SetPortionOfRing(portionOfRing);
}//FindAreaAndPortionOfRing()
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
  xIntersect = AliRICHParam::PcSizeX();
  yIntersect = m*(xIntersect - x0) + y0;
  if (yIntersect >= 0 && yIntersect <= AliRICHParam::PcSizeY() && xIntersect >= x1 && xIntersect <= x2)
    {
      SetIntersectionX(xIntersect);
      SetIntersectionY(yIntersect);
      return;
    }
  //
  xIntersect = 0;
  yIntersect = m*(xIntersect - x0) + y0;
  if (yIntersect >= 0 && yIntersect <= AliRICHParam::PcSizeY() && xIntersect >= x1 && xIntersect <= x2)
    {
      SetIntersectionX(xIntersect);
      SetIntersectionY(yIntersect);
      return;
    }
  //
  yIntersect = AliRICHParam::PcSizeY();
  xIntersect = (yIntersect - y0)/m + x0;
  if (xIntersect >= 0 && xIntersect <= AliRICHParam::PcSizeX() && yIntersect >= y1 && yIntersect <= y2)
    {
      SetIntersectionX(xIntersect);
      SetIntersectionY(yIntersect);
      return;
    }
  //
  yIntersect = 0;
  xIntersect = (yIntersect - y0)/m + x0;
  if (xIntersect >= 0 && xIntersect <= AliRICHParam::PcSizeX() && yIntersect >= y1 && yIntersect <= y2)
    {
      SetIntersectionX(xIntersect);
      SetIntersectionY(yIntersect);
      return;
    }
}

//__________________________________________________________________________________________________
Int_t AliRICHRecon::CheckDetectorAcceptance() const
{
  // check for the acceptance

  // crosses X -2.6 2.6 cm
  // crosses Y -1 1 cm

  Float_t xcoord = GetDetectorWhereX();
  Float_t ycoord = GetDetectorWhereY();

  if(xcoord > AliRICHParam::PcSizeX())
    {
      if(ycoord > AliRICHParam::PcSizeY()) return 2;
      if(ycoord > 0 && ycoord < AliRICHParam::PcSizeY()) return 3;
      if(ycoord < 0) return 4;
    }
  if(xcoord < 0)
    {
      if(ycoord > AliRICHParam::PcSizeY()) return 8;
      if(ycoord > 0 && ycoord < AliRICHParam::PcSizeY()) return 7;
      if(ycoord < 0) return 6;
    }
  if(xcoord > 0 && xcoord < AliRICHParam::PcSizeX())
    {
      if(ycoord > AliRICHParam::PcSizeY()) return 1;
      if(ycoord > 0 && ycoord < AliRICHParam::PcSizeY()) return 0;
      if(ycoord < 0) return 5;
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
// Trace current photon from emission point somewhere in radiator to photocathode
// Arguments: none
//   Returns:    

  Float_t nfreon, nquartz, ngas; 

  //fParam->Print();

  nfreon  = fParam->IdxC6F14(fParam->EckovMean());
  nquartz = fParam->IdxSiO2(fParam->EckovMean());
  ngas    = fParam->IdxCH4(fParam->EckovMean());

  Float_t trackTheta = GetTrackTheta();
  Float_t trackPhi = GetTrackPhi();
  Float_t lengthOfEmissionPoint = GetEmissionPoint();

  Float_t theta = GetThetaPhotonInDRS();
  Float_t phi   = GetPhiPhotonInDRS();

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

  Float_t xw = (AliRICHParam::RadThick() - lengthOfEmissionPoint)*cos(phi)*tan(theta);
  Float_t xq = AliRICHParam::WinThick()*cos(phi)*tan(thetaquar);
  Float_t xg = AliRICHParam::Pc2Win()*cos(phi)*tan(thetagap);
  Float_t yw = (AliRICHParam::RadThick() - lengthOfEmissionPoint)*sin(phi)*tan(theta);
  Float_t yq = AliRICHParam::WinThick()*sin(phi)*tan(thetaquar);
  Float_t yg = AliRICHParam::Pc2Win()*sin(phi)*tan(thetagap);


  Float_t xtot = xemiss + xw + xq + xg;
  Float_t ytot = yemiss + yw + yq + yg;

  SetXPointOnCathode(xtot);
  SetYPointOnCathode(ytot);


  Float_t distanceFromEntrance = TMath::Sqrt(TMath::Power(fPhotonLimitX,2)+TMath::Power(fPhotonLimitY,2)); 

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
Double_t AliRICHRecon::HoughResponse()
{
//
//
//       
  Int_t nChannels = (Int_t)(fThetaMax/fDTheta+0.5);
  TH1F *phots   = new TH1F("phots"  ,"phots"  ,nChannels,0,fThetaMax);
  TH1F *photsw  = new TH1F("photsw" ,"photsw" ,nChannels,0,fThetaMax);
  TH1F *resultw = new TH1F("resultw","resultw",nChannels,0,fThetaMax);
  Int_t nBin = (Int_t)(fThetaMax/fDTheta);
  Int_t nCorrBand = (Int_t)(fWindowWidth/(2*fDTheta));
  AliDebug(1,Form("Ring reconstruction for track with theta %f",GetTrackTheta()*TMath::RadToDeg()));	      
  for (Int_t kPhot=0; kPhot< GetPhotonsNumber(); kPhot++){
    SetPhotonIndex(kPhot);
    Double_t angle = GetPhotonEta();
    if(angle<0||angle>fThetaMax) continue;
    phots->Fill(angle);
    Int_t bin = (Int_t)(0.5+angle/(fDTheta));
    Double_t weight=1.;
    if(fIsWEIGHT){
      Double_t lowerlimit = ((Float_t)bin)*fDTheta - 0.5*fDTheta;
      SetThetaOfRing(lowerlimit);
      FindAreaAndPortionOfRing();
      Float_t area1 = GetAreaOfRing();
      Double_t upperlimit = ((Float_t)bin)*fDTheta + 0.5*fDTheta;
      SetThetaOfRing(upperlimit);
      FindAreaAndPortionOfRing();
      Float_t area2 = GetAreaOfRing();
      AliDebug(1,Form("lowerlimit %f  area %f ; upperlimit %f area %f",lowerlimit,area1,upperlimit,area2));	    
      Float_t diffarea = area2 - area1;
      if(diffarea>0){weight = 1./(area2-area1);}else{weight = 1.;}
    }
    AliDebug(1,Form("Calculated weight %f",weight));	    
    photsw->Fill(angle,weight);
    SetPhotonWeight(weight);
  }  
  for (Int_t i=1; i<=nBin;i++){
    Int_t bin1= i-nCorrBand;
    Int_t bin2= i+nCorrBand;
    if(bin1<1) bin1=1;
    if(bin2>nBin)bin2=nBin;
    Double_t sumPhots=phots->Integral(bin1,bin2);
    if(sumPhots<fMinNumPhots) continue; // cut on minimum n. of photons per ring
    Double_t sumPhotsw=photsw->Integral(bin1,bin2);
    resultw->Fill((Float_t)((i+0.5)*fDTheta),sumPhotsw);
  } 
// evaluate the "BEST" theta ckov as the maximum value of histogramm
  Float_t *pVec = resultw->GetArray();
  Int_t locMax = TMath::LocMax(nBin,pVec);
  phots->Delete();photsw->Delete();resultw->Delete(); // Reset and delete objects
  
  return (Double_t)(locMax*fDTheta+0.5*fDTheta); //final most probable track theta ckov   
}//HoughResponse
//__________________________________________________________________________________________________
void AliRICHRecon::FindThetaCerenkov()
{
// Loops on all Ckov candidates and estimates the best Theta Ckov for a ring formed by those candidates. Also estimates an error for that Theat Ckov
// collecting errors for all single Ckov candidates thetas. (Assuming they are independent)  
// Arguments: none
//    Return: none    

  Float_t wei = 0.;
  Float_t weightThetaCerenkov = 0.;

  Double_t etaMin=9999.,etaMax=0.;
  Double_t sigma2 = 0;   //to collect error squared for this ring
  
  for(Int_t i=0;i<GetPhotonsNumber();i++){
    SetPhotonIndex(i);
    if(GetPhotonFlag() == 2){
      Float_t photonEta = GetPhotonEta();
      if(photonEta<etaMin) etaMin=photonEta;
      if(photonEta>etaMax) etaMax=photonEta;
      Float_t photonWeight = GetPhotonWeight();
      weightThetaCerenkov += photonEta*photonWeight;
      wei += photonWeight;      
      //here comes sigma of the reconstructed ring
      
     //Double_t phiref=(GetPhiPoint()-GetTrackPhi());
       if(GetPhotonEta()<=0) continue;//?????????????????Flag photos = 2 may imply CkovEta = 0?????????????? 
                                      //???????????  Look at SetPhoton Flag method    
      Double_t phiref=GetTrackPhi();
  
      Double_t beta = 1./(TMath::Cos(GetPhotonEta())*fParam->IdxC6F14(AliRICHParam::EckovMean()));
      sigma2 += 1./AliRICHParam::SigmaSinglePhotonFormula(GetPhotonEta(),GetPhiPoint(),GetTrackTheta(),phiref,beta);
    }
  }
  
  if(sigma2>0) SetRingSigma2(1./sigma2);
  else         SetRingSigma2(1e10);  
  
  if(wei != 0.) weightThetaCerenkov /= wei; else weightThetaCerenkov = 0.;  
  SetThetaCerenkov(weightThetaCerenkov);

  // estimate of the n. of bkg photons
  SetThetaOfRing(etaMin); FindAreaAndPortionOfRing(); Double_t internalArea = GetAreaOfRing();
  SetThetaOfRing(etaMax); FindAreaAndPortionOfRing(); Double_t externalArea = GetAreaOfRing();

  Double_t effArea = (AliRICHParam::PcSizeX()-AliRICHParam::DeadZone())*(AliRICHParam::PcSizeY()-2*AliRICHParam::DeadZone());
  Double_t nPhotBKG = (externalArea-internalArea)/effArea*fClusters->GetEntries();
  if(nPhotBKG<0) nPhotBKG=0; //just protection from funny angles...
  SetPhotBKG(nPhotBKG);
  
  AliDebug(1,Form(" thetac weighted -> %f",weightThetaCerenkov));
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t AliRICHRecon::FlagPhotons(Double_t thetaCkovHough)
{
// flag photon candidates if their individual theta ckov inside the window  around theta ckov of Hough transform 
// Arguments: thetaCkovHough- value of most probable theta ckov for track as returned by HoughResponse()
//   Returns: number of photon candidates happened to be inside the window

  Int_t steps = (Int_t)((thetaCkovHough - fThetaMin)/ fDTheta); //how many times we need to have fDTheta to fill the distance betwee fThetaMin and thetaCkovHough

  Float_t tmin = fThetaMin + (Float_t)(steps - 1)*fDTheta;
  Float_t tmax = fThetaMin + (Float_t)(steps)*fDTheta;
  Float_t tavg = 0.5*(tmin+tmax);

  tmin = tavg - 0.5*fWindowWidth;  tmax = tavg + 0.5*fWindowWidth;

  Int_t iInsideCnt = 0; //count photons which theta inside prdefined window
  for(Int_t i=0;i<GetPhotonsNumber();i++){//photon candidates loop
    SetPhotonIndex(i);  Float_t photonEta = GetPhotonEta();
    if(photonEta == -999.) continue;
    if(photonEta >= tmin && photonEta <= tmax)	{ 
      SetPhotonFlag(2);	  
      iInsideCnt++;
    }
  }
  return iInsideCnt;
}//FlagPhotons
