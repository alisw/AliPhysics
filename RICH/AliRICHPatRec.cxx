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
  Revision 1.9  2001/02/13 20:38:48  jbarbosa
  Changes to make it work with new IO.

  Revision 1.8  2000/11/01 15:37:18  jbarbosa
  Updated to use its own rec. point object.

  Revision 1.7  2000/10/03 21:44:09  morsch
  Use AliSegmentation and AliHit abstract base classes.

  Revision 1.6  2000/10/02 21:28:12  fca
  Removal of useless dependecies via forward declarations

  Revision 1.5  2000/10/02 15:50:25  jbarbosa
  Fixed forward declarations.

  Revision 1.4  2000/06/30 16:33:43  dibari
  Several changes (ring drawing, fiducial selection, etc.)

  Revision 1.3  2000/06/15 15:47:12  jbarbosa
  Corrected compilation errors on HP-UX (replaced pow with TMath::Power)

  Revision 1.2  2000/06/12 15:26:09  jbarbosa
  Cleaned up version.

  Revision 1.1  2000/06/09 14:53:01  jbarbosa
  Bari's pattern recognition algorithm

*/

#include "AliRICHHit.h"
#include "AliRICHCerenkov.h"
#include "AliRICHSDigit.h"
#include "AliRICHDigit.h"
#include "AliRICHRawCluster.h"
#include "AliRICHRecHit1D.h"
#include "AliRun.h"
#include "AliDetector.h"
#include "AliRICH.h"
#include "AliRICHPoints.h"
#include "AliSegmentation.h"
#include "AliRICHPatRec.h"
#include "AliRICH.h"
#include "AliRICHConst.h"
#include "AliRICHPoints.h"
#include "AliConst.h"
#include "AliHitMap.h"

#include <TParticle.h>
#include <TMath.h>
#include <TRandom.h>
#include <TCanvas.h>
#include <TH2.h>
#include <TTree.h>


ClassImp(AliRICHPatRec)
//___________________________________________
AliRICHPatRec::AliRICHPatRec() : TObject()
{
  // Default constructor
  
    //fChambers = 0;
}
//___________________________________________
AliRICHPatRec::AliRICHPatRec(const char *name, const char *title)
    : TObject()
{
    //Constructor for Bari's pattern recogniton method object
}

void AliRICHPatRec::PatRec()
{

// Pattern recognition algorithm

  AliRICHChamber*       iChamber;
  AliSegmentation*  segmentation;
  	
  Int_t ntracks, ndigits[kNCH];
  Int_t itr, ich, i;
  Int_t goodPhotons;
  Int_t x,y,q;
  Float_t rx,ry,rz;
  Int_t nent,status;
  Int_t padsUsedX[100];
  Int_t padsUsedY[100];

  Float_t rechit[7];

  printf("PatRec started\n");

  AliRICH *pRICH  = (AliRICH*)gAlice->GetDetector("RICH");
  TTree *treeH = gAlice->TreeH();

  ntracks =(Int_t) treeH->GetEntries();
  //  ntracks = 1;
  for (itr=0; itr<ntracks; itr++) {
 
    status = TrackParam(itr,ich);    
    if(status==1) continue;
    //printf(" theta %f phi %f track \n",fTrackTheta,fTrackPhi);
    //    ring->Fill(fTrackLoc[0],fTrackLoc[1],100.);

    iChamber = &(pRICH->Chamber(ich));
    segmentation=iChamber->GetSegmentationModel();

    nent=(Int_t)gAlice->TreeD()->GetEntries();
    gAlice->TreeD()->GetEvent(1);
    TClonesArray *pDigitss = pRICH->DigitsAddress(ich);
    ndigits[ich] = pDigitss->GetEntriesFast();
    printf("Digits in chamber %d: %d\n",ich,ndigits[ich]);
    AliRICHDigit *padI = 0;

    goodPhotons = 0;

    for (Int_t dig=0;dig<ndigits[ich];dig++) {
      padI=(AliRICHDigit*) pDigitss->UncheckedAt(dig);
      x=padI->fPadX;
      y=padI->fPadY;
      q=padI->fSignal;
      segmentation->GetPadC(x,y,rx,ry,rz);      

      //printf("Pad coordinates x:%d, Real coordinates x:%f\n",x,rx);
      //printf("Pad coordinates y:%d, Real coordinates y:%f\n",y,ry);

      fXpad = rx-fXshift;
      fYpad = ry-fYshift;
      fQpad = q;
    
      fCerenkovAnglePad = PhotonCerenkovAngle();
      if(fCerenkovAnglePad==-999) continue;

      if(!PhotonInBand()) continue;

      Int_t xpad;
      Int_t ypad;

      segmentation->GetPadI(fXpad,fYpad,0,xpad,ypad);

      padsUsedX[goodPhotons]=xpad;
      padsUsedY[goodPhotons]=ypad;
      
      goodPhotons++;
      fEtaPhotons[goodPhotons-1] = fCerenkovAnglePad;
    }
    fNumEtaPhotons = goodPhotons;

    BackgroundEstimation();

    HoughResponse();
    //CerenkovRingDrawing();

    rechit[0] = 0;
    rechit[1]   = 0;
    rechit[2] = fThetaCerenkov;
    rechit[3] = fXshift + fTrackLoc[0];
    rechit[4] = fYshift + fTrackLoc[1];
    rechit[5] = fEmissPoint;
    rechit[6] = goodPhotons;

    //printf("Center coordinates:%f %f\n",rechit[3],rechit[4]);
    
    pRICH->AddRecHit1D(ich,rechit,fEtaPhotons,padsUsedX,padsUsedY);
    
  }    

  gAlice->TreeR()->Fill();
  TClonesArray *fRec;
  for (i=0;i<kNCH;i++) {
    fRec=pRICH->RecHitsAddress1D(i);
    int ndig=fRec->GetEntriesFast();
    printf ("Chamber %d, rings %d\n",i,ndig);
  }
  pRICH->ResetRecHits1D();
  
}     


Int_t AliRICHPatRec::TrackParam(Int_t itr, Int_t &ich)
{
  // Get Local coordinates of track impact  

  AliRICHChamber*       iChamber;
  AliSegmentation*      segmentation;

  Float_t trackglob[3];
  Float_t trackloc[3];
  Float_t thetatr;  
  Float_t phitr;
  Float_t iloss;
  Float_t part;
  Float_t pX, pY, pZ;

  //printf("Calling TrackParam\n");

    gAlice->ResetHits();
    TTree *treeH = gAlice->TreeH();
    treeH->GetEvent(itr);
 
    AliRICH *pRICH  = (AliRICH*)gAlice->GetDetector("RICH");
    AliRICHHit* mHit=(AliRICHHit*)pRICH->FirstHit(-1);
    if(mHit==0) return 1;
    ich = mHit->fChamber-1;
    trackglob[0] = mHit->X();
    trackglob[1] = mHit->Y();
    trackglob[2] = mHit->Z();
    pX = mHit->fMomX;
    pY = mHit->fMomY;
    pZ = mHit->fMomZ;
    fTrackMom = sqrt(TMath::Power(pX,2)+TMath::Power(pY,2)+TMath::Power(pZ,2));
    thetatr = (mHit->fTheta)*(Float_t)kDegrad;
    phitr = mHit->fPhi*(Float_t)kDegrad;
    iloss = mHit->fLoss;
    part  = mHit->fParticle;

    iChamber = &(pRICH->Chamber(ich));
    iChamber->GlobaltoLocal(trackglob,trackloc);

    segmentation=iChamber->GetSegmentationModel();

    // retrieve geometrical params

    AliRICHGeometry* fGeometry=iChamber->GetGeometryModel();   
    
    fRw   = fGeometry->GetFreonThickness();
    fQw   = fGeometry->GetQuartzThickness();
    fTgap = fGeometry->GetGapThickness(); 
    Float_t radiatorToPads= fGeometry->GetRadiatorToPads(); 
      //+ fGeometry->GetProximityGapThickness();

    //printf("Distance to pads. From geometry:%f, From calculations:%f\n",radiatorToPads,fRw + fQw + fTgap); 

    //Float_t apar = (fRw + fQw + fTgap)*tan(thetatr);
    Float_t apar = radiatorToPads*tan(thetatr);
    fTrackLoc[0] = apar*cos(phitr);
    fTrackLoc[1] = apar*sin(phitr);
    //fTrackLoc[2] = fRw + fQw + fTgap;
    fTrackLoc[2] = radiatorToPads;
    fTrackTheta  = thetatr;
    fTrackPhi    = phitr;
   
    fXshift = trackloc[0] - fTrackLoc[0];    
    fYshift = trackloc[2] - fTrackLoc[1];
     
    return 0;
}

Float_t AliRICHPatRec::EstimationAtLimits(Float_t lim, Float_t radius,
                                                 Float_t phiphot)
{

// Estimation of emission point

  Float_t nquartz = 1.585;
  Float_t ngas    = 1.;
  Float_t nfreon  = 1.295;
  Float_t value;

  //  printf("Calling EstimationLimits\n");

  Float_t apar = (fRw -fEmissPoint + fQw + fTgap)*tan(fTrackTheta);
  Float_t b1 = (fRw-fEmissPoint)*tan(lim);
  Float_t b2 = fQw / sqrt(TMath::Power(nquartz,2)-TMath::Power(nfreon*sin(lim),2));
  Float_t b3 = fTgap / sqrt(TMath::Power(ngas,2)-TMath::Power(nfreon*sin(lim),2));
  Float_t bpar = b1 + nfreon*sin(lim)*(b2+b3);
  value = TMath::Power(radius,2)
    -TMath::Power((apar*cos(fTrackPhi)-bpar*cos(phiphot)),2)
    -TMath::Power((apar*sin(fTrackPhi)-bpar*sin(phiphot)),2);
  return value;
}


Float_t AliRICHPatRec::PhotonCerenkovAngle()
{
  // Cherenkov pad angle reconstruction
     
  Float_t radius;
  Float_t cherMin = 0;
  Float_t cherMax = 0.8;
  Float_t phiphot;
  Float_t eps = 0.0001;
  Int_t   niterEmiss = 0;
  Int_t   niterEmissMax = 0;
  Float_t x1,x2,x3=0,p1,p2,p3;
  Float_t argY,argX;
  Int_t niterFun;

  //  printf("Calling PhotonCerenkovAngle\n");

  radius = sqrt(TMath::Power(fTrackLoc[0]-fXpad,2)+TMath::Power(fTrackLoc[1]-fYpad,2));
  fEmissPoint = fRw/2.;  //Start value of EmissionPoint
  
  while(niterEmiss<=niterEmissMax) {
 
    niterFun = 0;
    argY = fYpad - fEmissPoint*tan(fTrackTheta)*sin(fTrackPhi);
    argX = fXpad - fEmissPoint*tan(fTrackTheta)*cos(fTrackPhi);
    phiphot = atan2(argY,argX); 
    p1 = EstimationAtLimits(cherMin,radius,phiphot);
    p2 = EstimationAtLimits(cherMax,radius,phiphot);
    if(p1*p2>0)
      {
	//        printf("PhotonCerenkovAngle failed\n");
        return -999;
      }

    //start to find the Cherenkov pad angle 
    x1 = cherMin;
    x2 = cherMax;
    x3 = (x1+x2)/2.;
    p3 = EstimationAtLimits(x3,radius,phiphot);
    while(TMath::Abs(p3)>eps){
      if(p1*p3<0) x2 = x3;
      if(p1*p3>0) {
        x1 = x3;
        p1 = EstimationAtLimits(x1,radius,phiphot);
      }
      x3 = (x1+x2)/2.;
      p3 = EstimationAtLimits(x3,radius,phiphot);
      niterFun++;

      if(niterFun>=1000) {
	//	printf(" max iterations in PhotonCerenkovAngle\n");
	return x3;
      }
    }
    //    printf("niterFun %i \n",niterFun);
    niterEmiss++;
    if (niterEmiss != niterEmissMax+1) EmissionPoint();
  } 
  /* 
   printf(" phiphot %f fXpad %f fYpad %f fEmiss %f \n",
                         phiphot,fXpad,fYpad,fEmissPoint);
  */

  return x3;

}


void AliRICHPatRec::EmissionPoint()
{ 

// Find emission point

  Float_t absorbtionLength=7.83*fRw; //absorption length in the freon (cm)
  // 7.83 = -1/ln(T0) where 
  // T0->Trasmission freon at 180nm = 0.88 (Eph=6.85eV)
  Float_t photonLength, photonLengthMin, photonLengthMax;

  photonLength=exp(-fRw/(absorbtionLength*cos(fCerenkovAnglePad)));
  photonLengthMin=fRw*photonLength/(1.-photonLength);
  photonLengthMax=absorbtionLength*cos(fCerenkovAnglePad);
  fEmissPoint = fRw + photonLengthMin - photonLengthMax;

}

void AliRICHPatRec::PhotonSelection(Int_t track, Int_t &nphot, Float_t &thetamean)
{

// not implemented yet

  printf("Calling PhotonSelection\n");
}

void AliRICHPatRec::BackgroundEstimation()
{

// estimate background noise

  Float_t stepEta   = 0.001;  
  Float_t etaMinBkg = 0.72;
  Float_t etaMaxBkg = 0.75;
  Float_t etaMin    = 0.;
  Float_t etaMax    = 0.75;
  Float_t ngas      = 1.;
  Float_t nfreon    = 1.295;

  Float_t etaStepMin,etaStepMax,etaStepAvg;
  Int_t i,ip,nstep;
  Int_t numPhotBkg, numPhotonStep;
  Float_t funBkg,areaBkg,normBkg;
  Float_t densityBkg,storeBkg,numStore;
  Float_t thetaSig;
  
  numPhotBkg = 0;
  areaBkg = 0.;

  nstep = (int)((etaMaxBkg-etaMinBkg)/stepEta);

  for (i=0;i<fNumEtaPhotons;i++) {

    if(fEtaPhotons[i]>etaMinBkg && fEtaPhotons[i]<etaMaxBkg) {
      numPhotBkg++;
    }    
  }
  if (numPhotBkg == 0) {
     for (i=0;i<fNumEtaPhotons;i++) {
        fWeightPhotons[i] = 1.;
     }
    return;
  }

  //  printf(" numPhotBkg %i ",numPhotBkg);

  for (i=0;i<nstep;i++) {
    etaStepMin = etaMinBkg + (Float_t)(i)*stepEta;
    etaStepMax = etaMinBkg + (Float_t)(i+1)*stepEta;    
    etaStepAvg = 0.5*(etaStepMax + etaStepMin);
    /*
    funBkg = tan(etaStepAvg)*TMath::Power((1.+TMath::Power(tan(etaStepAvg),2)),
				  5.52)-7.803 + 22.02*tan(etaStepAvg);
    */
    
    //printf("etaStepAvg: %f, etaStepMax: %f, etaStepMin: %f", etaStepAvg,etaStepMax,etaStepMin);

    thetaSig = TMath::ASin(nfreon/ngas*TMath::Sin(etaStepAvg));
    funBkg = tan(thetaSig)*(1.+TMath::Power(tan(thetaSig),2))*nfreon
       /ngas*cos(etaStepAvg)/cos(thetaSig);
    areaBkg += stepEta*funBkg;
  }

  densityBkg = 0.95*(Float_t)(numPhotBkg)/areaBkg;
  //  printf(" densityBkg %f \n",densityBkg);
  
  nstep = (int)((etaMax-etaMin)/stepEta); 
  storeBkg = 0.;
  numStore = 0;
  for (i=0;i<nstep;i++) {
    etaStepMin = etaMinBkg + (Float_t)(i)*stepEta;
    etaStepMax = etaMinBkg + (Float_t)(i+1)*stepEta;    
    etaStepAvg = 0.5*(etaStepMax + etaStepMin);
    /*
    funBkg = tan(etaStepAvg)*TMath::Power((1.+TMath::Power(tan(etaStepAvg),2)),
				  5.52)-7.803 + 22.02*tan(etaStepAvg);
    */

    thetaSig = asin(nfreon/ngas*sin(etaStepAvg));
    funBkg = tan(thetaSig)*(1.+TMath::Power(tan(thetaSig),2))*nfreon
       /ngas*cos(etaStepAvg)/cos(thetaSig);

    areaBkg = stepEta*funBkg;
    normBkg = densityBkg*areaBkg; 
    numPhotonStep = 0;
    for (ip=0;ip<fNumEtaPhotons;ip++) {
      if(fEtaPhotons[ip]>etaStepMin && fEtaPhotons[ip]<etaStepMax) {
        numPhotonStep++;
      }
    }
    if (numPhotonStep == 0) {
      storeBkg += normBkg;
      numStore++;
      if (numStore>50) {
        numStore = 0;
        storeBkg = 0.;
      }
    }
    if (numPhotonStep == 0) continue;
    for (ip=0;ip<fNumEtaPhotons;ip++) {
      if(fEtaPhotons[ip]>etaStepMin && fEtaPhotons[ip]<etaStepMax) {
        normBkg +=storeBkg;
        storeBkg = 0;
        numStore = 0;
        fWeightPhotons[ip] = 1. - normBkg/(Float_t)(numPhotonStep);
	/*
        printf(" normBkg %f numPhotonStep %i fW %f \n",
	       normBkg, numPhotonStep, fWeightPhotons[ip]);
	*/
        if(fWeightPhotons[ip]<0) fWeightPhotons[ip] = 0.;
      }
    }
  }
}


void AliRICHPatRec::FlagPhotons(Int_t track, Float_t theta)
{

// not implemented yet

  printf("Calling FlagPhotons\n");
}


//////////////////////////////////////////





Int_t AliRICHPatRec::PhotonInBand()
{ 
  //0=label for parameters giving internal band ellipse
  //1=label for parameters giving external band ellipse  

  Float_t imp[2], mass[2], energy[2], beta[2]; 
  Float_t emissPointLength[2];
  Float_t e1, e2, f1, f2;
  Float_t nfreon[2], nquartz[2]; 
  Int_t times;
  Float_t pointsOnCathode[3];

  Float_t phpad, thetacer[2]; 
  Float_t bandradius[2], padradius;

  imp[0] = 5.0; //threshold momentum for the proton Cherenkov emission
  imp[1] = 1.2;
 
  mass[0] = 0.938; //proton mass 
  mass[1] = 0.139; //pion mass

  emissPointLength[0] = fRw-0.0001; //at the beginning of the radiator
  emissPointLength[1] = 0.;//at the end of radiator
  
  //parameters to calculate freon window refractive index vs. energy
  Float_t a = 1.177;
  Float_t b = 0.0172;
 
  //parameters to calculate quartz window refractive index vs. energy
  /*
  Energ[0]  = 5.6;
  Energ[1]  = 7.7;
  */
  energy[0]  = 5.0;
  energy[1]  = 8.0;
  e1 = 10.666;
  e2  = 18.125;
  f1  = 46.411;
  f2  = 228.71;

  phpad = PhiPad();  

  for (times=0; times<=1; times++) {
  
    nfreon[times]   = a+b*energy[times];
    //nfreon[times]   = 1;

    nquartz[times] = sqrt(1+(f1/(TMath::Power(e1,2)-TMath::Power(energy[times],2)))+
			  (f2/(TMath::Power(e2,2)-TMath::Power(energy[times],2))));

    beta[times]  = imp[times]/sqrt(TMath::Power(imp[times],2)+TMath::Power(mass[times],2));
   
    thetacer[times] =  CherenkovAngle( nfreon[times], beta[times]);

    bandradius[times] = DistanceFromMip( nfreon[times], nquartz[times],
    				emissPointLength[times], 
                                      thetacer[times], phpad, pointsOnCathode);
    //printf(" ppp %f %f %f \n",pointsOnCathode);  
}

  bandradius[0] -= 1.6;
  bandradius[1] += 1.6;
  padradius = sqrt(TMath::Power(fXpad,2)+TMath::Power(fYpad,2));
  //  printf(" rmin %f r %f rmax %f \n",bandradius[0],padradius,bandradius[1]);

  if(padradius>=bandradius[0] && padradius<=bandradius[1]) return 1;
  return 0; 
}

Float_t AliRICHPatRec::DistanceFromMip(Float_t nfreon, Float_t nquartz, 
		       Float_t emissPointLength, Float_t thetacer, 
		       Float_t phpad, Float_t pointsOnCathode[3])
{ 

// Find the distance to MIP impact

  Float_t distanceValue;

  TVector3 radExitPhot(1,1,1);//photon impact at the radiator exit with respect
  //to local reference sistem with the origin in the MIP entrance 
   
  TVector3 vectEmissPointLength(1,1,1);
  Float_t magEmissPointLenght;

  TVector3 radExitPhot2(1,1,1);//photon impact at the radiator exit with respect
  Float_t magRadExitPhot2;
  //to a reference sistem with origin in the photon emission point and  
  //axes parallel to the MIP reference sistem

  TVector3 quarExitPhot(1,1,1);//photon impact at the quartz exit with respect
  Float_t magQuarExitPhot;
  // 
  TVector3 gapExitPhot(1,1,1) ;
  Float_t magGapExitPhot;
  //
  TVector3 PhotocatExitPhot(1,1,1);
  Double_t theta2;
  Double_t thetarad , phirad ;
  Double_t thetaquar, phiquar;
  Double_t thetagap , phigap ;

  Float_t ngas    = 1.;

  magEmissPointLenght =  emissPointLength/cos(fTrackTheta);

  vectEmissPointLength.SetMag(magEmissPointLenght);
  vectEmissPointLength.SetTheta(fTrackTheta);
  vectEmissPointLength.SetPhi(fTrackPhi);


  radExitPhot2.SetTheta(thetacer);  
  radExitPhot2.SetPhi(phpad); 


  TRotation r1;
  TRotation r2;
  TRotation r;

  r1. RotateY(fTrackTheta);
  r2. RotateZ(fTrackPhi);
  


  r = r2 * r1;//rotation about the y axis by MIP theta incidence angle
  //following by a rotation about the z axis by MIP phi incidence angle; 


  radExitPhot2    = r * radExitPhot2;
  theta2          = radExitPhot2.Theta();
  magRadExitPhot2 = (fRw -  vectEmissPointLength(2))/cos(theta2);
  radExitPhot2.SetMag(magRadExitPhot2);


  radExitPhot = vectEmissPointLength + radExitPhot2;
  thetarad    = radExitPhot.Theta();
  phirad  =  radExitPhot.Phi(); //check on the original file //

  thetaquar   = SnellAngle( nfreon, nquartz, theta2); 
  phiquar     = radExitPhot2.Phi(); 
  if(thetaquar == 999.) return thetaquar;
  magQuarExitPhot    = fQw/cos(thetaquar);
  quarExitPhot.SetMag( magQuarExitPhot);
  quarExitPhot.SetTheta(thetaquar);
  quarExitPhot.SetPhi(phiquar);

  thetagap = SnellAngle( nquartz, ngas, thetaquar); 
  phigap   = phiquar; 
  if(thetagap == 999.) return thetagap;
  magGapExitPhot    = fTgap/cos(thetagap);
  gapExitPhot.SetMag( magGapExitPhot);
  gapExitPhot.SetTheta(thetagap);
  gapExitPhot.SetPhi(phigap);

  PhotocatExitPhot =  radExitPhot + quarExitPhot + gapExitPhot; 

  distanceValue = sqrt(TMath::Power(PhotocatExitPhot(0),2)
                           +TMath::Power(PhotocatExitPhot(1),2)); 
  pointsOnCathode[0] = (Float_t) PhotocatExitPhot(0) + fXshift - fTrackLoc[0];
  pointsOnCathode[1] = (Float_t) PhotocatExitPhot(1) + fYshift - fTrackLoc[1];
  pointsOnCathode[2] = (Float_t) PhotocatExitPhot(2);

  //printf(" point in Distance.2. %f %f %f \n",pointsOnCathode[0],pointsOnCathode[1],pointsOnCathode[2]);

  return distanceValue;

 }

Float_t AliRICHPatRec::PhiPad()
{

// ??

  Float_t zpad;
  Float_t thetapad, phipad;
  Float_t thetarot, phirot;

  zpad = fRw + fQw + fTgap;

  TVector3 photonPad(fXpad, fYpad, zpad);
  thetapad = photonPad.Theta();
  phipad = photonPad.Phi();

  TRotation r1;
  TRotation r2;
  TRotation r;

  thetarot = - fTrackTheta;
  phirot   = - fTrackPhi;
  r1. RotateZ(phirot);
  r2. RotateY(thetarot);

  r = r2 * r1;//rotation about the z axis by MIP -phi incidence angle
  //following by a rotation about the y axis by MIP -theta incidence angle; 

  photonPad  = r * photonPad;

  phipad = photonPad.Phi(); 

  return phipad;
}

Float_t AliRICHPatRec:: SnellAngle(Float_t n1, Float_t n2, Float_t theta1)
{ 

// Compute the Snell angle

  Float_t sinrefractangle;
  Float_t refractangle;

  sinrefractangle = (n1/n2)*sin(theta1);

  if(sinrefractangle>1.) {
     refractangle = 999.;
     return refractangle;
  }
  
  refractangle = asin(sinrefractangle);  
  return refractangle;
}

Float_t AliRICHPatRec::CherenkovAngle(Float_t n, Float_t beta)
{ 

// Compute the cerenkov angle

  Float_t thetacer;  
      
  if((n*beta)<1.) {
    thetacer = 999.;
    return thetacer;
  }

  thetacer = acos (1./(n*beta));
  return thetacer;
}

Float_t AliRICHPatRec::BetaCerenkov(Float_t n, Float_t theta)
{ 

// Find beta

  Float_t beta;  
      
  beta = 1./(n*cos(theta));
  return beta;
}




void AliRICHPatRec::HoughResponse()

{	

// Implement Hough response pat. rec. method

  int 		bin=0;
  int           bin1=0;
  int           bin2=0;
  int           i, j, k, nCorrBand;
  int           etaBin = 750;
  float         hcs[750];
  float         angle, thetaCerMean;

  float         etaPeak[30];
  float         etaMin = 0.00;
  float         etaMax = 0.75;
  float         stepEta = 0.001;
  float         windowEta = 0.040;

  int           nBin;

  float etaPeakPos  = -1;
  Int_t   etaPeakCount = -1;
  
  thetaCerMean   = 0.;
  fThetaCerenkov = 0.;    
    
  nBin = (int)(0.5+etaMax/(stepEta));
  nCorrBand = (int)(0.5+ windowEta/(2 * stepEta)); 
  memset ((void *)hcs, 0, etaBin*sizeof(int));

  for (k=0; k< fNumEtaPhotons; k++) {

    angle = fEtaPhotons[k];

    if (angle>=etaMin && angle<= etaMax) {
      bin = (int)(0.5+angle/(stepEta));
      bin1= bin-nCorrBand;
      bin2= bin+nCorrBand;
      if (bin1<0)    bin1=0;
      if (bin2>nBin) bin2=nBin;
      
      for (j=bin1; j<bin2; j++) {
        hcs[j] += fWeightPhotons[k]; 
      }

      thetaCerMean += angle;
    }
  }
 
 thetaCerMean /= fNumEtaPhotons; 
 
  HoughFiltering(hcs);

  for (bin=0; bin <nBin; bin++) {
    angle = (bin+0.5) * (stepEta);
    if (hcs[bin] && hcs[bin] > etaPeakPos) {
      etaPeakCount = 0;
      etaPeakPos = hcs[bin];
      etaPeak[0]=angle;
    }
    else { 
      if (hcs[bin] == etaPeakPos) {
	etaPeak[++etaPeakCount] = angle;
      }
    }
  } 

  for (i=0; i<etaPeakCount+1; i++) {
    fThetaCerenkov += etaPeak[i];
  }
  if (etaPeakCount>=0) {
    fThetaCerenkov /= etaPeakCount+1;
    fThetaPeakPos = etaPeakPos;
  }
}


void AliRICHPatRec::HoughFiltering(float hcs[])
{

// hough filtering

   float hcsFilt[750];
   float k[5] = {0.05, 0.25, 0.4, 0.25, 0.05};
   int nx, i, nxDx;
   int sizeHCS;
   int nBin;

   int   etaBin = 750;
   float etaMax = 0.75;
   float stepEta = 0.001;

   nBin =  (int)(1+etaMax/stepEta); 
   sizeHCS = etaBin*sizeof(float);

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

/*void AliRICHPatRec::CerenkovRingDrawing()
{

//to draw Cherenkov ring by known Cherenkov angle

    Int_t nmaxdegrees;
    Int_t Nphpad;
    Float_t phpad;
    Float_t nfreonave, nquartzave;
    Float_t aveEnerg;
    Float_t energy[2];
    Float_t e1, e2, f1, f2;
    Float_t bandradius;

//parameters to calculate freon window refractive index vs. energy

    Float_t a = 1.177;
    Float_t b = 0.0172;
    
//parameters to calculate quartz window refractive index vs. energy

    energy[0]  = 5.0;
    energy[1]  = 8.0;
    e1  = 10.666;
    e2  = 18.125;
    f1  = 46.411;
    f2  = 228.71;
   

   nmaxdegrees = 36;
   
   for (Nphpad=0; Nphpad<nmaxdegrees;Nphpad++) { 

       phpad = (360./(Float_t)nmaxdegrees)*(Float_t)Nphpad;
      
       aveEnerg =  (energy[0]+energy[1])/2.;
       
       nfreonave  = a+b*aveEnerg;
       nquartzave = sqrt(1+(f1/(TMath::Power(e1,2)-TMath::Power(aveEnerg,2)))+
			 (f2/(TMath::Power(e2,2)-TMath::Power(aveEnerg,2))));
       
       bandradius = DistanceFromMip(nfreonave, nquartzave,
				   fEmissPoint,fThetaCerenkov, phpad);

       fCoordEllipse[0][Nphpad] = fOnCathode[0];
       fCoordEllipse[1][Nphpad] = fOnCathode[1];
       printf(" values %f %f \n",fOnCathode[0],fOnCathode[1]);
       
   }

}*/


