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
*/

#include "DataStructures.h"
#include "AliRun.h"
#include "AliDetector.h"
#include "AliRICH.h"
#include "AliRICHPoints.h"
#include "AliRICHSegResV0.h"
#include "AliRICHPatRec.h"
#include "AliRICH.h"
#include "AliRICHConst.h"
#include "AliRICHPoints.h"
#include "AliConst.h"
#include "TParticle.h"
#include "TMath.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TH2.h"


ClassImp(AliRICHPatRec)
//___________________________________________
AliRICHPatRec::AliRICHPatRec() : TObject()
{
    //fChambers = 0;
}
//___________________________________________
AliRICHPatRec::AliRICHPatRec(const char *name, const char *title)
    : TObject()
{
    
}

void AliRICHPatRec::PatRec()
{

  AliRICHChamber*       iChamber;
  AliRICHSegmentation*  segmentation;
  	
  Int_t ntracks, ndigits[7];
  Int_t itr, ich, i;
  Int_t GoodPhotons;
  Int_t x,y,q;
  Float_t rx,ry;
  Int_t nent,status;

  Float_t gamma,MassCer,BetaCer;

  Float_t rechit[5];

  printf("PatRec started\n");

  TCanvas *c1    = new TCanvas("c1","Alice RICH pad hits",50,10,700,700);

  TH2F *ring     = new TH2F("ring","map",90,-30.,30.,90,-30.,30.);
  TH2F *ringband = new TH2F("ringband","map",90,-65.,65.,90,-65.,65.);
  TH1F *cerangle = new TH1F("cerangle","phot",300,0.45,0.75);
  TH1F *ceranglew= new TH1F("ceranglew","phot",300,0.45,0.75);
  TH1F *hough    = new TH1F("hough","hough",75,0.45,0.75);
  TH1F *mass     = new TH1F("mass","mass",100,50.,1050.); 
 
  AliRICH *RICH  = (AliRICH*)gAlice->GetDetector("RICH");
  TTree *TH = gAlice->TreeH();

  ntracks =(Int_t) TH->GetEntries();
  //  ntracks = 1;
  for (itr=0; itr<ntracks; itr++) {
 
    status = TrackParam(itr,ich);    
    if(status==1) continue;
    printf(" theta %f phi %f track \n",fTrackTheta,fTrackPhi);
    //    ring->Fill(fTrackLoc[0],fTrackLoc[1],100.);

    iChamber = &(RICH->Chamber(ich));
    segmentation=iChamber->GetSegmentationModel();

    nent=(Int_t)gAlice->TreeD()->GetEntries();
    gAlice->TreeD()->GetEvent(nent-1);
    TClonesArray *Digits = RICH->DigitsAddress(ich);
    ndigits[ich] = Digits->GetEntriesFast();
    printf("ndigits %d in chamber %d\n",ndigits[ich],ich);
    AliRICHDigit *padI = 0;

    GoodPhotons = 0;

    for (Int_t dig=0;dig<ndigits[ich];dig++) {
      padI=(AliRICHDigit*) Digits->UncheckedAt(dig);
      x=padI->fPadX;
      y=padI->fPadY;
      q=padI->fSignal;
      segmentation->GetPadCxy(x,y,rx,ry);      

      fXpad = rx-fXshift;
      fYpad = ry-fYshift;
      fQpad = q;
    
      ringband->Fill(x,y,1.);
      fCerenkovAnglePad = PhotonCerenkovAngle();
      if(fCerenkovAnglePad==-999) continue;

      if(!PhotonInBand()) continue;

      ring->Fill(fXpad,fYpad,1.);
      cerangle->Fill(fCerenkovAnglePad,1.);

      GoodPhotons++;
      fEtaPhotons[GoodPhotons] = fCerenkovAnglePad;
    }
    fNumEtaPhotons = GoodPhotons;

    BackgroundEstimation();

    for(i=0;i<GoodPhotons;i++) {
      ceranglew->Fill(fEtaPhotons[i],fWeightPhotons[i]);
      //      printf(" Eta %f weight %f \n",fEtaPhotons[i],fWeightPhotons[i]);
    }

    HoughResponse();

    rechit[0] = 0;
    rechit[1]   = 0;
    rechit[2] = fThetaCerenkov;
    rechit[3] = 0;
    rechit[4] = 0;


    hough->Fill(fThetaCerenkov,1.);
    
    RICH->AddRecHit(ich,rechit);
    
    BetaCer = BetaCerenkov(1.29,fThetaCerenkov);
    gamma  = 1./sqrt(1.-pow(BetaCer,2));
    MassCer = fTrackMom/(BetaCer*gamma);
    //    printf(" mass %f \n",MassCer);
    mass->Fill(MassCer*1000,1.);
  }    

  gAlice->TreeR()->Fill();
  TClonesArray *fRec;
  for (i=0;i<7;i++) {
    fRec=RICH->RecHitsAddress(i);
    int ndig=fRec->GetEntriesFast();
    printf ("Chamber %d, rings %d\n",i,ndig);
  }
  RICH->ResetRecHits();
  

  c1->Divide(2,2);
  c1->cd(1);
  ring->Draw("box");
  c1->cd(2);
  ringband->Draw("box");
  c1->cd(3);
  cerangle->Draw();
  c1->cd(4);
  hough->Draw();

}     


Int_t AliRICHPatRec::TrackParam(Int_t itr, Int_t &ich)
{
  // Get Local coordinates of track impact  

  AliRICHChamber*       iChamber;
  AliRICHSegmentation*  segmentation;

  Float_t trackglob[3];
  Float_t trackloc[3];
  Float_t thetatr;  
  Float_t phitr;
  Float_t iloss;
  Float_t part;
  Float_t pX, pY, pZ;

  printf("Calling TrackParam\n");

    gAlice->ResetHits();
    TTree *TH = gAlice->TreeH();
    TH->GetEvent(itr);
 
    AliRICH *RICH  = (AliRICH*)gAlice->GetDetector("RICH");
    AliRICHHit* mHit=(AliRICHHit*)RICH->FirstHit(-1);
    if(mHit==0) return 1;
    ich = mHit->fChamber-1;
    trackglob[0] = mHit->fX;
    trackglob[1] = mHit->fY;
    trackglob[2] = mHit->fZ;
    pX = mHit->fMomX;
    pY = mHit->fMomY;
    pZ = mHit->fMomZ;
    fTrackMom = sqrt(pow(pX,2)+pow(pY,2)+pow(pZ,2));
    thetatr = (180 - mHit->fTheta)*(Float_t)kDegrad;
    phitr = mHit->fPhi*(Float_t)kDegrad;
    iloss = mHit->fLoss;
    part  = mHit->fParticle;

    iChamber = &(RICH->Chamber(ich));
    iChamber->GlobaltoLocal(trackglob,trackloc);

    segmentation=iChamber->GetSegmentationModel();

    // retrieve geometrical params

    AliRICHGeometry* fGeometry=iChamber->GetGeometryModel();   
    
    fRw   = fGeometry->GetFreonThickness();
    fQw   = fGeometry->GetQuartzThickness();
    fTgap = fGeometry->GetGapThickness() 
            + fGeometry->GetProximityGapThickness();  

    Float_t apar = (fRw + fQw + fTgap)*tan(thetatr);  
    fTrackLoc[0] = apar*cos(phitr);
    fTrackLoc[1] = apar*sin(phitr);
    fTrackLoc[2] = fRw + fQw + fTgap;
    fTrackTheta  = thetatr;
    fTrackPhi    = phitr;
   
    fXshift = trackloc[0] - fTrackLoc[0];    
    fYshift = trackloc[2] - fTrackLoc[1];
     
    return 0;
}

Float_t AliRICHPatRec::EstimationAtLimits(Float_t lim, Float_t radius,
                                                 Float_t phiphot)
{
  Float_t nquartz = 1.585;
  Float_t ngas    = 1.;
  Float_t nfreon  = 1.295;
  Float_t value;

  //  printf("Calling EstimationLimits\n");

  Float_t apar = (fRw -fEmissPoint + fQw + fTgap)*tan(fTrackTheta);
  Float_t b1 = (fRw-fEmissPoint)*tan(lim);
  Float_t b2 = fQw / sqrt(pow(nquartz,2)-pow(nfreon*sin(lim),2));
  Float_t b3 = fTgap / sqrt(pow(ngas,2)-pow(nfreon*sin(lim),2));
  Float_t bpar = b1 + nfreon*sin(lim)*(b2+b3);
  value = pow(radius,2)
    -pow((apar*cos(fTrackPhi)-bpar*cos(phiphot)),2)
    -pow((apar*sin(fTrackPhi)-bpar*sin(phiphot)),2);
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
  Float_t x1,x2,x3,p1,p2,p3;
  Float_t argY,argX;
  Int_t niterFun;

  //  printf("Calling PhotonCerenkovAngle\n");

  radius = sqrt(pow(fTrackLoc[0]-fXpad,2)+pow(fTrackLoc[1]-fYpad,2));
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
  Float_t AbsorLength=7.83*fRw; //absorption length in the freon (cm)
  // 7.83 = -1/ln(T0) where 
  // T0->Trasmission freon at 180nm = 0.88 (Eph=6.85eV)
  Float_t PhotonLength, PhotonLengthMin, PhotonLengthMax;

  PhotonLength=exp(-fRw/(AbsorLength*cos(fCerenkovAnglePad)));
  PhotonLengthMin=fRw*PhotonLength/(1.-PhotonLength);
  PhotonLengthMax=AbsorLength*cos(fCerenkovAnglePad);
  fEmissPoint = fRw + PhotonLengthMin - PhotonLengthMax;

}

void AliRICHPatRec::PhotonSelection(Int_t track, Int_t &nphot, Float_t &thetamean)
{
  printf("Calling PhotonSelection\n");
}

void AliRICHPatRec::BackgroundEstimation()
{
  Float_t StepEta   = 0.001;  
  Float_t EtaMinBkg = 0.72;
  Float_t EtaMaxBkg = 0.75;
  Float_t EtaMin    = 0.;
  Float_t EtaMax    = 0.75;
  Float_t ngas      = 1.;
  Float_t nfreon    = 1.295;

  Float_t EtaStepMin,EtaStepMax,EtaStepAvg;
  Int_t i,ip,nstep;
  Int_t NumPhotBkg, NumPhotonStep;
  Float_t FunBkg,AreaBkg,NormBkg;
  Float_t DensityBkg,StoreBkg,NumStore;
  Float_t ThetaSig;
  
  NumPhotBkg = 0;
  AreaBkg = 0.;

  nstep = (int)((EtaMaxBkg-EtaMinBkg)/StepEta);

  for (i=0;i<fNumEtaPhotons;i++) {

    if(fEtaPhotons[i]>EtaMinBkg && fEtaPhotons[i]<EtaMaxBkg) {
      NumPhotBkg++;
    }    
  }
  if (NumPhotBkg == 0) {
     for (i=0;i<fNumEtaPhotons;i++) {
        fWeightPhotons[i] = 1.;
     }
    return;
  }

  //  printf(" NumPhotBkg %i ",NumPhotBkg);

  for (i=0;i<nstep;i++) {
    EtaStepMin = EtaMinBkg + (Float_t)(i)*StepEta;
    EtaStepMax = EtaMinBkg + (Float_t)(i+1)*StepEta;    
    EtaStepAvg = 0.5*(EtaStepMax + EtaStepMin);
    /*
    FunBkg = tan(EtaStepAvg)*pow((1.+pow(tan(EtaStepAvg),2)),
				  5.52)-7.803 + 22.02*tan(EtaStepAvg);
    */
    ThetaSig = asin(nfreon/ngas*sin(EtaStepAvg));
    FunBkg = tan(ThetaSig)*(1.+pow(tan(ThetaSig),2))*nfreon
       /ngas*cos(EtaStepAvg)/cos(ThetaSig);
    AreaBkg += StepEta*FunBkg;
  }

  DensityBkg = 0.95*(Float_t)(NumPhotBkg)/AreaBkg;
  //  printf(" DensityBkg %f \n",DensityBkg);
  
  nstep = (int)((EtaMax-EtaMin)/StepEta); 
  StoreBkg = 0.;
  NumStore = 0;
  for (i=0;i<nstep;i++) {
    EtaStepMin = EtaMinBkg + (Float_t)(i)*StepEta;
    EtaStepMax = EtaMinBkg + (Float_t)(i+1)*StepEta;    
    EtaStepAvg = 0.5*(EtaStepMax + EtaStepMin);
    /*
    FunBkg = tan(EtaStepAvg)*pow((1.+pow(tan(EtaStepAvg),2)),
				  5.52)-7.803 + 22.02*tan(EtaStepAvg);
    */

    ThetaSig = asin(nfreon/ngas*sin(EtaStepAvg));
    FunBkg = tan(ThetaSig)*(1.+pow(tan(ThetaSig),2))*nfreon
       /ngas*cos(EtaStepAvg)/cos(ThetaSig);

    AreaBkg = StepEta*FunBkg;
    NormBkg = DensityBkg*AreaBkg; 
    NumPhotonStep = 0;
    for (ip=0;ip<fNumEtaPhotons;ip++) {
      if(fEtaPhotons[ip]>EtaStepMin && fEtaPhotons[ip]<EtaStepMax) {
        NumPhotonStep++;
      }
    }
    if (NumPhotonStep == 0) {
      StoreBkg += NormBkg;
      NumStore++;
      if (NumStore>50) {
        NumStore = 0;
        StoreBkg = 0.;
      }
    }
    if (NumPhotonStep == 0) continue;
    for (ip=0;ip<fNumEtaPhotons;ip++) {
      if(fEtaPhotons[ip]>EtaStepMin && fEtaPhotons[ip]<EtaStepMax) {
        NormBkg +=StoreBkg;
        StoreBkg = 0;
        NumStore = 0;
        fWeightPhotons[ip] = 1. - NormBkg/(Float_t)(NumPhotonStep);
	/*
        printf(" NormBkg %f NumPhotonStep %i fW %f \n",
	       NormBkg, NumPhotonStep, fWeightPhotons[ip]);
	*/
        if(fWeightPhotons[ip]<0) fWeightPhotons[ip] = 0.;
      }
    }
  }
}


void AliRICHPatRec::FlagPhotons(Int_t track, Float_t theta)
{
  printf("Calling FlagPhotons\n");
}


//////////////////////////////////////////





Int_t AliRICHPatRec::PhotonInBand()
{ 
  //0=label for parameters giving internal band ellipse
  //1=label for parameters giving external band ellipse  

  Float_t imp[2], mass[2], Energ[2], beta[2]; 
  Float_t EmissPointLength[2];
  Float_t E1, E2, F1, F2;
  Float_t nfreon[2], nquartz[2]; 
  Int_t times;


  Float_t phpad, thetacer[2]; 
  Float_t bandradius[2], padradius;

  imp[0] = 5.0; //threshold momentum for the proton Cherenkov emission
  imp[1] = 1.2;
 
  mass[0] = 0.938; //proton mass 
  mass[1] = 0.139; //pion mass

  EmissPointLength[0] = fRw-0.0001; //at the beginning of the radiator
  EmissPointLength[1] = 0.;//at the end of radiator
  
  //parameters to calculate freon window refractive index vs. energy
  Float_t a = 1.177;
  Float_t b = 0.0172;
 
  //parameters to calculate quartz window refractive index vs. energy
  /*
  Energ[0]  = 5.6;
  Energ[1]  = 7.7;
  */
  Energ[0]  = 5.0;
  Energ[1]  = 8.0;
  E1 = 10.666;
  E2  = 18.125;
  F1  = 46.411;
  F2  = 228.71;


  phpad = PhiPad();  

  for (times=0; times<=1; times++) {
  
    nfreon[times]   = a+b*Energ[times];

    nquartz[times] = sqrt(1+(F1/(pow(E1,2)-pow(Energ[times],2)))+
			  (F2/(pow(E2,2)-pow(Energ[times],2))));

    beta[times]  = imp[times]/sqrt(pow(imp[times],2)+pow(mass[times],2));
   
    thetacer[times] =  CherenkovAngle( nfreon[times], beta[times]);

    bandradius[times] = DistanceFromMip( nfreon[times], nquartz[times],
					EmissPointLength[times], 
                                        thetacer[times], phpad);
  }

  bandradius[0] -= 1.6;
  bandradius[1] += 1.6;
  padradius = sqrt(pow(fXpad,2)+pow(fYpad,2));
  //  printf(" rmin %f r %f rmax %f \n",bandradius[0],padradius,bandradius[1]);

  if(padradius>=bandradius[0] && padradius<=bandradius[1]) return 1;
  return 0; 
}

Float_t AliRICHPatRec::DistanceFromMip(Float_t nfreon, Float_t nquartz, 
		       Float_t EmissPointLength, Float_t thetacer, 
		       Float_t phpad)
{ 
  Float_t DistanceValue;

  TVector3 RadExitPhot(1,1,1);//photon impact at the radiator exit with respect
  //to local reference sistem with the origin in the MIP entrance 
   
  TVector3 VectEmissPointLength(1,1,1);
  Float_t MagEmissPointLenght;

  TVector3 RadExitPhot2(1,1,1);//photon impact at the radiator exit with respect
  Float_t MagRadExitPhot2;
  //to a reference sistem with origin in the photon emission point and  
  //axes parallel to the MIP reference sistem

  TVector3 QuarExitPhot(1,1,1);//photon impact at the quartz exit with respect
  Float_t MagQuarExitPhot;
  // 
  TVector3 GapExitPhot(1,1,1) ;
  Float_t MagGapExitPhot;
  //
  TVector3 fPhotocatExitPhot(1,1,1);
  Double_t theta2;
  Double_t thetarad , phirad ;
  Double_t thetaquar, phiquar;
  Double_t thetagap , phigap ;

  Float_t ngas    = 1.;

  MagEmissPointLenght =  EmissPointLength/cos(fTrackTheta);

  VectEmissPointLength.SetMag(MagEmissPointLenght);
  VectEmissPointLength.SetTheta(fTrackTheta);
  VectEmissPointLength.SetPhi(fTrackPhi);


  RadExitPhot2.SetTheta(thetacer);  
  RadExitPhot2.SetPhi(phpad); 


  TRotation r1;
  TRotation r2;
  TRotation r;

  r1. RotateY(fTrackTheta);
  r2. RotateZ(fTrackPhi);
  


  r = r2 * r1;//rotation about the y axis by MIP theta incidence angle
  //following by a rotation about the z axis by MIP phi incidence angle; 


  RadExitPhot2    = r * RadExitPhot2;
  theta2          = RadExitPhot2.Theta();
  MagRadExitPhot2 = (fRw -  VectEmissPointLength(2))/cos(theta2);
  RadExitPhot2.SetMag(MagRadExitPhot2);


  RadExitPhot = VectEmissPointLength + RadExitPhot2;
  thetarad    = RadExitPhot.Theta();

  phirad  =  RadExitPhot.Phi(); //check on the original file //

  thetaquar   = SnellAngle( nfreon, nquartz, theta2); 
  phiquar     = RadExitPhot2.Phi(); 
  if(thetaquar == 999.) return thetaquar;
  MagQuarExitPhot    = fQw/cos(thetaquar);
  QuarExitPhot.SetMag( MagQuarExitPhot);
  QuarExitPhot.SetTheta(thetaquar);
  QuarExitPhot.SetPhi(phiquar);

  thetagap = SnellAngle( nquartz, ngas, thetaquar); 
  phigap   = phiquar; 
  if(thetagap == 999.) return thetagap;
  MagGapExitPhot    = fTgap/cos(thetagap);
  GapExitPhot.SetMag( MagGapExitPhot);
  GapExitPhot.SetTheta(thetagap);
  GapExitPhot.SetPhi(phigap);

  fPhotocatExitPhot =  RadExitPhot + QuarExitPhot + GapExitPhot; 

  DistanceValue = sqrt(pow(fPhotocatExitPhot(0),2)
                           +pow(fPhotocatExitPhot(1),2)); 
  return  DistanceValue ;
}

Float_t AliRICHPatRec::PhiPad()
{
  Float_t zpad;
  Float_t thetapad, phipad;
  Float_t thetarot, phirot;

  zpad = fRw + fQw + fTgap;

  TVector3 PhotonPad(fXpad, fYpad, zpad);
  thetapad = PhotonPad.Theta();
  phipad = PhotonPad.Phi();

  TRotation r1;
  TRotation r2;
  TRotation r;

  thetarot = - fTrackTheta;
  phirot   = - fTrackPhi;
  r1. RotateZ(phirot);
  r2. RotateY(thetarot);

  r = r2 * r1;//rotation about the z axis by MIP -phi incidence angle
  //following by a rotation about the y axis by MIP -theta incidence angle; 

  PhotonPad  = r * PhotonPad;

  phipad = PhotonPad.Phi(); 

  return phipad;
}

Float_t AliRICHPatRec:: SnellAngle(Float_t n1, Float_t n2, Float_t theta1)
{ 
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
  Float_t beta;  
      
  beta = 1./(n*cos(theta));
  return beta;
}




void AliRICHPatRec::HoughResponse()

{	
  int 		bin=0;
  int           bin1=0;
  int           bin2=0;
  int           i, j, k, NcorrBand;
  int           EtaBin = 750;
  float         HCS[750];
  float         angle, ThetaCerMean;

  float         EtaPeak[30];
  float         EtaMin = 0.00;
  float         EtaMax = 0.75;
  float         StepEta = 0.001;
  float         WindowEta = 0.040;

  int           Nbin;

  float EtaPeakPos  = -1;
  Int_t   EtaPeakCount = -1;
  
  ThetaCerMean   = 0.;
  fThetaCerenkov = 0.;    
    
  Nbin = (int)(0.5+EtaMax/(StepEta));
  NcorrBand = (int)(0.5+ WindowEta/(2 * StepEta)); 
  memset ((void *)HCS, 0, EtaBin*sizeof(int));

  for (k=0; k< fNumEtaPhotons; k++) {

    angle = fEtaPhotons[k];

    if (angle>=EtaMin && angle<= EtaMax) {
      bin = (int)(0.5+angle/(StepEta));
      bin1= bin-NcorrBand;
      bin2= bin+NcorrBand;
      if (bin1<0)    bin1=0;
      if (bin2>Nbin) bin2=Nbin;
      
      for (j=bin1; j<bin2; j++) {
        HCS[j] += fWeightPhotons[k]; 
      }

      ThetaCerMean += angle;
    }
  }
 
 ThetaCerMean /= fNumEtaPhotons; 
 
  HoughFiltering(HCS);

  for (bin=0; bin <Nbin; bin++) {
    angle = (bin+0.5) * (StepEta);
    if (HCS[bin] && HCS[bin] > EtaPeakPos) {
      EtaPeakCount = 0;
      EtaPeakPos = HCS[bin];
      EtaPeak[0]=angle;
    }
    else { 
      if (HCS[bin] == EtaPeakPos) {
	EtaPeak[++EtaPeakCount] = angle;
      }
    }
  } 

  for (i=0; i<EtaPeakCount+1; i++) {
    fThetaCerenkov += EtaPeak[i];
  }
  if (EtaPeakCount>=0) {
    fThetaCerenkov /= EtaPeakCount+1;
    fThetaPeakPos = EtaPeakPos;
  }
}


void AliRICHPatRec::HoughFiltering(float HCS[])
{
   float HCS_filt[750];
   float K[5] = {0.05, 0.25, 0.4, 0.25, 0.05};
   int nx, i, nx_dx;
   int sizeHCS;
   int Nbin;

   int   EtaBin = 750;
   float EtaMax = 0.75;
   float StepEta = 0.001;

   Nbin =  (int)(1+EtaMax/StepEta); 
   sizeHCS = EtaBin*sizeof(float);

   memset ((void *)HCS_filt, 0, sizeHCS); 

   for (nx = 0; nx < Nbin; nx++) {
      for (i = 0; i < 5; i++)	{
        nx_dx = nx + (i-2);
	if (nx_dx> -1 && nx_dx<Nbin)
             HCS_filt[nx] +=  HCS[nx_dx] * K[i];
      }      
   }
     
   for (nx = 0; nx < Nbin; nx++) {
     HCS[nx] = HCS_filt[nx];
   }
}

Float_t AliRICHPatRec::CherenkovRingDrawing(Float_t fixedthetacer)

{

//to draw Cherenkov ring by known Cherenkov angle

   Int_t nmaxdegrees, nstepdegrees;
   Float_t phpad, thetacer;
   Float_t nfreonave, nquartzave;
   Float_t AveEnerg;
   Float_t Energ[2];
   Float_t E1, E2, F1, F2;
   Float_t bandradius;
   Float_t CoordPadRing;

//parameters to calculate freon window refractive index vs. energy
   Float_t a = 1.177;
   Float_t b = 0.0172;

//parameters to calculate quartz window refractive index vs. energy
/*
   Energ[0]  = 5.6;
   Energ[1]  = 7.7;
*/	
   Energ[0]  = 5.0;
   Energ[1]  = 8.0;
   E1  = 10.666;
   E2  = 18.125;
   F1  = 46.411;
   F2  = 228.71;


   nmaxdegrees = 360;

   nstepdegrees = 36;

   for (phpad=0; phpad<nmaxdegrees;phpad++) { 
      
     AveEnerg =  (Energ[0]+Energ[1])/2.;

     nfreonave  = a+b*AveEnerg;
     nquartzave = sqrt(1+(F1/(pow(E1,2)-pow(AveEnerg,2)))+
			 (F2/(pow(E2,2)-pow(AveEnerg,2))));

     thetacer =  fixedthetacer;

     bandradius = DistanceFromMip(nfreonave, nquartzave,
				   fEmissPoint,thetacer, phpad); 

     CoordPadRing=fPhotocatExitPhot;

     phpad = (nmaxdegrees/nstepdegrees)*phpad;

     return CoordPadRing;
										    }
 }
