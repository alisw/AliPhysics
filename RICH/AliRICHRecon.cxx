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


#include "AliRICH.h"
#include "AliRICHRecon.h"
#include "AliRICHParam.h"
#include <AliLoader.h>
#include <AliRun.h>
#include <AliStack.h>
#include <Riostream.h>
#include <TParticle.h>
#include <TH2.h>
#include <TMath.h>
#include <TRandom.h>
#include <TMinuit.h>
#include <TNtuple.h>
#include <TMath.h>
#include <TRotation.h>
#include <TVector3.h>
#include <TCanvas.h>

#define NPointsOfRing 201

TMinuit *gAliRICHminuit ;

void fcnrecon(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
//__________________________________________________________________________________________________
AliRICHRecon::AliRICHRecon(const char*name, const char*title)
             :TTask(name,title)
{
  // main ctor
  fThetaBin=750; fThetaMin = 0.0; fThetaMax = 0.75; 
  fDTheta       = 0.001;   fWindowWidth  = 0.060;       
  fNpadX         = AliRICHParam::NpadsY();
  fNpadY         = AliRICHParam::NpadsX();
  fPadSizeX      = AliRICHParam::PadSizeY();
  fPadSizeY      = AliRICHParam::PadSizeX();
  fRadiatorWidth = AliRICHParam::FreonThickness();
  fQuartzWidth   = AliRICHParam::QuartzThickness();
  fGapWidth      = AliRICHParam::RadiatorToPads();
  fXmin         = -AliRICHParam::PcSizeY()/2.;
  fXmax         =  AliRICHParam::PcSizeY()/2.;
  fYmin         = -AliRICHParam::PcSizeX()/2.;
  fYmax         =  AliRICHParam::PcSizeX()/2.;  
  fRich = (AliRICH*)gAlice->GetDetector("RICH");
  fOutFile=new TFile("Anal.root","RECREATE","My Analysis histos"); 
  if(fIsDISPLAY) fDisplay = new TCanvas("Display","RICH Display",0,0,1200,750);      
  fNtuple=new TNtuple("hn","ntuple",
"Run:Trig:VertZ:Pmod:Pt:Eta:TrackTheta:TrackPhi:TrackThetaFit:TrackPhiFit:Charge::NPhotons:NPhotonsFit:InRing:MassOfParticle:HoughArea:Multiplicity:TPCLastZ");
}
//__________________________________________________________________________________________________
void AliRICHRecon::StartProcessEvent()
{
  //start to process for pattern recognition
  
  Float_t trackThetaStored    = 0;
  Float_t trackPhiStored      = 0;
  Float_t thetaCerenkovStored = 0;
  Int_t houghPhotonsStored    = 0;
  
  SetFreonScaleFactor(0.994);

  if(fIsDISPLAY) 
    {
      DrawEvent(0);
//      Waiting();
    }

    Rich()->GetLoader()->LoadHits();
    Rich()->GetLoader()->LoadRecPoints();
    Rich()->GetLoader()->LoadDigits();
    gAlice->GetRunLoader()->LoadHeader();
    gAlice->GetRunLoader()->LoadKinematics();    
    
    Rich()->GetLoader()->TreeR()->GetEntry(0);

    Float_t clusX[7][500],clusY[7][500];
    Int_t clusQ[7][500],clusMul[7][500];    
    Int_t nClusters[7];
    
    for (Int_t ich=0;ich<7;ich++) {
      nClusters[ich] = Rich()->Clusters(ich+1)->GetEntries();    
      for(Int_t k=0;k<nClusters[ich];k++) {
        AliRICHcluster *pCluster = (AliRICHcluster *)Rich()->Clusters(ich+1)->At(k);
        clusX[ich][k] = pCluster->X();
        clusY[ich][k] = pCluster->Y();
        clusQ[ich][k] = pCluster->Q();
        clusMul[ich][k] = pCluster->Size();
        pCluster->Print();
      }
    }
        
    Int_t nPrimaries = (Int_t)Rich()->GetLoader()->TreeH()->GetEntries();

    cout << " N. primaries " << nPrimaries << endl;
        
    for(Int_t i=0;i<nPrimaries;i++){
      
      Rich()->GetLoader()->TreeH()->GetEntry(i);

      Rich()->Hits()->Print();
      Int_t iPrim = 0;

      AliRICHhit* pHit=0;
      
      for(Int_t j=0;j<Rich()->Hits()->GetEntries();j++) {

        pHit = (AliRICHhit*)Rich()->Hits()->At(j);
        if(pHit->GetTrack() < nPrimaries) break;
        iPrim++;
      }

      cout << " iPrim " << iPrim << " pHit " << pHit << endl;
      
      if (!pHit) return;
      
      pHit->Print();
      
      TParticle *pParticle = gAlice->GetRunLoader()->Stack()->Particle(pHit->GetTrack());
      Float_t pmod     = pParticle->P();
      Float_t pt       = pParticle->Pt();
      Float_t trackEta = pParticle->Eta();
      Int_t q          = (Int_t)TMath::Sign(1.,pParticle->GetPDG()->Charge());        

      pParticle->Print();
      
      cout << " pmod " << pmod << " pt " << pt << " Eta " << trackEta << " charge " << q << endl;
      
      SetTrackMomentum(pmod); 
      SetTrackPt(pt);
      SetTrackEta(trackEta);
      SetTrackCharge(q);

      TVector3 pGlob(pHit->MomFreoX(),pHit->MomFreoY(),pHit->MomFreoZ());
      TVector3 pLocal = Rich()->C(pHit->Chamber())->Glob2Loc(pGlob,1);
      
      Float_t primGlobalX = pHit->X();
      Float_t primGlobalY = pHit->Y();
      Float_t primGlobalZ = pHit->Z();
      TVector3 primGlobal(primGlobalX,primGlobalY,primGlobalZ);
      TVector3 primLocal = Rich()->C(pHit->Chamber())->Glob2Loc(primGlobal);
      
//      Float_t pmodFreo = pLocal.Mag();
      Float_t trackTheta = pLocal.Theta();
      Float_t trackPhi = pLocal.Phi();

//      cout << " trackTheta " << trackTheta << " trackPhi " << trackPhi << endl;
      
      SetTrackTheta(trackTheta);
      SetTrackPhi(trackPhi);
 
      Int_t maxInd = 0;
      Float_t minDist =  999.;

//      cout << " n Clusters " << nClusters[pHit->Chamber()-1] << " for chamber n. " << pHit->Chamber() << endl;
      
      for(Int_t j=0;j<nClusters[pHit->Chamber()-1];j++)
	{
	  Float_t diffx = primLocal.X() - clusX[pHit->Chamber()-1][j];
	  Float_t diffy = primLocal.Y() - clusY[pHit->Chamber()-1][j];

          
          Float_t diff = sqrt(diffx*diffx + diffy*diffy);

	  if(diff < minDist)
	    {
	      minDist = diff;
	      maxInd = j;
	    }

	}

      Float_t diffx = primLocal.X() - clusX[pHit->Chamber()-1][maxInd];
      Float_t diffy = primLocal.Y() - clusY[pHit->Chamber()-1][maxInd];

      cout << " diffx " << diffx << " diffy " << diffy << endl;
      

      SetMipIndex(maxInd);
      SetTrackIndex(i);

      Float_t shiftX = primLocal.X()/primLocal.Z()*(fRadiatorWidth+fQuartzWidth+fGapWidth) + primLocal.X();
      Float_t shiftY = primLocal.Y()/primLocal.Z()*(fRadiatorWidth+fQuartzWidth+fGapWidth) + primLocal.Y();
      
      SetShiftX(shiftX);
      SetShiftY(shiftY);

      Float_t *pclusX = &clusX[pHit->Chamber()-1][0];
      Float_t *pclusY = &clusY[pHit->Chamber()-1][0];
      
      SetCandidatePhotonX(pclusX);
      SetCandidatePhotonY(pclusY);
      SetCandidatePhotonsNumber(nClusters[pHit->Chamber()-1]);

      Int_t qch = clusQ[pHit->Chamber()-1][maxInd];

       
      if(minDist < 3.0 && qch > 120 && maxInd !=0) 
	{
	  
	  if(fIsBACKGROUND)
	    {
	      
	      Float_t xrndm = fXmin + (fXmax-fXmin)*gRandom->Rndm(280964);
	      Float_t yrndm = fYmin + (fYmax-fYmin)*gRandom->Rndm(280964);
	      SetShiftX(xrndm);
	      SetShiftY(yrndm);
	      
	    }

	  PatRec();

	  trackThetaStored = GetTrackTheta();
	  trackPhiStored = GetTrackPhi();
	  thetaCerenkovStored = GetThetaCerenkov();
	  houghPhotonsStored = GetHoughPhotons();
	  
          Int_t diffNPhotons = 999;
          Int_t nsteps = 0;
          Float_t diffTrackTheta = 999.;
          Float_t diffTrackPhi   = 999.;

	  while(fIsMINIMIZER && GetHoughPhotons() > 2 
                            && diffNPhotons !=0 
                            && diffTrackTheta > 0.0001
                            && nsteps < 10)
	    {

	      Int_t   houghPhotonsBefore  = GetHoughPhotons();

	      Float_t trackThetaBefore = GetTrackTheta();
	      Float_t trackPhiBefore   = GetTrackPhi();
	  
	      Minimization(); 

              PatRec();
 
              diffNPhotons = TMath::Abs(houghPhotonsBefore - GetHoughPhotons()); 

	      Float_t trackThetaAfter = GetTrackTheta();
	      Float_t trackPhiAfter   = GetTrackPhi();

              diffTrackTheta = TMath::Abs(trackThetaAfter - trackThetaBefore);
              diffTrackPhi   = TMath::Abs(trackPhiAfter - trackPhiBefore);

              if(fDebug)
              cout << " houghPhotonsBefore " << houghPhotonsBefore
                   << " GetHoughPhotons()  " << GetHoughPhotons();

              nsteps++;
	    }

	  SetFittedThetaCerenkov(GetThetaCerenkov());
	  SetFittedHoughPhotons(GetHoughPhotons());

	  SetTrackTheta(trackThetaStored);
	  SetTrackPhi(trackPhiStored);
	  SetThetaCerenkov(thetaCerenkovStored);
	  SetHoughPhotons(houghPhotonsStored);

          SetMinDist(minDist);

	  FillHistograms();
      
	  if(fIsDISPLAY) DrawEvent(1);

	  Waiting();

	}
    }
  if(fIsDISPLAY) fDisplay->Print("display.ps");
}//StartProcessEvent()
//__________________________________________________________________________________________________
void AliRICHRecon::EndProcessEvent()
{
  // function called at the end of the event loop

  fOutFile->Write();
  fOutFile->Close();                                                     
}
//__________________________________________________________________________________________________
void AliRICHRecon::PatRec()
{
  //pattern recognition method based on Hough transform

  
  Float_t trackTheta = GetTrackTheta();
  Float_t trackPhi   = GetTrackPhi();
  Float_t pmod       = GetTrackMomentum();
  Int_t iMipIndex   = GetMipIndex();

  Bool_t kPatRec = kFALSE;  

  Int_t candidatePhotons = 0;

  Float_t shiftX = GetShiftX();
  Float_t shiftY = GetShiftY();

  Float_t* candidatePhotonX = GetCandidatePhotonX();
  Float_t* candidatePhotonY = GetCandidatePhotonY();

  Int_t candidatePhotonsNumber = GetCandidatePhotonsNumber();

  if(fDebug) cout << " n " << candidatePhotonsNumber << endl;

  SetThetaCerenkov(999.);
  SetHoughPhotons(0);
  SetHoughPhotonsNorm(0);
  SetHoughRMS(999.);

  for (Int_t j=0; j < candidatePhotonsNumber; j++)
    {

      SetPhotonIndex(j);

      SetPhotonFlag(0);
      SetPhotonEta(-999.);
      SetPhotonWeight(0.);

      if (j == iMipIndex) continue;

        
      if(candidatePhotonX[j] < -64.) continue; /* avoid artificial clusters from edge uesd by Yale.... */

      Float_t xtoentr = candidatePhotonX[j] - shiftX;
      Float_t ytoentr = candidatePhotonY[j] - shiftY;

      //      Float_t chargehit = fHits_charge[j]; 
      //      if(chargehit > 150) continue;

      SetEntranceX(xtoentr);
      SetEntranceY(ytoentr);

      FindPhiPoint();

      Int_t photonStatus = PhotonInBand();
 
      if(fDebug)
         {
            cout << " Photon n. " << j << " Status " << photonStatus << " accepted " << endl;
            cout << " CandidatePhotonX[j] " << candidatePhotonX[j] << " CandidatePhotonY[j] " << candidatePhotonY[j] << endl;
         }
    
      if(photonStatus == 0) continue;

      SetPhotonFlag(1);

      FindThetaPhotonCerenkov();

      Float_t thetaPhotonCerenkov = GetThetaPhotonCerenkov();

      if(fDebug) cout << " theta photon " << thetaPhotonCerenkov << endl;

      SetPhotonEta(thetaPhotonCerenkov);

      candidatePhotons++;

      
    }

  if(candidatePhotons >= 1) kPatRec = kTRUE;

  if(!kPatRec) return;
    {
       SetThetaCerenkov(999.);
       SetHoughPhotons(0);
    }
  SetPhotonsNumber(candidatePhotonsNumber);

  HoughResponse();
  
  fNrings++;

  FlagPhotons();
  Int_t nPhotonHough = GetHoughPhotons();
 
  if(nPhotonHough < 1) 
    {
      SetThetaCerenkov(999.);
      SetHoughPhotonsNorm(0.);
      return;
    }

  if(fIsWEIGHT) FindWeightThetaCerenkov();

  Float_t thetaCerenkov = GetThetaCerenkov();

  SetThetaOfRing(thetaCerenkov);
  FindAreaAndPortionOfRing();

  Float_t nPhotonHoughNorm = ((Float_t)nPhotonHough)/GetPortionOfRing();
  SetHoughPhotonsNorm(nPhotonHoughNorm);

  // Calculate the area where the photon are accepted...

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

  if(fDebug)
    {
      cout << " ----- SUMMARY OF RECONSTRUCTION ----- " << endl; 
      cout << " Rings found " << fNrings << " with thetac " << thetaCerenkov << endl;
      
      
      cout << " Nphotons " << GetPhotonsNumber() 
	   << " Hough    " << nPhotonHough 
	   << " norm     " << nPhotonHoughNorm << endl;
      
      cout << " In PatRec:p " << pmod << " theta " << trackTheta << " phi " << trackPhi << endl;
      cout << " ------------------------------------- " << endl; 
    }

  Int_t nPhotons = GetPhotonsNumber();

  Float_t xmean = 0.;
  Float_t x2mean = 0.;
  Int_t nev = 0;

  for (Int_t j=0; j < nPhotons;j++)
    {
      SetPhotonIndex(j);

      Float_t eta = GetPhotonEta();

      if(eta != -999.) 
	{
	  if(GetPhotonFlag() == 2) 
	    {


	      xmean += eta;
	      x2mean += eta*eta;
	      nev++;
	    }
	}
    }

  if(nev > 0)
    {
      xmean /=(Float_t)nev;
      x2mean /=(Float_t)nev;
    } else {
      xmean = 0.;
      x2mean = 0.;
    }

  Float_t vRMS = sqrt(x2mean - xmean*xmean);

  SetHoughRMS(vRMS);

  if(fDebug) cout << " RMS " << vRMS << endl;

}

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
}


Int_t AliRICHRecon::PhotonInBand()
{
  //search band fro photon candidates

  //  Float_t massOfParticle;
  Float_t beta;
  Float_t nfreon;

  Float_t thetacer;

  Float_t xtoentr = GetEntranceX();
  Float_t ytoentr = GetEntranceY();

  Float_t innerRadius;
  Float_t outerRadius;

  Float_t phpad = GetPhiPoint();

  //  Float_t pmod = GetTrackMomentum();
  //  Float_t trackTheta = GetTrackTheta();
  //  Float_t trackPhi = GetTrackPhi();

  // inner radius //
  SetPhotonEnergy(5.6);
  SetEmissionPoint(fRadiatorWidth -0.0001);
  SetMassHypotesis(0.93828);

  SetBetaOfParticle();
  SetFreonRefractiveIndex();

  beta   = GetBetaOfParticle();
  nfreon = GetFreonRefractiveIndex();

  thetacer = Cerenkovangle(nfreon,beta);

  thetacer = 0.;

  if(fDebug) cout << " thetacer in photoninband min " << thetacer << endl;

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
  SetMassHypotesis(0.);

  SetBetaOfParticle();
  SetFreonRefractiveIndex();

  beta   = GetBetaOfParticle();
  nfreon = GetFreonRefractiveIndex();

  thetacer = Cerenkovangle(nfreon,beta);

  //  thetacer = 0.75;

  if(fDebug) cout << " thetacer in photoninband max " << thetacer << endl;

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
  
  if(fDebug) printf(" rmin %f r %f rmax %f \n",innerRadius,padradius,outerRadius);

  if(padradius>=innerRadius && padradius<=outerRadius) return 1;
  return 0;
}

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

    if(fDebug) cout << " Theta sol unique " << thetaCerenkov << endl;  

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

  if(fDebug)
    {
      cout << " trackTheta    " << trackTheta    << endl;
      cout << " TrackPhi      " << trackPhi      << endl;
      cout << " PhiPoint      " << phiPoint      << endl;
      cout << " ThetaCerenkov " << thetaCerenkov << endl;
      cout << " den b c " << den << " b " << b << " c " << c << endl;
    }

  if(underSqrt < 0) {
    if(fDebug) cout << " sqrt negative !!!!" << underSqrt << endl;
    SetThetaAtQuartz(999.);
    return;
  }

  Double_t sol1 = (1+TMath::Sqrt(underSqrt))/(b-c);
  Double_t sol2 = (1-TMath::Sqrt(underSqrt))/(b-c);

  Double_t thetaSol1 = 2*TMath::ATan(sol1);
  Double_t thetaSol2 = 2*TMath::ATan(sol2);

  if(fDebug) cout << " Theta sol 1 " << thetaSol1 
		  << " Theta sol 2 " << thetaSol2 << endl;  

  if(thetaSol1>0 && thetaSol1 < TMath::Pi()) thetaAtQuartz = (Float_t)thetaSol1;
  if(thetaSol2>0 && thetaSol2 < TMath::Pi()) thetaAtQuartz = (Float_t)thetaSol2;

  SetThetaAtQuartz(thetaAtQuartz);
}

void AliRICHRecon::FindThetaPhotonCerenkov()
{
  //find theta cerenkov of ring

  Float_t thetaCerMin = 0.;
  Float_t thetaCerMax = 0.75;
  Float_t thetaCerMean;

  Float_t radiusMin, radiusMax, radiusMean;
  Int_t nIteration = 0;

  const Float_t kTollerance = 0.05;

  //  Float_t pmod = GetTrackMomentum();
  //  Float_t trackTheta = GetTrackTheta();
  //  Float_t trackPhi = GetTrackPhi();

  Float_t phiPoint = GetPhiPoint();

  SetPhotonEnergy(6.85);
  SetEmissionPoint(fRadiatorWidth/2);

  Float_t xPoint = GetEntranceX();
  Float_t yPoint = GetEntranceY();
  Float_t distPoint = sqrt(xPoint*xPoint + yPoint*yPoint);

  if(fDebug) cout << " DistPoint " << distPoint << endl;

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

  if(fDebug) cout << " r1 " << radiusMin << " rmean " 
		  << radiusMean << " r2 " << radiusMax << endl;

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
	if(fDebug) printf(" max iterations in FindPhotonCerenkov\n");
	SetThetaPhotonCerenkov(999.);
	return;
      }
    }

  SetThetaPhotonCerenkov(thetaCerMean);

}

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

  //  Float_t pmod = GetTrackMomentum();
  //  Float_t trackTheta = GetTrackTheta();
  //  Float_t trackPhi = GetTrackPhi();

  SetPhotonEnergy(6.85);
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
Float_t AliRICHRecon::PhotonPositionOnCathode()
{ 
  // find the photon position on the CsI
  //  Float_t massOfParticle;
  Float_t beta;
  Float_t nfreon;


  SetPhotonEnergy(6.85);
  SetEmissionPoint(fRadiatorWidth/2.);
  SetMassHypotesis(0.139567);

  SetBetaOfParticle();
  SetFreonRefractiveIndex();

  beta   = GetBetaOfParticle();   
  nfreon = GetFreonRefractiveIndex();


  Float_t radius = FromEmissionToCathode();
  if (radius == 999.) return 999.;

  return 0;
}

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
  Float_t phipad = atan2(argY,argX); 

  SetPhiPoint(phipad);

}

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

void AliRICHRecon::FindWeightThetaCerenkov()
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

  cout << " thetac weighted -> " << weightThetaCerenkov << endl;
}


void AliRICHRecon::FlagPhotons()
{
  // flag photons

  Int_t nPhotonHough = 0;

  Float_t thetaCerenkov = GetThetaCerenkov();
  if(fDebug) cout << " fThetaCerenkov " << thetaCerenkov << endl;

  Float_t thetaDist= thetaCerenkov - fThetaMin;
  Int_t steps = (Int_t)(thetaDist / fDTheta);

  Float_t tmin = fThetaMin + (Float_t)(steps - 1)*fDTheta;
  Float_t tmax = fThetaMin + (Float_t)(steps)*fDTheta;
  Float_t tavg = 0.5*(tmin+tmax);

  tmin = tavg - 0.5*fWindowWidth;
  tmax = tavg + 0.5*fWindowWidth;

  if(fDebug) cout << " tmin " << tmin << " tmax " << tmax << endl;
  if(fDebug) cout << " thetac " << thetaCerenkov << endl;

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

void AliRICHRecon::DrawEvent(Int_t flag) const
{
  // draw event with rings

  flag=1; // dummy to be removed...
}

Float_t  AliRICHRecon::FindMassOfParticle()
{
  // find mass of the particle from theta cerenkov

  Float_t pmod = GetTrackMomentum();

  SetPhotonEnergy(6.85);
  SetFreonRefractiveIndex();

  Float_t thetaCerenkov = GetThetaCerenkov();
  FindBetaFromTheta(thetaCerenkov);

  Double_t beta = (Double_t)(GetBetaOfParticle());
  Double_t den = 1. - beta*beta;
  if(den<=0.) return 999.;

  Double_t gamma = 1./TMath::Sqrt(den);

  Float_t mass = pmod/(beta*(Float_t)gamma);

  return mass;
}


void AliRICHRecon::FillHistograms()
{
  // fill histograms..

  Float_t fittedTrackTheta, fittedTrackPhi;

  Float_t thetaCerenkov    = GetThetaCerenkov();
  if(thetaCerenkov == 999.) return;

  Float_t vertZ = GetEventVertexZ();

  Float_t trackTheta = GetTrackTheta();
  Float_t trackPhi   = GetTrackPhi();
  Float_t pmod       = GetTrackMomentum();
  Float_t pt         = GetTrackPt();
  Float_t trackEta   = GetTrackEta();
  Int_t q            = GetTrackCharge();
  Float_t tPCLastZ   = GetTrackTPCLastZ(); 
  Float_t minDist    = GetMinDist(); 

  fittedTrackTheta = GetFittedTrackTheta();
  fittedTrackPhi   = GetFittedTrackPhi();
  Int_t fittednPhotonHough = GetFittedHoughPhotons();
  
  if(fDebug)
    {
      cout << " p " << pmod  << " ThetaC " << thetaCerenkov 
	   << " rings " << fNrings << endl;
    }

  Int_t nPhotonHough     = GetHoughPhotons();
//  Float_t nPhotonHoughNorm = GetHoughPhotonsNorm();
  Float_t inRing = GetPortionOfRing();

  Float_t massOfParticle = FindMassOfParticle();

  Float_t houghArea = GetHoughArea();
  Float_t multiplicity = GetEventMultiplicity();


  Float_t var[20];

  var[0] = 0; 
  var[1] = 0;
  var[2] = vertZ;
  var[3] = pmod;
  var[4] = pt;
  var[5] = trackEta;
  var[6] = trackTheta;
  var[7] = trackPhi;
  var[8] = fittedTrackTheta;
  var[9] = fittedTrackPhi;
  var[10] = q;
  var[11] = thetaCerenkov;
  var[12] = (Float_t)nPhotonHough;
  var[13] = (Float_t)fittednPhotonHough;
  var[14] = inRing;
  var[15] = massOfParticle;
  var[16] = houghArea;
  var[17] = multiplicity;
  var[18] = tPCLastZ;
  var[19] = minDist;

  fNtuple->Fill(var);


  fittedTrackTheta = GetFittedTrackTheta();
  fittedTrackPhi = GetFittedTrackPhi();



  if(thetaCerenkov > 0.505 && thetaCerenkov < 0.605) {
      SetPhotonEnergy(6.85);
      SetFreonRefractiveIndex();
  }

  Int_t nPhotons = GetPhotonsNumber();

  for (Int_t j=0; j < nPhotons;j++)
    SetPhotonIndex(j);
}//FillHistograms()
//__________________________________________________________________________________________________
void AliRICHRecon::Minimization()
{
  // minimization to find the best theta and phi of the track

  Double_t arglist;
  Int_t ierflag = 0;

  static Double_t vstart[2];
  static Double_t lower[2], upper[2];
  static Double_t step[2]={0.001,0.001};

  Double_t trackThetaNew,trackPhiNew;
  TString chname;
  Double_t eps, b1, b2;
  Int_t ierflg;

  gAliRICHminuit = new TMinuit(2);
  gAliRICHminuit->SetObjectFit((TObject *)this);
  gAliRICHminuit->SetFCN(fcnrecon);
  gAliRICHminuit->mninit(5,10,7);

  vstart[0] = (Double_t)GetTrackTheta();
  vstart[1] = (Double_t)GetTrackPhi();

  lower[0] = vstart[0] - 0.03;
  if(lower[0] < 0) lower[0] = 0.;
  upper[0] = vstart[0] + 0.03;
  lower[1] = vstart[1] - 0.03;
  upper[1] = vstart[1] + 0.03;


  gAliRICHminuit->mnparm(0,"theta",vstart[0],step[0],lower[0],upper[0],ierflag);
  gAliRICHminuit->mnparm(1," phi ",vstart[1],step[1],lower[1],upper[1],ierflag);

  arglist = -1;

  //  gAliRICHminuit->FixParameter(0);

  gAliRICHminuit->SetPrintLevel(-1);
//  gAliRICHminuit->mnexcm("SET PRI",&arglist, 1, ierflag);
  gAliRICHminuit->mnexcm("SET NOGR",&arglist, 1, ierflag);
  gAliRICHminuit->mnexcm("SET NOW",&arglist, 1, ierflag);
  arglist = 1;
  gAliRICHminuit->mnexcm("SET ERR", &arglist, 1,ierflg);
  arglist = -1;

  //  gAliRICHminuit->mnscan();

//  gAliRICHminuit->mnexcm("SIMPLEX",&arglist, 0, ierflag);
  gAliRICHminuit->mnexcm("MIGRAD",&arglist, 0, ierflag);
  gAliRICHminuit->mnexcm("EXIT" ,&arglist, 0, ierflag);
  
  gAliRICHminuit->mnpout(0,chname, trackThetaNew, eps , b1, b2, ierflg);
  gAliRICHminuit->mnpout(1,chname, trackPhiNew, eps , b1, b2, ierflg);

  //values after the fit...
  SetFittedTrackTheta((Float_t)trackThetaNew);
  SetFittedTrackPhi((Float_t)trackPhiNew);

  delete gAliRICHminuit;

}

void AliRICHRecon::EstimationOfTheta()
{
  // theta estimate

  Int_t nPhotons = 0;

  Float_t shiftX = GetShiftX();
  Float_t shiftY = GetShiftY();

  Float_t *candidatePhotonX = GetCandidatePhotonX();
  Float_t *candidatePhotonY = GetCandidatePhotonY();

  Int_t nPhotonsCandidates = GetCandidatePhotonsNumber();

  //  cout << "MINIM: Nphotons " << nPhotonsCandidates << endl;

  for (Int_t j=0; j < nPhotonsCandidates; j++)
    {

      SetPhotonIndex(j);

      if(!GetPhotonFlag()) continue;

      Float_t xtoentr = candidatePhotonX[j] - shiftX;
      Float_t ytoentr = candidatePhotonY[j] - shiftY;

      SetEntranceX(xtoentr);
      SetEntranceY(ytoentr);

      FindPhiPoint();

      FindThetaPhotonCerenkov();

      Float_t thetaPhotonCerenkov = GetThetaPhotonCerenkov();

      //      cout << " ACCEPTED!!! " << thetaPhotonCerenkov << endl;

      SetPhotonEta(thetaPhotonCerenkov);

      nPhotons++;

    }

  Float_t xmean = 0.;
  Float_t x2mean = 0.;
  Int_t nev = 0;

  for (Int_t j=0; j < nPhotonsCandidates;j++)
    {
      SetPhotonIndex(j);

      Float_t eta = GetPhotonEta();

      if(eta != -999.) 
	{
	  if(GetPhotonFlag() == 2) 
	    {
	      xmean += eta;
	      x2mean += eta*eta;
	      nev++;
	    }
	}
    }

  if(nev > 0)
    {
      xmean /=(Float_t)nev;
      x2mean /=(Float_t)nev;
    } else {
      xmean = 0.;
      x2mean = 0.;
    }

  Float_t vRMS = sqrt(x2mean - xmean*xmean);

  //  cout << " RMS " << vRMS;

  SetEstimationOfTheta(xmean);
  SetEstimationOfThetaRMS(vRMS);
}

void fcnrecon(Int_t& /*npar*/, Double_t* /*gin*/, Double_t &f, Double_t *par, Int_t)
{
  // function to be minimized
  AliRICHRecon *gMyRecon = (AliRICHRecon*)gAliRICHminuit->GetObjectFit();

  Float_t p0 = (Float_t)par[0];
  Float_t p1 = (Float_t)par[1];

  gMyRecon->SetTrackTheta(p0);
  gMyRecon->SetTrackPhi(p1);

  gMyRecon->EstimationOfTheta();
  Float_t vRMS = gMyRecon->GetEstimationOfThetaRMS();

  Int_t houghPhotons = gMyRecon->GetHoughPhotons();


  f = (Double_t)(1000*vRMS/(Float_t)houghPhotons);

//   if(fDebug) cout << "   f   " << f
// 		  << " theta " << par[0] << " phi " << par[1] 
//                   << " HoughPhotons " << houghPhotons << endl;
//   
//   if(fDebug&&iflag == 3)
//     {
//             cout << " --- end convergence...summary --- " << endl;
//             cout << " theta " << par[0] << endl;
//             cout << "  phi  " << par[1] << endl;
//     }
}

void AliRICHRecon::Waiting()
{
  // wait, wait....
  if(!fIsDISPLAY) return;
  cout << " Press any key to continue...";

//  gSystem->ProcessEvents();
  getchar(); 

  cout << endl;

  return;
}

