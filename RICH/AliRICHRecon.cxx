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

#include "AliRICH.h"
#include "AliRICHRecon.h"
#include "AliRICHParam.h"
#include <AliLoader.h>
#include <AliStack.h>
#include <Riostream.h>
#include <TCanvas.h>
#include <TParticle.h>
#include <TStyle.h>
#include <TF1.h>
#include <TH2.h>
#include <TMath.h>
#include <TRandom.h>
#include <TMinuit.h>
#include <TNtuple.h>
#include <TMath.h>
#include <TGraph.h>
#include <TLine.h>
#include <TPolyLine.h>
#include <TMarker.h>
#include <TText.h>
#include <TProfile.h>
#include <TRotation.h>
#include <TSystem.h>
#include <TVector3.h>
#include <TEventList.h>  

#define NPointsOfRing 201

// Geometry of the RICH at Star...

static const Int_t nPadX      = AliRICHParam::NpadsY();
static const Int_t nPadY      = AliRICHParam::NpadsX();
static const Float_t PadSizeX = AliRICHParam::PadSizeY();
static const Float_t PadSizeY = AliRICHParam::PadSizeX();
static const Float_t spacer   = AliRICHParam::DeadZone();
static const Float_t degree = 180/3.1415926535;

static const Float_t pi = TMath::Pi();

static const Float_t RadiatorWidth = AliRICHParam::FreonThickness();
static const Float_t QuartzWidth   = AliRICHParam::QuartzThickness();
static const Float_t GapWidth      = AliRICHParam::RadiatorToPads();

static const Float_t fDTheta = 0.001;          // Step for sliding window
//static const Float_t fWindowWidth = 0.040;     // Hough width of sliding window
static const Float_t fWindowWidth = 0.060;     // Hough width of sliding window

static const Int_t fThetaBin = 750;            // Hough bins
static const Float_t fThetaMin = 0.0;          // Theta band min
static const Float_t fThetaMax = 0.75;         // Theta band max

static const Float_t Xmin = -AliRICHParam::PcSizeY()/2.;
static const Float_t Xmax =  AliRICHParam::PcSizeY()/2.;
static const Float_t Ymin = -AliRICHParam::PcSizeX()/2.;
static const Float_t Ymax =  AliRICHParam::PcSizeX()/2.;


// Global variables...

Bool_t fDebug      = kFALSE;
Bool_t kDISPLAY    = kFALSE;
Bool_t kWEIGHT     = kFALSE;
Bool_t kBACKGROUND = kFALSE;
Bool_t kMINIMIZER  = kFALSE;
//

Int_t TotEvents = 0;

static Float_t xGraph[3000],yGraph[3000];

static Int_t NRings = 0;
static Int_t NevTOT = 0;

TMinuit *gMyMinuit ;

void fcnrecon(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

// Float_t fEmissionPoint;
// Float_t fTrackTheta;
// Float_t fTrackPhi;
// Float_t fXtoentr;
// Float_t fYtoentr;

//

TFile *outputfile;

TH1F *h1_photons,*h1_photacc,*h1_hough;
TH2F *h2_tvsppos, *h2_tvspneg,*h2_func;

TH2F *h2_disp;

TH2F *h2_test1, *h2_test2, *h2_test4, *h2_testmap; 
TH2F *h2_dist_p;

TH1F *h1_photons1, *h1_photons2;
TH1F *h1_houghpos, *h1_houghneg;
TH1F *h1_mass;
TH2F *h2_mvsp;

TH1F *h1_hcs, *h1_hcsw;

TH1F *h1_nprotons;

TProfile *hp_1pos, *hp_1neg;
TProfile *hp_1posnorm, *hp_1negnorm;
TH2F *h2_1pos, *h2_1neg;
TH2F *h2_1posnorm, *h2_1negnorm;
TH2F *h2_mvst;

TH1F *h1_deltap, *h1_deltapop;
TH1F *h1_diffTrackTheta, *h1_diffTrackPhi;
TH1F *h1_photaccspread;

TH2F *h2_diffpos, *h2_diffneg;
TH2F *h2_map, *h2_mapw;

TH1F *photris;

TNtuple *hn;

TCanvas *StarCanvas,*Display,*Displayhcs;
TGraph *gra;
TLine *line;
TPolyLine *poll;
TPolyMarker *polm;
TMarker *Point, *TrackPoints, *Photon, *PhotonAcc;
TText *text;

AliRICHRecon::AliRICHRecon(const char*, const char*)
{

  fRich = (AliRICH*)gAlice->GetDetector("RICH");

}

void AliRICHRecon::InitRecon()
{

   outputfile = new TFile("Anal.root","RECREATE","My Analysis histos"); 
   if(kDISPLAY) Display = new TCanvas("Display","RICH Display",0,0,1200,750);      

   h1_photons = new TH1F("h1_photons","photons",750,0.,0.75);
   h1_photacc = new TH1F("h1_photacc","photons",750,0.,0.75);
   h1_hough   = new TH1F("h1_hough","hough",750,0.,0.75);
   h1_houghpos= new TH1F("h1_houghpos","hough",750,0.,0.75);
   h1_houghneg= new TH1F("h1_houghneg","hough",750,0.,0.75);

   h2_tvsppos = new TH2F("h2_tvsppos","thetac vs p",100,0.,5.,750,0.,0.75);
   h2_tvspneg = new TH2F("h2_tvspneg","thetac vs p",100,0.,5.,750,0.,0.75);
   h2_func    = new TH2F("h2_func"," func ",800,0.,0.8,100,-100.,100.);
   h2_mvsp    = new TH2F("h2_mvsp","mass vs p",100,0.,5.,200,0.,2.);
   h2_mvst    = new TH2F("h2_mvst","mass vs t",750,0.,0.75,200,0.,2.);
   h2_map     = new TH2F("h2_map","h2_map",160,0.,160.,96,0.,96.);
   h2_mapw    = new TH2F("h2_mapw","h2_mapw",160,0.,160.,96,0.,96.);

   h2_dist_p = new TH2F("h2_dist_p","h2_dist_p",100,0.,5.,100,0.,5.);
   //

   h2_disp = new TH2F("h2_disp","STAR-RICH Event Display",165,Xmin,Xmax,100,Ymin,Ymax);

   //   h2_test1  = new TH2F("h2_test1","test1 map",165,-64.,64.,100,-42.,42.);
   h2_test2  = new TH2F("h2_test2","test2 map",165,-64.,64.,100,-42.,42.);
   //   h2_test4  = new TH2F("h2_test4","test4 map",165,-64.,64.,100,-42.,42.);
   h2_testmap= new TH2F("h2_testmap","test map",165,-64.,64.,100,-42.,42.);

   //
   h1_photons1 = new TH1F("h1_photons1","photons",750,0.,0.75);
   h1_photons2 = new TH1F("h1_photons2","photons",750,0.,0.75);
   //
   h1_hcs  = new TH1F("h1_hcs","hcs",750,0.,750.);
   h1_hcsw = new TH1F("h1_hcsw","hcsw",750,0.,750.);
   //
   h1_nprotons = new TH1F("h1_nprotons","n prot",30,0.,30.);
   //
   hp_1pos = new TProfile("hp_1pos","Nphot vs thetac pos",250,0.,0.75); 
   hp_1neg = new TProfile("hp_1neg","Nphot vs thetac neg",250,0.,0.75); 
   hp_1posnorm = new TProfile("hp_1posnorm","Nphot vs thetac pos norm",250,0.,0.75); 
   hp_1negnorm = new TProfile("hp_1negnorm","Nphot vs thetac neg norm",250,0.,0.75); 
   //
   h2_1pos     = new TH2F("h2_1pos","Nphot vs p pos",100,0.,5.,30,0.,30.); 
   h2_1neg     = new TH2F("h2_1neg","Nphot vs p neg",100,0.,5.,30,0.,30.); 
   h2_1posnorm = new TH2F("h2_1posnorm","Nphot vs p pos norm",100,0.,5.,30,0.,30.); 
   h2_1negnorm = new TH2F("h2_1negnorm","Nphot vs p neg norm",100,0.,5.,30,0.,30.); 

   h1_deltap = new TH1F("h1_deltap","delta_p",200,-0.5,0.5);
   h1_deltapop = new TH1F("h1_deltapop","deltapop",200,-1.,1.);
   h1_diffTrackTheta = new TH1F("h1_diffTrackTheta","delta theta",200,-0.25,0.25);
   h1_diffTrackPhi   = new TH1F("h1_diffTrackPhi","delta phi",200,-0.25,0.25);

   h1_photaccspread = new TH1F("h1_photaccspread","photons spread",200,-0.1,0.1);

   //

   h1_mass = new TH1F("h1_mass","mass",200,0.,2.);
   photris = new TH1F("photris","photris",1000,0.,1.);
   h2_diffneg = new TH2F("h2_diffneg","diff neg",100,-2.5,2.5,100,-2.5,2.5);
   h2_diffpos = new TH2F("h2_diffpos","diff pos",100,-2.5,2.5,100,-2.5,2.5);

   hn = new TNtuple("hn","ntuple",
"Run:Trig:VertZ:Pmod:Pt:Eta:TrackTheta:TrackPhi:TrackThetaFit:TrackPhiFit:Charge:ThetaCerenkov:NPhotons:NPhotonsFit:InRing:MassOfParticle:HoughArea:Multiplicity:TPCLastZ");
}

void AliRICHRecon::StartProcessEvent()
{
  
  Float_t TrackThetaStored    = 0;
  Float_t TrackPhiStored      = 0;
  Float_t ThetaCerenkovStored = 0;
  Int_t HoughPhotonsStored    = 0;
  
  SetFreonScaleFactor(0.994);

  InitRecon();
  
  if(kDISPLAY) 
    {
      DrawEvent(0);
//      waiting();
    }

    Rich()->GetLoader()->LoadHits();
    Rich()->GetLoader()->LoadRecPoints();
    Rich()->GetLoader()->LoadDigits();
    gAlice->GetRunLoader()->LoadHeader();
    gAlice->GetRunLoader()->LoadKinematics();    
    
    Rich()->GetLoader()->TreeR()->GetEntry(0);

    Float_t clusX[7][500],clusY[7][500];
    Int_t clusQ[7][500];    
    Int_t nClusters[7];
    
    for (Int_t ich=0;ich<7;ich++) {
      nClusters[ich] = Rich()->ClustersOld(ich+1)->GetEntries();    
      for(Int_t k=0;k<nClusters[ich];k++) {
        AliRICHRawCluster *pCluster = (AliRICHRawCluster *)Rich()->ClustersOld(ich+1)->At(k);
        clusX[ich][k] = pCluster->fX;
        clusY[ich][k] = pCluster->fY;
        clusQ[ich][k] = pCluster->fQ;
//        pCluster->Print();
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

      cout << " iPrim " << iPrim << endl;
//      if(iPrim==0) return;
//      if(iPrim>1) Fatal("StartProcessEvent"," problems with prim to hit!!! = %3i", iPrim);
      
      TParticle *pParticle = gAlice->GetRunLoader()->Stack()->Particle(pHit->GetTrack());
      Float_t pmod     = pParticle->P();
      Float_t pt       = pParticle->Pt();
      Float_t TrackEta = pParticle->Eta();
      Int_t q          = (Int_t)TMath::Sign(1.,pParticle->GetPDG()->Charge());        

      pParticle->Print();
      
      cout << " pmod " << pmod << " pt " << pt << " Eta " << TrackEta << " charge " << q << endl;
      
      SetTrackMomentum(pmod); 
      SetTrackPt(pt);
      SetTrackEta(TrackEta);
      SetTrackCharge(q);

      TVector3 pGlob(pHit->MomFreoX(),pHit->MomFreoY(),pHit->MomFreoZ());
      TVector3 pLocal = Rich()->C(pHit->Chamber())->Global2Local(pGlob,1);
      
      Float_t primGlobalX = pHit->X();
      Float_t primGlobalY = pHit->Y();
      Float_t primGlobalZ = pHit->Z();
      TVector3 primGlobal(primGlobalX,primGlobalY,primGlobalZ);
      TVector3 primLocal = Rich()->C(pHit->Chamber())->Global2Local(primGlobal);
      
//      Float_t pmodFreo = pLocal.Mag();
      Float_t TrackTheta = pLocal.Theta();
      Float_t TrackPhi = pLocal.Phi();

      cout << " TrackTheta " << TrackTheta << " TrackPhi " << TrackPhi << endl;
      
      SetTrackTheta(TrackTheta);
      SetTrackPhi(TrackPhi);
 
      Int_t MaxInd = 0;
      Float_t MinDist =  999.;

      cout << " n Clusters " << nClusters[pHit->Chamber()-1] << " for chamber n. " << pHit->Chamber() << endl;
      
      for(Int_t j=0;j<nClusters[pHit->Chamber()-1];j++)
	{
	  Float_t diffx = primLocal.X() - clusX[pHit->Chamber()-1][j];
	  Float_t diffy = primLocal.Y() - clusY[pHit->Chamber()-1][j];

          cout << " cluster x " << clusX[pHit->Chamber()-1][j] << " hit track x " << primLocal.X();
          cout << " cluster y " << clusY[pHit->Chamber()-1][j] << " hit track y " << primLocal.Y() << endl;
          
          Float_t diff = sqrt(diffx*diffx + diffy*diffy);

	  if(diff < MinDist)
	    {
	      MinDist = diff;
	      MaxInd = j;
	    }

	}

      Float_t diffx = primLocal.X() - clusX[pHit->Chamber()-1][MaxInd];
      Float_t diffy = primLocal.Y() - clusY[pHit->Chamber()-1][MaxInd];

      cout << " diffx " << diffx << " diffy " << diffy << endl;
      
      if(q>0)
      {
         h2_diffpos->Fill(diffx,diffy);
      } else {
         h2_diffneg->Fill(diffx,diffy);
      }

      SetMipIndex(MaxInd);
      SetTrackIndex(i);

      Float_t ShiftX = primLocal.X()/primLocal.Z()*(RadiatorWidth+QuartzWidth+GapWidth) + primLocal.X();
      Float_t ShiftY = primLocal.Y()/primLocal.Z()*(RadiatorWidth+QuartzWidth+GapWidth) + primLocal.Y();
      
      SetShiftX(ShiftX);
      SetShiftY(ShiftY);

      Float_t *pclusX = &clusX[pHit->Chamber()-1][0];
      Float_t *pclusY = &clusY[pHit->Chamber()-1][0];
      
      SetCandidatePhotonX(pclusX);
      SetCandidatePhotonY(pclusY);
      SetCandidatePhotonsNumber(nClusters[pHit->Chamber()-1]);

      Int_t qch = clusQ[pHit->Chamber()-1][MaxInd];

      if(MinDist < 3.0 && qch > 120 && MaxInd !=0) 
	{
	  
	  if(kBACKGROUND)
	    {
	      
	      Float_t Xrndm = Xmin + (Xmax-Xmin)*gRandom->Rndm(280964);
	      Float_t Yrndm = Ymin + (Ymax-Ymin)*gRandom->Rndm(280964);

	      cout << " Xrndm " << Xrndm << " Yrndm " << Yrndm << endl;

	      SetShiftX(Xrndm);
	      SetShiftY(Yrndm);
	      
	    }

	  PatRec();

	  TrackThetaStored = GetTrackTheta();
	  TrackPhiStored = GetTrackPhi();
	  ThetaCerenkovStored = GetThetaCerenkov();
	  HoughPhotonsStored = GetHoughPhotons();
	  
          Int_t DiffNPhotons = 999;
          Int_t Nsteps = 0;
          Float_t DiffTrackTheta = 999.;
          Float_t DiffTrackPhi   = 999.;

	  while( kMINIMIZER && GetHoughPhotons() > 2 
                            && DiffNPhotons !=0 
                            && DiffTrackTheta > 0.0001
                            && Nsteps < 10)
	    {

	      Int_t   HoughPhotonsBefore  = GetHoughPhotons();

	      Float_t TrackThetaBefore = GetTrackTheta();
	      Float_t TrackPhiBefore   = GetTrackPhi();
	  
	      Minimization(); 

              PatRec();
 
              DiffNPhotons = abs(HoughPhotonsBefore - GetHoughPhotons()); 

	      Float_t TrackThetaAfter = GetTrackTheta();
	      Float_t TrackPhiAfter   = GetTrackPhi();

              DiffTrackTheta = TMath::Abs(TrackThetaAfter - TrackThetaBefore);
              DiffTrackPhi   = TMath::Abs(TrackPhiAfter - TrackPhiBefore);

              if(fDebug)
              cout << " HoughPhotonsBefore " << HoughPhotonsBefore
                   << " GetHoughPhotons()  " << GetHoughPhotons();

              Nsteps++;
	    }

	  SetFittedThetaCerenkov(GetThetaCerenkov());
	  SetFittedHoughPhotons(GetHoughPhotons());

	  SetTrackTheta(TrackThetaStored);
	  SetTrackPhi(TrackPhiStored);
	  SetThetaCerenkov(ThetaCerenkovStored);
	  SetHoughPhotons(HoughPhotonsStored);

          SetMinDist(MinDist);

	  FillHistograms();
      
	  if(kDISPLAY) DrawEvent(1);

	  waiting();

	}
    }
  //
  if(kDISPLAY) Display->Print("display.ps");
}


void AliRICHRecon::EndProcessEvent()
{
// function called at the end of the event loop

  printf("Processed events: %d Total events: %d \n",TotEvents,NevTOT); 

  outputfile->Write();
  outputfile->Close();                                                     
}

void AliRICHRecon::PatRec()
{

  Float_t TrackTheta = GetTrackTheta();
  Float_t TrackPhi   = GetTrackPhi();
  Float_t pmod       = GetTrackMomentum();
  //  Int_t q            = GetTrackCharge();

  //  Int_t TrackIndex = GetTrackIndex();
  Int_t MipIndex   = GetMipIndex();

  Bool_t kPatRec = kFALSE;  

  Int_t CandidatePhotons = 0;

  Float_t ShiftX = GetShiftX();
  Float_t ShiftY = GetShiftY();

  Float_t* CandidatePhotonX = GetCandidatePhotonX();
  Float_t* CandidatePhotonY = GetCandidatePhotonY();

  Int_t CandidatePhotonsNumber = GetCandidatePhotonsNumber();

  if(fDebug) cout << " n " << CandidatePhotonsNumber << endl;

  SetThetaCerenkov(999.);
  SetHoughPhotons(0);
  SetHoughPhotonsNorm(0);
  SetHoughRMS(999.);

  for (Int_t j=0; j < CandidatePhotonsNumber; j++)
    {

      SetPhotonIndex(j);

      SetPhotonFlag(0);
      SetPhotonEta(-999.);
      SetPhotonWeight(0.);

      if (j == MipIndex) continue;

      //      h2_test1->Fill(CandidatePhotonX[j],CandidatePhotonY[j]);
        
      if(CandidatePhotonX[j] < -64.) continue; /* avoid artificial clusters from edge uesd by Yale.... */

      Float_t Xtoentr = CandidatePhotonX[j] - ShiftX;
      Float_t Ytoentr = CandidatePhotonY[j] - ShiftY;

      //      Float_t chargehit = fHits_charge[j]; 
      //      if(chargehit > 150) continue;

      SetEntranceX(Xtoentr);
      SetEntranceY(Ytoentr);

      FindPhiPoint();

      Int_t PhotonStatus = PhotonInBand();
 
      if(fDebug)
         {
            cout << " Photon n. " << j << " Status " << PhotonStatus << " accepted " << endl;
            cout << " CandidatePhotonX[j] " << CandidatePhotonX[j] << " CandidatePhotonY[j] " << CandidatePhotonY[j] << endl;
         }
    
      if(PhotonStatus == 0) continue;

      SetPhotonFlag(1);

      FindThetaPhotonCerenkov();

      Float_t ThetaPhotonCerenkov = GetThetaPhotonCerenkov();

      if(fDebug) cout << " theta photon " << ThetaPhotonCerenkov << endl;

      SetPhotonEta(ThetaPhotonCerenkov);

      CandidatePhotons++;

      // fill histograms

      //      h2_test4->Fill(CandidatePhotonX[j],CandidatePhotonY[j]);

      //      if(kDISPLAY) h1_photons->Fill(ThetaPhotonCerenkov);
      
    }

  if(CandidatePhotons >= 1) kPatRec = kTRUE;

  if(!kPatRec) return;
    {
       SetThetaCerenkov(999.);
       SetHoughPhotons(0);
    }
  SetPhotonsNumber(CandidatePhotonsNumber);

  HoughResponse();
  
  NRings++;

  FlagPhotons();
  Int_t NPhotonHough = GetHoughPhotons();
 
  if(NPhotonHough < 1) 
    {
      SetThetaCerenkov(999.);
      SetHoughPhotonsNorm(0.);
      return;
    }

  if(kWEIGHT) FindWeightThetaCerenkov();

  Float_t ThetaCerenkov = GetThetaCerenkov();

  SetThetaOfRing(ThetaCerenkov);
  FindAreaAndPortionOfRing();

  Float_t NPhotonHoughNorm = ((Float_t)NPhotonHough)/GetPortionOfRing();
  SetHoughPhotonsNorm(NPhotonHoughNorm);

  // Calculate the area where the photon are accepted...

  Float_t ThetaInternal = ThetaCerenkov - 0.5*fWindowWidth; 
  SetThetaOfRing(ThetaInternal);
  FindAreaAndPortionOfRing();
  Float_t InternalArea = GetAreaOfRing();

  Float_t ThetaExternal = ThetaCerenkov + 0.5*fWindowWidth; 
  SetThetaOfRing(ThetaExternal);
  FindAreaAndPortionOfRing();
  Float_t ExternalArea = GetAreaOfRing();

  Float_t HoughArea = ExternalArea - InternalArea;

  SetHoughArea(HoughArea);

  if(fDebug)
    {
      cout << " ----- SUMMARY OF RECONSTRUCTION ----- " << endl; 
      cout << " Rings found " << NRings << " with thetac " << ThetaCerenkov << endl;
      
      h1_hough->Fill(ThetaCerenkov,1.);
      
      cout << " Nphotons " << GetPhotonsNumber() 
	   << " Hough    " << NPhotonHough 
	   << " norm     " << NPhotonHoughNorm << endl;
      
      cout << " In PatRec:p " << pmod << " theta " << TrackTheta << " phi " << TrackPhi << endl;
      cout << " ------------------------------------- " << endl; 
    }

  Int_t NPhotons = GetPhotonsNumber();

  Float_t xmean = 0.;
  Float_t x2mean = 0.;
  Int_t nev = 0;

  for (Int_t j=0; j < NPhotons;j++)
    {
      SetPhotonIndex(j);

      Float_t eta = GetPhotonEta();

      if(eta != -999.) 
	{
	  if(GetPhotonFlag() == 2) 
	    {

	      if(pmod>2.5&&ThetaCerenkov>0.65) photris->Fill(eta);

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

  Float_t RMS = sqrt(x2mean - xmean*xmean);

  SetHoughRMS(RMS);

  if(fDebug) cout << " RMS " << RMS << endl;

}

void AliRICHRecon::FindEmissionPoint()
{ 

// Find emission point

  Float_t absorbtionLenght=7.83*RadiatorWidth; //absorption length in the freon (cm)
  // 7.83 = -1/ln(T0) where 
  // T0->Trasmission freon at 180nm = 0.88 (Eph=6.85eV)
  Float_t photonLenght, photonLenghtMin, photonLenghtMax;

  photonLenght=exp(-RadiatorWidth/(absorbtionLenght*cos(fCerenkovAnglePad)));
  photonLenghtMin=RadiatorWidth*photonLenght/(1.-photonLenght);
  photonLenghtMax=absorbtionLenght*cos(fCerenkovAnglePad);
  Float_t EmissionPoint = RadiatorWidth + photonLenghtMin - photonLenghtMax;

  SetEmissionPoint(EmissionPoint);
}


Int_t AliRICHRecon::PhotonInBand()
{

  //  Float_t MassOfParticle;
  Float_t beta;
  Float_t nfreon;

  Float_t thetacer;

  Float_t Xtoentr = GetEntranceX();
  Float_t Ytoentr = GetEntranceY();

  Float_t InnerRadius;
  Float_t OuterRadius;

  Float_t phpad = GetPhiPoint();

  //  Float_t pmod = GetTrackMomentum();
  //  Float_t TrackTheta = GetTrackTheta();
  //  Float_t TrackPhi = GetTrackPhi();

  // inner radius //
  SetPhotonEnergy(5.6);
  SetEmissionPoint(RadiatorWidth -0.0001);
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
      InnerRadius = -999.;
      SetXInnerRing(-999.);
      SetYInnerRing(-999.);
      SetRadiusInnerRing(-999.);
    }
  else
    {
      SetThetaPhotonInDRS(GetThetaAtQuartz());
      SetPhiPhotonInDRS(phpad);

      InnerRadius = FromEmissionToCathode();
       if(InnerRadius == 999.) InnerRadius = -999.;
      
      SetXInnerRing(GetXPointOnCathode());
      SetYInnerRing(GetYPointOnCathode());
      SetRadiusInnerRing(InnerRadius);
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
      OuterRadius = 999.;
      SetXOuterRing(999.);
      SetYOuterRing(999.);
      SetRadiusOuterRing(999.);
    }
  else
    {
      SetThetaPhotonInDRS(GetThetaAtQuartz());
      SetPhiPhotonInDRS(phpad);

      OuterRadius = FromEmissionToCathode();
//      cout << " OuterRadius " << OuterRadius << endl;
      SetXOuterRing(GetXPointOnCathode());
      SetYOuterRing(GetYPointOnCathode());
      SetRadiusOuterRing(OuterRadius);
    }

  Float_t padradius = sqrt(TMath::Power(Xtoentr,2)+TMath::Power(Ytoentr,2));
  
  if(fDebug) printf(" rmin %f r %f rmax %f \n",InnerRadius,padradius,OuterRadius);

  if(padradius>=InnerRadius && padradius<=OuterRadius) return 1;
  return 0;
}

void AliRICHRecon::FindThetaAtQuartz(Float_t ThetaCerenkov)
{

  if(ThetaCerenkov == 999.) 
    {
      SetThetaAtQuartz(999.);
      return;
    }

  Float_t ThetaAtQuartz = 999.;

  Float_t TrackTheta = GetTrackTheta();

  if(TrackTheta == 0) {

    if(fDebug) cout << " Theta sol unique " << ThetaCerenkov << endl;  

    ThetaAtQuartz = ThetaCerenkov;
    SetThetaAtQuartz(ThetaAtQuartz);
    return;
  }

  Float_t TrackPhi   = GetTrackPhi();
  Float_t PhiPoint = GetPhiPoint();

  Double_t den = TMath::Sin((Double_t)TrackTheta)
    *TMath::Cos((Double_t)TrackPhi)
    *TMath::Cos((Double_t)PhiPoint) +
    TMath::Sin((Double_t)TrackTheta)
    *TMath::Sin((Double_t)TrackPhi)
    *TMath::Sin((Double_t)PhiPoint); 
  Double_t b = TMath::Cos((Double_t)TrackTheta)/den;
  Double_t c = -TMath::Cos((Double_t)ThetaCerenkov)/den;

  Double_t UnderSqrt = 1 + b*b - c*c;

  if(fDebug)
    {
      cout << " TrackTheta    " << TrackTheta    << endl;
      cout << " TrackPhi      " << TrackPhi      << endl;
      cout << " PhiPoint      " << PhiPoint      << endl;
      cout << " ThetaCerenkov " << ThetaCerenkov << endl;
      cout << " den b c " << den << " b " << b << " c " << c << endl;
    }

  if(UnderSqrt < 0) {
    if(fDebug) cout << " sqrt negative !!!!" << UnderSqrt << endl;
    SetThetaAtQuartz(999.);
    return;
  }

  Double_t sol1 = (1+TMath::Sqrt(UnderSqrt))/(b-c);
  Double_t sol2 = (1-TMath::Sqrt(UnderSqrt))/(b-c);

  Double_t ThetaSol1 = 2*TMath::ATan(sol1);
  Double_t ThetaSol2 = 2*TMath::ATan(sol2);

  if(fDebug) cout << " Theta sol 1 " << ThetaSol1 
		  << " Theta sol 2 " << ThetaSol2 << endl;  

  if(ThetaSol1>0 && ThetaSol1 < pi) ThetaAtQuartz = (Float_t)ThetaSol1;
  if(ThetaSol2>0 && ThetaSol2 < pi) ThetaAtQuartz = (Float_t)ThetaSol2;

  SetThetaAtQuartz(ThetaAtQuartz);
}

void AliRICHRecon::FindThetaPhotonCerenkov()
{

  Float_t ThetaCerMin = 0.;
  Float_t ThetaCerMax = 0.75;
  Float_t ThetaCerMean;

  Float_t RadiusMin, RadiusMax, RadiusMean;
  Int_t nIteration = 0;

  const Float_t Tollerance = 0.05;

  //  Float_t pmod = GetTrackMomentum();
  //  Float_t TrackTheta = GetTrackTheta();
  //  Float_t TrackPhi = GetTrackPhi();

  Float_t PhiPoint = GetPhiPoint();

  SetPhotonEnergy(6.85);
  SetEmissionPoint(RadiatorWidth/2);

  Float_t XPoint = GetEntranceX();
  Float_t YPoint = GetEntranceY();
  Float_t DistPoint = sqrt(XPoint*XPoint + YPoint*YPoint);

  if(fDebug) cout << " DistPoint " << DistPoint << endl;

  // Star minimization...

  // First value...

  FindThetaAtQuartz(ThetaCerMin);
  
  if(GetThetaAtQuartz() == 999.)
    {
      RadiusMin = -999.;
    }
  else
    {
      SetThetaPhotonInDRS(GetThetaAtQuartz());
      SetPhiPhotonInDRS(PhiPoint);
      
      RadiusMin = FromEmissionToCathode();
    }

  // Second value...

  FindThetaAtQuartz(ThetaCerMax);
  if(GetThetaAtQuartz() == 999.)
    {
      RadiusMax = 999.;
    }
  else
    {
      SetThetaPhotonInDRS(GetThetaAtQuartz());
      SetPhiPhotonInDRS(PhiPoint);
      
      RadiusMax = FromEmissionToCathode();
    }
  // Mean value...

  ThetaCerMean = (ThetaCerMax + ThetaCerMin)/2;

  FindThetaAtQuartz(ThetaCerMean);
  if(GetThetaAtQuartz() == 999.)
    {
      RadiusMean = 999.;
    }
  else
    {
      SetThetaPhotonInDRS(GetThetaAtQuartz());
      SetPhiPhotonInDRS(PhiPoint);
      
      RadiusMean = FromEmissionToCathode();
    }

  if(fDebug) cout << " r1 " << RadiusMin << " rmean " 
		  << RadiusMean << " r2 " << RadiusMax << endl;

  while (TMath::Abs(RadiusMean-DistPoint) > Tollerance)
    {

      if((RadiusMin-DistPoint)*(RadiusMean-DistPoint) < 0) ThetaCerMax = ThetaCerMean;
      if((RadiusMin-DistPoint)*(RadiusMean-DistPoint) > 0) {

	ThetaCerMin = ThetaCerMean;

	FindThetaAtQuartz(ThetaCerMin);
	SetThetaPhotonInDRS(GetThetaAtQuartz());
	SetPhiPhotonInDRS(PhiPoint);

	RadiusMin =FromEmissionToCathode();
      }

      ThetaCerMean = (ThetaCerMax + ThetaCerMin)/2;

      FindThetaAtQuartz(ThetaCerMean);
      SetThetaPhotonInDRS(GetThetaAtQuartz());
      SetPhiPhotonInDRS(PhiPoint);

      RadiusMean = FromEmissionToCathode();

      nIteration++;
      if(nIteration>=50) {
	if(fDebug) printf(" max iterations in FindPhotonCerenkov\n");
	SetThetaPhotonCerenkov(999.);
	return;
      }
    }

  SetThetaPhotonCerenkov(ThetaCerMean);

}

void AliRICHRecon::FindAreaAndPortionOfRing()
{

  Float_t XPoint[NPointsOfRing], YPoint[NPointsOfRing];

  //  Float_t Xtoentr = GetEntranceX();
  //  Float_t Ytoentr = GetEntranceY();
  Float_t ShiftX = GetShiftX();
  Float_t ShiftY = GetShiftY();

  Float_t XEmiss = GetXCoordOfEmission(); 
  Float_t YEmiss = GetYCoordOfEmission(); 

  Float_t x0 = XEmiss + ShiftX;
  Float_t y0 = YEmiss + ShiftY;

  //  Float_t pmod = GetTrackMomentum();
  //  Float_t TrackTheta = GetTrackTheta();
  //  Float_t TrackPhi = GetTrackPhi();

  SetPhotonEnergy(6.85);
  SetFreonRefractiveIndex();

  SetEmissionPoint(RadiatorWidth/2.);

  Float_t Theta = GetThetaOfRing();
  
  Int_t nPoints = 0;
  Int_t NPsiAccepted = 0;
  Int_t NPsiTotal = 0;

  for(Int_t i=0;i<NPointsOfRing-1;i++)
    {

      Float_t Psi = 2*TMath::Pi()*i/NPointsOfRing;
      
      SetThetaPhotonInTRS(Theta);
      SetPhiPhotonInTRS(Psi);
      FindPhotonAnglesInDRS();
      
      Float_t Radius = FromEmissionToCathode();
      if (Radius == 999.) continue;
      
      NPsiTotal++;

      Float_t XPointRing = GetXPointOnCathode() + ShiftX;
      Float_t YPointRing = GetYPointOnCathode() + ShiftY;
      
      SetDetectorWhereX(XPointRing);
      SetDetectorWhereY(YPointRing);
      
      Int_t Zone = CheckDetectorAcceptance();

//      cout << " XPointing " << XPointRing << " YPointing " << YPointRing << " Zone " << Zone << endl;
//      cout << " ShiftX " << ShiftX << " ShiftY " << ShiftY << endl;
//      cout << " GetXPointOnCathode() " << GetXPointOnCathode() << endl;
//      cout << " GetYPointOnCathode() " << GetYPointOnCathode() << endl;

      if (Zone != 0) 
	{
	  FindIntersectionWithDetector();
	  XPoint[nPoints] = GetIntersectionX();
	  YPoint[nPoints] = GetIntersectionY();
	}
      else
	{
	  XPoint[nPoints] = XPointRing;
	  YPoint[nPoints] = YPointRing;
	  NPsiAccepted++;
	}

      nPoints++;

    }

  XPoint[nPoints] = XPoint[0];
  YPoint[nPoints] = YPoint[0];
  
  // find area...

  Float_t Area = 0;

  for (Int_t i = 0; i < nPoints; i++)
    {
      Area += TMath::Abs((XPoint[i]-x0)*(YPoint[i+1]-y0) - (XPoint[i+1]-x0)*(YPoint[i]-y0));
    }
  
  Area *= 0.5;
  
  Float_t PortionOfRing = ((Float_t)NPsiAccepted)/((Float_t)(NPsiTotal));

  //  cout << " Area " << Area << " Portion of ring " << PortionOfRing << endl;

  SetAreaOfRing(Area);
  SetPortionOfRing(PortionOfRing);
}

void AliRICHRecon::FindIntersectionWithDetector()
{

  Float_t XIntersect, YIntersect;
  Float_t x1, x2, y1, y2;

  Float_t ShiftX = GetShiftX();
  Float_t ShiftY = GetShiftY();

  Float_t XPoint = GetXPointOnCathode() + ShiftX;
  Float_t YPoint = GetYPointOnCathode() + ShiftY;

  Float_t XEmiss = GetXCoordOfEmission(); 
  Float_t YEmiss = GetYCoordOfEmission(); 

  Float_t Phi = GetPhiPhotonInDRS();
  Float_t m = tan(Phi);

  Float_t x0 = XEmiss + ShiftX;
  Float_t y0 = YEmiss + ShiftY;

  if(XPoint > x0)
    {
      x1 = x0;
      x2 = XPoint;
    }
  else
    {
      x2 = x0;
      x1 = XPoint;
    }
  if(YPoint > y0)
    {
      y1 = y0;
      y2 = YPoint;
    }
  else
    {
      y2 = y0;
      y1 = YPoint;
    }
  //
  XIntersect = Xmax;
  YIntersect = m*(XIntersect - x0) + y0;
  if (YIntersect >= Ymin && YIntersect <= Ymax && XIntersect >= x1 && XIntersect <= x2)
    {
      SetIntersectionX(XIntersect);
      SetIntersectionY(YIntersect);
      return;
    }
  //
  XIntersect = Xmin;
  YIntersect = m*(XIntersect - x0) + y0;
  if (YIntersect >= Ymin && YIntersect <= Ymax && XIntersect >= x1 && XIntersect <= x2)
    {
      SetIntersectionX(XIntersect);
      SetIntersectionY(YIntersect);
      return;
    }
  //
  YIntersect = Ymax;
  XIntersect = (YIntersect - y0)/m + x0;
  if (XIntersect >= Xmin && XIntersect <= Xmax && YIntersect >= y1 && YIntersect <= y2)
    {
      SetIntersectionX(XIntersect);
      SetIntersectionY(YIntersect);
      return;
    }
  //
  YIntersect = Ymin;
  XIntersect = (YIntersect - y0)/m + x0;
  if (XIntersect >= Xmin && XIntersect <= Xmax && YIntersect >= y1 && YIntersect <= y2)
    {
      SetIntersectionX(XIntersect);
      SetIntersectionY(YIntersect);
      return;
    }
  
  cout << " sono fuori!!!!!!" << endl;
//  cout << " x1 " << x1 << " x2 " << x2 << endl;
//  cout << " y1 " << y1 << " y2 " << y2 << endl;
//  cout << " Xmin " << Xmin << " Xmax " << Xmax << endl;
//  cout << " Ymin " << Ymin << " Ymax " << Ymax << endl;
  
}

Int_t AliRICHRecon::CheckDetectorAcceptance()
{

  // crosses X -2.6 2.6 cm
  // crosses Y -1 1 cm

  Float_t Xcoord = GetDetectorWhereX();
  Float_t Ycoord = GetDetectorWhereY();

//  cout << " Xcoord " << Xcoord << " Ycoord " << Ycoord << endl;
  if(Xcoord > Xmax)
    {
      if(Ycoord > Ymax) return 2;
      if(Ycoord > Ymin && Ycoord < Ymax) return 3;
      if(Ycoord < Ymin) return 4;
    }
  if(Xcoord < Xmin)
    {
      if(Ycoord > Ymax) return 8;
      if(Ycoord > Ymin && Ycoord < Ymax) return 7;
      if(Ycoord < Ymin) return 6;
    }
  if(Xcoord > Xmin && Xcoord < Xmax)
    {
      if(Ycoord > Ymax) return 1;
      if(Ycoord > Ymin && Ycoord < Ymax) return 0;
      if(Ycoord < Ymin) return 5;
    }
  return 999;
}

void AliRICHRecon::DrawRing()
{

  //  Float_t xGraph[1000],yGraph[1000];

  Float_t type;
  //  Float_t MassOfParticle;
  Float_t beta;
  Float_t nfreon;

  Float_t ThetaCerenkov;

  //  Float_t Xtoentr = GetEntranceX();
  //  Float_t Ytoentr = GetEntranceY();

  //  Float_t pmod = GetTrackMomentum();
  //  Float_t TrackTheta = GetTrackTheta();
  //  Float_t TrackPhi = GetTrackPhi();

  SetPhotonEnergy(6.85);
  SetFreonRefractiveIndex();

  SetEmissionPoint(RadiatorWidth/2.);

  type = 1;

  if(type == 1)
    {
      SetMassHypotesis(0.139567);
      SetBetaOfParticle();
      
      beta   = GetBetaOfParticle();   
      
    }
  else if(type == 2)
    {
      ThetaCerenkov = GetThetaCerenkov();
      FindBetaFromTheta(ThetaCerenkov);
    }
  
  nfreon = GetFreonRefractiveIndex();
  
  Float_t thetacer = Cerenkovangle(nfreon,beta);

  if(fDebug) cout << " TetaCer in DrawRing " << thetacer << endl;

  Int_t nPoints = 100;

  Int_t nPointsToDraw = 0;
  for(Int_t i=0;i<nPoints;i++)
    {
      Float_t phpad = 2*TMath::Pi()*i/nPoints;
      SetThetaPhotonInTRS(thetacer);
      SetPhiPhotonInTRS(phpad);
      FindPhotonAnglesInDRS();
      Float_t Radius = FromEmissionToCathode();
      if (Radius == 999.) continue;
      xGraph[nPointsToDraw] = GetXPointOnCathode() + GetShiftX();
      yGraph[nPointsToDraw] = GetYPointOnCathode() + GetShiftY();
      //      cout << " get shift X " << GetShiftX() << endl;
      //      cout << " get shift Y " << GetShiftY() << endl;
      nPointsToDraw++;
    }


  if(fDebug) cout << " Drawing the Ring... with " << nPointsToDraw << " points " << endl;

  //  pol = new TPolyLine(nPointsToDraw,xGraph,yGraph);
  //  pol->Draw("same");
  gra = new TGraph(nPointsToDraw,xGraph,yGraph);
  gra->Draw("AC");
  StarCanvas->Update();
 
}

Float_t AliRICHRecon::PhotonPositionOnCathode()
{ 
  //  Float_t MassOfParticle;
  Float_t beta;
  Float_t nfreon;

  //  Float_t pmod = GetTrackMomentum();
  //  Float_t TrackTheta = GetTrackTheta();
  //  Float_t TrackPhi = GetTrackPhi();

  //  Float_t phpad = GetPhiPoint();

  SetPhotonEnergy(6.85);
  SetEmissionPoint(RadiatorWidth/2.);
  SetMassHypotesis(0.139567);

  SetBetaOfParticle();
  SetFreonRefractiveIndex();

  beta   = GetBetaOfParticle();   
  nfreon = GetFreonRefractiveIndex();

  //  Float_t thetacer = Cerenkovangle(nfreon,beta);

  //  cout << " FromEmissionToCathode: thetacer " << thetacer << " phpad " << phpad << endl;

  Float_t Radius = FromEmissionToCathode();
  if (Radius == 999.) return 999.;

  //  Float_t Xphoton = GetXPointOnCathode();
  //  Float_t Yphoton = GetYPointOnCathode();
  //  cout << " PhotonPositionOnCathode: Xphoton " << Xphoton << " Yphoton " << Yphoton <<
  //  " Radius for photon " << Radius << endl;
  return 0;
}

void AliRICHRecon::FindPhotonAnglesInDRS()
{
  // Setup the rotation matrix of the track...

  TRotation Mtheta;
  TRotation Mphi;
  TRotation Minv;
  TRotation Mrot;
  
  Float_t TrackTheta = GetTrackTheta();
  Float_t TrackPhi = GetTrackPhi();

  Mtheta.RotateY(TrackTheta);
  Mphi.RotateZ(TrackPhi);
  
  Mrot = Mphi * Mtheta;
  //  Minv = Mrot.Inverse();

  TVector3 PhotonInRadiator(1,1,1);

  Float_t ThetaCerenkov = GetThetaPhotonInTRS();
  Float_t PhiCerenkov   = GetPhiPhotonInTRS();

  PhotonInRadiator.SetTheta(ThetaCerenkov);
  PhotonInRadiator.SetPhi(PhiCerenkov);
  PhotonInRadiator = Mrot * PhotonInRadiator;
  Float_t Theta = PhotonInRadiator.Theta();
  Float_t Phi = PhotonInRadiator.Phi();
  SetThetaPhotonInDRS(Theta);
  SetPhiPhotonInDRS(Phi);

}

Float_t AliRICHRecon::FromEmissionToCathode()
{

  Float_t nfreon, nquartz, ngas; 

  SetFreonRefractiveIndex();
  SetQuartzRefractiveIndex();
  SetGasRefractiveIndex();

  nfreon  = GetFreonRefractiveIndex();
  nquartz = GetQuartzRefractiveIndex();
  ngas    = GetGasRefractiveIndex();

  Float_t TrackTheta = GetTrackTheta();
  Float_t TrackPhi = GetTrackPhi();
  Float_t LengthOfEmissionPoint = GetEmissionPoint();

  Float_t Theta = GetThetaPhotonInDRS();
  Float_t Phi   = GetPhiPhotonInDRS();

//   cout << " Theta " << Theta << " Phi " << Phi << endl;

  Float_t xEmiss = LengthOfEmissionPoint*tan(TrackTheta)*cos(TrackPhi);
  Float_t yEmiss = LengthOfEmissionPoint*tan(TrackTheta)*sin(TrackPhi);

  SetXCoordOfEmission(xEmiss);
  SetYCoordOfEmission(yEmiss);
  
  Float_t thetaquar = SnellAngle(nfreon, nquartz, Theta);

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

  Float_t xw = (RadiatorWidth - LengthOfEmissionPoint)*cos(Phi)*tan(Theta);
  Float_t xq = QuartzWidth*cos(Phi)*tan(thetaquar);
  Float_t xg = GapWidth*cos(Phi)*tan(thetagap);
  Float_t yw = (RadiatorWidth - LengthOfEmissionPoint)*sin(Phi)*tan(Theta);
  Float_t yq = QuartzWidth*sin(Phi)*tan(thetaquar);
  Float_t yg = GapWidth*sin(Phi)*tan(thetagap);

//  Float_t xtot = x1 + xw + xq + xg;
//  Float_t ytot = y1 + yw + yq + yg;

  Float_t xtot = xEmiss + xw + xq + xg;
  Float_t ytot = yEmiss + yw + yq + yg;

  SetXPointOnCathode(xtot);
  SetYPointOnCathode(ytot);

//  cout << " xtot " << xtot << " ytot " << ytot << endl;

  Float_t DistanceFromEntrance = sqrt(TMath::Power(fPhotonLimitX,2)
				    +TMath::Power(fPhotonLimitY,2)); 

  return DistanceFromEntrance;

}


void AliRICHRecon::FindPhiPoint()
{

  Float_t Xtoentr = GetEntranceX();
  Float_t Ytoentr = GetEntranceY();

  Float_t TrackTheta = GetTrackTheta();
  Float_t TrackPhi = GetTrackPhi();

  Float_t EmissionPoint = GetEmissionPoint();

  Float_t argY = Ytoentr - EmissionPoint*tan(TrackTheta)*sin(TrackPhi);
  Float_t argX = Xtoentr - EmissionPoint*tan(TrackTheta)*cos(TrackPhi);
  Float_t phipad = atan2(argY,argX); 

  SetPhiPoint(phipad);

}

Float_t AliRICHRecon::Cerenkovangle(Float_t n, Float_t beta)
{

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

// Implement Hough response pat. rec. method

  Float_t *HCSspace;

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
  
  Float_t ThetaCerenkov = 0.;
    
  nBin = (int)(0.5+fThetaMax/(fDTheta));
  nCorrBand = (int)(0.5+ fWindowWidth/(2 * fDTheta)); 

  memset ((void *)hcs, 0, fThetaBin*sizeof(float));
  memset ((void *)hcsw, 0, fThetaBin*sizeof(float));

  Int_t NPhotons = GetPhotonsNumber();

  Int_t WeightFlag = 0;

  for (k=0; k< NPhotons; k++) {

    SetPhotonIndex(k);

    angle = GetPhotonEta();

    if(angle == -999.) continue;

    if (angle>=fThetaMin && angle<= fThetaMax) 

      {

	bin = (int)(0.5+angle/(fDTheta));

	bin1= bin-nCorrBand;
	bin2= bin+nCorrBand;

	// calculate weights

	if(kWEIGHT)
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
                WeightFlag = 1;
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

	h1_photons1->Fill(angle);
	h1_photons2->Fill(angle,weight);

	if (bin1<0)    bin1=0;
	if (bin2>nBin) bin2=nBin;
      
	for (j=bin1; j<bin2; j++) 
	  {
	    hcs[j] += 1; 
	    hcsw[j] += weight;
	  }
      }
  }
  
//   if(kDISPLAY)
//     {
//       for(Int_t j=0;j<750;j++)
// 	{
// 	  h1_hcs->Fill(((Float_t)j),hcs[j]);
// 	  h1_hcsw->Fill(((Float_t)j),hcsw[j]);
// 	}
//     }

  if(WeightFlag == 0) 
    {
      HCSspace = hcsw;
    }
  else
    {
      HCSspace = hcs;
      //      cout << " probems with weight...normal procedure adopted " << endl;
    }

  HoughFiltering(HCSspace);

  for (bin=0; bin <nBin; bin++) {
    angle = (bin+0.5) * (fDTheta);
    if (HCSspace[bin] && HCSspace[bin] > etaPeakPos) {
      etaPeakCount = 0;
      etaPeakPos = HCSspace[bin];
      etaPeak[0]=angle;
    }
    else { 
      if (HCSspace[bin] == etaPeakPos) {
	etaPeak[++etaPeakCount] = angle;
      }
    }
  } 

  for (i=0; i<etaPeakCount+1; i++) {
    ThetaCerenkov += etaPeak[i];
  }
  if (etaPeakCount>=0) {
    ThetaCerenkov /= etaPeakCount+1;
    fThetaPeakPos = etaPeakPos;
  }

  SetThetaCerenkov(ThetaCerenkov);
}


void AliRICHRecon::HoughFiltering(float hcs[])
{

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

  Float_t wei = 0.;
  Float_t WeightThetaCerenkov = 0.;

  Int_t NPhotons = GetPhotonsNumber();
  for(Int_t i=0;i<NPhotons;i++)
    {
      SetPhotonIndex(i);

      if(GetPhotonFlag() == 2)
	{
	  Float_t PhotonEta = GetPhotonEta();
	  Float_t PhotonWeight = GetPhotonWeight();
	  WeightThetaCerenkov += PhotonEta*PhotonWeight;
	  wei += PhotonWeight;
	}
    }

  if(wei != 0.) 
    {
      WeightThetaCerenkov /= wei;
    }
  else
    {
      WeightThetaCerenkov = 0.;
    }
  
  SetThetaCerenkov(WeightThetaCerenkov);

  cout << " thetac weighted -> " << WeightThetaCerenkov << endl;
}


void AliRICHRecon::FlagPhotons()
{

  Int_t NPhotonHough = 0;

  Float_t ThetaCerenkov = GetThetaCerenkov();
  if(fDebug) cout << " fThetaCerenkov " << ThetaCerenkov << endl;

  Float_t ThetaDist= ThetaCerenkov - fThetaMin;
  Int_t steps = (Int_t)(ThetaDist / fDTheta);

  Float_t tmin = fThetaMin + (Float_t)(steps - 1)*fDTheta;
  Float_t tmax = fThetaMin + (Float_t)(steps)*fDTheta;
  Float_t tavg = 0.5*(tmin+tmax);

  tmin = tavg - 0.5*fWindowWidth;
  tmax = tavg + 0.5*fWindowWidth;

  if(fDebug) cout << " tmin " << tmin << " tmax " << tmax << endl;
  if(fDebug) cout << " thetac " << ThetaCerenkov << endl;

  //  Int_t CandidatePhotonsNumber = GetCandidatePhotonsNumber();

  Int_t NPhotons = GetPhotonsNumber();

  //  for(Int_t i=0;i<CandidatePhotonsNumber;i++)

  for(Int_t i=0;i<NPhotons;i++)
    {
      SetPhotonIndex(i);

      Float_t PhotonEta = GetPhotonEta();

      if(PhotonEta == -999.) continue;

      if(PhotonEta >= tmin && PhotonEta <= tmax)
	{
	  SetPhotonFlag(2);
	  NPhotonHough++;
	}
    }
  SetHoughPhotons(NPhotonHough);
}

void AliRICHRecon::DrawEvent(Int_t flag)
{

  flag=1; // dummy to be removed...
/*
  Float_t xGraph[3000],yGraph[3000];

  Float_t ThetaCerenkov;

  // Display event...

  gStyle->SetPalette(1,0);

  if(flag == 0) 
    {

      //      Display = new TCanvas("Display","Star Display",0,0,1200,750);      
      
      Display->ToggleEventStatus();
      Display->Modified()
      
      text = new TText(0,0,"");
      text->SetTextFont(61);
      text->SetTextSize(0.03);
      text->SetTextAlign(22);                                                       
      
      Display->Resize();

      h2_disp->Reset();
      
      for(Int_t j=1;j<=nPixels;j++)
	{
	  Float_t xpad = fPixels_localX[j-1];
	  Float_t ypad = fPixels_localY[j-1];
	  h2_disp->Fill(xpad,ypad,fPixels_charge[j-1]);
	}

      h2_disp->SetMaximum(200);
      //      h2_disp->SetMaximum(1);
      h2_disp->SetStats(0);
      h2_disp->Draw("colz");
      
      for(Int_t i=0; i<nRichPrimaries;i++)
	
	{
	  
	  TrackPoints = new TMarker(fRichPrimaries_localPadX[i],
				    fRichPrimaries_localPadY[i],3);

	  TrackPoints->SetMarkerSize(1.5);
	  
	  Float_t pmod = sqrt(fRichPrimaries_localPadPx[i] * fRichPrimaries_localPadPx[i] +
			      fRichPrimaries_localPadPy[i] * fRichPrimaries_localPadPy[i] +
			      fRichPrimaries_localPadPz[i] * fRichPrimaries_localPadPz[i]); 
	  
	  if(pmod < 1) TrackPoints->SetMarkerColor(kBlue);
	  if(pmod > 1 && pmod < 2) TrackPoints->SetMarkerColor(kGreen);
	  if(pmod > 2) TrackPoints->SetMarkerColor(kRed);
	  
	  TrackPoints->Draw();

	  line = new TLine(-0.13,-42.,-0.13,42.);
	  line->Draw();
	  line = new TLine(0.13,-42.,0.13,42.);
	  line->Draw();
	  line = new TLine(-64.,-0.13,64.,-0.13);
	  line->Draw();
	  line = new TLine(-64.,0.13,64.,0.13);
	  line->Draw();                        

	}
      
      return;

    }
  
  //

  // Draw rings...

  //

  //  Float_t Xtoentr = GetEntranceX();
  //  Float_t Ytoentr = GetEntranceY();

  //  Float_t pmod = GetTrackMomentum();
  //  Float_t TrackTheta = GetTrackTheta();
  //  Float_t TrackPhi = GetTrackPhi();

  SetPhotonEnergy(6.85);
  SetFreonRefractiveIndex();

  SetEmissionPoint(RadiatorWidth/2.);

  ThetaCerenkov = GetThetaCerenkov();

  if (ThetaCerenkov == 999.) return;

  Int_t nPointsToDraw = 0;

  for(Int_t i=0;i<99;i++)
    {
      Float_t phpad = 2*TMath::Pi()*i/99;
      SetThetaPhotonInTRS(ThetaCerenkov);
      SetPhiPhotonInTRS(phpad);
      FindPhotonAnglesInDRS();
      Float_t Radius = FromEmissionToCathode();
      
      if (Radius == 999.) continue;
      
      Float_t ShiftX = GetShiftX();
      Float_t ShiftY = GetShiftY();
      
      Float_t XPointRing = GetXPointOnCathode() + ShiftX;
      Float_t YPointRing = GetYPointOnCathode() + ShiftY;
      
      SetDetectorWhereX(XPointRing);
      SetDetectorWhereY(YPointRing);
      
      Int_t Zone = CheckDetectorAcceptance();
      
      if (Zone != 0) 
	{
	  FindIntersectionWithDetector();
	  xGraph[nPointsToDraw] = GetIntersectionX();
	  yGraph[nPointsToDraw] = GetIntersectionY();
	  nPointsToDraw++;
	}
      else
	{
	  xGraph[nPointsToDraw] = GetXPointOnCathode() + GetShiftX();
	  yGraph[nPointsToDraw] = GetYPointOnCathode() + GetShiftY();
	  nPointsToDraw++;
	}
    }
  
  xGraph[nPointsToDraw] = xGraph[0];  
  yGraph[nPointsToDraw] = yGraph[0];  

  poll = new TPolyLine(nPointsToDraw+1,xGraph,yGraph);
  poll->SetLineColor(2);
  poll->SetLineWidth(3);

  Display->Update();

  //  waiting();
  poll->Draw();

  for(Int_t j=0;j<nHits;j++)
    {
      
      Float_t xhit = fHits_localX[j];
      Float_t yhit = fHits_localY[j];

      SetPhotonIndex(j);
      Int_t FlagPhoton = GetPhotonFlag();

//       if(FlagPhoton >= 1) 
// 	{

// 	  Photon = new TMarker(xhit,yhit,4);
// 	  Photon->SetMarkerSize(1.5);
// 	  Photon->Draw("same");

// 	}


      if(FlagPhoton == 2) 
	{
	  
	  PhotonAcc = new TMarker(xhit,yhit,30);
	  PhotonAcc->SetMarkerSize(1.5);
	  PhotonAcc->SetMarkerColor(50);
	  PhotonAcc->Draw("same");
	  
	}
    }  

  Display->Update();

//   waiting();
//   h1_photons->Draw();
//   Display->Update();

//   waiting();
//   h1_photacc->Draw();
//   Display->Update();

//   waiting();

//   Display->Update();

//   h1_photons->Reset();
//   h1_photacc->Reset();

*/
}

Float_t  AliRICHRecon::FindMassOfParticle()
{

  Float_t pmod = GetTrackMomentum();

  SetPhotonEnergy(6.85);
  SetFreonRefractiveIndex();

  Float_t ThetaCerenkov = GetThetaCerenkov();
  FindBetaFromTheta(ThetaCerenkov);

  Double_t beta = (Double_t)(GetBetaOfParticle());
  Double_t den = 1. - beta*beta;
  if(den<=0.) return 999.;

  Double_t gamma = 1./TMath::Sqrt(den);

  Float_t mass = pmod/(beta*(Float_t)gamma);

  return mass;
}


void AliRICHRecon::FillHistograms()
{

  Float_t FittedTrackTheta, FittedTrackPhi;

  Float_t ThetaCerenkov    = GetThetaCerenkov();
  if(ThetaCerenkov == 999.) return;

  Float_t VertZ = GetEventVertexZ();

  Float_t TrackTheta = GetTrackTheta();
  Float_t TrackPhi   = GetTrackPhi();
  Float_t pmod       = GetTrackMomentum();
  Float_t pt         = GetTrackPt();
  Float_t TrackEta   = GetTrackEta();
  Int_t q            = GetTrackCharge();
  Float_t TPCLastZ   = GetTrackTPCLastZ(); 
  Float_t MinDist    = GetMinDist(); 

  FittedTrackTheta = GetFittedTrackTheta();
  FittedTrackPhi   = GetFittedTrackPhi();
  Int_t FittedNPhotonHough = GetFittedHoughPhotons();
  
  if(fDebug)
    {
      cout << " p " << pmod  << " ThetaC " << ThetaCerenkov 
	   << " rings " << NRings << endl;
    }

  Int_t NPhotonHough     = GetHoughPhotons();
  Float_t NPhotonHoughNorm = GetHoughPhotonsNorm();
  Float_t InRing = GetPortionOfRing();

  Float_t MassOfParticle = FindMassOfParticle();

  Float_t HoughArea = GetHoughArea();
  Float_t Multiplicity = GetEventMultiplicity();

//  cout << " area " << HoughArea << " mult " << Multiplicity << endl;

  Float_t var[20];

//  var[0] = (Float_t)runID; 
//  var[1] = (Float_t)evID;
  var[0] = 0; 
  var[1] = 0;
  var[2] = VertZ;
  var[3] = pmod;
  var[4] = pt;
  var[5] = TrackEta;
  var[6] = TrackTheta;
  var[7] = TrackPhi;
  var[8] = FittedTrackTheta;
  var[9] = FittedTrackPhi;
  var[10] = q;
  var[11] = ThetaCerenkov;
  var[12] = (Float_t)NPhotonHough;
  var[13] = (Float_t)FittedNPhotonHough;
  var[14] = InRing;
  var[15] = MassOfParticle;
  var[16] = HoughArea;
  var[17] = Multiplicity;
  var[18] = TPCLastZ;
  var[19] = MinDist;

  hn->Fill(var);

  h1_mass->Fill(MassOfParticle);
  h2_mvsp->Fill(pmod,MassOfParticle);
  h2_mvst->Fill(ThetaCerenkov,MassOfParticle);

  FittedTrackTheta = GetFittedTrackTheta();
  FittedTrackPhi = GetFittedTrackPhi();

  Float_t DiffTheta = FittedTrackTheta - TrackTheta;
  Float_t DiffPhi = FittedTrackPhi - TrackPhi;

  h1_diffTrackTheta -> Fill(DiffTheta);
  h1_diffTrackPhi -> Fill(DiffPhi);

  if(ThetaCerenkov > 0.505 && ThetaCerenkov < 0.605) 
    {
      SetPhotonEnergy(6.85);
      SetFreonRefractiveIndex();

      Float_t pmom = GetTrackMomentum();
      Float_t beta = 1./(cos(ThetaCerenkov)*GetFreonRefractiveIndex());
      Float_t gamma = 1./sqrt(1.-beta*beta);

      Float_t pmomnew = 0.93828*beta*gamma;
      Float_t deltap = pmomnew - pmom;
      h1_deltap->Fill(deltap);
      Float_t deltapop = deltap/pmom;
      h1_deltapop->Fill(deltapop);

      h1_nprotons->Fill((Float_t)NPhotonHoughNorm);
    }

  if(q > 0)
    {
      h2_tvsppos->Fill(pmod,ThetaCerenkov);
      hp_1pos->Fill(ThetaCerenkov,(Float_t)NPhotonHough);
      hp_1posnorm->Fill(ThetaCerenkov,(Float_t)NPhotonHoughNorm);
      h2_1pos->Fill(pmod,(Float_t)NPhotonHough);
      h2_1posnorm->Fill(pmod,(Float_t)NPhotonHoughNorm);
      h1_houghpos->Fill(ThetaCerenkov);
    }
else
  {
      h2_tvspneg->Fill(pmod,ThetaCerenkov);
      hp_1neg->Fill(ThetaCerenkov,(Float_t)NPhotonHough);
      hp_1negnorm->Fill(ThetaCerenkov,(Float_t)NPhotonHoughNorm);
      h2_1neg->Fill(pmod,(Float_t)NPhotonHough);
      h2_1negnorm->Fill(pmod,(Float_t)NPhotonHoughNorm);
      h1_houghneg->Fill(ThetaCerenkov);
  }

  Int_t NPhotons = GetPhotonsNumber();

  for (Int_t j=0; j < NPhotons;j++)

    {
      SetPhotonIndex(j);

      Float_t eta = GetPhotonEta();

      if(GetPhotonFlag() == 2) 
	{
	  h1_photacc->Fill(eta);
	  Float_t photaccspread = eta - ThetaCerenkov;
	  h1_photaccspread->Fill(photaccspread);
	}

    }
}

void AliRICHRecon::Minimization()
{

  Double_t arglist;
  Int_t ierflag = 0;

  static Double_t vstart[2];
  static Double_t lower[2], upper[2];
  static Double_t step[2]={0.001,0.001};

  Double_t TrackThetaNew,TrackPhiNew;
  TString chname;
  Double_t eps, b1, b2;
  Int_t ierflg;

  gMyMinuit = new TMinuit(2);
  gMyMinuit->SetObjectFit((TObject *)this);
  gMyMinuit->SetFCN(fcnrecon);
  gMyMinuit->mninit(5,10,7);

  vstart[0] = (Double_t)GetTrackTheta();
  vstart[1] = (Double_t)GetTrackPhi();

  lower[0] = vstart[0] - 0.03;
  if(lower[0] < 0) lower[0] = 0.;
  upper[0] = vstart[0] + 0.03;
  lower[1] = vstart[1] - 0.03;
  upper[1] = vstart[1] + 0.03;


  gMyMinuit->mnparm(0,"theta",vstart[0],step[0],lower[0],upper[0],ierflag);
  gMyMinuit->mnparm(1," phi ",vstart[1],step[1],lower[1],upper[1],ierflag);

  arglist = -1;

  //  gMyMinuit->FixParameter(0);

  gMyMinuit->SetPrintLevel(-1);
//  gMyMinuit->mnexcm("SET PRI",&arglist, 1, ierflag);
  gMyMinuit->mnexcm("SET NOGR",&arglist, 1, ierflag);
  gMyMinuit->mnexcm("SET NOW",&arglist, 1, ierflag);
  arglist = 1;
  gMyMinuit->mnexcm("SET ERR", &arglist, 1,ierflg);
  arglist = -1;

  //  gMyMinuit->mnscan();

//  gMyMinuit->mnexcm("SIMPLEX",&arglist, 0, ierflag);
  gMyMinuit->mnexcm("MIGRAD",&arglist, 0, ierflag);
  gMyMinuit->mnexcm("EXIT" ,&arglist, 0, ierflag);
  
  gMyMinuit->mnpout(0,chname, TrackThetaNew, eps , b1, b2, ierflg);
  gMyMinuit->mnpout(1,chname, TrackPhiNew, eps , b1, b2, ierflg);

  //values after the fit...
  SetFittedTrackTheta((Float_t)TrackThetaNew);
  SetFittedTrackPhi((Float_t)TrackPhiNew);

  delete gMyMinuit;

}

void AliRICHRecon::EstimationOfTheta()
{

  Int_t NPhotons = 0;

  Float_t ShiftX = GetShiftX();
  Float_t ShiftY = GetShiftY();

  Float_t *CandidatePhotonX = GetCandidatePhotonX();
  Float_t *CandidatePhotonY = GetCandidatePhotonY();

  Int_t NPhotonsCandidates = GetCandidatePhotonsNumber();

  //  cout << "MINIM: Nphotons " << NPhotonsCandidates << endl;

  for (Int_t j=0; j < NPhotonsCandidates; j++)
    {

      SetPhotonIndex(j);

      if(!GetPhotonFlag()) continue;

      Float_t Xtoentr = CandidatePhotonX[j] - ShiftX;
      Float_t Ytoentr = CandidatePhotonY[j] - ShiftY;

      SetEntranceX(Xtoentr);
      SetEntranceY(Ytoentr);

      FindPhiPoint();

      FindThetaPhotonCerenkov();

      Float_t ThetaPhotonCerenkov = GetThetaPhotonCerenkov();

      //      cout << " ACCEPTED!!! " << ThetaPhotonCerenkov << endl;

      SetPhotonEta(ThetaPhotonCerenkov);

      NPhotons++;

    }

  Float_t xmean = 0.;
  Float_t x2mean = 0.;
  Int_t nev = 0;

  for (Int_t j=0; j < NPhotonsCandidates;j++)
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

  Float_t RMS = sqrt(x2mean - xmean*xmean);

  //  cout << " RMS " << RMS;

  SetEstimationOfTheta(xmean);
  SetEstimationOfThetaRMS(RMS);
}

void fcnrecon(Int_t& /*npar*/, Double_t* /*gin*/, Double_t &f, Double_t *par, Int_t iflag)
{
  AliRICHRecon *gMyRecon = (AliRICHRecon*)gMyMinuit->GetObjectFit();

  Float_t p0 = (Float_t)par[0];
  Float_t p1 = (Float_t)par[1];

  gMyRecon->SetTrackTheta(p0);
  gMyRecon->SetTrackPhi(p1);

  gMyRecon->EstimationOfTheta();
  Float_t RMS = gMyRecon->GetEstimationOfThetaRMS();

  Int_t HoughPhotons = gMyRecon->GetHoughPhotons();


  f = (Double_t)(1000*RMS/(Float_t)HoughPhotons);

  if(fDebug) cout << "   f   " << f
		  << " theta " << par[0] << " phi " << par[1] 
                  << " HoughPhotons " << HoughPhotons << endl;
  
  if(fDebug&&iflag == 3)
    {
            cout << " --- end convergence...summary --- " << endl;
            cout << " theta " << par[0] << endl;
            cout << "  phi  " << par[1] << endl;
    }
}

void AliRICHRecon::waiting()
{
  if(!kDISPLAY) return;
  cout << " Press any key to continue...";

//  gSystem->ProcessEvents();
  getchar(); 

  cout << endl;

  return;
}

/*
void ~AliRICHRecon()
{
}
*/
