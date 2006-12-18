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

/* $Id$ */

// Class to generate correlated Heavy Flavor hadron pairs (one or several pairs
// per event) using paramtrized kinematics of quark pairs from some generator
// and quark fragmentation functions.
// Is a generalisation of AliGenParam class for correlated pairs of hadrons.
// In this version quark pairs and fragmentation functions are obtained from
// Pythia6.124 using 100K events generated with kCharmppMNRwmi & kBeautyppMNRwmi 
// in pp collisions at 14 TeV.
// Decays are performed by Pythia. Used AliRoot version: v4-04-Release
// Author: S. Grigoryan, LPC Clermont-Fd & YerPhI, Smbat.Grigoryan@cern.ch
//
//-------------------------------------------------------------------------
// How it works (for the given flavor):
//
// 1) Reads QQbar kinematical grid from the Input file and generates
// quark pairs according to the weights of the cells.
// It is a 5D grid in y1,y2,pt1,pt2 and deltaphi, with occupancy weights
// of the cells obtained from Pythia (see details in GetQuarkPair).
// 2) Reads "soft" and "hard" fragmentation functions (12 2D-histograms each,
// for 12 pt bins) from the Input file, applies to quarks and produces hadrons
// (only lower states, with proportions of species obtained from Pythia).
// Fragmentation functions are the same for all hadron species and depend
// on 2 variables - light cone energy-momentum fractions:
//     z1=(E_H + Pz_H)/(E_Q + Pz_Q),  z2=(E_H - Pz_H)/(E_Q - Pz_Q).
// "soft" & "hard" FFs correspond to "slower" & "faster" quark of a pair 
// (see details in GetHadronPair).
// 3) Decays the hadrons and saves all the particles in the event stack in the
// following order: HF hadron from Q, then its decay products, then HF hadron
// from Qbar, then its decay productes, then next HF hadon pair (if any) 
// in the same way, etc... 
// 4) It is fast, e.g., generates the same number of events with a beauty pair 
//  ~15 times faster than AliGenPythia with kBeautyppMNRwmi (w/o tracking)
//
// An Input file for each quark flavor is included in EVGEN/dataCorrHF/
// One can use also user-defined Input files.
//
// More details could be found in my presentation at DiMuonNet Workshop, Dec 2006: 
// http://www-dapnia.cea.fr/Sphn/Alice/DiMuonNet
// and will be published in an Internal Note. 
//
//-------------------------------------------------------------------------
// How to use it:
//
// add the following typical lines in Config.C
/*
  if (!strcmp(option,"corr")) {
    // Example for correlated charm or beauty hadron pair production 

    // AliGenCorrHF *gener = new AliGenCorrHF(1, 4);  // for charm, 1 pair per event
    AliGenCorrHF *gener = new AliGenCorrHF(1, 5);  // for beauty, 1 pair per event

    gener->SetMomentumRange(0,9999);
    gener->SetCutOnChild(0);        // 1/0 means cuts on children enable/disable
    gener->SetChildThetaRange(171.0,178.0);
    gener->SetOrigin(0,0,0);          //vertex position    
    gener->SetSigma(0,0,0);           //Sigma in (X,Y,Z) (cm) on IP position
    gener->SetForceDecay(kSemiMuonic);
    gener->SetTrackingFlag(0);
    gener->Init();
}
*/
// and in aliroot do e.g. gAlice->Run(10,"Config.C") to produce 10 events.
// One can include AliGenCorrHF in an AliGenCocktail generator.
//--------------------------------------------------------------------------

#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TMath.h>
#include <TRandom.h>
#include <TROOT.h>
#include <TLorentzVector.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TDatabasePDG.h>
#include <TVirtualMC.h>
#include <TCanvas.h>
#include <Riostream.h>

#include "AliGenCorrHF.h"
#include "AliLog.h"
#include "AliConst.h"
#include "AliDecayer.h"
#include "AliMC.h"
#include "AliRun.h"

ClassImp(AliGenCorrHF)

  //Begin_Html
  /*
    <img src="picts/AliGenCorrHF.gif">
  */
  //End_Html

Double_t AliGenCorrHF::fgdph[19] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180};
Double_t AliGenCorrHF::fgy[31] = {-10,-7, -6.5, -6, -5.5, -5, -4.5, -4, -3.5, -3, -2.5, -2,- 1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 10};
Double_t AliGenCorrHF::fgpt[33] = {0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.6, 7.2, 7.8, 8.4, 9, 9.6, 10.3, 11.1, 12, 13.1, 14.3, 15.6, 17.1, 19, 21, 24, 28, 35, 50, 100};
Int_t AliGenCorrHF::fgnptbins = 12;
Double_t AliGenCorrHF::fgptbmin[12] = {0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 9};
Double_t AliGenCorrHF::fgptbmax[12] = {0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 9, 100};

Double_t* AliGenCorrHF::fgIntegral = 0;

//____________________________________________________________
    AliGenCorrHF::AliGenCorrHF():
	fFileName(0),
	fFile(0),
	fQuark(0),
	fBias(0.),
	fTrials(0),
	fDecayer(0)
{
// Default constructor
}

//____________________________________________________________
AliGenCorrHF::AliGenCorrHF(Int_t npart, Int_t param):
    AliGenMC(npart),
    fFileName(0),
    fFile(0),
    fQuark(param),
    fBias(0.),
    fTrials(0),
    //    fDecayer(new AliDecayerPythia())
    fDecayer(0)
{
// Constructor using number of particles, quark type & default InputFile
//
    if (fQuark != 5) fQuark = 4;
    fFileName = "$ALICE_ROOT/EVGEN/dataCorrHF/CharmppMNRwmiCorr100K.root";
    if (fQuark == 5) fFileName = "$ALICE_ROOT/EVGEN/dataCorrHF/BeautyppMNRwmiCorr100K.root";

    fName = "Default";
    fTitle= "Generator for correlated pairs of HF hadrons";
      
    fChildSelect.Set(5);
    for (Int_t i=0; i<5; i++) fChildSelect[i]=0;
    SetForceDecay();
    SetCutOnChild();
    SetChildMomentumRange();
    SetChildPtRange();
    SetChildPhiRange();
    SetChildThetaRange(); 
}

//___________________________________________________________________
AliGenCorrHF::AliGenCorrHF(char* tname, Int_t npart, Int_t param):
    AliGenMC(npart),
    fFileName(tname),
    fFile(0),
    fQuark(param),
    fBias(0.),
    fTrials(0),
    //    fDecayer(new AliDecayerPythia())
    fDecayer(0)
{
// Constructor using number of particles, quark type & user-defined InputFile
//
    if (fQuark != 5) fQuark = 4;
    fName = "UserDefined";
    fTitle= "Generator for correlated pairs of HF hadrons";
      
    fChildSelect.Set(5);
    for (Int_t i=0; i<5; i++) fChildSelect[i]=0;
    SetForceDecay();
    SetCutOnChild();
    SetChildMomentumRange();
    SetChildPtRange();
    SetChildPhiRange();
    SetChildThetaRange(); 
}

//____________________________________________________________
AliGenCorrHF::AliGenCorrHF(const AliGenCorrHF & CorrHF)
    :AliGenMC(CorrHF),	
     fFileName(0),
     fFile(0),
     fQuark(0),
     fBias(0.),
     fTrials(0),
     fDecayer(0)
{
// Copy constructor
    CorrHF.Copy(*this);
}

//____________________________________________________________
AliGenCorrHF::~AliGenCorrHF()
{
// Destructor
  delete fFile;
}

//____________________________________________________________
void AliGenCorrHF::Init()
{
// Initialisation

  AliInfo(Form(" QQbar kinematics and fragm. functions from:  %s",fFileName.Data() )); 
    fFile = TFile::Open(fFileName.Data());
    if(!fFile->IsOpen()){
      AliError(Form("Could not open file %s",fFileName.Data() ));
    }

    ComputeIntegral(fFile);

    fParentWeight = 1./fNpart;   // fNpart is number of HF-hadron pairs

// particle decay related initialization

    if (gMC) fDecayer = gMC->GetDecayer();
    fDecayer->SetForceDecay(fForceDecay);
    fDecayer->Init();

//
    AliGenMC::Init();
}

//____________________________________________________________
void AliGenCorrHF::Generate()
{
//
// Generate fNpart of correlated HF hadron pairs per event
// in the the desired theta and momentum windows (phi = 0 - 2pi). 
// Gaussian smearing on the vertex is done if selected. 
// The decay of heavy hadrons is done using lujet, 
//    and the childern particle are tracked by GEANT
// However, light mesons are directly tracked by GEANT 
// setting fForceDecay = nodecay (SetForceDecay(nodecay)) 
//


  Float_t polar[3]= {0,0,0};  // Polarisation of the parent particle (for GEANT tracking)
  Float_t origin0[3];         // Origin of the generated parent particle (for GEANT tracking)
  Float_t pt, pl, ptot;       // Transverse, logitudinal and total momenta of the parent particle
  Float_t phi, theta;         // Phi and theta spherical angles of the parent particle momentum
  Float_t p[3], pc[3], och[3];// Momentum, polarisation and origin of the children particles from lujet


  Double_t dphi=0, ptq[2], yq[2], pth[2], plh[2], ph[2], phih[2];
  Int_t i, j, ipair, ihadron[2];
  for (i=0; i<2; i++) { 
    ptq[i] =0; 
    yq[i]  =0; 
    pth[i] =0; 
    plh[i] =0; 
    ihadron[i] =0; 
  }

  static TClonesArray *particles;
  //
  if(!particles) particles = new TClonesArray("TParticle",1000);
  
  TDatabasePDG* pDataBase = TDatabasePDG::Instance();
 
// Calculating vertex position per event
  for (j=0;j<3;j++) origin0[j]=fOrigin[j];
  if(fVertexSmear==kPerEvent) {
      Vertex();
      for (j=0;j<3;j++) origin0[j]=fVertex[j];
  }
  
  Float_t wgtp, wgtch, random[6];
  Int_t ipap = 0;
  Int_t nt   = 0;

// Generating fNpart HF-hadron pairs per event
  while(ipap<fNpart) {

    while(1) {

      GetQuarkPair(fFile, fgIntegral, yq[0], yq[1], ptq[0], ptq[1], dphi);

      GetHadronPair(fFile, fQuark, yq[0], yq[1], ptq[0], ptq[1], ihadron[0], ihadron[1], plh[0], plh[1], pth[0], pth[1]);

// Here we assume that  |phi_H1 - phi_H2| = |phi_Q1 - phi_Q2| = dphi
// which is a good approximation for heavy flavors in Pythia

      /* // doesn't work if PhiMax < k2PI or PhiMin > 0, since dphi = 0 - 180 
      phih[0] = fPhiMin + Rndm()*(fPhiMax-fPhiMin);
      phih[1] = phih[0] + dphi*kDegrad;
      if (phih[0] > fPhiMax/2.) phih[1] = phih[0] - dphi*kDegrad;
      */
      phih[0] = Rndm()*k2PI;
      phih[1] = phih[0] + dphi*kDegrad; 
      if (phih[0] > TMath::Pi()) phih[1] = phih[0] - dphi*kDegrad;

// Cut on theta
      theta=TMath::ATan2(pth[0],plh[0]);
      if(theta<fThetaMin || theta>fThetaMax) continue;
      theta=TMath::ATan2(pth[1],plh[1]);
      if(theta<fThetaMin || theta>fThetaMax) continue;

// Cut on momentum
      ph[0]=TMath::Sqrt(pth[0]*pth[0]+plh[0]*plh[0]);
      if (ph[0]<fPMin || ph[0]>fPMax) continue;
      ph[1]=TMath::Sqrt(pth[1]*pth[1]+plh[1]*plh[1]);
      if (ph[1]<fPMin || ph[1]>fPMax) continue;

// Common origin for particles of the HF-hadron pair
      if(fVertexSmear==kPerTrack) {
	  Rndm(random,6);
	  for (j=0;j<3;j++) {
	      origin0[j]=
		fOrigin[j]+fOsigma[j]*TMath::Cos(2*random[2*j]*TMath::Pi())*
		TMath::Sqrt(-2*TMath::Log(random[2*j+1]));
	  }
      }

      Int_t np1=0, kf1[100], select1[100], iparent1[100], trackIt1[100];
      Float_t  wgtch1=0, p1[3], pc1[100][3], och1[100][3];

      for (j=0; j<3; j++)   p1[j] = 0;
      for (i=0; i<100; i++) {
	kf1[i]      = 0;  
	select1[i]  = 0;  
	iparent1[i] = 0;  
	trackIt1[i] = 0;
	for (j=0; j<3; j++) {
	  pc1[i][j]  = 0;
	  och1[i][j] = 0;
	}
      }

//
// Loop over particles of the HF-hadron pair
      Int_t nhadron = 0;
      for (ipair=0;ipair<2;ipair++) {
	phi  = phih[ipair];
	pl   = plh[ipair];
	pt   = pth[ipair];
	ptot = ph[ipair];
//
// particle type 
	  Int_t iPart = ihadron[ipair];
	  Float_t am  = pDataBase->GetParticle(iPart)->Mass();
          fChildWeight=(fDecayer->GetPartialBranchingRatio(iPart))*fParentWeight;

          wgtp  = fParentWeight;
	  wgtch = fChildWeight;

//
	  p[0]=pt*TMath::Cos(phi);
	  p[1]=pt*TMath::Sin(phi);
	  p[2]=pl;
	  
// Looking at fForceDecay : 
// if fForceDecay != none Primary particle decays using 
// AliPythia and children are tracked by GEANT
//
// if fForceDecay == none Primary particle is tracked by GEANT 
// (In the latest, make sure that GEANT actually does all the decays you want)	  
//

	  if (fForceDecay != kNoDecay) {
// Using lujet to decay particle
	      Float_t energy=TMath::Sqrt(ptot*ptot+am*am);
	      TLorentzVector pmom(p[0], p[1], p[2], energy);
   	      fDecayer->Decay(iPart,&pmom);
//
// select decay particles
	      Int_t np=fDecayer->ImportParticles(particles);

//  Selecting  GeometryAcceptance for particles fPdgCodeParticleforAcceptanceCut;
	      if (fGeometryAcceptance)
		if (!CheckAcceptanceGeometry(np,particles)) break;
	      Int_t ncsel=0;
	      Int_t* pParent    = new Int_t[np];
	      Int_t* pSelected  = new Int_t[np];
	      Int_t* trackIt    = new Int_t[np];

	      for (i=0; i<np; i++) {
		  pSelected[i] =  0;
		  pParent[i]   = -1;
	      }
	      
	      if (np >1) {
		  TParticle* iparticle =  (TParticle *) particles->At(0);
		  for (i=1; i<np; i++) {
		      trackIt[i] = 1;
		      iparticle = (TParticle *) particles->At(i);
		      Int_t kf = iparticle->GetPdgCode();
		      Int_t ks = iparticle->GetStatusCode();

// particles with long life-time (c tau > .3 mum)
		      if (ks != 1) { 
			  Double_t lifeTime = fDecayer->GetLifetime(kf);
			  if (lifeTime <= (Double_t) fMaxLifeTime) {
			      trackIt[i]     = 0;
			      pSelected[i]   = 1;
			  }
		      } // ks==1 ?
//
// children, discard neutrinos
		      if (TMath::Abs(kf) == 12 || TMath::Abs(kf) == 14) continue;
		      if (trackIt[i])
		      {
			  if (fCutOnChild) {
			      pc[0]=iparticle->Px();
			      pc[1]=iparticle->Py();
			      pc[2]=iparticle->Pz();
			      Bool_t  childok = KinematicSelection(iparticle, 1);
			      if(childok) {
				  pSelected[i]  = 1;
				  ncsel++;
			      } else {
				  ncsel=-1;
				  break;
			      } // child kine cuts
			  } else {
			      pSelected[i]  = 1;
			      ncsel++;
			  } // if child selection
		      } // select muon
		  } // decay particle loop
	      } // if decay products
	      
	      Int_t iparent;
	      if ((fCutOnChild && ncsel >0) || !fCutOnChild){

         	  nhadron++;
//
// Parents and Decay Products
		  if (ipair == 0) {
		      np1    = np;
		      wgtch1 = wgtch;
		      p1[0] = p[0]; p1[1] = p[1]; p1[2] = p[2];
		  } else {
		      ipap++;
		      PushTrack(0, -1, ihadron[0], p1, origin0, polar, 0, 
				kPPrimary, nt, wgtp);
		      KeepTrack(nt);
		      for (i = 1; i < np1; i++) {
			if (select1[i]) {
			  for (j=0; j<3; j++) {
			      och[j] = och1[i][j];
			      pc[j]  = pc1[i][j];
			  }
			  PushTrack(fTrackIt*trackIt1[i], iparent1[i], kf1[i], pc, och,
				    polar, 0, kPDecay, nt, wgtch1);
			  KeepTrack(nt);
			  }
		      }
		      PushTrack(0, -1, iPart, p, origin0, polar, 0, kPPrimary, nt, wgtp);
		      KeepTrack(nt);
		  } 
		  pParent[0] = nt;
//
// Decay Products
		  Int_t ntcount = 0;
		  for (i = 1; i < np; i++) {
		    if (pSelected[i]) {
			  TParticle* iparticle = (TParticle *) particles->At(i);
			  Int_t kf  = iparticle->GetPdgCode();
			  Int_t ipa = iparticle->GetFirstMother()-1;

			  och[0] = origin0[0]+iparticle->Vx()/10;
			  och[1] = origin0[1]+iparticle->Vy()/10;
			  och[2] = origin0[2]+iparticle->Vz()/10;
			  pc[0]  = iparticle->Px();
			  pc[1]  = iparticle->Py();
			  pc[2]  = iparticle->Pz();
			  
			  if (ipa > -1) {
			      iparent = pParent[ipa];
			  } else {
			      iparent = -1;
			  }

			  if (ipair == 0) {
			      kf1[i]      = kf;  
			      select1[i]  = pSelected[i];
			      iparent1[i] = iparent;
			      trackIt1[i] = trackIt[i];
			      for (j=0; j<3; j++) {
				  och1[i][j] = och[j];
				  pc1[i][j]  = pc[j];
			      }
			      ntcount++;
			  } else {
			    PushTrack(fTrackIt*trackIt[i], iparent, kf, pc, och,
				      polar, 0, kPDecay, nt, wgtch);
			    KeepTrack(nt); 
			  } 
			  pParent[i] = nt + ntcount;
		    } // Selected
		  } // Particle loop 
	      }  // Decays by Lujet
	      particles->Clear();
	      if (pParent)    delete[] pParent;
	      if (pSelected)  delete[] pSelected;	   
	      if (trackIt)    delete[] trackIt;
	  } // kinematic selection
	  else  // nodecay option, so parent will be tracked by GEANT
	  {
	    nhadron++;
	    if (ipair == 0) {
	        p1[0] = p[0]; p1[1] = p[1]; p1[2] = p[2];
	    } else {
	        ipap++;
		gAlice->GetMCApp()->
		PushTrack(fTrackIt,-1,ihadron[0],p1,origin0,polar,0,kPPrimary,nt,wgtp);
		gAlice->GetMCApp()->
		PushTrack(fTrackIt,-1,iPart,p,origin0,polar,0,kPPrimary,nt,wgtp);
	    }
	  }
	  if (nhadron == 0) break;
      } // ipair loop
      if (nhadron != 2) continue;
      break;
    } // while(1)
    nt++;
  } // while(ipa<fNpart) --> event loop

  SetHighWaterMark(nt);
}

//____________________________________________________________________________________
AliGenCorrHF& AliGenCorrHF::operator=(const  AliGenCorrHF& rhs)
{
// Assignment operator
    rhs.Copy(*this);
    return *this;
}

//____________________________________________________________________________________
Int_t AliGenCorrHF::IpCharm(TRandom* ran)
{  
// Composition of lower state charm hadrons, containing a c-quark
    Float_t random;
    Int_t ip;            // +- 411,421,431,4122,4132,4232,4332
    random = ran->Rndm();
//  Rates from Pythia6.214 using 100Kevents with kPyCharmppMNRwmi at 14 TeV.   
  
    if (random < 0.6027) {                       
        ip=421;
    } else if (random < 0.7962) {
        ip=411;
    } else if (random < 0.9127) {
        ip=431;        
    } else if (random < 0.9899) {
        ip=4122;        
    } else if (random < 0.9948) {
        ip=4132;        
    } else if (random < 0.9999) {
        ip=4232;        
    } else {
        ip=4332;
    }
    
    return ip;
}

Int_t AliGenCorrHF::IpBeauty(TRandom* ran)
{  
// Composition of lower state beauty hadrons, containing a b-quark
    Float_t random;
    Int_t ip;            // +- 511,521,531,5122,5132,5232,5332
    random = ran->Rndm(); 
//  Rates from Pythia6.214 using 100Kevents with kPyBeautyppMNRwmi at 14 TeV.   
                        // B-Bbar mixing will be done by Pythia at the decay point
 if (random < 0.3965) {                       
        ip=-511;
    } else if (random < 0.7930) {
        ip=-521;
    } else if (random < 0.9112) {
        ip=-531;        
    } else if (random < 0.9887) {
        ip=5122;        
    } else if (random < 0.9943) {
        ip=5132;        
    } else if (random < 0.9999) {
        ip=5232;        
    } else {
        ip=5332;
    }
    
  return ip;
}

//____________________________________________________________________________________
Double_t AliGenCorrHF::ComputeIntegral(TFile* fG)       // needed by GetQuarkPair
{
   // Read QQbar kinematical 5D grid's cell occupancy weights
   Int_t* cell = new Int_t[6];           // cell[6]={wght,iy1,iy2,ipt1,ipt2,idph}
   TTree* tG = (TTree*) fG->Get("tGqq");
   tG->GetBranch("cell")->SetAddress(cell);
   Int_t nbins = tG->GetEntries();

   //   delete previously computed integral (if any)
   if(fgIntegral) delete [] fgIntegral;

   fgIntegral = new Double_t[nbins+1];
   fgIntegral[0] = 0;
   Int_t bin;
   for(bin=0;bin<nbins;bin++) {
     tG->GetEvent(bin);
     fgIntegral[bin+1] = fgIntegral[bin] + cell[0];
   }
   //   Normalize integral to 1
   if (fgIntegral[nbins] == 0 ) {
      return 0;
   }
   for (bin=1;bin<=nbins;bin++)  fgIntegral[bin] /= fgIntegral[nbins];

   return fgIntegral[nbins];
}

//____________________________________________________________________________________
void AliGenCorrHF::GetQuarkPair(TFile* fG, Double_t* fInt, Double_t &y1, Double_t &y2, Double_t &pt1, Double_t &pt2, Double_t &dphi)              
                                 // modification of ROOT's TH3::GetRandom3 for 5D
{
   // Read QQbar kinematical 5D grid's cell coordinates
   Int_t* cell = new Int_t[6];           // cell[6]={wght,iy1,iy2,ipt1,ipt2,idph}
   TTree* tG = (TTree*) fG->Get("tGqq");
   tG->GetBranch("cell")->SetAddress(cell);
   Int_t nbins = tG->GetEntries();
   Double_t rand[6];
   gRandom->RndmArray(6,rand);
   Int_t ibin = TMath::BinarySearch(nbins,fInt,rand[0]);
   tG->GetEvent(ibin);
   y1   = fgy[cell[1]]  + (fgy[cell[1]+1]-fgy[cell[1]])*rand[1];
   y2   = fgy[cell[2]]  + (fgy[cell[2]+1]-fgy[cell[2]])*rand[2];
   pt1  = fgpt[cell[3]] + (fgpt[cell[3]+1]-fgpt[cell[3]])*rand[3];
   pt2  = fgpt[cell[4]] + (fgpt[cell[4]+1]-fgpt[cell[4]])*rand[4];
   dphi = fgdph[cell[5]]+ (fgdph[cell[5]+1]-fgdph[cell[5]])*rand[5];
}

//____________________________________________________________________________________
void AliGenCorrHF::GetHadronPair(TFile* fG, Int_t idq, Double_t y1, Double_t y2, Double_t pt1, Double_t pt2, Int_t &id3, Int_t &id4, Double_t &pz3, Double_t &pz4, Double_t &pt3, Double_t &pt4) 
{
    // Generate a hadron pair
    Int_t (*fIpParaFunc )(TRandom*);//Pointer to particle type parametrisation function
    fIpParaFunc = IpCharm;
    Double_t mq = 1.2;              // c & b quark masses (used in AliPythia)
    if (idq == 5) {
      fIpParaFunc = IpBeauty;
      mq = 4.75;
    }
    Double_t z11, z12, z21, z22, pz1, pz2, e1, e2, mh, ptemp, rand[2];
    char tag[100]; 
    TH2F *h2h[12], *h2s[12];      // hard & soft Fragmentation Functions
    for (Int_t ipt = 0; ipt<fgnptbins; ipt++) { 
      sprintf(tag,"h2h_pt%d",ipt); 
      h2h[ipt] = (TH2F*) fG->Get(tag); 
      sprintf(tag,"h2s_pt%d",ipt); 
      h2s[ipt] = (TH2F*) fG->Get(tag); 
    }

       if (y1*y2 < 0) {
	 for (Int_t ipt = 0; ipt<fgnptbins; ipt++) { 
	   if(pt1 >= fgptbmin[ipt] && pt1 < fgptbmax[ipt]) 
      	     h2h[ipt]->GetRandom2(z11, z21);
	   if(pt2 >= fgptbmin[ipt] && pt2 < fgptbmax[ipt]) 
      	     h2h[ipt]->GetRandom2(z12, z22); 
	 }
       }
       else {
	 if (TMath::Abs(y1) > TMath::Abs(y2)) {
	   for (Int_t ipt = 0; ipt<fgnptbins; ipt++) { 
	     if(pt1 >= fgptbmin[ipt] && pt1 < fgptbmax[ipt]) 
	       h2h[ipt]->GetRandom2(z11, z21);
	     if(pt2 >= fgptbmin[ipt] && pt2 < fgptbmax[ipt]) 
	       h2s[ipt]->GetRandom2(z12, z22); 
	   }
	 }
	 else {
	   for (Int_t ipt = 0; ipt<fgnptbins; ipt++) { 
	     if(pt1 >= fgptbmin[ipt] && pt1 < fgptbmax[ipt]) 
	       h2s[ipt]->GetRandom2(z11, z21);
	     if(pt2 >= fgptbmin[ipt] && pt2 < fgptbmax[ipt]) 
	       h2h[ipt]->GetRandom2(z12, z22); 
	   }
	 }
       }
      gRandom->RndmArray(2,rand);
      ptemp = TMath::Sqrt(pt1*pt1 + mq*mq);
      pz1   = ptemp*TMath::SinH(y1); 
      e1    = ptemp*TMath::CosH(y1); 
      ptemp = TMath::Sqrt(pt2*pt2 + mq*mq);
      pz2   = ptemp*TMath::SinH(y2); 
      e2    = ptemp*TMath::CosH(y2); 

      id3   = fIpParaFunc(gRandom);
      mh    = TDatabasePDG::Instance()->GetParticle(id3)->Mass();
      ptemp = z11*z21*(e1*e1-pz1*pz1) - mh*mh;
      pt3   = (idq-3)*rand[0];                // some smearing at low pt, try better
      if (ptemp > 0) pt3 = TMath::Sqrt(ptemp);
      if (pz1 > 0)   pz3 = (z11*(e1 + pz1) - z21*(e1 - pz1)) / 2;
      else           pz3 = (z21*(e1 + pz1) - z11*(e1 - pz1)) / 2;
      e1 = TMath::Sqrt(pz3*pz3 + pt3*pt3 + mh*mh);

      id4   = - fIpParaFunc(gRandom);
      mh    = TDatabasePDG::Instance()->GetParticle(id4)->Mass();
      ptemp = z12*z22*(e2*e2-pz2*pz2) - mh*mh;
      pt4   = (idq-3)*rand[1];                // some smearing at low pt, try better
      if (ptemp > 0) pt4 = TMath::Sqrt(ptemp);
      if (pz2 > 0)   pz4 = (z12*(e2 + pz2) - z22*(e2 - pz2)) / 2;
      else           pz4 = (z22*(e2 + pz2) - z12*(e2 - pz2)) / 2;
      e2 = TMath::Sqrt(pz4*pz4 + pt4*pt4 + mh*mh);

      // small corr. instead of using Frag. Func. depending on yQ (in addition to ptQ)
      Float_t ycorr = 0.2, y3, y4;
      gRandom->RndmArray(2,rand);
      y3 = 0.5 * TMath::Log((e1 + pz3 + 1.e-13)/(e1 - pz3 + 1.e-13));
      y4 = 0.5 * TMath::Log((e2 + pz4 + 1.e-13)/(e2 - pz4 + 1.e-13));
      if(TMath::Abs(y3)<ycorr && TMath::Abs(y4)<ycorr && rand[0]>0.5) {
	ptemp = TMath::Sqrt(e1*e1 - pz3*pz3);
	y3  = 4*(1 - 2*rand[1]);
	pz3 = ptemp*TMath::SinH(y3);
	pz4 = pz3;
      }
}
