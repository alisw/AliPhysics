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
// ~2.10^6 Pythia6.214 events generated with kCharmppMNRwmi & kBeautyppMNRwmi, 
// CTEQ5L PDF and Pt_hard = 2.76 GeV/c for p-p collisions at 7, 10 and 14 TeV,
// and with kCharmppMNR (Pt_hard = 2.10 GeV/c) & kBeautyppMNR (Pt_hard = 2.75 GeV/c), 
// CTEQ4L PDF for Pb-Pb at 3.94 TeV, for p-Pb & Pb-p at 8.8 TeV. 
// Decays are performed by Pythia.
// Author: S. Grigoryan, LPC Clermont-Fd & YerPhI, Smbat.Grigoryan@cern.ch
// July 07: added quarks in the stack (B. Vulpescu)
// April 09: added energy choice between 10 and 14 TeV (S. Grigoryan)
// Sept 09: added hadron pair composition probabilities via 2D histo (X.M. Zhang)
// Oct 09: added energy choice between 7, 10, 14 TeV (for p-p), 4 TeV (for Pb-Pb),
// 9 TeV (for p-Pb) and -9 TeV (for Pb-p) (S. Grigoryan)
// April 10: removed "static" from definition of some variables (B. Vulpescu)
// May 11: added Flag for transportation of background particles while using 
// SetForceDecay() function (L. Manceau)
// June 11: added modifications allowing the setting of cuts on HF-hadron children.
// Quarks, hadrons and decay particles are loaded in the stack outside the loop
// of HF-hadrons, when the cuts on their children are satisfied (L. Manceau)
// Oct 11: added Pb-Pb at 2.76 TeV (S. Grigoryan)
// 
//-------------------------------------------------------------------------
// How it works (for the given flavor and p-p energy):
//
// 1) Reads QQbar kinematical grid (TTree) from the Input file and generates
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
// (see details in GetHadronPair). Fragmentation does not depend on p-p energy.
// 3) Decays the hadrons and saves all the particles in the event stack in the
// following order: HF hadron from Q, then its decay products, then HF hadron
// from Qbar, then its decay productes, then next HF hadon pair (if any) 
// in the same way, etc... 
// 4) It is fast, e.g., generates the same number of events with a beauty pair 
//  ~15 times faster than AliGenPythia with kBeautyppMNRwmi (w/o tracking)
//
// An Input file for each quark flavor and p-p energy is in EVGEN/dataCorrHF/
// One can use also user-defined Input files.
//
// More details could be found in my presentation at DiMuonNet Workshop, Dec 2006: 
// http://www-dapnia.cea.fr/Sphn/Alice/DiMuonNet.
//
//-------------------------------------------------------------------------
// How to use it:
//
// add the following typical lines in Config.C
/*
  if (!strcmp(option,"corr")) {
    // An example for correlated charm or beauty hadron pair production at 14 TeV

    // AliGenCorrHF *gener = new AliGenCorrHF(1, 4, 14);  // for charm, 1 pair per event
    AliGenCorrHF *gener = new AliGenCorrHF(1, 5, 14);  // for beauty, 1 pair per event

    gener->SetMomentumRange(0,9999);
    gener->SetCutOnChild(0);          // 1/0 means cuts on children enable/disable
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

#include <Riostream.h>
#include <TCanvas.h>
#include <TClonesArray.h>
#include <TDatabasePDG.h>
#include <TFile.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TROOT.h>
#include <TRandom.h>
#include <TTree.h>
#include <TVirtualMC.h>
#include <TVector3.h>

#include "AliGenCorrHF.h"
#include "AliLog.h"
#include "AliConst.h"
#include "AliDecayer.h"
#include "AliMC.h"
#include "AliRun.h"
#include "AliGenEventHeader.h"

ClassImp(AliGenCorrHF)

  //Begin_Html
  /*
    <img src="picts/AliGenCorrHF.gif">
  */
  //End_Html

Double_t AliGenCorrHF::fgdph[19] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180};
Double_t AliGenCorrHF::fgy[31] = {-10,-7, -6.5, -6, -5.5, -5, -4.5, -4, -3.5, -3, -2.5, -2,- 1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 10};
Double_t AliGenCorrHF::fgpt[51] = {0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.6, 7.2, 7.8, 8.4, 9, 9.6, 10.3, 11.1, 12, 13, 14, 15, 16, 17, 18, 19, 20.1, 21.5, 23, 24.5, 26, 27.5, 29.1, 31, 33, 35, 37, 39.2, 42, 45, 48, 51, 55.2, 60, 65, 71, 81, 100};
Int_t AliGenCorrHF::fgnptbins = 12;
Double_t AliGenCorrHF::fgptbmin[12] = {0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 9};
Double_t AliGenCorrHF::fgptbmax[12] = {0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 9, 100};

//____________________________________________________________
    AliGenCorrHF::AliGenCorrHF():
	fFileName(0),
	fFile(0),
	fQuark(0),
	fEnergy(0),
	fBias(0.),
	fTrials(0),
	fSelectAll(kFALSE),
	fDecayer(0),
	fgIntegral(0)
{
// Default constructor
}

//____________________________________________________________
AliGenCorrHF::AliGenCorrHF(Int_t npart, Int_t idquark, Int_t energy):
    AliGenMC(npart),
    fFileName(0),
    fFile(0),
    fQuark(idquark),
    fEnergy(energy),
    fBias(0.),
    fTrials(0),
    fSelectAll(kFALSE),
    fDecayer(0),
    fgIntegral(0)
{
// Constructor using particle number, quark type, energy & default InputFile
//
    if (fQuark == 5) {
      if (fEnergy == 7)
           fFileName = "$ALICE_ROOT/EVGEN/dataCorrHF/BeautyPP7PythiaMNRwmi.root";
      else if (fEnergy == 10)
           fFileName = "$ALICE_ROOT/EVGEN/dataCorrHF/BeautyPP10PythiaMNRwmi.root";
      else if (fEnergy == 14)
           fFileName = "$ALICE_ROOT/EVGEN/dataCorrHF/BeautyPP14PythiaMNRwmi.root";
      else if (fEnergy == 3)
	   fFileName = "$ALICE_ROOT/EVGEN/dataCorrHF/BeautyPbPb276PythiaMNR.root";
      else if (fEnergy == 4)
	   fFileName = "$ALICE_ROOT/EVGEN/dataCorrHF/BeautyPbPb394PythiaMNR.root";
      else if (fEnergy == 9 || fEnergy == -9)
	   fFileName = "$ALICE_ROOT/EVGEN/dataCorrHF/BeautyPPb88PythiaMNR.root";
      else fFileName = "$ALICE_ROOT/EVGEN/dataCorrHF/BeautyPbPb394PythiaMNR.root";
    }
    else {
      fQuark = 4;
      if (fEnergy == 7)
           fFileName = "$ALICE_ROOT/EVGEN/dataCorrHF/CharmPP7PythiaMNRwmi.root";
      else if (fEnergy == 10)
           fFileName = "$ALICE_ROOT/EVGEN/dataCorrHF/CharmPP10PythiaMNRwmi.root";
      else if (fEnergy == 14)
           fFileName = "$ALICE_ROOT/EVGEN/dataCorrHF/CharmPP14PythiaMNRwmi.root";
      else if (fEnergy == 3)
           fFileName = "$ALICE_ROOT/EVGEN/dataCorrHF/CharmPbPb276PythiaMNR.root";
      else if (fEnergy == 4)
           fFileName = "$ALICE_ROOT/EVGEN/dataCorrHF/CharmPbPb394PythiaMNR.root";
      else if (fEnergy == 9 || fEnergy == -9)
           fFileName = "$ALICE_ROOT/EVGEN/dataCorrHF/CharmPPb88PythiaMNR.root";
      else fFileName = "$ALICE_ROOT/EVGEN/dataCorrHF/CharmPbPb394PythiaMNR.root";
    }
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
AliGenCorrHF::AliGenCorrHF(char* tname, Int_t npart, Int_t idquark, Int_t energy):
    AliGenMC(npart),
    fFileName(tname),
    fFile(0),
    fQuark(idquark),
    fEnergy(energy),
    fBias(0.),
    fTrials(0),
    fSelectAll(kFALSE),
    fDecayer(0),
    fgIntegral(0)
{
// Constructor using particle number, quark type, energy & user-defined InputFile
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
AliGenCorrHF::~AliGenCorrHF()
{
// Destructor
  delete fFile;
}

//____________________________________________________________
void AliGenCorrHF::Init()
{
// Initialisation
  AliInfo(Form("Number of HF-hadron pairs = %d",fNpart)); 
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
//

//  Reinitialize decayer

  fDecayer->SetForceDecay(fForceDecay);
  fDecayer->Init();

  Float_t polar[2][3];        // Polarisation of the parent particle (for GEANT tracking)
  Float_t origin0[2][3];      // Origin of the generated parent particle (for GEANT tracking)
  Float_t pt, pl, ptot;       // Transverse, logitudinal and total momenta of the parent particle
  Float_t phi, theta;         // Phi and theta spherical angles of the parent particle momentum
  Float_t p[2][3];            // Momenta
  Int_t nt, i, j, ihad, ipa, ipa0, ipa1, ihadron[2], iquark[2];
  Float_t  wgtp[2], wgtch[2], random[6];
  Float_t pq[2][3], pc[3];    // Momenta of the two quarks
  Double_t tanhy2, qm = 0;
  Int_t np[2];
  Double_t dphi=0, ptq[2], yq[2], pth[2], plh[2], ph[2], phih[2], phiq[2];
  Int_t ncsel[2];
  Int_t** pSelected = new Int_t* [2];
  Int_t** trackIt = new Int_t* [2];

  for (i=0; i<2; i++) { 
    ptq[i]     =0; 
    yq[i]      =0; 
    pth[i]     =0; 
    plh[i]     =0;
    phih[i]    =0; 
    phiq[i]    =0;
    ihadron[i] =0; 
    iquark[i]  =0;
    for (j=0; j<3; j++) polar[i][j]=0;
  }

  // same quarks mass as in the fragmentation functions
  if (fQuark == 4) qm = 1.20;
  else             qm = 4.75;
  
  TClonesArray *particleshad1 = new TClonesArray("TParticle",1000);
  TClonesArray *particleshad2 = new TClonesArray("TParticle",1000);
  
  TList *particleslist = new TList();
  particleslist->Add(particleshad1);
  particleslist->Add(particleshad2);
  
  TDatabasePDG *pDataBase = TDatabasePDG::Instance();

  // Calculating vertex position per event
  for (i=0;i<2;i++){
    for (j=0;j<3;j++) origin0[i][j]=fOrigin[j];
    if (fVertexSmear==kPerEvent) {
      Vertex();
      for (j=0;j<3;j++) origin0[i][j]=fVertex[j];
    }
  }
  
  ipa  = 0;
  ipa1 = 0;
  ipa0 = 0;
  
  // Generating fNpart HF-hadron pairs
  fNprimaries = 0;
 
  while (ipa<2*fNpart) {

    GetQuarkPair(fFile, fgIntegral, yq[0], yq[1], ptq[0], ptq[1], dphi);
    
    GetHadronPair(fFile, fQuark, yq[0], yq[1], ptq[0], ptq[1], ihadron[0], ihadron[1], plh[0], plh[1], pth[0], pth[1]);
    
    if (fEnergy == 9 || fEnergy == -9) {      // boost particles from c.m.s. to ALICE lab frame
      Double_t dyBoost = 0.47;
      Double_t beta  = TMath::TanH(dyBoost);
      Double_t gamma = 1./TMath::Sqrt((1.-beta)*(1.+beta));
      Double_t gb    = gamma * beta;
      yq[0] += dyBoost;
      yq[1] += dyBoost;
      plh[0] = gb * TMath::Sqrt(plh[0]*plh[0] + pth[0]*pth[0]) + gamma * plh[0];
      plh[1] = gb * TMath::Sqrt(plh[1]*plh[1] + pth[1]*pth[1]) + gamma * plh[1];
      if (fEnergy == 9) {
	yq[0] *= -1;
	yq[1] *= -1;
	plh[0] *= -1;
	plh[1] *= -1;
      }
    }      
    
    // Cuts from AliGenerator
    
    // Cut on theta
    theta=TMath::ATan2(pth[0],plh[0]);
    if (theta<fThetaMin || theta>fThetaMax) continue;
    theta=TMath::ATan2(pth[1],plh[1]);
    if (theta<fThetaMin || theta>fThetaMax) continue;
    
    // Cut on momentum
    ph[0]=TMath::Sqrt(pth[0]*pth[0]+plh[0]*plh[0]);
    if (ph[0]<fPMin || ph[0]>fPMax) continue;
    ph[1]=TMath::Sqrt(pth[1]*pth[1]+plh[1]*plh[1]);
    if (ph[1]<fPMin || ph[1]>fPMax) continue;
    
    // Add the quarks in the stack
    
    phiq[0] = Rndm()*k2PI;
    if (Rndm() < 0.5) {
      phiq[1] = phiq[0] + dphi*kDegrad; 
    } else {
      phiq[1] = phiq[0] - dphi*kDegrad; 
    }    
    if (phiq[1] > k2PI) phiq[1] -= k2PI;
    if (phiq[1] < 0   ) phiq[1] += k2PI;
    
    // quarks pdg
    iquark[0] = +fQuark;
    iquark[1] = -fQuark;
    
    // px and py
    TVector2 qvect1 = TVector2();
    TVector2 qvect2 = TVector2();
    qvect1.SetMagPhi(ptq[0],phiq[0]);
    qvect2.SetMagPhi(ptq[1],phiq[1]);
    pq[0][0] = qvect1.Px();
    pq[0][1] = qvect1.Py();
    pq[1][0] = qvect2.Px();
    pq[1][1] = qvect2.Py();
    
    // pz
    tanhy2 = TMath::TanH(yq[0]);
    tanhy2 *= tanhy2;
    pq[0][2] = TMath::Sqrt((ptq[0]*ptq[0]+qm*qm)*tanhy2/(1-tanhy2));
    pq[0][2] = TMath::Sign((Double_t)pq[0][2],yq[0]);
    tanhy2 = TMath::TanH(yq[1]);
    tanhy2 *= tanhy2;
    pq[1][2] = TMath::Sqrt((ptq[1]*ptq[1]+qm*qm)*tanhy2/(1-tanhy2));
    pq[1][2] = TMath::Sign((Double_t)pq[1][2],yq[1]);
    
    // Here we assume that  |phi_H1 - phi_H2| = |phi_Q1 - phi_Q2| = dphi
    // which is a good approximation for heavy flavors in Pythia
    // ... moreover, same phi angles as for the quarks ...
    
    phih[0] = phiq[0];    
    phih[1] = phiq[1];    
    
    ipa1 = 0;

    for (ihad = 0; ihad < 2; ihad++) {
      while(1) {
	
	ipa0=ipa1;
	
	// particle type 
	fChildWeight=(fDecayer->GetPartialBranchingRatio(ihadron[ihad]))*fParentWeight;
	wgtp[ihad]=fParentWeight;
	wgtch[ihad]=fChildWeight;
	TParticlePDG *particle = pDataBase->GetParticle(ihadron[ihad]);
	Float_t am = particle->Mass();
	phi = phih[ihad];
	pt  = pth[ihad];
	pl  = plh[ihad];
	ptot=TMath::Sqrt(pt*pt+pl*pl);
	
	p[ihad][0]=pt*TMath::Cos(phi);
	p[ihad][1]=pt*TMath::Sin(phi);
	p[ihad][2]=pl;
	
	if(fVertexSmear==kPerTrack) {
	  Rndm(random,6);
	  for (j=0;j<3;j++) {
	    origin0[ihad][j]=
	      fOrigin[j]+fOsigma[j]*TMath::Cos(2*random[2*j]*TMath::Pi())*
	      TMath::Sqrt(-2*TMath::Log(random[2*j+1]));
	  }
	}
	
	// Looking at fForceDecay : 
	// if fForceDecay != none Primary particle decays using 
	// AliPythia and children are tracked by GEANT
	//
	// if fForceDecay == none Primary particle is tracked by GEANT 
	// (In the latest, make sure that GEANT actually does all the decays you want)	  

	if (fForceDecay != kNoDecay) {
	  // Using lujet to decay particle
	  Float_t energy=TMath::Sqrt(ptot*ptot+am*am);
	  TLorentzVector pmom(p[ihad][0], p[ihad][1], p[ihad][2], energy);
	  fDecayer->Decay(ihadron[ihad],&pmom);
	  
	  // select decay particles
	 
	  np[ihad]=fDecayer->ImportParticles((TClonesArray *)particleslist->At(ihad));

	  //  Selecting  GeometryAcceptance for particles fPdgCodeParticleforAcceptanceCut;
	  
	  if (fGeometryAcceptance) 
	    if (!CheckAcceptanceGeometry(np[ihad],(TClonesArray*)particleslist->At(ihad))) continue;
	  
	  trackIt[ihad]     = new Int_t [np[ihad]];
	  pSelected[ihad]   = new Int_t [np[ihad]];
	  Int_t* pFlag      = new Int_t [np[ihad]];
	  
	  for (i=0; i<np[ihad]; i++) {
	    pFlag[i]     =  0;
	    pSelected[ihad][i] =  0;
	  }

	  if (np[ihad] >1) {
	    TParticle* iparticle = 0;
	    Int_t ipF, ipL;
	    
	    for (i = 1; i<np[ihad] ; i++) {
	      trackIt[ihad][i] = 1;
	      iparticle = 
		(TParticle *) ((TClonesArray *) particleslist->At(ihad))->At(i);
	      Int_t kf = iparticle->GetPdgCode();
	      Int_t ks = iparticle->GetStatusCode();
	      // flagged particle
	      if (pFlag[i] == 1) {
		ipF = iparticle->GetFirstDaughter();
		ipL = iparticle->GetLastDaughter();	
		if (ipF > 0) for (j=ipF-1; j<ipL; j++) pFlag[j]=1;
		continue;
	      }
	      
	      // flag decay products of particles with long life-time (c tau > .3 mum)
	      if (ks != 1) { 
		Double_t lifeTime = fDecayer->GetLifetime(kf);
		if (lifeTime > (Double_t) fMaxLifeTime) {
		  ipF = iparticle->GetFirstDaughter();
		  ipL = iparticle->GetLastDaughter();	
		  if (ipF > 0) for (j=ipF-1; j<ipL; j++) pFlag[j]=1;
		} else {
		  trackIt[ihad][i]     = 0;
		  pSelected[ihad][i]   = 1;
		}
	      } // ks==1 ?
	      //
	      // children
	      if ((ChildSelected(TMath::Abs(kf)) || fForceDecay == kAll || fSelectAll) && trackIt[ihad][i])
		{      
		  if (fCutOnChild) {
		    pc[0]=iparticle->Px();
		    pc[1]=iparticle->Py();
		    pc[2]=iparticle->Pz();
		    //printf("px %f py %f pz %f\n",pc[0],pc[1],pc[2]);
		    Bool_t  childok = KinematicSelection(iparticle, 1);
		    if(childok) {
		      pSelected[ihad][i]  = 1;
		      ncsel[ihad]++;
		    } else {
		      ncsel[ihad]=-1;
		      break;
		    } // child kine cuts
		  } else {
		    pSelected[ihad][i]  = 1;
		    ncsel[ihad]++;
		  } // if child selection
		} // select muon
	    } // decay particle loop
	  } // if decay products

	  if ((fCutOnChild && ncsel[ihad] >0) || !fCutOnChild) ipa1++;

	  if (pFlag) delete[] pFlag;

	} // kinematic selection
	else  // nodecay option, so parent will be tracked by GEANT (pions, kaons, eta, omegas, baryons)
	  {
	    gAlice->GetMCApp()->
	      PushTrack(fTrackIt,-1,ihadron[ihad],p[ihad],origin0[ihad],polar[ihad],0,kPPrimary,nt,wgtp[ihad]);
	    ipa1++; 
	    fNprimaries++;
	    
	  }
	break;
      } // while(1) loop
      if (ipa1<ipa0+1){
	ipa1=0; 
	if (pSelected[ihad]) delete pSelected[ihad];
	if (trackIt[ihad])   delete trackIt[ihad];
	particleshad1->Clear();
	particleshad2->Clear();
	break;
      }//go out of loop and generate new pair if at least one hadron is rejected
    } // hadron pair loop
    if(ipa1==2){ 
 
      ipa=ipa+ipa1;
 
      if(fForceDecay != kNoDecay){
	for(ihad=0;ihad<2;ihad++){

	  //load tracks in the stack if both hadrons of the pair accepted
	  LoadTracks(iquark[ihad],pq[ihad],ihadron[ihad],p[ihad],np[ihad],
		     (TClonesArray *)particleslist->At(ihad),origin0[ihad],
		     polar[ihad],wgtp[ihad],wgtch[ihad],nt,ncsel[ihad],
		     pSelected[ihad],trackIt[ihad]);

	  if (pSelected[ihad]) delete pSelected[ihad];
	  if (trackIt[ihad])   delete trackIt[ihad];

	}
	  particleshad1->Clear();
	  particleshad2->Clear();
      }
    }
  }   // while (ipa<2*fNpart) loop
  
  SetHighWaterMark(nt);
  
  AliGenEventHeader* header = new AliGenEventHeader("CorrHF");
  header->SetPrimaryVertex(fVertex);
  header->SetNProduced(fNprimaries);
  AddHeader(header);
  
  
  delete particleshad1;
  delete particleshad2;
  delete particleslist;
 
  delete[] pSelected;
  delete[] trackIt;
}
//____________________________________________________________________________________
void AliGenCorrHF::IpCharm(TH2F *hProbHH, Int_t &pdg3, Int_t &pdg4)
{  
// Composition of a lower state charm hadron pair from a ccbar quark pair
   Int_t pdgH[] = {411, 421, 431, 4122, 4132, 4232, 4332};

   Double_t id3, id4;
   hProbHH->GetRandom2(id3, id4);
   pdg3 = pdgH[(Int_t)TMath::Floor(id3)];
   pdg4 = -1*pdgH[(Int_t)TMath::Floor(id4)];

   return;
}

void AliGenCorrHF::IpBeauty(TH2F *hProbHH, Int_t &pdg3, Int_t &pdg4)
{  
// Composition of a lower state beauty hadron pair from a bbbar quark pair
   // B-Bbar mixing will be done by Pythia at their decay point
   Int_t pdgH[] = {511, 521, 531, 5122, 5132, 5232, 5332};

   Double_t id3, id4;
   hProbHH->GetRandom2(id3, id4);
   pdg3 = pdgH[(Int_t)TMath::Floor(id3)];
   pdg4 = -1*pdgH[(Int_t)TMath::Floor(id4)];

   if ( (pdg3== 511) || (pdg3== 521) || (pdg3== 531) ) pdg3 *= -1;
   if ( (pdg4==-511) || (pdg4==-521) || (pdg4==-531) ) pdg4 *= -1;

   return;
}

//____________________________________________________________________________________
Double_t AliGenCorrHF::ComputeIntegral(TFile* fG)       // needed by GetQuarkPair
{
   // Read QQbar kinematical 5D grid's cell occupancy weights
   Int_t cell[6];           // cell[6]={wght,iy1,iy2,ipt1,ipt2,idph}
   TTree* tG = (TTree*) fG->Get("tGqq");
   tG->GetBranch("cell")->SetAddress(&cell);
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
   Int_t cell[6];           // cell[6]={wght,iy1,iy2,ipt1,ipt2,idph}
   TTree* tG = (TTree*) fG->Get("tGqq");
   tG->GetBranch("cell")->SetAddress(&cell);
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
    void (*fIpParaFunc)(TH2F *, Int_t &, Int_t &);//Pointer to hadron pair composition function
    fIpParaFunc = IpCharm;
    Double_t mq = 1.2;              // c & b quark masses (used in AliPythia)
    if (idq == 5) {
      fIpParaFunc = IpBeauty;
      mq = 4.75;
    }
    Double_t z11 = 0.;
    Double_t z12 = 0.;
    Double_t z21 = 0.;
    Double_t z22 = 0.;
    Double_t pz1, pz2, e1, e2, mh, ptemp, rand[2];
    char tag[100]; 
    TH2F *h2h[12], *h2s[12], *hProbHH; // hard & soft fragmentation and HH-probability functions
    for (Int_t ipt = 0; ipt<fgnptbins; ipt++) { 
      snprintf(tag,100, "h2h_pt%d",ipt); 
      h2h[ipt] = (TH2F*) fG->Get(tag); 
      snprintf(tag,100, "h2s_pt%d",ipt); 
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

      hProbHH = (TH2F*)fG->Get("hProbHH");
      fIpParaFunc(hProbHH, id3, id4);
      mh    = TDatabasePDG::Instance()->GetParticle(id3)->Mass();
      ptemp = z11*z21*(e1*e1-pz1*pz1) - mh*mh;
      if (idq==5) pt3   = pt1;                // an approximation at low pt, try better
      else        pt3   = rand[0];            // pt3=pt1 gives less D-hadrons at low pt 
      if (ptemp > 0) pt3 = TMath::Sqrt(ptemp);
      if (pz1 > 0)   pz3 = (z11*(e1 + pz1) - z21*(e1 - pz1)) / 2;
      else           pz3 = (z21*(e1 + pz1) - z11*(e1 - pz1)) / 2;
      e1 = TMath::Sqrt(pz3*pz3 + pt3*pt3 + mh*mh);

      mh    = TDatabasePDG::Instance()->GetParticle(id4)->Mass();
      ptemp = z12*z22*(e2*e2-pz2*pz2) - mh*mh;
      if (idq==5) pt4   = pt2;                // an approximation at low pt, try better
      else        pt4   = rand[1];
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
	ptemp = TMath::Sqrt((e1-pz3)*(e1+pz3));
	y3  = 4*(1 - 2*rand[1]);
	pz3 = ptemp*TMath::SinH(y3);
	pz4 = pz3;
      }
}

		
//____________________________________________________________________________________
void AliGenCorrHF::LoadTracks(Int_t iquark, Float_t *pq, 
			      Int_t iPart, Float_t *p, 
			      Int_t np, TClonesArray *particles,
			      Float_t *origin0, Float_t *polar, 
			      Float_t wgtp, Float_t wgtch,
			      Int_t &nt, Int_t ncsel, Int_t *pSelected, 
			      Int_t *trackIt){
  Int_t i; 
  Int_t ntq=-1;
  Int_t* pParent = new Int_t[np];
  Float_t pc[3], och[3];
  Int_t iparent;

  for(i=0;i<np;i++) pParent[i]=-1;

  if ((fCutOnChild && ncsel >0) || !fCutOnChild){  
    // Parents
    // quark
    PushTrack(0, -1, iquark, pq, origin0, polar, 0, kPPrimary, nt, wgtp);
    KeepTrack(nt);
    ntq = nt;
    // hadron
    PushTrack(0, ntq, iPart, p, origin0, polar, 0, kPDecay, nt, wgtp);
    pParent[0] = nt;
    KeepTrack(nt); 
    fNprimaries++;

    // Decay Products  
    for (i = 1; i < np; i++) {
      if (pSelected[i]) {

	TParticle* iparticle = (TParticle *) particles->At(i);
	Int_t kf  = iparticle->GetPdgCode();
	Int_t jpa = iparticle->GetFirstMother()-1;

	och[0] = origin0[0]+iparticle->Vx()/10;
	och[1] = origin0[1]+iparticle->Vy()/10;
	och[2] = origin0[2]+iparticle->Vz()/10;
	pc[0]  = iparticle->Px();
	pc[1]  = iparticle->Py();
	pc[2]  = iparticle->Pz();
	
	if (jpa > -1) {
	  iparent = pParent[jpa];
	} else {
	  iparent = -1;
	}
	
	PushTrack(fTrackIt*trackIt[i], iparent, kf,
		  pc, och, polar,
		  0, kPDecay, nt, wgtch);
	pParent[i] = nt;
	KeepTrack(nt); 
	fNprimaries++;

      } // Selected
    } // Particle loop
  }
  if (pParent) delete[] pParent;
 
  return;
}

