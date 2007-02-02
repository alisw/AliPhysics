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

// Class to generate particles from using paramtrized pT and y distributions.
// Distributions are obtained from pointer to object of type AliGenLib.
// (For example AliGenMUONlib)
// Decays are performed using Pythia.
// andreas.morsch@cern.ch

#include <TCanvas.h>
#include <TClonesArray.h>
#include <TDatabasePDG.h>
#include <TF1.h>
#include <TH1F.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TROOT.h>
#include <TVirtualMC.h>

#include "AliDecayer.h"
#include "AliGenMUONlib.h"
#include "AliGenParam.h"
#include "AliMC.h"
#include "AliRun.h"

ClassImp(AliGenParam)

//------------------------------------------------------------

  //Begin_Html
  /*
    <img src="picts/AliGenParam.gif">
  */
  //End_Html

//____________________________________________________________
    AliGenParam::AliGenParam():
	fPtParaFunc(0),
	fYParaFunc(0),
	fIpParaFunc(0),
	fPtPara(0),
	fYPara(0),
	fParam(0),
	fdNdy0(0.),
	fYWgt(0.),
	fPtWgt(0.),
	fBias(0.),
	fTrials(0),
	fDeltaPt(0.01),
	fDecayer(0)
{
// Default constructor
}
//____________________________________________________________
AliGenParam::AliGenParam(Int_t npart, AliGenLib * Library,  Int_t param, char* tname)
    :AliGenMC(npart),
     fPtParaFunc(Library->GetPt(param, tname)),
     fYParaFunc (Library->GetY (param, tname)),
     fIpParaFunc(Library->GetIp(param, tname)),
     fPtPara(0),
     fYPara(0),
     fParam(param),
     fdNdy0(0.),
     fYWgt(0.),
     fPtWgt(0.),
     fBias(0.),
     fTrials(0),
     fDeltaPt(0.01),
     fDecayer(0)
{
// Constructor using number of particles parameterisation id and library
    fName = "Param";
    fTitle= "Particle Generator using pT and y parameterisation";
    fAnalog = kAnalog;
    SetForceDecay();
}
//____________________________________________________________
AliGenParam::AliGenParam(Int_t npart, Int_t param, const char* tname, const char* name):
    AliGenMC(npart),
    fPtParaFunc(0),
    fYParaFunc (0),
    fIpParaFunc(0),
    fPtPara(0),
    fYPara(0),
    fParam(param),
    fdNdy0(0.),
    fYWgt(0.),
    fPtWgt(0.),
    fBias(0.),
    fTrials(0),
    fDeltaPt(0.01),
    fDecayer(0)
{
// Constructor using parameterisation id and number of particles
//
    fName = name;
    fTitle= "Particle Generator using pT and y parameterisation";
      
    AliGenLib* pLibrary = new AliGenMUONlib();
    fPtParaFunc = pLibrary->GetPt(param, tname);
    fYParaFunc  = pLibrary->GetY (param, tname);
    fIpParaFunc = pLibrary->GetIp(param, tname);
    
    fAnalog = kAnalog;
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

AliGenParam::AliGenParam(Int_t npart, Int_t param,
                         Double_t (*PtPara) (Double_t*, Double_t*),
                         Double_t (*YPara ) (Double_t* ,Double_t*),
		         Int_t    (*IpPara) (TRandom *))		 
    :AliGenMC(npart),
     
     fPtParaFunc(PtPara),
     fYParaFunc(YPara),
     fIpParaFunc(IpPara),
     fPtPara(0),
     fYPara(0),
     fParam(param),
     fdNdy0(0.),
     fYWgt(0.),
     fPtWgt(0.),
     fBias(0.),
     fTrials(0),
     fDeltaPt(0.01),
     fDecayer(0)
{
// Constructor
// Gines Martinez 1/10/99 
    fName = "Param";
    fTitle= "Particle Generator using pT and y parameterisation";

    fAnalog = kAnalog;
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
AliGenParam::~AliGenParam()
{
// Destructor
    delete  fPtPara;
    delete  fYPara;
}

//____________________________________________________________
void AliGenParam::Init()
{
// Initialisation

    if (gMC) fDecayer = gMC->GetDecayer();
  //Begin_Html
  /*
    <img src="picts/AliGenParam.gif">
  */
  //End_Html
    char name[256];
    sprintf(name, "pt-parameterisation for %s", GetName());
    
    if (fPtPara) fPtPara->Delete();
    fPtPara = new TF1(name, fPtParaFunc, fPtMin, fPtMax,0);
    gROOT->GetListOfFunctions()->Remove(fPtPara);
//  Set representation precision to 10 MeV
    Int_t npx= Int_t((fPtMax - fPtMin) / fDeltaPt);
    
    fPtPara->SetNpx(npx);

    sprintf(name, "y-parameterisation  for %s", GetName());
    if (fYPara) fYPara->Delete();
    fYPara  = new TF1(name, fYParaFunc, fYMin, fYMax, 0);
    gROOT->GetListOfFunctions()->Remove(fYPara);

    
    sprintf(name, "pt-for-%s", GetName());
    TF1 ptPara(name ,fPtParaFunc, 0, 15, 0);
    sprintf(name, "y-for-%s", GetName());
    TF1 yPara(name, fYParaFunc, -6, 6, 0);

//
// dN/dy| y=0
    Double_t y1=0;
    Double_t y2=0;
    
    fdNdy0=fYParaFunc(&y1,&y2);
//
// Integral over generation region
    Float_t intYS  = yPara.Integral(fYMin, fYMax,(Double_t*) 0x0,1.e-6);
    Float_t intPt0 = ptPara.Integral(0,15,(Double_t *) 0x0,1.e-6);
    Float_t intPtS = ptPara.Integral(fPtMin,fPtMax,(Double_t*) 0x0,1.e-6);
    Float_t phiWgt=(fPhiMax-fPhiMin)/2./TMath::Pi();
    if (fAnalog == kAnalog) {
	fYWgt  = intYS/fdNdy0;
	fPtWgt = intPtS/intPt0;
	fParentWeight = fYWgt*fPtWgt*phiWgt/fNpart;
    } else {
	fYWgt = intYS/fdNdy0;
	fPtWgt = (fPtMax-fPtMin)/intPt0;
	fParentWeight = fYWgt*fPtWgt*phiWgt/fNpart;
    }
//
// particle decay related initialization
    fDecayer->SetForceDecay(fForceDecay);
    fDecayer->Init();

//
    AliGenMC::Init();
}

//____________________________________________________________
void AliGenParam::Generate()
{
//
// Generate 'npart' of light and heavy mesons (J/Psi, upsilon or phi, Pion,
// Kaons, Etas, Omegas) and Baryons (proton, antiprotons, neutrons and 
// antineutrons in the the desired theta, phi and momentum windows; 
// Gaussian smearing on the vertex is done if selected. 
// The decay of heavy mesons is done using lujet, 
//    and the childern particle are tracked by GEANT
// However, light mesons are directly tracked by GEANT 
// setting fForceDecay = nodecay (SetForceDecay(nodecay)) 
//
//
//  Reinitialize decayer
  fDecayer->SetForceDecay(fForceDecay);
  fDecayer->Init();

//
  Float_t polar[3]= {0,0,0};  // Polarisation of the parent particle (for GEANT tracking)
  Float_t origin0[3];         // Origin of the generated parent particle (for GEANT tracking)
  Float_t pt, pl, ptot;       // Transverse, logitudinal and total momenta of the parent particle
  Float_t phi, theta;         // Phi and theta spherical angles of the parent particle momentum
  Float_t p[3], pc[3], 
          och[3];             // Momentum, polarisation and origin of the children particles from lujet
  Double_t ty, xmt;
  Int_t nt, i, j;
  Float_t  wgtp, wgtch;
  Double_t dummy;
  static TClonesArray *particles;
  //
  if(!particles) particles = new TClonesArray("TParticle",1000);
  
  TDatabasePDG *pDataBase = TDatabasePDG::Instance();
  //
  Float_t random[6];
 
// Calculating vertex position per event
  for (j=0;j<3;j++) origin0[j]=fOrigin[j];
  if(fVertexSmear==kPerEvent) {
      Vertex();
      for (j=0;j<3;j++) origin0[j]=fVertex[j];
  }
  
  Int_t ipa=0;
  
// Generating fNpart particles
  while (ipa<fNpart) {
      while(1) {
//
// particle type 
	  Int_t iPart = fIpParaFunc(fRandom);
	  fChildWeight=(fDecayer->GetPartialBranchingRatio(iPart))*fParentWeight;	   
	  TParticlePDG *particle = pDataBase->GetParticle(iPart);
	  Float_t am = particle->Mass();
	  
	  Rndm(random,2);
//
// phi
	  phi=fPhiMin+random[0]*(fPhiMax-fPhiMin);
//
// y
	  ty = TMath::TanH(fYPara->GetRandom());
//
// pT
	  if (fAnalog == kAnalog) {
	      pt=fPtPara->GetRandom();
	      wgtp=fParentWeight;
	      wgtch=fChildWeight;
	  } else {
	      pt=fPtMin+random[1]*(fPtMax-fPtMin);
	      Double_t ptd=pt;
	      wgtp=fParentWeight*fPtParaFunc(& ptd, &dummy);
	      wgtch=fChildWeight*fPtParaFunc(& ptd, &dummy);
	  }
	  xmt=sqrt(pt*pt+am*am);
	  if (TMath::Abs(ty)==1.) {
	      ty=0.;
	      Fatal("AliGenParam", 
		    "Division by 0: Please check you rapidity range !");
	  }
	  
	  pl=xmt*ty/sqrt(1.-ty*ty);
	  theta=TMath::ATan2(pt,pl);
// Cut on theta
	  if(theta<fThetaMin || theta>fThetaMax) continue;
	  ptot=TMath::Sqrt(pt*pt+pl*pl);
// Cut on momentum
	  if(ptot<fPMin || ptot>fPMax) continue;
//
	  p[0]=pt*TMath::Cos(phi);
	  p[1]=pt*TMath::Sin(phi);
	  p[2]=pl;
	  if(fVertexSmear==kPerTrack) {
	      Rndm(random,6);
	      for (j=0;j<3;j++) {
		  origin0[j]=
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
		if (!CheckAcceptanceGeometry(np,particles)) continue;
	      Int_t ncsel=0;
	      Int_t* pFlag      = new Int_t[np];
	      Int_t* pParent    = new Int_t[np];
	      Int_t* pSelected  = new Int_t[np];
	      Int_t* trackIt    = new Int_t[np];

	      for (i=0; i<np; i++) {
		  pFlag[i]     =  0;
		  pSelected[i] =  0;
		  pParent[i]   = -1;
	      }
	      
	      if (np >1) {
		  TParticle* iparticle =  (TParticle *) particles->At(0);
		  Int_t ipF, ipL;
		  for (i = 1; i<np ; i++) {
		      trackIt[i] = 1;
		      iparticle = (TParticle *) particles->At(i);
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
//			  TParticlePDG *particle = pDataBase->GetParticle(kf);
			  
			  Double_t lifeTime = fDecayer->GetLifetime(kf);
//			  Double_t mass     = particle->Mass();
//			  Double_t width    = particle->Width();
			  if (lifeTime > (Double_t) fMaxLifeTime) {
			      ipF = iparticle->GetFirstDaughter();
			      ipL = iparticle->GetLastDaughter();	
			      if (ipF > 0) for (j=ipF-1; j<ipL; j++) pFlag[j]=1;
			  } else{
			      trackIt[i]     = 0;
			      pSelected[i]   = 1;
			  }
		      } // ks==1 ?
//
// children
		      
		      if (ChildSelected(TMath::Abs(kf)) || fForceDecay == kAll && trackIt[i])
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
		  ipa++;
//
// Parent
		  PushTrack(0, -1, iPart, p, origin0, polar, 0, kPPrimary, nt, wgtp);
		  pParent[0] = nt;
		  KeepTrack(nt); 
//
// Decay Products
//		  
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
			 
			  PushTrack(fTrackIt*trackIt[i], iparent, kf,
					   pc, och, polar,
					   0, kPDecay, nt, wgtch);
			  pParent[i] = nt;
			  KeepTrack(nt); 
		      } // Selected
		  } // Particle loop 
	      }  // Decays by Lujet
	      particles->Clear();
	      if (pFlag)      delete[] pFlag;
	      if (pParent)    delete[] pParent;
	      if (pSelected)  delete[] pSelected;	   
	      if (trackIt)    delete[] trackIt;
	  } // kinematic selection
	  else  // nodecay option, so parent will be tracked by GEANT (pions, kaons, eta, omegas, baryons)
	  {
	    gAlice->GetMCApp()->
		PushTrack(fTrackIt,-1,iPart,p,origin0,polar,0,kPPrimary,nt,wgtp);
            ipa++; 
	  }
	  break;
    } // while
  } // event loop
  SetHighWaterMark(nt);
}
//____________________________________________________________________________________
Float_t AliGenParam::GetRelativeArea(Float_t ptMin, Float_t ptMax, Float_t yMin, Float_t yMax, Float_t phiMin, Float_t phiMax)
{
//
// Normalisation for selected kinematic region
//
  Float_t ratio =  
    fPtPara->Integral(ptMin,ptMax,(Double_t *)0,1.e-6) / fPtPara->Integral( fPtPara->GetXmin(), fPtPara->GetXmax(),(Double_t *)0,1.e-6) *
    fYPara->Integral(yMin,yMax,(Double_t *)0,1.e-6)/fYPara->Integral(fYPara->GetXmin(),fYPara->GetXmax(),(Double_t *)0,1.e-6)   *
    (phiMax-phiMin)/360.;
  return TMath::Abs(ratio);
}

//____________________________________________________________________________________

void AliGenParam::Draw( const char * /*opt*/)
{
    //
    // Draw the pT and y Distributions
    //
     TCanvas *c0 = new TCanvas("c0","Canvas 0",400,10,600,700);
     c0->Divide(2,1);
     c0->cd(1);
     fPtPara->Draw();
     fPtPara->GetHistogram()->SetXTitle("p_{T} (GeV)");     
     c0->cd(2);
     fYPara->Draw();
     fYPara->GetHistogram()->SetXTitle("y");     
}




