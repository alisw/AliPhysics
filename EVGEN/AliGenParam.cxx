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
Revision 1.42  2003/03/15 19:48:01  morsch
AliDecayerPythia replaced by AliDecayer

Revision 1.41  2003/01/09 17:38:47  morsch
Draw() method added.

Revision 1.40  2002/10/14 14:55:35  hristov
Merging the VirtualMC branch to the main development branch (HEAD)

Revision 1.36.6.3  2002/10/10 16:40:08  hristov
Updating VirtualMC to v3-09-02

Revision 1.39  2002/09/16 08:21:16  morsch
Use TDatabasePDG::Instance();

Revision 1.38  2002/05/30 14:59:12  morsch
Check geometrical acceptance. (G. Martinez)

Revision 1.37  2002/04/17 10:20:44  morsch
Coding Rule violations corrected.

Revision 1.36  2002/02/08 16:50:50  morsch
Add name and title in constructor.

Revision 1.35  2002/01/21 10:02:40  morsch
ty is Double_t
Abort if too high rapidity causes numerical paroblem. User has to specify
meaningful y-range.

Revision 1.34  2001/11/27 13:13:07  morsch
Maximum lifetime for long-lived particles to be put on the stack is parameter.
It can be set via SetMaximumLifetime(..).

Revision 1.33  2001/10/21 18:35:56  hristov
Several pointers were set to zero in the default constructors to avoid memory management problems

Revision 1.32  2001/07/27 17:09:36  morsch
Use local SetTrack, KeepTrack and SetHighWaterMark methods
to delegate either to local stack or to stack owned by AliRun.
(Piotr Skowronski, A.M.)

Revision 1.31  2001/07/13 10:58:54  morsch
- Some coded moved to AliGenMC
- Improved handling of secondary vertices.

Revision 1.30  2001/06/15 07:55:04  morsch
Put only first generation decay products on the stack.

Revision 1.29  2001/03/27 10:58:41  morsch
Initialize decayer before generation. Important if run inside cocktail.

Revision 1.28  2001/03/09 13:01:41  morsch
- enum constants for paramterisation type (particle family) moved to AliGen*lib.h
- use AliGenGSIlib::kUpsilon, AliGenPHOSlib::kEtaPrime to access the constants

Revision 1.27  2001/02/02 15:21:10  morsch
Set high water mark after last particle.
Use Vertex() method for Vertex.

Revision 1.26  2000/12/21 16:24:06  morsch
Coding convention clean-up

Revision 1.25  2000/11/30 07:12:50  alibrary
Introducing new Rndm and QA classes

Revision 1.24  2000/10/18 19:11:27  hristov
Division by zero fixed

Revision 1.23  2000/10/02 21:28:06  fca
Removal of useless dependecies via forward declarations

Revision 1.22  2000/09/12 14:14:55  morsch
Call fDecayer->ForceDecay() at the beginning of Generate().

Revision 1.21  2000/09/08 15:39:01  morsch
Handle the case fForceDecay=all during the generation, i.e. select all secondaries.

Revision 1.20  2000/09/06 14:35:44  morsch
Use AliDecayerPythia for particle decays.

Revision 1.19  2000/07/11 18:24:56  fca
Coding convention corrections + few minor bug fixes

Revision 1.18  2000/06/29 21:08:27  morsch
All paramatrisation libraries derive from the pure virtual base class AliGenLib.
This allows to pass a pointer to a library directly to AliGenParam and avoids the
use of function pointers in Config.C.

Revision 1.17  2000/06/09 20:33:30  morsch
All coding rule violations except RS3 corrected

Revision 1.16  2000/05/02 07:51:31  morsch
- Control precision of pT sampling TF1::SetNpx(..)
- Correct initialisation of child-cuts in all constructors.
- Most coding rule violations corrected.

Revision 1.15  2000/04/03 15:42:12  morsch
Cuts on primary particles are separated from those on the decay products. Methods
SetChildMomentumRange, SetChildPtRange, SetChildPhiRange, SetChildThetaRange added.

Revision 1.14  1999/11/09 07:38:48  fca
Changes for compatibility with version 2.23 of ROOT

Revision 1.13  1999/11/04 11:30:31  fca
Correct the logics for SetForceDecay

Revision 1.12  1999/11/03 17:43:20  fca
New version from G.Martinez & A.Morsch

Revision 1.11  1999/09/29 09:24:14  fca
Introduction of the Copyright and cvs Log

*/



// Class to generate particles from using paramtrized pT and y distributions.
// Distributions are obtained from pointer to object of type AliGenLib.
// (For example AliGenMUONlib)
// Decays are performed using Pythia.
// andreas.morsch@cern.ch

#include "AliGenParam.h"
#include "AliDecayer.h"
#include "AliGenMUONlib.h"
#include "AliRun.h"
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>

#include <TF1.h>
#include <TCanvas.h>
#include <TH1.h>

ClassImp(AliGenParam)

//------------------------------------------------------------

  //Begin_Html
  /*
    <img src="picts/AliGenParam.gif">
  */
  //End_Html

//____________________________________________________________
//____________________________________________________________
AliGenParam::AliGenParam()
{
// Deafault constructor
    fPtPara = 0;
    fYPara  = 0;
    fParam  = 0;
    fAnalog = kAnalog;
    SetDeltaPt();
    fDecayer = 0;


}

AliGenParam::AliGenParam(Int_t npart, AliGenLib * Library,  Int_t param, char* tname):AliGenMC(npart)
{
// Constructor using number of particles parameterisation id and library
    fName = "Param";
    fTitle= "Particle Generator using pT and y parameterisation";
    
    fPtParaFunc = Library->GetPt(param, tname);
    fYParaFunc  = Library->GetY (param, tname);
    fIpParaFunc = Library->GetIp(param, tname);
    
    fPtPara = 0;
    fYPara  = 0;
    fParam  = param;
    fAnalog = kAnalog;
    SetForceDecay();
    SetDeltaPt(); 
}

//____________________________________________________________

AliGenParam::AliGenParam(Int_t npart, Int_t param, char* tname):AliGenMC(npart)
{
// Constructor using parameterisation id and number of particles
//
    fName = "Param";
    fTitle= "Particle Generator using pT and y parameterisation";
      
    AliGenLib* pLibrary = new AliGenMUONlib();
 
    fPtParaFunc = pLibrary->GetPt(param, tname);
    fYParaFunc  = pLibrary->GetY (param, tname);
    fIpParaFunc = pLibrary->GetIp(param, tname);
    
    fPtPara = 0;
    fYPara  = 0;
    fParam  = param;
    fAnalog = kAnalog;
    fChildSelect.Set(5);
    for (Int_t i=0; i<5; i++) fChildSelect[i]=0;
    SetForceDecay();
    SetCutOnChild();
    SetChildMomentumRange();
    SetChildPtRange();
    SetChildPhiRange();
    SetChildThetaRange(); 
    SetDeltaPt(); 
}

AliGenParam::AliGenParam(Int_t npart, Int_t param,
                         Double_t (*PtPara) (Double_t*, Double_t*),
                         Double_t (*YPara ) (Double_t* ,Double_t*),
		         Int_t    (*IpPara) (TRandom *))		 
    :AliGenMC(npart)
{
// Constructor
// Gines Martinez 1/10/99 
    fName = "Param";
    fTitle= "Particle Generator using pT and y parameterisation";

    fPtParaFunc = PtPara; 
    fYParaFunc  = YPara;  
    fIpParaFunc = IpPara;
//  
    fPtPara = 0;
    fYPara  = 0;
    fParam  = param;
    fAnalog = kAnalog;
    fChildSelect.Set(5);
    for (Int_t i=0; i<5; i++) fChildSelect[i]=0;
    SetForceDecay();
    SetCutOnChild();
    SetChildMomentumRange();
    SetChildPtRange();
    SetChildPhiRange();
    SetChildThetaRange();  
    SetDeltaPt();
}


AliGenParam::AliGenParam(const AliGenParam & Param)
{
// copy constructor
    Param.Copy(*this);
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
 
    fPtPara = new TF1("Pt-Parametrization",fPtParaFunc,fPtMin,fPtMax,0);
//  Set representation precision to 10 MeV
    Int_t npx= Int_t((fPtMax-fPtMin)/fDeltaPt);
    
    fPtPara->SetNpx(npx);
    
    fYPara  = new TF1("Y -Parametrization",fYParaFunc,fYMin,fYMax,0);
    TF1* ptPara = new TF1("Pt-Parametrization",fPtParaFunc,0,15,0);
    TF1* yPara  = new TF1("Y -Parametrization",fYParaFunc,-6,6,0);

//
// dN/dy| y=0
    Double_t y1=0;
    Double_t y2=0;
    
    fdNdy0=fYParaFunc(&y1,&y2);
//
// Integral over generation region
    Float_t intYS  = yPara ->Integral(fYMin, fYMax);
    Float_t intPt0 = ptPara->Integral(0,15);
    Float_t intPtS = ptPara->Integral(fPtMin,fPtMax);
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
		  SetTrack(0, -1, iPart, p, origin0, polar, 0, kPPrimary, nt, wgtp);
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
			 
			  SetTrack(fTrackIt*trackIt[i], iparent, kf,
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
	    gAlice->
		SetTrack(fTrackIt,-1,iPart,p,origin0,polar,0,kPPrimary,nt,wgtp);
            ipa++; 
	  }
	  break;
    } // while
  } // event loop
  SetHighWaterMark(nt);
}

void AliGenParam::Draw()
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

AliGenParam& AliGenParam::operator=(const  AliGenParam& rhs)
{
// Assignment operator
    return *this;
}



