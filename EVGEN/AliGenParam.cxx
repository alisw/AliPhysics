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

#include "AliGenParam.h"
#include "AliDecayerPythia.h"
#include "AliGenMUONlib.h"
#include "AliRun.h"
#include "AliPythia.h"
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>

#include <TF1.h>

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
                 :AliGenerator()
{
// Deafault constructor
    fPtPara = 0;
    fYPara  = 0;
    fParam  = jpsi_p;
    fAnalog = analog;
    SetCutOnChild();
    SetChildMomentumRange();
    SetChildPtRange();
    SetChildPhiRange();
    SetChildThetaRange();  
    SetDeltaPt();
}

AliGenParam::AliGenParam(Int_t npart, AliGenLib * Library,  Param_t param, char* tname):AliGenerator(npart)
{
// Constructor using number of particles parameterisation id and library
    
    fPtParaFunc = Library->GetPt(param, tname);
    fYParaFunc  = Library->GetY (param, tname);
    fIpParaFunc = Library->GetIp(param, tname);
    
    fPtPara = 0;
    fYPara  = 0;
    fParam  = param;
    fAnalog = analog;
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

//____________________________________________________________

AliGenParam::AliGenParam(Int_t npart, Param_t param, char* tname):AliGenerator(npart)
{
// Constructor using parameterisation id and number of particles
//      
    AliGenLib* Library = new AliGenMUONlib();
 
    fPtParaFunc = Library->GetPt(param, tname);
    fYParaFunc  = Library->GetY (param, tname);
    fIpParaFunc = Library->GetIp(param, tname);
    
    fPtPara = 0;
    fYPara  = 0;
    fParam  = param;
    fAnalog = analog;
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

AliGenParam::AliGenParam(Int_t npart, Param_t param,
                         Double_t (*PtPara) (Double_t*, Double_t*),
                         Double_t (*YPara ) (Double_t* ,Double_t*),
		         Int_t    (*IpPara) ())			 
    :AliGenerator(npart)
{
// Constructor
// Gines Martinez 1/10/99 
    fPtParaFunc = PtPara; 
    fYParaFunc  = YPara;  
    fIpParaFunc = IpPara;
//  
    fPtPara = 0;
    fYPara  = 0;
    fParam  = param;
    fAnalog = analog;
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


AliGenParam::AliGenParam(const AliGenParam & Paramd)
{
// copy constructor
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

    fDecayer = new AliDecayerPythia();
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
    if (fAnalog == analog) {
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
    fDecayer->Init();
// semimuonic decays of charm and beauty
    fDecayer->ForceDecay(fForceDecay);
//
    switch (fForceDecay) 
    {
    case semielectronic:
    case dielectron:
    case b_jpsi_dielectron:
    case b_psip_dielectron:
	fChildSelect[0]=11;	
	break;
    case semimuonic:
    case dimuon:
    case b_jpsi_dimuon:
    case b_psip_dimuon:
	fChildSelect[0]=13;
	break;
    case pitomu:
	fChildSelect[0]=13;
	break;
    case katomu:
	fChildSelect[0]=13;
	break;
    case nodecay:
	break;
    case all:
	break;
    }
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


  Float_t polar[3]= {0,0,0};  // Polarisation of the parent particle (for GEANT tracking)
  Float_t origin0[3];         // Origin of the generated parent particle (for GEANT tracking)
  Float_t pt, pl, ptot;       // Transverse, logitudinal and total momenta of the parent particle
  Float_t phi, theta;         // Phi and theta spherical angles of the parent particle momentum
  Float_t p[3], pc[3], 
          och[3], pch[10][3]; // Momentum, polarisation and origin of the children particles from lujet
  Float_t ty, xmt;
  Int_t nt, i, j, kfch[10];
  Float_t  wgtp, wgtch;
  Double_t dummy;
  static TClonesArray *particles;
  //
  if(!particles) particles=new TClonesArray("TParticle",1000);

  static TDatabasePDG *DataBase = new TDatabasePDG();
  if(!DataBase) DataBase = new TDatabasePDG();
  //
  Float_t random[6];
 
// Calculating vertex position per event
  for (j=0;j<3;j++) origin0[j]=fOrigin[j];
  if(fVertexSmear==kPerEvent) {
      gMC->Rndm(random,6);
      for (j=0;j<3;j++) {
	  origin0[j]+=fOsigma[j]*TMath::Cos(2*random[2*j]*TMath::Pi())*
	      TMath::Sqrt(-2*TMath::Log(random[2*j+1]));
      }
  }
  Int_t ipa=0;
// Generating fNpart particles
  while (ipa<fNpart) {
      while(1) {
//
// particle type
	  Int_t iPart = fIpParaFunc();
	  fChildWeight=(fDecayer->GetPartialBranchingRatio(iPart))*fParentWeight;	   
	  TParticlePDG *particle = DataBase->GetParticle(iPart);
	  Float_t am = particle->Mass();

	  gMC->Rndm(random,2);
//
// phi
	  phi=fPhiMin+random[0]*(fPhiMax-fPhiMin);
//
// y
	  ty=Float_t(TMath::TanH(fYPara->GetRandom()));
//
// pT
	  if (fAnalog == analog) {
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
	  pl=xmt*ty/sqrt(1.-ty*ty);
	  theta=TMath::ATan2(pt,pl);
// Cut on theta
	  if(theta<fThetaMin || theta>fThetaMax) continue;
	  ptot=TMath::Sqrt(pt*pt+pl*pl);
// Cut on momentum
	  if(ptot<fPMin || ptot>fPMax) continue;
	  p[0]=pt*TMath::Cos(phi);
	  p[1]=pt*TMath::Sin(phi);
	  p[2]=pl;
	  if(fVertexSmear==kPerTrack) {
	      gMC->Rndm(random,6);
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
	  if (fForceDecay != nodecay) {
// Using lujet to decay particle
	      Float_t energy=TMath::Sqrt(ptot*ptot+am*am);
	      TLorentzVector pmom(p[0], p[1], p[2], energy);
	      fDecayer->Decay(iPart,&pmom);
//
// select decay particles
	      Int_t np=fDecayer->ImportParticles(particles);
	      Int_t ncsel=0;
	      for (i = 1; i<np; i++) {
		  TParticle *  iparticle = (TParticle *) particles->At(i);
		  Int_t kf = iparticle->GetPdgCode();
//
// children
		  if (ChildSelected(TMath::Abs(kf)))
		  {
		      pc[0]=iparticle->Px();
		      pc[1]=iparticle->Py();
		      pc[2]=iparticle->Pz();
		      och[0]=origin0[0]+iparticle->Vx()/10;
		      och[1]=origin0[1]+iparticle->Vy()/10;
		      och[2]=origin0[2]+iparticle->Vz()/10;
		      if (fCutOnChild) {
			  Float_t ptChild=TMath::Sqrt(pc[0]*pc[0]+pc[1]*pc[1]);
			  Float_t pChild=TMath::Sqrt(ptChild*ptChild+pc[2]*pc[2]);
			  Float_t thetaChild=TMath::ATan2(ptChild,pc[2]);
			  Float_t phiChild=TMath::ATan2(pc[1],pc[0]);
			  Bool_t childok = 
			      ((ptChild    > fChildPtMin    && ptChild    <fChildPtMax)      &&
			       (pChild     > fChildPMin     && pChild     <fChildPMax)       &&
			       (thetaChild > fChildThetaMin && thetaChild <fChildThetaMax)   &&
			       (phiChild   > fChildPhiMin   && phiChild   <fChildPhiMax));
			  if(childok)
			  {
			      pch[ncsel][0]=pc[0];
			      pch[ncsel][1]=pc[1];
			      pch[ncsel][2]=pc[2];
			      kfch[ncsel]=kf;
			      ncsel++;
			  } else {
			      ncsel=-1;
			      break;
			  } // child kine cuts
		      } else {
			  pch[ncsel][0]=pc[0];
			  pch[ncsel][1]=pc[1];
			  pch[ncsel][2]=pc[2];
			  kfch[ncsel]=kf;
			  ncsel++;
		      } // if child selection
		  } // select muon
	      } // decay particle loop
	      Int_t iparent;
	      if ((fCutOnChild && ncsel >0) || !fCutOnChild){
		  ipa++;
//
// parent
		  gAlice->
		      SetTrack(0,-1,iPart,p,origin0,polar,0,"Primary",nt,wgtp);
		  iparent=nt;
		  gAlice->KeepTrack(nt); 
		  for (i=0; i< ncsel; i++) {
		      gAlice->SetTrack(fTrackIt,iparent,kfch[i],
				       &pch[i][0],och,polar,
				       0,"Decay",nt,wgtch);
		      gAlice->KeepTrack(nt); 
		  }
	      }  // Decays by Lujet
	  } // kinematic selection
	  else  // nodecay option, so parent will be tracked by GEANT (pions, kaons, eta, omegas, baryons)
	  {
	    gAlice->
		SetTrack(fTrackIt,-1,iPart,p,origin0,polar,0,"Primary",nt,wgtp);
            ipa++; 
	  }
	  break;
    } // while
  } // event loop
}

Bool_t AliGenParam::ChildSelected(Int_t ip)
{
// True if particle is in list of selected children
    for (Int_t i=0; i<5; i++)
    {
	if (fChildSelect[i]==ip) return kTRUE;
    }
    return kFALSE;
}

Bool_t AliGenParam::KinematicSelection(TParticle *particle)
{
// Perform kinematic cuts
    Float_t px=particle->Px();
    Float_t py=particle->Py();
    Float_t pz=particle->Pz();
//
// momentum cut
    Float_t p=TMath::Sqrt(px*px+py*py+pz*pz);
    if (p > fPMax || p < fPMin) 
    {
//	printf("\n failed p cut %f %f %f \n",p,fPMin,fPMax);
	return kFALSE;
    }
    Float_t pt=TMath::Sqrt(px*px+py*py);
    
//
// theta cut
    Float_t  theta = Float_t(TMath::ATan2(Double_t(pt),Double_t(p)));
    if (theta > fThetaMax || theta < fThetaMin) 
    {
//	printf("\n failed theta cut %f %f %f \n",theta,fThetaMin,fThetaMax);
	return kFALSE;
    }

    return kTRUE;
}


AliGenParam& AliGenParam::operator=(const  AliGenParam& rhs)
{
// Assignment operator
    return *this;
}



