#include "AliGenerator.h"
#include "AliGenPythia.h"
#include "TGeant3.h"
#include "AliRun.h"
#include "AliPythia.h"
#include <TDirectory.h>
#include <TFile.h>
#include <TTree.h>
#include <stdlib.h>
#include <AliPythia.h>
#include <TMCParticle.h>
#include <GParticle.h>
 ClassImp(AliGenPythia)

AliGenPythia::AliGenPythia()
                 :AliGenerator()
{
}

AliGenPythia::AliGenPythia(Int_t npart)
                 :AliGenerator(npart)
{
// default charm production at 5. 5 TeV
// semimuonic decay
// structure function GRVHO
//
    fXsection  = 0.;
    fParentSelect.Set(5);
    fChildSelect.Set(5);
    for (Int_t i=0; i<5; i++) fParentSelect[i]=fChildSelect[i]=0;
    SetProcess();
    SetStrucFunc();
    ForceDecay();
    SetPtHard();
    SetEnergyCMS();
}

AliGenPythia::~AliGenPythia()
{
}

void AliGenPythia::Init()
{
    SetMC(new AliPythia());
    fPythia=(AliPythia*) fgMCEvGen;
//
    fParentWeight=1./Float_t(fNpart);
//
//  Forward Paramters to the AliPythia object
    fPythia->DefineParticles();
    fPythia->SetCKIN(3,fPtHardMin);
    fPythia->SetCKIN(4,fPtHardMax);    
    fPythia->ProcInit(fProcess,fEnergyCMS,fStrucFunc);
    fPythia->ForceDecay(fForceDecay);
    fPythia->LuList(0);
    fPythia->PyStat(2);
//  Parent and Children Selection
    switch (fProcess) 
    {
    case charm:

	fParentSelect[0]=411;
	fParentSelect[1]=421;
	fParentSelect[2]=431;
	fParentSelect[3]=4122;	
	break;
    case charm_unforced:

	fParentSelect[0]=411;
	fParentSelect[1]=421;
	fParentSelect[2]=431;
	fParentSelect[3]=4122;	
	break;
    case beauty:
	fParentSelect[0]=511;
	fParentSelect[1]=521;
	fParentSelect[2]=531;
	fParentSelect[3]=5122;	
	break;
    case beauty_unforced:
	fParentSelect[0]=511;
	fParentSelect[1]=521;
	fParentSelect[2]=531;
	fParentSelect[3]=5122;	
	break;
    case jpsi_chi:
    case jpsi:
	fParentSelect[0]=443;
	break;
    case mb:
	break;
    }

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
    case pitomu:
    case katomu:
	fChildSelect[0]=13;
	break;
    }
}

void AliGenPythia::Generate()
{
    AliMC* pMC = AliMC::GetMC();

    Float_t polar[3] =   {0,0,0};
    Float_t origin[3]=   {0,0,0};
    Float_t origin_p[3]= {0,0,0};
    Float_t origin0[3]=  {0,0,0};
    Float_t p[3], p_p[4], random[6];
//  converts from mm/c to s
    const Float_t kconv=0.001/2.999792458e8;
    
    
//
    Int_t nt=0;
    Int_t nt_p=0;
    Int_t jev=0;
    Int_t j;
    
    fTrials=0;
    for (j=0;j<3;j++) origin0[j]=fOrigin[j];
    if(fVertexSmear==perEvent) {
	pMC->Rndm(random,6);
	for (j=0;j<3;j++) {
	    origin0[j]+=fOsigma[j]*TMath::Cos(2*random[2*j]*TMath::Pi())*
		TMath::Sqrt(-2*TMath::Log(random[2*j+1]));
	    fPythia->SetMSTP(151,0);
	}
    } else if (fVertexSmear==perTrack) {
	fPythia->SetMSTP(151,0);
	for (j=0;j<3;j++) {
	    fPythia->SetPARP(151+j, fOsigma[j]*10.);
	}
    }
    while(1)
    {
	fPythia->PyEvnt();
	fTrials++;
	TObjArray* particles = fPythia->GetPrimaries() ;
	Int_t np = particles->GetEntriesFast();
	printf("\n **************************************************%d\n",np);
	Int_t nc=0;
	if (np == 0 ) continue;
	if (fProcess != mb) {
	    for (Int_t i = 0; i<np; i++) {
		TMCParticle *  iparticle = (TMCParticle *) particles->At(i);
		Int_t kf = iparticle->GetKF();
		fChildWeight=(fPythia->GetBraPart(kf))*fParentWeight;	  
//
// Parent
		if (ParentSelected(TMath::Abs(kf))) {
		    if (KinematicSelection(iparticle)) {
			if (nc==0) {
//
// Store information concerning the hard scattering process
//
			    Float_t mass_p = fPythia->GetPARI(13);
			    Float_t   pt_p = fPythia->GetPARI(17);
			    Float_t    y_p = fPythia->GetPARI(37);
			    Float_t  xmt_p = sqrt(pt_p*pt_p+mass_p*mass_p);
			    Float_t     ty = Float_t(TMath::TanH(y_p));
			    p_p[0] = pt_p;
			    p_p[1] = 0;
			    p_p[2] = xmt_p*ty/sqrt(1.-ty*ty);
			    p_p[3] = mass_p;
			    gAlice->SetTrack(0,-1,-1,
					     p_p,origin_p,polar,
					     0,"Hard Scat.",nt_p,fParentWeight);
			    gAlice->KeepTrack(nt_p);
			}
			nc++;
//
// store parent track information
			p[0]=iparticle->GetPx();
			p[1]=iparticle->GetPy();
			p[2]=iparticle->GetPz();
			origin[0]=origin0[0]+iparticle->GetVx()/10;
			origin[1]=origin0[1]+iparticle->GetVy()/10;
			origin[2]=origin0[2]+iparticle->GetVz()/10;

			Int_t ifch=iparticle->GetFirstChild();
			Int_t ilch=iparticle->GetLastChild();	
			if (ifch !=0 && ilch !=0) {
			    gAlice->SetTrack(0,nt_p,fPythia->GetGeantCode(kf),
					     p,origin,polar,
					     0,"Primary",nt,fParentWeight);
			    gAlice->KeepTrack(nt);
			    Int_t iparent = nt;
//
// Children	    

			    for (Int_t j=ifch; j<=ilch; j++)
			    {
				TMCParticle *  ichild = 
				    (TMCParticle *) particles->At(j-1);
				Int_t kf = ichild->GetKF();
//
// 
				if (ChildSelected(TMath::Abs(kf))) {
				    Int_t kg=fPythia->GetGeantCode(kf);
				    origin[0]=ichild->GetVx();
				    origin[1]=ichild->GetVy();
				    origin[2]=ichild->GetVz();		
				    p[0]=ichild->GetPx();
				    p[1]=ichild->GetPy();
				    p[2]=ichild->GetPz();
				    Float_t tof=kconv*ichild->GetTime();
				    gAlice->SetTrack(fTrackIt, iparent, kg,
						     p,origin,polar,
						     tof,"Decay",nt,fChildWeight);
				    gAlice->KeepTrack(nt);
				} // select child
			    } // child loop
			}
		    } // kinematic selection
		} // select particle
	    } // particle loop
	} else {
	    for (Int_t i = 0; i<np; i++) {
		TMCParticle *  iparticle = (TMCParticle *) particles->At(i);
		Int_t kf = iparticle->GetKF();
		Int_t ks = iparticle->GetKS();
		Int_t gc = fPythia->GetGeantCode(kf);
		if (ks==1 && gc!=0 && KinematicSelection(iparticle)) {
			nc++;
//
// store track information
			p[0]=iparticle->GetPx();
			p[1]=iparticle->GetPy();
			p[2]=iparticle->GetPz();
			origin[0]=origin0[0]+iparticle->GetVx()/10;
			origin[1]=origin0[1]+iparticle->GetVy()/10;
			origin[2]=origin0[2]+iparticle->GetVz()/10;
			Float_t tof=kconv*iparticle->GetTime();
			gAlice->SetTrack(fTrackIt,-1,gc,p,origin,polar,
					 tof,"Primary",nt);
			gAlice->KeepTrack(nt);
		} // select particle
	    } // particle loop 
	    printf("\n I've put %i particles on the stack \n",nc);
	} // mb ?
	if (nc > 0) {
	    jev+=nc;
	    if (jev >= fNpart) {
		fKineBias=Float_t(fNpart)/Float_t(fTrials);
		printf("\n Trials: %i\n",fTrials);
		break;
	    }
	}
    } // event loop
//  adjust weight due to kinematic selection
    AdjustWeights();
//  get cross-section
    fXsection=fPythia->GetPARI(1);
}

Bool_t AliGenPythia::ParentSelected(Int_t ip)
{
    for (Int_t i=0; i<5; i++)
    {
	if (fParentSelect[i]==ip) return kTRUE;
    }
    return kFALSE;
}

Bool_t AliGenPythia::ChildSelected(Int_t ip)
{
    for (Int_t i=0; i<5; i++)
    {
	if (fChildSelect[i]==ip) return kTRUE;
    }
    return kFALSE;
}

Bool_t AliGenPythia::KinematicSelection(TMCParticle *particle)
{
    Float_t px=particle->GetPx();
    Float_t py=particle->GetPy();
    Float_t pz=particle->GetPz();
    Float_t  e=particle->GetEnergy();

//
//  transverse momentum cut    
    Float_t pt=TMath::Sqrt(px*px+py*py);
    if (pt > fPtMax || pt < fPtMin) 
    {
//	printf("\n failed pt cut %f %f %f \n",pt,fPtMin,fPtMax);
	return kFALSE;
    }
//
// momentum cut
    Float_t p=TMath::Sqrt(px*px+py*py+pz*pz);
    if (p > fPMax || p < fPMin) 
    {
//	printf("\n failed p cut %f %f %f \n",p,fPMin,fPMax);
	return kFALSE;
    }
    
//
// theta cut
    Float_t  theta = Float_t(TMath::ATan2(Double_t(pt),Double_t(pz)));
    if (theta > fThetaMax || theta < fThetaMin) 
    {
//	printf("\n failed theta cut %f %f %f \n",theta,fThetaMin,fThetaMax);
	return kFALSE;
    }

//
// rapidity cut
    Float_t y = 0.5*TMath::Log((e+pz)/(e-pz));
    if (y > fYMax || y < fYMin)
    {
//	printf("\n failed y cut %f %f %f \n",y,fYMin,fYMax);
	return kFALSE;
    }

//
// phi cut
    Float_t phi=Float_t(TMath::ATan2(Double_t(py),Double_t(px)))+TMath::Pi();
    if (phi > fPhiMax || phi < fPhiMin)
    {
//	printf("\n failed phi cut %f %f %f \n",phi,fPhiMin,fPhiMax);
	return kFALSE;
    }

    return kTRUE;
}
void AliGenPythia::AdjustWeights()
{
    TClonesArray *PartArray = gAlice->Particles();
    GParticle *Part;
    Int_t ntrack=gAlice->GetNtrack();
    for (Int_t i=0; i<ntrack; i++) {
	Part= (GParticle*) PartArray->UncheckedAt(i);
	Part->SetWgt(Part->GetWgt()*fKineBias);
    }
}









