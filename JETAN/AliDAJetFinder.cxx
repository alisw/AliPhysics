
//  **************************************************************************
//  * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
//  *                                                                        *
//  * Author: The ALICE Off-line Project.                                    *
//  * Contributors are mentioned in the code where appropriate.              *
//  *                                                                        *
//  * Permission to use, copy, modify and distribute this software and its   *
//  * documentation strictly for non-commercial purposes is hereby granted   *
//  * without fee, provided that the above copyright notice appears in all   *
//  * copies and that both the copyright notice and this permission notice   *
//  * appear in the supporting documentation. The authors make no claims     *
//  * about the suitability of this software for any purpose. It is          *
//  * provided "as is" without express or implied warranty.                  *
//  **************************************************************************

//-----------------------------------------------------------------------------------
// Jet finder based on Deterministic Annealing
// For further informations about the DA working features see:
// Phys.Lett. B601 (2004) 56-63 (http://arxiv.org/abs/hep-ph/0407214)
// Author: Davide Perrino (davide.perrino@ba.infn.it, davide.perrino@cern.ch)
//-----------------------------------------------------------------------------------

#include <TMath.h>
#include <TRandom.h>
#include <TClonesArray.h>
#include "AliJetReader.h"
#include "AliDAJetHeader.h"
#include "AliDAJetFinder.h"


ClassImp(AliDAJetFinder)


//-----------------------------------------------------------------------------------
AliDAJetFinder::AliDAJetFinder():
	fAlpha(1.01),
	fDelta(1e-8),
	fAvDist(1e-6),
	fEps(1e-4),
	fEpsMax(1e-2),
	fNloopMax(100),
	fBeta(0.1),
	fNclustMax(0),
	fPyx(0x0),
	fY(0x0),
	fPx(0x0),
	fPy(0x0),
	fXEta(0x0),
	fXPhi(0x0),
	fNin(0)
{
	// Constructor
}

//-----------------------------------------------------------------------------------
AliDAJetFinder::~AliDAJetFinder()
{
	// Destructor
	delete fPyx;
	delete fY;
	delete fPx;
	delete fPy;
	delete [] fXEta;
	delete [] fXPhi;
}

//-----------------------------------------------------------------------------------
void AliDAJetFinder::FindJets()  
{
// Find the jets in current event
// 
	Float_t betaStop=100.;

	Double_t dEtSum=0;
	InitDetAnn(dEtSum);
	if (!fNin) return;

	Int_t nc=1,nk;
	DoubleClusters(nc,nk);
	do{					//loop over beta
		fBeta*=fAlpha;
		Annealing(nk);
		NumCl(nc,nk);
	}while((fBeta<betaStop || nc<4) && nc<fNclustMax);

	Int_t *xx=new Int_t[fNin];
	EndDetAnn(nk,xx,dEtSum);
	StoreJets(nk,xx);
	delete [] xx;

}

//-----------------------------------------------------------------------------------
void AliDAJetFinder::InitDetAnn(Double_t &dEtSum)
{
//Initialise the variables used by the algorithm
	fBeta=0.1;
	TClonesArray *lvArray = fReader->GetMomentumArray();
	fNin = lvArray->GetEntries();
	fNclustMax= ((AliDAJetHeader*)fHeader)->GetFixedCl() ? 
	    ((AliDAJetHeader*)fHeader)->GetNclustMax() 
	    : 
	    TMath::Max((Int_t)TMath::Sqrt(fNin),5);
	fXEta=new Double_t[fNin]; fXPhi=new Double_t[fNin];
	fPx = new TVectorD(fNin);
	for (Int_t iIn=0; iIn<fNin; iIn++){
		TLorentzVector *lv=(TLorentzVector*)lvArray->At(iIn);
		fXEta[iIn] = lv->Eta();
		fXPhi[iIn] = lv->Phi()<0 ? lv->Phi() + 2*TMath::Pi() : lv->Phi();
		(*fPx)(iIn)=lv->Pt();
		dEtSum+=(*fPx)(iIn);
	}
	for (Int_t iIn=0; iIn<fNin; iIn++) (*fPx)(iIn)=(*fPx)(iIn)/dEtSum;

	Int_t njdim=2*fNclustMax+1;
	fPyx = new TMatrixD(fNin,njdim);
	fY = new TMatrixD(4,njdim);
	fPy= new TVectorD(njdim);
	fY->Zero();fPyx->Zero();fPy->Zero();
	(*fPy)(0)=1;
	TMatrixDColumn(*fPyx,0)=1;
	Double_t ypos=0,xpos=0;
	for (Int_t iIn=0; iIn<fNin; iIn++){
		(*fY)(0,0)+=(*fPx)(iIn)*fXEta[iIn];
		ypos+=(*fPx)(iIn)*TMath::Sin(fXPhi[iIn]);
		xpos+=(*fPx)(iIn)*TMath::Cos(fXPhi[iIn]);
	}
	(*fY)(1,0)=(atan2(ypos,xpos)>0) ? atan2(ypos,xpos) : atan2(ypos,xpos)+2*TMath::Pi();
	lvArray->Delete();
}

//-----------------------------------------------------------------------------------
void AliDAJetFinder::DoubleClusters(Int_t nc,Int_t &nk)
{
	for(Int_t iClust=0; iClust<nc; iClust++){
		(*fPy)(iClust)=(*fPy)(iClust)/2;
		(*fPy)(nc+iClust)=(*fPy)(iClust);
		for(Int_t iComp=0; iComp<3; iComp++) (*fY)(iComp,nc+iClust)=(*fY)(iComp,iClust);
	}
	nk=2*nc;
}

//-----------------------------------------------------------------------------------
void AliDAJetFinder::Annealing(Int_t nk)
{
// Main part of the algorithm
	const Double_t pi=TMath::Pi();
	TVectorD *py = new TVectorD(nk);
	TVectorD *p  = new TVectorD(nk);
	TMatrixD *y  = new TMatrixD(4,nk);
	TMatrixD *y1 = new TMatrixD(4,nk);
	TMatrixD *ry = new TMatrixD(2,nk);
	Double_t Dist(TVectorD,TVectorD);

	Double_t df[2]={((AliDAJetHeader*)fHeader)->GetEtaCut(),pi};
	TVectorD vPart(2);
	Double_t *m = new Double_t[nk];
	Double_t chi,chi1;
	do{
		Int_t nloop=0;
		for (Int_t iClust=0; iClust<nk; iClust++){
			for (Int_t i=0; i<3; i++)(*y1)(i,iClust)=(*fY)(i,iClust);
			(*py)(iClust)=(*fPy)(iClust);
		}
	//perturbation of codevectors
		Double_t seed=1000000*gRandom->Rndm(24);
		ry->Randomize(-0.5,0.5,seed);
		for (Int_t i=0; i<2; i++){
			for (Int_t iClust=0; iClust<nk/2; iClust++)
				(*y1)(i,iClust)+=((*ry)(i,iClust)+TMath::Sign(0.5,(*ry)(i,iClust)))*fDelta*df[i];
			for (Int_t iClust=nk/2; iClust<nk; iClust++)
				(*y1)(i,iClust)-=((*ry)(i,iClust-nk/2)+TMath::Sign(0.5,(*ry)(i,iClust-nk/2)))*fDelta*df[i];
		}
		do{
	//recalculate conditional probabilities
			nloop++;
			for (Int_t iIn=0; iIn<fNin; iIn++){
				vPart(0)=fXEta[iIn]; vPart(1)=fXPhi[iIn];
				for(Int_t iClust=0; iClust<nk; iClust++){
					(*fPyx)(iIn,iClust)=-log((*py)(iClust))+fBeta*Dist(vPart,TMatrixDColumn(*y1,iClust));
					m[iClust]=(*fPyx)(iIn,iClust);
				}
				Double_t pyxNorm=0;
				Double_t minPyx=TMath::MinElement(nk,m);
				for (Int_t iClust=0; iClust<nk; iClust++){
					(*fPyx)(iIn,iClust)-=minPyx;
					(*fPyx)(iIn,iClust)=exp(-(*fPyx)(iIn,iClust));
					pyxNorm+=(*fPyx)(iIn,iClust);
				}
				for (Int_t iClust=0; iClust<nk; iClust++) (*fPyx)(iIn,iClust)/=pyxNorm;
			}
			p->Zero();
			y->Zero();
	//recalculate codevectors
			for (Int_t iClust=0; iClust<nk; iClust++){
				Double_t xpos=0,ypos=0,pxy;
				for (Int_t iIn=0; iIn<fNin; iIn++) (*p)(iClust)+=(*fPx)(iIn)*(*fPyx)(iIn,iClust);
				for (Int_t iIn=0; iIn<fNin; iIn++){
					pxy=(*fPx)(iIn)*(*fPyx)(iIn,iClust)/(*p)(iClust);
					ypos+=pxy*TMath::Sin(fXPhi[iIn]);
					xpos+=pxy*TMath::Cos(fXPhi[iIn]);
					(*y)(0,iClust)+=pxy*fXEta[iIn];
				}
				(*y)(1,iClust)=(atan2(ypos,xpos)>0) ? atan2(ypos,xpos) : atan2(ypos,xpos)+2*pi;
			}
	//verify codevectors' stability
			chi=0;
			for (Int_t iClust=0; iClust<nk; iClust++){
				chi1=TMath::CosH((*y1)(0,iClust)-(*y)(0,iClust))-TMath::Cos((*y1)(1,iClust)-(*y)(1,iClust));
				chi1/=(2*TMath::CosH((*y1)(0,iClust))*TMath::CosH((*y)(0,iClust)));
				chi1*=chi1;
				if (chi1>chi) chi=chi1;
			}
			chi=TMath::Sqrt(chi);
			for (Int_t iClust=0; iClust<nk; iClust++){
				for (Int_t i=0; i<2; i++) (*y1)(i,iClust)=(*y)(i,iClust);
				(*py)(iClust)=(*p)(iClust);
			}
			if (nloop>fNloopMax){
				if (chi<fEpsMax || nloop>500) break;
			}
		}while (chi>fEps);
	}while (chi>fEpsMax);
	for (Int_t iClust=0; iClust<nk; iClust++){				//set codevectors and probability equal to those calculated
		for (Int_t i=0; i<2; i++) (*fY)(i,iClust)=(*y)(i,iClust);
		(*fPy)(iClust)=(*p)(iClust);
	}
    delete py;
    delete p;
    delete y;
    delete y1;
    delete ry;
	delete [] m;
}

//-----------------------------------------------------------------------------------
void AliDAJetFinder::NumCl(Int_t &nc,Int_t &nk)
{
	static Bool_t growcl=true;
	
	if (nk==2) growcl=true;
	if (growcl){
//verify if two codevectors are equal within fAvDist
		Int_t *nSame = new Int_t[nk];
		Int_t **iSame = new Int_t*[nk];
		Int_t **cont = new Int_t*[nk];
		for (Int_t iClust=0; iClust<nk; iClust++) cont[iClust]=new Int_t[nk],iSame[iClust]=new Int_t[nk];
		for (Int_t iClust=0; iClust<nk; iClust++){
			iSame[iClust][iClust]=1;
			for (Int_t iClust1=iClust+1; iClust1<nk; iClust1++){
				Double_t eta  = (*fY)(0,iClust) ; Double_t phi  = (*fY)(1,iClust);
				Double_t eta1 = (*fY)(0,iClust1); Double_t phi1 = (*fY)(1,iClust1);
				Double_t distCl=(TMath::CosH(eta-eta1)-TMath::Cos(phi-phi1))/(2*TMath::CosH(eta)*TMath::CosH(eta1));
				if (distCl < fAvDist) iSame[iClust][iClust1]=iSame[iClust1][iClust]=1;
			}
		}
		ReduceClusters(iSame,nk,nc,cont,nSame);
		if (nc >= fNclustMax) growcl=false;
//recalculate the nc distinct codevectors
		TMatrixD *pyx = new TMatrixD(fNin,2*nc);
		TVectorD *py = new TVectorD(nk);
		TMatrixD *y1  = new TMatrixD(3,nk);
		for (Int_t iClust=0; iClust<nc; iClust++){
			for(Int_t j=0; j<nSame[iClust]; j++){
				Int_t iClust1 = cont[iClust][j];
				for (Int_t iIn=0; iIn<fNin; iIn++) (*pyx)(iIn,iClust)+=(*fPyx)(iIn,iClust1);
				(*py)(iClust)+=(*fPy)(iClust1);
				for (Int_t i=0; i<2; i++) (*y1)(i,iClust)+=(*fY)(i,iClust1);
			}
			for (Int_t i=0; i<2; i++) (*y1)(i,iClust)/=nSame[iClust];
		}
		if (nc > nk/2){
			for (Int_t iClust=0; iClust<nc; iClust++){
				for (Int_t iIn=0; iIn<fNin; iIn++) (*fPyx)(iIn,iClust)=(*pyx)(iIn,iClust);
				for (Int_t iComp=0; iComp<2; iComp++) (*fY)(iComp,iClust)=(*y1)(iComp,iClust);
				(*fPy)(iClust)=(*py)(iClust);
			}
			nk=nc;
			if (growcl) DoubleClusters(nc,nk);
		}
		delete [] nSame;
		delete [] iSame;
		delete [] cont;
		delete pyx;
		delete py;
		delete y1;
	}

}

//-----------------------------------------------------------------------------------
void AliDAJetFinder::ReduceClusters(Int_t **iSame,Int_t nc,Int_t &ncout,Int_t **cont,Int_t *nSameOut)
{
	Int_t *nSame = new Int_t[nc];
	Int_t *iperm = new Int_t[nc];
	Int_t *go = new Int_t[nc];
	for (Int_t iCl=0; iCl<nc; iCl++){
		nSame[iCl]=0;
		for (Int_t jCl=0; jCl<nc; jCl++) nSame[iCl]+=iSame[iCl][jCl];
		iperm[iCl]=iCl;
		go[iCl]=1;
	}
	TMath::Sort(nc,nSame,iperm,true);
	Int_t l=0;
	for (Int_t iCl=0; iCl<nc; iCl++){
		Int_t k=iperm[iCl];
		if (go[k] == 1){
			Int_t m=0;
			for (Int_t jCl=0; jCl<nc; jCl++){
				if (iSame[k][jCl] == 1){
					cont[l][m]=jCl;
					go[jCl]=0;
					m++;
				}
			}
			nSameOut[l]=m;
			l++;
		}
	}
	ncout=l;
	delete [] nSame;
	delete [] iperm;
	delete [] go;
}

//-----------------------------------------------------------------------------------
void AliDAJetFinder::EndDetAnn(Int_t &nk,Int_t *xx,Double_t etx)
{
//now assign each particle to only one cluster
	Double_t *clusters=new Double_t[nk];
	for (Int_t iIn=0; iIn<fNin; iIn++){
		for (Int_t iClust=0; iClust<nk; iClust++) clusters[iClust]=(*fPyx)(iIn,iClust);
		xx[iIn]=TMath::LocMax(nk,clusters);
	}
    delete [] clusters;
	
//recalculate codevectors, having all p(y|x)=0 or 1
	fY->Zero();
	fPyx->Zero();
	fPy->Zero();
	for (Int_t iIn=0; iIn<fNin; iIn++){
		Int_t iClust=xx[iIn];
		(*fPyx)(iIn,iClust)=1;
		(*fPy)(iClust)+=(*fPx)(iIn);
		(*fY)(0,iClust)+=(*fPx)(iIn)*fXEta[iIn];
		(*fY)(3,iClust)+=(*fPx)(iIn)*etx;
	}
	Int_t k=0;
	for (Int_t iClust=0; iClust<nk; iClust++){
		if ((*fPy)(iClust)>0){
			Double_t xpos=0,ypos=0,pxy;
			for (Int_t iIn=0; iIn<fNin; iIn++){
				pxy=(*fPx)(iIn)*(*fPyx)(iIn,iClust)/(*fPy)(iClust);
				ypos+=pxy*TMath::Sin(fXPhi[iIn]);
				xpos+=pxy*TMath::Cos(fXPhi[iIn]);
				if (xx[iIn]==iClust) xx[iIn]=k;
			}
			(*fY)(0,k)=(*fY)(0,iClust)/(*fPy)(iClust);
			(*fY)(1,k)=(atan2(ypos,xpos)>0) ? atan2(ypos,xpos) : atan2(ypos,xpos)+2*TMath::Pi();
			(*fY)(3,k)=(*fY)(3,iClust);
			k++;
		}
	}
	nk=k;
}

//-----------------------------------------------------------------------------------
void AliDAJetFinder::StoreJets(Int_t nk,Int_t *xx)
{
//evaluate significant clusters properties
	const Double_t pi=TMath::Pi();
	Double_t dMeanDist=TMath::Sqrt(4*((AliDAJetHeader*)fHeader)->GetEtaCut()*pi/fNin);
	Bool_t   *isJet = new Bool_t[nk];
	Double_t *etNoBg= new Double_t[nk];
	Double_t *dDeltaEta=new Double_t[nk];
	Double_t *dDeltaPhi=new Double_t[nk];
	Double_t *surf  = new Double_t[nk];
	Double_t *etDens= new Double_t[nk];
	for (Int_t iClust=0; iClust<nk; iClust++){									//clusters loop
		isJet[iClust]=false;
		Double_t dEtaMin=10.,dEtaMax=-10.,dPhiMin=10.,dPhiMax=0.;
		for (Int_t iIn=0; iIn<fNin; iIn++){
			if (xx[iIn]!=iClust) continue;
			if (fXEta[iIn] < dEtaMin) dEtaMin=fXEta[iIn];
			if (fXEta[iIn] > dEtaMax) dEtaMax=fXEta[iIn];
			Double_t dPhi=fXPhi[iIn]-(*fY)(1,iClust);
			if      (dPhi > pi     ) dPhi-=2*pi;
			else if (dPhi < (-1)*pi) dPhi+=2*pi;
			if      (dPhi < dPhiMin) dPhiMin=dPhi;
			else if (dPhi > dPhiMax) dPhiMax=dPhi;
		}
		dDeltaEta[iClust]=dEtaMax-dEtaMin+dMeanDist;
		dDeltaPhi[iClust]=dPhiMax-dPhiMin+dMeanDist;
		surf[iClust]=0.25*pi*dDeltaEta[iClust]*dDeltaPhi[iClust];
		etDens[iClust]=(*fY)(3,iClust)/surf[iClust];
	}

	if (((AliDAJetHeader*)fHeader)->GetSelJets()){
		for (Int_t iClust=0; iClust<nk; iClust++){
			if (!isJet[iClust]){
				Double_t etDensMed=0.;
				Double_t etDensSqr=0.;
				Int_t norm=0;
				for (Int_t iClust1=0; iClust1<nk; iClust1++){
					if(iClust1!=iClust){
						etDensMed+=etDens[iClust1];
						etDensSqr+=TMath::Power(etDens[iClust1],2);
						norm++;
					}
				}
				etDensMed/=TMath::Max(norm,1);
				etDensSqr/=TMath::Max(norm,1);
				Double_t deltaEtDens=TMath::Sqrt(etDensSqr-TMath::Power(etDensMed,2));
				if ((*fY)(3,iClust) > (etDensMed+deltaEtDens)*surf[iClust]) isJet[iClust]=kTRUE;
				etNoBg[iClust]=(*fY)(3,iClust)-etDensMed*surf[iClust];
			}
		}
		for (Int_t iClust=0; iClust<nk; iClust++){
			if (isJet[iClust]){
				Double_t etDensMed=0;
				Double_t extSurf=4*((AliDAJetHeader*)fHeader)->GetEtaCut()*pi;
				for (Int_t iClust1=0; iClust1<nk; iClust1++){
					if (!isJet[iClust1]) etDensMed+=(*fY)(3,iClust1);
					else extSurf-=surf[iClust1];
				}
				etDensMed/=extSurf;
				etNoBg[iClust]=(*fY)(3,iClust)-etDensMed*surf[iClust];
				if (etNoBg[iClust]<((AliDAJetHeader*)fHeader)->GetEtMin()){
					isJet[iClust]=kFALSE;
					iClust=-1;
				}
			}
		}
	} else {
		for (Int_t iClust=0; iClust<nk; iClust++) isJet[iClust]=true;
	}
	delete [] etDens;
	delete [] surf;
	
//now add selected jets to the list
	Int_t *inJet = new Int_t[fNin];
	TRefArray *refs = 0;
	Bool_t fromAod = !strcmp(fReader->ClassName(),"AliJetAODReader");
	if (fromAod) refs = fReader->GetReferences();
	for (Int_t iClust=0; iClust<nk; iClust++){									//clusters loop
		if (isJet[iClust]){														//choose cluster
			Float_t px,py,pz,en;
			px = (*fY)(3,iClust)*TMath::Cos((*fY)(1,iClust));
			py = (*fY)(3,iClust)*TMath::Sin((*fY)(1,iClust));
			pz = (*fY)(3,iClust)/TMath::Tan(2.0 * TMath::ATan(TMath::Exp(-(*fY)(0,iClust))));
			en = TMath::Sqrt(px * px + py * py + pz * pz);
			AliAODJet jet(px, py, pz, en);
			if (fromAod) 
			    for (Int_t iIn=0; iIn<fNin; iIn++) if (xx[iIn]==iClust) jet.AddTrack(refs->At(iIn));
			AddJet(jet);
			printf("jet %d, Eta: %f, Phi: %f, Et: %f\n",iClust,jet.Eta(),jet.Phi(),jet.Pt());
		}
	}
	delete [] dDeltaEta; delete [] dDeltaPhi;
	delete [] etNoBg;
	delete [] isJet;
	delete [] inJet;
}

//-----------------------------------------------------------------------------------
Double_t Dist(TVectorD x,TVectorD y)
{
// Squared distance
	const Double_t pi=TMath::Pi();
	Double_t dphi=TMath::Abs(x(1)-y(1));
	if (dphi > pi) dphi=2*pi-dphi;
	Double_t dist=pow(x(0)-y(0),2)+pow(dphi,2);
	return dist;
}
