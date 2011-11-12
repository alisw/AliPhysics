
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
#include <TRandom2.h>
#include <TClonesArray.h>
#include "AliJetReaderHeader.h"
#include "AliJetReader.h"
#include "AliDAJetHeader.h"
#include "AliDAJetFinder.h"

ClassImp(AliDAJetFinder)


//-----------------------------------------------------------------------------------
AliDAJetFinder::AliDAJetFinder():
	AliJetFinder(),
	fAlpha(1.01),
	fDelta(1e-8),
	fAvDist(1e-6),
	fEps(1e-4),
	fEpsMax(1e-2),
	fNloopMax(100),
	fBeta(0.1),
	fNclustMax(0),
	fNin(0),
	fNeff(0)
{
	// Constructor
}

//-----------------------------------------------------------------------------------
AliDAJetFinder::~AliDAJetFinder()
{
	// Destructor
}

//-----------------------------------------------------------------------------------
void AliDAJetFinder::FindJets()  
{
// Find the jets in current event
// 
	Float_t betaStop=100.;
	fDebug = fHeader->GetDebug();

	Double_t dEtSum=0;
	Double_t *xData[2];
	TVectorD *vPx = new TVectorD();
	TVectorD *vPy = new TVectorD();
	TMatrixD *mPyx= new TMatrixD();
	TMatrixD *mY  = new TMatrixD();
	InitDetAnn(dEtSum,xData,vPx,vPy,mPyx,mY);
	if (fNin < fNclustMax){
	  delete [] xData[0], delete [] xData[1];
	  delete vPx;
	  delete vPy;
	  delete mPyx;
	  delete mY;
	  return;
	}
	Int_t nc=1, nk;
	DoubleClusters(nc,nk,vPy,mY);
	do{					//loop over beta
		fBeta*=fAlpha;
		Annealing(nk,xData,vPx,vPy,mPyx,mY);
		NumCl(nc,nk,vPy,mPyx,mY);
	}while((fBeta<betaStop || nc<4) && nc<fNclustMax);

	Int_t *xx=new Int_t[fNeff];
	for (Int_t i = 0; i < fNeff; i++) xx[i] = 0;
	
	EndDetAnn(nk,xData,xx,dEtSum,vPx,vPy,mPyx,mY);
	StoreJets(nk,xData,xx,mY);
	delete [] xx;

	delete [] xData[0], delete [] xData[1];
	delete mPyx;
	delete mY;
	delete vPx;
	delete vPy;

}

//-----------------------------------------------------------------------------------
void AliDAJetFinder::InitDetAnn(Double_t &dEtSum,Double_t **xData,TVectorD *vPx,TVectorD *vPy,TMatrixD *mPyx,TMatrixD *mY)
{
//Initialise the variables used by the algorithm
	fBeta=0.1;
	fNclustMax = ((AliDAJetHeader*)fHeader)->GetFixedCl() ? 
	    ((AliDAJetHeader*)fHeader)->GetNclustMax() : 
	    TMath::Max((Int_t)TMath::Sqrt(fNin),5);
	Float_t etaEff = ((AliDAJetHeader*)fHeader)->GetEtaEff();
	TClonesArray *lvArray = fReader->GetMomentumArray();
	Int_t nEntr = lvArray->GetEntries();
	fNin=0;
	for (Int_t iEn=0; iEn<nEntr; iEn++) if (fReader->GetCutFlag(iEn)==1) fNin++;

	fNeff = ((AliDAJetHeader*)fHeader)->GetNeff();
	fNeff = TMath::Max(fNeff,fNin);
	Double_t *xEta = new Double_t[fNeff];
	Double_t *xPhi = new Double_t[fNeff];
	xData[0]=xEta; xData[1]=xPhi;
	vPx->ResizeTo(fNeff);
	Int_t iIn=0;
	for (Int_t iEn=0; iEn<nEntr; iEn++){
		if (fReader->GetCutFlag(iEn)==0) continue;
		TLorentzVector *lv=(TLorentzVector*)lvArray->At(iEn);
		xEta[iIn] = lv->Eta();
		xPhi[iIn] = lv->Phi()<0 ? lv->Phi() + 2*TMath::Pi() : lv->Phi();
		(*vPx)(iIn)=lv->Pt();
		dEtSum+=(*vPx)(iIn);
		iIn++;
	}
	TRandom2 r;
	r.SetSeed(0);
	for (iIn=fNin; iIn<fNeff; iIn++){
		xEta[iIn]=r.Uniform(-1*etaEff,etaEff);
		xPhi[iIn]=r.Uniform(0.,2*TMath::Pi());
		(*vPx)(iIn)=r.Uniform(0.01,0.02);
		dEtSum+=(*vPx)(iIn);
	}
	for (iIn=0; iIn<fNeff; iIn++) (*vPx)(iIn)=(*vPx)(iIn)/dEtSum;

	Int_t njdim=2*fNclustMax+1;
	mPyx->ResizeTo(fNeff,njdim);
	mY->ResizeTo(4,njdim);
	vPy->ResizeTo(njdim);
	mY->Zero();mPyx->Zero();vPy->Zero();
	(*vPy)(0)=1;
	TMatrixDColumn(*mPyx,0)=1;
	Double_t ypos=0,xpos=0;
	for (iIn=0; iIn<fNeff; iIn++){
		(*mY)(0,0)+=(*vPx)(iIn)*xEta[iIn];
		ypos+=(*vPx)(iIn)*TMath::Sin(xPhi[iIn]);
		xpos+=(*vPx)(iIn)*TMath::Cos(xPhi[iIn]);
	}
	(*mY)(1,0)=(atan2(ypos,xpos)>0) ? atan2(ypos,xpos) : atan2(ypos,xpos)+2*TMath::Pi();
}

//-----------------------------------------------------------------------------------
void AliDAJetFinder::DoubleClusters(Int_t nc,Int_t &nk,  TVectorD *vPy, TMatrixD *mY) const
{
// Return double clusters
	for(Int_t iClust=0; iClust<nc; iClust++){
		(*vPy)(iClust)=(*vPy)(iClust)/2;
		(*vPy)(nc+iClust)=(*vPy)(iClust);
		for(Int_t iComp=0; iComp<3; iComp++) (*mY)(iComp,nc+iClust)=(*mY)(iComp,iClust);
	}
	nk=2*nc;
}

//-----------------------------------------------------------------------------------
void AliDAJetFinder::Annealing(Int_t nk,Double_t **xData,  TVectorD *vPx,  TVectorD *vPy,  TMatrixD *mPyx,  TMatrixD *mY)
{
// Main part of the algorithm
	const Double_t pi=TMath::Pi();
	TVectorD *py = new TVectorD(nk);
	TVectorD *p  = new TVectorD(nk);
	TMatrixD *y  = new TMatrixD(4,nk);
	TMatrixD *y1 = new TMatrixD(4,nk);
	TMatrixD *ry = new TMatrixD(2,nk);
	Double_t *xEta = xData[0];
	Double_t *xPhi = xData[1];
	Double_t Dist(TVectorD,TVectorD);

	Double_t df[2]={fReader->GetReaderHeader()->GetFiducialEtaMax(),pi};
	TVectorD vPart(2);
	Double_t *m = new Double_t[nk];
	Double_t chi,chi1;
	do{
		Int_t nloop=0;
		for (Int_t iClust=0; iClust<nk; iClust++){
			for (Int_t i=0; i<3; i++)(*y1)(i,iClust)=(*mY)(i,iClust);
			(*py)(iClust)=(*vPy)(iClust);
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
			for (Int_t iIn=0; iIn<fNeff; iIn++){
				vPart(0)=xEta[iIn]; vPart(1)=xPhi[iIn];
				for(Int_t iClust=0; iClust<nk; iClust++){
					(*mPyx)(iIn,iClust)=-log((*py)(iClust))+fBeta*Dist(vPart,TMatrixDColumn(*y1,iClust));
					m[iClust]=(*mPyx)(iIn,iClust);
				}
				Double_t pyxNorm=0;
				Double_t minPyx=TMath::MinElement(nk,m);
				for (Int_t iClust=0; iClust<nk; iClust++){
					(*mPyx)(iIn,iClust)-=minPyx;
					(*mPyx)(iIn,iClust)=exp(-(*mPyx)(iIn,iClust));
					pyxNorm+=(*mPyx)(iIn,iClust);
				}
				for (Int_t iClust=0; iClust<nk; iClust++) (*mPyx)(iIn,iClust)/=pyxNorm;
			}
			p->Zero();
			y->Zero();
	//recalculate codevectors
			for (Int_t iClust=0; iClust<nk; iClust++){
				Double_t xpos=0,ypos=0,pxy;
				for (Int_t iIn=0; iIn<fNeff; iIn++) (*p)(iClust)+=(*vPx)(iIn)*(*mPyx)(iIn,iClust);
				for (Int_t iIn=0; iIn<fNeff; iIn++){
					pxy=(*vPx)(iIn)*(*mPyx)(iIn,iClust)/(*p)(iClust);
					ypos+=pxy*TMath::Sin(xPhi[iIn]);
					xpos+=pxy*TMath::Cos(xPhi[iIn]);
					(*y)(0,iClust)+=pxy*xEta[iIn];
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
		for (Int_t i=0; i<2; i++) (*mY)(i,iClust)=(*y)(i,iClust);
		(*vPy)(iClust)=(*p)(iClust);
	}
    delete py;
    delete p;
    delete y;
    delete y1;
    delete ry;
    delete [] m;
}

//-----------------------------------------------------------------------------------
void AliDAJetFinder::NumCl(Int_t &nc,Int_t &nk,TVectorD *vPy,  TMatrixD *mPyx,TMatrixD *mY)
{
    // Number of clusters
	static Bool_t growcl=true;
	
	if (nk==2) growcl=true;
	if (growcl){
//verify if two codevectors are equal within fAvDist
		Int_t *nSame = new Int_t[nk];
		Int_t **iSame = new Int_t*[nk];
		Int_t **cont = new Int_t*[nk];
		for (Int_t iClust=0; iClust<nk; iClust++) {
		    cont[iClust] =new Int_t[nk];
		    iSame[iClust]=new Int_t[nk];
		}
		
		for (Int_t iClust=0; iClust<nk; iClust++){
			iSame[iClust][iClust]=1;
			for (Int_t iClust1=iClust+1; iClust1<nk; iClust1++){
				Double_t eta  = (*mY)(0,iClust) ; Double_t phi  = (*mY)(1,iClust);
				Double_t eta1 = (*mY)(0,iClust1); Double_t phi1 = (*mY)(1,iClust1);
				Double_t distCl=(TMath::CosH(eta-eta1)-TMath::Cos(phi-phi1))/(2*TMath::CosH(eta)*TMath::CosH(eta1));
				if (distCl < fAvDist) iSame[iClust][iClust1]=iSame[iClust1][iClust]=1;
				else iSame[iClust][iClust1]=iSame[iClust1][iClust]=0;
			}
		}
		ReduceClusters(iSame,nk,nc,cont,nSame);
		if (nc >= fNclustMax) growcl=false;
//recalculate the nc distinct codevectors
		TMatrixD *pyx = new TMatrixD(fNeff,2*nc);
		TVectorD *py = new TVectorD(nk);
		TMatrixD *y1  = new TMatrixD(3,nk);
		for (Int_t iClust=0; iClust<nc; iClust++){
			for(Int_t j=0; j<nSame[iClust]; j++){
				Int_t iClust1 = cont[iClust][j];
				for (Int_t iIn=0; iIn<fNeff; iIn++) (*pyx)(iIn,iClust)+=(*mPyx)(iIn,iClust1);
				(*py)(iClust)+=(*vPy)(iClust1);
				for (Int_t i=0; i<2; i++) (*y1)(i,iClust)+=(*mY)(i,iClust1);
			}
			for (Int_t i=0; i<2; i++) (*y1)(i,iClust)/=nSame[iClust];
		}
		for (Int_t iClust=0; iClust<nk; iClust++) delete [] cont[iClust], delete [] iSame[iClust];
		delete [] iSame;
		delete [] cont;
		delete [] nSame;
		if (nc > nk/2){
			for (Int_t iClust=0; iClust<nc; iClust++){
				for (Int_t iIn=0; iIn<fNeff; iIn++) (*mPyx)(iIn,iClust)=(*pyx)(iIn,iClust);
				for (Int_t iComp=0; iComp<2; iComp++) (*mY)(iComp,iClust)=(*y1)(iComp,iClust);
				(*vPy)(iClust)=(*py)(iClust);
			}
			nk=nc;
			if (growcl) DoubleClusters(nc,nk,vPy,mY);
		}
		delete pyx;
		delete py;
		delete y1;
	}

}

//-----------------------------------------------------------------------------------
void AliDAJetFinder::ReduceClusters(Int_t **iSame,Int_t nc,Int_t &ncout,Int_t **cont,Int_t *nSameOut) const
{
// Reduction step
	Int_t *nSame = new Int_t[nc];
	Int_t *iperm = new Int_t[nc];
	Int_t *go = new Int_t[nc];
	for (Int_t iCl=0; iCl<nc; iCl++){
		nSame[iCl]=0;
		for (Int_t jCl=0; jCl<nc; jCl++) nSame[iCl]+=iSame[iCl][jCl], cont[iCl][jCl]=0;
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
void AliDAJetFinder::EndDetAnn(Int_t &nk,Double_t **xData,Int_t *xx,Double_t etx,TVectorD *vPx,TVectorD *vPy,TMatrixD *mPyx,TMatrixD *mY)
{
//now assign each particle to only one cluster
	Double_t *clusters=new Double_t[nk];
	for (Int_t iIn=0; iIn<fNeff; iIn++){
		for (Int_t iClust=0; iClust<nk; iClust++) clusters[iClust]=(*mPyx)(iIn,iClust);
		xx[iIn]=TMath::LocMax(nk,clusters);
	}
	delete [] clusters;
	
//recalculate codevectors, having all p(y|x)=0 or 1
	Double_t *xEta = xData[0];
	Double_t *xPhi = xData[1];
	mY->Zero();
	mPyx->Zero();
	vPy->Zero();
	for (Int_t iIn=0; iIn<fNin; iIn++){
		Int_t iClust=xx[iIn];
		(*mPyx)(iIn,iClust)=1;
		(*vPy)(iClust)+=(*vPx)(iIn);
		(*mY)(0,iClust)+=(*vPx)(iIn)*xEta[iIn];
		(*mY)(3,iClust)+=(*vPx)(iIn)*etx;
	}
	Int_t k=0;
	for (Int_t iClust=0; iClust<nk; iClust++){
		if ((*vPy)(iClust)>0){
			Double_t xpos=0,ypos=0,pxy;
			for (Int_t iIn=0; iIn<fNin; iIn++){
				pxy=(*vPx)(iIn)*(*mPyx)(iIn,iClust)/(*vPy)(iClust);
				ypos+=pxy*TMath::Sin(xPhi[iIn]);
				xpos+=pxy*TMath::Cos(xPhi[iIn]);
				if (xx[iIn]==iClust) xx[iIn]=k;
			}
			(*mY)(0,k)=(*mY)(0,iClust)/(*vPy)(iClust);
			(*mY)(1,k)=(atan2(ypos,xpos)>0) ? atan2(ypos,xpos) : atan2(ypos,xpos)+2*TMath::Pi();
			(*mY)(3,k)=(*mY)(3,iClust);
			k++;
		}
	}
	nk=k;
}

//-----------------------------------------------------------------------------------
void AliDAJetFinder::StoreJets(Int_t nk, Double_t **xData, Int_t *xx, TMatrixD *mY)
{
//evaluate significant clusters properties
	const Double_t pi=TMath::Pi();
	AliJetReaderHeader *rHeader=fReader->GetReaderHeader();
	Float_t dFidEtaMax = rHeader->GetFiducialEtaMax();
	Float_t dFidEtaMin = rHeader->GetFiducialEtaMin();
	Float_t dFiducialEta= dFidEtaMax - dFidEtaMin;
	Double_t *xEta = xData[0];
	Double_t *xPhi = xData[1];
	Int_t nEff = 0;
	for (Int_t i=0; i<fNeff; i++) if (xEta[i]<dFidEtaMax && xEta[i]>dFidEtaMin) nEff++;
	Double_t dMeanDist=TMath::Sqrt(2*dFiducialEta*pi/nEff);
	Bool_t   *isJet = new Bool_t[nk];
	Double_t *etNoBg= new Double_t[nk];
	Double_t *dDeltaEta=new Double_t[nk];
	Double_t *dDeltaPhi=new Double_t[nk];
	Double_t *surf  = new Double_t[nk];
	Double_t *etDens= new Double_t[nk];
	for (Int_t iClust=0; iClust<nk; iClust++){
		isJet[iClust]=false;
		Double_t dEtaMin=10.,dEtaMax=-10.,dPhiMin=10.,dPhiMax=0.;
		for (Int_t iIn=0; iIn<fNeff; iIn++){
			if (xx[iIn]!=iClust || xEta[iIn]>dFidEtaMax || xEta[iIn]<dFidEtaMin) continue;
			if (xEta[iIn] < dEtaMin) dEtaMin=xEta[iIn];
			if (xEta[iIn] > dEtaMax) dEtaMax=xEta[iIn];
			Double_t dPhi=xPhi[iIn]-(*mY)(1,iClust);
			if      (dPhi > pi     ) dPhi-=2*pi;
			else if (dPhi < (-1)*pi) dPhi+=2*pi;
			if      (dPhi < dPhiMin) dPhiMin=dPhi;
			else if (dPhi > dPhiMax) dPhiMax=dPhi;
		}
		dDeltaEta[iClust]=dEtaMax-dEtaMin+dMeanDist;
		dDeltaPhi[iClust]=dPhiMax-dPhiMin+dMeanDist;
		surf[iClust]=0.25*pi*dDeltaEta[iClust]*dDeltaPhi[iClust];
		etDens[iClust]=(*mY)(3,iClust)/surf[iClust];
	}

	if (((AliDAJetHeader*)fHeader)->GetSelJets()){
		for (Int_t iClust=0; iClust<nk; iClust++){
			if (!isJet[iClust] && (*mY)(0,iClust)<dFidEtaMax && (*mY)(0,iClust)>dFidEtaMin){
				Double_t etDensMed=0.;
				Double_t etDensSqr=0.;
				Int_t norm=0;
				for (Int_t iClust1=0; iClust1<nk; iClust1++){
					if(iClust1!=iClust && (*mY)(0,iClust)<dFidEtaMax && (*mY)(0,iClust)>dFidEtaMin){
						etDensMed+=etDens[iClust1];
						etDensSqr+=TMath::Power(etDens[iClust1],2);
						norm++;
					}
				}
				etDensMed/=TMath::Max(norm,1);
				etDensSqr/=TMath::Max(norm,1);
				Double_t deltaEtDens=TMath::Sqrt(etDensSqr-TMath::Power(etDensMed,2));
				if ((*mY)(3,iClust) > (etDensMed+deltaEtDens)*surf[iClust]) isJet[iClust]=kTRUE;
				etNoBg[iClust]=(*mY)(3,iClust)-etDensMed*surf[iClust];
			}
		}
		for (Int_t iClust=0; iClust<nk; iClust++){
			if (isJet[iClust]){
				Double_t etDensMed=0;
				Double_t extSurf=2*dFiducialEta*pi;
				for (Int_t iClust1=0; iClust1<nk; iClust1++){
					if (!isJet[iClust1]) etDensMed+=(*mY)(3,iClust1);
					else extSurf-=surf[iClust1];
				}
				etDensMed/=extSurf;
				etNoBg[iClust]=(*mY)(3,iClust)-etDensMed*surf[iClust];
				if (etNoBg[iClust]<((AliDAJetHeader*)fHeader)->GetEtMin()){
					isJet[iClust]=kFALSE;
					iClust=-1;
				}
			}
		}
	} else {
		for (Int_t iClust=0; iClust<nk; iClust++){
			isJet[iClust]=true;
			etNoBg[iClust]=(*mY)(3,iClust);
		}
	}
	delete [] etDens;
	delete [] surf;
	
//now add selected jets to the list
	Int_t *iSort = new Int_t[nk];
	TMath::Sort(nk,etNoBg,iSort,true);
	Int_t iCl = 0;
	TRefArray *refs = 0;
	Bool_t fromAod = !strcmp(fReader->ClassName(),"AliJetAODReader");
	if (fromAod) refs = fReader->GetReferences();
	for (Int_t iClust=0; iClust<nk; iClust++){									//clusters loop
		iCl=iSort[iClust];
		if (isJet[iCl]){														//choose cluster
			Float_t px,py,pz,en;
			px = (*mY)(3,iCl)*TMath::Cos((*mY)(1,iCl));
			py = (*mY)(3,iCl)*TMath::Sin((*mY)(1,iCl));
			pz = (*mY)(3,iCl)/TMath::Tan(2.0 * TMath::ATan(TMath::Exp(-(*mY)(0,iCl))));
			en = TMath::Sqrt(px * px + py * py + pz * pz);
			AliAODJet jet(px, py, pz, en);
			if (fromAod){
				Int_t iIn=0;
				Int_t nEntr = fReader->GetMomentumArray()->GetEntries();
				for (Int_t iEn=0; iEn<nEntr; iEn++){
					if (fReader->GetCutFlag(iEn)==0) continue;
					if (xx[iIn]==iCl) jet.AddTrack(refs->At(iEn));
					iIn++;
				}
			}
			AddJet(jet);
			if (fDebug > 0) printf("jet %d, Eta: %f, Phi: %f, Et: %f\n",iCl,jet.Eta(),jet.Phi(),jet.Pt());
		}
	}
	delete [] dDeltaEta; delete [] dDeltaPhi;
	delete [] etNoBg;
	delete [] isJet;
	delete [] iSort;
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
