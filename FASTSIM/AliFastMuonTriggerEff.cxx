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

#include "AliFastMuonTriggerEff.h"
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <stdlib.h>

// Debugging flag
//#define MYTRIGDEBUG

ClassImp(AliFastMuonTriggerEff)

AliFastMuonTriggerEff::AliFastMuonTriggerEff():
    AliFastResponse("Efficiency", "Muon Trigger Efficiency")
{
//
// Default constructor
//
    SetCut(kLow);
    fZones=0;
}


void AliFastMuonTriggerEff::Init()
{
//
//  Initialization
//
    printf("Initializing %s / %s", fName.Data(), fTitle.Data());
//

    fPhiMin   = -90. ;
    fPhiMax   =  90. ;
    fDphi     =   9. ;
    fThetaMin =   2. ;
    fThetaMax =   9. ;
    fDtheta   =   0.7;
    fDpt      =   0.1;

    InitTree();
}

void AliFastMuonTriggerEff::InitTree()
{
//
//  Initialize tables from tree
//
    TTree          *chain;   // pointer to the analyzed TTree or TChain
    Int_t           nSim, nZona;
    Double_t        pmin, pmax, tmin, tmax, phimin, phimax, bkg;
    Double_t        len50, hen50, leff, heff;
    Double_t        vLPt[50];
    Double_t        vHPt[50];
    Char_t file[100]="$(ALICE_ROOT)/FASTSIM/data/vettorpara.root";
//
//  Avoid memory leak in case of reinitialization
    if(fZones!=0) {
        printf("\nWarning: reinitialization of an object of class: AliFastMuonTriggerEff\n");
        for (Int_t i=0; i<fZones; i++) {
            if(fEffLow [i])delete[] fEffLow[i];
	    if(fEffHigh[i])delete[] fEffHigh[i];
        }
	if(fEffLow) {
	    delete[] fEffLow;
	    fEffLow=0;
        }
	if(fEffHigh) {
	    delete[] fEffHigh;
	    fEffHigh=0;
        }
    }
    printf("AliFastMuonTriggerEff: Initialization\n");
    TFile *f = new TFile(file);
    if(f->IsZombie()) {
        printf("Cannot open file: %s\n",file);
        return;
    }
    f->ls();
    
    TTree* tree = (TTree*)gDirectory->Get("fitPar");

//
//
    chain = tree;
    chain->SetMakeClass(1);
    
    chain->SetBranchAddress("nSim",&nSim);
    chain->SetBranchAddress("nZona",&nZona);
    chain->SetBranchAddress("ptmin",&pmin);
    chain->SetBranchAddress("ptmax",&pmax);
    chain->SetBranchAddress("Thetamin",&tmin);
    chain->SetBranchAddress("Thetamax",&tmax);
    chain->SetBranchAddress("Phimin",&phimin);
    chain->SetBranchAddress("Phimax",&phimax);
    chain->SetBranchAddress("Bkg",&bkg);
    chain->SetBranchAddress("EffLPt",vLPt);
    chain->SetBranchAddress("EffHPt",vHPt);
    chain->SetBranchAddress("Pt0.5Low",&len50);
    chain->SetBranchAddress("Pt0.5High",&hen50);
    chain->SetBranchAddress("EffLowMax",&leff);
    chain->SetBranchAddress("EffHighMax",&heff);
// 
//
    Int_t nentries = Int_t(chain->GetEntries());
    Int_t nbytes = 0, nb = 0;

// Count the number of zones of the parametrization
    Int_t nzone0=0, nzone1=0;
    for (Int_t jentry=0; jentry<nentries; jentry++) {
	nb = chain->GetEntry(jentry);
	if(bkg==0.) {
	     if(nSim==0)nzone0++;
	     if(nSim==1)nzone1++;
	}
    }
    
    printf("Trigger parametrisation for %d zones for Pt: 0. 5. GeV/c\n",nzone0);
    printf("and %d zones extended to 10. GeV/c\n",nzone1);
    fZones=nzone0+nzone1;
//    printf("Ciao\n");
    if(fZones<=0){
        printf("Number of zones must be greater than 0\n");
	exit(6);
    }
    
    fEffLow =new Float_t*[fZones];
    fEffHigh=new Float_t*[fZones];
    for (Int_t i=0; i<fZones; i++) {
        fEffLow [i]=new Float_t[50];
	fEffHigh[i]=new Float_t[50];
    }
    
//  Initialize look-up table to standard values
    Int_t isim, itheta, iphi;
    for (isim=0; isim<2; isim++) {
        for (itheta=0; itheta<10; itheta++) {
	    for (iphi=0; iphi<20; iphi++) {
	        fLook[isim][itheta][iphi]=0;
            }
	}
    }

//  Loading Trigger efficiencies
    Int_t myzone=0;
#ifdef MYTRIGDEBUG
            printf("Filling nSim nZona pmin pmax tmin tmax phimin phimax: ....\n");
#endif
    for (Int_t jentry=0; jentry<nentries; jentry++) {
//	Int_t ientry = LoadTree(jentry); 
	nb = chain->GetEntry(jentry);   
	nbytes += nb;
#ifdef MYTRIGDEBUG
       printf("Getting entry %d... ",jentry);
#endif
// For the time being it works with background 0
	if ((nSim == 0 || nSim == 1)&&bkg==0.) {
#ifdef MYTRIGDEBUG
            printf("Filling %d %d %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f: ",
	    nSim,nZona,pmin,pmax,tmin,tmax,phimin,phimax);
#endif
	    for (Int_t k = 0; k < 50; k++) {
		fEffLow [myzone][k] = vLPt[k];
		fEffHigh[myzone][k] = vHPt[k];
#ifdef MYTRIGDEBUG
		if(k<15)printf(" %5.3f",vLPt[k]);
#endif
	    }
#ifdef MYTRIGDEBUG
            printf("\n");
#endif
	    myzone++;
	    iphi=Int_t(nZona/10)-10;
	    itheta=nZona%10;
	    if(iphi<0||iphi>19||itheta<0||itheta>9) {
	        printf("The format of data file is not consistent\nwith this version of the code\n");
	        printf("This should never happen: iphi %d, itheta: %d\n",iphi,itheta);
		exit(7);
	    }
	    fLook[nSim][itheta][iphi]=myzone;
	} else {
	    printf("Skipping entry with nSim=%d, bkg=%f\n",nSim,bkg);
	}
    }
#ifdef MYTRIGDEBUG
    printf("This is the content of the LUT after first step of initialization\n");
    for(isim=0; isim<2; isim++) {
        printf("isim=%d\n",isim);
        for(iphi=0; iphi<20; iphi++) {
	    printf("iphi: %2d:",iphi);
	    for(itheta=0; itheta<10; itheta++) {
	        printf(" %4d",fLook[isim][itheta][iphi]);
            }
	    printf("\n");
	}
    }
#endif
// Filling look=up table for the zones where the extended simulation does
// not exists
    for(iphi=0; iphi<20; iphi++) {
    	for(itheta=0; itheta<10; itheta++) {
	    if(fLook[0][itheta][iphi]==0) {
	        printf("Missing entry isim=%d itheta=%d iphi=%d\n",isim,itheta,iphi);
		exit(8);
	    }
	    if(fLook[0][itheta][iphi]<0||fLook[0][itheta][iphi]>fZones) {
	        printf("Problem with entry isim=%d itheta=%d iphi=%d\n",isim,itheta,iphi);
		exit(9);
	    }
	    if(fLook[1][itheta][iphi]==0) {
                 fLook[1][itheta][iphi]=-fLook[0][itheta][iphi];
	    }
	}
    }
#ifdef MYTRIGDEBUG
    for(isim=0; isim<2; isim++) {
        printf("isim=%d\n",isim);
        for(iphi=0; iphi<20; iphi++) {
	    printf("iphi: %2d:",iphi);
	    for(itheta=0; itheta<10; itheta++) {
	        printf(" %4d",fLook[isim][itheta][iphi]);
            }
	    printf("\n");
	}
    }
        for(iphi=0; iphi<20; iphi++) {
	    for(itheta=0; itheta<10; itheta++) {
    for(isim=0; isim<2; isim++) {
	        printf("%d %2d %2d %4d:",isim,iphi,itheta,fLook[isim][itheta][iphi]);
		if(fLook[isim][itheta][iphi]>0) {
		    myzone=fLook[isim][itheta][iphi]-1;
		    for(Int_t ipt=0; ipt<20; ipt++) {
		        //printf(" %5.3f",fEffLow[myzone][ipt]);
			printf(" %5.3f",fEffHigh[myzone][ipt]);
		    }
		    printf(" ...");
		    for(Int_t ipt=40; ipt<50; ipt++) {
		        //printf(" %5.3f",fEffLow[myzone][ipt]);
			printf(" %5.3f",fEffHigh[myzone][ipt]);
		    }		    
		}
		printf("\n");
            }
	}
    }    
#endif
    f->Close();
}

void AliFastMuonTriggerEff::Evaluate(Float_t charge, Float_t pt,
                Float_t theta, Float_t phi,Float_t& effLow, Float_t& effHigh,
		Float_t& /*eff*/)
{
//
//  Trigger efficiency for pt, theta, phi (low and high cut)
//
    effLow=0.;
    effHigh=0.;
    if(fZones==0) {
        printf("Call to uninitialized object of class: AliFastMuonTriggerEff\n");
	return;
    }
    if(pt<0) {
        printf("Warning: pt: %f < 0. GeV/c\n",pt);
	return;	
    }
    Int_t iPt   = Int_t(pt/fDpt)%50;
    Int_t iSim  = Int_t(pt/fDpt)/50;
    Int_t iPhi  = Int_t((phi-fPhiMin)/fDphi);
    if(phi<fPhiMin)iPhi=iPhi-1;
    Int_t iTheta = Int_t((theta-fThetaMin)/fDtheta);
#ifdef MYTRIGDEBUG
    printf("iSim iPt iTheta iPhi: %d %d %d %d\n",iSim,iPt,iTheta,iPhi);
#endif
    iPhi=iPhi-40*(iPhi/40);
#ifdef MYTRIGDEBUG
    printf("1: iPhi converted to: %d for angle equivalence\n",iPhi);
#endif
    if(iPhi<0)iPhi=-iPhi-1;
    if(iPhi>19)iPhi=39-iPhi;
#ifdef MYTRIGDEBUG
    printf("2: iPhi converted to: %d for the symmetry of the spectrometer\n",iPhi);
#endif
    if(charge==1.){
    } else if(charge==-1.) {
    iPhi=19-iPhi;
#ifdef MYTRIGDEBUG
    printf("3: iPhi converted to: %d for charge symmetry\n",iPhi);
#endif
    } else {
        printf("Warning: not understand charge: %f\n",charge);
        return;
    }
    if(iTheta<0||iTheta>9) {
        printf("Warning: theta: %f outside acceptance\n",theta);
        return;
    }
    if(iPt<0) {
        printf("Warning: what do you mean with pt: %f <0?\n",pt);
	return;
    }
    if(iSim>=fSim) {
        iSim=fSim-1;
	iPt=49;
#ifdef MYTRIGDEBUG
    printf("4: iSim iPt converted to: %d %d (last zone)\n",iSim,iPt);
#endif
    }
    Int_t iLook=fLook[iSim][iTheta][iPhi];
    if(iLook<0) {
#ifdef MYTRIGDEBUG
    printf("5: iLook iPt: %d %d converted to: ",iLook,iPt);
#endif
        iLook=-iLook-1;
	iPt=49;
#ifdef MYTRIGDEBUG
    printf("%d %d from look up table contents\n",iLook,iPt);
#endif
    } else {
        iLook=iLook-1;
    }
    effLow=fEffLow[iLook][iPt];
    effHigh=fEffHigh[iLook][iPt];
#ifdef MYTRIGDEBUG
    printf("6: charge, iSim, iTheta, iPhi, iPt: %f %d %d %d %d effLow: %f, effHigh: %f\n",
               charge,iSim,iTheta,iPhi,iPt,effLow,effHigh);
#endif
    
    //fEffLow [iPhi][iTheta][iPt];
    //fEffHigh[iPhi][iTheta][iPt];    
}

Float_t AliFastMuonTriggerEff::Evaluate(Float_t charge, Float_t pt,
                   Float_t theta, Float_t phi)
{
    //
    // Trigger efficiency for pt, theta, phi depending of fCut
    // 
    if(fZones==0) {
        printf("Call to uninitialized object of class: AliFastMuonTriggerEff\n");
	return 0.;
    }
    Float_t eff;
    Float_t effLow, effHigh;
    
    Evaluate(charge,pt,theta,phi,effLow,effHigh);
    if (fCut == kLow) 
	eff  = effLow;
    else
	eff  = effHigh;

    return eff;
}
