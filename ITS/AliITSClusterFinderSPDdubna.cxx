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

//  Cluster Finder for
//  Silicon
//  Pixel
//  developed by the Dubna group
//
#include <Riostream.h>
#include "AliITSClusterFinderSPDdubna.h"
#include "AliITSMapA1.h"
#include "AliITS.h"
#include "AliITSdigitSPD.h"
#include "AliITSRawClusterSPD.h"
#include "AliITSRecPoint.h"
#include "AliITSsegmentation.h"
#include "AliRun.h"

//#define DEBUG

ClassImp(AliITSClusterFinderSPDdubna)

//----------------------------------------------------------------------
AliITSClusterFinderSPDdubna::AliITSClusterFinderSPDdubna(
           AliITSsegmentation *seg, TClonesArray *digits, TClonesArray *recp){
    // constructor

    fSegmentation = seg;
    fDigits       = digits;
    fClusters     = recp;
    fNclusters    = 0;
    fMap          = new AliITSMapA1(fSegmentation,fDigits);
    SetDx();
    SetDz();
    SetNCells();
}
//_______________________________________________________________________
AliITSClusterFinderSPDdubna::AliITSClusterFinderSPDdubna(){
    // constructor

    fSegmentation = 0;
    fDigits       = 0;
    fClusters     = 0;
    fNclusters    = 0;
    fMap          = 0;
    SetDx();
    SetDz();
    SetNCells();
}
//______________________________________________________________________
AliITSClusterFinderSPDdubna &AliITSClusterFinderSPDdubna::operator=(const
                                          AliITSClusterFinderSPDdubna &s){
    // The = operator for the class AliITSClusterFinderSPDdugna
    // Inputs:
    //    AliITSClusterFinderSPDdubna @s  The AliITSClusterFinderSPDdubna
    //                                    class to be copy constructed.
    // Outputs:
    //    none.
    // Returned:
    //    this, a copy of the class s.

    this->fNclusters = s.fNclusters;
    this->fDz        = s.fDz;
    this->fDx        = s.fDx;
    this->fMinNCells = s.fMinNCells;
    this->fClusters  = new TClonesArray(*(s.fClusters));
    return *this;
}
//______________________________________________________________________
AliITSClusterFinderSPDdubna::AliITSClusterFinderSPDdubna(
      const AliITSClusterFinderSPDdubna &s): AliITSClusterFinder(s){
    // The Copy constructortor the class AliITSClusterFinderSPDdugna
    // It calles the = operator.
    // Inputs:
    //    AliITSClusterFinderSPDdubna @s  The AliITSClusterFinderSPDdubna
    //                                    class to be copy constructed.
    // Outputs:
    //    none.
    // Returned:
    //    none.

    *this = s;
    return;
}
//______________________________________________________________________
AliITSClusterFinderSPDdubna::~AliITSClusterFinderSPDdubna(){
    // destructor

    if (fMap) delete fMap;
}
//______________________________________________________________________
void AliITSClusterFinderSPDdubna::Find1DClusters(Int_t mod){
    // Find one dimensional clusters, i.e.
    // in r*phi(x) direction for each colunm in z direction

    static AliITS *iTS=(AliITS*)gAlice->GetModule("ITS");
  
    // retrieve the parameters 
    Int_t fNofPixels     = fSegmentation->Npz(); 
    Int_t fMaxNofSamples = fSegmentation->Npx();

    // read in digits -> do not apply threshold 
    // signal in fired pixels is always 1
    fMap->FillMap();

    Int_t nofFoundClusters = 0;

    Int_t k,it,m;
    for(k=0;k<fNofPixels;k++) {
	Int_t mmax = 10;  // a size of the window for the cluster finding
	for(it=0;it<fMaxNofSamples;it++) {
	    Int_t lclx   = 0;
	    Int_t xstart = 0;
	    Int_t xstop  = 0;
	    Int_t id     = 0;
	    Int_t ilcl   = 0;
	    for(m=0;m<mmax;m++) {  // find the cluster inside the window
		id = it+m;
		if(id >= fMaxNofSamples) break;    // ! no possible for the fadc
		if(fMap->TestHit(k,id) == kUnused) {   // start of the cluster
		    lclx += 1;
		    if(lclx == 1) xstart = id;
		} // end if fMap
		if(lclx > 0 && fMap->TestHit(k,id) == kEmpty) {  
		    // end of cluster if a gap exists
		    xstop = id-1;
		    ilcl = 1;
		    break;
		} // end if lclx >0 ...
	    }   //  end of m-loop
	    if(lclx == 0 && ilcl == 0) it = id; // no cluster in the window,
	                                        // continue the "it" loop
	    if(id >= fMaxNofSamples && lclx == 0) break; // the x row finished
	    if(id < fMaxNofSamples && ilcl == 0 && lclx > 0) {  
		// cluster end is outside of the window,
		mmax += 5;     // increase mmax and repeat the cluster finding
		it -= 1;
	    } // end if id<
	    if(id >= fMaxNofSamples && lclx > 0) {  // the x row finished but
		xstop = fMaxNofSamples - 1;           // the end cluster exists
		ilcl = 1;
	    }  // end if id>=      
	    // ---  Calculate z and x coordinates for one dimensional clusters
	    if(ilcl == 1) {         // new cluster exists
		it = id;
		mmax = 10;
		nofFoundClusters++;
		Float_t clusterCharge = 0.;
		Float_t zpitch        = fSegmentation->Dpz(k+1);
		Int_t   dummy         = 0;
		Float_t clusterZ, dummyX; 
		fSegmentation->GetPadCxz(dummy,k,dummyX,clusterZ);
		Float_t zstart        = clusterZ - 0.5*zpitch;
		Float_t zstop         = clusterZ + 0.5*zpitch;
		Float_t clusterX      = 0.;
		Int_t   clusterSizeX  = lclx;
		Int_t   clusterSizeZ  = 1;
		Int_t   its;
		for(its=xstart; its<=xstop; its++) {
		    Int_t firedpixel=0;
		    if (fMap->GetHitIndex(k,its)>=0) firedpixel=1; 
		    clusterCharge += firedpixel;
		    clusterX      += its + 0.5;
		} // end for its
		Float_t fRphiPitch = fSegmentation->Dpx(dummy);
		clusterX /= (clusterSizeX/fRphiPitch); // center of gravity for x
		// Write the points (coordinates and some cluster information) 
		// to the AliITSRawClusterSPD object

		AliITSRawClusterSPD clust(clusterZ,clusterX,clusterCharge,
					  clusterSizeZ,clusterSizeX,
					  xstart,xstop,zstart,zstop,k,mod);
		iTS->AddCluster(0,&clust);
	    }    // new cluster (ilcl=1)
	} // X direction loop (it)
    } // Z direction loop (k)
    //fMap->ClearMap();
    return;
}
//______________________________________________________________________
void  AliITSClusterFinderSPDdubna::GroupClusters(){
    // Find two dimensional clusters, i.e. group one dimensional clusters
    // into two dimensional ones (go both in x and z directions).
    // get number of clusters for this module

    Int_t nofClusters = fClusters->GetEntriesFast();
    nofClusters -= fNclusters;
    AliITSRawClusterSPD *clusterI;
    AliITSRawClusterSPD *clusterJ;
    Int_t *label=new Int_t[nofClusters];  
    Int_t i,j;
    for(i=0; i<nofClusters; i++) label[i] = 0;
    for(i=0; i<nofClusters; i++) {
	if(label[i] != 0) continue;
	for(j=i+1; j<nofClusters; j++) { 
	    if(label[j] != 0) continue;
	    clusterI = (AliITSRawClusterSPD*) fClusters->At(i);
	    clusterJ = (AliITSRawClusterSPD*) fClusters->At(j);
	    Bool_t pair = clusterI->Brother(clusterJ,fDz,fDx);
	    if(pair) {
#ifndef DEBUG
		cout << "clusters " << i << "," << j << " before grouping" << endl;
		clusterI->PrintInfo();
		clusterJ->PrintInfo();
#endif
		clusterI->Add(clusterJ);
#ifndef DEBUG
		cout << "remove cluster " << j << endl;
#endif
		label[j] = 1;
		fClusters->RemoveAt(j);
#ifndef DEBUG
		cout << "cluster  " << i << " after grouping" << endl;
		clusterI->PrintInfo();
#endif
	    }  // pair
	} // J clusters  
	label[i] = 1;
    } // I clusters
    fClusters->Compress();
#ifndef DEBUG
    Int_t totalNofClusters = fClusters->GetEntriesFast();
    cout << " Nomber of clusters at the group end ="<< totalNofClusters<<endl;
#endif
    delete [] label;
    return;
}
//______________________________________________________________________
void AliITSClusterFinderSPDdubna::TracksInCluster(){
    // Find tracks creating one cluster
    // get number of clusters for this module
    Int_t nofClusters = fClusters->GetEntriesFast();
    Int_t i, ix, iz, jx, jz, xstart, xstop, zstart, zstop, nclx, nclz;
    const Int_t kTrmax = 100;
    Int_t cltracks[kTrmax], itr, tracki, ii, is, js, ie, ntr, tr0, tr1, tr2;

    nofClusters -= fNclusters;
    for(i=0; i<nofClusters; i++) { 
	ii = 0;
	memset(cltracks,-1,sizeof(int)*kTrmax);
	tr0=tr1=tr2=-1;
	AliITSRawClusterSPD *clusterI = (AliITSRawClusterSPD*) fClusters->At(i);
	nclx = clusterI->NclX();
	nclz = clusterI->NclZ();
	xstart = clusterI->XStartf();
	xstop = clusterI->XStopf();
	zstart = clusterI->Zend()-nclz+1;
	zstop = clusterI->Zend();
	Int_t ind;
	for(iz=0; iz<nclz; iz++) { 
	    jz = zstart + iz;
	    for(ix=0; ix<nclx; ix++) { 
		jx = xstart + ix;
		ind = fMap->GetHitIndex(jz,jx);
		if(ind < 0) continue;
		if(ind == 0 && iz >= 0 && ix > 0) continue;
		if(ind == 0 && iz > 0 && ix >= 0) continue;
		if(ind == 0 && iz == 0 && ix == 0 && i > 0) continue;
		AliITSdigitSPD *dig = (AliITSdigitSPD*)fMap->GetHit(jz,jx);
		for(itr=0; itr<3; itr++) { 
		    tracki = dig->GetTrack(i);
		    if(tracki >= 0) {
			ii += 1;
			if(ii < 99) cltracks[ii-1] = tracki;
		    } // end if tracki >-0
		} // end for itr
	    } // ix pixel
	}  // iz pixel
	for(is=0; is<kTrmax; is++) { 
	    if(cltracks[is]<0) continue;
	    for(js=is+1; js<kTrmax; js++) { 
		if(cltracks[js]<0) continue;
		if(cltracks[js]==cltracks[is]) cltracks[js]=-5;
	    } // end for js
	} // end for is
	ntr = 0;
	for(ie=0; ie<kTrmax; ie++) { 
	    if(cltracks[ie] >= 0) {
		ntr=ntr+1;
		if(ntr==1) tr0=cltracks[ie];
		if(ntr==2) tr1=cltracks[ie];
		if(ntr==3) tr2=cltracks[ie];
	    } // end if cltracks[ie] >-=
	} // end for ie
	// if delta ray only
	if(ntr == 0) ntr = 1;
	clusterI->SetNTracks(ntr);
	clusterI->SetTracks(tr0,tr1,tr2);
    } // I cluster
}
//______________________________________________________________________
void AliITSClusterFinderSPDdubna::GetRecPoints(){
    // get rec points
    static AliITS *iTS = (AliITS*)gAlice->GetModule("ITS");
    // get number of clusters for this module
    Int_t nofClusters   = fClusters->GetEntriesFast();
    const Float_t kconv = 1.0e-4;
    const Float_t kRMSx = 12.0*kconv; // microns -> cm ITS TDR Table 1.3
    const Float_t kRMSz = 70.0*kconv; // microns -> cm ITS TDR Table 1.3
    Float_t spdLength   = fSegmentation->Dz();
    Float_t spdWidth    = fSegmentation->Dx();
    Int_t   dummy       = 0;
    Float_t xpitch      = fSegmentation->Dpx(dummy);
    Int_t i;
    Int_t track0, track1, track2;

    nofClusters -= fNclusters;
    for(i=0; i<nofClusters; i++) {
	AliITSRawClusterSPD *clusterI = (AliITSRawClusterSPD*) fClusters->At(i);
	Int_t clustersizex = clusterI->NclX();
	Int_t clustersizez = clusterI->NclZ();
	Int_t xstart = clusterI->XStartf();
	Float_t clusterx = 0;
	Float_t clusterz = clusterI->Z();
	for(Int_t ii=0; ii<clustersizex; ii++) clusterx += (xstart+0.5+ii)*xpitch;
	clusterx /= clustersizex;
	clusterz /= clustersizez;   
	clusterI->GetTracks(track0, track1, track2); 
	AliITSRecPoint rnew;
	rnew.SetX((clusterx - spdWidth/2)*kconv);
	rnew.SetZ((clusterz - spdLength/2)*kconv);
	rnew.SetQ(1.);
	rnew.SetdEdX(0.);
	rnew.SetSigmaX2(kRMSx*kRMSx);
	rnew.SetSigmaZ2(kRMSz*kRMSz);
	rnew.fTracks[0]=track0;
	rnew.fTracks[1]=track1;
	rnew.fTracks[2]=track2;
	iTS->AddRecPoint(rnew);
    } // I clusters
    fMap->ClearMap();
}
//______________________________________________________________________
void AliITSClusterFinderSPDdubna::FindRawClusters(Int_t mod){
    // find raw clusters

#ifndef DEBUG
    cout<<"SPDdubna: module ="<<mod<<endl;
#endif
    Find1DClusters(mod);
    GroupClusters();
    TracksInCluster();
    GetRecPoints();
}
