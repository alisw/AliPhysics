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
Revision 1.10  2000/10/03 13:51:57  egangler
Removal of useless dependencies via forward declarations

Revision 1.9  2000/10/02 16:58:29  egangler
Cleaning of the code :
-> coding conventions
-> void Streamers
-> some useless includes removed or replaced by "class" statement

Revision 1.8  2000/07/03 11:54:57  morsch
AliMUONSegmentation and AliMUONHitMap have been replaced by AliSegmentation and AliHitMap in STEER
The methods GetPadIxy and GetPadXxy of AliMUONSegmentation have changed name to GetPadI and GetPadC.

Revision 1.7  2000/06/28 15:16:35  morsch
(1) Client code adapted to new method signatures in AliMUONSegmentation (see comments there)
to allow development of slat-muon chamber simulation and reconstruction code in the MUON
framework. The changes should have no side effects (mostly dummy arguments).
(2) Hit disintegration uses 3-dim hit coordinates to allow simulation
of chambers with overlapping modules (MakePadHits, Disintegration).

Revision 1.6  2000/06/28 12:19:18  morsch
More consequent seperation of global input data services (AliMUONClusterInput singleton) and the
cluster and hit reconstruction algorithms in AliMUONClusterFinderVS.
AliMUONClusterFinderVS becomes the base class for clustering and hit reconstruction.
It requires two cathode planes. Small modifications in the code will make it usable for
one cathode plane and, hence, more general (for test beam data).
AliMUONClusterFinder is now obsolete.

Revision 1.5  2000/06/28 08:06:10  morsch
Avoid global variables in AliMUONClusterFinderVS by seperating the input data for the fit from the
algorithmic part of the class. Input data resides inside the AliMUONClusterInput singleton.
It also naturally takes care of the TMinuit instance.

Revision 1.4  2000/06/27 16:18:47  gosset
Finally correct implementation of xm, ym, ixm, iym sizes
when at least three local maxima on cathode 1 or on cathode 2

Revision 1.3  2000/06/22 14:02:45  morsch
Parameterised size of xm[], ym[], ixm[], iym[] correctly implemented (PH)
Some HP scope problems corrected (PH)

Revision 1.2  2000/06/15 07:58:48  morsch
Code from MUON-dev joined

Revision 1.1.2.3  2000/06/09 21:58:33  morsch
Most coding rule violations corrected.

Revision 1.1.2.2  2000/02/15 08:33:52  morsch
Error in calculation of contribution map for double clusters (Split method) corrected   (A.M.)
Error in determination of track list for double cluster (FillCluster method) corrected  (A.M.)
Revised and extended SplitByLocalMaxima method (Isabelle Chevrot):
	- For clusters with more than 2 maxima on one of the cathode planes all valid
	combinations of maxima on the two cathodes are preserved. The position of the maxima is
	taken as the hit position.
	- New FillCluster method with 2 arguments to find tracks associated to the clusters
	defined above added. (Method destinction by argument list not very elegant in this case,
	should be revides (A.M.)
	- Bug in if-statement to handle maximum 1 maximum per plane corrected
	- Two cluster per cathode but only 1 combination valid is handled.
	- More rigerous treatment of 1-2 and 2-1 combinations of maxima.

*/

#include "AliMUONClusterFinderVS.h"
#include "AliMUONDigit.h"
#include "AliMUONRawCluster.h"
#include "AliSegmentation.h"
#include "AliMUONResponse.h"
#include "AliMUONClusterInput.h"
#include "AliMUONHitMapA1.h"
#include "AliRun.h"
#include "AliMUON.h"

#include <TTree.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TPad.h>
#include <TGraph.h> 
#include <TPostScript.h> 
#include <TMinuit.h> 
#include <TF1.h>

#include <stdio.h>
#include <iostream.h>

//_____________________________________________________________________
// This function is minimized in the double-Mathieson fit
void fcnS2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
void fcnS1(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
void fcnCombiS1(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
void fcnCombiS2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

ClassImp(AliMUONClusterFinderVS)

    AliMUONClusterFinderVS::AliMUONClusterFinderVS()
{
// Default constructor
    fInput=AliMUONClusterInput::Instance();
    fHitMap[0] = 0;
    fHitMap[1] = 0;
    fTrack[0]=fTrack[1]=-1;
}

AliMUONClusterFinderVS::AliMUONClusterFinderVS(
    const AliMUONClusterFinderVS & clusterFinder)
{
// Dummy copy Constructor
    ;
}

void AliMUONClusterFinderVS::Decluster(AliMUONRawCluster *cluster)
{
// Decluster by local maxima
    SplitByLocalMaxima(cluster);
}

void AliMUONClusterFinderVS::SplitByLocalMaxima(AliMUONRawCluster *c)
{
// Split complex cluster by local maxima 
    Int_t cath, i;

    fInput->SetCluster(c);

    fMul[0]=c->fMultiplicity[0];
    fMul[1]=c->fMultiplicity[1];

//
//  dump digit information into arrays
//

    Float_t qtot;
    
    for (cath=0; cath<2; cath++) {
	qtot=0;
	for (i=0; i<fMul[cath]; i++)
	{
	    // pointer to digit
	    fDig[i][cath]=fInput->Digit(cath, c->fIndexMap[i][cath]);
	    // pad coordinates
	    fIx[i][cath]= fDig[i][cath]->fPadX;
	    fIy[i][cath]= fDig[i][cath]->fPadY;
	    // pad charge
	    fQ[i][cath] = fDig[i][cath]->fSignal;
	    // pad centre coordinates
	    fSeg[cath]->
		GetPadC(fIx[i][cath], fIy[i][cath], fX[i][cath], fY[i][cath], fZ[i][cath]);
	} // loop over cluster digits
    }  // loop over cathodes


    FindLocalMaxima(c);

//
//  Initialise and perform mathieson fits
    Float_t chi2, oldchi2;
//  ++++++++++++++++++*************+++++++++++++++++++++
//  (1) No more than one local maximum per cathode plane 
//  +++++++++++++++++++++++++++++++*************++++++++
    if ((fNLocal[0]==1 && (fNLocal[1]==0 ||  fNLocal[1]==1)) || 
	(fNLocal[0]==0 && fNLocal[1]==1)) {

// Perform combined single Mathieson fit
// Initial values for coordinates (x,y) 

	// One local maximum on cathodes 1 and 2 (X->cathode 2, Y->cathode 1)
	if (fNLocal[0]==1 &&  fNLocal[1]==1) {
	    fXInit[0]=c->fX[1];
	    fYInit[0]=c->fY[0];
	    // One local maximum on cathode 1 (X,Y->cathode 1)
	} else if (fNLocal[0]==1) {
	    fXInit[0]=c->fX[0];
	    fYInit[0]=c->fY[0];
	    // One local maximum on cathode 2  (X,Y->cathode 2)
	} else {
	    fXInit[0]=c->fX[1];
	    fYInit[0]=c->fY[1];
	}
	fprintf(stderr,"\n cas (1) CombiSingleMathiesonFit(c)\n");
	chi2=CombiSingleMathiesonFit(c);
// 	Int_t ndf = fgNbins[0]+fgNbins[1]-2;
// 	Float_t prob = TMath::Prob(Double_t(chi2),ndf);
// 	prob1->Fill(prob);
// 	chi2_1->Fill(chi2);
	oldchi2=chi2;
	fprintf(stderr," chi2 %f ",chi2);

	c->fX[0]=fXFit[0];
	c->fY[0]=fYFit[0];

	c->fX[1]=fXFit[0];
	c->fY[1]=fYFit[0];
	c->fChi2[0]=chi2;
	c->fChi2[1]=chi2;
	c->fX[0]=fSeg[0]->GetAnod(c->fX[0]);
	c->fX[1]=fSeg[1]->GetAnod(c->fX[1]);
	
// If reasonable chi^2 add result to the list of rawclusters
	//	if (chi2 < 50) {
	if (chi2 < 0.3) {
	    AddRawCluster(*c);
// If not try combined double Mathieson Fit
	} else {
	    fprintf(stderr," MAUVAIS CHI2 !!!\n");
	    if (fNLocal[0]==1 &&  fNLocal[1]==1) {
		fXInit[0]=fX[fIndLocal[0][1]][1];
		fYInit[0]=fY[fIndLocal[0][0]][0];
		fXInit[1]=fX[fIndLocal[0][1]][1];
		fYInit[1]=fY[fIndLocal[0][0]][0];
	    } else if (fNLocal[0]==1) {
		fXInit[0]=fX[fIndLocal[0][0]][0];
		fYInit[0]=fY[fIndLocal[0][0]][0];
		fXInit[1]=fX[fIndLocal[0][0]][0];
		fYInit[1]=fY[fIndLocal[0][0]][0];
	    } else {
		fXInit[0]=fX[fIndLocal[0][1]][1];
		fYInit[0]=fY[fIndLocal[0][1]][1];
		fXInit[1]=fX[fIndLocal[0][1]][1];
		fYInit[1]=fY[fIndLocal[0][1]][1];
	    }
	    
//  Initial value for charge ratios
	    fQrInit[0]=0.5;
	    fQrInit[1]=0.5;
	    fprintf(stderr,"\n cas (1) CombiDoubleMathiesonFit(c)\n");
	    chi2=CombiDoubleMathiesonFit(c);
// 	    Int_t ndf = fgNbins[0]+fgNbins[1]-6;
// 	    Float_t prob = TMath::Prob(chi2,ndf);
// 	    prob2->Fill(prob);
//	    chi2_2->Fill(chi2);
	    
// Was this any better ??
	    fprintf(stderr," Old and new chi2 %f %f ", oldchi2, chi2);
	    if (fFitStat!=0 && chi2>0 && (2.*chi2 < oldchi2)) {
		fprintf(stderr," Split\n");
		// Split cluster into two according to fit result
		Split(c);
	    } else {
		fprintf(stderr," Don't Split\n");
		// Don't split
		AddRawCluster(*c);
	    }
	}

//  +++++++++++++++++++++++++++++++++++++++
//  (2) Two local maxima per cathode plane 
//  +++++++++++++++++++++++++++++++++++++++
    } else if (fNLocal[0]==2 &&  fNLocal[1]==2) {
//
//  Let's look for ghosts first 
//
	Float_t xm[4][2], ym[4][2];
	Float_t dpx, dpy, dx, dy;
	Int_t ixm[4][2], iym[4][2];
	Int_t isec, im1, im2, ico;
//
//  Form the 2x2 combinations
//  0-0, 0-1, 1-0, 1-1	
        ico=0;
	for (im1=0; im1<2; im1++) {
	    for (im2=0; im2<2; im2++) {	    
		xm[ico][0]=fX[fIndLocal[im1][0]][0];
		ym[ico][0]=fY[fIndLocal[im1][0]][0];
		xm[ico][1]=fX[fIndLocal[im2][1]][1];
		ym[ico][1]=fY[fIndLocal[im2][1]][1];

		ixm[ico][0]=fIx[fIndLocal[im1][0]][0];
		iym[ico][0]=fIy[fIndLocal[im1][0]][0];
		ixm[ico][1]=fIx[fIndLocal[im2][1]][1];
		iym[ico][1]=fIy[fIndLocal[im2][1]][1];
		ico++;
	    }
	}
// ico = 0 : first local maximum on cathodes 1 and 2
// ico = 1 : fisrt local maximum on cathode 1 and second on cathode 2
// ico = 2 : second local maximum on cathode 1 and first on cathode 1
// ico = 3 : second local maximum on cathodes 1 and 2

// Analyse the combinations and keep those that are possible !
// For each combination check consistency in x and y	
	Int_t iacc;
	Bool_t accepted[4];
	iacc=0;
	
	for (ico=0; ico<4; ico++) {
	    accepted[ico]=kFALSE;
// cathode one: x-coordinate
	    isec=fSeg[0]->Sector(ixm[ico][0], iym[ico][0]);
	    dpx=fSeg[0]->Dpx(isec)/2.;
	    dx=TMath::Abs(xm[ico][0]-xm[ico][1]);
// cathode two: y-coordinate
	    isec=fSeg[1]->Sector(ixm[ico][1], iym[ico][1]);
	    dpy=fSeg[1]->Dpy(isec)/2.;
	    dy=TMath::Abs(ym[ico][0]-ym[ico][1]);
//	    printf("\n %i %f %f %f %f \n", ico, ym[ico][0], ym[ico][1], dy, dpy );
	    if ((dx <= dpx) && (dy <= dpy)) {
		// consistent
		accepted[ico]=kTRUE;
		iacc++;
	    } else {
		// reject
		accepted[ico]=kFALSE;
	    }
	}

	if (iacc==2) {
	    fprintf(stderr,"\n iacc=2: No problem ! \n");
	} else if (iacc==4) {
	    fprintf(stderr,"\n iacc=4: Ok, but ghost problem !!! \n");
	} else if (iacc==0) {
	    fprintf(stderr,"\n iacc=0: I don't know what to do with this !!!!!!!!! \n");
	}

//  Initial value for charge ratios
	fQrInit[0]=Float_t(fQ[fIndLocal[0][0]][0])/
	    Float_t(fQ[fIndLocal[0][0]][0]+fQ[fIndLocal[1][0]][0]);
	fQrInit[1]=Float_t(fQ[fIndLocal[0][1]][1])/
	    Float_t(fQ[fIndLocal[0][1]][1]+fQ[fIndLocal[1][1]][1]);
	
// ******* iacc = 0 *******
// No combinations found between the 2 cathodes
// We keep the center of gravity of the cluster
	if (iacc==0) {
	    AddRawCluster(*c);
	}

// ******* iacc = 1 *******
// Only one combination found between the 2 cathodes
	if (iacc==1) {

// Initial values for the 2 maxima (x,y)

// 1 maximum is initialised with the maximum of the combination found (X->cathode 2, Y->cathode 1)
// 1 maximum is initialised with the other maximum of the first cathode  
	    if (accepted[0]){
		fprintf(stderr,"ico=0\n");
		fXInit[0]=xm[0][1];
		fYInit[0]=ym[0][0];
		fXInit[1]=xm[3][0];
		fYInit[1]=ym[3][0];
	    } else if (accepted[1]){
		fprintf(stderr,"ico=1\n");
		fXInit[0]=xm[1][1];
		fYInit[0]=ym[1][0];
		fXInit[1]=xm[2][0];
		fYInit[1]=ym[2][0];
	    } else if (accepted[2]){
		fprintf(stderr,"ico=2\n");
		fXInit[0]=xm[2][1];
		fYInit[0]=ym[2][0];
		fXInit[1]=xm[1][0];
		fYInit[1]=ym[1][0];
	    } else if (accepted[3]){
		fprintf(stderr,"ico=3\n");
		fXInit[0]=xm[3][1];
		fYInit[0]=ym[3][0];
		fXInit[1]=xm[0][0];
		fYInit[1]=ym[0][0];
	    }
	    fprintf(stderr,"\n cas (2) CombiDoubleMathiesonFit(c)\n");
	    chi2=CombiDoubleMathiesonFit(c);
// 	    Int_t ndf = fgNbins[0]+fgNbins[1]-6;
// 	    Float_t prob = TMath::Prob(chi2,ndf);
// 	    prob2->Fill(prob);
// 	    chi2_2->Fill(chi2);
	    fprintf(stderr," chi2 %f\n",chi2);

// If reasonable chi^2 add result to the list of rawclusters
	    if (chi2<10) {
		Split(c);

	    } else {
// 1 maximum is initialised with the maximum of the combination found (X->cathode 2, Y->cathode 1)
// 1 maximum is initialised with the other maximum of the second cathode  
		if (accepted[0]){
		    fprintf(stderr,"ico=0\n");
		    fXInit[0]=xm[0][1];
		    fYInit[0]=ym[0][0];
		    fXInit[1]=xm[3][1];
		    fYInit[1]=ym[3][1];
		} else if (accepted[1]){
		    fprintf(stderr,"ico=1\n");
		    fXInit[0]=xm[1][1];
		    fYInit[0]=ym[1][0];
		    fXInit[1]=xm[2][1];
		    fYInit[1]=ym[2][1];
		} else if (accepted[2]){
		    fprintf(stderr,"ico=2\n");
		    fXInit[0]=xm[2][1];
		    fYInit[0]=ym[2][0];
		    fXInit[1]=xm[1][1];
		    fYInit[1]=ym[1][1];
		} else if (accepted[3]){
		    fprintf(stderr,"ico=3\n");
		    fXInit[0]=xm[3][1];
		    fYInit[0]=ym[3][0];
		    fXInit[1]=xm[0][1];
		    fYInit[1]=ym[0][1];
		}
		fprintf(stderr,"\n cas (2) CombiDoubleMathiesonFit(c)\n");
		chi2=CombiDoubleMathiesonFit(c);
// 		Int_t ndf = fgNbins[0]+fgNbins[1]-6;
// 		Float_t prob = TMath::Prob(chi2,ndf);
// 		prob2->Fill(prob);
// 		chi2_2->Fill(chi2);
		fprintf(stderr," chi2 %f\n",chi2);

// If reasonable chi^2 add result to the list of rawclusters
		if (chi2<10) {
		    Split(c);
		} else {
//We keep only the combination found (X->cathode 2, Y->cathode 1)
		    for (Int_t ico=0; ico<2; ico++) {
			if (accepted[ico]) {
			    AliMUONRawCluster cnew;
			    Int_t cath;    
			    for (cath=0; cath<2; cath++) {
				cnew.fX[cath]=Float_t(xm[ico][1]);
				cnew.fY[cath]=Float_t(ym[ico][0]);
				cnew.fMultiplicity[cath]=c->fMultiplicity[cath];
			      	for (i=0; i<fMul[cath]; i++) {
				    cnew.fIndexMap[i][cath]=c->fIndexMap[i][cath];
				    fSeg[cath]->SetPad(fIx[i][cath], fIy[i][cath]);
				}
				fprintf(stderr,"\nRawCluster %d cath %d\n",ico,cath);
				fprintf(stderr,"mult_av %d\n",c->fMultiplicity[cath]);
				FillCluster(&cnew,cath);
			    } 
			    cnew.fClusterType=cnew.PhysicsContribution();
			    AddRawCluster(cnew);
			    fNPeaks++;
			}
		    }
		}
	    }
	}
	
// ******* iacc = 2 *******
// Two combinations found between the 2 cathodes
	if (iacc==2) {

// Was the same maximum taken twice
	    if ((accepted[0]&&accepted[1]) || (accepted[2]&&accepted[3])) {
		fprintf(stderr,"\n Maximum taken twice !!!\n");

// Have a try !! with that 
		if (accepted[0]&&accepted[3]) {
		    fXInit[0]=xm[0][1];
		    fYInit[0]=ym[0][0];
		    fXInit[1]=xm[1][1];
		    fYInit[1]=ym[1][0];
		} else {
		    fXInit[0]=xm[2][1];
		    fYInit[0]=ym[2][0];
		    fXInit[1]=xm[3][1];
		    fYInit[1]=ym[3][0];
		}
		fprintf(stderr,"\n cas (2) CombiDoubleMathiesonFit(c)\n");
		chi2=CombiDoubleMathiesonFit(c);
// 		    Int_t ndf = fgNbins[0]+fgNbins[1]-6;
// 		    Float_t prob = TMath::Prob(chi2,ndf);
// 		    prob2->Fill(prob);
// 		    chi2_2->Fill(chi2);
		Split(c);
		
	    } else {
// No ghosts ! No Problems ! -  Perform one fit only !
		if (accepted[0]&&accepted[3]) {
		    fXInit[0]=xm[0][1];
		    fYInit[0]=ym[0][0];
		    fXInit[1]=xm[3][1];
		    fYInit[1]=ym[3][0];
		} else {
		    fXInit[0]=xm[1][1];
		    fYInit[0]=ym[1][0];
		    fXInit[1]=xm[2][1];
		    fYInit[1]=ym[2][0];
		}
		fprintf(stderr,"\n cas (2) CombiDoubleMathiesonFit(c)\n");
		chi2=CombiDoubleMathiesonFit(c);
// 		    Int_t ndf = fgNbins[0]+fgNbins[1]-6;
// 		    Float_t prob = TMath::Prob(chi2,ndf);
// 		    prob2->Fill(prob);
// 		    chi2_2->Fill(chi2);
		fprintf(stderr," chi2 %f\n",chi2);
		Split(c);
	    }
	    
// ******* iacc = 4 *******
// Four combinations found between the 2 cathodes
// Ghost !!
	} else if (iacc==4) {
// Perform fits for the two possibilities !!	
	    fXInit[0]=xm[0][1];
	    fYInit[0]=ym[0][0];
	    fXInit[1]=xm[3][1];
	    fYInit[1]=ym[3][0];
	    fprintf(stderr,"\n cas (2) CombiDoubleMathiesonFit(c)\n");
	    chi2=CombiDoubleMathiesonFit(c);
// 		Int_t ndf = fgNbins[0]+fgNbins[1]-6;
// 		Float_t prob = TMath::Prob(chi2,ndf);
// 		prob2->Fill(prob);
// 		chi2_2->Fill(chi2);
	    fprintf(stderr," chi2 %f\n",chi2);
	    Split(c);
	    fXInit[0]=xm[1][1];
	    fYInit[0]=ym[1][0];
	    fXInit[1]=xm[2][1];
	    fYInit[1]=ym[2][0];
	    fprintf(stderr,"\n cas (2) CombiDoubleMathiesonFit(c)\n");
	    chi2=CombiDoubleMathiesonFit(c);
// 		ndf = fgNbins[0]+fgNbins[1]-6;
// 		prob = TMath::Prob(chi2,ndf);
// 		prob2->Fill(prob);
// 		chi2_2->Fill(chi2);
	    fprintf(stderr," chi2 %f\n",chi2);
	    Split(c);
	}

    } else if (fNLocal[0]==2 &&  fNLocal[1]==1) {
//  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//  (3) Two local maxima on cathode 1 and one maximum on cathode 2 
//  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
	Float_t xm[4][2], ym[4][2];
	Float_t dpx, dpy, dx, dy;
	Int_t ixm[4][2], iym[4][2];
	Int_t isec, im1, ico;
//
//  Form the 2x2 combinations
//  0-0, 0-1, 1-0, 1-1	
        ico=0;
	for (im1=0; im1<2; im1++) {
	    xm[ico][0]=fX[fIndLocal[im1][0]][0];
	    ym[ico][0]=fY[fIndLocal[im1][0]][0];
	    xm[ico][1]=fX[fIndLocal[0][1]][1];
	    ym[ico][1]=fY[fIndLocal[0][1]][1];
	    
	    ixm[ico][0]=fIx[fIndLocal[im1][0]][0];
	    iym[ico][0]=fIy[fIndLocal[im1][0]][0];
	    ixm[ico][1]=fIx[fIndLocal[0][1]][1];
	    iym[ico][1]=fIy[fIndLocal[0][1]][1];
	    ico++;
	}
// ico = 0 : first local maximum on cathodes 1 and 2
// ico = 1 : second local maximum on cathode 1 and first on cathode 2

// Analyse the combinations and keep those that are possible !
// For each combination check consistency in x and y	
	Int_t iacc;
	Bool_t accepted[4];
	iacc=0;
	
	for (ico=0; ico<2; ico++) {
	    accepted[ico]=kFALSE;
	    isec=fSeg[0]->Sector(ixm[ico][0], iym[ico][0]);
	    dpx=fSeg[0]->Dpx(isec)/2.;
	    dx=TMath::Abs(xm[ico][0]-xm[ico][1]);
	    isec=fSeg[1]->Sector(ixm[ico][1], iym[ico][1]);
	    dpy=fSeg[1]->Dpy(isec)/2.;
	    dy=TMath::Abs(ym[ico][0]-ym[ico][1]);
//	    printf("\n %i %f %f %f %f \n", ico, ym[ico][0], ym[ico][1], dy, dpy );
	    if ((dx <= dpx) && (dy <= dpy)) {
		// consistent
		accepted[ico]=kTRUE;
		iacc++;
	    } else {
		// reject
		accepted[ico]=kFALSE;
	    }
	}
	
	Float_t chi21 = 100;
	Float_t chi22 = 100;
	
	if (accepted[0]) {
	    fXInit[0]=xm[0][1];
	    fYInit[0]=ym[0][0];
	    fXInit[1]=xm[1][0];
	    fYInit[1]=ym[1][0];
	    chi21=CombiDoubleMathiesonFit(c);
// 	    Int_t ndf = fgNbins[0]+fgNbins[1]-6;
// 	    Float_t prob = TMath::Prob(chi2,ndf);
// 	    prob2->Fill(prob);
// 	    chi2_2->Fill(chi21);
	    fprintf(stderr," chi2 %f\n",chi21);
	    if (chi21<10) Split(c);
	} else if (accepted[1]) {
	    fXInit[0]=xm[1][1];
	    fYInit[0]=ym[1][0];
	    fXInit[1]=xm[0][0];
	    fYInit[1]=ym[0][0];
	    chi22=CombiDoubleMathiesonFit(c);
// 	    Int_t ndf = fgNbins[0]+fgNbins[1]-6;
// 	    Float_t prob = TMath::Prob(chi2,ndf);
// 	    prob2->Fill(prob);
// 	    chi2_2->Fill(chi22);
	    fprintf(stderr," chi2 %f\n",chi22);
	    if (chi22<10) Split(c);
	}

	if (chi21 > 10 && chi22 > 10) {
// We keep only the combination found (X->cathode 2, Y->cathode 1)
	    for (Int_t ico=0; ico<2; ico++) {
		if (accepted[ico]) {
		    AliMUONRawCluster cnew;
		    Int_t cath;    
		    for (cath=0; cath<2; cath++) {
			cnew.fX[cath]=Float_t(xm[ico][1]);
			cnew.fY[cath]=Float_t(ym[ico][0]);
			cnew.fMultiplicity[cath]=c->fMultiplicity[cath];
			for (i=0; i<fMul[cath]; i++) {
			    cnew.fIndexMap[i][cath]=c->fIndexMap[i][cath];
			    fSeg[cath]->SetPad(fIx[i][cath], fIy[i][cath]);
			}
			fprintf(stderr,"\nRawCluster %d cath %d\n",ico,cath);
			fprintf(stderr,"mult_av %d\n",c->fMultiplicity[cath]);
			FillCluster(&cnew,cath);
		    } 
		    cnew.fClusterType=cnew.PhysicsContribution();
		    AddRawCluster(cnew);
		    fNPeaks++;
		}
	    }
	}
	
//  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//  (3') One local maximum on cathode 1 and two maxima on cathode 2 
//  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    } else if (fNLocal[0]==1 && fNLocal[1]==2) {
	
	Float_t xm[4][2], ym[4][2];
	Float_t dpx, dpy, dx, dy;
	Int_t ixm[4][2], iym[4][2];
	Int_t isec, im1, ico;
//
//  Form the 2x2 combinations
//  0-0, 0-1, 1-0, 1-1	
        ico=0;
	for (im1=0; im1<2; im1++) {
	    xm[ico][0]=fX[fIndLocal[0][0]][0];
	    ym[ico][0]=fY[fIndLocal[0][0]][0];
	    xm[ico][1]=fX[fIndLocal[im1][1]][1];
	    ym[ico][1]=fY[fIndLocal[im1][1]][1];
	    
	    ixm[ico][0]=fIx[fIndLocal[0][0]][0];
	    iym[ico][0]=fIy[fIndLocal[0][0]][0];
	    ixm[ico][1]=fIx[fIndLocal[im1][1]][1];
	    iym[ico][1]=fIy[fIndLocal[im1][1]][1];
	    ico++;
	}
// ico = 0 : first local maximum on cathodes 1 and 2
// ico = 1 : first local maximum on cathode 1 and second on cathode 2

// Analyse the combinations and keep those that are possible !
// For each combination check consistency in x and y	
	Int_t iacc;
	Bool_t accepted[4];
	iacc=0;
	
	for (ico=0; ico<2; ico++) {
	    accepted[ico]=kFALSE;
	    isec=fSeg[0]->Sector(ixm[ico][0], iym[ico][0]);
	    dpx=fSeg[0]->Dpx(isec)/2.;
	    dx=TMath::Abs(xm[ico][0]-xm[ico][1]);
	    isec=fSeg[1]->Sector(ixm[ico][1], iym[ico][1]);
	    dpy=fSeg[1]->Dpy(isec)/2.;
	    dy=TMath::Abs(ym[ico][0]-ym[ico][1]);
//	    printf("\n %i %f %f %f %f \n", ico, ym[ico][0], ym[ico][1], dy, dpy );
	    if ((dx <= dpx) && (dy <= dpy)) {
		// consistent
		accepted[ico]=kTRUE;
		fprintf(stderr,"ico %d\n",ico);
		iacc++;
	    } else {
		// reject
		accepted[ico]=kFALSE;
	    }
	}

	Float_t chi21 = 100;
	Float_t chi22 = 100;

	if (accepted[0]) {
	    fXInit[0]=xm[0][0];
	    fYInit[0]=ym[0][1];
	    fXInit[1]=xm[1][1];
	    fYInit[1]=ym[1][1];
	    chi21=CombiDoubleMathiesonFit(c);
// 	    Int_t ndf = fgNbins[0]+fgNbins[1]-6;
// 	    Float_t prob = TMath::Prob(chi2,ndf);
// 	    prob2->Fill(prob);
// 	    chi2_2->Fill(chi21);
	    fprintf(stderr," chi2 %f\n",chi21);
	    if (chi21<10) Split(c);
	} else if (accepted[1]) {
	    fXInit[0]=xm[1][0];
	    fYInit[0]=ym[1][1];
	    fXInit[1]=xm[0][1];
	    fYInit[1]=ym[0][1];
	    chi22=CombiDoubleMathiesonFit(c);
// 	    Int_t ndf = fgNbins[0]+fgNbins[1]-6;
// 	    Float_t prob = TMath::Prob(chi2,ndf);
// 	    prob2->Fill(prob);
// 	    chi2_2->Fill(chi22);
	    fprintf(stderr," chi2 %f\n",chi22);
	    if (chi22<10) Split(c);
	}

	if (chi21 > 10 && chi22 > 10) {
//We keep only the combination found (X->cathode 2, Y->cathode 1)
	    for (Int_t ico=0; ico<2; ico++) {
		if (accepted[ico]) {
		    AliMUONRawCluster cnew;
		    Int_t cath;    
		    for (cath=0; cath<2; cath++) {
			cnew.fX[cath]=Float_t(xm[ico][1]);
			cnew.fY[cath]=Float_t(ym[ico][0]);
			cnew.fMultiplicity[cath]=c->fMultiplicity[cath];
			for (i=0; i<fMul[cath]; i++) {
			    cnew.fIndexMap[i][cath]=c->fIndexMap[i][cath];
			    fSeg[cath]->SetPad(fIx[i][cath], fIy[i][cath]);
			}
			fprintf(stderr,"\nRawCluster %d cath %d\n",ico,cath);
			fprintf(stderr,"mult_av %d\n",c->fMultiplicity[cath]);
			FillCluster(&cnew,cath);
		    } 
		    cnew.fClusterType=cnew.PhysicsContribution();
		    AddRawCluster(cnew);
		    fNPeaks++;
		}
	    }
	}

//  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//  (4) At least three local maxima on cathode 1 or on cathode 2 
//  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    } else if (fNLocal[0]>2 || fNLocal[1]>2) {
	
	Int_t param = fNLocal[0]*fNLocal[1];
	Int_t ii;

	Float_t ** xm = new Float_t * [param];
	for (ii=0; ii<param; ii++) xm[ii]=new Float_t [2];
	Float_t ** ym = new Float_t * [param];
	for (ii=0; ii<param; ii++) ym[ii]=new Float_t [2];
	Int_t ** ixm = new Int_t * [param];
	for (ii=0; ii<param; ii++) ixm[ii]=new Int_t [2];
	Int_t ** iym = new Int_t * [param];
	for (ii=0; ii<param; ii++) iym[ii]=new Int_t [2];
	
	Int_t isec, ico;
	Float_t dpx, dpy, dx, dy;

        ico=0;
	for (Int_t im1=0; im1<fNLocal[0]; im1++) {
	    for (Int_t im2=0; im2<fNLocal[1]; im2++) {
		xm[ico][0]=fX[fIndLocal[im1][0]][0];
		ym[ico][0]=fY[fIndLocal[im1][0]][0];
		xm[ico][1]=fX[fIndLocal[im2][1]][1];
		ym[ico][1]=fY[fIndLocal[im2][1]][1];

		ixm[ico][0]=fIx[fIndLocal[im1][0]][0];
		iym[ico][0]=fIy[fIndLocal[im1][0]][0];
		ixm[ico][1]=fIx[fIndLocal[im2][1]][1];
		iym[ico][1]=fIy[fIndLocal[im2][1]][1];
		ico++;
	    }
	}
	
	Int_t nIco = ico;
	
	fprintf(stderr,"nIco %d\n",nIco);
	for (ico=0; ico<nIco; ico++) {
	    fprintf(stderr,"ico = %d\n",ico);
	    isec=fSeg[0]->Sector(ixm[ico][0], iym[ico][0]);
	    dpx=fSeg[0]->Dpx(isec)/2.;
	    dx=TMath::Abs(xm[ico][0]-xm[ico][1]);
	    isec=fSeg[1]->Sector(ixm[ico][1], iym[ico][1]);
	    dpy=fSeg[1]->Dpy(isec)/2.;
	    dy=TMath::Abs(ym[ico][0]-ym[ico][1]);

	    fprintf(stderr,"dx %f dpx %f dy %f dpy %f\n",dx,dpx,dy,dpy);
	    fprintf(stderr,"  X %f Y %f\n",xm[ico][1],ym[ico][0]);
	    if ((dx <= dpx) && (dy <= dpy)) {
		fprintf(stderr,"ok\n");
		Int_t cath;    
		AliMUONRawCluster cnew;
		for (cath=0; cath<2; cath++) {
		    cnew.fX[cath]=Float_t(xm[ico][1]);
		    cnew.fY[cath]=Float_t(ym[ico][0]);
		    cnew.fMultiplicity[cath]=c->fMultiplicity[cath];
		    for (i=0; i<fMul[cath]; i++) {
			cnew.fIndexMap[i][cath]=c->fIndexMap[i][cath];
			fSeg[cath]->SetPad(fIx[i][cath], fIy[i][cath]);
		    }
		    FillCluster(&cnew,cath);
		} 
		cnew.fClusterType=cnew.PhysicsContribution();
		AddRawCluster(cnew);
		fNPeaks++;
	    }
	}
	delete [] xm;
	delete [] ym;
	delete [] ixm;
	delete [] iym;
    }
}

void AliMUONClusterFinderVS::FindLocalMaxima(AliMUONRawCluster* c)
{
// Find all local maxima of a cluster
    printf("\n Find Local maxima  !");
    
    AliMUONDigit* digt;
    
    Int_t cath, cath1; // loops over cathodes
    Int_t i;           // loops over digits
    Int_t j;           // loops over cathodes
//
//  Find local maxima
//
//  counters for number of local maxima
    fNLocal[0]=fNLocal[1]=0;
//  flags digits as local maximum
    Bool_t isLocal[100][2];
    for (i=0; i<100;i++) {
	isLocal[i][0]=isLocal[i][1]=kFALSE;
    }
//  number of next neighbours and arrays to store them 
    Int_t nn;
    Int_t x[10], y[10];
// loop over cathodes
    for (cath=0; cath<2; cath++) {
// loop over cluster digits
	for (i=0; i<fMul[cath]; i++) {
// get neighbours for that digit and assume that it is local maximum	    
	    fSeg[cath]->Neighbours(fIx[i][cath], fIy[i][cath], &nn, x, y);
	    isLocal[i][cath]=kTRUE;
	    Int_t isec= fSeg[cath]->Sector(fIx[i][cath], fIy[i][cath]);
	    Float_t a0 = fSeg[cath]->Dpx(isec)*fSeg[cath]->Dpy(isec);
// loop over next neighbours, if at least one neighbour has higher charger assumption
// digit is not local maximum 
	    for (j=0; j<nn; j++) {
		if (fHitMap[cath]->TestHit(x[j], y[j])==kEmpty) continue;
		digt=(AliMUONDigit*) fHitMap[cath]->GetHit(x[j], y[j]);
		isec=fSeg[cath]->Sector(x[j], y[j]);
		Float_t a1 = fSeg[cath]->Dpx(isec)*fSeg[cath]->Dpy(isec);
		if (digt->fSignal/a1 > fQ[i][cath]/a0) {
		    isLocal[i][cath]=kFALSE;
		    break;
//
// handle special case of neighbouring pads with equal signal
		} else if (digt->fSignal == fQ[i][cath]) {
		    if (fNLocal[cath]>0) {
			for (Int_t k=0; k<fNLocal[cath]; k++) {
			    if (x[j]==fIx[fIndLocal[k][cath]][cath] 
				&& y[j]==fIy[fIndLocal[k][cath]][cath])
			    {
				isLocal[i][cath]=kFALSE;
			    } 
			} // loop over local maxima
		    } // are there already local maxima
		} // same charge ? 
	    } // loop over next neighbours
	    if (isLocal[i][cath]) {
		fIndLocal[fNLocal[cath]][cath]=i;
		fNLocal[cath]++;
	    } 
	} // loop over all digits
    } // loop over cathodes
    
    printf("\n Found %d %d %d %d local Maxima\n",
	   fNLocal[0], fNLocal[1], fMul[0], fMul[1]);
    fprintf(stderr,"\n Cathode 1 local Maxima %d Multiplicite %d\n",fNLocal[0], fMul[0]);
    fprintf(stderr," Cathode 2 local Maxima %d Multiplicite %d\n",fNLocal[1], fMul[1]);
    Int_t ix, iy, isec;
    Float_t dpx, dpy;
    
    
    if (fNLocal[1]==2 &&  (fNLocal[0]==1 || fNLocal[0]==0)) {
	Int_t iback=fNLocal[0];
	
//  Two local maxima on cathode 2 and one maximum on cathode 1 
//  Look for local maxima considering up and down neighbours on the 1st cathode only
//
//  Loop over cluster digits
	cath=0;
	cath1=1;
	
	for (i=0; i<fMul[cath]; i++) {
	    isec=fSeg[cath]->Sector(fIx[i][cath],fIy[i][cath]);
	    dpy=fSeg[cath]->Dpy(isec);
	    dpx=fSeg[cath]->Dpx(isec);
	    if (isLocal[i][cath]) continue;
// Pad position should be consistent with position of local maxima on the opposite cathode
	    if ((TMath::Abs(fX[i][cath]-fX[fIndLocal[0][cath1]][cath1]) > dpx/2.) && 
		(TMath::Abs(fX[i][cath]-fX[fIndLocal[1][cath1]][cath1]) > dpx/2.))
 		continue;

// get neighbours for that digit and assume that it is local maximum	    
	    isLocal[i][cath]=kTRUE;
// compare signal to that on the two neighbours on the left and on the right
// iNN counts the number of neighbours with signal, it should be 1 or 2
	    Int_t iNN=0;

 	    for (fSeg[cath]
		     ->FirstPad(fX[i][cath], fY[i][cath], fZPlane, 0., dpy);
		 fSeg[cath]
		     ->MorePads();
		 fSeg[cath]
		     ->NextPad())
	    {
		ix = fSeg[cath]->Ix();
		iy = fSeg[cath]->Iy();
		// skip the current pad
		if (iy == fIy[i][cath]) continue;
		
		if (fHitMap[cath]->TestHit(ix, iy)!=kEmpty) {
		    iNN++;
		    digt=(AliMUONDigit*) fHitMap[cath]->GetHit(ix,iy);
		    if (digt->fSignal > fQ[i][cath]) isLocal[i][cath]=kFALSE;
		}
	    } // Loop over pad neighbours in y
	    if (isLocal[i][cath] && iNN>0) {
		fIndLocal[fNLocal[cath]][cath]=i;
		fNLocal[cath]++;
	    } 
	} // loop over all digits
// if one additional maximum has been found we are happy 
// if more maxima have been found restore the previous situation
	fprintf(stderr,
		"\n New search gives %d local maxima for cathode 1 \n",
		fNLocal[0]);
	fprintf(stderr,
		"                  %d local maxima for cathode 2 \n",
		fNLocal[1]);
	if (fNLocal[cath]>2) {
	    fNLocal[cath]=iback;
	}
	
    } // 1,2 local maxima
    
    if (fNLocal[0]==2 &&  (fNLocal[1]==1 || fNLocal[1]==0)) {
	Int_t iback=fNLocal[1];
	
//  Two local maxima on cathode 1 and one maximum on cathode 2 
//  Look for local maxima considering left and right neighbours on the 2nd cathode only
	cath=1;
	Int_t cath1=0;
//
//  Loop over cluster digits
	for (i=0; i<fMul[cath]; i++) {
	    isec=fSeg[cath]->Sector(fIx[i][cath],fIy[i][cath]);
	    dpx=fSeg[cath]->Dpx(isec);
	    dpy=fSeg[cath]->Dpy(isec);
	    if (isLocal[i][cath]) continue;
// Pad position should be consistent with position of local maxima on the opposite cathode
	    if ((TMath::Abs(fY[i][cath]-fY[fIndLocal[0][cath1]][cath1]) > dpy/2.) && 
		(TMath::Abs(fY[i][cath]-fY[fIndLocal[1][cath1]][cath1]) > dpy/2.))
 		continue;
//
// get neighbours for that digit and assume that it is local maximum	    
	    isLocal[i][cath]=kTRUE;
// compare signal to that on the two neighbours on the left and on the right

// iNN counts the number of neighbours with signal, it should be 1 or 2
	    Int_t iNN=0;
 	    for (fSeg[cath]
		     ->FirstPad(fX[i][cath], fY[i][cath], fZPlane, 0., dpx);
		 fSeg[cath]
		     ->MorePads();
		 fSeg[cath]
		     ->NextPad())
	    {
		ix = fSeg[cath]->Ix();
		iy = fSeg[cath]->Iy();
		
		// skip the current pad
		if (ix == fIx[i][cath]) continue;
		
		if (fHitMap[cath]->TestHit(ix, iy)!=kEmpty) {
		    iNN++;
		    digt=(AliMUONDigit*) fHitMap[cath]->GetHit(ix,iy);
		    if (digt->fSignal > fQ[i][cath]) isLocal[i][cath]=kFALSE;
		}
	    } // Loop over pad neighbours in x
	    if (isLocal[i][cath] && iNN>0) {
		fIndLocal[fNLocal[cath]][cath]=i;
		fNLocal[cath]++;
	    } 
	} // loop over all digits
// if one additional maximum has been found we are happy 
// if more maxima have been found restore the previous situation
	fprintf(stderr,"\n New search gives %d local maxima for cathode 1 \n",fNLocal[0]);
	fprintf(stderr,"\n                  %d local maxima for cathode 2 \n",fNLocal[1]);
//	printf("\n New search gives %d %d \n",fNLocal[0],fNLocal[1]);
	if (fNLocal[cath]>2) {
	    fNLocal[cath]=iback;
	}
    } // 2,1 local maxima
}


void  AliMUONClusterFinderVS::FillCluster(AliMUONRawCluster* c, Int_t flag, Int_t cath) 
{
//
//  Completes cluster information starting from list of digits
//
    AliMUONDigit* dig;
    Float_t x, y, z;
    Int_t  ix, iy;
    
    if (cath==1) {
	c->fPeakSignal[cath]=c->fPeakSignal[0];	
    } else {
	c->fPeakSignal[cath]=0;
    }
    
    
    if (flag) {
	c->fX[cath]=0;
	c->fY[cath]=0;
	c->fQ[cath]=0;
    }

//    fprintf(stderr,"\n fPeakSignal %d\n",c->fPeakSignal[cath]);
    for (Int_t i=0; i<c->fMultiplicity[cath]; i++)
    {
	dig= fInput->Digit(cath,c->fIndexMap[i][cath]);
	ix=dig->fPadX+c->fOffsetMap[i][cath];
	iy=dig->fPadY;
	Int_t q=dig->fSignal;
	if (!flag) q=Int_t(q*c->fContMap[i][cath]);
//	fprintf(stderr,"q %d c->fPeakSignal[ %d ] %d\n",q,cath,c->fPeakSignal[cath]);
	if (dig->fPhysics >= dig->fSignal) {
	    c->fPhysicsMap[i]=2;
	} else if (dig->fPhysics == 0) {
	    c->fPhysicsMap[i]=0;
	} else  c->fPhysicsMap[i]=1;
//
// 
//	fprintf(stderr,"q %d c->fPeakSignal[cath] %d\n",q,c->fPeakSignal[cath]);
// peak signal and track list
	if (q>c->fPeakSignal[cath]) {
	    c->fPeakSignal[cath]=q;
	    c->fTracks[0]=dig->fHit;
	    c->fTracks[1]=dig->fTracks[0];
	    c->fTracks[2]=dig->fTracks[1];
//	    fprintf(stderr," c->fTracks[0] %d c->fTracks[1] %d\n",dig->fHit,dig->fTracks[0]);
	}
//
	if (flag) {
	    fSeg[cath]->GetPadC(ix, iy, x, y, z);
	    c->fX[cath] += q*x;
	    c->fY[cath] += q*y;
	    c->fQ[cath] += q;
	}
    } // loop over digits
//    fprintf(stderr," fin du cluster c\n");


    if (flag) {
    	c->fX[cath]/=c->fQ[cath];
     	c->fX[cath]=fSeg[cath]->GetAnod(c->fX[cath]);
     	c->fY[cath]/=c->fQ[cath]; 
//
//  apply correction to the coordinate along the anode wire
//
     	x=c->fX[cath];   
     	y=c->fY[cath];
     	fSeg[cath]->GetPadI(x, y, fZPlane, ix, iy);
     	fSeg[cath]->GetPadC(ix, iy, x, y, z);
     	Int_t isec=fSeg[cath]->Sector(ix,iy);
    	TF1* cogCorr = fSeg[cath]->CorrFunc(isec-1);
	
     	if (cogCorr) {
	    Float_t yOnPad=(c->fY[cath]-y)/fSeg[cath]->Dpy(isec);
	    c->fY[cath]=c->fY[cath]-cogCorr->Eval(yOnPad, 0, 0);
     	}
    }
}

void  AliMUONClusterFinderVS::FillCluster(AliMUONRawCluster* c, Int_t cath) 
{
//
//  Completes cluster information starting from list of digits
//
    static Float_t dr0;

    AliMUONDigit* dig;

    if (cath==0) {
	dr0 = 10000;
    }
    
    Float_t xpad, ypad, zpad;
    Float_t dx, dy, dr;

    for (Int_t i=0; i<c->fMultiplicity[cath]; i++)
    {
	dig = fInput->Digit(cath,c->fIndexMap[i][cath]);
	fSeg[cath]->
	GetPadC(dig->fPadX,dig->fPadY,xpad,ypad, zpad);
	fprintf(stderr,"x %f y %f cx %f cy %f\n",xpad,ypad,c->fX[0],c->fY[0]);
	dx = xpad - c->fX[0];
	dy = ypad - c->fY[0];
	dr = TMath::Sqrt(dx*dx+dy*dy);

	if (dr < dr0) {
	    dr0 = dr;
	    fprintf(stderr," dr %f\n",dr);
	    Int_t q=dig->fSignal;
	    if (dig->fPhysics >= dig->fSignal) {
		c->fPhysicsMap[i]=2;
	    } else if (dig->fPhysics == 0) {
		c->fPhysicsMap[i]=0;
	    } else  c->fPhysicsMap[i]=1;
	    c->fPeakSignal[cath]=q;
	    c->fTracks[0]=dig->fHit;
	    c->fTracks[1]=dig->fTracks[0];
	    c->fTracks[2]=dig->fTracks[1];
	    fprintf(stderr," c->fTracks[0] %d c->fTracks[1] %d\n",dig->fHit,dig->fTracks[0]);
	}
//
    } // loop over digits

//  apply correction to the coordinate along the anode wire
    c->fX[cath]=fSeg[cath]->GetAnod(c->fX[cath]);
}

void  AliMUONClusterFinderVS::FindCluster(Int_t i, Int_t j, Int_t cath, AliMUONRawCluster &c){


//
//  Find a super cluster on both cathodes
//
//
//  Add i,j as element of the cluster
//
    
    Int_t idx = fHitMap[cath]->GetHitIndex(i,j);
    AliMUONDigit* dig = (AliMUONDigit*) fHitMap[cath]->GetHit(i,j);
    Int_t q=dig->fSignal;
    Int_t theX=dig->fPadX;
    Int_t theY=dig->fPadY; 
   
    if (q > TMath::Abs(c.fPeakSignal[0]) && q > TMath::Abs(c.fPeakSignal[1])) {
	c.fPeakSignal[cath]=q;
	c.fTracks[0]=dig->fHit;
	c.fTracks[1]=dig->fTracks[0];
	c.fTracks[2]=dig->fTracks[1];
    }

//
//  Make sure that list of digits is ordered 
// 
    Int_t mu=c.fMultiplicity[cath];
    c.fIndexMap[mu][cath]=idx;
    
    if (dig->fPhysics >= dig->fSignal) {
        c.fPhysicsMap[mu]=2;
    } else if (dig->fPhysics == 0) {
        c.fPhysicsMap[mu]=0;
    } else  c.fPhysicsMap[mu]=1;

    
    if (mu > 0) {
	for (Int_t ind = mu-1; ind >= 0; ind--) {
	    Int_t ist=(c.fIndexMap)[ind][cath];
	    Int_t ql=fInput->Digit(cath, ist)->fSignal;
	    Int_t ix=fInput->Digit(cath, ist)->fPadX;
	    Int_t iy=fInput->Digit(cath, ist)->fPadY;
	    
	    if (q>ql || (q==ql && theX > ix && theY < iy)) {
		c.fIndexMap[ind][cath]=idx;
		c.fIndexMap[ind+1][cath]=ist;
	    } else {
		
		break;
	    }
	}
    }

    c.fMultiplicity[cath]++;
    if (c.fMultiplicity[cath] >= 50 ) {
	printf("FindCluster - multiplicity >50  %d \n",c.fMultiplicity[0]);
	c.fMultiplicity[cath]=49;
    }

// Prepare center of gravity calculation
    Float_t x, y, z;
    fSeg[cath]->GetPadC(i, j, x, y, z);
    
    c.fX[cath] += q*x;
    c.fY[cath] += q*y;
    c.fQ[cath] += q;
//
// Flag hit as "taken"  
    fHitMap[cath]->FlagHit(i,j);
//
//  Now look recursively for all neighbours and pad hit on opposite cathode
//
//  Loop over neighbours
    Int_t ix,iy;
    ix=iy=0;
    Int_t nn;
    Int_t xList[10], yList[10];
    fSeg[cath]->Neighbours(i,j,&nn,xList,yList);
    for (Int_t in=0; in<nn; in++) {
	ix=xList[in];
	iy=yList[in];
	
	if (fHitMap[cath]->TestHit(ix,iy)==kUnused) {
//	    printf("\n Neighbours %d %d %d", cath, ix, iy);
	    FindCluster(ix, iy, cath, c);
	}
	
   }
    Int_t nOpp=0;
    Int_t iXopp[50], iYopp[50];
    
//  Neighbours on opposite cathode 
//  Take into account that several pads can overlap with the present pad
    Int_t isec=fSeg[cath]->Sector(i,j);    
    Int_t iop;
    Float_t dx, dy;

    if (cath==0) {
	iop = 1;
	dx  = (fSeg[cath]->Dpx(isec))/2.;
	dy  = 0.;
    } else {
	iop = 0;
	dx  = 0.;
	dy  = (fSeg[cath]->Dpy(isec))/2;
    }
// loop over pad neighbours on opposite cathode
    for (fSeg[iop]->FirstPad(x, y, fZPlane, dx, dy);
	 fSeg[iop]->MorePads();
	 fSeg[iop]->NextPad())
    {
	
	ix = fSeg[iop]->Ix(); iy = fSeg[iop]->Iy();
//	    printf("\n ix, iy: %f %f %f %d %d %d", x,y,z,ix, iy, fSector);
	if (fHitMap[iop]->TestHit(ix,iy)==kUnused){
	    iXopp[nOpp]=ix;
	    iYopp[nOpp++]=iy;
//	    printf("\n Opposite %d %d %d", iop, ix, iy);
	}
	
    } // Loop over pad neighbours
//  This had to go outside the loop since recursive calls inside the iterator are not possible
//
    Int_t jopp;
    for (jopp=0; jopp<nOpp; jopp++) {
	if (fHitMap[iop]->TestHit(iXopp[jopp],iYopp[jopp]) == kUnused) 
	    FindCluster(iXopp[jopp], iYopp[jopp], iop, c);
    }
}

//_____________________________________________________________________________

void AliMUONClusterFinderVS::FindRawClusters()
{
  //
  // MUON cluster finder from digits -- finds neighbours on both cathodes and 
  // fills the tree with raw clusters
  //

//  Return if no input datad available
    if (!fInput->NDigits(0) && !fInput->NDigits(1)) return;

    fSeg[0] = fInput->Segmentation(0);
    fSeg[1] = fInput->Segmentation(1);

    fHitMap[0]  = new AliMUONHitMapA1(fSeg[0], fInput->Digits(0));
    fHitMap[1]  = new AliMUONHitMapA1(fSeg[1], fInput->Digits(1));

 
    AliMUONDigit *dig;

    Int_t ndig, cath;
    Int_t nskip=0;
    Int_t ncls=0;
    fHitMap[0]->FillHits();
    fHitMap[1]->FillHits();
//
//  Outer Loop over Cathodes
    for (cath=0; cath<2; cath++) {
	for (ndig=0; ndig<fInput->NDigits(cath); ndig++) {
	    dig = fInput->Digit(cath, ndig);
	    Int_t i=dig->fPadX;
	    Int_t j=dig->fPadY;
	    if (fHitMap[cath]->TestHit(i,j)==kUsed ||fHitMap[0]->TestHit(i,j)==kEmpty) {
		nskip++;
		continue;
	    }
	    fprintf(stderr,"\n CATHODE %d CLUSTER %d\n",cath,ncls);
	    AliMUONRawCluster c;
	    c.fMultiplicity[0]=0;
	    c.fMultiplicity[1]=0;
	    c.fPeakSignal[cath]=dig->fSignal;
	    c.fTracks[0]=dig->fHit;
	    c.fTracks[1]=dig->fTracks[0];
	    c.fTracks[2]=dig->fTracks[1];
	    // tag the beginning of cluster list in a raw cluster
	    c.fNcluster[0]=-1;
	    Float_t xcu, ycu;
	    fSeg[cath]->GetPadC(i,j,xcu, ycu, fZPlane);
	    fSector= fSeg[cath]->Sector(i,j)/100;
//	    printf("\n New Seed %d %d ", i,j);
	    
	    FindCluster(i,j,cath,c);
//          ^^^^^^^^^^^^^^^^^^^^^^^^
	    // center of gravity
	    c.fX[0] /= c.fQ[0];
	    c.fX[0]=fSeg[0]->GetAnod(c.fX[0]);
	    c.fY[0] /= c.fQ[0];
	    c.fX[1] /= c.fQ[1];
	    c.fX[1]=fSeg[0]->GetAnod(c.fX[1]);
	    c.fY[1] /= c.fQ[1];
	    fprintf(stderr,"\n Cathode 1 multiplicite %d X(CG) %f Y(CG) %f\n",
		    c.fMultiplicity[0],c.fX[0],c.fY[0]);
	    fprintf(stderr," Cathode 2 multiplicite %d X(CG) %f Y(CG) %f\n",
		    c.fMultiplicity[1],c.fX[1],c.fY[1]);
//
//      Analyse cluster and decluster if necessary
//	
	ncls++;
	c.fNcluster[1]=fNRawClusters;
	c.fClusterType=c.PhysicsContribution();

	fNPeaks=0;
//
//
	Decluster(&c);
//
//      reset Cluster object
	{ // begin local scope
	    for (int k=0;k<c.fMultiplicity[0];k++) c.fIndexMap[k][0]=0;
	} // end local scope

	{ // begin local scope
	    for (int k=0;k<c.fMultiplicity[1];k++) c.fIndexMap[k][1]=0;
	} // end local scope
	
	c.fMultiplicity[0]=c.fMultiplicity[0]=0;

	
	} // end loop ndig
    } // end loop cathodes
    delete fHitMap[0];
    delete fHitMap[1];
}

Float_t AliMUONClusterFinderVS::SingleMathiesonFit(AliMUONRawCluster *c, Int_t cath)
{
// Performs a single Mathieson fit on one cathode
// 
    AliMUONClusterInput& clusterInput = *(AliMUONClusterInput::Instance());
    
    clusterInput.Fitter()->SetFCN(fcnS1);
    clusterInput.Fitter()->mninit(2,10,7);
    Double_t arglist[20];
    Int_t ierflag=0;
    arglist[0]=1;
// Set starting values 
    static Double_t vstart[2];
    vstart[0]=c->fX[1];
    vstart[1]=c->fY[0];
    
    
// lower and upper limits
    static Double_t lower[2], upper[2];
    Int_t ix,iy;
    fSeg[cath]->GetPadI(c->fX[cath], c->fY[cath], fZPlane, ix, iy);
    Int_t isec=fSeg[cath]->Sector(ix, iy);
    lower[0]=vstart[0]-fSeg[cath]->Dpx(isec)/2;
    lower[1]=vstart[1]-fSeg[cath]->Dpy(isec)/2;
    
    upper[0]=lower[0]+fSeg[cath]->Dpx(isec);
    upper[1]=lower[1]+fSeg[cath]->Dpy(isec);
    
// step sizes
    static Double_t step[2]={0.0005, 0.0005};
    
    clusterInput.Fitter()->mnparm(0,"x1",vstart[0],step[0],lower[0],upper[0],ierflag);
    clusterInput.Fitter()->mnparm(1,"y1",vstart[1],step[1],lower[1],upper[1],ierflag);
// ready for minimisation	
    clusterInput.Fitter()->SetPrintLevel(1);
    clusterInput.Fitter()->mnexcm("SET OUT", arglist, 0, ierflag);
    arglist[0]= -1;
    arglist[1]= 0;
    
    clusterInput.Fitter()->mnexcm("SET NOGR", arglist, 0, ierflag);
    clusterInput.Fitter()->mnexcm("MIGRAD", arglist, 0, ierflag);
    clusterInput.Fitter()->mnexcm("EXIT" , arglist, 0, ierflag);
    Double_t fmin, fedm, errdef;
    Int_t   npari, nparx, istat;
      
    clusterInput.Fitter()->mnstat(fmin, fedm, errdef, npari, nparx, istat);  
    fFitStat=istat;
    
// Print results
// Get fitted parameters
    Double_t xrec, yrec;
    TString chname;
    Double_t epxz, b1, b2;
    Int_t ierflg;
    clusterInput.Fitter()->mnpout(0, chname, xrec, epxz, b1, b2, ierflg);	
    clusterInput.Fitter()->mnpout(1, chname, yrec, epxz, b1, b2, ierflg);	
    fXFit[cath]=xrec;
    fYFit[cath]=yrec;
    return fmin;
}

Float_t AliMUONClusterFinderVS::CombiSingleMathiesonFit(AliMUONRawCluster *c)
{
// Perform combined Mathieson fit on both cathode planes
//
    AliMUONClusterInput& clusterInput = *(AliMUONClusterInput::Instance());
    clusterInput.Fitter()->SetFCN(fcnCombiS1);
    clusterInput.Fitter()->mninit(2,10,7);
    Double_t arglist[20];
    Int_t ierflag=0;
    arglist[0]=1;
    static Double_t vstart[2];
    vstart[0]=fXInit[0];
    vstart[1]=fYInit[0];
    
    
// lower and upper limits
    static Float_t lower[2], upper[2];
    Int_t ix,iy,isec;
    fSeg[0]->GetPadI(fXInit[0], fYInit[0], fZPlane, ix, iy);
    isec=fSeg[0]->Sector(ix, iy);
    Float_t dpy=fSeg[0]->Dpy(isec);
    fSeg[1]->GetPadI(fXInit[0], fYInit[0], fZPlane, ix, iy);
    isec=fSeg[1]->Sector(ix, iy);
    Float_t dpx=fSeg[1]->Dpx(isec);

    Int_t icount;
    Float_t xdum, ydum, zdum;

//  Find save upper and lower limits    
    
    icount = 0;
    
    for (fSeg[1]->FirstPad(fXInit[0], fYInit[0], fZPlane, dpx, 0.); 
	 fSeg[1]->MorePads(); fSeg[1]->NextPad())
    {
	ix=fSeg[1]->Ix(); iy=fSeg[1]->Iy();
	fSeg[1]->GetPadC(ix,iy, upper[0], ydum, zdum);	
	if (icount ==0) lower[0]=upper[0];
	icount++;
    }

    if (lower[0]>upper[0]) {xdum=lower[0]; lower[0]=upper[0]; upper[0]=xdum;}
	
    icount=0;
    
    for (fSeg[0]->FirstPad(fXInit[0], fYInit[0], fZPlane, 0., dpy); 
	 fSeg[0]->MorePads(); fSeg[0]->NextPad())
    {
	ix=fSeg[0]->Ix(); iy=fSeg[0]->Iy();
	fSeg[0]->GetPadC(ix,iy,xdum,upper[1],zdum);	
	if (icount ==0) lower[1]=upper[1];
	icount++;
    }
    
    if (lower[1]>upper[1]) {xdum=lower[1]; lower[1]=upper[1]; upper[1]=xdum;}

// step sizes
    static Double_t step[2]={0.00001, 0.0001};
    
    clusterInput.Fitter()->mnparm(0,"x1",vstart[0],step[0],lower[0],upper[0],ierflag);
    clusterInput.Fitter()->mnparm(1,"y1",vstart[1],step[1],lower[1],upper[1],ierflag);
// ready for minimisation	
    clusterInput.Fitter()->SetPrintLevel(1);
    clusterInput.Fitter()->mnexcm("SET OUT", arglist, 0, ierflag);
    arglist[0]= -1;
    arglist[1]= 0;
    
    clusterInput.Fitter()->mnexcm("SET NOGR", arglist, 0, ierflag);
    clusterInput.Fitter()->mnexcm("MIGRAD", arglist, 0, ierflag);
    clusterInput.Fitter()->mnexcm("EXIT" , arglist, 0, ierflag);
    Double_t fmin, fedm, errdef;
    Int_t   npari, nparx, istat;
      
    clusterInput.Fitter()->mnstat(fmin, fedm, errdef, npari, nparx, istat);  
    fFitStat=istat;
    
// Print results
// Get fitted parameters
    Double_t xrec, yrec;
    TString chname;
    Double_t epxz, b1, b2;
    Int_t ierflg;
    clusterInput.Fitter()->mnpout(0, chname, xrec, epxz, b1, b2, ierflg);	
    clusterInput.Fitter()->mnpout(1, chname, yrec, epxz, b1, b2, ierflg);	
    fXFit[0]=xrec;
    fYFit[0]=yrec;
    return fmin;
}

Bool_t AliMUONClusterFinderVS::DoubleMathiesonFit(AliMUONRawCluster *c, Int_t cath)
{
// Performs a double Mathieson fit on one cathode
// 

//
//  Initialise global variables for fit
    AliMUONClusterInput& clusterInput = *(AliMUONClusterInput::Instance());
    clusterInput.Fitter()->SetFCN(fcnS2);
    clusterInput.Fitter()->mninit(5,10,7);
    Double_t arglist[20];
    Int_t ierflag=0;
    arglist[0]=1;
// Set starting values 
    static Double_t vstart[5];
    vstart[0]=fX[fIndLocal[0][cath]][cath];
    vstart[1]=fY[fIndLocal[0][cath]][cath];	
    vstart[2]=fX[fIndLocal[1][cath]][cath];
    vstart[3]=fY[fIndLocal[1][cath]][cath];	
    vstart[4]=Float_t(fQ[fIndLocal[0][cath]][cath])/
	Float_t(fQ[fIndLocal[0][cath]][cath]+fQ[fIndLocal[1][cath]][cath]);
// lower and upper limits
    static Float_t lower[5], upper[5];
    Int_t isec=fSeg[cath]->Sector(fIx[fIndLocal[0][cath]][cath], fIy[fIndLocal[0][cath]][cath]);
    lower[0]=vstart[0]-fSeg[cath]->Dpx(isec);
    lower[1]=vstart[1]-fSeg[cath]->Dpy(isec);
    
    upper[0]=lower[0]+2.*fSeg[cath]->Dpx(isec);
    upper[1]=lower[1]+2.*fSeg[cath]->Dpy(isec);
    
    isec=fSeg[cath]->Sector(fIx[fIndLocal[1][cath]][cath], fIy[fIndLocal[1][cath]][cath]);
    lower[2]=vstart[2]-fSeg[cath]->Dpx(isec)/2;
    lower[3]=vstart[3]-fSeg[cath]->Dpy(isec)/2;
    
    upper[2]=lower[2]+fSeg[cath]->Dpx(isec);
    upper[3]=lower[3]+fSeg[cath]->Dpy(isec);
    
    lower[4]=0.;
    upper[4]=1.;
// step sizes
    static Double_t step[5]={0.0005, 0.0005, 0.0005, 0.0005, 0.0001};
    
    clusterInput.Fitter()->mnparm(0,"x1",vstart[0],step[0],lower[0],upper[0],ierflag);
    clusterInput.Fitter()->mnparm(1,"y1",vstart[1],step[1],lower[1],upper[1],ierflag);
    clusterInput.Fitter()->mnparm(2,"x2",vstart[2],step[2],lower[2],upper[2],ierflag);
    clusterInput.Fitter()->mnparm(3,"y2",vstart[3],step[3],lower[3],upper[3],ierflag);
    clusterInput.Fitter()->mnparm(4,"a0",vstart[4],step[4],lower[4],upper[4],ierflag);
// ready for minimisation	
    clusterInput.Fitter()->SetPrintLevel(-1);
    clusterInput.Fitter()->mnexcm("SET OUT", arglist, 0, ierflag);
    arglist[0]= -1;
    arglist[1]= 0;
    
    clusterInput.Fitter()->mnexcm("SET NOGR", arglist, 0, ierflag);
    clusterInput.Fitter()->mnexcm("MIGRAD", arglist, 0, ierflag);
    clusterInput.Fitter()->mnexcm("EXIT" , arglist, 0, ierflag);
// Get fitted parameters
    Double_t xrec[2], yrec[2], qfrac;
    TString chname;
    Double_t epxz, b1, b2;
    Int_t ierflg;
    clusterInput.Fitter()->mnpout(0, chname, xrec[0], epxz, b1, b2, ierflg);	
    clusterInput.Fitter()->mnpout(1, chname, yrec[0], epxz, b1, b2, ierflg);	
    clusterInput.Fitter()->mnpout(2, chname, xrec[1], epxz, b1, b2, ierflg);	
    clusterInput.Fitter()->mnpout(3, chname, yrec[1], epxz, b1, b2, ierflg);	
    clusterInput.Fitter()->mnpout(4, chname, qfrac,   epxz, b1, b2, ierflg);	

    Double_t fmin, fedm, errdef;
    Int_t   npari, nparx, istat;
      
    clusterInput.Fitter()->mnstat(fmin, fedm, errdef, npari, nparx, istat);  
    fFitStat=istat;
    return kTRUE;
}

Float_t AliMUONClusterFinderVS::CombiDoubleMathiesonFit(AliMUONRawCluster *c)
{
//
// Perform combined double Mathieson fit on both cathode planes
//
    AliMUONClusterInput& clusterInput = *(AliMUONClusterInput::Instance());
    clusterInput.Fitter()->SetFCN(fcnCombiS2);
    clusterInput.Fitter()->mninit(6,10,7);
    Double_t arglist[20];
    Int_t ierflag=0;
    arglist[0]=1;
// Set starting values 
    static Double_t vstart[6];
    vstart[0]=fXInit[0];
    vstart[1]=fYInit[0];
    vstart[2]=fXInit[1];
    vstart[3]=fYInit[1];
    vstart[4]=fQrInit[0];
    vstart[5]=fQrInit[1];
// lower and upper limits
    static Float_t lower[6], upper[6];
    Int_t ix,iy,isec;
    Float_t dpx, dpy;
    
    fSeg[1]->GetPadI(fXInit[0], fYInit[0], fZPlane, ix, iy);
    isec=fSeg[1]->Sector(ix, iy);
    dpx=fSeg[1]->Dpx(isec);

    fSeg[0]->GetPadI(fXInit[0], fYInit[0], fZPlane, ix, iy);
    isec=fSeg[0]->Sector(ix, iy);
    dpy=fSeg[0]->Dpy(isec);


    Int_t icount;
    Float_t xdum, ydum, zdum;
//    printf("\n Cluster Finder: %f %f %f %f  ", fXInit[0], fXInit[1],fYInit[0], fYInit[1] );
    
//  Find save upper and lower limits    
    icount = 0;
    
    for (fSeg[1]->FirstPad(fXInit[0], fYInit[0], fZPlane, dpx, 0.); 
	 fSeg[1]->MorePads(); fSeg[1]->NextPad())
    {
	ix=fSeg[1]->Ix(); iy=fSeg[1]->Iy();
	fSeg[1]->GetPadC(ix,iy,upper[0],ydum,zdum);	
	if (icount ==0) lower[0]=upper[0];
	icount++;
    }
    if (lower[0]>upper[0]) {xdum=lower[0]; lower[0]=upper[0]; upper[0]=xdum;}    
    icount=0;
    
    for (fSeg[0]->FirstPad(fXInit[0], fYInit[0], fZPlane, 0., dpy); 
	 fSeg[0]->MorePads(); fSeg[0]->NextPad())
    {
	ix=fSeg[0]->Ix(); iy=fSeg[0]->Iy();
	fSeg[0]->GetPadC(ix,iy,xdum,upper[1],zdum);	
	if (icount ==0) lower[1]=upper[1];
	icount++;
    }
    if (lower[1]>upper[1]) {xdum=lower[1]; lower[1]=upper[1]; upper[1]=xdum;}    

    fSeg[1]->GetPadI(fXInit[1], fYInit[1], fZPlane, ix, iy);
    isec=fSeg[1]->Sector(ix, iy);
    dpx=fSeg[1]->Dpx(isec);
    fSeg[0]->GetPadI(fXInit[1], fYInit[1], fZPlane, ix, iy);
    isec=fSeg[0]->Sector(ix, iy);
    dpy=fSeg[0]->Dpy(isec);


//  Find save upper and lower limits    

    icount=0;
    
    for (fSeg[1]->FirstPad(fXInit[1], fYInit[1], fZPlane, dpx, 0); 
	 fSeg[1]->MorePads(); fSeg[1]->NextPad())
    {
	ix=fSeg[1]->Ix(); iy=fSeg[1]->Iy();
	fSeg[1]->GetPadC(ix,iy,upper[2],ydum,zdum);	
	if (icount ==0) lower[2]=upper[2];
	icount++;
    }
    if (lower[2]>upper[2]) {xdum=lower[2]; lower[2]=upper[2]; upper[2]=xdum;}    

    icount=0;
    
    for (fSeg[0]->FirstPad(fXInit[1], fYInit[1], fZPlane, 0, dpy); 
	 fSeg[0]-> MorePads(); fSeg[0]->NextPad())
    {
	ix=fSeg[0]->Ix(); iy=fSeg[0]->Iy();
	fSeg[0]->GetPadC(ix,iy,xdum,upper[3],zdum);	
	if (icount ==0) lower[3]=upper[3];
	icount++;
    }
    if (lower[3]>upper[3]) {xdum=lower[3]; lower[3]=upper[3]; upper[3]=xdum;}    

    lower[4]=0.;
    upper[4]=1.;
    lower[5]=0.;
    upper[5]=1.;

// step sizes
    static Double_t step[6]={0.0005, 0.0005, 0.0005, 0.0005, 0.001, 0.001};
    clusterInput.Fitter()->mnparm(0,"x1",vstart[0],step[0],lower[0],upper[0],ierflag);
    clusterInput.Fitter()->mnparm(1,"y1",vstart[1],step[1],lower[1],upper[1],ierflag);
    clusterInput.Fitter()->mnparm(2,"x2",vstart[2],step[2],lower[2],upper[2],ierflag);
    clusterInput.Fitter()->mnparm(3,"y2",vstart[3],step[3],lower[3],upper[3],ierflag);
    clusterInput.Fitter()->mnparm(4,"a0",vstart[4],step[4],lower[4],upper[4],ierflag);
    clusterInput.Fitter()->mnparm(5,"a1",vstart[5],step[5],lower[5],upper[5],ierflag);
// ready for minimisation	
    clusterInput.Fitter()->SetPrintLevel(-1);
    clusterInput.Fitter()->mnexcm("SET OUT", arglist, 0, ierflag);
    arglist[0]= -1;
    arglist[1]= 0;
    
    clusterInput.Fitter()->mnexcm("SET NOGR", arglist, 0, ierflag);
    clusterInput.Fitter()->mnexcm("MIGRAD", arglist, 0, ierflag);
    clusterInput.Fitter()->mnexcm("EXIT" , arglist, 0, ierflag);
// Get fitted parameters
    TString chname;
    Double_t epxz, b1, b2;
    Int_t ierflg;
    clusterInput.Fitter()->mnpout(0, chname, fXFit[0],  epxz, b1, b2, ierflg);	
    clusterInput.Fitter()->mnpout(1, chname, fYFit[0],  epxz, b1, b2, ierflg);	
    clusterInput.Fitter()->mnpout(2, chname, fXFit[1],  epxz, b1, b2, ierflg);	
    clusterInput.Fitter()->mnpout(3, chname, fYFit[1],  epxz, b1, b2, ierflg);	
    clusterInput.Fitter()->mnpout(4, chname, fQrFit[0], epxz, b1, b2, ierflg);	
    clusterInput.Fitter()->mnpout(5, chname, fQrFit[1], epxz, b1, b2, ierflg);	

    Double_t fmin, fedm, errdef;
    Int_t   npari, nparx, istat;
      
    clusterInput.Fitter()->mnstat(fmin, fedm, errdef, npari, nparx, istat);  
    fFitStat=istat;
    
    fChi2[0]=fmin;
    fChi2[1]=fmin;
    return fmin;
}

void AliMUONClusterFinderVS::Split(AliMUONRawCluster* c)
{
//
// One cluster for each maximum
//
    Int_t i, j, cath;
    AliMUONClusterInput& clusterInput = *(AliMUONClusterInput::Instance());
    for (j=0; j<2; j++) {
	AliMUONRawCluster cnew;
	for (cath=0; cath<2; cath++) {
	    cnew.fChi2[cath]=fChi2[0];
	    
	    if (fNPeaks == 0) {
		cnew.fNcluster[0]=-1;
		cnew.fNcluster[1]=fNRawClusters;
	    } else {
		cnew.fNcluster[0]=fNPeaks;
		cnew.fNcluster[1]=0;
	    }
	    cnew.fMultiplicity[cath]=0;
	    cnew.fX[cath]=Float_t(fXFit[j]);
	    cnew.fY[cath]=Float_t(fYFit[j]);
	    if (j==0) {
		cnew.fQ[cath]=Int_t(clusterInput.TotalCharge(cath)*fQrFit[cath]);
	    } else {
		cnew.fQ[cath]=Int_t(clusterInput.TotalCharge(cath)*(1-fQrFit[cath]));
	    }
	    fSeg[cath]->SetHit(fXFit[j],fYFit[j],fZPlane);
	    for (i=0; i<fMul[cath]; i++) {
		cnew.fIndexMap[cnew.fMultiplicity[cath]][cath]=
		    c->fIndexMap[i][cath];
		fSeg[cath]->SetPad(fIx[i][cath], fIy[i][cath]);
		Float_t q1=fInput->Response()->IntXY(fSeg[cath]);
		cnew.fContMap[i][cath]
		    =(q1*Float_t(cnew.fQ[cath]))/Float_t(fQ[i][cath]);
		cnew.fMultiplicity[cath]++;
	    }
	    FillCluster(&cnew,0,cath);
	} // cathode loop
	
	cnew.fClusterType=cnew.PhysicsContribution();
	if (cnew.fQ[0]>0 && cnew.fQ[1]>0) AddRawCluster(cnew);
	fNPeaks++;
    }
}


//
// Minimisation functions
// Single Mathieson
void fcnS1(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
    AliMUONClusterInput& clusterInput = *(AliMUONClusterInput::Instance());    
    Int_t i;
    Float_t delta;
    Float_t chisq=0;
    Float_t qcont=0;
    Float_t qtot=0;

    for (i=0; i<clusterInput.Nmul(0); i++) {
	Float_t q0=clusterInput.Charge(i,0);
	Float_t q1=clusterInput.DiscrChargeS1(i,par);
	delta=(q0-q1)/q0;
	chisq+=delta*delta;
	qcont+=q1;
	qtot+=q0;
    }
    f=chisq;
}

void fcnCombiS1(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
    AliMUONClusterInput& clusterInput = *(AliMUONClusterInput::Instance());    
    Int_t i, cath;
    Float_t delta;
    Float_t chisq=0;
    Float_t qcont=0;
    Float_t qtot=0;

    for (cath=0; cath<2; cath++) {
	for (i=0; i<clusterInput.Nmul(cath); i++) {
	    Float_t q0=clusterInput.Charge(i,cath);
	    Float_t q1=clusterInput.DiscrChargeCombiS1(i,par,cath);
	    delta=(q0-q1)/q0;
	    chisq+=delta*delta;
	    qcont+=q1;
	    qtot+=q0;
	}
    }
    f=chisq;
}

// Double Mathieson
void fcnS2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
    AliMUONClusterInput& clusterInput = *(AliMUONClusterInput::Instance());    
    Int_t i;
    Float_t delta;
    Float_t chisq=0;
    Float_t qcont=0;
    Float_t qtot=0;
    
    for (i=0; i<clusterInput.Nmul(0); i++) {

	Float_t q0=clusterInput.Charge(i,0);
	Float_t q1=clusterInput.DiscrChargeS2(i,par);
	delta=(q0-q1)/q0;
	chisq+=delta*delta;
	qcont+=q1;
	qtot+=q0;
    }
    f=chisq;
}

// Double Mathieson
void fcnCombiS2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
    AliMUONClusterInput& clusterInput = *(AliMUONClusterInput::Instance());    
    Int_t i, cath;
    Float_t delta;
    Float_t chisq=0;
    Float_t qcont=0;
    Float_t qtot=0;
    for (cath=0; cath<2; cath++) {
	for (i=0; i<clusterInput.Nmul(cath); i++) {
	    Float_t q0=clusterInput.Charge(i,cath);
	    Float_t q1=clusterInput.DiscrChargeCombiS2(i,par,cath);
	    delta=(q0-q1)/q0;
	    chisq+=delta*delta;
	    qcont+=q1;
	    qtot+=q0;
	}
    }
    f=chisq;
}

void AliMUONClusterFinderVS::AddRawCluster(const AliMUONRawCluster c)
{
  //
  // Add a raw cluster copy to the list
  //
    AliMUON *pMUON=(AliMUON*)gAlice->GetModule("MUON");
    pMUON->AddRawCluster(fInput->Chamber(),c); 
    fNRawClusters++;
    fprintf(stderr,"\nfNRawClusters %d\n",fNRawClusters);
}

Bool_t AliMUONClusterFinderVS::TestTrack(Int_t t) {
    if (fTrack[0]==-1 || fTrack[1]==-1) {
	return kTRUE;
    } else if (t==fTrack[0] || t==fTrack[1]) {
	return kTRUE;
    } else {
	return kFALSE;
    }
}

AliMUONClusterFinderVS& AliMUONClusterFinderVS
::operator = (const AliMUONClusterFinderVS& rhs)
{
// Dummy assignment operator
    return *this;
}









