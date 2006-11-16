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

// -------------------------------
// Class AliMUONClusterFinderVS
// -------------------------------
// Class for clustering and reconstruction of space points
// (Not used by default)

#include "AliMUONClusterFinderVS.h"
#include "AliMUONDigit.h"
#include "AliMUONRawCluster.h"
#include "AliMUONGeometrySegmentation.h"
#include "AliMUONMathieson.h"
#include "AliMUONClusterInput.h"
#include "AliMUONDigitMapA1.h"

#include "AliLog.h"

#include <TMinuit.h> 
#include <TF1.h>
#include <TMinuit.h> 
#include <Riostream.h>


//_____________________________________________________________________
// This function is minimized in the double-Mathieson fit
void fcnS2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
void fcnS1(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
void fcnCombiS1(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
void fcnCombiS2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

/// \cond CLASSIMP
ClassImp(AliMUONClusterFinderVS)
/// \endcond

AliMUONClusterFinderVS::AliMUONClusterFinderVS()
  : TObject(),
    fInput(AliMUONClusterInput::Instance()),
    fDeclusterFlag(0),
    fClusterSize(0),
    fNperMax(0),
    fGhostChi2Cut(1e6),
    fNPeaks(0),
    fNRawClusters(0),
    fRawClusters(0x0),
    fZPlane(0.),
    fSector(0),
    fFitStat(0)
{
/// Default constructor
    fDigitMap[0] = 0;
    fDigitMap[1] = 0;
    fTrack[0]=fTrack[1]=-1;
    fSeg2[0]    = 0;
    fSeg2[1]    = 0;

    for(Int_t i=0; i<100; i++) {
      for (Int_t j=0; j<2; j++) {
        fDig[i][j] = 0;
      }
    } 
    fRawClusters = new TClonesArray("AliMUONRawCluster",1000);
}
 //____________________________________________________________________________
AliMUONClusterFinderVS::~AliMUONClusterFinderVS()
{
/// Destructor

  // Reset tracks information
   fNRawClusters = 0;
   if (fRawClusters) {
     fRawClusters->Delete();
     delete fRawClusters;
   }
}

//____________________________________________________________________________
void AliMUONClusterFinderVS::ResetRawClusters()
{
/// Reset tracks information
  fNRawClusters = 0;
  if (fRawClusters) fRawClusters->Clear();
}
//____________________________________________________________________________
void AliMUONClusterFinderVS::Decluster(AliMUONRawCluster *cluster)
{
/// Decluster by local maxima
    SplitByLocalMaxima(cluster);
}
//____________________________________________________________________________
void AliMUONClusterFinderVS::SplitByLocalMaxima(AliMUONRawCluster *c)
{
/// Split complex cluster by local maxima 
    Int_t cath, i;

    fInput->SetCluster(c);

    fMul[0]=c->GetMultiplicity(0);
    fMul[1]=c->GetMultiplicity(1);

//
//  dump digit information into arrays
//

    Float_t qtot;
    
    for (cath=0; cath<2; cath++) {
      qtot=0;

      for (i=0; i<fMul[cath]; i++) {
	// pointer to digit
	fDig[i][cath]=fInput->Digit(cath, c->GetIndex(i, cath));
	// pad coordinates
	fIx[i][cath]= fDig[i][cath]->PadX();
	fIy[i][cath]= fDig[i][cath]->PadY();
	// pad charge
	fQ[i][cath] = fDig[i][cath]->Signal();
	// pad centre coordinates
	  fSeg2[cath]->
	    GetPadC(fInput->DetElemId(), fIx[i][cath], fIy[i][cath], fX[i][cath], fY[i][cath], fZ[i][cath]);
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
	    fXInit[0]=c->GetX(1);
	    fYInit[0]=c->GetY(0);
	    // One local maximum on cathode 1 (X,Y->cathode 1)
	} else if (fNLocal[0]==1) {
	    fXInit[0]=c->GetX(0);
	    fYInit[0]=c->GetY(0);
	    // One local maximum on cathode 2  (X,Y->cathode 2)
	} else {
	    fXInit[0]=c->GetX(1);
	    fYInit[0]=c->GetY(1);
	}
	AliDebug(1,"cas (1) CombiSingleMathiesonFit(c)");
	chi2=CombiSingleMathiesonFit(c);
// 	Int_t ndf = fgNbins[0]+fgNbins[1]-2;
// 	Float_t prob = TMath::Prob(Double_t(chi2),ndf);
// 	prob1->Fill(prob);
// 	chi2_1->Fill(chi2);
	oldchi2=chi2;
	AliDebug(1,Form(" chi2 %f ",chi2));	   

	c->SetX(0, fXFit[0]);
	c->SetY(0, fYFit[0]);

	c->SetX(1,fXFit[0]);
	c->SetY(1,fYFit[0]);
	c->SetChi2(0,chi2);
	c->SetChi2(1,chi2);
        // Force on anod

	c->SetX(0, fSeg2[0]->GetAnod(fInput->DetElemId(), c->GetX(0)));
	c->SetX(1, fSeg2[1]->GetAnod(fInput->DetElemId(), c->GetX(1)));

	//	c->SetDetElemId(fInput->DetElemId());
	// If reasonable chi^2 add result to the list of rawclusters
	if (chi2 < 0.3) {
	    AddRawCluster(*c);
	    // If not try combined double Mathieson Fit
	} else {
	  	AliDebug(1," MAUVAIS CHI2 !!!\n");
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
	    AliDebug(1,"\n cas (1) CombiDoubleMathiesonFit(c)\n");
	    chi2=CombiDoubleMathiesonFit(c);
// 	    Int_t ndf = fgNbins[0]+fgNbins[1]-6;
// 	    Float_t prob = TMath::Prob(chi2,ndf);
// 	    prob2->Fill(prob);
//	    chi2_2->Fill(chi2);
	    
// Was this any better ??
	    AliDebug(1,Form(" Old and new chi2 %f %f ", oldchi2, chi2));
	    if (fFitStat!=0 && chi2>0 && (2.*chi2 < oldchi2)) {
	      AliDebug(1,"Split");
		// Split cluster into two according to fit result
		Split(c);
	    } else {
	      AliDebug(1,"Do not Split");
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
	Int_t   iacc;
	Bool_t  accepted[4];
	Float_t dr[4] = {1.e4, 1.e4, 1.e4, 1.e4};
	iacc=0;

// In case of staggering maxima are displaced by exactly half the pad-size in y. 
// We have to take into account the numerical precision in the consistency check; 	
	Float_t eps = 1.e-5;
//
	for (ico=0; ico<4; ico++) {
	    accepted[ico]=kFALSE;
// cathode one: x-coordinate
	    isec=fSeg2[0]->Sector(fInput->DetElemId(), ixm[ico][0], iym[ico][0]);
	    dpx=fSeg2[0]->Dpx(fInput->DetElemId(), isec)/2.;
	   
	    dx=TMath::Abs(xm[ico][0]-xm[ico][1]);
// cathode two: y-coordinate

	    isec=fSeg2[1]->Sector(fInput->DetElemId(), ixm[ico][1], iym[ico][1]);
	    dpy=fSeg2[1]->Dpy(fInput->DetElemId(), isec)/2.;
	    
	    dy=TMath::Abs(ym[ico][0]-ym[ico][1]);
	    AliDebug(2,Form("\n %i %f %f %f %f %f %f \n", ico, ym[ico][0], ym[ico][1], dy, dpy, dx, dpx ));
	    if ((dx <= dpx) && (dy <= dpy+eps)) {
		// consistent
		accepted[ico]=kTRUE;
		dr[ico] = TMath::Sqrt(dx*dx+dy*dy);
		iacc++;
	    } else {
		// reject
		accepted[ico]=kFALSE;
	    }
	}
	AliDebug(1,Form("\n iacc= %d:\n", iacc));
	if (iacc == 3) {
	    if (accepted[0] && accepted[1]) {
		if (dr[0] >= dr[1]) {
		    accepted[0]=kFALSE;
		} else {
		    accepted[1]=kFALSE;
		}
	    }

	    if (accepted[2] && accepted[3]) {
		if (dr[2] >= dr[3]) {
		    accepted[2]=kFALSE;
		} else {
		    accepted[3]=kFALSE;
		}
	    }
/*	    
// eliminate one candidate
	    Float_t drmax = 0;
	    Int_t icobad = -1;

	    for (ico=0; ico<4; ico++) {
		if (accepted[ico] && dr[ico] > drmax) {
		    icobad = ico;
		    drmax  = dr[ico];
		}
	    }
	    
	    accepted[icobad] = kFALSE;
*/
	    iacc = 2;
	}
	
	
  	AliDebug(1,Form("\n iacc= %d:\n", iacc));
	if (iacc==2) {
		AliDebug(1,"\n iacc=2: No problem ! \n");
	} else if (iacc==4) {
		AliDebug(1,"\n iacc=4: Ok, but ghost problem !!! \n");
	} else if (iacc==0) {
		AliDebug(1,"\n iacc=0: I don't know what to do with this !!!!!!!!! \n");
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
		AliDebug(1,"ico=0");
		fXInit[0]=xm[0][1];
		fYInit[0]=ym[0][0];
		fXInit[1]=xm[3][0];
		fYInit[1]=ym[3][0];
	    } else if (accepted[1]){
		AliDebug(1,"ico=1");
		fXInit[0]=xm[1][1];
		fYInit[0]=ym[1][0];
		fXInit[1]=xm[2][0];
		fYInit[1]=ym[2][0];
	    } else if (accepted[2]){
		AliDebug(1,"ico=2");
		fXInit[0]=xm[2][1];
		fYInit[0]=ym[2][0];
		fXInit[1]=xm[1][0];
		fYInit[1]=ym[1][0];
	    } else if (accepted[3]){
		AliDebug(1,"ico=3");
		fXInit[0]=xm[3][1];
		fYInit[0]=ym[3][0];
		fXInit[1]=xm[0][0];
		fYInit[1]=ym[0][0];
	    }
		AliDebug(1,"cas (2) CombiDoubleMathiesonFit(c)");
	    chi2=CombiDoubleMathiesonFit(c);
// 	    Int_t ndf = fgNbins[0]+fgNbins[1]-6;
// 	    Float_t prob = TMath::Prob(chi2,ndf);
// 	    prob2->Fill(prob);
// 	    chi2_2->Fill(chi2);
	    AliDebug(1,Form(" chi2 %f\n",chi2));

// If reasonable chi^2 add result to the list of rawclusters
	    if (chi2<10) {
		Split(c);

	    } else {
// 1 maximum is initialised with the maximum of the combination found (X->cathode 2, Y->cathode 1)
// 1 maximum is initialised with the other maximum of the second cathode  
		if (accepted[0]){
			AliDebug(1,"ico=0");
		    fXInit[0]=xm[0][1];
		    fYInit[0]=ym[0][0];
		    fXInit[1]=xm[3][1];
		    fYInit[1]=ym[3][1];
		} else if (accepted[1]){
			AliDebug(1,"ico=1");
		    fXInit[0]=xm[1][1];
		    fYInit[0]=ym[1][0];
		    fXInit[1]=xm[2][1];
		    fYInit[1]=ym[2][1];
		} else if (accepted[2]){
			AliDebug(1,"ico=2");
		    fXInit[0]=xm[2][1];
		    fYInit[0]=ym[2][0];
		    fXInit[1]=xm[1][1];
		    fYInit[1]=ym[1][1];
		} else if (accepted[3]){
			AliDebug(1,"ico=3");
		    fXInit[0]=xm[3][1];
		    fYInit[0]=ym[3][0];
		    fXInit[1]=xm[0][1];
		    fYInit[1]=ym[0][1];
		}
		AliDebug(1,"\n cas (2) CombiDoubleMathiesonFit(c)\n");
		chi2=CombiDoubleMathiesonFit(c);
// 		Int_t ndf = fgNbins[0]+fgNbins[1]-6;
// 		Float_t prob = TMath::Prob(chi2,ndf);
// 		prob2->Fill(prob);
// 		chi2_2->Fill(chi2);
		AliDebug(1,Form(" chi2 %f\n",chi2));

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
				cnew.SetX(cath, Float_t(xm[ico][1]));
				cnew.SetY(cath, Float_t(ym[ico][0]));
				cnew.SetZ(cath, fZPlane);
				cnew.SetMultiplicity(cath,c->GetMultiplicity(cath));
			      	for (i=0; i<fMul[cath]; i++) {
				  cnew.SetIndex(i, cath, c->GetIndex(i,cath));
				  fSeg2[cath]->SetPad(fInput->DetElemId(), fIx[i][cath], fIy[i][cath]);
				}
				AliDebug(1,Form("\nRawCluster %d cath %d\n",ico,cath));
				AliDebug(1,Form("mult_av %d\n",c->GetMultiplicity(cath)));
				FillCluster(&cnew,cath);
			    } 
			    cnew.SetClusterType(cnew.PhysicsContribution());
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
		AliDebug(1,"\n Maximum taken twice !!!\n");

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
		AliDebug(1,"\n cas (2) CombiDoubleMathiesonFit(c)\n");
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
		AliDebug(1,"\n cas (2) CombiDoubleMathiesonFit(c)\n");
		chi2=CombiDoubleMathiesonFit(c);
// 		    Int_t ndf = fgNbins[0]+fgNbins[1]-6;
// 		    Float_t prob = TMath::Prob(chi2,ndf);
// 		    prob2->Fill(prob);
// 		    chi2_2->Fill(chi2);
		AliDebug(1,Form(" chi2 %f\n",chi2));
		Split(c);
	    }
	    
// ******* iacc = 4 *******
// Four combinations found between the 2 cathodes
// Ghost !!
	} else if (iacc==4) {
// Perform fits for the two possibilities !!	
// Accept if charges are compatible on both cathodes
// If none are compatible, keep everything
	    fXInit[0]=xm[0][1];
	    fYInit[0]=ym[0][0];
	    fXInit[1]=xm[3][1];
	    fYInit[1]=ym[3][0];
	    AliDebug(1,"\n cas (2) CombiDoubleMathiesonFit(c)\n");
	    chi2=CombiDoubleMathiesonFit(c);
// 		Int_t ndf = fgNbins[0]+fgNbins[1]-6;
// 		Float_t prob = TMath::Prob(chi2,ndf);
// 		prob2->Fill(prob);
// 		chi2_2->Fill(chi2);
	    AliDebug(1,Form(" chi2 %f\n",chi2));
	    // store results of fit and postpone decision
	    Double_t sXFit[2],sYFit[2],sQrFit[2];
	    Float_t sChi2[2];
	    for (Int_t i=0;i<2;i++) {
		sXFit[i]=fXFit[i];
		sYFit[i]=fYFit[i];
		sQrFit[i]=fQrFit[i];
		sChi2[i]=fChi2[i];
	    }
	    fXInit[0]=xm[1][1];
	    fYInit[0]=ym[1][0];
	    fXInit[1]=xm[2][1];
	    fYInit[1]=ym[2][0];
	    AliDebug(1,"\n cas (2) CombiDoubleMathiesonFit(c)\n");
	    chi2=CombiDoubleMathiesonFit(c);
// 		ndf = fgNbins[0]+fgNbins[1]-6;
// 		prob = TMath::Prob(chi2,ndf);
// 		prob2->Fill(prob);
// 		chi2_2->Fill(chi2);
	    AliDebug(1,Form(" chi2 %f\n",chi2));
	    // We have all informations to perform the decision
	    // Compute the chi2 for the 2 possibilities
	    Float_t chi2fi,chi2si,chi2f,chi2s;

	    chi2f = (TMath::Log(fInput->TotalCharge(0)*fQrFit[0]
		  /  (fInput->TotalCharge(1)*fQrFit[1]) )
		  / fInput->ChargeCorrel() );
	    chi2f *=chi2f;
	    chi2fi = (TMath::Log(fInput->TotalCharge(0)*(1-fQrFit[0])
		  /  (fInput->TotalCharge(1)*(1-fQrFit[1])) )
		  / fInput->ChargeCorrel() );
	    chi2f += chi2fi*chi2fi;

	    chi2s = (TMath::Log(fInput->TotalCharge(0)*sQrFit[0]
		  /  (fInput->TotalCharge(1)*sQrFit[1]) )
		  / fInput->ChargeCorrel() );
	    chi2s *=chi2s;
	    chi2si = (TMath::Log(fInput->TotalCharge(0)*(1-sQrFit[0])
		  /  (fInput->TotalCharge(1)*(1-sQrFit[1])) )
		  / fInput->ChargeCorrel() );
	    chi2s += chi2si*chi2si;

	    // usefull to store the charge matching chi2 in the cluster
	    // fChi2[0]=sChi2[1]=chi2f;
	    // fChi2[1]=sChi2[0]=chi2s;

	    if (chi2f<=fGhostChi2Cut && chi2s<=fGhostChi2Cut)
		c->SetGhost(1);
	    if	 (chi2f>fGhostChi2Cut && chi2s>fGhostChi2Cut) {
		// we keep the ghost
		c->SetGhost(2);
		chi2s=-1;
		chi2f=-1;
	    }
	    if (chi2f<=fGhostChi2Cut)
		Split(c);
	    if (chi2s<=fGhostChi2Cut) {
		// retreive saved values
		for (Int_t i=0;i<2;i++) {
		    fXFit[i]=sXFit[i];
		    fYFit[i]=sYFit[i];
		    fQrFit[i]=sQrFit[i];
		    fChi2[i]=sChi2[i];
		}
		Split(c);
  	    }
	    c->SetGhost(0);
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
	// In case of staggering maxima are displaced by exactly half the pad-size in y. 
        // We have to take into account the numerical precision in the consistency check;
	
	Float_t eps = 1.e-5;

	for (ico=0; ico<2; ico++) {
	    isec=fSeg2[0]->Sector(fInput->DetElemId(), ixm[ico][0], iym[ico][0]);
	    dpx=fSeg2[0]->Dpx(fInput->DetElemId(), isec)/2.;
	    
	    dx=TMath::Abs(xm[ico][0]-xm[ico][1]);
	    isec=fSeg2[1]->Sector(fInput->DetElemId(), ixm[ico][1], iym[ico][1]);
	    dpy=fSeg2[1]->Dpy(fInput->DetElemId(), isec)/2.;
	   
	    dy=TMath::Abs(ym[ico][0]-ym[ico][1]);
	    AliDebug(2,Form("\n %i %f %f %f %f \n", ico, ym[ico][0], ym[ico][1], dy, dpy ));
	    if ((dx <= dpx) && (dy <= dpy+eps)) {
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
	Float_t chi23 = 100;

	//  Initial value for charge ratios
	fQrInit[0]=Float_t(fQ[fIndLocal[0][0]][0])/
	    Float_t(fQ[fIndLocal[0][0]][0]+fQ[fIndLocal[1][0]][0]);
	fQrInit[1]=fQrInit[0];
	
	if (accepted[0] && accepted[1]) {
	    
	    fXInit[0]=0.5*(xm[0][1]+xm[0][0]);
	    fYInit[0]=ym[0][0];
	    fXInit[1]=0.5*(xm[0][1]+xm[1][0]);
	    fYInit[1]=ym[1][0];
	    fQrInit[0]=0.5;
	    fQrInit[1]=0.5;
	    chi23=CombiDoubleMathiesonFit(c);
	    if (chi23<10) {
		Split(c);
		Float_t yst;
		yst = fYFit[0];
		fYFit[0] = fYFit[1];
		fYFit[1] = yst;
		Split(c);
	    }
	} else if (accepted[0]) {
	    fXInit[0]=xm[0][1];
	    fYInit[0]=ym[0][0];
	    fXInit[1]=xm[1][0];
	    fYInit[1]=ym[1][0];
	    chi21=CombiDoubleMathiesonFit(c);
// 	    Int_t ndf = fgNbins[0]+fgNbins[1]-6;
// 	    Float_t prob = TMath::Prob(chi2,ndf);
// 	    prob2->Fill(prob);
// 	    chi2_2->Fill(chi21);
	    AliDebug(1,Form(" chi2 %f\n",chi21));
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
	    AliDebug(1,Form(" chi2 %f\n",chi22));
	    if (chi22<10) Split(c);
	}

	if (chi21 > 10 && chi22 > 10 && chi23 > 10) {
// We keep only the combination found (X->cathode 2, Y->cathode 1)
	    for (Int_t ico=0; ico<2; ico++) {
		if (accepted[ico]) {
		    AliMUONRawCluster cnew;
		    Int_t cath;    
		    for (cath=0; cath<2; cath++) {
			cnew.SetX(cath, Float_t(xm[ico][1]));
			cnew.SetY(cath, Float_t(ym[ico][0]));
			cnew.SetZ(cath, fZPlane);
			cnew.SetMultiplicity(cath, c->GetMultiplicity(cath));
			for (i=0; i<fMul[cath]; i++) {
			    cnew.SetIndex(i, cath, c->GetIndex(i, cath));
			    fSeg2[cath]->SetPad(fInput->DetElemId(), fIx[i][cath], fIy[i][cath]);

			}
			AliDebug(1,Form("\nRawCluster %d cath %d\n",ico,cath));
			AliDebug(1,Form("mult_av %d\n",c->GetMultiplicity(cath)));
			
			FillCluster(&cnew,cath);
		    } 
		    cnew.SetClusterType(cnew.PhysicsContribution());
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
        // In case of staggering maxima are displaced by exactly half the pad-size in y. 
        // We have to take into account the numerical precision in the consistency check; 	
	Float_t eps = 1.e-5;

	
	for (ico=0; ico<2; ico++) {
	    accepted[ico]=kFALSE;
	    isec=fSeg2[0]->Sector(fInput->DetElemId(), ixm[ico][0], iym[ico][0]);
	    dpx=fSeg2[0]->Dpx(fInput->DetElemId(), isec)/2.;
	   
	    dx=TMath::Abs(xm[ico][0]-xm[ico][1]);
	    isec=fSeg2[1]->Sector(fInput->DetElemId(), ixm[ico][1], iym[ico][1]);
	    dpy=fSeg2[1]->Dpy(fInput->DetElemId(), isec)/2.;
	    
	    dy=TMath::Abs(ym[ico][0]-ym[ico][1]);
	    AliDebug(1,Form("\n %i %f %f %f %f \n", ico, ym[ico][0], ym[ico][1], dy, dpy ));
	    if ((dx <= dpx) && (dy <= dpy+eps)) {
		// consistent
		accepted[ico]=kTRUE;
		AliDebug(1,Form("ico %d\n",ico));
		iacc++;
	    } else {
		// reject
		accepted[ico]=kFALSE;
	    }
	}

	Float_t chi21 = 100;
	Float_t chi22 = 100;
	Float_t chi23 = 100;

	fQrInit[1]=Float_t(fQ[fIndLocal[0][1]][1])/
	    Float_t(fQ[fIndLocal[0][1]][1]+fQ[fIndLocal[1][1]][1]);
	
	fQrInit[0]=fQrInit[1];

	
	if (accepted[0] && accepted[1]) {
	    fXInit[0]=xm[0][1];
	    fYInit[0]=0.5*(ym[0][0]+ym[0][1]);
	    fXInit[1]=xm[1][1];
	    fYInit[1]=0.5*(ym[0][0]+ym[1][1]);
	    fQrInit[0]=0.5;
	    fQrInit[1]=0.5;
	    chi23=CombiDoubleMathiesonFit(c);
	    if (chi23<10) {
		Split(c);
		Float_t yst;
		yst = fYFit[0];
		fYFit[0] = fYFit[1];
		fYFit[1] = yst;
		Split(c);
	    }
	} else if (accepted[0]) {
	    fXInit[0]=xm[0][0];
	    fYInit[0]=ym[0][1];
	    fXInit[1]=xm[1][1];
	    fYInit[1]=ym[1][1];
	    chi21=CombiDoubleMathiesonFit(c);
// 	    Int_t ndf = fgNbins[0]+fgNbins[1]-6;
// 	    Float_t prob = TMath::Prob(chi2,ndf);
// 	    prob2->Fill(prob);
// 	    chi2_2->Fill(chi21);
	    AliDebug(1,Form(" chi2 %f\n",chi21));
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
	    AliDebug(1,Form(" chi2 %f\n",chi22));
	    if (chi22<10) Split(c);
	}

	if (chi21 > 10 && chi22 > 10 && chi23 > 10) {
//We keep only the combination found (X->cathode 2, Y->cathode 1)
	    for (Int_t ico=0; ico<2; ico++) {
		if (accepted[ico]) {
		    AliMUONRawCluster cnew;
		    Int_t cath;    
		    for (cath=0; cath<2; cath++) {
			cnew.SetX(cath, Float_t(xm[ico][1]));
			cnew.SetY(cath, Float_t(ym[ico][0]));
			cnew.SetZ(cath, fZPlane);
			cnew.SetMultiplicity(cath, c->GetMultiplicity(cath));
			for (i=0; i<fMul[cath]; i++) {
			  cnew.SetIndex(i, cath, c->GetIndex(i, cath));
			  fSeg2[cath]->SetPad(fInput->DetElemId(), fIx[i][cath], fIy[i][cath]);
			}
			AliDebug(1,Form("\nRawCluster %d cath %d\n",ico,cath));
			AliDebug(1,Form("mult_av %d\n",c->GetMultiplicity(cath)));
			FillCluster(&cnew,cath);
		    } 
		    cnew.SetClusterType(cnew.PhysicsContribution());
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
	AliDebug(1,Form("nIco %d\n",nIco));
	for (ico=0; ico<nIco; ico++) {
	    AliDebug(1,Form("ico = %d\n",ico));
	    isec=fSeg2[0]->Sector(fInput->DetElemId(), ixm[ico][0], iym[ico][0]);
	    dpx=fSeg2[0]->Dpx(fInput->DetElemId(), isec)/2.;
	    
	    dx=TMath::Abs(xm[ico][0]-xm[ico][1]);
	    isec=fSeg2[1]->Sector(fInput->DetElemId(), ixm[ico][1], iym[ico][1]);
	    dpy=fSeg2[1]->Dpy(fInput->DetElemId(), isec)/2.;
	    
	    dy=TMath::Abs(ym[ico][0]-ym[ico][1]);
		AliDebug(1,Form("dx %f dpx %f dy %f dpy %f\n",dx,dpx,dy,dpy));
		AliDebug(1,Form("  X %f Y %f\n",xm[ico][1],ym[ico][0]));
	    if ((dx <= dpx) && (dy <= dpy)) {
			AliDebug(1,"ok\n");
		Int_t cath;    
		AliMUONRawCluster cnew;
		for (cath=0; cath<2; cath++) {
		    cnew.SetX(cath, Float_t(xm[ico][1]));
		    cnew.SetY(cath, Float_t(ym[ico][0]));
		    cnew.SetZ(cath, fZPlane);
		    cnew.SetMultiplicity(cath, c->GetMultiplicity(cath));
		    for (i=0; i<fMul[cath]; i++) {
			cnew.SetIndex(i, cath, c->GetIndex(i, cath));
			fSeg2[cath]->SetPad(fInput->DetElemId(), fIx[i][cath], fIy[i][cath]);
		    }
		    FillCluster(&cnew,cath);
		} 
		cnew.SetClusterType(cnew.PhysicsContribution());
		//		cnew.SetDetElemId(fInput->DetElemId());
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

void AliMUONClusterFinderVS::FindLocalMaxima(AliMUONRawCluster* /*c*/)
{
/// Find all local maxima of a cluster
    AliDebug(1,"\n Find Local maxima  !");
    
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
	Int_t isec;
	Float_t a0;

	fSeg2[cath]->Neighbours(fInput->DetElemId(), fIx[i][cath], fIy[i][cath], &nn, x, y);
	  
	isLocal[i][cath]=kTRUE;
	isec = fSeg2[cath]->Sector(fInput->DetElemId(), fIx[i][cath], fIy[i][cath]);
	a0   = fSeg2[cath]->Dpx(fInput->DetElemId(), isec)*fSeg2[cath]->Dpy(fInput->DetElemId(), isec);
	
	// loop over next neighbours, if at least one neighbour has higher charger assumption
	// digit is not local maximum 
	for (j=0; j<nn; j++) {
	  if (fDigitMap[cath]->TestHit(x[j], y[j])==kEmpty) continue;
	  digt=(AliMUONDigit*) fDigitMap[cath]->GetHit(x[j], y[j]);
	  Float_t a1;
	  isec=fSeg2[cath]->Sector(fInput->DetElemId(), x[j], y[j]);
	  a1 = fSeg2[cath]->Dpx(fInput->DetElemId(),isec)*fSeg2[cath]->Dpy(fInput->DetElemId(), isec);
	 
	  if (digt->Signal()/a1 > fQ[i][cath]/a0) {
	    isLocal[i][cath]=kFALSE;
	    break;
	    //
	    // handle special case of neighbouring pads with equal signal
	  } else if (digt->Signal() == fQ[i][cath]) {
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

    AliDebug(1,Form("\n Found %d %d %d %d local Maxima\n",
	       fNLocal[0], fNLocal[1], fMul[0], fMul[1]));
	AliDebug(1,Form("\n Cathode 1 local Maxima %d Multiplicite %d\n",fNLocal[0], fMul[0]));
	AliDebug(1,Form(" Cathode 2 local Maxima %d Multiplicite %d\n",fNLocal[1], fMul[1]));
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
	  isec=fSeg2[cath]->Sector(fInput->DetElemId(), fIx[i][cath],fIy[i][cath]);
	  dpy=fSeg2[cath]->Dpy(fInput->DetElemId(), isec);
	  dpx=fSeg2[cath]->Dpx(fInput->DetElemId(), isec);
	  
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


	    for (fSeg2[cath]->FirstPad(fInput->DetElemId(), fX[i][cath], fY[i][cath], fZPlane, 0., dpy);
		   fSeg2[cath]->MorePads(fInput->DetElemId());
		   fSeg2[cath]->NextPad(fInput->DetElemId()))
		{
		  ix = fSeg2[cath]->Ix();
		  iy = fSeg2[cath]->Iy();
		  // skip the current pad
		  if (iy == fIy[i][cath]) continue;
		
		  if (fDigitMap[cath]->TestHit(ix, iy)!=kEmpty) {
		    iNN++;
		    digt=(AliMUONDigit*) fDigitMap[cath]->GetHit(ix,iy);
		    if (digt->Signal() > fQ[i][cath]) isLocal[i][cath]=kFALSE;
		  }
		} // Loop over pad neighbours in y
	    
	    if (isLocal[i][cath] && iNN>0) {
		fIndLocal[fNLocal[cath]][cath]=i;
		fNLocal[cath]++;
	    } 
	} // loop over all digits
// if one additional maximum has been found we are happy 
// if more maxima have been found restore the previous situation
	AliDebug(1,Form("\n New search gives %d local maxima for cathode 1 \n",
		    fNLocal[0]));
	AliDebug(1,Form("                  %d local maxima for cathode 2 \n",
		    fNLocal[1]));
	if (fNLocal[cath]>2) {
	    fNLocal[cath]=iback;
	}
	
    } // 1,2 local maxima
    
    if (fNLocal[0]==2 &&  (fNLocal[1]==1 || fNLocal[1]==0)) {
	Int_t iback=fNLocal[1];
	
//  Two local maxima on cathode 1 and one maximum on cathode 2 
//  Look for local maxima considering left and right neighbours on the 2nd cathode only
	cath=1;
	Int_t cath1 = 0;
	Float_t eps = 1.e-5;
	
//
//  Loop over cluster digits
	for (i=0; i<fMul[cath]; i++) {
	    isec=fSeg2[cath]->Sector(fInput->DetElemId(), fIx[i][cath],fIy[i][cath]);
	    dpx=fSeg2[cath]->Dpx(fInput->DetElemId(), isec);
	    dpy=fSeg2[cath]->Dpy(fInput->DetElemId(), isec);
	  
	
	    if (isLocal[i][cath]) continue;
// Pad position should be consistent with position of local maxima on the opposite cathode
	    if ((TMath::Abs(fY[i][cath]-fY[fIndLocal[0][cath1]][cath1]) > dpy/2.+eps) && 
		(TMath::Abs(fY[i][cath]-fY[fIndLocal[1][cath1]][cath1]) > dpy/2.+eps))
 		continue;
	    
//
// get neighbours for that digit and assume that it is local maximum	    
	    isLocal[i][cath]=kTRUE;
// compare signal to that on the two neighbours on the left and on the right

// iNN counts the number of neighbours with signal, it should be 1 or 2
	    Int_t iNN=0;
	    for (fSeg2[cath]->FirstPad(fInput->DetElemId(), fX[i][cath], fY[i][cath], fZPlane, dpx, 0.);
		   fSeg2[cath]->MorePads(fInput->DetElemId());
		   fSeg2[cath]->NextPad(fInput->DetElemId()))
		{

		  ix = fSeg2[cath]->Ix();
		  iy = fSeg2[cath]->Iy();

		  // skip the current pad
		  if (ix == fIx[i][cath]) continue;
		
		  if (fDigitMap[cath]->TestHit(ix, iy)!=kEmpty) {
		    iNN++;
		    digt=(AliMUONDigit*) fDigitMap[cath]->GetHit(ix,iy);
		    if (digt->Signal() > fQ[i][cath]) isLocal[i][cath]=kFALSE;
		  }
		} // Loop over pad neighbours in x
	    
	    if (isLocal[i][cath] && iNN>0) {
		fIndLocal[fNLocal[cath]][cath]=i;
		fNLocal[cath]++;
	    } 
	} // loop over all digits
// if one additional maximum has been found we are happy 
// if more maxima have been found restore the previous situation
	AliDebug(1,Form("\n New search gives %d local maxima for cathode 1 \n",fNLocal[0]));
	AliDebug(1,Form("\n                  %d local maxima for cathode 2 \n",fNLocal[1]));
	AliDebug(1,Form("\n New search gives %d %d \n",fNLocal[0],fNLocal[1]));
	if (fNLocal[cath]>2) {
	    fNLocal[cath]=iback;
	}
    } // 2,1 local maxima
}


void  AliMUONClusterFinderVS::FillCluster(AliMUONRawCluster* c, Int_t flag, Int_t cath) 
{
///  Completes cluster information starting from list of digits

  AliMUONDigit* dig;
  Float_t x, y, z;
  Int_t  ix, iy;
  
  if (cath==1) {
    c->SetPeakSignal(cath,c->GetPeakSignal(0));	
  } else {
    c->SetPeakSignal(cath,0);
  }
  
  
  if (flag) {
    c->SetX(cath,0.);
    c->SetY(cath,0.);
    c->SetCharge(cath,0);
  }
  
  AliDebug(1,Form("\n fPeakSignal %d\n",c->GetPeakSignal(cath)));
  for (Int_t i=0; i<c->GetMultiplicity(cath); i++)
    {
      dig= fInput->Digit(cath,c->GetIndex(i,cath));
      ix=dig->PadX()+c->GetOffset(i,cath);
      iy=dig->PadY();
      Float_t q=dig->Signal();
      if (!flag) q *= c->GetContrib(i,cath);
      //	fprintf(stderr,"q %d c->fPeakSignal[ %d ] %d\n",q,cath,c->fPeakSignal[cath]);
      if (dig->Physics() >= dig->Signal()) {
	c->SetPhysics(i,2);
      } else if (dig->Physics() == 0) {
	c->SetPhysics(i,0);
      } else  c->SetPhysics(i,1);
      //
      // 
      AliDebug(2,Form("q %d c->fPeakSignal[cath] %d\n",q,c->GetPeakSignal(cath)));
      // peak signal and track list
      if (q>c->GetPeakSignal(cath)) {
	c->SetPeakSignal(cath, q);
	c->SetTrack(0,dig->Hit());
	c->SetTrack(1,dig->Track(0));
	c->SetTrack(2,dig->Track(1));
	//	    fprintf(stderr," c->fTracks[0] %d c->fTracks[1] %d\n",dig->fHit,dig->fTracks[0]);
      }
      //
      if (flag) {
	fSeg2[cath]->GetPadC(fInput->DetElemId(), ix, iy, x, y, z);
	
	c->AddX(cath, q*x);
	c->AddY(cath, q*y);
	c->AddCharge(cath, q);
      }
    } // loop over digits
  AliDebug(1," fin du cluster c\n");
  
  
  if (flag) {
    c->SetX(cath, c->GetX(cath)/c->GetCharge(cath));
    // Force on anod
    c->SetX(cath, fSeg2[cath]->GetAnod(fInput->DetElemId(), c->GetX(cath)));
    c->SetY(cath, c->GetY(cath)/c->GetCharge(cath)); 
    //
    //  apply correction to the coordinate along the anode wire
    //
    x=c->GetX(cath);   
    y=c->GetY(cath);
    TF1* cogCorr;
    Int_t isec;
    fSeg2[cath]->GetPadI(fInput->DetElemId(), x, y, fZPlane, ix, iy);
    fSeg2[cath]->GetPadC(fInput->DetElemId(), ix, iy, x, y, z);
    isec=fSeg2[cath]->Sector(fInput->DetElemId(), ix,iy);
    cogCorr = fSeg2[cath]->CorrFunc(fInput->DetElemId(), isec-1);
    
    
    if (cogCorr) {
      Float_t yOnPad;
      yOnPad=(c->GetY(cath)-y)/fSeg2[cath]->Dpy(fInput->DetElemId(), isec);
      
      c->SetY(cath, c->GetY(cath)-cogCorr->Eval(yOnPad, 0, 0));
      // slat ID from digit
      
    }
  }
}

void  AliMUONClusterFinderVS::FillCluster(AliMUONRawCluster* c, Int_t cath) 
{
///  Completes cluster information starting from list of digits

  static Float_t dr0;
  
  AliMUONDigit* dig;
  
  if (cath==0) {
    dr0 = 10000;
  }
  
  Float_t xpad, ypad, zpad;
  Float_t dx, dy, dr;
  
  for (Int_t i=0; i<c->GetMultiplicity(cath); i++)
    {
      dig = fInput->Digit(cath,c->GetIndex(i,cath));
      fSeg2[cath]->
	  GetPadC(fInput->DetElemId(),dig->PadX(),dig->PadY(),xpad,ypad, zpad);
      AliDebug(1,Form("x %f y %f cx %f cy %f\n",xpad,ypad,c->GetX(0),c->GetY(0)));
      dx = xpad - c->GetX(0);
      dy = ypad - c->GetY(0);
      dr = TMath::Sqrt(dx*dx+dy*dy);
      
      if (dr < dr0) {
	dr0 = dr;
	AliDebug(1,Form(" dr %f\n",dr));
	Float_t q=dig->Signal();
	if (dig->Physics() >= dig->Signal()) {
	  c->SetPhysics(i,2);
	} else if (dig->Physics() == 0) {
	  c->SetPhysics(i,0);
	} else  c->SetPhysics(i,1);
	c->SetPeakSignal(cath,q);
	c->SetTrack(0,dig->Hit());
	c->SetTrack(1,dig->Track(0));
	c->SetTrack(2,dig->Track(1));
	
	AliDebug(1,Form(" c->fTracks[0] %d c->fTracks[1] %d\n",dig->Hit(),
			dig->Track(0)));
      }
      //
    } // loop over digits
  
  //  apply correction to the coordinate along the anode wire
  // Force on anod
    c->SetX(cath,fSeg2[cath]->GetAnod(fInput->DetElemId(), c->GetX(cath)));
}

void  AliMUONClusterFinderVS::FindCluster(Int_t i, Int_t j, Int_t cath, AliMUONRawCluster &c)
{
///  Find a super cluster on both cathodes
///  Add i,j as element of the cluster
    
    Int_t idx = fDigitMap[cath]->GetHitIndex(i,j);
    AliMUONDigit* dig = (AliMUONDigit*) fDigitMap[cath]->GetHit(i,j);
    Float_t q=dig->Signal();
    Int_t theX=dig->PadX();
    Int_t theY=dig->PadY(); 
   
    if (q > TMath::Abs(c.GetPeakSignal(0)) && q > TMath::Abs(c.GetPeakSignal(1))) {
	c.SetPeakSignal(cath,q);
	c.SetTrack(0,dig->Hit());
	c.SetTrack(1,dig->Track(0));
	c.SetTrack(2,dig->Track(1));
    }

//
//  Make sure that list of digits is ordered 
// 
    Int_t mu=c.GetMultiplicity(cath);
    c.SetIndex(mu, cath, idx);
    
    if (dig->Physics() >= dig->Signal()) {
        c.SetPhysics(mu,2);
    } else if (dig->Physics() == 0) {
        c.SetPhysics(mu,0);
    } else  c.SetPhysics(mu,1);

    
    if (mu > 0) {
	for (Int_t ind = mu-1; ind >= 0; ind--) {
	    Int_t ist=c.GetIndex(ind,cath);
	    Float_t ql=fInput->Digit(cath, ist)->Signal();
	    Int_t ix=fInput->Digit(cath, ist)->PadX();
	    Int_t iy=fInput->Digit(cath, ist)->PadY();
	    
	    if (q>ql || (q==ql && theX > ix && theY < iy)) {
		c.SetIndex(ind, cath, idx);
		c.SetIndex(ind+1, cath, ist);
	    } else {
		
		break;
	    }
	}
    }

    c.SetMultiplicity(cath, c.GetMultiplicity(cath)+1);
    if (c.GetMultiplicity(cath) >= 50 ) {
      AliDebug(1,Form("FindCluster - multiplicity >50  %d \n",c.GetMultiplicity(0)));
	c.SetMultiplicity(cath, 49);
    }

// Prepare center of gravity calculation
    Float_t x, y, z;
    fSeg2[cath]->GetPadC(fInput->DetElemId(), i, j, x, y, z);
    c.AddX(cath,q*x);
    c.AddY(cath,q*y);
    c.AddCharge(cath,q);
//
// Flag hit as "taken"  
    fDigitMap[cath]->FlagHit(i,j);
//
//  Now look recursively for all neighbours and pad hit on opposite cathode
//
//  Loop over neighbours
    Int_t ix,iy;
    ix=iy=0;
    Int_t nn;
    Int_t xList[10], yList[10];
    fSeg2[cath]->Neighbours(fInput->DetElemId(), i,j,&nn,xList,yList);
    for (Int_t in=0; in<nn; in++) {
	ix=xList[in];
	iy=yList[in];
	
	if (fDigitMap[cath]->TestHit(ix,iy)==kUnused) {
	    AliDebug(2,Form("\n Neighbours %d %d %d", cath, ix, iy));
	    FindCluster(ix, iy, cath, c);
	}
	
   }
    Int_t nOpp=0;
    Int_t iXopp[50], iYopp[50];
    
//  Neighbours on opposite cathode 
//  Take into account that several pads can overlap with the present pad
    Int_t isec;
    isec=fSeg2[cath]->Sector(fInput->DetElemId(), i,j);    

    Int_t iop;
    Float_t dx, dy;
  
    if (cath==0) {
      iop = 1;
      dx  = (fSeg2[cath]->Dpx(fInput->DetElemId(), isec))/2.;
      dy  = 0.;
    } else {
      iop = 0;
      dx  = 0.;
      dy  = (fSeg2[cath]->Dpy(fInput->DetElemId(), isec))/2;
    }
   

    
    // loop over pad neighbours on opposite cathode
    for (fSeg2[iop]->FirstPad(fInput->DetElemId(), x, y, fZPlane, dx, dy);
	 fSeg2[iop]->MorePads(fInput->DetElemId());
	 fSeg2[iop]->NextPad(fInput->DetElemId()))
      {
	
	ix = fSeg2[iop]->Ix(); iy = fSeg2[iop]->Iy();
	AliDebug(2,Form("\n ix, iy: %f %f %f %d %d %d", x,y,z,ix, iy, fSector));
	if (fDigitMap[iop]->TestHit(ix,iy)==kUnused){
	  iXopp[nOpp]=ix;
	  iYopp[nOpp++]=iy;
	  AliDebug(2,Form("\n Opposite %d %d %d", iop, ix, iy));
	}
	
      } // Loop over pad neighbours
    //  This had to go outside the loop since recursive calls inside the iterator are not possible
    //
    Int_t jopp;
    for (jopp=0; jopp<nOpp; jopp++) {
      if (fDigitMap[iop]->TestHit(iXopp[jopp],iYopp[jopp]) == kUnused) 
	FindCluster(iXopp[jopp], iYopp[jopp], iop, c);
    }

}

//_____________________________________________________________________________

void AliMUONClusterFinderVS::FindRawClusters()
{
/// MUON cluster finder from digits -- finds neighbours on both cathodes and 
/// fills the tree with raw clusters

    ResetRawClusters();
//  Return if no input datad available
    if (!fInput->NDigits(0) && !fInput->NDigits(1)) return;

    fSeg2[0] = fInput->Segmentation2(0);
    fSeg2[1] = fInput->Segmentation2(1);
    
    Int_t detElemId = fInput->DetElemId();
    
    Int_t npx0  = fSeg2[0]->Npx(detElemId)+1;
    Int_t npy0  = fSeg2[0]->Npy(detElemId)+1;
    fDigitMap[0]  = new AliMUONDigitMapA1(detElemId, npx0, npy0);

    Int_t npx1  = fSeg2[0]->Npx(detElemId)+1;
    Int_t npy1  = fSeg2[0]->Npy(detElemId)+1;
    fDigitMap[1]  = new AliMUONDigitMapA1(detElemId, npx1, npy1);
    
    AliMUONDigit *dig;

    Int_t ndig, cath;
    Int_t nskip=0;
    Int_t ncls=0;

    fDigitMap[0]->FillHits(fInput->Digits(0));
    fDigitMap[1]->FillHits(fInput->Digits(1));
//
//  Outer Loop over Cathodes
    for (cath = 0; cath < 2; cath++) {
      
	for (ndig=0; ndig<fInput->NDigits(cath); ndig++) {
	  dig = fInput->Digit(cath, ndig);
	  Int_t padx = dig->PadX();
	  Int_t pady = dig->PadY();
	  if (fDigitMap[cath]->TestHit(padx,pady)==kUsed ||fDigitMap[0]->TestHit(padx,pady)==kEmpty) {
	    nskip++;
	    continue;
	  }
	  AliDebug(1,Form("\n CATHODE %d CLUSTER %d\n",cath,ncls));
	  AliMUONRawCluster clus;
	  clus.SetMultiplicity(0, 0);
	  clus.SetMultiplicity(1, 0);
	  clus.SetPeakSignal(cath,dig->Signal());
	  clus.SetTrack(0, dig->Hit());
	  clus.SetTrack(1, dig->Track(0));
	  clus.SetTrack(2, dig->Track(1));
	  
	  AliDebug(1,Form("idDE %d Padx %d Pady %d", fInput->DetElemId(), padx, pady));
	  
	  // tag the beginning of cluster list in a raw cluster
	  clus.SetNcluster(0,-1);
	  Float_t xcu, ycu;
	  fSeg2[cath]->GetPadC(fInput->DetElemId(), padx, pady, xcu, ycu, fZPlane);
	  fSector= fSeg2[cath]->Sector(fInput->DetElemId(), padx, pady)/100;
	  
	  
	  
	  
	  FindCluster(padx,pady,cath,clus);
	  //^^^^^^^^^^^^^^^^^^^^^^^^
	  // center of gravity
	  if (clus.GetX(0)!=0.) clus.SetX(0, clus.GetX(0)/clus.GetCharge(0)); // clus.fX[0] /= clus.fQ[0];
	  
	  // Force on anod
	  clus.SetX(0,fSeg2[0]->GetAnod(fInput->DetElemId(), clus.GetX(0)));
	  if (clus.GetY(0)!=0.) clus.SetY(0, clus.GetY(0)/clus.GetCharge(0)); // clus.fY[0] /= clus.fQ[0];
	  
	  if(clus.GetCharge(1)!=0.) clus.SetX(1, clus.GetX(1)/clus.GetCharge(1));  // clus.fX[1] /= clus.fQ[1];
	  
	  // Force on anod
	  clus.SetX(1, fSeg2[0]->GetAnod(fInput->DetElemId(),clus.GetX(1)));
	  if(clus.GetCharge(1)!=0.) clus.SetY(1, clus.GetY(1)/clus.GetCharge(1));// clus.fY[1] /= clus.fQ[1];
	  
	  clus.SetZ(0, fZPlane);
	  clus.SetZ(1, fZPlane);	    
	  
	  AliDebug(1,Form("\n Cathode 1 multiplicite %d X(CG) %f Y(CG) %f\n",
			  clus.GetMultiplicity(0),clus.GetX(0),clus.GetY(0)));
	  AliDebug(1,Form(" Cathode 2 multiplicite %d X(CG) %f Y(CG) %f\n",
			  clus.GetMultiplicity(1),clus.GetX(1),clus.GetY(1)));
	  //      Analyse cluster and decluster if necessary
	  //	
	  ncls++;
	  clus.SetNcluster(1,fNRawClusters);
	  clus.SetClusterType(clus.PhysicsContribution());
	  
	  fNPeaks=0;
	  //
	  //
	  Decluster(&clus);
	  //
	  //      reset Cluster object
	  { // begin local scope
	    for (int k=0;k<clus.GetMultiplicity(0);k++) clus.SetIndex(k, 0, 0);
	  } // end local scope
	  
	  { // begin local scope
	    for (int k=0;k<clus.GetMultiplicity(1);k++) clus.SetIndex(k, 1, 0);
	  } // end local scope
	  
	  clus.SetMultiplicity(0,0);
	  clus.SetMultiplicity(1,0);
	  
	
	} // end loop ndig
    } // end loop cathodes
    delete fDigitMap[0];
    delete fDigitMap[1];
}

Float_t AliMUONClusterFinderVS::SingleMathiesonFit(AliMUONRawCluster *c, Int_t cath)
{
/// Performs a single Mathieson fit on one cathode

    Double_t arglist[20];
    Int_t ierflag=0;
    AliMUONClusterInput& clusterInput = *(AliMUONClusterInput::Instance());
    
    clusterInput.Fitter()->SetFCN(fcnS1);
    clusterInput.Fitter()->mninit(2,10,7);
    clusterInput.Fitter()->SetPrintLevel(-1 + AliLog::GetGlobalDebugLevel());
    arglist[0]=-1;
    clusterInput.Fitter()->mnexcm("SET NOW", arglist, 0, ierflag);
// Set starting values 
    static Double_t vstart[2];
    vstart[0]=c->GetX(1);
    vstart[1]=c->GetY(0);
    
    
// lower and upper limits
    static Double_t lower[2], upper[2];
    Int_t ix,iy, isec;
    fSeg2[cath]->GetPadI(fInput->DetElemId(), c->GetX(cath), c->GetY(cath), fZPlane, ix, iy);
    isec=fSeg2[cath]->Sector(fInput->DetElemId(), ix, iy);

    lower[0]=vstart[0]-fSeg2[cath]->Dpx(fInput->DetElemId(), isec)/2;
    lower[1]=vstart[1]-fSeg2[cath]->Dpy(fInput->DetElemId(), isec)/2;
    
    upper[0]=lower[0]+fSeg2[cath]->Dpx(fInput->DetElemId(), isec);
    upper[1]=lower[1]+fSeg2[cath]->Dpy(fInput->DetElemId(), isec);
    

// step sizes
    static Double_t step[2]={0.0005, 0.0005};
    
    clusterInput.Fitter()->mnparm(0,"x1",vstart[0],step[0],lower[0],upper[0],ierflag);
    clusterInput.Fitter()->mnparm(1,"y1",vstart[1],step[1],lower[1],upper[1],ierflag);
// ready for minimisation	
    arglist[0]= -1;
    arglist[1]= 0;
    
    clusterInput.Fitter()->mnexcm("SET NOGR", arglist, 0, ierflag);
    clusterInput.Fitter()->mnexcm("MIGRAD", arglist, 0, ierflag);
    //    clusterInput.Fitter()->mnexcm("EXIT" , arglist, 0, ierflag);
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

Float_t AliMUONClusterFinderVS::CombiSingleMathiesonFit(AliMUONRawCluster * /*c*/)
{
/// Perform combined Mathieson fit on both cathode planes

    Double_t arglist[20];
    Int_t ierflag=0;
    AliMUONClusterInput& clusterInput = *(AliMUONClusterInput::Instance());
    clusterInput.Fitter()->SetFCN(fcnCombiS1);
    clusterInput.Fitter()->mninit(2,10,7);
    clusterInput.Fitter()->SetPrintLevel(-1 + AliLog::GetGlobalDebugLevel());
    arglist[0]=-1;
    clusterInput.Fitter()->mnexcm("SET NOW", arglist, 0, ierflag);
    static Double_t vstart[2];
    vstart[0]=fXInit[0];
    vstart[1]=fYInit[0];
    
    
// lower and upper limits
    static Float_t lower[2], upper[2];
    Int_t ix,iy,isec;
    Float_t dpy, dpx;

    fSeg2[0]->GetPadI(fInput->DetElemId(), fXInit[0], fYInit[0], fZPlane, ix, iy);
    isec=fSeg2[0]->Sector(fInput->DetElemId(), ix, iy);
    dpy=fSeg2[0]->Dpy(fInput->DetElemId(), isec);
    fSeg2[1]->GetPadI(fInput->DetElemId(), fXInit[0], fYInit[0], fZPlane, ix, iy);
    isec=fSeg2[1]->Sector(fInput->DetElemId(), ix, iy);
    dpx=fSeg2[1]->Dpx(fInput->DetElemId(), isec);
      
    Int_t icount;
    Float_t xdum, ydum, zdum;

//  Find save upper and lower limits    
    
    icount = 0;
    for (fSeg2[1]->FirstPad(fInput->DetElemId(),fXInit[0], fYInit[0], fZPlane, dpx, 0.); 
	 fSeg2[1]->MorePads(fInput->DetElemId()); 
	 fSeg2[1]->NextPad(fInput->DetElemId()))
	{
	  ix=fSeg2[1]->Ix(); iy=fSeg2[1]->Iy();
	  fSeg2[1]->GetPadC(fInput->DetElemId(), ix,iy, upper[0], ydum, zdum);	
	  if (icount ==0) lower[0]=upper[0];
	  icount++;
	}
    
    if (lower[0]>upper[0]) {xdum=lower[0]; lower[0]=upper[0]; upper[0]=xdum;}
	
    icount=0;
    AliDebug(1,Form("\n single y %f %f", fXInit[0], fYInit[0]));
    
    for (fSeg2[0]->FirstPad(fInput->DetElemId(), fXInit[0], fYInit[0], fZPlane, 0., dpy); 
	 fSeg2[0]->MorePads(fInput->DetElemId()); 
	 fSeg2[0]->NextPad(fInput->DetElemId()))
	{
	  ix=fSeg2[0]->Ix(); iy=fSeg2[0]->Iy();
	  fSeg2[0]->GetPadC(fInput->DetElemId(), ix,iy,xdum,upper[1],zdum);	
	  if (icount ==0) lower[1]=upper[1];
	  icount++;
	  AliDebug(1,Form("\n upper lower %d %f %f", icount, upper[1], lower[1]));
	}
    
    if (lower[1]>upper[1]) {xdum=lower[1]; lower[1]=upper[1]; upper[1]=xdum;}

// step sizes
    static Double_t step[2]={0.00001, 0.0001};
    
    clusterInput.Fitter()->mnparm(0,"x1",vstart[0],step[0],lower[0],upper[0],ierflag);
    clusterInput.Fitter()->mnparm(1,"y1",vstart[1],step[1],lower[1],upper[1],ierflag);
// ready for minimisation	
    arglist[0]= -1;
    arglist[1]= 0;
    
    clusterInput.Fitter()->mnexcm("SET NOGR", arglist, 0, ierflag);
    clusterInput.Fitter()->mnexcm("MIGRAD", arglist, 0, ierflag);
    //    clusterInput.Fitter()->mnexcm("EXIT" , arglist, 0, ierflag);
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

Bool_t AliMUONClusterFinderVS::DoubleMathiesonFit(AliMUONRawCluster * /*c*/, Int_t cath)
{
/// Performs a double Mathieson fit on one cathode 

//
//  Initialise global variables for fit
    Double_t arglist[20];
    Int_t ierflag=0;
    AliMUONClusterInput& clusterInput = *(AliMUONClusterInput::Instance());
    clusterInput.Fitter()->SetFCN(fcnS2);
    clusterInput.Fitter()->mninit(5,10,7);
    clusterInput.Fitter()->SetPrintLevel(-1 + AliLog::GetGlobalDebugLevel());
    arglist[0]=-1;
    clusterInput.Fitter()->mnexcm("SET NOW", arglist, 0, ierflag);
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
    Int_t isec;

    isec=fSeg2[cath]->Sector(fInput->DetElemId(),fIx[fIndLocal[0][cath]][cath], 
			     fIy[fIndLocal[0][cath]][cath]);
    lower[0]=vstart[0]-fSeg2[cath]->Dpx(fInput->DetElemId(),isec);
    lower[1]=vstart[1]-fSeg2[cath]->Dpy(fInput->DetElemId(),isec);
    
    upper[0]=lower[0]+2.*fSeg2[cath]->Dpx(fInput->DetElemId(),isec);
    upper[1]=lower[1]+2.*fSeg2[cath]->Dpy(fInput->DetElemId(),isec);
    
    isec=fSeg2[cath]->Sector(fInput->DetElemId(),fIx[fIndLocal[1][cath]][cath], 
			     fIy[fIndLocal[1][cath]][cath]);
    lower[2]=vstart[2]-fSeg2[cath]->Dpx(fInput->DetElemId(),isec)/2;
    lower[3]=vstart[3]-fSeg2[cath]->Dpy(fInput->DetElemId(),isec)/2;
    
    upper[2]=lower[2]+fSeg2[cath]->Dpx(fInput->DetElemId(),isec);
    upper[1]=lower[1]+2.*fSeg2[cath]->Dpy(fInput->DetElemId(),isec);

    

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
    arglist[0]= -1;
    arglist[1]= 0;
    
    clusterInput.Fitter()->mnexcm("SET NOGR", arglist, 0, ierflag);
    clusterInput.Fitter()->mnexcm("MIGRAD", arglist, 0, ierflag);
    //    clusterInput.Fitter()->mnexcm("EXIT" , arglist, 0, ierflag);
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

Float_t AliMUONClusterFinderVS::CombiDoubleMathiesonFit(AliMUONRawCluster * /*c*/)
{
/// Perform combined double Mathieson fit on both cathode planes

    Double_t arglist[20];
    Int_t ierflag=0;
    AliMUONClusterInput& clusterInput = *(AliMUONClusterInput::Instance());
    clusterInput.Fitter()->SetFCN(fcnCombiS2);
    clusterInput.Fitter()->mninit(6,10,7);
    clusterInput.Fitter()->SetPrintLevel(-1 + AliLog::GetGlobalDebugLevel());
    arglist[0]=-1;
    clusterInput.Fitter()->mnexcm("SET NOW", arglist, 0, ierflag);
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

    fSeg2[1]->GetPadI(fInput->DetElemId(),fXInit[0], fYInit[0], fZPlane, ix, iy);
    isec=fSeg2[1]->Sector(fInput->DetElemId(),ix, iy);
    dpx=fSeg2[1]->Dpx(fInput->DetElemId(), isec);

    fSeg2[0]->GetPadI(fInput->DetElemId(), fXInit[0], fYInit[0], fZPlane, ix, iy);
    isec=fSeg2[0]->Sector(fInput->DetElemId(), ix, iy);
    dpy=fSeg2[0]->Dpy(fInput->DetElemId(), isec);

  

    Int_t icount;
    Float_t xdum, ydum, zdum;
    AliDebug(1,Form("\n Cluster Finder: %f %f %f %f  ", fXInit[0], fXInit[1],fYInit[0], fYInit[1] ));

    //  Find save upper and lower limits    
    icount = 0;
    
    for (fSeg2[1]->FirstPad(fInput->DetElemId(),fXInit[0], fYInit[0], fZPlane, dpx, 0.); 
	 fSeg2[1]->MorePads(fInput->DetElemId()); 
	 fSeg2[1]->NextPad(fInput->DetElemId()))
      {
	ix=fSeg2[1]->Ix(); iy=fSeg2[1]->Iy();
	//	if (fDigitMap[1]->TestHit(ix, iy) == kEmpty) continue;
	fSeg2[1]->GetPadC(fInput->DetElemId(),ix,iy,upper[0],ydum,zdum);	
	if (icount ==0) lower[0]=upper[0];
	icount++;
      }
    if (lower[0]>upper[0]) {xdum=lower[0]; lower[0]=upper[0]; upper[0]=xdum;}    
    //    vstart[0] = 0.5*(lower[0]+upper[0]);

    
    icount=0;
    
    for (fSeg2[0]->FirstPad(fInput->DetElemId(),fXInit[0], fYInit[0], fZPlane, 0., dpy); 
	 fSeg2[0]->MorePads(fInput->DetElemId()); 
	 fSeg2[0]->NextPad(fInput->DetElemId()))
      {
	ix=fSeg2[0]->Ix(); iy=fSeg2[0]->Iy();
	//	if (fDigitMap[0]->TestHit(ix, iy) == kEmpty) continue;
	fSeg2[0]->GetPadC(fInput->DetElemId(),ix,iy,xdum,upper[1],zdum);	
	if (icount ==0) lower[1]=upper[1];
	icount++;
      }
    
    if (lower[1]>upper[1]) {xdum=lower[1]; lower[1]=upper[1]; upper[1]=xdum;}    
    //     vstart[1] = 0.5*(lower[1]+upper[1]);


    fSeg2[1]->GetPadI(fInput->DetElemId(),fXInit[1], fYInit[1], fZPlane, ix, iy);
    isec=fSeg2[1]->Sector(fInput->DetElemId(),ix, iy);
    dpx=fSeg2[1]->Dpx(fInput->DetElemId(),isec);
    fSeg2[0]->GetPadI(fInput->DetElemId(),fXInit[1], fYInit[1], fZPlane, ix, iy);
    isec=fSeg2[0]->Sector(fInput->DetElemId(),ix, iy);
    dpy=fSeg2[0]->Dpy(fInput->DetElemId(),isec);


    //  Find save upper and lower limits    

    icount=0;
    
    for (fSeg2[1]->FirstPad(fInput->DetElemId(),fXInit[1], fYInit[1], fZPlane, dpx, 0); 
	 fSeg2[1]->MorePads(fInput->DetElemId()); 
	 fSeg2[1]->NextPad(fInput->DetElemId()))
      {
	ix=fSeg2[1]->Ix(); iy=fSeg2[1]->Iy();
	//	if (fDigitMap[1]->TestHit(ix, iy) == kEmpty) continue;
	fSeg2[1]->GetPadC(fInput->DetElemId(),ix,iy,upper[2],ydum,zdum);	
	if (icount ==0) lower[2]=upper[2];
	icount++;
      }
    if (lower[2]>upper[2]) {xdum=lower[2]; lower[2]=upper[2]; upper[2]=xdum;}    
    //    vstart[2] = 0.5*(lower[2]+upper[2]);

    icount=0;
    
    for (fSeg2[0]->FirstPad(fInput->DetElemId(),fXInit[1], fYInit[1], fZPlane, 0, dpy); 
	 fSeg2[0]-> MorePads(fInput->DetElemId()); 
	 fSeg2[0]->NextPad(fInput->DetElemId()))
      {
	ix=fSeg2[0]->Ix(); iy=fSeg2[0]->Iy();
	//	if (fDigitMap[0]->TestHit(ix, iy) != kEmpty) continue;
	
	fSeg2[0]->GetPadC(fInput->DetElemId(),ix,iy,xdum,upper[3],zdum);	
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
    arglist[0]= -1;
    arglist[1]= 0;
    
    clusterInput.Fitter()->mnexcm("SET NOGR", arglist, 0, ierflag);
    clusterInput.Fitter()->mnexcm("MIGRAD", arglist, 0, ierflag);
    //    clusterInput.Fitter()->mnexcm("EXIT" , arglist, 0, ierflag);
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
/// One cluster for each maximum

    Int_t i, j, cath;
    AliMUONClusterInput& clusterInput = *(AliMUONClusterInput::Instance());
    for (j=0; j<2; j++) {
	AliMUONRawCluster cnew;
	cnew.SetGhost(c->GetGhost());
	for (cath=0; cath<2; cath++) {
	    cnew.SetChi2(cath,fChi2[0]);
	    // ?? why not cnew.fChi2[cath]=fChi2[cath];
	    
	    if (fNPeaks == 0) {
		cnew.SetNcluster(0,-1);
		cnew.SetNcluster(1,fNRawClusters);
	    } else {
		cnew.SetNcluster(0,fNPeaks);
		cnew.SetNcluster(1,0);
	    }
	    cnew.SetMultiplicity(cath,0);
	    cnew.SetX(cath, Float_t(fXFit[j]));
	    cnew.SetY(cath, Float_t(fYFit[j]));
	    cnew.SetZ(cath, fZPlane);
	    if (j==0) {
		cnew.SetCharge(cath, Int_t(clusterInput.TotalCharge(cath)*fQrFit[cath]));
	    } else {
		cnew.SetCharge(cath, Int_t(clusterInput.TotalCharge(cath)*(1-fQrFit[cath])));
	    }
	    fSeg2[cath]->SetHit(fInput->DetElemId(), fXFit[j],fYFit[j],fZPlane);

	    for (i=0; i<fMul[cath]; i++) {
	      Float_t q1;
		cnew.SetIndex(cnew.GetMultiplicity(cath), cath, c->GetIndex(i,cath));

		fSeg2[cath]->SetPad(fInput->DetElemId(),fIx[i][cath], fIy[i][cath]);
		q1 = fInput->Mathieson()->IntXY(fInput->DetElemId(),fSeg2[cath]);
		
		cnew.SetContrib(i, cath, q1*Float_t(cnew.GetCharge(cath))/Float_t(fQ[i][cath]));
		cnew.SetMultiplicity(cath, cnew.GetMultiplicity(cath)+1 );
	    }
	    FillCluster(&cnew,0,cath);
	} // cathode loop
	cnew.SetClusterType(cnew.PhysicsContribution());
	if (cnew.GetCharge(0)>0 && cnew.GetCharge(1)>0) AddRawCluster(cnew);
	fNPeaks++;
    }
}
void AliMUONClusterFinderVS::AddRawCluster(AliMUONRawCluster& c)
{
/// Add a raw cluster copy to the list

  //     AliMUON *pMUON=(AliMUON*)gAlice->GetModule("MUON");
  //     pMUON->GetMUONData()->AddRawCluster(fInput->Chamber(),c); 
  //     fNRawClusters++;
  
  // Setting detection element in raw cluster for alignment
  // BB 19/05/05
  c.SetDetElemId(fInput->DetElemId());
  
  TClonesArray &lrawcl = *fRawClusters;
  new(lrawcl[fNRawClusters++]) AliMUONRawCluster(c);
  AliDebug(1,Form("\nfNRawClusters %d\n",fNRawClusters));
}

AliMUONClusterFinderVS& AliMUONClusterFinderVS
::operator = (const AliMUONClusterFinderVS& rhs)
{
// Protected assignement operator

  if (this == &rhs) return *this;

  AliFatal("Not implemented.");
    
  return *this;  
}

//
// Minimisation functions
// Single Mathieson
void fcnS1(Int_t & /*npar*/, Double_t * /*gin*/, Double_t &f, Double_t *par, Int_t /*iflag*/)
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

void fcnCombiS1(Int_t & /*npar*/, Double_t * /*gin*/, Double_t &f, Double_t *par, Int_t /*iflag*/)
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
void fcnS2(Int_t & /*npar*/, Double_t * /*gin*/, Double_t &f, Double_t *par, Int_t /*iflag*/)
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
void fcnCombiS2(Int_t & /*npar*/, Double_t * /*gin*/, Double_t &f, Double_t *par, Int_t /*iflag*/)
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
