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

//-------------------------------------------------------------------------
//               Implementation of the V0 vertexer class
//                  reads tracks writes out V0 vertices
//                      fills the ESD with the V0s       
//     Origin: Iouri Belikov, IPHC, Strasbourg, Jouri.Belikov@cern.ch
//-------------------------------------------------------------------------


#include "AliESDEvent.h"
#include "AliESDv0.h"
#include "AliV0vertexer.h"
#include "AliV0HypSel.h"
#include "AliLog.h"

ClassImp(AliV0vertexer)


//A set of very loose cuts
Double_t AliV0vertexer::fgEtaMax =5.0; // max eta
Double_t AliV0vertexer::fgChi2max=33.; //max chi2
Double_t AliV0vertexer::fgDNmin=0.05;  //min imp parameter for the 1st daughter
Double_t AliV0vertexer::fgDPmin=0.05;  //min imp parameter for the 2nd daughter
Double_t AliV0vertexer::fgDCAmax=1.5;  //max DCA between the daughter tracks
Double_t AliV0vertexer::fgCPAmin=0.9;  //min cosine of V0's pointing angle
Double_t AliV0vertexer::fgRmin=0.2;    //min radius of the fiducial volume
Double_t AliV0vertexer::fgRmax=200.;   //max radius of the fiducial volume

Int_t AliV0vertexer::Tracks2V0vertices(AliESDEvent *event) {
  //--------------------------------------------------------------------
  //This function reconstructs V0 vertices
  //--------------------------------------------------------------------
  
  const AliESDVertex *vtxT3D=event->GetPrimaryVertex();
  
  Double_t xPrimaryVertex=vtxT3D->GetX();
  Double_t yPrimaryVertex=vtxT3D->GetY();
  Double_t zPrimaryVertex=vtxT3D->GetZ();
  
  Int_t nentr=event->GetNumberOfTracks();
  Double_t b=event->GetMagneticField();
  
  if (nentr<2) return 0;
  
  AliV0HypSel::AccountBField(b);
  
  TArrayI neg(nentr);
  TArrayI pos(nentr);
  
  Int_t nneg=0, npos=0, nvtx=0;
  
  Int_t i;
  for (i=0; i<nentr; i++) {
    AliESDtrack *esdTrack=event->GetTrack(i);
    ULong64_t status=esdTrack->GetStatus();
    
    //if ((status&AliESDtrack::kITSrefit)==0)//not to accept the ITS SA tracks
    if ((status&AliESDtrack::kTPCrefit)==0) continue;
    
    Double_t d=0.0;
    if (esdTrack->GetSign() < 0.){
      d=esdTrack->GetD(xPrimaryVertex,yPrimaryVertex,b);
      if (TMath::Abs(d)<fDNmin) continue;
      if (TMath::Abs(d)>fRmax) continue;
      neg[nneg++]=i;
    }else{
      d=esdTrack->GetD(xPrimaryVertex,yPrimaryVertex,b);
      if (TMath::Abs(d)<fDPmin) continue;
      if (TMath::Abs(d)>fRmax) continue;
      pos[npos++]=i;
    }
  }
  
  int nHypSel = fV0HypSelArray ? fV0HypSelArray->GetEntriesFast() : 0;
  
  for (i=0; i<nneg; i++) {
    Int_t nidx=neg[i];
    AliESDtrack *ntrk=event->GetTrack(nidx);
    
    for (Int_t k=0; k<npos; k++) {
      Int_t pidx=pos[k];
      AliESDtrack *ptrk=event->GetTrack(pidx);
      
      //This is fully redundant with respect to the first loop, avoid for speed -> OK ?
      //if (TMath::Abs(ntrk->GetD(xPrimaryVertex,yPrimaryVertex,b))<fDNmin)
      //if (TMath::Abs(ptrk->GetD(xPrimaryVertex,yPrimaryVertex,b))<fDNmin) continue
      
      AliExternalTrackParam nt(*ntrk), pt(*ptrk), *ntp=&nt, *ptp=&pt;
      
      Double_t xn, xp;
      if(fUseImprovedFinding && Preoptimize(ntp,ptp,&xn,&xp,b)){
        //Move tracks to a better position if that helps
        nt.PropagateTo(xn,b); pt.PropagateTo(xp,b);
      }
      Double_t dca=ntp->GetDCA(ptp,b,xn,xp);
      if (dca > fDCAmax) continue;
      if ((xn+xp) > 2*fRmax) continue;
      if ((xn+xp) < 2*fRmin) continue;
      
      /* check if this can be removed with Iouri ?
       this does absolutely nothing, only repeats GetDCA call
       Bool_t corrected=kFALSE;
       if ((nt.GetX() > 3.) && (xn < 3.)) {
       //correct for the beam pipe material
       corrected=kTRUE;
       }
       if ((pt.GetX() > 3.) && (xp < 3.)) {
       //correct for the beam pipe material
       corrected=kTRUE;
       }
       if (corrected) {
       dca=nt.GetDCA(&pt,b,xn,xp);
       if (dca > fDCAmax) continue;
       if ((xn+xp) > 2*fRmax) continue;
       if ((xn+xp) < 2*fRmin) continue;
       }
       */
      
      nt.PropagateTo(xn,b); pt.PropagateTo(xp,b);
      
      AliESDv0 vertex(nt,nidx,pt,pidx);
      if(fUseImprovedFinding) vertex.Refit(); //imp pos + cov mat
      
      if (vertex.GetChi2V0() > fChi2max) continue;
      
      Double_t x=vertex.Xv(), y=vertex.Yv();
      Double_t r2=x*x + y*y;
      if (r2 < fRmin*fRmin) continue;
      if (r2 > fRmax*fRmax) continue;
      
      Float_t cpa=vertex.GetV0CosineOfPointingAngle(xPrimaryVertex,yPrimaryVertex,zPrimaryVertex);
      const Double_t pThr=1.5;
      Double_t pv0=vertex.P();
      
      //Warning; this needs to be removed when loosening cuts ... FIXME!
      if (pv0<pThr) {
        //Below the threshold "pThr", try a momentum dependent cos(PA) cut
        const Double_t bend=0.03; // approximate Xi bending angle
        const Double_t qt=0.211;  // max Lambda pT in Omega decay
        const Double_t cpaThr=TMath::Cos(TMath::ATan(qt/pThr) + bend);
        Double_t
        cpaCut=(fCPAmin/cpaThr)*TMath::Cos(TMath::ATan(qt/pv0) + bend);
        if (cpa < cpaCut) continue;
      } else
      if (cpa < fCPAmin) continue;
      
      if (nHypSel) { // do we select particular hypthesese?
        Bool_t reject = kTRUE;
        float pt = vertex.Pt();
        for (int ih=0;ih<nHypSel;ih++) {
          const AliV0HypSel* hyp = (const AliV0HypSel*)(*fV0HypSelArray)[ih];
          double m = vertex.GetEffMassExplicit(hyp->GetM0(),hyp->GetM1());
          if (TMath::Abs(m - hyp->GetMass())<hyp->GetMassMargin(pt)) {
            reject = kFALSE;
            break;
          }
        }
        if (reject) continue;
      }
      
      vertex.SetDcaV0Daughters(dca);
      vertex.SetV0CosineOfPointingAngle(cpa);
      vertex.ChangeMassHypothesis(kK0Short);
      
      event->AddV0(&vertex);
      
      nvtx++;
    }
  }
  
  Info("Tracks2V0vertices","Number of reconstructed V0 vertices: %d",nvtx);
  
  return nvtx;
}


//________________________________________________
void AliV0vertexer::SetV0HypSel(const TObjArray* selArr)
{
  if (!selArr || !selArr->GetEntriesFast()) {
    AliInfo("No V0 hypothesis selection will be performed");
    return;
  }
  fV0HypSelArray = selArr;
}

//________________________________________________
void AliV0vertexer::SetUseImprovedFinding(const Bool_t lInput)
{
  fUseImprovedFinding = lInput;
}

//________________________________________________
void AliV0vertexer::GetHelixCenter(const AliExternalTrackParam *track,Double_t center[2], const Double_t b)
{
  // Copied from AliV0ReaderV1::GetHelixCenter
  // Get Center of the helix track parametrization
  
  Int_t charge=track->Charge();
  
  Double_t    helix[6];
  track->GetHelixParameters(helix,b);
  
  Double_t xpos =    helix[5];
  Double_t ypos =    helix[0];
  Double_t radius = TMath::Abs(1./helix[4]);
  Double_t phi = helix[2];
  if(phi < 0){
    phi = phi + 2*TMath::Pi();
  }
  phi -= TMath::Pi()/2.;
  Double_t xpoint =    radius * TMath::Cos(phi);
  Double_t ypoint =    radius * TMath::Sin(phi);
  if(b<0&&charge > 0){
    xpoint = - xpoint;
    ypoint = - ypoint;
  }
  if(b>0 && charge < 0){
    xpoint = - xpoint;
    ypoint = - ypoint;
  }
  center[0] =    xpos + xpoint;
  center[1] =    ypos + ypoint;
  return;
}

//________________________________________________
Bool_t AliV0vertexer::Preoptimize(const AliExternalTrackParam *nt, AliExternalTrackParam *pt, Double_t *lPreprocessxn, Double_t *lPreprocessxp, const Double_t b)
//This function pre-optimizes a two-track pair in the XY plane
//and provides two X values for the tracks if successful
{
  Double_t lMinimumX = -3.0 ;
  Double_t lMaximumX = 300.0;
  
  Double_t nhelix[6], phelix[6];
  nt->GetHelixParameters(nhelix,b);
  pt->GetHelixParameters(phelix,b);
  Double_t lNegCenterR[2], lPosCenterR[2];
  
  //Negative track parameters in XY
  GetHelixCenter( nt , lNegCenterR, b);
  Double_t xNegCenter = lNegCenterR[0];
  Double_t yNegCenter = lNegCenterR[1];
  Double_t NegRadius = TMath::Abs(1./nhelix[4]);
  
  //Positive track parameters in XY
  GetHelixCenter( pt , lPosCenterR, b );
  Double_t xPosCenter = lPosCenterR[0];
  Double_t yPosCenter = lPosCenterR[1];
  Double_t PosRadius = TMath::Abs(1./phelix[4]);
  
  //Define convenient coordinate system
  //Logical zero: position of negative center
  Double_t ux = xPosCenter - xNegCenter;
  Double_t uy = yPosCenter - yNegCenter;
  
  //Check center-to-center distance
  Double_t lDist = TMath::Sqrt(
                               TMath::Power( xNegCenter - xPosCenter , 2) +
                               TMath::Power( yNegCenter - yPosCenter , 2)
                               );
  //Normalize ux, uz to unit vector
  ux /= lDist; uy /= lDist;
  
  //Calculate perpendicular vector (normalized)
  Double_t vx = -uy;
  Double_t vy = +ux;
  
  Double_t lPreprocessDCAxy = 1e+3; //define outside scope
  *lPreprocessxp = pt->GetX(); //start at current location
  *lPreprocessxn = nt->GetX(); //start at current location
  
  //============================================================
  //Pre-optimization in the XY plane: cases considered here
  //============================================================
  //
  //  Case 1: Circles do not touch, centers far away
  //          (D > R1 + R2)
  //
  //  Case 2: Circles touch, centers at reasonable distance wrt D
  //          (D < R1 + R2) && (D > |R1-R2|)
  //
  //  Case 3: Circles do not touch, one inside the other
  //          (D < |R1-R2|)
  //
  //  Cases 1 and 2 are treated. Case 3 is not treated (unlikely
  //  to be a problem with unlike-sign charged tracks): brute
  //  force minimization takes place in any case
  //
  //============================================================
  
  //______________________
  //CASE 1
  if( (lDist > NegRadius + PosRadius) ){
    //================================================================
    //Case 1: distance bigger than sum of radii ("gamma-like")
    //        re-position tracks along the center-to-center axis
    //Re-position negative track
    Double_t xNegOptPosition = xNegCenter + NegRadius*ux;
    Double_t yNegOptPosition = yNegCenter + NegRadius*uy;
    Double_t csNeg=TMath::Cos(nt->GetAlpha());
    Double_t snNeg=TMath::Sin(nt->GetAlpha());
    Double_t xThisNeg=xNegOptPosition*csNeg + yNegOptPosition*snNeg;
    
    //Re-position positive track
    Double_t xPosOptPosition = xPosCenter - PosRadius*ux;
    Double_t yPosOptPosition = yPosCenter - PosRadius*uy;
    Double_t csPos=TMath::Cos(pt->GetAlpha());
    Double_t snPos=TMath::Sin(pt->GetAlpha());
    Double_t xThisPos=xPosOptPosition*csPos + yPosOptPosition*snPos;
    
    if( xThisNeg < lMaximumX && xThisPos < lMaximumX && xThisNeg > lMinimumX && xThisPos > lMinimumX) {
      Double_t lCase1NegR[3]={0.},lCase1PosR[3]={0.};
      if (nt->GetXYZAt(xThisNeg,b, lCase1NegR) && pt->GetXYZAt(xThisPos,b, lCase1PosR)) {
	lPreprocessDCAxy = TMath::Sqrt(
				       TMath::Power(lCase1NegR[0]-lCase1PosR[0],2)+
				       TMath::Power(lCase1NegR[1]-lCase1PosR[1],2)+
				       TMath::Power(lCase1NegR[2]-lCase1PosR[2],2)
				       );
	//Pass coordinates
	if( lPreprocessDCAxy<999){
	  *lPreprocessxp = xThisPos;
	  *lPreprocessxn = xThisNeg;
	}
      }
    }
    //================================================================
  }
  
  //______________________
  //CASE 2
  if( (lDist > TMath::Abs(NegRadius-PosRadius)) && (lDist < NegRadius + PosRadius) ){
    //================================================================
    //Case 2: distance smaller than sum of radii (cowboy/sailor configs)
    
    //Calculate coordinate for radical line
    Double_t lRadical = (lDist*lDist - PosRadius*PosRadius + NegRadius*NegRadius) / (2*lDist);
    
    //Calculate absolute displacement from center-to-center axis
    Double_t lDisplace = (0.5/lDist) * TMath::Sqrt(
                                                   (-lDist + PosRadius - NegRadius) *
                                                   (-lDist - PosRadius + NegRadius) *
                                                   (-lDist + PosRadius + NegRadius) *
                                                   ( lDist + PosRadius + NegRadius)
                                                   );
    
    //3D distances in the two cases studied (prefer smallest)
    Double_t lCase2aDCA = 1e+3;
    Double_t lCase2bDCA = 1e+3;
    
    //2 cases: positive and negative displacement
    Double_t xNegOptPosition[2], yNegOptPosition[2], xPosOptPosition[2], yPosOptPosition[2];
    Double_t csNeg, snNeg, csPos, snPos;
    Double_t xThisNeg[2], xThisPos[2];
    
    csNeg=TMath::Cos(nt->GetAlpha());
    snNeg=TMath::Sin(nt->GetAlpha());
    csPos=TMath::Cos(pt->GetAlpha());
    snPos=TMath::Sin(pt->GetAlpha());
    
    //Case 2a: Positive displacement along v vector
    //Re-position negative track
    xNegOptPosition[0] = xNegCenter + lRadical*ux + lDisplace*vx;
    yNegOptPosition[0] = yNegCenter + lRadical*uy + lDisplace*vy;
    xThisNeg[0] = xNegOptPosition[0]*csNeg + yNegOptPosition[0]*snNeg;
    //Re-position positive track
    xPosOptPosition[0] = xNegCenter + lRadical*ux + lDisplace*vx;
    yPosOptPosition[0] = yNegCenter + lRadical*uy + lDisplace*vy;
    xThisPos[0] = xPosOptPosition[0]*csPos + yPosOptPosition[0]*snPos;
    
    //Case 2b: Negative displacement along v vector
    //Re-position negative track
    xNegOptPosition[1] = xNegCenter + lRadical*ux - lDisplace*vx;
    yNegOptPosition[1] = yNegCenter + lRadical*uy - lDisplace*vy;
    xThisNeg[1] = xNegOptPosition[1]*csNeg + yNegOptPosition[1]*snNeg;
    //Re-position positive track
    xPosOptPosition[1] = xNegCenter + lRadical*ux - lDisplace*vx;
    yPosOptPosition[1] = yNegCenter + lRadical*uy - lDisplace*vy;
    xThisPos[1] = xPosOptPosition[1]*csPos + yPosOptPosition[1]*snPos;
    
    //Test the two cases, please
    
    //Case 2a
    if( xThisNeg[0] < lMaximumX && xThisPos[0] < lMaximumX && xThisNeg[0] > lMinimumX && xThisPos[0] > lMinimumX ){
      Double_t lCase2aNegR[3]={0.}, lCase2aPosR[3]={0.};
      if (nt->GetXYZAt(xThisNeg[0],b, lCase2aNegR) && pt->GetXYZAt(xThisPos[0],b, lCase2aPosR)) {
	lCase2aDCA = TMath::Sqrt(
				 TMath::Power(lCase2aNegR[0]-lCase2aPosR[0],2)+
				 TMath::Power(lCase2aNegR[1]-lCase2aPosR[1],2)+
				 TMath::Power(lCase2aNegR[2]-lCase2aPosR[2],2)
				 );
      }
    }
    
    //Case 2b
    if( xThisNeg[1] < lMaximumX && xThisPos[1] < lMaximumX && xThisNeg[1] > lMinimumX && xThisPos[1] > lMinimumX ){      
      Double_t lCase2bNegR[3]={0.}, lCase2bPosR[3]={0.};
      if (nt->GetXYZAt(xThisNeg[1],b, lCase2bNegR) && pt->GetXYZAt(xThisPos[1],b, lCase2bPosR)) {
	lCase2bDCA = TMath::Sqrt(
				 TMath::Power(lCase2bNegR[0]-lCase2bPosR[0],2)+
				 TMath::Power(lCase2bNegR[1]-lCase2bPosR[1],2)+
				 TMath::Power(lCase2bNegR[2]-lCase2bPosR[2],2)
				 );
      }
    }
    //Minor detail: all things being equal, prefer closest X
    Double_t lCase2aSumX = xThisPos[0]+xThisNeg[0];
    Double_t lCase2bSumX = xThisPos[1]+xThisNeg[1];
    
    Double_t lDCAxySmallestR = lCase2aDCA;
    Double_t lxpSmallestR = xThisPos[0];
    Double_t lxnSmallestR = xThisNeg[0];
    
    Double_t lDCAxyLargestR = lCase2bDCA;
    Double_t lxpLargestR = xThisPos[1];
    Double_t lxnLargestR = xThisNeg[1];
    
    if( lCase2bSumX+1e-6 < lCase2aSumX ){
      lDCAxySmallestR = lCase2bDCA;
      lxpSmallestR = xThisPos[1];
      lxnSmallestR = xThisNeg[1];
      lDCAxyLargestR = lCase2aDCA;
      lxpLargestR = xThisPos[0];
      lxnLargestR = xThisNeg[0];
    }
    
    //Pass conclusion to lPreprocess variables, please
    lPreprocessDCAxy = lDCAxySmallestR;
    *lPreprocessxp = lxpSmallestR;
    *lPreprocessxn = lxnSmallestR;
    if( lDCAxyLargestR+1e-6 < lDCAxySmallestR ){ //beware epsilon: numerical calculations are unstable here
      lPreprocessDCAxy = lDCAxyLargestR;
      *lPreprocessxp = lxpLargestR;
      *lPreprocessxn = lxnLargestR;
    }
    //Protection against something too crazy, please
    if( lPreprocessDCAxy>999){
      *lPreprocessxp = pt->GetX(); //start at current location
      *lPreprocessxn = nt->GetX(); //start at current location
    }
    
  }
  //End of preprocessing stage!
  //at this point lPreprocessxp, lPreprocessxn are already good starting points: update helixparams
  Bool_t lWorked = kFALSE;
  if( lPreprocessDCAxy < 999 ) lWorked = kTRUE;
  return lWorked;
}
