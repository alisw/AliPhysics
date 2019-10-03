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
//                Adapted for more customizability and
//                   lightweight operation in Pb-Pb
//
//          This is still being tested! Use at your own risk!
//-------------------------------------------------------------------------

#include "AliESDEvent.h"
#include "AliESDv0.h"
#include "AliLightV0vertexer.h"

ClassImp(AliLightV0vertexer)


//A set of very loose cuts
Double_t AliLightV0vertexer::fgChi2max=33.; //max chi2
Double_t AliLightV0vertexer::fgDNmin=0.05;  //min imp parameter for the 1st daughter
Double_t AliLightV0vertexer::fgDPmin=0.05;  //min imp parameter for the 2nd daughter
Double_t AliLightV0vertexer::fgDCAmax=1.5;  //max DCA between the daughter tracks
Double_t AliLightV0vertexer::fgCPAmin=0.9;  //min cosine of V0's pointing angle
Double_t AliLightV0vertexer::fgRmin=0.2;    //min radius of the fiducial volume
Double_t AliLightV0vertexer::fgRmax=200.;   //max radius of the fiducial volume

Double_t AliLightV0vertexer::fgMaxEta=0.8;        //max |eta|
Double_t AliLightV0vertexer::fgMinClusters=70;   //min clusters (>=)

Int_t AliLightV0vertexer::Tracks2V0vertices(AliESDEvent *event) {
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
    
    TArrayI neg(nentr);
    TArrayI pos(nentr);
    
    Int_t nneg=0, npos=0, nvtx=0;
    
    Int_t i;
    for (i=0; i<nentr; i++) {
        AliESDtrack *esdTrack=event->GetTrack(i);
        ULong_t status=esdTrack->GetStatus();
        
        //if ((status&AliESDtrack::kITSrefit)==0)//not to accept the ITS SA tracks
        if ((status&AliESDtrack::kTPCrefit)==0) continue;
        
        //Track pre-selection: clusters
        if (esdTrack->GetTPCNcls() < fMinClusters ) continue;
        
        Double_t d=esdTrack->GetD(xPrimaryVertex,yPrimaryVertex,b);
        if (TMath::Abs(d)<fDPmin) continue;
        if (TMath::Abs(d)>fRmax) continue;
        
        if (esdTrack->GetSign() < 0.) neg[nneg++]=i;
        else pos[npos++]=i;
    }
    
    
    for (i=0; i<nneg; i++) {
        Int_t nidx=neg[i];
        AliESDtrack *ntrk=event->GetTrack(nidx);
        
        for (Int_t k=0; k<npos; k++) {
            Int_t pidx=pos[k];
            AliESDtrack *ptrk=event->GetTrack(pidx);
            
            //Track pre-selection: clusters
            if (ptrk->GetTPCNcls() < fMinClusters ) continue;
            
            if (TMath::Abs(ntrk->GetD(xPrimaryVertex,yPrimaryVertex,b))<fDNmin)
                if (TMath::Abs(ptrk->GetD(xPrimaryVertex,yPrimaryVertex,b))<fDNmin) continue;
            
            Double_t xn, xp, dca=ntrk->GetDCA(ptrk,b,xn,xp);
            if (dca > fDCAmax) continue;
            if ((xn+xp) > 2*fRmax) continue;
            if ((xn+xp) < 2*fRmin) continue;
            
            AliExternalTrackParam nt(*ntrk), pt(*ptrk);
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
            
            nt.PropagateTo(xn,b); pt.PropagateTo(xp,b);
            
            //select maximum eta range (after propagation)
            if (TMath::Abs(nt.Eta())>fMaxEta) continue;
            if (TMath::Abs(pt.Eta())>fMaxEta) continue;
            
            AliESDv0 vertex(nt,nidx,pt,pidx);
            
            //Experimental: refit V0 if asked to do so 
            if( fkDoRefit ) vertex.Refit();
            
            //No selection: it was not previously applied, don't  apply now. 
            //if (vertex.GetChi2V0() > fChi2max) continue;
            
            Double_t x=vertex.Xv(), y=vertex.Yv();
            Double_t r2=x*x + y*y;
            if (r2 < fRmin*fRmin) continue;
            if (r2 > fRmax*fRmax) continue;
            
            Float_t cpa=vertex.GetV0CosineOfPointingAngle(xPrimaryVertex,yPrimaryVertex,zPrimaryVertex);
            
            //Simple cosine cut (no pt dependence for now)
            if (cpa < fCPAmin) continue;
            
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














