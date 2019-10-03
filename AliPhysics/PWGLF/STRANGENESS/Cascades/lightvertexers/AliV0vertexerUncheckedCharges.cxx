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
//                Further adapted to allow for like-sign
//                     daughter track combination
//
//          This is still being tested! Use at your own risk!
//-------------------------------------------------------------------------

#include "AliESDEvent.h"
#include "AliESDv0.h"
#include "AliV0vertexerUncheckedCharges.h"

ClassImp(AliV0vertexerUncheckedCharges)


//A set of very loose cuts
Double_t AliV0vertexerUncheckedCharges::fgChi2max=33.; //max chi2
Double_t AliV0vertexerUncheckedCharges::fgDNmin=0.05;  //min imp parameter for the 1st daughter
Double_t AliV0vertexerUncheckedCharges::fgDPmin=0.05;  //min imp parameter for the 2nd daughter
Double_t AliV0vertexerUncheckedCharges::fgDCAmax=1.5;  //max DCA between the daughter tracks
Double_t AliV0vertexerUncheckedCharges::fgCPAmin=0.9;  //min cosine of V0's pointing angle
Double_t AliV0vertexerUncheckedCharges::fgRmin=0.2;    //min radius of the fiducial volume
Double_t AliV0vertexerUncheckedCharges::fgRmax=200.;   //max radius of the fiducial volume

Double_t AliV0vertexerUncheckedCharges::fgMaxEta=0.8;        //max |eta|
Double_t AliV0vertexerUncheckedCharges::fgMinClusters=70;   //min clusters (>=)

Int_t AliV0vertexerUncheckedCharges::Tracks2V0vertices(AliESDEvent *event) {
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
    
    TArrayI trackarray(nentr);
    
    Int_t ntracks=0, nvtx=0;
    
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
        
        //Disregard charges
        trackarray[ntracks++]=i;
    }
    
    
    for (i=0; i<ntracks; i++) {
        //originally: negative (now track 1)
        Int_t  idx1=trackarray[i];
        AliESDtrack *trk1=event->GetTrack(idx1);
        
        for (Int_t k=0; k<ntracks; k++) {
            if( i==k ) continue; //don't combine a track with itself, please
            
            //originally: positive (now track 2) 
            Int_t idx2=trackarray[k];
            AliESDtrack *trk2=event->GetTrack(idx2);
            
            
            if (TMath::Abs(trk1->GetD(xPrimaryVertex,yPrimaryVertex,b))<fDNmin)
                if (TMath::Abs(trk2->GetD(xPrimaryVertex,yPrimaryVertex,b))<fDNmin) continue;
            
            Double_t xn, xp, dca=trk1->GetDCA(trk2,b,xn,xp);
            if (dca > fDCAmax) continue;
            if ((xn+xp) > 2*fRmax) continue;
            if ((xn+xp) < 2*fRmin) continue;
            
            AliExternalTrackParam t1(*trk1), t2(*trk2);
            Bool_t corrected=kFALSE;
            if ((t1.GetX() > 3.) && (xn < 3.)) {
                //correct for the beam pipe material
                corrected=kTRUE;
            }
            if ((t2.GetX() > 3.) && (xp < 3.)) {
                //correct for the beam pipe material
                corrected=kTRUE;
            }
            if (corrected) {
                dca=t1.GetDCA(&t2,b,xn,xp);
                if (dca > fDCAmax) continue;
                if ((xn+xp) > 2*fRmax) continue;
                if ((xn+xp) < 2*fRmin) continue;
            }
            
            t1.PropagateTo(xn,b); t2.PropagateTo(xp,b);
            
            //select maximum eta range (after propagation)
            if (TMath::Abs(t1.Eta())>fMaxEta) continue;
            if (TMath::Abs(t2.Eta())>fMaxEta) continue;
            
            AliESDv0 vertex(t1,idx1,t2,idx2);
            
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














