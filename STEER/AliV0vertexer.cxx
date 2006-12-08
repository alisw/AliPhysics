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
//     Origin: Iouri Belikov, IReS, Strasbourg, Jouri.Belikov@cern.ch
//-------------------------------------------------------------------------


#include <TObjArray.h>
#include <TTree.h>

#include "AliESD.h"
#include "AliESDv0.h"
#include "AliESDtrack.h"
#include "AliV0vertexer.h"

ClassImp(AliV0vertexer)


//A set of very loose cuts 
Double_t AliV0vertexer::fgChi2max=33.; //max chi2
Double_t AliV0vertexer::fgDNmin=0.05;  //min imp parameter for the 1st daughter
Double_t AliV0vertexer::fgDPmin=0.05;  //min imp parameter for the 2nd daughter
Double_t AliV0vertexer::fgDCAmax=0.5;  //max DCA between the daughter tracks
Double_t AliV0vertexer::fgCPAmax=0.99; //max cosine of V0's pointing angle
Double_t AliV0vertexer::fgRmin=0.2;    //min radius of the fiducial volume
Double_t AliV0vertexer::fgRmax=100.;   //max radius of the fiducial volume

Int_t AliV0vertexer::Tracks2V0vertices(AliESD *event) {
  //--------------------------------------------------------------------
  //This function reconstructs V0 vertices
  //--------------------------------------------------------------------
   const AliESDVertex *vtx=event->GetVertex();
   Double_t xv=vtx->GetXv(), yv=vtx->GetYv(), zv=vtx->GetZv();

   Int_t nentr=event->GetNumberOfTracks();
   Double_t b=event->GetMagneticField();

   if (nentr<2) return 0; 

   TArrayI neg(nentr/2);
   TArrayI pos(nentr/2);

   Int_t nneg=0, npos=0, nvtx=0;

   Int_t i;
   for (i=0; i<nentr; i++) {
     AliESDtrack *esdTrack=event->GetTrack(i);
     UInt_t status=esdTrack->GetStatus();
     UInt_t flags=AliESDtrack::kITSin|AliESDtrack::kTPCin|
                  AliESDtrack::kTPCpid|AliESDtrack::kESDpid;

     if ((status&AliESDtrack::kITSrefit)==0)
        if (flags!=status) continue;

     Double_t d=esdTrack->GetD(xv,yv,b);
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

         if (TMath::Abs(ntrk->GetD(xv,yv,b))<fDNmin)
	   if (TMath::Abs(ptrk->GetD(xv,yv,b))<fDNmin) continue;

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

         AliESDv0 vertex(nt,nidx,pt,pidx);
         if (vertex.GetChi2V0() > fChi2max) continue;
	 
	 Float_t cpa=vertex.GetV0CosineOfPointingAngle(xv,yv,zv);
	 if (cpa < fCPAmax) continue;
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














