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
#include "AliESDtrack.h"
#include "AliITStrackV2.h"
#include "AliV0vertex.h"
#include "AliV0vertexer.h"

ClassImp(AliV0vertexer)

Int_t AliV0vertexer::Tracks2V0vertices(AliESD *event) {
  //--------------------------------------------------------------------
  //This function reconstructs V0 vertices
  //--------------------------------------------------------------------

   Int_t nentr=event->GetNumberOfTracks();

   TObjArray negtrks(nentr/2);
   TObjArray postrks(nentr/2);

   Int_t nneg=0, npos=0, nvtx=0;

   Int_t i;
   for (i=0; i<nentr; i++) {
     AliESDtrack *esd=event->GetTrack(i);
     UInt_t status=esd->GetStatus();
     UInt_t flags=AliESDtrack::kITSin|AliESDtrack::kTPCin|
                  AliESDtrack::kTPCpid|AliESDtrack::kESDpid;

     if ((status&AliESDtrack::kITSrefit)==0)
        if (flags!=status) continue;

     AliITStrackV2 *iotrack=new AliITStrackV2(*esd);
     iotrack->SetLabel(i);  // now it is the index in array of ESD tracks
     if ((status&AliESDtrack::kITSrefit)==0)   //correction for the beam pipe
        if (!iotrack->PropagateTo(3.,0.0023,65.19)) {
	  delete iotrack;
	  continue;
	}
     if (!iotrack->PropagateTo(2.5,0.,0.)) {
       delete iotrack;
       continue;
     }

     if (iotrack->Get1Pt() < 0.) {nneg++; negtrks.AddLast(iotrack);}
     else {npos++; postrks.AddLast(iotrack);}
   }   


   for (i=0; i<nneg; i++) {
      //if (i%10==0) cerr<<nneg-i<<'\r';
      AliITStrackV2 *ntrk=(AliITStrackV2 *)negtrks.UncheckedAt(i);

      if (TMath::Abs(ntrk->GetD(fX,fY))<fDPmin) continue;
      if (TMath::Abs(ntrk->GetD(fX,fY))>fRmax) continue;

      for (Int_t k=0; k<npos; k++) {
         AliITStrackV2 *ptrk=(AliITStrackV2 *)postrks.UncheckedAt(k);

         if (TMath::Abs(ptrk->GetD(fX,fY))<fDPmin) continue;
         if (TMath::Abs(ptrk->GetD(fX,fY))>fRmax) continue;

         if (TMath::Abs(ntrk->GetD(fX,fY))<fDNmin)
         if (TMath::Abs(ptrk->GetD(fX,fY))<fDNmin) continue;


         AliITStrackV2 nt(*ntrk), pt(*ptrk), *pnt=&nt, *ppt=&pt;

         Double_t dca=PropagateToDCA(pnt,ppt);
         if (dca > fDCAmax) continue;

         AliV0vertex vertex(*pnt,*ppt);
         if (vertex.GetChi2() > fChi2max) continue;
	 
         /*  Think of something better here ! 
         nt.PropagateToVertex(); if (TMath::Abs(nt.GetZ())<0.04) continue;
         pt.PropagateToVertex(); if (TMath::Abs(pt.GetZ())<0.04) continue;
	 */

         Double_t x,y,z; vertex.GetXYZ(x,y,z); 
         Double_t r2=x*x + y*y; 
         if (r2 > fRmax*fRmax) continue;
         if (r2 < fRmin*fRmin) continue;

         Double_t px,py,pz; vertex.GetPxPyPz(px,py,pz);
         Double_t p2=px*px+py*py+pz*pz;
         Double_t cost=((x-fX)*px + (y-fY)*py + (z-fZ)*pz)/
               TMath::Sqrt(p2*((x-fX)*(x-fX) + (y-fY)*(y-fY) + (z-fZ)*(z-fZ)));

        //if (cost < (5*fCPAmax-0.9-TMath::Sqrt(r2)*(fCPAmax-1))/4.1) continue;
         if (cost < fCPAmax) continue;
	 vertex.SetDcaDaughters(dca);
         //vertex.ChangeMassHypothesis(); //default is Lambda0 

         event->AddV0(&vertex);

         nvtx++;
      }
   }

   Info("Tracks2V0vertices","Number of reconstructed V0 vertices: %d",nvtx);

   negtrks.Delete();
   postrks.Delete();

   return 0;
}

Int_t AliV0vertexer::Tracks2V0vertices(TTree *tTree, TTree *vTree) {
  //--------------------------------------------------------------------
  //This function reconstructs V0 vertices
  //--------------------------------------------------------------------
  Warning("Tracks2V0vertices(TTree*,TTree*)",
	  "Will be removed soon !  Use Tracks2V0vertices(AliESD*) instead");

   TBranch *branch=tTree->GetBranch("tracks");
   if (!branch) {
      Error("Tracks2V0vertices","Can't get the branch !");
      return 1;
   }
   Int_t nentr=(Int_t)tTree->GetEntries();

   TObjArray negtrks(nentr/2);
   TObjArray postrks(nentr/2);

   Int_t nneg=0, npos=0, nvtx=0;

   Int_t i;
   for (i=0; i<nentr; i++) {
       AliITStrackV2 *iotrack=new AliITStrackV2;
       branch->SetAddress(&iotrack);
       tTree->GetEvent(i);

       if (!iotrack->PropagateTo(3.,0.0023,65.19)) {
	 delete iotrack;
	 continue; 
       }
       if (!iotrack->PropagateTo(2.5,0.,0.)) {
	 delete iotrack;
	 continue;
       }

       if (iotrack->Get1Pt() > 0.) {nneg++; negtrks.AddLast(iotrack);}
       else {npos++; postrks.AddLast(iotrack);}
   }   


   AliV0vertex *ioVertex=0;
   branch=vTree->GetBranch("vertices");
   if (!branch) vTree->Branch("vertices","AliV0vertex",&ioVertex,32000,3);
   else branch->SetAddress(&ioVertex);


   for (i=0; i<nneg; i++) {
      //if (i%10==0) cerr<<nneg-i<<'\r';
      AliITStrackV2 *ntrk=(AliITStrackV2 *)negtrks.UncheckedAt(i);

      if (TMath::Abs(ntrk->GetD(fX,fY))<fDPmin) continue;
      if (TMath::Abs(ntrk->GetD(fX,fY))>fRmax) continue;

      for (Int_t k=0; k<npos; k++) {
         AliITStrackV2 *ptrk=(AliITStrackV2 *)postrks.UncheckedAt(k);

         if (TMath::Abs(ptrk->GetD(fX,fY))<fDPmin) continue;
         if (TMath::Abs(ptrk->GetD(fX,fY))>fRmax) continue;

         if (TMath::Abs(ntrk->GetD(fX,fY))<fDNmin)
         if (TMath::Abs(ptrk->GetD(fX,fY))<fDNmin) continue;


         AliITStrackV2 nt(*ntrk), pt(*ptrk), *pnt=&nt, *ppt=&pt;

         Double_t dca=PropagateToDCA(pnt,ppt);
         if (dca > fDCAmax) continue;

         AliV0vertex vertex(*pnt,*ppt);
         if (vertex.GetChi2() > fChi2max) continue;
	 
         /*  Think of something better here ! 
         nt.PropagateToVertex(); if (TMath::Abs(nt.GetZ())<0.04) continue;
         pt.PropagateToVertex(); if (TMath::Abs(pt.GetZ())<0.04) continue;
	 */

         Double_t x,y,z; vertex.GetXYZ(x,y,z); 
         Double_t r2=x*x + y*y; 
         if (r2 > fRmax*fRmax) continue;
         if (r2 < fRmin*fRmin) continue;

         Double_t px,py,pz; vertex.GetPxPyPz(px,py,pz);
         Double_t p2=px*px+py*py+pz*pz;
         Double_t cost=((x-fX)*px + (y-fY)*py + (z-fZ)*pz)/
               TMath::Sqrt(p2*((x-fX)*(x-fX) + (y-fY)*(y-fY) + (z-fZ)*(z-fZ)));

        //if (cost < (5*fCPAmax-0.9-TMath::Sqrt(r2)*(fCPAmax-1))/4.1) continue;
         if (cost < fCPAmax) continue;
	 vertex.SetDcaDaughters(dca);
         //vertex.ChangeMassHypothesis(); //default is Lambda0 

         ioVertex=&vertex; vTree->Fill();

         nvtx++;
      }
   }

   Info("Tracks2V0vertices","Number of reconstructed V0 vertices: %d",nvtx);

   negtrks.Delete();
   postrks.Delete();

   return 0;
}

Double_t 
AliV0vertexer::PropagateToDCA(AliITStrackV2 *n, AliITStrackV2 *p) const {
  //--------------------------------------------------------------------
  // This function returns the DCA between two tracks
  // The tracks will be moved to the point of DCA ! 
  //--------------------------------------------------------------------
  return n->PropagateToDCA(p);
}
















