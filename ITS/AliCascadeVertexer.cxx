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
//               Implementation of the cascade vertexer class
//          Reads V0s and tracks, writes out cascade vertices
//                     Fills the ESD with the cascades 
//    Origin: Christian Kuhn, IReS, Strasbourg, christian.kuhn@ires.in2p3.fr
//-------------------------------------------------------------------------
#include <TObjArray.h>
#include <TTree.h>

#include "AliESD.h"
#include "AliESDv0.h"
#include "AliCascadeVertex.h"
#include "AliCascadeVertexer.h"
#include "AliITStrackV2.h"
#include "AliV0vertex.h"

ClassImp(AliCascadeVertexer)

Int_t AliCascadeVertexer::V0sTracks2CascadeVertices(AliESD *event) {
  //--------------------------------------------------------------------
  // This function reconstructs cascade vertices
  //      Adapted to the ESD by I.Belikov (Jouri.Belikov@cern.ch)
  //--------------------------------------------------------------------

   Int_t nV0=(Int_t)event->GetNumberOfV0s();
   TObjArray vtcs(nV0);
   Int_t i;
   for (i=0; i<nV0; i++) {
       const AliESDv0 *esdV0=event->GetV0(i);
       vtcs.AddLast(new AliV0vertex(*esdV0));
   }


   Int_t ntr=(Int_t)event->GetNumberOfTracks();
   TObjArray trks(ntr);
   for (i=0; i<ntr; i++) {
       AliESDtrack *esdtr=event->GetTrack(i);
       UInt_t status=esdtr->GetStatus();
       UInt_t flags=AliESDtrack::kITSin|AliESDtrack::kTPCin;

       if ((status&AliESDtrack::kITSrefit)==0)
          if ((status&flags)!=status) continue;

       AliITStrackV2 *iotrack=new AliITStrackV2(*esdtr);
       iotrack->SetLabel(i);  // now it is the index in array of ESD tracks
       if ((status&AliESDtrack::kITSrefit)==0)   //correction for the beam pipe
          if (!iotrack->PropagateTo(3.,0.0023,65.19)) continue; 
       if (!iotrack->PropagateTo(2.5,0.,0.)) continue;
       trks.AddLast(iotrack);
   }   
   ntr=trks.GetEntriesFast();

   Double_t massLambda=1.11568;
   Int_t ncasc=0;

   // Looking for the cascades...
   for (i=0; i<nV0; i++) {
      AliV0vertex *v=(AliV0vertex*)vtcs.UncheckedAt(i);
      v->ChangeMassHypothesis(kLambda0); // the v0 must be Lambda 
      if (TMath::Abs(v->GetEffMass()-massLambda)>fMassWin) continue; 
      if (v->GetD(0,0,0)<fDV0min) continue;
      for (Int_t j=0; j<ntr; j++) {
	 AliITStrackV2 *b=(AliITStrackV2*)trks.UncheckedAt(j);

         if (TMath::Abs(b->GetD())<fDBachMin) continue;
         if (b->Get1Pt()<0.) continue;  // bachelor's charge 
          
	 AliV0vertex v0(*v), *pv0=&v0;
         AliITStrackV2 bt(*b), *pbt=&bt;

         Double_t dca=PropagateToDCA(pv0,pbt);
         if (dca > fDCAmax) continue;

         AliCascadeVertex cascade(*pv0,*pbt);
         if (cascade.GetChi2() > fChi2max) continue;

	 Double_t x,y,z; cascade.GetXYZ(x,y,z); 
         Double_t r2=x*x + y*y; 
         if (r2 > fRmax*fRmax) continue;   // condition on fiducial zone
         if (r2 < fRmin*fRmin) continue;

         {
         Double_t x1,y1,z1; pv0->GetXYZ(x1,y1,z1);
         if (r2 > (x1*x1+y1*y1)) continue;
         if (z*z > z1*z1) continue;
         }

	 Double_t px,py,pz; cascade.GetPxPyPz(px,py,pz);
         Double_t p2=px*px+py*py+pz*pz;
         Double_t cost=(x*px+y*py+z*pz)/TMath::Sqrt(p2*(r2+z*z));

         if (cost<fCPAmax) continue; //condition on the cascade pointing angle 
         //cascade.ChangeMassHypothesis(); //default is Xi

         event->AddCascade(&cascade);

         ncasc++;

      }
   }

   // Looking for the anti-cascades...
   for (i=0; i<nV0; i++) {
      AliV0vertex *v=(AliV0vertex*)vtcs.UncheckedAt(i);
      v->ChangeMassHypothesis(kLambda0Bar); //the v0 must be anti-Lambda 
      if (TMath::Abs(v->GetEffMass()-massLambda)>fMassWin) continue; 
      if (v->GetD(0,0,0)<fDV0min) continue;
      for (Int_t j=0; j<ntr; j++) {
	 AliITStrackV2 *b=(AliITStrackV2*)trks.UncheckedAt(j);

         if (TMath::Abs(b->GetD())<fDBachMin) continue;
         if (b->Get1Pt()>0.) continue;  // bachelor's charge 
          
	 AliV0vertex v0(*v), *pv0=&v0;
         AliITStrackV2 bt(*b), *pbt=&bt;

         Double_t dca=PropagateToDCA(pv0,pbt);
         if (dca > fDCAmax) continue;

         AliCascadeVertex cascade(*pv0,*pbt);
         if (cascade.GetChi2() > fChi2max) continue;

	 Double_t x,y,z; cascade.GetXYZ(x,y,z); 
         Double_t r2=x*x + y*y; 
         if (r2 > fRmax*fRmax) continue;   // condition on fiducial zone
         if (r2 < fRmin*fRmin) continue;

         {
         Double_t x1,y1,z1; pv0->GetXYZ(x1,y1,z1);
         if (r2 > (x1*x1+y1*y1)) continue;
         if (z*z > z1*z1) continue;
         }

	 Double_t px,py,pz; cascade.GetPxPyPz(px,py,pz);
         Double_t p2=px*px+py*py+pz*pz;
         Double_t cost=(x*px+y*py+z*pz)/TMath::Sqrt(p2*(r2+z*z));

         if (cost<fCPAmax) continue; //condition on the cascade pointing angle 
         //cascade.ChangeMassHypothesis(); //default is Xi

         event->AddCascade(&cascade);

         ncasc++;

      }
   }

Info("V0sTracks2CascadeVertices","Number of reconstructed cascades: %d",ncasc);

   trks.Delete();
   vtcs.Delete();

   return 0;
}

Int_t AliCascadeVertexer::
V0sTracks2CascadeVertices(TTree *vTree,TTree *tTree, TTree *xTree) {
  //--------------------------------------------------------------------
  // This function reconstructs cascade vertices
  //--------------------------------------------------------------------
  Warning("V0sTracks2CascadeVertices(TTree*,TTree*,TTree*)",
  "Will be removed soon !  Use V0sTracks2CascadeVertices(AliESD*) instead");

   TBranch *branch=vTree->GetBranch("vertices");
   if (!branch) {
      Error("V0sTracks2CascadeVertices","Can't get the V0 branch !");
      return 1;
   }   
   Int_t nentrV0=(Int_t)vTree->GetEntries();
   
   TObjArray vtxV0(nentrV0);

   // fill TObjArray vtxV0 with vertices

   Int_t i, nV0=0;
   for (i=0; i<nentrV0; i++) {

       AliV0vertex *ioVertex=new AliV0vertex;
       branch->SetAddress(&ioVertex);
       vTree->GetEvent(i);
       nV0++; 
       vtxV0.AddLast(ioVertex);
       
   }

   branch=tTree->GetBranch("tracks");
   if (!branch) {
      Error("V0sTracks2CascadeVertices","Can't get the track branch !");
      return 2;
   }
   Int_t nentr=(Int_t)tTree->GetEntries();

   TObjArray trks(nentr);

   // fill TObjArray trks with tracks

   Int_t ntrack=0;

   for (i=0; i<nentr; i++) {

       AliITStrackV2 *iotrack=new AliITStrackV2;
       branch->SetAddress(&iotrack);
       tTree->GetEvent(i);

       if (!iotrack->PropagateTo(3.,0.0023,65.19)) continue; 
       if (!iotrack->PropagateTo(2.5,0.,0.)) continue;

       ntrack++; trks.AddLast(iotrack);
       
   }   

  AliCascadeVertex *ioCascade=0;
  branch=xTree->GetBranch("cascades");
  if (!branch) xTree->Branch("cascades","AliCascadeVertex",&ioCascade,32000,3);
  else branch->SetAddress(&ioCascade); 

   // loop on all vertices

   Double_t massLambda=1.11568;
   Int_t ncasc=0;

   for (i=0; i<nV0; i++) {

       AliV0vertex *lV0ver=(AliV0vertex *)vtxV0.UncheckedAt(i);

       lV0ver->ChangeMassHypothesis(kLambda0); //I.B.

       if (lV0ver->GetEffMass()<massLambda-fMassWin ||       // condition of the V0 mass window (cut fMassWin)
           lV0ver->GetEffMass()>massLambda+fMassWin) continue; 

       if (lV0ver->GetD(0,0,0)<fDV0min) continue;          // condition of minimum impact parameter of the V0 (cut fDV0min) 
                                                          // here why not cuting on pointing angle ???

   // for each vertex in the good mass range, loop on all tracks (= bachelor candidates)

       for (Int_t k=0; k<ntrack; k++) {

	  AliITStrackV2 *bachtrk=(AliITStrackV2 *)trks.UncheckedAt(k);

          if (TMath::Abs(bachtrk->GetD())<fDBachMin) continue;        // eliminate to small impact parameters

          if (lV0ver->GetPdgCode()==kLambda0 && bachtrk->Get1Pt()<0.) continue;     // condition on V0 label 
          if (lV0ver->GetPdgCode()==kLambda0Bar && bachtrk->Get1Pt()>0.) continue;  // + good sign for bachelor
          
	  AliV0vertex lV0(*lV0ver), *pV0=&lV0;
          AliITStrackV2 bt(*bachtrk), *pbt=&bt;

   // calculation of the distance of closest approach between the V0 and the bachelor

          Double_t dca=PropagateToDCA(pV0,pbt);
          if (dca > fDCAmax) continue;                         // cut on dca

   // construction of a cascade object

          AliCascadeVertex cascade(*pV0,*pbt);
          if (cascade.GetChi2() > fChi2max) continue;

   // get cascade decay position (V0, bachelor)
          
	 Double_t x,y,z; cascade.GetXYZ(x,y,z); 
         Double_t r2=x*x + y*y; 
         if (r2 > fRmax*fRmax) continue;   // condition on fiducial zone
         if (r2 < fRmin*fRmin) continue;

         {
   //I.B.
         Double_t x1,y1,z1; lV0ver->GetXYZ(x1,y1,z1);
         if (r2 > (x1*x1+y1*y1)) continue;
         if (z*z > z1*z1) continue;
         }

   // get cascade momentum
         
	 Double_t px,py,pz; cascade.GetPxPyPz(px,py,pz);
         Double_t p2=px*px+py*py+pz*pz;
         Double_t cost=(x*px+y*py+z*pz)/TMath::Sqrt(p2*(r2+z*z));

         if (cost<fCPAmax) continue; // condition on the cascade pointing angle 

         //cascade.ChangeMassHypothesis(); //default is Xi


         ioCascade=&cascade; xTree->Fill();

         ncasc++;

      }
   }

Info("V0sTracks2CascadeVertices","Number of reconstructed cascades: %d",ncasc);

   trks.Delete();
   vtxV0.Delete();

   return 0;
}

inline Double_t det(Double_t a00, Double_t a01, Double_t a10, Double_t a11){
  // determinant 2x2
  return a00*a11 - a01*a10;
}

inline Double_t det (Double_t a00,Double_t a01,Double_t a02,
                     Double_t a10,Double_t a11,Double_t a12,
                     Double_t a20,Double_t a21,Double_t a22) {
  // determinant 3x3
  return 
  a00*det(a11,a12,a21,a22)-a01*det(a10,a12,a20,a22)+a02*det(a10,a11,a20,a21);
}




Double_t 
AliCascadeVertexer::PropagateToDCA(AliV0vertex *v, AliITStrackV2 *t) {
  //--------------------------------------------------------------------
  // This function returns the DCA between the V0 and the track
  //--------------------------------------------------------------------

  Double_t phi, x, par[5];
  Double_t alpha, cs1, sn1;

  t->GetExternalParameters(x,par); alpha=t->GetAlpha();
  phi=TMath::ASin(par[2]) + alpha;  
  Double_t px1=TMath::Cos(phi), py1=TMath::Sin(phi), pz1=par[3];

 
  cs1=TMath::Cos(alpha); sn1=TMath::Sin(alpha);
  Double_t x1=x*cs1 - par[0]*sn1;
  Double_t y1=x*sn1 + par[0]*cs1;
  Double_t z1=par[1];

  
  Double_t x2,y2,z2;     // position and momentum of V0
  Double_t px2,py2,pz2;
  
  v->GetXYZ(x2,y2,z2);
  v->GetPxPyPz(px2,py2,pz2);
 
// calculation dca
   
  Double_t dd= det(x2-x1,y2-y1,z2-z1,px1,py1,pz1,px2,py2,pz2);
  Double_t ax= det(py1,pz1,py2,pz2);
  Double_t ay=-det(px1,pz1,px2,pz2);
  Double_t az= det(px1,py1,px2,py2);

  Double_t dca=TMath::Abs(dd)/TMath::Sqrt(ax*ax + ay*ay + az*az);

//points of the DCA
  Double_t t1 = det(x2-x1,y2-y1,z2-z1,px2,py2,pz2,ax,ay,az)/
                det(px1,py1,pz1,px2,py2,pz2,ax,ay,az);
  
  x1 += px1*t1; y1 += py1*t1; //z1 += pz1*t1;
  

  //propagate track to the points of DCA

  x1=x1*cs1 + y1*sn1;
  if (!t->PropagateTo(x1,0.,0.)) {
    Error("PropagateToDCA","Propagation failed !");
    return 1.e+33;
  }  

  return dca;
}











