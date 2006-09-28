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

//modified by R. Vernet 30/6/2006 : daughter label
//modified by R. Vernet  3/7/2006 : causality


#include <TObjArray.h>
#include <TTree.h>

#include "AliESD.h"
#include "AliESDv0.h"
#include "AliESDcascade.h"
#include "AliCascadeVertexer.h"

ClassImp(AliCascadeVertexer)

Int_t AliCascadeVertexer::V0sTracks2CascadeVertices(AliESD *event) {
  //--------------------------------------------------------------------
  // This function reconstructs cascade vertices
  //      Adapted to the ESD by I.Belikov (Jouri.Belikov@cern.ch)
  //--------------------------------------------------------------------
   Double_t b=event->GetMagneticField();
   Int_t nV0=(Int_t)event->GetNumberOfV0s();

   //stores relevant V0s in an array
   TObjArray vtcs(nV0);
   Int_t i;
   for (i=0; i<nV0; i++) {
       AliESDv0 *v=event->GetV0(i);
       if (v->GetOnFlyStatus()) continue;
       if (v->GetD(fX,fY,fZ)<fDV0min) continue;
       vtcs.AddLast(v);
   }
   nV0=vtcs.GetEntriesFast();

   // stores relevant tracks in another array
   Int_t nentr=(Int_t)event->GetNumberOfTracks();
   TArrayI trk(nentr); Int_t ntr=0;
   for (i=0; i<nentr; i++) {
       AliESDtrack *esdtr=event->GetTrack(i);
       UInt_t status=esdtr->GetStatus();
       UInt_t flags=AliESDtrack::kITSin|AliESDtrack::kTPCin|
                    AliESDtrack::kTPCpid|AliESDtrack::kESDpid;

       if ((status&AliESDtrack::kITSrefit)==0)
          if (flags!=status) continue;

       if (TMath::Abs(esdtr->GetD(fX,fY,b))<fDBachMin) continue;

       trk[ntr++]=i;
   }   

   Double_t massLambda=1.11568;
   Int_t ncasc=0;

   // Looking for the cascades...

   for (i=0; i<nV0; i++) { //loop on V0s

      AliESDv0 *v=(AliESDv0*)vtcs.UncheckedAt(i);
      v->ChangeMassHypothesis(kLambda0); // the v0 must be Lambda 
      if (TMath::Abs(v->GetEffMass()-massLambda)>fMassWin) continue; 

      for (Int_t j=0; j<ntr; j++) {//loop on tracks
	 Int_t bidx=trk[j];
 	 if (bidx==v->GetNindex()) continue; //bachelor and v0's negative tracks must be different
	 AliESDtrack *btrk=event->GetTrack(bidx);
         if (btrk->GetSign()>0) continue;  // bachelor's charge 
          
    	 AliESDv0 v0(*v), *pv0=&v0;
         AliExternalTrackParam bt(*btrk), *pbt=&bt;

         Double_t dca=PropagateToDCA(pv0,pbt,b);
         if (dca > fDCAmax) continue;

         AliESDcascade cascade(*pv0,*pbt,bidx);//constucts a cascade candidate
         if (cascade.GetChi2Xi() > fChi2max) continue;

	 Double_t x,y,z; cascade.GetXYZ(x,y,z); 
         Double_t r2=x*x + y*y; 
         if (r2 > fRmax*fRmax) continue;   // condition on fiducial zone
         if (r2 < fRmin*fRmin) continue;

	 Double_t pxV0,pyV0,pzV0;
	 pv0->GetPxPyPz(pxV0,pyV0,pzV0);
	 if (x*pxV0+y*pyV0+z*pzV0 < 0) continue; //causality

         Double_t x1,y1,z1; pv0->GetXYZ(x1,y1,z1);
         if (r2 > (x1*x1+y1*y1)) continue;

  	 if (cascade.GetCascadeCosineOfPointingAngle(fX,fY,fZ) <fCPAmax) continue; //condition on the cascade pointing angle 
	 
	 event->AddCascade(&cascade);
         ncasc++;
      } // end loop tracks
   } // end loop V0s

   // Looking for the anti-cascades...

   for (i=0; i<nV0; i++) { //loop on V0s
      AliESDv0 *v=(AliESDv0*)vtcs.UncheckedAt(i);
      v->ChangeMassHypothesis(kLambda0Bar); //the v0 must be anti-Lambda 
      if (TMath::Abs(v->GetEffMass()-massLambda)>fMassWin) continue; 

      for (Int_t j=0; j<ntr; j++) {//loop on tracks
	 Int_t bidx=trk[j];
 	 if (bidx==v->GetPindex()) continue; //bachelor and v0's positive tracks must be different
	 AliESDtrack *btrk=event->GetTrack(bidx);
         if (btrk->GetSign()<0) continue;  // bachelor's charge 
          
	 AliESDv0 v0(*v), *pv0=&v0;
         AliESDtrack bt(*btrk), *pbt=&bt;

         Double_t dca=PropagateToDCA(pv0,pbt,b);
         if (dca > fDCAmax) continue;

         AliESDcascade cascade(*pv0,*pbt,bidx); //constucts a cascade candidate
         if (cascade.GetChi2Xi() > fChi2max) continue;

	 Double_t x,y,z; cascade.GetXYZ(x,y,z); 
         Double_t r2=x*x + y*y; 
         if (r2 > fRmax*fRmax) continue;   // condition on fiducial zone
         if (r2 < fRmin*fRmin) continue;

	 Double_t pxV0,pyV0,pzV0;
	 pv0->GetPxPyPz(pxV0,pyV0,pzV0);
	 if (x*pxV0+y*pyV0+z*pzV0 < 0) continue; //causality

         Double_t x1,y1,z1; pv0->GetXYZ(x1,y1,z1);
         if (r2 > (x1*x1+y1*y1)) continue;
         if (z*z > z1*z1) continue;

	 if (cascade.GetCascadeCosineOfPointingAngle(fX,fY,fZ) < fCPAmax) continue; //condition on the cascade pointing angle 
	 event->AddCascade(&cascade);
         ncasc++;

      } // end loop tracks
   } // end loop V0s

Info("V0sTracks2CascadeVertices","Number of reconstructed cascades: %d",ncasc);

   return 0;
}


Double_t det(Double_t a00, Double_t a01, Double_t a10, Double_t a11){
  // determinant 2x2
  return a00*a11 - a01*a10;
}

Double_t det (Double_t a00,Double_t a01,Double_t a02,
                     Double_t a10,Double_t a11,Double_t a12,
                     Double_t a20,Double_t a21,Double_t a22) {
  // determinant 3x3
  return 
  a00*det(a11,a12,a21,a22)-a01*det(a10,a12,a20,a22)+a02*det(a10,a11,a20,a21);
}




Double_t AliCascadeVertexer::
PropagateToDCA(AliESDv0 *v, AliExternalTrackParam *t, Double_t b) {
  //--------------------------------------------------------------------
  // This function returns the DCA between the V0 and the track
  //--------------------------------------------------------------------
  Double_t alpha=t->GetAlpha(), cs1=TMath::Cos(alpha), sn1=TMath::Sin(alpha);
  Double_t r[3]; t->GetXYZ(r);
  Double_t x1=r[0], y1=r[1], z1=r[2];
  Double_t p[3]; t->GetPxPyPz(p);
  Double_t px1=p[0], py1=p[1], pz1=p[2];
  
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
  if (!t->PropagateTo(x1,b)) {
    Error("PropagateToDCA","Propagation failed !");
    return 1.e+33;
  }  

  return dca;
}











