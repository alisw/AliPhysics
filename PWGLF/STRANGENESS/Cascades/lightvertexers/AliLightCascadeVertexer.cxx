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
//modified by I. Belikov 24/11/2006 : static setter for the default cuts

#include "AliESDEvent.h"
#include "AliESDcascade.h"
#include "AliLightCascadeVertexer.h"

ClassImp(AliLightCascadeVertexer)

//A set of loose cuts
Double_t 
  AliLightCascadeVertexer::fgChi2max=33.;   //maximal allowed chi2 
Double_t 
  AliLightCascadeVertexer::fgDV0min=0.01;   //min V0 impact parameter
Double_t 
  AliLightCascadeVertexer::fgMassWin=0.008; //"window" around the Lambda mass
Double_t 
  AliLightCascadeVertexer::fgDBachMin=0.01; //min bachelor impact parameter
Double_t 
  AliLightCascadeVertexer::fgDCAmax=2.0;    //max DCA between the V0 and the track 
Double_t 
  AliLightCascadeVertexer::fgCPAmin=0.98; //min cosine of the cascade pointing angle
Double_t 
  AliLightCascadeVertexer::fgRmin=0.2;      //min radius of the fiducial volume
Double_t 
  AliLightCascadeVertexer::fgRmax=100.;     //max radius of the fiducial volume

Double_t AliLightCascadeVertexer::fgMaxEta=0.8;        //max |eta|
Double_t AliLightCascadeVertexer::fgMinClusters=70;   //min clusters (>=)
Bool_t AliLightCascadeVertexer::fgSwitchCharges=kFALSE;   //min clusters (>=)

Int_t AliLightCascadeVertexer::V0sTracks2CascadeVertices(AliESDEvent *event) {
  //--------------------------------------------------------------------
  // This function reconstructs cascade vertices
  //      Adapted to the ESD by I.Belikov (Jouri.Belikov@cern.ch)
  //--------------------------------------------------------------------
   const AliESDVertex *vtxT3D=event->GetPrimaryVertex();

   Double_t xPrimaryVertex=vtxT3D->GetX();
   Double_t yPrimaryVertex=vtxT3D->GetY();
   Double_t zPrimaryVertex=vtxT3D->GetZ();

   Double_t b=event->GetMagneticField();
   Int_t nV0=(Int_t)event->GetNumberOfV0s();

   //stores relevant V0s in an array
   TObjArray vtcs(nV0);
   Int_t i;
   for (i=0; i<nV0; i++) {
       AliESDv0 *v=event->GetV0(i);
       if (v->GetOnFlyStatus()) continue;
       if (v->GetD(xPrimaryVertex,yPrimaryVertex,zPrimaryVertex)<fDV0min) continue;
       vtcs.AddLast(v);
   }
   nV0=vtcs.GetEntriesFast();

   // stores relevant tracks in another array
   Int_t nentr=(Int_t)event->GetNumberOfTracks();
   TArrayI trk(nentr); Int_t ntr=0;
   for (i=0; i<nentr; i++) {
       AliESDtrack *esdtr=event->GetTrack(i);
       ULong_t status=esdtr->GetStatus();

       if ((status&AliESDtrack::kITSrefit)==0)
	  if ((status&AliESDtrack::kTPCrefit)==0) continue;
       
       //Track pre-selection: clusters
       if (esdtr->GetTPCNcls() < fMinClusters ) continue;

       if (TMath::Abs(esdtr->GetD(xPrimaryVertex,yPrimaryVertex,b))<fDBachMin) continue;

       trk[ntr++]=i;
   }   

   Double_t massLambda=1.11568;
   Int_t ncasc=0;

   // Looking for the cascades...

   for (i=0; i<nV0; i++) { //loop on V0s

      AliESDv0 *v=(AliESDv0*)vtcs.UncheckedAt(i);
      AliESDv0 v0(*v);
      v0.ChangeMassHypothesis(kLambda0); // the v0 must be Lambda 
      if (TMath::Abs(v0.GetEffMass()-massLambda)>fMassWin) continue; 

      for (Int_t j=0; j<ntr; j++) {//loop on tracks
	 Int_t bidx=trk[j];
 	 //Bo:   if (bidx==v->GetNindex()) continue; //bachelor and v0's negative tracks must be different
          if (!fSwitchCharges && bidx==v0.GetIndex(0)) continue; //Bo:  consistency 0 for neg
          if ( fSwitchCharges && bidx==v0.GetIndex(1)) continue; //Bo:  consistency 0 for neg
          
          AliESDtrack *btrk=event->GetTrack(bidx);
          
         if (!fSwitchCharges && btrk->GetSign()>0) continue;  // bachelor's charge
         if ( fSwitchCharges && btrk->GetSign()<0) continue;  // bachelor's charge
          
    	 AliESDv0 *pv0=&v0;
         AliExternalTrackParam bt(*btrk), *pbt=&bt;

         Double_t dca=PropagateToDCA(pv0,pbt,b);
         if (dca > fDCAmax) continue;
          
          //eta cut - test
            if (TMath::Abs(pbt->Eta())>fMaxEta) continue;

         AliESDcascade cascade(*pv0,*pbt,bidx);//constucts a cascade candidate
	 //PH        if (cascade.GetChi2Xi() > fChi2max) continue;

	 Double_t x,y,z; cascade.GetXYZcascade(x,y,z); // Bo: bug correction
         Double_t r2=x*x + y*y; 
         if (r2 > fRmax*fRmax) continue;   // condition on fiducial zone
         if (r2 < fRmin*fRmin) continue;

	 Double_t pxV0,pyV0,pzV0;
	 pv0->GetPxPyPz(pxV0,pyV0,pzV0);
	 if (x*pxV0+y*pyV0+z*pzV0 < 0) continue; //causality

         Double_t x1,y1,z1; pv0->GetXYZ(x1,y1,z1);
         if (r2 > (x1*x1+y1*y1)) continue;

  	 if (cascade.GetCascadeCosineOfPointingAngle(xPrimaryVertex,yPrimaryVertex,zPrimaryVertex) <fCPAmin) continue; //condition on the cascade pointing angle 
	 
         cascade.SetDcaXiDaughters(dca);
	 event->AddCascade(&cascade);
         ncasc++;
      } // end loop tracks
   } // end loop V0s

   // Looking for the anti-cascades...

   for (i=0; i<nV0; i++) { //loop on V0s
      AliESDv0 *v=(AliESDv0*)vtcs.UncheckedAt(i);
      AliESDv0 v0(*v);
      v0.ChangeMassHypothesis(kLambda0Bar); //the v0 must be anti-Lambda 
      if (TMath::Abs(v0.GetEffMass()-massLambda)>fMassWin) continue; 

      for (Int_t j=0; j<ntr; j++) {//loop on tracks
	 Int_t bidx=trk[j];
 	 //Bo:   if (bidx==v->GetPindex()) continue; //bachelor and v0's positive tracks must be different
         if (!fSwitchCharges && bidx==v0.GetIndex(1)) continue; //Bo:  consistency 1 for pos
         if ( fSwitchCharges && bidx==v0.GetIndex(0)) continue; //Bo:  consistency 1 for pos
          
          AliESDtrack *btrk=event->GetTrack(bidx);
          
         if (!fSwitchCharges && btrk->GetSign()<0) continue;  // bachelor's charge
         if ( fSwitchCharges && btrk->GetSign()>0) continue;  // bachelor's charge
          
	 AliESDv0 *pv0=&v0;
         AliExternalTrackParam bt(*btrk), *pbt=&bt;

         Double_t dca=PropagateToDCA(pv0,pbt,b);
         if (dca > fDCAmax) continue;

          //eta cut - test
          if (TMath::Abs(pbt->Eta())>fMaxEta) continue;
          
         AliESDcascade cascade(*pv0,*pbt,bidx); //constucts a cascade candidate
	 //PH         if (cascade.GetChi2Xi() > fChi2max) continue;

	 Double_t x,y,z; cascade.GetXYZcascade(x,y,z); // Bo: bug correction
         Double_t r2=x*x + y*y; 
         if (r2 > fRmax*fRmax) continue;   // condition on fiducial zone
         if (r2 < fRmin*fRmin) continue;

	 Double_t pxV0,pyV0,pzV0;
	 pv0->GetPxPyPz(pxV0,pyV0,pzV0);
	 if (x*pxV0+y*pyV0+z*pzV0 < 0) continue; //causality

         Double_t x1,y1,z1; pv0->GetXYZ(x1,y1,z1);
         if (r2 > (x1*x1+y1*y1)) continue;

	 if (cascade.GetCascadeCosineOfPointingAngle(xPrimaryVertex,yPrimaryVertex,zPrimaryVertex) < fCPAmin) continue; //condition on the cascade pointing angle 

         cascade.SetDcaXiDaughters(dca);
	 event->AddCascade(&cascade);
         ncasc++;

      } // end loop tracks
   } // end loop V0s

Info("V0sTracks2CascadeVertices","Number of reconstructed cascades: %d",ncasc);

   return 0;
}


Double_t AliLightCascadeVertexer::Det(Double_t a00, Double_t a01, Double_t a10, Double_t a11) const {
  //--------------------------------------------------------------------
  // This function calculates locally a 2x2 determinant
  //--------------------------------------------------------------------
  return a00*a11 - a01*a10;
}

Double_t AliLightCascadeVertexer::Det(Double_t a00,Double_t a01,Double_t a02,
				 Double_t a10,Double_t a11,Double_t a12,
				 Double_t a20,Double_t a21,Double_t a22) const {
  //--------------------------------------------------------------------
  // This function calculates locally a 3x3 determinant
  //--------------------------------------------------------------------
  return  a00*Det(a11,a12,a21,a22)-a01*Det(a10,a12,a20,a22)+a02*Det(a10,a11,a20,a21);
}




Double_t AliLightCascadeVertexer::PropagateToDCA(AliESDv0 *v, AliExternalTrackParam *t, Double_t b) {
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
   
  Double_t dd= Det(x2-x1,y2-y1,z2-z1,px1,py1,pz1,px2,py2,pz2);
  Double_t ax= Det(py1,pz1,py2,pz2);
  Double_t ay=-Det(px1,pz1,px2,pz2);
  Double_t az= Det(px1,py1,px2,py2);

  Double_t dca=TMath::Abs(dd)/TMath::Sqrt(ax*ax + ay*ay + az*az);

//points of the DCA
  Double_t t1 = Det(x2-x1,y2-y1,z2-z1,px2,py2,pz2,ax,ay,az)/
                Det(px1,py1,pz1,px2,py2,pz2,ax,ay,az);
  
  x1 += px1*t1; y1 += py1*t1; //z1 += pz1*t1;
  

  //propagate track to the points of DCA

  x1=x1*cs1 + y1*sn1;
  if (!t->PropagateTo(x1,b)) {
    Error("PropagateToDCA","Propagation failed !");
    return 1.e+33;
  }  

  return dca;
}











