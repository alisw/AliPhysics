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
//
//    Origin: Christian Kuhn, IReS, Strasbourg, christian.kuhn@ires.in2p3.fr
//-------------------------------------------------------------------------
#include <TFile.h>
#include <TTree.h>
#include <TObjArray.h>
#include <Riostream.h>

#include "AliCascadeVertex.h"
#include "AliCascadeVertexer.h"
#include "AliV0vertex.h"
#include "AliITStrackV2.h"

ClassImp(AliCascadeVertexer)

Int_t 
AliCascadeVertexer::V0sTracks2CascadeVertices(const TFile *inp, TFile *out) {

  //--------------------------------------------------------------------
  //This function reconstructs cascade vertices
  //--------------------------------------------------------------------
  TFile *in=(TFile*)inp;
  TDirectory *savedir=gDirectory;

   // Tree for vertices(V0's)

   TTree *vtxTree=(TTree*)gDirectory->Get("TreeV0");
   TBranch *branch=vtxTree->GetBranch("vertices");
   Int_t nentrV0=(Int_t)vtxTree->GetEntries();
   
   TObjArray vtxV0(nentrV0);

   // fill TObjArray vtxV0 with vertices

   Int_t i, nV0=0;
   for (i=0; i<nentrV0; i++) {

       AliV0vertex *ioVertex=new AliV0vertex;
       branch->SetAddress(&ioVertex);
       vtxTree->GetEvent(i);
       nV0++; 
       vtxV0.AddLast(ioVertex);
       
   }


   in->cd();

  // Tree for tracks

   TTree *trkTree=(TTree*)in->Get("TreeT_ITS_0");
   branch=trkTree->GetBranch("tracks");
   Int_t nentr=(Int_t)trkTree->GetEntries();

   TObjArray trks(nentr);

   // fill TObjArray trks with tracks

   Int_t ntrack=0;

   for (i=0; i<nentr; i++) {

       AliITStrackV2 *iotrack=new AliITStrackV2;
       branch->SetAddress(&iotrack);
       trkTree->GetEvent(i);

       iotrack->PropagateTo(3.,0.0023,65.19); iotrack->PropagateTo(2.5,0.,0.);

       ntrack++; trks.AddLast(iotrack);
       
   }   

   // create Tree for cascades 

   out->cd();

   TTree cascTree("TreeCasc","Tree with cascades");
   AliCascadeVertex *ioCascade=0;
   cascTree.Branch("cascades","AliCascadeVertex",&ioCascade,32000,0);

   // loop on all vertices

   Double_t massLambda=1.11568;
   Int_t ncasc=0;

   for (i=0; i<nV0; i++) {

       AliV0vertex *V0ver=(AliV0vertex *)vtxV0.UncheckedAt(i);

       V0ver->ChangeMassHypothesis(kLambda0); //I.B.

       if (V0ver->GetEffMass()<massLambda-fMassWin ||       // condition of the V0 mass window (cut fMassWin)
           V0ver->GetEffMass()>massLambda+fMassWin) continue; 

       if (V0ver->GetD(0,0,0)<fDV0min) continue;          // condition of minimum impact parameter of the V0 (cut fDV0min) 
                                                          // here why not cuting on pointing angle ???

   // for each vertex in the good mass range, loop on all tracks (= bachelor candidates)

       for (Int_t k=0; k<ntrack; k++) {

	  AliITStrackV2 *bachtrk=(AliITStrackV2 *)trks.UncheckedAt(k);

          if (TMath::Abs(bachtrk->GetD())<fDBachMin) continue;        // eliminate to small impact parameters

          if (V0ver->GetPdgCode()==kLambda0 && bachtrk->Get1Pt()<0.) continue;     // condition on V0 label 
          if (V0ver->GetPdgCode()==kLambda0Bar && bachtrk->Get1Pt()>0.) continue;  // + good sign for bachelor
          
	  AliV0vertex V0(*V0ver), *pV0=&V0;
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
         Double_t x1,y1,z1; V0ver->GetXYZ(x1,y1,z1);
         if (r2 > (x1*x1+y1*y1)) continue;
         if (z*z > z1*z1) continue;
         }

   // get cascade momentum
         
	 Double_t px,py,pz; cascade.GetPxPyPz(px,py,pz);
         Double_t p2=px*px+py*py+pz*pz;
         Double_t cost=(x*px+y*py+z*pz)/TMath::Sqrt(p2*(r2+z*z));

         if (cost<fCPAmax) continue; // condition on the cascade pointing angle 

         //cascade.ChangeMassHypothesis(); //default is Xi


         ioCascade=&cascade; cascTree.Fill();

         ncasc++;

      }
   }

   cerr<<"Number of reconstructed cascades: "<<ncasc<<endl;

   cascTree.Write();

   trks.Delete();
   vtxV0.Delete();

   savedir->cd();

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
    cerr<<"AliV0vertexer::PropagateToDCA: propagation failed !\n";
    return 1.e+33;
  }  

  return dca;
}











