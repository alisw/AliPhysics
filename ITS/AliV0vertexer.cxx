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
//
//     Origin: Iouri Belikov, IReS, Strasbourg, Jouri.Belikov@cern.ch
//-------------------------------------------------------------------------
#include <TFile.h>
#include <TTree.h>
#include <TObjArray.h>
#include <Riostream.h>

#include "AliV0vertex.h"
#include "AliV0vertexer.h"
#include "AliITStrackV2.h"

ClassImp(AliV0vertexer)

Int_t AliV0vertexer::Tracks2V0vertices(const TFile *inp, TFile *out) {
  //--------------------------------------------------------------------
  //This function reconstructs V0 vertices
  //--------------------------------------------------------------------
   TFile *in=(TFile*)inp;
   TDirectory *savedir=gDirectory; 

   if (!in->IsOpen()) {
     cerr<<"AliV0vertexer::Tracks2V0vertices(): ";
     cerr<<"file with ITS tracks has not been open !\n";
     return 1;
   }

   if (!out->IsOpen()) {
     cerr<<"AliV0vertexer::Tracks2V0vertices(): ";
     cerr<<"file for V0 vertices has not been open !\n";
     return 2;
   }

   in->cd();

   TTree *trkTree=(TTree*)in->Get("TreeT_ITS_0");
   TBranch *branch=trkTree->GetBranch("tracks");
   Int_t nentr=(Int_t)trkTree->GetEntries();

   TObjArray negtrks(nentr/2);
   TObjArray postrks(nentr/2);

   Int_t nneg=0, npos=0, nvtx=0;

   Int_t i;
   for (i=0; i<nentr; i++) {
       AliITStrackV2 *iotrack=new AliITStrackV2;
       branch->SetAddress(&iotrack);
       trkTree->GetEvent(i);

       iotrack->PropagateTo(3.,0.0023,65.19); iotrack->PropagateTo(2.5,0.,0.);

       if (iotrack->Get1Pt() > 0.) {nneg++; negtrks.AddLast(iotrack);}
       else {npos++; postrks.AddLast(iotrack);}
   }   


   out->cd();
   TTree vtxTree("TreeV","Tree with V0 vertices");
   AliV0vertex *ioVertex=0;
   vtxTree.Branch("vertices","AliV0vertex",&ioVertex,32000,0);


   for (i=0; i<nneg; i++) {
      if (i%10==0) cerr<<nneg-i<<'\r';
      AliITStrackV2 *ntrk=(AliITStrackV2 *)negtrks.UncheckedAt(i);

      if (TMath::Abs(ntrk->GetD())<fDPmin) continue;
      if (TMath::Abs(ntrk->GetD())>fRmax) continue;

      for (Int_t k=0; k<npos; k++) {
         AliITStrackV2 *ptrk=(AliITStrackV2 *)postrks.UncheckedAt(k);

         if (TMath::Abs(ptrk->GetD())<fDPmin) continue;
         if (TMath::Abs(ptrk->GetD())>fRmax) continue;

         if (TMath::Abs(ntrk->GetD())<fDNmin)
         if (TMath::Abs(ptrk->GetD())<fDNmin) continue;


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
         Double_t cost=(x*px+y*py+z*pz)/TMath::Sqrt(p2*(r2+z*z));

         if (cost < (5*fCPAmax-0.9-TMath::Sqrt(r2)*(fCPAmax-1))/4.1) continue;
         //if (cost < fCPAmax) continue;

         //vertex.ChangeMassHypothesis(); //default is Lambda0 

         ioVertex=&vertex; vtxTree.Fill();

         nvtx++;
      }
   }

   cerr<<"Number of reconstructed V0 vertices: "<<nvtx<<endl;

   vtxTree.Write();

   negtrks.Delete();
   postrks.Delete();

   savedir->cd();
   
   delete trkTree;

   return 0;
}


static void External2Helix(const AliITStrackV2 *t, Double_t helix[6]) { 
  //--------------------------------------------------------------------
  // External track parameters -> helix parameters 
  //--------------------------------------------------------------------
  Double_t alpha,x,cs,sn;
  t->GetExternalParameters(x,helix); alpha=t->GetAlpha();

  cs=TMath::Cos(alpha); sn=TMath::Sin(alpha);
  helix[5]=x*cs - helix[0]*sn;            // x0
  helix[0]=x*sn + helix[0]*cs;            // y0
//helix[1]=                               // z0
  helix[2]=TMath::ASin(helix[2]) + alpha; // phi0
//helix[3]=                               // tgl
  helix[4]=helix[4]/t->GetConvConst();    // C
}

static void Evaluate(const Double_t *h, Double_t t,
                     Double_t r[3],  //radius vector
                     Double_t g[3],  //first defivatives
                     Double_t gg[3]) //second derivatives
{
  //--------------------------------------------------------------------
  // Calculate position of a point on a track and some derivatives
  //--------------------------------------------------------------------
  Double_t phase=h[4]*t+h[2];
  Double_t sn=TMath::Sin(phase), cs=TMath::Cos(phase);

  r[0] = h[5] + (sn - h[6])/h[4];
  r[1] = h[0] - (cs - h[7])/h[4];  
  r[2] = h[1] + h[3]*t;

  g[0] = cs; g[1]=sn; g[2]=h[3];
  
  gg[0]=-h[4]*sn; gg[1]=h[4]*cs; gg[2]=0.;
}

Double_t AliV0vertexer::PropagateToDCA(AliITStrackV2 *n, AliITStrackV2 *p) {
  //--------------------------------------------------------------------
  // This function returns the DCA between two tracks
  // The tracks will be moved to the point of DCA ! 
  //--------------------------------------------------------------------
  Double_t dy2=n->GetSigmaY2() + p->GetSigmaY2();
  Double_t dz2=n->GetSigmaZ2() + p->GetSigmaZ2();
  Double_t dx2=dy2; 

  //dx2=dy2=dz2=1.;

  Double_t p1[8]; External2Helix(n,p1);
  p1[6]=TMath::Sin(p1[2]); p1[7]=TMath::Cos(p1[2]);
  Double_t p2[8]; External2Helix(p,p2);
  p2[6]=TMath::Sin(p2[2]); p2[7]=TMath::Cos(p2[2]);


  Double_t r1[3],g1[3],gg1[3]; Double_t t1=0.;
  Evaluate(p1,t1,r1,g1,gg1);
  Double_t r2[3],g2[3],gg2[3]; Double_t t2=0.;
  Evaluate(p2,t2,r2,g2,gg2);

  Double_t dx=r2[0]-r1[0], dy=r2[1]-r1[1], dz=r2[2]-r1[2];
  Double_t dm=dx*dx/dx2 + dy*dy/dy2 + dz*dz/dz2;

  Int_t max=27;
  while (max--) {
     Double_t gt1=-(dx*g1[0]/dx2 + dy*g1[1]/dy2 + dz*g1[2]/dz2);
     Double_t gt2=+(dx*g2[0]/dx2 + dy*g2[1]/dy2 + dz*g2[2]/dz2);
     Double_t h11=(g1[0]*g1[0] - dx*gg1[0])/dx2 + 
                  (g1[1]*g1[1] - dy*gg1[1])/dy2 +
                  (g1[2]*g1[2] - dz*gg1[2])/dz2;
     Double_t h22=(g2[0]*g2[0] + dx*gg2[0])/dx2 + 
                  (g2[1]*g2[1] + dy*gg2[1])/dy2 +
                  (g2[2]*g2[2] + dz*gg2[2])/dz2;
     Double_t h12=-(g1[0]*g2[0]/dx2 + g1[1]*g2[1]/dy2 + g1[2]*g2[2]/dz2);

     Double_t det=h11*h22-h12*h12;

     Double_t dt1,dt2;
     if (TMath::Abs(det)<1.e-33) {
        //(quasi)singular Hessian
        dt1=-gt1; dt2=-gt2;
     } else {
        dt1=-(gt1*h22 - gt2*h12)/det; 
        dt2=-(h11*gt2 - h12*gt1)/det;
     }

     if ((dt1*gt1+dt2*gt2)>0) {dt1=-dt1; dt2=-dt2;}

     //check delta(phase1) ?
     //check delta(phase2) ?

     if (TMath::Abs(dt1)/(TMath::Abs(t1)+1.e-3) < 1.e-4)
     if (TMath::Abs(dt2)/(TMath::Abs(t2)+1.e-3) < 1.e-4) {
        if ((gt1*gt1+gt2*gt2) > 1.e-4/dy2/dy2) 
	   cerr<<"AliV0vertexer::PropagateToDCA:"
                 " stopped at not a stationary point !\n";
        Double_t lmb=h11+h22; lmb=lmb-TMath::Sqrt(lmb*lmb-4*det);
        if (lmb < 0.) 
	   cerr<<"AliV0vertexer::PropagateToDCA:"
                 " stopped at not a minimum !\n";
        break;
     }

     Double_t dd=dm;
     for (Int_t div=1 ; ; div*=2) {
        Evaluate(p1,t1+dt1,r1,g1,gg1);
        Evaluate(p2,t2+dt2,r2,g2,gg2);
        dx=r2[0]-r1[0]; dy=r2[1]-r1[1]; dz=r2[2]-r1[2];
        dd=dx*dx/dx2 + dy*dy/dy2 + dz*dz/dz2;
	if (dd<dm) break;
        dt1*=0.5; dt2*=0.5;
        if (div>512) {
           cerr<<"AliV0vertexer::PropagateToDCA: overshoot !\n"; break;
        }   
     }
     dm=dd;

     t1+=dt1;
     t2+=dt2;

  }

  if (max<=0) cerr<<"AliV0vertexer::PropagateToDCA: too many iterations !\n";  
  
  //propagate tracks to the points of DCA
  Double_t cs=TMath::Cos(n->GetAlpha());
  Double_t sn=TMath::Sin(n->GetAlpha());
  Double_t x=r1[0]*cs + r1[1]*sn;
  if (!n->PropagateTo(x,0.,0.)) {
    //cerr<<"AliV0vertexer::PropagateToDCA: propagation failed !\n";
    return 1.e+33;
  }  

  cs=TMath::Cos(p->GetAlpha());
  sn=TMath::Sin(p->GetAlpha());
  x=r2[0]*cs + r2[1]*sn;
  if (!p->PropagateTo(x,0.,0.)) {
    //cerr<<"AliV0vertexer::PropagateToDCA: propagation failed !\n";
    return 1.e+33;
  }  

  return TMath::Sqrt(dm*TMath::Sqrt(dy2*dz2));
  //return TMath::Sqrt(dx*dx + dy*dy + dz*dz);
}
















