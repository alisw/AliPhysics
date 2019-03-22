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
#include "AliCascadeVertexer.h"

ClassImp(AliCascadeVertexer)

//A set of loose cuts
Double_t 
AliCascadeVertexer::fgChi2max=33.;   //maximal allowed chi2
Double_t 
AliCascadeVertexer::fgDV0min=0.01;   //min V0 impact parameter
Double_t 
AliCascadeVertexer::fgMassWin=0.008; //"window" around the Lambda mass
Double_t 
AliCascadeVertexer::fgDBachMin=0.01; //min bachelor impact parameter
Double_t 
AliCascadeVertexer::fgDCAmax=2.0;    //max DCA between the V0 and the track
Double_t 
AliCascadeVertexer::fgCPAmin=0.98; //min cosine of the cascade pointing angle
Double_t 
AliCascadeVertexer::fgRmin=0.2;      //min radius of the fiducial volume
Double_t 
AliCascadeVertexer::fgRmax=100.;     //max radius of the fiducial volume


Int_t AliCascadeVertexer::V0sTracks2CascadeVertices(AliESDEvent *event) {
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
  int trk[nentr], ntr=0;
  for (i=0; i<nentr; i++) {
    AliESDtrack *esdtr=event->GetTrack(i);
    ULong64_t status=esdtr->GetStatus();
    if (status&AliESDtrack::kITSpureSA) continue;
    if ((status&AliESDtrack::kITSrefit)==0)
      if ((status&AliESDtrack::kTPCrefit)==0) continue;
    
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
      if (bidx==v0.GetIndex(0)) continue; //Bo:  consistency 0 for neg
      AliESDtrack *btrk=event->GetTrack(bidx);
      if (btrk->GetSign()>0) continue;  // bachelor's charge
      
      AliESDv0 *pv0=&v0;
      AliExternalTrackParam bt(*btrk), *pbt=&bt;
      
      Double_t dca=100;
      if(!fUseImprovedFinding){
        dca = PropagateToDCA(pv0,pbt,b);
      }else{
        dca = PropagateToDCACurvedBachelor(pv0,pbt,b);
      }
      
      if (dca > fDCAmax) continue;
      
      AliESDcascade cascade(*pv0,*pbt,bidx);//constucts a cascade candidate
      //PH        if (cascade.GetChi2Xi() > fChi2max) continue;
      
      if(fUseImprovedFinding) cascade.RefitCascade(pbt); //imp pos + cov mat
      
      Double_t x,y,z; cascade.GetXYZcascade(x,y,z); // Bo: bug correction
      Double_t r2=x*x + y*y;
      if (r2 > fRmax2) continue;   // condition on fiducial zone
      if (r2 < fRmin2) continue;
      
      Double_t pxV0,pyV0,pzV0;
      pv0->GetPxPyPz(pxV0,pyV0,pzV0);
      if (x*pxV0+y*pyV0+z*pzV0 < 0) continue; //causality
      
      Double_t x1,y1,z1; pv0->GetXYZ(x1,y1,z1);
      if (r2 > (x1*x1+y1*y1)) continue;
      
      if (cascade.GetCascadeCosineOfPointingAngle(xPrimaryVertex,yPrimaryVertex,zPrimaryVertex) <fCPAmin) continue; //condition on the cascade pointing angle
      
      cascade.SetDcaXiDaughters(dca);
      event->AddCascade(&cascade);
      v->SetUsedByCascade(kTRUE);
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
      if (bidx==v0.GetIndex(1)) continue; //Bo:  consistency 1 for pos
      AliESDtrack *btrk=event->GetTrack(bidx);
      if (btrk->GetSign()<0) continue;  // bachelor's charge
      
      AliESDv0 *pv0=&v0;
      AliExternalTrackParam bt(*btrk), *pbt=&bt;
      
      Double_t dca=100;
      if(!fUseImprovedFinding){
        dca = PropagateToDCA(pv0,pbt,b);
      }else{
        dca = PropagateToDCACurvedBachelor(pv0,pbt,b);
      }
      if (dca > fDCAmax) continue;
      
      AliESDcascade cascade(*pv0,*pbt,bidx); //constucts a cascade candidate
      //PH         if (cascade.GetChi2Xi() > fChi2max) continue;
      
      if(fUseImprovedFinding) cascade.RefitCascade(pbt); //imp pos + cov mat
      
      Double_t x,y,z; cascade.GetXYZcascade(x,y,z); // Bo: bug correction
      Double_t r2=x*x + y*y;
      if (r2 > fRmax2) continue;   // condition on fiducial zone
      if (r2 < fRmin2) continue;
      
      Double_t pxV0,pyV0,pzV0;
      pv0->GetPxPyPz(pxV0,pyV0,pzV0);
      if (x*pxV0+y*pyV0+z*pzV0 < 0) continue; //causality
      
      Double_t x1,y1,z1; pv0->GetXYZ(x1,y1,z1);
      if (r2 > (x1*x1+y1*y1)) continue;
      
      if (cascade.GetCascadeCosineOfPointingAngle(xPrimaryVertex,yPrimaryVertex,zPrimaryVertex) < fCPAmin) continue; //condition on the cascade pointing angle
      
      cascade.SetDcaXiDaughters(dca);
      event->AddCascade(&cascade);
      v->SetUsedByCascade(kTRUE);
      ncasc++;
      
    } // end loop tracks
  } // end loop V0s
  
  Info("V0sTracks2CascadeVertices","Number of reconstructed cascades: %d",ncasc);
  
  return 0;
}


Double_t AliCascadeVertexer::Det(Double_t a00, Double_t a01, Double_t a10, Double_t a11) const {
  //--------------------------------------------------------------------
  // This function calculates locally a 2x2 determinant
  //--------------------------------------------------------------------
  return a00*a11 - a01*a10;
}

Double_t AliCascadeVertexer::Det(Double_t a00,Double_t a01,Double_t a02,
                                 Double_t a10,Double_t a11,Double_t a12,
                                 Double_t a20,Double_t a21,Double_t a22) const {
  //--------------------------------------------------------------------
  // This function calculates locally a 3x3 determinant
  //--------------------------------------------------------------------
  return  a00*Det(a11,a12,a21,a22)-a01*Det(a10,a12,a20,a22)+a02*Det(a10,a11,a20,a21);
}




Double_t AliCascadeVertexer::PropagateToDCA(AliESDv0 *v, AliExternalTrackParam *t, Double_t b) {
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
  if (dca > fDCAmax) return 1.e+33;
  
  //points of the DCA
  Double_t t1 = Det(x2-x1,y2-y1,z2-z1,px2,py2,pz2,ax,ay,az)/
  Det(px1,py1,pz1,px2,py2,pz2,ax,ay,az);
  
  x1 += px1*t1; y1 += py1*t1; //z1 += pz1*t1;
  
  if (x1*x1+y1*y1 > fRmaxMargin2) return 1.e+33;
  
  //propagate track to the points of DCA
  
  x1=x1*cs1 + y1*sn1;
  
  if (!t->PropagateTo(x1,b)) {
    AliError("Propagation failed");
    //    AliErrorF("Propagation failed for X=%f | V0: %f %f %f",x1,x2,y2,z2);
    //    t->Print();
    //
    return 1.e+33;
  }
  
  return dca;
}

Double_t AliCascadeVertexer::PropagateToDCACurvedBachelor(AliESDv0 *v, AliExternalTrackParam *t, Double_t b)
{
  //--------------------------------------------------------------------
  // This function returns the DCA between the V0 and the track
  // assumes that bachelor track is not straight
  // algorithm based on AliExternalTrackParam::GetDCA with zero curvature track
  //--------------------------------------------------------------------
  // Double_t alpha=t->GetAlpha(), cs1=TMath::Cos(alpha), sn1=TMath::Sin(alpha);
  Double_t r[3]; t->GetXYZ(r);
  // Double_t x1=r[0], y1=r[1], z1=r[2];
  Double_t p[3]; t->GetPxPyPz(p);
  // Double_t px1=p[0], py1=p[1], pz1=p[2];

  Double_t x2,y2,z2;     // position and momentum of V0
  Double_t px2,py2,pz2;
  
  v->GetXYZ(x2,y2,z2);
  v->GetPxPyPz(px2,py2,pz2);
  
  Double_t dca = 1e+33;
  Double_t dy2=1e-10;
  Double_t dz2=1e-10;
  Double_t dx2=1e-10;
  
  //Create dummy V0 track
  //V0 properties to get started
  Double_t xyz[3], pxpypz[3], cv[21];
  for(Int_t ii=0;ii<21;ii++) cv[ii]=0.0; //something small
  
  v->GetXYZ(xyz[0],xyz[1],xyz[2]);
  v->GetPxPyPz( pxpypz[0],pxpypz[1],pxpypz[2] );
  
  //Mockup track for V0 trajectory (no covariance)
  //AliExternalTrackParam *hV0Traj = new AliExternalTrackParam(xyz,pxpypz,cv,+1);
  AliExternalTrackParam lV0TrajObject(xyz,pxpypz,cv,+1), *hV0Traj = &lV0TrajObject;
  hV0Traj->ResetCovariance(1); //won't use
  
  //Re-acquire helix parameters for bachelor (necessary!)
  Double_t p1[8]; t->GetHelixParameters(p1,b);
  p1[6]=TMath::Sin(p1[2]);
  p1[7]=TMath::Cos(p1[2]);
  
  Double_t p2[8]; hV0Traj->GetHelixParameters(p2,0.0); //p2[4]=0 -> no curvature (fine, predicted in Evaluate)
  p2[6]=TMath::Sin(p2[2]); p2[7]=TMath::Cos(p2[2]);
  
  Double_t r1[3],g1[3],gg1[3]; Double_t t1=0.;
  Evaluate(p1,t1,r1,g1,gg1);
  Double_t r2[3],g2[3],gg2[3]; Double_t t2=0.;
  Evaluate(p2,t2,r2,g2,gg2);
  
  Double_t dx=r2[0]-r1[0], dy=r2[1]-r1[1], dz=r2[2]-r1[2];
  Double_t dm=dx*dx/dx2 + dy*dy/dy2 + dz*dz/dz2;
  
  Int_t max=27; //standard in AliExternalTrackParam::GetDCA, good performance
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
        if ((gt1*gt1+gt2*gt2) > 1.e-4/dy2/dy2){
          AliDebug(1," stopped at not a stationary point !");
        }
        Double_t lmb=h11+h22; lmb=lmb-TMath::Sqrt(lmb*lmb-4*det);
        if (lmb < 0.){
          AliDebug(1," stopped at not a minimum !");
        }
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
        AliDebug(1," overshoot !"); break;
      }
    }
    dm=dd;
    t1+=dt1;
    t2+=dt2;
  }
  if (max<=0){
    AliDebug(1," too many iterations !");
  }
  Double_t cs=TMath::Cos(t->GetAlpha());
  Double_t sn=TMath::Sin(t->GetAlpha());
  Double_t xthis=r1[0]*cs + r1[1]*sn;
  
  //Propagate bachelor to the point of DCA
  if (!t->PropagateTo(xthis,b)) {
    AliDebug(1," propagation failed !");
    return 1e+33;
  }
  
  //V0 distance to bachelor: the desired distance
  Double_t rBachDCAPt[3]; t->GetXYZ(rBachDCAPt);
  dca = v->GetD(rBachDCAPt[0],rBachDCAPt[1],rBachDCAPt[2]);
  return dca;
}

//________________________________________________________________________
void AliCascadeVertexer::Evaluate(const Double_t *h, Double_t t,
                                  Double_t r[3],  //radius vector
                                  Double_t g[3],  //first defivatives
                                  Double_t gg[3]) //second derivatives
{
  //--------------------------------------------------------------------
  // Calculate position of a point on a track and some derivatives
  //--------------------------------------------------------------------
  Double_t phase=h[4]*t+h[2];
  Double_t sn=TMath::Sin(phase), cs=TMath::Cos(phase);
  
  r[0] = h[5];
  r[1] = h[0];
  if (TMath::Abs(h[4])>kAlmost0) {
    r[0] += (sn - h[6])/h[4];
    r[1] -= (cs - h[7])/h[4];
  } else {
    r[0] += t*cs;
    r[1] -= -t*sn;
  }
  r[2] = h[1] + h[3]*t;
  
  g[0] = cs; g[1]=sn; g[2]=h[3];
  
  gg[0]=-h[4]*sn; gg[1]=h[4]*cs; gg[2]=0.;
}

//________________________________________________
void AliCascadeVertexer::SetUseImprovedFinding(const Bool_t lInput)
{
  fUseImprovedFinding = lInput;
}
