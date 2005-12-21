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

/* $Id$ */

//-------------------------------------------------------------------------
//                Implementation of the AliKalmanTrack class
//   that is the base for AliTPCtrack, AliITStrackV2 and AliTRDtrack
//        Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//-------------------------------------------------------------------------

#include "AliKalmanTrack.h"
#include "TGeoManager.h"

ClassImp(AliKalmanTrack)

const AliMagF *AliKalmanTrack::fgkFieldMap=0;
Double_t AliKalmanTrack::fgConvConst=0.;

//_______________________________________________________________________
AliKalmanTrack::AliKalmanTrack():
  fLab(-3141593),
  fFakeRatio(0),
  fChi2(0),
  fMass(AliPID::ParticleMass(AliPID::kPion)),
  fN(0),
  fLocalConvConst(0),
  fStartTimeIntegral(kFALSE),
  fIntegratedLength(0)
{
  //
  // Default constructor
  //
    if (fgkFieldMap==0) {
      AliFatal("The magnetic field has not been set!");
    }

    for(Int_t i=0; i<AliPID::kSPECIES; i++) fIntegratedTime[i] = 0;
}

//_______________________________________________________________________
AliKalmanTrack::AliKalmanTrack(const AliKalmanTrack &t):
  TObject(t),
  fLab(t.fLab),
  fFakeRatio(t.fFakeRatio),
  fChi2(t.fChi2),
  fMass(t.fMass),
  fN(t.fN),
  fLocalConvConst(t.fLocalConvConst),
  fStartTimeIntegral(t.fStartTimeIntegral),
  fIntegratedLength(t.fIntegratedLength)
{
  //
  // Copy constructor
  //
  if (fgkFieldMap==0) {
    AliFatal("The magnetic field has not been set!");
  }
  
  for (Int_t i=0; i<AliPID::kSPECIES; i++)
      fIntegratedTime[i] = t.fIntegratedTime[i];
}

//_______________________________________________________________________
void AliKalmanTrack::StartTimeIntegral() 
{
  // Sylwester Radomski, GSI
  // S.Radomski@gsi.de
  //
  // Start time integration
  // To be called at Vertex by ITS tracker
  //
  
  //if (fStartTimeIntegral) 
  //  AliWarning("Reseting Recorded Time.");

  fStartTimeIntegral = kTRUE;
  for(Int_t i=0; i<AliPID::kSPECIES; i++) fIntegratedTime[i] = 0;  
  fIntegratedLength = 0;
}
//_______________________________________________________________________
void AliKalmanTrack:: AddTimeStep(Double_t length) 
{
  // 
  // Add step to integrated time
  // this method should be called by a sublasses at the end
  // of the PropagateTo function or by a tracker
  // each time step is made.
  //
  // If integration not started function does nothing
  //
  // Formula
  // dt = dl * sqrt(p^2 + m^2) / p
  // p = pT * (1 + tg^2 (lambda) )
  //
  // pt = 1/external parameter [4]
  // tg lambda = external parameter [3]
  //
  //
  // Sylwester Radomski, GSI
  // S.Radomski@gsi.de
  // 
  
  static const Double_t kcc = 2.99792458e-2;

  if (!fStartTimeIntegral) return;
  
  fIntegratedLength += length;

  Double_t xr, param[5];
  Double_t pt, tgl;
  
  GetExternalParameters(xr, param);
  pt =  1/param[4] ;
  tgl = param[3];

  Double_t p = TMath::Abs(pt * TMath::Sqrt(1+tgl*tgl));

  if (length > 100) return;

  for (Int_t i=0; i<AliPID::kSPECIES; i++) {
    
    Double_t mass = AliPID::ParticleMass(i);
    Double_t correction = TMath::Sqrt( pt*pt * (1 + tgl*tgl) + mass * mass ) / p;
    Double_t time = length * correction / kcc;

    fIntegratedTime[i] += time;
  }
}

//_______________________________________________________________________

Double_t AliKalmanTrack::GetIntegratedTime(Int_t pdg) const 
{
  // Sylwester Radomski, GSI
  // S.Radomski@gsi.de
  //
  // Return integrated time hypothesis for a given particle
  // type assumption.
  //
  // Input parameter:
  // pdg - Pdg code of a particle type
  //


  if (!fStartTimeIntegral) {
    AliWarning("Time integration not started");
    return 0.;
  }

  for (Int_t i=0; i<AliPID::kSPECIES; i++)
    if (AliPID::ParticleCode(i) == TMath::Abs(pdg)) return fIntegratedTime[i];

  AliWarning(Form("Particle type [%d] not found", pdg));
  return 0;
}

void AliKalmanTrack::GetIntegratedTimes(Double_t *times) const {
  for (Int_t i=0; i<AliPID::kSPECIES; i++) times[i]=fIntegratedTime[i];
}

void AliKalmanTrack::SetIntegratedTimes(const Double_t *times) {
  for (Int_t i=0; i<AliPID::kSPECIES; i++) fIntegratedTime[i]=times[i];
}

//_______________________________________________________________________

void AliKalmanTrack::PrintTime() const
{
  // Sylwester Radomski, GSI
  // S.Radomski@gsi.de
  //
  // For testing
  // Prints time for all hypothesis
  //

  for (Int_t i=0; i<AliPID::kSPECIES; i++)
    printf("%d: %.2f  ", AliPID::ParticleCode(i), fIntegratedTime[i]);
  printf("\n");  
}

void AliKalmanTrack::External2Helix(Double_t helix[6]) const { 
  //--------------------------------------------------------------------
  // External track parameters -> helix parameters 
  //--------------------------------------------------------------------
  Double_t alpha,x,cs,sn;
  GetExternalParameters(x,helix); alpha=GetAlpha();

  cs=TMath::Cos(alpha); sn=TMath::Sin(alpha);
  helix[5]=x*cs - helix[0]*sn;            // x0
  helix[0]=x*sn + helix[0]*cs;            // y0
//helix[1]=                               // z0
  helix[2]=TMath::ASin(helix[2]) + alpha; // phi0
//helix[3]=                               // tgl
  helix[4]=helix[4]/GetLocalConvConst();  // C
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

Double_t AliKalmanTrack::
GetDCA(const AliKalmanTrack *p, Double_t &xthis, Double_t &xp) const {
  //------------------------------------------------------------
  // Returns the (weighed !) distance of closest approach between 
  // this track and the track passed as the argument.
  // Other returned values:
  //   xthis, xt - coordinates of tracks' reference planes at the DCA 
  //-----------------------------------------------------------
  Double_t dy2=GetSigmaY2() + p->GetSigmaY2();
  Double_t dz2=GetSigmaZ2() + p->GetSigmaZ2();
  Double_t dx2=dy2; 

  //dx2=dy2=dz2=1.;

  Double_t p1[8]; External2Helix(p1);
  p1[6]=TMath::Sin(p1[2]); p1[7]=TMath::Cos(p1[2]);
  Double_t p2[8]; p->External2Helix(p2);
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
	  AliWarning(" stopped at not a stationary point !");
        Double_t lmb=h11+h22; lmb=lmb-TMath::Sqrt(lmb*lmb-4*det);
        if (lmb < 0.) 
	  AliWarning(" stopped at not a minimum !");
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
           AliWarning(" overshoot !"); break;
        }   
     }
     dm=dd;

     t1+=dt1;
     t2+=dt2;

  }

  if (max<=0) AliWarning(" too many iterations !");

  Double_t cs=TMath::Cos(GetAlpha());
  Double_t sn=TMath::Sin(GetAlpha());
  xthis=r1[0]*cs + r1[1]*sn;

  cs=TMath::Cos(p->GetAlpha());
  sn=TMath::Sin(p->GetAlpha());
  xp=r2[0]*cs + r2[1]*sn;

  return TMath::Sqrt(dm*TMath::Sqrt(dy2*dz2));
}

Double_t AliKalmanTrack::
PropagateToDCA(AliKalmanTrack *p, Double_t d, Double_t x0) {
  //--------------------------------------------------------------
  // Propagates this track and the argument track to the position of the
  // distance of closest approach. 
  // Returns the (weighed !) distance of closest approach.
  //--------------------------------------------------------------
  Double_t xthis,xp;
  Double_t dca=GetDCA(p,xthis,xp);

  if (!PropagateTo(xthis,d,x0)) {
    //AliWarning(" propagation failed !");
    return 1e+33;
  }  

  if (!p->PropagateTo(xp,d,x0)) {
    //AliWarning(" propagation failed !";
    return 1e+33;
  }  

  return dca;
}





Double_t AliKalmanTrack::MeanMaterialBudget(Double_t *start, Double_t *end, Double_t *mparam)
{
  //
  // calculate mean material budget and material properties beween point start and end
  // mparam - returns parameters used for dEdx and multiple scatering
  //
  // mparam[0] - density mean 
  // mparam[1] - rad length
  // mparam[2] - A mean
  // mparam[3] - Z mean
  // mparam[4] - length
  // mparam[5] - Z/A mean
  // mparam[6] - number of boundary crosses
  //
    mparam[0]=0; mparam[1]=1; mparam[2] =0; mparam[3] =0, mparam[4]=0, mparam[5]=0; mparam[6]=0; 
  //
  Double_t bparam[6], lparam[6];          // bparam - total param - lparam - local parameters
  for (Int_t i=0;i<6;i++) bparam[i]=0;    // 

  if (!gGeoManager) {
    printf("ERROR: no TGeo\n");
    return 0.;
  }
  //
  Double_t length;
  Double_t dir[3];
  length = TMath::Sqrt((end[0]-start[0])*(end[0]-start[0])+
		       (end[1]-start[1])*(end[1]-start[1])+
		       (end[2]-start[2])*(end[2]-start[2]));
  mparam[4]=length;
  if (length<TGeoShape::Tolerance()) return 0.0;
  Double_t invlen = 1./length;
  dir[0] = (end[0]-start[0])*invlen;
  dir[1] = (end[1]-start[1])*invlen;
  dir[2] = (end[2]-start[2])*invlen;
  // Initialize start point and direction
  TGeoNode *currentnode = 0;
  TGeoNode *startnode = gGeoManager->InitTrack(start, dir);
  //  printf("%s length=%f\n",gGeoManager->GetPath(),length);
  if (!startnode) {
    printf("ERROR: start point out of geometry\n");
    return 0.0;
  }
  TGeoMaterial *material = startnode->GetVolume()->GetMedium()->GetMaterial();
  lparam[0] = material->GetDensity();
  lparam[1]   = material->GetRadLen();
  lparam[2]   = material->GetA();
  lparam[3]   = material->GetZ();
  lparam[5]   = lparam[3]/lparam[2];
  if (material->IsMixture()) {
    lparam[1]*=lparam[0];  // different normalization in the modeler for mixture
    TGeoMixture * mixture = (TGeoMixture*)material;
    lparam[5] =0;
    Double_t sum =0;
    for (Int_t iel=0;iel<mixture->GetNelements();iel++){
      sum  += mixture->GetWmixt()[iel];
      lparam[5]+= mixture->GetZmixt()[iel]*mixture->GetWmixt()[iel]/mixture->GetAmixt()[iel];
    }
    lparam[5]/=sum;
  }
  gGeoManager->FindNextBoundary(length);
  Double_t snext = gGeoManager->GetStep();
  Double_t step = 0.0;
  // If no boundary within proposed length, return current density
  if (snext>=length) {
    for (Int_t ip=0;ip<5;ip++) mparam[ip] = lparam[ip];  
    return lparam[0];
  }
  // Try to cross the boundary and see what is next
  while (length>TGeoShape::Tolerance()) {
    mparam[6]+=1.;
    currentnode = gGeoManager->Step();
    step += snext+1.E-6;
    bparam[1]    += snext*lparam[1];
    bparam[2]    += snext*lparam[2];
    bparam[3]    += snext*lparam[3];
    bparam[5]    += snext*lparam[5];    
    bparam[0]    += snext*lparam[0];

    if (snext>=length) break;
    if (!currentnode) break;
    //    printf("%s snext=%f  density=%f bparam[0]=%f\n", gGeoManager->GetPath(),snext,density,bparam[0]);
    if (!gGeoManager->IsEntering()) {
      gGeoManager->SetStep(1.E-3);
      currentnode = gGeoManager->Step();
      if (!gGeoManager->IsEntering() || !currentnode) {
	//	printf("ERROR: cannot cross boundary\n"); 
	mparam[0] = bparam[0]/step;
	mparam[1] = bparam[1]/step;
	mparam[2] = bparam[2]/step;
	mparam[3] = bparam[3]/step;	
	mparam[5] = bparam[5]/step;	
	mparam[4] = step;
	mparam[0] = 0.;             // if crash of navigation take mean density 0 
	mparam[1] = 1000000;        // and infinite rad length
        return bparam[0]/step;
      }
      step += 1.E-3;
      snext += 1.E-3;
      bparam[0] += lparam[0]*1.E-3;
      bparam[1]    += lparam[1]*1.E-3;
      bparam[2]    += lparam[2]*1.E-3;
      bparam[3]    += lparam[3]*1.E-3;
      bparam[5]    += lparam[5]*1.E-3;
    }
    length -= snext;
    material = currentnode->GetVolume()->GetMedium()->GetMaterial();
    lparam[0] = material->GetDensity();
    lparam[1]  = material->GetRadLen();
    lparam[2]  = material->GetA();
    lparam[3]  = material->GetZ();
    lparam[5]   = lparam[3]/lparam[2];
    if (material->IsMixture()) {
      lparam[1]*=lparam[0];
      TGeoMixture * mixture = (TGeoMixture*)material;
      lparam[5]=0; 
      Double_t sum =0;
      for (Int_t iel=0;iel<mixture->GetNelements();iel++){
	sum+= mixture->GetWmixt()[iel];
	lparam[5]+= mixture->GetZmixt()[iel]*mixture->GetWmixt()[iel]/mixture->GetAmixt()[iel];
      }
      lparam[5]/=sum;
    }
    gGeoManager->FindNextBoundary(length);
    snext = gGeoManager->GetStep();
  }   
  mparam[0] = bparam[0]/step;
  mparam[1] = bparam[1]/step;
  mparam[2] = bparam[2]/step;
  mparam[3] = bparam[3]/step;  
  mparam[5] = bparam[5]/step;  
  return bparam[0]/step;   
  
}
