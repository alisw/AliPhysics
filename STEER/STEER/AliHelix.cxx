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
//                Implementation of the AliHelix class
//        Origin: Marian Ivanov, CERN, marian.ivanov@cern.ch
//-------------------------------------------------------------------------


#include "AliHelix.h"
#include "AliKalmanTrack.h"
#include "AliTracker.h"
#include "TMath.h"
ClassImp(AliHelix)


//_______________________________________________________________________
AliHelix::AliHelix()
{
  //
  // Default constructor
  //
  for (Int_t i =0;i<9;i++) fHelix[i]=0;
}

//_______________________________________________________________________
AliHelix::AliHelix(const AliHelix &t):TObject(t){
  //
  //
  for (Int_t i=0;i<9;i++) 
    fHelix[i]=t.fHelix[i];
}

AliHelix::AliHelix(const AliKalmanTrack &t)
{
  //
  // 
  Double_t alpha,x,cs,sn;
  t.GetExternalParameters(x,fHelix); 
  alpha=t.GetAlpha();
  //
  //circle parameters
  //PH Sometimes fP4 and fHelix[4] are very big and the calculation
  //PH of the Sqrt cannot be done. To be investigated...
  fHelix[4]=fHelix[4]/(-1000/0.299792458/AliTracker::GetBz());    // C
  cs=TMath::Cos(alpha); sn=TMath::Sin(alpha);

  Double_t xc, yc, rc;
  rc  =  1/fHelix[4];
  xc  =  x-fHelix[2]*rc;
  Double_t dummy = 1-(x-xc)*(x-xc)*fHelix[4]*fHelix[4];
  if (dummy<0) {
    AliError(Form("The argument of the Sqrt is %f => set to 0\n",dummy));
    dummy = 0;
  }
  yc  =  fHelix[0]+TMath::Sqrt(dummy)/fHelix[4];
  
  fHelix[6] = xc*cs - yc*sn;
  fHelix[7] = xc*sn + yc*cs;
  fHelix[8] =  TMath::Abs(rc);
  //
  //
  fHelix[5]=x*cs - fHelix[0]*sn;            // x0
  fHelix[0]=x*sn + fHelix[0]*cs;            // y0
  //fHelix[1]=                               // z0
  fHelix[2]=TMath::ATan2(-(fHelix[5]-fHelix[6]),fHelix[0]-fHelix[7]); // phi0
  if (fHelix[4]>0) fHelix[2]-=TMath::Pi();

  //fHelix[3]=                               // tgl
  //
  //
  fHelix[5]   = fHelix[6];
  fHelix[0]   = fHelix[7];
}


AliHelix::AliHelix(const AliExternalTrackParam &t)
{
  //
  // 
  Double_t alpha,x,cs,sn;
  const Double_t *param =t.GetParameter(); 
  for (Int_t i=0;i<5;i++) fHelix[i]=param[i]; 
  x = t.GetX();
  alpha=t.GetAlpha();
  //
  //circle parameters
  //PH Sometimes fP4 and fHelix[4] are very big and the calculation
  //PH of the Sqrt cannot be done. To be investigated...
  fHelix[4]=fHelix[4]/(-1000/0.299792458/AliTracker::GetBz());    // C
  cs=TMath::Cos(alpha); sn=TMath::Sin(alpha);

  Double_t xc, yc, rc;
  rc  =  1/fHelix[4];
  xc  =  x-fHelix[2]*rc;
  Double_t dummy = 1-(x-xc)*(x-xc)*fHelix[4]*fHelix[4];
  if (dummy<0) {
    AliError(Form("The argument of the Sqrt is %f => set to 0\n",dummy));
    dummy = 0;
  }
  yc  =  fHelix[0]+TMath::Sqrt(dummy)/fHelix[4];
  
  fHelix[6] = xc*cs - yc*sn;
  fHelix[7] = xc*sn + yc*cs;
  fHelix[8] =  TMath::Abs(rc);
  //
  //
  fHelix[5]=x*cs - fHelix[0]*sn;            // x0
  fHelix[0]=x*sn + fHelix[0]*cs;            // y0
  //fHelix[1]=                               // z0
  fHelix[2]=TMath::ASin(fHelix[2]) + alpha; // phi0
  //fHelix[3]=                               // tgl
  //
  //
  fHelix[5]   = fHelix[6];
  fHelix[0]   = fHelix[7];
}

AliHelix::AliHelix(Double_t x[3], Double_t p[3], Double_t charge, Double_t conversion)
{
  //
  //
  //
  Double_t pt = TMath::Sqrt(p[0]*p[0]+p[1]*p[1]);
  if (TMath::Abs(conversion)<0.00000001) 
    conversion = -1000/0.299792458/AliTracker::GetBz();
  //
  //  
  fHelix[4] = charge/(conversion*pt); // C
  fHelix[3] = p[2]/pt;    // tgl
  //  
  Double_t xc, yc, rc;
  rc  =  1/fHelix[4];
  xc  =  x[0]  -rc*p[1]/pt;
  yc  =  x[1]  +rc*p[0]/pt; 
  //
  fHelix[5] = x[0];   // x0
  fHelix[0] = x[1];   // y0
  fHelix[1] = x[2];   // z0
  //
  fHelix[6] = xc;
  fHelix[7] = yc;
  fHelix[8] = TMath::Abs(rc);
  //
  fHelix[5]=xc; 
  fHelix[0]=yc; 
  //
  if (TMath::Abs(p[1])<TMath::Abs(p[0])){     
    fHelix[2]=TMath::ASin(p[1]/pt);
    if (charge*yc<charge*x[1])  fHelix[2] = TMath::Pi()-fHelix[2];
  }
  else{
    fHelix[2]=TMath::ACos(p[0]/pt);
    if (charge*xc>charge*x[0])  fHelix[2] = -fHelix[2];
  }

}

void  AliHelix::GetMomentum(Double_t phase, Double_t p[4],Double_t conversion, Double_t *xr)
{
  // return  momentum at given phase
  Double_t x[3],g[3],gg[3];
  Evaluate(phase,x,g,gg);
  //  if (TMath::Abs(conversion)<0.0001) conversion = -1000/0.299792458/AliTracker::GetBz();
  if (TMath::Abs(conversion)<0.0001) conversion = TMath::Abs(1./kB2C/AliTracker::GetBz());

  Double_t mt = TMath::Sqrt(g[0]*g[0]+g[1]*g[1]);
  p[0] = fHelix[8]*g[0]/(mt*conversion);
  p[1] = fHelix[8]*g[1]/(mt*conversion);
  p[2] = fHelix[8]*g[2]/(mt*conversion);
  if (xr){
    xr[0] = x[0]; xr[1] = x[1]; xr[2] = x[2];
  }
}

void   AliHelix::GetAngle(Double_t t1, AliHelix &h, Double_t t2, Double_t angle[3])
{
  //
  //
  //
  Double_t x1[3],g1[3],gg1[3];
  Double_t x2[3],g2[3],gg2[3];
  Evaluate(t1,x1,g1,gg1);
  h.Evaluate(t2,x2,g2,gg2);

  //
  Double_t norm1r = g1[0]*g1[0]+g1[1]*g1[1];
  Double_t norm1  = TMath::Sqrt(norm1r+g1[2]*g1[2]);
  norm1r         = TMath::Sqrt(norm1r);
  //
  Double_t norm2r = g2[0]*g2[0]+g2[1]*g2[1];
  Double_t norm2  = TMath::Sqrt(norm2r+g2[2]*g2[2]);
  norm2r         = TMath::Sqrt(norm2r);
  //
  angle[0]  = (g1[0]*g2[0]+g1[1]*g2[1])/(norm1r*norm2r);   // angle in phi projection
  if (TMath::Abs(angle[0])<1.) angle[0] = TMath::ACos(angle[0]);
  else{ 
    if (angle[0]>0) angle[0] = 0;
    if (angle[0]<0) angle[0] = TMath::Pi();
  }
  //
  angle[1]  = ((norm1r*norm2r)+g1[2]*g2[2])/(norm1*norm2); // angle in rz  projection
  if (TMath::Abs(angle[1])<1.) angle[1] = TMath::ACos(angle[1]);
  else 
    angle[1]=0;

  angle[2]  = (g1[0]*g2[0]+g1[1]*g2[1]+g1[2]*g2[2])/(norm1*norm2); //3D angle
  if (TMath::Abs(angle[2])<1.) angle[2] = TMath::ACos(angle[2]);
  else 
    angle[2]=0;

  
  

}


void AliHelix::Evaluate(Double_t t,
                     Double_t r[3],  //radius vector
                     Double_t g[3],  //first defivatives
                     Double_t gg[3]) //second derivatives
{
  //--------------------------------------------------------------------
  // Calculate position of a point on a track and some derivatives at given phase
  //--------------------------------------------------------------------
  Double_t phase=fHelix[4]*t+fHelix[2];
  Double_t sn=TMath::Sin(phase), cs=TMath::Cos(phase);

  r[0] = fHelix[5] + sn/fHelix[4];
  r[1] = fHelix[0] - cs/fHelix[4];  
  r[2] = fHelix[1] + fHelix[3]*t;

  g[0] = cs; g[1]=sn; g[2]=fHelix[3];
  
  gg[0]=-fHelix[4]*sn; gg[1]=fHelix[4]*cs; gg[2]=0.;
}

Int_t     AliHelix::GetClosestPhases(AliHelix &h, Double_t phase[2][2])
{
  //
  // get phases to minimize distances
  //
  Double_t xyz0[3];
  Double_t xyz1[3];

  for (Int_t i=0;i<2;i++){  
    Evaluate(phase[i][0]  ,xyz0);	
    h.Evaluate(phase[i][1],xyz1);
    Double_t mindist = TMath::Sqrt((xyz0[0]-xyz1[0])*(xyz0[0]-xyz1[0])+
				   (xyz0[1]-xyz1[1])*(xyz0[1]-xyz1[1])+
				   (xyz0[2]-xyz1[2])*(xyz0[2]-xyz1[2]));  
    Double_t tbest[2]={phase[i][0],phase[i][1]};
    for (Int_t i0=-1;i0<=1;i0++){
      Double_t t0 = ((phase[i][0]*fHelix[4])+i0*2.*TMath::Pi())/fHelix[4];
      Evaluate(t0,xyz0);
      for (Int_t i1=-1;i1<=1;i1++){
	Double_t t1 = ((phase[i][1]*h.fHelix[4])+i1*2.*TMath::Pi())/h.fHelix[4];    
	h.Evaluate(t1,xyz1);    
	Double_t dist = TMath::Sqrt((xyz0[0]-xyz1[0])*(xyz0[0]-xyz1[0])+
				    (xyz0[1]-xyz1[1])*(xyz0[1]-xyz1[1])+
				    (xyz0[2]-xyz1[2])*(xyz0[2]-xyz1[2])); 
	if (dist<=mindist){
	  tbest[0] = t0;
	  tbest[1] = t1;
	  mindist=dist;
	}
      }
    }
    phase[i][0] = tbest[0];
    phase[i][1] = tbest[1];
  } 
  return 1;
}

Double_t  AliHelix::GetPointAngle(AliHelix &h, Double_t phase[2], const Float_t * vertex)
{
  //
  // get point angle bettwen two helixes
  // 
  Double_t r0[3],p0[4];
  Double_t r1[3],p1[4];
  GetMomentum(phase[0],p0,1,r0);
  h.GetMomentum(phase[1],p1,1,r1);
  //
  Double_t r[3] = {(r0[0]+r1[0])*0.5-vertex[0],(r0[1]+r1[1])*0.5-vertex[1],(r0[2]+r1[2])*0.5-vertex[2]};
  //intersection point - relative to the prim vertex
  Double_t p[3] = { p0[0]+p1[0], p0[1]+p1[1],p0[2]+p1[2]};
  // derivation vector
  Double_t normr = TMath::Sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
  Double_t normp = TMath::Sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);   
  Double_t pointAngle = (r[0]*p[0]+r[1]*p[1]+r[2]*p[2])/(normr*normp);
  return pointAngle;
}

Double_t  AliHelix::GetPhase(Double_t x, Double_t y )
			
{
  //
  //calculate helix param at given x,y  point
  //
  //Double_t phase2 = TMath::ATan2((y-fHelix[0]), (x-fHelix[5]))- TMath::Pi()/2.;
  Double_t phase2 = TMath::ATan2(-(x-fHelix[5]),(y-fHelix[0]));
  Int_t sign = (fHelix[4]>0)? 1:-1;
  if (sign>0) phase2 = phase2-TMath::Pi();
  //
  Float_t delta = TMath::Nint((phase2-fHelix[2])/(2.*TMath::Pi()));
  phase2-= 2*TMath::Pi()*delta;
  if ( (phase2-fHelix[2])>TMath::Pi()) phase2 -=2.*TMath::Pi();
  if ( (phase2-fHelix[2])<-TMath::Pi()) phase2+=2.*TMath::Pi();

  Double_t t     = (phase2-fHelix[2]);
  t/=fHelix[4];
  return t;
}

Int_t AliHelix::GetPhase(Double_t /*r0*/, Double_t * /*t[2]*/) 
{
  //
  //calculate helix param at given r  point - return nearest point ()
  //
  // not implemented yet
  

  return 0;
}


Double_t  AliHelix::GetPhaseZ(Double_t z0)
{
  //
  //
  return (z0-fHelix[1])/fHelix[3];
}


Int_t    AliHelix::GetRPHIintersections(AliHelix &h, Double_t phase[2][2], Double_t ri[2], Double_t cut)
{
  //--------------------------------------------------------------------
  // This function returns  phase vectors with intesection between helix (0, 1 or 2)
  // in x-y plane projection  
  //--------------------------------------------------------------------
  //    
  //  Double_t * c1 = &fHelix[6];
  //Double_t * c2 = &(h.fHelix[6]);
  //  Double_t  c1[3] = {fHelix[5],fHelix[0],fHelix[8]};

  // PH initiaziation in case of return
  phase[0][0]=phase[0][1]=phase[1][0]=phase[1][1]=0;
  ri[0]=ri[1]=1000000;

  Double_t  c1[3] = {0,0,fHelix[8]};
  Double_t  c2[3] = {h.fHelix[5]-fHelix[5],h.fHelix[0]-fHelix[0],h.fHelix[8]};

  Double_t d  = TMath::Sqrt(c2[0]*c2[0]+c2[1]*c2[1]); 
  if (d<0.000000000001) return 0;
  //
  Double_t x0[2];
  Double_t y0[2];
  //  
  if ( d>=(c1[2]+c2[2])){
    if (d>=(c1[2]+c2[2]+cut)) return 0;
    x0[0] = (d+c1[2]-c2[2])*c2[0]/(2*d)+ fHelix[5];
    y0[0] = (d+c1[2]-c2[2])*c2[1]/(2*d)+ fHelix[0];
    //    return 0;
    phase[1][0] = phase[0][0] = GetPhase(x0[0],y0[0]);
    phase[1][1] = phase[0][1] = h.GetPhase(x0[0],y0[0]);
    ri[1] = ri[0] = x0[0]*x0[0]+y0[0]*y0[0];
    return 1;
  }
  if ( (d+c2[2])<c1[2]){
    if ( (d+c2[2])+cut<c1[2]) return 0;
    //
    Double_t xx = c2[0]+ c2[0]*c2[2]/d+ fHelix[5];
    Double_t yy = c2[1]+ c2[1]*c2[2]/d+ fHelix[0]; 
    phase[1][1] = phase[0][1] = h.GetPhase(xx,yy);
    //
    Double_t xx2 = c2[0]*c1[2]/d+ fHelix[5];
    Double_t yy2 = c2[1]*c1[2]/d+ fHelix[0]; 
    phase[1][0] = phase[0][0] = GetPhase(xx2,yy2);
    ri[1] = ri[0] = xx*xx+yy*yy;
    return 1;
  }

  if ( (d+c1[2])<c2[2]){
    if ( (d+c1[2])+cut<c2[2]) return 0;
    //
    Double_t xx = -c2[0]*c1[2]/d+ fHelix[5];
    Double_t yy = -c2[1]*c1[2]/d+ fHelix[0]; 
    phase[1][1] = phase[0][1] = GetPhase(xx,yy);
    //
    Double_t xx2 = c2[0]- c2[0]*c2[2]/d+ fHelix[5];
    Double_t yy2 = c2[1]- c2[1]*c2[2]/d+ fHelix[0]; 
    phase[1][0] = phase[0][0] = h.GetPhase(xx2,yy2);
    ri[1] = ri[0] = xx*xx+yy*yy;
    return 1;
  }

  Double_t d1 = (d*d+c1[2]*c1[2]-c2[2]*c2[2])/(2.*d);
  Double_t v1 = c1[2]*c1[2]-d1*d1;
  if (v1<0) return 0;
  v1 = TMath::Sqrt(v1);
  //
  x0[0] = (c2[0]*d1+c2[1]*v1)/d + fHelix[5];
  y0[0] = (c2[1]*d1-c2[0]*v1)/d + fHelix[0];            
  //
  x0[1] = (c2[0]*d1-c2[1]*v1)/d + fHelix[5];
  y0[1] = (c2[1]*d1+c2[0]*v1)/d + fHelix[0];      
  //
  for (Int_t i=0;i<2;i++){
    phase[i][0] = GetPhase(x0[i],y0[i]);
    phase[i][1] = h.GetPhase(x0[i],y0[i]);
    ri[i] = x0[i]*x0[i]+y0[i]*y0[i];    
  }      
  return 2;
} 



Int_t   AliHelix::LinearDCA(AliHelix &h, Double_t &t1, Double_t &t2, 
		      Double_t &R, Double_t &dist)
{
  //
  //
  // find intersection using linear approximation
  Double_t r1[3],g1[3],gg1[3];
  Double_t r2[3],g2[3],gg2[3];
  //
  Evaluate(t1,r1,g1,gg1);
  h.Evaluate(t2,r2,g2,gg2);
  // 
  Double_t g1_2 = g1[0]*g1[0] +g1[1]*g1[1] +g1[2]*g1[2];
  Double_t g2_2 = g2[0]*g2[0] +g2[1]*g2[1] +g2[2]*g2[2];
  Double_t g1x2 = g1[0]*g2[0] +g1[1]*g2[1] +g1[2]*g2[2];  
  Double_t det  = g1_2*g2_2   - g1x2*g1x2;
  //  
  if (TMath::Abs(det)>0){
    //
    Double_t r1g1 = r1[0]*g1[0] +r1[1]*g1[1] +r1[2]*g1[2];    
    Double_t r2g1 = r2[0]*g1[0] +r2[1]*g1[1] +r2[2]*g1[2];      
    Double_t r1g2 = r1[0]*g2[0] +r1[1]*g2[1] +r1[2]*g2[2];
    Double_t r2g2 = r2[0]*g2[0] +r2[1]*g2[1] +r2[2]*g2[2];
    //    
    Double_t dt    = - ( g2_2*(r1g1-r2g1) - g1x2*(r1g2-r2g2)) / det;      
    Double_t dp    = - ( g1_2*(r2g2-r1g2) - g1x2*(r2g1-r1g1)) / det;
    //
    t1+=dt;
    t2+=dp;
    Evaluate(t1,r1);
    h.Evaluate(t2,r2);
    //
    dist = (r1[0]-r2[0])*(r1[0]-r2[0])+
				  (r1[1]-r2[1])*(r1[1]-r2[1])+
				  (r1[2]-r2[2])*(r1[2]-r2[2]);    
    R = ((r1[0]+r2[0])*(r1[0]+r2[0])+(r1[1]+r2[1])*(r1[1]+r2[1]))/4.;
  }     
  return 0;
}


Int_t  AliHelix::ParabolicDCA(AliHelix&h,  //helixes
			       Double_t &t1, Double_t &t2, 
			       Double_t &R, Double_t &dist, Int_t iter)
{
  //
  //
  // find intersection using linear fit
  Double_t r1[3],g1[3],gg1[3];
  Double_t r2[3],g2[3],gg2[3];
  //
  Evaluate(t1,r1,g1,gg1);
  h.Evaluate(t2,r2,g2,gg2);

  //
  Double_t dx2=1.;
  Double_t dy2=1.;
  Double_t dz2=1.;
  //
  Double_t dx=r2[0]-r1[0], dy=r2[1]-r1[1], dz=r2[2]-r1[2];
  Double_t dm=dx*dx/dx2 + dy*dy/dy2 + dz*dz/dz2;
  //

 iter++;
 while (iter--) {
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
    
    //if (TMath::Abs(dt1)/(TMath::Abs(t1)+1.e-3) < 1.e-4)
    //  if (TMath::Abs(dt2)/(TMath::Abs(t2)+1.e-3) < 1.e-4) {
    //	break;
    //  }
    
    Double_t dd=dm;
    for (Int_t div=1 ; div<512 ; div*=2) {
      Evaluate(t1+dt1,r1,g1,gg1);
      h.Evaluate(t2+dt2,r2,g2,gg2);
      dx=r2[0]-r1[0]; dy=r2[1]-r1[1]; dz=r2[2]-r1[2];
      dd=dx*dx/dx2 + dy*dy/dy2 + dz*dz/dz2;
      if (dd<dm) break;
      dt1*=0.5; dt2*=0.5;
  //     if (div==0){
// 	div =1;
//       }
//       if (div>512) {	  
// 	break;
//       }   
    }
    dm=dd;
    t1+=dt1;
    t2+=dt2;
 }
 Evaluate(t1,r1,g1,gg1);
 h.Evaluate(t2,r2,g2,gg2);
 //
 dist = (r1[0]-r2[0])*(r1[0]-r2[0])+
   (r1[1]-r2[1])*(r1[1]-r2[1])+
   (r1[2]-r2[2])*(r1[2]-r2[2]);    
 
 R = ((r1[0]+r2[0])*(r1[0]+r2[0])+(r1[1]+r2[1])*(r1[1]+r2[1]))/4;
 return 0;
 
}


Int_t  AliHelix::ParabolicDCA2(AliHelix&h,  //helixes
			       Double_t &t1, Double_t &t2, 
			       Double_t &R, Double_t &dist,  Double_t err[3], Int_t iter)
{
  //
  //
  // find intersection using linear fit
  Double_t r1[3],g1[3],gg1[3];
  Double_t r2[3],g2[3],gg2[3];
  //
  Evaluate(t1,r1,g1,gg1);
  h.Evaluate(t2,r2,g2,gg2);

  //
  Double_t dx2=err[0];
  Double_t dy2=err[1];
  Double_t dz2=err[2];
  //
  Double_t dx=r2[0]-r1[0], dy=r2[1]-r1[1], dz=r2[2]-r1[2];
  Double_t dm=dx*dx/dx2 + dy*dy/dy2 + dz*dz/dz2;
  //

 iter++;
 while (iter--) {
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
    
    //if (TMath::Abs(dt1)/(TMath::Abs(t1)+1.e-3) < 1.e-4)
    //  if (TMath::Abs(dt2)/(TMath::Abs(t2)+1.e-3) < 1.e-4) {
    //	break;
    //  }
    
    Double_t dd=dm;
    for (Int_t div=1 ; div<512 ; div*=2) {
      Evaluate(t1+dt1,r1,g1,gg1);
      h.Evaluate(t2+dt2,r2,g2,gg2);
      dx=r2[0]-r1[0]; dy=r2[1]-r1[1]; dz=r2[2]-r1[2];
      dd=dx*dx/dx2 + dy*dy/dy2 + dz*dz/dz2;
      if (dd<dm) break;
      dt1*=0.5; dt2*=0.5;
      if (div==0){
	div =1;
      }
      if (div>512) {	  
	break;
      }   
    }
    dm=dd;
    t1+=dt1;
    t2+=dt2;
 }
 Evaluate(t1,r1,g1,gg1);
 h.Evaluate(t2,r2,g2,gg2);
 //
 dist = (r1[0]-r2[0])*(r1[0]-r2[0])+
   (r1[1]-r2[1])*(r1[1]-r2[1])+
   (r1[2]-r2[2])*(r1[2]-r2[2]);    
 
 R = ((r1[0]+r2[0])*(r1[0]+r2[0])+(r1[1]+r2[1])*(r1[1]+r2[1]))/4;
 return 0;
 
}

