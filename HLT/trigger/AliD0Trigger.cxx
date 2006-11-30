#include "AliD0Trigger.h"
#include <TMath.h>
#include "AliITStrackV2.h"
#include "AliHLTTransform.h"
#include <TVector3.h>
#include <iostream.h>

ClassImp(AliD0Trigger)

AliD0Trigger::AliD0Trigger()
{
  posTrack=0;
  negTrack=0;
}
AliD0Trigger::AliD0Trigger(double c[7],double b,double pv[3])
{
  posTrack=0;
  negTrack=0;
  cutV0low=c[0]; 
  cutV0high=c[1];
  cutInvMass=c[2]; 
  cutPointAngle=c[3];
  cutd0d0=c[4];
  cutCosThetaStar=c[5];
  cutpTchild=c[6];
  Bfield=b;
  primaryVertex[0]=pv[0];
  primaryVertex[1]=pv[1];
  primaryVertex[2]=pv[2];
}
AliD0Trigger::AliD0Trigger(double c[7],double b,double pv[3],AliITStrackV2 * posT, AliITStrackV2 * negT)
{
  posTrack=posT;
  negTrack=negT;
  cutV0low=c[0]; 
  cutV0high=c[1];
  cutInvMass=c[2]; 
  cutPointAngle=c[3];
  cutd0d0=c[4];
  cutCosThetaStar=c[5];
  cutpTchild=c[6];
  Bfield=b;
  primaryVertex[0]=pv[0];
  primaryVertex[1]=pv[1];
  primaryVertex[2]=pv[2];
}
AliD0Trigger::AliD0Trigger(AliITStrackV2 * posT, AliITStrackV2 * negT)
{ 
  posTrack=posT;
  negTrack=negT;
}


AliD0Trigger::~AliD0Trigger()
{
 
}

void AliD0Trigger::SetTracks(AliITStrackV2 * posT, AliITStrackV2 * negT){
  posTrack=posT;
  negTrack=negT;
}
bool AliD0Trigger::FindInvMass() {
   
  double mD0 = 1.8645;
  double mPi = 0.13957;
  double mKa = 0.49368;
  
  double mom2[2],momTot2,energy[2];
  
  mom2[0]=pow(momenta[0],2)+pow(momenta[1],2)+pow(momenta[2],2);
  mom2[1]=pow(momenta[3],2)+pow(momenta[4],2)+pow(momenta[5],2);
  
  momTot2 = pow(momenta[0]+momenta[3],2)+pow(momenta[1]+momenta[4],2)+pow(momenta[2]+momenta[5],2);
  
  // D0 -> K- Pi+
  energy[1] = sqrt(pow(mKa,2)+mom2[1]);
  energy[0] = sqrt(pow(mPi,2)+mom2[0]);
  
  double minvD0 = sqrt((energy[0]+energy[1])*(energy[0]+energy[1])-momTot2);
   
  // D0bar -> K+ Pi-
  energy[0] = sqrt(mKa*mKa+mom2[0]);
  energy[1] = sqrt(mPi*mPi+mom2[1]);
  
  double minvD0bar = sqrt((energy[0]+energy[1])*(energy[0]+energy[1])-momTot2);
  
  if(fabs(minvD0-mD0) < cutInvMass)    return true;
  if(fabs(minvD0bar-mD0) < cutInvMass) return true;
  return false;
  
}
bool AliD0Trigger::FindV0(){

  bool goodV0=false;

  double Gxpos=0,Gypos=0,Gxneg=0,Gyneg=0;
  double r1=0,r2=0,a1=0,a2=0,b1=0,b2=0;
  double q=0,p=0,t=0,Fa=0,Fb=0,Fc=0;
  double y1=0,y2=0,x1=0,x2=0,x3=0,x4=0;
  double V01=0,V02=0,V03=0,V04=0;
  double V0[4]={0,0,0,0};
  bestV0[0]=0;  bestV0[1]=0;  bestV0[2]=0;
  
  //Do not get right values
  Gxpos=posTrack->GetX()*cos(posTrack->GetAlpha())-posTrack->GetY()*sin(posTrack->GetAlpha());
  Gypos=posTrack->GetX()*sin(posTrack->GetAlpha())+posTrack->GetY()*cos(posTrack->GetAlpha());
  Gxneg=negTrack->GetX()*cos(negTrack->GetAlpha())-negTrack->GetY()*sin(negTrack->GetAlpha());
  Gyneg=negTrack->GetX()*sin(negTrack->GetAlpha())+negTrack->GetY()*cos(negTrack->GetAlpha());
    
  r1=fabs(1/(AliHLTTransform::GetBFact()*Bfield*posTrack->Get1Pt()));
  r2=fabs(1/(AliHLTTransform::GetBFact()*Bfield*negTrack->Get1Pt()));
  
  a1=Gxpos-(r1*cos(posTrack->GetAlpha()+asin(posTrack->GetSnp())+TMath::PiOver2()));
  a2=Gxneg-(r2*cos(negTrack->GetAlpha()+asin(negTrack->GetSnp())-TMath::PiOver2()));
  b1=Gypos-(r1*sin(posTrack->GetAlpha()+asin(posTrack->GetSnp())+TMath::PiOver2()));
  b2=Gyneg-(r2*sin(negTrack->GetAlpha()+asin(negTrack->GetSnp())-TMath::PiOver2()));

  //double a3=a1-a2;
  //double b3=b1-b2;
  //double r3=sqrt(pow(a3,2)+pow(b3,2));
  
  //if(r3<r1+r2){
  q=a1-a2;
  p=2*(b1-b2);
  t=pow(r2,2)-pow(r2,2)+pow(b1,2)-pow(b2,2)-pow(q,2);
  
  Fa=pow(p,2)+(4*pow(q,2));
  Fb=-(2*t*p+8*pow(q,2)*b1);
  Fc=pow(t,2)-(4*pow(q,2)*pow(r1,2))+(4*pow(q,2)*pow(b1,2));
  
  if(pow(Fb,2)-(4*Fa*Fc)>=0){
    y1=(-Fb+(sqrt(pow(Fb,2)-(4*Fa*Fc))))/(2*Fa);  
    y2=(-Fb-(sqrt(pow(Fb,2)-(4*Fa*Fc))))/(2*Fa);  
    
    x1=sqrt(pow(r1,2)-pow((y1-b1),2))+a1;    
    x2=sqrt(pow(r1,2)-pow((y2-b1),2))+a1;    
    x3=-sqrt(pow(r1,2)-pow((y1-b1),2))+a1;
    x4=-sqrt(pow(r1,2)-pow((y2-b1),2))+a1;
    
    V01=sqrt(pow(x1-primaryVertex[0],2)+pow(y1-primaryVertex[1],2));
    V02=sqrt(pow(x2-primaryVertex[0],2)+pow(y2-primaryVertex[1],2));
    V03=sqrt(pow(x3-primaryVertex[0],2)+pow(y1-primaryVertex[1],2));  
    V04=sqrt(pow(x4-primaryVertex[0],2)+pow(y2-primaryVertex[1],2));  
    
    V0[0]=V01;V0[1]=V02;V0[2]=V03;V0[3]=V04;

    int index=0;
    
    double nearestV0=V0[0];

    for(int i=0;i<4;i++){                        
      if(nearestV0>V0[i]){
	nearestV0=V0[i];
	index=i;
      }
    }
    cout<<"index: "<<index<<endl;
    if(nearestV0<cutV0high){
      if(nearestV0>cutV0low){  
	goodV0=true;
	if(index==0){
	  bestV0[0]=x1;   bestV0[1]=y1;   bestV0[2]=0.;
	}
	if(index==1){
	  bestV0[0]=x2;   bestV0[1]=y2;   bestV0[2]=0.;
	}
	if(index==2){
	  bestV0[0]=x3;   bestV0[1]=y1;   bestV0[2]=0.;
	}
	if(index==3){
	  bestV0[0]=x4;   bestV0[1]=y2;   bestV0[2]=0.;
	}
      }
    }
    //}
  }
  
  cout<<"My V0: x: "<<bestV0[0]<<" y: "<<bestV0[1]<<" z: "<<bestV0[2]<<endl;

  return goodV0;	
}
void AliD0Trigger::FindMomentaAtVertex(){

  //This method moves the momenta to the secondary vertex

  double r1=fabs(1/(AliHLTTransform::GetBFact()*Bfield*posTrack->Get1Pt()));
  double r2=fabs(1/(AliHLTTransform::GetBFact()*Bfield*negTrack->Get1Pt()));
  
  double Gx1=posTrack->GetX()*cos(posTrack->GetAlpha())-posTrack->GetY()*sin(posTrack->GetAlpha());
  double Gy1=posTrack->GetX()*sin(posTrack->GetAlpha())+posTrack->GetY()*cos(posTrack->GetAlpha());
  double Gx2=negTrack->GetX()*cos(negTrack->GetAlpha())-negTrack->GetY()*sin(negTrack->GetAlpha());
  double Gy2=negTrack->GetX()*sin(negTrack->GetAlpha())+negTrack->GetY()*cos(negTrack->GetAlpha());
  
  double centerx1=Gx1-(r1*cos(posTrack->GetAlpha()+asin(posTrack->GetSnp())+TMath::PiOver2()));
  double centerx2=Gx2-(r2*cos(negTrack->GetAlpha()+asin(negTrack->GetSnp())-TMath::PiOver2()));
  double centery1=Gy1-(r1*sin(posTrack->GetAlpha()+asin(posTrack->GetSnp())+TMath::PiOver2()));
  double centery2=Gy2-(r2*sin(negTrack->GetAlpha()+asin(negTrack->GetSnp())-TMath::PiOver2()));
  
  double a1=sqrt(pow(Gx1-bestV0[0],2)+pow(Gy1-bestV0[1],2));
  double b1=sqrt(pow(bestV0[0]-centerx1,2)+pow(bestV0[1]-centery1,2));
  double c1=sqrt(pow(Gx1-centerx1,2)+pow(Gy1-centery1,2));
  double a2=sqrt(pow(Gx2-bestV0[0],2)+pow(Gy2-bestV0[1],2));
  double b2=sqrt(pow(bestV0[0]-centerx2,2)+pow(bestV0[1]-centery2,2));
  double c2=sqrt(pow(Gx2-centerx2,2)+pow(Gy2-centery2,2));
  
  double alpha1=acos((pow(b1,2)+pow(c1,2)-pow(a1,2))/(2*b1*c1));
  double alpha2=acos((pow(b2,2)+pow(c2,2)-pow(a2,2))/(2*b2*c2));

  double sign1=0;
  double sign2=0;
  
  if((sqrt(pow(Gx1,2)+pow(Gy1,2)))>(sqrt(pow(bestV0[0],2)+pow(bestV0[1],2))))
    sign1=-1;
  else
    sign1=1;
  if((sqrt(pow(Gx2,2)+pow(Gy2,2)))>(sqrt(pow(bestV0[0],2)+pow(bestV0[1],2))))
    sign2=1;
  else
    sign2=-1;
  
  double psi1=(posTrack->GetAlpha()+asin(posTrack->GetSnp()))+(sign1*alpha1);
  double psi2=(negTrack->GetAlpha()+asin(negTrack->GetSnp()))+(sign2*alpha2);
  
  double ptP = 1./fabs(posTrack->Get1Pt());
  momenta[0] = ptP*cos(psi1); 
  momenta[1] = ptP*sin(psi1);
  momenta[2] = ptP*posTrack->GetTgl();
  
  double ptN = 1./fabs(negTrack->Get1Pt());
  momenta[3] = ptN*cos(psi2); 
  momenta[4] = ptN*sin(psi2);
  momenta[5] = ptN*negTrack->GetTgl();
  
}
Bool_t AliD0Trigger::PointingAngle(){

  TVector3 mom(Px(),Py(),Pz());
  TVector3 flight(bestV0[0]-primaryVertex[0],bestV0[1]-primaryVertex[1],bestV0[2]-primaryVertex[2]);

  double pta = mom.Angle(flight);

  if(cos(pta)>cutPointAngle)
    return true;
  else
    return false;
}
void AliD0Trigger::SetMomenta(double m[6])
{
  momenta[0] = m[0]; 
  momenta[1] = m[1]; 
  momenta[2] = m[2]; 
  momenta[3] = m[3]; 
  momenta[4] = m[4]; 
  momenta[5] = m[5]; 
}
void AliD0Trigger::FindMomentaOffline()
{
  double ptP,alphaP,phiP,ptN,alphaN,phiN;
  // momenta of the tracks at the vertex
  ptP = 1./TMath::Abs(posTrack->Get1Pt());
  alphaP = posTrack->GetAlpha();
  phiP = alphaP+TMath::ASin(posTrack->GetSnp());
  momenta[0] = ptP*TMath::Cos(phiP); 
  momenta[1] = ptP*TMath::Sin(phiP);
  momenta[2] = ptP*posTrack->GetTgl();
  
  ptN = 1./TMath::Abs(negTrack->Get1Pt());
  alphaN = negTrack->GetAlpha();
  phiN = alphaN+TMath::ASin(negTrack->GetSnp());
  momenta[3] = ptN*TMath::Cos(phiN); 
  momenta[4] = ptN*TMath::Sin(phiN);
  momenta[5] = ptN*negTrack->GetTgl();  
}
bool AliD0Trigger::FindV0offline(double v[3])
{

  bool goodV0=false;
  double r=0;
  
  bestV0[0]=v[0];
  bestV0[1]=v[1];
  bestV0[2]=v[2];

  r=sqrt(pow(bestV0[0]-primaryVertex[0],2)+pow(bestV0[1]-primaryVertex[1],2)+pow(bestV0[2]-primaryVertex[2],2));

  if(r<cutV0high){
    if(r>cutV0low){
      goodV0=true;
    }
  }

  return goodV0;

}
bool AliD0Trigger::d0d0()
{
  bool goodd0=false;
  double d00=0, d01=0;

  d00 =  10000.*posTrack->GetD(primaryVertex[0],primaryVertex[1]);
  d01 = -10000.*negTrack->GetD(primaryVertex[0],primaryVertex[1]);

  if(d00*d01<cutd0d0){
    goodd0=true;
  }

  return goodd0;
}
double AliD0Trigger::Energy()
{
  double kMD0 = 1.8645;  // D0  mass
  return sqrt(P()*P()+kMD0*kMD0);
}
bool AliD0Trigger::CosThetaStar()
{
  bool goodtheta=false;
  double kMD0 = 1.8645;  // D0  mass
  double kMK  = 0.49368; // K+  mass
  double kMPi = 0.13957; // Pi+ mass
  double qL=0,qL2=0,ctsD0=0,ctsD0bar=0;

  Double_t pStar = TMath::Sqrt(TMath::Power(kMD0*kMD0-kMK*kMK-kMPi*kMPi,2.)-4.*kMK*kMK*kMPi*kMPi)/(2.*kMD0);

  Double_t beta = P()/Energy();
  Double_t gamma = Energy()/kMD0;

  TVector3 mom(momenta[0],momenta[1],momenta[2]);
  TVector3 mom2(momenta[3],momenta[4],momenta[5]);
  TVector3 momD(Px(),Py(),Pz());

  qL = mom.Dot(momD)/momD.Mag();

  ctsD0 = (qL/gamma-beta*TMath::Sqrt(pStar*pStar+kMK*kMK))/pStar;

  qL2 = mom.Dot(momD)/momD.Mag();

  ctsD0bar = (qL2/gamma-beta*TMath::Sqrt(pStar*pStar+kMK*kMK))/pStar;

  if(fabs(ctsD0) > cutCosThetaStar || fabs(ctsD0bar) > cutCosThetaStar){
      goodtheta=true;
    }
  return goodtheta;
}
bool AliD0Trigger::pTchild()
{
  bool goodpT=false;
  double pT1=0,pT2=0;
  
  pT1=sqrt(pow(momenta[0],2)+pow(momenta[1],2));
  pT2=sqrt(pow(momenta[3],2)+pow(momenta[4],2));

  if(pT1>cutpTchild){
    if(pT2>cutpTchild){
      goodpT=true;
    }
  }
  return goodpT;
}
void AliD0Trigger::SetV0(double v[3])
{
  bestV0[0]=v[0];
  bestV0[1]=v[1];
  bestV0[2]=v[2];
}
