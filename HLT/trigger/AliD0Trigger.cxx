#include "AliD0Trigger.h"
#include <TMath.h>
#include "AliITStrackV2.h"
#include "AliL3Transform.h"
#include <TVector3.h>

ClassImp(AliD0Trigger)

AliD0Trigger::AliD0Trigger()
{
  posTrack=0;
  negTrack=0;
}
AliD0Trigger::AliD0Trigger(double c[4],double b,double pv[3])
{
  posTrack=0;
  negTrack=0;
  cutV0low=c[0]; 
  cutV0high=c[1];
  cutInvMass=c[2]; 
  cutPointAngle=c[3];
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
  
  double Gxpos=posTrack->GetX()*cos(posTrack->GetAlpha())-posTrack->GetY()*sin(posTrack->GetAlpha());
  double Gypos=posTrack->GetX()*sin(posTrack->GetAlpha())+posTrack->GetY()*cos(posTrack->GetAlpha());
  double Gxneg=negTrack->GetX()*cos(negTrack->GetAlpha())-negTrack->GetY()*sin(negTrack->GetAlpha());
  double Gyneg=negTrack->GetX()*sin(negTrack->GetAlpha())+negTrack->GetY()*cos(negTrack->GetAlpha());
    
  double r1=fabs(1/(AliL3Transform::GetBFact()*Bfield*posTrack->Get1Pt()));
  double r2=fabs(1/(AliL3Transform::GetBFact()*Bfield*negTrack->Get1Pt()));
  //double a1=posTrack->GetX()-(r1*cos(asin(posTrack->GetSnp())+TMath::PiOver2()));
  //double a2=negTrack->GetX()-(r1*cos(asin(negTrack->GetSnp())-TMath::PiOver2()));
  //double b1=posTrack->GetY()-(r1*sin(asin(posTrack->GetSnp())+TMath::PiOver2()));
  //double b2=negTrack->GetY()-(r1*sin(asin(negTrack->GetSnp())-TMath::PiOver2()));
  
  double a1=Gxpos-(r1*cos(posTrack->GetAlpha()+asin(posTrack->GetSnp())+TMath::PiOver2()));
  double a2=Gxneg-(r2*cos(negTrack->GetAlpha()+asin(negTrack->GetSnp())-TMath::PiOver2()));
  double b1=Gypos-(r1*sin(posTrack->GetAlpha()+asin(posTrack->GetSnp())+TMath::PiOver2()));
  double b2=Gyneg-(r2*sin(negTrack->GetAlpha()+asin(negTrack->GetSnp())-TMath::PiOver2()));

  double a3=a1-a2;
  double b3=b1-b2;
  double r3=sqrt(pow(a3,2)+pow(b3,2));
  
  if(r3<r1+r2){
    double q=a1-a2;
    double p=2*(b1-b2);
    double t=pow(r2,2)-pow(r2,2)+pow(b1,2)-pow(b2,2)-pow(q,2);
 
    double Fa=pow(p,2)+(4*pow(q,2));
    double Fb=-(2*t*p+8*pow(q,2)*b1);
    double Fc=pow(t,2)-(4*pow(q,2)*pow(r1,2))+(4*pow(q,2)*pow(b1,2));
 
    if(pow(Fb,2)-(4*Fa*Fc)>=0){
      double y1=(-Fb+(sqrt(pow(Fb,2)-(4*Fa*Fc))))/(2*Fa);  //noe feil her. floating point
      double y2=(-Fb-(sqrt(pow(Fb,2)-(4*Fa*Fc))))/(2*Fa);  //trolig negativ under rot 
      
      double x1=sqrt(pow(r1,2)-pow((y1-b1),2))+a1;    
      double x2=sqrt(pow(r1,2)-pow((y2-b1),2))+a1;    
      double x3=-sqrt(pow(r1,2)-pow((y1-b1),2))+a1;
      double x4=-sqrt(pow(r1,2)-pow((y2-b1),2))+a1;
      
      double V01=sqrt(pow(x1,2)+pow(y1,2));
      double V02=sqrt(pow(x2,2)+pow(y2,2));
      double V03=sqrt(pow(x3,2)+pow(y1,2));
      double V04=sqrt(pow(x4,2)+pow(y2,2));
      
      double V0[4]={V01,V02,V03,V04};
      int index=0;
      
      double nearestV0=V0[0];
      for(int i=1;i<4;i++){
	if(nearestV0>V0[i]){
	  nearestV0=V0[i];
	  index=i;
	}
      }
      if(nearestV0<cutV0high && nearestV0>cutV0low){  
      	goodV0=true;
	switch (index){
	case 0 : bestV0[0]=x1;bestV0[1]=y1;bestV0[0]=0.;break;
	case 1 : bestV0[0]=x2;bestV0[1]=y2;bestV0[0]=0.;break;
	case 2 : bestV0[0]=x3;bestV0[1]=y1;bestV0[0]=0.;break;
	case 3 : bestV0[0]=x4;bestV0[1]=y2;bestV0[0]=0.;break;
	default : printf("Can't set V0");
	}
      }
    }
  }
  return goodV0;	
}
void AliD0Trigger::FindMomentaAtVertex(){

  double r1=fabs(1/(AliL3Transform::GetBFact()*Bfield*posTrack->Get1Pt()));
  double r2=fabs(1/(AliL3Transform::GetBFact()*Bfield*negTrack->Get1Pt()));
  
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

  TVector3 mom(momenta[0]+momenta[3],momenta[1]+momenta[4],momenta[2]+momenta[5]);
  TVector3 flight(bestV0[0]-primaryVertex[0],bestV0[1]-primaryVertex[1],bestV0[2]-primaryVertex[2]);

  double pta = mom.Angle(flight);

  if(cos(pta)<cutPointAngle)
    return true;
  else
    return false;
}
