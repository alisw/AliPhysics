// Author: Dariusz Miskowiec <mailto:d.miskowiec@gsi.de> 2002

#ifndef ALIDPAIR_H
#define ALIDPAIR_H

#include <TObject.h>
#include "TLorentzVector.h"

//=============================================================================
class AliDPair : public TObject {
   
 protected:
   TLorentzVector p0;             // LAB four-momentum of track 0
   TLorentzVector p1;             // LAB four-momentum of track 1
   TLorentzVector p;              // LAB total four-momentum
   TLorentzVector q;              // LAB four-momentum difference
   TVector3 beta;                 // LAB pair velocity
   TVector3 betat;                // LAB pair velocity transverse
   TVector3 betaz;                // LAB pair velocity along z
   TVector3 ubeta;                // LAB pair velocity direction
   TVector3 ubetat;               // LAB pair velocity transverse direction
   TVector3 ubetaz;               // LAB pair velocity along z direction (hm)
   TLorentzVector CMp;            // CM total four-momentum
   TLorentzVector CMq;            // CM four-momentum difference
   TLorentzVector buf;            // dummy buffer for swapping

 public:
   AliDPair();                       // constructor
   virtual ~AliDPair()               {printf("AliDPair object deleted\n");}
   void Set0(double m,double p,double theta,double phi)  {p0.SetXYZM(0,0,p,m); p0.SetTheta(theta); p0.SetPhi(phi);}
   void Set1(double m,double p,double theta,double phi)  {p1.SetXYZM(0,0,p,m); p1.SetTheta(theta); p1.SetPhi(phi);}
   void SetMXYZ0(double m,double px,double py,double pz) {p0.SetXYZM(px,py,pz,m);}
   void SetMXYZ1(double m,double px,double py,double pz) {p1.SetXYZM(px,py,pz,m);}
   void CalcLAB()                 {p=p0+p1; q=p1-p0; beta=p.BoostVector(); 
                                   betaz.SetXYZ(0,0,beta.Z()); betat=beta; betat.SetZ(0); 
                                   ubeta=beta.Unit(); ubetat=betat.Unit(); ubetaz=betaz.Unit();}
   double Rapidity()              {return p.Rapidity();}
   double Pt()                    {return p.Pt();}
   double Phi()                   {return p.Phi();}
   double DTheta()                {return p1.Theta()-p0.Theta();}
   double DPhi()                  {return TVector2::Phi_mpi_pi(p1.Phi()-p0.Phi());}
   void CalcPairCM()              {CMp=p; CMp.Boost(-beta);  CMq=q; CMq.Boost(-beta);}
   void CalcLcmsCM()              {CMp=p; CMp.Boost(-betaz); CMq=q; CMq.Boost(-betaz);}
   void Swap()                    {buf=p0; p0=p1; p1=buf; q=-q; CMq=-CMq;}
   double Minv()                  {return p.M();}
   double Qinv2()                 {return -q.M2();}
   double QCM()                   {return CMq.Vect().Mag();}
   double QCMpar()                {return CMq.Vect()*ubeta;}
   double QCMper()                {return TMath::Sqrt(TMath::Abs(QCM()*QCM()-QCMpar()*QCMpar()));}
   double QCMout()                {return +CMq.Vect().X()*ubetat.X()+CMq.Vect().Y()*ubetat.Y();}
   double QCMside()               {return -CMq.Vect().X()*ubetat.Y()+CMq.Vect().Y()*ubetat.X();}
   double QCMlong()               {return CMq.Vect().Z();}
   double QCMTheta()              {return CMq.Theta();}
   double QCMPhi()                {return CMq.Phi();}
   double QCMPhiOut()             {return TMath::ATan2(QCMside(),QCMout());} // phi w.r.t. out

   ClassDef(AliDPair,1)
};
//=============================================================================
#endif
