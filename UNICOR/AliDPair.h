// Author: Dariusz Miskowiec <mailto:d.miskowiec@gsi.de> 2002

//=============================================================================
// particle (track) pair
//=============================================================================

#ifndef ALIDPAIR_H
#define ALIDPAIR_H

#include <cmath>
#include <TObject.h>
#include <TLorentzVector.h>

//=============================================================================
class AliDPair : public TObject {
   
 public:
   AliDPair();                       // constructor
   virtual ~AliDPair()               {printf("AliDPair object deleted\n");}
   void Set0(double m,double p,double theta,double phi)  {fP0.SetXYZM(0,0,p,m); fP0.SetTheta(theta); fP0.SetPhi(phi);}
   void Set1(double m,double p,double theta,double phi)  {fP1.SetXYZM(0,0,p,m); fP1.SetTheta(theta); fP1.SetPhi(phi);}
   void SetMXYZ0(double m,double px,double py,double pz) {fP0.SetXYZM(px,py,pz,m);}
   void SetMXYZ1(double m,double px,double py,double pz) {fP1.SetXYZM(px,py,pz,m);}
   void CalcLAB()                 {fP=fP0+fP1; fQ=fP1-fP0; fBeta=fP.BoostVector(); 
                                   fBetaz.SetXYZ(0,0,fBeta.Z()); fBetat=fBeta; fBetat.SetZ(0); 
                                   fUbeta=fBeta.Unit(); fUbetat=fBetat.Unit(); fUbetaz=fBetaz.Unit();}
   double Rapidity()        const {return fP.Rapidity();}
   double Pt()              const {return fP.Pt();}
   double Phi()             const {return fP.Phi();}
   double DTheta()          const {return fP1.Theta()-fP0.Theta();}
   double DPhi()            const {return TVector2::Phi_mpi_pi(fP1.Phi()-fP0.Phi());}
   void CalcPairCM()              {fCMp=fP; fCMp.Boost(-fBeta);  fCMq=fQ; fCMq.Boost(-fBeta);}
   void CalcLcmsCM()              {fCMp=fP; fCMp.Boost(-fBetaz); fCMq=fQ; fCMq.Boost(-fBetaz);}
   void Swap()                    {fBuf=fP0; fP0=fP1; fP1=fBuf; fQ=-fQ; fCMq=-fCMq;}
   double Minv()            const {return fP.M();}
   double Qinv2()           const {return -fQ.M2();}
   double QCM()             const {return fCMq.Vect().Mag();}
   double QCMpar()          const {return fCMq.Vect()*fUbeta;}
   double QCMper()          const {return sqrt(fabs(QCM()*QCM()-QCMpar()*QCMpar()));}
   double QCMout()          const {return +fCMq.Vect().X()*fUbetat.X()+fCMq.Vect().Y()*fUbetat.Y();}
   double QCMside()         const {return -fCMq.Vect().X()*fUbetat.Y()+fCMq.Vect().Y()*fUbetat.X();}
   double QCMlong()         const {return fCMq.Vect().Z();}
   double QCMTheta()        const {return fCMq.Theta();}
   double QCMPhi()          const {return fCMq.Phi();}
   double QCMPhiOut()       const {return atan2(QCMside(),QCMout());} // phi w.r.t. out

 protected:
   TLorentzVector fP0;            // LAB four-momentum of track 0
   TLorentzVector fP1;            // LAB four-momentum of track 1
   TLorentzVector fP;             // LAB total four-momentum
   TLorentzVector fQ;             // LAB four-momentum difference
   TVector3 fBeta;                // LAB pair velocity
   TVector3 fBetat;               // LAB pair velocity transverse
   TVector3 fBetaz;               // LAB pair velocity along z
   TVector3 fUbeta;               // LAB pair velocity direction
   TVector3 fUbetat;              // LAB pair velocity transverse direction
   TVector3 fUbetaz;              // LAB pair velocity along z direction (hm)
   TLorentzVector fCMp;           // CM total four-momentum
   TLorentzVector fCMq;           // CM four-momentum difference
   TLorentzVector fBuf;           // dummy buffer for swapping

   ClassDef(AliDPair,1)
};
//=============================================================================
#endif
