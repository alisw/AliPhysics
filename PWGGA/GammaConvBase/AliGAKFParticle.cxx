#include "AliLog.h"
#include "AliGAKFParticle.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"

#ifdef PWGGAUSEKFPARTICLE
ClassImp(AliGAKFParticle);


AliGAKFParticle::AliGAKFParticle( const AliVTrack &track, int PID ) : KFParticle() {
  //  AliError("THIS CLASS MEMBER NEEDS TO BE IMPLEMENTED");
 
  double myfP_t[8];
  double myfC_t[36]; 
  int    myfQ;    

  track.GetXYZ(myfP_t);
  track.PxPyPz(myfP_t+3);
  myfQ = int(track.Charge());
  track.GetCovarianceXYZPxPyPz( myfC_t );

  float myfP[8]; for (int i = 0; i < 8; ++i) myfP[i] = float(myfP_t[i]);
  float myfC[36]; for (int i = 0; i < 36; ++i) myfC[i] = float(myfC_t[i]);

  TParticlePDG* particlePDG = TDatabasePDG::Instance()->GetParticle(PID);
  float mass = (particlePDG) ? particlePDG->Mass() :0.13957;

  Create(myfP,myfC,myfQ,mass);
}

AliGAKFParticle::AliGAKFParticle( const AliExternalTrackParam &track, double Mass, int Charge ) : KFParticle() {
  // AliError("THIS CLASS MEMBER NEEDS TO BE IMPLEMENTED");
  double myfP_t[8];
  double myfC_t[36]; 
  int    myfQ;    
  
    track.GetXYZ(myfP_t);
    track.PxPyPz(myfP_t+3);
    float myfP[8]; for (int i = 0; i < 8; ++i) myfP[i] = float(myfP_t[i]);
    myfQ = int(track.Charge()*TMath::Abs(Charge));
    myfP[3] *= TMath::Abs(Charge);
    myfP[4] *= TMath::Abs(Charge);
    myfP[5] *= TMath::Abs(Charge);
    
    Double_t pt=1./TMath::Abs(track.GetParameter()[4]) * TMath::Abs(Charge);
    Double_t cs=TMath::Cos(track.GetAlpha()), sn=TMath::Sin(track.GetAlpha());
    Double_t r=TMath::Sqrt((1.-track.GetParameter()[2])*(1.+track.GetParameter()[2]));
    
    Double_t m00=-sn, m10=cs;
    Double_t m23=-pt*(sn + track.GetParameter()[2]*cs/r), m43=-pt*pt*(r*cs - track.GetParameter()[2]*sn);
    Double_t m24= pt*(cs - track.GetParameter()[2]*sn/r), m44=-pt*pt*(r*sn + track.GetParameter()[2]*cs);
    Double_t m35=pt, m45=-pt*pt*track.GetParameter()[3];
    
    m43*=track.GetSign();
    m44*=track.GetSign();
    m45*=track.GetSign();
    
    const Double_t *cTr = track.GetCovariance();
    
    myfC_t[0 ] = cTr[0]*m00*m00;
    myfC_t[1 ] = cTr[0]*m00*m10; 
    myfC_t[2 ] = cTr[0]*m10*m10;
    myfC_t[3 ] = cTr[1]*m00; 
    myfC_t[4 ] = cTr[1]*m10; 
    myfC_t[5 ] = cTr[2];
    myfC_t[6 ] = m00*(cTr[3]*m23 + cTr[10]*m43); 
    myfC_t[7 ] = m10*(cTr[3]*m23 + cTr[10]*m43); 
    myfC_t[8 ] = cTr[4]*m23 + cTr[11]*m43; 
    myfC_t[9 ] = m23*(cTr[5]*m23 + cTr[12]*m43)  +  m43*(cTr[12]*m23 + cTr[14]*m43);
    myfC_t[10] = m00*(cTr[3]*m24 + cTr[10]*m44); 
    myfC_t[11] = m10*(cTr[3]*m24 + cTr[10]*m44); 
    myfC_t[12] = cTr[4]*m24 + cTr[11]*m44; 
    myfC_t[13] = m23*(cTr[5]*m24 + cTr[12]*m44)  +  m43*(cTr[12]*m24 + cTr[14]*m44);
    myfC_t[14] = m24*(cTr[5]*m24 + cTr[12]*m44)  +  m44*(cTr[12]*m24 + cTr[14]*m44);
    myfC_t[15] = m00*(cTr[6]*m35 + cTr[10]*m45); 
    myfC_t[16] = m10*(cTr[6]*m35 + cTr[10]*m45); 
    myfC_t[17] = cTr[7]*m35 + cTr[11]*m45; 
    myfC_t[18] = m23*(cTr[8]*m35 + cTr[12]*m45)  +  m43*(cTr[13]*m35 + cTr[14]*m45);
    myfC_t[19] = m24*(cTr[8]*m35 + cTr[12]*m45)  +  m44*(cTr[13]*m35 + cTr[14]*m45); 
    myfC_t[20] = m35*(cTr[9]*m35 + cTr[13]*m45)  +  m45*(cTr[13]*m35 + cTr[14]*m45);

    float myfC[36]; for (int i = 0; i < 36; ++i) myfC[i] = float(myfC_t[i]);
    Create(myfP,myfC,fQ,Mass);

}



void AliGAKFParticle::GetDStoParticle(KFParticle & particle, double &  dS, double & dSp ) const {
  //AliError("THIS CLASS MEMBER NEEDS TO BE IMPLEMENTED");
  float my_dS[2] = {0.0,0.0}; // default values to define later
  float my_dsdr[4][6] = {{0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0}}; // default values to define later

  // In AliKF ONLY dS and dSp are calculated, in KF particle both dS, dSp, and the derivatives; uses an array to store distance to particle and particle p
  my_dS[0] = dS;
  my_dS[1] = dSp;
  KFParticle::GetDStoParticle(particle,my_dS, my_dsdr);
  dS  = my_dS[0];
  dSp = my_dS[1];
}


void AliGAKFParticle::TransportToPoint(double xyz[3]) {
  // AliError("THIS CLASS MEMBER NEEDS TO BE IMPLEMENTED");
    float my_store[3] = {float(xyz[0]),float(xyz[1]),float(xyz[2])};
    KFParticle::TransportToPoint(my_store); 

}

void AliGAKFParticle::ConstructGamma(const KFParticle&, const KFParticle&) {
    AliError("THIS CLASS MEMBER NEEDS TO BE IMPLEMENTED");
}

 void AliGAKFParticle::GetArmenterosPodolanski( KFParticle& positive,  KFParticle& negative, double QtAlpha[2]) {
   //    AliErrorClass("THIS CLASS MEMBER NEEDS TO BE IMPLEMENTED");
    float my_store[2] = {float(QtAlpha[0]),float(QtAlpha[1])};
    KFParticle::GetArmenterosPodolanski(positive, negative, my_store);
    QtAlpha[0] = double(my_store[0]);
    QtAlpha[1] = double(my_store[1]);

}

#endif // PWGGAUSEKFPARTICLE
