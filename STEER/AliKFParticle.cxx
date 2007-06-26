//----------------------------------------------------------------------------
// Implementation of the AliKFParticle class
// .
// @author  S.Gorbunov, I.Kisel
// @version 1.0
// @since   13.05.07
// 
// Class to reconstruct and store the decayed particle parameters.
// The method is described in CBM-SOFT note 2007-003, 
// ``Reconstruction of decayed particles based on the Kalman filter'', 
// http://www.gsi.de/documents/DOC-2007-May-14-1.pdf
//
// This class is ALICE interface to general mathematics in AliKFParticleCore
// 
//  -= Copyright &copy ALICE HLT Group =-
//____________________________________________________________________________


#include "AliKFParticle.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "AliExternalTrackParam.h"
#include "AliHelix.h"
//#include "TMath.h"

ClassImp(AliKFParticle);

Double_t AliKFParticle::fgBz = -5.;  //* Bz compoment of the magnetic field

AliKFParticle::AliKFParticle( const AliExternalTrackParam &track, Int_t PID )
{
  // Constructor from ALICE track, PID hypothesis should be provided
 
  TParticlePDG* particlePDG = TDatabasePDG::Instance()->GetParticle(PID);
  Double_t mass = (particlePDG) ? particlePDG->Mass() :0.13957;

  track.GetXYZ(fP);
  track.GetPxPyPz(fP+3);
  Double_t energy = TMath::Sqrt( mass*mass + fP[3]*fP[3] + fP[4]*fP[4] + fP[5]*fP[5]);
  fP[6] = energy;
  fP[7] = 0;
  fQ = (track.Get1Pt() >0 ) ?1 :-1;
  fNDF = 0;
  fChi2 = 0;
  fAtProductionVertex = 0;
  fIsLinearized = 0;
  fSFromDecay = 0;

  Double_t energyInv = 1./energy;
  Double_t 
    h0 = fP[3]*energyInv,
    h1 = fP[4]*energyInv,
    h2 = fP[5]*energyInv;

  track.GetCovarianceXYZPxPyPz( fC );

  fC[21] = h0*fC[ 6] + h1*fC[10] + h2*fC[15];
  fC[22] = h0*fC[ 7] + h1*fC[11] + h2*fC[16];
  fC[23] = h0*fC[ 8] + h1*fC[12] + h2*fC[17];
  fC[24] = h0*fC[ 9] + h1*fC[13] + h2*fC[18];
  fC[25] = h0*fC[13] + h1*fC[14] + h2*fC[19];
  fC[26] = h0*fC[18] + h1*fC[19] + h2*fC[20];
  fC[27] = h0*h0*fC[ 9] + h1*h1*fC[14] + h2*h2*fC[20] 
    + 2*(h0*h1*fC[13] + h0*h2*fC[18] + h1*h2*fC[19] );
  for( int i=28; i<36; i++ ) fC[i] = 0;
  fC[35] = 1.;
}

AliKFParticle::AliKFParticle( const AliESDVertex &vertex )
{
  // Constructor from ALICE vertex

  vertex.GetXYZ( fP );
  vertex.GetCovMatrix( fC );
  fChi2 = vertex.GetChi2();
  fNDF = 2*vertex.GetNContributors() - 3;
  fQ = 0;
  fAtProductionVertex = 0;
  fIsLinearized = 0;
  fSFromDecay = 0;
}


void AliKFParticle::GetExternalTrackParam( const AliKFParticleBase &p, Double_t &X, Double_t &Alpha, Double_t P[5] ) 
{
  Double_t cosA = p.GetPx(), sinA = p.GetPy(); 
  Double_t pt = TMath::Sqrt(cosA*cosA + sinA*sinA);
  Double_t pti = 0;
  if( pt<1.e-4 ){
    cosA = 1;
    sinA = 0;
  } else {
    pti = 1./pt;
    cosA*=pti;
    sinA*=pti;
  }
  Alpha = TMath::ATan2(sinA,cosA);  
  X   = p.GetX()*cosA + p.GetY()*sinA;   
  P[0]= p.GetY()*cosA - p.GetX()*sinA;
  P[1]= p.GetZ();
  P[2]= 0;
  P[3]= p.GetPz()*pti;
  P[4]= p.GetQ()*pti;
}


void AliKFParticle::GetDStoParticleALICE( const AliKFParticleBase &p, 
					  Double_t &DS, Double_t &DS1 ) 
  const
{ 
  DS = DS1 = 0;   
  Double_t x1, a1, x2, a2;
  Double_t par1[5], par2[5], cov[15];
  for(int i=0; i<15; i++) cov[i] = 0;
  cov[0] = cov[2] = cov[5] = cov[9] = cov[14] = .001;

  GetExternalTrackParam( *this, x1, a1, par1 );
  GetExternalTrackParam( p, x2, a2, par2 );

  AliExternalTrackParam t1(x1,a1, par1, cov);
  AliExternalTrackParam t2(x2,a2, par2, cov);

  Double_t xe1=0, xe2=0;
  t1.GetDCA( &t2, -GetFieldAlice(), xe1, xe2 );
  t1.PropagateTo( xe1, -GetFieldAlice() );
  t2.PropagateTo( xe2, -GetFieldAlice() );

  Double_t xyz1[3], xyz2[3];
  t1.GetXYZ( xyz1 );
  t2.GetXYZ( xyz2 );
  
  DS = GetDStoPoint( xyz1 );
  DS1 = p.GetDStoPoint( xyz2 );

  return;
}
