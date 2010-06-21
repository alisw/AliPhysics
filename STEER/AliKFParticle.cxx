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
#include "AliVTrack.h"
#include "AliVVertex.h"

ClassImp(AliKFParticle)

Double_t AliKFParticle::fgBz = -5.;  //* Bz compoment of the magnetic field

void AliKFParticle::Create( const Double_t Param[], const Double_t Cov[], Int_t Charge, Int_t PID )
{
  // Constructor from "cartesian" track, PID hypothesis should be provided
  //
  // Param[6] = { X, Y, Z, Px, Py, Pz } - position and momentum
  // Cov [21] = lower-triangular part of the covariance matrix:
  //
  //                (  0  .  .  .  .  . )
  //                (  1  2  .  .  .  . )
  //  Cov. matrix = (  3  4  5  .  .  . ) - numbering of covariance elements in Cov[]
  //                (  6  7  8  9  .  . )
  //                ( 10 11 12 13 14  . )
  //                ( 15 16 17 18 19 20 )
  Double_t C[21];
  for( int i=0; i<21; i++ ) C[i] = Cov[i];
  
  TParticlePDG* particlePDG = TDatabasePDG::Instance()->GetParticle(PID);
  Double_t mass = (particlePDG) ? particlePDG->Mass() :0.13957;
  
  AliKFParticleBase::Initialize( Param, C, Charge, mass );
}

AliKFParticle::AliKFParticle( const AliVTrack &track, Int_t PID )
{
  // Constructor from ALICE track, PID hypothesis should be provided

  track.XvYvZv(fP);
  track.PxPyPz(fP+3);
  fQ = track.Charge();
  track.GetCovarianceXYZPxPyPz( fC );
  Create(fP,fC,fQ,PID);
}

AliKFParticle::AliKFParticle( const AliVVertex &vertex )
{
  // Constructor from ALICE vertex

  vertex.GetXYZ( fP );
  vertex.GetCovarianceMatrix( fC );  
  fChi2 = vertex.GetChi2();
  fNDF = 2*vertex.GetNContributors() - 3;
  fQ = 0;
  fAtProductionVertex = 0;
  fIsLinearized = 0;
  fSFromDecay = 0;
}

void AliKFParticle::GetExternalTrackParam( const AliKFParticleBase &p, Double_t &X, Double_t &Alpha, Double_t P[5] ) 
{
  // Conversion to AliExternalTrackParam parameterization

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

Bool_t AliKFParticle::GetDistanceFromVertexXY( const Double_t vtx[], const Double_t Cv[], Double_t &val, Double_t &err ) const
{
  //* Calculate DCA distance from vertex (transverse impact parameter) in XY
  //* v = [xy], Cv=[Cxx,Cxy,Cyy ]-covariance matrix

  Bool_t ret = 0;
  
  Double_t mP[8];
  Double_t mC[36];
  
  Transport( GetDStoPoint(vtx), mP, mC );  

  Double_t dx = mP[0] - vtx[0];
  Double_t dy = mP[1] - vtx[1];
  Double_t px = mP[3];
  Double_t py = mP[4];
  Double_t pt = TMath::Sqrt(px*px + py*py);
  Double_t ex=0, ey=0;
  if( pt<1.e-4 ){
    ret = 1;
    pt = 1.;
    val = 1.e4;
  } else{
    ex = px/pt;
    ey = py/pt;
    val = dy*ex - dx*ey;
  }

  Double_t h0 = -ey;
  Double_t h1 = ex;
  Double_t h3 = (dy*ey + dx*ex)*ey/pt;
  Double_t h4 = -(dy*ey + dx*ex)*ex/pt;
  
  err = 
    h0*(h0*GetCovariance(0,0) + h1*GetCovariance(0,1) + h3*GetCovariance(0,3) + h4*GetCovariance(0,4) ) +
    h1*(h0*GetCovariance(1,0) + h1*GetCovariance(1,1) + h3*GetCovariance(1,3) + h4*GetCovariance(1,4) ) +
    h3*(h0*GetCovariance(3,0) + h1*GetCovariance(3,1) + h3*GetCovariance(3,3) + h4*GetCovariance(3,4) ) +
    h4*(h0*GetCovariance(4,0) + h1*GetCovariance(4,1) + h3*GetCovariance(4,3) + h4*GetCovariance(4,4) );

  if( Cv ){
    err+= h0*(h0*Cv[0] + h1*Cv[1] ) + h1*(h0*Cv[1] + h1*Cv[2] ); 
  }

  err = TMath::Sqrt(TMath::Abs(err));

  return ret;
}

Bool_t AliKFParticle::GetDistanceFromVertexXY( const Double_t vtx[], Double_t &val, Double_t &err ) const
{
  return GetDistanceFromVertexXY( vtx, 0, val, err );
}


Bool_t AliKFParticle::GetDistanceFromVertexXY( const AliKFParticle &Vtx, Double_t &val, Double_t &err ) const 
{
  //* Calculate distance from vertex [cm] in XY-plane

  return GetDistanceFromVertexXY( Vtx.fP, Vtx.fC, val, err );
}

Bool_t AliKFParticle::GetDistanceFromVertexXY( const AliVVertex &Vtx, Double_t &val, Double_t &err ) const 
{
  //* Calculate distance from vertex [cm] in XY-plane

  return GetDistanceFromVertexXY( AliKFParticle(Vtx), val, err );
}

Double_t AliKFParticle::GetDistanceFromVertexXY( const Double_t vtx[] ) const
{
  //* Calculate distance from vertex [cm] in XY-plane
  Double_t val, err;
  GetDistanceFromVertexXY( vtx, 0, val, err );
  return val;
}

Double_t AliKFParticle::GetDistanceFromVertexXY( const AliKFParticle &Vtx ) const 
{
  //* Calculate distance from vertex [cm] in XY-plane

  return GetDistanceFromVertexXY( Vtx.fP );
}

Double_t AliKFParticle::GetDistanceFromVertexXY( const AliVVertex &Vtx ) const 
{
  //* Calculate distance from vertex [cm] in XY-plane

  return GetDistanceFromVertexXY( AliKFParticle(Vtx).fP );
}

Double_t AliKFParticle::GetDistanceFromParticleXY( const AliKFParticle &p ) const 
{
  //* Calculate distance to other particle [cm]

  Double_t dS, dS1;
  GetDStoParticleXY( p, dS, dS1 );   
  Double_t mP[8], mC[36], mP1[8], mC1[36];
  Transport( dS, mP, mC ); 
  p.Transport( dS1, mP1, mC1 ); 
  Double_t dx = mP[0]-mP1[0]; 
  Double_t dy = mP[1]-mP1[1]; 
  return TMath::Sqrt(dx*dx+dy*dy);
}

Double_t AliKFParticle::GetDeviationFromParticleXY( const AliKFParticle &p ) const 
{
  //* Calculate sqrt(Chi2/ndf) deviation from other particle

  Double_t dS, dS1;
  GetDStoParticleXY( p, dS, dS1 );   
  Double_t mP1[8], mC1[36];
  p.Transport( dS1, mP1, mC1 ); 

  Double_t d[2]={ fP[0]-mP1[0], fP[1]-mP1[1] };

  Double_t sigmaS = .1+10.*TMath::Sqrt( (d[0]*d[0]+d[1]*d[1] )/
					(mP1[3]*mP1[3]+mP1[4]*mP1[4] )  );

  Double_t h[2] = { mP1[3]*sigmaS, mP1[4]*sigmaS };       
  
  mC1[0] +=h[0]*h[0];
  mC1[1] +=h[1]*h[0]; 
  mC1[2] +=h[1]*h[1]; 

  return GetDeviationFromVertexXY( mP1, mC1 )*TMath::Sqrt(2./1.);
}


Double_t AliKFParticle::GetDeviationFromVertexXY( const Double_t vtx[], const Double_t Cv[] ) const 
{
  //* Calculate sqrt(Chi2/ndf) deviation from vertex
  //* v = [xyz], Cv=[Cxx,Cxy,Cyy,Cxz,Cyz,Czz]-covariance matrix

  Double_t val, err;
  Bool_t problem = GetDistanceFromVertexXY( vtx, Cv, val, err );
  if( problem || err<1.e-20 ) return 1.e4;
  else return val/err;
}


Double_t AliKFParticle::GetDeviationFromVertexXY( const AliKFParticle &Vtx ) const  
{
  //* Calculate sqrt(Chi2/ndf) deviation from vertex
  //* v = [xyz], Cv=[Cxx,Cxy,Cyy,Cxz,Cyz,Czz]-covariance matrix

  return GetDeviationFromVertexXY( Vtx.fP, Vtx.fC );
}

Double_t AliKFParticle::GetDeviationFromVertexXY( const AliVVertex &Vtx ) const 
{
  //* Calculate sqrt(Chi2/ndf) deviation from vertex
  //* v = [xyz], Cv=[Cxx,Cxy,Cyy,Cxz,Cyz,Czz]-covariance matrix

  AliKFParticle v(Vtx);
  return GetDeviationFromVertexXY( v.fP, v.fC );
}

Double_t AliKFParticle::GetAngle  ( const AliKFParticle &p ) const 
{
  //* Calculate the opening angle between two particles

  Double_t dS, dS1;
  GetDStoParticle( p, dS, dS1 );   
  Double_t mP[8], mC[36], mP1[8], mC1[36];
  Transport( dS, mP, mC ); 
  p.Transport( dS1, mP1, mC1 ); 
  Double_t n = TMath::Sqrt( mP[3]*mP[3] + mP[4]*mP[4] + mP[5]*mP[5] );
  Double_t n1= TMath::Sqrt( mP1[3]*mP1[3] + mP1[4]*mP1[4] + mP1[5]*mP1[5] );
  n*=n1;
  Double_t a = 0;
  if( n>1.e-8 ) a = ( mP[3]*mP1[3] + mP[4]*mP1[4] + mP[5]*mP1[5] )/n;
  if (TMath::Abs(a)<1.) a = TMath::ACos(a);
  else a = (a>=0) ?0 :TMath::Pi();
  return a;
}

Double_t AliKFParticle::GetAngleXY( const AliKFParticle &p ) const 
{
  //* Calculate the opening angle between two particles in XY plane

  Double_t dS, dS1;
  GetDStoParticleXY( p, dS, dS1 );   
  Double_t mP[8], mC[36], mP1[8], mC1[36];
  Transport( dS, mP, mC ); 
  p.Transport( dS1, mP1, mC1 ); 
  Double_t n = TMath::Sqrt( mP[3]*mP[3] + mP[4]*mP[4] );
  Double_t n1= TMath::Sqrt( mP1[3]*mP1[3] + mP1[4]*mP1[4] );
  n*=n1;
  Double_t a = 0;
  if( n>1.e-8 ) a = ( mP[3]*mP1[3] + mP[4]*mP1[4] )/n;
  if (TMath::Abs(a)<1.) a = TMath::ACos(a);
  else a = (a>=0) ?0 :TMath::Pi();
  return a;
}

Double_t AliKFParticle::GetAngleRZ( const AliKFParticle &p ) const 
{
  //* Calculate the opening angle between two particles in RZ plane

  Double_t dS, dS1;
  GetDStoParticle( p, dS, dS1 );   
  Double_t mP[8], mC[36], mP1[8], mC1[36];
  Transport( dS, mP, mC ); 
  p.Transport( dS1, mP1, mC1 ); 
  Double_t nr = TMath::Sqrt( mP[3]*mP[3] + mP[4]*mP[4] );
  Double_t n1r= TMath::Sqrt( mP1[3]*mP1[3] + mP1[4]*mP1[4]  );
  Double_t n = TMath::Sqrt( nr*nr + mP[5]*mP[5] );
  Double_t n1= TMath::Sqrt( n1r*n1r + mP1[5]*mP1[5] );
  n*=n1;
  Double_t a = 0;
  if( n>1.e-8 ) a = ( nr*n1r +mP[5]*mP1[5])/n; 
  if (TMath::Abs(a)<1.) a = TMath::ACos(a);
  else a = (a>=0) ?0 :TMath::Pi();
  return a;
}


/*

#include "AliExternalTrackParam.h"

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
*/
