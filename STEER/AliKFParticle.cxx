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
//#include "TMath.h"

ClassImp(AliKFParticle)


  AliKFParticle::AliKFParticle( const AliExternalTrackParam &track, Double_t bz, Int_t PID ) : fBz(bz)
{
  // Constructor from ALICE track, PID hypothesis can be provided

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

AliKFParticle::AliKFParticle( const AliESDVertex &vertex, Double_t bz): fBz(bz)
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
