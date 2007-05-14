//----------------------------------------------------------------------------
// Implementation of the AliKFVertex class
// .
// @author  S.Gorbunov, I.Kisel
// @version 1.0
// @since   13.05.07
// 
// Class to reconstruct and store primary and secondary vertices
// The method is described in CBM-SOFT note 2007-003, 
// ``Reconstruction of decayed particles based on the Kalman filter'', 
// http://www.gsi.de/documents/DOC-2007-May-14-1.pdf
//
// This class is ALICE interface to general mathematics in AliKFParticleCore
// 
//  -= Copyright &copy ALICE HLT Group =-
//____________________________________________________________________________


#include "AliKFVertex.h"


ClassImp(AliKFVertex)


AliKFVertex::AliKFVertex( const AliESDVertex &vertex )
{
  // Constructor from ALICE ESD vertex

  vertex.GetXYZ( fP );
  vertex.GetCovMatrix( fC );
  fChi2 = vertex.GetChi2();
  fNDF = 2*vertex.GetNContributors() - 3;
  fQ = 0;
  fAtProductionVertex = 0;
  fIsLinearized = 0;
  fSFromDecay = 0;
}

void AliKFVertex::ConstructPrimaryVertex( const AliKFParticle *vDaughters[], 
					  int NDaughters, Double_t ChiCut  )
{
  //* Primary vertex finder with simple rejection of outliers

  if( NDaughters<2 ) return;
  Construct( &(vDaughters[0]), NDaughters );
    
  Int_t nt=NDaughters;
  for( Int_t it=0; it<NDaughters; it++){
    if( nt<3) return;
    const AliKFParticle &p = *(vDaughters[it]);
    AliKFVertex tmp = *this - p;
    Double_t d = p.GetDeviationFromVertex( tmp );
    if( d>ChiCut ){  
      *this = tmp;
      nt--;
    }
  }   
}
