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

/*****************************************************************
  AliStarTrack: Event container for flow analysis

  origin:   Mikolaj Krzewicki  (mikolaj.krzewicki@cern.ch)
*****************************************************************/

#include <string.h>
#include <TObject.h>
#include <TString.h>
#include "AliStarTrack.h"

ClassImp(AliStarTrack)

//______________________________________________________________________________
AliStarTrack::AliStarTrack():
  TObject(),
  fParams()
{
  //ctor
}

//______________________________________________________________________________
AliStarTrack::AliStarTrack( const Float_t* params ):
  TObject(),
  fParams()
{
  //ctor
  memcpy(fParams,params,fgkNparams*sizeof(Float_t));
}

//______________________________________________________________________________
AliStarTrack::AliStarTrack( const AliStarTrack& track ):
  TObject(),
  fParams()
{
  //copy ctor
  memcpy(fParams,track.fParams,fgkNparams*sizeof(Float_t));
}

//______________________________________________________________________________
AliStarTrack& AliStarTrack::operator=( const AliStarTrack& track )
{
  //assignment
  if (this == &track) return *this;
  TObject::operator=(track);
  memcpy(fParams,track.fParams,fgkNparams*sizeof(Float_t));
  return *this;
}

//______________________________________________________________________________
void AliStarTrack::SetParams( const Float_t* params )
{
  //set params
  memcpy(fParams,params,fgkNparams*sizeof(Float_t));
}

//______________________________________________________________________________
AliStarTrack* AliStarTrack::Clone( const char* /*option*/) const
{
  //clone "constructor"
  return new AliStarTrack(*this);
}

//______________________________________________________________________________
void AliStarTrack::Print( Option_t* option ) const
{
  // print info
  // TNtuple* track: names are documented in the next 2 lines
  // tracks = new TNtuple("tracks","tracks",
  //   "ID:Charge:Eta:Phi:Pt:Dca:nHits:nHitsFit:nHitsPoss:nHitsDedx:dEdx:nSigElect:nSigPi:nSigK:nSigProton" ) ;
  //
  TString optionstr(option);
  if ( optionstr.Contains("legend") )
  {
    printf(
	  "    id charge     eta     phi      pt     dca  nHits  nFit nPoss ndEdx   dEdx nSigElec nSigPi  nSigK nSigPr\n") ;
    return;
  }

  Int_t   id             = (Int_t)   fParams[0]   ;  // id - a unique integer for each track in this event 
  Int_t   charge         = (Int_t)   fParams[1]   ;  // +1 or -1 
  Float_t eta            = (Float_t) fParams[2]   ;  // Pseudo-rapidity at the vertex
  Float_t phi            = (Float_t) fParams[3]   ;  // Phi emission angle at the vertexcd 
  Float_t pt             = (Float_t) fParams[4]   ;  // Pt of the track at the vertex
  Float_t dca            = (Float_t) fParams[5]   ;  // Magnitude of 3D DCA to vertex
  Int_t   nHits          = (Int_t)   fParams[6]   ;  // Number of clusters available to the Kalman Filter
  Int_t   nHitsFit       = (Int_t)   fParams[7]   ;  // Number of clusters used in the fit (after cuts)
  Int_t   nHitsPoss      = (Int_t)   fParams[8]   ;  // Number of possible cluster on track (theoretical max) 
  Int_t   nHitsDedx      = (Int_t)   fParams[9]   ;  // Number of clusters used in the fit (after dEdx cuts)
  Float_t dEdx           = 1.e6*(Float_t)fParams[10]  ;  // Measured dEdx (Note: GeV/cm so convert to keV/cm!!)
  Float_t nSigmaElectron = (Float_t) fParams[11]  ;  // Number of sigma from electron Bethe-Bloch curve
  Float_t nSigmaPion     = (Float_t) fParams[12]  ;  // Number of sigma from pion Bethe-Bloch curve
  Float_t nSigmaKaon     = (Float_t) fParams[13]  ;  // Number of sigma from kaon Bethe-Bloch curve
  Float_t nSigmaProton   = (Float_t) fParams[14]  ;  // Number of sigma from proton Bethe-Bloch curve
  
  // Alternative way to access the data 
  // nHitsPoss      = (Int_t) ( fTracks->GetLeaf("nHitsPoss")->GetValue() ) ;  // Note alternative method to retrieve data
  // Using the definition of the original NTuple
  // TrackTuple      = new TNtuple("NTtracks","NTtracks",
  // "ID:Charge:Eta:Phi:Pt:Dca:nHits:nHitsFit:nHitsPoss:nHitsDedx:dEdx:nSigElect:nSigPi:nSigK:nSigProton" ) 
  
  printf("%6d %4d   %7.3f %7.3f %7.3f %7.4f %6d %5d %5d %5d %6.2f   %6.2f %6.2f %6.2f %6.2f \n",
	 id, charge, eta, phi, pt, dca, nHits, nHitsFit, nHitsPoss, nHitsDedx, dEdx,
	 nSigmaElectron, nSigmaPion, nSigmaKaon, nSigmaProton ) ;
}

//______________________________________________________________________________
Int_t AliStarTrack::PID() const
{
  // Note: This is a very simple PID selection scheme.  More elaborate methods (with multiple cuts) may be required.
  // When you *are* using dEdx information, you must chose a finite number of good Dedx hits ... but the limit should
  // be about 2/3 of nHitsMin.  This is because some clusters do not form good dEdx hits due to track
  // merging, etc., and so nHitsDedx is always less than nHitsFit.  A rule of thumb says ~2/3 ratio.

  Int_t pid = 0 ;

  const Int_t   nHitDedxMin =    15  ;       // 10 to 20 is typical.  nHitDedxMin is often chosen to be about 2/3 of nHitMin.
  const Float_t nSigmaPID   =    2.0 ;       // Number of Sigma cut to apply to PID bands

  // Test on Number of dE/dx hits required, return 0 if not enough hits
  if ( GetNHitsDedx() <  nHitDedxMin ) return pid;

  // Begin PID

  if ( TMath::Abs( GetNSigElect() ) >= nSigmaPID )
  {
    if ( TMath::Abs( GetNSigK()  ) <= nSigmaPID )
    {
      pid = 321  ;
    }
    if ( TMath::Abs( GetNSigProton()  ) <= nSigmaPID )
    {
      pid = 2212 ;
    }
    if ( TMath::Abs( GetNSigPi()  ) <= nSigmaPID )
    {
      pid = 211  ;
    }
  }

  // Pion is the default in case of ambiguity because it is most abundent. Don't re-arrange order, above.

  return pid ;
}

