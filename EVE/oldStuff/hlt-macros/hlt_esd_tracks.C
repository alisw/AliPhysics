// $Id: esd_tracks.C 35148 2009-10-01 11:21:06Z mtadel $
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TString.h"
#include "TMath.h"
#include "TGListTree.h"
#include "TEveVSDStructs.h"
#include "TEveManager.h"
#include "TEveTrackPropagator.h"

#include "AliESDEvent.h"
#include "AliESDtrackCuts.h"
#include "AliESDtrack.h"
#include "AliExternalTrackParam.h"

#include "EVE/EveBase/AliEveTrack.h"
#include "EVE/EveBase/AliEveMagField.h"
#include "EVE/EveBase/AliEveEventManager.h"

#include "AliHLTTPCCATrackParam.h"
#include "AliHLTTPCCATrackConvertor.h"

#endif

AliEveTrack* hlt_esd_make_track(AliESDtrack *at, TEveTrackList* cont)
{
  // Make a standard track representation and put it into given container.

  // Choose which parameters to use a track's starting point.
  // If gkFixFailedITSExtr is TRUE (FALSE by default) and
  // if ITS refit failed, take track parameters at inner TPC radius.

  const double kCLight = 0.000299792458;
  double bz = - kCLight*10.*( cont->GetPropagator()->GetMagField(0,0,0).fZ);

  Bool_t innerTaken = kFALSE;
  if ( ! at->IsOn(AliESDtrack::kITSrefit) && g_esd_tracks_use_ip_on_failed_its_refit)
  {
    //tp = at->GetInnerParam();
    innerTaken = kTRUE;
  }

  // Add inner/outer track parameters as path-marks.

  Double_t     pbuf[3], vbuf[3];

  AliExternalTrackParam trackParam = *at;

  // take parameters constrained to vertex (if they are)

  if( at->GetConstrainedParam() ){
    trackParam = *at->GetConstrainedParam();
  }
  else if( at->GetInnerParam() ){
    trackParam = *(at->GetInnerParam());
  }
  if( at->GetStatus()&AliESDtrack::kTRDin ){
    // transport to TRD in
    trackParam = *at;
    trackParam.PropagateTo( 290.45, -10.*( cont->GetPropagator()->GetMagField(0,0,0).fZ) );
  }

  TEveRecTrack rt;
  {
    rt.fLabel  = at->GetLabel();
    rt.fIndex  = (Int_t) at->GetID();
    rt.fStatus = (Int_t) at->GetStatus();
    rt.fSign   = (Int_t) trackParam.GetSign();  
    trackParam.GetXYZ(vbuf);
    trackParam.GetPxPyPz(pbuf);    
    rt.fV.Set(vbuf);
    rt.fP.Set(pbuf);
    Double_t ep = at->GetP(), mc = at->GetMass();
    rt.fBeta = ep/TMath::Sqrt(ep*ep + mc*mc);
  }

  AliEveTrack* track = new AliEveTrack(&rt, cont->GetPropagator());
  track->SetAttLineAttMarker(cont);
  track->SetName(Form("AliEveTrack %d", at->GetID()));
  track->SetElementTitle(esd_track_title(at));
  track->SetSourceObject(at);


  // Set reference points along the trajectory
  // and the last point

  { 
    TEvePathMark startPoint(TEvePathMark::kReference);
    trackParam.GetXYZ(vbuf);
    trackParam.GetPxPyPz(pbuf);    
    startPoint.fV.Set(vbuf);
    startPoint.fP.Set(pbuf);
    rt.fV.Set(vbuf);
    rt.fP.Set(pbuf);
    Double_t ep = at->GetP(), mc = at->GetMass();
    rt.fBeta = ep/TMath::Sqrt(ep*ep + mc*mc);

    track->AddPathMark( startPoint );    
  }


  if( at->GetTPCPoints(2)>80 ){
  
    //
    // use AliHLTTPCCATrackParam propagator 
    // since AliExternalTrackParam:PropagateTo()
    // has an offset at big distances
    //
    
    AliHLTTPCCATrackParam t;
    AliHLTTPCCATrackConvertor::SetExtParam( t, trackParam );
    
    Double_t x0 = trackParam.GetX();
    Double_t dx = at->GetTPCPoints(2) - x0;
    
    //
    // set a reference at the half of trajectory for better drawing
    //
    
    for( double dxx=dx/2; TMath::Abs(dxx)>=1.; dxx*=.9 ){
      if( !t.TransportToX(x0+dxx, bz, .99 ) ) continue;
      AliHLTTPCCATrackConvertor::GetExtParam( t, trackParam, trackParam.GetAlpha() ); 
      trackParam.GetXYZ(vbuf);
      trackParam.GetPxPyPz(pbuf);
      TEvePathMark midPoint(TEvePathMark::kReference);
      midPoint.fV.Set(vbuf);
      midPoint.fP.Set(pbuf);    
      track->AddPathMark( midPoint );
      break;
    }
    
    //
    // Set a reference at the end of the trajectory
    // and a "decay point", to let the event display know where the track ends
    //
    
    for( ; TMath::Abs(dx)>=1.; dx*=.9 ){
      if( !t.TransportToX(x0+dx, bz, .99 ) ) continue;
      AliHLTTPCCATrackConvertor::GetExtParam( t, trackParam, trackParam.GetAlpha() ); 
      trackParam.GetXYZ(vbuf);
      trackParam.GetPxPyPz(pbuf);
      TEvePathMark endPoint(TEvePathMark::kReference);
      TEvePathMark decPoint(TEvePathMark::kDecay);
      endPoint.fV.Set(vbuf);
      endPoint.fP.Set(pbuf);
      decPoint.fV.Set(vbuf);
      decPoint.fP.Set(pbuf);
      track->AddPathMark( endPoint );
      track->AddPathMark( decPoint );
      break;
    }  
  }

  if (at->IsOn(AliESDtrack::kTPCrefit))
  {
    if ( ! innerTaken)
    {
      esd_track_add_param(track, at->GetInnerParam());
    }
    esd_track_add_param(track, at->GetOuterParam());
  }


  return track;
}

