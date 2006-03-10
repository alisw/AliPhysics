/*******************************************************************************
 * Copyright(c) 2003, IceCube Experiment at the South Pole. All rights reserved.
 *
 * Author: The IceCube RALICE-based Offline Project.
 * Contributors are mentioned in the code where appropriate.
 *
 * Permission to use, copy, modify and distribute this software and its
 * documentation strictly for non-commercial purposes is hereby granted
 * without fee, provided that the above copyright notice appears in all
 * copies and that both the copyright notice and this permission notice
 * appear in the supporting documentation.
 * The authors make no claims about the suitability of this software for
 * any purpose. It is provided "as is" without express or implied warranty.
 *******************************************************************************/

// $Id$

///////////////////////////////////////////////////////////////////////////
// Class IceLinefit
// TTask derived class to perform a linefit track reconstruction.
// The procedure is based on the method described in the Amanda publication
// in Nuclear Instruments and Methods A524 (2004) 179-180.
// To prevent waisting CPU time in trying to reconstruct (high-energy) cascade
// events, or to select specifically reconstruction of low multiplicity events,
// the user may invoke the memberfunctions SetMaxModA() and SetMinModA.
// This allows selection of events for processing with a certain maximum
// and/or minimum number of good Amanda OMs firing.
// By default the minimum and maximum are set to 0 and 999, respectively,
// in the constructor, which implies no multiplicity selection. 
//
// The reconstructed track is stored in the IceEvent structure with as
// default "IceLinefit" as the name of the track.
// This track name identifier can be modified by the user via the
// SetTrackName() memberfunction. This will allow unique identification
// of tracks which are produced when re-processing existing data with
// different criteria.
// The track 3-momentum is set to the reconstructed velocity, normalised
// to 1 GeV. The mass and charge of the track are left 0.
// The r0 and t0 can be obtained from the reference point of the track,
// whereas the t0 ia also available from the track timestamp .
//
// For further details the user is referred to NIM A524 (2004) 169.
//
// Note : This algorithm works best on data which has been calibrated
//        (IceCalibrate), cross talk corrected (IceXtalk) and cleaned
//        from noise hits etc. (IceCleanHits).
//
//--- Author: Nick van Eijndhoven 10-mar-2006 Utrecht University
//- Modified: NvE $Date$ Utrecht University
///////////////////////////////////////////////////////////////////////////
 
#include "IceLinefit.h"
#include "Riostream.h"

ClassImp(IceLinefit) // Class implementation to enable ROOT I/O

IceLinefit::IceLinefit(const char* name,const char* title) : TTask(name,title)
{
// Default constructor.
 fMaxmodA=999;
 fMinmodA=0;
 fTrackname="IceLinefit";
}
///////////////////////////////////////////////////////////////////////////
IceLinefit::~IceLinefit()
{
// Default destructor.
}
///////////////////////////////////////////////////////////////////////////
void IceLinefit::SetMaxModA(Int_t nmax)
{
// Set the maximum number of good Amanda modules that may have fired
// in order to process this event.
// This allows suppression of processing (high-energy) cascade events
// with this linefit tracking to prevent waisting cpu time for cases
// in which tracking doesn't make sense anyhow.
// Furthermore it allows selection of low multiplicity events for processing.
// By default the maximum number of Amanda modules is set to 999 in the ctor,
// which implies no selection on maximum module multiplicity.
// See also the memberfunction SetMinModA().
 fMaxmodA=nmax;
}
///////////////////////////////////////////////////////////////////////////
void IceLinefit::SetMinModA(Int_t nmin)
{
// Set the minimum number of good Amanda modules that must have fired
// in order to process this event.
// This allows selection of a minimal multiplicity for events to be processed.
// By default the minimum number of Amanda modules is set to 0 in the ctor,
// which implies no selection on minimum module multiplicity.
// See also the memberfunction SetMaxModA().
 fMinmodA=nmin;
}
///////////////////////////////////////////////////////////////////////////
void IceLinefit::SetTrackName(TString s)
{
// Set (alternative) name identifier for the produced first guess tracks.
// This allows unique identification of (newly) produced linefit tracks
// in case of re-processing of existing data with different criteria.
// By default the produced first guess tracks have the name "IceLinefit"
// which is set in the constructor of this class.
 fTrackname=s;
}
///////////////////////////////////////////////////////////////////////////
void IceLinefit::Exec(Option_t* opt)
{
// Implementation of the linefit reconstruction.

 TString name=opt;
 AliJob* parent=(AliJob*)(gROOT->GetListOfTasks()->FindObject(name.Data()));

 if (!parent) return;

 IceEvent* evt=(IceEvent*)parent->GetObject("IceEvent");
 if (!evt) return;

 // Fetch all fired Amanda OMs for this event
 TObjArray* aoms=evt->GetDevices("IceAOM");
 Int_t naoms=aoms->GetEntries();
 if (!naoms) return;

 // Check for the minimum and/or maximum number of good fired Amanda OMs
 Int_t ngood=0;
 for (Int_t iom=0; iom<naoms; iom++)
 {
  IceGOM* omx=(IceGOM*)aoms->At(iom);
  if (!omx) continue;
  if (omx->GetDeadValue("ADC") || omx->GetDeadValue("LE") || omx->GetDeadValue("TOT")) continue;
  ngood++;
 } 
 if (ngood<fMinmodA || ngood>fMaxmodA) return;

 AliSignal* sx=0;
 Ali3Vector rom,sumr;
 Ali3Vector rt,sumrt;
 Float_t thit;
 Float_t sumt=0,sumt2=0;
 TObjArray hits;

 // Loop over all OMs and hits to determine the linefit parameters.
 // Also all the used hits are recorded for association with the track.
 for (Int_t iom=0; iom<naoms; iom++)
 {
  IceGOM* omx=(IceGOM*)aoms->At(iom);
  if (!omx) continue;
  if (omx->GetDeadValue("LE")) continue;
  rom=(Ali3Vector)omx->GetPosition();
  // Use all the good hits of this OM
  for (Int_t ih=1; ih<=omx->GetNhits(); ih++)
  {
   sx=omx->GetHit(ih);
   if (!sx) continue;
   if (sx->GetDeadValue("ADC") || sx->GetDeadValue("LE") || sx->GetDeadValue("TOT")) continue;

   thit=sx->GetSignal("LE",7);
   rt=rom*thit;
   sumr+=rom;
   sumrt+=rt;
   sumt+=thit;
   sumt2+=thit*thit;

   // Record this hit for association with the track
   hits.Add(sx);
  }
 }

 Int_t nused=hits.GetEntries();
 if (!nused) return;

 sumr/=float(nused);
 sumrt/=float(nused);
 sumt/=float(nused);
 sumt2/=float(nused);

 Ali3Vector v;
 Ali3Vector temp;
 temp=sumr*sumt;
 v=sumrt-temp;
 Float_t dum=sumt2-(sumt*sumt);
 if (dum) v/=dum;

 Ali3Vector r;
 temp=v*sumt;
 r=sumr-temp;

 AliTrack t; 
 t.SetNameTitle(fTrackname.Data(),"IceLinefit linefit track");
 evt->AddTrack(t);
 AliTrack* trk=evt->GetTrack(evt->GetNtracks());
 if (!trk) return;

 Ali3Vector p;
 Float_t vec[3];
 v.GetVector(vec,"sph");
 vec[0]=1;
 p.SetVector(vec,"sph");

 AliPosition r0;
 r0.SetPosition(r);
 r0.SetTimestamp((AliTimestamp&)*evt);
 AliTimestamp* t0=r0.GetTimestamp();
 t0->Add(0,0,(int)sumt);

 trk->Set3Momentum(p);
 trk->SetReferencePoint(r0);
 trk->SetTimestamp(*t0);

 // Link the used hits to the track (and vice versa)
 for (Int_t i=0; i<nused; i++)
 {
  sx=(AliSignal*)hits.At(i);
  if (sx) sx->AddLink(trk);
 } 
}
///////////////////////////////////////////////////////////////////////////
