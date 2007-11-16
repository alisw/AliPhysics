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
//
// In case an event has been rejected by an AliEventSelector (based) processor,
// this task (and its sub-tasks) is not executed.
//
// Note : Amanda OMs and InIce DOMs are treated seperately, which means that
//        for events with both OMs and DOMs firing, 2 linefit tracks will
//        be produced. The 2 linefit tracks can be distinguished on basis of
//        their name as explained below. 
//
// The procedure is based on the method described in the Amanda publication
// in Nuclear Instruments and Methods A524 (2004) 179-180.
// To prevent waisting CPU time in trying to reconstruct (high-energy) cascade
// events, or to select specifically reconstruction of low multiplicity events,
// the user may invoke the memberfunctions SetMaxMod() and SetMinMod.
// This allows selection of events for processing with a certain maximum
// and/or minimum number of good (D)OMs firing.
// By default the minimum and maximum are set to 0 and 999999, respectively,
// in the constructor, which implies no multiplicity selection. 
// The maximum number of good hits per (D)OM to be used for the reconstruction
// can be specified via the memberfunction SetMaxHits().
// By default all good hits of each (D)OM are used but the user may want
// to restrict this number to the first n hits of each (D)OM to account
// for possible noise and/or afterpulse signals that are not recognised by the
// hit cleaning procedure.
//
// Information about the actual parameter settings can be found in the event
// structure itself via the device named "IceLinefit".
//
// The reconstructed track is stored in the IceEvent structure with as
// default the Classname of the producing processor as the name of the track.
// A suffix "A" for an Amanda (OM) track or a suffix "I" for an InIce (DOM) track
// will be added to the name automatically.
// This track name identifier can be modified by the user via the
// SetTrackName() memberfunction. This will allow unique identification
// of tracks which are produced when re-processing existing data with
// different criteria.
// Note that a suffix "A" or "I" will always be generated automatically.
// The track 3-momentum is set to the reconstructed velocity, normalised
// to 1 GeV. The mass and charge of the track are left 0.
// The r0 and t0 can be obtained from the reference point of the track,
// whereas the t0 ia also available from the track timestamp .
// By default the charge of the produced tracks is set to 0, since
// no distinction can be made between positive or negative tracks.
// However, the user can define the track charge by invokation
// of the memberfunction SetCharge().
// This facility may be used to distinguish tracks produced by the
// various reconstruction algorithms in a (3D) colour display
// (see the class AliHelix for further details).
// The value of beta=v/c for the reconstructed velocity is available
// from the fitdetails as stored for the reconstructed track. 
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
 fEvt=0;
 fMaxmodA=999999;
 fMinmodA=0;
 fMaxhitsA=0;
 fMaxmodI=999999;
 fMinmodI=0;
 fMaxhitsI=0;
 fTrackname="";
 fCharge=0;
}
///////////////////////////////////////////////////////////////////////////
IceLinefit::~IceLinefit()
{
// Default destructor.
}
///////////////////////////////////////////////////////////////////////////
void IceLinefit::SetMaxMod(Int_t nmax,TString s)
{
// Set the maximum number of good (D)OMs that may have fired
// in order to process this event.
// This allows suppression of processing (high-energy) cascade events
// with this linefit tracking to prevent waisting cpu time for cases
// in which tracking doesn't make sense anyhow.
// Furthermore it allows selection of low multiplicity events for processing.
// By default the maximum number of (D)OMs is set to 999999 in the ctor,
// which implies no selection on maximum module multiplicity.
// See also the memberfunction SetMinMod().
//
// The input argument "s" allows for detector specification.
//
// s = "A" --> Amanda OMs
//     "I" --> InIce DOMs
//
// The default is s="A" for backward compatibility.

 if (s=="A") fMaxmodA=nmax;
 if (s=="I") fMaxmodI=nmax;
}
///////////////////////////////////////////////////////////////////////////
void IceLinefit::SetMinMod(Int_t nmin,TString s)
{
// Set the minimum number of good (D)OMs that must have fired
// in order to process this event.
// This allows selection of a minimal multiplicity for events to be processed.
// By default the minimum number of (D)OMs is set to 0 in the ctor,
// which implies no selection on minimum module multiplicity.
// See also the memberfunction SetMaxMod().
//
// The input argument "s" allows for detector specification.
//
// s = "A" --> Amanda OMs
//     "I" --> InIce DOMs
//
// The default is s="A" for backward compatibility.

 if (s=="A") fMinmodA=nmin;
 if (s=="I") fMinmodI=nmin;
}
///////////////////////////////////////////////////////////////////////////
void IceLinefit::SetMaxHits(Int_t nmax,TString s)
{
// Set the maximum number of good hits per Amanda module to be processed.
//
// Special values :
// nmax = 0 : No maximum limit set; all good hits will be processed
//      < 0 : No hits will be processed
//
// In case the user selects a maximum number of good hits per module, all the
// hits of each module will be ordered w.r.t. increasing hit time (LE).
// This allows selection of processing e.g. only the first good hits etc...
// By default the maximum number of hits per (D)OM is set to 0 in the ctor,
// which implies just processing all good hits without any maximum limit.
//
// The input argument "s" allows for detector specification.
//
// s = "A" --> Amanda OMs
//     "I" --> InIce DOMs
//
// The default is s="A" for backward compatibility.

 if (s=="A") fMaxhitsA=nmax;
 if (s=="I") fMaxhitsI=nmax;
}
///////////////////////////////////////////////////////////////////////////
void IceLinefit::SetTrackName(TString s)
{
// Set (alternative) name identifier for the produced first guess tracks.
// This allows unique identification of (newly) produced linefit tracks
// in case of re-processing of existing data with different criteria.
// By default the produced first guess tracks have the name of the class
// by which they were produced.
 fTrackname=s;
}
///////////////////////////////////////////////////////////////////////////
void IceLinefit::SetCharge(Float_t charge)
{
// Set user defined charge for the produced first guess tracks.
// This allows identification of these tracks on color displays.
// By default the produced first guess tracks have charge=0
// which is set in the constructor of this class.
 fCharge=charge;
}
///////////////////////////////////////////////////////////////////////////
void IceLinefit::Exec(Option_t* opt)
{
// Implementation of the linefit reconstruction.

 TString name=opt;
 AliJob* parent=(AliJob*)(gROOT->GetListOfTasks()->FindObject(name.Data()));

 if (!parent) return;

 fEvt=(IceEvent*)parent->GetObject("IceEvent");
 if (!fEvt) return;

 // Only process accepted events
 AliDevice* seldev=(AliDevice*)fEvt->GetDevice("AliEventSelector");
 if (seldev)
 {
  if (seldev->GetSignal("Select") < 0.1) return;
 }

 // Enter the reco parameters as a device in the event
 AliSignal params;
 params.SetNameTitle("IceLinefit","IceLinefit reco parameters");
 params.AddNamedSlot("MaxmodA");
 params.AddNamedSlot("MinmodA");
 params.AddNamedSlot("MaxhitsA");
 params.AddNamedSlot("MaxmodI");
 params.AddNamedSlot("MinmodI");
 params.AddNamedSlot("MaxhitsI");

 params.SetSignal(fMaxmodA,"MaxmodA");
 params.SetSignal(fMinmodA,"MinmodA");
 params.SetSignal(fMaxhitsA,"MaxhitsA");
 params.SetSignal(fMaxmodI,"MaxmodI");
 params.SetSignal(fMinmodI,"MinmodI");
 params.SetSignal(fMaxhitsI,"MaxhitsI");

 fEvt->AddDevice(params);

 // Perform linefit reconstruction for the various hits
 Amanda();
 InIce();
}
///////////////////////////////////////////////////////////////////////////
void IceLinefit::Amanda()
{
// Implementation of the linefit reconstruction for Amanda OMs.

 if (fMaxhitsA<0) return;

 // Fetch all fired Amanda OMs for this event
 TObjArray* aoms=fEvt->GetDevices("IceAOM");
 if (!aoms) return;
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

 const Float_t c=0.299792; // Light speed in vacuum in meters per ns

 AliSignal* sx=0;
 Ali3Vector rom,sumr;
 Ali3Vector rt,sumrt;
 Float_t thit;
 Float_t sumt=0,sumt2=0;
 TObjArray hits;
 TObjArray* ordered;
 Int_t nh;

 // Loop over all OMs and hits to determine the linefit parameters.
 // Also all the used hits are recorded for association with the track.
 for (Int_t iom=0; iom<naoms; iom++)
 {
  IceGOM* omx=(IceGOM*)aoms->At(iom);
  if (!omx) continue;
  if (omx->GetDeadValue("LE")) continue;
  rom=(Ali3Vector)omx->GetPosition();
  // Use the specified good hits of this OM
  ordered=0;
  if (fMaxhitsA>0 && omx->GetNhits()>fMaxhitsA) ordered=omx->SortHits("LE",1,0,7);
  nh=0;
  for (Int_t ih=1; ih<=omx->GetNhits(); ih++)
  {
   if (ordered)
   {
    if (nh>=fMaxhitsA) break;
    sx=(AliSignal*)ordered->At(ih-1);
   }
   else
   {
    sx=omx->GetHit(ih);
   }
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
   nh++;
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

 Float_t beta=v.GetNorm()/c;
 AliSignal fitstats;
 fitstats.SetNameTitle("Fitstats","Fit stats for IceLinefit");
 fitstats.SetSlotName("Beta",1);
 fitstats.SetSignal(beta,1);

 Ali3Vector r;
 temp=v*sumt;
 r=sumr-temp;

 TString name=fTrackname;
 if (name=="") name=ClassName();
 name+="A";
 TString title=ClassName();
 title+=" Amanda track";
 AliTrack t; 
 t.SetNameTitle(name.Data(),title.Data());
 t.SetCharge(fCharge);
 fEvt->AddTrack(t);
 AliTrack* trk=fEvt->GetTrack(fEvt->GetNtracks());
 if (!trk) return;

 trk->SetId(fEvt->GetNtracks(1)+1);

 Ali3Vector p;
 Float_t vec[3];
 v.GetVector(vec,"sph");
 vec[0]=1;
 p.SetVector(vec,"sph");

 AliPosition r0;
 r0.SetPosition(r);
 r0.SetTimestamp((AliTimestamp&)*fEvt);
 AliTimestamp* t0=r0.GetTimestamp();
 t0->Add(0,0,(int)sumt);

 trk->Set3Momentum(p);
 trk->SetReferencePoint(r0);
 trk->SetTimestamp(*t0);
 trk->SetFitDetails(fitstats);

 // Link the used hits to the track (and vice versa)
 for (Int_t i=0; i<nused; i++)
 {
  sx=(AliSignal*)hits.At(i);
  if (sx) sx->AddTrack(*trk);
 } 
}
///////////////////////////////////////////////////////////////////////////
void IceLinefit::InIce()
{
// Implementation of the linefit reconstruction for InIce DOMs.

 if (fMaxhitsI<0) return;

 // Fetch all fired InIce DOMs for this event
 TObjArray* idoms=fEvt->GetDevices("IceIDOM");
 if (!idoms) return;
 Int_t nidoms=idoms->GetEntries();
 if (!nidoms) return;

 // Check for the minimum and/or maximum number of good fired InIce DOMs
 Int_t ngood=0;
 for (Int_t idom=0; idom<nidoms; idom++)
 {
  IceGOM* omx=(IceGOM*)idoms->At(idom);
  if (!omx) continue;
  if (omx->GetDeadValue("ADC") || omx->GetDeadValue("LE") || omx->GetDeadValue("TOT")) continue;
  ngood++;
 } 
 if (ngood<fMinmodI || ngood>fMaxmodI) return;

 const Float_t c=0.299792; // Light speed in vacuum in meters per ns

 AliSignal* sx=0;
 Ali3Vector rom,sumr;
 Ali3Vector rt,sumrt;
 Float_t thit;
 Float_t sumt=0,sumt2=0;
 TObjArray hits;
 TObjArray* ordered;
 Int_t nh;

 // Loop over all DOMs and hits to determine the linefit parameters.
 // Also all the used hits are recorded for association with the track.
 for (Int_t idom=0; idom<nidoms; idom++)
 {
  IceGOM* omx=(IceGOM*)idoms->At(idom);
  if (!omx) continue;
  if (omx->GetDeadValue("LE")) continue;
  rom=(Ali3Vector)omx->GetPosition();
  // Use the specified good hits of this DOM
  ordered=0;
  if (fMaxhitsI>0 && omx->GetNhits()>fMaxhitsI) ordered=omx->SortHits("LE",1,0,7);
  nh=0;
  for (Int_t ih=1; ih<=omx->GetNhits(); ih++)
  {
   if (ordered)
   {
    if (nh>=fMaxhitsI) break;
    sx=(AliSignal*)ordered->At(ih-1);
   }
   else
   {
    sx=omx->GetHit(ih);
   }
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
   nh++;
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

 Float_t beta=v.GetNorm()/c;
 AliSignal fitstats;
 fitstats.SetNameTitle("Fitstats","Fit stats for IceLinefit");
 fitstats.SetSlotName("Beta",1);
 fitstats.SetSignal(beta,1);

 Ali3Vector r;
 temp=v*sumt;
 r=sumr-temp;

 TString name=fTrackname;
 if (name=="") name=ClassName();
 name+="I";
 TString title=ClassName();
 title+=" InIce track";
 AliTrack t; 
 t.SetNameTitle(name.Data(),title.Data());
 t.SetCharge(fCharge);
 fEvt->AddTrack(t);
 AliTrack* trk=fEvt->GetTrack(fEvt->GetNtracks());
 if (!trk) return;

 trk->SetId(fEvt->GetNtracks(1)+1);

 Ali3Vector p;
 Float_t vec[3];
 v.GetVector(vec,"sph");
 vec[0]=1;
 p.SetVector(vec,"sph");

 AliPosition r0;
 r0.SetPosition(r);
 r0.SetTimestamp((AliTimestamp&)*fEvt);
 AliTimestamp* t0=r0.GetTimestamp();
 t0->Add(0,0,(int)sumt);

 trk->Set3Momentum(p);
 trk->SetReferencePoint(r0);
 trk->SetTimestamp(*t0);
 trk->SetFitDetails(fitstats);

 // Link the used hits to the track (and vice versa)
 for (Int_t i=0; i<nused; i++)
 {
  sx=(AliSignal*)hits.At(i);
  if (sx) sx->AddTrack(*trk);
 } 
}
///////////////////////////////////////////////////////////////////////////
