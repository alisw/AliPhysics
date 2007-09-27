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
// Class IceDwalk
// TTask derived class to perform direct walk track reconstruction.
//
// In case an event has been rejected by an AliEventSelector (based) processor,
// this task (and its sub-tasks) is not executed.
//
// The procedure is based on the method described in the Amanda publication
// in Nuclear Instruments and Methods A524 (2004) 179-180.
// However, the Amanda method has been extended with the intention to
// take also multiple (muon) tracks within 1 event into account.
// This will not only provide a means to reconstruct muon bundles and
// multiple track events in IceCube, but will also allow to reduce the
// background of faked upgoing muons as a result of multiple downgoing
// muons hitting the top and bottom parts of the detector.
// A further extension of the original Amanda method is the separate treatment
// of the phase and group velocities as introduced in collaboration with
// George Japaridze (Clark Atlanta University, USA) which will provide more
// accurate time residuals due to the different velocities of the Cerenkov
// wave front (v_phase) and the actually detected photons (v_group).
// This distinction between v_phase and v_group can be (de)activated via the
// memberfunction SetVgroupUsage(). By default the distinction between v_phase
// and v_group is activated in the constructor of this class.
// To prevent waisting CPU time in trying to reconstruct (high-energy) cascade
// events, or to select specifically reconstruction of low multiplicity events,
// the user may invoke the memberfunctions SetMaxModA() and SetMinModA().
// This allows selection of events for processing with a certain maximum and/or
// minimum number of good Amanda OMs firing.
// By default the minimum and maximum are set to 0 and 999, respectively,
// in the constructor, which implies no multiplicity selection. 
// The maximum number of good hits per Amanda OM to be used for the reconstruction
// can be specified via the memberfunction SetMaxHitsA().
// By default only the first good hit of each Amanda OM is used.
// Note that when all the good hits of an OM are used, this may lead to large
// processing time in case many noise and/or afterpulse signals are not
// recognised by the hit cleaning procedure.
//
// Information about the actual parameter settings can be found in the event
// structure itself via the device named "IceDwalk".
//
// The various reconstruction steps are summarised as follows :
//
// 1) Construction of track elements (TE's).
//    A track element is a straight line connecting two hits that
//    appeared at some minimum distance d and within some maximum
//    time difference dt, according to eq. (20) of the NIM article.
//    The default value for d is 75 meter, but this can be modified
//    via the memberfunction SetDmin().
//    By default dt=(hit distance)/c but an additional time margin
//    may be specified via the memberfunction SetDtmarg().    
//    The reference point r0 of the TE is taken as the center between
//    the two hit positions and the TE timestamp t0 at the position r0
//    is taken as the IceEvent timestamp increased by the average of the
//    two hit times. So, all timestamps contain the overall IceEvent
//    timestamp as a basis. This means that time differences can be
//    obtained via the AliTimestamp facilities (supporting upto picosecond
//    precision when available).
//    The TE direction is given by the relative position of the two hits.
//
// 2) Each TE will obtain so called associated hits.
//    A hit is associated to a TE when it fulfills both the conditions
//
//      -30 < tres < 300 ns
//      dhit/lambda < F
//
//    tres   : time residual
//             Difference between the observed hit time and the time expected
//             for a direct photon hit.     
//    dhit   : Distance traveled by the cherenkov photon from the track to the hit position
//    lambda : Photon scattering length in ice
//
//    By default F is set to 3.07126, but this can be modified via the memberfunction
//    SetMaxDhit(). 
//
// 3) Construction of track candidates (TC's).
//    These are TE's that fulfill both the conditions
//
//     nah >= 1
//     qtc >= 0.8*qtcmax
//
//    where we have defined :
//
//    nah    : Number of associated hits for the specific TE.
//    qtc    : The track quality number (see hereafter).
//    qtcmax : Maximum quality number encountered for the TE's.
//
//    The track quality number qtc is defined as follows :
//
//     qtc=nah*(term1+term2)-term3-term4-term5
//
//    here we have defined :
//
//    term1=2*spread/span
//    term2=2*spreadL/spanL
//    term3=|spread-expspread|/spread
//    term4=|spreadL-expspreadL|/spreadL
//    term5=|medianT|/spreadT
//
//    The central observables here are the projected positions X on the track
//    of the various associated hits w.r.t. the track reference point r0.
//    Note that X can be negative as well as positive.
//    Therefore we also introduce XL=|X|. 
//
//    span       : max(X)-min(X)
//    spanL      : max(XL)-min(XL)
//    Xmedian    : median of X
//    XmedianL   : median of XL
//    spread     : < |X-Xmedian| >
//    spreadL    : < |XL-XmedianL| >
//    expspread  : expected spread in X for a flat distribution of nah hits over span
//    expspreadL : expected spread in XL for a flat distribution of nah hits over spanL
//    medianT    : median of tres
//    spreadT    : < |tres-medianT| >
//
//    However, if |Xmedian| > span/2 we set qtc=0 in order to always require
//    projected hits to appear on both sides of r0 on the track.
//
//    Note : The qtc quality number is used to define the norm of the momentum
//           of the track candidate. As such it serves as a weight for the jet
//           momentum (direction) after clustering of the TC's and lateron
//           merging of the jets (see hereafter).
//
// 4) The remaining track candidates are clustered into jets when their directions
//    are within a certain maximum opening angle.
//    In addition a track candidate must within a certain maximum distance
//    of the jet starting TC in order to get clustered. 
//    The latter criterion prevents clustering of (nearly) parallel track candidates
//    crossing the detector a very different locations (e.g. muon bundles).
//    The default maximum track opening angle is 15 degrees, but can be modified
//    via the SetTangmax memberfunction.
//    The default maximum track distance is 20 meters, but can be modified
//    via the SetTdistmax memberfunction. This memberfunction also allows to
//    specify whether the distance is determined within the detector volume or not.
//
//    The average of all the r0 and t0 values of the constituent TC's
//    of the jet will provide the r0 and t0 (i.e. reference point) of the jet.
//
//    The jet total momentum consists of the vector sum of the momenta of the
//    constituent TC's. This implies that the qtc quality numbers of the various
//    TC's define a weight for each track in the construction of the jet direction.
//    In addition it means that the total jet momentum represents the sum of the
//    qtc quality numbers of the constituent TC's weighted by the opening angles
//    between the various TC's.  
//    As such each jet is given an absolute quality number defined as :
//
//      qtcjet=|jet momentum|/ntracks 
//
//    This jet quality number is refined on basis of the number of hits
//    associated to the jet as :
//
//      qtcjet=qtcjet+0.2*(nah-nahmax)
//
//    where we have defined :
//
//    nah    : Number of associated hits for the specific jet.
//    nahmax : Maximum number of associated hits encountered for the jets.
//
//    This qtcjet value is then used to order the various jets w.r.t.
//    decreasing qtcjet quality number.
//
//    Note : The qtcjet value is stored as "energy" of the jet, such that
//           it is always available for each jet and can also be used for
//           ordering the jets according to this value using the generic
//           AliEvent::SortJets() facility. 
//
// 5) The jets (after having been ordered w.r.t. decreasing qtcjet value)
//    are merged when their directions are within a certain maximum opening angle.
//    In addition a jet must within a certain maximum distance of the starting jet
//    in order to get merged. 
//    The latter criterion prevents merging of (nearly) parallel tracks/jets
//    crossing the detector a very different locations (e.g. muon bundles).
//    The jet ordering before the merging process is essential, since the starting jet
//    will "eat up" the jets that will be merged into it. 
//    The jet ordering ensures that the jet with the highest quality number will
//    always initiate the merging process.
//    The default maximum opening angle is half the TC maximum opening angle,
//    but can be modified via the SetJangmax memberfunction. This memberfunction
//    also allows to specify whether jet merging will be performed iteratively or not.
//    In case iteration has been activated, the jet ordering is performed after each
//    iteration step. This has to be done because since the quality numbers of the
//    resulting merged jets have been automatically updated in the merging process.
//    
//    The default maximum jet distance is 30 meters, but can be modified
//    via the SetJdistmax memberfunction. This memberfunction also allows to
//    specify whether the distance is determined within the detector volume or not.
//
//    Note : Setting the maximum jet opening angle to <=0 will prevent
//           the merging of jets.
//
//    The average of all the r0 and t0 values of the merged jets will provide
//    the r0 and t0 (i.e. reference point) of the final jet.
//
// 6) The remaining (merged) jets are ordered w.r.t. decreasing jet quality number.
//    As such the jet with the highest quality number will be the first one
//    in the list, which will result in the fact that the final tracks are also
//    ordered w.r.t. decreasing quality number, as outlined hereafter.
//    Each remaining jet will provide the parameters (e.g. direction)
//    for a reconstructed track.
//    The track 3-momentum is set to the total jet 3-momentum, normalised
//    to 1 GeV. The mass and charge of the track are left 0.
//    The reference point data of the jet will provide the r0 and t0
//    (i.e. reference point) of the track.
//
//    All these tracks will be stored in the IceEvent structure with as
//    default "IceDwalk" as the name of the track.
//    This track name identifier can be modified by the user via the
//    SetTrackName() memberfunction. This will allow unique identification
//    of tracks which are produced when re-processing existing data with
//    different criteria.
//    By default the charge of the produced tracks is set to 0, since
//    no distinction can be made between positive or negative tracks.
//    However, the user can define the track charge by invokation
//    of the memberfunction SetCharge().
//    This facility may be used to distinguish tracks produced by the
//    various reconstruction algorithms in a (3D) colour display
//    (see the class AliHelix for further details).  
//
//    Note : In case the maximum jet opening angle was specified <0,
//           only the jet with the highest quality number will appear
//           as a reconstructed track in the IceEvent structure.
//           This will allow comparison with the standard Sieglinde
//           single track direct walk reconstruction results. 
//    
// For further details the user is referred to NIM A524 (2004) 169.
//
// Note : This algorithm works best on data which has been calibrated
//        (IceCalibrate), cross talk corrected (IceXtalk) and cleaned
//        from noise hits etc. (IceCleanHits).
//
//--- Author: Nick van Eijndhoven 07-oct-2005 Utrecht University
//- Modified: NvE $Date$ Utrecht University
///////////////////////////////////////////////////////////////////////////
 
#include "IceDwalk.h"
#include "Riostream.h"

ClassImp(IceDwalk) // Class implementation to enable ROOT I/O

IceDwalk::IceDwalk(const char* name,const char* title) : TTask(name,title)
{
// Default constructor.
// The various reconstruction parameters are initialised to the values
// as mentioned in NIM A524 (2004) 179-180.
// The newly introduced angular separation parameter for jet merging
// is initialised as half the value of the angular separation parameter
// for track candidate clustering.    
 fEvt=0;
 fDmin=75;
 fDtmarg=0;
 fMaxdhit=3.07126;
 fTangmax=15;
 fTdistmax=20;
 fTinvol=1;
 fJangmax=fTangmax/2.;
 fJiterate=1;
 fJdistmax=30;
 fJinvol=1;
 fMaxmodA=999;
 fMinmodA=0;
 fMaxhitsA=1;
 fVgroup=1;
 fTrackname="";
 fCharge=0;
}
///////////////////////////////////////////////////////////////////////////
IceDwalk::~IceDwalk()
{
// Default destructor.
}
///////////////////////////////////////////////////////////////////////////
void IceDwalk::SetDmin(Float_t d)
{
// Set minimum hit distance (in m) to form a track element.
// In the constructor the default has been set to 75 meter.
 fDmin=d;
}
///////////////////////////////////////////////////////////////////////////
void IceDwalk::SetDtmarg(Int_t dt)
{
// Set maximum hit time difference margin (in ns) for track elements.
// In the constructor the default has been set to 0 ns.
 fDtmarg=dt;
}
///////////////////////////////////////////////////////////////////////////
void IceDwalk::SetMaxDhit(Float_t d)
{
// Set maximum distance (in scattering length) for a hit to get associated.
// In the constructor the default has been set to 2 lambda_scat.
 fMaxdhit=d;
}
///////////////////////////////////////////////////////////////////////////
void IceDwalk::SetTangmax(Float_t ang)
{
// Set maximum angular separation (in deg) for track candidate clustering
// into jets.
// In the constructor the default has been set to 15 deg, in accordance
// to NIM A524 (2004) 180.
//
// Note : This function also sets automatically the value of the maximum
//        angular separation for jet merging into 1 single track to ang/2.
//        In order to specify a different max. jet merging separation angle,
//        one has to invoke the memberfunction SetJangmax afterwards.
 
 fTangmax=ang;
 fJangmax=ang/2.;
}
///////////////////////////////////////////////////////////////////////////
void IceDwalk::SetTdistmax(Float_t d,Int_t invol)
{
// Set maximum distance (in m) of the two track candidates in the track
// clustering process.
// The distance between the two tracks can be determined restricted to the
// detector volume (invol=1) or in the overall space (invol=0).  
// The former will prevent clustering of (nearly) parallel tracks which cross
// the detector volume at very different locations, whereas the latter will
// enable clustering of tracks with a common location of origin (e.g. muon
// bundles from an air shower) even if they cross the detector volume at
// very different locations. 
// At invokation of this memberfunction the default is invol=1.
// In the constructor the default has been set to 20 meter with invol=1.
 
 fTdistmax=d;
 fTinvol=invol;
}
///////////////////////////////////////////////////////////////////////////
void IceDwalk::SetJangmax(Float_t ang,Int_t iter)
{
// Set angular separation (in deg) within which jets are merged into 1
// single track.
// The merging process is a dynamic procedure and can be carried out by
// iteration (iter=1) until no further merging of the various jets occurs anymore.
// However, by specification of iter=0 the user can also select to go only
// once through all the jet combinations to check for mergers.
// For large events the latter will in general result in more track candidates.  
// At invokation of this memberfunction the default is iter=1.
// In the constructor the default angle has been set 7.5 deg, being half
// of the value of the default track candidate clustering separation angle.
// The iteration flag was set to 1 in the constructor.
//
// Notes :
// -------
// 1)  Setting ang=0 will prevent jet merging.
//     Consequently, every jet will appear as a separate track in the
//     reconstruction result.  
// 2)  Setting ang<0 will prevent jet merging.
//     In addition, only the jet with the maximum number of tracks will
//     appear as a track in the reconstruction result.
//     This situation resembles the standard Sieglinde direct walk processing
//     and as such can be used to perform comparison studies.

 fJangmax=ang;
 fJiterate=iter;
}
///////////////////////////////////////////////////////////////////////////
void IceDwalk::SetJdistmax(Float_t d,Int_t invol)
{
// Set maximum distance (in m) of the two jets in the jet merging process.
// The distance between the two jets can be determined restricted to the
// detector volume (invol=1) or in the overall space (invol=0).  
// The former will prevent clustering of (nearly) parallel tracks which cross
// the detector volume at very different locations, whereas the latter will
// enable clustering of tracks with a common location of origin (e.g. muon
// bundles from an air shower) even if they cross the detector volume at
// very different locations. 
// At invokation of this memberfunction the default is invol=1.
// In the constructor the default has been set to 30 meter with invol=1.
 
 fJdistmax=d;
 fJinvol=invol;
}
///////////////////////////////////////////////////////////////////////////
void IceDwalk::SetMaxModA(Int_t nmax)
{
// Set the maximum number of good Amanda modules that may have fired
// in order to process this event.
// This allows suppression of processing (high-energy) cascade events
// with this direct walk tracking to prevent waisting cpu time for cases
// in which tracking doesn't make sense anyhow.
// Furthermore it allows selection of low multiplicity events for processing.
// By default the maximum number of Amanda modules is set to 999 in the ctor,
// which implies no selection on maximum module multiplicity.
// See also the memberfunction SetMinModA().
 fMaxmodA=nmax;
}
///////////////////////////////////////////////////////////////////////////
void IceDwalk::SetMinModA(Int_t nmin)
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
void IceDwalk::SetMaxHitsA(Int_t nmax)
{
// Set the maximum number of good hits per Amanda module to be processed.
//
// Special values :
// nmax = 0 : No maximum limit set; all good hits will be processed
//      < 0 : No hits will be processed
//
// In case the user selects a maximum number of good hits per module, all the
// hits of each module will be ordered w.r.t. increasing hit time (LE).
// This allows selection of processing e.g. only the first hits etc...
// By default the maximum number of hits per Amanda modules is set to 1 in the ctor,
// which implies processing only the first good hit of each Amanda OM.
 fMaxhitsA=nmax;
}
///////////////////////////////////////////////////////////////////////////
void IceDwalk::SetVgroupUsage(Int_t flag)
{
// (De)activate the distinction between v_phase and v_group of the Cherenkov light.
//
// flag = 0 : No distinction between v_phase and v_group
//      = 1 : Separate treatment of v_phase and v_group
//
// By default the distinction between v_phase and v_group is activated
// in the constructor of this class.
 fVgroup=flag;
}
///////////////////////////////////////////////////////////////////////////
void IceDwalk::SetTrackName(TString s)
{
// Set (alternative) name identifier for the produced first guess tracks.
// This allows unique identification of (newly) produced direct walk tracks
// in case of re-processing of existing data with different criteria.
// By default the produced first guess tracks have the name of the class
// by which they were produced.
 fTrackname=s;
}
///////////////////////////////////////////////////////////////////////////
void IceDwalk::SetCharge(Float_t charge)
{
// Set user defined charge for the produced first guess tracks.
// This allows identification of these tracks on color displays.
// By default the produced first guess tracks have charge=0
// which is set in the constructor of this class.
 fCharge=charge;
}
///////////////////////////////////////////////////////////////////////////
void IceDwalk::Exec(Option_t* opt)
{
// Implementation of the direct walk track reconstruction.

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
 params.SetNameTitle(ClassName(),"Reco parameters");
 params.SetSlotName("Dmin",1);
 params.SetSlotName("Dtmarg",2);
 params.SetSlotName("Maxdhit",3);
 params.SetSlotName("Tangmax",4);
 params.SetSlotName("Tdistmax",5);
 params.SetSlotName("Tinvol",6);
 params.SetSlotName("Jangmax",7);
 params.SetSlotName("Jiterate",8);
 params.SetSlotName("Jdistmax",9);
 params.SetSlotName("Jinvol",10);
 params.SetSlotName("MaxmodA",11);
 params.SetSlotName("MinmodA",12);
 params.SetSlotName("MaxhitsA",13);
 params.SetSlotName("Vgroup",14);

 params.SetSignal(fDmin,1);
 params.SetSignal(fDtmarg,2);
 params.SetSignal(fMaxdhit,3);
 params.SetSignal(fTangmax,4);
 params.SetSignal(fTdistmax,5);
 params.SetSignal(fTinvol,6);
 params.SetSignal(fJangmax,7);
 params.SetSignal(fJiterate,8);
 params.SetSignal(fJdistmax,9);
 params.SetSignal(fJinvol,10);
 params.SetSignal(fMaxmodA,11);
 params.SetSignal(fMinmodA,12);
 params.SetSignal(fMaxhitsA,13);
 params.SetSignal(fVgroup,14);

 fEvt->AddDevice(params);

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

 const Float_t c=0.299792458; // Light speed in vacuum in meters per ns

 // Storage of track elements.
 TObjArray tes;
 tes.SetOwner();

 AliPosition r1;
 AliPosition r2;
 Ali3Vector r12;
 Ali3Vector rsum;
 AliPosition r0;
 TObjArray hits1;
 TObjArray hits2;
 Int_t nh1,nh2;
 AliSignal* sx1=0;
 AliSignal* sx2=0;
 Float_t dist=0;
 Float_t t1,t2,dt,t0;
 Float_t dtmax;
 TObjArray hits;
 TObjArray* ordered;

 // Check the hits of Amanda OM pairs for possible track elements.
 // Also all the good hits are stored in the meantime (to save CPU time)
 // for hit association with the various track elements lateron.
 AliTrack* te=0;
 for (Int_t i1=0; i1<naoms; i1++) // First OM of the pair
 {
  IceGOM* omx1=(IceGOM*)aoms->At(i1);
  if (!omx1) continue;
  if (omx1->GetDeadValue("LE")) continue;
  r1=omx1->GetPosition();
  // Select all the good hits of this first OM
  hits1.Clear();
  // Determine the max. number of hits to be processed for this OM
  ordered=0;
  if (fMaxhitsA>0 && omx1->GetNhits()>fMaxhitsA) ordered=omx1->SortHits("LE",1,0,7);
  nh1=0;
  for (Int_t j1=1; j1<=omx1->GetNhits(); j1++)
  {
   if (ordered)
   {
    if (nh1>=fMaxhitsA) break;
    sx1=(AliSignal*)ordered->At(j1-1);
   }
   else
   {
    sx1=omx1->GetHit(j1);
   }
   if (!sx1) continue;
   if (sx1->GetDeadValue("ADC") || sx1->GetDeadValue("LE") || sx1->GetDeadValue("TOT")) continue;
   hits1.Add(sx1);
   // Also store all good hits in the total hit array
   hits.Add(sx1);
   nh1++;
  }

  // No further pair to be formed with the last OM in the list 
  if (i1==(naoms-1)) break;

  nh1=hits1.GetEntries();
  if (!nh1) continue;

  for (Int_t i2=i1+1; i2<naoms; i2++) // Second OM of the pair
  {
   IceGOM* omx2=(IceGOM*)aoms->At(i2);
   if (!omx2) continue;
   if (omx2->GetDeadValue("LE")) continue;
   r2=omx2->GetPosition();
   r12=r2-r1;
   dist=r12.GetNorm();

   if (dist<=fDmin) continue;

   // Select all the good hits of this second OM
   hits2.Clear();
   // Determine the max. number of hits to be processed for this OM
   ordered=0;
   if (fMaxhitsA>0 && omx2->GetNhits()>fMaxhitsA) ordered=omx2->SortHits("LE",1,0,7);
   nh2=0;
   for (Int_t j2=1; j2<=omx2->GetNhits(); j2++)
   {
    if (ordered)
    {
     if (nh2>=fMaxhitsA) break;
     sx2=(AliSignal*)ordered->At(j2-1);
    }
    else
    {
     sx2=omx2->GetHit(j2);
    }
    if (!sx2) continue;
    if (sx2->GetDeadValue("ADC") || sx2->GetDeadValue("LE") || sx2->GetDeadValue("TOT")) continue;
    hits2.Add(sx2);
    nh2++;
   }
 
   nh2=hits2.GetEntries();
   if (!nh2) continue;

   // Position r0 in between the two OMs and normalised relative direction r12
   rsum=(r1+r2)/2.;
   r0.SetPosition((Ali3Vector&)rsum);
   r12/=dist;

   // Check all hit pair combinations of these two OMs for possible track elements  
   dtmax=dist/c+float(fDtmarg);
   for (Int_t ih1=0; ih1<nh1; ih1++) // Hits of first OM
   {
    sx1=(AliSignal*)hits1.At(ih1);
    if (!sx1) continue;
    for (Int_t ih2=0; ih2<nh2; ih2++) // Hits of second OM
    {
     sx2=(AliSignal*)hits2.At(ih2);
     if (!sx2) continue;
     t1=sx1->GetSignal("LE",7);
     t2=sx2->GetSignal("LE",7);
     dt=t2-t1;
     t0=(t1+t2)/2.;

     if (fabs(dt)>=dtmax) continue;

     te=new AliTrack();
     tes.Add(te);
     if (dt<0) r12*=-1.;
     r0.SetTimestamp((AliTimestamp&)*fEvt);
     AliTimestamp* tsx=r0.GetTimestamp();
     tsx->Add(0,0,(int)t0);
     te->SetReferencePoint(r0);
     te->Set3Momentum(r12);
    }
   }
  } // end of loop over the second OM of the pair
 } // end of loop over first OM of the pair

 // Association of hits to the various track elements.
 Float_t qmax=0;
 Int_t nahmax=0;
 AssociateHits(tes,hits,qmax,nahmax);

 // Selection on quality (Q value) in case of multiple track candidates
 SelectQvalue(tes,qmax);

 Int_t nte=tes.GetEntries();
 if (!nte) return;

 // Clustering of track candidates into jets
 TObjArray jets;
 jets.SetOwner();
 ClusterTracks(tes,jets,qmax);

 Int_t njets=jets.GetEntries();
 if (!njets) return;

 // Order the jets w.r.t. decreasing quality value
 ordered=fEvt->SortJets(-2,&jets);
 TObjArray jets2(*ordered);

 // Merging f jets
 MergeJets(jets2);

 // Production and storage of the final tracks
 StoreTracks(jets2);
}
///////////////////////////////////////////////////////////////////////////
void IceDwalk::AssociateHits(TObjArray& tes,TObjArray& hits,Float_t& qmax,Int_t& nahmax)
{
 // Association of hits to the various track elements.

 const Float_t pi=acos(-1.);
 const Float_t c=0.299792458;         // Light speed in vacuum in meters per ns
 const Float_t npice=1.31768387;      // Phase refractive index (c/v_phase) of ice
 const Float_t ngice=1.35075806;      // Group refractive index (c/v_group) of ice
 const Float_t lambda=33.3;           // Light scattering length in ice
 const Float_t thetac=acos(1./npice); // Cherenkov angle (in radians)

 // Angular reduction of complement of thetac due to v_phase and v_group difference
 Float_t alphac=0;
 if (fVgroup) alphac=atan((1.-npice/ngice)/sqrt(npice*npice-1.));

 Int_t nte=tes.GetEntries();
 Int_t nh=hits.GetEntries();
 Float_t d=0;
 Ali3Vector p;
 AliPosition r1;
 AliPosition r2;
 Ali3Vector r12;
 Float_t t1;
 Float_t dist,t0,tgeo,tres;
 AliSample levers;      // Statistics of the assoc. hit lever arms
 levers.SetStoreMode(1);// Enable median calculation
 AliSample hprojs;      // Statistics of the assoc. hit position projections on the track w.r.t. r0
 hprojs.SetStoreMode(1);// Enable median calculation
 AliSample times;       // Statistics of the time residuals of the associated hits
 times.SetStoreMode(1); // Enable median calculation
 AliSignal fit;         // Storage of Q value etc... for each track candidate
 fit.AddNamedSlot("QTC");
 fit.AddNamedSlot("SpanL");
 fit.AddNamedSlot("MedianL");
 fit.AddNamedSlot("MeanL");
 fit.AddNamedSlot("SigmaL");
 fit.AddNamedSlot("SpreadL");
 fit.AddNamedSlot("ExpSpreadL");
 fit.AddNamedSlot("Span");
 fit.AddNamedSlot("Median");
 fit.AddNamedSlot("Mean");
 fit.AddNamedSlot("Sigma");
 fit.AddNamedSlot("Spread");
 fit.AddNamedSlot("ExpSpread");
 fit.AddNamedSlot("MedianT");
 fit.AddNamedSlot("MeanT");
 fit.AddNamedSlot("SigmaT");
 fit.AddNamedSlot("SpreadT");
 fit.AddNamedSlot("term1");
 fit.AddNamedSlot("term2");
 fit.AddNamedSlot("term3");
 fit.AddNamedSlot("term4");
 fit.AddNamedSlot("term5");
 Float_t qtc=0;
 Int_t nah; // Number of associated hits for a certain TE
 Float_t lmin,lmax,spanl,medianl,meanl,sigmal,spreadl,expspreadl;
 Float_t hproj,hprojmin,hprojmax,span,median,mean,sigma,spread,expspread;
 Float_t mediant,meant,sigmat,spreadt;
 Float_t term1,term2,term3,term4,term5;
 qmax=0;
 nahmax=0;
 for (Int_t jte=0; jte<nte; jte++)
 {
  AliTrack* te=(AliTrack*)tes.At(jte);
  if (!te) continue;
  AliPosition* tr0=te->GetReferencePoint();
  AliTimestamp* tt0=tr0->GetTimestamp();
  t0=fEvt->GetDifference(tt0,"ns");
  p=te->Get3Momentum();
  levers.Reset();
  hprojs.Reset();
  times.Reset();
  for (Int_t jh=0; jh<nh; jh++)
  {
   AliSignal* sx1=(AliSignal*)hits.At(jh);
   if (!sx1) continue;
   IceGOM* omx=(IceGOM*)sx1->GetDevice();
   if (!omx) continue;
   r1=omx->GetPosition();
   d=te->GetDistance(r1);
   r12=r1-(*tr0);
   hproj=p.Dot(r12);
   dist=hproj+d/tan(pi/2.-thetac-alphac);
   tgeo=t0+dist/c;
   t1=sx1->GetSignal("LE",7);
   tres=t1-tgeo;

   d=d/sin(thetac); // The distance traveled by a cherenkov photon

   if (tres<-30 || tres>300 || d>fMaxdhit*lambda) continue;

   // Associate this hit to the TE
   te->AddSignal(*sx1);
   levers.Enter(fabs(hproj));
   hprojs.Enter(hproj);
   times.Enter(tres);
  }

  // Determine the Q quality of the various TE's.
  // Good quality TE's will be called track candidates (TC's)
  nah=te->GetNsignals();
  if (nah>nahmax) nahmax=nah;
  lmin=levers.GetMinimum(1);
  lmax=levers.GetMaximum(1);
  spanl=lmax-lmin;
  medianl=levers.GetMedian(1);
  meanl=levers.GetMean(1);
  sigmal=levers.GetSigma(1);
  spreadl=levers.GetSpread(1);
  // Expected spread for a flat distribution
  expspreadl=0;
  if (spanl>0) expspreadl=(0.5*pow(lmin,2)+0.5*pow(lmax,2)+pow(medianl,2)-medianl*(lmin+lmax))/spanl;
  hprojmin=hprojs.GetMinimum(1);
  hprojmax=hprojs.GetMaximum(1);
  span=hprojmax-hprojmin;
  median=hprojs.GetMedian(1);
  mean=hprojs.GetMean(1);
  sigma=hprojs.GetSigma(1);
  spread=hprojs.GetSpread(1);
  // Expected spread for a flat distribution
  expspread=0;
  if (span>0) expspread=(0.5*pow(hprojmin,2)+0.5*pow(hprojmax,2)+pow(median,2)-median*(hprojmin+hprojmax))/span;
  mediant=times.GetMedian(1);
  meant=times.GetMean(1);
  sigmat=times.GetSigma(1);
  spreadt=times.GetSpread(1);

  term1=0;
  if (span>0) term1=2.*spread/span;

  term2=0;
  if (spanl>0) term2=2.*spreadl/spanl;

  term3=0;
  if (spread>0) term3=fabs(spread-expspread)/spread;

  term4=0;
  if (spreadl>0) term4=fabs(spreadl-expspreadl)/spreadl;

  term5=0;
  if (spreadt>0) term5=fabs(mediant)/spreadt;

  qtc=float(nah)*(term1+term2)-term3-term4-term5;
  if (fabs(median)>span/2.) qtc=0; // Require projected hits on both sides of r0

  fit.SetSignal(qtc,"QTC");
  fit.SetSignal(spanl,"SpanL");
  fit.SetSignal(medianl,"MedianL");
  fit.SetSignal(meanl,"MeanL");
  fit.SetSignal(sigmal,"SigmaL");
  fit.SetSignal(spreadl,"SpreadL");
  fit.SetSignal(expspreadl,"ExpSpreadL");
  fit.SetSignal(span,"Span");
  fit.SetSignal(median,"Median");
  fit.SetSignal(mean,"Mean");
  fit.SetSignal(sigma,"Sigma");
  fit.SetSignal(spread,"Spread");
  fit.SetSignal(expspread,"ExpSpread");
  fit.SetSignal(mediant,"MedianT");
  fit.SetSignal(meant,"MeanT");
  fit.SetSignal(sigmat,"SigmaT");
  fit.SetSignal(spreadt,"SpreadT");
  fit.SetSignal(term1,"term1");
  fit.SetSignal(term2,"term2");
  fit.SetSignal(term3,"term3");
  fit.SetSignal(term4,"term4");
  fit.SetSignal(term5,"term5");
  te->SetFitDetails(fit);
  if (qtc>qmax) qmax=qtc;
 }
}
///////////////////////////////////////////////////////////////////////////
void IceDwalk::SelectQvalue(TObjArray& tes,Float_t qmax)
{
 // Perform selection on Q value in case of multiple track candidates

 Int_t nte=tes.GetEntries();
 Int_t nah;
 Float_t qtc;
 Ali3Vector p;
 for (Int_t jtc=0; jtc<nte; jtc++)
 {
  AliTrack* te=(AliTrack*)tes.At(jtc);
  if (!te) continue;
  nah=te->GetNsignals();
  AliSignal* sx1=(AliSignal*)te->GetFitDetails();
  qtc=-1;
  if (sx1) qtc=sx1->GetSignal("QTC");
  if (!nah || qtc<0.8*qmax)
  {
   tes.RemoveAt(jtc);
   delete te;
  }
  else // Set Q value as momentum to provide a weight for jet clustering
  {
   if (qtc>0)
   {
    p=te->Get3Momentum();
    p*=qtc;
    te->Set3Momentum(p);
   }
  }
 } 
 tes.Compress();
}
///////////////////////////////////////////////////////////////////////////
void IceDwalk::ClusterTracks(TObjArray& tes,TObjArray& jets,Float_t qmax)
{
 // Cluster track candidates within a certain opening angle into jets.
 // Also the track should be within a certain maximum distance of the
 // starting track in order to get clustered.
 // The latter prevents clustering of (nearly) parallel track candidates
 // crossing the detector a very different locations (e.g. muon bundles).
 // The average r0 and t0 of the constituent tracks will be taken as the
 // jet reference point. 

 Int_t nte=tes.GetEntries();
 Float_t ang=0;
 AliSample pos;
 AliSample time;
 Float_t vec[3],err[3];
 AliPosition r0;
 Float_t t0,dist,dist2;
 Int_t nah=0,nahmax=0; // Determine the max. number of associated hits for the jets
 Float_t qtc;
 for (Int_t jtc1=0; jtc1<nte; jtc1++)
 {
  AliTrack* te=(AliTrack*)tes.At(jtc1);
  if (!te) continue;
  AliPosition* x1=te->GetReferencePoint();
  if (!x1) continue;
  AliTimestamp* ts1=x1->GetTimestamp();
  if (!ts1) continue;
  AliJet* jx=new AliJet();
  jx->AddTrack(te);
  pos.Reset();
  time.Reset();
  x1->GetPosition(vec,"car");
  pos.Enter(vec[0],vec[1],vec[2]);
  t0=fEvt->GetDifference(ts1,"ns");
  time.Enter(t0);
  for (Int_t jtc2=0; jtc2<nte; jtc2++)
  {
   if (jtc2==jtc1) continue;
   AliTrack* te2=(AliTrack*)tes.At(jtc2);
   if (!te2) continue;
   ang=te->GetOpeningAngle(*te2,"deg");
   if (ang<=fTangmax)
   {
    AliPosition* x2=te2->GetReferencePoint();
    if (!x2) continue;
    AliTimestamp* ts2=x2->GetTimestamp();
    if (!ts2) continue;
    if (!fTinvol)
    {
     dist=te->GetDistance(te2);
    }
    else
    {
     dist=te->GetDistance(x2);
     dist2=te2->GetDistance(x1);
     if (dist2<dist) dist=dist2;
    }
    if (dist<=fTdistmax)
    {
     x2->GetPosition(vec,"car");
     pos.Enter(vec[0],vec[1],vec[2]);
     t0=fEvt->GetDifference(ts2,"ns");
     time.Enter(t0);
     jx->AddTrack(te2);
    }
   }
  }

  // Set the reference point data for this jet
  for (Int_t j=1; j<=3; j++)
  {
   vec[j-1]=pos.GetMean(j);
   err[j-1]=pos.GetSigma(j);
  }
  r0.SetPosition(vec,"car");
  r0.SetPositionErrors(err,"car");
  r0.SetTimestamp((AliTimestamp&)*fEvt);
  AliTimestamp* jt0=r0.GetTimestamp();
  t0=time.GetMean(1);
  jt0->Add(0,0,(int)t0);
  jx->SetReferencePoint(r0);

  // Store this jet for further processing if ntracks>1
  if (jx->GetNtracks() > 1 || fTangmax<=0)
  {
   jets.Add(jx);
   nah=jx->GetNsignals();
   if (nah>nahmax) nahmax=nah;
  }
  else // Only keep single-track jets which have qtc=qmax 
  {
   AliSignal* sx1=(AliSignal*)te->GetFitDetails();
   qtc=-1;
   if (sx1) qtc=sx1->GetSignal("QTC");
   if (qtc>=(qmax-1.e-10))
   {
    jets.Add(jx);
    nah=jx->GetNsignals();
    if (nah>nahmax) nahmax=nah;
   }
   else
   {
    delete jx;
   }
  }
 }

 Int_t njets=jets.GetEntries();
 if (!njets) return;

 // The sum of 0.15*(nah-nahmax) and average qtc value per track for each jet
 // will be stored as the jet energy to enable sorting on this value lateron
 Float_t sortval=0;
 Int_t ntk=0;
 for (Int_t ijet=0; ijet<njets; ijet++)
 {
  AliJet* jx=(AliJet*)jets.At(ijet);
  if (!jx) continue;
  nah=jx->GetNsignals();
  ntk=jx->GetNtracks();
  sortval=0.15*float(nah-nahmax);
  if (ntk) sortval+=jx->GetMomentum()/float(ntk);
  jx->SetScalar(sortval);
 }
}
///////////////////////////////////////////////////////////////////////////
void IceDwalk::MergeJets(TObjArray& jets2)
{
 // Merge jets within a certain opening to provide the final track(s).
 // Also the jet should be within a certain maximum distance of the
 // starting jet in order to get merged.
 // The latter prevents merging of (nearly) parallel jets/tracks
 // crossing the detector a very different locations (e.g. muon bundles).
 // The average r0 and t0 of the constituent jets will be taken as the
 // final reference point. 

 Int_t njets=jets2.GetEntries();
 AliJet* jx1=0;
 AliJet* jx2=0;
 Int_t merged=1;
 Int_t ntk,nah,nahmax;
 Float_t ang,dist,dist2,t0;
 AliSample pos;
 AliSample time;
 AliPosition r0;
 Float_t vec[3],err[3];
 Float_t sortval;
 if (fJangmax>=0)
 {
  while (merged)
  {
   merged=0;
   nahmax=0;
   for (Int_t jet1=0; jet1<njets; jet1++)
   {
    jx1=(AliJet*)jets2.At(jet1);
    if (!jx1) continue;
    AliPosition* x1=jx1->GetReferencePoint();
    if (!x1) continue;
    AliTimestamp* ts1=x1->GetTimestamp();
    if (!ts1) continue;
    pos.Reset();
    time.Reset();
    x1->GetPosition(vec,"car");
    pos.Enter(vec[0],vec[1],vec[2]);
    t0=fEvt->GetDifference(ts1,"ns");
    time.Enter(t0);
    for (Int_t jet2=0; jet2<njets; jet2++)
    {
     jx2=(AliJet*)jets2.At(jet2);
     if (!jx2 || jet2==jet1) continue;
     AliPosition* x2=jx2->GetReferencePoint();
     if (!x2) continue;
     AliTimestamp* ts2=x2->GetTimestamp();
     if (!ts2) continue;
     ang=jx1->GetOpeningAngle(*jx2,"deg");
     if (ang<=fJangmax)
     {
      if (!fJinvol)
      {
       dist=jx1->GetDistance(jx2);
      }
      else
      {
       dist=jx1->GetDistance(x2);
       dist2=jx2->GetDistance(x1);
       if (dist2<dist) dist=dist2;
      }
      if (dist<=fJdistmax)
      {
       x2->GetPosition(vec,"car");
       pos.Enter(vec[0],vec[1],vec[2]);
       t0=fEvt->GetDifference(ts2,"ns");
       time.Enter(t0);
       for (Int_t jtk=1; jtk<=jx2->GetNtracks(); jtk++)
       {
        AliTrack* te=jx2->GetTrack(jtk);
        if (!te) continue;
        jx1->AddTrack(te);
       }
       jets2.RemoveAt(jet2);
       if (fJiterate) merged=1;
      }
     }
    } // End of jet2 loop

    // Set the reference point data for this jet
    for (Int_t k=1; k<=3; k++)
    {
     vec[k-1]=pos.GetMean(k);
     err[k-1]=pos.GetSigma(k);
    }
    r0.SetPosition(vec,"car");
    r0.SetPositionErrors(err,"car");
    r0.SetTimestamp((AliTimestamp&)*fEvt);
    AliTimestamp* jt0=r0.GetTimestamp();
    t0=time.GetMean(1);
    jt0->Add(0,0,(int)t0);
    jx1->SetReferencePoint(r0);

    nah=jx1->GetNsignals();
    if (nah>nahmax) nahmax=nah;
   } // End of jet1 loop

   jets2.Compress();

   // The sum of 0.15*(nah-nahmax) and average qtc value per track for each jet
   // will be stored as the jet energy to enable sorting on this value
   for (Int_t jjet=0; jjet<njets; jjet++)
   {
    AliJet* jx=(AliJet*)jets2.At(jjet);
    if (!jx) continue;
    nah=jx->GetNsignals();
    ntk=jx->GetNtracks();
    sortval=0.15*float(nah-nahmax);
    if (ntk) sortval+=jx->GetMomentum()/float(ntk);
    jx->SetScalar(sortval);
   }

   // Order the jets w.r.t. decreasing quality value
   TObjArray* ordered=fEvt->SortJets(-2,&jets2);
   njets=ordered->GetEntries();
   jets2.Clear();
   for (Int_t icopy=0; icopy<njets; icopy++)
   {
    jets2.Add(ordered->At(icopy));
   }
  } // End of iterative while loop
 }
}
///////////////////////////////////////////////////////////////////////////
void IceDwalk::StoreTracks(TObjArray& jets2)
{
 // Store every jet as a reconstructed track in the event structure.
 // The jet 3-momentum (normalised to 1) and reference point
 // (i.e.the average r0 and t0 of the constituent tracks) will make up
 // the final track parameters.
 // All the associated hits of all the constituent tracks of the jet
 // will be associated to the final track.
 // In case the jet angular separation was set <0, only the jet with
 // the maximum number of tracks (i.e. the first one in the array)
 // will be used to form a track. This will allow comparison with
 // the standard Sieglinde processing.

 Int_t njets=jets2.GetEntries();

 if (fTrackname=="") fTrackname=ClassName();
 TString title=ClassName();
 title+=" reco track";
 AliTrack t; 
 t.SetNameTitle(fTrackname.Data(),title.Data());
 t.SetCharge(fCharge);
 Ali3Vector p;
 for (Int_t jet=0; jet<njets; jet++)
 {
  AliJet* jx=(AliJet*)jets2.At(jet);
  if (!jx) continue;
  AliPosition* ref=jx->GetReferencePoint();
  if (!ref) continue;
  fEvt->AddTrack(t);
  AliTrack* trk=fEvt->GetTrack(fEvt->GetNtracks());
  if (!trk) continue;
  trk->SetId(fEvt->GetNtracks(1)+1);
  p=jx->Get3Momentum();
  p/=p.GetNorm();
  trk->Set3Momentum(p);
  trk->SetReferencePoint(*ref);
  AliTimestamp* tt0=ref->GetTimestamp();
  if (tt0) trk->SetTimestamp(*tt0);
  for (Int_t jt=1; jt<=jx->GetNtracks(); jt++)
  {
   AliTrack* tx=jx->GetTrack(jt);
   if (!tx) continue;
   for (Int_t is=1; is<=tx->GetNsignals(); is++)
   {
    AliSignal* sx1=tx->GetSignal(is);
    if (sx1) sx1->AddTrack(*trk);
   }
  }

  // Only take the jet with the highest quality number
  // (i.e. the first jet in the list) when the user had selected
  // this reconstruction mode.
  if (fJangmax<0) break;
 }
}
///////////////////////////////////////////////////////////////////////////
