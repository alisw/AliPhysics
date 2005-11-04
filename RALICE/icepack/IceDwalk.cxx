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
// The procedure is based on the method described in the Amanda publication
// in Nuclear Instruments and Methods A524 (2004) 179-180.
// However, the Amanda method has been extended with the intention to
// take also multiple (muon) tracks within 1 event into account.
// This will not only provide a means to reconstruct muon bundles and
// multiple track events in IceCube, but will also allow to reduce the
// background of faked upgoing muons as a result of multiple downgoing
// muons hitting the top and bottom parts of the detector. 
// The various reconstruction steps are summarised as follows :
//
// 1) Construction of track elements (TE's).
//    A track element is a straight line connecting two hits that
//    appeared at some minimum distance d and within some maximum
//    time difference dt.
//    The default values for d and dt are given in eq. (20) of the
//    NIM article, but can be modified by the appropriate Set functions.
//    For dt a default margin of 30 ns is used (according to eq. (20)),
//    but also this margin may be modified via the appropriate Set function.    
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
//      dhit < 25*(tres+30)^(1/4) meter
//
//    tres : time residual
//           Difference between the observed hit time and the time expected
//           for a direct photon hit.     
//    dhit : Distance between the hit and the TE
//
// 3) Construction of track candidates (TC's).
//    These are TE's that fulfill both the conditions
//
//      nah >= 10
//      sigmal >= 20 meter
// 
//    nah    : Number of associated hits.
//    sigmal : rms variance of the distances between r0 and the point on
//             the track which is closest to the various associated hits. 
//
// 4) Only track candidates are kept which fulfill the quality criterion
//    (see eq. (21) in the NIM article)
//
//     qtc >= 0.7*qtcmax
//
//     qtc=min(nah,0.3*sigmal+7)
//     qtcmax=max(qtc)
//
// 5) The surviving track candidates are clustered into jets when
//    their directions are within a certain maximum opening angle.
//    The default maximum opening angle is 15 degrees, but can be modified
//    via the SetTangsep memberfunction.
//
// 6) The jets are merged when their directions are within
//    a certain maximum opening angle. 
//    The default maximum opening angle is half the TC maximum opening angle,
//    but can be modified via the SetJangsep memberfunction.
//    Note : Setting the maximum jet opening angle to <=0 will prevent
//           the merging of jets.
//
// 7) The remaining jets are ordered w.r.t. decreasing number of tracks.
//    Each remaining jet will provide the parameters (e.g. direction)
//    for a reconstructed track.
//    The track 3-momentum is set to the total jet 3-momentum, normalised
//    to 1 GeV. The mass and charge of the track are left 0.
//    The average of all the r0 and t0 values of the constituent TC's
//    of the jet will provide the r0 and t0 of the track.
//    All these tracks will be stored in the IceEvent structure with "IceDwalk"
//    as the name of the track.
//    Note : In case the maximum jet opening angle was specified <0,
//           only the jet with the maximum number of tracks will appear
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
 fDmin=50;
 fDtmarg=30;
 fTangsep=15;
 fJangsep=fTangsep/2.;
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
// In the constructor the default has been set to 50 meter, in accordance
// to eq.(20) of NIM A524 (2004) 179.
 fDmin=d;
}
///////////////////////////////////////////////////////////////////////////
void IceDwalk::SetDtmarg(Int_t dt)
{
// Set maximum hit time difference margin (in ns) for track elements.
// In the constructor the default has been set to 30 ns, in accordance
// to eq.(20) of NIM A524 (2004) 179.
 fDtmarg=dt;
}
///////////////////////////////////////////////////////////////////////////
void IceDwalk::SetTangsep(Float_t ang)
{
// Set angular separation (in deg) within which track candidates are
// clustered into jets.
// In the constructor the default has been set to 15 deg, in accordance
// to NIM A524 (2004) 180.
//
// Note : This function also sets automatically the value of the angular
//        separation within which jets are merged into 1 single track
//        to ang/2.
//        In order to specify a different jet merging separation angle,
//        one has to invoke the memberfunction SetJangsep afterwards.
 
 fTangsep=ang;
 fJangsep=ang/2.;
}
///////////////////////////////////////////////////////////////////////////
void IceDwalk::SetJangsep(Float_t ang)
{
// Set angular separation (in deg) within which jets are merged into 1
// single track.
// In the constructor the default has been set 7.5 deg, being half of the
// value of the default track candidate clustering separation angle. 
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

 fJangsep=ang;
}
///////////////////////////////////////////////////////////////////////////
void IceDwalk::Exec(Option_t* opt)
{
// Implementation of the direct walk track reconstruction.

 TString name=opt;
 AliJob* parent=(AliJob*)(gROOT->GetListOfTasks()->FindObject(name.Data()));

 if (!parent) return;

 IceEvent* evt=(IceEvent*)parent->GetObject("IceEvent");
 if (!evt) return;

 Float_t c=0.3;                // Light speed in vacuum in meters per ns
 Float_t nice=1.33;            // Refractive index of ice
 Float_t thetac=acos(1./nice); // Cherenkov angle (in radians)

 // Storage of track elements with various time difference margins.
 // temap(i,j) holds the i-th track element (TE) with a time difference margin
 // of j*3 nanoseconds. Currently we use a maximum margin of 30 ns.
 TObjArray tes;
 tes.SetOwner();
 AliObjMatrix temap;

 Int_t* ntes=new Int_t[fDtmarg/3]; // Counter of TEs for each 3 ns margin slot
 for (Int_t i=0; i<fDtmarg/3; i++)
 {
  ntes[i]=0;
 }  

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
 Float_t dtmax,dttest;
 TObjArray hits;

 // Check the hits of Amanda OM pairs for posible track elements.
 // Also all the good hits are stored in the meantime (to save CPU time)
 // for hit association with the various track elements lateron.
 TObjArray* aoms=evt->GetDevices("IceAOM");
 Int_t naoms=aoms->GetEntries();
 AliTrack* te=0;
 Int_t ite=0;
 for (Int_t i1=0; i1<naoms; i1++) // First OM of the pair
 {
  IceGOM* omx1=(IceGOM*)aoms->At(i1);
  if (!omx1) continue;
  if (omx1->GetDeadValue("LE")) continue;
  r1=omx1->GetPosition();
  // Select all the good hits of this first OM
  hits1.Clear();
  for (Int_t j1=1; j1<=omx1->GetNhits(); j1++)
  {
   sx1=omx1->GetHit(j1);
   if (!sx1) continue;
   if (sx1->GetDeadValue("ADC") || sx1->GetDeadValue("LE") || sx1->GetDeadValue("TOT")) continue;
   hits1.Add(sx1);
   // Also store all good hits in the total hit array
   hits.Add(sx1);
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
   dtmax=dist/c+float(fDtmarg);
   if (dist<=fDmin) continue;
   // Select all the good hits of this second OM
   hits2.Clear();
   for (Int_t j2=1; j2<=omx2->GetNhits(); j2++)
   {
    sx2=omx2->GetHit(j2);
    if (!sx2) continue;
    if (sx2->GetDeadValue("ADC") || sx2->GetDeadValue("LE") || sx2->GetDeadValue("TOT")) continue;
    hits2.Add(sx2);
   }
 
   nh2=hits2.GetEntries();
   if (!nh2) continue;

   // Check all hit pair combinations of these two OMs for possible track elements  
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
     if (dt && fabs(dt)<dtmax)
     {
      te=new AliTrack();
      tes.Add(te);
      ite++;
      if (dt<0) r12*=-1.;
      rsum=(r1+r2)/2.;
      r0.SetPosition((Ali3Vector&)rsum);
      te->SetReferencePoint(r0);
      te->SetTimestamp((AliTimestamp&)*evt);
      AliTimestamp* tsx=te->GetTimestamp();
      tsx->Add(0,0,(int)t0);
      r12/=r12.GetNorm();
      te->Set3Momentum(r12);
      dttest=dtmax;
      for (Int_t jt=fDtmarg/3; jt>0; jt--)
      {
       if (fabs(dt)>=dttest) break;
       temap.EnterObject(ite,jt,te);
       ntes[jt-1]++;
       dttest-=3.;
      }
     }
    }
   }
  } // end of loop over the second OM of the pair
 } // end of loop over first OM of the pair

 // Association of hits to the various track elements
 // For the time being all track elements will be treated,
 // but in a later stage one could select only the TE's of a certain
 // 3 ns margin slot in the TE map to save CPU time.
 Int_t nte=tes.GetEntries();
 Int_t nh=hits.GetEntries();
 Float_t d=0;
 Ali3Vector p;
 Float_t tgeo,tres;
 AliSample levers;  // Statistics of the assoc. hit lever arms
 AliSignal fit;     // Storage of Q value etc... for each track candidate
 fit.SetSlotName("QTC",1);
 fit.SetSlotName("SIGMAL",2);
 Float_t qtc=0,qmax=0;
 Int_t nah;      // Number of associated hits for a certain TE
 Float_t sigmal; // The mean lever arm of the various associated hits
 for (Int_t jte=0; jte<nte; jte++)
 {
  te=(AliTrack*)tes.At(jte);
  if (!te) continue;
  p=te->Get3Momentum();
  AliTimestamp* tt0=te->GetTimestamp();
  t0=evt->GetDifference(tt0,"ns");
  AliPosition* tr0=te->GetReferencePoint();
  levers.Reset();
  for (Int_t jh=0; jh<nh; jh++)
  {
   sx1=(AliSignal*)hits.At(jh);
   if (!sx1) continue;
   IceGOM* omx=(IceGOM*)sx1->GetDevice();
   if (!omx) continue;
   r1=omx->GetPosition();
   d=tr0->GetDistance(r1);
   d*=sin(thetac);
   r12=r1-(*tr0);
   dist=p.Dot(r12)+d*tan(thetac);
   tgeo=t0+dist/c;
   t1=sx1->GetSignal("LE",7);
   tres=t1-tgeo;

   if (tres<-30 || tres>300 || d>25.*pow(tres+30.,0.25)) continue;

   // Associate this hit to the TE
   te->AddSignal(*sx1);
   levers.Enter(d/tan(thetac));
  }
  // Quality check of the various TE's.
  // Survivors will be called track candidates (TC's)
  // and their Q quality value will be determined.
  nah=te->GetNsignals();
  sigmal=levers.GetSigma(1);
  if (nah<10 || sigmal<20)  // Remove the TE's of poor quality
  {
   temap.RemoveObjects(te);
   tes.RemoveAt(jte);
   delete te;
  }
  else // Specify the Q factor for this TC
  {
   qtc=0.3*sigmal+7.;
   if (qtc>nah) qtc=nah;
   fit.SetSignal(qtc,"QTC");
   fit.SetSignal(sigmal,"SIGMAL");
   te->SetFitDetails(fit);
   if (qtc>qmax) qmax=qtc;
  }
 }
 tes.Compress();
 nte=tes.GetEntries();

 // Perform selection on Q value in case of multiple track candidates
 for (Int_t jtc=0; jtc<nte; jtc++)
 {
  te=(AliTrack*)tes.At(jtc);
  if (!te) continue;
  sx1=(AliSignal*)te->GetFitDetails();
  if (!sx1) continue;
  qtc=sx1->GetSignal("QTC");
  if (qtc<0.7*qmax)
  {
   temap.RemoveObjects(te);
   tes.RemoveAt(jtc);
   delete te;
  }
 } 
 tes.Compress();
 nte=tes.GetEntries();

 // Exit in case no track candidates are left
 if (!nte)
 {
  if (ntes) delete [] ntes;
  return;
 }

 // Cluster track candidates within a certain opening angle into jets. 
 TObjArray jets;
 jets.SetOwner();
 AliTrack* te2=0;
 Float_t ang=0;
 for (Int_t jtc1=0; jtc1<nte; jtc1++)
 {
  te=(AliTrack*)tes.At(jtc1);
  if (!te) continue;
  AliJet* jx=new AliJet();
  jx->AddTrack(te);
  jets.Add(jx);
  for (Int_t jtc2=0; jtc2<nte; jtc2++)
  {
   te2=(AliTrack*)tes.At(jtc2);
   if (!te2) continue;
   ang=te->GetOpeningAngle(*te2,"deg");
   if (ang<fTangsep) jx->AddTrack(te2);
  }
 }

 // Order the jets w.r.t. decreasing number of tracks
 TObjArray* ordered=evt->SortJets(-1,&jets);
 TObjArray jets2(*ordered);
 Int_t njets=jets2.GetEntries();

 // Merge jets within a certain opening to provide the final track(s).
 AliJet* jx1=0;
 AliJet* jx2=0;
 if (fJangsep>0)
 {
  for (Int_t jet1=0; jet1<njets; jet1++)
  {
   jx1=(AliJet*)jets2.At(jet1);
   if (!jx1) continue;
   for (Int_t jet2=jet1+1; jet2<njets; jet2++)
   {
    jx2=(AliJet*)jets2.At(jet2);
    if (!jx2) continue;
    ang=jx1->GetOpeningAngle(*jx2,"deg");
    if (ang<fJangsep)
    {
     for (Int_t jtk=1; jtk<=jx2->GetNtracks(); jtk++)
     {
      te=jx2->GetTrack(jtk);
      if (!te) continue;
      jx1->AddTrack(te);
     }
     jets2.RemoveAt(jet2);    
    }
   }
  }
  jets2.Compress();
  njets=jets2.GetEntries();
 }

 // Store every jet as a reconstructed track in the event structure.
 // The jet 3-momentum (normalised to 1) and the average r0 and t0
 // of the constituent tracks will make up the final track parameters.
 // All the associated hits of all the constituent tracks of the jet
 // will be associated to the final track.
 // In case the jet angular separation was set <0, only the jet with
 // the maximum number of tracks (i.e. the first one in the array)
 // will be used to form a track. This will allow comparison with
 // the standard Sieglinde processing.
 AliSample pos;
 AliSample time;
 AliPosition* ref=0;
 AliTrack t; 
 t.SetNameTitle("IceDwalk","Direct walk track");
 t.SetTimestamp((AliTimestamp&)*evt);
 Float_t vec[3],err[3];
 for (Int_t jet=0; jet<njets; jet++)
 {
  AliJet* jx=(AliJet*)jets2.At(jet);
  if (!jx) continue;
  evt->AddTrack(t);
  AliTrack* trk=evt->GetTrack(evt->GetNtracks());
  if (!trk) continue;
  trk->SetId(evt->GetNtracks(1)+1);
  p=jx->Get3Momentum();
  p/=p.GetNorm();
  trk->Set3Momentum(p);
  pos.Reset();
  time.Reset();
  for (Int_t jt=1; jt<=jx->GetNtracks(); jt++)
  {
   AliTrack* tx=jx->GetTrack(jt);
   if (!tx) continue;
   AliTimestamp* tsx=tx->GetTimestamp();
   t0=evt->GetDifference(tsx,"ns");
   time.Enter(t0);
   ref=tx->GetReferencePoint();
   if (ref)
   {
    ref->GetPosition(vec,"car");
    pos.Enter(vec[0],vec[1],vec[2]);
   }
   for (Int_t is=1; is<=tx->GetNsignals(); is++)
   {
    sx1=tx->GetSignal(is);
    if (sx1) sx1->AddLink(trk);
   }
  }
  for (Int_t k=1; k<=3; k++)
  {
   vec[k-1]=pos.GetMean(k);
   err[k-1]=pos.GetSigma(k);
  }
  r0.SetPosition(vec,"car");
  r0.SetPositionErrors(err,"car");
  t0=time.GetMean(1);
  AliTimestamp* tt0=trk->GetTimestamp();
  tt0->Add(0,0,(int)t0);
  trk->SetReferencePoint(r0);

  // Only take the jet with the maximum number of tracks
  // (i.e. the first jet in the list) when the user had selected
  // this reconstruction mode.
  if (fJangsep<0) break;
 }

 if (ntes) delete [] ntes;
}
///////////////////////////////////////////////////////////////////////////
