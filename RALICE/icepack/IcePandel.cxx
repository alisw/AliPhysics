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
// Class IcePandel
// TTask derived class to perform Pandel fitting.
//
// The code in this processor is based on the algorithms as developed
// by Oladipo Fadiran and George Japaridze (Clark Atlanta University, USA).
//
// Use the UseTracks memberfunction to specify the first guess tracks
// to be processed by the minimiser.
// By default only the first encountered IceDwalk track will be processed.
//
// Use the SelectHits memberfunction to specify the hits to be used.
// By default only the hits associated to the first guess track are used.
//
// The fit processor (TFitter which is basically Minuit) printlevel
// can be selected via the memberfunction SetPrintLevel.
// By default all printout is suppressed (i.e. level=-2).
//
// An example of how to invoke this processor after Xtalk, hit cleaning
// and a direct walk first guess estimate can be found in the ROOT example
// macro icepandel.cc which resides in the /macros subdirectory.
// 
// The minimisation results are stored in the IceEvent structure as
// tracks with name "IcePandel" (just like the first guess results
// of e.g. IceDwalk).
// A pointer to the first guess track which was used as input is available
// via the GetParentTrack facility of these "IcePandel" tracks.
// Furthermore, all the hits that were used in the minisation are available
// via the GetSignal facility of a certain track.
// The statistics of the TFitter result are stored as an AliSignal object
// in the track, which can be obtained via the GetFitDetails memberfunction.
// Currently no overall probability has yet been defined to indicate the
// quality of a certain fit result. So, for the time being the user has
// to judge the fit quality his/herself by means of the various TFitter stats
// as accessible via the GetFitDetails facility.
//
// An example of how the various data can be accessed is given below,
// where "evt" indicates the pointer to the IceEvent structure.
//
// Example for accessing data :
// ----------------------------
// TObjArray* tracks=evt->GetTracks("IcePandel");
// if (!tracks) return;
// AliPosition* r0=0;
// Float_t fcn=0;
// for (Int_t jtk=0; jtk<tracks->GetEntries(); jtk++)
// {
//  AliTrack* tx=(AliTrack*)tracks->At(jtk);
//  if (!tx) continue;
//  tx->Data("sph");
//  r0=tx->GetReferencePoint();
//  if (r0) r0->Data();
//  sx=(AliSignal*)tx->GetFitDetails();
//  if (sx) fcn=sx->GetSignal("FCN");
//  AliTrack* tx2=tx->GetParentTrack();
//  if (!tx2) continue;
//  tx2->Data("sph");
//  r0=tx2->GetReferencePoint();
//  if (r0) r0->Data();
// }
//
// Notes :
// -------
// 1) This processor only works properly on data which are Time and ADC
//    calibrated and contain tracks from first guess algorithms like
//    e.g. IceDwalk.
// 2) In view of the usage of TFitter/Minuit minimisation, a global pointer
//    to the instance of this class (gIcePandel) and a global static
//    wrapper function (IcePandelFCN) have been introduced, to allow the
//    actual minimisation to be performed via the memberfunction FitFCN.
//    This implies that in a certain processing job only 1 instance of
//    this IcePandel class may occur.
//
//--- Author: Nick van Eijndhoven 09-feb-2006 Utrecht University
//- Modified: NvE $Date$ Utrecht University
///////////////////////////////////////////////////////////////////////////
 
#include "IcePandel.h"
#include "Riostream.h"

// Global pointer to the instance of this object
 IcePandel* gIcePandel=0;

// TFitter/Minuit interface to IcePandel::FitFCN
 void IcePandelFCN(Int_t& npar,Double_t* gin,Double_t& f,Double_t* u,Int_t flag)
 {
  if (gIcePandel) gIcePandel->FitFCN(npar,gin,f,u,flag);
 }

ClassImp(IcePandel) // Class implementation to enable ROOT I/O

IcePandel::IcePandel(const char* name,const char* title) : TTask(name,title)
{
// Default constructor.
 fFirst=1;
 fPrint=-2;
 fSelhits=1;
 fEvt=0;
 fUseNames=0;
 fUseNtk=0;
 fTrack=0;
 fHits=0;
 fFitter=0;

 // Set the global pointer to this instance
 gIcePandel=this;
}
///////////////////////////////////////////////////////////////////////////
IcePandel::~IcePandel()
{
// Default destructor.
 if (fUseNames)
 {
  delete fUseNames;
  fUseNames=0;
 }
 if (fUseNtk)
 {
  delete fUseNtk;
  fUseNtk=0;
 }
 if (fHits)
 {
  delete fHits;
  fHits=0;
 }
 if (fFitter)
 {
  delete fFitter;
  fFitter=0;
 }
}
///////////////////////////////////////////////////////////////////////////
void IcePandel::Exec(Option_t* opt)
{
// Implementation of the hit fitting procedure.

 TString name=opt;
 AliJob* parent=(AliJob*)(gROOT->GetListOfTasks()->FindObject(name.Data()));

 if (!parent) return;

 fEvt=(IceEvent*)parent->GetObject("IceEvent");
 if (!fEvt) return;

 if (!fUseNames) UseTracks("IceDwalk",1);

 Int_t nclasses=fUseNames->GetEntries(); // Number of first guess classes to be processed
 Int_t ntkmax=0; // Max. number of tracks for a certain class
 TObjString* strx=0;
 TString str;

 if (fFirst)
 {
  cout << " *IcePandel* First guess selections to be processed (-1=all)." << endl;
  for (Int_t i=0; i<nclasses; i++)
  {
   strx=(TObjString*)fUseNames->At(i);
   if (!strx) continue;
   str=strx->GetString();
   ntkmax=fUseNtk->At(i);
   cout << " Maximally " << ntkmax << " track(s) per event for procedure : " << str.Data() << endl;
  }
  cout << " *IcePandel* Hit selection mode : " << fSelhits << endl;
  cout << endl;
  fFirst=0;
 }

 Double_t pi=acos(-1.);

 // Initialisation of the minimisation processor
 Double_t arglist[100];
 if (!fFitter) fFitter=new TFitter();

 // The number of reconstructed tracks already present in the event
 Int_t ntkreco=fEvt->GetNtracks(1);

 if (!fHits)
 {
  fHits=new TObjArray();
 }
 else
 {
  fHits->Clear();
 }

 // If selected, use all the good quality hits of the complete event
 if (fSelhits==0)
 {
  TObjArray* hits=fEvt->GetHits("IceGOM");
  for (Int_t ih=0; ih<hits->GetEntries(); ih++)
  {
   AliSignal* sx=(AliSignal*)hits->At(ih);
   if (!sx) continue;
   if (sx->GetDeadValue("ADC") || sx->GetDeadValue("LE") || sx->GetDeadValue("TOT")) continue;
   fHits->Add(sx);
  }
 }

 // Track by track processing of the selected first guess classes
 Float_t vec[3];
 Float_t err[3];
 Float_t margin,minmarg=25; // (minimal) margin for the fitter r0 search area
 Ali3Vector p;
 AliPosition pos;
 AliTrack tkfit;
 tkfit.SetNameTitle("IcePandel","Pandel fit result");
 AliSignal fitstats;
 fitstats.SetNameTitle("Fitstats","TFitter stats for Pandel fit");
 fitstats.SetSlotName("IER",1);
 fitstats.SetSlotName("FCN",2);
 fitstats.SetSlotName("EDM",3);
 fitstats.SetSlotName("NVARS",4);
 Float_t x,y,z,r,theta,phi;
 Float_t xmin,xmax,ymin,ymax,zmin,zmax,rmin,rmax,thetamin,thetamax,phimin,phimax;
 Double_t amin,edm,errdef; // Minimisation stats
 Int_t ier,nvpar,nparx;    // Minimisation stats
 Int_t ntk=0;
 Int_t nsig=0;
 for (Int_t iclass=0; iclass<nclasses; iclass++) // Loop over first guess classes
 {
  strx=(TObjString*)fUseNames->At(iclass);
  if (!strx) continue;
  str=strx->GetString();
  ntkmax=fUseNtk->At(iclass);
  TObjArray* tracks=fEvt->GetTracks(str);
  ntk=tracks->GetEntries();
  if (ntkmax>0 && ntk>ntkmax) ntk=ntkmax;

  for (Int_t jtk=0; jtk<ntk; jtk++) // Loop over tracks of a certain class
  {
   fTrack=(AliTrack*)tracks->At(jtk);
   if (!fTrack) continue;

   AliPosition* r0=fTrack->GetReferencePoint();
   if (!r0) continue;

   // If selected, use only the first guess track associated hits
   if (fSelhits==1)
   {
    fHits->Clear();
    nsig=fTrack->GetNsignals();
    for (Int_t is=1; is<=nsig; is++)
    {
     AliSignal* sx=fTrack->GetSignal(is);
     if (!sx) continue;
     if (!sx->GetDevice()->InheritsFrom("IceGOM")) continue;
     if (sx->GetDeadValue("ADC") || sx->GetDeadValue("LE") || sx->GetDeadValue("TOT")) continue;
     fHits->Add(sx);
    }
   }
   
   if (!fHits->GetEntries()) continue;

   r0->GetVector(vec,"car");
   r0->GetErrors(err,"car");

   x=vec[0];
   y=vec[1];
   z=vec[2];
   margin=5.*err[0];
   if (margin<minmarg) margin=minmarg;
   xmin=x-margin;
   xmax=x+margin;
   margin=5.*err[1];
   if (margin<minmarg) margin=minmarg;
   ymin=y-margin;
   ymax=y+margin;
   margin=5.*err[2];
   if (margin<minmarg) margin=minmarg;
   zmin=z-margin;
   zmax=z+margin;

   p=fTrack->Get3Momentum();
   p.GetVector(vec,"sph");

   theta=vec[1];
   phi=vec[2];
   thetamin=theta-(pi/3.);
   if (thetamin<0) thetamin=0;
   thetamax=theta+(pi/3.);
   if (thetamax>pi) thetamax=pi;
   phimin=phi-(pi/2.);
   phimax=phi+(pi/2.);

   // Process this first guess track with its associated hits
   fFitter->Clear();

   arglist[0]=fPrint;
   if (fPrint==-2) arglist[0]=-1;
   fFitter->ExecuteCommand("SET PRINT",arglist,1);
   if (fPrint==-2) fFitter->ExecuteCommand("SET NOWARNINGS",arglist,0);

   fFitter->SetFitMethod("loglikelihood");

   fFitter->SetParameter(0,"r0x",x,0.001,xmin,xmax);
   fFitter->SetParameter(1,"r0y",y,0.001,ymin,ymax);
   fFitter->SetParameter(2,"r0z",z,0.001,zmin,zmax);
   fFitter->SetParameter(3,"theta",theta,0.001,thetamin,thetamax);
   fFitter->SetParameter(4,"phi",phi,0.001,phimin,phimax);

   fFitter->SetFCN(IcePandelFCN);

   arglist[0]=0;
   ier=fFitter->ExecuteCommand("MINIMIZE",arglist,0);

   fFitter->GetStats(amin,edm,errdef,nvpar,nparx);
   fitstats.Reset();
   fitstats.SetSignal(ier,1);
   fitstats.SetSignal(amin,2);
   fitstats.SetSignal(edm,3);
   fitstats.SetSignal(nvpar,4);

   // Resulting parameters after minimisation
   vec[0]=fFitter->GetParameter(0);
   vec[1]=fFitter->GetParameter(1);
   vec[2]=fFitter->GetParameter(2);
   err[0]=fFitter->GetParError(0);
   err[1]=fFitter->GetParError(1);
   err[2]=fFitter->GetParError(2);
   pos.SetPosition(vec,"car");
   pos.SetPositionErrors(err,"car");

   vec[0]=1.;
   vec[1]=fFitter->GetParameter(3);
   vec[2]=fFitter->GetParameter(4);
   err[0]=0.;
   err[1]=fFitter->GetParError(3);
   err[2]=fFitter->GetParError(4);
   p.SetVector(vec,"sph");
   p.SetErrors(err,"sph");

   // Enter the fit result as a track in the event structure
   tkfit.Reset();
   ntkreco++;
   tkfit.SetId(ntkreco);
   tkfit.SetParentTrack(fTrack);
   AliTimestamp* tt0=r0->GetTimestamp();
   pos.SetTimestamp(*tt0);
   tkfit.SetTimestamp(*tt0);
   tkfit.SetReferencePoint(pos);
   tkfit.Set3Momentum(p);
   tkfit.SetFitDetails(fitstats);
   for (Int_t ihit=0; ihit<fHits->GetEntries(); ihit++)
   {
    AliSignal* sx=(AliSignal*)fHits->At(ihit);
    if (sx) tkfit.AddSignal(*sx);
   }
   fEvt->AddTrack(tkfit);
  } // End loop over tracks
 } // End loop over first guess classes

}
///////////////////////////////////////////////////////////////////////////
void IcePandel::SetPrintLevel(Int_t level)
{
// Set the fitter (Minuit) print level.
// See the TFitter and TMinuit docs for details.
//
// Note : level=-2 suppresses also all fit processor warnings.
//        
// The default in the constructor is level=-2. 

 fPrint=level;
}
///////////////////////////////////////////////////////////////////////////
void IcePandel::UseTracks(TString classname,Int_t n)
{
// Specification of the first guess tracks to be used.
//
// classname : Specifies the first guess algorithm (e.g. "IceDwalk");
// n : Specifies the max. number of these tracks to be used
//
// Note : n<0 will use all the existing tracks of the specified classname
//
// The default is n=-1.
//
// Consecutive invokations of this memberfunction with different classnames
// will result in an incremental effect.
//
// Example :
// ---------
// UseTracks("IceDwalk",5);
// UseTracks("IceLfit",2);
// UseTracks("IceJams");
//
// This will use the first 5 IceDwalk, the first 2 IceLfit and all the
// IceJams tracks which are encountered in the event structure.

 if (!fUseNames)
 {
  fUseNames=new TObjArray();
  fUseNames->SetOwner();
 }
 
 if (!fUseNtk) fUseNtk=new TArrayI();

 // Check if this classname has already been specified before 
 TString s;
 Int_t nen=fUseNames->GetEntries();
 for (Int_t i=0; i<nen; i++)
 {
  TObjString* sx=(TObjString*)fUseNames->At(i);
  if (!sx) continue;
  s=sx->GetString();
  if (s==classname) return;
 }

 // New classname to be added into the storage
 if (nen >= fUseNames->GetSize()) fUseNames->Expand(nen+1);
 if (nen >= fUseNtk->GetSize()) fUseNtk->Set(nen+1);
 
 TObjString* name=new TObjString();
 name->SetString(classname);
 fUseNames->Add(name);
 fUseNtk->AddAt(n,nen);
}
///////////////////////////////////////////////////////////////////////////
void IcePandel::SelectHits(Int_t mode)
{
// Specification of the hits to be used in the minimisation.
//
// mode = 0 : All hit cleaning survived hits of the complete event are used
//        1 : Only the associated hits are used for each first guess track
//
// The default is mode=1.

 if (mode==0 || mode==1) fSelhits=mode;
}
///////////////////////////////////////////////////////////////////////////
void IcePandel::FitFCN(Int_t&,Double_t*,Double_t& f,Double_t* x,Int_t)
{
// The Pandel function used for the minimisation process.

 const Float_t c=0.299792;           // Light speed in vacuum in meters per ns
 const Float_t nice=1.35634;         // Refractive index of ice
 const Float_t thetac=acos(1./nice); // Cherenkov angle (in radians)
 const Float_t lambda=33.3;          // Light scattering length in ice
 const Float_t labs=98;              // Light absorbtion length in ice
 const Float_t cice=c/nice;          // Light speed in ice in meters per ns
 const Float_t tau=557;
 const Double_t rho=((1./tau)+(cice/labs));

 f=0;

 // The first original guess track parameters
 AliPosition* tr0=fTrack->GetReferencePoint();
 AliTimestamp* tt0=tr0->GetTimestamp();
 Float_t t0=fEvt->GetDifference(tt0,"ns");
 Ali3Vector tp=fTrack->Get3Momentum();

 // The new r0 and p vectors from the minimisation
 Float_t vec[3];
 
 AliPosition r0;
 vec[0]=x[0];
 vec[1]=x[1];
 vec[2]=x[2];
 r0.SetPosition(vec,"car");

 Ali3Vector p;
 vec[0]=1;
 vec[1]=x[3];
 vec[2]=x[4];
 p.SetVector(vec,"sph");

 Int_t nhits=fHits->GetEntries();
 AliPosition rhit;
 Ali3Vector r12;
 Float_t d,dist,thit,tgeo;
 Double_t tres,zeta,pandel;
 for (Int_t i=0; i<nhits; i++)
 {
  AliSignal* sx=(AliSignal*)fHits->At(i);
  if (!sx) continue;
  IceGOM* omx=(IceGOM*)sx->GetDevice();
  if (!omx) continue;
  rhit=omx->GetPosition();
  d=fTrack->GetDistance(rhit);
  zeta=d/lambda;
  d*=sin(thetac);
  r12=rhit-r0;
  dist=p.Dot(r12)+d*tan(thetac);
  tgeo=t0+dist/c;
  thit=sx->GetSignal("LE",7);
  tres=thit-tgeo;

  // The Pandel function evaluation
  // Avoid minimiser problem for infinite derivative on Pandel surface
  // This problem will be absent when using a smooth convoluted Pandel
  if (tres>0 && zeta>=1)
  {
   pandel=pow(rho,zeta)*pow(tres,(zeta-1.))*exp(-rho*tres)/TMath::Gamma(zeta);
  }
  else
  {
   pandel=1.e-6;
  }

  // Use 10*log10 expression to obtain intuitive decibel scale
  f-=10.*log10(pandel);
 }
}
///////////////////////////////////////////////////////////////////////////
