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

// $Id$

///////////////////////////////////////////////////////////////////////////
// Class AliHelix
// Representation and extrapolation of AliTracks in a magnetic field.
//
// This class is meant to provide a means to display and extrapolate
// AliTrack objects in the presence of a constant homogeneous magnetic field. 
//
// For track/event displays the line width, colour etc... can be set using the
// standard facilities (see TAttLine).
// By default the linewith is set to 2 and the colour set to -1 in the constructor.
// The latter results in an automatic colour coding according to the track charge
// with the convention positive=red neutral=green negative=blue.
//
// To indicate the track starting point, the memberfunction SetMarker()
// may be used.
// By default no marker will be displayed. 
//
// Examples :
// ==========
//
// Display and extrapolation of individual tracks 
// ----------------------------------------------
// Float_t vec[3];
// AliPosition r1;
// Ali3Vector p;
// AliTrack t;
//
// vec[0]=0;
// vec[1]=0;
// vec[2]=0;
// r1.SetVector(vec,"car");
//
// vec[0]=1;
// vec[1]=0;
// vec[2]=0.3;
// p.SetVector(vec,"car");
//
// t.Set3Momentum(p);
// t.SetBeginPoint(r1);
// t.SetCharge(-1);
// t.SetMass(0.139);
//
// // The magnetic field vector in Tesla
// Ali3Vector b;
// vec[0]=0;
// vec[1]=0;
// vec[2]=1;
// b.SetVector(vec,"car");
//
// AliHelix* helix=new AliHelix();
// helix->SetB(b);
// helix->SetTofmax(1e-7);
//
// TCanvas* c1=new TCanvas("c1","c1");
// TView* view=new TView(1);
// view->SetRange(-1000,-1000,-1000,1000,1000,1000);
// view->ShowAxis();
//
// // Track displays 
// Double_t range[2]={0,600};
// helix->Display(&t,range,3);
// t.SetCharge(-t.GetCharge());
// helix->Display(&t);
//
// // Track extrapolation
// Double_t pars[3]={550,0.001,3};
// AliPosition* rext=helix->Extrapolate(&t,pars);
// if (rext) rext->Data();
// ======================================================================
//
// Online display of events generated via AliCollider
// -------------------------------------------------- 
// Int_t nevents=5;   // Number of events to be generated
// Int_t jrun=1;      // The run number of this batch of generated events
//
// cout << " ***" << endl;
// cout << " *** AliCollider run for " << nevents << " events." << endl; 
// cout << " ***" << endl;
//
// AliCollider* gen=new AliCollider();
//
// gen->OpenFortranFile(6,"dump.log");
//
// gen->SetVertexMode(2);
// gen->SetResolution(1e-4);
//
// gen->SetRunNumber(jrun);
// gen->SetPrintFreq(1);
//
// gen->SetSpectatorPmin(0.01);
//
// Int_t zp=1;
// Int_t ap=1;
// Int_t zt=2;
// Int_t at=4;
//
// gen->Init("fixt",zp,ap,zt,at,158);
//
// AliHelix* helix=new AliHelix();
// Float_t vec[3]={0,2,0};
// Ali3Vector b;
// b.SetVector(vec,"car");
// helix->SetB(b);
//
// helix->Refresh(-1); // Refresh display after each event
//
// TCanvas* c1=new TCanvas("c1","c1");
// TView* view=new TView(1);
// view->SetRange(-200,-200,-200,200,200,200);
// view->ShowAxis();
//
// // Prepare random number sequence for this run
// // to obtain the number of participants for each event
// AliRandom rndm(abs(jrun));
// Float_t* rans=new Float_t[nevents];
// rndm.Uniform(rans,nevents,2,ap+at);
// Int_t npart=0;
// Int_t ntk=0;
// for (Int_t i=0; i<nevents; i++)
// {
//  npart=rans[i];
//  gen->MakeEvent(npart);
//  AliEvent* evt=gen->GetEvent();
//  if (evt)
//  {
//   helix->Display(evt);
//   c1->Update();
//   gSystem->Sleep(5000); // Some delay to keep the display on screen
//  }
// }
// ======================================================================
//
//--- Author: Nick van Eijndhoven 17-jun-2004 Utrecht University
//- Modified: NvE $Date$ Utrecht University
///////////////////////////////////////////////////////////////////////////

#include "AliHelix.h"
#include "Riostream.h"
 
ClassImp(AliHelix) // Class implementation to enable ROOT I/O
 
AliHelix::AliHelix() : THelix()
{
// Default constructor
 fRefresh=0;
 fCurves=0;
 fExt=0;
 fTofmax=1e-8;
 fMstyle=-1;
 fMsize=0;
 fMcol=0;
 fEnduse=1;

 fLineWidth=2;
 fLineColor=-1;
}
///////////////////////////////////////////////////////////////////////////
AliHelix::~AliHelix()
{
// Destructor to delete dynamically allocated memory.
 if (fCurves)
 {
  delete fCurves;
  fCurves=0;
 }
 if (fExt)
 {
  delete fExt;
  fExt=0;
 }
}
///////////////////////////////////////////////////////////////////////////
AliHelix::AliHelix(const AliHelix& h) : THelix(h)
{
// Copy constructor
 fB=h.fB;
 fRefresh=h.fRefresh;
 fTofmax=h.fTofmax;
 fMstyle=h.fMstyle;
 fMsize=h.fMsize;
 fMcol=h.fMcol;
 fEnduse=h.fEnduse;
}
///////////////////////////////////////////////////////////////////////////
void AliHelix::SetB(Ali3Vector& b)
{
// Set the magnetic field vector in Tesla.
 fB=b;

 if (fB.GetNorm()>0)
 {
  Double_t axis[3];
  fB.GetVector(axis,"car");
  SetAxis(axis);
 }
}
///////////////////////////////////////////////////////////////////////////
Ali3Vector& AliHelix::GetB()
{
// Provide the magnetic field vector in Tesla.
 return fB;
}
///////////////////////////////////////////////////////////////////////////
void AliHelix::SetTofmax(Float_t tof)
{
// Set the maximum time of flight for straight tracks in seconds.
// This maximum tof will be used for drawing etc... in case no begin
// and endpoints can be determined from the track info.
// Notes :
// -------
// 1) In case the user specifies an explicit range, it will override
//    the maximum tof limit.
// 2) By default the tofmax is set to 10 ns in the AliHelix constructor.
 fTofmax=tof;
}
///////////////////////////////////////////////////////////////////////////
Float_t AliHelix::GetTofmax() const
{
// Provide the maximum time of flight for straight tracks in seconds.
 return fTofmax;
}
///////////////////////////////////////////////////////////////////////////
void AliHelix::SetMarker(Int_t style,Float_t size,Int_t col)
{
// Specify the marker (style, size and colour) to indicate the starting point
// of a track in a display.
// In case col<0 the marker will have the same color as the track itself.
// 
// Defaults are style=8, size=0.2 and col=-1.
 
 fMstyle=style;
 fMsize=size;
 fMcol=col;
}
///////////////////////////////////////////////////////////////////////////
void AliHelix::UseEndPoint(Int_t mode)
{
// Select usage of track endpoint in drawing and extrapolation.
// This allows correct event displays even for very long tracks.
//
// mode = 0 : Do not use the track endpoint
//        1 : Use the track endpoint
// 
// The default value is mode=1 (which is also set in the constructor).

 if (mode==0 || mode==1) fEnduse=mode; 
}
///////////////////////////////////////////////////////////////////////////
void AliHelix::MakeCurve(AliTrack* t,Double_t* range,Int_t iaxis,Double_t scale)
{
// Make the helix curve for the specified AliTrack.
// Detailed information of all the helix points can be obtained via the
// GetN() and GetP() memberfunctions of TPolyLine3D.
// In case one wants to display or extrapolate an AliTrack it is preferable
// to use the Display() or Extrapolate() memberfunctions.
// It is assumed that the track charge is stored in elementary units
// (i.e. charge=1 for a proton).
// The input argument "scale" specifies the unit scale for the various
// locations where scale=0.01 indicates unit scales in cm etc...
// In case scale<=0, the unit scale for locations is determined from the
// begin, reference or endpoint of the track. If neither of these
// positions is present, all locations are assumed to be given in meter.
// The lower and upper bounds for the range are specified by range[0] and
// range[1] and the argument "iaxis" indicates along which axis this range
// is specified.
// The range can be specified either in the LAB frame or in the Helix frame.
// The latter is the frame in which the Z axis points in the B direction.
//
// The conventions for the "iaxis" argument are the following :
// iaxis = 1 ==> X axis in the LAB frame
//         2 ==> Y axis in the LAB frame
//         3 ==> Z axis in the LAB frame
//        -1 ==> X axis in the Helix frame
//        -2 ==> Y axis in the Helix frame
//        -3 ==> Z axis in the Helix frame
//
// In case range=0 the begin/end/reference points of the AliTrack and the
// maximum time of flight (see the SetTofmax() memberfunction) will be used
// and an appropriate choice for the iaxis parameter will be made automatically
// based on the track kinematics.
// In case the reference point is not present, the begin or endpoint will be used
// as reference point for the 3-momentum specification. If neither of these positions
// is present, (0,0,0) will be taken as the reference point.
// 
// The default values are range=0, iaxis=3 and scale=-1.

 SetPolyLine(0); // Reset the polyline data points

 if (!t || (range && !iaxis)) return;

 Double_t energy=t->GetEnergy(1); // Track energy in GeV
 Double_t betanorm=t->GetBeta();

 if (energy<=0 || betanorm<=0) return;

 AliPosition* rbeg=t->GetBeginPoint();
 AliPosition* rend=0;
 if (fEnduse) rend=t->GetEndPoint();
 AliPosition* rref=t->GetReferencePoint();

 // Magnetic field vector or default Z-direction
 Double_t bvec[3]={0,0,1};
 if (fB.GetNorm()>0) fB.GetVector(bvec,"car");

 // The unit scale for locations if not specified by the user
 if (scale<=0)
 {
  scale=1; // Set default to meter
  if (rbeg)
  {
   scale=rbeg->GetUnitScale();
  }
  else if (rend)
  {
   scale=rend->GetUnitScale();
  }
  else if (rref)
  {
   scale=rref->GetUnitScale();
  }
 }

 Double_t c=2.99792458e8/scale; // Lightspeed in the selected unit scale

 // The helix angular frequency
 Double_t w=9e7*(t->GetCharge()*fB.GetNorm())/energy;

 // The particle velocity in the LAB frame
 Ali3Vector beta=t->GetBetaVector();
 Ali3Vector v=beta*c;
 Double_t vel[3];
 v.GetVector(vel,"car");

 // The particle velocity in the Helix frame
 Ali3Vector betaprim=beta.GetPrimed(fRotMat);
 v=v.GetPrimed(fRotMat);
 Double_t velprim[3];
 v.GetVector(velprim,"car");

 // Check compatibility of velocity and range specification.
 if (range)
 {
  Double_t betavec[3];
  if (iaxis>0) beta.GetVector(betavec,"car");
  if (iaxis<0) betaprim.GetVector(betavec,"car");
  if (fabs(betavec[abs(iaxis)-1])/betanorm<1e-10) return;
 }

 // The LAB location in which the velocity of the particle is defined
 Double_t loc[3]={0,0,0};
 Ali3Vector rx;
 Double_t scalex=-1;
 if (rref)
 {
  rx=(Ali3Vector)(*rref);
  scalex=rref->GetUnitScale();
 }
 else if (rbeg)
 {
  rx=(Ali3Vector)(*rbeg);
  scalex=rbeg->GetUnitScale();
 }
 else if (rend)
 {
  rx=(Ali3Vector)(*rend);
  scalex=rend->GetUnitScale();
 }

 if (scalex>0 && (scalex/scale>1.1 || scale/scalex>1.1)) rx*=scalex/scale;
 rx.GetVector(loc,"car");

 // Initialisation of Helix kinematics
 SetHelix(loc,vel,w,0,kUnchanged,bvec);

 Int_t bend=0;
 if (fabs(w)>0 && fabs(fVt)>0) bend=1;

 // Flight time boundaries.
 // The time origin t=0 is chosen to indicate the position in which
 // the particle velocity was defined.
 // The total flight time is initialised to the (user specified) tofmax.
 Double_t tmin=0,tmax=0;
 Double_t tof=fTofmax;
 Double_t dum=0;

 // The trajectory begin and end points
 Double_t vec1[3]={0,0,0};
 Double_t vec2[3]={0,0,0};
 Ali3Vector r1;
 Ali3Vector r2;
 Double_t scale1=1;
 Double_t scale2=1;

 if (!bend)
 {
  ////////////////////////////////////////
  // Treatment of straight trajectories //
  ////////////////////////////////////////
  Ali3Vector r;
  if (range) // Specified range allows for exact flight time boundaries
  {
   if (iaxis>0)
   {
    tmin=(range[0]-loc[iaxis-1])/vel[iaxis-1];
    tmax=(range[1]-loc[iaxis-1])/vel[iaxis-1];
   }
   else
   {
    loc[0]=fX0;
    loc[1]=fY0;
    loc[2]=fZ0;
    tmin=(range[0]-loc[abs(iaxis)-1])/velprim[abs(iaxis)-1];
    tmax=(range[1]-loc[abs(iaxis)-1])/velprim[abs(iaxis)-1];
   }
   if (tmax<tmin)
   {
    dum=tmin;
    tmin=tmax;
    tmax=dum;
   }
   // Make the 'curve' in the LAB frame and exit.
   // Use the parametrisation : r(t)=r0+t*v
   // using the range based flight time boundaries.
   // An additional point in the middle of the trajectory is
   // generated in view of accuracy in the case of extrapolations.
   tof=tmax-tmin;
   v=beta*c;
   r1=rx;
   r=v*tmin;
   r1=r1+r;
   r1.GetVector(vec1,"car");
   SetNextPoint(float(vec1[0]),float(vec1[1]),float(vec1[2]));
   r=v*(tof/2.);
   r2=r1+r;
   r2.GetVector(vec2,"car");
   SetNextPoint(float(vec2[0]),float(vec2[1]),float(vec2[2]));
   r=v*tof;
   r2=r1+r;
   r2.GetVector(vec2,"car");
   SetNextPoint(float(vec2[0]),float(vec2[1]),float(vec2[2]));
  }
  else // Automatic range determination
  {
   // Initially the point with Z=0 in the Helix frame is taken as a starting point.
   // In case this point can't be reached, the point in which the particle velocity
   // was defined is taken as the starting point.
   // The endpoint is initially obtained by applying the tofmax from the start point.
   tmin=0;
   if (fabs(fVz)>0) tmin=-fZ0/fVz;
   v=beta*c;
   r1=rx;
   r=v*tmin;
   r1=r1+r;

   // Override the initial begin and endpoint settings by the track data
   if (rbeg)
   {
    r1=(Ali3Vector)(*rbeg); 
    scale1=rbeg->GetUnitScale();
    // All coordinates in the selected unit scale
    if (scale1/scale>1.1 || scale/scale1>1.1) r1*=scale1/scale;
   }

   r=v*fTofmax;
   r2=r1+r;
   if (rend)
   {
    r2=(Ali3Vector)(*rend); 
    scale2=rend->GetUnitScale();
    // All coordinates in the selected unit scale
    if (scale2/scale>1.1 || scale/scale2>1.1) r2*=scale2/scale;
   }
   
   r1.GetVector(vec1,"car");
   r2.GetVector(vec2,"car");

   // Make the 'curve' in the LAB frame and exit.
   SetNextPoint(float(vec1[0]),float(vec1[1]),float(vec1[2]));
   SetNextPoint(float(vec2[0]),float(vec2[1]),float(vec2[2]));
  }
 }
 else
 {
  //////////////////////////////////////
  // Treatment of curved trajectories //
  //////////////////////////////////////

  // Initialisation of the flight time boundaries.
  // Based on the constant motion of the particle along the Helix Z-axis,
  // the parametrisation z(t)=z0+fVz*t in the Helix frame is used. 
  // If possible the point with Z=0 in the Helix frame is taken as a starting point.
  // In case this point can't be reached, the point in which the particle velocity
  // was defined is taken as the starting point.
  tmin=0;
  if (fabs(fVz)>0) tmin=-fZ0/fVz;
  tmax=tmin+fTofmax;

  if (tmax<tmin)
  {
   dum=tmin;
   tmin=tmax;
   tmax=dum;
  }

  // Determination of the range in the helix frame

  if (!range) // Automatic range determination
  {
   scale1=1;
   scale2=1;
   if (rbeg)
   {
    r1=rbeg->GetPrimed(fRotMat);
    scale1=rbeg->GetUnitScale();
    // All coordinates in the selected unit scale
    if (scale1/scale>1.1 || scale/scale1>1.1) r1*=scale1/scale;
    // Re-calculate the tmin for this new starting point
    r1.GetVector(vec1,"car");
    if (fabs(fVz)>0) tmin=(vec1[2]-fZ0)/fVz;
    tmax=tmin+fTofmax;
   }
   if (rend)
   {
    r2=rend->GetPrimed(fRotMat);
    scale2=rend->GetUnitScale();
    // All coordinates in the selected unit scale
    if (scale2/scale>1.1 || scale/scale2>1.1) r2*=scale2/scale;
    r2.GetVector(vec2,"car");
    if (fabs(fVz)>0) tmax=(vec2[2]-fZ0)/fVz;
   }
   // Make the curve on basis of the flight time boundaries and exit
   if (tmax<tmin)
   {
    dum=tmin;
    tmin=tmax;
    tmax=dum;
   }
   SetRange(tmin,tmax,kHelixT);
  }
  else // User explicitly specified range
  {
   vec1[abs(iaxis)-1]=range[0];
   vec2[abs(iaxis)-1]=range[1];
   r1.SetVector(vec1,"car");
   r2.SetVector(vec2,"car");
   if (iaxis>0) // Range specified in LAB frame
   {
    r1=r1.GetPrimed(fRotMat);
    r1.GetVector(vec1,"car");
    r2=r2.GetPrimed(fRotMat);
    r2.GetVector(vec2,"car");
   } 
   // Determination of the axis component with the
   // largest range difference
   Double_t dmax=0;
   Int_t imax=0;
   Double_t test=0;
   for (Int_t i=0; i<3; i++)
   {
    test=fabs(vec1[i]-vec2[i]);
    if (test>dmax)
    {
     dmax=test;
     imax=i;
    }
   }

   Double_t rmin=vec1[imax];
   Double_t rmax=vec2[imax];
   if (rmax<rmin)
   {
    dum=rmin;
    rmin=rmax;
    rmax=dum;
   }

   // The kinematic range boundaries in the helix frame
   Double_t xmin=fX0-fVt/fW;
   Double_t xmax=fX0+fVt/fW;
   Double_t ymin=fY0-fVt/fW;
   Double_t ymax=fY0+fVt/fW;

   if (xmax<xmin)
   {
    dum=xmin;
    xmin=xmax;
    xmax=dum;
   }
   if (ymax<ymin)
   {
    dum=ymin;
    ymin=ymax;
    ymax=dum;
   }

   // Set the range for the helix
   if (imax==2 && dmax>0) SetRange(rmin,rmax,kHelixZ);
   if (imax==1)
   {
    // Limit range to kinematic boundaries if needed
    if (rmin<=ymin) rmin=ymin+1e-6*dmax;
    if (rmax>=ymax) rmax=ymax-1e-6*dmax;
    if (rmin<rmax) SetRange(rmin,rmax,kHelixY);
   }
   if (imax==0)
   {
    // Limit range to kinematic boundaries if needed
    if (rmin<=xmin) rmin=xmin+1e-6*dmax;
    if (rmax>=xmax) rmax=xmax-1e-6*dmax;
    if (rmin<rmax) SetRange(rmin,rmax,kHelixX);
   }
  }
 }
 return;
}
///////////////////////////////////////////////////////////////////////////
void AliHelix::Display(AliTrack* t,Double_t* range,Int_t iaxis,Double_t scale)
{
// Display the helix curve of an AliTrack.
// Various curves can be displayed together or individually; please refer to
// the memberfunction Refresh() for further details.
// It is assumed that the track charge is stored in elementary units
// (i.e. charge=1 for a proton).
// The input argument "scale" specifies the unit scale for the various
// locations where scale=0.01 indicates unit scales in cm etc...
// In case scale<=0, the unit scale for locations is determined from the
// begin, reference or endpoint of the track. If neither of these
// positions is present, all locations are assumed to be given in meter.
// The lower and upper bounds for the range are specified by range[0] and
// range[1] and the argument "iaxis" indicates along which axis this range
// is specified.
// The range can be specified either in the LAB frame or in the Helix frame.
// The latter is the frame in which the Z axis points in the B direction.
//
// The conventions for the "iaxis" argument are the following :
// iaxis = 1 ==> X axis in the LAB frame
//         2 ==> Y axis in the LAB frame
//         3 ==> Z axis in the LAB frame
//        -1 ==> X axis in the Helix frame
//        -2 ==> Y axis in the Helix frame
//        -3 ==> Z axis in the Helix frame
//
// In case range=0 the begin/end/reference points of the AliTrack and the
// maximum time of flight (see the SetTofmax() memberfunction) will be used
// and an appropriate choice for the iaxis parameter will be made automatically
// based on the track kinematics.
// In case the reference point is not present, the begin or endpoint will be used
// as reference point for the 3-momentum specification. If neither of these positions
// is present, (0,0,0) will be taken as the reference point.
// 
// The default values are range=0, iaxis=3 and scale=-1.
//
// Note :
// ------
// Before any display activity, a TCanvas and a TView have to be initiated
// first by the user like for instance
// 
// TCanvas* c1=new TCanvas("c1","c1");
// TView* view=new TView(1);
// view->SetRange(-1000,-1000,-1000,1000,1000,1000);
// view->ShowAxis();
//
// The user can also use the 3D viewing facilities from the TCanvas menu
// to open an appropriate view. 

 if (!t || (range && !iaxis)) return;

 MakeCurve(t,range,iaxis,scale);

 if (fRefresh>0) Refresh(fRefresh);

 Int_t np=GetLastPoint()+1;
 if (!np) return;

 Float_t* points=GetP();
 TPolyLine3D* curve=new TPolyLine3D(np,points);

 curve->SetLineWidth(fLineWidth);
 if (fLineColor<0)
 {
  Float_t q=t->GetCharge();
  curve->SetLineColor(kGreen);
  if (q>0) curve->SetLineColor(kRed);
  if (q<0) curve->SetLineColor(kBlue);
 }
 else
 {
  curve->SetLineColor(fLineColor);
 }
 curve->Draw();

 if (!fCurves)
 {
  fCurves=new TObjArray();
  fCurves->SetOwner();
 }
 fCurves->Add(curve);

 // Display the marker for the track starting point
 if (fMstyle>0)
 {
  TPolyMarker3D* m=new TPolyMarker3D();
  m->SetPoint(0,points[0],points[1],points[2]);
  m->SetMarkerStyle(fMstyle);
  m->SetMarkerSize(fMsize);
  Int_t col=curve->GetLineColor();
  if (fMcol>0) col=fMcol;
  m->SetMarkerColor(col);
  m->Draw();
  fCurves->Add(m);
 }
}
///////////////////////////////////////////////////////////////////////////
void AliHelix::Refresh(Int_t mode)
{
// Refresh the display screen before showing the next curve.
//
// mode = 0 : refreshing fully under user control.
//        1 : the display screen will be refreshed automatically
//            at each individual track display.
//       -1 : the display screen will be refreshed automatically
//            at each event display.
//
// The default is mode=0.

 if (abs(mode)<2) fRefresh=mode;
 if (fCurves) fCurves->Delete();
}
///////////////////////////////////////////////////////////////////////////
void AliHelix::Display(AliEvent* evt,Double_t* range,Int_t iaxis,Double_t scale)
{
// Display the helix curves of all tracks of the specified event.
// Various events can be displayed together or individually; please refer to
// the memberfunction Refresh() for further details.
// Please refer to the track display memberfunction for further details
// on the input arguments.
// 
// The default values are range=0, iaxis=3 and scale=-1.
//
// Note :
// ------
// Before any display activity, a TCanvas and a TView have to be initiated
// first by the user like for instance
// 
// TCanvas* c1=new TCanvas("c1","c1");
// TView* view=new TView(1);
// view->SetRange(-1000,-1000,-1000,1000,1000,1000);
// view->ShowAxis();
//
// The user can also use the 3D viewing facilities from the TCanvas menu
// to open an appropriate view. 

 if (!evt) return;

 if (fRefresh<0) Refresh(fRefresh);

 Int_t ntk=evt->GetNtracks();
 for (Int_t jtk=1; jtk<=ntk; jtk++)
 {
  AliTrack* tx=evt->GetTrack(jtk);
  if (tx) Display(tx,range,iaxis,scale);
 }
}
///////////////////////////////////////////////////////////////////////////
void AliHelix::Display(TObjArray* arr,Double_t* range,Int_t iaxis,Double_t scale)
{
// Display the helix curves of all tracks in the specified array.
// A convenient way to obtain an array with selected tracks from e.g. an AliEvent
// is to make use of its GetTracks() selection facility.
// Various arrays can be displayed together or individually; please refer to
// the memberfunction Refresh() for further details.
// Please refer to the track display memberfunction for further details
// on the input arguments.
// 
// The default values are range=0, iaxis=3 and scale=-1.
//
// Note :
// ------
// Before any display activity, a TCanvas and a TView have to be initiated
// first by the user like for instance
// 
// TCanvas* c1=new TCanvas("c1","c1");
// TView* view=new TView(1);
// view->SetRange(-1000,-1000,-1000,1000,1000,1000);
// view->ShowAxis();
//
// The user can also use the 3D viewing facilities from the TCanvas menu
// to open an appropriate view. 

 if (!arr) return;

 Int_t ntk=arr->GetEntries();
 for (Int_t jtk=0; jtk<ntk; jtk++)
 {
  TObject* obj=arr->At(jtk);
  if (!obj) continue;
  if (!(obj->InheritsFrom("AliTrack"))) continue;
  AliTrack* tx=(AliTrack*)obj;
  Display(tx,range,iaxis,scale);
 }
}
///////////////////////////////////////////////////////////////////////////
AliPosition* AliHelix::Extrapolate(AliTrack* t,Double_t* pars,Double_t scale)
{
// Extrapolate an AliTrack according to the corresponding helix curve
// and provide a pointer to the impact position w.r.t. a specified plane.
// In case the track can never reach the specified plane, the returned
// position pointer is zero.
// Detailed information of all the helix points used in the extrapolation
// can be obtained via the GetN() and GetP() memberfunctions of TPolyLine3D.
// It is assumed that the track charge is stored in elementary units
// (i.e. charge=1 for a proton).
// The input argument "scale" specifies the unit scale for the various
// locations where scale=0.01 indicates unit scales in cm etc...
// In case scale<=0, the unit scale for locations is determined from the
// begin, reference or endpoint of the track. If neither of these
// positions is present, all locations are assumed to be given in meter.
// The extrapolation parameters for the impact plane and required accuracy
// are specified by pars[0], pars[1] and pars[2], respectively.
// pars[0] = coordinate value of the plane for the impact point
// pars[1] = required accuracy on the specified impact plane coordinate
// pars[2] = the axis along which the value of par[0] is specified
//
// The parameters can be specified either w.r.t. the LAB frame or the Helix frame.
// The latter is the frame in which the Z axis points in the B direction.
//
// The conventions for the par[2] argument are the following :
// par[2] = 1 ==> X axis in the LAB frame
//          2 ==> Y axis in the LAB frame
//          3 ==> Z axis in the LAB frame
//         -1 ==> X axis in the Helix frame
//         -2 ==> Y axis in the Helix frame
//         -3 ==> Z axis in the Helix frame
//
// Example :
// ---------
// To obtain an extrapolation to the plane Z=0 in the LAB frame
// with an accuracy of 0.001 cm the input arguments would be
// pars[0]=0  pars[1]=0.001  pars[2]=3  scale=0.01
//
// Note : The default value for the scale is -1.

 if (fExt)
 {
  delete fExt;
  fExt=0;
 }

 if (!t || !pars) return fExt;

 AliPosition* rbeg=t->GetBeginPoint();
 AliPosition* rend=t->GetEndPoint();
 AliPosition* rref=t->GetReferencePoint();

 // The unit scale for locations if not specified by the user
 if (scale<=0)
 {
  scale=1; // Set default to meter
  if (rbeg)
  {
   scale=rbeg->GetUnitScale();
  }
  else if (rend)
  {
   scale=rend->GetUnitScale();
  }
  else if (rref)
  {
   scale=rref->GetUnitScale();
  }
 }

 Double_t range[2];
 range[0]=pars[0]-fabs(pars[1])/2.;
 range[1]=pars[0]+fabs(pars[1])/2.;

 Int_t iaxis=int(pars[2]);

 MakeCurve(t,range,iaxis,scale);

 Int_t np=GetLastPoint()+1;
 if (!np) return fExt;

 Float_t* points=GetP();

 // First point of the curve around the impact
 Int_t ip=0;
 Float_t first[3]={points[3*ip],points[3*ip+1],points[3*ip+2]};

 // Last point of the curve around the impact
 ip=GetLastPoint();
 Float_t last[3]={points[3*ip],points[3*ip+1],points[3*ip+2]};

 // The accuracy on the impact point 
 Float_t err[3];
 err[0]=fabs(first[0]-last[0]);
 err[1]=fabs(first[1]-last[1]);
 err[2]=fabs(first[2]-last[2]);

 // Take the middle point as impact location
 ip=np/2;
 Float_t imp[3]={points[3*ip],points[3*ip+1],points[3*ip+2]};

 fExt=new AliPosition();
 fExt->SetUnitScale(scale);
 fExt->SetPosition(imp,"car");
 fExt->SetPositionErrors(err,"car");

 return fExt;
}
///////////////////////////////////////////////////////////////////////////
