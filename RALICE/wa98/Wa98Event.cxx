// $Id$

///////////////////////////////////////////////////////////////////////////
// Class Wa98Event
// Creation and investigation of a Wa98 physics event.
// This event class is derived from AliEvent and has some Wa98 specific
// extensions like e.g. the information from the trigger calorimeters
// and some LEDA specific functions.
//
//--- Author: Nick van Eijndhoven 24-apr-2002 UU-SAP Utrecht
//- Modified: NvE $Date$ UU-SAP Utrecht
///////////////////////////////////////////////////////////////////////////

#include "Wa98Event.h"
#include "Riostream.h"
 
ClassImp(Wa98Event) // Class implementation to enable ROOT I/O
 
Wa98Event::Wa98Event() : AliEvent()
{
// Default constructor.
// All variables initialised to default values.
 Reset();
}
///////////////////////////////////////////////////////////////////////////
Wa98Event::Wa98Event(Int_t n) : AliEvent(n)
{
// Create an event to hold initially a maximum of n tracks
// All variables initialised to default values
 Reset();
}
///////////////////////////////////////////////////////////////////////////
Wa98Event::~Wa98Event()
{
// Default destructor
}
///////////////////////////////////////////////////////////////////////////
Wa98Event::Wa98Event(Wa98Event& evt) : AliEvent(evt)
{
// Copy constructor.
 fTrig=evt.fTrig;
 fWeight=evt.fWeight;
 fZdc=evt.fZdc;
 fEmir=evt.fEmir;
 fEmire=evt.fEmire;
 fEmirh=evt.fEmirh;
 fEtm=evt.fEtm;
 fEtme=evt.fEtme;
 fEtmh=evt.fEtmh;
}
///////////////////////////////////////////////////////////////////////////
void Wa98Event::Reset()
{
// Reset all variables to default values
// The max. number of tracks is set to the initial value again
// The max. number of vertices is set to the default value again
 fTrig=0;
 fWeight=0;
 fZdc=0;
 fEmir=0;
 fEmire=0;
 fEmirh=0;
 fEtm=0;
 fEtme=0;
 fEtmh=0;

 AliEvent::Reset();
}
///////////////////////////////////////////////////////////////////////////
void Wa98Event::SetTrig(Int_t trig)
{
// Set the trigger class.
// Trigger classes : 1=nsc 3=cen 5=per 6=mbias 7=beam 8=inbeam ped
 fTrig=trig;
}
///////////////////////////////////////////////////////////////////////////
Int_t Wa98Event::GetTrig()
{
// Provide the trigger class.
// Trigger classes : 1=nsc 3=cen 5=per 6=mbias 7=beam 8=inbeam ped
 return fTrig;
}
///////////////////////////////////////////////////////////////////////////
void Wa98Event::SetWeight(Int_t w)
{
// Set the event weight to account for the downscale factor.
 fWeight=w;
}
///////////////////////////////////////////////////////////////////////////
Int_t Wa98Event::GetWeight()
{
// Provide the event weight factor to account for the DAQ downscaling.
 return fWeight;
}
///////////////////////////////////////////////////////////////////////////
void Wa98Event::SetZdc(Float_t zdc)
{
// Set the ZDC signal in GeV.
 fZdc=zdc;
}
///////////////////////////////////////////////////////////////////////////
Float_t Wa98Event::GetZdc()
{
// Provide the ZDC signal in GeV.
 return fZdc;
}
///////////////////////////////////////////////////////////////////////////
void Wa98Event::SetMiracE(Float_t tot,Float_t em,Float_t had)
{
// Set the total, EM and hadronic signals for MIRAC in GeV.
 fEmir=tot;
 fEmire=em;
 fEmirh=had;
}
///////////////////////////////////////////////////////////////////////////
Float_t Wa98Event::GetEmir()
{
// Provide the total MIRAC signal in GeV.
 return fEmir;
}
///////////////////////////////////////////////////////////////////////////
Float_t Wa98Event::GetEmire()
{
// Provide the MIRAC EM signal in GeV.
 return fEmire;
}
///////////////////////////////////////////////////////////////////////////
Float_t Wa98Event::GetEmirh()
{
// Provide the MIRAC hadronic signal in GeV.
 return fEmirh;
}
///////////////////////////////////////////////////////////////////////////
void Wa98Event::SetMiracEt(Float_t tot,Float_t em,Float_t had)
{
// Set the total, EM and hadronic Et signals for MIRAC in GeV.
 fEtm=tot;
 fEtme=em;
 fEtmh=had;
}
///////////////////////////////////////////////////////////////////////////
Float_t Wa98Event::GetEtm()
{
// Provide the total MIRAC Et signal in GeV.
 return fEtm;
}
///////////////////////////////////////////////////////////////////////////
Float_t Wa98Event::GetEtme()
{
// Provide the MIRAC EM Et signal in GeV.
 return fEtme;
}
///////////////////////////////////////////////////////////////////////////
Float_t Wa98Event::GetEtmh()
{
// Provide the MIRAC hadronic Et signal in GeV.
 return fEtmh;
}
///////////////////////////////////////////////////////////////////////////
void Wa98Event::InitLeda(AliCalorimeter* cal)
{
// Set module positions and flag bad modules for this LEDA part. 
 if (cal)
 {
  Int_t nrows=cal->GetNrows();
  if (nrows==40)
  {
   SetPositionsLedalw(cal);
   SetBadModulesLedalw(cal);
  } 
  if (nrows==44)
  {
   SetPositionsLedaup(cal);
   SetBadModulesLedaup(cal);
  } 
  if (nrows!=40 && nrows!=44)
  {
   cout << " *Wa98Event::InitLeda* Not a LEDA configuration. nrows = " << nrows << endl;
  } 
 }
 else
 {
  cout << " *Wa98Event::InitLeda* Calorimeter pointer was zero. " << endl;
 }
}
///////////////////////////////////////////////////////////////////////////
void Wa98Event::SetBadModulesLedaup(AliCalorimeter* cal)
{
// Marking of the bad modules the upper LEDA 

 Int_t nr=cal->GetNrows();
 Int_t nc=cal->GetNcolumns();

// Declare the non-existing modules as dead
 Int_t row0=0;
 for (Int_t col=1; col<=nc; col++)
 {
  if (col<19 || col>126) row0=41;
  if (col>36 && col<109) row0=41;
  if (col<13 || col>132) row0=37;
  if (col>42 && col<103) row0=37;
  if (col< 7 || col>138) row0=33;
  if (col>48 && col< 97) row0=33;

  for (Int_t row=row0; row<=nr; row++)
  {
   cal->SetDead(row,col);
  }
 }

 // The bad area as seen in the first 50 events of Pb96 run 9066
 for (Int_t i=24; i<=33; i++)
 {
  for (Int_t j=1; j<=20; j++)
  {
   cal->SetDead(i,j);
  }
 }
}
//////////////////////////////////////////////////////////////////////////
void Wa98Event::SetBadModulesLedalw(AliCalorimeter* cal)
{
// Marking of the bad modules the lower LEDA 

 Int_t nr=cal->GetNrows();
 Int_t nc=cal->GetNcolumns();

// Declare the non-existing modules as dead
 Int_t row0=0;
 for (Int_t col=1; col<=nc; col++)
 {
  if (col<19 || col>126) row0=37;
  if (col>36 && col<109) row0=37;
  if (col<13 || col>132) row0=33;
  if (col>42 && col<103) row0=33;
  if (col< 7 || col>138) row0=29;
  if (col>48 && col< 97) row0=29;

  for (Int_t row=row0; row<=nr; row++)
  {
   cal->SetDead(row,col);
  }
 }
}
//////////////////////////////////////////////////////////////////////////
void Wa98Event::SetPositionsLedaup(AliCalorimeter* cal)
{
// Determination of the lab. position of each module of the upper LEDA 

 Float_t sx=4.105;  // X-dimension of a module in cm
 Float_t sy=4.085;  // Y-dimension of a module in cm
 Float_t tilt=8.31; // Tilt angle in degrees

 Float_t pi=acos(-1.);
 Float_t tiltr=tilt*pi/180.; // Tilt angle in radians

 Float_t dx=sx;            // X-displacement skipping 1 module horizontal
 Float_t dy=sy*cos(tiltr); // Y-displacement skipping 1 module vertical
 Float_t dz=sy*sin(tiltr); // Z-displacement skipping 1 module vertical

 // Position of the left upper module (1,1) looking downstream
 Float_t x0=293.;  
 Float_t y0=345.73;
 Float_t z0=2180.12;

 // Determine and store position of each module
 Int_t nr=cal->GetNrows();
 Int_t nc=cal->GetNcolumns();
 Float_t pos[3];
 for (Int_t i=1; i<=nr; i++)
 {
  for (Int_t j=1; j<=nc; j++)
  {
   pos[0]=x0-float(j-1)*dx;
   pos[1]=y0-float(i-1)*dy;
   pos[2]=z0+float(i-1)*dz;
   cal->SetPosition(i,j,pos,"car");
  }
 }
}
//////////////////////////////////////////////////////////////////////////
void Wa98Event::SetPositionsLedalw(AliCalorimeter* cal)
{
// Determination of the lab. position of each module of the lower LEDA 

 Float_t sx=4.105;  // X-dimension of a module in cm
 Float_t sy=4.085;  // Y-dimension of a module in cm
 Float_t tilt=8.31; // Tilt angle in degrees

 Float_t pi=acos(-1.);
 Float_t tiltr=tilt*pi/180.; // Tilt angle in radians

 Float_t dx=sx;            // X-displacement skipping 1 module horizontal
 Float_t dy=sy*cos(tiltr); // Y-displacement skipping 1 module vertical
 Float_t dz=sy*sin(tiltr); // Z-displacement skipping 1 module vertical

 // Position of the left lower module (1,1) looking downstream
 Float_t x0=291.91;  
 Float_t y0=-331.14;
 Float_t z0=2183.14;

 // Determine and store position of each module
 Int_t nr=cal->GetNrows();
 Int_t nc=cal->GetNcolumns();
 Float_t pos[3];
 for (Int_t i=1; i<=nr; i++)
 {
  for (Int_t j=1; j<=nc; j++)
  {
   pos[0]=x0-float(j-1)*dx;
   pos[1]=y0+float(i-1)*dy;
   pos[2]=z0+float(i-1)*dz;
   cal->SetPosition(i,j,pos,"car");
  }
 }
}
//////////////////////////////////////////////////////////////////////////
void Wa98Event::ClusterLeda(AliCalorimeter* cal,Int_t n,Int_t mode)
{
// Group LEDA modules into clusters.
// The parameter n indicates the number of rings for the grouping process
// (default n=2) and the parameter mode indicates the sorting algorithm
// (default mode=1). See AliCalorimeter::Group() for further details.
// This function invokes AliCalorimeter::Group(n) and automatically
// sets the uncertainties on the cluster positions. 
// The precision of a cluster position in the X-Y plane has been seen
// to be about half the size of a module, so dx=dy=2cm.

 Float_t err[3]={2,2,0};
 AliCalcluster* c=0;

 if (cal)
 {
  cal->Group(n,mode);
  for (Int_t i=1; i<=cal->GetNclusters(); i++)
  {
   c=cal->GetCluster(i);
   if (c) c->SetPositionErrors(err,"car");
  }
 }
 else
 {
  cout << " *Wa98Event::ClusterLeda* Calorimeter pointer was zero. " << endl;
 }
}
//////////////////////////////////////////////////////////////////////////
void Wa98Event::VetoLeda(AliCalorimeter* cal,Float_t dtheta,Float_t dphi)
{
// Perform LEDA-Veto cluster association.
// An association is only made if the Veto signal appears into a small
// cone around the LEDA cluster position.
// The cone dimensions are given by the parameters "dtheta" and "dphi"
// (in degrees) and the default values are dtheta=0.2 and dphi=1.
// The automatic straight line extrapolation of the Veto signal position
// onto the LEDA surface is used in the association procedure.

 Float_t pi=acos(-1.);

 AliCalcluster* c;
 AliSignal* v;
 Float_t posc[3],posv[3];
 Float_t dth,dph;

 Int_t nc=cal->GetNclusters();
 Int_t nv=cal->GetNvetos();
 if (nv)
 {
  for (Int_t i=1; i<=nc; i++)
  {
   c=cal->GetCluster(i);
   c->GetPosition(posc,"sph");

   for (Int_t j=1; j<=nv; j++)
   {
    v=cal->GetVetoSignal(j);
    v->GetPosition(posv,"sph");
    // Associate a close-by veto to this LEDA cluster
    dth=(posc[1]-posv[1])*180./pi;
    dph=(posc[2]-posv[2])*180./pi;
    if (fabs(dth)<dtheta && fabs(dph)<dphi) c->AddVetoSignal(v);
   }
  }
 }
}
///////////////////////////////////////////////////////////////////////////
