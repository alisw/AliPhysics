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
// Class AliEventSelector
// TTask based processor to perform generic event selection.
// This class is derived from AliAstrolab in order to also provide event
// selection based on space and time matching with external (astrophysical)
// objects and phenomena.
//
// After having applied the various selection criteria, this processor
// introduces an AliDevice with the name "AliEventSelector" into the event.
// This device contains named signal slots to indicate the settings
// of the various selection parameters.
// One of the slots has the name "Select" and the signal value of this
// slot indicates the final selection result.
//
// value : -1 ==> Event rejected
//          0 ==> Decision unknown (incomplete selection parameters)
//          1 ==> Event selected
//
// Event selection may be performed based on various selection types,
// e.g. individual track observables, total event observables or
// space and time matching with external objects.
// These types can be (de)activated via the SetSelector member function.
//
// The specific selection criteria for each selection type may be
// specified via the SetRange memberfunction.
//
// The logic to be used in the selection process with the various criteria
// is set via the memberfunction SetLogic.
// Obviously, matching of tracks with various external objects is always
// performed in logical "or".
//
// For investigation of individual track observables and/or matching with
// external objects, the user may define a restricted set of tracks to
// be used in the evaluation procedures. The definition of such a restricted
// track set is performed via the memberfunction UseTracks.  
//
// Example :
// ---------
// gSystem->Load("ralice");
//
// AliEventSelector sel;
// sel.SetLabPosition(0.5,-90,"deg");    // The South Pole Neutrino Observatory
// sel.SetLocalFrame(90,180,90,270,0,0); // Local frame has X-axis to the North
//
// sel.SetSelector("astro");
// sel.SetSelector("track");
// sel.SetSelector("event");
// sel.SetLogic("and");
// sel.UseTracks("*");
//
// sel.SetRange("track","p",0.5,1e20);   // Require at least 500 MeV/c track momentum
// sel.SetRange("event","ntrk","*",1,3); // Only low multiplicity events 
//
// // Match with Astrophysical objects within 5 degrees and 10 seconds
// Float_t da=5;
// Float_t dt=10;
// sel.SetAstroMatch(da,dt,"from");
//
// // Some observed event to be investigated
// AliEvent evt;
// evt.SetUT(1989,7,30,8,14,23,738504,0);
//
// Float_t vec[3]={1,23.8,118.65};
// Ali3Vector r;
// r.SetVector(vec,"sph","deg");
//
// AliTrack t;
// t.SetNameTitle("SomeTrack","Just a dummy test track");
// r*=-1; // Let track originate from specified location
// t.Set3Momentum(r);
//
// evt.AddTrack(t);
//
// // Enter some external (astrophysical) reference signals
// Float_t alpha=194818.0;
// Float_t delta=84400.;
// sel.SetSignal(alpha,delta,"B",1950,"M",-1,"Altair");
// alpha=124900.0;
// delta=272400.;
// sel.SetSignal(alpha,delta,"B",1950,"M",-1,"NGP");
// alpha=64508.917;
// delta=-164258.02;
// sel.SetSignal(alpha,delta,"J",2000,"M",-1,"Sirius");
// alpha=23149.08;
// delta=891550.8;
// sel.SetSignal(alpha,delta,"J",2000,"M",-1,"Polaris");
// alpha=43600.;
// delta=163100.;
// sel.SetSignal(alpha,delta,"J",2000,"M",-1,"Aldebaran");
// Float_t l=327.531;
// Float_t b=-35.8903;
// Float_t pos[3]={1,90.-b,l};
// r.SetVector(pos,"sph","deg");
// sel.SetUT(1989,7,30,8,14,16,0,0);
// sel.SetSignal(&r,"gal","M",0,-1,"GRB890730"); // This matches our track
//
// // List all stored reference objects
// sel.ListSignals("equ","T",5);
//
// // Let's see what the event selection makes of it
// AliJob q;
// q.Add(&sel);
// q.ListEnvironment();
// q.ProcessObject(&evt);
//
// AliDevice* evtsel=(AliDevice*)evt.GetDevice("AliEventSelector");
// if (evtsel) evtsel->Data();
//
//--- Author: Nick van Eijndhoven 17-sep-2007 Utrecht University
//- Modified: NvE $Date$ Utrecht University
///////////////////////////////////////////////////////////////////////////
 
#include "AliEventSelector.h"
#include "Riostream.h"

ClassImp(AliEventSelector) // Class implementation to enable ROOT I/O

AliEventSelector::AliEventSelector(const char* name,const char* title) : AliAstrolab(name,title)
{
// Default constructor.

 fFirst=1;
 fEvt=0;
 fParams=0;
 fSelect=0;
 fTrackflag=0;
 fEventflag=0;
 fAstroflag=0;
 fLogic=0;
 fUseNames=0;
 fUseNtk=0;
 fAstroDa=-1;
 fAstroDt=-1;
 fAstroDir=0;

 for (Int_t i=0; i<=8; i=i+2)
 {
  fEventTracks[i]=0;
  fEventTracks[i+1]=-1;
  if (i<=4)
  {
   fTrackMomenta[i]=0;
   fTrackMomenta[i+1]=-1;
   fTrackEnergies[i]=0;
   fTrackEnergies[i+1]=-1;
   fEventMomenta[i]=0;
   fEventMomenta[i+1]=-1;
   fEventEnergies[i]=0;
   fEventEnergies[i+1]=-1;
  }
  if (i<=2)
  {
   fTrackRapidities[i]=0;
   fTrackRapidities[i+1]=-1;
  }
  if (i==0)
  {
   fTrackMasses[i]=0;
   fTrackMasses[i+1]=-1;
   fEventMasses[i]=0;
   fEventMasses[i+1]=-1;
   fTrackCharges[i]=0;
   fTrackCharges[i+1]=-1;
   fEventCharges[i]=0;
   fEventCharges[i+1]=-1;
   fTrackDevices[i]=0;
   fTrackDevices[i+1]=-1;
   fEventDevices[i]=0;
   fEventDevices[i+1]=-1;
  }
 }
 fTrackDevClass="";
 fEventDevClass="";
 fEventTrkName="";
}
///////////////////////////////////////////////////////////////////////////
AliEventSelector::~AliEventSelector()
{
// Default destructor.

 if (fParams)
 {
  delete fParams;
  fParams=0;
 }

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
}
///////////////////////////////////////////////////////////////////////////
AliEventSelector::AliEventSelector(const AliEventSelector& q) : AliAstrolab(q)
{
// Copy constructor

 fFirst=1;
 fEvt=0;
 fParams=new AliDevice(*(q.fParams));
 fSelect=q.fSelect;
 fTrackflag=q.fTrackflag;
 fEventflag=q.fEventflag;
 fAstroflag=q.fAstroflag;
 fLogic=q.fLogic;
 fAstroDa=q.fAstroDa;
 fAstroDt=q.fAstroDt;
 fAstroDir=q.fAstroDir;

 for (Int_t i=0; i<10; i++)
 {
  fEventTracks[i]=q.fEventTracks[i];

  if (i<6)
  {
   fTrackMomenta[i]=q.fTrackMomenta[i];
   fTrackEnergies[i]=q.fTrackEnergies[i];
   fEventMomenta[i]=q.fEventMomenta[i];
   fEventEnergies[i]=q.fEventEnergies[i];
  }

  if (i<4) fTrackRapidities[i]=q.fTrackRapidities[i];

  if (i<2)
  {
   fTrackMasses[i]=q.fTrackMasses[i];
   fEventMasses[i]=q.fEventMasses[i];
   fTrackCharges[i]=q.fTrackCharges[i];
   fEventCharges[i]=q.fEventCharges[i];
   fTrackDevices[i]=q.fTrackDevices[i];
   fEventDevices[i]=q.fEventDevices[i];
  }
 }
 fTrackDevClass=q.fTrackDevClass;
 fEventDevClass=q.fEventDevClass;
 fEventTrkName=q.fEventTrkName;
}
///////////////////////////////////////////////////////////////////////////
void AliEventSelector::SetSelector(TString type,Int_t flag)
{
// Specify the selection types to be used.
// The various types my be selected in a cumulative way by specification
// of the input argument "type".
// The various possibilities are :
//
// type = "track" ==> Selection based on individual track observables (e.g. Pt)
//        "event" ==> Selection based on total event observables (e.g. Invmass)
//        "astro" ==> Selection based on correlation with external objects
//
// The specified selection types can be (de)activated via the input
// argument "flag".
//
// flag = 0 ==> Don't use the specified selection type
//        1 ==> Use the specified selection type
//
// For type="astro" the flag>0 value specifies further selections as follows :
//
// flag = 1 ==> Match individual track momentum directions with external (astrophysical) objects
//        2 ==> Match event total momentum direction with external (astrophysical) objects
//        3 ==> Match event position with external (astrophysical) objects
//
// For further details see memberfunction SetAstroMatch.
//
// The default value is flag=1.
//
// Note : In the default constructor all selection types are de-activated.

 if (type=="track") fTrackflag=flag;
 if (type=="event") fEventflag=flag;
 if (type=="astro") fAstroflag=flag;
}
///////////////////////////////////////////////////////////////////////////
void AliEventSelector::SetLogic(TString type)
{
// Set type of the decision logic.
//
// type = "and"  ==> Event selection based on logical "and"
//        "or"   ==> Event selection based on logical "or"
//        "nand" ==> Event selection based on logical "nand"
//        "nor"  ==> Event selection based on logical "nor"
//
// Note : In the default constructor the decision logic is set to "unknown".

 if (type=="and") fLogic=1;
 if (type=="or") fLogic=2;
 if (type=="nand") fLogic=-1;
 if (type=="nor") fLogic=-2;
}
///////////////////////////////////////////////////////////////////////////
void AliEventSelector::UseTracks(TString name,Int_t n)
{
// Specification of the track names to be used for the investigation
// of individual track observables and matching with external objects.
//
// name : Specifies the track name (e.g. "IceDwalk")
//        In case name="*" all track names will be accepted.

// n : Specifies the max. number of these tracks to be used
//
// Note : n<0 will use all the existing tracks of the specified name
//
// The default is n=-1.
//
// Consecutive invokations of this memberfunction with different names
// will result in an incremental effect.
//
// Example :
// ---------
// UseTracks("IceDwalk",5);
// UseTracks("IceLinefit",2);
// UseTracks("Pythia");
//
// This will use the first 5 IceDwalk, the first 2 IceLinefit and all the
// Pythia tracks which are encountered in the event structure.

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
  if (s==name) return;
 }

 // New name to be added into the storage
 if (nen >= fUseNames->GetSize()) fUseNames->Expand(nen+1);
 if (nen >= fUseNtk->GetSize()) fUseNtk->Set(nen+1);
 
 TObjString* namex=new TObjString();
 namex->SetString(name);
 fUseNames->Add(namex);
 fUseNtk->AddAt(n,nen);
}
///////////////////////////////////////////////////////////////////////////
void AliEventSelector::SetAstroMatch(Double_t da,Double_t dt,TString dir)
{
// Set the parameters for the matching of reference objects.
//
// da  : Maximum angular difference in degrees
// dt  : Maximum absolute time difference in seconds
// dir : "to"   ==> Check the location the track (or event) points to
//       "from" ==> Check the location the track (or event) originates from

 fAstroDa=fabs(da);
 fAstroDt=fabs(dt);
 fAstroDir=0;
 if (dir=="to") fAstroDir=1;
 if (dir=="from") fAstroDir=-1;

 if (!fParams) fParams=new AliDevice();

 fParams->AddNamedSlot("AstroDa");
 fParams->AddNamedSlot("AstroDt");
 fParams->AddNamedSlot("AstroDir");
 fParams->SetSignal(fAstroDa,"AstroDa");
 fParams->SetSignal(fAstroDt,"AstroDt");
 fParams->SetSignal(fAstroDir,"AstroDir");
}
///////////////////////////////////////////////////////////////////////////
void AliEventSelector::SetRange(TString type,TString obs,Double_t low,Double_t up)
{
// Set range for the specified observable.
//
// type : Selection type specifier (e.g. "track" or "event").
// obs  : Observable specification.
// low  : Lower bound of acceptance range
// up   : Upper bound of acceptance range
//
// The various observables that are available for selection criteria are :
//
// obs : "p"   ==> Momentum value in GeV/c
//       "pt"  ==> Transverse momentum value in GeV/c
//       "pl"  ==> Longitudinal momentum value in GeV/c
//       "e"   ==> Energy value in GeV
//       "et"  ==> Transverse momentum value in GeV
//       "el"  ==> Longitudinal momentum value in GeV
//       "m"   ==> (Invariant) mass in GeV/c^2
//       "q"   ==> Charge (electron charge is defined as -1)
//       "y"   ==> Rapidity (only for "track")
//       "eta" ==> Pseudo-rapidity (only for "track")
//
// Note : When up<low the specified observable will not be used for selection.
//
// In the default constructor all observables are de-activated for selection.

 if (!fParams) fParams=new AliDevice();

 if (type=="track") // Individual track observables
 {
  if (obs=="p")
  {
   fTrackMomenta[0]=low;
   fTrackMomenta[1]=up;
   fParams->AddNamedSlot("TrackMinP");
   fParams->AddNamedSlot("TrackMaxP");
   fParams->SetSignal(low,"TrackMinP");
   fParams->SetSignal(up,"TrackMaxP");
  }
  if (obs=="pt")
  {
   fTrackMomenta[2]=low;
   fTrackMomenta[3]=up;
   fParams->AddNamedSlot("TrackMinPt");
   fParams->AddNamedSlot("TrackMaxPt");
   fParams->SetSignal(low,"TrackMinPt");
   fParams->SetSignal(up,"TrackMaxPt");
  }
  if (obs=="pl")
  {
   fTrackMomenta[4]=low;
   fTrackMomenta[5]=up;
   fParams->AddNamedSlot("TrackMinPl");
   fParams->AddNamedSlot("TrackMaxPl");
   fParams->SetSignal(low,"TrackMinPl");
   fParams->SetSignal(up,"TrackMaxPl");
  }
  if (obs=="e")
  {
   fTrackEnergies[0]=low;
   fTrackEnergies[1]=up;
   fParams->AddNamedSlot("TrackMinE");
   fParams->AddNamedSlot("TrackMaxE");
   fParams->SetSignal(low,"TrackMinE");
   fParams->SetSignal(up,"TrackMaxE");
  }
  if (obs=="et")
  {
   fTrackEnergies[2]=low;
   fTrackEnergies[3]=up;
   fParams->AddNamedSlot("TrackMinEt");
   fParams->AddNamedSlot("TrackMaxEt");
   fParams->SetSignal(low,"TrackMinEt");
   fParams->SetSignal(up,"TrackMaxEt");
  }
  if (obs=="el")
  {
   fTrackEnergies[4]=low;
   fTrackEnergies[5]=up;
   fParams->AddNamedSlot("TrackMinEl");
   fParams->AddNamedSlot("TrackMaxEl");
   fParams->SetSignal(low,"TrackMinEl");
   fParams->SetSignal(up,"TrackMaxEl");
  }
  if (obs=="m")
  {
   fTrackMasses[0]=low;
   fTrackMasses[1]=up;
   fParams->AddNamedSlot("TrackMinM");
   fParams->AddNamedSlot("TrackMaxM");
   fParams->SetSignal(low,"TrackMinM");
   fParams->SetSignal(up,"TrackMaxM");
  }
  if (obs=="q")
  {
   fTrackCharges[0]=low;
   fTrackCharges[1]=up;
   fParams->AddNamedSlot("TrackMinQ");
   fParams->AddNamedSlot("TrackMaxQ");
   fParams->SetSignal(low,"TrackMinQ");
   fParams->SetSignal(up,"TrackMaxQ");
  }
  if (obs=="y")
  {
   fTrackRapidities[0]=low;
   fTrackRapidities[1]=up;
   fParams->AddNamedSlot("TrackMinY");
   fParams->AddNamedSlot("TrackMaxY");
   fParams->SetSignal(low,"TrackMinY");
   fParams->SetSignal(up,"TrackMaxY");
  }
  if (obs=="eta")
  {
   fTrackRapidities[2]=low;
   fTrackRapidities[3]=up;
   fParams->AddNamedSlot("TrackMinEta");
   fParams->AddNamedSlot("TrackMaxEta");
   fParams->SetSignal(low,"TrackMinEta");
   fParams->SetSignal(up,"TrackMaxEta");
  }
 }

 if (type=="event") // Total event observables
 {
  if (obs=="p")
  {
   fEventMomenta[0]=low;
   fEventMomenta[1]=up;
   fParams->AddNamedSlot("EventMinP");
   fParams->AddNamedSlot("EventMaxP");
   fParams->SetSignal(low,"EventMinP");
   fParams->SetSignal(up,"EventMaxP");
  }
  if (obs=="pt")
  {
   fEventMomenta[2]=low;
   fEventMomenta[3]=up;
   fParams->AddNamedSlot("EventMinPt");
   fParams->AddNamedSlot("EventMaxPt");
   fParams->SetSignal(low,"EventMinPt");
   fParams->SetSignal(up,"EventMaxPt");
  }
  if (obs=="pl")
  {
   fEventMomenta[4]=low;
   fEventMomenta[5]=up;
   fParams->AddNamedSlot("EventMinPl");
   fParams->AddNamedSlot("EventMaxPl");
   fParams->SetSignal(low,"EventMinPl");
   fParams->SetSignal(up,"EventMaxPl");
  }
  if (obs=="e")
  {
   fEventEnergies[0]=low;
   fEventEnergies[1]=up;
   fParams->AddNamedSlot("EventMinE");
   fParams->AddNamedSlot("EventMaxE");
   fParams->SetSignal(low,"EventMinE");
   fParams->SetSignal(up,"EventMaxE");
  }
  if (obs=="et")
  {
   fEventEnergies[2]=low;
   fEventEnergies[3]=up;
   fParams->AddNamedSlot("EventMinEt");
   fParams->AddNamedSlot("EventMaxEt");
   fParams->SetSignal(low,"EventMinEt");
   fParams->SetSignal(up,"EventMaxEt");
  }
  if (obs=="el")
  {
   fEventEnergies[4]=low;
   fEventEnergies[5]=up;
   fParams->AddNamedSlot("EventMinEl");
   fParams->AddNamedSlot("EventMaxEl");
   fParams->SetSignal(low,"EventMinEl");
   fParams->SetSignal(up,"EventMaxEl");
  }
  if (obs=="m")
  {
   fEventMasses[0]=low;
   fEventMasses[1]=up;
   fParams->AddNamedSlot("EventMinM");
   fParams->AddNamedSlot("EventMaxM");
   fParams->SetSignal(low,"EventMinM");
   fParams->SetSignal(up,"EventMaxM");
  }
  if (obs=="q")
  {
   fEventCharges[0]=low;
   fEventCharges[1]=up;
   fParams->AddNamedSlot("EventMinQ");
   fParams->AddNamedSlot("EventMaxQ");
   fParams->SetSignal(low,"EventMinQ");
   fParams->SetSignal(up,"EventMaxQ");
  }
 }
}
///////////////////////////////////////////////////////////////////////////
void AliEventSelector::SetRange(TString type,TString obs,TString name,Int_t nlow,Int_t nup)
{
// Set range for the specified observable.
//
// type  : Selection type specifier (e.g. "track" or "event").
// obs   : Observable specification.
// name  : (Class) name of the objects to be searched for
// nlow  : Lower bound of acceptance range
// nup   : Upper bound of acceptance range
//
// The various observables that are available for selection criteria are :
//
// obs : "ndev" ==> Number of associated devices of the specified (derived) class name
//       "ntrk" ==> Number of tracks with the specified name (name="*" ==> all tracks)
//       "ntkc" ==> Total number of charged tracks (no name selection)
//       "ntk0" ==> Total number of neutral tracks (no name selection)
//       "ntk+" ==> Total number of positive tracks (no name selection)
//       "ntk-" ==> Total number of negative tracks (no name selection)
//
// Note : For a certain (type,obs) combination only one (class) name can be specified.  

 if (!fParams) fParams=new AliDevice();

 if (type=="track") // Individual track observables
 {
  if (obs=="ndev")
  {
   fTrackDevices[0]=nlow;
   fTrackDevices[1]=nup;
   fTrackDevClass=name;
   fParams->AddNamedSlot("TrackMinNdev");
   fParams->AddNamedSlot("TrackMaxNdev");
   fParams->SetSignal(float(nlow),"TrackMinNdev");
   fParams->SetSignal(float(nup),"TrackMaxNdev");
  }
 }

 if (type=="event") // Total event observables
 {
  if (obs=="ndev")
  {
   fEventDevices[0]=nlow;
   fEventDevices[1]=nup;
   fEventDevClass=name;
   fParams->AddNamedSlot("EventMinNdev");
   fParams->AddNamedSlot("EventMaxNdev");
   fParams->SetSignal(float(nlow),"EventMinNdev");
   fParams->SetSignal(float(nup),"EventMaxNdev");
  }
  if (obs=="ntrk")
  {
   fEventTracks[0]=nlow;
   fEventTracks[1]=nup;
   fEventTrkName=name;
   fParams->AddNamedSlot("EventMinNtrk");
   fParams->AddNamedSlot("EventMaxNtrk");
   fParams->SetSignal(float(nlow),"EventMinNtrk");
   fParams->SetSignal(float(nup),"EventMaxNtrk");
  }
  if (obs=="ntkc")
  {
   fEventTracks[2]=nlow;
   fEventTracks[3]=nup;
   fParams->AddNamedSlot("EventMinNtkc");
   fParams->AddNamedSlot("EventMaxNtkc");
   fParams->SetSignal(float(nlow),"EventMinNtkc");
   fParams->SetSignal(float(nup),"EventMaxNtkc");
  }
  if (obs=="ntk0")
  {
   fEventTracks[4]=nlow;
   fEventTracks[5]=nup;
   fParams->AddNamedSlot("EventMinNtk0");
   fParams->AddNamedSlot("EventMaxNtk0");
   fParams->SetSignal(float(nlow),"EventMinNtk0");
   fParams->SetSignal(float(nup),"EventMaxNtk0");
  }
  if (obs=="ntk+")
  {
   fEventTracks[6]=nlow;
   fEventTracks[7]=nup;
   fParams->AddNamedSlot("EventMinNtk+");
   fParams->AddNamedSlot("EventMaxNtk+");
   fParams->SetSignal(float(nlow),"EventMinNtk+");
   fParams->SetSignal(float(nup),"EventMaxNtk+");
  }
  if (obs=="ntk-")
  {
   fEventTracks[8]=nlow;
   fEventTracks[9]=nup;
   fParams->AddNamedSlot("EventMinNtk-");
   fParams->AddNamedSlot("EventMaxNtk-");
   fParams->SetSignal(float(nlow),"EventMinNtk-");
   fParams->SetSignal(float(nup),"EventMaxNtk-");
  }
 }
}
///////////////////////////////////////////////////////////////////////////
void AliEventSelector::Exec(Option_t* opt)
{
// Implementation of the event selection procedures.

 TString name=opt;
 AliJob* parent=(AliJob*)(gROOT->GetListOfTasks()->FindObject(name.Data()));

 if (!parent) return;

 fEvt=(AliEvent*)parent->GetObject("AliEvent");
 if (!fEvt) return;

 Int_t ntkmax=0; // Max. number of tracks for a certain name
 Int_t nnames=0; // Number of track names to be processed
 if (fUseNames) nnames=fUseNames->GetEntries();
 TObjString* strx=0;
 TString str;

 if (fFirst)
 {
  cout << " *AliEventSelector* Selection parameters." << endl;
  cout << " Selection types in use :";
  if (fTrackflag) cout << " track";
  if (fEventflag) cout << " event";
  if (fAstroflag) cout << " astro";
  if (!fTrackflag && !fEventflag && !fAstroflag) cout << " none";
  cout << endl;
  cout << " Selection logic in use :";
  if (fLogic==1) cout << " and";
  if (fLogic==2) cout << " or";
  if (fLogic==-1) cout << " nand";
  if (fLogic==-2) cout << " nor";
  if (!fLogic) cout << " unknown";
  cout << endl;
  if (nnames) cout << " Track name selections to be processed (-1=all)." << endl;
  for (Int_t i=0; i<nnames; i++)
  {
   strx=(TObjString*)fUseNames->At(i);
   if (!strx) continue;
   str=strx->GetString();
   ntkmax=fUseNtk->At(i);
   cout << " Maximally " << ntkmax << " track(s) per event of name : " << str.Data() << endl;
  }
  cout << endl;

  fFirst=0;
 }

 // Storage of the used parameters in the AliEventSelector device
 if (!fParams) fParams=new AliDevice();
 fParams->SetNameTitle("AliEventSelector","AliEventSelector processor parameters");

 fSelect=0;
 if (fLogic)
 {
  if (fEventflag) Event();  // Check criteria for total event observables
  if (fTrackflag) Track(0); // Check criteria for total track observables
  if (fAstroflag) Astro();  // Check for matches with external objects
 }

 if (fLogic<0) fSelect*=-1; // In case of "nand"/"nor" logic

 fParams->AddNamedSlot("Logic");
 fParams->SetSignal(float(fLogic),"Logic");
 fParams->AddNamedSlot("Eventflag");
 fParams->SetSignal(float(fEventflag),"Eventflag");
 fParams->AddNamedSlot("Trackflag");
 fParams->SetSignal(float(fTrackflag),"Trackflag");
 fParams->AddNamedSlot("Astroflag");
 fParams->SetSignal(float(fAstroflag),"Astroflag");
 fParams->AddNamedSlot("Select");
 fParams->SetSignal(float(fSelect),"Select");

 AliDevice* dx=(AliDevice*)fEvt->GetDevice("AliEventSelector");
 if (dx) fEvt->RemoveDevice(dx);
 fEvt->AddDevice(fParams);
}
///////////////////////////////////////////////////////////////////////////
void AliEventSelector::Track(Int_t mode)
{
// Check criteria for individual track observables.
// This memberfunction serves also the track direction checking
// for external (astrophysical) objects.
// mode = 0 : Track observables (e.g. P, Pt etc...) are checked
//        1 : Track direction is checked w.r.t. external (astrophysical) objects

 if (abs(fLogic)==1) fSelect=-1; // Selections are made in logical "(n)and"

 if (fSelect>0) return; // Event is already flagged as select

 if (!fUseNames) return;

 if (mode==1 && !fAstroDir) return;

 Int_t nnames=fUseNames->GetEntries(); // Number of track names to be processed
 Int_t ntkmax=0; // Max. number of tracks for a certain name
 TObjString* strx=0;
 TString str;
 Int_t ntk;

 AliTrack* track=0;
 AliTimestamp* ts=0;
 Ali3Vector p;
 TObjArray* tracks=0;
 Double_t val;
 Int_t ival;
 for (Int_t i=0; i<nnames; i++) // Loop over selected track names
 {
  strx=(TObjString*)fUseNames->At(i);
  if (!strx) continue;
  str=strx->GetString();
  ntkmax=fUseNtk->At(i);
  if (str=="*")
  {
   tracks=fEvt->GetTracks();
  }
  else
  {
   tracks=fEvt->GetTracks(str);
  }
  ntk=0;
  if (tracks) ntk=tracks->GetEntries();
  if (ntkmax>0 && ntk>ntkmax) ntk=ntkmax;

  for (Int_t jtk=0; jtk<ntk; jtk++) // Loop over tracks of a certain name
  {
   track=(AliTrack*)tracks->At(jtk);
   if (!track) continue;

   if (!mode) // Check track observables
   {
    if (fTrackMomenta[1]>fTrackMomenta[0]) // Selection on P
    {
     if (abs(fLogic)==1) fSelect=-1; // Selections are made in logical "(n)and"
     val=track->GetMomentum(1);
     if (val>=fTrackMomenta[0] && val<=fTrackMomenta[1])
     {
      fSelect=1;
      if (abs(fLogic)==2) return; // Selections are made in logical "(n)or"
     }
    }
    if (fTrackMomenta[3]>fTrackMomenta[2]) // Selection on Pt
    {
     if (abs(fLogic)==1) fSelect=-1; // Selections are made in logical "(n)and"
     val=track->GetPt(1);
     if (val>=fTrackMomenta[2] && val<=fTrackMomenta[3])
     {
      fSelect=1;
      if (abs(fLogic)==2) return; // Selections are made in logical "(n)or"
     }
    }
    if (fTrackMomenta[5]>fTrackMomenta[4]) // Selection on Pl
    {
     if (abs(fLogic)==1) fSelect=-1; // Selections are made in logical "(n)and"
     val=track->GetPl(1);
     if (val>=fTrackMomenta[4] && val<=fTrackMomenta[5])
     {
      fSelect=1;
      if (abs(fLogic)==2) return; // Selections are made in logical "(n)or"
     }
    }
    if (fTrackEnergies[1]>fTrackEnergies[0]) // Selection on E
    {
     if (abs(fLogic)==1) fSelect=-1; // Selections are made in logical "(n)and"
     val=track->GetEnergy(1);
     if (val>=fTrackEnergies[0] && val<=fTrackEnergies[1])
     {
      fSelect=1;
      if (abs(fLogic)==2) return; // Selections are made in logical "(n)or"
     }
    }
    if (fTrackEnergies[3]>fTrackEnergies[2]) // Selection on Et
    {
     if (abs(fLogic)==1) fSelect=-1; // Selections are made in logical "(n)and"
     val=track->GetEt(1);
     if (val>=fTrackEnergies[2] && val<=fTrackEnergies[3])
     {
      fSelect=1;
      if (abs(fLogic)==2) return; // Selections are made in logical "(n)or"
     }
    }
    if (fTrackEnergies[5]>fTrackEnergies[4]) // Selection on El
    {
     if (abs(fLogic)==1) fSelect=-1; // Selections are made in logical "(n)and"
     val=track->GetEl(1);
     if (val>=fTrackEnergies[4] && val<=fTrackEnergies[5])
     {
      fSelect=1;
      if (abs(fLogic)==2) return; // Selections are made in logical "(n)or"
     }
    }
    if (fTrackMasses[1]>fTrackMasses[0]) // Selection on M
    {
     if (abs(fLogic)==1) fSelect=-1; // Selections are made in logical "(n)and"
     val=track->GetMass(1);
     if (val>=fTrackMasses[0] && val<=fTrackMasses[1])
     {
      fSelect=1;
      if (abs(fLogic)==2) return; // Selections are made in logical "(n)or"
     }
    }
    if (fTrackCharges[1]>fTrackCharges[0]) // Selection on Q
    {
     if (abs(fLogic)==1) fSelect=-1; // Selections are made in logical "(n)and"
     val=track->GetCharge();
     if (val>=fTrackCharges[0] && val<=fTrackCharges[1])
     {
      fSelect=1;
      if (abs(fLogic)==2) return; // Selections are made in logical "(n)or"
     }
    }
    if (fTrackRapidities[1]>fTrackRapidities[0]) // Selection on Y
    {
     if (abs(fLogic)==1) fSelect=-1; // Selections are made in logical "(n)and"
     val=track->GetRapidity();
     if (val>=fTrackRapidities[0] && val<=fTrackRapidities[1])
     {
      fSelect=1;
      if (abs(fLogic)==2) return; // Selections are made in logical "(n)or"
     }
    }
    if (fTrackRapidities[3]>fTrackRapidities[2]) // Selection on Eta
    {
     if (abs(fLogic)==1) fSelect=-1; // Selections are made in logical "(n)and"
     val=track->GetPseudoRapidity();
     if (val>=fTrackRapidities[2] && val<=fTrackRapidities[3])
     {
      fSelect=1;
      if (abs(fLogic)==2) return; // Selections are made in logical "(n)or"
     }
    }
    if (fTrackDevices[1]>fTrackDevices[0]) // Selection on Ndev
    {
     if (abs(fLogic)==1) fSelect=-1; // Selections are made in logical "(n)and"
     ival=track->GetNsignals(fTrackDevClass.Data());
     if (ival>=fTrackDevices[0] && ival<=fTrackDevices[1])
     {
      fSelect=1;
      if (abs(fLogic)==2) return; // Selections are made in logical "(n)or"
     }
    }
   }
   if (mode==1) // Check track direction w.r.t. external (astrophysical) objects
   {
    p=track->Get3Momentum();
    if (fAstroDir<0) p*=-1;
    ts=track->GetTimestamp();
    if (!ts) ts=(AliTimestamp*)fEvt;

    SetSignal(&p,"loc","T",ts,0,"Track");
    TArrayI* arr=MatchRefSignal(fAstroDa,"deg",fAstroDt,"s");
    if (arr)
    {
     fSelect=1;
     return;
    }
   }
  } // End of loop over tracks of a certain name
 } // End of loop over selected track names 
}
///////////////////////////////////////////////////////////////////////////
void AliEventSelector::Event()
{
// Check criteria for total event observables.

 if (abs(fLogic)==1) fSelect=-1; // Selections are made in logical "(n)and"

 if (fSelect>0) return; // Event is already flagged as select

 Double_t val=0;
 Int_t ival=0;
 if (fEventMomenta[1]>fEventMomenta[0]) // Selection on P
 {
  if (abs(fLogic)==1) fSelect=-1; // Selections are made in logical "(n)and"
  val=fEvt->GetMomentum(1);
  if (val>=fEventMomenta[0] && val<=fEventMomenta[1])
  {
   fSelect=1;
   if (abs(fLogic)==2) return; // Selections are made in logical "(n)or"
  }
 }
 if (fEventMomenta[3]>fEventMomenta[2]) // Selection on Pt
 {
  if (abs(fLogic)==1) fSelect=-1; // Selections are made in logical "(n)and"
  val=fEvt->GetPt(1);
  if (val>=fEventMomenta[2] && val<=fEventMomenta[3])
  {
   fSelect=1;
   if (abs(fLogic)==2) return; // Selections are made in logical "(n)or"
  }
 }
 if (fEventMomenta[5]>fEventMomenta[4]) // Selection on Pl
 {
  if (abs(fLogic)==1) fSelect=-1; // Selections are made in logical "(n)and"
  val=fEvt->GetPl(1);
  if (val>=fEventMomenta[4] && val<=fEventMomenta[5])
  {
   fSelect=1;
   if (abs(fLogic)==2) return; // Selections are made in logical "(n)or"
  }
 }
 if (fEventEnergies[1]>fEventEnergies[0]) // Selection on E
 {
  if (abs(fLogic)==1) fSelect=-1; // Selections are made in logical "(n)and"
  val=fEvt->GetEnergy(1);
  if (val>=fEventEnergies[0] && val<=fEventEnergies[1])
  {
   fSelect=1;
   if (abs(fLogic)==2) return; // Selections are made in logical "(n)or"
  }
 }
 if (fEventEnergies[3]>fEventEnergies[2]) // Selection on Et
 {
  if (abs(fLogic)==1) fSelect=-1; // Selections are made in logical "(n)and"
  val=fEvt->GetEt(1);
  if (val>=fEventEnergies[2] && val<=fEventEnergies[3])
  {
   fSelect=1;
   if (abs(fLogic)==2) return; // Selections are made in logical "(n)or"
  }
 }
 if (fEventEnergies[5]>fEventEnergies[4]) // Selection on El
 {
  if (abs(fLogic)==1) fSelect=-1; // Selections are made in logical "(n)and"
  val=fEvt->GetEl(1);
  if (val>=fEventEnergies[4] && val<=fEventEnergies[5])
  {
   fSelect=1;
   if (abs(fLogic)==2) return; // Selections are made in logical "(n)or"
  }
 }
 if (fEventMasses[1]>fEventMasses[0]) // Selection on Minv
 {
  if (abs(fLogic)==1) fSelect=-1; // Selections are made in logical "(n)and"
  val=fEvt->GetInvmass(1);
  if (val>=fEventMasses[0] && val<=fEventMasses[1])
  {
   fSelect=1;
   if (abs(fLogic)==2) return; // Selections are made in logical "(n)or"
  }
 }
 if (fEventCharges[1]>fEventCharges[0]) // Selection on Q
 {
  if (abs(fLogic)==1) fSelect=-1; // Selections are made in logical "(n)and"
  val=fEvt->GetCharge();
  if (val>=fEventCharges[0] && val<=fEventCharges[1])
  {
   fSelect=1;
   if (abs(fLogic)==2) return; // Selections are made in logical "(n)or"
  }
 }
 if (fEventDevices[1]>fEventDevices[0]) // Selection on Ndev
 {
  if (abs(fLogic)==1) fSelect=-1; // Selections are made in logical "(n)and"
  ival=fEvt->GetNdevices(fEventDevClass.Data());
  if (ival>=fEventDevices[0] && ival<=fEventDevices[1])
  {
   fSelect=1;
   if (abs(fLogic)==2) return; // Selections are made in logical "(n)or"
  }
 }
 if (fEventTracks[1]>fEventTracks[0]) // Selection on Ntrk
 {
  if (abs(fLogic)==1) fSelect=-1; // Selections are made in logical "(n)and"
  if (fEventTrkName=="*")
  {
   ival=fEvt->GetNtracks();
  }
  else
  {
   ival=fEvt->GetNtracks(fEventTrkName);
  }
  if (ival>=fEventTracks[0] && ival<=fEventTracks[1])
  {
   fSelect=1;
   if (abs(fLogic)==2) return; // Selections are made in logical "(n)or"
  }
 }
 if (fEventTracks[3]>fEventTracks[2]) // Selection on Ntkc
 {
  if (abs(fLogic)==1) fSelect=-1; // Selections are made in logical "(n)and"
  ival=fEvt->GetNtracks(0,3,0);
  if (ival>=fEventTracks[2] && ival<=fEventTracks[3])
  {
   fSelect=1;
   if (abs(fLogic)==2) return; // Selections are made in logical "(n)or"
  }
 }
 if (fEventTracks[5]>fEventTracks[4]) // Selection on Ntk0
 {
  if (abs(fLogic)==1) fSelect=-1; // Selections are made in logical "(n)and"
  ival=fEvt->GetNtracks(0,0,0);
  if (ival>=fEventTracks[4] && ival<=fEventTracks[5])
  {
   fSelect=1;
   if (abs(fLogic)==2) return; // Selections are made in logical "(n)or"
  }
 }
 if (fEventTracks[7]>fEventTracks[6]) // Selection on Ntk+
 {
  if (abs(fLogic)==1) fSelect=-1; // Selections are made in logical "(n)and"
  ival=fEvt->GetNtracks(0,1,0);
  if (ival>=fEventTracks[6] && ival<=fEventTracks[7])
  {
   fSelect=1;
   if (abs(fLogic)==2) return; // Selections are made in logical "(n)or"
  }
 }
 if (fEventTracks[9]>fEventTracks[8]) // Selection on Ntk-
 {
  if (abs(fLogic)==1) fSelect=-1; // Selections are made in logical "(n)and"
  ival=fEvt->GetNtracks(0,-1,0);
  if (ival>=fEventTracks[8] && ival<=fEventTracks[9])
  {
   fSelect=1;
   if (abs(fLogic)==2) return; // Selections are made in logical "(n)or"
  }
 }
}
///////////////////////////////////////////////////////////////////////////
void AliEventSelector::Astro()
{
// Check for matches with external objects.

 if (abs(fLogic)==1) fSelect=-1; // Selections are made in logical "(n)and"

 if (fSelect>0) return; // Event is already flagged as select

 // Check track directions w.r.t. external (astrophysical) objects
 if (fAstroflag==1)
 {
  Track(1);
  return;
 }

 // Check total event momentum direction w.r.t. external (astrophysical) objects
 if (fAstroflag==2)
 {
  Ali3Vector p;
  p=fEvt->Get3Momentum();
  if (fAstroDir<0) p*=-1;
  SetSignal(&p,"loc","T",(AliTimestamp*)fEvt,0,"Event");
  TArrayI* arr=MatchRefSignal(fAstroDa,"deg",fAstroDt,"s");
  if (arr && fLogic) fSelect=1;
  return;
 }

 // Check event position w.r.t. external (astrophysical) objects
 if (fAstroflag==3)
 {
  SetSignal((Ali3Vector*)fEvt,"loc","T",(AliTimestamp*)fEvt,0,"Event");
  TArrayI* arr=MatchRefSignal(fAstroDa,"deg",fAstroDt,"s");
  if (arr && fLogic) fSelect=1;
  return;
 }
}
///////////////////////////////////////////////////////////////////////////
TObject* AliEventSelector::Clone(const char* name) const
{
// Make a deep copy of the current object and provide the pointer to the copy.
// This memberfunction enables automatic creation of new objects of the
// correct type depending on the object type, a feature which may be very useful
// for containers when adding objects in case the container owns the objects.

 AliEventSelector* sel=new AliEventSelector(*this);
 if (name)
 {
  if (strlen(name)) sel->SetName(name);
 }
 return sel;
}
///////////////////////////////////////////////////////////////////////////
