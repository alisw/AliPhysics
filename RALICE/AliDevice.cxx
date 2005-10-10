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
// Class AliDevice
// Signal (Hit) handling of a generic device.
// Basically this class provides a user interface to group and handle
// various instances of AliSignal objects, called generically "hits".
// An AliDevice object itself has (in addition to hit storage) also the
// complete functionality of the class AliSignal.
//
// Example :
// =========
//
// AliDevice m;
// // Set user defined status word to indicate e.g. readout electronics version 
// m.SetStatus(100201);
// m.SetHitCopy(1);
// m.SetName("OM123");
//
// Float_t pos[3]={1,2,3};
// m.SetPosition(pos,"car");
//
// AliSignal s;
//
// s.Reset(1);
// s.SetName("OM123 Hit 1");
// s.SetSlotName("ADC");
// s.SetSignal(10);
// s.SetSlotName("LE",2);
// s.SetSignal(-100,2);
// s.SetSlotName("TOT",3);
// s.SetSignal(-1000,3);
// m.AddHit(s);
//
// s.Reset(1);
// s.SetName("OM123 Hit 2");
// s.SetSlotName("ADC");
// s.SetSignal(11);
// s.SetSlotName("LE",2);
// s.SetSignal(-101,2);
// s.SetSlotName("TOT",3);
// s.SetSignal(1001,3);
// m.AddHit(s);
//
// s.Reset(1);
// s.SetName("OM123 Hit 3");
// s.SetSlotName("ADC");
// s.SetSignal(12);
// s.SetSlotName("LE",2);
// s.SetSignal(-102,2);
// s.SetSlotName("TOT",3);
// s.SetSignal(-1002,3);
// m.AddHit(s);
//
// TObjArray* ordered=m.SortHits("TOT");
// nhits=ordered->GetEntries();
// for (Int_t i=0; i<nhits; i++)
// {
//  AliSignal* sx=(AliSignal*)ordered->At(i);
//  if (sx) sx->Data();
// }
//
//--- Author: Nick van Eijndhoven 23-jun-2004 Utrecht University
//- Modified: NvE $Date$ Utrecht University
///////////////////////////////////////////////////////////////////////////

#include "AliDevice.h"
#include "Riostream.h"
 
ClassImp(AliDevice) // Class implementation to enable ROOT I/O
 
AliDevice::AliDevice() : AliSignal()
{
// Default constructor.
// The user definable status word is set to zero.
// By default private copies of the recorded hits will be made.
// This implies that by default the device will own the registered hits.
// See the SetHitCopy() memberfunction for further details.
 fStatus=0;
 fHitCopy=1;
 fHits=0;
 fOrdered=0;
 fMarkers=0;
}
///////////////////////////////////////////////////////////////////////////
AliDevice::~AliDevice()
{
// Default destructor.

 // Remove backward links to this device from the hits
 // which were not owned by it.
 if (!fHitCopy)
 {
  for (Int_t ih=1; ih<=GetNhits(); ih++)
  {
   AliSignal* sx=GetHit(ih);
   if (sx) sx->ResetLinks(this);
  }
 }

 if (fHits)
 {
  delete fHits;
  fHits=0;
 }

 if (fOrdered)
 {
  delete fOrdered;
  fOrdered=0;
 }

 if (fMarkers)
 {
  delete fMarkers;
  fMarkers=0;
 }
}
///////////////////////////////////////////////////////////////////////////
AliDevice::AliDevice(const AliDevice& dev) : AliSignal(dev)
{
// Copy constructor.

 fHits=0;
 fOrdered=0;
 fMarkers=0;

 fStatus=dev.GetStatus();
 fHitCopy=dev.GetHitCopy();

 Int_t nhits=dev.GetNhits();
 if (nhits)
 {
  fHits=new TObjArray(nhits);
  if (fHitCopy) fHits->SetOwner();
  for (Int_t ih=1; ih<=nhits; ih++)
  {
   AliSignal* sx=dev.GetHit(ih);
   if (fHitCopy)
   {
    fHits->Add(sx->Clone());
    AliSignal* s=(AliSignal*)fHits->Last();
    s->ResetLinks((AliDevice*)&dev);
    s->SetDevice(this);
   }
   else
   {
    sx->AddLink(this);
    fHits->Add(sx);
   }
  }
 }
}
///////////////////////////////////////////////////////////////////////////
void AliDevice::Reset(Int_t mode)
{
// Reset registered hits and AliSignal attributes.
// Note : The status word and HitCopy flag are NOT modified.
//        Use SetStatus() and SetHitCopy() to modify these parameters. 
// See AliSignal::Reset() for further details.
 RemoveHits();
 AliSignal::Reset(mode);
}
///////////////////////////////////////////////////////////////////////////
void AliDevice::SetHitCopy(Int_t j)
{
// (De)activate the creation of private copies of the AliSignals added as hits.
// j=0 ==> No private copies are made; pointers of original hits are stored.
// j=1 ==> Private copies of the hits are made and these pointers are stored.
//
// Note : Once the storage contains pointer(s) to hit(s) one cannot
//        change the HitCopy mode anymore.
//        To change the HitCopy mode for an existing AliDevice containing
//        hits one first has to invoke either RemoveHits() or Reset().
 if (!fHits)
 {
  if (j==0 || j==1)
  {
   fHitCopy=j;
  }
  else
  {
   cout << "*AliDevice::SetHitCopy* Invalid argument : " << j << endl;
  }
 }
 else
 {
  cout << "*AliDevice::SetHitCopy* Storage already contained hits."
       << "  ==> HitCopy mode not changed." << endl; 
 }
}
///////////////////////////////////////////////////////////////////////////
Int_t AliDevice::GetHitCopy() const
{
// Provide value of the HitCopy mode.
// 0 ==> No private copies are made; pointers of original hits are stored.
// 1 ==> Private copies of the hits are made and these pointers are stored.
 return fHitCopy;
}
///////////////////////////////////////////////////////////////////////////
void AliDevice::SetStatus(Int_t word)
{
// Set a user defined status word for this device.
 fStatus=word;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliDevice::GetStatus() const
{
// Provide the user defined status word for this device.
 return fStatus;
}
///////////////////////////////////////////////////////////////////////////
void AliDevice::AddHit(AliSignal& s)
{
// Register an AliSignal object as a hit to this device.
// Note : In case this device owns the AliSignal object, the pointer to
//        this device will be stored in the special owning device
//        pointer of the AliSignal object.
//        In case this device does not own the AliSignal object, a (backward)
//        link to this device is added to the first slot of the AliSignal
//        if there was no link to this device already present.

 if (!fHits)
 {
  fHits=new TObjArray(1);
  if (fHitCopy) fHits->SetOwner();
 }

 // Check if this signal is already stored for this device.
 Int_t nhits=GetNhits();
 for (Int_t i=0; i<nhits; i++)
 {
  if (&s==fHits->At(i)) return; 
 }

 // Check for existing (backward) link to this device.
 Int_t nlinks=s.GetNlinks(this);

 if (fHitCopy)
 {
  fHits->Add(s.Clone());
  // Remove unnecessary backward link(s) from the various slots
  // and set the owning link to this device
  AliSignal* sx=(AliSignal*)fHits->Last();
  if (nlinks) sx->ResetLinks(this);
  sx->SetDevice(this);
 }
 else
 {
  fHits->Add(&s);
  // Set (backward) link to the this device
  if (!nlinks) s.AddLink(this);
 }
}
///////////////////////////////////////////////////////////////////////////
void AliDevice::RemoveHit(AliSignal& s)
{
// Remove AliSignal object registered as a hit from this device.
 if (fHits)
 {
  AliSignal* test=(AliSignal*)fHits->Remove(&s);
  if (test)
  {
   fHits->Compress();
   if (fHitCopy) delete test;
  }
 }
 if (fOrdered)
 {
  AliSignal* test=(AliSignal*)fOrdered->Remove(&s);
  if (test) fOrdered->Compress();
 }
}
///////////////////////////////////////////////////////////////////////////
void AliDevice::RemoveHits()
{
// Remove all AliSignal objects registered as hits from this device.
 if (fHits)
 {
  delete fHits;
  fHits=0;
 }
 if (fOrdered)
 {
  delete fOrdered;
  fOrdered=0;
 }
 if (fMarkers)
 {
  delete fMarkers;
  fMarkers=0;
 }
}
///////////////////////////////////////////////////////////////////////////
Int_t AliDevice::GetNhits() const
{
// Provide the number of registered hits for this device.
 Int_t nhits=0;
 if (fHits) nhits=fHits->GetEntries();
 return nhits;
}
///////////////////////////////////////////////////////////////////////////
AliSignal* AliDevice::GetHit(Int_t j) const
{
// Provide the AliSignal object registered as hit number j.
// Note : j=1 denotes the first hit.
 if (!fHits) return 0;

 if ((j >= 1) && (j <= GetNhits()))
 {
  return (AliSignal*)fHits->At(j-1);
 }
 else
 {
  return 0;
 }
}
///////////////////////////////////////////////////////////////////////////
AliSignal* AliDevice::GetIdHit(Int_t id) const
{
// Return the hit with unique identifier "id".
 if (!fHits || id<0) return 0;

 AliSignal* sx=0;
 Int_t sid=0;
 for (Int_t i=0; i<GetNhits(); i++)
 {
  sx=(AliSignal*)fHits->At(i);
  if (sx)
  {
   sid=sx->GetUniqueID();
   if (id==sid) return sx;
  }
 }
 return 0; // No matching id found
}
///////////////////////////////////////////////////////////////////////////
TObjArray* AliDevice::GetHits()
{
// Provide the references to all the registered hits.
 return fHits;
}
///////////////////////////////////////////////////////////////////////////
void AliDevice::ShowHit(Int_t j) const
{
// Show data of the registered j-th hit.
// If j=0 all associated hits will be shown.
// The default is j=0.
 if (!j)
 {
  Int_t nhits=GetNhits();
  for (Int_t ih=1; ih<=nhits; ih++)
  {
   AliSignal* sx=GetHit(ih);
   if (sx) sx->Data();
  }
 }
 else
 {
  AliSignal* s=GetHit(j);
  if (s) s->Data();
 }
}
///////////////////////////////////////////////////////////////////////////
void AliDevice::Data(TString f) const
{
// Print the device and all registered hit info according to the specified
// coordinate frame.
 AliSignal::Data(f);
 Int_t nhits=GetNhits();
 if (nhits)
 {
  cout << " The following " << nhits << " hits are registered : " << endl;
  ShowHit();
 }
 else
 {
  cout << " No hits have been registered for this device." << endl;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliDevice::GetExtremes(Float_t& vmin,Float_t& vmax,Int_t idx,TObjArray* hits,Int_t mode) const
{
// Provide the min. and max. signal values of an array of hits.
// The input argument "idx" denotes the index of the signal slots to be investigated.
// The default is idx=1;
// In case hits=0 (default), the registered hits of the current device are used. 
// Signals which were declared as "Dead" will be rejected.
// The gain etc... corrected signals will be used in the process as specified
// by the  "mode" argument. The definition of this "mode" parameter corresponds to
// the description provided in the GetSignal memberfunction of class AliSignal.
// The default is mode=1 (for backward compatibility reasons).

 vmin=0;
 vmax=0;

 if (!hits) hits=fHits;
 
 if (idx<=0 || !hits) return;

 Int_t nhits=hits->GetEntries();

 Float_t sig=0;
 for (Int_t i=0; i<nhits; i++)
 {
  AliSignal* sx=(AliSignal*)hits->At(i);

  if (!sx) continue;
  if (idx > sx->GetNvalues()) continue; // User specified slotindex out of range for this signal
  if (sx->GetDeadValue(idx)) continue;  // Only take alive signals

  sig=sx->GetSignal(idx,mode);
  if (i==0)
  {
   vmin=sig;
   vmax=sig;
  }
  else
  {
   if (sig<vmin) vmin=sig;
   if (sig>vmax) vmax=sig;
  }
 }
}
///////////////////////////////////////////////////////////////////////////
void AliDevice::GetExtremes(Float_t& vmin,Float_t& vmax,TString name,TObjArray* hits,Int_t mode) const
{
// Provide the min. and max. signal values of an array of hits.
// The input argument "name" denotes the name of the signal slots to be investigated.
// In case hits=0 (default), the registered hits of the current device are used. 
// Signals which were declared as "Dead" will be rejected.
// The gain etc... corrected signals will be used in the process as specified
// by the  "mode" argument. The definition of this "mode" parameter corresponds to
// the description provided in the GetSignal memberfunction of class AliSignal.
// The default is mode=1 (for backward compatibility reasons).

 vmin=0;
 vmax=0;

 if (!hits) hits=fHits;
 
 if (!hits) return;

 Int_t nhits=hits->GetEntries();

 Int_t idx=0; // The signal slotindex to perform the sorting on

 Float_t sig=0;
 for (Int_t i=0; i<nhits; i++)
 {
  AliSignal* sx=(AliSignal*)hits->At(i);

  if (!sx) continue;

  // Obtain the slotindex corresponding to the user selection
  idx=sx->GetSlotIndex(name);
  if (!idx) continue;

  if (sx->GetDeadValue(idx)) continue; // Only take alive signals

  sig=sx->GetSignal(idx,mode);
  if (i==0)
  {
   vmin=sig;
   vmax=sig;
  }
  else
  {
   if (sig<vmin) vmin=sig;
   if (sig>vmax) vmax=sig;
  }
 }
}
///////////////////////////////////////////////////////////////////////////
TObjArray* AliDevice::SortHits(Int_t idx,Int_t mode,TObjArray* hits,Int_t mcal)
{
// Order the references to an array of hits by looping over the input array "hits"
// and checking the signal value. The ordered array is returned as a TObjArray.
// In case hits=0 (default), the registered hits of the current device are used. 
// Note that the original hit array is not modified.
// A "hit" represents an abstract object which is derived from AliSignal.
// The user can specify the index of the signal slot to perform the sorting on.
// By default the slotindex will be 1.
// Via the "mode" argument the user can specify ordering in decreasing
// order (mode=-1) or ordering in increasing order (mode=1).
// The default is mode=-1.
// Signals which were declared as "Dead" will be rejected.
// The gain etc... corrected signals will be used in the ordering process as
// specified by the "mcal" argument. The definition of this "mcal" parameter
// corresponds to the signal correction mode described in the GetSignal
// memberfunction of class AliSignal.
// The default is mcal=1 (for backward compatibility reasons).

 if (fOrdered)
 {
  delete fOrdered;
  fOrdered=0;
 }

 if (!hits) hits=fHits;
 
 if (idx<=0 || abs(mode)!=1 || !hits) return fOrdered;

 Int_t nhits=hits->GetEntries();
 if (!nhits)
 {
  return fOrdered;
 }
 else
 {
  fOrdered=new TObjArray(nhits);
 }

 Int_t nord=0;
 for (Int_t i=0; i<nhits; i++) // Loop over all hits of the array
 {
  AliSignal* s=(AliSignal*)hits->At(i);

  if (!s) continue;

  if (idx > s->GetNvalues()) continue; // User specified slotindex out of range for this signal
  if (s->GetDeadValue(idx)) continue;  // Only take alive signals
 
  if (nord == 0) // store the first hit with a signal at the first ordered position
  {
   nord++;
   fOrdered->AddAt(s,nord-1);
   continue;
  }
 
  for (Int_t j=0; j<=nord; j++) // put hit in the right ordered position
  {
   if (j == nord) // module has smallest (mode=-1) or largest (mode=1) signal seen so far
   {
    nord++;
    fOrdered->AddAt(s,j); // add hit at the end
    break; // go for next hit
   }
 
   if (mode==-1 && s->GetSignal(idx,mcal) < ((AliSignal*)fOrdered->At(j))->GetSignal(idx,mcal)) continue;
   if (mode==1 && s->GetSignal(idx,mcal) > ((AliSignal*)fOrdered->At(j))->GetSignal(idx,mcal)) continue;
 
   nord++;
   for (Int_t k=nord-1; k>j; k--) // create empty position
   {
    fOrdered->AddAt(fOrdered->At(k-1),k);
   }
   fOrdered->AddAt(s,j); // put hit at empty position
   break; // go for next matrix module
  }
 }
 return fOrdered;
}
///////////////////////////////////////////////////////////////////////////
TObjArray* AliDevice::SortHits(TString name,Int_t mode,TObjArray* hits,Int_t mcal)
{
// Order the references to an array of hits by looping over the input array "hits"
// and checking the signal value. The ordered array is returned as a TObjArray.
// In case hits=0 (default), the registered hits of the current device are used. 
// Note that the input array is not modified.
// A "hit" represents an abstract object which is derived from AliSignal.
// The user can specify the name of the signal slot to perform the sorting on.
// In case no matching slotname is found, the signal will be skipped.
// Via the "mode" argument the user can specify ordering in decreasing
// order (mode=-1) or ordering in increasing order (mode=1).
// The default is mode=-1.
// Signals which were declared as "Dead" will be rejected.
// The gain etc... corrected signals will be used in the ordering process as
// specified by the "mcal" argument. The definition of this "mcal" parameter
// corresponds to the signal correction mode described in the GetSignal
// memberfunction of class AliSignal.
// The default is mcal=1 (for backward compatibility reasons).

 if (fOrdered)
 {
  delete fOrdered;
  fOrdered=0;
 }

 if (!hits) hits=fHits;
 
 if (abs(mode)!=1 || !hits) return fOrdered;

 Int_t nhits=hits->GetEntries();
 if (!nhits)
 {
  return fOrdered;
 }
 else
 {
  fOrdered=new TObjArray(nhits);
 }

 Int_t idx=0; // The signal slotindex to perform the sorting on

 Int_t nord=0;
 for (Int_t i=0; i<nhits; i++) // loop over all hits of the array
 {
  AliSignal* s=(AliSignal*)hits->At(i);

  if (!s) continue;

  // Obtain the slotindex corresponding to the user selection
  idx=s->GetSlotIndex(name);
  if (!idx) continue;

  if (s->GetDeadValue(idx)) continue; // only take alive signals
 
  if (nord == 0) // store the first hit with a signal at the first ordered position
  {
   nord++;
   fOrdered->AddAt(s,nord-1);
   continue;
  }
 
  for (Int_t j=0; j<=nord; j++) // put hit in the right ordered position
  {
   if (j == nord) // module has smallest (mode=-1) or largest (mode=1) signal seen so far
   {
    nord++;
    fOrdered->AddAt(s,j); // add hit at the end
    break; // go for next hit
   }
 
   if (mode==-1 && s->GetSignal(idx,mcal) < ((AliSignal*)fOrdered->At(j))->GetSignal(idx,mcal)) continue;
   if (mode==1 && s->GetSignal(idx,mcal) > ((AliSignal*)fOrdered->At(j))->GetSignal(idx,mcal)) continue;
 
   nord++;
   for (Int_t k=nord-1; k>j; k--) // create empty position
   {
    fOrdered->AddAt(fOrdered->At(k-1),k);
   }
   fOrdered->AddAt(s,j); // put hit at empty position
   break; // go for next matrix module
  }
 }
 return fOrdered;
}
///////////////////////////////////////////////////////////////////////////
void AliDevice::DisplayHits(Int_t idx,Float_t scale,TObjArray* hits,Int_t dp,Int_t mode,Int_t mcol)
{
// 3D color display of an array hits.
// The user can specify the index (default=1) of the signal slot to perform the display for.
// The marker size will indicate the absolute value of the signal (specified by the slotindex)
// as a percentage of the input argument "scale".
// In case scale<0 the maximum absolute signal value encountered in the hit array will be used
// to define the 100% scale. The default is scale=-1.
// In case hits=0 (default), the registered hits of the current device are used. 
// Note that the input array is not modified.
// In case dp=1 the device position will be used, otherwise the hit position will
// be used in the display. The default is dp=0.
// Via the "mcol" argument the user can specify the marker color (see TPolyMarker3D).
// The default is mcol=blue.
// Signals which were declared as "Dead" will not be displayed.
// The gain etc... corrected signals will be used to determine the marker size.
// The gain correction is performed according to "mode" argument. The definition of this
// "mode" parameter corresponds to the description provided in the GetSignal
// memberfunction of class AliSignal.
// The default is mode=1 (for backward compatibility reasons).
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

 Int_t thisdev=0; // Indicate whether this is the owning device or not 
 if (!hits)
 {
  hits=fHits;
  thisdev=1;
 }
 
 if (idx<=0 || !hits) return;

 Int_t nhits=hits->GetEntries();
 if (!nhits) return;

 Float_t sigmax=fabs(scale);
 if (scale<0)
 {
  Float_t vmin,vmax;
  GetExtremes(vmin,vmax,idx,hits);
  sigmax=fabs(vmax);
  if (fabs(vmin)>sigmax) sigmax=fabs(vmin);
 }

 if (sigmax <=0) return;

 if (fMarkers)
 {
  delete fMarkers;
  fMarkers=0;
 }
 fMarkers=new TObjArray(nhits);
 fMarkers->SetOwner();

 Float_t pos[3];
 GetPosition(pos,"car");

 Float_t sig=0;
 for (Int_t ih=0; ih<nhits; ih++)
 {
  AliSignal* sx=(AliSignal*)hits->At(ih);
  if (!sx) continue;
  if (!dp)
  {
   sx->GetPosition(pos,"car");
  }
  else
  {
   if (!thisdev)
   {
    AliDevice* dev=sx->GetDevice();
    if (dev) dev->GetPosition(pos,"car");
   }
  }
  sig=sx->GetSignal(idx,mode);
  TPolyMarker3D* m=new TPolyMarker3D();
  m->SetMarkerStyle(8);
  m->SetMarkerColor(mcol);
  m->SetMarkerSize(100.*fabs(sig)/sigmax);
  m->SetPoint(0,pos[0],pos[1],pos[2]);
  fMarkers->Add(m);
  m->Draw();
 }
}
///////////////////////////////////////////////////////////////////////////
void AliDevice::DisplayHits(TString name,Float_t scale,TObjArray* hits,Int_t dp,Int_t mode,Int_t mcol)
{
// 3D color display of an array hits.
// The user can specify the name of the signal slot to perform the display for.
// The marker size will indicate the absolute value of the signal (specified by the slotname)
// as a percentage of the input argument "scale".
// In case scale<0 the maximum absolute signal value encountered in the hit array will be used
// to define the 100% scale. The default is scale=-1.
// In case hits=0 (default), the registered hits of the current device are used. 
// Note that the input array is not modified.
// In case dp=1 the device position will be used, otherwise the hit position will
// be used in the display. The default is dp=0.
// The marker size will indicate the percentage of the maximum encountered value
// of the absolute value of the name-specified input signal slots.
// Via the "mcol" argument the user can specify the marker color (see TPolyMarker3D).
// The default is mcol=blue.
// Signals which were declared as "Dead" will not be displayed.
// The gain etc... corrected signals will be used to determine the marker size.
// The gain correction is performed according to "mode" argument. The definition of this
// "mode" parameter corresponds to the description provided in the GetSignal
// memberfunction of class AliSignal.
// The default is mode=1 (for backward compatibility reasons).
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

 Int_t thisdev=0; // Indicate whether this is the owning device or not 
 if (!hits)
 {
  hits=fHits;
  thisdev=1;
 }
 
 if (!hits) return;

 Int_t nhits=hits->GetEntries();

 if (!nhits) return;

 Float_t sigmax=fabs(scale);
 if (scale<0)
 {
  Float_t vmin,vmax;
  GetExtremes(vmin,vmax,name,hits);
  sigmax=fabs(vmax);
  if (fabs(vmin)>sigmax) sigmax=fabs(vmin);
 }

 if (sigmax <=0) return;

 if (fMarkers)
 {
  delete fMarkers;
  fMarkers=0;
 }
 fMarkers=new TObjArray(nhits);
 fMarkers->SetOwner();

 Float_t pos[3];
 GetPosition(pos,"car");

 Int_t idx=0; // The slot index corresponding to the user specified name
 Float_t sig=0;
 for (Int_t ih=0; ih<nhits; ih++)
 {
  AliSignal* sx=(AliSignal*)hits->At(ih);
  if (!sx) continue;
  idx=sx->GetSlotIndex(name);
  if (!idx) continue;
  if (!dp)
  {
   sx->GetPosition(pos,"car");
  }
  else
  {
   if (!thisdev)
   {
    AliDevice* dev=sx->GetDevice();
    if (dev) dev->GetPosition(pos,"car");
   }
  }
  sig=sx->GetSignal(idx,mode);
  TPolyMarker3D* m=new TPolyMarker3D();
  m->SetMarkerStyle(8);
  m->SetMarkerColor(mcol);
  m->SetMarkerSize(100.*fabs(sig)/sigmax);
  m->SetPoint(0,pos[0],pos[1],pos[2]);
  fMarkers->Add(m);
  m->Draw();
 }
}
///////////////////////////////////////////////////////////////////////////
TObject* AliDevice::Clone(const char* name) const
{
// Make a deep copy of the current object and provide the pointer to the copy.
// This memberfunction enables automatic creation of new objects of the
// correct type depending on the object type, a feature which may be very useful
// for containers like AliEvent when adding objects in case the
// container owns the objects. This feature allows e.g. AliEvent
// to store either AliDevice objects or objects derived from AliDevice
// via tha AddDevice memberfunction, provided these derived classes also have
// a proper Clone memberfunction. 

 AliDevice* dev=new AliDevice(*this);
 if (name)
 {
  if (strlen(name)) dev->SetName(name);
 }
 return dev;
}
///////////////////////////////////////////////////////////////////////////
