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
// Class AliSignal
// Generic handling of (extrapolated) detector signals.
//
// Note :
// ------
// Signal positions (r) and reference frames (f) are specified via
// SetPosition(r,f) under the following conventions :
//
// f="car" ==> r is Cartesian   (x,y,z)
// f="sph" ==> r is Spherical   (r,theta,phi)
// f="cyl" ==> r is Cylindrical (rho,phi,z)
//
// The same holds for SetPositionErrors().
//
// All angles are in radians.
//
// Example :
// ---------
//
// AliSignal s;
// s.SetName("Start counter");
// Float_t pos[3]={-1,25,7};
// Float_t err[3]={0.03,0.7,0.18};
// Float_t signal=120.8;
// Float_t error=1.73;
// Float_t offset=-12.78;
// Float_t gain=250;
// s.SetPosition(pos,"car");
// s.SetPositionErrors(err,"car");
// s.SetSignal(signal);
// s.SetSignalError(error);
// s.SetOffset(offset);
// s.SetGain(gain);
// Float_t loc[3],dr[3],sigma;
// s.GetPosition(loc,"sph");
// s.GetPositionErrors(dr,"sph");
// Float_t adc=s.GetSignal();
// Float_t sigma=s.GetSignalError();
//
// AliSignal q;    // In the example below a signal contains the
//                 // following data : timing, ADC and dE/dx
// q.SetName("TOF hit");
// q.SetPosition(pos,"car");
// q.SetPositionErrors(err,"car");
// signal=82.5; // e.g. signal time in ns
// error=2.01;
// offset=0.003;
// q.SetSignal(signal,1);
// q.SetSignalError(error,1);
// q.SetOffset(offset,1);
// signal=268.1; // e.g. ADC value of signal
// error=3.75;
// gain=120.78;
// q.SetSignal(signal,2);
// q.SetSignalError(error,2);
// q.SetGain(gain,2);
// signal=23.7; // e.g. corresponding dE/dx value
// error=0.48;
// offset=0.2;
// gain=150;
// q.SetSignal(signal,3);
// q.SetSignalError(error,3);
// q.SetOffset(offset,3);
// q.SetGain(gain,3);
//
//--- Author: Nick van Eijndhoven 23-jan-1999 UU-SAP Utrecht
//- Modified: NvE $Date$ UU-SAP Utrecht
///////////////////////////////////////////////////////////////////////////

#include "AliSignal.h"
#include "AliTrack.h"
#include "Riostream.h"
 
ClassImp(AliSignal) // Class implementation to enable ROOT I/O
 
AliSignal::AliSignal() : TNamed(),AliPosition(),AliAttrib()
{
// Creation of an AliSignal object and initialisation of parameters.
// Several signal values (with errors) can be stored in different slots.
// If needed, the storage for values (and errors) will be expanded automatically
// when entering values and/or errors.
 fSignals=0;
 fDsignals=0;
 fWaveforms=0;
 fLinks=0;
}
///////////////////////////////////////////////////////////////////////////
AliSignal::~AliSignal()
{
// Destructor to delete dynamically allocated memory
 if (fSignals)
 {
  delete fSignals;
  fSignals=0;
 }
 if (fDsignals)
 {
  delete fDsignals;
  fDsignals=0;
 }
 if (fWaveforms)
 {
  delete fWaveforms;
  fWaveforms=0;
 }
 if (fLinks)
 {
  delete fLinks;
  fLinks=0;
 }
}
///////////////////////////////////////////////////////////////////////////
AliSignal::AliSignal(AliSignal& s) : TNamed(s),AliPosition(s),AliAttrib(s)
{
// Copy constructor
 fSignals=0;
 fDsignals=0;
 fWaveforms=0;
 fLinks=0;

 Int_t n=s.GetNvalues();
 Double_t val;
 for (Int_t i=1; i<=n; i++)
 {
  val=s.GetSignal(i);
  SetSignal(val,i);
 } 

 n=s.GetNerrors();
 for (Int_t j=1; j<=n; j++)
 {
  val=s.GetSignalError(j);
  SetSignalError(val,j);
 }

 n=s.GetNwaveforms();
 for (Int_t k=1; k<=n; k++)
 {
  TH1F* hist=s.GetWaveform(k);
  if (hist) SetWaveform(hist,k); 
 }

 n=s.GetNlinks();
 for (Int_t il=1; il<=n; il++)
 {
  TObject* obj=s.GetLink(il);
  if (obj) SetLink(obj,il); 
 }
}
///////////////////////////////////////////////////////////////////////////
void AliSignal::Reset(Int_t mode)
{
// Reset all signal and position values and errors to 0.
//
// mode = 0 Reset position and all signal values and their errors to 0.
//          The waveform histograms are reset, but the calibration
//          constants (i.e. gains and offsets) are kept.
//        1 Reset position and delete the signal and error storage arrays.
//          Also the waveform histograms, gains and offset arrays are deleted.
//
// The default when invoking Reset() corresponds to mode=0.
//
// Note : In all cases the storage of the various links will be cleared
//        and the container itself will be deleted to recover the memory.
//
// The usage of mode=0 allows to re-use the allocated memory for new
// signal (and error) values. This behaviour is preferable (i.e. faster)
// in case the various signals always contain the same number of values
// and have the same calibration constants.
// The usage of mode=1 is slower, but allows a more efficient memory
// occupation (and smaller output file size) in case the different
// signals have a variable number of values.
//
// For more specific actions see ResetPosition(), ResetSignals(),
// DeleteSignals(), ResetGain(), ResetOffset() and DeleteCalibrations().
//

 if (mode<0 || mode>1)
 {
  cout << " *AliSignal::Reset* Invalid argument mode = " << mode << endl;
  cout << " Default mode=0 will be used." << endl;
  mode=0;
 }

 ResetPosition();
 if (!mode)
 {
  ResetSignals();
 }
 else
 {
  DeleteSignals();
  DeleteCalibrations();
 }

 if (fLinks)
 {
  delete fLinks;
  fLinks=0;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliSignal::ResetSignals(Int_t mode)
{
// Reset various signal data according to user selection.
//
// mode = 0 Reset all signal values and their errors to 0.
//        1 Reset only signal values
//        2 Reset only signal errors
//
// The default when invoking ResetSignals() corresponds to mode=0.
//
// Irrespective of the mode, the waveform histograms are reset.

 if (mode<0 || mode>2)
 {
  cout << " *AliSignal::ResetSignals* Invalid argument mode = " << mode << endl;
  cout << " Default mode=0 will be used." << endl;
  mode=0;
 }

 if (fSignals && (mode==0 || mode==1))
 {
  for (Int_t i=0; i<fSignals->GetSize(); i++)
  {
   fSignals->AddAt(0,i);
  }
 }

 if (fDsignals && (mode==0 || mode==2))
 {
  for (Int_t j=0; j<fDsignals->GetSize(); j++)
  {
   fDsignals->AddAt(0,j);
  }
 }

 ResetWaveform(0);
}
///////////////////////////////////////////////////////////////////////////
void AliSignal::DeleteSignals(Int_t mode)
{
// Delete storage arrays of various signal data according to user selection.
//
// mode = 0 Delete arrays of both signal values and their errors.
//        1 Delete only signal values array
//        2 Delete only signal errors array
//
// The default when invoking DeleteSignals() corresponds to mode=0.
//
// Irrespective of the mode, the waveform histograms are deleted.

 if (mode<0 || mode>2)
 {
  cout << " *AliSignal::DeleteSignals* Invalid argument mode = " << mode << endl;
  cout << " Default mode=0 will be used." << endl;
  mode=0;
 }

 if (fSignals && (mode==0 || mode==1))
 {
  delete fSignals;
  fSignals=0;
 }

 if (fDsignals && (mode==0 || mode==2))
 {
  delete fDsignals;
  fDsignals=0;
 }

 DeleteWaveform(0);
}
///////////////////////////////////////////////////////////////////////////
void AliSignal::SetSignal(Double_t sig,Int_t j)
{
// Store value in the j-th (default j=1) signal slot.
// Note : The first signal slot is at j=1.
// In case the value of the index j exceeds the maximum number of reserved
// slots for signal values, the number of reserved slots for the
// signal values is increased automatically.

 if (!fSignals)
 {
  fSignals=new TArrayF(j);
  ResetSignals(1);
 }

 Int_t size=fSignals->GetSize();

 if (j>size)
 {
  fSignals->Set(j);
 }

 fSignals->AddAt(float(sig),j-1);
}
///////////////////////////////////////////////////////////////////////////
void AliSignal::AddSignal(Double_t sig,Int_t j)
{
// Add value to the j-th (default j=1) signal slot.
// Note : The first signal slot is at j=1.
// In case the value of the index j exceeds the maximum number of reserved
// slots for signal values, the number of reserved slots for the
// signal values is increased automatically.

 if (!fSignals)
 {
  fSignals=new TArrayF(j);
  ResetSignals(1);
 }

 Int_t size=fSignals->GetSize();

 if (j>size)
 {
  fSignals->Set(j);
 }

 Float_t sum=(fSignals->At(j-1))+sig;
 fSignals->AddAt(sum,j-1);
}
///////////////////////////////////////////////////////////////////////////
Float_t AliSignal::GetSignal(Int_t j,Int_t mode)
{
// Provide value of the j-th (default j=1) signal slot.
// Note : The first signal slot is at j=1.
// In case no signal is present or the argument j is invalid, 0 is returned.
// The parameter "mode" allows for automatic gain etc... correction of the signal.
//
// mode = 0 : Just the j-th signal is returned.
//        1 : The j-th signal is corrected for the gain, offset, dead flag etc...
//            In case the gain value was not set, gain=1 will be assumed.
//            In case the gain value was 0, a signal value of 0 is returned.
//            In case the offset value was not set, offset=0 will be assumed.
//            In case the j-th slot was marked dead, 0 is returned.
//
// The corrected signal (sigc) is determined as follows :
//
//              sigc=(signal/gain)-offset 
//
// The default is mode=0.

 Float_t sig=0;
 Float_t gain=1;
 Float_t offset=0;
 if (fSignals)
 {
  if (j>0 && j<=(fSignals->GetSize()))
  {
   sig=fSignals->At(j-1);

   if (mode==0) return sig;

   // Correct the signal for the gain, offset, dead flag etc...
   if (GetDeadValue(j)) return 0;

   if (GetGainFlag(j)) gain=GetGain(j);
   if (GetOffsetFlag(j)) offset=GetOffset(j);

   if (fabs(gain)>0.)
   {
    sig=(sig/gain)-offset;
   }
   else
   {
    sig=0;
   }
  }
  else
  {
   cout << " *AliSignal::GetSignal* Index j = " << j << " invalid." << endl;
  } 
 }
 return sig;
}
///////////////////////////////////////////////////////////////////////////
void AliSignal::SetSignalError(Double_t dsig,Int_t j)
{
// Store error for the j-th (default j=1) signal slot.
// Note : The first signal slot is at j=1.
// In case the value of the index j exceeds the maximum number of reserved
// slots for signal error values, the number of reserved slots for the
// signal errors is increased automatically.

 if (!fDsignals)
 {
  fDsignals=new TArrayF(j);
  ResetSignals(2);
 }

 Int_t size=fDsignals->GetSize();

 if (j>size)
 {
  fDsignals->Set(j);
 }

 fDsignals->AddAt(float(dsig),j-1);
}
///////////////////////////////////////////////////////////////////////////
Float_t AliSignal::GetSignalError(Int_t j)
{
// Provide error of the j-th (default j=1) signal slot.
// Note : The first signal slot is at j=1.
// In case no signal is present or the argument j is invalid, 0 is returned.
 Float_t err=0;
 if (fDsignals)
 {
  if (j>0 && j<=(fDsignals->GetSize()))
  {
   err=fDsignals->At(j-1);
  }
  else
  {
   cout << " *AliSignal::GetSignalError* Index j = " << j << " invalid." << endl;
  } 
 }
 return err;
}
///////////////////////////////////////////////////////////////////////////
void AliSignal::Data(TString f)
{
// Provide all signal information within the coordinate frame f.

 const char* name=GetName();
 const char* title=GetTitle();

 cout << " *" << ClassName() << "::Data*";
 if (strlen(name))  cout << " Name : " << name;
 if (strlen(title)) cout << " Title : " << title;
 cout << endl;
 cout << " Position";
 Ali3Vector::Data(f);

 List(-1);
} 
///////////////////////////////////////////////////////////////////////////
void AliSignal::List(Int_t j)
{
// Provide signal information for the j-th slot.
// The first slot is at j=1.
// In case j=0 (default) the data of all slots will be listed.
// In case j=-1 the data of all slots will be listed, but the header
// information will be suppressed.

 if (j<-1) 
 {
  cout << " *AliSignal::List* Invalid argument j = " << j << endl;
  return;
 }

 if (j != -1)
 {
  const char* name=GetName();
  const char* title=GetTitle();

  cout << " *" << ClassName() << "::Data*";
  if (strlen(name))  cout << " Name : " << name;
  if (strlen(title)) cout << " Title : " << title;
  cout << endl;
 }

 Int_t nvalues=GetNvalues();
 Int_t nerrors=GetNerrors();
 Int_t nwforms=GetNwaveforms();
 Int_t nlinks=GetNlinks();
 Int_t n=nvalues;
 if (nerrors>n) n=nerrors;
 if (nwforms>n) n=nwforms;
 if (nlinks>n) n=nlinks;

 TObject* obj=0;

 if (j<=0)
 {
  for (Int_t i=1; i<=n; i++)
  {
   cout << "   Signal";
   if (i<=nvalues) cout << " value : " << GetSignal(i);
   if (i<=nerrors) cout << " error : " << GetSignalError(i);
   AliAttrib::List(i);
   cout << endl;
   obj=GetWaveform(i);
   if (obj)
   {
    const char* wfname=obj->GetName();
    const char* wftitle=obj->GetTitle();
    cout << "    Waveform : " << obj->ClassName();
    if (strlen(wfname))  cout << " Name : " << wfname;
    if (strlen(wftitle)) cout << " Title : " << wftitle;
    cout << endl;
   }
   obj=GetLink(i);
   if (obj)
   {
    cout << "    Link to : " << obj->ClassName();
    if (obj->InheritsFrom("TNamed"))
    {
     const char* lname=obj->GetName();
     const char* ltitle=obj->GetTitle();
     if (strlen(lname))  cout << " Name : " << lname;
     if (strlen(ltitle)) cout << " Title : " << ltitle;
    }
    cout << endl;
   }
  }
 }
 else
 {
  if (j<=n)
  {
   cout << "   Signal";
   if (j<=nvalues) cout << " value : " << GetSignal(j);
   if (j<=nerrors) cout << " error : " << GetSignalError(j);
   AliAttrib::List(j);
   cout << endl;
   obj=GetWaveform(j);
   if (obj)
   {
    const char* wfnamej=obj->GetName();
    const char* wftitlej=obj->GetTitle();
    cout << "    Waveform : " << obj->ClassName();
    if (strlen(wfnamej))  cout << " Name : " << wfnamej;
    if (strlen(wftitlej)) cout << " Title : " << wftitlej;
    cout << endl;
   }
   obj=GetLink(j);
   if (obj)
   {
    cout << "    Link to : " << obj->ClassName();
    if (obj->InheritsFrom("TNamed"))
    {
     const char* lnamej=obj->GetName();
     const char* ltitlej=obj->GetTitle();
     if (strlen(lnamej))  cout << " Name : " << lnamej;
     if (strlen(ltitlej)) cout << " Title : " << ltitlej;
    }
    cout << endl;
   }
  }
 }
} 
///////////////////////////////////////////////////////////////////////////
Int_t AliSignal::GetNvalues()
{
// Provide the number of values for this signal.
 Int_t n=0;
 if (fSignals) n=fSignals->GetSize();
 return n;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliSignal::GetNerrors()
{
// Provide the number specified errors on the values for this signal.
 Int_t n=0;
 if (fDsignals) n=fDsignals->GetSize();
 return n;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliSignal::GetNwaveforms()
{
// Provide the number specified waveforms for this signal.
 Int_t n=0;
 if (fWaveforms) n=fWaveforms->GetSize();
 return n;
}
///////////////////////////////////////////////////////////////////////////
TH1F* AliSignal::GetWaveform(Int_t j)
{
// Provide pointer to the j-th waveform histogram.
 TH1F* waveform=0;
 if (j <= GetNwaveforms()) waveform=(TH1F*)fWaveforms->At(j-1);
 return waveform;
}
///////////////////////////////////////////////////////////////////////////
void AliSignal::SetWaveform(TH1F* waveform,Int_t j)
{
// Set the 1D waveform histogram corresponding to the j-th signal value.
//
// Notes :
//  The waveform of the first signal value is at j=1.
//  j=1 is the default value.
//
// In case the value of the index j exceeds the maximum number of reserved
// slots for the waveforms, the number of reserved slots for the waveforms
// is increased automatically.
//
// In case the histo pointer argument has the same value as the current waveform
// histogram pointer value, no action is taken since the user has already
// modified the actual histogram.
//
// In case the histo pointer argument is zero, the current waveform histogram
// is deleted and the pointer set to zero.
//
// In all other cases the current waveform histogram is deleted and a new
// copy of the input histogram is created which becomes the current waveform
// histogram.

 if (!fWaveforms)
 {
  fWaveforms=new TObjArray(j);
  fWaveforms->SetOwner();
 }

 if (j > fWaveforms->GetSize()) fWaveforms->Expand(j);

 TH1F* hcur=(TH1F*)fWaveforms->At(j-1);
 if (waveform != hcur)
 {
  if (hcur)
  {
   fWaveforms->Remove(hcur);
   delete hcur;
   hcur=0;
  }
  if (waveform)
  {
   hcur=new TH1F(*waveform);
   fWaveforms->AddAt(hcur,j-1);
  }
 } 
}
///////////////////////////////////////////////////////////////////////////
void AliSignal::ResetWaveform(Int_t j)
{
// Reset the waveform of the j-th (default j=1) signal value.
// This memberfunction invokes TH1F::Reset() for the corresponding waveform(s).
// To actually delete the histograms from memory, use DeleteWaveform().
// Notes : The first signal value is at j=1.
//         j=0 ==> All waveforms will be reset.
 
 if (!fWaveforms) return;

 Int_t size=fWaveforms->GetSize();

 if ((j>=0) && (j<=size))
 {
  if (j)
  {
   TH1F* hwave=(TH1F*)fWaveforms->At(j-1);
   if (hwave) hwave->Reset();
  }
  else
  {
   for (Int_t i=0; i<size; i++)
   {
    TH1F* hwave=(TH1F*)fWaveforms->At(i);
    if (hwave) hwave->Reset();
   }
  }
 }
 else
 {
  cout << " *AliSignal::ResetWaveform* Index j = " << j << " invalid." << endl;
  return;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliSignal::DeleteWaveform(Int_t j)
{
// Delete the waveform of the j-th (default j=1) signal value.
// Notes : The first signal value is at j=1.
//         j=0 ==> All waveforms will be deleted.
 
 if (!fWaveforms) return;

 Int_t size=fWaveforms->GetSize();

 if ((j>=0) && (j<=size))
 {
  if (j)
  {
   TH1F* hwave=(TH1F*)fWaveforms->At(j-1);
   if (hwave)
   {
    fWaveforms->Remove(hwave);
    delete hwave;
   }
  }
  else
  {
   delete fWaveforms;
   fWaveforms=0;
  }
 }
 else
 {
  cout << " *AliSignal::DeleteWaveform* Index j = " << j << " invalid." << endl;
  return;
 }
}
///////////////////////////////////////////////////////////////////////////
Int_t AliSignal::GetNlinks()
{
// Provide the highest slot number with a specified link for this signal.
 Int_t n=0;
 if (fLinks) n=fLinks->GetSize();
 return n;
}
///////////////////////////////////////////////////////////////////////////
TObject* AliSignal::GetLink(Int_t j)
{
// Provide pointer of the object linked to the j-th slot.
 TObject* obj=0;
 if (j <= GetNlinks()) obj=fLinks->At(j-1);
 return obj;
}
///////////////////////////////////////////////////////////////////////////
void AliSignal::SetLink(TObject* obj,Int_t j)
{
// Introduce a link (=pointer) to an object for the j-th slot.
// Only the pointer values are stored for (backward) reference, meaning
// that the objects of which the pointers are stored are NOT owned
// by the AliSignal object.
//
// Notes :
//  The link of the first slot is at j=1.
//  j=1 is the default value.
//
// In case the value of the index j exceeds the maximum number of reserved
// slots for the links, the number of reserved slots for the links
// is increased automatically (unless the pointer argument is zero).
//
// In case the pointer argument is zero, indeed a value of zero will be
// stored for the specified slot (unless j exceeds the current maximum).
//
// In principle any object derived from TObject can be referred to by this
// mechanism.
// However, this "linking back" facility was introduced to enable AliSignal slots
// to refer directly to the various AliTracks to which the AliSignal object itself
// is related (see AliTrack::AddSignal).
// Therefore, in case the input argument "obj" points to an AliTrack (or derived)
// object, the current signal is automatically related to this AliTrack
// (or derived) object.
// 
// Please also have a look at the docs of the memberfunction ResetLink()
// to prevent the situation of stored pointers to non-existent object. 

 if (!fLinks && obj)
 {
  fLinks=new TObjArray(j);
 }

 if (j>fLinks->GetSize() && obj) fLinks->Expand(j);

 if (j<=fLinks->GetSize())
 {
  fLinks->AddAt(obj,j-1);
  if (obj) 
  {
   if (obj->InheritsFrom("AliTrack"))
   {
    AliTrack* t=(AliTrack*)obj;
    t->AddSignal(*this);
   }
  }
 }
}
///////////////////////////////////////////////////////////////////////////
void AliSignal::ResetLink(Int_t j)
{
// Reset the link of the j-th (default j=1) slot.
// Notes : The first link position is at j=1.
//         j=0 ==> All links will be reset and the storage array deleted.
//
// In general the user should take care of properly clearing the corresponding
// pointer here when the referred object is deleted.
// However, this "linking back" facility was introduced to enable AliSignal slots
// to refer directly to the various AliTracks to which the AliSignal object itself
// is related (see AliTrack::AddSignal).
// As such, the AliTrack destructor already takes care of clearing the corresponding
// links from the various AliSignal slots for all the AliSignal objects that were
// related to that AliTrack. 
// So, in case the link introduced via SetLink() is the pointer of an AliTrack object,
// the user doesn't have to worry about clearing the corresponding AliTrack link from
// the AliSignal object when the corresponding AliTrack object is deleted.
 
 if (!fLinks || j<0) return;

 TObject* obj=0;

 if (j)
 {
  SetLink(obj,j);
 }
 else
 {
  delete fLinks;
  fLinks=0;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliSignal::ResetLink(TObject* obj)
{
// Reset the link of the all the slots referring to the specified object.
//
// In general the user should take care of properly clearing the corresponding
// pointer here when the referred object is deleted.
// However, this "linking back" facility was introduced to enable AliSignal slots
// to refer directly to the various AliTracks to which the AliSignal object itself
// is related (see AliTrack::AddSignal).
// As such, the AliTrack destructor already takes care of clearing the corresponding
// links from the various AliSignal slots for all the AliSignal objects that were
// related to that AliTrack. 
// So, in case the link introduced via SetLink() is the pointer of an AliTrack object,
// the user doesn't have to worry about clearing the corresponding AliTrack link from
// the AliSignal object when the corresponding AliTrack object is deleted.
 
 if (!fLinks || !obj) return;

 Int_t nlinks=GetNlinks();
 for (Int_t i=1; i<=nlinks; i++)
 {
  TObject* obj2=GetLink(i);
  if (obj2==obj) ResetLink(i);
 }
}
///////////////////////////////////////////////////////////////////////////
TObject* AliSignal::Clone(const char* name)
{
// Make a deep copy of the current object and provide the pointer to the copy.
// This memberfunction enables automatic creation of new objects of the
// correct type depending on the object type, a feature which may be very useful
// for containers when adding objects in case the container owns the objects.
// This feature allows e.g. AliTrack to store either AliSignal objects or
// objects derived from AliSignal via the AddSignal memberfunction, provided
// these derived classes also have a proper Clone memberfunction. 

 AliSignal* sig=new AliSignal(*this);
 if (name)
 {
  if (strlen(name)) sig->SetName(name);
 }
 return sig;
}
///////////////////////////////////////////////////////////////////////////
