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
// The user can decide to store either calibrated or uncalibrated signals.
// Via the specification of a gain and offset or/and an explicit
// (de)calibration function both calibrated and uncalibrated signals
// can always be obtained. For details see the documentation of the
// memberfunction GetSignal() and the class AliAttrib.
// The explicit specification of a (de)calibration function offers the
// maximum flexibility and also allows automatic indication whether
// calibrated or uncalibrated data has been stored.
// The latter can be achieved by only specifying a calibration function
// (and no de-calibration function) in case uncalibrated data is stored,
// whereas in case of stored calibrated data the user should only
// provide a de-calibration function (and no calibration function).
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
// q.SetNameTitle("Hybrid","Test for multiple signal data");
// q.SetPosition(pos,"car");
// q.SetPositionErrors(err,"car");
// signal=82.5; // e.g. signal time in ns
// error=2.01;
// offset=0.003;
// q.SetSlotName("TOF",1);
// q.SetSignal(signal,1);
// q.SetSignalError(error,1);
// q.SetOffset(offset,1);
// signal=268.1; // e.g. ADC value of signal
// error=3.75;
// gain=120.78;
// offset=1.5732;
// // Addressing via name specification instead of index 
// q.SetSlotName("ADC",2);
// q.SetSignal(signal,"ADC");
// q.SetSignalError(error,"ADC");
// q.SetGain(gain,"ADC");
// q.SetOffset(offset,"ADC");
// signal=23.7; // e.g. corresponding dE/dx value
// error=0.48;
// TF1 f=("calib","[0]*pow(x,2)+[1]"); // dE/dx calib. function
// f.SetParameter(0,3.285);
// f.SetParameter(1,-18.67);
// q.SetSlotName("dE/dx",3);
// q.SetCalFunction(&f,"dE/dx");
// q.SetSignal(signal,"dE/dx");
// q.SetSignalError(error,"dE/dx");
//
// // Signal retrieval with various (de)calibration modes
// Float_t tof=q.GetSignal("TOF");
// Float_t adc=q.GetSignal("ADC",1);
// Float_t dedx=q.GetSignal("dE/dx",3);
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
 fSigflags=0;
 fWaveforms=0;
 fLinks=0;
 fDevice=0;
 fTracks=0;
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
 if (fSigflags)
 {
  delete fSigflags;
  fSigflags=0;
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
 if (fTracks)
 {
  // Remove this signal from all related tracks
  for (Int_t i=1; i<=GetNtracks(); i++)
  {
   AliTrack* tx=GetTrack(i);
   if (tx) tx->RemoveSignal(*this,0);
  }
  delete fTracks;
  fTracks=0;
 }
}
///////////////////////////////////////////////////////////////////////////
AliSignal::AliSignal(const AliSignal& s) : TNamed(s),AliPosition(s),AliAttrib(s)
{
// Copy constructor
 fSignals=0;
 fDsignals=0;
 fSigflags=0;
 fWaveforms=0;
 fLinks=0;
 fTracks=0;

 // Don't copy the owning device pointer for the copy
 fDevice=0;

 Int_t n=s.GetNvalues();
 Double_t val;
 Int_t lock;
 for (Int_t i=1; i<=n; i++)
 {
  if (s.GetSignalFlag(i))
  {
   val=s.GetSignal(i);
   lock=s.GetLockValue();
   if (lock) Unlock();
   SetSignal(val,i);
   if (lock) Lock();
  }
 } 

 n=s.GetNerrors();
 for (Int_t j=1; j<=n; j++)
 {
  if (s.GetErrorFlag(j))
  {
   val=s.GetSignalError(j);
   lock=s.GetLockValue();
   if (lock) Unlock();
   SetSignalError(val,j);
   if (lock) Lock();
  }
 }

 n=s.GetNwaveforms();
 for (Int_t k=1; k<=n; k++)
 {
  TH1F* hist=s.GetWaveform(k);
  if (hist) SetWaveform(hist,k); 
 }

 TArrayI slotarr;
 TArrayI posarr;
 TObject* dum=0;
 n=s.GetIndices(dum,slotarr,posarr);
 Int_t slot,pos;
 for (Int_t idx=0; idx<n; idx++)
 {
  slot=slotarr.At(idx);
  pos=posarr.At(idx);
  TObject* obj=s.GetLink(slot,pos);
  if (obj) SetLink(obj,slot,pos); 
 }

 Int_t ntk=s.GetNtracks();
 if (ntk)
 {
  fTracks=new TObjArray(ntk);
  for (Int_t it=1; it<=ntk; it++)
  {
   AliTrack* tx=s.GetTrack(it);
   fTracks->Add(tx);
  }
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
// Note : In all cases the storage of the various links will be reset.
//        The UniqueID, name and title will NOT be reset.
//        In case the user wants to reset these attributes, this has to
//        be done explicitly via the SET facilities. 
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
// DeleteSignals(), ResetGain(), ResetOffset(), ResetLink(), ResetWaveform(),
// DeleteWaveform() and DeleteCalibrations().
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

 if (fLinks) fLinks->Reset();
 fDevice=0;

 if (fTracks)
 {
  delete fTracks;
  fTracks=0;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliSignal::ResetSignals(Int_t mode)
{
// Reset various signal data according to user selection.
//
// mode = 0 Reset all signal values, their errors and all waveform histos.
//        1 Reset only signal values and waveform histos.
//        2 Reset only signal errors and waveform histos.
//       -1 Reset only signal values.
//       -2 Reset only signal errors.
//
// The default when invoking ResetSignals() corresponds to mode=0.

 if (abs(mode)>2)
 {
  cout << " *AliSignal::ResetSignals* Invalid argument mode = " << mode << endl;
  cout << " Default mode=0 will be used." << endl;
  mode=0;
 }

 Int_t sflag=0;
 Int_t eflag=0;

 if (fSignals && (abs(mode)==0 || abs(mode)==1))
 {
  for (Int_t i=1; i<=fSignals->GetSize(); i++)
  {
   fSignals->AddAt(0,i-1);
   eflag=GetErrorFlag(i);
   SetSigFlags(0,eflag,i);
  }
 }

 if (fDsignals && (abs(mode)==0 || abs(mode)==2))
 {
  for (Int_t j=1; j<=fDsignals->GetSize(); j++)
  {
   fDsignals->AddAt(0,j-1);
   sflag=GetSignalFlag(j);
   SetSigFlags(sflag,0,j);
  }
 }

 if (mode>=0) ResetWaveform(0);
}
///////////////////////////////////////////////////////////////////////////
void AliSignal::DeleteSignals(Int_t mode)
{
// Delete storage arrays of various signal data according to user selection.
//
// mode = 0 Delete arrays of signal values, their errors and all waveform histos.
//        1 Delete only signal values array and waveform histos.
//        2 Delete only signal errors array and waveform histos.
//       -1 Delete only signal values array.
//       -2 Delete only signal errors array.
//
// The default when invoking DeleteSignals() corresponds to mode=0.

 if (abs(mode)>2)
 {
  cout << " *AliSignal::DeleteSignals* Invalid argument mode = " << mode << endl;
  cout << " Default mode=0 will be used." << endl;
  mode=0;
 }

 if (fSignals && (abs(mode)==0 || abs(mode)==1))
 {
  delete fSignals;
  fSignals=0;
 }

 if (fDsignals && (abs(mode)==0 || abs(mode)==2))
 {
  delete fDsignals;
  fDsignals=0;
 }

 Int_t sflag=0;
 Int_t eflag=0;

 if (mode==0)
 {
  delete fSigflags;
  fSigflags=0;
 }
 else if (abs(mode)==1)
 {
  for (Int_t i=1; i<=fSigflags->GetSize(); i++)
  {
   eflag=GetErrorFlag(i);
   SetSigFlags(0,eflag,i);
  }
 }
 else if (abs(mode)==2)
 {
  for (Int_t j=1; j<=fSigflags->GetSize(); j++)
  {
   sflag=GetSignalFlag(j);
   SetSigFlags(sflag,0,j);
  }
 }

 if (mode>=0) DeleteWaveform(0);
}
///////////////////////////////////////////////////////////////////////////
void AliSignal::SetSignal(Double_t sig,Int_t j)
{
// Store signal value for the j-th (default j=1) slot.
// Notes :
// -------
// 1) The first signal slot is at j=1.
// 2) In case the 'lock' flag was set for the specified slot, the new
//    signal value will not be stored.
//    One has to unlock the specified slot first in case one really wants
//    to overwite the signal value.
// In case the value of the index j exceeds the maximum number of reserved
// slots for signal values, the number of reserved slots for the
// signal values is increased automatically.

 if (GetLockValue(j)) return;

 if (!fSignals)
 {
  fSignals=new TArrayF(j);
  ResetSignals(-1);
 }

 Int_t size=fSignals->GetSize();

 if (j>size)
 {
  fSignals->Set(j);
 }

 fSignals->AddAt(float(sig),j-1);

 Int_t eflag=GetErrorFlag(j);
 SetSigFlags(1,eflag,j);
}
///////////////////////////////////////////////////////////////////////////
void AliSignal::SetSignal(Double_t sig,TString name)
{
// Store signal value for the name-specified slot.
//
// Note :
// ------
// In case the 'lock' flag was set for the specified slot, the new
// signal value will not be stored.
// One has to unlock the specified slot first in case one really wants
// to overwite the signal value.
//
// This procedure involves a slot-index search based on the specified name
// at each invokation. This may become slow in case many slots have been
// defined and/or when this procedure is invoked many times.
// In such cases it is preferable to use indexed addressing in the user code
// either directly or via a few invokations of GetSlotIndex().

 Int_t j=GetSlotIndex(name);
 if (j>0)
 {
  if (!GetLockValue(j)) SetSignal(sig,j);
 }
}
///////////////////////////////////////////////////////////////////////////
void AliSignal::AddSignal(Double_t sig,Int_t j)
{
// Add value to the signal of the j-th (default j=1) slot.
// Notes :
// -------
// 1) The first signal slot is at j=1.
// 2) In case the 'lock' flag was set for the specified slot, the new
//    signal value will not be stored.
//    One has to unlock the specified slot first in case one really wants
//    to overwite the signal value.
// In case the value of the index j exceeds the maximum number of reserved
// slots for signal values, the number of reserved slots for the
// signal values is increased automatically.

 if (GetLockValue(j)) return;

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

 Int_t eflag=GetErrorFlag(j);
 SetSigFlags(1,eflag,j);
}
///////////////////////////////////////////////////////////////////////////
void AliSignal::AddSignal(Double_t sig,TString name)
{
// Add value to the signal of the name-specified slot.
//
// Note :
// ------
// In case the 'lock' flag was set for the specified slot, the new
// signal value will not be stored.
// One has to unlock the specified slot first in case one really wants
// to overwite the signal value.
//
// This procedure involves a slot-index search based on the specified name
// at each invokation. This may become slow in case many slots have been
// defined and/or when this procedure is invoked many times.
// In such cases it is preferable to use indexed addressing in the user code
// either directly or via a few invokations of GetSlotIndex().

 Int_t j=GetSlotIndex(name);
 if (j>0)
 {
  if (!GetLockValue(j)) AddSignal(sig,j);
 }
}
///////////////////////////////////////////////////////////////////////////
Float_t AliSignal::GetSignal(Int_t j,Int_t mode) const
{
// Provide signal value of the j-th (default j=1) slot.
// Note : The first signal slot is at j=1.
// In case no signal is present or the input argument "j" or "mode" is invalid,
// the value 0 is returned.
// The parameter "mode" allows for automatic (de)calibration of the signal
// (e.g. gain etc... correction or via explicit (de)calibration functions).
//
// mode = 0 : Just the j-th signal is returned.
//        1 : The j-th signal is corrected for the gain, offset, dead flag etc...
//            In case the j-th slot was marked dead, 0 is returned.
//            In case the gain value was not set, gain=1 will be assumed.
//            In case the gain value was 0, a signal value of 0 is returned.
//            In case the offset value was not set, offset=0 will be assumed.
//        2 : Same as mode=1 but gain, offset dead flag etc... are taken from
//            the AliDevice which owns this AliSignal object.
//            The corresponding AliDevice slot is obtained via matching of
//            the slotnames. In case this fails, the slotindex "j" of the
//            input argument will be used. 
//            In case this AliSignal object has no parent AliDevice, just
//            the j-th signal is returned (like with mode=0).
//        3 : The j-th signal is corrected using the corresponding calibration
//            function.
//            In case the j-th slot was marked dead, 0 is returned.
//            In case no calibration function is present, just the j-th signal
//            is returned (like with mode=0).
//        4 : Same as mode=3 but the calibration function and dead flag are
//            taken from the AliDevice which owns this AliSignal object.
//            The corresponding AliDevice slot is obtained via matching of
//            the slotnames. In case this fails, the slotindex "j" of the
//            input argument will be used. 
//        5 : Same as mode=2 but in case no parent AliDevice is present
//            an automatic switch to mode=1 will be made.
//        6 : Same as mode=4 but in case no parent AliDevice is present
//            an automatic switch to mode=3 will be made.
//        7 : Same as mode=3 but in case no calibration function is present
//            an automatic switch to mode=4 will be made.
//        8 : Same as mode=7 but also the corresponding dead flag of the
//            parent device (if any) will be checked.
//            If either the dead flag of the requested signal slot of this
//            AliSignal object or the corresponding parent device slot is
//            set, 0 is returned.
//
//       <0 : The corresponding de-correction or de-calibration is performed
//
// The corrected signal (sigc) is determined as follows :
//
//              sigc=(signal/gain)-offset 
//
// The de-corrected signal is determined as follows :
//
//              signal=(sigc+offset)*gain 
//
// The default is mode=0.

 if (abs(mode)>8) return 0;

 Int_t jcal=j;
 Float_t sig=0;
 Float_t gain=1;
 Float_t offset=0;

 // Get the corresponding slot index (and dead flag) of the parent device
 Int_t pj=0;
 Int_t pdead=0;
 AliSignal* parent=(AliSignal*)GetDevice();
 if ((abs(mode)==2 || abs(mode)>=4) && parent)
 {
   TString name=GetSlotName(j);
   if (strlen(name.Data())) pj=parent->GetSlotIndex(name);
   if (abs(mode)==8 && pj) pdead=parent->GetDeadValue(pj);
 }
 if (mode==8) mode=7;
 if (mode==-8) mode=-7;

 AliSignal* sx=(AliSignal*)this;

 TF1* f=0;
 if (mode==7)
 {
  f=sx->GetCalFunction(jcal);
  if (f)
  {
   mode=3;
  }
  else
  {
   mode=4;
  } 
 }
 if (mode==-7)
 {
  f=sx->GetDecalFunction(jcal);
  if (f)
  {
   mode=-3;
  }
  else
  {
   mode=-4;
  } 
 }

 if (abs(mode)==2 || abs(mode)>=4)
 {
  sx=(AliSignal*)GetDevice();
  if (pj) jcal=pj;
 }
 if (!sx && abs(mode)>=5) sx=(AliSignal*)this;
 if (mode==5) mode=2;
 if (mode==-5) mode=-2;
 if (mode==6) mode=3;
 if (mode==-6) mode=-3;

 if (fSignals)
 {
  if (j>0 && j<=(fSignals->GetSize()))
  {
   sig=fSignals->At(j-1);

   if (mode==0 || !sx) return sig;

   // Check for the dead flag setting
   if (sx->GetDeadValue(jcal) || pdead) return 0;

   // (De)correct the signal for the gain and offset
   if (abs(mode)==1 || abs(mode)==2)
   {
    if (sx->GetGainFlag(jcal)) gain=sx->GetGain(jcal);
    if (sx->GetOffsetFlag(jcal)) offset=sx->GetOffset(jcal);

    if (fabs(gain)>0.)
    {
     if (mode>0) sig=(sig/gain)-offset; // Gain and offset correction
     if (mode<0) sig=(sig+offset)*gain; // Gain and offset de-correction
    }
    else
    {
     sig=0;
    }
    return sig;
   }

   // (De)calibrate the signal with the corresponding (de)calibration function
   if (abs(mode)==3 || abs(mode)==4)
   {
    f=sx->GetCalFunction(jcal);
    if (mode<0) f=sx->GetDecalFunction(jcal);
    if (f) sig=f->Eval(sig);
    return sig;
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
Float_t AliSignal::GetSignal(TString name,Int_t mode) const
{
// Provide signal value of the name-specified slot.
// In case no signal is present, 0 is returned.
// The parameter "mode" allows for automatic (de)calibration of the signal
// (e.g. gain etc... correction or via explicit (de)calibration functions).
// For further details about the (de)calibration modes, please refer to the
// corresponding slot-index based memberfunction. 
//
// The default is mode=0.
//
// This procedure involves a slot-index search based on the specified name
// at each invokation. This may become slow in case many slots have been
// defined and/or when this procedure is invoked many times.
// In such cases it is preferable to use indexed addressing in the user code
// either directly or via a few invokations of GetSlotIndex().

 Int_t j=GetSlotIndex(name);
 Float_t val=0;
 if (j>0) val=GetSignal(j,mode);
 return val;
}
///////////////////////////////////////////////////////////////////////////
void AliSignal::SetSignalError(Double_t dsig,Int_t j)
{
// Store error on the signal for the j-th (default j=1) slot.
// Notes :
// -------
// 1) The first signal slot is at j=1.
// 2) In case the 'lock' flag was set for the specified slot, the new
//    signal error value will not be stored.
//    One has to unlock the specified slot first in case one really wants
//    to overwite the signal error value.
// In case the value of the index j exceeds the maximum number of reserved
// slots for signal error values, the number of reserved slots for the
// signal errors is increased automatically.

 if (GetLockValue(j)) return;

 if (!fDsignals)
 {
  fDsignals=new TArrayF(j);
  ResetSignals(-2);
 }

 Int_t size=fDsignals->GetSize();

 if (j>size)
 {
  fDsignals->Set(j);
 }

 fDsignals->AddAt(float(dsig),j-1);

 Int_t sflag=GetSignalFlag(j);
 SetSigFlags(sflag,1,j);
}
///////////////////////////////////////////////////////////////////////////
void AliSignal::SetSignalError(Double_t dsig,TString name)
{
// Store error on the signal for the name-specified slot.
//
// Note :
// ------
// In case the 'lock' flag was set for the specified slot, the new
// signal error value will not be stored.
// One has to unlock the specified slot first in case one really wants
// to overwite the signal error value.
//
// This procedure involves a slot-index search based on the specified name
// at each invokation. This may become slow in case many slots have been
// defined and/or when this procedure is invoked many times.
// In such cases it is preferable to use indexed addressing in the user code
// either directly or via a few invokations of GetSlotIndex().

 Int_t j=GetSlotIndex(name);
 if (j>0)
 {
  if (!GetLockValue(j)) SetSignalError(dsig,j);
 }
}
///////////////////////////////////////////////////////////////////////////
Float_t AliSignal::GetSignalError(Int_t j) const
{
// Provide error on the signal of the j-th (default j=1) slot.
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
Float_t AliSignal::GetSignalError(TString name) const
{
// Provide error on the signal of the name-specified slot.
//
// This procedure involves a slot-index search based on the specified name
// at each invokation. This may become slow in case many slots have been
// defined and/or when this procedure is invoked many times.
// In such cases it is preferable to use indexed addressing in the user code
// either directly or via a few invokations of GetSlotIndex().

 Int_t j=GetSlotIndex(name);
 Float_t val=0;
 if (j>0) val=GetSignalError(j);
 return val;
}
///////////////////////////////////////////////////////////////////////////
void AliSignal::Data(TString f,TString u) const
{
// Provide all signal information within the coordinate frame f.
//
// The string argument "u" allows to choose between different angular units
// in case e.g. a spherical frame is selected.
// u = "rad" : angles provided in radians
//     "deg" : angles provided in degrees
//
// The defaults are f="car" and u="rad".

 const char* name=GetName();
 const char* title=GetTitle();

 cout << " *" << ClassName() << "::Data* Id : " << GetUniqueID();
 if (strlen(name))  cout << " Name : " << name;
 if (strlen(title)) cout << " Title : " << title;
 cout << endl;
 cout << "   Position";
 AliPosition::Data(f,u);
 if (fDevice)
 {
  const char* devname=fDevice->GetName();
  const char* devtitle=fDevice->GetTitle();
  cout << "   Owned by device : " << fDevice->ClassName()
       << " Id : " << fDevice->GetUniqueID();
  if (strlen(devname))  cout << " Name : " << devname;
  if (strlen(devtitle)) cout << " Title : " << devtitle;
  cout << endl;
 }

 // Provide an overview of the stored waveforms
 ListWaveform(-1);

 // Provide an overview of the associated tracks
 ListTrack(-1);

 // Provide an overview of all the data and attribute slots
 List(-1);
} 
///////////////////////////////////////////////////////////////////////////
void AliSignal::List(Int_t j) const
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

  cout << " *" << ClassName() << "::Data* Id :" << GetUniqueID();
  if (strlen(name))  cout << " Name : " << name;
  if (strlen(title)) cout << " Title : " << title;
  cout << endl;
  if (fDevice)
  {
   const char* devname=fDevice->GetName();
   const char* devtitle=fDevice->GetTitle();
   cout << "   Owned by device : " << fDevice->ClassName();
   if (strlen(devname))  cout << " Name : " << devname;
   if (strlen(devtitle)) cout << " Title : " << devtitle;
   cout << endl;
  }
 }

 Int_t n=GetNslots();
 Int_t nlinkslots=0;
 if (GetNlinks()) nlinkslots=fLinks->GetMaxColumn();
 if (nlinkslots>n) n=nlinkslots;
 
 TObject* obj=0;
 Int_t nrefs=0;
 TArrayI posarr;
 Int_t pos;

 if (j<=0)
 {
  for (Int_t i=1; i<=n; i++)
  {
   obj=0;
   nrefs=GetIndices(obj,i,posarr);

   if (GetSignalFlag(i) || GetErrorFlag(i) || GetCalFunction(i) || GetDecalFunction(i) || GetCalWord(i) || nrefs)
   {
    cout << "   Slot : " << i;
    if (GetSignalFlag(i)) cout << " Signal value : " << GetSignal(i);
    if (GetErrorFlag(i))  cout << " error : " << GetSignalError(i);
    AliAttrib::List(i);
    cout << endl;

    for (Int_t k=0; k<nrefs; k++)
    {
     pos=posarr.At(k);
     obj=GetLink(i,pos);
     if (obj)
     {
      cout << "    Link at position " << pos << " to : " << obj->ClassName();
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
  }
 }
 else
 {
  if (j<=n)
  {
   obj=0;
   nrefs=GetIndices(obj,j,posarr);

   if (GetSignalFlag(j) || GetErrorFlag(j) || GetCalFunction(j) || GetDecalFunction(j) || GetCalWord(j) || nrefs)
   {
    cout << "   Slot : " << j;
    if (GetSignalFlag(j)) cout << " Signal value : " << GetSignal(j);
    if (GetErrorFlag(j))  cout << " error : " << GetSignalError(j);
    AliAttrib::List(j);
    cout << endl;

    for (Int_t kj=0; kj<nrefs; kj++)
    {
     pos=posarr.At(kj);
     obj=GetLink(j,pos);
     if (obj)
     {
      cout << "    Link at position " << pos << " to : " << obj->ClassName();
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
 }
} 
///////////////////////////////////////////////////////////////////////////
void AliSignal::List(TString name) const
{
// Provide signal information for the name-specified slot.
//
// This procedure involves a slot-index search based on the specified name
// at each invokation. This may become slow in case many slots have been
// defined and/or when this procedure is invoked many times.
// In such cases it is preferable to use indexed addressing in the user code
// either directly or via a few invokations of GetSlotIndex().

 Int_t j=GetSlotIndex(name);
 if (j>0) List(j);
}
///////////////////////////////////////////////////////////////////////////
void AliSignal::ListWaveform(Int_t j) const
{
// Provide information for the j-th waveform.
// The first waveform is at j=1.
// In case j=0 (default) the info of all waveforms will be listed.
// In case j=-1 the info of all waveforms will be listed, but the header
// information will be suppressed.

 if (j<-1) 
 {
  cout << " *AliSignal::ListWaveform* Invalid argument j = " << j << endl;
  return;
 }

 if (j != -1)
 {
  const char* name=GetName();
  const char* title=GetTitle();

  cout << " *" << ClassName() << "::Data* Id :" << GetUniqueID();
  if (strlen(name))  cout << " Name : " << name;
  if (strlen(title)) cout << " Title : " << title;
  cout << endl;
  if (fDevice)
  {
   const char* devname=fDevice->GetName();
   const char* devtitle=fDevice->GetTitle();
   cout << "   Owned by device : " << fDevice->ClassName();
   if (strlen(devname))  cout << " Name : " << devname;
   if (strlen(devtitle)) cout << " Title : " << devtitle;
   cout << endl;
  }
 }

 Int_t n=GetNwaveforms();
 TObject* obj=0;

 if (j<=0)
 {
  for (Int_t i=1; i<=n; i++)
  {
   obj=GetWaveform(i);
   if (obj)
   {
    const char* wfname=obj->GetName();
    const char* wftitle=obj->GetTitle();
    cout << "    Waveform " << i << " : " << obj->ClassName();
    if (strlen(wfname))  cout << " Name : " << wfname;
    if (strlen(wftitle)) cout << " Title : " << wftitle;
    cout << endl;
   }
  }
 }
 else
 {
  if (j<=n)
  {
   obj=GetWaveform(j);
   if (obj)
   {
    const char* wfnamej=obj->GetName();
    const char* wftitlej=obj->GetTitle();
    cout << "    Waveform " << j << " : " << obj->ClassName();
    if (strlen(wfnamej))  cout << " Name : " << wfnamej;
    if (strlen(wftitlej)) cout << " Title : " << wftitlej;
    cout << endl;
   }
  }
 }
}
///////////////////////////////////////////////////////////////////////////
void AliSignal::ListTrack(Int_t j) const
{
// Provide information for the j-th associated track.
// The first associated track is at j=1.
// In case j=0 (default) the info of all associated tracks will be listed.
// In case j=-1 the info of all tracks will be listed, but the header
// information will be suppressed.

 if (j<-1) 
 {
  cout << " *AliSignal::ListTrack* Invalid argument j = " << j << endl;
  return;
 }

 if (j != -1)
 {
  const char* name=GetName();
  const char* title=GetTitle();

  cout << " *" << ClassName() << "::Data* Id :" << GetUniqueID();
  if (strlen(name))  cout << " Name : " << name;
  if (strlen(title)) cout << " Title : " << title;
  cout << endl;
  if (fDevice)
  {
   const char* devname=fDevice->GetName();
   const char* devtitle=fDevice->GetTitle();
   cout << "   Owned by device : " << fDevice->ClassName();
   if (strlen(devname))  cout << " Name : " << devname;
   if (strlen(devtitle)) cout << " Title : " << devtitle;
   cout << endl;
  }
 }

 Int_t n=GetNtracks();
 AliTrack* tx=0;

 if (j<=0)
 {
  for (Int_t i=1; i<=n; i++)
  {
   tx=GetTrack(i);
   if (tx)
   {
    const char* txname=tx->GetName();
    const char* txtitle=tx->GetTitle();
    cout << "    Track " << i << " : " << tx->ClassName() << " Id : " << tx->GetId();
    if (strlen(txname))  cout << " Name : " << txname;
    if (strlen(txtitle)) cout << " Title : " << txtitle;
    cout << endl;
   }
  }
 }
 else
 {
  if (j<=n)
  {
   tx=GetTrack(j);
   if (tx)
   {
    const char* txnamej=tx->GetName();
    const char* txtitlej=tx->GetTitle();
    cout << "    Track " << j << " : " << tx->ClassName() << " Id : " << tx->GetId();
    if (strlen(txnamej))  cout << " Name : " << txnamej;
    if (strlen(txtitlej)) cout << " Title : " << txtitlej;
    cout << endl;
   }
  }
 }
}
///////////////////////////////////////////////////////////////////////////
Int_t AliSignal::GetNvalues() const
{
// Provide the number of values for this signal.
 
 if (!fSignals) return 0;

 Int_t n=0;
 for (Int_t i=1; i<=fSigflags->GetSize(); i++)
 {
  if (GetSignalFlag(i)) n=i;
 }

 return n;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliSignal::GetNerrors() const
{
// Provide the number specified errors on the values for this signal.
 
 if (!fDsignals) return 0;

 Int_t n=0;
 for (Int_t i=1; i<=fSigflags->GetSize(); i++)
 {
  if (GetErrorFlag(i)) n=i;
 }

 return n;
}
///////////////////////////////////////////////////////////////////////////
void AliSignal::SetSigFlags(Int_t is,Int_t ie,Int_t j)
{
// Store signal and/or error value flags of the j-th (default j=1) slot.
// Note : The first slot is at j=1.
// In case the value of the index j exceeds the maximum number of reserved
// slots for the flags, the number of reserved slots for the flags is
// increased automatically.
// The value stored is : 10*signalflag + errorflag.

 if (j<1) 
 {
  cout << " *AliSignal::SetSigFlags* Invalid argument j = " << j << endl;
  return;
 }

 if (!fSigflags)
 {
  fSigflags=new TArrayI(j);
 }

 Int_t size=fSigflags->GetSize();

 if (j>size)
 {
  fSigflags->Set(j);
 }

 Int_t word=10*is+ie;
 
 fSigflags->AddAt(word,j-1);
}
///////////////////////////////////////////////////////////////////////////
Int_t AliSignal::GetSignalFlag(Int_t j) const
{
// Provide signal value flag of the j-th (default j=1) slot.
//
// flag = 1 : Signal value was set
//        0 : Signal value was not set
//
// Note : The first attribute slot is at j=1.
// In case j is invalid, 0 is returned.

 if (j<1) 
 {
  cout << " *AliSignal::GetSignalFlag* Invalid argument j = " << j << endl;
  return 0;
 }
 Int_t flag=0;
 if (fSigflags)
 {
  if (j>0 && j<=(fSigflags->GetSize()))
  {
   Int_t word=fSigflags->At(j-1);
   flag=word/10;
  }
 }
 return flag;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliSignal::GetSignalFlag(TString name) const
{
// Provide signal value flag of the name-specified slot.
//
// flag = 1 : Signal value was set
//        0 : Signal value was not set
//
//
// This procedure involves a slot-index search based on the specified name
// at each invokation. This may become slow in case many slots have been
// defined and/or when this procedure is invoked many times.
// In such cases it is preferable to use indexed addressing in the user code
// either directly or via a few invokations of GetSlotIndex().

 Int_t j=GetSlotIndex(name);
 Int_t flag=0;
 if (j>0) flag=GetSignalFlag(j);
 return flag;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliSignal::GetErrorFlag(Int_t j) const
{
// Provide error value flag of the j-th (default j=1) slot.
//
// flag = 1 : Error value was set
//        0 : Error value was not set
//
// Note : The first attribute slot is at j=1.
// In case j is invalid, 0 is returned.

 if (j<1) 
 {
  cout << " *AliSignal::GetErrorFlag* Invalid argument j = " << j << endl;
  return 0;
 }
 Int_t flag=0;
 if (fSigflags)
 {
  if (j>0 && j<=(fSigflags->GetSize()))
  {
   Int_t word=fSigflags->At(j-1);
   flag=word%10;
  }
 }
 return flag;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliSignal::GetErrorFlag(TString name) const
{
// Provide error value flag of the name-specified slot.
//
// flag = 1 : Error value was set
//        0 : Error value was not set
//
//
// This procedure involves a slot-index search based on the specified name
// at each invokation. This may become slow in case many slots have been
// defined and/or when this procedure is invoked many times.
// In such cases it is preferable to use indexed addressing in the user code
// either directly or via a few invokations of GetSlotIndex().

 Int_t j=GetSlotIndex(name);
 Int_t flag=0;
 if (j>0) flag=GetErrorFlag(j);
 return flag;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliSignal::GetNslots() const
{
// Provide the number of existing slots.

 Int_t n=AliAttrib::GetNslots();

 if (!fSigflags) return n;

 Int_t nflags=0;
 for (Int_t i=0; i<fSigflags->GetSize(); i++)
 {
  if (fSigflags->At(i)) nflags=i+1;
 }

 if (n<nflags) n=nflags;

 return n;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliSignal::GetNwaveforms() const
{
// Provide the number of specified waveforms for this signal.
// Actually the return value is the highest index of the stored waveforms.
// This allows an index dependent meaning of waveform info (e.g. waveforms
// with various gain values).
// So, when all waveforms are stored in consequetive positions (e.g. 1,2,3),
// this memberfunction returns 3, being both the highest filled position
// and the actual number of waveforms.
// In case only waveforms are stored at positions 1,2,5,7 this memberfunction
// returns a value 7 whereas only 4 actual waveforms are present.
// This implies that when looping over the various waveform slots, one
// always has to check whether the returned pointer value is non-zero
// (which is a good practice anyhow).
 Int_t n=-1;
 if (fWaveforms) n=fWaveforms->GetLast();
 return (n+1);
}
///////////////////////////////////////////////////////////////////////////
TH1F* AliSignal::GetWaveform(Int_t j) const
{
// Provide pointer to the j-th waveform histogram.
 TH1F* waveform=0;
 if (j <= GetNwaveforms()) waveform=(TH1F*)fWaveforms->At(j-1);
 return waveform;
}
///////////////////////////////////////////////////////////////////////////
TH1F* AliSignal::GetWaveform(TString name) const
{
// Provide pointer to the waveform histogram with the specified name.
// In case no match is found, zero is returned.
 Int_t n=GetNwaveforms();
 TString str;
 for (Int_t i=1; i<=n; i++)
 {
  TH1F* waveform=GetWaveform(i);
  if (waveform)
  {
   str=waveform->GetName();
   if (str == name) return waveform;
  }
 }
 return 0; // No match found
}
///////////////////////////////////////////////////////////////////////////
Int_t AliSignal::GetWaveformIndex(TString name) const
{
// Provide index to the waveform histogram with the specified name.
// In case no match is found, zero is returned.
 Int_t n=GetNwaveforms();
 TString str;
 for (Int_t i=1; i<=n; i++)
 {
  TH1F* waveform=GetWaveform(i);
  if (waveform)
  {
   str=waveform->GetName();
   if (str == name) return i;
  }
 }
 return 0; // No match found
}
///////////////////////////////////////////////////////////////////////////
void AliSignal::SetWaveform(TH1F* waveform,Int_t j)
{
// Set the 1D waveform histogram for the j-th waveform.
//
// Notes :
//  The first waveform position at j=1.
//  j=1 is the default value.
//
// In case the value of the index j exceeds the maximum number of reserved
// positions for the waveforms, the number of reserved positions for the waveforms
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

 if (j<1) return;

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
// Reset the histogram of the j-th (default j=1) waveform.
// This memberfunction invokes TH1F::Reset() for the corresponding waveform(s).
// To actually delete the histograms from memory, use DeleteWaveform().
// Notes : The first position is at j=1.
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
void AliSignal::ResetWaveform(TString name)
{
// Reset the waveform with the specified name.
 Int_t j=GetWaveformIndex(name);
 if (j>0) ResetWaveform(j);
}
///////////////////////////////////////////////////////////////////////////
void AliSignal::DeleteWaveform(Int_t j)
{
// Delete the histogram of the j-th (default j=1) waveform.
// Notes : The first position is at j=1.
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
void AliSignal::DeleteWaveform(TString name)
{
// Delete the waveform with the specified name.
 Int_t j=GetWaveformIndex(name);
 if (j>0) DeleteWaveform(j);
}
///////////////////////////////////////////////////////////////////////////
Int_t AliSignal::GetNlinks(TObject* obj,Int_t j) const
{
// Provide the number of links to the specified object for the j-th slot.
// If j=0 (default) all slots will be scanned for the specified object.
// If obj=0 (default) all encountered objects for the specified slot will be counted.
// So, invokation of the default GetNlinks() will return the total number of
// all references to all sorts of stored objects.
 if (j<0)
 {
  cout << " *AliSignal::GetNlinks* Index j = " << j << " invalid." << endl;
  return 0;
 }

 if (!fLinks) return 0;

 Int_t n=0;
 if (!j)
 {
  n=fLinks->GetNrefs(obj);
 }
 else
 {
  TArrayI posarr;
  n=GetIndices(obj,j,posarr);
 }
 return n;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliSignal::GetNlinks(TObject* obj,TString name) const
{
// Provide the number of links to the specified object for the name-spec. slot.
// If obj=0 all encountered objects for the specified slot will be counted.
//
// This procedure involves a slot-index search based on the specified name
// at each invokation. This may become slow in case many slots have been
// defined and/or when this procedure is invoked many times.
// In such cases it is preferable to use indexed addressing in the user code
// either directly or via a few invokations of GetSlotIndex().

 Int_t j=GetSlotIndex(name);
 Int_t n=0;
 if (j>0) n=GetNlinks(obj,j);
 return n;
}
///////////////////////////////////////////////////////////////////////////
TObject* AliSignal::GetLink(Int_t j,Int_t k) const
{
// Provide pointer of the object linked to the j-th slot at position k.

 TObject* obj=0;
 // Note : In the internal storage matrix slots=columns positions=rows 
 if (fLinks) obj=fLinks->GetObject(k,j);
 return obj;
}
///////////////////////////////////////////////////////////////////////////
TObject* AliSignal::GetLink(TString name,Int_t k) const
{
// Provide pointer of the object linked to the name-spec. slot at position k.
//
// This procedure involves a slot-index search based on the specified name
// at each invokation. This may become slow in case many slots have been
// defined and/or when this procedure is invoked many times.
// In such cases it is preferable to use indexed addressing in the user code
// either directly or via a few invokations of GetSlotIndex().

 Int_t j=GetSlotIndex(name);
 TObject* obj=0;
 if (j>0) obj=GetLink(j,k);
 return obj;
}
///////////////////////////////////////////////////////////////////////////
void AliSignal::SetLink(TObject* obj,Int_t j,Int_t k)
{
// Introduce a link (=pointer) to an object for the j-th slot at position k.
// Only the pointer values are stored for (backward) reference, meaning
// that the objects of which the pointers are stored are NOT owned
// by the AliSignal object.
//
// Notes :
//  The first slot is at j=1 and the first position is at k=1.
//  j=1 and k=1 are the default values.
//
// If needed, the storage area for the links is increased automatically.
//
// In case the pointer argument is zero, indeed a value of zero will be
// stored at the specified position (k) for the specified slot (j).
//
// In principle any object derived from TObject can be referred to by this
// mechanism.
// However, this "linking back" facility was introduced to enable AliSignal slots
// to refer directly to the various AliTracks to which the AliSignal object itself
// is related (see AliTrack::AddSignal).
// Therefore, in case the input argument "obj" points to an AliTrack (or derived)
// object, the current signal is automatically related to this AliTrack
// (or derived) object.
// Also a global link to this AliTrack (or derived) object will be stored
// via the AddTrack() facility.
//
// IMPORTANT NOTE :
// ----------------
// In case one just wants to relate the current AliSignal to a certain AliTrack
// without a specific signal slot association, it is much more efficient
// (both memory and CPU wise) to use the memberfunction AddTrack() instead.
// 
// Please also have a look at the docs of the memberfunction ResetLink()
// to prevent the situation of stored pointers to non-existent object. 

 if (!fLinks && obj) fLinks=new AliObjMatrix();

 if (!fLinks) return;

 // Note : In the internal storage matrix slots=columns positions=rows 
 fLinks->EnterObject(k,j,obj);
 if (obj) 
 {
  if (obj->InheritsFrom("AliTrack"))
  {
   AliTrack* t=(AliTrack*)obj;
   AddTrack(*t,1);
  }
 }
}
///////////////////////////////////////////////////////////////////////////
void AliSignal::SetLink(TObject* obj,TString name,Int_t k)
{
// Introduce a link (=pointer) to an object for the name-spec. slot at position k.
// Only the pointer values are stored for (backward) reference, meaning
// that the objects of which the pointers are stored are NOT owned
// by the AliSignal object.
//
// This procedure involves a slot-index search based on the specified name
// at each invokation. This may become slow in case many slots have been
// defined and/or when this procedure is invoked many times.
// In such cases it is preferable to use indexed addressing in the user code
// either directly or via a few invokations of GetSlotIndex().
//
// In case the pointer argument is zero, indeed a value of zero will be
// stored at the specified position (k) for the specified slotname.
//
// In principle any object derived from TObject can be referred to by this
// mechanism.
// However, this "linking back" facility was introduced to enable AliSignal slots
// to refer directly to the various AliTracks to which the AliSignal object itself
// is related (see AliTrack::AddSignal).
// Therefore, in case the input argument "obj" points to an AliTrack (or derived)
// object, the current signal is automatically related to this AliTrack
// (or derived) object.
// Also a global link to this AliTrack (or derived) object will be stored
// via the AddTrack() facility.
//
// IMPORTANT NOTE :
// ----------------
// In case one just wants to relate the current AliSignal to a certain AliTrack
// without a specific signal slot association, it is much more efficient
// (both memory and CPU wise) to use the memberfunction AddTrack() instead.
// 
// Please also have a look at the docs of the memberfunction ResetLink()
// to prevent the situation of stored pointers to non-existent object. 

 Int_t j=GetSlotIndex(name);
 if (j>0) SetLink(obj,j,k);
}
///////////////////////////////////////////////////////////////////////////
void AliSignal::AddLink(TObject* obj,Int_t j)
{
// Introduce a link (=pointer) to an object for the j-th slot at the first
// free position.
// Only the pointer values are stored for (backward) reference, meaning
// that the objects of which the pointers are stored are NOT owned
// by the AliSignal object.
//
// Notes :
//  The first slot is at j=1 and the first position is at k=1.
//  j=1 is the default value.
//
// If needed, the storage area for the links is increased automatically.
//
// In case the pointer argument is zero, no link will be added.
//
// In principle any object derived from TObject can be referred to by this
// mechanism.
// However, this "linking back" facility was introduced to enable AliSignal slots
// to refer directly to the various AliTracks to which the AliSignal object itself
// is related (see AliTrack::AddSignal).
// Therefore, in case the input argument "obj" points to an AliTrack (or derived)
// object, the current signal is automatically related to this AliTrack
// (or derived) object.
// Also a global link to this AliTrack (or derived) object will be stored
// via the AddTrack() facility.
//
// IMPORTANT NOTE :
// ----------------
// In case one just wants to relate the current AliSignal to a certain AliTrack
// without a specific signal slot association, it is much more efficient
// (both memory and CPU wise) to use the memberfunction AddTrack() instead.
// 
// Please also have a look at the docs of the memberfunction ResetLink()
// to prevent the situation of stored pointers to non-existent object. 

 if (!obj || j<=0) return;

 if (!fLinks) fLinks=new AliObjMatrix();

 TObject* dum=0;
 Int_t n=GetNlinks(dum,j);
 Int_t pos=1;
 for (Int_t k=1; k<=n; k++)
 {
  dum=GetLink(j,k);
  if (!dum) break;
  pos++;
 }

 SetLink(obj,j,pos);
}
///////////////////////////////////////////////////////////////////////////
void AliSignal::AddLink(TObject* obj,TString name)
{
// Introduce a link (=pointer) to an object for the name-spec slot at the first
// free position.
// Only the pointer values are stored for (backward) reference, meaning
// that the objects of which the pointers are stored are NOT owned
// by the AliSignal object.
//
// This procedure involves a slot-index search based on the specified name
// at each invokation. This may become slow in case many slots have been
// defined and/or when this procedure is invoked many times.
// In such cases it is preferable to use indexed addressing in the user code
// either directly or via a few invokations of GetSlotIndex().
//
// In case the pointer argument is zero, indeed a value of zero will be
// stored at the first free position of the specified slotname.
//
// In principle any object derived from TObject can be referred to by this
// mechanism.
// However, this "linking back" facility was introduced to enable AliSignal slots
// to refer directly to the various AliTracks to which the AliSignal object itself
// is related (see AliTrack::AddSignal).
// Therefore, in case the input argument "obj" points to an AliTrack (or derived)
// object, the current signal is automatically related to this AliTrack
// (or derived) object.
// Also a global link to this AliTrack (or derived) object will be stored
// via the AddTrack() facility.
//
// IMPORTANT NOTE :
// ----------------
// In case one just wants to relate the current AliSignal to a certain AliTrack
// without a specific signal slot association, it is much more efficient
// (both memory and CPU wise) to use the memberfunction AddTrack() instead.
// 
// Please also have a look at the docs of the memberfunction ResetLink()
// to prevent the situation of stored pointers to non-existent object. 

 Int_t j=GetSlotIndex(name);
 if (j>0) AddLink(obj,j);
}
///////////////////////////////////////////////////////////////////////////
void AliSignal::ResetLink(Int_t j,Int_t k)
{
// Reset the link of the j-th slot at position k.
//
// Notes :
//  The first slot is at j=1 and the first position is at k=1.
//  j=1 and k=1 are the default values.
//
//  This memberfunction is intended to reset only 1 specified link location.
//  For extended functionality, please refer to the memberfuction ResetLinks().
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
 
 // Note : In the internal storage matrix slots=columns positions=rows 
 if (fLinks) fLinks->RemoveObject(k,j);
}
///////////////////////////////////////////////////////////////////////////
void AliSignal::ResetLink(TString name,Int_t k)
{
// Reset the link of the name-specified slot at position k.
//
// This memberfunction is intended to reset only 1 specified link location.
// For extended functionality, please refer to the memberfuction ResetLinks().
//
// This procedure involves a slot-index search based on the specified name
// at each invokation. This may become slow in case many slots have been
// defined and/or when this procedure is invoked many times.
// In such cases it is preferable to use indexed addressing in the user code
// either directly or via a few invokations of GetSlotIndex().

 Int_t j=GetSlotIndex(name);
 if (j>0) ResetLink(j,k);
}
///////////////////////////////////////////////////////////////////////////
void AliSignal::ResetLinks(TObject* obj,Int_t j,Int_t k)
{
// Reset single or multiple slot link(s) according to user specified selections.
//
// IMPORTANT NOTE :
// ----------------
// This facility only acts on the slot related links.
// The global track reference list will not be affected.
// To remove all references to AliTrack (or derived) objects, please
// use the RemoveTrack() of RemoveTracks() memberfunctions.
//
// A link is only reset if the stored reference matches the argument "obj".
// In case obj=0 no check on the matching of the stored reference is performed
// and the stored link is always reset in accordance with the other
// selection criteria.
//
// In case the slot argument "j" is specified, only the links from that
// specified slot will be deleted.
// In case j=0 (default) no checking on the slot index is performed.
//
// In case the position argument "k" is specified, only the links from that
// specified position will be deleted.
// In case k=0 (default) no checking on the position index is performed.
//
// So, invokation of ResetLinks(obj) will remove all references to the
// object "obj" from the total AliSignal, whereas ResetLinks(obj,j)
// will remove all references to the object "obj" only from slot "j".
//
// Notes :
// -------
// The first slot is indicated as j=1, whereas the first position is at k=1.
//
// Invokation of ResetLinks(0,row,col) is equivalent to invoking the
// memberfunction ResetLink(row,col).
// Invoking the latter directly is slightly faster.
//
// Invokation of ResetLinks(0) will reset all stored references in this AliSignal.
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
 
 if (!fLinks) return;

 if (!obj && !j && !k)
 {
  fLinks->Reset();
 }
 else
 {
  // Note : In the internal storage matrix slots=columns positions=rows 
  fLinks->RemoveObjects(obj,k,j);
 }
}
///////////////////////////////////////////////////////////////////////////
void AliSignal::ResetLinks(TObject* obj,TString name,Int_t k)
{
// Reset single or multiple slot link(s) according to user specified selections.
//
// IMPORTANT NOTE :
// ----------------
// This facility only acts on the slot related links.
// The global track reference list will not be affected.
// To remove all references to AliTrack (or derived) objects, please
// use the RemoveTrack() of RemoveTracks() memberfunctions.
//
// A link is only reset if the stored reference matches the argument "obj".
// In case obj=0 no check on the matching of the stored reference is performed
// and the stored link is always reset in accordance with the other
// selection criteria.
//
// In case the position argument "k" is specified, only the links from that
// specified position will be deleted.
// In case k=0 (default) no checking on the position index is performed.
//
// This procedure involves a slot-index search based on the specified name
// at each invokation. This may become slow in case many slots have been
// defined and/or when this procedure is invoked many times.
// In such cases it is preferable to use indexed addressing in the user code
// either directly or via a few invokations of GetSlotIndex().
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

 Int_t j=GetSlotIndex(name);
 if (j>0) ResetLinks(obj,j,k);
}
///////////////////////////////////////////////////////////////////////////
Int_t AliSignal::GetIndices(TObject* obj,TArrayI& js,TArrayI& ks) const
{
// Provide the slot and position indices of all the storage locations
// of the specified object.
// The slot (j) and pos. (k) indices are returned in the two separate TArrayI arrays
// from which the (j,k) pairs can be obtained from the corresponding
// array indices like (j,k)=(js.At(i),ks.At(i)).
// The integer return argument represents the number of (j,k) pairs which
// were encountered for the specified object.
//
// If obj=0 no object selection is performed and all (j,k) indices
// of the stored references for all objects are returned.
//
// Notes :
// -------
// As usual the convention is that slot and position numbering starts at 1.
// 
// This memberfunction always resets the two TArrayI arrays at the start.
//
// This memberfunction can only be used to obtain the (j,k) indices
// of the object as stored via the SetLink() or AddLink() memberfunction.
// This means that in case the user has entered a TObjArray as object
// (to increase the dimension of the resulting structure), the (j,k)
// indices of that TObjArray are obtained and NOT the indices of the
// actual objects contained in that TObjArray structure.
//
 Int_t nrefs=0;
 js.Reset();
 ks.Reset();
 // Note : In the internal storage matrix slots=columns positions=rows 
 if (fLinks) nrefs=fLinks->GetIndices(obj,ks,js);
 return nrefs;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliSignal::GetIndices(TObject* obj,Int_t j,TArrayI& ks) const
{
// Provide the position indices of all the storage locations of the
// specified object in the j-th slot of this AliSignal.
// The position indices are returned in the TArrayI array.
// The integer return argument represents the number of storage locations which
// were encountered for the specified object in the j-th slot.
//
// If obj=0 no object selection is performed and all position indices
// of the stored references for all objects of the j-th slot are returned.
//
// If j=0 all slots will be scanned and all position indices matching the
// object selection are returned.
// Note that in this case multiple appearances of the same position index
// will only be recorded once in the returned TArrayI array.
//
// Notes :
// -------
// As usual the convention is that slot and position numbering starts at 1.
// 
// This memberfunction always resets the TArrayI array at the start.
//
// This memberfunction can only be used to obtain the position indices
// of the object as stored via the SetLink() or AddLink() memberfunction.
// This means that in case the user has entered a TObjArray as object
// (to increase the dimension of the resulting structure), the position
// indices of that TObjArray are obtained and NOT the indices of the
// actual objects contained in that TObjArray structure.
//
 Int_t nrefs=0;
 ks.Reset();
 // Note : In the internal storage matrix slots=columns positions=rows 
 if (fLinks) nrefs=fLinks->GetIndices(obj,ks,j);
 return nrefs;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliSignal::GetIndices(TObject* obj,TString name,TArrayI& ks) const
{
// Provide the position indices of all the storage locations of the
// specified object in the name-specified slot of this AliSignal.
// The position indices are returned in the TArrayI array.
// The integer return argument represents the number of storage locations which
// were encountered for the specified object in the j-th slot.
//
// If obj=0 no object selection is performed and all position indices
// of the stored references for all objects of the j-th slot are returned.
//
// This procedure involves a slot-index search based on the specified name
// at each invokation. This may become slow in case many slots have been
// defined and/or when this procedure is invoked many times.
// In such cases it is preferable to use indexed addressing in the user code
// either directly or via a few invokations of GetSlotIndex().

 Int_t j=GetSlotIndex(name);
 Int_t n=0;
 if (j>0) n=GetIndices(obj,j,ks);
 return n;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliSignal::GetIndices(TObject* obj,TArrayI& js,Int_t k) const
{
// Provide the slot indices of all the storage locations of the
// specified object for the k-th position in this AliSignal.
// The slot indices are returned in the TArrayI array.
// The integer return argument represents the number of storage locations which
// were encountered for the specified object in the k-th position.
//
// If obj=0 no object selection is performed and all slot indices
// of the stored references for all objects in the k-th position are returned.
//
// If k=0 all positions will be scanned and all slot indices matching the
// object selection are returned.
// Note that in this case multiple appearances of the same slot index
// will only be recorded once in the returned TArrayI array.
//
// Notes :
// -------
// As usual the convention is that slot and position numbering starts at 1.
// 
// This memberfunction always resets the TArrayI array at the start.
//
// This memberfunction can only be used to obtain the slot indices
// of the object as stored via the SetLink() or AddLink() memberfunction.
// This means that in case the user has entered a TObjArray as object
// (to increase the dimension of the resulting structure), the slot
// indices of that TObjArray are obtained and NOT the indices of the
// actual objects contained in that TObjArray structure.
//
 Int_t nrefs=0;
 js.Reset();
 // Note : In the internal storage matrix slots=columns positions=rows 
 if (fLinks) nrefs=fLinks->GetIndices(obj,k,js);
 return nrefs;
}
///////////////////////////////////////////////////////////////////////////
void AliSignal::SetSwapMode(Int_t swap)
{
// Set swapmode flag for the internal link storage.
// In case for the stored links the maximum slot number differs considerably
// from the maximum position number, it might be more efficient
// (w.r.t. memory usage and/or output file size) to internally store the
// link reference matrix with the rows and colums swapped.
// This swapping is only related with the internal storage and as such
// is completely hidden for the user.
// At invokation of this memberfunction the default argument is swap=1.
//
// Note : The swap mode can only be set as long as no links are
//        stored in the AliSignal (i.e. a new instance of AliSignal
//        or after invokation of the Reset() or ResetLinks() function).
 
 if (!fLinks) fLinks=new AliObjMatrix();
 fLinks->SetSwapMode(swap);
}
///////////////////////////////////////////////////////////////////////////
Int_t AliSignal::GetSwapMode() const
{
// Provide swapmode flag of the link storage.
 Int_t swap=0; 
 if (fLinks) swap=fLinks->GetSwapMode();
 return swap;
}
///////////////////////////////////////////////////////////////////////////
void AliSignal::SetDevice(TObject* dev)
{
// Store the pointer to the device which owns this AliSignal object.
// This memberfunction is meant for internal use in AliDevice.
 fDevice=dev;
}
///////////////////////////////////////////////////////////////////////////
AliDevice* AliSignal::GetDevice() const
{
// Provide the pointer to the device which owns this AliSignal object.
 return (AliDevice*)fDevice;
}
///////////////////////////////////////////////////////////////////////////
void AliSignal::AddTrack(AliTrack& t,Int_t mode)
{
// Relate an AliTrack object to this signal.
// Only the pointer values are stored for (backward) reference, meaning
// that the tracks of which the pointers are stored are NOT owned
// by the AliSignal object.
//
// mode = 0 : Only the reference to the specified track is stored in
//            the current signal, without storing the (backward) reference
//            to this signal into the AliTrack structure. 
//        1 : The (backward) reference to the current signal is also automatically
//            stored into the AliTrack (or derived) object specified in the
//            input argument.
//
// The default is mode=1.

 if (!fTracks) fTracks=new TObjArray(1);

 // Check if this track is already stored for this signal
 Int_t ntk=GetNtracks();
 for (Int_t i=0; i<ntk; i++)
 {
  if (&t==fTracks->At(i)) return; 
 }

 fTracks->Add(&t);
 if (mode==1) t.AddSignal(*this,0);
}
///////////////////////////////////////////////////////////////////////////
void AliSignal::RemoveTrack(AliTrack& t,Int_t mode)
{
// Remove related AliTrack object from this signal.
// Also all references (if any) to this track in the slot links area
// are removed.
//
// mode = 0 : All references to the specified track are removed from
//            the current signal, without removing the (backward) reference
//            to this signal from the AliTrack structure. 
//        1 : The (backward) reference to the current signal is also automatically
//            removed from the AliTrack (or derived) object specified in the
//            input argument.
//
// The default is mode=1.

 if (fTracks)
 {
  AliTrack* test=(AliTrack*)fTracks->Remove(&t);
  if (test) fTracks->Compress();
 }

 ResetLinks(&t);

 if (mode==1) t.RemoveSignal(*this,0);
}
///////////////////////////////////////////////////////////////////////////
void AliSignal::RemoveTracks(Int_t mode)
{
// Remove all related AliTrack objects from this signal.
// Also all references (if any) to the related tracks in the slot links area
// are removed.
//
// mode = 0 : All track references are removed from the current signal,
//            without removing the (backward) references to this signal from
//            the corresponding AliTrack objects. 
//        1 : The (backward) references to the current signal are also automatically
//            removed from the corresponding AliTrack (or derived) objects.
//
// The default is mode=1.
 
 if (!fTracks) return;

 Int_t ntk=GetNtracks();
 for (Int_t i=0; i<ntk; i++)
 {
  AliTrack* tx=(AliTrack*)fTracks->At(i);
  if (tx)
  {
   ResetLinks(tx);
   if (mode==1) tx->RemoveSignal(*this,0);
  }
 }

 delete fTracks;
 fTracks=0;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliSignal::GetNtracks(AliTrack* t) const
{
// Provide the number of related AliTracks.
// In case an AliTrack pointer is specified as input argument,
// the number returned will be the number of occurrences (i.e. 0 or 1)
// for that specified track.
// By default t=0, which implies that just the number of all associated
// tracks will be returned.

 if (!fTracks) return 0;

 Int_t ntk=fTracks->GetEntries();

 if (!t) return ntk;

 for (Int_t i=0; i<ntk; i++)
 {
  AliTrack* tx=(AliTrack*)fTracks->At(i);
  if (tx==t) return 1;
 }

 return 0;
}
///////////////////////////////////////////////////////////////////////////
AliTrack* AliSignal::GetTrack(Int_t j) const
{
// Provide the related AliTrack number j.
// Note : j=1 denotes the first track.

 if (!fTracks) return 0;

 if ((j >= 1) && (j <= GetNtracks()))
 {
  return (AliTrack*)fTracks->At(j-1);
 }
 else
 {
  cout << " *AliSignal* track number : " << j << " out of range."
       << " Ntk = " << GetNtracks() << endl;
  return 0;
 }
}
///////////////////////////////////////////////////////////////////////////
AliTrack* AliSignal::GetIdTrack(Int_t id) const
{
// Return the track with user identifier "id" of this signal
 if (!fTracks) return 0;

 AliTrack* tx=0;
 for (Int_t i=0; i<GetNtracks(); i++)
 {
  tx=(AliTrack*)fTracks->At(i);
  if (id == tx->GetId()) return tx;
 }
 return 0; // No matching id found
}
///////////////////////////////////////////////////////////////////////////
TObject* AliSignal::Clone(const char* name) const
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
