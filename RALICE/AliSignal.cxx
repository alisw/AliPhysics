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
// s.SetPosition(pos,"car");
// s.SetPositionErrors(err,"car");
// s.SetSignal(signal);
// s.SetSignalError(error);
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
// q.SetSignal(signal,1);
// q.SetSignalError(error,1);
// signal=268.1; // e.g. ADC value of signal
// error=3.75;
// q.SetSignal(signal,2);
// q.SetSignalError(error,2);
// signal=23.7; // e.g. corresponding dE/dx value
// error=0.48;
// q.SetSignal(signal,3);
// q.SetSignalError(error,3);
//
//--- Author: Nick van Eijndhoven 23-jan-1999 UU-SAP Utrecht
//- Modified: NvE $Date$ UU-SAP Utrecht
///////////////////////////////////////////////////////////////////////////

#include "AliSignal.h"
#include "Riostream.h"
 
ClassImp(AliSignal) // Class implementation to enable ROOT I/O
 
AliSignal::AliSignal() : TObject(),AliPosition()
{
// Creation of an AliSignal object and initialisation of parameters.
// Several values (with errors) can be stored.
// If needed, the storage for values (and errors) will be expanded automatically
// when entering values and/or errors.
 fSignal=0;
 fDsignal=0;
 fName="Unspecified";
 fHwaveform=0;
}
///////////////////////////////////////////////////////////////////////////
AliSignal::~AliSignal()
{
// Destructor to delete dynamically allocated memory
 if (fSignal)
 {
  delete fSignal;
  fSignal=0;
 }
 if (fDsignal)
 {
  delete fDsignal;
  fDsignal=0;
 }
 if (fHwaveform)
 {
  delete fHwaveform;
  fHwaveform=0;
 }
}
///////////////////////////////////////////////////////////////////////////
AliSignal::AliSignal(AliSignal& s) : TObject(s),AliPosition(s)
{
// Copy constructor
 fSignal=0;
 fDsignal=0;
 fName=s.fName;
 fHwaveform=0;

 Int_t nvalues=s.GetNvalues();
 Double_t sig;
 for (Int_t i=1; i<=nvalues; i++)
 {
  sig=s.GetSignal(i);
  SetSignal(sig,i);
 } 

 Int_t nerrors=s.GetNerrors();
 Double_t err;
 for (Int_t j=1; j<=nerrors; j++)
 {
  err=s.GetSignalError(j);
  SetSignalError(err,j);
 }

 TH1F* hist=s.GetWaveform();
 if (hist) fHwaveform=new TH1F(*hist); 
}
///////////////////////////////////////////////////////////////////////////
void AliSignal::Reset(Int_t mode)
{
// Reset all signal and position values and errors to 0.
//
// mode = 0 Reset position and all signal values and their errors to 0.
//          The waveform histogram is reset.
//        1 Reset position and delete the signal and error storage arrays.
//          The waveform histogram is deleted.
//
// The default when invoking Reset() corresponds to mode=0.
//
// The usage of mode=0 allows to re-use the allocated memory for new
// signal (and error) values. This behaviour is preferable (i.e. faster)
// in case the various signals always contain the same number of values.
// The usage of mode=1 is slower, but allows a more efficient memory
// occupation (and smaller output file size) in case the different
// signals have a variable number of values.
//
// For more specific actions see ResetPosition(), ResetSignals()
// and DeleteSignals().
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
// Irrespective of the mode, the waveform histogram is reset.

 if (mode<0 || mode>2)
 {
  cout << " *AliSignal::ResetSignals* Invalid argument mode = " << mode << endl;
  cout << " Default mode=0 will be used." << endl;
  mode=0;
 }

 if (fSignal && (mode==0 || mode==1))
 {
  for (Int_t i=0; i<fSignal->GetSize(); i++)
  {
   fSignal->AddAt(0,i);
  }
 }

 if (fDsignal && (mode==0 || mode==2))
 {
  for (Int_t j=0; j<fDsignal->GetSize(); j++)
  {
   fDsignal->AddAt(0,j);
  }
 }

 if (fHwaveform) fHwaveform->Reset();
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
// Irrespective of the mode, the waveform histogram is deleted.

 if (mode<0 || mode>2)
 {
  cout << " *AliSignal::DeleteSignals* Invalid argument mode = " << mode << endl;
  cout << " Default mode=0 will be used." << endl;
  mode=0;
 }

 if (fSignal && (mode==0 || mode==1))
 {
  delete fSignal;
  fSignal=0;
 }

 if (fDsignal && (mode==0 || mode==2))
 {
  delete fDsignal;
  fDsignal=0;
 }

 if (fHwaveform)
 {
  delete fHwaveform;
  fHwaveform=0;
 }
}
///////////////////////////////////////////////////////////////////////////
void AliSignal::ResetPosition()
{
// Reset the position and corresponding errors to 0.
 Double_t r[3]={0,0,0};
 SetPosition(r,"sph");
 SetErrors(r,"car");
}
///////////////////////////////////////////////////////////////////////////
void AliSignal::SetSignal(Double_t sig,Int_t j)
{
// Store j-th (default j=1) signal value.
// Note : The first signal value is at j=1.
// In case the value of the index j exceeds the maximum number of reserved
// slots for signal values, the number of reserved slots for the
// signal values is increased automatically.

 if (!fSignal)
 {
  fSignal=new TArrayF(j);
  ResetSignals(1);
 }

 Int_t size=fSignal->GetSize();

 if (j>size)
 {
  fSignal->Set(j);
 }

 fSignal->AddAt(float(sig),j-1);
}
///////////////////////////////////////////////////////////////////////////
void AliSignal::AddSignal(Double_t sig,Int_t j)
{
// Add value to j-th (default j=1) signal value.
// Note : The first signal value is at j=1.
// In case the value of the index j exceeds the maximum number of reserved
// slots for signal values, the number of reserved slots for the
// signal values is increased automatically.

 if (!fSignal)
 {
  fSignal=new TArrayF(j);
  ResetSignals(1);
 }

 Int_t size=fSignal->GetSize();

 if (j>size)
 {
  fSignal->Set(j);
 }

 Float_t sum=(fSignal->At(j-1))+sig;
 fSignal->AddAt(sum,j-1);
}
///////////////////////////////////////////////////////////////////////////
Float_t AliSignal::GetSignal(Int_t j)
{
// Provide j-th (default j=1) signal value.
// Note : The first signal value is at j=1.
// In case no signal is present or the argument j is invalid, 0 is returned.
 Float_t sig=0;
 if (fSignal)
 {
  if (j>0 && j<=(fSignal->GetSize()))
  {
   sig=fSignal->At(j-1);
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
// Store error on j-th (default j=1) signal value.
// Note : The error on the first signal value is at j=1.
// In case the value of the index j exceeds the maximum number of reserved
// slots for signal error values, the number of reserved slots for the
// signal errors is increased automatically.

 if (!fDsignal)
 {
  fDsignal=new TArrayF(j);
  ResetSignals(2);
 }

 Int_t size=fDsignal->GetSize();

 if (j>size)
 {
  fDsignal->Set(j);
 }

 fDsignal->AddAt(float(dsig),j-1);
}
///////////////////////////////////////////////////////////////////////////
Float_t AliSignal::GetSignalError(Int_t j)
{
// Provide error on the j-th (default j=1) signal value.
// Note : The error on the first signal value is at j=1.
// In case no signal is present or the argument j is invalid, 0 is returned.
 Float_t err=0;
 if (fDsignal)
 {
  if (j>0 && j<=(fDsignal->GetSize()))
  {
   err=fDsignal->At(j-1);
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
// Provide signal information within the coordinate frame f
 cout << " *AliSignal::Data* Signal of kind : " << fName.Data() << endl;
 cout << " Position";
 Ali3Vector::Data(f);

 Int_t nvalues=GetNvalues();
 Int_t nerrors=GetNerrors();

 if (fSignal)
 {
  for (Int_t i=0; i<nvalues; i++)
  {
   cout << "   Signal value : " << fSignal->At(i);
   if (fDsignal && i<nerrors) cout << " error : " << fDsignal->At(i);
   cout << endl;
  }
 }
} 
///////////////////////////////////////////////////////////////////////////
void AliSignal::SetName(TString name)
{
// Set the name tag to indicate the kind of signal.
 fName=name;
}
///////////////////////////////////////////////////////////////////////////
TString AliSignal::GetName()
{
// Provide the name tag indicating the kind of signal.
 return fName;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliSignal::GetNvalues()
{
// Provide the number of values for this signal.
 Int_t n=0;
 if (fSignal) n=fSignal->GetSize();
 return n;
}
///////////////////////////////////////////////////////////////////////////
Int_t AliSignal::GetNerrors()
{
// Provide the number specified errors on the values for this signal.
 Int_t n=0;
 if (fDsignal) n=fDsignal->GetSize();
 return n;
}
///////////////////////////////////////////////////////////////////////////
TH1F* AliSignal::GetWaveform()
{
// Provide pointer to the 1D waveform histogram.
 return fHwaveform;
}
///////////////////////////////////////////////////////////////////////////
void AliSignal::SetWaveform(TH1F* waveform)
{
// Set the 1D waveform histogram.
//
// In case the input argument has the same value as the current waveform
// histogram pointer value, no action is taken since the user has already
// modified the actual histogram.
//
// In case the input argument is zero, the current waveform histogram
// is deleted and the pointer set to zero.
//
// In all other cases the current waveform histogram is deleted and a new
// copy of the input histogram is created which becomes the current waveform
// histogram.

 if (waveform != fHwaveform)
 {
  if (fHwaveform)
  {
   delete fHwaveform;
   fHwaveform=0;
  }
  if (waveform) fHwaveform=new TH1F(*waveform);
 } 
}
///////////////////////////////////////////////////////////////////////////
AliSignal* AliSignal::MakeCopy(AliSignal& s)
{
// Make a deep copy of the input object and provide the pointer to the copy.
// This memberfunction enables automatic creation of new objects of the
// correct type depending on the argument type, a feature which may be very useful
// for containers when adding objects in case the container owns the objects.
// This feature allows e.g. AliTrack to store either AliSignal objects or
// objects derived from AliSignal via the AddSignal memberfunction, provided
// these derived classes also have a proper MakeCopy memberfunction. 

 AliSignal* sig=new AliSignal(s);
 return sig;
}
///////////////////////////////////////////////////////////////////////////
