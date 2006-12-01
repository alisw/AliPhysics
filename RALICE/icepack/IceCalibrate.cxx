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
// Class IceCalibrate
// TTask derived class to perform the various calibrations.
//
// This task takes the current event in memory and uses the attached
// OM database to access the various calibration functions.
// A specific OM database may be attached by means of the SetOMdbase()
// or SetCalibFile() memberfunctions.
// Further details about the OM database can be found in the docs
// of IceCal2Root.
//
// In the calibration procedure, all event data in memory is scanned and
// replaced by calibrated data if a calibration function is present.
// When data is succesfully calibrated, the corresponding de-calibration
// function is stored in the event data at the appropriate place to allow
// access to uncalibrated data as well (see AliSignal::GetSignal for 
// further details).
// When the input event in memory already contained calibrated data
// (i.e. de-calibration functions are present in the event data), the event
// data is first de-calibrated (using the corresponding de-calibration functions
// contained in the event data) before the new calibration is performed.
// In case no corresponding calibration function is present, the calibration
// of those specific data will not be performed.
// This implies that running this task on calibrated data without having
// attached an OM database, will result in fully de-calibrated data. 
// In case an OM slot was flagged as bad in the OM database, this flag
// will be copied into the event data for the corresponding OM.
//
// Information about the actual parameter settings can be found in the event
// structure itself via the device named "IceCalibrate".
//
//--- Author: Nick van Eijndhoven 18-sep-2005 Utrecht University
//- Modified: NvE $Date$ Utrecht University
///////////////////////////////////////////////////////////////////////////
 
#include "IceCalibrate.h"
#include "Riostream.h"

ClassImp(IceCalibrate) // Class implementation to enable ROOT I/O

IceCalibrate::IceCalibrate(const char* name,const char* title) : TTask(name,title)
{
// Default constructor.
 fCalfile=0;
 fOmdb=0;
}
///////////////////////////////////////////////////////////////////////////
IceCalibrate::~IceCalibrate()
{
// Default destructor.
 if (fCalfile)
 {
  delete fCalfile;
  fCalfile=0;
 }
}
///////////////////////////////////////////////////////////////////////////
void IceCalibrate::SetOMdbase(AliObjMatrix* omdb)
{
// Set the pointer to the OM database.
// Note : this will overrule a previously attached database. 
 fOmdb=omdb;
}
///////////////////////////////////////////////////////////////////////////
void IceCalibrate::SetCalibFile(TString name)
{
// Set the calibration ROOT file as created with IceCal2Root.
// Note : this will overrule a previously attached database. 
 fCalfile=new TFile(name.Data());
 if (fCalfile)
 {
  fOmdb=(AliObjMatrix*)fCalfile->Get("Cal-OMDBASE");
 }
}
///////////////////////////////////////////////////////////////////////////
void IceCalibrate::Exec(Option_t* opt)
{
// Implementation of the various calibration procedures.

 TString name=opt;
 AliJob* parent=(AliJob*)(gROOT->GetListOfTasks()->FindObject(name.Data()));

 if (!parent) return;

 IceEvent* evt=(IceEvent*)parent->GetObject("IceEvent");
 if (!evt) return;

 // Storage of the used parameters in the IceCalibrate device
 AliSignal params;
 params.SetNameTitle("IceCalibrate","IceCalibrate processor parameters");
 params.SetSlotName("Omdb",1);

 if (fOmdb) params.SetSignal(1,1);

 evt->AddDevice(params);

 // All OMs with a signal
 TObjArray* mods=evt->GetDevices("IceGOM");

 Int_t nmods=mods->GetEntries();
 if (!nmods) return;

 IceGOM* ome=0; // The event OM pointer
 IceGOM* omd=0; // The database OM pointer
 Int_t id=0;
 Float_t adc=0,le=0,tot=0;    // Uncalibrated values
 Float_t cadc=0,cle=0,ctot=0; // Calibrated values
 TF1* fcal=0;
 TF1* fdecal=0;
 TString slotname;

 for (Int_t imod=0; imod<nmods; imod++)
 {
  ome=(IceGOM*)mods->At(imod);
  if (!ome) continue;

  id=ome->GetUniqueID();

  omd=0;
  if (fOmdb) omd=(IceGOM*)fOmdb->GetObject(id,1);

  // Set global OM constants
  if (omd)
  {
   ome->SetPosition((Ali3Vector&)omd->GetPosition());
   for (Int_t ind=1; ind<=omd->GetNnames(); ind++)
   {
    slotname=omd->GetSlotName(ind);
    ome->AddNamedSlot(slotname);
   }
   for (Int_t isd=1; isd<=omd->GetNvalues(); isd++)
   {
    ome->SetSignal(omd->GetSignal(isd),omd->GetSlotName(isd));
   }
  }

  // Make signals of bad modules available
  ome->SetAlive("ADC");
  ome->SetAlive("LE");
  ome->SetAlive("TOT");

  // (De)calibrate all hit signals for this OM
  for (Int_t ithit=1; ithit<=ome->GetNhits(); ithit++)
  {
   AliSignal* sx=ome->GetHit(ithit);
   if (!sx) continue;

   // ADC (de)calibration
   adc=sx->GetSignal("ADC",-7); // Get uncalibrated signal
   fcal=0;
   fdecal=0;
   if (omd)
   {
    fcal=omd->GetCalFunction("ADC");
    fdecal=omd->GetDecalFunction("ADC");
   }
   if (fcal) // Store calibrated signal
   {
    cadc=fcal->Eval(adc);
    sx->SetSignal(cadc,"ADC");
   }
   else // Store uncalibrated signal
   {
    sx->SetSignal(adc,"ADC");
   }

   // LE (TDC) (de)calibration
   le=sx->GetSignal("LE",-7); // Get uncalibrated signal
   fcal=0;
   fdecal=0;
   if (omd)
   {
    fcal=omd->GetCalFunction("LE");
    fdecal=omd->GetDecalFunction("LE");
   }
   // Store the hit-specific ADC dependent (de)calibration function in the hit itself
   sx->SetCalFunction(fcal,"LE");
   sx->SetDecalFunction(fdecal,"LE");
   fcal=sx->GetCalFunction("LE");
   fdecal=sx->GetDecalFunction("LE");
   if (adc>0)
   {
    if (fcal) fcal->SetParameter(3,adc);
    if (fdecal) fdecal->SetParameter(3,adc);
   }
   else
   {
    if (fcal) fcal->SetParameter(3,1.e20);
    if (fdecal) fdecal->SetParameter(3,1.e20);
   }
   if (fcal) // Store calibrated signal
   {
    cle=fcal->Eval(le);
    sx->SetSignal(cle,"LE");
    sx->SetCalFunction(0,"LE");
    sx->SetDecalFunction(fdecal,"LE");
   }
   else // Store uncalibrated signal
   {
    sx->SetSignal(le,"LE");
    sx->SetCalFunction(fcal,"LE");
    sx->SetDecalFunction(0,"LE");
   }

   // TOT (de)calibration
   tot=sx->GetSignal("TOT",-7); // Get uncalibrated signal
   fcal=0;
   fdecal=0;
   if (omd)
   {
    fcal=omd->GetCalFunction("TOT");
    fdecal=omd->GetDecalFunction("TOT");
   }
   if (fcal) // Store calibrated signal
   {
    ctot=fcal->Eval(tot);
    sx->SetSignal(ctot,"TOT");
   }
   else // Store uncalibrated signal
   {
    sx->SetSignal(tot,"TOT");
   }
  } // End of loop over hits of the OM

  // Set bad OM flags according to dbase info
  if (omd)
  {  
   if (omd->GetDeadValue("ADC")) ome->SetDead("ADC");
   if (omd->GetDeadValue("LE")) ome->SetDead("LE");
   if (omd->GetDeadValue("TOT")) ome->SetDead("TOT");
  }

  // Store ADC (de)calibration function in this OM according to dbase info
  fcal=0;
  fdecal=0;
  if (omd)
  {
   fcal=omd->GetCalFunction("ADC");
   fdecal=omd->GetDecalFunction("ADC");
  }
  if (fcal) // Calibrated ADC signals were stored
  {
   ome->SetCalFunction(0,"ADC");
   ome->SetDecalFunction(fdecal,"ADC");
  }
  else // Uncalibrated ADC signals were stored
  {
   ome->SetCalFunction(fcal,"ADC");
   ome->SetDecalFunction(0,"ADC");
  }

  // Store ADC independent LE (de)calibration function in this OM according to dbase info
  fcal=0;
  fdecal=0;
  if (omd)
  {
   fcal=omd->GetCalFunction("LE");
   fdecal=omd->GetDecalFunction("LE");
  }
  if (fcal) // Calibrated LE signals were stored
  {
   ome->SetCalFunction(0,"LE");
   ome->SetDecalFunction(fdecal,"LE");
  }
  else // Uncalibrated LE signals were stored
  {
   ome->SetCalFunction(fcal,"LE");
   ome->SetDecalFunction(0,"LE");
  }

  // Store TOT (de)calibration function in this OM according to dbase info
  fcal=0;
  fdecal=0;
  if (omd)
  {
   fcal=omd->GetCalFunction("TOT");
   fdecal=omd->GetDecalFunction("TOT");
  }
  if (fcal) // Calibrated TOT signals were stored
  {
   ome->SetCalFunction(0,"TOT");
   ome->SetDecalFunction(fdecal,"TOT");
  }
  else // Uncalibrated TOT signals were stored
  {
   ome->SetCalFunction(fcal,"TOT");
   ome->SetDecalFunction(0,"TOT");
  }
 } // End of loop over OM's of the event
}
///////////////////////////////////////////////////////////////////////////
