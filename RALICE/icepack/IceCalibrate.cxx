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

 IceGOM* ome=0; // The event OM pointer
 IceGOM* omd=0; // The database OM pointer
 Int_t id=0;
 Float_t sig=0;
 Float_t sigc=0;
 TF1* fcal=0;
 TF1* fdecal=0;

 Int_t ndev=evt->GetNdevices();
 for (Int_t idev=1; idev<=ndev; idev++)
 {
  ome=(IceGOM*)evt->GetDevice(idev);
  if (!ome) continue;

  id=ome->GetUniqueID();

  omd=0;
  if (fOmdb) omd=(IceGOM*)fOmdb->GetObject(id,1);

  // Set global OM constants
  if (omd)
  {
   ome->SetPosition((Ali3Vector)omd->GetPosition());
   for (Int_t isd=4; isd<17; isd++)
   {
    ome->SetSignal(omd->GetSignal(isd),isd);
   }
  }
  else
  {
   for (Int_t ise=8; ise<17; ise++)
   {
    ome->SetSignal(0,ise);
   }
   ome->SetSignal(1,8);
   ome->SetSignal(1,12);
   ome->SetSignal(1,15);
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
   sig=sx->GetSignal("ADC",-7); // Get uncalibrated signal
   fcal=0;
   fdecal=0;
   if (omd)
   {
    fcal=omd->GetCalFunction("ADC");
    fdecal=omd->GetDecalFunction("ADC");
   }
   if (fcal) // Store calibrated signal
   {
    sigc=fcal->Eval(sig);
    sx->SetSignal(sigc,"ADC");
    ome->SetCalFunction(0,"ADC");
    ome->SetDecalFunction(fdecal,"ADC");
   }
   else // Store uncalibrated signal
   {
    sx->SetSignal(sig,"ADC");
    ome->SetCalFunction(fcal,"ADC");
    ome->SetDecalFunction(0,"ADC");
   }

   // LE (TDC) (de)calibration
   sig=sx->GetSignal("LE",-7); // Get uncalibrated signal
   fcal=0;
   fdecal=0;
   if (omd)
   {
    fcal=omd->GetCalFunction("LE");
    fdecal=omd->GetDecalFunction("LE");
   }
   // Store the ADC independent (de)calibration function in the OM
   if (fcal)
   {
    ome->SetCalFunction(0,"LE");
    ome->SetDecalFunction(fdecal,"LE");
   }
   else
   {
    ome->SetCalFunction(fcal,"LE");
    ome->SetDecalFunction(0,"LE");
   }
   // Store the hit-specific ADC dependent (de)calibration function in the hit itself
   sx->SetCalFunction(fcal,"LE");
   sx->SetDecalFunction(fdecal,"LE");
   fcal=sx->GetCalFunction("LE");
   fdecal=sx->GetDecalFunction("LE");
   sigc=sx->GetSignal("ADC",-7);
   if (sigc>0)
   {
    if (fcal) fcal->SetParameter(3,sigc);
    if (fdecal) fdecal->SetParameter(3,sigc);
   }
   else
   {
    if (fcal) fcal->SetParameter(3,1.e20);
    if (fdecal) fdecal->SetParameter(3,1.e20);
   }
   if (fcal) // Store calibrated signal
   {
    sigc=fcal->Eval(sig);
    sx->SetSignal(sigc,"LE");
    sx->SetCalFunction(0,"LE");
    sx->SetDecalFunction(fdecal,"LE");
   }
   else // Store uncalibrated signal
   {
    sx->SetSignal(sig,"LE");
    sx->SetCalFunction(fcal,"LE");
    sx->SetDecalFunction(0,"LE");
   }

   // TOT (de)calibration
   sig=sx->GetSignal("TOT",-7); // Get uncalibrated signal
   fcal=0;
   fdecal=0;
   if (omd)
   {
    fcal=omd->GetCalFunction("TOT");
    fdecal=omd->GetDecalFunction("TOT");
   }
   if (fcal) // Store calibrated signal
   {
    sigc=fcal->Eval(sig);
    sx->SetSignal(sigc,"TOT");
    ome->SetCalFunction(0,"TOT");
    ome->SetDecalFunction(fdecal,"TOT");
   }
   else // Store uncalibrated signal
   {
    sx->SetSignal(sig,"TOT");
    ome->SetCalFunction(fcal,"TOT");
    ome->SetDecalFunction(0,"TOT");
   }
  }

  if (omd)
  {  
   if (omd->GetDeadValue("ADC")) ome->SetDead("ADC");
   if (omd->GetDeadValue("LE")) ome->SetDead("LE");
   if (omd->GetDeadValue("TOT")) ome->SetDead("TOT");
  }
 }
}
///////////////////////////////////////////////////////////////////////////
