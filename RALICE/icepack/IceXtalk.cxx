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
// Class IceXtalk
// TTask derived class to perform cross talk hit correction.
//
// This task takes the current event in memory and uses the attached
// OM database to identify pairs of OMs which might induce cross talk.
// For each particular transmitter and receiver pair within the event
// a probability for cross talk induction is determined on basis of
// the hit data and the cross talk probability functions from the OM database.
// If this probability is above a certain threshold "pmin" and the difference
// in LE between transmitter and receiver signal is within the bounds as
// given by the Xtalk calibration data, the occurrence of a cross talk signal
// of "pe" photo-electrons in the receiver is assumed and this signal is
// lateron subtracted from the receiver ADC signal.
// The subtraction of the Xtalk signal from the various hits is delayed
// until the very end, since receiver OMs may also act as transmitter OMs
// for other pair combinations.
// In case a transmitter OM was flagged as a bad module, the check on the
// LE time difference will not be performed and a cross talk induction
// decision will be solely based on the corresponding probability, since
// the latter only depends on the un-calibrated signal amplitude of the
// transmitter hit.
// So, hits which have induced cross talk are not completely removed
// but their signal is corrected for the cross talk.
// This will prevent severe losses when studying UHE events, where cross
// talk effects are expected to have negligible effects on the observed
// large module signals.
//
// Note : The amount with which the ADC value was corrected is stored
//        in the signal as an ADC offset. This will allow easy investigation
//        of Xtalk corrected signals and also enables successive
//        applications of this Xtalk processor to investigate (systematic)
//        effects of various criteria.
//        Example : In case an amount of 1 pe was subtracted from the ADC
//                  value, the ADC offset will be set to 1. 
//
// The values of "pmin" and "pe" can be set by the user via the
// SetMinProb() and SetXtalkPE() memberfunctions.
// The default values are pmin=0.5 and pe=1.
//
// Information about the actual parameter settings can be found in the event
// structure itself via the device named "IceXtalk".
//
// In case this processor is followed by an ADC threshold hit cleaning
// procedure, hits which only contained cross talk can be efficiently
// skipped or removed from the event.
// In case one would like to always remove a hit which containes induced
// cross talk, one could set the "pe" parameter to some very large value.
// This will result in cross talk induced hits to get large negative ADC
// signals and can as such easily be skipped/removed afterwards. 
//
// Note : This processor only works properly on Time and ADC calibrated data.
//        In case no OM database has been specified for this processor,
//        no cross talk hit correction will be performed.
//
//--- Author: Nick van Eijndhoven 11-aug-2005 Utrecht University
//- Modified: NvE $Date$ Utrecht University
///////////////////////////////////////////////////////////////////////////
 
#include "IceXtalk.h"
#include "Riostream.h"

ClassImp(IceXtalk) // Class implementation to enable ROOT I/O

IceXtalk::IceXtalk(const char* name,const char* title) : TTask(name,title)
{
// Default constructor.
 fCalfile=0;
 fOmdb=0;
 fPmin=0.5;
 fPe=1;
}
///////////////////////////////////////////////////////////////////////////
IceXtalk::~IceXtalk()
{
// Default destructor.
 if (fCalfile)
 {
  delete fCalfile;
  fCalfile=0;
 }
}
///////////////////////////////////////////////////////////////////////////
void IceXtalk::SetOMdbase(AliObjMatrix* omdb)
{
// Set the pointer to the OM database.
// Note : this will overrule a previously attached database. 
 fOmdb=omdb;
}
///////////////////////////////////////////////////////////////////////////
void IceXtalk::SetCalibFile(TString name)
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
void IceXtalk::SetMinProb(Float_t pmin)
{
// Set the minimal probability for cross talk induction.
 fPmin=pmin;
}
///////////////////////////////////////////////////////////////////////////
void IceXtalk::SetXtalkPE(Float_t pe)
{
// Set the nominal Xtalk signal in photo-electron equivalent.
 fPe=pe;
}
///////////////////////////////////////////////////////////////////////////
void IceXtalk::Exec(Option_t* opt)
{
// Implementation of cross talk hit correction.

 if (!fOmdb) return;

 TString name=opt;
 AliJob* parent=(AliJob*)(gROOT->GetListOfTasks()->FindObject(name.Data()));

 if (!parent) return;

 IceEvent* evt=(IceEvent*)parent->GetObject("IceEvent");
 if (!evt) return;

 // Storage of the used parameters in the IceXtalk device
 AliSignal params;
 params.SetNameTitle("IceXtalk","IceXtalk processor parameters");
 params.SetSlotName("Pmin",1);
 params.SetSlotName("Pe",2);

 params.SetSignal(fPmin,1);
 params.SetSignal(fPe,2);

 evt->AddDevice(params);

 // All Amanda OMs with a signal
 TObjArray* mods=evt->GetDevices("IceAOM");

 Int_t nmods=mods->GetEntries();
 if (!nmods) return;

 TObjArray xhits; // Array with pointers to Xtalk hits to be corrected

 IceAOM* omt=0; // Transmitter OM
 IceAOM* omr=0; // Receiver OM
 TF1* fxtalk=0; // The Xtalk probability function
 Int_t idtrans=0;
 Int_t idrec=0;
 Float_t adct=0,let=0;
 Float_t p=0;
 Float_t adcr=0,ler=0;
 Int_t ibad=0;
 Float_t cpar=0,bpar=0,dlemin=0,dlemax=0,test=0;
 Float_t dle=0;
 Float_t sigcor=0;
 for (Int_t imod=0; imod<nmods; imod++)
 {
  omt=(IceAOM*)mods->At(imod);
  if (!omt) continue;
  idtrans=omt->GetUniqueID();

  // Also take bad transmitter modules into account
  ibad=omt->GetDeadValue("ADC");
  if (ibad) omt->SetAlive("ADC");

  // Check for corresponding receiver modules  
  for (Int_t jmod=0; jmod<nmods; jmod++)
  {
   // No Xtalk from a module to itself 
   if (jmod==imod) continue;

   omr=(IceAOM*)mods->At(jmod);
   if (!omr) continue;
   idrec=omr->GetUniqueID();

   fxtalk=(TF1*)fOmdb->GetObject(idtrans,idrec+1);
   if (!fxtalk) continue;

   // Determine Xtalk probability for each transmitter hit
   for (Int_t ithit=1; ithit<=omt->GetNhits(); ithit++)
   {
    AliSignal* st=omt->GetHit(ithit);
    if (!st) continue;

    // Eliminate possible previous Xtalk correction
    sigcor=st->GetOffset("ADC");
    st->AddSignal(sigcor,"ADC");
    st->ResetOffset("ADC");
    adct=st->GetSignal("ADC",-4); // Get uncalibrated amplitude
    dlemin=fxtalk->GetParameter("dLE-min");
    dlemax=fxtalk->GetParameter("dLE-max");
    // Protection against overflow in Xtalk prob. function
    cpar=fxtalk->GetParameter("C");
    bpar=fxtalk->GetParameter("B");
    test=(adct-cpar)/bpar;
    if (test>100)
    {
     p=1;
    }
    else if (test<-100)
    {
     p=0;
    }
    else
    {
     p=fxtalk->Eval(adct);
    }
    if (p<fPmin) continue;

    // Xtalk flagging for each receiver hit
    for (Int_t irhit=1; irhit<=omr->GetNhits(); irhit++)
    {
     AliSignal* sr=omr->GetHit(irhit);
     if (!sr) continue;

     // Check calibrated LE time differences
     if (!ibad)
     {
      let=st->GetSignal("LE",4);
      ler=sr->GetSignal("LE",4);
      dle=ler-let;
      if (dle<dlemin || dle>dlemax) continue;
     }
     // Register this receiver hit to be corrected
     xhits.Add(sr); 
    }
   } // end of transmitter hit loop
  } // end of receiver OM loop

  // Restore the original transmitter dead flag value
  if (ibad) omt->SetDead("ADC");
 } // end of transmitter OM loop
 
 // Perform the Xtalk signal correction to the registered receiver hits
 // The correction value is stored in the signal data as an offset value
 for (Int_t ix=0; ix<xhits.GetEntries(); ix++)
 {
  AliSignal* sx=(AliSignal*)xhits.At(ix);
  if (!sx) continue;

  // Eliminate possible previous Xtalk correction
  sigcor=sx->GetOffset("ADC");
  adcr=sx->GetSignal("ADC")+sigcor;
  sigcor=fPe;
  // If stored ADC data is un-calibrated, convert fPe to raw ADC value 
  AliSignal* parent=(AliSignal*)sx->GetDevice();
  if (parent)
  {
   if (parent->GetCalFunction("ADC"));
   {
    idrec=parent->GetUniqueID();
    omr=(IceAOM*)fOmdb->GetObject(idrec,1);
    TF1* fdecal=0;
    if (omr) fdecal=omr->GetDecalFunction("ADC");
    if (fdecal) sigcor=fdecal->Eval(fPe);
   }
  }
  adcr=adcr-sigcor;
  sx->SetSignal(adcr,"ADC");
  sx->SetOffset(sigcor,"ADC");
 }
}
///////////////////////////////////////////////////////////////////////////
