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
// Class IceCleanHits
// TTask derived class to perform hit cleaning.
//
// The code in this processor is based on the algorithms as developed
// by Oladipo Fadiran and George Japaridze (Clark Atlanta University, USA).
//
// Criteria applied for Amanda modules :
// -------------------------------------
// 1) ADC within [min,max]  Default : [0.3,999999] PE
// 2) TOT within [min,max]  Default : electrical [125,2000] optical [20,2000] ns
// 3) abs(LE-Ttrig)<=win    Default : win=2250 TDC counts
//    where : LE=uncalibrated hit LE (i.e. TDC counts)   Ttrig=trigger pulse LE in TDC counts
// 4) At least one other hit within radius R and time difference dt
//    to remove isolated hits. Defaults : R=70 m  dt=500 ns
//
// The defaults of the various parameters can be changed by the corresponding
// Set memberfunctions.
//
// Information about the actual parameter settings can be found in the event
// structure itself via the device named "IceCleanHits".
//
// Concerning the trigger time :
// -----------------------------
// By default the trigger time is obtained automatically from the IceEvent structure
// via the device called "Trigger".
// The uncalibrated LE (i.e. TDC counts) of a specified trigger pulse is used.
// The user can impose a specific trigger name or time to be used
// by invokation of the memberfunctions SetTnameA or SetTtimeA, respectively.
// Specification of a negative trigger time will result in the automatic
// trigger time setting corresponding to the "main" trigger.
// By default the trigger time of the "main" trigger will be used.
//
// The hits which do not fullfill the criteria are flagged "dead" for the
// corresponding signal slot. This means they are still present in the
// IceEvent structure and are as such still accessible.
// It is left to the user to decide (based on the various "dead" flag settings)
// whether or not to use these hits in his/her reconstruction or analysis.
//
// Note : This processor only works properly on Time and ADC calibrated data.
//
//--- Author: Nick van Eijndhoven 13-oct-2005 Utrecht University
//- Modified: NvE $Date$ Utrecht University
///////////////////////////////////////////////////////////////////////////
 
#include "IceCleanHits.h"
#include "Riostream.h"

ClassImp(IceCleanHits) // Class implementation to enable ROOT I/O

IceCleanHits::IceCleanHits(const char* name,const char* title) : TTask(name,title)
{
// Default constructor.
 fEvt=0;
 fAdcminA=0.3;
 fAdcmaxA=999999;
 fTotminAE=125;
 fTotmaxAE=2000;
 fTotminAO=20;
 fTotmaxAO=2000;
 fRmaxA=70;
 fDtmaxA=500;
 fTwinA=2250;
 fTtimA=-1;
 fTnamA="main";
}
///////////////////////////////////////////////////////////////////////////
IceCleanHits::~IceCleanHits()
{
// Default destructor.
}
///////////////////////////////////////////////////////////////////////////
void IceCleanHits::SetAdcRangeA(Float_t min,Float_t max)
{
// Set Amanda ADC range in PE.
// The default for the maximum is 999999.
 fAdcminA=min;
 fAdcmaxA=max;
}
///////////////////////////////////////////////////////////////////////////
void IceCleanHits::SetTotRangeAE(Float_t min,Float_t max)
{
// Set Amanda electrical TOT range in ns.
// The default for the maximum is 2000.
 fTotminAE=min;
 fTotmaxAE=max;
}
///////////////////////////////////////////////////////////////////////////
void IceCleanHits::SetTotRangeAO(Float_t min,Float_t max)
{
// Set Amanda optical TOT range in ns.
// The default for the maximum is 2000.
 fTotminAO=min;
 fTotmaxAO=max;
}
///////////////////////////////////////////////////////////////////////////
void IceCleanHits::SetIsolationA(Float_t rmax,Float_t dtmax)
{
// Set Amanda isolation radius (in m) and time difference (in ns).
 fRmaxA=rmax;
 fDtmaxA=dtmax;
}
///////////////////////////////////////////////////////////////////////////
void IceCleanHits::SetTwindowA(Float_t dtmax)
{
// Set Amanda maximal trigger window (in TDC counts).
// Only hits which occur in [T-dtmax,T+dtmax] will be kept,
// where T indicates the trigger time in TDC counts.
// For the Amanda DAQ hardware, 1 TDC corresponds to about 1.04 ns. 
 fTwinA=dtmax;
}
///////////////////////////////////////////////////////////////////////////
void IceCleanHits::SetTtimeA(Float_t t)
{
// Set Amanda trigger time (in TDC counts).
// A negative value will induce automatic trigger time setting based
// on "main" trigger as recorded in the IceEvent structure.
 fTtimA=t;
 fTnamA="user";
 if (t<0) fTnamA="main";
}
///////////////////////////////////////////////////////////////////////////
void IceCleanHits::SetTnameA(TString name)
{
// Set Amanda trigger name.
// Specification of a non-existing trigger name will result in a trigger time
// value of 0.
 fTtimA=0;
 fTnamA=name;
}
///////////////////////////////////////////////////////////////////////////
void IceCleanHits::Exec(Option_t* opt)
{
// Implementation of the hit cleaning procedures.

 TString name=opt;
 AliJob* parent=(AliJob*)(gROOT->GetListOfTasks()->FindObject(name.Data()));

 if (!parent) return;

 fEvt=(IceEvent*)parent->GetObject("IceEvent");
 if (!fEvt) return;

 // Storage of the used parameters in the IceCleanHits device
 AliSignal params;
 params.SetNameTitle("IceCleanHits","IceCleanHits processor parameters");
 params.SetSlotName("AdcminA",1);
 params.SetSlotName("AdcmaxA",2);
 params.SetSlotName("TotminAE",3);
 params.SetSlotName("TotmaxAE",4);
 params.SetSlotName("TotminAO",5);
 params.SetSlotName("TotmaxAO",6);
 params.SetSlotName("RmaxA",7);
 params.SetSlotName("DtmaxA",8);
 params.SetSlotName("TwinA",9);
 params.SetSlotName("TtimA",10);

 params.SetSignal(fAdcminA,1);
 params.SetSignal(fAdcmaxA,2);
 params.SetSignal(fTotminAE,3);
 params.SetSignal(fTotmaxAE,4);
 params.SetSignal(fTotminAO,5);
 params.SetSignal(fTotmaxAO,6);
 params.SetSignal(fRmaxA,7);
 params.SetSignal(fDtmaxA,8);
 params.SetSignal(fTwinA,9);
 params.SetSignal(fTtimA,10);

 fEvt->AddDevice(params);

 Amanda();
 InIce();
 IceTop();
}
///////////////////////////////////////////////////////////////////////////
void IceCleanHits::Amanda()
{
// Hit cleaning for Amanda modules.

 // Trigger time setting according to the user selection.
 // By default the "main" trigger of the event will be used.
 if (fTnamA != "user")
 {
  AliDevice* tdev=(AliDevice*)fEvt->GetDevice("Trigger");
  if (tdev)
  {
   AliSignal* trig=tdev->GetHit(fTnamA);
   if (trig) fTtimA=trig->GetSignal("trig_pulse_le");
  }
 }

 // All Amanda OMs with a signal
 TObjArray* aoms=fEvt->GetDevices("IceAOM");

 // Local OM array with bad/dead OMs (as indicated via IceCalibrate) discarded
 TObjArray oms;
 IceAOM* omx=0;
 for (Int_t i=0; i<aoms->GetEntries(); i++)
 {
  omx=(IceAOM*)aoms->At(i);
  if (!omx) continue;
  if (omx->GetDeadValue("ADC") || omx->GetDeadValue("LE") || omx->GetDeadValue("TOT")) continue;
  oms.Add(omx);
 }

 // Local array with clean hits
 TObjArray hits;
 Int_t omid=0;
 Int_t clean=1;
 AliSignal* sx=0;
 Float_t adc,le,tot;
 for (Int_t iom=0; iom<oms.GetEntries(); iom++)
 {
  omx=(IceAOM*)oms.At(iom);
  if (!omx) continue;
  omid=omx->GetUniqueID();
  for (Int_t ih=1; ih<=omx->GetNhits(); ih++)
  {
   sx=omx->GetHit(ih);
   if (!sx) continue;
   adc=sx->GetSignal("ADC",7);
   le=sx->GetSignal("LE",7);
   tot=sx->GetSignal("TOT",7);

   // Remove hits with an ADC value outside the range
   if (adc<fAdcminA || adc>fAdcmaxA)
   {
    sx->SetDead("ADC");
    clean=0;
   }
   // Remove hits with a TOT value outside the range
   // Note : Different ranges for electrical and optical modules
   if (omid<303) // Electrical OMs
   {
    if (tot<fTotminAE || tot>fTotmaxAE)
    {
     sx->SetDead("TOT");
     clean=0;
    }
   }
   else // Optical OMs
   {
    if (tot<fTotminAO || tot>fTotmaxAO)
    {
     sx->SetDead("TOT");
     clean=0;
    }
   }
   // Remove hits that are outside the trigger time window.
   // Since the trigger time was determined from uncalibrated LE's
   // (to include cable length effects) the uncalibrated LE of each
   // hit should be used here as well. 
   le=sx->GetSignal("LE",-7);
   if (fabs(le-fTtimA)>fTwinA)
   {
     sx->SetDead("LE");
     clean=0;
   }
   // Store only the current clean hits in our local hit array
   // This will save CPU time for the isolation criterion 
   if (clean) hits.Add(sx);
  }
 }
 
 // Isolation cut
 // Only retain hits that have at least one other hit within a certain
 // radius and within a certain time window
 Int_t nhits=hits.GetEntries();
 AliSignal* sx1=0;
 AliSignal* sx2=0;
 Float_t t1,t2;
 IceAOM* omx1=0;
 IceAOM* omx2=0;
 AliPosition r1;
 AliPosition r2;
 Float_t dt,dr;
 Int_t iso;
 for (Int_t jh1=0; jh1<nhits; jh1++)
 {
  sx1=(AliSignal*)hits.At(jh1);
  if (!sx1) continue;
  iso=1;
  for (Int_t jh2=0; jh2<nhits; jh2++)
  {
   if (jh1==jh2)
   {
    iso=0;
    continue;
   }
   sx2=(AliSignal*)hits.At(jh2);
   if (!sx2) continue;
   t1=sx1->GetSignal("LE",7);
   t2=sx2->GetSignal("LE",7);
   dt=fabs(t2-t1);
   if (dt>fDtmaxA) continue;
   omx1=(IceAOM*)sx1->GetDevice();
   omx2=(IceAOM*)sx2->GetDevice();
   if (omx1 && omx2)
   {
    r1=omx1->GetPosition();
    r2=omx2->GetPosition();
    dr=r1.GetDistance(r2);
    if (dr>fRmaxA) continue;
    iso=0;
   }
  }
  if (iso) sx1->SetDead("LE");
 }   
}
///////////////////////////////////////////////////////////////////////////
void IceCleanHits::InIce()
{
// Hit cleaning for IceCube InIce DOMs.
}
///////////////////////////////////////////////////////////////////////////
void IceCleanHits::IceTop()
{
// Hit cleaning for IceTop DOMs.
}
///////////////////////////////////////////////////////////////////////////
