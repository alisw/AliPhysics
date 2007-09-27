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
// In case an event has been rejected by an AliEventSelector (based) processor,
// this task (and its sub-tasks) is not executed.
//
// The code in this processor is based on the algorithms as developed
// by Oladipo Fadiran and George Japaridze (Clark Atlanta University, USA).
//
// Criteria applied for Amanda MuDaq data :
// ----------------------------------------
// 1) ADC within [min,max]  Default : [0.3,999999] PE
// 2) TOT within [min,max]  Default : electrical [125,2000] optical [20,2000] ns
// 3) abs(LE-Ttrig)<=win    Default : win=2250 TDC counts
//    where : LE=uncalibrated hit LE (i.e. TDC counts)   Ttrig=trigger pulse LE in TDC counts
// 4) At least one hit in another OM within radius R and time difference dt
//    to remove isolated hits. Defaults : R=70 m  dt=500 ns
//
// Criteria applied for Amanda TWRDaq data :
// ----------------------------------------
// 1) ADC within [min,max]  Default : [0.3,999999] PE
// 2) TOT within [min,max]  Default : electrical [125,2000] optical [20,2000] ns
// 3) abs(LE-Ttrig)<=win    Default : win=3000 ns
//    where : LE=uncalibrated hit LE    Ttrig=trigger pulse LE
// 4) At least one hit in another OM within radius R and time difference dt
//    to remove isolated hits. Defaults : R=70 m  dt=500 ns
//
// The actual DAQ system is obtained automatically from the IceEvent structure
// via the device called "Daq".
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
// The uncalibrated LE of a specified trigger pulse is used.
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
 fAdcminAM=0.3;
 fAdcmaxAM=999999;
 fAdcminAT=0.3;
 fAdcmaxAT=999999;
 fTotminAEM=125;
 fTotmaxAEM=2000;
 fTotminAOM=20;
 fTotmaxAOM=2000;
 fTotminAET=125;
 fTotmaxAET=2000;
 fTotminAOT=20;
 fTotmaxAOT=2000;
 fRmaxA=70;
 fDtmaxA=500;
 fTwinAM=2250;
 fTwinAT=3000;
 fTtimAM=-1;
 fTnamAM="main";
 fTtimAT=-1;
 fTnamAT="main";
}
///////////////////////////////////////////////////////////////////////////
IceCleanHits::~IceCleanHits()
{
// Default destructor.
}
///////////////////////////////////////////////////////////////////////////
void IceCleanHits::SetAdcRangeA(Float_t min,Float_t max,TString s)
{
// Set Amanda ADC range in PE.
// The default for the maximum is 999999.
// The argument "s" specifies the DAQ system (i.e. "MuDaq" or "TWRDaq").
// For backward compatibility the default is s="MuDaq".

 if (s=="MuDaq")
 {
  fAdcminAM=min;
  fAdcmaxAM=max;
 }
 if (s=="TWRDaq")
 {
  fAdcminAT=min;
  fAdcmaxAT=max;
 }
}
///////////////////////////////////////////////////////////////////////////
void IceCleanHits::SetTotRangeAE(Float_t min,Float_t max,TString s)
{
// Set Amanda electrical TOT range in ns.
// The default for the maximum is 2000.
// The argument "s" specifies the DAQ system (i.e. "MuDaq" or "TWRDaq").
// For backward compatibility the default is s="MuDaq".

 if (s=="MuDaq")
 {
  fTotminAEM=min;
  fTotmaxAEM=max;
 }
 if (s=="TWRDaq")
 {
  fTotminAET=min;
  fTotmaxAET=max;
 }
}
///////////////////////////////////////////////////////////////////////////
void IceCleanHits::SetTotRangeAO(Float_t min,Float_t max,TString s)
{
// Set Amanda optical TOT range in ns.
// The default for the maximum is 2000.
// The argument "s" specifies the DAQ system (i.e. "MuDaq" or "TWRDaq").
// For backward compatibility the default is s="MuDaq".

 if (s=="MuDaq")
 {
  fTotminAOM=min;
  fTotmaxAOM=max;
 }
 if (s=="TWRDaq")
 {
  fTotminAOT=min;
  fTotmaxAOT=max;
 }
}
///////////////////////////////////////////////////////////////////////////
void IceCleanHits::SetIsolationA(Float_t rmax,Float_t dtmax)
{
// Set Amanda isolation radius (in m) and time difference (in ns).
 fRmaxA=rmax;
 fDtmaxA=dtmax;
}
///////////////////////////////////////////////////////////////////////////
void IceCleanHits::SetTwindowA(Float_t dtmax,TString s)
{
// Set Amanda maximal trigger window.
// Only hits which occur in [T-dtmax,T+dtmax] will be kept,
// where T indicates the trigger time.
// For the Amanda MuDaq hardware, the times are all in TDC counts,
// where 1 TDC corresponds to about 1.04 ns. 
// For the TWRDaq, the times are all in nanoseconds.
// The argument "s" specifies the DAQ system (i.e. "MuDaq" or "TWRDaq").
// For backward compatibility the default is s="MuDaq".

 if (s=="MuDaq") fTwinAM=dtmax;
 if (s=="TWRDaq") fTwinAT=dtmax;
}
///////////////////////////////////////////////////////////////////////////
void IceCleanHits::SetTtimeA(Float_t t,TString s)
{
// Set Amanda trigger time.
// For the Amanda MuDaq hardware, the times are all in TDC counts,
// where 1 TDC corresponds to about 1.04 ns. 
// For the TWRDaq, the times are all in nanoseconds.
// The argument "s" specifies the DAQ system (i.e. "MuDaq" or "TWRDaq").
// For backward compatibility the default is s="MuDaq".
// A negative value will induce automatic trigger time setting based
// on "main" trigger as recorded in the IceEvent structure.

 if (s=="MuDaq")
 {
  fTtimAM=t;
  fTnamAM="user";
  if (t<0) fTnamAM="main";
 }
 if (s=="TWRDaq")
 {
  fTtimAT=t;
  fTnamAT="user";
  if (t<0) fTnamAT="main";
 }
}
///////////////////////////////////////////////////////////////////////////
void IceCleanHits::SetTnameA(TString name,TString s)
{
// Set Amanda trigger name.
// The argument "s" specifies the DAQ system (i.e. "MuDaq" or "TWRDaq").
// For backward compatibility the default is s="MuDaq".
// Specification of a non-existing trigger name will result in a trigger time
// value of 0.

 if (s=="MuDaq")
 {
  fTtimAM=0;
  fTnamAM=name;
 }
 if (s=="TWRDaq")
 {
  fTtimAT=0;
  fTnamAT=name;
 }
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

 // Only process accepted events
 AliDevice* seldev=(AliDevice*)fEvt->GetDevice("AliEventSelector");
 if (seldev)
 {
  if (seldev->GetSignal("Select") < 0.1) return;
 }

 // Storage of the used parameters in the IceCleanHits device
 AliSignal params;
 params.SetNameTitle("IceCleanHits","IceCleanHits processor parameters");
 params.SetSlotName("AdcminAM",1);
 params.SetSlotName("AdcmaxAM",2);
 params.SetSlotName("TotminAEM",3);
 params.SetSlotName("TotmaxAEM",4);
 params.SetSlotName("TotminAOM",5);
 params.SetSlotName("TotmaxAOM",6);
 params.SetSlotName("RmaxA",7);
 params.SetSlotName("DtmaxA",8);
 params.SetSlotName("TwinAM",9);
 params.SetSlotName("TtimAM",10);
 params.SetSlotName("AdcminAT",11);
 params.SetSlotName("AdcmaxAT",12);
 params.SetSlotName("TotminAET",13);
 params.SetSlotName("TotmaxAET",14);
 params.SetSlotName("TotminAOT",15);
 params.SetSlotName("TotmaxAOT",16);
 params.SetSlotName("TwinAT",17);
 params.SetSlotName("TtimAT",18);

 params.SetSignal(fAdcminAM,1);
 params.SetSignal(fAdcmaxAM,2);
 params.SetSignal(fTotminAEM,3);
 params.SetSignal(fTotmaxAEM,4);
 params.SetSignal(fTotminAOM,5);
 params.SetSignal(fTotmaxAOM,6);
 params.SetSignal(fRmaxA,7);
 params.SetSignal(fDtmaxA,8);
 params.SetSignal(fTwinAM,9);
 params.SetSignal(fTtimAM,10);
 params.SetSignal(fAdcminAT,11);
 params.SetSignal(fAdcmaxAT,12);
 params.SetSignal(fTotminAET,13);
 params.SetSignal(fTotmaxAET,14);
 params.SetSignal(fTotminAOT,15);
 params.SetSignal(fTotmaxAOT,16);
 params.SetSignal(fTwinAT,17);
 params.SetSignal(fTtimAT,18);

 fEvt->AddDevice(params);

 Amanda();
 InIce();
 IceTop();
}
///////////////////////////////////////////////////////////////////////////
void IceCleanHits::Amanda()
{
// Hit cleaning for Amanda modules.

 AliDevice* daq=(AliDevice*)fEvt->GetDevice("Daq");
 if (!daq) return;

 if (daq->GetSignal("Muon")) MuDaq();
 if (daq->GetSignal("TWR")) TWRDaq();
}
///////////////////////////////////////////////////////////////////////////
void IceCleanHits::MuDaq()
{
// Hit cleaning for Amanda MuDaq data.

 // Trigger time setting according to the user selection.
 // By default the "main" trigger of the event will be used.
 if (fTnamAM != "user")
 {
  AliDevice* tdev=(AliDevice*)fEvt->GetDevice("Trigger");
  if (tdev)
  {
   AliSignal* trig=tdev->GetHit(fTnamAM);
   if (trig) fTtimAM=trig->GetSignal("trig_pulse_le");
  }
 }

 // All Amanda OMs with a signal
 TObjArray* aoms=fEvt->GetDevices("IceAOM");
 if (!aoms) return;

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
 Int_t readout;
 for (Int_t iom=0; iom<oms.GetEntries(); iom++)
 {
  omx=(IceAOM*)oms.At(iom);
  if (!omx) continue;
  omid=omx->GetUniqueID();
  readout=int(omx->GetSignal("READOUT"));
  // General readout setting in case info was missing 
  if (!readout)
  {
   readout=1;
   if (omid>=303) readout=2; // Optical OMs
  }
  for (Int_t ih=1; ih<=omx->GetNhits(); ih++)
  {
   sx=omx->GetHit(ih);
   if (!sx) continue;
   adc=sx->GetSignal("ADC",7);
   le=sx->GetSignal("LE",7);
   tot=sx->GetSignal("TOT",7);

   // Remove hits with an ADC value outside the range
   if (adc<fAdcminAM || adc>fAdcmaxAM)
   {
    sx->SetDead("ADC");
    clean=0;
   }
   // Remove hits with a TOT value outside the range
   // Note : Different ranges for electrical and optical modules
   if (readout==1) // Electrical OMs
   {
    if (tot<fTotminAEM || tot>fTotmaxAEM)
    {
     sx->SetDead("TOT");
     clean=0;
    }
   }
   else // Optical OMs
   {
    if (tot<fTotminAOM || tot>fTotmaxAOM)
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
   if (fabs(le-fTtimAM)>fTwinAM)
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
 // Only retain hits that have at least one hit of another OM within a certain
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
   sx2=(AliSignal*)hits.At(jh2);
   if (!sx2) continue;

   omx1=(IceAOM*)sx1->GetDevice();
   omx2=(IceAOM*)sx2->GetDevice();
   if (omx1==omx2) continue;

   t1=sx1->GetSignal("LE",7);
   t2=sx2->GetSignal("LE",7);
   dt=fabs(t2-t1);
   if (dt>fDtmaxA) continue;

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
void IceCleanHits::TWRDaq()
{
// Hit cleaning for Amanda TWRDaq data.

 // Trigger time setting according to the user selection.
 // By default the "main" trigger of the event will be used.
 if (fTnamAT != "user")
 {
  AliDevice* tdev=(AliDevice*)fEvt->GetDevice("Trigger");
  if (tdev)
  {
   AliSignal* trig=tdev->GetHit(fTnamAT);
   if (trig) fTtimAT=trig->GetSignal("trig_pulse_le");
  }
 }

 // All Amanda OMs with a signal
 TObjArray* aoms=fEvt->GetDevices("IceAOM");
 if (!aoms) return;

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
 Int_t readout;
 for (Int_t iom=0; iom<oms.GetEntries(); iom++)
 {
  omx=(IceAOM*)oms.At(iom);
  if (!omx) continue;
  omid=omx->GetUniqueID();
  readout=int(omx->GetSignal("READOUT"));
  // General readout setting in case info was missing 
  if (!readout)
  {
   readout=1;
   if (omid>=303) readout=2; // Optical OMs
  }
  for (Int_t ih=1; ih<=omx->GetNhits(); ih++)
  {
   sx=omx->GetHit(ih);
   if (!sx) continue;
   adc=sx->GetSignal("ADC",7);
   le=sx->GetSignal("LE",7);
   tot=sx->GetSignal("TOT",7);

   // Remove hits with an ADC value outside the range
   if (adc<fAdcminAT || adc>fAdcmaxAT)
   {
    sx->SetDead("ADC");
    clean=0;
   }
   // Remove hits with a TOT value outside the range
   // Note : Different ranges for electrical and optical modules
   if (readout==1) // Electrical OMs
   {
    if (tot<fTotminAET || tot>fTotmaxAET)
    {
     sx->SetDead("TOT");
     clean=0;
    }
   }
   else // Optical OMs
   {
    if (tot<fTotminAOT || tot>fTotmaxAOT)
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
   if (fabs(le-fTtimAT)>fTwinAT)
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
 // Only retain hits that have at least one hit of another OM within a certain
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
   sx2=(AliSignal*)hits.At(jh2);
   if (!sx2) continue;

   omx1=(IceAOM*)sx1->GetDevice();
   omx2=(IceAOM*)sx2->GetDevice();
   if (omx1==omx2) continue;

   t1=sx1->GetSignal("LE",7);
   t2=sx2->GetSignal("LE",7);
   dt=fabs(t2-t1);
   if (dt>fDtmaxA) continue;

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
