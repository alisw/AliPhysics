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
// Class EvtAna
// TTask derived class to demonstrate the concept of analysis via Tasks.
// This example task investigates an IceEvent structure and is invoked
// in the example macros icextalk.cc and icecalib.cc.
//
//--- Author: Nick van Eijndhoven 07-may-2005 Utrecht University
//- Modified: NvE $Date$ Utrecht University
///////////////////////////////////////////////////////////////////////////
 
#include "EvtAna.h"
#include "Riostream.h"

ClassImp(EvtAna) // Class implementation to enable ROOT I/O

EvtAna::EvtAna(const char* name,const char* title) : TTask(name,title)
{
// Default constructor.
}
///////////////////////////////////////////////////////////////////////////
EvtAna::~EvtAna()
{
// Default destructor.
}
///////////////////////////////////////////////////////////////////////////
void EvtAna::Exec(Option_t* opt)
{
// Implementation of the event analysis to be performed.

 TString name=opt;
 AliJob* parent=(AliJob*)(gROOT->GetListOfTasks()->FindObject(name.Data()));

 if (!parent) return;

 IceEvent* evt=(IceEvent*)parent->GetObject("IceEvent");
 if (!evt) return;
 
 evt->Data();
 evt->ShowDevices();

  IceAOM* om=(IceAOM*)evt->GetIdDevice(125);
  if (om)
  {
   om->Data();
   AliSignal* sx=om->GetHit(1);
   if (sx)
   {
    cout << " Calibrated ADC : " << sx->GetSignal("ADC",7) << endl;
    cout << " De-calibrated ADC : " << sx->GetSignal("ADC",-7) << endl;
    cout << " Calibrated LE : " << sx->GetSignal("LE",7) << endl;
    cout << " De-calibrated LE : " << sx->GetSignal("LE",-7) << endl;
    cout << " Calibrated TOT : " << sx->GetSignal("TOT",7) << endl;
    cout << " De-calibrated TOT : " << sx->GetSignal("TOT",-7) << endl;
   }
  }
}
///////////////////////////////////////////////////////////////////////////
