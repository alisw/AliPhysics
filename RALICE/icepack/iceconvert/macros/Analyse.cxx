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
// Class Analyse
// TTask derived class to demonstrate the concept of analysis via Tasks.
// This example task investigates an OM database structure and is invoked
// in the example macro IceCal2Root.cc.
//
//--- Author: Nick van Eijndhoven 08-aug-2005 Utrecht University
//- Modified: NvE $Date$ Utrecht University
///////////////////////////////////////////////////////////////////////////
 
#include "Analyse.h"
#include "Riostream.h"

ClassImp(Analyse) // Class implementation to enable ROOT I/O

Analyse::Analyse(const char* name,const char* title) : TTask(name,title)
{
// Default constructor.
}
///////////////////////////////////////////////////////////////////////////
Analyse::~Analyse()
{
// Default destructor.
}
///////////////////////////////////////////////////////////////////////////
void Analyse::Exec(Option_t* opt)
{
// Implementation of the event analysis to be performed.

 TString name=opt;
 AliJob* parent=(AliJob*)(gROOT->GetListOfTasks()->FindObject(name.Data()));

 if (!parent) return;

 AliObjMatrix* omdb=(AliObjMatrix*)parent->GetMainObject();
 if (!omdb) return;

 IceAOM* omx=0;
 for (Int_t i=1; i<=10; i++)
 {
  omx=(IceAOM*)omdb->GetObject(i,1);
  if (omx) omx->Data();
 }

 TF1* fxtalk=0;
 fxtalk=(TF1*)omdb->GetObject(90,113+1);
 if (fxtalk) fxtalk->Print();
 fxtalk=(TF1*)omdb->GetObject(121,102+1);
 if (fxtalk) fxtalk->Print();
}
///////////////////////////////////////////////////////////////////////////
