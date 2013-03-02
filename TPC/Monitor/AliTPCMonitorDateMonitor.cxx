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

/*
$Log$
Revision 1.2  2007/09/17 16:34:54  cvetan
The package was overwriting the rootcint flags. This was fixed by applying the necessary changes in the DATE-dependent parts of the code

Revision 1.1  2007/09/17 10:23:31  cvetan
New TPC monitoring package from Stefan Kniege. The monitoring package can be started by running TPCMonitor.C macro located in macros folder.

*/ 

////////////////////////////////////////////////////////////////////////
//
// AliTPCMonitorDateMonitor class
// 
// Monitoring wrapper class for DATE raw data monitoring used by the TPC 
// raw data monitor.
// Online monitoring is only possible if DATE is installed on the machine.
// If not, raw data can be read using the AliTPCMonitorDateFile 
// (can be choosen from 'Sel.Format' in the TPCMonitor gui")
//
// Authors: Roland Bramm, 
//          Stefan Kniege, IKF, Frankfurt
//       
/////////////////////////////////////////////////////////////////////////



#include <iostream>
#include "AliTPCMonitorDateMonitor.h"
#include "event.h"
#include "monitor.h"
#include <cstdlib>

ClassImp(AliTPCMonitorDateMonitor)

//_____________________________________________________________________________
AliTPCMonitorDateMonitor::AliTPCMonitorDateMonitor():
fPointer(0)
{
  // Constructor
}


//_____________________________________________________________________________
AliTPCMonitorDateMonitor::AliTPCMonitorDateMonitor(const AliTPCMonitorDateMonitor &datemon) :
  TNamed(datemon.GetName(),datemon.GetTitle()),
  fPointer(datemon.fPointer)
{
  // copy constructor  
}

//_____________________________________________________________________________
AliTPCMonitorDateMonitor &AliTPCMonitorDateMonitor::operator =(const AliTPCMonitorDateMonitor& datemon)
{
  // assignement operator
  if(this!=&datemon){ 
    fPointer=datemon.fPointer;
  }

  return *this;
}

//_____________________________________________________________________________
AliTPCMonitorDateMonitor::~AliTPCMonitorDateMonitor(){
  // Destructor
}

//_____________________________________________________________________________
Int_t AliTPCMonitorDateMonitor::OpenMonitoring(string name){
  // Set data source for the monitor
	return monitorSetDataSource((char*)name.c_str());
}

//_____________________________________________________________________________
const char* AliTPCMonitorDateMonitor::DecodeError(Int_t error  ){
  // Return decoded error string 
  return monitorDecodeError(error);
}

//_____________________________________________________________________________
Int_t AliTPCMonitorDateMonitor::DeclareMonitor(string name){
  // Declare monitorn Map
  return monitorDeclareMp((char*)name.c_str());
}

//_____________________________________________________________________________
Int_t AliTPCMonitorDateMonitor::GetEvent(){
  // Get next event
  return monitorGetEventDynamic( (void**)&fPointer );
}
 
//_____________________________________________________________________________
char *AliTPCMonitorDateMonitor::GetEventPointerasChar(){
  // Return pointer to data mamory
  return fPointer;
}

//_____________________________________________________________________________
Int_t AliTPCMonitorDateMonitor::Logout(){
  // Logout the monitor
  return monitorLogout();
}
//_____________________________________________________________________________
int AliTPCMonitorDateMonitor::FlushEvents(){
  // Flush events 
  return monitorFlushEvents();
}

//_____________________________________________________________________________
void  AliTPCMonitorDateMonitor::Free(){
  // Free fPointer 
  if(fPointer)free(fPointer);
  fPointer = 0;
}
