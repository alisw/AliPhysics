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
*/ 

#include <iostream>
#include "AliTPCMonitorDateMonitor.h"
ClassImp(AliTPCMonitorDateMonitor)

//_____________________________________________________________________________
AliTPCMonitorDateMonitor::AliTPCMonitorDateMonitor(){
  // Constructor
  fPointer = 0;
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
char* AliTPCMonitorDateMonitor::DecodeError(Int_t error  ){
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
