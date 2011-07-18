#ifndef ALITPCMONITORDATEMONITOR_H
#define ALITPCMONITORDATEMONITOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


////////////////////////////////////////////////////////////////////////
////
//// AliTPCMonitorDateMonitor class
//// 
//// Monitoring wrapper class for DATE monitoring
//// 
//// Authors: Roland Bramm, 
////          Stefan Kniege, IKF, Frankfurt
////       
/////////////////////////////////////////////////////////////////////////

#include <string>
#include "TNamed.h"

using namespace std;

class AliTPCMonitorDateMonitor : public TNamed {
 
 public:
    AliTPCMonitorDateMonitor();
    AliTPCMonitorDateMonitor(const  AliTPCMonitorDateMonitor &config);
    AliTPCMonitorDateMonitor& operator= (const AliTPCMonitorDateMonitor& config);
    ~AliTPCMonitorDateMonitor();
    
    void    Free();
    Int_t   OpenMonitoring(string name);
    const Char_t* DecodeError(int error);
    Int_t   DeclareMonitor(string name); 
    Int_t   FlushEvents();
    Int_t   GetEvent();
    Char_t *GetEventPointerasChar();
    Int_t   Logout();
	
 private:
	
	Char_t *fPointer;                     // pointer to event
	
	ClassDef(AliTPCMonitorDateMonitor,1);
};
#endif
