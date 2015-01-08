//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 1998      Caltech, UCSB
//
// Module: EvtGen/EvtGen.hh
//
// Description:Main class to provide user interface to EvtGen.
//
// Modification history:
//
//    RYD     March 24, 1998     Module created
//
//    DJL     August 10, 1998    Additional Event member function added
//
//    RYD     December 25, 1999  Any application using EvtGen will need
//                               to instantiate an instance of this class
//                               and hold on to it untill done generating
//                               events. This class will now hold data used
//                               for the lifetime of the generator.
//    Lange   April 18, 2002 - split "status" info into own class

//------------------------------------------------------------------------

#ifndef EVTSTATUS_HH
#define EVTSTATUS_HH


class EvtStatus{

public:


  static void setRejectFlag() {int *temp=rejectFlag();  *temp=1; return;}
  static void initRejectFlag() {int *temp=rejectFlag();  *temp=0; return;}
  static int* rejectFlag() {static int rejectEvent=0; return &rejectEvent;}
  static int getRejectFlag() {int *temp=rejectFlag(); return *temp;}
  
};



#endif

