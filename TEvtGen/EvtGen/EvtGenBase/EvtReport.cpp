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
// Module: EvtReport.cc
//
// Description: definitions of global functions.
//
// Modification history:
//
//    Simon Patton   June 3, 1996           Module created
//
//------------------------------------------------------------------------
//
#include "EvtGenBase/EvtPatches.hh"

#include "EvtGenBase/EvtReport.hh"
using std::cerr;
using std::cout;
using std::endl;
using std::ostream;


//
// constants, enums and typedefs
//


ostream& report( Severity::Enum severity ,
                 const char* facility )
{
   int printNoFacility=1;

   if ( ( facility == 0 ) &&
        ( printNoFacility ==1) ) {
      cout << "There is no `facility' implemented in `report'"
                        << endl ;
      printNoFacility = 0 ;
   }
   if ( severity < Severity::Warning ) {
     if (facility[0]!=0){
       cerr<<facility<<":";
     }
     return ( cerr ) ;
   }
   if (facility[0]!=0){
     cout<<facility<<":";
   }    
   return cout;
}


