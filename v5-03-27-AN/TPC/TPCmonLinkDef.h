#ifdef __CINT__
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Log$
/* Revision 1.2  2007/09/17 16:34:54  cvetan
/* The package was overwriting the rootcint flags. This was fixed by applying the necessary changes in the DATE-dependent parts of the code
/*
/* Revision 1.1  2007/09/17 10:23:31  cvetan
/* New TPC monitoring package from Stefan Kniege. The monitoring package can be started by running TPCMonitor.C macro located in macros folder.
/* */


#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class  AliTPCMonitorConfig+; 
#pragma link C++ class  AliTPCMonitorAltro+; 
#pragma link C++ class  AliTPCMonitorFFT+; 
#pragma link C++ class  AliTPCMonitorMappingHandler+;

#ifdef ALI_DATE
#pragma link C++ class  AliTPCMonitorDateMonitor+;
#endif
#pragma link C++ class  AliTPCMonitorDateFile+;  
#pragma link C++ class  AliTPCMonitorDateFormat+;  
#pragma link C++ class  AliTPCMonitorDialog+;  
#pragma link C++ class  AliTPCMonitor+; 

#endif

