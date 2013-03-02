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

//  re - write for upgraded TPC --- not before

#pragma link C++ class  AliTPCMonitorConfig+;         // Config Conatainer
#pragma link C++ class  AliTPCMonitorAltro+;          // ALtro Mapping ...  to be checked if ituseds Base mapping
#pragma link C++ class  AliTPCMonitorFFT+;            // FFT of pad signals
#pragma link C++ class  AliTPCMonitorMappingHandler+; // TPC mapper ... duplication ... but no time to be spend

#ifdef ALI_DATE
#pragma link C++ class  AliTPCMonitorDateMonitor+;    // Used to read "Date Format" --- check if it is still used   
#endif
#pragma link C++ class  AliTPCMonitorDateFile+;       // Used to read "Date Format" --- check if it is still used   
#pragma link C++ class  AliTPCMonitorDateFormat+;     // Used to read "Date Format" --- check if it is still used   
#pragma link C++ class  AliTPCMonitorDialog+;         // dialog for opening data streams
#pragma link C++ class  AliTPCMonitor+;               // Main monitor class

#endif

