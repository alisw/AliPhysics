#if !defined(__CINT__) || defined(__MAKECINT__)

// Root include files
#include <Riostream.h>
#include <TString.h>
#include <TChain.h>
#include <TSystem.h>

// AliROOT include files
#include <AliLog.h>

#else
#endif

// Define global parameters 
// const TString kgESDTreeName     = "esdTree" ;
// const UInt_t  kgeventsToRead    = 1 ;
// const char *  kgPattern         = "Evt" ;

TChain * AliReadESDfromdisk(const UInt_t eventsToRead,
			       const TString dirName, 
			       const TString esdTreeName = "esdTree",
			       const char *  pattern     = "Evt") ;
TChain * AliReadESD(const UInt_t eventsToRead = 1,
		     TString fileName  = "",
		     const TString esdTreeName = "esdTree",
		     const char *  pattern     = "Evt" ) ;

// Needed for AliLog (return the macro name)
char * ClassName() { return "macro" ; } 
   
