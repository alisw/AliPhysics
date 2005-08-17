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
const TString kgESDTreeName     = "esdTree" ; 
const UInt_t  kgeventsToRead    = 1 ;
const char *  kgPattern         = "Evt" ;

TChain * AliReadESDfromdisk(const UInt_t eventsToRead, 
			       const TString dirName, 
			       const TString esdTreeName = kgESDTreeName, 
			       const char *  pattern     = kgPattern) ; 
TChain * AliReadESD(const UInt_t eventsToRead = kgeventsToRead,
		     TString fileName  = "", 
		     const TString esdTreeName = kgESDTreeName, 
		     const char *  pattern     = kgPattern ) ;  

// Needed for AliLog (return the macro name)
char * ClassName() { return "macro" ; } 
   
