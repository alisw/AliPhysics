/////////////////////////////////////////////////////////////////////
//                                                                 //
//  Class AliVZEROLoader                                           //
//                                                                 //
//  Base class for VZEROLoaders.                                   //
//  Loader provides base I/O facilities for standard data.         //
//  Each detector has a loader data member.                        //
//  Loader is accessible via folder structure as well.             //
//                                                                 //
/////////////////////////////////////////////////////////////////////

#include "AliVZEROLoader.h"

const TString AliVZEROLoader::fgkDefaultHitsFileName      = "VZERO.Hits.root";
const TString AliVZEROLoader::fgkDefaultDigitsFileName    = "VZERO.Digits.root";

ClassImp(AliVZEROLoader)

//_____________________________________________________________________________
AliVZEROLoader::AliVZEROLoader()
 { 
 // Default constructor
 }

//_____________________________________________________________________________
AliVZEROLoader::AliVZEROLoader(const Char_t *name,const Char_t *topfoldername)
 :AliLoader(name,topfoldername)
{
  if (GetDebug()) 
   Info("AliVZEROLoader"," name = %s; topfolder = %s",name,topfoldername);
}

//_____________________________________________________________________________
AliVZEROLoader::AliVZEROLoader(const Char_t *name,TFolder *topfolder)
 :AliLoader(name,topfolder)
 {
 }

//_______________________________________________________________________________
AliVZEROLoader::~AliVZEROLoader()
{
 // Destructor 
}
