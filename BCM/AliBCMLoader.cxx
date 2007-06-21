/////////////////////////////////////////////////////////////////////
//                                                                 //
//  Class AliBCMLoader                                             //
//                                                                 //
//  Base class for BCMLoaders.                                     //
//  Loader provides base I/O facilities for standard data.         //
//  Each detector has a loader data member.                        //
//  Loader is accessible via folder structure as well.             //
//                                                                 //
/////////////////////////////////////////////////////////////////////

#include "AliBCMLoader.h"
#include "AliLog.h"

const TString AliBCMLoader::fgkDefaultHitsFileName      = "BCM.Hits.root";
const TString AliBCMLoader::fgkDefaultDigitsFileName    = "BCM.Digits.root";

ClassImp(AliBCMLoader)

//_____________________________________________________________________________
AliBCMLoader::AliBCMLoader()
 { 
 // Default constructor
 }

//_____________________________________________________________________________
AliBCMLoader::AliBCMLoader(const Char_t *name,const Char_t *topfoldername)
 :AliLoader(name,topfoldername)
{
  AliDebug(1,Form("Name = %s; topfolder = %s",name,topfoldername));
}

//_____________________________________________________________________________
AliBCMLoader::AliBCMLoader(const Char_t *name,TFolder *topfolder)
 :AliLoader(name,topfolder)
 {
 }
