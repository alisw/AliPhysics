/////////////////////////////////////////////////////////////////////
//                                                                 //
//  Class AliADLoader                                              //
//                                                                 //
//  Base class for ADLoaders.                                      //
//  Loader provides base I/O facilities for standard data.         //
//  Each detector has a loader data member.                        //
//  Loader is accessible via folder structure as well.             //
//                                                                 //
/////////////////////////////////////////////////////////////////////

#include "AliADLoader.h"
#include "AliLog.h"

const TString AliADLoader::fgkDefaultHitsFileName      = "AD.Hits.root";
const TString AliADLoader::fgkDefaultDigitsFileName    = "AD.Digits.root";

ClassImp(AliADLoader);

//_____________________________________________________________________________
AliADLoader::AliADLoader()
{
}

//_____________________________________________________________________________
AliADLoader::AliADLoader(const Char_t *name,const Char_t *topfoldername)
  : AliLoader(name,topfoldername)
{
  AliDebugF(1, "Name = %s; topfolder = %s", name, topfoldername);
}

//_____________________________________________________________________________
AliADLoader::AliADLoader(const Char_t *name,TFolder *topfolder)
  : AliLoader(name,topfolder)
{
}
