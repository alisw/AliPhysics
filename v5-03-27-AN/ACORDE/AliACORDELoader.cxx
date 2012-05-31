/////////////////////////////////////////////////////////////////////
//                                                                 //
//  Class AliACORDELoader                                           //
//                                                                 //
//  Base class for ACORDELoaders.                                   //
//  Loader provides base I/O facilities for standard data.         //
//  Each detector has a loader data member.                        //
//  Loader is accessible via folder structure as well.             //
//                                                                 //
/////////////////////////////////////////////////////////////////////

#include "AliACORDELoader.h"
#include "AliLog.h"

const TString AliACORDELoader::fgkDefaultHitsFileName= "ACORDE.Hits.root";
const TString AliACORDELoader::fgkDefaultDigitsFileName= "ACORDE.Digits.root";

ClassImp(AliACORDELoader)

//_____________________________________________________________________________
AliACORDELoader::AliACORDELoader()
 { 
 // Default constructor
 }

//_____________________________________________________________________________
AliACORDELoader::AliACORDELoader(const Char_t *name,const Char_t *topfoldername)
 :AliLoader(name,topfoldername)
{
  AliDebug(1,Form("Name = %s; topfolder = %s",name,topfoldername));
}

//_____________________________________________________________________________
AliACORDELoader::AliACORDELoader(const Char_t *name,TFolder *topfolder)
 :AliLoader(name,topfolder)
 {
 }
