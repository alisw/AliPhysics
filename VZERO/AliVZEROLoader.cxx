#include "AliVZEROLoader.h"

const TString AliVZEROLoader::fgkDefaultHitsFileName      = "VZERO.Hits.root";
const TString AliVZEROLoader::fgkDefaultDigitsFileName    = "VZERO.Digits.root";

ClassImp(AliVZEROLoader)
AliVZEROLoader::AliVZEROLoader()
 {
 }
/*****************************************************************************/ 
AliVZEROLoader::AliVZEROLoader(const Char_t *name,const Char_t *topfoldername)
 :AliLoader(name,topfoldername)
{
  if (GetDebug()) 
   Info("AliVZEROLoader"," name = %s; topfolder = %s",name,topfoldername);
}
/*****************************************************************************/ 

AliVZEROLoader::AliVZEROLoader(const Char_t *name,TFolder *topfolder)
 :AliLoader(name,topfolder)
 {
 }

