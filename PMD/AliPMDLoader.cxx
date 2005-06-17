#include "AliPMDLoader.h"
#include "AliLog.h"

const TString AliPMDLoader::fgkDefaultHitsFileName      = "PMD.Hits.root";
const TString AliPMDLoader::fgkDefaultSDigitsFileName   = "PMD.SDigits.root";
const TString AliPMDLoader::fgkDefaultDigitsFileName    = "PMD.Digits.root";
const TString AliPMDLoader::fgkDefaultRecPointsFileName = "PMD.RecPoints.root";
const TString AliPMDLoader::fgkDefaultTracksFileName    = "PMD.Tracks.root";


ClassImp(AliPMDLoader)
AliPMDLoader::AliPMDLoader()
 {
 }
/*****************************************************************************/ 
AliPMDLoader::AliPMDLoader(const Char_t *name,const Char_t *topfoldername)
 :AliLoader(name,topfoldername)
{
  AliInfoClass(Form(" name = %s; topfolder = %s",name,topfoldername));
}
/*****************************************************************************/ 

AliPMDLoader::AliPMDLoader(const Char_t *name,TFolder *topfolder)
 :AliLoader(name,topfolder)
 {
 }

