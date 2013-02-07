/**
 * @file   GridDownload.C
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Wed Jan 23 21:51:40 2013
 * 
 * @brief  Helper script to download results from the Grid 
 * 
 * 
 * @ingroup pwglf_forward_trains_helper
 */
#ifndef __CINT__
# include <TString.h>
# include <TSystem.h>
# include <TGrid.h>
# include <TFile.h>
# include <TObjArray.h>
# include <TObjString.h>
# include <TError.h>
# include <TEnv.h>
#else
class TString;
#endif

/** 
 * Get one file 
 * 
 * @param base Base name 
 * @param dir  Directory 
 * 
 * @return true on success 
 */
Bool_t GetOne(const TString& base, const TString& dir)
{
  TString src(gSystem->ConcatFileName(base,dir));
  src = gSystem->ConcatFileName(src,"root_archive.zip");
  TString dest;
  dest.Form("root_archive_%s.zip",dir.Data());

  if (!TFile::Cp(src, dest)) {
    Error("GetOne","Failed to download %s -> %s",
          src.Data(), dest.Data());
    return false;
  }
  return true;
}

void GridDownload(const TString& base, const TString& runs)
{
  gEnv->SetValue("XSec.GSI.DelegProxy", "2");
  if (!TGrid::Connect("alien://")) {
    Error("Download","Failed to connect to AliEn");
    return;
  }

  TObjArray*  runArray = runs.Tokenize(" ");
  TObjString* run      = 0;
  TIter       next(runArray);
  while ((run = static_cast<TObjString*>(next()))) {
    GetOne(base, run->String());
  }
}
// EOF

