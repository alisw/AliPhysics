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
 * @param base   Base name 
 * @param dir    Directory 
 * @param unpack If true, also unzip the retrieved archive 
 * 
 * @return true on success 
 */
Bool_t GetOne(const TString& base, const TString& dir, Bool_t unpack)
{
  TString src(gSystem->ConcatFileName(base,dir));
  src = gSystem->ConcatFileName(src,"root_archive.zip");
  TString name;
  name.Form("root_archive_%s",dir.Data());
  TString dest;
  dest.Form("%s.zip",name.Data());

  if (!TFile::Cp(src, dest)) {
    Error("GetOne","Failed to download %s -> %s",
          src.Data(), dest.Data());
    return false;
  }
  if (!unpack) return true;
  gSystem->Exec(Form("mkdir -p %s && (cd %s && unzip -n ../%s)", 
		     name.Data(), name.Data(), dest.Data()));
  return true;
}

void GridDownload(const TString& base, const TString& runs, Bool_t unpack)
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
    GetOne(base, run->String(), unpack);
  }
}
// EOF

