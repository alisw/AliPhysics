#ifndef __CINT__
#include <TSystemDirectory.h>
#include <TFile.h>
#include <TList.h>
#include <AliCDBEntry.h>
#include <AliFMDMap.h>
#include <AliFMDCalibPedestal.h>
#include <AliFMDCalibGain.h>
#include <AliFMDCalibStripRange.h>
#include <AliFMDCalibSampleRate.h>
#include <TString.h>
#include <TSystem.h>
#include <TError.h>
#else
class AliFMDMap;
#endif

Bool_t 
CheckMap(const char* path, const AliFMDMap* map)
{
  if (!map) { 
    Warning("CheckFile", "No map in %s", path);
    return false;
  }
  if (!map->Ptr() || map->MaxIndex() <= 0) { 
    Warning("CheckFile", "Map %p (%d) has no data in %s", 
	    map->Ptr(), map->MaxIndex(), path);
    return false;
  }
  return true;
}

enum {
  kMap, 
  kPedestal, 
  kGain, 
  kRate, 
  kRange
};
  
Bool_t
CheckFile(const char* name, const char* dirName, Int_t which)
{
  TString path(gSystem->ConcatFileName(dirName, name));
  TFile* file = TFile::Open(path, "READ");
  if (!file) { 
    Warning("CheckFile", "Failed to open %s", path.Data());
    return false;
  }
  AliCDBEntry* entry = static_cast<AliCDBEntry*>(file->Get("AliCDBEntry"));
  if (!entry) { 
    Warning("CheckFile", "No entry in %s", path.Data());
    file->Close();
    return false;
  }
  TObject* object = entry->GetObject();
  if (!object) { 
    Warning("CheckFile", "Entry has no object in %s", path.Data());
    file->Close();
    return false;
  }

  const AliFMDMap* map = 0;
  if (which == kMap) map = static_cast<AliFMDMap*>(object);
  else if (which == kPedestal) 
    map = &(static_cast<AliFMDCalibPedestal*>(object)->Values());
  else if (which == kGain) 
    map = &(static_cast<AliFMDCalibGain*>(object)->Values());
  else if (which == kRate) 
    map = &(static_cast<AliFMDCalibSampleRate*>(object)->Rates());
  else if (which == kRange) 
    map = &(static_cast<AliFMDCalibStripRange*>(object)->Ranges()); 
  else {
    Warning("CheckFile", "Don't now how to deal with what=%d", which);
    file->Close();
    return false;
  }
  if (!CheckMap(path.Data(), map)) { 
    file->Close();
    return false;
  }
  Info("CheckFile", "Map OK in %s", path.Data());
  file->Close();
  return true;
}

    
void
CheckCalibData(const char* dirName)
{
  TString dirS(dirName);
  if (dirS.EndsWith("/")) dirS.Remove(dirS.Length()-1);
  dirS = gSystem->BaseName(dirS.Data());
  Int_t what = 0;
  if (dirS == "Dead" || dirS == "ZeroSuppression") what = kMap;
  else if (dirS == "Pedestal") what = kPedestal;
  else if (dirS == "PulseGain") what = kGain;
  else if (dirS == "SampleRate") what = kRate;
  else if (dirS == "StripRange") what = kRange;
  else {
    Error("CheckCalibData", "Don't know how to deal with %s in %s", 
	  dirS.Data(), dirName);
    return;
  }
    
  TSystemDirectory dir(dirName, dirName);
  TList* files(dir.GetListOfFiles());
  TIter next(files);
  TObject* obj = 0;
  
  Int_t nTotal = 0;
  Int_t nOk    = 0;
  while ((obj = next())) { 
    TString name(obj->GetName());
    if (!name.EndsWith(".root")) continue;
    nTotal++;
    if (CheckFile(name, dirName, what)) nOk++;
  }
  Info("CheckCalibData", "Total: %d, OK: %d, Bad: %d in %s ", 
       nTotal, nOk, nTotal - nOk, dirName);
}
