// -*- C++ -*-
#ifndef _ALI_VDM_DATA_H
#define _ALI_VDM_DATA_H

#include <TSystem.h>

#include "AliXMLEngine.h"
#include "AliTriggerBCMask.h"

class AliVdMData {
public:
  AliVdMData(const char* xmlFileName)
    : fXML(xmlFileName)
    , fNodeVdM(fXML.GetRootNode())
    , fFillNumber(std::stoi(fNodeVdM.GetAttr("fill").GetData()))
    , fTriggerBCMask(Form("BCMask fill %d", fFillNumber),
                     fNodeVdM.GetAttr("bcMask").GetData())
    , fNodeScans(fNodeVdM.GetChild("Scans")) {
    Printf("fill %d  bcs=%d", fFillNumber, fTriggerBCMask.GetNUnmaskedBCs());
  }

  AliXMLEngine::Node GetScansBegin() const { return fNodeScans.GetChildBegin(); }
  AliXMLEngine::Node GetScansEnd()   const { return fNodeScans.GetChildEnd(); }

  Int_t GetFillNumber() const { return fFillNumber; }
  const AliTriggerBCMask& GetTriggerBCMask() const { return fTriggerBCMask; }

  static TString GetScanName(const AliXMLEngine::Node& n) {
    const TString scanType(n.GetAttr("type").GetData());
    TString scanName("Scan");
    if (scanType == "X" || scanType == "Y")
      scanName += n.GetAttr("number").GetData();
    scanName += scanType;
    return scanName;
  }
  static int GetScanType(const AliXMLEngine::Node& n) {
    return (TString(n.GetAttr("type").GetData()).Contains("X") ? 0 : 1);
  }
  static TString GetFileName(const TString& path) {
    const TString localPath(gSystem->ExpandPathName("$HOME/cernbox/data/vdM/"+path));
    if (!gSystem->AccessPathName(localPath.Data()))
      return localPath;
    const TString eosBasePath("root://eosuser.cern.ch//eos/user/c/cmayer/data/vdM/");
    return eosBasePath+path;
  }
protected:
private:
  AliXMLEngine             fXML;           //!
  const AliXMLEngine::Node fNodeVdM;       //!
  const Int_t              fFillNumber;    //!
  const AliTriggerBCMask   fTriggerBCMask; //!
  const AliXMLEngine::Node fNodeScans;     //!
} ;

#endif // _ALI_VDM_DATA_H
