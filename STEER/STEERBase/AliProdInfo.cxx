/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


//-----------------------------------------------------------------
//        Base class for parsing alirootVersion                  //
//        TNamed object stored in ESD/AOD data                   //
//        and setup easy getters for end user                    //
//                                                               //
//                                                               //
//   Origin: Pietro Antonioli, INFN-BO Pietro.Antonioli@cern.ch  //
//                                                               //
//-----------------------------------------------------------------

#include <AliLog.h>
#include "AliProdInfo.h"

ClassImp(AliProdInfo);

AliProdInfo::AliProdInfo()
  :fMcFlag(kFALSE)
  ,fAlirootSvnVersion(0)
  ,fRootSvnVersion(0)
  ,fRecoPass(-1)
  ,fPeriod("")
  ,fProductionTag("")
  ,fAlirootVersion("")
  ,fRootVersion("")
{	
  //
  // default constructor
  //	
  AliLog::SetClassDebugLevel("AliProdInfo",10);
}

AliProdInfo::AliProdInfo(const TString& name,const TString& title)
  :TNamed(name,title)
  ,fMcFlag(kFALSE)
  ,fAlirootSvnVersion(0)
  ,fRootSvnVersion(0)
  ,fRecoPass(-1)
  ,fPeriod("")
  ,fProductionTag("")
  ,fAlirootVersion("")
  ,fRootVersion("")
{	
  //
  // default constructor
  //	
  AliLog::SetClassDebugLevel("AliProdInfo",10);
}

AliProdInfo::AliProdInfo(TList *userInfoList)
  :fMcFlag(kFALSE)
  ,fAlirootSvnVersion(0)
  ,fRootSvnVersion(0)
  ,fRecoPass(-1)
  ,fPeriod("")
  ,fProductionTag("")
  ,fAlirootVersion("")
  ,fRootVersion("")
{	
  //
  // default constructor & init
  //	
  AliLog::SetClassDebugLevel("AliProdInfo",10);
  Init(userInfoList);
}

//-------------------------------------------------------------------------------------------------	
AliProdInfo::~AliProdInfo() {

}

//-------------------------------------------------------------------------------------------------	
void AliProdInfo::Init(TList *userInfoList) 
{
  fPeriod="";
  fAlirootVersion="";
  fAlirootSvnVersion=0;
  fRootVersion="";
  fRootSvnVersion=0;
  fMcFlag=kFALSE;
  fRecoPass=-1;
  TNamed *prodInfo = (TNamed *)userInfoList->FindObject("alirootVersion");
  if (!prodInfo) {
    AliError("No alirootVersion named object found in user info");
    return;
  }
  ParseProdInfo(prodInfo);
 }

//-------------------------------------------------------------------------------------------------	
void AliProdInfo::ParseProdInfo(TNamed *prodInfoData) 
{
  // parse information
  const char* key[] = {
    "aliroot"
    ,"root"
    ,"OutputDir="
    ,"LPMRawPass="
    ,"LPMProductionType="
    ,"LPMProductionTag="
    ,"LPMAnchorProduction="
  };
  enum {kAliroot,kRoot,kOutDir,kPass,kProdType,kProdTag,kPeriod};
  TString aliroot="";
  TString root="";
  TString outDir="";
  TString metadata="";
  TString lpmRawPass="";
  TString lpmProdType="";
  TString lpmProdTag="";
  TString lpmPeriod="";
  TString tmpStr="";
  //
  TString str(prodInfoData->GetTitle());
  TObjArray *tokens = str.Tokenize(";");  
  //
  for (Int_t i=0;i<=tokens->GetLast();i++) {
    TObjString *stObj = (TObjString *)tokens->At(i);
    tmpStr = stObj->GetString().Strip(TString::kBoth,' '); // strip irrelevant spaces
    //
    if (tmpStr.BeginsWith( key[kAliroot] ))
      aliroot    = tmpStr.ReplaceAll( key[kAliroot] ,"").Strip(TString::kBoth,' ');
    //
    else if (tmpStr.BeginsWith( key[kRoot] ))          
      root       = tmpStr.ReplaceAll( key[kRoot] ,"").Strip(TString::kBoth,' ');
    //
    else if (tmpStr.BeginsWith( key[kOutDir] ))
      outDir     = tmpStr.ReplaceAll( key[kOutDir] ,"").Strip(TString::kBoth,' ');
    //
    else if (tmpStr.BeginsWith( key[kPass] ))   
      lpmRawPass = tmpStr.ReplaceAll( key[kPass] ,"").Strip(TString::kBoth,' ');
    //
    else if (tmpStr.BeginsWith( key[kProdType] )) 
      lpmProdType = tmpStr.ReplaceAll( key[kProdType] ,"").Strip(TString::kBoth,' ');
    //
    else if (tmpStr.BeginsWith( key[kProdTag] )) 
      lpmProdTag = tmpStr.ReplaceAll( key[kProdTag] ,"").Strip(TString::kBoth,' ');
    //
    else if (tmpStr.BeginsWith( key[kPeriod] )) 
      lpmPeriod = tmpStr.ReplaceAll( key[kPeriod] ,"").Strip(TString::kBoth,' ');
    //
  }  
  delete tokens;
  // now interpret ...
  //
  // extract ALIROOT version
  if (!aliroot.IsNull()) {
    TObjArray *tali = (TObjArray *)aliroot.Tokenize(":");
    TObjString *tos = (TObjString *)tali->At(0);
    fAlirootVersion = "";
    if (tos) {
      fAlirootVersion = tos->GetString().Strip(TString::kBoth,' ');
      if (fAlirootVersion.IsNull()) AliWarning("Cannot extract AliROOT version string. Might be git related.");
    }
    //
    tos = (TObjString*)tali->At(1);
    if (tos){
      tmpStr = tos->GetString().Strip(TString::kBoth,' ');
      if (tmpStr.IsDigit()){
	fAlirootSvnVersion = tmpStr.Atoi();
      } 
      else if (tmpStr.IsHex()) {
	if (tmpStr.Length()>6) tmpStr.Resize(6);
	sscanf(tmpStr.Data(),"%x",&fAlirootSvnVersion);
	AliWarningF("ALIROOT SVN version number not decimal, might be on git. Reading as hex %s -> %d",tmpStr.Data(),fAlirootSvnVersion);
      }
      else {
	fAlirootSvnVersion=65263;
	AliWarningF("AliRoot SVN version is not extracted, setting to %d",fAlirootSvnVersion);
      }
    }
    delete tali;
  }
  else AliWarningF("Failed to extract %s version information",key[kAliroot]);
  //
  // extract ROOT version
  if (!root.IsNull()) { 
    TObjArray *tali = root.Tokenize(":");
    TObjString *tos = (TObjString *)tali->At(0);
    fRootVersion = "";
    if (tos) {
      fRootVersion = tos->GetString().Strip(TString::kBoth,' ');
      if (fRootVersion.IsNull()) AliWarning("Cannot extract ROOT version string. Might be git related.");
    }
    //
    tos = (TObjString*)tali->At(1);
    if (tos){
      tmpStr = tos->GetString().Strip(TString::kBoth,' ');
      if (tmpStr.IsDigit()){
	fRootSvnVersion = tmpStr.Atoi();
      } 
      else if (tmpStr.IsHex()) {
	if (tmpStr.Length()>6) tmpStr.Resize(6);
	sscanf(tmpStr.Data(),"%x",&fRootSvnVersion);
	AliWarningF("ROOT SVN version number not decimal, might be on git. Reading as hex %s -> %d",
		    tmpStr.Data(),fRootSvnVersion);
      }
      else {
	fRootSvnVersion=0;
	AliWarningF("ROOT SVN version is not extracted, setting to %d",fRootSvnVersion);
      }
    }
    delete tali;
  }
  else AliWarningF("Failed to extract %s version information",key[kRoot]);
  //
  // extract PASS
  if (!lpmRawPass.IsNull() && lpmRawPass.IsDigit()) fRecoPass = lpmRawPass.Atoi();
  else {
    AliWarningF("No %s record found, attempting to extract pass from OutputDir",key[kPass]);
    tmpStr = "/pass";
    if (outDir.IsNull() || !outDir.Contains(tmpStr)
	|| !sscanf(outDir.Data()+outDir.Index(tmpStr)+tmpStr.Length(),"%d",&fRecoPass)) 
      AliWarningF("Failed to extract pass number, set to %d",fRecoPass);
  }
  //
  // extract production type (RAW/MC)
  if (!lpmProdType.IsNull()) fMcFlag = (lpmProdType=="MC") ? kTRUE:kFALSE;
  else {
    AliWarningF("No %s record found, attempting to extract production type from OutputDir",key[kProdType]);
    if (lpmProdType.Contains("/alice/sim")) fMcFlag = kTRUE;
  }
  //
  // extract production tag
  if (!lpmProdTag.IsNull()) fProductionTag = lpmProdTag;
  else {
    AliWarningF("No %s record found, attempting to extract production tag from OutputDir",key[kProdTag]);
    tmpStr = "/LHC";
    if (!outDir.IsNull() && outDir.Contains(tmpStr)) {
      fProductionTag = outDir.Data()+outDir.Index(tmpStr)+1;
      if (fProductionTag.Contains("/")) fProductionTag.Resize(tmpStr.Index("/"));
    }
  }
  //
  // extract (anchored) period
  if (!lpmPeriod.IsNull()) fPeriod = lpmPeriod;
  else {
    AliWarningF("No %s record found, for raw data production tag %s will be assigned",key[kPeriod],fProductionTag.Data());
    if (!fMcFlag) fPeriod = fProductionTag;
  }
  //
  SetParsed(kTRUE);
}

//-------------------------------------------------------------------------------------------------	
void AliProdInfo::List() const {

  if (IsParsed()) {
    AliInfo("ALICE Production Info found in UserInfo: ");
    AliInfo(Form("  ALIROOT Version: %s [SVN #: %d]",fAlirootVersion.Data(),fAlirootSvnVersion));
    AliInfo(Form("  ROOT Version: %s [SVN #: %d]",fRootVersion.Data(),fRootSvnVersion));
    AliInfo(Form("  Reconstruction Pass: %d",fRecoPass));
    AliInfo(Form("  LHC Period: %s",fPeriod.Data()));
    AliInfo(Form("  ProductionTag: %s",fProductionTag.Data()));
    AliInfo(Form("  MC Flag: %d",fMcFlag));
  } else {
    AliInfo("ALICE Production Info not available in UserInfo");
  }
}
