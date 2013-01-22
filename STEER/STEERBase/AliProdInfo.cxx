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

AliProdInfo::AliProdInfo():
  fPeriod(""),
  fAlirootVersion(""),
  fAlirootSvnVersion(0),
  fRootVersion(""),
  fRootSvnVersion(0),
  fMcFlag(kFALSE),
  fRecoPass(-1)
{	
  //
  // default constructor
  //	
  AliLog::SetClassDebugLevel("AliProdInfo",10);
}

AliProdInfo::AliProdInfo(const TString& name,const TString& title):
  TNamed(name,title),
  fPeriod(""),
  fAlirootVersion(""),
  fAlirootSvnVersion(0),
  fRootVersion(""),
  fRootSvnVersion(0),
  fMcFlag(kFALSE),
  fRecoPass(-1)
{	
  //
  // default constructor
  //	
  AliLog::SetClassDebugLevel("AliProdInfo",10);
}

AliProdInfo::AliProdInfo(TList *userInfoList):
  fPeriod(""),
  fAlirootVersion(""),
  fAlirootSvnVersion(0),
  fRootVersion(""),
  fRootSvnVersion(0),
  fMcFlag(kFALSE),
  fRecoPass(-1)
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
void AliProdInfo::Init(TList *userInfoList) {
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
void AliProdInfo::ParseProdInfo(TNamed *prodInfoData) {

  TString str(prodInfoData->GetTitle());
  TObjArray *tokens = str.Tokenize(";");

  for (Int_t i=0;i<=tokens->GetLast();i++) {
    TObjString *stObj = (TObjString *)tokens->At(i);
    
    if (stObj->GetString().Contains("aliroot") && (i==0) ) {  // aliroot version
      TObjArray *tali = (TObjArray *)stObj->GetString().Tokenize(":");
      TObjString *tos = (TObjString *)tali->At(0);
      TObjArray *tali2 = (TObjArray *)tos->GetString().Tokenize(" ");
      TObjString *av = (TObjString *)tali2->At(1);
      fAlirootVersion=av->GetString().Data();
      TObjString *ts = (TObjString*)tali->At(1); 
      fAlirootSvnVersion = ts->GetString().Atoi();
      delete tali;
      delete tali2;
    }
    else if (stObj->GetString().Contains("root") && (i==1) ) {  // root version
      TObjArray *tali = (TObjArray *)stObj->GetString().Tokenize(":");
      TObjString *tos = (TObjString *)tali->At(0);
      TObjArray *tali2 = (TObjArray *)tos->GetString().Tokenize(" ");
      TObjString *av = (TObjString *)tali2->At(1);
      fRootVersion=av->GetString().Data();
      TObjString *ts = (TObjString*)tali->At(1); 
      fRootSvnVersion = ts->GetString().Atoi();
      delete tali;
      delete tali2;
    }
    else if (stObj->GetString().Contains("OutputDir")) {
	  if (stObj->GetString().Contains("pass1") ) {
	    fRecoPass=1;
	  } else if (stObj->GetString().Contains("pass2") ) {
	    fRecoPass=2;
	  } else if (stObj->GetString().Contains("pass3") ) {
	    fRecoPass=3;
	  } else if (stObj->GetString().Contains("pass4") ) {
	    fRecoPass=4;
	  } else if (stObj->GetString().Contains("pass5") ) {
	    fRecoPass=5;
	  }
	  if (stObj->GetString().Contains("/alice/sim") ) fMcFlag = kTRUE;
	  else {
	    TObjArray *tit = (TObjArray *)stObj->GetString().Tokenize("=");
	    TObjString *tos = (TObjString*)tit->At(1);
	    tit=(TObjArray *)tos->GetString().Tokenize("/");
            tos=(TObjString*)tit->At(3);
            if (tos) {
              if (tos->GetString().Contains("_")) {
		tit=(TObjArray *)tos->GetString().Tokenize("_");
		tos=(TObjString*)tit->At(0);
	      }
	      if (tos) fPeriod=tos->GetString().Data();
	    }
	    delete tit;
	  }
    }
    //    else if (stObj->GetString().Contains("LPMProductionType=MC") ) {
    //      fMcFlag = kTRUE;
    //    }
    else if (stObj->GetString().Contains("LPMAnchorProduction") ) {
      TObjArray *tit = (TObjArray *)stObj->GetString().Tokenize("=");
      TObjString *tos = (TObjString *)tit->At(1);
      fPeriod=tos->GetString().Data();
      delete tit;
    }
    else if (stObj->GetString().Contains("LPMProductionTag")) {
      // for the moment we don't parse the tag, given unsafe standards...
      /*
      TObjArray *tit = (TObjArray *)stObj->GetString().Tokenize(" ");
      if (tit->GetLast()>1) {
	TObjString *tos = (TObjString *)tit->At(2);
	if (tos) fPeriod=tos->GetString().Data();
      }
      */
    }
  }
  delete tokens;
}

//-------------------------------------------------------------------------------------------------	
void AliProdInfo::List() const {

  if (fAlirootSvnVersion > 0) {
    AliInfo("ALICE Production Info found in UserInfo: ");
    AliInfo(Form("  ALIROOT Version: %s [SVN #: %d]",fAlirootVersion.Data(),fAlirootSvnVersion));
    AliInfo(Form("  ROOT Version: %s [SVN #: %d]",fRootVersion.Data(),fRootSvnVersion));
    if (!fMcFlag) AliInfo(Form("  Reconstruction Pass: %d",fRecoPass));
    AliInfo(Form("  LHC Period: %s",fPeriod.Data()));
    AliInfo(Form("  MC Flag: %d",fMcFlag));
  } else {
    AliInfo("ALICE Production Info not available in UserInfo");
  }
}
