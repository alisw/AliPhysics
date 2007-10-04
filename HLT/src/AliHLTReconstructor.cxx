// $Id$

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for HLT reconstruction                                              //
// <Cvetan.Cheshkov@cern.ch>                                                 //
// <loizides@ikf.uni-frankfurt.de>                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TSystem.h>
#include <TObjString.h>
#include "AliHLTReconstructor.h"
#include "AliLog.h"
#include "AliHLTSystem.h"

ClassImp(AliHLTReconstructor)

AliHLTReconstructor::AliHLTReconstructor()
  : 
  AliReconstructor(),
  fpSystem(NULL)
{ 
  //constructor
}

AliHLTReconstructor::~AliHLTReconstructor()
{ 
  //destructor
  if (fpSystem) {
    delete fpSystem;
  }
  fpSystem=NULL;
}

void AliHLTReconstructor::Init()
{
  // init the reconstructor
  if (!fpSystem) fpSystem=new AliHLTSystem;
  if (!fpSystem) {
    AliError("can not create AliHLTSystem object");
    return;
  }
  if (fpSystem->CheckStatus(AliHLTSystem::kError)) {
    AliError("HLT system in error state");
    return;
  }

  // the options scan has been moved to AliHLTSystem, the old code
  // here is kept to be able to run an older version of the HLT code
  // with newer AliRoot versions.
  TString libs("");
  TString option = GetOption();
  TObjArray* pTokens=option.Tokenize(" ");
  option="";
  if (pTokens) {
    int iEntries=pTokens->GetEntries();
    for (int i=0; i<iEntries; i++) {
      TString token=(((TObjString*)pTokens->At(i))->GetString());
      if (token.Contains("loglevel=")) {
	TString param=token.ReplaceAll("loglevel=", "");
	if (param.IsDigit()) {
	  fpSystem->SetGlobalLoggingLevel((AliHLTComponentLogSeverity)param.Atoi());
	} else if (param.BeginsWith("0x") &&
		   param.Replace(0,2,"",0).IsHex()) {
	  int severity=0;
	  sscanf(param.Data(),"%x", &severity);
	  fpSystem->SetGlobalLoggingLevel((AliHLTComponentLogSeverity)severity);
	} else {
	  AliWarning("wrong parameter for option \'loglevel=\', (hex) number expected");
	}
      } else if (token.Contains("alilog=off")) {
	fpSystem->SwitchAliLog(0);
      } else if (token.BeginsWith("lib") && token.EndsWith(".so")) {
	libs+=token;
	libs+=" ";
      } else {
	if (option.Length()>0) option+=" ";
	option+=token;
      }
    }
    delete pTokens;
  }
  
  if (!libs.IsNull() &&
      (!fpSystem->CheckStatus(AliHLTSystem::kLibrariesLoaded)) &&
      (fpSystem->LoadComponentLibraries(libs.Data())<0)) {
    AliError("error while loading HLT libraries");
    return;
  }

  if (!fpSystem->CheckStatus(AliHLTSystem::kReady)) {
    typedef int (*AliHLTSystemSetOptions)(AliHLTSystem* pInstance, const char* options);
    gSystem->Load("libHLTinterface.so");
    AliHLTSystemSetOptions pFunc=(AliHLTSystemSetOptions)(gSystem->DynFindSymbol("libHLTinterface.so", "AliHLTSystemSetOptions"));
    if (pFunc) {
      if ((pFunc)(fpSystem, option.Data())<0) {
      AliError("error setting options for HLT system");
      return;	
      }
    } else if (option.Length()>0) {
      AliError(Form("version of HLT system does not support the options \'%s\'", option.Data()));
      return;
    }
    if ((fpSystem->Configure())<0) {
      AliError("error during HLT system configuration");
      return;
    }
  }
}

void AliHLTReconstructor::Reconstruct(AliRawReader* /*rawReader*/, TTree* /*clustersTree*/) const 
{
  // reconstruction of real data without writing of ESD

  // all reconstruction has been moved to FillESD
//   int iResult=0;
//   if (fpSystem) {
//     if (fpSystem->CheckStatus(AliHLTSystem::kError)) {
//       AliError("HLT system in error state");
//       return;
//     }
//     if ((iResult=fpSystem->Reconstruct(1, NULL, rawReader))>=0) {
//     }
//   }
}

void AliHLTReconstructor::FillESD(AliRawReader* rawReader, TTree* /*clustersTree*/, 
				  AliESDEvent* esd) const
{
  // reconstruct real data and fill ESD
  int iResult=0;
  if (fpSystem) {
    if (fpSystem->CheckStatus(AliHLTSystem::kError)) {
      AliError("HLT system in error state");
      return;
    }
    if (!fpSystem->CheckStatus(AliHLTSystem::kReady)) {
      AliError("HLT system in wrong state");
      return;
    }
    if ((iResult=fpSystem->Reconstruct(1, NULL, rawReader))>=0) {
      fpSystem->FillESD(-1, NULL, esd);
    }
  }
}
