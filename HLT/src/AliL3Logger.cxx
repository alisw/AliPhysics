// @8#) $Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>, Constantin Loizides <mailto:loizides@ikf.uni-frankfurt.de>
//*-- Copyright &copy ALICE HLT Group

#include "AliL3StandardIncludes.h"

#include "AliL3Logging.h"
#include "AliL3Logger.h"

ClassImp(AliL3Logger)

#ifdef use_logging

Int_t AliL3Logger::kAll= AliL3Log::kAll;
Int_t AliL3Logger::kDebug = AliL3Log::kDebug;
Int_t AliL3Logger::kInformational = AliL3Log::kInformational;
Int_t AliL3Logger::kWarning = AliL3Log::kWarning;
Int_t AliL3Logger::kError = AliL3Log::kError;
Int_t AliL3Logger::kFatal = AliL3Log::kFatal;

AliL3Logger::AliL3Logger(){
  gLogLevel=AliL3Log::kAll;
  dn = so = se = sm =0;
  of = 0;
}
AliL3Logger::~AliL3Logger(){
  if(dn) {gLog.DelServer(dn);delete dn;}
  if(so) {gLog.DelServer(so);delete so;}
  if(se) {gLog.DelServer(se);delete se;}
  if(sm) {gLog.DelServer(sm);delete sm;}
  if(of) {of->close();delete of;}
}
void AliL3Logger::Set(Int_t l){gLogLevel |=l;}
void AliL3Logger::UnSet(Int_t l){gLogLevel &=(~l);}
void AliL3Logger::UseDevNull(){
  if(dn) return;
  dn = new AliL3DevNullLogServer();
  gLog.AddServer(dn);
}
void AliL3Logger::UseStdout(){
  if(so)return;
  so = new AliL3StdoutLogServer(); 
  gLog.AddServer(so);
}
void AliL3Logger::UseStderr(){
  if(se) return;
  se = new AliL3StderrLogServer();
  gLog.AddServer(se);
}

void AliL3Logger::UseStream(Char_t *name){
  if(sm) return;
  if(of) of->close();
  delete of;
  of = 0;
  of = new ofstream();
  of->open(name);
  sm = new  AliL3StreamLogServer(*of);
  gLog.AddServer(sm);
}
void AliL3Logger::NotUseDevNull(){
  if(dn) {gLog.DelServer(dn);delete dn;dn=0;}
}
void AliL3Logger::NotUseStdout(){
  if(so) {gLog.DelServer(so);delete so;so=0;}
}
void AliL3Logger::NotUseStderr(){
  if(se) {gLog.DelServer(se);delete se;se=0;}
}

void AliL3Logger::NotUseStream(){
  if(sm) {gLog.DelServer(sm);delete sm;sm=0;}
  if(of) {of->close();delete of;of=0;}
}
#else

Int_t AliL3Logger::kAll= AliL3Log::kAll;
Int_t AliL3Logger::kDebug = AliL3Log::kDebug;
Int_t AliL3Logger::kInformational = AliL3Log::kInformational;
Int_t AliL3Logger::kWarning = AliL3Log::kWarning;
Int_t AliL3Logger::kError = AliL3Log::kError;
Int_t AliL3Logger::kFatal = AliL3Log::kFatal;

AliL3Logger::AliL3Logger(){;}
AliL3Logger::~AliL3Logger(){;}
void AliL3Logger::Set(Int_t /*l*/){;}
void AliL3Logger::UnSet(Int_t /*l*/){;}
void AliL3Logger::UseDevNull(){;}
void AliL3Logger::UseStdout(){;}
void AliL3Logger::UseStderr(){;}
void AliL3Logger::UseStream(Char_t */*name*/){;}
void AliL3Logger::NotUseDevNull(){;}
void AliL3Logger::NotUseStdout(){;}
void AliL3Logger::NotUseStderr(){;}
void AliL3Logger::NotUseStream(){;}
#endif



