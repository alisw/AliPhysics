#include "AliL3Logging.h"
#ifdef use_logging
#include "AliL3Logger.h"
#include <fstream.h>

int AliL3Logger::kAll= AliL3Log::kAll;
int AliL3Logger::kDebug = AliL3Log::kDebug;
int AliL3Logger::kInformational = AliL3Log::kInformational;
int AliL3Logger::kWarning = AliL3Log::kWarning;
int AliL3Logger::kError = AliL3Log::kError;
int AliL3Logger::kFatal = AliL3Log::kFatal;

AliL3Logger::AliL3Logger(){
  gLogLevel=AliL3Log::kAll;
  dn = so = se = sm =0; 
}
AliL3Logger::~AliL3Logger(){
  if(dn) {gLog.DelServer(dn);delete dn;}
  if(so) {gLog.DelServer(so);delete so;}
  if(se) {gLog.DelServer(se);delete se;}
  if(sm) {gLog.DelServer(sm);delete sm;}
  if(of) {of->close();delete of;}
}
void AliL3Logger::Set(int l){gLogLevel |=l;}
void AliL3Logger::UnSet(int l){gLogLevel &=(~l);}
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

void AliL3Logger::UseStream(char *name){
  if(sm) return;
//  static ofstream of;
  if(of) of->close();
  delete of;
  of =0;
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

#endif

