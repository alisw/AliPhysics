// @8#) $Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>, Constantin Loizides <mailto:loizides@ikf.uni-frankfurt.de>
//*-- Copyright &copy ALICE HLT Group

#include "AliL3StandardIncludes.h"

#include "AliL3Logging.h"
#include "AliL3Logger.h"

ClassImp(AliL3Logger)

#ifdef use_logging

Int_t AliL3Logger::fgAll= AliL3Log::kAll;
Int_t AliL3Logger::fgDebug = AliL3Log::kDebug;
Int_t AliL3Logger::fgInformational = AliL3Log::kInformational;
Int_t AliL3Logger::fgWarning = AliL3Log::kWarning;
Int_t AliL3Logger::fgError = AliL3Log::kError;
Int_t AliL3Logger::fgFatal = AliL3Log::kFatal;

AliL3Logger::AliL3Logger()
{
  //constructor
 if (!gLogP) {
   //printf( "Setting gLogP...\n" );
   //printf( "&gLogP: 0x%08lX\n", (unsigned long)&gLogP );
   //printf( "gLogP: 0x%08lX\n", (unsigned long)gLogP );
   gLogP = &MLUCLog::gLog;
   //printf( "gLogP set\n" );
 }
 if (!gLogLevelP){
   //printf( "Setting gLogLevelP...\n" );
   gLogLevelP = &MLUCLog::gLogLevel;
   //printf( "gLogLevelP set...\n" );
 }
 (*gLogLevelP)=AliL3Log::kAll;
  fdn = fso = fse = fsm =0;
  fof = 0;
}

AliL3Logger::~AliL3Logger()
{
  //destructor
  if(fdn) {gLog->DelServer(fdn);delete fdn;}
  if(fso) {gLog->DelServer(fso);delete fso;}
  if(fse) {gLog->DelServer(fse);delete fse;}
  if(fsm) {gLog->DelServer(fsm);delete fsm;}
  if(fof) {fof->close();delete fof;}
}

void AliL3Logger::Set(Int_t l)
{ 
  //set logger
  (*gLogLevel) |=l;
}

void AliL3Logger::UnSet(Int_t l)
{ 
  //unset logger
  (*gLogLevel) &=(~l);
}

void AliL3Logger::UseDevNull()
{
  //use dev null
  if(fdn) return;
  fdn = new AliL3DevNullLogServer();
  gLog->AddServer(dn);
}
void AliL3Logger::UseStdout()
{
  //use stdout
  if(fso)return;
  fso = new AliL3StdoutLogServer(); 
  gLog->AddServer(fso);
}
void AliL3Logger::UseStderr()
{
  //use stderr
  if(fse) return;
  fse = new AliL3StderrLogServer();
  gLog->AddServer(fse);
}

void AliL3Logger::UseStream(Char_t *name)
{
  //use stream
  if(fsm) return;
  if(fof) fof->close();
  delete fof;
  fof = 0;
  fof = new ofstream();
  fof->open(name);
  fsm = new  AliL3StreamLogServer(*fof);
  gLog->AddServer(fsm);
}

void AliL3Logger::NotUseDevNull()
{
  //not dev null
  if(fdn) {gLog->DelServer(fdn);delete fdn;fdn=0;}
}

void AliL3Logger::NotUseStdout()
{
  //not stdout
  if(fso) {gLog->DelServer(fso);delete fso;fso=0;}
}

void AliL3Logger::NotUseStderr()
{
  //not stderr
  if(fse) {gLog->DelServer(fse);delete fse;fse=0;}
}

void AliL3Logger::NotUseStream()
{
  //not stream
  if(fsm) {gLog->DelServer(fsm);delete fsm;fsm=0;}
  if(fof) {fof->close();delete fof;fof=0;}
}

#else /*not use_logging*/

Int_t AliL3Logger::fgAll= AliL3Log::kAll;
Int_t AliL3Logger::fgDebug = AliL3Log::kDebug;
Int_t AliL3Logger::fgInformational = AliL3Log::kInformational;
Int_t AliL3Logger::fgWarning = AliL3Log::kWarning;
Int_t AliL3Logger::fgError = AliL3Log::kError;
Int_t AliL3Logger::fgFatal = AliL3Log::kFatal;

AliL3Logger::AliL3Logger()
{
  //
  ;
}

AliL3Logger::~AliL3Logger()
{
  //
  ;
}

void AliL3Logger::Set(Int_t /*l*/)
{
  //
  ;
}

void AliL3Logger::UnSet(Int_t /*l*/)
{
  //
  ;
}

void AliL3Logger::UseDevNull()
{
  //
  ;
}

void AliL3Logger::UseStdout()
{
  //
  ;
}

void AliL3Logger::UseStderr()
{
  //
  ;
}

void AliL3Logger::UseStream(Char_t */*name*/)
{
  //
  ;
}

void AliL3Logger::NotUseDevNull()
{
  //
  ;
}

void AliL3Logger::NotUseStdout()
{
  //
  ;
}

void AliL3Logger::NotUseStderr()
{
  //
  ;
}

void AliL3Logger::NotUseStream()
{
  //
  ;
}
#endif



