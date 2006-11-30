// @8#) $Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>, Constantin Loizides <mailto:loizides@ikf.uni-frankfurt.de>
//*-- Copyright &copy ALICE HLT Group

#include "AliHLTStandardIncludes.h"

#include "AliHLTLogging.h"
#include "AliHLTLogger.h"

ClassImp(AliHLTLogger)

#ifdef use_logging

Int_t AliHLTLogger::fgAll= AliHLTLog::kAll;
Int_t AliHLTLogger::fgDebug = AliHLTLog::kDebug;
Int_t AliHLTLogger::fgInformational = AliHLTLog::kInformational;
Int_t AliHLTLogger::fgWarning = AliHLTLog::kWarning;
Int_t AliHLTLogger::fgError = AliHLTLog::kError;
Int_t AliHLTLogger::fgFatal = AliHLTLog::kFatal;

AliHLTLogger::AliHLTLogger()
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
 (*gLogLevelP)=AliHLTLog::kAll;
  fdn = fso = fse = fsm =0;
  fof = 0;
}

AliHLTLogger::~AliHLTLogger()
{
  //destructor
  if(fdn) {gLogP->DelServer(fdn);delete fdn;}
  if(fso) {gLogP->DelServer(fso);delete fso;}
  if(fse) {gLogP->DelServer(fse);delete fse;}
  if(fsm) {gLogP->DelServer(fsm);delete fsm;}
  if(fof) {fof->close();delete fof;}
}

void AliHLTLogger::Set(Int_t l)
{ 
  //set logger
  (*gLogLevelP) |=l;
}

void AliHLTLogger::UnSet(Int_t l)
{ 
  //unset logger
  (*gLogLevelP) &=(~l);
}

void AliHLTLogger::UseDevNull()
{
  //use dev null
  if(fdn) return;
  fdn = new AliHLTDevNullLogServer();
  gLogP->AddServer(fdn);
}
void AliHLTLogger::UseStdout()
{
  //use stdout
  if(fso)return;
  fso = new AliHLTStdoutLogServer(); 
  gLogP->AddServer(fso);
}
void AliHLTLogger::UseStderr()
{
  //use stderr
  if(fse) return;
  fse = new AliHLTStderrLogServer();
  gLogP->AddServer(fse);
}

void AliHLTLogger::UseStream(Char_t *name)
{
  //use stream
  if(fsm) return;
  if(fof) fof->close();
  delete fof;
  fof = 0;
  fof = new ofstream();
  fof->open(name);
  fsm = new  AliHLTStreamLogServer(*fof);
  gLogP->AddServer(fsm);
}

void AliHLTLogger::NotUseDevNull()
{
  //not dev null
  if(fdn) {gLogP->DelServer(fdn);delete fdn;fdn=0;}
}

void AliHLTLogger::NotUseStdout()
{
  //not stdout
  if(fso) {gLogP->DelServer(fso);delete fso;fso=0;}
}

void AliHLTLogger::NotUseStderr()
{
  //not stderr
  if(fse) {gLogP->DelServer(fse);delete fse;fse=0;}
}

void AliHLTLogger::NotUseStream()
{
  //not stream
  if(fsm) {gLogP->DelServer(fsm);delete fsm;fsm=0;}
  if(fof) {fof->close();delete fof;fof=0;}
}

#else /*not use_logging*/

Int_t AliHLTLogger::fgAll= AliHLTLog::kAll;
Int_t AliHLTLogger::fgDebug = AliHLTLog::kDebug;
Int_t AliHLTLogger::fgInformational = AliHLTLog::kInformational;
Int_t AliHLTLogger::fgWarning = AliHLTLog::kWarning;
Int_t AliHLTLogger::fgError = AliHLTLog::kError;
Int_t AliHLTLogger::fgFatal = AliHLTLog::kFatal;

AliHLTLogger::AliHLTLogger()
{
  //
  ;
}

AliHLTLogger::~AliHLTLogger()
{
  //
  ;
}

void AliHLTLogger::Set(Int_t /*l*/)
{
  //
  ;
}

void AliHLTLogger::UnSet(Int_t /*l*/)
{
  //
  ;
}

void AliHLTLogger::UseDevNull()
{
  //
  ;
}

void AliHLTLogger::UseStdout()
{
  //
  ;
}

void AliHLTLogger::UseStderr()
{
  //
  ;
}

void AliHLTLogger::UseStream(Char_t */*name*/)
{
  //
  ;
}

void AliHLTLogger::NotUseDevNull()
{
  //
  ;
}

void AliHLTLogger::NotUseStdout()
{
  //
  ;
}

void AliHLTLogger::NotUseStderr()
{
  //
  ;
}

void AliHLTLogger::NotUseStream()
{
  //
  ;
}
#endif
