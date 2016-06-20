/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

//************************************************************************
//  Author: Mihai Niculescu (mihai.niculescu@cern.ch)
//  Jan 31 2012
//************************************************************************

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TList.h>
#include <TObjString.h>
#include <TString.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TEveUtil.h>
#endif

/* Compiles all AliEve macros
 * 
 * Usage: 
 *  1. cd $ALICE_ROOT/EVE/test-macros/
 *  2. Launch alieve from terminal: alieve
 *  3. from terminal execute this script: .x compileEveMacros.C
 *  3. Wait for compilation to finish
 * 
 * NOTE:
 *  At the end of compilation,you might get a list 
 * with all macros that failed compilation
 * 
 * Parameters:
 *  macDir  - absolute directory path where to find the macros
 *  opt     - compilation options (See: TSystem::CompileMacro)
 * 
 * Default:
 *  - compiles all macros from AliRoot installation dir ($ALICE_ROOT/EVE/macros)
 *  - options for compilation (- k f c):
 * The possible options are:
     k : keep the shared library after the session end.
     f : force recompilation.
     g : compile with debug symbol
     O : optimized the code (ignore if 'g' is specified)
     c : compile only, do not attempt to load the library.
     - : if buildir is set, use a flat structure (see buildir below)
 */

void compileEveMacros( const char * macDir="", Option_t *opt="")
{

   // solves current issue in ROOT when pre-compiling with ACLiC
   gSystem->Load("libCint");
   gSystem->Load("libTender");
   gSystem->Load("libPWGPP");

  if(macDir == "")
    macDir = Form("%s/EVE/macros", gSystem->Getenv("ALICE_ROOT") );
  
  if(opt == "")
    opt = "-kc"; // compilation options
  
  TObjString *mac;
  TList * listOfFailedMacros = new TList; // list of macros that failed compilation
    
  TSystemDirectory *curDir = new TSystemDirectory(macDir, macDir);
  TSystemFile *curFile = 0;
  
  TList* listOfMacros = curDir->GetListOfFiles();
  TListIter next(listOfMacros);
  
  TPMERegexp regex("\\.C$");
    
  printf("Directory: %s \tfiles:%d\n",  curDir->GetName(), listOfMacros->GetSize() );
  gSystem->cd(curDir->GetName());

  const char* incPath = gSystem->GetIncludePath();
  
  gSystem->SetIncludePath(Form("%s -I%s", incPath, "$ROOTSYS/include"));
  
  
  Int_t nMacros = 0;
  
  while( (curFile = static_cast<TSystemFile*>(next())) )
  {
    mac = new TObjString(curFile->GetName());
    
    if( regex.Match(mac->String().Data() ) > 0 )
    {
        nMacros++;
      
        printf("Macro %s\n", mac->String().Data() );
        TEveUtil::CheckMacro(mac->String().Data() );
        if(!gSystem->CompileMacro(mac->String().Data(), opt ) )
          listOfFailedMacros->Add(mac);
    }
    
  }
  
  printf("\n\nTotal Macros:%d \tFailed:%d \n", nMacros, listOfFailedMacros->GetSize() );
  
  printf("\nFollowing macros failed compilation: \n");
  TListIter failed(listOfFailedMacros);
  while( (mac = static_cast<TObjString*>(failed())) )
  {
    printf( "%s\n", mac->String().Data() );
  }
  
  delete listOfFailedMacros;
  delete listOfMacros;
  delete curDir;
}
