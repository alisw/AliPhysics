// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveMacroExecutor.h"
#include "AliEveMacro.h"
#include "AliEveEventManager.h"
#include "AliSysInfo.h"

#include <TEveUtil.h>
#include <TList.h>
#include <TROOT.h>

#include <TEveManager.h>
#include <TGFileDialog.h>
#include <TGMenu.h>

#include <TSystem.h>
#include <TPRegexp.h>
#include <RVersion.h>


//______________________________________________________________________________
//
// Contains a list of AliEveMacros.
// The macros are added via AddMacro() and are owned by the executor.
// The macros can be executed via ExecMacros().
// They are executed in order in which they are registered.

ClassImp(AliEveMacroExecutor)

//______________________________________________________________________________
AliEveMacroExecutor::AliEveMacroExecutor() :
TObject(),
fMacros(new TList)
{
    // Constructor.
    
    fMacros->SetOwner(kTRUE);
}

//______________________________________________________________________________
AliEveMacroExecutor::~AliEveMacroExecutor()
{
    // Destructor.
    
    delete fMacros;
}

/******************************************************************************/

void AliEveMacroExecutor::AddMacro(AliEveMacro* mac)
{
    // Add a new macro. Ownership transfered to the executor.
    
    static const TEveException kEH("AliEveMacroExecutor::AddMacro ");
    
    const TString mname = mac->GetMacro();
    if ( ! mname.IsNull() && TEveUtil::CheckMacro(mname) == kFALSE)
    {
        TEveUtil::LoadMacro(mname);
    }
    fMacros->Add(mac);
}

AliEveMacro* AliEveMacroExecutor::FindMacro(const TString& func)
{
    // Find macro with given function name (it is supposed to be unique).
    // Returns 0 if not found.
    
    TIter next(fMacros);
    AliEveMacro* mac;
    while ((mac = (AliEveMacro*) next()))
    {
        if (mac->GetFunc() == func)
            return mac;
    }
    return 0;
}

/******************************************************************************/

void AliEveMacroExecutor::RemoveMacros()
{
    fMacros->Clear();
}

/******************************************************************************/

#if ROOT_VERSION_CODE < ROOT_VERSION(5,99,0)
#include "Api.h"
#endif
#include "TInterpreter.h"

void AliEveMacroExecutor::ExecMacros()
{
    // Execute registered macros.
    
    TIter next(fMacros);
    AliEveMacro* mac;
    while ((mac = (AliEveMacro*) next()))
    {
        if (mac->GetFunc().IsNull())
        {
            continue;
        }
                
        TString cmd(mac->FormForExec());
        try
        {
            Long_t                   result = 0;
            TInterpreter::EErrorCode error  = TInterpreter::kNoError;
            
            AliSysInfo::AddStamp(Form("%s_%s_before",mac->GetMacro().Data(), mac->GetFunc().Data()));
            result = gInterpreter->ProcessLine(cmd, &error);
            AliSysInfo::AddStamp(Form("%s_%s_after",mac->GetMacro().Data(), mac->GetFunc().Data()));
            
#if ROOT_VERSION_CODE < ROOT_VERSION(5,99,0)
            // Try to fix broken cint state? Code taken form pyroot.
            if (G__get_return(0) > G__RETURN_NORMAL)
            {
                printf ("*** FIXING CINT STATE AFTER RETURN ***\n");
                G__security_recover(0);
            }
#endif
            if (error != TInterpreter::kNoError)
            {
                Error("ExecMacros", "Executing %s::%s, CINT error ... hopefully recovered.",
                      mac->GetMacro().Data(), cmd.Data());
            }
            else
            {
                TEveElement *el  = (TEveElement*) result;
                TObject     *obj = dynamic_cast<TObject*>(el);
                if (el != 0 && obj == 0)
                {
                    Warning("ExecMacros", "Executing %s::%s, returned TEveElement seems bad, setting it to 0.",
                            mac->GetMacro().Data(), cmd.Data());
                    el = 0;
                }
            }
        }
        catch(TEveException& exc)
        {
#if ROOT_VERSION_CODE < ROOT_VERSION(5,99,0)
            // Try to fix broken cint state? Code taken form pyroot.
            if (G__get_return(0) > G__RETURN_NORMAL)
            {
                printf ("*** FIXING CINT STATE AFTER EXCEPTION ***\n");
                G__security_recover(0);
            }
#endif
            Error("ExecMacros", "Executing %s::%s, caught exception: '%s'.",
                  mac->GetMacro().Data(), cmd.Data(), exc.Data());
        }
    }
}

/******************************************************************************/

#include <iostream>
#include <fstream>
using namespace std;

namespace
{
    const char *gMacroSaveAsTypes[] = {"CINT Macro", "*.C",
        0, 0};
}

//void AliEveMacroExecutor::SaveAddedMacros()
//{
//    
//    TGFileInfo fi;
//    fi.fFileTypes   = gMacroSaveAsTypes;
//    fi.fIniDir      = StrDup(""); // current directory
//    fi.fFileTypeIdx = 0;
//    fi.fOverwrite   = kTRUE;
//    new TGFileDialog(gClient->GetDefaultRoot(), gEve->GetMainWindow(), kFDSave, &fi);
//    if (!fi.fFilename) return;
//    
//    TPMERegexp filere(".*/([^/]+$)");
//    if (filere.Match(fi.fFilename) != 2)
//    {
//        Warning("AliEvePopupHandler", "file '%s' bad.", fi.fFilename);
//        return;
//    }
//    printf("Saving...\n");
//    
//    TString file(filere[1]);
//    TString file1;
//    if (!file.EndsWith(".C"))
//        file1 = file + ".C";
//    gSystem->ChangeDirectory(fi.fIniDir);
//    ofstream myfile;
//    myfile.open (file1);
//    
//    TIter next(fMacros);
//    AliEveMacro* mac;
//    
//    
//    myfile <<"//Macro generated automatically by AliEveMacroExecutor\n\n";
//    
//    myfile <<"void "<<file<<"(){\n\n";
//    myfile <<"  AliEveMacroExecutor *exec = AliEveEventManager::Instance()->GetExecutor();\n";
//    myfile <<"  exec->RemoveMacros();\n";
//    myfile <<"  TEveBrowser *browser = gEve->GetBrowser();\n";
//    myfile <<"  browser->ShowCloseTab(kFALSE);\n";
//    
//    while ((mac = (AliEveMacro*) next()))
//    {
//        myfile <<"  exec->AddMacro(new AliEveMacro("<<mac->GetSources()<<", "<<char(34)<<mac->GetTags()<<char(34)<<", "
//        <<char(34)<<mac->GetMacro()<<char(34)<<", "<<char(34)<<mac->GetFunc()<<char(34)<<", "<<char(34)<<mac->GetArgs()
//        <<char(34)<<", "<<mac->GetActive()<<"));\n\n";
//    }
//    
//    myfile <<"  TEveWindowSlot *slot = TEveWindow::CreateWindowInTab(browser->GetTabRight());\n";
//    myfile <<"  slot->StartEmbedding();\n";
//    myfile <<"  AliEveMacroExecutorWindow* exewin = new AliEveMacroExecutorWindow(exec);\n";
//    myfile <<"  slot->StopEmbedding("<<char(34)<<"DataSelection"<<char(34)<<");\n";
//    myfile <<"  exewin->PopulateMacros();\n\n";
//    
//    myfile <<"\n}";
//    myfile.close();
//    printf("Saved...\n");
//    
//}
