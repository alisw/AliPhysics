#include <TGLabel.h>
#include <TGComboBox.h>
#include <TGButton.h>
#include <TG3DLine.h>
#include <TGTextEntry.h>
#include <TMath.h>
#include <TSystem.h>

#include <TGFileDialog.h>

#include "AliEveUtil.h"
#include "AliEveFileDialog.h"

ClassImp(AliEveFileDialog)

AliEveFileDialog::AliEveFileDialog(const TGWindow* p,const TGWindow* main, EAliEveFileDialogMode mode)
    : TGTransientFrame(p, main, 100,50, kVerticalFrame),
    fESDFilesFrame(0),
    fAdvancedOptsButton(0),
    fAdvancedOptsFrame(0),
    fESDfriendFilesFrame(0),
    fAODFilesFrame(0),
    fAODfriendFilesFrame(0),
    fRawFilesFrame(0),
    fUrlFrame(0),
    fCDBFrame(0),
    fDialogButtonsFrame(0),
    fPathEntryESD(0),
    fPathEntryESDfriend(0),
    fPathEntryAOD(0),
    fPathEntryAODfriend(0),
    fPathEntryRawFile(0),
    fPathEntryUrl(0),
    fCDBPathCB(0),
    fIsAccepted(kFALSE),
    fMode(mode)
{
    SetCleanup(kDeepCleanup);

    /************************
   *** Local Files Frame ***
   *************************/
    //ESD Frame
    fESDFilesFrame = new TGHorizontalFrame(this, 100, 100);
    TGLabel* esdLabel = new TGLabel(fESDFilesFrame, "ESD File:");
    esdLabel->Resize(110, esdLabel->GetDefaultHeight());
    esdLabel->SetMargins(0,0,0,0);

    fPathEntryESD= new TGTextEntry(fESDFilesFrame);
    fPathEntryESD->Resize(250, fPathEntryESD->GetDefaultHeight());

    TGTextButton* browseButtonESD = new TGTextButton(fESDFilesFrame, "Browse...");
    browseButtonESD->Connect("Clicked()", "AliEveFileDialog", this, "onBrowseESDFile()");

    fESDFilesFrame->AddFrame(esdLabel, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 3, 3, 3, 3));
    fESDFilesFrame->AddFrame(fPathEntryESD, new TGLayoutHints(kLHintsExpandX, 3, 3, 3, 3));
    fESDFilesFrame->AddFrame(browseButtonESD, new TGLayoutHints(kLHintsNormal, 3, 3, 3, 3));

    AddFrame(fESDFilesFrame, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 3, 3, 3, 3));
    
    //Advanced Options Button
    fAdvancedOptsButton = new TGCheckButton(this, "Advanced Options");
    fAdvancedOptsButton->Connect("Toggled(Bool_t)", "AliEveFileDialog", this, "showAdvancedOpts(Bool_t)");

    AddFrame(fAdvancedOptsButton,new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 3, 3, 3, 3));

//Advanced Options Frame
    //ESDfriend Frame
    fAdvancedOptsFrame = new TGGroupFrame(this, "Advanced Options");
    
    fESDfriendFilesFrame = new TGHorizontalFrame(fAdvancedOptsFrame, 100, 100);
    TGLabel* esdFriendLabel = new TGLabel(fESDfriendFilesFrame, "ESDfriend File:");
    esdFriendLabel->Resize(110, esdFriendLabel->GetDefaultHeight());
    esdFriendLabel->SetMargins(0,0,0,0);

    fPathEntryESDfriend= new TGTextEntry(fESDfriendFilesFrame);
    fPathEntryESDfriend->Resize(250, fPathEntryESDfriend->GetDefaultHeight());

    TGTextButton* browseButtonESDfriend = new TGTextButton(fESDfriendFilesFrame, "Browse...");
    browseButtonESDfriend->Connect("Clicked()", "AliEveFileDialog", this, "onBrowseESDfriendFile()");

    fESDfriendFilesFrame->AddFrame(esdFriendLabel, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 3, 3, 3, 3));
    fESDfriendFilesFrame->AddFrame(fPathEntryESDfriend, new TGLayoutHints(kLHintsExpandX, 3, 3, 3, 3));
    fESDfriendFilesFrame->AddFrame(browseButtonESDfriend, new TGLayoutHints(kLHintsNormal, 3, 3, 3, 3));

    fAdvancedOptsFrame->AddFrame(fESDfriendFilesFrame, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 3, 3, 3, 3));
    
    // AOD Frame
    fAODFilesFrame = new TGHorizontalFrame(fAdvancedOptsFrame, 100, 100);
    TGLabel* aodLabel = new TGLabel(fAODFilesFrame, "AOD File:");
    aodLabel->Resize(110, aodLabel->GetDefaultHeight());

    fPathEntryAOD = new TGTextEntry(fAODFilesFrame);
    fPathEntryAOD->Resize(250, fPathEntryAOD->GetDefaultHeight());

    TGTextButton* browseButtonAOD = new TGTextButton(fAODFilesFrame, "Browse...");
    browseButtonAOD->Connect("Clicked()", "AliEveFileDialog", this, "onBrowseAODFile()");

    fAODFilesFrame->AddFrame(aodLabel, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 3, 3, 3, 3));
    fAODFilesFrame->AddFrame(fPathEntryAOD, new TGLayoutHints(kLHintsExpandX, 3, 3, 3, 3));
    fAODFilesFrame->AddFrame(browseButtonAOD, new TGLayoutHints(kLHintsNormal, 3, 3, 3, 3));

    fAdvancedOptsFrame->AddFrame(fAODFilesFrame, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 3, 3, 3, 3));

    // AODfriend Frame
    fAODfriendFilesFrame = new TGHorizontalFrame(fAdvancedOptsFrame, 100, 100);
    TGLabel* aodFriendLabel = new TGLabel(fAODfriendFilesFrame, "AODfriend File:");
    aodFriendLabel->Resize(110, aodFriendLabel->GetDefaultHeight());

    fPathEntryAODfriend = new TGTextEntry(fAODfriendFilesFrame);
    fPathEntryAODfriend->Resize(250, fPathEntryAODfriend->GetDefaultHeight());

    TGTextButton* browseButtonAODfriend = new TGTextButton(fAODfriendFilesFrame, "Browse...");
    browseButtonAODfriend->Connect("Clicked()", "AliEveFileDialog", this, "onBrowseAODfriendFile()");

    fAODfriendFilesFrame->AddFrame(aodFriendLabel, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 3, 3, 3, 3));
    fAODfriendFilesFrame->AddFrame(fPathEntryAODfriend, new TGLayoutHints(kLHintsExpandX, 3, 3, 3, 3));
    fAODfriendFilesFrame->AddFrame(browseButtonAODfriend, new TGLayoutHints(kLHintsNormal, 3, 3, 3, 3));

    fAdvancedOptsFrame->AddFrame(fAODfriendFilesFrame, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 3, 3, 3, 3));

    // Raw Frame
    fRawFilesFrame = new TGHorizontalFrame(fAdvancedOptsFrame, 100, 100);
    TGLabel* rawLabel = new TGLabel(fRawFilesFrame, "Raw File:");
    rawLabel->Resize(110, rawLabel->GetDefaultHeight());

    fPathEntryRawFile = new TGTextEntry(fRawFilesFrame);
    fPathEntryRawFile->Resize(250, fPathEntryRawFile->GetDefaultHeight());

    TGTextButton* browseButtonRawFile = new TGTextButton(fRawFilesFrame, "Browse...");
    browseButtonRawFile->Connect("Clicked()", "AliEveFileDialog", this, "onBrowseRawFile()");

    fRawFilesFrame->AddFrame(rawLabel, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 3, 3, 3, 3));
    fRawFilesFrame->AddFrame(fPathEntryRawFile, new TGLayoutHints(kLHintsExpandX, 3, 3, 3, 3));
    fRawFilesFrame->AddFrame(browseButtonRawFile, new TGLayoutHints(kLHintsNormal, 3, 3, 3, 3));

    fAdvancedOptsFrame->AddFrame(fRawFilesFrame, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 3, 3, 3, 3));

   AddFrame(fAdvancedOptsFrame, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 3, 3, 3, 3) );

   /****************
   *** Url Frame ***
   *****************/
    fUrlFrame = new TGHorizontalFrame(this, 100, 100);
    TGLabel* urlLabel = new TGLabel(fUrlFrame, "Url:");
    urlLabel->Resize(110, urlLabel->GetDefaultHeight());

    fPathEntryUrl= new TGTextEntry(fUrlFrame, "alien:///alice/");
    fPathEntryUrl->Resize(250, fPathEntryUrl->GetDefaultHeight());

    fUrlFrame->AddFrame(urlLabel, new TGLayoutHints(kLHintsLeft | kLHintsCenterY , 3, 3, 3, 3));
    fUrlFrame->AddFrame(fPathEntryUrl, new TGLayoutHints(kLHintsExpandX, 3, 3, 3, 3));

    AddFrame(fUrlFrame, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 3, 3, 3, 3));

    /****************
   *** CDB Frame ***
   *****************/
    fCDBFrame = new TGHorizontalFrame(this, 100, 100);
    TGLabel *cdbLabel = new TGLabel(fCDBFrame, "CDB Storage:");
    fCDBPathCB = new TGComboBox(fCDBFrame);
    fCDBPathCB->SetEditable(kFALSE);

    fCDBPathCB->AddEntry("local", AliEveUtil::CDB_LOCAL);
    fCDBPathCB->AddEntry("raw", AliEveUtil::CDB_RAW);
    fCDBPathCB->AddEntry("mcideal", AliEveUtil::CDB_MCIDEAL);
    fCDBPathCB->AddEntry("mcresidual", AliEveUtil::CDB_MCRESIDUAL);
    fCDBPathCB->AddEntry("mcfull", AliEveUtil::CDB_MCFULL);

    fCDBPathCB->Resize(100, fPathEntryESD->GetDefaultHeight());
    fCDBPathCB->Select(0, kTRUE);

    fCDBFrame->AddFrame(cdbLabel, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 3, 3, 3, 3) );
    fCDBFrame->AddFrame(fCDBPathCB,  new TGLayoutHints(kLHintsNormal | kLHintsCenterY, 3, 3, 3, 3) );

    AddFrame(fCDBFrame, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 3, 3, 3, 3) );

    /***************************
   *** Dialog Buttons Frame ***
   ****************************/
    fDialogButtonsFrame = new TGHorizontalFrame(this, 100, 100);

    TGTextButton* okButton = new TGTextButton(fDialogButtonsFrame, "OK");
    okButton->Associate(this);
    okButton->Resize(50, okButton->GetDefaultHeight());
    okButton->Connect("Clicked()", "AliEveFileDialog", this, "onAccept()");

    TGTextButton* cancelButton = new TGTextButton(fDialogButtonsFrame, "Cancel");
    okButton->Associate(this);
    cancelButton->Connect("Clicked()", "AliEveFileDialog", this, "onReject()");
    cancelButton->Resize(50, cancelButton->GetDefaultHeight());

    fDialogButtonsFrame->AddFrame(new TGHorizontal3DLine, new TGLayoutHints(kLHintsExpandX, 3, 3, 3, 3));
    fDialogButtonsFrame->AddFrame(okButton,    new TGLayoutHints(kLHintsExpandX | kLHintsCenterY, 3, 3, 3, 3));
    fDialogButtonsFrame->AddFrame(cancelButton,new TGLayoutHints(kLHintsExpandX | kLHintsCenterY, 3, 3, 3, 3));

    AddFrame(fDialogButtonsFrame, new TGLayoutHints(kLHintsNormal | kLHintsExpandX , 3, 3, 3, 3) );

    setMode(mode);

    TGDimension size = GetDefaultSize();

    Resize(size);
    SetWindowName("Open File...");

    Layout();
    HideFrame(this);
    UnmapWindow();
}

AliEveFileDialog::~AliEveFileDialog()
{

}

const TString AliEveFileDialog::GetPathESD() const
{
    return TString(fPathEntryESD->GetText());
}

const TString AliEveFileDialog::GetPathESDfriend() const
{
    return TString(fPathEntryESDfriend->GetText());
}

const TString AliEveFileDialog::GetPathAOD() const
{
    return TString(fPathEntryAOD->GetText());
}

const TString AliEveFileDialog::GetPathAODfriend() const
{
    return TString(fPathEntryAODfriend->GetText());
}

const TString AliEveFileDialog::GetPathRaw() const
{
    return TString(fPathEntryRawFile->GetText());
}

const TString AliEveFileDialog::GetUrl() const
{
    return TString(fPathEntryUrl->GetText());
}

const TString AliEveFileDialog::GetCDBStoragePath() const
{
    TString cdbSelected;

    switch(fCDBPathCB->GetSelected()){
    case AliEveUtil::CDB_LOCAL:
        cdbSelected = "local://";
        break;
    case AliEveUtil::CDB_MCFULL:
        cdbSelected = "mcfull://";
        break;
    case AliEveUtil::CDB_MCIDEAL:
        cdbSelected = "mcideal://";
        break;
    case AliEveUtil::CDB_MCRESIDUAL:
        cdbSelected = "mcresidual://";
        break;
    case AliEveUtil::CDB_RAW:
        cdbSelected = "raw://";
        break;
    default:
        cdbSelected= "local://";
    }


    return cdbSelected;
}

void AliEveFileDialog::onBrowseESDFile()
{
    TGFileInfo* fileInfo = new TGFileInfo;

    const char* types[] = { "ALICE ESD file", "*.root",
                            "ROOT Archive", "*.zip",
                            0, 0};

    fileInfo->fFileTypes = types;

    new TGFileDialog(GetParent(), GetMain(), kFDOpen, fileInfo);

    fPathEntryESD->SetText(fileInfo->fFilename);
    
    // look for the other files in the current directory
    
    

    delete fileInfo;

}

void AliEveFileDialog::onBrowseESDfriendFile()
{
    TGFileInfo* fileInfo = new TGFileInfo;

    const char* types[] = { "ALICE ESDfriends file", "*.root",
                            "ROOT Archive", "*.zip",
                            0, 0};

    fileInfo->fFileTypes = types;

    new TGFileDialog(GetParent(), GetMain(), kFDOpen, fileInfo);

    fPathEntryESDfriend->SetText(fileInfo->fFilename);

    delete fileInfo;

}

void AliEveFileDialog::onBrowseAODFile()
{
    TGFileInfo* fileInfo = new TGFileInfo;

    const char* types[] = { "ALICE AOD file", "*.root",
                            "ROOT Archive", "*.zip",
                            0, 0};

    fileInfo->fFileTypes = types;

    new TGFileDialog(GetParent(), GetMain(), kFDOpen, fileInfo);

    fPathEntryAOD->SetText(fileInfo->fFilename);

    delete fileInfo;

}

void AliEveFileDialog::onBrowseAODfriendFile()
{
    TGFileInfo* fileInfo = new TGFileInfo;

    const char* types[] = { "ALICE AODfriend file", "*.root",
                            "ROOT Archive", "*.zip",
                            0, 0};

    fileInfo->fFileTypes = types;

    new TGFileDialog(GetParent(), GetMain(), kFDOpen, fileInfo);

    fPathEntryAODfriend->SetText(fileInfo->fFilename);

    delete fileInfo;

}

void AliEveFileDialog::onBrowseRawFile()
{
    TGFileInfo* fileInfo = new TGFileInfo;

    const char* types[] = { "ALICE RAW file", "*.root",
                            "ROOT Archive", "*.zip",
                            0, 0};

    fileInfo->fFileTypes = types;

    new TGFileDialog(GetParent(), GetMain(), kFDOpen, fileInfo);

    fPathEntryRawFile->SetText(fileInfo->fFilename);

    delete fileInfo;

}

void AliEveFileDialog::onAccept()
{
    fIsAccepted = kTRUE;
    UnmapWindow();
}

void AliEveFileDialog::onReject()
{
    fIsAccepted = kFALSE;
    UnmapWindow();
}

void AliEveFileDialog::setMode(EAliEveFileDialogMode mode)
{
    fMode = mode;
}

void AliEveFileDialog::MapWindow()
{
    if(IsMapped()) return;

    MapSubwindows();
    TGTransientFrame::MapWindow();
    Layout();

    gClient->WaitFor(this);
}

void AliEveFileDialog::UnmapWindow()
{
    TGTransientFrame::UnmapWindow();
    gClient->ResetWaitFor(this);
}

void AliEveFileDialog::MapSubwindows()
{
    TGTransientFrame::MapSubwindows();

    // Show/Hide Widgets according to the current Mode
    if(fMode==kAliEveFDLocal){
        ShowFrame(fESDFilesFrame);
        ShowFrame(fAdvancedOptsButton);
        showAdvancedOpts(fAdvancedOptsButton->IsDown());

        HideFrame(fUrlFrame);
    }
    else{ // remote file
        HideFrame(fESDFilesFrame);
        HideFrame(fAdvancedOptsButton);
        HideFrame(fAdvancedOptsFrame);

        ShowFrame(fUrlFrame);
    }

}

void AliEveFileDialog::showAdvancedOpts(Bool_t shown)
{
    UInt_t w, h;
    if(shown){
        ShowFrame(fAdvancedOptsFrame);
        fAdvancedOptsFrame->MapWindow();
        Layout();

        h = fESDFilesFrame->GetSize().fHeight+fAdvancedOptsFrame->GetSize().fHeight+fAdvancedOptsButton->GetSize().fHeight+fCDBFrame->GetSize().fHeight+fDialogButtonsFrame->GetSize().fHeight;
    }else{
        HideFrame(fAdvancedOptsFrame);
        fAdvancedOptsFrame->UnmapWindow();
        Layout();

        h = fESDFilesFrame->GetSize().fHeight+fAdvancedOptsButton->GetSize().fHeight+fCDBFrame->GetSize().fHeight+fDialogButtonsFrame->GetSize().fHeight;
    }

    TGDimension size = GetSize();
    w = size.fWidth;

    Resize(w,h+30);

}

void AliEveFileDialog::CloseWindow()
{
    onReject();
}
