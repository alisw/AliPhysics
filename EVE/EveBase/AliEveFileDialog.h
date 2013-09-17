// Author: Mihai Niculescu 2013

/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveFileDialog_H
#define AliEveFileDialog_H

#include <TString.h>
#include <TGFrame.h>

class TGComboBox;
class TGCheckButton;
class TGLabel;
class TGFileInfo;
class TGTextButton;
class TGTextEntry;

//______________________________________________________________________________
// AliEveMainWindow
//
enum EAliEveFileDialogMode
{
    kAliEveFDLocal,
    kAliEveFDRemote
};

class AliEveFileDialog : public TGTransientFrame
{
public:
    AliEveFileDialog(const TGWindow* p = 0, const TGWindow* main = 0, EAliEveFileDialogMode mode= kAliEveFDLocal);
    ~AliEveFileDialog();

    //const TString GetDirectory(); // directory where the ESD resides
    const TString GetPathESD() const;
    const TString GetPathESDfriend() const;
    const TString GetPathAOD() const;
    const TString GetPathAODfriend() const;
    const TString GetPathRaw() const;
    const TString GetUrl() const;
    UInt_t GetMode() const { return fMode; }

    const TString GetCDBStoragePath() const;

    void onBrowseESDFile();
    void onBrowseESDfriendFile();
    void onBrowseAODFile();
    void onBrowseAODfriendFile();
    void onBrowseRawFile();

    Bool_t accepted() const { return fIsAccepted; } // true if clicked the OK button, false otherwise

    void onAccept();// ok button was clicked
    void onReject();// cancel button was clicked

    void setMode(EAliEveFileDialogMode mode);
    void showAdvancedOpts(Bool_t shown=kTRUE);

    void CloseWindow();
    void MapWindow();
    void UnmapWindow();
    void MapSubwindows();
private:
    AliEveFileDialog(const AliEveFileDialog&);              // not implemented
    AliEveFileDialog& operator=(const AliEveFileDialog&);   // not implemented

    TGHorizontalFrame* fESDFilesFrame;
    TGCheckButton* fAdvancedOptsButton;
    
    TGGroupFrame* fAdvancedOptsFrame;
    TGHorizontalFrame* fESDfriendFilesFrame;
    TGHorizontalFrame* fAODFilesFrame;
    TGHorizontalFrame* fAODfriendFilesFrame;
    TGHorizontalFrame* fRawFilesFrame;
    TGHorizontalFrame* fUrlFrame;
    TGHorizontalFrame* fCDBFrame;
    TGHorizontalFrame* fDialogButtonsFrame;
    TGTextEntry* fPathEntryESD;
    TGTextEntry* fPathEntryESDfriend;
    TGTextEntry* fPathEntryAOD;
    TGTextEntry* fPathEntryAODfriend;
    TGTextEntry* fPathEntryRawFile;
    TGTextEntry* fPathEntryUrl; // remote file
    TGComboBox* fCDBPathCB;

    Bool_t fIsAccepted; // true if clicked OK, false otherwise
    EAliEveFileDialogMode fMode;

public:

    ClassDef(AliEveFileDialog, 0)
};

#endif

