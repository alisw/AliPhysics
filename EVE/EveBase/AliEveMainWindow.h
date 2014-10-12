// Author: Mihai Niculescu 2013

/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveMainWindow_H
#define AliEveMainWindow_H

#include <TGFrame.h>

class TEveManager;
class TGMenuBar;
class TGPicturePool;
class TGPopupMenu;
class TGToolBar;
class AliEveMainWindow;
class AliEveFileDialog;

// AliEveMainWindow
//

class AliEveMainWindow : public TGMainFrame
{
public:
    AliEveMainWindow(const char* title, UInt_t width=800, UInt_t height=600);
    ~AliEveMainWindow();
// Items IDs for MenuBar
enum MENU_ITEM_IDS
{   // File Menu
    MENU_FILE_OPEN=0,
    MENU_FILE_OPEN_URL,
    MENU_FILE_OPEN_CONNECTION,
    MENU_FILE_EXPORT_VIEWS,
    MENU_FILE_EXIT,
    // Edit Menu
    MENU_EDIT_UNDO,
    MENU_EDIT_REDO,
    MENU_EDIT_CUT,
    MENU_EDIT_COPY,
    MENU_EDIT_PASTE,
    MENU_EDIT_DELETE,
    MENU_EDIT_PROP,
    // View Menu
    MENU_VIEW_TOOLBAR_MAIN,
    MENU_VIEW_TOOLBAR_NAV, // event navigation toolbar
    MENU_VIEW_TOOLBAR_PROP, // Properties Sidebar - event info, objects list, etc...
    MENU_VIEW_TOOLBAR_HIST, // History sidebar
    MENU_VIEW_RELOAD,
    MENU_VIEW_ZOOM_IN,
    MENU_VIEW_ZOOM_OUT,
    MENU_VIEW_ZOOM_RESET,
    // Go Menu
    MENU_GO_NEXT_EVENT,
    MENU_GO_PREV_EVENT,
    MENU_GO_FIRST_EVENT,
    MENU_GO_LAST_EVENT,
    MENU_GO_PLAY,
    // History Menu
    MENU_HIST_SHOW_ALL,
    MENU_HIST_CLEAR_RECENT,
    // Tools Menu
    MENU_TOOLS_MACROS,
    MENU_TOOLS_QA,
    // Help Menu
    MENU_HELP_CONTENTS,
    MENU_HELP_ABOUT
};

public: // SLOTS
    void onMenuFileItem(UInt_t id);
    void onMenuEditItem(UInt_t id);
    void onMenuViewItem(UInt_t id);
    void onMenuGoItem(UInt_t id);

protected:
    void setupMenus();
    void setupToolbars();

    void loadFiles();

private:
    AliEveMainWindow(const AliEveMainWindow& other);// Not implemented
    AliEveMainWindow& operator=(const AliEveMainWindow& other);


    // Menubar
    TGMenuBar* fMenuBar;
    TGPopupMenu *fMenuFile;
    TGPopupMenu *fMenuEdit;
    TGPopupMenu *fMenuView;
    TGPopupMenu *fMenuViewToolbars;
    TGPopupMenu *fMenuViewSidebars;
    TGPopupMenu *fMenuGo;
    TGPopupMenu *fMenuTools;
    TGPopupMenu *fMenuHelp;

    // Toolbar
    TGToolBar *fToolBar;

    TGPicturePool* fPicturePool;

//    TEveManager* fEve;
    AliEveFileDialog* fFileDialog;

    ClassDef(AliEveMainWindow, 0); // Short description.
};

#endif

