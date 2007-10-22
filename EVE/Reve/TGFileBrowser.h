#ifndef ROOT_TGFileBrowser
#define ROOT_TGFileBrowser

#ifndef ROOT_TGFrame
#include "TGFrame.h"
#endif

#ifndef ROOT_TBrowserImp
#include "TBrowserImp.h"
#endif

class TGCanvas;
class TGListTree;
class TGListTreeItem;
class TGPicture;
class TGLabel;
class TGComboBox;
class TGTextEntry;
class TGTextBuffer;
class TGTextView;
class TContextMenu;
class TRegexp;
class TGMenuBar;
class TGPopupMenu;
class TGStatusBar;
class TGPictureButton;
class TString;
class TGNewBrowser;

class TGFileBrowser : public TGMainFrame, public TBrowserImp {

protected:
   // TGNewBrowser      *fNewBrowser;
   TGHorizontalFrame *fTopFrame;
   TGHorizontalFrame *fBotFrame;
   TGCanvas          *fCanvas;
   TGListTree        *fListTree;
   TGListTreeItem    *fListLevel;         // current TGListTree level
   TGListTreeItem    *fCurrentDir;        // 
   TGListTreeItem    *fRootDir;           // 
   TGComboBox        *fDrawOption;        // draw options combobox
   TGComboBox        *fFileType;          // file type combobox
   TContextMenu      *fContextMenu;       // context menu pointer
   const TGPicture   *fRootIcon;
   const TGPicture   *fFileIcon;
   const TGPicture   *fCachedPic;         //
   TString            fCachedPicName;     //
   TRegexp           *fFilter;
   Int_t              fGroupSize;         // total number of items when icon box switched to "global view" mode
   Long_t             fNKeys, fCnt;
   Bool_t             fGrouped;           //
   Bool_t             fShowHidden;

   TGNewBrowser      *fNewBrowser;

   void CreateBrowser(const char *name);

public:
   TGFileBrowser(TBrowser* b=0, const char *name="ROOT Browser", UInt_t w=200, UInt_t h=400);
   TGFileBrowser(TBrowser* b, const char *name, Int_t x, Int_t y, UInt_t w, UInt_t h);
   virtual ~TGFileBrowser();

   virtual void Add(TObject *obj, const char *name = 0, Int_t check = -1);
   virtual void BrowseObj(TObject *obj);
   virtual void RecursiveRemove(TObject *obj);
   virtual void Refresh(Bool_t force = kFALSE);
   virtual void Show() { MapRaised(); }
   Option_t    *GetDrawOption() const;

   TGNewBrowser *GetNewBrowser() const          { return fNewBrowser; }
   void          SetNewBrowser(TGNewBrowser* b) { fNewBrowser = b;    }

   void        AddFSDirectory(const char* entry, const char* path=0);
   void        AddKey(TGListTreeItem *itm, TObject *obj, const char *name = 0);
   void        ApplyFilter(Int_t id);
   void        Chdir(TGListTreeItem *item);
   void        Clicked(TGListTreeItem *item, Int_t btn, Int_t x, Int_t y);
   TString     DirName(TGListTreeItem* item);
   void        DoubleClicked(TGListTreeItem *item, Int_t btn);
   Long_t      XXExecuteDefaultAction(TObject *obj);
   char       *FormatFileInfo(const char *fname, Long64_t size, Long_t modtime);
   void        GetObjPicture(const TGPicture **pic, TObject *obj);
   void        GotoDir(const char *path);
   
   // overridden from TGMainFrame
   void        ReallyDelete();

   ClassDef(TGFileBrowser, 0) // File browser.
};

#endif
