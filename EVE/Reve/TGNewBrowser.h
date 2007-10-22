#ifndef ROOT_TGNewBrowser
#define ROOT_TGNewBrowser

#ifndef ROOT_TGFrame
#include "TGFrame.h"
#endif

//#ifndef ROOT_TBrowserImp
//#include "TBrowserImp.h"
//#endif

class TGLayoutHints;
class TGTab;
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
class TGVSplitter;
class TGHSplitter;

class TGFileBrowser;

class TGNewBrowser : public TGMainFrame { //, public TBrowserImp {

protected:

   TGLayoutHints     *fLH0, *fLH1, *fLH2, *fLH3, *fLH4;
   TGLayoutHints     *fLH5, *fLH6, *fLH7;
   TGTab             *fTabLeft;
   TGTab             *fTabRight;
   TGTab             *fTabBottom;
   TGTab             *fEditTab;
   TGVerticalFrame   *fVf;
   TGHorizontalFrame *fHf;
   TGHorizontalFrame *fH1;
   TGHorizontalFrame *fH2;
   TGVerticalFrame   *fV1;
   TGVerticalFrame   *fV2;
   TGVSplitter       *fVSplitter;
   TGHSplitter       *fHSplitter;
   TGCompositeFrame  *fEditFrame;
   TGHorizontalFrame *fTopMenuFrame;
   TGHorizontalFrame *fPreMenuFrame;
   TGHorizontalFrame *fMenuFrame;
   TGHorizontalFrame *fToolbarFrame;
   TGMenuBar         *fMenuBar;
   TGPopupMenu       *fMenuFile;
   TGCompositeFrame  *fActMenuBar;
   TGStatusBar       *fStatusBar;
   Int_t              fNbTab[3];
   Int_t              fCrTab[3];
   Int_t              fPid;               // current process id

public:
   enum EInsertPosition {
      kLeft, kRight, kBottom
   };

   TGNewBrowser(const char *name = "ROOT Browser", UInt_t width = 800, UInt_t height = 500, Bool_t initshow=kTRUE);
   TGNewBrowser(const char *name, Int_t x, Int_t y, UInt_t width, UInt_t height, Bool_t initshow=kTRUE);
   virtual ~TGNewBrowser();

   TGFileBrowser*    MakeFileBrowser();
   void              InitPlugins();

   void              CreateBrowser(const char *name);
   void              CloseWindow();
   void              DoTab(Int_t id);
   TGFrame          *GetActFrame() const { return (TGFrame *)fEditFrame; }
   TGStatusBar      *GetStatusBar() const { return fStatusBar; }
   TGTab            *GetTabLeft() const { return fTabLeft; }
   TGTab            *GetTabRight() const { return fTabRight; }
   TGTab            *GetTabBottom() const { return fTabBottom; }
   TGTab            *GetTab(Int_t pos) const;
   void              SetTab(Int_t pos = kRight, Int_t subpos = -1);
   void              SetTabTitle(const char *title, Int_t pos = kRight, Int_t subpos = -1);
   void              HandleMenu(Int_t id);
   using             TGCompositeFrame::RemoveFrame;
   void              RecursiveReparent(TGPopupMenu *popup);
   void              RemoveFrame(Int_t pos, Int_t subpos);
   void              ShowMenu(TGCompositeFrame *menu);
   TGCompositeFrame *StartEmbedding(Int_t pos = kRight, Int_t subpos = -1);
   void              StopEmbedding(TGLayoutHints *layout=0);
   void              SwitchMenus(TGCompositeFrame *from);

   virtual void      ExecPlugin(const char *fname, Int_t pos = kRight, Int_t subpos = -1);
   virtual Bool_t    HandleKey(Event_t *event);

   // static TBrowserImp *NewBrowser(TBrowser *b = 0, const char *title = "ROOT Browser", UInt_t width = 800, UInt_t height = 500);
   // static TBrowserImp *NewBrowser(TBrowser *b, const char *title, Int_t x, Int_t y, UInt_t width, UInt_t height);

   // overridden from TGMainFrame
   // void              ReallyDelete();

   ClassDef(TGNewBrowser, 0)
};

#endif
