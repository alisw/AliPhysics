#ifndef ROOT_TGClassBrowser
#define ROOT_TGClassBrowser

#ifndef ROOT_TGFrame
#include "TGFrame.h"
#endif

class TGCanvas;
class TGListTree;
class TGListTreeItem;
class TGPicture;

class TGClassBrowser : public TGMainFrame {

protected:
   TGCanvas          *fCanvas;
   TGListTree        *fListTree;
   const TGPicture   *fClassIcon;
   const TGPicture   *fMemberIcon;
   const TGPicture   *fMethodIcon;

public:
   TGClassBrowser(const TGWindow *p, UInt_t w, UInt_t h);
   virtual ~TGClassBrowser();

   void     DisplayClass(TGListTreeItem *item, const TString &fname);
   void     DoubleClicked(TGListTreeItem *item, Int_t btn);

   ClassDef(TGClassBrowser, 0) // ROOT Classes browser
};

#endif

