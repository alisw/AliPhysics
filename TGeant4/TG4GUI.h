// $Id$
// Category: interfaces
//
// Author: D. Adamova
//
//======================================================
//
//------------TG4GUI.h--------------------------------//
//---------Main Window for the AG4 Geometry Browser---//
//
//=======================================================

#ifndef TG4_GUI_H
#define TG4_GUI_H

#include <TGFrame.h>

class TObject;
class TGListTreeItem;
class TGPicture;
class TGListTree;
class TGTab;
class TGMenuBar;
class TGPopupMenu;
class TGCanvas;

class TG4GUI : public TGMainFrame {

public:   
    
    TG4GUI(const TGWindow *p, UInt_t w, UInt_t h);
    virtual ~TG4GUI();
 
    virtual TGListTreeItem*
        AddItem(TObject* obj, TGListTreeItem* parent,const char* name,
                const TGPicture* open, const TGPicture* closed);    
 
     virtual void CloseWindow();
 
     virtual Bool_t ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2);    
      
private:
    TGTab*      fTab;           // Contains Tab entries
    TGCanvas*   fCanvasWindow;  // Canvas window for the list tree
    TGListTree* fLt;            // Volumes list tree 
   
    TGMenuBar*          fMenuBar;          // main menu bar   
    TGPopupMenu*        fMenuFile;         // popup for window manipulations
    TGPopupMenu*        fMenuTest;         // popup for test messages
    TGPopupMenu*        fMenuHelp;         // popup for help messages
    TGLayoutHints*      fMenuBarItemLayout;// layout left
    TGLayoutHints*      fMenuBarHelpLayout;// layout right
    TGLayoutHints*      fMenuBarLayout;    // main bar layout 
  
   TG4GUI(const TG4GUI& gm) 
    : TGMainFrame( (const TGMainFrame&) gm) {}
  virtual TG4GUI& operator=(const TG4GUI& gm) {return *this;}
    

    ClassDef(TG4GUI,0)    
  };
  
#endif
