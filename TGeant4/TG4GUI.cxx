// $Id$
// Category: interfaces
//
// Author: D. Adamova
//
//========================================================
//
//------------TG4GUI.cxx--------------------------------//
//---------Main Window for the AG4 Geometry Browser---//
//
//========================================================= 

#include "TG4GUI.h"
#include "TG4Editor.h"
 
#include <TGListTree.h>
#include <TGCanvas.h>
#include <TGTab.h>
#include <TGMenu.h>
#include <TGCanvas.h>
#include <TApplication.h>
#include <TGMsgBox.h>
#include <TGTextBuffer.h>
 

ClassImp(TG4GUI)

TG4GUI::TG4GUI(const TGWindow* p, UInt_t w, UInt_t h)
    : TGMainFrame(p, w, h)
{
    // Create a main frame 

  fMenuBarLayout = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX, 0, 0, 1, 1);
  fMenuBarItemLayout = new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0);
  fMenuBarHelpLayout = new TGLayoutHints(kLHintsTop | kLHintsRight);

  fMenuFile = new TGPopupMenu(gClient->GetRoot());
  fMenuFile->AddLabel("CloseWindow/Exit");
  fMenuFile->AddEntry("&Close", 1);
  fMenuFile->AddEntry("&Exit", 2);
  
  fMenuTest = new TGPopupMenu(this);
  fMenuTest->AddEntry("&Message", 3);
    
  fMenuHelp = new TGPopupMenu(gClient->GetRoot());
  fMenuHelp->AddEntry("&About", 4);

   fMenuFile->Associate(this);
   fMenuTest->Associate(this);
   fMenuHelp->Associate(this);
 
   fMenuBar = new TGMenuBar(this, 1, 1, kHorizontalFrame);
   fMenuBar->AddPopup("&Quit??", fMenuFile, fMenuBarItemLayout);
   fMenuBar->AddPopup("&Draw Control", fMenuTest, fMenuBarItemLayout);
   fMenuBar->AddPopup("&Report", fMenuHelp, fMenuBarHelpLayout);

   AddFrame(fMenuBar, fMenuBarLayout);
 
//------>Volumes

   fTab = new TGTab(this, 400, 400);
   TGCompositeFrame *tf = fTab->AddTab("Volumes");
   TGLayoutHints *lTab = new TGLayoutHints(kLHintsBottom | kLHintsExpandX |
                                          kLHintsExpandY, 2, 2, 5, 1);
   AddFrame(fTab, lTab);

//------>TGCanvas and a canvas container  
   fCanvasWindow = new TGCanvas(tf, 400, 240);
 
//----->List
   fLt = new TGListTree(fCanvasWindow->GetViewPort(), 10, 10, kHorizontalFrame,
                        fgWhitePixel);
   fLt->Associate(this);
   fCanvasWindow->SetContainer(fLt);

 
   TGLayoutHints *lo = new TGLayoutHints(kLHintsExpandX | kLHintsExpandY);
   tf->AddFrame(fCanvasWindow, lo);
 
//----------------------------------------------------------------------------- 
 
//----->Window name and final mapping

   SetWindowName("AG4 Geometry Browser");
   MapSubwindows();
   
 
   Layout();
   MapWindow();    
 

}

TG4GUI::~TG4GUI()
{

  // Delete created widgets.
    
   delete fCanvasWindow;

   delete fMenuBarLayout;
   delete fMenuBarItemLayout;
   delete fMenuBarHelpLayout;

   delete fMenuFile;
   delete fMenuTest;
   delete fMenuHelp;

}

TGListTreeItem*  TG4GUI::
AddItem(TObject* obj, TGListTreeItem* parent, const char* name, 
                       const TGPicture* open, const TGPicture* closed)
{
//----->Add item to the list tree
    return fLt->AddItem(parent, name, obj, open, closed);
}

void TG4GUI::CloseWindow()
{
   //----->Close the window & exit root    

   TGMainFrame::CloseWindow();
   gApplication->Terminate(0);
}

Bool_t TG4GUI::ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2)
{

/*************************************************************/
//----->Message Box

   int  retval;
   EMsgBoxIcon icontype = kMBIconExclamation;
   int buttons = kMBOk;
   TGTextBuffer* fTbtitle = new TGTextBuffer(100);
   TGTextBuffer* fTbmsg   = new TGTextBuffer(100);
   fTbtitle->AddText(0, "MsgBox");
   fTbmsg->AddText(0, "Volumes drawing & materials coming soon!");

/*************************************************************/

//----->Editor window text
const char *editortxt =
"This is a to-be AG4 Geometry Browser. \n"
"The volumes drawing and materials coming soon.\n" ;


/***********************************************************/

//----->Process messages to widgets
    switch (GET_MSG(msg)) {

//---->case Handle Popup menus    
    case kC_COMMAND:
        switch (GET_SUBMSG(msg)) {
            case kCM_MENU:
               switch (parm1) {

                   case 1:
                     printf("CLOSING MAIN WINDOW\n");
		     TGMainFrame::CloseWindow();
                     break;  

                  case 2:
                     CloseWindow();   //-->and exit root  
                     break;
		     
		   case 3:
		      //-->popup Message 
		      new TGMsgBox(fClient->GetRoot(), this,
                                  fTbtitle->GetString(), fTbmsg->GetString(),
                                  icontype, buttons, &retval);
		      break;
		      
		  case 4:
		     //-->editor window
		     {
		      TG4Editor *ed = new TG4Editor(this, 400, 150);
                      ed->LoadBuffer(editortxt);
                      ed->Popup(); 
		     };  

                   default:
                     break;
               }
            default:
               break;
         }
	 break;	   

//----->case Handle volumes listing
    case kC_LISTTREE:
	switch (GET_SUBMSG(msg)) {

//----->Cases to Handle mouse click
   //-->case 1 
	case kCT_ITEMCLICK:    
	    break;
   //-->case 2	    
	case kCT_ITEMDBLCLICK:
	    if (parm1 == kButton1) {
		if (fLt->GetSelected() != 0) {
		    gClient->NeedRedraw(fLt);
		}
	    }
	    break;
   //-->default for GET_SUBMSG	    
	default:
	    break;
	}
	break;
//---->default for GET_MSG	
    default:
	break;
    }
    return kTRUE;
}
