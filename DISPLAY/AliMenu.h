#ifndef ALIMENU_H
#define ALIMENU_H
/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/////////////////////////////////////////////////////////////////////////
// ALICE MENU CLASS                                                    //
// Author: Mayeul   ROUSSELET                                          //
// e-mail: Mayeul.Rousselet@cern.ch                                    //
// Last update:26/08/2003                                              //
/////////////////////////////////////////////////////////////////////////

#include <Rtypes.h>
#include <RQ_OBJECT.h>
#include <TGToolBar.h>

class TGMenuBar;
class TGLayoutHints;
class TGPopupMenu;
class TGButton;

class AliMenu {
  //This class implement both the menu and the toolbar

public:

 AliMenu(TGCompositeFrame *p, UInt_t w, UInt_t h, UInt_t options);
 virtual ~AliMenu();
 
 void 		      		AddPictureButton( const char *fname, const char *tiptext,UInt_t id, UInt_t spacing);//add a picture button to the toolbar
 //slots
 void		      		DoMenu(Int_t id=0);
 void		       		DoToolBar(Int_t id=0);

private:

 TGMenuBar			*fMenuBar; // Menu bar
 TGToolBar		       	*fToolBar; // Tool bar
 TGLayoutHints			*fMenuBarLayout; // Menu bar layout
 TGLayoutHints			*fMenuBarItemLayout; // Menu bar item layout
 TGPopupMenu		       	*fMenuFile; // Menu file
 TGPopupMenu		       	*fMenuOptions; // Menu options
 TGPopupMenu		       	*fMenuHelp; // Menu help
 TGPopupMenu		       	*fMenuView; // Menu view
 TGLayoutHints			*fToolBarLayout; // Tool bar layout
 ToolBarData_t 			*fTBD; // Tool bar data
 TGButton                       *fButton; // Button

 RQ_OBJECT("AliMenu")

 ClassDef(AliMenu,0);
};

#endif
