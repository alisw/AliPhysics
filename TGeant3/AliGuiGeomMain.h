#ifndef ALIGUIGEOMMAIN_H
#define ALIGUIGEOMMAIN_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "TClonesArray.h"
#include "TGFrame.h"
#include "TGListTree.h"
#include "TGComboBox.h"


class TGTab;
class TGMenuBar;
class TGPopupMenu;
class TGTextBuffer;
class TGTextEntry;
class TGLabel;
class TGTextButton;
class AliNode;
class TObjArray;
class TFolder;

class AliG3Material;
class AliG3Medium;
class AliGuiGeomDialog;
class AliG3Volume;

class AliGuiGeomMain : public TGMainFrame {
 public:
    AliGuiGeomMain(const TGWindow *p, UInt_t w, UInt_t h);
    virtual ~AliGuiGeomMain();
    // Destroy the main window
    virtual void CloseWindow();
    // Add item to ListTree
    virtual TGListTreeItem *
	AddItem(TObject *obj, TGListTreeItem* parent,
		const char* name,
		const TGPicture* open, const TGPicture* closed);
    // Add Material to ComboBox
    virtual void AddMaterial(AliG3Material *Material, Int_t i);
    // Add Medium to ComboBox
    virtual void AddMedium(AliG3Medium *Medium, Int_t i);
    // Process messages from this window
    virtual Bool_t ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2);
    // Update widgets
    virtual void Update();
    // Update ComboBoxes
    virtual void UpdateCombo();
    virtual void UpdateListBox();
    // Relate objects to ComboEntries
    // Currently ComboBox Entries are strings only, hence we need this construction
    virtual void SetMaterialComboEntries(TClonesArray *entries);
    virtual void SetMediaComboEntries(TClonesArray *entries);
    virtual void AddFoldersRecursively(TFolder* folder=0, TGListTreeItem* parent=NULL);
    virtual void Plot();
private:
    TGTab              *fTab;           // Contains Tab entries: volumes, materials..
    TGCanvas           *fCanvasWindow;  // Canvas window for list tree
    TGCompositeFrame   *fF2, *fF21, *fF3, *fF31, *fF4, *fF5;      // Frames for combos
    TGCompositeFrame   *fF6, *fF61, *fF62, *fF63;                 // Frames for combos
    TGListTree         *fLt;                                      // Volumes list tree
    TGMenuBar          *fMenuBar;                                 // Menu bar: File, Draw Control ...
    TGPopupMenu        *fMenuFile, *fMenuTest, *fMenuHelp;        // Pop-up menus
    TGLayoutHints      *fMenuBarItemLayout, *fMenuBarHelpLayout,  // Lay-out hints
	               *fMenuBarLayout, fLTab;                    // Lay-out hints
    TGLayoutHints      *fL2;                                      // Lay-out hints
    AliGuiGeomDialog   *fDialog;                                  //! no output please
    TGComboBox         *fMaterialCombo;                           // Material  combo box
    TGComboBox         *fMechanismCombo;                          // Mechanism combo box
    TGComboBox         *fMediaCombo, *fParticleCombo;             // Media and particle combo boxes
    TGListBox          *fProcessLB, *fCutsLB;                     // List boxes for cuts and processes
    TClonesArray       *fComboMaterialEntries;                    // List of materials
    TClonesArray       *fComboMediaEntries;                       // List of media
    TGHorizontalFrame  *fHframe[6],*fHframeM[8];                  // sub frames 
    TGTextBuffer       *fTbh[6], *fTbhM[8], *fTbh61, *fTbh62, *fTbh63;  // text frames
    TGTextEntry        *fTeh[6], *fTehM[8], *fTeh61, *fTeh62, *fTeh63;  // text entries
    TGLabel            *fLabel[6], *fLabelM[8], *fSLabel61;             // labels
    TGTextButton       *fPlotButton;                                    // Plot-Button
    Float_t            fEmin;         // minimum energy for de/dx plot
    Float_t            fEmax;         // maximum energy for de/dx plot
    Int_t              fNbins;        // number of bins for de/dx plot
  AliGuiGeomMain(const AliGuiGeomMain &gm) 
    : TGMainFrame((const TGMainFrame&)gm) {}
  virtual AliGuiGeomMain & operator=(const AliGuiGeomMain &) {return *this;}
    

    ClassDef(AliGuiGeomMain,1)  // MainFrame for Geometry Browser
};

R__EXTERN AliG3Material  *gCurrentMaterial;
R__EXTERN AliG3Medium    *gCurrentMedium;

#endif
