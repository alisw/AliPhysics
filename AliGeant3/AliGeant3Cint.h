/********************************************************************
* AliGeant3Cint.h
********************************************************************/
#ifdef __CINT__
#error AliGeant3Cint.h/C is only for compilation. Abort cint.
#endif
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define G__ANSIHEADER
#define G__DICTIONARY
#include "G__ci.h"
extern "C" {
extern void G__cpp_setup_tagtableAliGeant3Cint();
extern void G__cpp_setup_inheritanceAliGeant3Cint();
extern void G__cpp_setup_typetableAliGeant3Cint();
extern void G__cpp_setup_memvarAliGeant3Cint();
extern void G__cpp_setup_globalAliGeant3Cint();
extern void G__cpp_setup_memfuncAliGeant3Cint();
extern void G__cpp_setup_funcAliGeant3Cint();
extern void G__set_cpp_environmentAliGeant3Cint();
}


#include "TROOT.h"
#include "TMemberInspector.h"
#include "AliG3Medium.h"
#include "AliG3Material.h"
#include "AliG3Volume.h"
#include "AliGUISliders.h"
#include "AliGuiGeomDialog.h"
#include "AliGuiGeomMain.h"
#include "AliGeant3GeometryGUI.h"
#include "AliNode.h"
#include "AliG3toRoot.h"

#ifndef G__MEMFUNCBODY
#endif

extern G__linked_taginfo G__AliGeant3CintLN_TClass;
extern G__linked_taginfo G__AliGeant3CintLN_TList;
extern G__linked_taginfo G__AliGeant3CintLN_TObjArray;
extern G__linked_taginfo G__AliGeant3CintLN_TObject;
extern G__linked_taginfo G__AliGeant3CintLN_TNamed;
extern G__linked_taginfo G__AliGeant3CintLN_TFolder;
extern G__linked_taginfo G__AliGeant3CintLN_AliG3Medium;
extern G__linked_taginfo G__AliGeant3CintLN_TAttFill;
extern G__linked_taginfo G__AliGeant3CintLN_TMaterial;
extern G__linked_taginfo G__AliGeant3CintLN_AliG3Material;
extern G__linked_taginfo G__AliGeant3CintLN_TGObject;
extern G__linked_taginfo G__AliGeant3CintLN_TGWindow;
extern G__linked_taginfo G__AliGeant3CintLN_TAttLine;
extern G__linked_taginfo G__AliGeant3CintLN_TQObject;
extern G__linked_taginfo G__AliGeant3CintLN_TGFrame;
extern G__linked_taginfo G__AliGeant3CintLN_TGCompositeFrame;
extern G__linked_taginfo G__AliGeant3CintLN_TGLayoutHints;
extern G__linked_taginfo G__AliGeant3CintLN_TGHorizontalFrame;
extern G__linked_taginfo G__AliGeant3CintLN_TGMainFrame;
extern G__linked_taginfo G__AliGeant3CintLN_TGCanvas;
extern G__linked_taginfo G__AliGeant3CintLN_TGListTreeItem;
extern G__linked_taginfo G__AliGeant3CintLN_TGListTree;
extern G__linked_taginfo G__AliGeant3CintLN_TClonesArray;
extern G__linked_taginfo G__AliGeant3CintLN_Decay_t;
extern G__linked_taginfo G__AliGeant3CintLN_Quest_t;
extern G__linked_taginfo G__AliGeant3CintLN_Gcbank_t;
extern G__linked_taginfo G__AliGeant3CintLN_Gclink_t;
extern G__linked_taginfo G__AliGeant3CintLN_Gcflag_t;
extern G__linked_taginfo G__AliGeant3CintLN_Gckine_t;
extern G__linked_taginfo G__AliGeant3CintLN_Gcking_t;
extern G__linked_taginfo G__AliGeant3CintLN_Gckin2_t;
extern G__linked_taginfo G__AliGeant3CintLN_Gckin3_t;
extern G__linked_taginfo G__AliGeant3CintLN_Gcmate_t;
extern G__linked_taginfo G__AliGeant3CintLN_Gctmed_t;
extern G__linked_taginfo G__AliGeant3CintLN_Gctrak_t;
extern G__linked_taginfo G__AliGeant3CintLN_Gcvolu_t;
extern G__linked_taginfo G__AliGeant3CintLN_Gcsets_t;
extern G__linked_taginfo G__AliGeant3CintLN_Gcnum_t;
extern G__linked_taginfo G__AliGeant3CintLN_Gccuts_t;
extern G__linked_taginfo G__AliGeant3CintLN_Gcmulo_t;
extern G__linked_taginfo G__AliGeant3CintLN_Gcphys_t;
extern G__linked_taginfo G__AliGeant3CintLN_Gcphlt_t;
extern G__linked_taginfo G__AliGeant3CintLN_Gcopti_t;
extern G__linked_taginfo G__AliGeant3CintLN_Gctlit_t;
extern G__linked_taginfo G__AliGeant3CintLN_Gcvdma_t;
extern G__linked_taginfo G__AliGeant3CintLN_Gctpol_t;
extern G__linked_taginfo G__AliGeant3CintLN_Ertrio_t;
extern G__linked_taginfo G__AliGeant3CintLN_Eropts_t;
extern G__linked_taginfo G__AliGeant3CintLN_Eroptc_t;
extern G__linked_taginfo G__AliGeant3CintLN_Erwork_t;
extern G__linked_taginfo G__AliGeant3CintLN_TArrayF;
extern G__linked_taginfo G__AliGeant3CintLN_AliG3Volume;
extern G__linked_taginfo G__AliGeant3CintLN_TGTextBuffer;
extern G__linked_taginfo G__AliGeant3CintLN_TGTextEntry;
extern G__linked_taginfo G__AliGeant3CintLN_TGLabel;
extern G__linked_taginfo G__AliGeant3CintLN_TGComboBox;
extern G__linked_taginfo G__AliGeant3CintLN_TGTab;
extern G__linked_taginfo G__AliGeant3CintLN_AliGuiGeomDialog;
extern G__linked_taginfo G__AliGeant3CintLN_TGTextButton;
extern G__linked_taginfo G__AliGeant3CintLN_TGListBox;
extern G__linked_taginfo G__AliGeant3CintLN_TGMenuBar;
extern G__linked_taginfo G__AliGeant3CintLN_TGPopupMenu;
extern G__linked_taginfo G__AliGeant3CintLN_AliNode;
extern G__linked_taginfo G__AliGeant3CintLN_AliGuiGeomMain;
extern G__linked_taginfo G__AliGeant3CintLN_AliGeant3GeometryGUI;
extern G__linked_taginfo G__AliGeant3CintLN_TAtt3D;
extern G__linked_taginfo G__AliGeant3CintLN__x3d_data_;
extern G__linked_taginfo G__AliGeant3CintLN__x3d_sizeof_;
extern G__linked_taginfo G__AliGeant3CintLN_TNode;
extern G__linked_taginfo G__AliGeant3CintLN_TGeometry;
extern G__linked_taginfo G__AliGeant3CintLN_AliG3toRoot;

/* STUB derived class for protected member access */
