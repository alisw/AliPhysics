#ifndef ROOT_guitest
#define ROOT_guitest
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TROOT.h>
#include <TVirtualX.h>
#include <TGListBox.h>
#include <TGListTree.h>
#include <TGClient.h>
#include <TGFrame.h>
#include <TGIcon.h>
#include <TGLabel.h>
#include <TGButton.h>
#include <TGTextEntry.h>
#include <TGMsgBox.h>
#include <TGMenu.h>
#include <TGCanvas.h>
#include <TGComboBox.h>
#include <TGTab.h>
#include <TGSlider.h>
#include <TGDoubleSlider.h>
#include <TGFileDialog.h>
#include <TRootEmbeddedCanvas.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TClonesArray.h>
#include <TGeant3.h>

class AliGuiGeomDialog;
class AliGUISliders;
class AliGuiGeomMain;
class AliDrawVolume;


class AliGeant3GeometryGUI : public TObject {
 private:
    AliGuiGeomMain *fPanel;  // the main gui panel
    Int_t          fNstack;  // number of volumes
    TClonesArray   *fVolumes;
    Int_t*    fZlq;
    Float_t*  fZq;    
    Gclink_t* fGclink;
    
 public:
    AliGeant3GeometryGUI();
    void  ReadGeometryTree();
 private:
    virtual AliDrawVolume* Volume(Int_t id)
	{return (AliDrawVolume *) (fVolumes->UncheckedAt(id));}
	
    Int_t NChildren(Int_t idvol);
    Int_t Child(Int_t idvol, Int_t idc);
    ClassDef(AliGeant3GeometryGUI,1)  // GUI for Geant3 geometry visualisation
};


class AliGuiGeomMain : public TGMainFrame {

private:
   TGCanvas           *fCanvasWindow;
   TGCompositeFrame   *fContainer;
   TGListTree         *fLt;
   TGMenuBar          *fMenuBar;
   TGPopupMenu        *fMenuFile, *fMenuTest, *fMenuHelp;
   TGLayoutHints      *fMenuBarItemLayout, *fMenuBarHelpLayout, *fMenuBarLayout;
   AliGuiGeomDialog   *fDialog; //! no output please
   
public:
   AliGuiGeomMain(const TGWindow *p, UInt_t w, UInt_t h);
   virtual ~AliGuiGeomMain();

   virtual void CloseWindow();
   virtual TGListTreeItem *
       AddItem(TObject *, TGListTreeItem* ,
	       const char*, const TGPicture*, const TGPicture*);
   virtual Bool_t ProcessMessage(Long_t msg, Long_t parm1, Long_t);
   virtual void Update();
   ClassDef(AliGuiGeomMain,1)  // MainFrame for Geometry Browser
};


class AliGuiGeomDialog : public TGTransientFrame {

private:
    AliGUISliders       *fF1;
    TGCompositeFrame    *fFrame1, *fF2, *fF3, *fF4;
    TGButton            *fOkButton, *fCancelButton;
    TGButton            *fChk1, *fChk2, *fChk3;
    TGComboBox          *fCombo, *fCombo2;
    TGLabel             *fLabel1, *fLabel2;
    TGTab               *fTab;
    TGLayoutHints       *fL1, *fL2, *fL3, *fL4, *fBly, *fBfly1;
    TGHorizontalFrame   *fHSframe1, *fHSframe2, *fHSframe3;
    TGTextBuffer        *fTbh11, *fTbh12, *fTbh21, *fTbh22, *fTbh31, *fTbh32;
    TGTextEntry         *fTeh11, *fTeh12, *fTeh21, *fTeh22, *fTeh31, *fTeh32;
    TGDoubleHSlider     *fDslider1, *fDslider2, *fDslider3;
    TGLabel             *fSLabel1,  *fSLabel2,  *fSLabel3;
public:
   AliGuiGeomDialog(const TGWindow *p, const TGWindow *main, UInt_t w, UInt_t h,
               UInt_t options = kMainFrame | kVerticalFrame);
   virtual ~AliGuiGeomDialog();

   virtual void CloseWindow();
   virtual Bool_t ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2);
   virtual void Update();
//   ClassDef(AliGuiGeomDialog,1)  // Just a test 
};

class AliGUISliders : public  TGCompositeFrame {

private:
//
    TGHorizontalFrame *fHframe[8];
    TGLayoutHints     *fBly, *fBfly1;
    TGHSlider         *fHslider[8];
    TGTextEntry       *fTeh[8];
    TGTextBuffer      *fTbh[8];
    TGLabel           *fLabel[8];
    Text_t            fLabelText[8];
public:
   AliGUISliders(const TGWindow *p, const TGWindow *main, UInt_t w, UInt_t h);
   virtual ~AliGUISliders();

   virtual void CloseWindow();
   virtual Bool_t ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2);
   virtual void Update();
      
   //   ClassDef(AliGUISliders,1)  // Just a test 
};

class AliDrawVolume : public TObject 
{
public:
    AliDrawVolume(char* name);
    virtual ~AliDrawVolume(){;}
    virtual void    Draw(Option_t * =0);
    virtual void    DrawSpec();
    virtual char*   Name();
    virtual void    SetParam(Int_t, Float_t);
    virtual Float_t GetParam(Int_t);

    virtual void  SetIdVolume(Int_t id) {fIdVolume = id;}
    virtual void  SetIdCopy(Int_t id)   {fIdCopy = id;}    
    virtual Int_t GetIdVolume()         {return fIdVolume;}
    virtual Int_t GetIdCopy()           {return fIdCopy;}
    virtual void  AddCopy()             {fIdCopy ++;}
    virtual void  SetItem(TGListTreeItem *item) {fItem = item;}


    virtual TGListTreeItem* GetItem() {return fItem;}
	    
private:
    char*   fName;        // name of the volume 
    Float_t fTheta;       // theta-angle for drawing
    Float_t fPhi;         // phi-angle   for drawing
    Float_t fPsi;         // psi-angle   for drawing 
    Float_t fU;           // u-position
    Float_t fV;           // v-position
    Float_t fUscale;      // u-scaling factor
    Float_t fVscale;      // v-scaling factor
    Bool_t  fHide;        // hide flag
    Bool_t  fShadow;      // shadow flag
    Int_t   fFill;        // fill option 1-6
    Int_t   fSeen;        // seen option -2 - 1
    Bool_t  fClip;        // clipping flag
    Float_t fClipXmin;    // clip box range xmin
    Float_t fClipXmax;    // clip box range xmax
    Float_t fClipYmin;    // clip box range ymin
    Float_t fClipYmax;    // clip box range ymax
    Float_t fClipZmin;    // clip box range zmin
    Float_t fClipZmax;    // clip box range zmax
    Int_t   fIdVolume;    // geant volume id
    Int_t   fIdCopy;      // copy flag
    TGListTreeItem        *fItem;
    ClassDef(AliDrawVolume,1) // Volume Object for Drawing 
};
#endif







