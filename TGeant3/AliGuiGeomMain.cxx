/* *************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/*
$Log$
Revision 1.8  2000/07/12 08:56:32  fca
Coding convention correction and warning removal

Revision 1.7  2000/06/28 21:27:45  morsch
Most coding rule violations corrected.
Still to do: Split the file (on file per class) ? Avoid the global variables.
Copy constructors and assignment operators (dummy ?)

Revision 1.6  2000/04/14 11:07:46  morsch
Correct volume to medium assignment in case several media are asigned to the
same material.

Revision 1.5  2000/03/20 15:11:03  fca
Mods to make the code compile on HP

Revision 1.4  2000/01/18 16:12:08  morsch
Bug in calculation of number of volume divisions and number of positionings corrected
Browser for Material and Media properties added

Revision 1.3  1999/11/14 14:31:14  fca
Correct small error and remove compilation warnings on HP

Revision 1.2  1999/11/10 16:53:35  fca
The new geometry viewer from A.Morsch

*/

/* 
 *  Version: 0
 *  Written by Andreas Morsch
 *  
 * 
 *
 * For questions critics and suggestions to this part of the code
 * contact andreas.morsch@cern.ch
 * 
 **************************************************************************/

#include <stdlib.h>

#include "TGMsgBox.h"
#include "TGMenu.h"
#include "TGTab.h"
#include "TGFrame.h"
#include "TGTextBuffer.h"
#include "TGTextEntry.h"
#include "TGLabel.h"
#include "TGButton.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TApplication.h"
#include "TGFileDialog.h"

#include "AliGuiGeomMain.h"
#include "AliGUIMaterial.h"
#include "AliGuiGeomDialog.h"
#include "AliDrawVolume.h"
#include "AliGUIMedium.h"
#include "TGeant3.h"

ClassImp(AliGuiGeomMain)

static Int_t           gCurrentParticle = 1;
static Int_t           gCurrentProcess  = 1;

 const Text_t* kLabelTextP[19]  = 
{"PAIR  ", "COMP  ", "PHOT  ", "PFIS  ", "DRAY  ", "ANNI  ", "BREM  ", 
 "HADR  ", "MUNU  ", "DCAY  ", "LOSS  ", "MULS  ", "GHCOR1", "BIRK1 ", 
 "BIRK2 ", "BIRK3 ", "LABS  ", "SYNC  ", "STRA  "};


 const Text_t* kLabelTextC[10]  = 
 {"CUTGAM", "CUTELE", "CUTNEU", "CUTHAD", "CUTMUO", "BCUTE", "BCUTM",
  "DCUTE ", "DCUTM ", "PPCUTM"};

const Text_t* kLabelTextPart[24]  = 
{"Photon", "Positron", "Electron", "Neutrino", "Muon+", "Muon-", 
 "Pi0", "Pi+", "Pi-", "Kaon_L", "Kaon+", "Kaon-", "Neutron", "Proton", 
 "Anti Proton", "Kaon_S", "Eta", "Lambda", "Sigma+", "Sigma0", "Sigma-",
 "Xi0", "Xi-", "Omega-"};

const Text_t* kLabelTextMechanism[24]  = 
{"HADF", "INEF", "ELAF", "FISF", "CAPF",
 "HADG", "INEG", "ELAG", "FISG", "CAPG",
 "LOSS", "PHOT", "ANNI", "COMP", "BREM",
 "PAIR", "DRAY", "PFIS", "RAYL", "HADG",
 "MUNU", "RANG", "STEP", "MUON"};




enum ETestCommandIdentifiers {
   kFileOpen,
   kFileSave,
   kFileSaveAs,
   kFileExit,

   kTestDlg,

   kHelpContents,
   kHelpSearch,
   kHelpAbout,


   kVId1,
   kHId1,
   kVId2,
   kHId2,

   kVSId1,
   kHSId1,
   kVSId2,
   kHSId2
};


Int_t mbButtonId[9] = { kMBYes, kMBNo, kMBOk, kMBApply,
                          kMBRetry, kMBIgnore, kMBCancel,
                          kMBClose, kMBDismiss };

EMsgBoxIcon mbIcon[4] = { kMBIconStop, kMBIconQuestion,
                           kMBIconExclamation, kMBIconAsterisk };

const char *kFileTypes[] = { "All files",     "*",
                            "ROOT files",    "*.root",
                            "ROOT macros",   "*.C",
                            0,               0 };




TGListTreeItem*  AliGuiGeomMain::
AddItem(TObject * obj, TGListTreeItem *parent, const char* name, const TGPicture *open, const TGPicture *closed)
{
// Add item to the list tree
    return fLt->AddItem(parent, name, obj, open, closed);
}

AliGuiGeomMain::AliGuiGeomMain(const TGWindow *p, UInt_t w, UInt_t h)
      : TGMainFrame(p, w, h)
{
    // Create test main frame. A TGMainFrame is a top level window.
    // Create menubar and popup menus. The hint objects are used to place
    // and group the different menu widgets with respect to eachother.
    
    fDialog=0;
    fMenuBarLayout = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX,
				       0, 0, 1, 1);
    fMenuBarItemLayout = new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0);
    fMenuBarHelpLayout = new TGLayoutHints(kLHintsTop | kLHintsRight);
    
    fMenuFile = new TGPopupMenu(gClient->GetRoot());
    fMenuFile->AddEntry("&Open...", kFileOpen);
    fMenuFile->AddEntry("&Save", kFileSave);
    fMenuFile->AddEntry("S&ave as...", kFileSaveAs);
    fMenuFile->AddEntry("&Close", -1);
    fMenuFile->AddSeparator();
    fMenuFile->AddEntry("E&xit", kFileExit);
    
    fMenuFile->DisableEntry(kFileSaveAs);
    fMenuFile->DisableEntry(kFileOpen);
    fMenuFile->DisableEntry(kFileSave);
    


   fMenuTest = new TGPopupMenu(this);
   fMenuTest->AddLabel("Draw");
   fMenuTest->AddSeparator();
   fMenuTest->AddEntry("&Volume Draw Control", kTestDlg);
   fMenuTest->AddSeparator();
   fMenuTest->Associate(this);
   
   
   fMenuHelp = new TGPopupMenu(gClient->GetRoot());
   fMenuHelp->AddEntry("&Contents", kHelpContents);
   fMenuHelp->AddEntry("&Search...", kHelpSearch);
   fMenuHelp->AddSeparator();
   fMenuHelp->AddEntry("&About", kHelpAbout);

   fMenuFile->DisableEntry(kHelpContents);
   fMenuFile->DisableEntry(kHelpSearch);
   fMenuFile->DisableEntry(kHelpAbout);
   // Menu button messages are handled by the main frame (i.e. "this")
   // ProcessMessage() method.
   fMenuFile->Associate(this);
   fMenuTest->Associate(this);
   fMenuHelp->Associate(this);

   fMenuBar = new TGMenuBar(this, 1, 1, kHorizontalFrame);
   fMenuBar->AddPopup("&File", fMenuFile, fMenuBarItemLayout);
   fMenuBar->AddPopup("&Draw Control", fMenuTest, fMenuBarItemLayout);
   fMenuBar->AddPopup("&Help", fMenuHelp, fMenuBarHelpLayout);

   AddFrame(fMenuBar, fMenuBarLayout);

// 
// Volumes
//
   fTab = new TGTab(this, 400, 400);
   TGCompositeFrame *tf = fTab->AddTab("Volumes");
   TGLayoutHints *lTab = new TGLayoutHints(kLHintsBottom | kLHintsExpandX |
                                          kLHintsExpandY, 2, 2, 5, 1);
   AddFrame(fTab, lTab);

// Create TGCanvas and a canvas container which uses a tile layout manager
   fCanvasWindow = new TGCanvas(tf, 400, 240);
// Create TreeList
   fLt = new TGListTree(fCanvasWindow->GetViewPort(), 10, 10, kHorizontalFrame,
                        fgWhitePixel);
   fLt->Associate(this);
   fCanvasWindow->SetContainer(fLt);

    
   TGLayoutHints *lo = new TGLayoutHints(kLHintsExpandX | kLHintsExpandY);
   tf->AddFrame(fCanvasWindow, lo);
//
// Materials
//
   tf = fTab->AddTab("Materials");
   fF2 = new TGCompositeFrame(tf, 60, 20, kHorizontalFrame);
   fL2 = new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 5);
// ComboBox for materials
   fMaterialCombo = new TGComboBox(fF2, 1);
   fF2->AddFrame(fMaterialCombo, fL2);
// 
// text labels with material properties 
//
   Text_t* labelText[6]  = 
   {"Material Number  ", 
    "Atomic Weight    ",
    "Atomic Number    ", 
    "Density          ",
    "Radiation Length ", 
    "Absorption Length"};

   TGLayoutHints* bly   = 
       new TGLayoutHints(kLHintsTop | kLHintsExpandY, 5, 5, 5, 5);
   TGLayoutHints* bFly1 = 
       new TGLayoutHints(kLHintsLeft | kLHintsExpandX );
   fF21 = new TGCompositeFrame(fF2, 60, 20, kVerticalFrame);
   fF2->AddFrame(fF21,fL2);
   { //Begin local scope for i
     for (Int_t i=0; i<6; i++) {
       Int_t idT=i+1;       
       fHframe[i] = new TGHorizontalFrame(fF21, 400, 100, kFixedWidth);
       fF21->AddFrame(fHframe[i], bly);
       fTbh[i] = new TGTextBuffer(10);
       fTeh[i] = new TGTextEntry(fHframe[i], fTbh[i],idT);
       fTeh[i]->Associate(this);
       fLabel[i] = new TGLabel(fHframe[i], labelText[i]);
       fHframe[i]->AddFrame(fLabel[i], bFly1);
       fHframe[i]->AddFrame(fTeh[i], bFly1);
     }
   } //End local scope for i
   tf->AddFrame(fF2, fL2);
   fMaterialCombo->Resize(200, 20);
   fMaterialCombo->Associate(this);
   

// Media Combo
//   
   tf = fTab->AddTab("Media");
   fF3 = new TGCompositeFrame(tf, 60, 20, kHorizontalFrame);
// ComboBox for tracking media
   fMediaCombo = new TGComboBox(fF3, 2);
   fF3->AddFrame(fMediaCombo, fL2);
// 
// text labels with material properties 
//
   Text_t* labelTextM[8]  = 
   {"Sensitivity Flag      ", 
    "Magnetic Field Flag   ",
    "Maximum Field         ", 
    "Max. Ang. Deviation   ",
    "Maximum Step          ", 
    "Max. Frac. Energy Loss",
    "Crossing Precission   ",
    "Minimum Step Size     "};

   fF31 = new TGCompositeFrame(fF3, 60, 20, kVerticalFrame);
   fF3->AddFrame(fF31,fL2);
   { //Begin local scope for i
     for (Int_t i=0; i<8; i++) {
       Int_t idT=i+1;       
       fHframeM[i] = new TGHorizontalFrame(fF31, 400, 100, kFixedWidth);
       fF31->AddFrame(fHframeM[i], bly);
       fTbhM[i] = new TGTextBuffer(10);
       fTehM[i] = new TGTextEntry(fHframeM[i], fTbhM[i],idT);
       fTehM[i]->Associate(this);
       fLabelM[i] = new TGLabel(fHframeM[i], labelTextM[i]);
       fHframeM[i]->AddFrame(fLabelM[i], bFly1);
       fHframeM[i]->AddFrame(fTehM[i], bFly1);
     }
   } //End local scope for i
   tf->AddFrame(fF3, fL2);
   fMediaCombo->Resize(200, 20);
   fMediaCombo->Associate(this);
//
// Processes
   tf = fTab->AddTab("Processes");
   fF4 = new TGCompositeFrame(tf, 60, 20, kHorizontalFrame);
   fProcessLB = new TGListBox(fF4, 1);
   fF4->AddFrame(fProcessLB, fL2);
   tf->AddFrame(fF4, fL2);
   fProcessLB->Resize(150, 350);
   fProcessLB->Associate(this);
   
//
// Cuts
   tf = fTab->AddTab("Cuts");
   fF5 = new TGCompositeFrame(tf, 60, 20, kHorizontalFrame);
   fCutsLB = new TGListBox(fF5, 1);
   fF5->AddFrame(fCutsLB, fL2);
   tf->AddFrame(fF5, fL2);
   fCutsLB->Resize(150, 350);

//
// de/dx and cross-sections
   tf = fTab->AddTab("DE/DX and X-Sections");
   fF6 = new TGCompositeFrame(tf, 60, 20, kHorizontalFrame);

// ComboBox for particles
   fF61 = new TGCompositeFrame(fF6, 60, 20, kVerticalFrame);
   fF6->AddFrame(fF61, fL2);
   fParticleCombo = new TGComboBox(fF61, 3);
   fF61->AddFrame(fParticleCombo, fL2);
   { //Begin local scope for i
     for (Int_t i = 0; i < 24; i++) {
       char tmp[20];
       sprintf(tmp, "%s", kLabelTextPart[i]);
       fParticleCombo->AddEntry(tmp, i+1);
     }
   } //End local scope for i
   fParticleCombo->Select(1);
   fParticleCombo->Resize(100, 20);
   fParticleCombo->Associate(this);

// ComboBox for mechanisms
   fF63 = new TGCompositeFrame(fF6, 60, 20, kVerticalFrame);
   fF6->AddFrame(fF63, fL2);
   fMechanismCombo = new TGComboBox(fF63, 4);
   fF63->AddFrame(fMechanismCombo, fL2);
   { //Begin local scope for i
     for (Int_t i = 0; i < 24; i++) {
       char tmp[20];
       sprintf(tmp, "%s", kLabelTextMechanism[i]);
       fMechanismCombo->AddEntry(tmp, i+1);
     }
   } //End local scope for i
   fMechanismCombo->Select(1);
   fMechanismCombo->Resize(100, 20);
   fMechanismCombo->Associate(this);

//
// Energy Range
//
   fTbh61 = new TGTextBuffer(10);
   fTeh61 = new TGTextEntry(fF61, fTbh61,10);
   fTbh61->AddText(0, "  100");
   fTeh61->Associate(this);

   fTbh62 = new TGTextBuffer(10);
   fTeh62 = new TGTextEntry(fF61, fTbh62,11);
   fTbh62->AddText(0, "0.001");
   fTeh62->Associate(this);
    
   fTbh63 = new TGTextBuffer(10);
   fTeh63 = new TGTextEntry(fF61, fTbh63,12);
   fTbh63->AddText(0, "10.");
   fTeh63->Associate(this);

   fEmin=0.001;
   fEmax=10.;
   fNbins=100;
   
   fSLabel61 = new TGLabel(fF61, "Nbins-Emin-Emax");
   bFly1 = new TGLayoutHints(kLHintsLeft | kLHintsExpandX );
   fF61->AddFrame(fSLabel61, bFly1);
   fF61->AddFrame(fTeh61, bFly1);
   fF61->AddFrame(fTeh62, bFly1);
   fF61->AddFrame(fTeh63, bFly1);
//
// Plot Button
   fF62 = new TGCompositeFrame(fF6, 60, 20, kHorizontalFrame);
   fF6->AddFrame(fF62, fL2);
   fPlotButton = new TGTextButton(fF62, "Plot", 1);
   fPlotButton -> Associate(this);
   fF62->AddFrame(fPlotButton);

   tf->AddFrame(fF6, fL2);
// Window name and final mapping
//
   SetWindowName("AliRoot Geometry Browser");
   MapSubwindows();
   // We need to use GetDefault...() to initialize the layout algorithm...
   Resize(GetDefaultSize());
   MapWindow();
}

AliGuiGeomMain::~AliGuiGeomMain()
{
   // Delete all created widgets.

   delete fCanvasWindow;

   delete fMenuBarLayout;
   delete fMenuBarItemLayout;
   delete fMenuBarHelpLayout;

   delete fMenuFile;
   delete fMenuTest;
   delete fMenuHelp;
}

void AliGuiGeomMain::Streamer(TBuffer &)
{
// Dummy streamer
;
}

void AliGuiGeomMain::Plot()
{
// plot de/dx or cross-sections
    const Float_t kAvo=0.60221367;
    Float_t *tkin  = new Float_t[fNbins];
    Float_t *value = new Float_t[fNbins];
    Float_t *pcut  = new Float_t[fNbins];
    Int_t ixst;
    Int_t imate=gCurrentMaterial->Id();
    Float_t z=gCurrentMaterial->Z();
    Float_t a=gCurrentMaterial->A();
    Float_t density=gCurrentMaterial->Density();
    
    Int_t ipart=gCurrentParticle;
    const char *kChMeca= kLabelTextMechanism[gCurrentProcess-1];
    char* tmp;
    tmp = new char[5];
    strncpy(tmp, kChMeca, 4);
    tmp[4]='\0';
    Int_t kdim=fNbins;
    Float_t de=(fEmax-fEmin)/fNbins;
    { //Begin local scope for i
      for (Int_t i=0; i<kdim; i++) {
	tkin[i]=fEmin+Float_t(i)*de;
	value[i]=0.;
      }
    } //End local scope for i
    if (kChMeca!="MUON") {
	((TGeant3*) gMC)->Gftmat(imate, ipart, tmp, kdim, tkin, value, pcut, ixst);
    } else {
	for (Int_t i=0; i<kdim; i++) {
	    Float_t ekin=tkin[i];
	    value[i]+=((TGeant3*) gMC)->Gbrelm(z,ekin,1.e10);
	    value[i]+=((TGeant3*) gMC)->Gprelm(z,ekin,1.e10);
	    value[i]*=1000.*kAvo*density/a;
	}
    }
    
    
    printf("\n %d %d %s %d \n",imate, ipart, kChMeca, kdim);
    if (ixst) {
	TGraph *graph = new TGraph(kdim,tkin,value);
	TCanvas *c1 = new TCanvas("c1"," ",400,10,600,700);
	c1->Divide(1,1);
	c1->cd(1);
	
	graph->SetFillColor(42);
	graph->SetMarkerColor(4);
	graph->SetMarkerStyle(21);
	graph->Draw("AC");
	graph->GetHistogram()->SetXTitle("Energy (GeV)");
	if (kChMeca == "RANG" || kChMeca == "STEP") {
	    graph->GetHistogram()->SetYTitle
		("Distance (cm)");   
	} else if (kChMeca == "LOSS" || kChMeca == "MUON") {
	    graph->GetHistogram()->SetYTitle("dE/dx (MeV/cm)");   
	} else {
	    graph->GetHistogram()->SetYTitle
		("Macroscopic X-Section (1/cm)"); 
	}
    }
    
    delete tkin;
    delete value;
    delete pcut;
    
    
}

void AliGuiGeomMain::Update()
{
// Update widgets
    if (fDialog) {
	fDialog->Update();
    }

    Int_t imat=gCurrentVolume->Material();
    Int_t nmat=fComboEntries->GetEntriesFast();
    Int_t i=0;
    for (i=0; i<nmat; i++) {
	gCurrentMaterial = (AliGUIMaterial*) 
	    (fComboEntries->UncheckedAt(i));
	if (gCurrentMaterial->Id()==imat) break;
    }
    Int_t imed = gCurrentVolume->Medium();
    Int_t nmed=fComboMediaEntries->GetEntriesFast();
    for (i=0; i<nmed; i++) {
	gCurrentMedium = (AliGUIMedium*) 
	    (fComboMediaEntries->UncheckedAt(i));
	if (gCurrentMedium->Id()==imed) break;
    }

    UpdateCombo();
    UpdateListBox();
    
}

void AliGuiGeomMain::UpdateCombo()
{
// Update combos

    Int_t   imat = gCurrentMaterial->Id();
    Float_t    a = gCurrentMaterial->A();
    Float_t    z = gCurrentMaterial->Z();    
    Float_t dens = gCurrentMaterial->Density();
    Float_t radl = gCurrentMaterial->RadiationLength();
    Float_t absl = gCurrentMaterial->AbsorptionLength();
    Int_t entry=gCurrentMaterial->ItemId();
//    
//
    fMaterialCombo->Select(entry);
    fMediaCombo->Select(entry);    
    char buf[10];
     
    sprintf(buf, "%10i", imat);
    fTbh[0]->Clear();
    fTbh[0]->AddText(0, buf);
    gClient->NeedRedraw(fTeh[0]);
    sprintf(buf, "%10.2e", a);
    fTbh[1]->Clear();
    fTbh[1]->AddText(0, buf);
    gClient->NeedRedraw(fTeh[1]);

    sprintf(buf, "%10.2e", z);
    fTbh[2]->Clear();
    fTbh[2]->AddText(0, buf);
    gClient->NeedRedraw(fTeh[2]);

    sprintf(buf, "%10.2e", dens);
    fTbh[3]->Clear();
    fTbh[3]->AddText(0, buf);
    gClient->NeedRedraw(fTeh[3]);

    sprintf(buf, "%10.2e", radl);
    fTbh[4]->Clear();
    fTbh[4]->AddText(0, buf);
    gClient->NeedRedraw(fTeh[4]);

    sprintf(buf, "%10.2e", absl);
    fTbh[5]->Clear();
    fTbh[5]->AddText(0, buf);
    gClient->NeedRedraw(fTeh[5]);

//  Media Combo
    sprintf(buf, "%10i", gCurrentMedium->Isvol());
    fTbhM[0]->Clear();
    fTbhM[0]->AddText(0, buf);
    gClient->NeedRedraw(fTehM[0]);


    sprintf(buf, "%10i", gCurrentMedium->Ifield());
    fTbhM[1]->Clear();
    fTbhM[1]->AddText(0, buf);
    gClient->NeedRedraw(fTehM[1]);

    sprintf(buf, "%10.2e", gCurrentMedium->Fieldm());
    fTbhM[2]->Clear();
    fTbhM[2]->AddText(0, buf);
    gClient->NeedRedraw(fTehM[2]);

    sprintf(buf, "%10.2e", gCurrentMedium->Tmaxfd());
    fTbhM[3]->Clear();
    fTbhM[3]->AddText(0, buf);
    gClient->NeedRedraw(fTehM[3]);

    sprintf(buf, "%10.2e", gCurrentMedium->Stemax());
    fTbhM[4]->Clear();
    fTbhM[4]->AddText(0, buf);
    gClient->NeedRedraw(fTehM[4]);

    sprintf(buf, "%10.2e", gCurrentMedium->Deemax());
    fTbhM[5]->Clear();
    fTbhM[5]->AddText(0, buf);
    gClient->NeedRedraw(fTehM[5]);

    sprintf(buf, "%10.2e", gCurrentMedium->Epsil());
    fTbhM[6]->Clear();
    fTbhM[6]->AddText(0, buf);
    gClient->NeedRedraw(fTehM[6]);

    sprintf(buf, "%10.2e", gCurrentMedium->Stmin());
    fTbhM[7]->Clear();
    fTbhM[7]->AddText(0, buf);
    gClient->NeedRedraw(fTehM[7]);
}

void AliGuiGeomMain::UpdateListBox()
{
// Update the list box
    Int_t i;
    fProcessLB->RemoveEntries(1,19);
    for (i=11; i < 30; i++) {
	Float_t p=gCurrentMedium->GetPar(i);
	char tmp[20];
	sprintf(tmp, "%6s%5d", kLabelTextP[i-11], Int_t(p));
	fProcessLB->AddEntry(tmp, i-10);
    }
    fProcessLB->MapSubwindows();
    fProcessLB->Layout();

    fCutsLB->RemoveEntries(1,10);
    for (i=1; i < 11; i++) {
	Float_t p=gCurrentMedium->GetPar(i);
	char tmp[20];
	sprintf(tmp, "%6s%10.3e", kLabelTextC[i-1], p);
	fCutsLB->AddEntry(tmp,i);
    }
    fCutsLB->MapSubwindows();
    fCutsLB->Layout();
}


void AliGuiGeomMain::CloseWindow()
{
   // Got close message for this MainFrame. Calls parent CloseWindow()
   // (which destroys the window) and terminate the application.
   // The close message is generated by the window manager when its close
   // window menu item is selected.

   TGMainFrame::CloseWindow();
   gApplication->Terminate(0);
}

Bool_t AliGuiGeomMain::ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2)
{
// Process messages to widgets
    switch (GET_MSG(msg)) {
//
//  Text entries for binning of cross-section plots
    case kC_TEXTENTRY:
	switch (GET_SUBMSG(msg)) {
	case kTE_TEXTCHANGED:
	    switch (parm1) {
	    case 10:
		fNbins=(Int_t) atof(fTbh61->GetString());
		break;
	    case 11:
		fEmin= atof(fTbh62->GetString());
		break;
	    case 12:
		fEmax= atof(fTbh63->GetString());
		break;
	    }
	    printf("\n %d %f %f binning", fNbins, fEmin, fEmax);
	}
	break;
//
// ListTree for volumes
    case kC_LISTTREE:
	switch (GET_SUBMSG(msg)) {
//
// Handle mouse click 
	case kCT_ITEMCLICK:
//
// Button 1: Select volume
	    if (parm1 == kButton1) {
	      TGListTreeItem *item;
	      if ((item = fLt->GetSelected())) 
		{
		  gCurrentVolume=((AliDrawVolume *) item->GetUserData());
		  Update();
		}
	    }
//
// Button 2: Draw volume specifications

	    if (parm1 == kButton2) {
		TGListTreeItem *item;
		if ((item = fLt->GetSelected())) 
		{

		    ((AliDrawVolume *) item->GetUserData())->DrawSpec();

		    gCurrentVolume=((AliDrawVolume *) item->GetUserData());
		    Update();
		}
	    }
//
// Button 3: Draw Volume
	    if (parm1 == kButton3) {
		TGListTreeItem *item;
		if ((item = fLt->GetSelected())) 
		{
		    ((AliDrawVolume *) item->GetUserData())->Draw();
		    gCurrentVolume=((AliDrawVolume *) item->GetUserData());
		    Update();
		}
	    }
	    
	    
	    break;
	case kCT_ITEMDBLCLICK:
	    if (parm1 == kButton1) {
		if (fLt->GetSelected() != 0) {
		    gClient->NeedRedraw(fLt);
		}
	    }
	    break;
	default:
	    break;
	}
	break;
    case kC_COMMAND:
	switch (GET_SUBMSG(msg)) {
	case kCM_BUTTON:
	    switch(parm1) {
	    case 1:
		Plot();
		break;
	    }
	    break;
	case kCM_COMBOBOX:
//
// Combo box bindings
	    switch(parm1) {
//
// Material Combo
	    case 1:
		gCurrentMaterial=(AliGUIMaterial*) 
		    (fComboEntries->UncheckedAt(Int_t(parm2-1)));
		gCurrentMedium=(AliGUIMedium*) 
		    (fComboMediaEntries->UncheckedAt(Int_t(parm2-1)));
		UpdateCombo();
		UpdateListBox();
		break;
//
// Media Combo
	    case 2:
		gCurrentMedium=(AliGUIMedium*) 
		    (fComboMediaEntries->UncheckedAt(Int_t(parm2-1)));
		gCurrentMaterial=(AliGUIMaterial*) 
		    (fComboEntries->UncheckedAt(Int_t(parm2-1)));
		UpdateCombo();
		UpdateListBox();
		break;
//
// Particle Combo
	    case 3:
		gCurrentParticle=Int_t(parm2);
		break;
//
// Mechanism Combo
	    case 4:
		gCurrentProcess=Int_t(parm2);
		printf("\n Process %d", gCurrentProcess);
		break;
	    }
	    break;
	case kCM_MENUSELECT:
	    break;
	    
	case kCM_MENU:
	    switch (parm1) {
		
	    case kFileOpen:
	    {
		TGFileInfo fi;
		fi.fFileTypes = (char **)kFileTypes;
		new TGFileDialog(gClient->GetRoot(), this, kFDOpen,&fi);
	    }
	    break;
	    
	    case kTestDlg:
		fDialog = new AliGuiGeomDialog
		    (gClient->GetRoot(), this, 400, 200);
		break;
		
	    case kFileSave:
		printf("kFileSave\n");
		break;
		
	    case kFileExit:
		CloseWindow();   // this also terminates theApp
		break;
	    default:
		break;
	    }
	default:
	    break;
	}
    default:
	break;
    }
    return kTRUE;
}

void AliGuiGeomMain::AddMaterial(AliGUIMaterial *Material, Int_t i)
{
// Add material to material combo
    char* tmp;
    tmp=Material->Name();
    Material->SetItemId(i);
    fMaterialCombo->AddEntry(tmp, i);
    fMaterialCombo->Select(i);
    fMaterialCombo->Resize(200, 20);
}

void AliGuiGeomMain::AddMedium(AliGUIMedium *Medium, Int_t i)
{
// Add medium to medium combo
    char* tmp;
    tmp=Medium->Name();
    Medium->SetItemId(i);
    
    fMediaCombo->AddEntry(tmp, i);
    fMediaCombo->Select(i);
    fMediaCombo->Resize(200, 20);
}


