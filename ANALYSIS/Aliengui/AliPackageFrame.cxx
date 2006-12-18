/**************************************************************************
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

/* $Id$ */

//-----------------------------------------------------------------
//           AliPackageFrame class
//   The class that deals with the package tab of the GUI
//   Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-----------------------------------------------------------------

#include "TGLabel.h"
#include "TGTextEntry.h"
#include "TGButton.h"
#include "TGFileDialog.h"

#include "TSystem.h"
#include "TObjString.h"

#include "AliAnalysisGUI.h"
#include "AliPackageFrame.h"

ClassImp(AliPackageFrame)

//___________________________________________________________________________
  AliPackageFrame::AliPackageFrame (const TGWindow *main, UInt_t w, UInt_t h, AliAnalysisGUI * v): 
    TGHorizontalFrame(main, w, h), 
    fVFrame1(0), fVFrame2(0),
    fLabel1(0), fTextPackage(0), fButtonSelect(0), fButtonBuild(0),
    fAliAnalysisGUI(v) {
   // ctor.
  
   fVFrame1 = new TGVerticalFrame(this, 100, 100);
   AddFrame(fVFrame1, new TGLayoutHints(kLHintsLeft, 5,5,5,5));

   fLabel1 = new TGLabel(fVFrame1, new TGString("PAR file"));
   fVFrame1->AddFrame(fLabel1, new TGLayoutHints(kLHintsTop, 5,5,10,5));
   
   fTextPackage = new TGTextEntry(fVFrame1, new TGTextBuffer(50));
   fVFrame1->AddFrame(fTextPackage, new TGLayoutHints(kLHintsBottom, 5,5,5,5));
   fTextPackage->SetEnabled(false);

   fVFrame2 = new TGVerticalFrame(this, 100, 200);
   AddFrame(fVFrame2, new TGLayoutHints(kLHintsLeft, 5,5,5,5));

   fButtonSelect = new TGTextButton(fVFrame2, "Select...", 1);
   fVFrame2->AddFrame(fButtonSelect,  new TGLayoutHints(kLHintsExpandX | kLHintsTop, 5,5,5,5));
   fButtonSelect->Connect("Clicked()", "AliPackageFrame", this, "OnSelect()");

   fButtonBuild = new TGTextButton(fVFrame2, "Build...", 2);
   fVFrame2->AddFrame(fButtonBuild,  new TGLayoutHints(kLHintsExpandX | kLHintsBottom, 5,5,5,5));

   fButtonBuild->Connect("Clicked()", "AliPackageFrame", this, "OnBuild()");

   MapWindow();
   Resize();
   MapSubwindows();
}

//___________________________________________________________________________
void AliPackageFrame::OnSelect() {
   // When Select button is pressed.
   const char *filetypes[] = { 
     "PAR files",    "*.par",
     0,               0 
   };


   static TString dir(".");
   TGFileInfo fi;
   fi.fFileTypes = filetypes;
   fi.fIniDir    = StrDup(dir);

   new TGFileDialog(gClient->GetRoot(), fAliAnalysisGUI, kFDOpen, &fi);

   fTextPackage->SetText(fi.fFilename);

}

//___________________________________________________________________________
void AliPackageFrame::OnBuild() {

   // When Build button is pressed.  

   TString fname (fTextPackage->GetText());

   TObjArray *a = fname.Tokenize("/");
   TObjString * ostr = (TObjString*)a->At(a->GetEntries()-1);

   CreatePARFile(ostr->GetString().Data());

}

//___________________________________________________________________________
void AliPackageFrame::CreatePARFile(const char* pararchivename) {
   // Create PAR file

   if (pararchivename) {
     char processline[1024];
     //     sprintf(processline,".! tar xvzf ../tags/%s.par",pararchivename);
     sprintf(processline,".! tar xvzf %s",pararchivename);
     gROOT->ProcessLine(processline);

     //    const char* ocwd = gSystem->WorkingDirectory();
     gSystem->ChangeDirectory("ESD");

     // check for BUILD.sh and execute
     if (!gSystem->AccessPathName("PROOF-INF/BUILD.sh")) {
       printf("*******************************\n");
       printf("*** Building PAR archive    ***\n");
       printf("*******************************\n");
       if (gSystem->Exec("PROOF-INF/BUILD.sh")) {
 	 Error("batchSelector","Cannot Build the PAR Archive! - Abort!");
	 return;
       }           

       if (!gSystem->AccessPathName("PROOF-INF/SETUP.C")) {

 	 printf("*******************************\n");
	 printf("*** Setup PAR archive       ***\n");
	 printf("*******************************\n");

	 gROOT->Macro("PROOF-INF/SETUP.C");
       }

       gSystem->ChangeDirectory("../");

     }

   }
}
