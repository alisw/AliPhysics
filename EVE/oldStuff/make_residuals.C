/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
 
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TCanvas.h>
#include <TGButton.h>
#include <TGButtonGroup.h>
#include <TFile.h>
#include <TGLabel.h>
#include <TGNumberEntry.h>
#include <TG3DLine.h>
#include <TApplication.h>
#include <TGComboBox.h>
#include <TLatex.h>
#include <TTree.h>
#include <TEveUtil.h>

#include <AliCDBManager.h>
#include <AliESDEvent.h>
#include <AliESDfriendTrack.h>
#include <AliGeomManager.h>
#include <AliEveEventManager.h>

/* Not sure which ConfigCalibTrain.C macro ? 
 * From ANALYSIS or from PWGPP?
 */
#include <ANALYSIS/macros/ConfigCalibTrain.C>
#endif

class ButtonWindow : public TGMainFrame {

protected:
   TGComboBox *option1;
   TGComboBox *option2;
   TGComboBox *option3;
   TGComboBox *option4;
   TGComboBox *option5;
   TGComboBox *option6;
   TGComboBox *option7;
   TGComboBox *option8;
   TGComboBox *option9;
   TGComboBox *option10;
   TGTextEntry *cut1;
   TGTextEntry *cut2;
   TGTextEntry *cut3;
   TGTextEntry *cut4;
   TGTextEntry *customCutSelection;
   TGTextEntry *customDrawSelection;
   TGNumberEntry *nEntries;
   TGNumberEntry *firstEntry;

public:
   ButtonWindow();
   void DrawResiduals();
   
   ClassDef(ButtonWindow, 0)
};

//________________________________________________

ButtonWindow::ButtonWindow() : TGMainFrame(gClient->GetRoot(), 10, 10, kHorizontalFrame)
{
   // Main test window.

   SetCleanup(kDeepCleanup);

   // Controls on right
   TGVerticalFrame *controls = new TGVerticalFrame(this);
   AddFrame(controls, new TGLayoutHints(kLHintsRight | kLHintsExpandY, 5, 5, 5, 5));

   // control margins of the text

   TGLabel *label1 = 0;
   TGLabel *label2 = 0;
   TGLabel *label3 = 0;
  
   TGHorizontal3DLine *separator = 0;

   TGGroupFrame *margins = new TGGroupFrame(controls, "Residuals Drawing Options");
   margins->SetTitlePos(TGGroupFrame::kCenter);

//==========================

   separator = new TGHorizontal3DLine(margins);

   margins->AddFrame(separator, new TGLayoutHints(kLHintsExpandX));

   label1 = new TGLabel(margins, "Axes");

   margins->AddFrame(label1, new TGLayoutHints(kLHintsExpandX, 5, 5, 3, 4));

   TGHorizontalFrame *hframe1 = new TGHorizontalFrame(margins, 400, 20, kFixedWidth);
   
   label1 = new TGLabel(hframe1, "X axis: ");

   option1 = new TGComboBox(hframe1,"Draw option");
   option1->AddEntry("fX",1);
   option1->AddEntry("fY",2);
   option1->AddEntry("fZ",3);
   option1->AddEntry("fR",4);
   option1->Resize(120,20);

   label3 = new TGLabel(hframe1, "Y axis: ");

   option2 = new TGComboBox(hframe1,"Draw option");
   option2->AddEntry("fPx",1);
   option2->AddEntry("fPy",2);
   option2->AddEntry("fPz",3);
   option2->AddEntry("fPt",4);
   option2->AddEntry("fP",5);
   option2->Resize(120,20);

   hframe1->AddFrame(label1, new TGLayoutHints(kLHintsExpandX, 5, 5, 3, 4));
   hframe1->AddFrame(option1, new TGLayoutHints(kLHintsExpandX, 5, 5, 3, 4));
   hframe1->AddFrame(label3, new TGLayoutHints(kLHintsExpandX, 5, 5, 3, 4));
   hframe1->AddFrame(option2, new TGLayoutHints(kLHintsExpandX, 5, 5, 3, 4));

   margins->AddFrame(hframe1, new TGLayoutHints(kLHintsExpandX, 5, 5, 3, 4));

//==========================

   separator = new TGHorizontal3DLine(margins);

   margins->AddFrame(separator, new TGLayoutHints(kLHintsExpandX));

   label1 = new TGLabel(margins, "Cut Selection");

   margins->AddFrame(label1, new TGLayoutHints(kLHintsExpandX, 5, 5, 3, 4));

//===============================

   TGHorizontalFrame *hframe3 = new TGHorizontalFrame(margins, 400, 20, kFixedWidth);
   
   label1 = new TGLabel(hframe3, " 1. ");

   option3 = new TGComboBox(hframe3,"Cut");
   option3->AddEntry("fPx",1);
   option3->AddEntry("fPy",2);
   option3->AddEntry("fPz",3);
   option3->AddEntry("fPt",4);
   option3->AddEntry("fP",5);
   option3->Resize(100,20);

   option4 = new TGComboBox(hframe3,"-");
   option4->AddEntry("(no cut)",0);
   option4->AddEntry("<",1);
   option4->AddEntry(">",2);
   option4->Resize(100,20);

   cut1 = new TGTextEntry(hframe3);

   hframe3->AddFrame(label1, new TGLayoutHints(kLHintsExpandX, 5, 5, 1, 1));
   hframe3->AddFrame(option3, new TGLayoutHints(kLHintsExpandX, 5, 5, 1, 1));
   hframe3->AddFrame(option4, new TGLayoutHints(kLHintsExpandX, 5, 5, 1, 1));
   hframe3->AddFrame(cut1, new TGLayoutHints(kLHintsExpandX, 5, 5, 1, 1));

   margins->AddFrame(hframe3, new TGLayoutHints(kLHintsExpandX, 5, 5, 1, 1));

//=====================================

   TGHorizontalFrame *hframe4 = new TGHorizontalFrame(margins, 400, 20, kFixedWidth);
   
   label1 = new TGLabel(hframe4, "2. ");

   option5 = new TGComboBox(hframe4,"Cut");
   option5->AddEntry("fPx",1);
   option5->AddEntry("fPy",2);
   option5->AddEntry("fPz",3);
   option5->AddEntry("fPt",4);
   option5->AddEntry("fP",5);
   option5->Resize(100,20);

   option6 = new TGComboBox(hframe4,"-");
   option6->AddEntry("(no cut)",0);
   option6->AddEntry("<",1);
   option6->AddEntry(">",2);
   option6->Resize(100,20);

   cut2 = new TGTextEntry(hframe4);
//   cut2 = new TGNumberEntryField(hframe4, 40, 20, kFixedWidth);

   hframe4->AddFrame(label1, new TGLayoutHints(kLHintsExpandX, 5, 5, 1, 1));
   hframe4->AddFrame(option5, new TGLayoutHints(kLHintsExpandX, 5, 5, 1, 1));
   hframe4->AddFrame(option6, new TGLayoutHints(kLHintsExpandX, 5, 5, 1, 1));
   hframe4->AddFrame(cut2, new TGLayoutHints(kLHintsExpandX, 5, 5, 1, 1));

   margins->AddFrame(hframe4, new TGLayoutHints(kLHintsExpandX, 5, 5, 1, 1));

   separator = new TGHorizontal3DLine(margins);

   margins->AddFrame(separator, new TGLayoutHints(kLHintsExpandX));

   label1 = new TGLabel(margins, "Custom Cut Selection");

   margins->AddFrame(label1, new TGLayoutHints(kLHintsExpandX, 5, 5, 3, 4));

   customCutSelection = new TGTextEntry(margins);

   margins->AddFrame(customCutSelection, new TGLayoutHints(kLHintsExpandX, 5, 5, 3, 4));

//========================================

   separator = new TGHorizontal3DLine(margins);

   margins->AddFrame(separator, new TGLayoutHints(kLHintsExpandX));

   label1 = new TGLabel(margins, "Draw Selection");

   margins->AddFrame(label1, new TGLayoutHints(kLHintsExpandX, 5, 5, 3, 4));

//========================================

   TGHorizontalFrame *hframe5 = new TGHorizontalFrame(margins, 400, 20, kFixedWidth);
   
   label1 = new TGLabel(hframe5, "1. ");

   option7 = new TGComboBox(hframe5,"Cut");
   option7->AddEntry("fPx",1);
   option7->AddEntry("fPy",2);
   option7->AddEntry("fPz",3);
   option7->AddEntry("fPt",4);
   option7->AddEntry("fP",5);
   option7->Resize(100,20);

   option8 = new TGComboBox(hframe5,"-");
   option8->AddEntry("(no cut)",0);
   option8->AddEntry("<",1);
   option8->AddEntry(">",2);
   option8->Resize(100,20);

   cut3 = new TGTextEntry(hframe5);
//   cut3 = new TGNumberEntryField(hframe5, 40, 20, kFixedWidth);

   hframe5->AddFrame(label1, new TGLayoutHints(kLHintsExpandX, 5, 5, 1, 1));
   hframe5->AddFrame(option7, new TGLayoutHints(kLHintsExpandX, 5, 5, 1, 1));
   hframe5->AddFrame(option8, new TGLayoutHints(kLHintsExpandX, 5, 5, 1, 1));
   hframe5->AddFrame(cut3, new TGLayoutHints(kLHintsExpandX, 5, 5, 1, 1));

   margins->AddFrame(hframe5, new TGLayoutHints(kLHintsExpandX, 5, 5, 1, 1));

//=======================================

   TGHorizontalFrame *hframe6 = new TGHorizontalFrame(margins, 400, 20, kFixedWidth);
   
   label1 = new TGLabel(hframe6, "2. ");

   option9 = new TGComboBox(hframe6,"Cut");
   option9->AddEntry("fPx",1);
   option9->AddEntry("fPy",2);
   option9->AddEntry("fPz",3);
   option9->AddEntry("fPt",4);
   option9->AddEntry("fP",5);
   option9->Resize(100,20);

   option10 = new TGComboBox(hframe6,"-");
   option10->AddEntry("(no cut)",0);
   option10->AddEntry("<",1);
   option10->AddEntry(">",2);
   option10->Resize(100,20);

   cut4 = new TGTextEntry(hframe6);
//   cut4 = new TGNumberEntryField(hframe6, 40, 20, kFixedWidth);

   hframe6->AddFrame(label1, new TGLayoutHints(kLHintsExpandX, 5, 5, 1, 1));
   hframe6->AddFrame(option9, new TGLayoutHints(kLHintsExpandX, 5, 5, 1, 1));
   hframe6->AddFrame(option10, new TGLayoutHints(kLHintsExpandX, 5, 5, 1, 1));
   hframe6->AddFrame(cut4, new TGLayoutHints(kLHintsExpandX, 5, 5, 1, 1));

   margins->AddFrame(hframe6, new TGLayoutHints(kLHintsExpandX, 5, 5, 1, 1));

   separator = new TGHorizontal3DLine(margins);

   margins->AddFrame(separator, new TGLayoutHints(kLHintsExpandX));

   label1 = new TGLabel(margins, "Custom Draw Selection");

   margins->AddFrame(label1, new TGLayoutHints(kLHintsExpandX, 5, 5, 3, 4));

   customDrawSelection = new TGTextEntry(margins);

   margins->AddFrame(customDrawSelection, new TGLayoutHints(kLHintsExpandX, 5, 5, 3, 4));

//==========================

   separator = new TGHorizontal3DLine(margins);

   margins->AddFrame(separator, new TGLayoutHints(kLHintsExpandX));

   label1 = new TGLabel(margins, "Entries Selection");

   margins->AddFrame(label1, new TGLayoutHints(kLHintsExpandX, 5, 5, 3, 4));

//========================================

   TGHorizontalFrame *hframe7 = new TGHorizontalFrame(margins, 400, 20, kFixedWidth);
   
   label1 = new TGLabel(hframe7, "nEntries");

   nEntries = new TGNumberEntry(hframe7, 40, 20, kFixedWidth);

   label2 = new TGLabel(hframe7, "firstEntry");

   firstEntry = new TGNumberEntry(hframe7, 40, 20, kFixedWidth);

   hframe7->AddFrame(label1, new TGLayoutHints(kLHintsExpandX, 5, 5, 1, 1));
   hframe7->AddFrame(nEntries, new TGLayoutHints(kLHintsExpandX, 5, 5, 1, 1));
   hframe7->AddFrame(label2, new TGLayoutHints(kLHintsExpandX, 5, 5, 1, 1));
   hframe7->AddFrame(firstEntry, new TGLayoutHints(kLHintsExpandX, 5, 5, 1, 1));

   margins->AddFrame(hframe7, new TGLayoutHints(kLHintsExpandX, 5, 5, 1, 1));

//==========================

   separator = new TGHorizontal3DLine(margins);

   margins->AddFrame(separator, new TGLayoutHints(kLHintsExpandX));

   const TGFont *font = gClient->GetFont("-*-times-bold-r-*-*-16-*-*-*-*-*-*-*");
//   const TGFont *font = gClient->GetFont("-*-symbol-medium-r-normal-*-16-*-*-*-*-*-*-*");
   FontStruct_t buttonFont = font->GetFontStruct();
   ULong_t buttonRedColor;
   gClient->GetColorByName("red", buttonRedColor);
   TGTextButton *draw = new TGTextButton(margins,"Draw Residuals");
   draw->SetTextColor(buttonRedColor);
   draw->SetFont(buttonFont);
   draw->Connect("Clicked()", "ButtonWindow", this, "DrawResiduals()");
   margins->AddFrame(draw, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));

   controls->AddFrame(margins, new TGLayoutHints(kLHintsExpandX));

   MapSubwindows();
   Resize();

   SetWMSizeHints(GetDefaultWidth(), GetDefaultHeight(), 1000, 1000, 0 ,0);
   SetWindowName("Residuals");
   MapRaised();
}

//______________________________________________________________________________
void ButtonWindow::DrawResiduals()
{

      TString selection1, selection2, selection3, selection4, selection5, selection6, selection7, selection8, selection9, selection10;

      switch(option1->GetSelected())
      {
         case 1:
            selection1 = "fX";
            break;
         case 2:
            selection1 = "fY";
            break;
         case 3:
            selection1 = "fZ";
            break;
         default:
            selection1 = "fX";
            break;
      }
    
      switch(option2->GetSelected())
      {
         case 1:
            selection2 = "fPx";
            break;
         case 2:
            selection2 = "fPy";
            break;
         case 3:
            selection2 = "fPz";
            break;
         case 4:
            selection2 = "fPt";
            break;
         case 5:
            selection2 = "fP";
            break;
         default:
            selection2 = "fP";
            break;
      }

      switch(option3->GetSelected())
      {
         case 1:
            selection3 = "fPx";
            break;
         case 2:
            selection3 = "fPy";
            break;
         case 3:
            selection3 = "fPz";
            break;
         case 4:
            selection3 = "fPt";
            break;
         case 5:
            selection3 = "fP";
            break;
         default:
            selection3 = "fP";
            break;
      }

      switch(option4->GetSelected())
      {
         case 1:
            selection4 = "<";
            break;
         case 2:
            selection4 = ">";
            break;
         default:
            selection4 = "<";
            break;
      }

      switch(option5->GetSelected())
      {
         case 1:
            selection5 = "fPx";
            break;
         case 2:
            selection5 = "fPy";
            break;
         case 3:
            selection5 = "fPz";
            break;
         case 4:
            selection5 = "fPt";
            break;
         case 5:
            selection5 = "fP";
            break;
         default:
            selection5 = "fP";
            break;
      }

      switch(option6->GetSelected())
      {
         case 1:
            selection6 = "<";
            break;
         case 2:
            selection6 = ">";
            break;
         default:
            selection6 = "<";
            break;
      }

      switch(option7->GetSelected())
      {
         case 1:
            selection7 = "fPx";
            break;
         case 2:
            selection7 = "fPy";
            break;
         case 3:
            selection7 = "fPz";
            break;
         case 4:
            selection7 = "fPt";
            break;
         case 5:
            selection7 = "fP";
            break;
         default:
            selection7 = "fP";
            break;
      }

      switch(option8->GetSelected())
      {
         case 1:
            selection8 = "<";
            break;
         case 2:
            selection8 = ">";
            break;
         default:
            selection8 = "<";
            break;
      }

      switch(option9->GetSelected())
      {
         case 1:
            selection9 = "fPx";
            break;
         case 2:
            selection9 = "fPy";
            break;
         case 3:
            selection9 = "fPz";
            break;
         case 4:
            selection9 = "fPt";
            break;
         case 5:
            selection9 = "fP";
            break;
         default:
            selection3 = "fP";
            break;
      }

      switch(option10->GetSelected())
      {
         case 1:
            selection10 = "<";
            break;
         case 2:
            selection10 = ">";
            break;
         default:
            selection10 = "<";
            break;
      }

      TString cutSelectionMerged;// = ("abs(ESDfriend.fTracks[].fTPCOut."+selection3+"[4])"+selection4+cut1->GetText());
      TString drawSelectionMerged1 = ("track[]."+selection2+"[0]:track[]."+selection1);
      TString drawSelectionMerged2;// = ("abs(track[]."+selection7+"[0])"+selection8+cut3->GetText());

      if(customCutSelection->GetText())
        cutSelectionMerged = customCutSelection->GetText();

      if(customDrawSelection->GetText())
        drawSelectionMerged2 = customDrawSelection->GetText();


      Info("make_residuals::DrawResiduals", "%s abs(ESDfriend.fTracks[].fTPCOut.fP[4])<0.5", cutSelectionMerged.Data());
      Info("make_residuals::DrawResiduals", "%s track[].fP[0]:track[].fX", drawSelectionMerged1.Data());
      Info("make_residuals::DrawResiduals", "%s abs(track[].fP[0])<20", drawSelectionMerged2.Data());
      Info("make_residuals::DrawResiduals", "nEntries: %f", nEntries->GetNumber());
      Info("make_residuals::DrawResiduals", "firstEntry: %f", firstEntry->GetNumber());
    
      AliESDEvent *esd = AliEveEventManager::Instance()->AssertESD();
    
        // OCDB
        
//        printf("setting run to %d\n",esd->GetRunNumber());
//        AliCDBManager::Instance()->SetDefaultStorage("raw://");
//        AliCDBManager::Instance()->SetRun(esd->GetRunNumber());
    
        // geometry
//        AliGeomManager::LoadGeometry();
//        AliGeomManager::ApplyAlignObjsFromCDB("GRP ITS TPC");
    
    

      /* OBSOLETE CODE - No function members defined
       * 
      AliESDfriendTrack::MakeDebugStreamer(kTRUE);

      TTreeSRedirector *pcstreamRes = AliESDfriendTrack::GetDebugStreamer();
      */
      
      TFile fRes("AliESDfriends.root");

      TTree *treeRes = (TTree*)fRes.Get("esdFriendTree");     

      TCanvas* cnv = new TCanvas();
      
      
      treeRes->Draw("ESDfriend.fTracks[].MakeResidualGraph(1)",cutSelectionMerged," ",nEntries->GetNumber(),firstEntry->GetNumber());
//      treeRes->Draw("ESDfriend.fTracks[].MakeResidualGraph(1)","abs(ESDfriend.fTracks[].fTPCOut.fP[4])<0.5"," ",nEntries,firstEntry);

      //delete pcstreamRes;

      /*  OBSOLETE CODE -dump not defined

      TFile fRes("residual.root");

      dump.Draw(drawSelectionMerged1,drawSelectionMerged2);
      */
//      dump.Draw("track[].fP[0]:track[].fX","abs(track[].fP[0])<20");

}

//_____________________________________________________________________________
void make_residuals()
{

      new ButtonWindow();

}
