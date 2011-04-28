//
// Author: Michele Floris APR2011
// Derived by root's the statusBar.C tutorial
// This macro allows to browse the OADB of the spectra


#include <TApplication.h>
#include <TGClient.h>
#include <TGButton.h>
#include <TGFrame.h>
#include <TFrame.h>
#include <TRootEmbeddedCanvas.h>
#include <TGStatusBar.h>
#include <TGComboBox.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TRandom.h>
#include <TGraph.h>
#include <TAxis.h>
#include <AliOADBPWG2Spectra.h>
#include "TSystem.h"
#include "AliOADBContainer.h"
#include <TGTextEntry.h>
#include "TFile.h"
#include <iostream>
#include "TH1D.h"
#include "TInterpreter.h"

using namespace std;


class MyMainFrame : public TGMainFrame {

private:
  TRootEmbeddedCanvas  *fEcan; // embedded canvas
  TGStatusBar          *fStatusBar;// status bar
  TGComboBox           *fComboDetector;//detectors
  TGComboBox           *fComboParticle;//particle
  TGComboBox           *fComboCharge;//charge
  TGComboBox           *fComboPID;//PID
  TGComboBox           *fComboContainer;//container
  TGTextEntry          *fRunNumber; // run number 
  TGComboBox           *fComboCentr; // centr tag
  TGTextEntry          *fTextCentrBin; // centrality bin
  TGTextEntry          *fTextDrawOpt; // draw opt
  TGTextEntry          *fTextSaveName; // save name
  AliOADBPWG2Spectra   *fOADBSpectra; // spectra OADB
  AliOADBContainer     *fOADBContainer; // OADB container
   
  static const char * fkContainers[];
  static const Int_t fkNContainers ;
  static const char * fkCentrTags[];
  static const Int_t fkNCentrTags ;

public:
  MyMainFrame(const TGWindow *p, UInt_t w, UInt_t h);
  virtual ~MyMainFrame();
  void DoExit();
  void DoDraw();
  void DoRatio();
  void DoLegend();
  void DoLoad();
  void DoClear();
  void DoPrint();
  void DoSelectedDetector(Int_t id);
  void SetStatusText(const char *txt, Int_t pi);
  void EventInfo(Int_t event, Int_t px, Int_t py, TObject *selected);
   
  ClassDef(MyMainFrame, 0)
};

const char * MyMainFrame::fkContainers[] = {"Corrected", "Raw"};
const Int_t MyMainFrame::fkNContainers = 2;
const char * MyMainFrame::fkCentrTags[] = {"MB", "Ntracks", "SPD2"};
const Int_t MyMainFrame::fkNCentrTags = 3;

void MyMainFrame::DoSelectedDetector(Int_t id) {
  // Change the analysis type based on the detector
  cout << "ID " << id << endl;
  if (id == AliOADBPWG2Spectra::kTOFTPC) {
    fComboPID->Select(AliOADBPWG2Spectra::kNSigma);
  }
  else {
    fComboPID->Select(AliOADBPWG2Spectra::kGaussFit);
  }

}

void MyMainFrame::DoPrint()
{
  // print to file
  TCanvas *c1 = fEcan->GetCanvas();
  c1->Print(fTextSaveName->GetText());
}
void MyMainFrame::DoClear()
{
  // clear canvas
  TCanvas *c1 = fEcan->GetCanvas();
  c1->Clear();
  c1->Update();
  c1->Modified();
  c1->Update(); 

}

void MyMainFrame::DoLoad()
{
  // Load file
  cout << "Getting " << fkContainers[fComboContainer->GetSelected()] << endl;
  cout << fRunNumber->GetText() << endl;
  
  TString fileName = AliOADBPWG2Spectra::GetOADBPWG2SpectraFileName();
  TFile * f = new TFile (fileName);
  fOADBContainer = (AliOADBContainer*) f->Get(fkContainers[fComboContainer->GetSelected()]);
  f->Close();
  fOADBSpectra = (AliOADBPWG2Spectra*) fOADBContainer->GetObject(atoi(fRunNumber->GetText() ));
  if(!fOADBSpectra) fOADBContainer->List();
  fOADBSpectra->Print();
}

void MyMainFrame::DoLegend()
{
  // legend
  gInterpreter->ProcessLine("NewLegend(\"\", \"lpf\",0,1,0);");
  TCanvas *c1 = fEcan->GetCanvas();
  c1->Update();
  c1->Modified();
  c1->Update(); 
}
void MyMainFrame::DoRatio() {
  // Divide 2 histos on canvas
  gInterpreter->ProcessLine("Divide2HistosOnCanvas();");
}


void MyMainFrame::DoDraw()
{
   // Draw something in the canvas

   Printf("Slot DoDraw()");
   cout << "DET " << fComboDetector->GetSelected() << endl;
   if(!fOADBSpectra) {
     cout << "spectra not loaded" << endl;
     DoLoad();
     if(!fOADBSpectra)      {
       cout << "ERROR: Cannot load spectra" << endl;       
       return;
     }
   }

   TCanvas *c1 = fEcan->GetCanvas();
   c1->cd();
   TH1D * h = 0;
   if(strcmp(fkCentrTags[fComboCentr->GetSelected()],"")){
       h = fOADBSpectra->GetHisto(fComboDetector->GetSelected(), 
				  fComboPID->GetSelected(), 
				  fComboParticle->GetSelected(), 
				  fComboCharge->GetSelected(), 
				  fkCentrTags[fComboCentr->GetSelected()],
				  atoi(fTextCentrBin->GetText()));

     }
     else {
       h = fOADBSpectra->GetHisto(fComboDetector->GetSelected(), 
				  fComboPID->GetSelected(), 
				  fComboParticle->GetSelected(), 
				  fComboCharge->GetSelected());

     }
   // Draw the selected histogram

   if(!h) {
     cout << "Cannot get pointer to histo" << endl;
     
     cout << fkCentrTags[fComboCentr->GetSelected()] << " " << 
   
     fOADBSpectra->GetHistoName(fComboDetector->GetSelected(), 
				fComboPID->GetSelected(), 
				fComboParticle->GetSelected(), 
				fComboCharge->GetSelected(), 
				fkCentrTags[fComboCentr->GetSelected()],
				atoi(fTextCentrBin->GetText()))
	<< endl;

     return;
   }
   TString opt = fTextDrawOpt->GetText();
   if(opt=="auto") {
     c1->GetListOfPrimitives()->Print();
     if (c1->GetListOfPrimitives()->GetEntries()>0) opt = "same";
     else opt = "";
   }
   
   h->SetXTitle("p_{T} (GeV/c)");
   //   h->SetXTitle("dN/dp_{T}");
   h->Draw(opt);
   
   // TCanvas::Update() draws the frame, after which it can be changed
   c1->Update();
   c1->Modified();
   c1->Update();
}

void MyMainFrame::DoExit()
{
   printf("Exit application...");
   gApplication->Terminate(0);
}

void MyMainFrame::SetStatusText(const char *txt, Int_t pi)
{
   // Set text in status bar.
   fStatusBar->SetText(txt,pi);
}

void MyMainFrame::EventInfo(Int_t event, Int_t px, Int_t py, TObject *selected)
{
//  Writes the event status in the status bar parts

   const char *text0, *text1, *text3;
   char text2[50];
   text0 = selected->GetTitle();
   SetStatusText(text0,0);
   text1 = selected->GetName();
   SetStatusText(text1,1);
   if (event == kKeyPress)
      sprintf(text2, "%c", (char) px);
   else
      sprintf(text2, "%d,%d", px, py);
   SetStatusText(text2,2);
   text3 = selected->GetObjectInfo(px,py);
   SetStatusText(text3,3);
}

MyMainFrame::MyMainFrame(const TGWindow *p, UInt_t w, UInt_t h) :
  TGMainFrame(p, w, h), fOADBContainer(0), fOADBSpectra(0)
{
  // Create a horizontal frame (canvas and status bar on the left, buttons on the right)
   TGHorizontalFrame *hframeMain = new TGHorizontalFrame(this, 200, 40);

   // Create a vertical frame for the canvas and for the status bar
   TGVerticalFrame *vframeCanvas = new TGVerticalFrame(this, 200, 40);

   // Create the embedded canvas
   fEcan = new TRootEmbeddedCanvas(0,this,500,400);
   Int_t wid = fEcan->GetCanvasWindowId();
   TCanvas *myc = new TCanvas("MyCanvas", 10,10,wid);
   cout << myc << endl;
   
   fEcan->AdoptCanvas(myc);
   myc->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)","MyMainFrame",this, 
               "EventInfo(Int_t,Int_t,Int_t,TObject*)");

   vframeCanvas->AddFrame(fEcan, new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX  | kLHintsExpandY,0,0,1,1));
   // status bar
   Int_t parts[] = {45, 15, 10, 30};
   fStatusBar = new TGStatusBar(this, 50, 10, kVerticalFrame);
   fStatusBar->SetParts(parts, 4);
   fStatusBar->Draw3DCorner(kFALSE);
   vframeCanvas->AddFrame(fStatusBar, new TGLayoutHints(kLHintsExpandX, 0, 0, 10, 0));
   
   hframeMain->AddFrame(vframeCanvas, new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX  | kLHintsExpandY, 2, 2, 2, 2));

   // Create a vertical frame containing  buttons and controls
   TGVerticalFrame *vframeButtons = new TGVerticalFrame(this, 200, 40);
  
   // Container
   fComboContainer = new TGComboBox(vframeButtons);
   vframeButtons->AddFrame(fComboContainer, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));
   for(Int_t icont = 0; icont < fkNContainers; icont++){
     fComboContainer->AddEntry(fkContainers[icont],icont);
   }
   fComboContainer->Select(0);
   fComboContainer->Resize(120,20);
   fRunNumber = new TGTextEntry (vframeButtons, "116562");
   vframeButtons->AddFrame(fRunNumber, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));
   

   TGTextButton *load = new TGTextButton(vframeButtons, "&Load");
   load->Connect("Clicked()", "MyMainFrame", this, "DoLoad()");
   vframeButtons->AddFrame(load, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

   // Combos
   fComboDetector = new TGComboBox(vframeButtons);
   vframeButtons->AddFrame(fComboDetector, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));   
   for(Int_t idet = 0; idet < AliOADBPWG2Spectra::kNDetectors; idet++){
     fComboDetector->AddEntry(AliOADBPWG2Spectra::GetDetectorName(idet),idet);
   }
   fComboDetector->Resize(120,20);
   fComboDetector->Select(1);
   fComboDetector->Connect("Selected(Int_t)", "MyMainFrame", this, "DoSelectedDetector(Int_t)");

   fComboPID = new TGComboBox(vframeButtons);
   vframeButtons->AddFrame(fComboPID, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));   
   for(Int_t idet = 0; idet < AliOADBPWG2Spectra::kNPIDTypes; idet++){
     fComboPID->AddEntry(AliOADBPWG2Spectra::GetPIDName(idet),idet);
   }
   fComboPID->Resize(120,20);
   fComboPID->Select(0);

   fComboCharge = new TGComboBox(vframeButtons);
   vframeButtons->AddFrame(fComboCharge, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));   
   for(Int_t idet = 0; idet < AliOADBPWG2Spectra::kNCharge; idet++){
     fComboCharge->AddEntry(AliOADBPWG2Spectra::GetChargeName(idet),idet);
   }
   fComboCharge->Resize(120,20);
   fComboCharge->Select(0);

   fComboParticle = new TGComboBox(vframeButtons);
   vframeButtons->AddFrame(fComboParticle, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));   
   for(Int_t idet = 0; idet < AliOADBPWG2Spectra::kNParticle; idet++){
     fComboParticle->AddEntry(AliOADBPWG2Spectra::GetParticleName(idet),idet);
   }
   fComboParticle->Resize(120,20);
   fComboParticle->Select(0);

   fComboCentr = new TGComboBox(vframeButtons);
   vframeButtons->AddFrame(fComboCentr, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));   
   for(Int_t icentr = 0; icentr < fkNCentrTags; icentr++){
     fComboCentr->AddEntry(fkCentrTags[icentr],icentr);
   }
   fComboCentr->Resize(120,20);
   fComboCentr->Select(0);

   // Text fields
   fRunNumber; // run number 
   fTextCentrBin = new TGTextEntry (vframeButtons, "1");
   fTextCentrBin->Resize(120,20);
   vframeButtons->AddFrame(fTextCentrBin, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));   
   fTextDrawOpt = new TGTextEntry (vframeButtons, "auto");
   vframeButtons->AddFrame(fTextDrawOpt, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));   
   fTextDrawOpt->Resize(120,20);
   fTextSaveName = new TGTextEntry (vframeButtons, "spectra.png");
   vframeButtons->AddFrame(fTextSaveName, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));   
   fTextSaveName->Resize(120,20);

   // buttons
   TGTextButton *draw = new TGTextButton(vframeButtons, "&Draw");
   draw->Connect("Clicked()", "MyMainFrame", this, "DoDraw()");
   vframeButtons->AddFrame(draw, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));
   TGTextButton *legend = new TGTextButton(vframeButtons, "&Legend");
   legend->Connect("Clicked()", "MyMainFrame", this, "DoLegend()");
   vframeButtons->AddFrame(legend, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));
   TGTextButton *ratio = new TGTextButton(vframeButtons, "&Ratio");
   ratio->Connect("Clicked()", "MyMainFrame", this, "DoRatio()");
   vframeButtons->AddFrame(ratio, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));
   TGTextButton *clear = new TGTextButton(vframeButtons, "&Clear");
   clear->Connect("Clicked()", "MyMainFrame", this, "DoClear()");
   vframeButtons->AddFrame(clear, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));
   TGTextButton *print = new TGTextButton(vframeButtons, "&Print");
   print->Connect("Clicked()", "MyMainFrame", this, "DoPrint()");
   vframeButtons->AddFrame(print, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));
   TGTextButton *exit = new TGTextButton(vframeButtons, "&Exit ");
   exit->Connect("Pressed()", "MyMainFrame", this, "DoExit()");
   vframeButtons->AddFrame(exit, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));




   hframeMain->AddFrame(vframeButtons, new TGLayoutHints(kLHintsRight, 2, 2, 2, 2));

   AddFrame(hframeMain,new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX  | kLHintsExpandY, 2, 2, 2, 2));

   // Set a name to the main frame   
   SetWindowName("Spectra OADB Browser");
   MapSubwindows();

   // Initialize the layout algorithm via Resize()
   Resize(GetDefaultSize());

   // Map main frame
   MapWindow();
}


MyMainFrame::~MyMainFrame()
{
   // Clean up main frame...
   Cleanup();
   delete fEcan;
}


void spectraOADBGUI()
{
  // Popup the GUI...
  
  new MyMainFrame(gClient->GetRoot(), 200, 200);
}
