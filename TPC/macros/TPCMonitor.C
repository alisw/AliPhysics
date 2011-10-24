////////////////////////////////////////////////////////////////////////////////
// 
// Macro to launch Monitoring program for TPC raw data
// 
// Author: Stefan Kniege, IKF, Frankfurt
//
////////////////////////////////////////////////////////////////////////////////

#include "iostream"

#include "TGFrame.h"
#include "TGButton.h"
#include "TGTextEntry.h"
#include "TGLabel.h"
#include "TGFileDialog.h"
#include "TGClient.h"
#include "TGApplication.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TDirectory.h"
#include "TFileMerger.h"
#include "TGComboBox.h"

#include "AliTPCMonitor.h"
#include "AliTPCMonitorMappingHandler.h"
//#include "AliTPCMonitorEditor.h"
#include "AliTPCMonitorDialog.h"

using namespace std;


static    AliTPCMonitor*               fMon               = 0;
static    AliTPCMonitorMappingHandler* fMapHand           = 0;

static    TGMainFrame*                 fFrameMain         = 0;
static    TGCheckButton*               fFrameChDisFit     = 0;
static    TGCheckButton*               fFrameCh10bit      = 0;
static    TGCheckButton*               fFrameCheckVerb    = 0;
static    TGCheckButton*               fFrameCheckProcOne = 0;
static    TGCheckButton*               fFrameCheckPed     = 0;
static    TGCheckButton*               fFrameCalcBSL      = 0;           

static    TGComboBox*                  fTableField        = 0;

static    Int_t                        fVerb              = 0;
static    TGTextEntry*                 fTextEvId          = 0;

void      WriteChannel();
void      WriteHistos();  
void      SetConfig();
void      DisableFit();
//void      ProcessSector(char* fdata,char* ffil, Int_t side, Int_t sector);
void      ProcessSector(Int_t sid, Int_t sector);
//Int_t     DrawFit();
Int_t     DrawRMSMap();
void      ShowSelected();
//void      SetSize() ;
void      ResizeCanv();
void      ReadMe();
void      OpenDir();
void      SetWrite10Bit();
void      SetCheckVerb();
void      SetProcOne();
void      SetPedestalRun(Int_t val);
void      InitDialog(Int_t id);
//void      Resize(Int_t update,Int_t doit , Int_t side);
void      MonitorGui(AliTPCMonitor* fMon);
void      SetStyle();
 

//_________________________________________________________________________
void TPCMonitor()
{
  // Initialize the monitor 
  SetStyle();
  
  // MappingHandler
  Char_t fglobalmap[256] ; sprintf(fglobalmap,"%s/TPC/mapping/MappingGlobal.txt", gSystem->Getenv("ALICE_ROOT"));
  Char_t frowmap[256]    ; sprintf(frowmap,   "%s/TPC/mapping/MappingRow.txt",    gSystem->Getenv("ALICE_ROOT"));
  Char_t ffecmap[256]    ; sprintf(ffecmap,   "%s/TPC/mapping/MappingCards.txt",  gSystem->Getenv("ALICE_ROOT"));
  Char_t fnameconf[256]  ; sprintf(fnameconf, "%s/TPC/AliTPCMonitorConfig.txt",   gSystem->Getenv("ALICE_ROOT"));
  
  cout << " ALICE_ROOT : " << gSystem->Getenv("ALICE_ROOT") << endl;
  cout << " DATE_ROOT  : " << gSystem->Getenv("DATE_ROOT")  << endl;

  fMapHand = new AliTPCMonitorMappingHandler("maphand","maphand");
  fMapHand->ReadMapping(fglobalmap);
  fMapHand->ReadRowMappingGlob(frowmap);
  fMapHand->ReadFECMapping(    ffecmap);  
  if (gDirectory) { gDirectory->Append(fMapHand); }
  
  // Monitor
  fMon = new AliTPCMonitor("monitor","monitor");
  fMon->ReadConfig(fnameconf);
  fMon->SetMappingHandler(fMapHand);
  MonitorGui(fMon);
  
}

//_________________________________________________________________________
void MonitorGui(AliTPCMonitor */*fMon*/)
{
  // Display the main Window 

  Float_t xsize   = fMon->GetButtonXSize();
  Float_t ysize   = fMon->GetButtonYSize();
  Float_t xfirst1 = fMon->GetButtonXFirst1();
  Float_t xfirst2 = fMon->GetButtonXFirst2();
  Float_t yfirst  = fMon->GetButtonYFirst();
  Float_t mainx   = fMon->GetMainXSize();
  Float_t mainy   = fMon->GetMainYSize();
  
  // main frame
  fFrameMain = new TGMainFrame(gClient->GetRoot(),10,10,kMainFrame | kVerticalFrame);
  fFrameMain->SetLayoutBroken(kTRUE);

//  TGTextButton*  fFrameSelForm    = new TGTextButton(fFrameMain,  "Sel. Format "         );
  TGTextButton*  fFrameSelFil     = new TGTextButton(fFrameMain,  "Sel. File/Stream"     );
                 fTableField      = new TGComboBox(fFrameMain                            );
  TGTextButton*  fFrameFFT        = new TGTextButton(fFrameMain,  "FFT"                  );
  TGTextButton*  fFrameWRITE      = new TGTextButton(fFrameMain,  "Write Channel"        );
  fFrameChDisFit                  = new TGCheckButton(fFrameMain, "Disable G4-fit"       );
  TGTextButton*  fFrameRMS        = new TGTextButton(fFrameMain,  "Show RMS map"         );
  TGTextButton*  fFrameSetConf    = new TGTextButton(fFrameMain,  "Conf. Ranges"         );
  fFrameCh10bit                   = new TGCheckButton(fFrameMain, "Write 10bit      "    );
  fFrameCheckVerb                 = new TGCheckButton(fFrameMain, "Set Verbose      "    );
  fFrameCheckProcOne              = new TGCheckButton(fFrameMain, "Proc one Sector  "    );
  fFrameCalcBSL                   = new TGCheckButton(fFrameMain, "Calc BSL (onl)   "    );
  fFrameCheckPed                  = new TGCheckButton(fFrameMain, "No BSL sub.      "    );
  TGTextButton*  fFrameSideA      = new TGTextButton(fFrameMain,  "Side A"               );
  TGTextButton*  fFrameSideB      = new TGTextButton(fFrameMain,  "Side C"               );
  TGTextButton*  fFrameNextEvent  = new TGTextButton(fFrameMain,  "Next Event"           );
  TGTextButton*  fFramesel        = new TGTextButton(fFrameMain,  "Show Component"       );
  TGTextButton*  fFrameres        = new TGTextButton(fFrameMain,  "Resize Canvases"      );
  TGTextButton*  fFramewrite      = new TGTextButton(fFrameMain,  "Write Histos"         );
  TGTextButton*  fFrameQuit       = new TGTextButton(fFrameMain,  "Quit"                 );
 
//  fFrameSelForm->SetCommand(   "InitDialog(1)"     );
  fFrameSelFil->SetCommand(    "OpenDir()"        );
  fFrameFFT->SetCommand(       "FFT()"             );
  fFrameWRITE->SetCommand(     "WriteChannel()"    );
  fFrameCheckPed->SetCommand(  "SetPedestalRun(0)" );
  fFrameCh10bit->SetCommand(   "SetWrite10Bit()"   );
  fFrameCalcBSL->SetCommand(   "SetPedestalRun(1)" );
  fFrameChDisFit->SetCommand(  "DisableFit()"      );
  fFrameRMS->SetCommand(       "DrawRMSMap()"      );
  fFrameSetConf->SetCommand(   "SetConfig()"       );
  fFrameCheckVerb->SetCommand( "SetCheckVerb()"    );
  fFrameCheckProcOne->SetCommand("SetProcOne()"      );
  fFrameNextEvent->SetCommand( "NextEvent()"       );  fFrameNextEvent->SetTextColor(200);
  fFramesel->SetCommand(       "ShowSelected()"    );
  fFrameres->SetCommand(       "ResizeCanv()"      );
  fFramewrite->SetCommand(     "WriteHistos()"     );
  fFrameQuit->SetCommand(      "Quit()"            );

  
  fTableField->EnableTextInput(kTRUE);
  Int_t id=0;
  fTableField->RemoveAll();
  fTableField->AddEntry("PHY,Y,73,*",id++);
  fTableField->AddEntry("PHY,Y,*,*",id++);
  TString fields=gSystem->GetFromPipe("cat ~/TPCMonTable");
  TObjArray *arr=fields.Tokenize("\n");
  for (Int_t ifield=0; ifield<arr->GetEntries(); ++ifield) fTableField->AddEntry(arr->At(ifield)->GetName(),id++);
  fTableField->Select(0,kFALSE);
//   fTableField->Connect("ReturnPressed()", "", this, "DoCustomDraw()");
//   fTableField->Connect("Selected(Int_t)", "", this, "DoCustomDraw()");

  fFrameCalcBSL->SetDown(  1);
  fFrameCheckPed->SetDown( 1);
  fFrameChDisFit->SetDown( 1);
  fFrameCh10bit->SetDown(  0);
  fFrameCheckVerb->SetDown(0);
  fFrameCheckProcOne->SetDown(0);
  SetPedestalRun(0); // !!
  
  TGLayoutHints* fLayout = new TGLayoutHints(kLHintsLeft | kLHintsTop,5,5,5,5);
  fFrameMain->AddFrame(fFrameSideA     , fLayout);
  fFrameMain->AddFrame(fFrameSideB     , fLayout);
  fFrameMain->AddFrame(fFrameQuit      , fLayout);
  fFrameMain->AddFrame(fFrameCh10bit   , fLayout);
  fFrameMain->AddFrame(fFrameRMS       , fLayout);
  fFrameMain->AddFrame(fFrameSetConf   , fLayout);
  fFrameMain->AddFrame(fFrameSelFil    , fLayout);
  fFrameMain->AddFrame(fTableField     , fLayout);
//  fFrameMain->AddFrame(fFrameSelForm   , fLayout);
  fFrameMain->AddFrame(fFrameFFT       , fLayout);
//  Int_t step = (Int_t)(ysize/2.);
  Int_t start = 5; 
  
//  fFrameSelForm->MoveResize(     10, (Int_t)(start+ 1.0*ysize)        ,(UInt_t)(mainx-20)  ,(UInt_t)ysize);
  fFrameSelFil->MoveResize(      10, (Int_t)(start+ 1.0*ysize)        ,(UInt_t)(mainx-20)  ,(UInt_t)ysize);
 
  fFrameSetConf->MoveResize(     10, (Int_t)(start+ 2.0*ysize)        ,(UInt_t)(mainx-20)  ,(UInt_t)ysize);
  fTableField->MoveResize(       10, (Int_t)(start+ 3.0*ysize)        ,(UInt_t)(mainx-20)  ,(UInt_t)ysize+5);
 
  fFrameCalcBSL->MoveResize(     10, (Int_t)(start+   5*ysize)        ,(UInt_t)(mainx-20)  ,(UInt_t)ysize);
  fFrameCheckPed->MoveResize(    10, (Int_t)(start+   6*ysize)        ,(UInt_t)(mainx-20)  ,(UInt_t)ysize);
  fFrameChDisFit->MoveResize(    10, (Int_t)(start+   7*ysize)        ,(UInt_t)(mainx-20)  ,(UInt_t)ysize);
  fFrameCh10bit->MoveResize(     10, (Int_t)(start+   8*ysize)        ,(UInt_t)(mainx-20)  ,(UInt_t)ysize);
  fFrameCheckVerb->MoveResize(   10, (Int_t)(start+   9*ysize)        ,(UInt_t)(mainx-20)  ,(UInt_t)ysize);
  
  
  yfirst = start+ 14*ysize;
  fFrameSideA->MoveResize(  (Int_t)xfirst1, (Int_t)yfirst-30               ,(UInt_t)xsize     ,(UInt_t)ysize);
  fFrameSideB->MoveResize(  (Int_t)xfirst2, (Int_t)yfirst-30               ,(UInt_t)xsize     ,(UInt_t)ysize);
  
  // sector buttons 
  TObjArray     * fFrameArr   = new TObjArray();
  TGTextButton  * fTextButton = 0;
  Int_t side   = 0;  
  Int_t sector = 0;
  for(Int_t i = 0; i<36; i++)
    {
      char nameb[256] ; 
      if(i<18)sprintf(nameb,"Sector %i",i);
      else    sprintf(nameb,"Sector %i",i-18);
      fTextButton = new TGTextButton(fFrameMain,nameb);
      fFrameMain->AddFrame(fTextButton, new TGLayoutHints(kLHintsLeft | kLHintsTop,5,5,5,5));
      if(i<18)fTextButton->MoveResize((Int_t)xfirst1,(Int_t)(yfirst     +i*ysize),(UInt_t)xsize,(UInt_t)ysize);
      else    fTextButton->MoveResize((Int_t)xfirst2,(Int_t)(yfirst+(i-18)*ysize),(UInt_t)xsize,(UInt_t)ysize);
      if(i<18){ side = 0; sector = i;   }
      else    { side = 1; sector = i-18;}
      char bef[50]; sprintf(bef,"ProcessSector(%i,%i)",side,sector);
      fTextButton->SetCommand(bef);
      fFrameArr->Add(fTextButton);
    }
  
  TGLabel*      flab  = new TGLabel(fFrameMain, new TGHotString("Next Ev.ID"));
  TGTextBuffer* ftbuf = new TGTextBuffer(10); ftbuf->AddText(0, "1");
  fTextEvId = new TGTextEntry(fFrameMain, ftbuf);
  fTextEvId->SetTextColor(200);
  
  fFramesel->MoveResize(         10         , (Int_t)(mainy- 15.0*ysize)   ,(UInt_t)mainx-20  ,(UInt_t)ysize);
  fFrameRMS->MoveResize(         10         , (Int_t)(mainy- 14.0*ysize)   ,(UInt_t)mainx-20  ,(UInt_t)ysize);

  fFrameFFT->MoveResize(         10         , (Int_t)(mainy- 12.0*ysize)   ,(UInt_t)mainx-20  ,(UInt_t)ysize);
  fFramewrite->MoveResize(       10         , (Int_t)(mainy- 11.0*ysize)   ,(UInt_t)mainx-20  ,(UInt_t)ysize);
 
  fFrameWRITE->MoveResize(       10         , (Int_t)(mainy-  9.0*ysize)   ,(UInt_t)mainx-20  ,(UInt_t)ysize);
  fFrameres->MoveResize(         10         , (Int_t)(mainy-  8.0*ysize)   ,(UInt_t)mainx-20  ,(UInt_t)ysize);
  

  fFrameCheckProcOne->MoveResize(10         , (Int_t)(mainy-  6.5*ysize)   ,(UInt_t)mainx-20  ,(UInt_t)ysize);
  
  flab->MoveResize(             10          , (Int_t)(mainy-  5.0*ysize)   ,(UInt_t)xsize+5   ,(UInt_t)ysize);
  fTextEvId->MoveResize(      (Int_t)(mainx/2 +10)   , (Int_t)(mainy-  5.0*ysize)   ,(UInt_t)xsize-10  ,(UInt_t)ysize);
  
  fFrameNextEvent->MoveResize(   10         , (Int_t)(mainy-  3.5*ysize)   ,(UInt_t)mainx-20  ,(UInt_t)ysize);

  fFrameQuit->MoveResize(        30         , (Int_t)(mainy-  1.5*ysize)   ,(UInt_t)mainx-60  ,(UInt_t)ysize);
  
  fFrameMain->MapSubwindows();
  fFrameMain->MapWindow();
  fFrameMain->SetWindowName("OM");
  fFrameMain->MoveResize(0,0,(UInt_t)mainx,(UInt_t)mainy);
}  

//_________________________________________________________________________
void WriteHistos()  
{ 
  // Write Monitor Histos to file 
  fMon->WriteHistos() ;
} 

//_________________________________________________________________________
void WriteChannel()
{
  // Write 10 bit words to file for current channel
  fMon->Write10bitChannel();
}

//_________________________________________________________________________
void FFT()
{
  // Make Fourier Transformation for current channel
  fMon->ExecTransform();
}

//_________________________________________________________________________
void NextEvent()
{
  // Process next event
  fMon->SetProcNextEvent(1);
  fMon->SetupMonitoringTable(fTableField->GetTextEntry()->GetText());
  Int_t eventid =0;
  TString s1 =  fTextEvId->GetDisplayText();
  if(!s1.IsDigit()) { cout << " Invalid EventID  " << endl; return ;}
  else eventid = s1.Atoi();
  fMon->SetEventID(eventid);
  fMon->ResizeCanv();
  fMon->ProcessEvent();
  fMon->SetProcNextEvent(0);
  char tenttext[20] ; sprintf(tenttext,"%i",fMon->GetEventID()+1);
  fTextEvId->SetText(tenttext);
}

//_________________________________________________________________________
void SetCheckVerb()
{
  // Change verbose mode
  if(fFrameCheckVerb->IsDown())     {     fVerb = 1; fMon->SetVerbose(1); }
  else                              {     fVerb = 0; fMon->SetVerbose(0); }
}
//_________________________________________________________________________
void SetPedestalRun(Int_t val)
{
  // Set Pedestal calculation mode 
  if(val==0)
    {
      // check pedestal run button
      if(fFrameCheckPed->IsDown())  {	  fMon->SetPedestals(0);	  fFrameCalcBSL->SetDown(0);	}
      else	                    { 	  fMon->SetPedestals(1);	  fFrameCalcBSL->SetDown(1);	}
    }
  else  
    {
      // check online calc button
      if(fFrameCalcBSL->IsDown())   { 	  fMon->SetPedestals(1);	  fFrameCheckPed->SetDown(0);	}
      else	                    { 	  fMon->SetPedestals(2);	  fFrameCheckPed->SetDown(0);	}
    } 
}

//_________________________________________________________________________
void SetWrite10Bit()
{
  // Set Write 10Bit mode
  // 10 bit words will be written in AliTPCMonitorAltro for each equipment
  if(fMon->GetWrite10Bit()==1) fMon->SetWrite10Bit(0);
  else                         fMon->SetWrite10Bit(1);
}

//_________________________________________________________________________
void DisableFit()
{
  // Disable Gamma4 fit to maximum peak
  if(fMon->GetFitPulse())      fMon->SetFitPulse(0);
  else                         fMon->SetFitPulse(1);
}

//_________________________________________________________________________
Int_t DrawRMSMap()
{
  // Draw RMS map for IROC and OROC
  if(fMon->GetLastSector()==-1) cout << " no sector written yet " << endl;
  else fMon->DrawRMSMap();
  return 0;
}

//_________________________________________________________________________
void ProcessSector(Int_t sid, Int_t sec)
{
  // Process specific sector. Initiated by button
  if(fMon->GetSectorFilled(sec,sid))
    {
      fMon->SetLastSector(sec +sid*18);
      fMon->SetProcNextEvent(0);
      fMon->ProcessEvent();
    }
  else 
    {
      cout << " Sector not filled " << endl;
    }
}

//_________________________________________________________________________
void OpenDir()
{
  TString dir ;
  TString lastfile(fMon->GetLastProcFile());
  if ( lastfile == "" ) dir=gSystem->pwd();
  else dir=gSystem->DirName(lastfile);
  dir="/";
  TGFileInfo* fi = new TGFileInfo();
  fi->fIniDir    = StrDup(dir);
  fi->fFilename  = StrDup(lastfile);
  if(fVerb) printf("fIniDir = %s\n", fi->fIniDir);
  
  new TGFileDialog(gClient->GetRoot(), fFrameMain, kFDOpen, fi);
  if(!fi->fFilename)            return;
  if(!strcmp(fi->fFilename,"")) return;
  string fname(fi->fFilename);
  Int_t    ffirst   = fname.find_first_not_of("/",0);
  Int_t    fstart   = 0;  
  string   firsts   = fname.substr(ffirst,1);
  Int_t    firstsl  = strcmp(firsts.data(),":");
  Int_t    firstcol = strcmp(firsts.data(),"@");
  if ( (fname.find("mem:")!=string::npos) || (fname.find("rfio")!=string::npos) || (fname.find("http")!=string::npos)
       || !firstsl || !firstcol  )  fstart= ffirst  ;
  else                              fstart= ffirst-1; 
  
  if(fstart <0){ cout << " return : missing slash at beginning of file " << endl ; return ;}
  
  string fsubname =  fname.substr(fstart,fname.length()-fstart);
  fMon->SetFile((Char_t *)fsubname.data());

  if(fVerb) 
    {
      cout << "Monitor.C OpenDir:: dir    : " << fi->fIniDir       << endl; 
      cout << "Monitor.C OpenDir:: file   : " << fsubname.data()   << endl;
      cout << "Monitor.C OpenDir:: format : " << fMon->GetFormat() << endl;
    }
}

//_________________________________________________________________________
void InitDialog(Int_t id)
{
  if(id>0&&id<4)  new AliTPCMonitorDialog((TGWindow*)gClient->GetRoot(), fFrameMain, 400, 400,0,id,fMon);
  else cout << "Error : Invalid Id " << endl;
  return;
}

//_________________________________________________________________________
void ShowSelected()
{
  InitDialog(2);
  return;
}

//_________________________________________________________________________
void SetConfig()
{
  InitDialog(3);
}

//_________________________________________________________________________
void ResizeCanv()
{
  fMon->ResizeCanv();
  return;
}

//_________________________________________________________________________
void SetStyle() 
{
  gStyle->SetScreenFactor(1);
  gStyle->SetPadTopMargin(0.17);
  gStyle->SetPadBottomMargin(0.17);
  gStyle->SetPadLeftMargin(0.17);
  gStyle->SetPadRightMargin(0.19);
  gStyle->SetStatColor(kWhite);
  gStyle->SetPadColor(kWhite);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadBorderSize(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasBorderSize(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(1);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleOffset(1.3,"X");
  gStyle->SetTitleOffset(1.9,"Y");
  gStyle->SetTitleOffset(1.7,"Z");
  gStyle->SetPalette(1);

}

//_________________________________________________________________________
void ReadMe()
{
//  AliTPCMonitorEditor *ed = new AliTPCMonitorEditor(fFrameMain, 700, 400);
//  char nameread[256]; sprintf(nameread,"%s/TPC/AliTPCMonitorReadMe.txt",gSystem->Getenv("ALICE_ROOT"));
//  ed->LoadFile(nameread);
//  ed->Popup();
}

//_________________________________________________________________________
void SetProcOne()
{
  if(fFrameCheckProcOne->IsDown())  { fMon->SetProcOneSector(1); }
  else                              { fMon->SetProcOneSector(0); }
  
}
//_________________________________________________________________________
void Quit() 
{
  gApplication->Terminate(0);
}

