#include "TMCal.h"
#include "TRandom.h"
#include "TF1.h"

#include "AliRawEventHeaderBase.h"
#include "AliRawReaderRoot.h"

ClassImp(TMCal)
  
//_____________________________________________________________________________
TMCal::TMCal(const TGWindow *p, UInt_t w, UInt_t h, AliCaloCalibPedestal::kDetType detectorType) : 
  TMBaseModule(p, w, h),
  fDetType(detectorType),
  fEvents(0),
  fPedestals(0),
  fRawReader(0),
  fCaloRawStream(0),
  fLowGainMode(kFALSE),
  fSingleModule(kTRUE),
  fVisibleModule(2),
  fMaxPed(100.0),
  fMaxPeak(100.0),
  fMaxPeakDiff(100.0),
  fMaxPedDiff(100.0),
  fMaxPedRatio(2.0),
  fMaxPeakRatio(3.0)
{
  //PHOS/EMCAL MOOD Module; 
  // GUI related pointers (frames, canvases, buttons are initialized in the Construct* methods
  // called later in this ctor.
  
  //Create the object for making the histograms
  fPedestals = new AliCaloCalibPedestal( fDetType ); 
  // AliCaloCalibPedestal knows how many modules we have for PHOS or EMCAL
  fNumModules = fPedestals->GetModules(); 

  //Load the reference object (Note: a fixed path. This should be made a variable that can somehow be set other than via the constants file and recompiled..)
  if (fDetType == AliCaloCalibPedestal::kPhos) { 
    if (fPedestals->LoadReferenceCalib(fgkReferenceFilePhos, fgkReferenceObject) != kTRUE) {
      printf("Could not load PHOS reference calib.\n");
    }
  }
  else {
    if (fPedestals->LoadReferenceCalib(fgkReferenceFileEmCal, fgkReferenceObject) != kTRUE) {
      printf("Could not load EMCAL reference calib.\n");
    }
  }
   
  //Init the readers and streams
  InitRawReaderAndStream(fgkDefaultRunNo);
  
  // standard style options:  
  //Make all histos use a more sensible palette.
  gStyle->SetPalette(1);
  //Kill the titles etc...
  gStyle->SetOptStat(0);
  //gStyle->SetOptTitle(0);
  gStyle->SetTitleW(0.8);
  gStyle->SetTitleH(0.08);
  gStyle->SetTitleX(0.1);
  gROOT->ForceStyle(kFALSE);
  
  //Construct the various buttons and info 'menu' for the lower panel
  ConstructMenu(w);
  
  //Construct the various tabs and frames
  ConstructFullGUI();
  
  //Set the handler to handle the tab changing, so we can enforce the palette we like.
#ifdef fgkCustomDeadMapPalette
  GetTab()->Connect("Selected(Int_t)", "TMCal", this, "OnTabChange(Int_t)");
#endif
}

//_____________________________________________________________________________
void TMCal::InitRawReaderAndStream(unsigned int runNo)
{
  char file[256];  
  if (fDetType == AliCaloCalibPedestal::kPhos) { 
    snprintf(file, 256, fgkSourceFileTemplatePhos, runNo);
  }
  else {
    snprintf(file, 256, fgkSourceFileTemplateEmCal, runNo);
  }
  printf("\n TMCal::InitRawReaderAndStream file= %s \n", file);
 
  fEvents = 0;
  
  //Initialize the raw readers
  if (fCaloRawStream) {
    delete fCaloRawStream;
  }
  if (fRawReader) {
    delete fRawReader;
  }

  fRawReader = new AliRawReaderRoot(file);
  fCaloRawStream = new AliCaloRawStream(fRawReader, (fDetType==AliCaloCalibPedestal::kPhos) ? "PHOS" : "EMCAL");

  fCaloRawStream->SetOldRCUFormat(kTRUE); // should maybe have a switch for this later?
  
  //Not sure if this should be here... get to the first event, or never look at the first one?
  fRawReader->NextEvent();
 
  //Reset the actual analysis, since we changed the source file
  fPedestals->Reset();
 
  fPedestals->SetRunNumber(runNo);
 
  //Set the histo titles
  SetTitles();
  printf("Reference run #%i.\n", fPedestals->GetRefRunNumber());
}

//_____________________________________________________________________________
void TMCal::SetTitles()
{
  char title[256];
  for (int i = 0; i < fNumModules; i++) {
    sprintf(title, "Pedestals, high gain, run #%.4i, module %i", fPedestals->GetRunNumber(), i);
    fPedestals->GetPedProfileHighGain(i)->SetTitle(title);
    sprintf(title, "Peak-Pedestals, high gain, run #%.4i, module %i", fPedestals->GetRunNumber(), i);
    fPedestals->GetPeakProfileHighGain(i)->SetTitle(title);
    
    sprintf(title, "Pedestals, low gain, run #%.4i, module %i", fPedestals->GetRunNumber(), i);
    fPedestals->GetPedProfileLowGain(i)->SetTitle(title);
    sprintf(title, "Peak-Pedestals, low gain, run #%.4i, module %i", fPedestals->GetRunNumber(), i);
    fPedestals->GetPeakProfileLowGain(i)->SetTitle(title);
    
    sprintf(title, "Peak-Pedestals difference, high gain, run #%.4i, module %i", fPedestals->GetRunNumber(), i);
    fPedestals->GetPeakProfileHighGainDiff(i)->SetTitle(title);
    sprintf(title, "Pedestals difference, high gain, run #%.4i, module %i", fPedestals->GetRunNumber(), i);
    fPedestals->GetPedProfileHighGainDiff(i)->SetTitle(title);
    
    sprintf(title, "Peak-Pedestals difference, low gain, run #%.4i, module %i", fPedestals->GetRunNumber(), i);
    fPedestals->GetPeakProfileLowGainDiff(i)->SetTitle(title);
    sprintf(title, "Pedestals difference, low gain, run #%.4i, module %i", fPedestals->GetRunNumber(), i);
    fPedestals->GetPedProfileLowGainDiff(i)->SetTitle(title);
    
    sprintf(title, "Peak-Pedestals ratio, high gain, run #%.4i, module %i", fPedestals->GetRunNumber(), i);
    fPedestals->GetPeakProfileHighGainRatio(i)->SetTitle(title);
    sprintf(title, "Pedestals ratio, high gain, run #%.4i, module %i", fPedestals->GetRunNumber(), i);
    fPedestals->GetPedProfileHighGainRatio(i)->SetTitle(title);
    
    sprintf(title, "Peak-Pedestals ratio, low gain, run #%.4i, module %i", fPedestals->GetRunNumber(), i);
    fPedestals->GetPeakProfileLowGainRatio(i)->SetTitle(title);
    sprintf(title, "Pedestals ratio, low gain, run #%.4i, module %i", fPedestals->GetRunNumber(), i);
    fPedestals->GetPedProfileLowGainRatio(i)->SetTitle(title);

  } 

  // Note: Dead map is currently only for high gain.. and we try to use a special palette style for it
  gStyle->SetPalette(AliCaloCalibPedestal::kNumDeadMapStates, fDeadMapPalette);
  for (int i = 0; i < fNumModules; i++) {
    sprintf(title, "Dead map, high gain, run #%.4i, module %i", fPedestals->GetRunNumber(), i);
    fPedestals->GetDeadMap(i)->SetTitle(title);
    fPedestals->GetDeadMap(i)->UseCurrentStyle();
  }
  // back to standard palette
  gStyle->SetPalette(1);
  return;
}

//_____________________________________________________________________________
void TMCal::ConstructMenu(UInt_t w)
{//Constructs the buttons and info 'menu'

  //create frame for the High gain/Low gain -button
  TGHorizontalFrame *buttonFrame = new TGHorizontalFrame(GetMainViewFrame(),(UInt_t)w,(UInt_t)30);
  GetMainViewFrame()->AddFrame(buttonFrame,new TGLayoutHints(kLHintsTop|kLHintsExpandX,0,0,0,0));
  
  //Create the high/lowgain button
  fGainButton = new TGTextButton(buttonFrame,"Low gain");
  fGainButton->Connect("Clicked()","TMCal",this,"OnGainButton()");
  //fGainButton->AllowStayDown(kTRUE);
  buttonFrame->AddFrame(fGainButton,new TGLayoutHints(kLHintsCenterY,0,(UInt_t)5,0,0));
  
  fDeadCountText = new TGLabel(buttonFrame,"Run #0000: Dead: 000, 000%; 000 new, 0000 recov.; ref. run 0000");
  //We set generous margins for the label to make sure the text fits even with 4-digit numbers. 
  fDeadCountText->SetMargins(fDeadCountText->GetLeftMargin()+fgkDeadCountTextMargin, fDeadCountText->GetRightMargin()+fgkDeadCountTextMargin, fDeadCountText->GetTopMargin(), fDeadCountText->GetBottomMargin());
  buttonFrame->AddFrame(fDeadCountText, new TGLayoutHints(kLHintsCenterY, 0, (UInt_t)5, 0, 0));
  
  //Create the one histo/multi-histo -selection buttons
  TGButtonGroup * buttGr = new TGButtonGroup(buttonFrame, "Display module", kHorizontalFrame);
  TGRadioButton * radButts[fNumModules + 1];
  char txt[64];
  for (int i = 0; i < fNumModules; i++) {
    sprintf(txt, "%i", i);
    radButts[i] = new TGRadioButton(buttGr, new TGHotString(txt));
    sprintf(txt, "OnModeButton(=%i)", i);
    radButts[i]->Connect("Clicked()", "TMCal", this, txt);
  }
  
  radButts[fNumModules] = new TGRadioButton(buttGr, new TGHotString("All"));
  sprintf(txt, "OnModeButton(=%i)", fNumModules);
  radButts[fNumModules]->Connect("Clicked()", "TMCal", this, txt);
  
  radButts[fVisibleModule]->SetState(kButtonDown);
  buttonFrame->AddFrame(buttGr,new TGLayoutHints(kLHintsCenterY,0,(UInt_t)5,0,0));
  //buttGr->Show();
  
  //Add the run number entry
  TGLabel * runnoText = new TGLabel(buttonFrame,"Run #: ");;
  buttonFrame->AddFrame(runnoText, new TGLayoutHints(kLHintsCenterY, 0, (UInt_t)5, 0, 0));
  fRunNoEntry = new TGNumberEntry(buttonFrame, fgkDefaultRunNo, 6, -1,TGNumberFormat::kNESInteger, //style
  					TGNumberFormat::kNEAPositive//input value filter
				  );
  fRunNoEntry->Connect("ValueSet(Long_t)", "TMCal",this,"OnRunNoChange(Long_t)");
  buttonFrame->AddFrame(fRunNoEntry,new TGLayoutHints(kLHintsCenterY|kLHintsLeft,0,(UInt_t)5,0,0));
  
  return;
}

//_____________________________________________________________________________
void TMCal::ConstructFullGUI()
{//Constructs the GUI
  
  //We need to set the palette for drawing the dead map
  fDeadMapPalette[AliCaloCalibPedestal::kAlive] = fgkLiveTowerColor;
  fDeadMapPalette[AliCaloCalibPedestal::kDead] = fgkDeadTowerColor;
  fDeadMapPalette[AliCaloCalibPedestal::kResurrected] = fgkResurrectedTowerColor;
  fDeadMapPalette[AliCaloCalibPedestal::kRecentlyDeceased] = fgkRecentlyDeceasedTowerColor;
  
  // set up the TPaveText to have the info on the palette also
  fDeadMapPaveText = new TPaveText(0.9,0.7,0.995,0.9,"brNDC");
  /* alive ones are not shown really..  
      TText *tAlive = fDeadMapPaveText->AddText("OK");
      tAlive->SetTextColor(fgkLiveTowerColor);
  */
  TText *tDead = fDeadMapPaveText->AddText("Dead");
  TText *tResurrected = fDeadMapPaveText->AddText("New OK");
  TText *tRecentlyDeceased = fDeadMapPaveText->AddText("New Dead");
  tDead->SetTextColor(fgkDeadTowerColor);
  tResurrected->SetTextColor(fgkResurrectedTowerColor);
  tRecentlyDeceased->SetTextColor(fgkRecentlyDeceasedTowerColor);

  fPeaksTab = GetTab()->AddTab("Peak-Pedestal");
  fPeaksFrameSingle = MakeSingleModuleTab(fPeaksTab, &fPeaksCanvasSingle, &fPeaksRangeSingle, kPeaksSingle, AliCaloCalibPedestal::GetSampleMin(),
					  AliCaloCalibPedestal::GetSampleMax(), fMaxPeak);
  fPeaksFrame = MakeAllModulesTab(fPeaksTab, fPeaksCanvas, &fPeaksRange, kPeaks, AliCaloCalibPedestal::GetSampleMin(),
  				  AliCaloCalibPedestal::GetSampleMax(), fMaxPeak);
  fPeaksCanvas[fNumModules] = fPeaksCanvasSingle;
  
  fPedsTab = GetTab()->AddTab("Pedestal");
  fPedsFrameSingle = MakeSingleModuleTab(fPedsTab, &fPedsCanvasSingle, &fPedsRangeSingle, kPedsSingle, AliCaloCalibPedestal::GetSampleMin(),
  				  AliCaloCalibPedestal::GetSampleMax(), fMaxPed);
  fPedsFrame = MakeAllModulesTab(fPedsTab, fPedsCanvas, &fPedsRange, kPeds, AliCaloCalibPedestal::GetSampleMin(),
				 AliCaloCalibPedestal::GetSampleMax(), fMaxPed); //Create the all modules version of the tab
  fPedsCanvas[fNumModules] = fPedsCanvasSingle;
  
  
  fPeaksDiffTab = GetTab()->AddTab("Peak diff");
  fPeaksDiffFrameSingle = MakeSingleModuleTab(fPeaksDiffTab, &fPeaksDiffCanvasSingle, &fPeaksDiffRangeSingle, kPeaksDiffSingle, AliCaloCalibPedestal::GetSampleMin(),
					      AliCaloCalibPedestal::GetSampleMax(), fMaxPeakDiff);
  fPeaksDiffFrame = MakeAllModulesTab(fPeaksDiffTab, fPeaksDiffCanvas, &fPeaksDiffRange, kPeaksDiff, AliCaloCalibPedestal::GetSampleMin(),
				      AliCaloCalibPedestal::GetSampleMax(), fMaxPeakDiff);
  fPeaksDiffTab->MapSubwindows();
  fPeaksDiffCanvas[fNumModules] = fPeaksDiffCanvasSingle; //Copy the pointer to the array, so we can just loop over this one too
  
  
  fPedsDiffTab = GetTab()->AddTab("Pedestal diff");
  fPedsDiffFrameSingle = MakeSingleModuleTab(fPedsDiffTab, &fPedsDiffCanvasSingle, &fPedsDiffRangeSingle, kPedsDiffSingle, AliCaloCalibPedestal::GetSampleMin(),
					     AliCaloCalibPedestal::GetSampleMax(), fMaxPedDiff);
  fPedsDiffFrame = MakeAllModulesTab(fPedsDiffTab, fPedsDiffCanvas, &fPedsDiffRange, kPedsDiff, AliCaloCalibPedestal::GetSampleMin(),
  				     AliCaloCalibPedestal::GetSampleMax(), fMaxPedDiff);
  fPedsDiffCanvas[fNumModules] = fPedsDiffCanvasSingle; 
  
  fPeaksRatioTab = GetTab()->AddTab("Peak ratio");
  fPeaksRatioFrameSingle = MakeSingleModuleTab(fPeaksRatioTab, &fPeaksRatioCanvasSingle, &fPeaksRatioRangeSingle, kPeaksRatioSingle,
					       fgkMinMaxRatioExp*fgkRatioStepScale, fgkMaxMaxRatioExp*fgkRatioStepScale,
					       -(log(fMaxPeakRatio)/log(2))*fgkRatioStepScale + fgkMaxMaxRatioExp*fgkRatioStepScale);
  fPeaksRatioFrame = MakeAllModulesTab(fPeaksRatioTab, fPeaksRatioCanvas, &fPeaksRatioRange, kPeaksRatio,
  				       fgkMinMaxRatioExp*fgkRatioStepScale, fgkMaxMaxRatioExp*fgkRatioStepScale,
				       -(log(fMaxPeakRatio)/log(2))*fgkRatioStepScale + fgkMaxMaxRatioExp*fgkRatioStepScale);
  fPeaksRatioCanvas[fNumModules] = fPeaksRatioCanvasSingle;
  
  fPedsRatioTab = GetTab()->AddTab("Pedestal ratio");
  fPedsRatioFrameSingle = MakeSingleModuleTab(fPedsRatioTab, &fPedsRatioCanvasSingle, &fPedsRatioRangeSingle, kPedsRatioSingle, 
					      fgkMinMaxRatioExp*fgkRatioStepScale, fgkMaxMaxRatioExp*fgkRatioStepScale,
					      -(log(fMaxPedRatio)/log(2))*fgkRatioStepScale + fgkMaxMaxRatioExp*fgkRatioStepScale);
  fPedsRatioFrame = MakeAllModulesTab(fPedsRatioTab, fPedsRatioCanvas, &fPedsRatioRange, kPedsRatio, 
                                      fgkMinMaxRatioExp*fgkRatioStepScale, fgkMaxMaxRatioExp*fgkRatioStepScale, 
				      -(log(fMaxPedRatio)/log(2))*fgkRatioStepScale + fgkMaxMaxRatioExp*fgkRatioStepScale);
  fPedsRatioCanvas[fNumModules] = fPedsRatioCanvasSingle;
  
  fDeadMapTab = GetTab()->AddTab(fgkDeadMapTabName);
  fDeadMapFrameSingle = MakeSingleModuleTab(fDeadMapTab, &fDeadMapCanvasSingle);
  fDeadMapFrame = MakeAllModulesTab(fDeadMapTab, fDeadMapCanvas);
  fDeadMapCanvas[fNumModules] = fDeadMapCanvasSingle;
  
  DrawHistos();
  
  UpdateMonitors();//This is actually just to update the status text... which should then perhaps be done elsewhere, so we could call that separately,
  		   //but since it's needed in exactly one place, I'm too lazy...
  
  return;  
}

//_____________________________________________________________________________
TGCompositeFrame * TMCal::MakeAllModulesTab(TGCompositeFrame * tab,
					    TCanvas **canvArray, //Array of pointers to the canvases to which the histos have been drawn. MUST be initialized.
					    TGVSlider **slider, //The pointer to the pointer to the slider, which to put the pointer to the slider to
					    kTab tabId, //The tab Id which the slider should give as a parameter when connecting
					    double sliderMin, double sliderMax, double sliderPos //The range and position of the slider
					    )
{//The function to create a tab full of histos
  
  //The frames. The magic numbers here don't really matter, because the ExpandX and ExpandY
  //hints will override them anyway. If this behaviour is modified, then these numbers
  //should be properly declared as constants.
  TGHorizontalFrame * verFrame = new TGHorizontalFrame(tab, (UInt_t)5*400, (UInt_t)600);
  TGCompositeFrame * histoFrame = new TGCompositeFrame(verFrame,(UInt_t)40, (UInt_t)400);
  
  //Make the layout manager a matrix style
  histoFrame->SetLayoutManager(new TGMatrixLayout(histoFrame, 0, fgkHistosPerRow));
  
  //Add a vertical frame to stack the high/low -gain histos on top of each others
  tab->AddFrame(verFrame, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
  
  
  if (slider) {//Add the slider to the histoframe
    TGVSlider * range = new TGVSlider(verFrame, fgk2DHistoHeight, kSlider1);
    range->SetRange((int)sliderMin, (int)sliderMax);
    range->SetPosition( (int)(sliderMax - sliderPos) ); //We have to invert this by hand
    char str[64];
    sprintf(str, "OnRange(=%u)", (int)tabId);
    //range->Connect("PositionChanged(Int_t)", "TMCal",this,str);
    range->Connect("Released()", "TMCal",this,str);
    verFrame->AddFrame(range,new TGLayoutHints(0,0,fgk2DHistoMarginX,0,0));
    *slider = range;
  }
  
  //Add the horizontal frames which layout the individual histos for each module
  verFrame->AddFrame(histoFrame, new TGLayoutHints(kLHintsExpandX|kLHintsExpandY, 0, 0, 0, 0));
  
  //-------------------------------------
  //Create the canvases
    
  TString canvasName = tab->GetName();
  for (int i = 0; i < fNumModules; i++) {
  
    TRootEmbeddedCanvas *ECanvas = new TRootEmbeddedCanvas(canvasName + " " + (TString)i,histoFrame,fgk2DHistoWidth,fgk2DHistoHeight);
    
  
    histoFrame->AddFrame(ECanvas, new TGLayoutHints(0, fgk2DHistoMarginX,0, fgk2DHistoMarginY, 0));
    
    canvArray[i] = ECanvas->GetCanvas();
  } 
  
  return verFrame; 
}

//_____________________________________________________________________________
TGCompositeFrame * TMCal::MakeSingleModuleTab(TGCompositeFrame * tab,
					      TCanvas **canvas,
					      TGVSlider **slider, //The pointer to the pointer to the slider, which to put the pointer to the slider to
					      kTab tabId, //The tab Id which the slider should give as a parameter when connecting
					      double sliderMin, double sliderMax, double sliderPos //The range and position of the slider
					      )
{//The function to create a tab with a single big histo.
  
  //TODO: consider if we really need the vertical frame at all, or should just arrange everything to the horizontal frame...
  
  //The frames. The magic numbers here don't really matter, because the ExpandX and ExpandY
  //hints will override them anyway. If this behaviour is modified, then these numbers
  //should be properly declared as constants.
  TGVerticalFrame * verFrame = new TGVerticalFrame(tab, (UInt_t)5*400, (UInt_t)600);
  TGHorizontalFrame * histoFrame = new TGHorizontalFrame(verFrame,(UInt_t)400, (UInt_t)400);
  
  //Add a vertical frame to stack the high/low -gain histos on top of each others
  tab->AddFrame(verFrame, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
  
  //Add the horizontal frames which layout the individual histos for each module
  verFrame->AddFrame(histoFrame, new TGLayoutHints(kLHintsExpandX|kLHintsExpandY, 0, 0, 0, 0));
  
  if (slider) {//Add the slider to the histoframe
    TGVSlider * range = new TGVSlider(histoFrame, fgk2DSingleHistoHeight, kSlider1);
    range->SetRange((int)sliderMin, (int)sliderMax);
    range->SetPosition( (int)(sliderMax - sliderPos) ); //We have to invert this by hand
    char str[64];
    sprintf(str, "OnRange(=%u)", (int)tabId);
    //range->Connect("PositionChanged(Int_t)", "TMCal",this,str);
    range->Connect("Released()", "TMCal",this,str);
    histoFrame->AddFrame(range,new TGLayoutHints(0,0,0,0,0));
    *slider = range;
  }
  
  //-------------------------------------
  //Create the canvas  
  
  TString canvasName = tab->GetName();
  
  TRootEmbeddedCanvas *ECanvas = new TRootEmbeddedCanvas(canvasName + "  single",histoFrame,fgk2DSingleHistoWidth,fgk2DSingleHistoHeight);
  
  histoFrame->AddFrame(ECanvas, new TGLayoutHints(0, fgk2DSingleHistoMarginX,0, fgk2DSingleHistoMarginY, 0));
    
  *canvas = ECanvas->GetCanvas();
  
  return verFrame; 
}

//_____________________________________________________________________________
void TMCal::MonitorEvent(void)
{
// printf("Monitoring for the %ith time this evening.\n", i);
 //Process the next event
 if (fRawReader->NextEvent())
 {
  AliRawEventHeaderBase *aliHeader = (AliRawEventHeaderBase*) fRawReader->GetEventHeader();
  
  if (aliHeader->Get("Type") == AliRawEventHeaderBase::kPhysicsEvent)
  {
     fEvents++;
    //printf("A physics event.\n");
    fPedestals->ProcessEvent(fCaloRawStream);
  } //else
   // printf("A non-physics event.\n");
 }
 
}

//_____________________________________________________________________________
void TMCal::UpdateMonitors(void)
{
  //Just update all the histograms
  fPedestals->ComputeDiffAndRatio();//Re-compute the difference and the ratio
  //Update the dead cells count
  fPedestals->ComputeDeadTowers(fgkDeadThreshold);
  for (int i = 0; (i < fNumModules) && (!fSingleModule); i++) {
    fPeaksCanvas[i]->Paint();
    fPeaksCanvas[i]->Update();
    fPeaksCanvas[i]->Paint();
    
    fPedsCanvas[i]->Paint();
    fPedsCanvas[i]->Update();
    fPedsCanvas[i]->Paint();
    
    fPedsDiffCanvas[i]->Paint();
    fPedsDiffCanvas[i]->Update();
    fPedsDiffCanvas[i]->Paint();
    
    fPeaksDiffCanvas[i]->Paint();
    fPeaksDiffCanvas[i]->Update();
    fPeaksDiffCanvas[i]->Paint();
    
    fPedsRatioCanvas[i]->Paint();
    fPedsRatioCanvas[i]->Update();
    fPedsRatioCanvas[i]->Paint();
    
    fPeaksRatioCanvas[i]->Paint();
    fPeaksRatioCanvas[i]->Update();
    fPeaksRatioCanvas[i]->Paint();
    
    fPedestals->GetDeadMap(i)->UseCurrentStyle();
    fDeadMapCanvas[i]->Paint();
    fDeadMapCanvas[i]->Update();
    fDeadMapCanvas[i]->Paint();
  }

  if (fSingleModule) {
    fPeaksCanvasSingle->Paint();
    fPeaksCanvasSingle->Update();
    fPeaksCanvasSingle->Paint();
    
    fPeaksCanvasSingle->Paint();
    fPedsCanvasSingle->Update();
    fPedsCanvasSingle->Paint();
    
    fPeaksCanvasSingle->Paint();
    fPedsDiffCanvasSingle->Update();
    fPedsDiffCanvasSingle->Paint();
    
    fPeaksCanvasSingle->Paint();
    fPeaksDiffCanvasSingle->Update();
    fPeaksDiffCanvasSingle->Paint();
    
    fPeaksCanvasSingle->Paint();
    fPedsRatioCanvasSingle->Update();
    fPedsRatioCanvasSingle->Paint();
    
    fPeaksCanvasSingle->Paint();
    fPeaksRatioCanvasSingle->Update();
    fPeaksRatioCanvasSingle->Paint();
    
    fPedestals->GetDeadMap(fVisibleModule)->UseCurrentStyle();
    fDeadMapCanvasSingle->Paint();
    fDeadMapCanvasSingle->Update();
    fDeadMapCanvasSingle->Paint();  
 }
  
  char tmp[128];
  sprintf(tmp, "Run #%.4i: Dead: %i, %.1f%%; %i new, %i recov.; Ref. run %i", 
	  fPedestals->GetRunNumber(), (int)fPedestals->GetDeadTowerCount(), 
	  (float)(100.0*fPedestals->GetDeadTowerRatio()),
	  (int)fPedestals->GetDeadTowerNew(), 
	  (int)fPedestals->GetDeadTowerResurrected(), (int)fPedestals->GetRefRunNumber());

  fDeadCountText->SetText(tmp);
 
  return; 
}

//_____________________________________________________________________________
void TMCal::PostMonitor(void)
{
  fPedestals->ComputeDeadTowers(fgkDeadThreshold, fgkDeadMapName);//Write out the deadmap
  return;
}

//_____________________________________________________________________________
TMCal::~TMCal() 
{ // Destructor  
  if (fCaloRawStream) delete fCaloRawStream; 
  if (fRawReader) delete fRawReader;
  if (fPedestals) delete fPedestals;  
}

//_____________________________________________________________________________
void TMCal::DrawHistos() { // Gain button handler
  //Recalculate the differences and ratios (because this causes a screen update)
  fPedestals->ComputeDiffAndRatio();
 
  //Finally draw the deadmap. This has to be done outside the low/highgain iffing.
  //There is probably a smoother way to do this than by repeating the loop, though...
  for (int i = 0; (i < fNumModules) && (!fSingleModule); i++) {
    fDeadMapCanvas[i]->cd();
    //Draw the deadmap
    fPedestals->GetDeadMap(i)->Draw("col");
    // also a TPaveText with info on what the colors mean
    fDeadMapPaveText->Draw();
    fDeadMapCanvas[i]->Update();
    fDeadMapCanvas[i]->Paint();
  }
 
 if (fSingleModule) {
   fDeadMapCanvasSingle->cd();
   fPedestals->GetDeadMap(fVisibleModule)->Draw("col");
    // also a TPaveText with info on what the colors mean
   fDeadMapPaveText->Draw();
   fDeadMapCanvasSingle->Update();
   fDeadMapCanvasSingle->Paint();
 }
 
 if (!fLowGainMode) {
   for (int i = 0; (i < fNumModules) && (!fSingleModule); i++) {
     fPeaksCanvas[i]->cd();
     //Draw the high gain peak-minus pedestals
     fPedestals->GetPeakProfileHighGain((i < fNumModules) ? i : fVisibleModule)->Draw("colz");
     fPeaksCanvas[i]->Update();
     fPeaksCanvas[i]->Paint();
     
     fPedsCanvas[i]->cd();
     //Draw the high gain peak-minus pedestals
     fPedestals->GetPedProfileHighGain((i < fNumModules) ? i : fVisibleModule)->Draw("colz");
     fPedsCanvas[i]->Update();
     fPedsCanvas[i]->Paint();
   
     fPedsDiffCanvas[i]->cd();
     //Draw the high gain peak-minus pedestals
     fPedestals->GetPedProfileHighGainDiff((i < fNumModules) ? i : fVisibleModule)->Draw("colz");
     fPedsDiffCanvas[i]->Update();
     fPedsDiffCanvas[i]->Paint();
     
     fPeaksDiffCanvas[i]->cd();
     //Draw the high gain peak-minus pedestals
     fPedestals->GetPeakProfileHighGainDiff((i < fNumModules) ? i : fVisibleModule)->Draw("colz");
     fPeaksDiffCanvas[i]->Update();
     fPeaksDiffCanvas[i]->Paint();
     
     fPeaksRatioCanvas[i]->cd();
     //Draw the high gain peak-minus pedestals
     fPedestals->GetPeakProfileHighGainRatio((i < fNumModules) ? i : fVisibleModule)->Draw("colz");
     fPeaksRatioCanvas[i]->Update();
     fPeaksRatioCanvas[i]->Paint();
   
     fPedsRatioCanvas[i]->cd();
     //Draw the high gain peak-minus pedestals
     fPedestals->GetPedProfileHighGainRatio((i < fNumModules) ? i : fVisibleModule)->Draw("colz");
     fPedsRatioCanvas[i]->Update();
     fPedsRatioCanvas[i]->Paint();
   } // i = modules

   if (fSingleModule) {
     fPeaksCanvas[fNumModules]->cd();
     //Draw the high gain peak-minus pedestals
     fPedestals->GetPeakProfileHighGain(fVisibleModule)->Draw("colz");
     fPeaksCanvas[fNumModules]->Update();
     fPeaksCanvas[fNumModules]->Paint();
     
     fPedsCanvas[fNumModules]->cd();
     //Draw the high gain peak-minus pedestals
     fPedestals->GetPedProfileHighGain(fVisibleModule)->Draw("colz");
     fPedsCanvas[fNumModules]->Update();
     fPedsCanvas[fNumModules]->Paint();
     
     fPedsDiffCanvas[fNumModules]->cd();
     //Draw the high gain peak-minus pedestals
     fPedestals->GetPedProfileHighGainDiff(fVisibleModule)->Draw("colz");
     fPedsDiffCanvas[fNumModules]->Update();
     fPedsDiffCanvas[fNumModules]->Paint();
     
     fPeaksDiffCanvas[fNumModules]->cd();
     //Draw the high gain peak-minus pedestals
     fPedestals->GetPeakProfileHighGainDiff(fVisibleModule)->Draw("colz");
     fPeaksDiffCanvas[fNumModules]->Update();
     fPeaksDiffCanvas[fNumModules]->Paint();
     
     fPeaksRatioCanvas[fNumModules]->cd();
     //Draw the high gain peak-minus pedestals
     fPedestals->GetPeakProfileHighGainRatio(fVisibleModule)->Draw("colz");
     fPeaksRatioCanvas[fNumModules]->Update();
     fPeaksRatioCanvas[fNumModules]->Paint();
     
     fPedsRatioCanvas[fNumModules]->cd();
     //Draw the high gain peak-minus pedestals
     fPedestals->GetPedProfileHighGainRatio(fVisibleModule)->Draw("colz");
     fPedsRatioCanvas[fNumModules]->Update();
     fPedsRatioCanvas[fNumModules]->Paint();        
   }    
 } // not low gain
 else {
   for (int i = 0; (i < fNumModules) && (!fSingleModule); i++) {
     fPeaksCanvas[i]->cd();
     //Draw the high gain peak-minus pedestals
     fPedestals->GetPeakProfileLowGain((i < fNumModules) ? i : fVisibleModule)->Draw("colz");
     fPeaksCanvas[i]->Update();
     fPeaksCanvas[i]->Paint();
     
     fPedsCanvas[i]->cd();
     //Draw the high gain peak-minus pedestals
     fPedestals->GetPedProfileLowGain((i < fNumModules) ? i : fVisibleModule)->Draw("colz");
     fPedsCanvas[i]->Update();
     fPedsCanvas[i]->Paint();
     
     fPedsDiffCanvas[i]->cd();
     //Draw the high gain peak-minus pedestals
     fPedestals->GetPedProfileLowGainDiff((i < fNumModules) ? i : fVisibleModule)->Draw("colz");
     fPedsDiffCanvas[i]->Update();
     fPedsDiffCanvas[i]->Paint();
     
     fPeaksDiffCanvas[i]->cd();
     //Draw the high gain peak-minus pedestals
     fPedestals->GetPeakProfileLowGainDiff((i < fNumModules) ? i : fVisibleModule)->Draw("colz");
     fPeaksDiffCanvas[i]->Update();
     fPeaksDiffCanvas[i]->Paint();
     
     fPeaksRatioCanvas[i]->cd();
     //Draw the high gain peak-minus pedestals
     fPedestals->GetPeakProfileLowGainRatio((i < fNumModules) ? i : fVisibleModule)->Draw("colz");
     fPeaksRatioCanvas[i]->Update();
     fPeaksRatioCanvas[i]->Paint();
     
     fPedsRatioCanvas[i]->cd();
     //Draw the high gain peak-minus pedestals
     fPedestals->GetPedProfileLowGainRatio((i < fNumModules) ? i : fVisibleModule)->Draw("colz");
     fPedsRatioCanvas[i]->Update();
     fPedsRatioCanvas[i]->Paint();
   }

   if (fSingleModule) {
     fPeaksCanvas[fNumModules]->cd();
     //Draw the high gain peak-minus pedestals
     fPedestals->GetPeakProfileLowGain(fVisibleModule)->Draw("colz");
     fPeaksCanvas[fNumModules]->Update();
     fPeaksCanvas[fNumModules]->Paint();
     
     fPedsCanvas[fNumModules]->cd();
     //Draw the high gain peak-minus pedestals
     fPedestals->GetPedProfileLowGain(fVisibleModule)->Draw("colz");
     fPedsCanvas[fNumModules]->Update();
     fPedsCanvas[fNumModules]->Paint();
     
     fPedsDiffCanvas[fNumModules]->cd();
     //Draw the high gain peak-minus pedestals
     fPedestals->GetPedProfileLowGainDiff(fVisibleModule)->Draw("colz");
     fPedsDiffCanvas[fNumModules]->Update();
     fPedsDiffCanvas[fNumModules]->Paint();
     
     fPeaksDiffCanvas[fNumModules]->cd();
     //Draw the high gain peak-minus pedestals
     fPedestals->GetPeakProfileLowGainDiff(fVisibleModule)->Draw("colz");
     fPeaksDiffCanvas[fNumModules]->Update();
     fPeaksDiffCanvas[fNumModules]->Paint();
     
     fPeaksRatioCanvas[fNumModules]->cd();
     //Draw the high gain peak-minus pedestals
     fPedestals->GetPeakProfileLowGainRatio(fVisibleModule)->Draw("colz");
     fPeaksRatioCanvas[fNumModules]->Update();
     fPeaksRatioCanvas[fNumModules]->Paint();
     
     fPedsRatioCanvas[fNumModules]->cd();
     //Draw the high gain peak-minus pedestals
     fPedestals->GetPedProfileLowGainRatio(fVisibleModule)->Draw("colz");
     fPedsRatioCanvas[fNumModules]->Update();
     fPedsRatioCanvas[fNumModules]->Paint();
   }
 } // low gain 
 
 return;
}

//_____________________________________________________________________________
void TMCal::OnGainButton() 
{ // Gain button handler

  if (fLowGainMode) {
    fLowGainMode = kFALSE;
    fGainButton->SetText("Low gain");
  }
  else {
    fLowGainMode = kTRUE;
    fGainButton->SetText("High gain");
  }
 
  DrawHistos();
  
  return;
}

//_____________________________________________________________________________
void TMCal::OnRange(kTab tab)
{//Range change handler
  int i = 0;
  switch(tab) {

  case kPeds:
    fPedsRangeSingle->SetPosition(fPedsRange->GetPosition());
    fMaxPed = AliCaloCalibPedestal::GetSampleMax() - fPedsRange->GetPosition();
    //printf("Value changed to %u.\n", (unsigned int)fMaxPedDiff);
    for (i = 0; i < fNumModules; i++) {
      fPedestals->GetPedProfileHighGain(i)->GetZaxis()->SetRangeUser(-1.0*fMaxPed, fMaxPed);
      fPedestals->GetPedProfileLowGain(i)->GetZaxis()->SetRangeUser(-1.0*fMaxPed, fMaxPed);
      if (!fSingleModule) {//Don't update all the monitors if we're in single module mode, to speed up the redraw there
	fPedsCanvas[i]->Paint(); //For some reason we need paint-update-paint, just update-paint or paint-update won't do...
	fPedsCanvas[i]->Update();
	fPedsCanvas[i]->Paint();
      }
    }
    fPedsCanvasSingle->Paint(); //For some reason we need paint-update-paint, just update-paint or paint-update won't do...
    fPedsCanvasSingle->Update();
    fPedsCanvasSingle->Paint();
    break;
   
  case kPeaks:
    fPeaksRangeSingle->SetPosition(fPeaksRange->GetPosition());
    fMaxPeak = AliCaloCalibPedestal::GetSampleMax() - fPeaksRange->GetPosition();
    //printf("Value changed to %u.\n", (unsigned int)fMaxPedDiff);
    for (i = 0; i < fNumModules; i++) {
      fPedestals->GetPeakProfileHighGain(i)->GetZaxis()->SetRangeUser(-1.0*fMaxPeak, fMaxPeak);
      fPedestals->GetPeakProfileLowGain(i)->GetZaxis()->SetRangeUser(-1.0*fMaxPeak, fMaxPeak);
      if (!fSingleModule) {//Don't update all the monitors if we're in single module mode, to speed up the redraw there
	fPeaksCanvas[i]->Paint(); //For some reason we need paint-update-paint, just update-paint or paint-update won't do...
	fPeaksCanvas[i]->Update();
	fPeaksCanvas[i]->Paint();
      }
    }
    fPeaksCanvasSingle->Paint(); //For some reason we need paint-update-paint, just update-paint or paint-update won't do...
    fPeaksCanvasSingle->Update();
    fPeaksCanvasSingle->Paint();
    break;
  
  case kPedsDiff:
    fPedsDiffRangeSingle->SetPosition(fPedsDiffRange->GetPosition());
    fMaxPedDiff = AliCaloCalibPedestal::GetSampleMax() - fPedsDiffRange->GetPosition();
    //printf("Value changed to %u.\n", (unsigned int)fMaxPedDiff);
    for (i = 0; i < fNumModules; i++) {
      fPedestals->GetPedProfileHighGainDiff(i)->GetZaxis()->SetRangeUser(-1.0*fMaxPedDiff, fMaxPedDiff);
      fPedestals->GetPedProfileLowGainDiff(i)->GetZaxis()->SetRangeUser(-1.0*fMaxPedDiff, fMaxPedDiff);
      if (!fSingleModule) {//Don't update all the monitors if we're in single module mode, to speed up the redraw there
	fPedsDiffCanvas[i]->Paint(); //For some reason we need paint-update-paint, just update-paint or paint-update won't do...
	fPedsDiffCanvas[i]->Update();
	fPedsDiffCanvas[i]->Paint();
      }
    }
    fPedsDiffCanvasSingle->Paint(); //For some reason we need paint-update-paint, just update-paint or paint-update won't do...
    fPedsDiffCanvasSingle->Update();
    fPedsDiffCanvasSingle->Paint();
    break;
   
  case kPeaksDiff:
    fPeaksDiffRangeSingle->SetPosition(fPeaksDiffRange->GetPosition());
    fMaxPeakDiff = AliCaloCalibPedestal::GetSampleMax() - fPeaksDiffRange->GetPosition();
    //printf("Value changed to %u.\n", (unsigned int)fMaxPedDiff);
    for (i = 0; i < fNumModules; i++) {
      fPedestals->GetPeakProfileHighGainDiff(i)->GetZaxis()->SetRangeUser(-1.0*fMaxPeakDiff, fMaxPeakDiff);
      fPedestals->GetPeakProfileLowGainDiff(i)->GetZaxis()->SetRangeUser(-1.0*fMaxPeakDiff, fMaxPeakDiff);
      if (!fSingleModule) {//Don't update all the monitors if we're in single module mode, to speed up the redraw there
	fPeaksDiffCanvas[i]->Paint();
	fPeaksDiffCanvas[i]->Update();
	fPeaksDiffCanvas[i]->Paint();
      } 
    }
    fPeaksDiffCanvasSingle->Paint();
    fPeaksDiffCanvasSingle->Update();
    fPeaksDiffCanvasSingle->Paint();
    break;

  case kPedsRatio:
    fPedsRatioRangeSingle->SetPosition(fPedsRatioRange->GetPosition());
    //fMaxPedRatio = fgkMaxMaxRatio - fPedsDiffRange->GetPosition() + fgkMinMaxRatio;
    fMaxPedRatio = pow(2.0, -fPedsRatioRange->GetPosition()/fgkRatioStepScale);//Make the control logarithmic
    //printf("Pedsratio changed to %g, position is set at %i.\n", fMaxPedRatio, (int)fPedsRatioRange->GetPosition());
    for (i = 0; i < fNumModules; i++) {
      fPedestals->GetPedProfileHighGainRatio(i)->GetZaxis()->SetRangeUser(0.0, fMaxPedRatio);
      fPedestals->GetPedProfileLowGainRatio(i)->GetZaxis()->SetRangeUser(0.0, fMaxPedRatio);
      if (!fSingleModule) {//Don't update all the monitors if we're in single module mode, to speed up the redraw there
	fPedsRatioCanvas[i]->Paint();
	fPedsRatioCanvas[i]->Update();
	fPedsRatioCanvas[i]->Paint();
      } 
    }
    fPedsRatioCanvasSingle->Paint();
    fPedsRatioCanvasSingle->Update();
    fPedsRatioCanvasSingle->Paint();
    break;

  case kPeaksRatio:
    fPeaksRatioRangeSingle->SetPosition(fPeaksRatioRange->GetPosition());
    fMaxPeakRatio = pow(2.0, -fPeaksRatioRange->GetPosition()/fgkRatioStepScale);//Make the control logarithmic
    //printf("Pedsratio changed to %g, position is set at %i.\n", fMaxPeakRatio, (int)fPeaksRatioRange->GetPosition());
    for (i = 0; i < fNumModules; i++) {
      fPedestals->GetPeakProfileHighGainRatio(i)->GetZaxis()->SetRangeUser(0.0, fMaxPeakRatio);
      fPedestals->GetPeakProfileLowGainRatio(i)->GetZaxis()->SetRangeUser(0.0, fMaxPeakRatio);
      if (!fSingleModule) {//Don't update all the monitors if we're in single module mode, to speed up the redraw there
	fPeaksRatioCanvas[i]->Paint();
	fPeaksRatioCanvas[i]->Update();
	fPeaksRatioCanvas[i]->Paint();
      }
    }
    fPeaksRatioCanvasSingle->Paint();
    fPeaksRatioCanvasSingle->Update();
    fPeaksRatioCanvasSingle->Paint();
    break;
    
    //The single-histo-at-a-time sliders. We'll just set the other slider to match and call the normal handler :)
  case kPedsSingle:
    fPedsRange->SetPosition(fPedsRangeSingle->GetPosition());
    OnRange(kPeds);
    break;
  case kPeaksSingle:
    fPeaksRange->SetPosition(fPeaksRangeSingle->GetPosition());
    OnRange(kPeaks);
    break;
  case kPedsDiffSingle:
    fPedsDiffRange->SetPosition(fPedsDiffRangeSingle->GetPosition());
    OnRange(kPedsDiff);
    break;
  case kPeaksDiffSingle:
    fPeaksDiffRange->SetPosition(fPeaksDiffRangeSingle->GetPosition());
    OnRange(kPeaksDiff);
    break;
  case kPedsRatioSingle:
    fPedsRatioRange->SetPosition(fPedsRatioRangeSingle->GetPosition());
    OnRange(kPedsRatio);
    break;
  case kPeaksRatioSingle:
    fPeaksRatioRange->SetPosition(fPeaksRatioRangeSingle->GetPosition());
    OnRange(kPeaksRatio);
    break;
  default:
    printf("Unknown tab slider %u.\n", (int)tab);
    break;
  }

  return;
}

//_____________________________________________________________________________
void TMCal::OnModeButton(int button) 
{ // Module selection radio button handler
  if (button < fNumModules) {
    fVisibleModule = button;
    fSingleModule = true;
    
    fPeaksTab->HideFrame(fPeaksFrame);
    fPeaksTab->ShowFrame(fPeaksFrameSingle);
    
    fPedsTab->HideFrame(fPedsFrame);
    fPedsTab->ShowFrame(fPedsFrameSingle);
    
    fPeaksDiffTab->HideFrame(fPeaksDiffFrame);
    fPeaksDiffTab->ShowFrame(fPeaksDiffFrameSingle);
    
    fPedsDiffTab->HideFrame(fPedsDiffFrame);
    fPedsDiffTab->ShowFrame(fPedsDiffFrameSingle);
    
    fPeaksRatioTab->HideFrame(fPeaksRatioFrame);
    fPeaksRatioTab->ShowFrame(fPeaksRatioFrameSingle);
    
    fPedsRatioTab->HideFrame(fPedsRatioFrame);
    fPedsRatioTab->ShowFrame(fPedsRatioFrameSingle);
    
    fDeadMapTab->HideFrame(fDeadMapFrame);
    fDeadMapTab->ShowFrame(fDeadMapFrameSingle);
    
    DrawHistos();
  }
  else { // all
    fSingleModule = false;
    fPeaksTab->HideFrame(fPeaksFrameSingle);
    fPeaksTab->ShowFrame(fPeaksFrame);
    
    fPedsTab->HideFrame(fPedsFrameSingle);
    fPedsTab->ShowFrame(fPedsFrame);
    
    fPeaksDiffTab->HideFrame(fPeaksDiffFrameSingle);
    fPeaksDiffTab->ShowFrame(fPeaksDiffFrame);
    
    fPedsDiffTab->HideFrame(fPedsDiffFrameSingle);
    fPedsDiffTab->ShowFrame(fPedsDiffFrame);
    
    fPeaksRatioTab->HideFrame(fPeaksRatioFrameSingle);
    fPeaksRatioTab->ShowFrame(fPeaksRatioFrame);
    
    fPedsRatioTab->HideFrame(fPedsRatioFrameSingle);
    fPedsRatioTab->ShowFrame(fPedsRatioFrame);
    
    fDeadMapTab->HideFrame(fDeadMapFrameSingle);
    fDeadMapTab->ShowFrame(fDeadMapFrame);
    
    DrawHistos();//this is because we may not have drawn the histos to their canvases.
  }

  return;  
}

//_____________________________________________________________________________
void TMCal::OnRunNoChange(Long_t val) 
{
  //printf("Run number changed to %u.\n", (int)fRunNoEntry->GetIntNumber());
  InitRawReaderAndStream(fRunNoEntry->GetIntNumber());
  UpdateMonitors();
  //We still need to re-paint the tabs to make the change visible
  fPeaksFrame->Paint();
  return;
}

//_____________________________________________________________________________
void TMCal::OnTabChange(Int_t Id)
{
  //This sets the palette for the deadmap tab to the four colors we want there, and otherwise to
  //the standard no. 1 palette. To make this work, we need to do updatemonitors, which is slow,
  //so in case you're willing to sacrifice the custom palette for the deadmap in favor of speed,
  //just comment out/remove this function and the GetTab()->Connect(...) in the constructor.
  if (!strcmp(GetTab()->GetTabTab(Id)->GetString(), fgkDeadMapTabName)) {
    gStyle->SetPalette(AliCaloCalibPedestal::kNumDeadMapStates, fDeadMapPalette);
  }
  else {
    gStyle->SetPalette(1);
  }

  UpdateMonitors();
   
  return; 
}
