/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/////////////////////////////////////////////////////////////////////////
// ALICE EVENT DISPLAY CLASS                                           //
// Author: Mayeul   ROUSSELET                                          //
// e-mail: Mayeul.Rousselet@cern.ch                                    //
// Last update:26/08/2003                                              //
/////////////////////////////////////////////////////////////////////////

#ifndef AliDISPLAY2_H
#define AliDISPLAY2_H

#include <RQ_OBJECT.h>
#include <TGFrame.h>
#include <TGDoubleSlider.h>
#include <TRootEmbeddedCanvas.h>
#include <TGStatusBar.h>
#include <TGToolBar.h>
#include <TGDimension.h>
#include <TVirtualPad.h>

class TGButton;
class TGMenu;
class TGMenuBar;
class TGPopupMenu;
class TGPopupMenu;
class TGPopupMenu;
class TGTab;
class TGLayoutHints;
class TObjArray;
class TPolyMarker3D;
class TGNumberEntry;
class TGNumberEntryField;
class TGCheckButton;
class TGTextButton;
class TGShutter;
class TGShutterItem;
class TGWindow;
class TEnv;
class TGLabel;
class TCanvas;

enum AliDisplay2Id {
	kIdsVIEW,		//id of the view shutter
	kIdbSIDEVIEW,	//id of sideview button
	kIdbTOPVIEW,
	kIdbFRONTVIEW,
	kIdbALLVIEW,

	kIdsEVENT,
	kIdbPrevEVENT,
	kIdbNextEVENT,

	kIdbZoomIN,
	kIdbZoomOUT,
	kIdbZoomZONE,

	kIdiNBPARTICULES,
	kIdiNBEVENTS,

	kIdsRAPIDITY,
	kIdsMOMENTUM,
	kIdsOPTIONS,

	kIdsMomentumMIN,
	kIdsMomentumMAX,
	kIdsRapidityMIN,
	kIdsRapidityMAX,

	kIdbSelectALL,
	kIdbSelectINVERT,

	kIdsDETECTORS,
	kIdbCheckTPC,
	kIdbCheckMUON,
	kIdbCheckTOF,
	kIdbCheckTRD,
	kIdbCheckPHOS,
	kIdbCheckVZERO,
	kIdbCheckEMCAL,
	kIdbCheckSTART,
	kIdbCheckPMD,
	kIdbCheckFMD,
	kIdbCheckZDC,
	kIdbCheckRICH,
	kIdbCheckITS,
	
	kIdbCheckHITS,
	kIdbCheckCLUSTERS,
	kIdbCheckHLT,
	kIdbCheckTRACKS,

	kIdmOPEN,
	kIdmSAVEAS,
	kIdmCLOSE,
	kIdmPRINT,
	kIdmPRINTSETUP,
	kIdmEXIT,
	kIdmSETTINGS,
	kIdmSAVESETTINGS,
	kIdmHELP,
	kIdmABOUT,
	kIdmVIEWX3D,
	kIdmVIEWGL,

	kIdtZoomSTEP,
	kIdtSliderSTEP,
	kIdtSliderUPDATE,
	
	kmCUTS,		//tag for the update fucntion
	kmMODULES,
	kmPOINTS
};

const char helpTxt[] = "\tAliDisplay v2.0\n\t\tHelp\n\n\nWelcome in the AliDisplay help.\nHere is a list of useful subjects which discribes\nthe main functionnalities of the software\n \nEvent:Use the arrows to get the next or previous event\nView:Each button corresponds to a different view\nDetectors:Select the module you want to see\nOptions:Select the view mode\nSliders:Use the rapidity (or eta) slider to cut the set of hits\n\tAnd the momentum slider to cut with respect to the momentum\n";


//const char *filetypes[] = {"ROOT files","*.root","All files","*",0,0};


//Constants to set the display mode
const Int_t kHits = BIT(0);
const Int_t kClusters = BIT(1);
const Int_t kHLT = BIT(2);
const Int_t kTracks = BIT(3);

//***Class ModulesInfo***//
class AliModuleInfo{
 public:
  AliModuleInfo(int n);
  virtual ~AliModuleInfo();

  void     SetId(char* name,Int_t id);
  void     Add(const char * name,Int_t i);
  Int_t    Id(char *name);
  char*    Name(Int_t id);
  Bool_t   IsEnabled(Int_t id);
  Bool_t   IsEnabled(char *name){return IsEnabled(Id(name));};
  void     Disable(Int_t id);
  void     Disable(char* name){Disable(Id(name));};
  void     Enable(Int_t id);
  void     Enable(char *name){Enable(Id(name));};
  void     Print();
  
 private:
  //The purposes of this class is to link each module to its Id
  char     **fName;
  Int_t    *fId;
  Bool_t   *fEnabled;
  Int_t    fNb;

  ClassDef(AliModuleInfo,0);
};

//***Class DisplayCluster***//

class AliDisplayClusters{
  //This class is an interface to the clusters data

public:
 AliDisplayClusters();
 virtual ~AliDisplayClusters();

 void          LoadClusters(char * name,Int_t nevent);
 void          LoadITSClusters(Int_t nevent);
 void          LoadTPCClusters(Int_t nevent);
 void          Draw();
 Int_t         GetNbClusters();

private: 
 TPolyMarker3D *fPoints; //fPoints[i]=set of cluster coordinates in detector i;
 Int_t         fNb;
 char          **fName; //fName[i]=name of the detector i 

 RQ_OBJECT("AliDisplayClusters")

 ClassDef(AliDisplayClusters,0);
};

//***class AliDisplayHLT***//

class AliDisplayHLT{
  //This classes is an interface to the HLT data
  //For the moment only for TPC, for adding modules there is two choices:
  //1) add the function LoadHLT[module](Int_t) and update the function LoadHLT
  //2) or inherit your class from AliDisplayHLT and overload LoadHLT

 public:

  AliDisplayHLT();
  virtual ~AliDisplayHLT();

  virtual void  LoadHLT(char *name,Int_t e);//Load L3 datas whose belong to detector name and from the event e
  virtual void  LoadHLTTPC(Int_t nevent);
  virtual void  Draw();

 private:
  TPolyMarker3D *fPoints; //fPoints[i]=set of cluster coordinates in detector i;
  Int_t         fNb;
  char          **fName; //fName[i]=name of the detector i 

 ClassDef(AliDisplayHLT,0);
};

//***Class AliSliderFrame***//
class AliSliderFrame{
  //This class implements the cuts manager

public:

 AliSliderFrame(const TGWindow *p, UInt_t w, UInt_t h);
 virtual ~AliSliderFrame();

 //Setters
 void		      	SetMomentumRange(Float_t min, Float_t max){fMomentumSlider->SetRange(min,max);};
 void		       	SetRapidityRange(Float_t min, Float_t max){fRapiditySlider->SetRange(min,max);};
 
 Float_t       	       	GetMomentumMax(){return fMomentumSlider->GetMaxPosition();};
 Float_t       		GetMomentumMin(){return fMomentumSlider->GetMinPosition();};
 Float_t	      	GetRapidityMax(){return fRapiditySlider->GetMaxPosition();};
 Float_t       		GetRapidityMin(){return fRapiditySlider->GetMinPosition();};
 TGCompositeFrame*	GetSliderFrame(){return fMainFrame;};

 //Slots
 void		     	CloseWindow();
 void		       	DoSlider(Int_t pos=0);
 void		       	DoField(Long_t pos=0);
 void		       	DoReleased(Int_t pos=0);
 void                   DoPositionChanged(Int_t pos=0);

 //I/O
 void		       	SaveToRC();
 void		       	LoadFromRC();

private:
 
 TGCompositeFrame  	*fMainFrame;
 TGCompositeFrame	*fMomentumFrame;
 TGCompositeFrame	*fRapidityFrame;
 TGLayoutHints		*fLayout;//Layout of the frame
 TGLayoutHints		*fMomentumLayout;
 TGLayoutHints		*fRapidityLayout;
 TGDoubleHSlider       	*fMomentumSlider;
 TGDoubleHSlider       	*fRapiditySlider;
 TGLabel       		*fMomentumLabel;
 TGLabel	       	*fRapidityLabel;
 TGNumberEntry		*fMomentumMaxValue;
 TGNumberEntry		*fMomentumMinValue;
 TGNumberEntry		*fRapidityMaxValue;
 TGNumberEntry		*fRapidityMinValue;

 RQ_OBJECT("AliSliderFrame")

 ClassDef(AliSliderFrame,0);
};

//***Class AliDetectorFrame***//
class AliDetectorFrame{
  //This classe implements the frame which contains the list of enabled detectors
  //and allows the user to change the current status of the detector (enabled/disabled)

public:

 AliDetectorFrame(const TGWindow *p, Int_t w, Int_t h,UInt_t bgc);
 virtual ~AliDetectorFrame();
 
 TGCompositeFrame*	GetDetectorFrame(){return fMainFrame;};
 
 //Slots
 void		      	DoButton(Int_t pos=0);
 void		      	DoCheckButton(Int_t pos=0);
 void		       	DoSpecific();	

private:

 TGCompositeFrame	*fMainFrame;
 TGCompositeFrame	*fButtonFrame;
 TGCheckButton		**fCheckButton;
 TGTextButton		*fButtonAll;
 TGTextButton		*fButtonInvert;
 Bool_t		       	*fCheckedButton; //fEnabledButton[i]=kTRUE if button checked
 Int_t		       	*fCheckButtonId;
 Bool_t		       	fCheckedMode;
 static int	       	fgBaseId;

 RQ_OBJECT("AliDetectorFrame")

 ClassDef(AliDetectorFrame,0);
};

//***Class AliShutterItem**//
class AliShutterItem{
  //This class implements the shutter item, ie the base element of a shutter and provides functions to add button... in the shutter
public:
	
 AliShutterItem(TGShutter *s,char *text,UInt_t id);
 virtual ~AliShutterItem();

 //Getters
 TGShutterItem*		GetShutterItem(){return fShutterItem;};
 TGCompositeFrame*	GetShutterItemFrame(){return fMainFrame;};

 //Fill functions
 void	      		AddTextButton(char *text,char *tiptext,  UInt_t idb);
 void	       		AddPictureButton(char *file,char *tiptext,UInt_t idb);
 void	       		AddCheckButton(char *txt,Int_t idb);

 //Slot
 void	       		DoButton(Int_t pos=0);

private:

 TGCompositeFrame	*fMainFrame;
 TGShutterItem		*fShutterItem;
 TGButton		*fButton;

 RQ_OBJECT("AliShutterItem")

 ClassDef(AliShutterItem,0);
};


//***Class AliShutterFrame***//
class AliShutterFrame{
  //This class implements the shutter frame
public:

 AliShutterFrame(TGCompositeFrame *p, UInt_t w, UInt_t h);
 virtual ~AliShutterFrame();

 TGCompositeFrame*	GetShutterFrame(){return fMainFrame;};

private:

 TGCompositeFrame	*fMainFrame;
 TGLayoutHints		*fLayout;
 TGShutter	       	*fShutter;
 AliDetectorFrame	*fDetectorFrame;
 TGLayoutHints		*fDetectorFrameLayout;

 RQ_OBJECT("AliShutterFrame")

 ClassDef(AliShutterFrame,0);
};

//***Class AliDisplayFrame***//
class AliDisplayFrame{
  //This class implements the display of the event

public:

 AliDisplayFrame(const TGWindow *p, UInt_t w, UInt_t h);
 virtual ~AliDisplayFrame();
 
 //Getters
 TGCompositeFrame*		GetDisplayFrame(){return fMainFrame;};
 TCanvas*	       		GetMainCanvas(){return fMainCanvas;};
 Int_t		       		GetPreviousW(){return fPreviousW;};
 Int_t		       		GetPreviousH(){return fPreviousH;};
 TGDimension	       		GetFrameDimension(){return ((TGCanvas*)fMainEmbeddedCanvas)->GetViewPort()->GetDefaultSize();};
 Int_t		       		GetNbActivePoints();
 Int_t                          GetNbClusters(){return fClusters->GetNbClusters();};

 //Setters
 void		       		SetPreviousW(Int_t w){fPreviousW=w;};
 void	       			SetPreviousH(Int_t h){fPreviousH=h;};
 void	       			SetEditable(Bool_t b){gPad->SetEditable(b);};
 
 void	       			DoView(Int_t view);
 void	       			Draw(Float_t theta,Float_t phi,Float_t psi);
 void	       			DrawDetector(const char *name);
 void	       			DrawHits();
 void	       			DrawX3d();
 void	       			DrawGL();
 void	       			LoadEnabledModules();
 void                           LoadClusters(Int_t nevent);
 void                           LoadHLTClusters(Int_t nevent);
 void		       		LoadHits();
 void		       		ApplyCuts();
 void		       		EnableDetector(const char *name);
 void		       		DisableDetector(const char *name);
 void		       		ExecuteEvent(Int_t event, Int_t px,Int_t py,TObject *);
 void                           SavePadGIF(const char *file);

private:

 TGCompositeFrame		*fMainFrame,*fFrame1,*fFrame2;
 TGTab			       	*fMainTab;
 Bool_t			       	fAllViews;
 TRootEmbeddedCanvas		*fMainEmbeddedCanvas;//embedded Canvas which contains the main view(s)
 TRootEmbeddedCanvas		*fSelectionEmbeddedCanvas;
 TCanvas	       		*fMainCanvas;
 TCanvas       			*fSelectionCanvas;
 Float_t     			fClipMin,fClipMax;
 Int_t				fPreviousW,fPreviousH;
 Float_t	      		fRange;
 AliDisplayClusters             *fClusters;
 AliDisplayHLT                  *fHLT;
 TObjArray	      		*fPoints;
 TObjArray		       	*fPoints2;
 TObjArray		        *fModules;
 Int_t			       	fNbModules;
 Bool_t		       		*fActivePoints;
 TObjArray	       		*fPolyMarkers;//Array for TPolyMarker3D
 Float_t	       		*fClustersPos;
 Int_t		       		fNbClusters;
	
 RQ_OBJECT("AliDisplayFrame")

 ClassDef(AliDisplayFrame,0);
};

//***Class AliInfoFrame***//
class AliInfoFrame{
  //This class implements the info frame where the number of particles... are displayed

public:
	
 AliInfoFrame(TGCompositeFrame *p, UInt_t w, UInt_t h);
 virtual ~AliInfoFrame(void);
 
 void			AddLabel(char *text, UInt_t options);
 TGCompositeFrame	*GetInfoFrame(){return fMainFrame;};
 void 			Update();
 
private:

 TGCompositeFrame	*fMainFrame,*fTitleFrame,*fFiguresFrame;
 TGLabel       		*fNbParticuleLabel;
 TGLabel       	      	*fNbEventLabel;
 TGLabel       	       	*fNbHitsLabel;
 TGLabel                *fNbClustersLabel;        

 RQ_OBJECT("AliInfoFrame")

 ClassDef(AliInfoFrame,0);
};


//***Class AliSettingFrame***//
class AliSettingFrame:public TGTransientFrame{
  //This classe implement the setting frame where the different otption can be set

public:

 AliSettingFrame(const TGWindow *p, const TGWindow *main, UInt_t w, UInt_t h);
 virtual ~AliSettingFrame();

 //Slots
 void					DoSettings(Int_t id=0);

private:

 TGCompositeFrame		*fMainFrame;
 TGCompositeFrame		*fZoomStepFrame;
 TGLayoutHints			*fZoomStepLayout;
 TGNumberEntryField		*fZoomStepEntry;
 TGLabel       			*fZoomStepLabel;     
 TGCompositeFrame		*fSliderStepFrame;
 TGLayoutHints			*fSliderStepLayout;
 TGNumberEntryField		*fSliderStepEntry;
 TGLabel       			*fSliderStepLabel;
 TGCompositeFrame		*fSliderUpdateFrame;
 TGLayoutHints			*fSliderUpdateLayout;
 TGCheckButton                  *fSliderUpdateButton; 
 Bool_t                         fIsLoading;//Used when retrieving the state of the check button

 RQ_OBJECT("AliSettingFrame")

 ClassDef(AliSettingFrame,0);
};


//***Class AliMenu***//
class AliMenu {
  //This class implement both the menu and the toolbar

public:

 AliMenu(TGCompositeFrame *p, UInt_t w, UInt_t h, UInt_t options);
 virtual ~AliMenu();
 
 void 		      		AddPictureButton(char *fname,char *tiptext,UInt_t id, UInt_t spacing);//add a picture button to the toolbar
 //slots
 void		      		DoMenu(Int_t id=0);
 void		       		DoToolBar(Int_t id=0);

private:

 TGMenuBar			*fMenuBar;
 TGToolBar		       	*fToolBar;
 TGLayoutHints			*fMenuBarLayout;
 TGLayoutHints			*fMenuBarItemLayout;
 TGPopupMenu		       	*fMenuFile;
 TGPopupMenu		       	*fMenuOptions;
 TGPopupMenu		       	*fMenuHelp;
 TGPopupMenu		       	*fMenuView;
 TGLayoutHints			*fToolBarLayout;
 ToolBarData_t 			*fTBD;
 TGButton                       *fButton;

 RQ_OBJECT("AliMenu")

 ClassDef(AliMenu,0);
};

//***Class AliDisplay2***//
class AliDisplay2 : public TObject{
  //This classe is the main component of the new display. 
  //It aggregates all the subframes and manage the relationnships 
  //between the widgets and the internal representation

public:

 AliDisplay2(const TGWindow *p, UInt_t w, UInt_t h); //Constructor, load all the widgets
 virtual ~AliDisplay2(void);

 //Slots
 void				CloseWindow();//
 void				DoView(Int_t pos){fDisplayFrame->DoView(pos);};//change the current view
		
 void				DoSaveSettings();//Save settings to the ressource file
 void				LoadSettings();//load the settings from the ressource file

 //Getter
 Int_t				GetCurrentView(void){return fCurrentView;};//return current view
 TGMainFrame*		        GetMainFrame(){return fMainFrame;};//
 Float_t		    	GetZoomStep(){return fZoomStep;};//
 Float_t	       		GetZoomFactor(){return fZoomFactor;};//
 Float_t	      		GetRapidityMin(){return fSliderFrame->GetRapidityMin();};//
 Float_t       			GetRapidityMax(){return fSliderFrame->GetRapidityMax();};//
 Float_t	      		GetMomentumMin(){return fSliderFrame->GetMomentumMin();};//
 Float_t       			GetMomentumMax(){return fSliderFrame->GetMomentumMax();};//
 Int_t				GetNbParticles(){return fNbParticles;};//
 Int_t		       		GetEventNumber(){return fEventNumber;};//
 Int_t			       	GetNbHits(){return fNbHits;};//
 Int_t                          GetNbClusters(){return fDisplayFrame->GetNbClusters();};//
 Float_t		       	GetSliderStep(){return fSliderStep;};//
 Bool_t			       	GetZoomMode(){return fZoomMode;};//
 Int_t			       	GetMode(){return fMode;};//
 TObjArray*		       	GetModules(){return fModules;};//
 Bool_t*	       		GetEnabledModules(){return fEnabledModules;};//
 Int_t			       	GetNbModules(){return fNbModules;};//
 char*                          GetIconsPath(){return fIconsPath;};//
 AliModuleInfo*                    GetModuleInfo(){return fModuleInfo;};//
 Bool_t                         GetSliderUpdate(){return fSliderUpdate;};//
 char*                          GetRawDataPath(){return fRawDataPath;};

 //Setter	
 void                           SetSliderUpdate(Bool_t b){fSliderUpdate = b;};//
 void			       	SetMode(Int_t m){fMode=m;};//
 void                           Enable(Int_t m);//Enable the given mode 
 Bool_t                         IsEnabled(Int_t m);//Return if the given mode is enabled
 void                           Disable(Int_t m);//Disable the given mode
 void		      		SetZoomMode(Bool_t b=kTRUE){fZoomMode = b;};//
 void		       		SetCurrentView(Int_t id){fCurrentView = id;};//
 void		       		SetZoomStep(Float_t zs){fZoomStep=zs;};//
 void		       		SetZoomFactor(Float_t zf){fZoomFactor=zf;};//
 void		       		SetNbHits(Int_t hi){fNbHits = hi;};//
 void		       		SetNbParticles(Int_t nbp){fNbParticles = nbp;};//
 void		       		SetSliderStep(Float_t st){fSliderStep = st;};//
 void		       		SetEditable(Bool_t b){fDisplayFrame->SetEditable(b);};//
 void		       		SetStatusBar(const char *s,Int_t p){fStatusBar->SetText(s,p);};//Change text in the p part of the status bar
 void                           SetRawDataPath(char *path){fRawDataPath=new char[strlen(path)];strcpy(fRawDataPath,path);};
	
 void                           SavePadGIF(const char* file){fDisplayFrame->SavePadGIF(file);};//save current pad to gif file
 void	       			Draw(Option_t *options=0);//draw 
 void	       			DrawDetector(const char *name){fDisplayFrame->DrawDetector(name);};//draw detector specific view (not implemented)
 void	       			DrawX3d();//Draw current pad in x3d view
 void	       			DrawGL();//Draw current pad in OpenGL view
 void 	       			ShowNextEvent(Int_t delta);//Load and display the current+delta event, if it exists
 void	       			LoadClusters(Int_t nevent);//Load the clusters of the event 
 void	       			LoadHits();//Load the hits
 void                           LoadHLTClusters(Int_t nevent);//Load the hlt clusters
 void	       			LoadEnabledModules(){fDisplayFrame->LoadEnabledModules();};//Load enabled modules
 void	       			LoadEnabledHits(){fDisplayFrame->LoadHits();};//Load enabled hits
 void  				ApplyCuts(){fDisplayFrame->ApplyCuts();}//Apply cuts from the slider frame
 void  				EnableDetector(const char *name){fDisplayFrame->EnableDetector(name);};//Enable detector "name"
 void  				DisableDetector(const char *name){fDisplayFrame->DisableDetector(name);};//Disable detector "name"
 void  				Update(Int_t tag=-1);//Update the view, if loading only the modified data from the previous changes, the integer tag indicates the kind of modification
 void  				HandleMouseWheel(Event_t *event);//Handle mouve event, not working yet
 void  				HandleResize(Event_t *event);//Handle resize event
 void  				FindModules();//Find the modules used for the simulation

 //I/O
 void                           LoadFromRC();
 void                           SaveToRC();

private:

 Float_t 			fZoomStep;//Step of the zoom, ie the factor by which the zoom factor will be multiplied(divided) when pressing the zoom in (out) button
 Float_t		        fZoomFactor;//Zoom factor, ie 1=no zoom, >1 zoom in & <1 zoom out
 Bool_t				fZoomMode;//kTrue if zoom on zone is enabled
 Int_t				fCurrentView;//Current view, see the enum for value
 Int_t				fNbParticles;//Current number of displayed particles
 Int_t				fEventNumber;//Number of event
 Int_t				fNbHits;//Number of displayed hits
 Float_t			fSliderStep;//Step of the slider
 Int_t				fMode;//Display mode : each bit is associated to a mode : Hits=BIT(0), clusters=BIT(1)...
 Int_t				fNbModules;//Total number of modules
 TObjArray			*fModules;//Array of modules
 Bool_t				*fEnabledModules;//Array of bool for coding the enabled modules
 char                           *fIconsPath;//Icon Path (by default $ALICE_ROOT/DISPLAY/icons/)
 char                           *fRawDataPath;//Raw data path
 Bool_t                         fHitsLoaded,fClustersLoaded,fHLTLoaded,fTracksLoaded;//Implement if the type of datas was loaded, to avoid uneeded load
 Bool_t                         fSliderUpdate;//True if display frame update on slider move, for fast machine only (desactivated by default)
 AliModuleInfo                     *fModuleInfo;//Pointer to the class which map module name and internal module ID, necessary because of the dynamic load of the available modules
	
 TEnv				*fAliDisplay2rc;//ressources file (.alidisplay.rc)
	
 TGMainFrame			*fMainFrame;//Main frame
 TGCompositeFrame	        *fLeftFrame;//frame used for the layout
 TGCompositeFrame          	*fRightFrame;//frame used for the layout
 TGCompositeFrame       	*fSubFrame;//frame used for the layout	

 //Slider Frame
 TGLayoutHints	        	*fSliderFrameLayout;
 AliSliderFrame			*fSliderFrame;//Frame which contains the rapidity and momentum sliders

 //Shutter Frame
 TGLayoutHints		        *fShutterFrameLayout;
 AliShutterFrame	        *fShutterFrame;

 //Display Frame
 TGLayoutHints	        	*fDisplayFrameLayout;
 AliDisplayFrame	        	*fDisplayFrame;

 //Info Frame
 TGLayoutHints	        	*fInfoFrameLayout;
 AliInfoFrame			*fInfoFrame;

 //Detector Option Frame
 TGLayoutHints	        	*fDetectorFrameLayout;
 AliDetectorFrame	        *fDetectorFrame;

 //MenuBar
 AliMenu				*fMenu;

 //Status bar
 TGStatusBar			*fStatusBar;

 RQ_OBJECT("AliDisplay2")

 ClassDef(AliDisplay2,0);
};

#endif

