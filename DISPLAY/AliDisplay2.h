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

#include <TObject.h>
#include <RQ_OBJECT.h>
#include <TGStatusBar.h>

#include "AliDisplayFrame.h"
#include "AliSliderFrame.h"

class TObjArray;
class TGWindow;
class TEnv;

class AliModuleInfo;
class AliShutterFrame;
class AliInfoFrame;
class AliDetectorFrame;
class AliMenu;

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

//Constants to set the display mode
const Int_t kHits = BIT(0);
const Int_t kClusters = BIT(1);
const Int_t kHLT = BIT(2);
const Int_t kTracks = BIT(3);


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
 Int_t				GetCurrentView(void) const {return fCurrentView;};//return current view
 TGMainFrame*		        GetMainFrame() const {return fMainFrame;};//
 Float_t		    	GetZoomStep() const {return fZoomStep;};//
 Float_t	       		GetZoomFactor() const {return fZoomFactor;};//
 Float_t	      		GetRapidityMin() const {return fSliderFrame->GetRapidityMin();};//
 Float_t       			GetRapidityMax() const {return fSliderFrame->GetRapidityMax();};//
 Float_t	      		GetMomentumMin() const {return fSliderFrame->GetMomentumMin();};//
 Float_t       			GetMomentumMax() const {return fSliderFrame->GetMomentumMax();};//
 Int_t				GetNbParticles() const {return fNbParticles;};//
 Int_t		       		GetEventNumber() const {return fEventNumber;};//
 Int_t			       	GetNbHits() const {return fNbHits;};//
 Int_t                          GetNbClusters() const {return fDisplayFrame->GetNbClusters();};//
 Float_t		       	GetSliderStep() const {return fSliderStep;};//
 Bool_t			       	GetZoomMode() const {return fZoomMode;};//
 Int_t			       	GetMode() const {return fMode;};//
 TObjArray*		       	GetModules() const {return fModules;};//
 Bool_t*	       		GetEnabledModules() const {return fEnabledModules;};//
 Int_t			       	GetNbModules() const {return fNbModules;};//
 char*                          GetIconsPath() const {return fIconsPath;};//
 AliModuleInfo*                    GetModuleInfo() const {return fModuleInfo;};//
 Bool_t                         GetSliderUpdate() const {return fSliderUpdate;};//
 char*                          GetRawDataPath() const {return fRawDataPath;};

 //Setter	
 void                           SetSliderUpdate(Bool_t b){fSliderUpdate = b;};//
 void			       	SetMode(Int_t m){fMode=m;};//
 void                           Enable(Int_t m);//Enable the given mode 
 Bool_t                         IsEnabled(Int_t m) const;//Return if the given mode is enabled
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
 void                           SaveToRC() const;

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
 TGLayoutHints	        	*fSliderFrameLayout; // Slider Frame Layout
 AliSliderFrame			*fSliderFrame;//Frame which contains the rapidity and momentum sliders

 //Shutter Frame
 TGLayoutHints		        *fShutterFrameLayout; // Shutter Frame Layout
 AliShutterFrame	        *fShutterFrame;// Shutter Frame

 //Display Frame
 TGLayoutHints	        	*fDisplayFrameLayout; // Display Frame Layout
 AliDisplayFrame	        	*fDisplayFrame;// Display Frame

 //Info Frame
 TGLayoutHints	        	*fInfoFrameLayout; // Info Frame Layout
 AliInfoFrame			*fInfoFrame; // Info Frame

 //Detector Option Frame
 TGLayoutHints	        	*fDetectorFrameLayout;// Detector Frame Layout
 AliDetectorFrame	        *fDetectorFrame; // Detector Frame

 //MenuBar
 AliMenu				*fMenu; // Menu

 //Status bar
 TGStatusBar			*fStatusBar; // Status Bar

 RQ_OBJECT("AliDisplay2")

 ClassDef(AliDisplay2,0);
};

R__EXTERN  AliDisplay2* gAliDisplay2;

#endif

