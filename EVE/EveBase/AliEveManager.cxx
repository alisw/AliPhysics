#include <TInterpreter.h>
#include <TGeoManager.h>
#include <TEveBrowser.h>
#include <TEveGedEditor.h>
#include <TEveManager.h>
#include <TEveViewer.h>
#include <TEveSelection.h>
#include <TSystem.h>
#include <TString.h>
#include <TROOT.h>

#include <AliLog.h>
#include <AliEveConfigManager.h>
#include <AliEveMultiView.h>

#include "AliEveManager.h"

//______________________________________________________________________________
// AliEveManager
//
// Central aplication manager for AliEve.
// Manages environment, gEve, elements, GUI, GL scenes and GL viewers.
//
// ALICE_ROOT must be defined prior creating this object
ClassImp(AliEveManager);

AliEveManager::AliEveManager(UInt_t w, UInt_t h, Bool_t map_window, Option_t* opt)
	: TEveManager(w, h, map_window, opt)
{

	Init();
	
}

AliEveManager* AliEveManager::Create(Bool_t map_window, Option_t* opt)
{
	Int_t w = 1024;
	Int_t h =  768;
  
  if(gEve == 0){
      gEve = new AliEveManager(w, h, map_window, opt);
  }
	
	return (AliEveManager*)gEve;
}

AliEveManager::~AliEveManager()
{ 
      AliEveMultiView* mv = AliEveMultiView::Instance();
      
      delete mv;
      mv = 0;
}

void AliEveManager::Init()
{
	GetDefaultViewer()->SetElementName("3D View");
  GetSelection()->SetPickToSelect(TEveSelection::kPS_PableCompound);
  GetHighlight()->SetPickToSelect(TEveSelection::kPS_PableCompound);

	RegisterGeometryAlias("Default", Form("%s/EVE/alice-data/default_geo.root",  gSystem->Getenv("ALICE_ROOT")) );

	AliEveConfigManager::InitializeMaster(); // initializes menus
}

void AliEveManager::CloseEveWindow()
{
		// Close button haas been clicked on EVE main window (browser).
   // Cleanup and terminate application.

   TEveBrowser *eb = dynamic_cast<TEveBrowser*>( GetMainWindow() );
 	 eb->DontCallClose();
 		
	TEveGedEditor::DestroyEditors();

}

void AliEveManager::Terminate()
{

}
