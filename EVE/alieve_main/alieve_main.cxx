// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007
// Mihai Niculescu 2012

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
#include <TEveBrowser.h>

#include <AliEveApplication.h>
#include <AliEveManager.h>

int main(int argc, char **argv)
{
  AliEveApplication app("AliEve", &argc, argv, 0 , -1);
  gEve = AliEveManager::Create();
  
  // close application when closing EveWindow from X window button
	gEve->GetBrowser()->Connect("CloseWindow()", "TApplication", &app, "Terminate(=0)");
	
  // process command line arguments
	app.GetOptions(&argc, argv);
	
	app.Run();

	if(gEve){
	 delete gEve;
   gEve = 0;
	}
	
  return 0;
}
