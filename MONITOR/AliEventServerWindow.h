// Author: Mihai Niculesu 2013

/**************************************************************************
 * Copyright(c) 1998-2013, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEventServerWindow_H
#define AliEventServerWindow_H

#include <TObject.h>
#include <TString.h>
#include <TGFrame.h>
#include <TGLabel.h>

class TGTextButton;
class TGCheckButton;
class TGListBox;

class AliRecoServer;
class AliDimIntNotifier;

//______________________________________________________________________________
// Short description of AliEventServerWindow
//

class AliEventServerWindow : public TGMainFrame
{
public:
enum TOOLBUTTON{
  	TOOLBUTTON_START=1,
  	TOOLBUTTON_STOP,
  	TOOLBUTTON_PREFERENCES,
  	TOOLBUTTON_EXIT 	
  };
  
  AliEventServerWindow();
  virtual ~AliEventServerWindow();

  //------------------------------------------------------------------------------
  // Handlers of DIM signals.
  //------------------------------------------------------------------------------

  void StartOfRun(Int_t run);
  void EndOfRun(Int_t run);
  
  //------------------------------------------------------------------------------
  // Handlers of button signals.
  //------------------------------------------------------------------------------
	void onStartServer();
	void onStopServer();
	void onExit();

	void HandleToolBarAction(Int_t id=-1);

   
private:
  AliEventServerWindow(const AliEventServerWindow&);            // Not implemented
  AliEventServerWindow& operator=(const AliEventServerWindow&); // Not implemented
	void InitDIMListeners();
	void FillRunsFromDatabase();
	void SetupToolbar();
	
	void LaunchRecoServer();	
	void StartReco(Int_t run);
	bool StopRecoServer();
	
  // GUI components.
  TGListBox     *fRunList;    // List-box for listing current runs.
  TGTextButton  *fStartServButt;  // Start server for selected run.
  TGTextButton  *fStopServButt;   // Close server for selected run.
  TGTextButton  *fExitButt;   // Close server and do Exit.
  
  // DIM interface. Could do without members and just leak them ...
  AliDimIntNotifier *fDimSORListener[5]; // DIM listeners for SOR.
  AliDimIntNotifier *fDimEORListener[5]; // DIM listeners for EOR.

  // server state & process management
  Int_t fRunRunning;   // Run which is executed.
	AliRecoServer* fRecoServer;


  ClassDef(AliEventServerWindow, 0);
};

#endif
