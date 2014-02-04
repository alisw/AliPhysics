// Author: Mihai Niculescu 2013

/**************************************************************************
 * Copyright(c) 1998-2013, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
 
#include <TEnv.h>
#include <TObjString.h>
#include <TString.h>
#include <TSystem.h>

#include <TGButton.h>
#include <TGLabel.h>
#include <TGMsgBox.h>
#include <TGTab.h>
#include <TGTextEntry.h>

#include <AliReconstruction.h>

#include "AliEventServerUtil.h"
#include "AliEventServerPreferencesWindow.h"

//_________________________________________________________
ClassImp(AliEventServerPreferencesWindow)
AliEventServerPreferencesWindow::AliEventServerPreferencesWindow(const TGWindow* p, const char* window_title)
	: TGTransientFrame(gClient->GetRoot(), p),
	fTab(0),
	fEntryServerHost(0),
	fEntryServerPort(0),
	fEntryCDBDefaultStorage(0),
	fEntryCDBSpecStoragePath1(0),
	fEntryCDBSpecStorageValue1(0),
	fEntryCDBSpecStoragePath2(0),
	fEntryCDBSpecStorageValue2(0),
	fEntryCDBSpecStoragePath3(0),
	fEntryCDBSpecStorageValue3(0),
	fEntryRecoRunQA(0),
	fEntryRecoQARefDefStorage(0),
	fChkRecoRunGlobalQA(0),
	fChkRecoRunPlaneEff(0),
	fChkRecoWriteESDf(0),
	fChkRecoWriteAlignment(0),
	fChkRecoCleanESD(0),
	fDetectors(0),
	fEntryLogbookHost(0),
	fEntryLogbookPort(0),
	fEntryLogbookDB(0),
	fEntryLogbookUser(0),
	fEntryLogbookPass(0)
{
	SetCleanup(kDeepCleanup);
	SetWindowName(window_title);
	
	 // Tab Preferences
	fTab = new TGTab(this);
  TGCompositeFrame* tab1 = fTab->AddTab("Server");
  TGCompositeFrame* tab2 = fTab->AddTab("Reconstruction");
  TGCompositeFrame* tab3 = fTab->AddTab("Logbook");
	
	SetupServerTab(tab1);
	SetupRecoTab(tab2);
	SetupLogbookTab(tab3);
	
	// dialog buttons
	TGHorizontalFrame* hf = new TGHorizontalFrame(this);
	TGTextButton* btDefaults = new TGTextButton(hf, "Restore Defaults");
	TGTextButton* btCancel = new TGTextButton(hf, "Cancel");
	TGTextButton* btApply = new TGTextButton(hf, "Apply");
	
	btDefaults->Connect("Clicked()", "AliEventServerPreferencesWindow", this, "onRestoreDefaults()");
	btCancel->Connect("Clicked()", "AliEventServerPreferencesWindow", this, "onCancel()");
	btApply->Connect("Clicked()", "AliEventServerPreferencesWindow", this, "onApply()");
	
	hf->AddFrame(btDefaults, new  TGLayoutHints(kLHintsNormal) );
	hf->AddFrame(btApply, new  TGLayoutHints(kLHintsRight) );
	hf->AddFrame(btCancel, new  TGLayoutHints(kLHintsRight) );
	
	
	// Add Tab Widget and Dialog buttons to Main Window
	AddFrame(fTab, new TGLayoutHints(kLHintsNormal | kLHintsExpandX | kLHintsExpandY));
	AddFrame(hf,  new TGLayoutHints(kLHintsNormal | kLHintsExpandX));
	
	if(ReadSettings()!=0) Info("AliEventServerPreferencesWindow::AliEventServerPreferencesWindow", "Problems while reading settings from file!\n");
	
	SetMWMHints(kMWMDecorAll | kMWMDecorResizeH  | kMWMDecorMaximize |
                              kMWMDecorMinimize | kMWMDecorMenu,
               kMWMFuncAll  | kMWMFuncResize    | kMWMFuncMaximize |
                              kMWMFuncMinimize,
               kMWMInputModeless);
               
	MapSubwindows();
  Resize();
  MapWindow();
  
  gClient->WaitFor(this);
}

AliEventServerPreferencesWindow::~AliEventServerPreferencesWindow()
{
	delete fDetectors; 
	fDetectors=0;
}

/*********************/
/* Server Settings Tab  */
/*********************/
void AliEventServerPreferencesWindow::SetupServerTab(TGCompositeFrame* tab)
{
	TGHorizontalFrame *t2hf = new TGHorizontalFrame(tab, 1, 20);
	TGLabel* lbServerHost = new TGLabel(t2hf, "Host:");
	t2hf->AddFrame(lbServerHost, new TGLayoutHints(kLHintsCenterY));
	
	fEntryServerHost = new TGTextEntry(t2hf, "tcp://*");
	t2hf->AddFrame(fEntryServerHost, new TGLayoutHints(kLHintsNormal | kLHintsExpandX));
	
	TGLabel* lbServerPort = new TGLabel(t2hf, "Port:");
	t2hf->AddFrame(lbServerPort, new TGLayoutHints(kLHintsCenterY));
	
	fEntryServerPort = new TGTextEntry(t2hf, "5024");
	fEntryServerPort->SetMaxLength(5);
	t2hf->AddFrame(fEntryServerPort, new TGLayoutHints(kLHintsNormal));
	
	tab->AddFrame(t2hf, new TGLayoutHints(kLHintsNormal | kLHintsExpandX));
}
	
  /*******************/
	/* Reco Settings Tab */
	/*******************/
void AliEventServerPreferencesWindow::SetupRecoTab(TGCompositeFrame* tab)
{
	/* CDB Manager Group Frame */
	TGGroupFrame* grCDBMan = new TGGroupFrame(tab, "CDB Manager");
	// Default Storage Frame
	TGHorizontalFrame *hfCDBDefaultStorage = new TGHorizontalFrame(grCDBMan, 1, 20);
	TGLabel* lbCDBDefaultStorage = new TGLabel(hfCDBDefaultStorage, "Default Storage:");
	fEntryCDBDefaultStorage = new TGTextEntry(hfCDBDefaultStorage, "local:///local/cdb");
	
	hfCDBDefaultStorage->AddFrame(lbCDBDefaultStorage, new TGLayoutHints(kLHintsCenterY));
	hfCDBDefaultStorage->AddFrame(fEntryCDBDefaultStorage, new TGLayoutHints(kLHintsNormal | kLHintsExpandX));
	
	// Specific Storage 1 Frame
	TGHorizontalFrame *hfCDBSpecStorage1 = new TGHorizontalFrame(grCDBMan, 1, 20);
	TGLabel* lbCDBSpecStorage1 = new TGLabel(hfCDBSpecStorage1, "Specific Storage:");
	fEntryCDBSpecStoragePath1 = new TGTextEntry(hfCDBSpecStorage1, "GRP/GRP/Data");
	fEntryCDBSpecStorageValue1 = new TGTextEntry(hfCDBSpecStorage1, "");
	
	hfCDBSpecStorage1->AddFrame(lbCDBSpecStorage1, new TGLayoutHints(kLHintsCenterY));
	hfCDBSpecStorage1->AddFrame(fEntryCDBSpecStoragePath1, new TGLayoutHints(kLHintsNormal | kLHintsExpandX));
	hfCDBSpecStorage1->AddFrame(fEntryCDBSpecStorageValue1, new TGLayoutHints(kLHintsNormal | kLHintsExpandX));
	
	// Specific Storage 2 Frame
	TGHorizontalFrame *hfCDBSpecStorage2 = new TGHorizontalFrame(grCDBMan, 1, 20);
	TGLabel* lbCDBSpecStorage2 = new TGLabel(hfCDBSpecStorage2, "Specific Storage:");
	fEntryCDBSpecStoragePath2 = new TGTextEntry(hfCDBSpecStorage2, "GRP/CTP/Config");
  fEntryCDBSpecStorageValue2 = new TGTextEntry(hfCDBSpecStorage2, "");
	
	hfCDBSpecStorage2->AddFrame(lbCDBSpecStorage2, new TGLayoutHints(kLHintsCenterY));
	hfCDBSpecStorage2->AddFrame(fEntryCDBSpecStoragePath2, new TGLayoutHints(kLHintsNormal | kLHintsExpandX));
	hfCDBSpecStorage2->AddFrame(fEntryCDBSpecStorageValue2, new TGLayoutHints(kLHintsNormal | kLHintsExpandX));
	
	// Specific Storage 3 Frame
	TGHorizontalFrame *hfCDBSpecStorage3 = new TGHorizontalFrame(grCDBMan, 1, 20);
	TGLabel* lbCDBSpecStorage3 = new TGLabel(hfCDBSpecStorage3, "Specific Storage:");
	fEntryCDBSpecStoragePath3 = new TGTextEntry(hfCDBSpecStorage3, "");
	fEntryCDBSpecStorageValue3 = new TGTextEntry(hfCDBSpecStorage3, "");
	
	hfCDBSpecStorage3->AddFrame(lbCDBSpecStorage3, new TGLayoutHints(kLHintsCenterY));
	hfCDBSpecStorage3->AddFrame(fEntryCDBSpecStoragePath3, new TGLayoutHints(kLHintsNormal | kLHintsExpandX));
	hfCDBSpecStorage3->AddFrame(fEntryCDBSpecStorageValue3, new TGLayoutHints(kLHintsNormal | kLHintsExpandX));
	
	grCDBMan->AddFrame(hfCDBDefaultStorage, new TGLayoutHints(kLHintsNormal | kLHintsExpandX));
	grCDBMan->AddFrame(hfCDBSpecStorage1, new TGLayoutHints(kLHintsNormal | kLHintsExpandX));
	grCDBMan->AddFrame(hfCDBSpecStorage2, new TGLayoutHints(kLHintsNormal | kLHintsExpandX));
	grCDBMan->AddFrame(hfCDBSpecStorage3, new TGLayoutHints(kLHintsNormal | kLHintsExpandX));
	
	tab->AddFrame(grCDBMan, new TGLayoutHints(kLHintsNormal | kLHintsExpandX));
	grCDBMan->Resize();
	
	/* Reconstruction Settings Group Frame */
	TGGroupFrame* grRecoFrame = new TGGroupFrame(tab, "Reconstruction");
	// SetRunQA
	TGHorizontalFrame *hfRecoRunQA = new TGHorizontalFrame(grRecoFrame, 1, 20);
	TGLabel* lbRecoRunQA = new TGLabel(hfRecoRunQA, "Run QA:");
	fEntryRecoRunQA = new TGTextEntry(hfRecoRunQA, ":");
	
	hfRecoRunQA->AddFrame( lbRecoRunQA, new TGLayoutHints(kLHintsCenterY));
	hfRecoRunQA->AddFrame( fEntryRecoRunQA, new TGLayoutHints(kLHintsNormal | kLHintsExpandX));
	
	// QARef Default Storage
	TGHorizontalFrame *hfRecoQARefDefStorage = new TGHorizontalFrame(grRecoFrame, 1, 20);
	TGLabel* lbRecoQARefDefStorage= new TGLabel(hfRecoQARefDefStorage, "QARef Default Storage:");
	fEntryRecoQARefDefStorage = new TGTextEntry(hfRecoQARefDefStorage, "local://$ALICE_ROOT/QAref");
	
	hfRecoQARefDefStorage->AddFrame( lbRecoQARefDefStorage, new TGLayoutHints(kLHintsCenterY));
	hfRecoQARefDefStorage->AddFrame( fEntryRecoQARefDefStorage, new TGLayoutHints(kLHintsNormal | kLHintsExpandX));
	
	// RunGlobalQA
	fChkRecoRunGlobalQA = new TGCheckButton(grRecoFrame, "Run Global QA");
	
	// Plane Efficiency Evalution
	fChkRecoRunPlaneEff= new TGCheckButton(grRecoFrame, "Plane Efficiency Evaluation");
	
	// Write ESD Friend
	fChkRecoWriteESDf= new TGCheckButton(grRecoFrame, "Write ESD friend");
	
	// Write Alignement Data
	fChkRecoWriteAlignment= new TGCheckButton(grRecoFrame, "Write Alignment Data");
	
	// Clean ESD
	fChkRecoCleanESD= new TGCheckButton(grRecoFrame, "Clean ESD");
	
	// Participating Detectors
	TGGroupFrame* grRecoDetFrame = new TGGroupFrame(grRecoFrame, "Participating Detectors");
	TGHorizontalFrame *hfRecoDets = new TGHorizontalFrame(grRecoDetFrame, 1, 20);
	TGTextButton* btRecoCheckAll = new TGTextButton(hfRecoDets, "Check All");
	TGTextButton* btRecoUnCheckAll = new TGTextButton(hfRecoDets, "Uncheck All");
	
	btRecoCheckAll->Connect("Clicked()", "AliEventServerPreferencesWindow", this, "onCheckAllDetectors()");
	btRecoUnCheckAll->Connect("Clicked()", "AliEventServerPreferencesWindow", this, "onUnCheckAllDetectors()");
	
	hfRecoDets->AddFrame( btRecoCheckAll, new TGLayoutHints(kLHintsNormal));
	hfRecoDets->AddFrame( btRecoUnCheckAll, new TGLayoutHints(kLHintsNormal));
	
	grRecoDetFrame->AddFrame(hfRecoDets, new TGLayoutHints(kLHintsNormal | kLHintsExpandX));
	
	// get detectors names from AliReconstruction
	const char** detnames= AliReconstruction::GetDetectorNames();
	TGCompositeFrame* cfDetectors = new TGCompositeFrame(grRecoDetFrame);
	cfDetectors->SetLayoutManager(new TGMatrixLayout(cfDetectors, 0, 4, 10));
	
	fDetectors = new TObjArray;
	for(int i=0;  i< AliReconstruction::kNDetectors; i++){
		TGCheckButton* chkRecoDet = new TGCheckButton(cfDetectors, detnames[i]);
		cfDetectors->AddFrame(chkRecoDet, new TGLayoutHints(kLHintsNormal));
		fDetectors->Add(chkRecoDet);
	}
	
	grRecoDetFrame->AddFrame(cfDetectors,  new TGLayoutHints(kLHintsNormal | kLHintsExpandX));
	grRecoDetFrame->Resize();
	
	grRecoFrame->AddFrame(hfRecoRunQA, new TGLayoutHints(kLHintsNormal | kLHintsExpandX));
	grRecoFrame->AddFrame(hfRecoQARefDefStorage, new TGLayoutHints(kLHintsNormal | kLHintsExpandX));
	grRecoFrame->AddFrame(fChkRecoRunGlobalQA, new TGLayoutHints(kLHintsNormal | kLHintsExpandX));
	grRecoFrame->AddFrame(fChkRecoRunPlaneEff, new TGLayoutHints(kLHintsNormal | kLHintsExpandX));
	grRecoFrame->AddFrame(fChkRecoWriteESDf, new TGLayoutHints(kLHintsNormal | kLHintsExpandX));
	grRecoFrame->AddFrame(fChkRecoWriteAlignment, new TGLayoutHints(kLHintsNormal | kLHintsExpandX));
	grRecoFrame->AddFrame(fChkRecoCleanESD, new TGLayoutHints(kLHintsNormal | kLHintsExpandX));
	grRecoFrame->AddFrame(grRecoDetFrame, new TGLayoutHints(kLHintsNormal | kLHintsExpandX));
	
	tab->AddFrame(grRecoFrame, new TGLayoutHints(kLHintsNormal | kLHintsExpandX | kLHintsExpandY));
	grRecoFrame->Resize();
}

	
/***********************/
/* Logbook Settings Tab  */
/***********************/
void AliEventServerPreferencesWindow::SetupLogbookTab(TGCompositeFrame* tab)
{
	// host and port
	TGCompositeFrame* cfLogbook = new TGCompositeFrame(tab);
	cfLogbook->SetLayoutManager(new TGMatrixLayout(cfLogbook, 0, 2));
	
	TGLabel* lbLogbookHost = new TGLabel(cfLogbook, "Host:");
	fEntryLogbookHost = new TGTextEntry(cfLogbook, "pcaldbl501");
	fEntryLogbookHost->Resize(150,0);
	TGLabel* lbLogbookPort = new TGLabel(cfLogbook, "Port:");
	fEntryLogbookPort = new TGTextEntry(cfLogbook, "3306");
	fEntryLogbookPort->SetMaxLength(5);
	
	// database name
	TGLabel* lbLogbookDB = new TGLabel(cfLogbook, "Database:");
	fEntryLogbookDB = new TGTextEntry(cfLogbook, "logbook");
	fEntryLogbookDB->Resize(150,0);
	
  // username
	TGLabel* lbLogbookUser = new TGLabel(cfLogbook, "User:");
	fEntryLogbookUser = new TGTextEntry(cfLogbook, "dqm");
	fEntryLogbookUser->Resize(150,0);
			
	// password
	TGLabel* lbLogbookPass = new TGLabel(cfLogbook, "Password:");
	fEntryLogbookPass = new TGTextEntry(cfLogbook, "dqm123");
	fEntryLogbookPass->SetEchoMode(TGTextEntry::kPassword);
	fEntryLogbookPass->Resize(150,0);
	
	cfLogbook->AddFrame(lbLogbookHost, new TGLayoutHints(kLHintsCenterY));
	cfLogbook->AddFrame(fEntryLogbookHost, new TGLayoutHints(kLHintsNormal | kLHintsExpandX));
	cfLogbook->AddFrame(lbLogbookPort, new TGLayoutHints(kLHintsCenterY));
	cfLogbook->AddFrame(fEntryLogbookPort, new TGLayoutHints(kLHintsNormal));
	cfLogbook->AddFrame(lbLogbookDB, new TGLayoutHints(kLHintsCenterY));
	cfLogbook->AddFrame(fEntryLogbookDB, new TGLayoutHints(kLHintsNormal));
	cfLogbook->AddFrame(lbLogbookUser, new TGLayoutHints(kLHintsCenterY));
	cfLogbook->AddFrame(fEntryLogbookUser, new TGLayoutHints(kLHintsNormal));
	cfLogbook->AddFrame(lbLogbookPass, new TGLayoutHints(kLHintsCenterY));
	cfLogbook->AddFrame(fEntryLogbookPass, new TGLayoutHints(kLHintsNormal));
	
	tab->AddFrame(cfLogbook, new TGLayoutHints(kLHintsNormal));
}

void AliEventServerPreferencesWindow::onRestoreDefaults()
{
	int retCode;
	const char* messageInfo = "Are you sure you want to restore the default settings?\n";
 
	new TGMsgBox(gClient->GetRoot(), this, "Restore Default Settings?", messageInfo, kMBIconQuestion, kMBNo|kMBYes, &retCode);
  
  if(retCode==kMBNo) return;
	RestoreDefaults();
}

void AliEventServerPreferencesWindow::onCancel()
{
	SendCloseMessage();
}

void AliEventServerPreferencesWindow::onApply()
{
	if(WriteSettings()==0){// success!
  	const char* messageInfo = "The changes were saved.\n"
  	"Notice: They will take effect when a new run starts or\n"
  	"after the reconstruction is restarted (if it is already running)";
  	 
  	new TGMsgBox(gClient->GetRoot(), this, "New Settings Notice", messageInfo, kMBIconExclamation, kMBOk);
  	
  	SendCloseMessage();
	}
	else{
		const char* messageInfo = "The changes could not be saved!\n"
  	"Check if you have permissions to write to that location.";
  	 
  	new TGMsgBox(gClient->GetRoot(), this, "New Settings Notice", messageInfo, kMBIconExclamation, kMBOk);
	}
		
}

void AliEventServerPreferencesWindow::onCheckAllDetectors()
{
	for(int i=0; i<fDetectors->GetEntries(); ++i){
		TGCheckButton* bt = (TGCheckButton*)fDetectors->At(i);
		bt->SetOn(kTRUE);
	}
}

void AliEventServerPreferencesWindow::onUnCheckAllDetectors()
{
	for(int i=0; i<fDetectors->GetEntries(); ++i){
		TGCheckButton* bt = (TGCheckButton*)fDetectors->At(i);
		bt->SetOn(kFALSE);
	}
}

Int_t AliEventServerPreferencesWindow::ReadSettings()
{
	TEnv settings;
	int readStatus = settings.ReadFile(Form("%s/MONITOR/%s", gSystem->Getenv("ALICE_ROOT"), ALIEVENTSERVER_CONF), kEnvLocal);
	//check if there was an error reading the file
	if(readStatus!=0) return readStatus;
	
	// server settings
	fEntryServerHost->SetText( settings.GetValue("server.host", DEFAULT_SERVER_HOST), kFALSE);
	fEntryServerPort->SetText(  Form("%d", settings.GetValue("server.port", DEFAULT_SERVER_PORT)), kFALSE);
	
	// reco settings
	fEntryCDBDefaultStorage->SetText( settings.GetValue( "cdb.defaultStorage", DEFAULT_CDB_STORAGE), kFALSE);
	fEntryCDBSpecStoragePath1->SetText( settings.GetValue( "cdb.specificStoragePath1", DEFAULT_CDB_SPEC_STORAGE_PATH1), kFALSE);
	fEntryCDBSpecStorageValue1->SetText( settings.GetValue( "cdb.specificStorageValue1", DEFAULT_CDB_SPEC_STORAGE_VALUE1), kFALSE);
	fEntryCDBSpecStoragePath2->SetText( settings.GetValue( "cdb.specificStoragePath2", DEFAULT_CDB_SPEC_STORAGE_PATH2), kFALSE);
	fEntryCDBSpecStorageValue2->SetText(settings.GetValue( "cdb.specificStorageValue2", DEFAULT_CDB_SPEC_STORAGE_VALUE2), kFALSE);
	fEntryCDBSpecStoragePath3->SetText( settings.GetValue( "cdb.specificStoragePath3", DEFAULT_CDB_SPEC_STORAGE_PATH3), kFALSE);
	fEntryCDBSpecStorageValue3->SetText(settings.GetValue( "cdb.specificStorageValue3",  DEFAULT_CDB_SPEC_STORAGE_VALUE3), kFALSE);
	fEntryRecoRunQA->SetText( settings.GetValue( "qa.runDetectors", DEFAULT_QA_RUN), kFALSE );
	fEntryRecoQARefDefStorage->SetText(settings.GetValue( "qa.defaultStorage",DEFAULT_QAREF_STORAGE), kFALSE);
	fChkRecoRunGlobalQA->SetOn(settings.GetValue( "qa.runGlobal", DEFAULT_QA_RUN_GLOBAL), kFALSE);
	fChkRecoRunPlaneEff->SetOn(settings.GetValue( "reco.runPlaneEff", DEFAULT_RECO_RUN_PLANE_EFF), kFALSE);
	fChkRecoWriteESDf->SetOn(settings.GetValue( "reco.writeESDfriend", DEFAULT_RECO_WRITE_ESDF), kFALSE);
	fChkRecoWriteAlignment->SetOn(settings.GetValue( "reco.writeAlignment",DEFAULT_RECO_WRITE_ALIGN), kFALSE);
	fChkRecoCleanESD->SetOn(settings.GetValue( "reco.cleanESD",DEFAULT_RECO_CLEAN_ESD), kFALSE);
	
	// parse reco run detectors from string
	TString strRunDetectors(settings.GetValue( "reco.detectors", DEFAULT_RECO_DETECTORS) );
	TObjArray* arrRecoDets = strRunDetectors.Tokenize(" ");
	
	for(int i=0; i<arrRecoDets->GetEntries(); ++i){
		TObjString* objStr = (TObjString*)arrRecoDets->At(i);
		if(objStr->GetString().BeginsWith("-")){ // detector is disabled
			for(int j=0; j<fDetectors->GetEntries();++j){
				TGCheckButton* bt = (TGCheckButton*)fDetectors->At(j);
				TString btTitle(Form("-%s",bt->GetTitle()));
				if(btTitle.CompareTo(objStr->GetString())==0) {// match
					bt->SetOn(kFALSE); // uncheck detector button
					break;
				}
			}
		}
		else if(objStr->GetString().CompareTo("ALL")==0){
			onCheckAllDetectors();
		}
		else{ // detector is enabled
			for(int j=0; j<fDetectors->GetEntries();++j){
				TGCheckButton* bt = (TGCheckButton*)fDetectors->At(j);
				TString btTitle(bt->GetTitle());
				if(btTitle.CompareTo(objStr->GetString())==0) {// match
					bt->SetOn(kTRUE); // check detector
					break;
				}
			}
		}
	}
	
	// logbook settings
	fEntryLogbookHost->SetText( settings.GetValue("logbook.host", DEFAULT_LOGBOOK_HOST), kFALSE);
	fEntryLogbookPort->SetText( Form("%d", settings.GetValue("logbook.port", DEFAULT_LOGBOOK_PORT)), kFALSE);
	fEntryLogbookDB->SetText( settings.GetValue("logbook.db", DEFAULT_LOGBOOK_DB), kFALSE);
	fEntryLogbookUser->SetText( settings.GetValue("logbook.user", DEFAULT_LOGBOOK_USER), kFALSE);
	fEntryLogbookPass->SetText( settings.GetValue("logbook.pass", DEFAULT_LOGBOOK_PASS), kFALSE);

	return readStatus;
}

// write settings to a rootrc file
// returns 0 in case of success, -1 in case of error
Int_t AliEventServerPreferencesWindow::WriteSettings()
{
	TEnv settings;
	// server settings
	settings.SetValue("server.host", fEntryServerHost->GetText());
	settings.SetValue("server.port", TString(fEntryServerPort->GetText()).Atoi());
	
	// reco settings
	settings.SetValue( "cdb.defaultStorage", fEntryCDBDefaultStorage->GetText());
	settings.SetValue( "cdb.specificStoragePath1", fEntryCDBSpecStoragePath1->GetText());
	settings.SetValue( "cdb.specificStorageValue1", fEntryCDBSpecStorageValue1->GetText());
	settings.SetValue( "cdb.specificStoragePath2", fEntryCDBSpecStoragePath2->GetText());
	settings.SetValue( "cdb.specificStorageValue2", fEntryCDBSpecStorageValue2->GetText());
	settings.SetValue( "cdb.specificStoragePath3", fEntryCDBSpecStoragePath3->GetText());
	settings.SetValue( "cdb.specificStorageValue3", fEntryCDBSpecStorageValue3->GetText());
	settings.SetValue( "qa.runDetectors", fEntryRecoRunQA->GetText());
	settings.SetValue( "qa.defaultStorage", fEntryRecoQARefDefStorage->GetText());
	settings.SetValue( "qa.runGlobal", fChkRecoRunGlobalQA->IsOn());
	settings.SetValue( "reco.runPlaneEff", fChkRecoRunPlaneEff->IsOn());
	settings.SetValue( "reco.writeESDfriend", fChkRecoWriteESDf->IsOn());
	settings.SetValue( "reco.writeAlignment", fChkRecoWriteAlignment->IsOn());
	settings.SetValue( "reco.cleanESD", fChkRecoCleanESD->IsOn());
	
	// will write reco run detectors as a single string
	TObjArray checkedDetectors;
	TObjArray uncheckedDetectors;
	TString strRunDetectors;
	
	for(int i=0; i<fDetectors->GetEntries(); ++i){
		TGCheckButton* bt = (TGCheckButton*)fDetectors->At(i);
		if(bt->IsOn()){
			checkedDetectors.Add(bt);
		}
		else{
			uncheckedDetectors.Add(bt);
		}
	}
	
	int nChkDets = checkedDetectors.GetEntries();
	int nUnChkDets = uncheckedDetectors.GetEntries();
	if(nChkDets>=nUnChkDets){
		strRunDetectors="ALL ";
		for(int i=0; i<nUnChkDets; ++i ){
			strRunDetectors.Append( Form(" -%s",uncheckedDetectors.At(i)->GetTitle()) );
		}
	}
	else {
		for(int i=0; i<nChkDets; ++i ){
			strRunDetectors.Append( Form(" %s",checkedDetectors.At(i)->GetTitle()) );
		}
	}
	
	settings.SetValue( "reco.detectors", strRunDetectors.Data());
	
	// logbook settings
	settings.SetValue("logbook.host", fEntryLogbookHost->GetText());
	settings.SetValue("logbook.port", TString(fEntryLogbookPort->GetText()).Atoi());
	settings.SetValue("logbook.db", fEntryLogbookDB->GetText());
	settings.SetValue("logbook.user", fEntryLogbookUser->GetText());
	settings.SetValue("logbook.pass", fEntryLogbookPass->GetText());
	
	return settings.WriteFile(Form("%s/MONITOR/%s", gSystem->Getenv("ALICE_ROOT"), ALIEVENTSERVER_CONF));
}

void AliEventServerPreferencesWindow::RestoreDefaults()
{
	// server settings
	fEntryServerHost->SetText( DEFAULT_SERVER_HOST, kFALSE);
	fEntryServerPort->SetText(  Form("%d", DEFAULT_SERVER_PORT), kFALSE);
	
	// reco settings
	fEntryCDBDefaultStorage->SetText( DEFAULT_CDB_STORAGE, kFALSE);
	fEntryCDBSpecStoragePath1->SetText( DEFAULT_CDB_SPEC_STORAGE_PATH1, kFALSE);
	fEntryCDBSpecStorageValue1->SetText( DEFAULT_CDB_SPEC_STORAGE_VALUE1, kFALSE);
	fEntryCDBSpecStoragePath2->SetText( DEFAULT_CDB_SPEC_STORAGE_PATH2, kFALSE);
	fEntryCDBSpecStorageValue2->SetText(DEFAULT_CDB_SPEC_STORAGE_VALUE2, kFALSE);
	fEntryCDBSpecStoragePath3->SetText( DEFAULT_CDB_SPEC_STORAGE_PATH3, kFALSE);
	fEntryCDBSpecStorageValue3->SetText( DEFAULT_CDB_SPEC_STORAGE_VALUE3, kFALSE);
	fEntryRecoRunQA->SetText( DEFAULT_QA_RUN, kFALSE );
	fEntryRecoQARefDefStorage->SetText(DEFAULT_QAREF_STORAGE, kFALSE);
	fChkRecoRunGlobalQA->SetOn(DEFAULT_QA_RUN_GLOBAL, kFALSE);
	fChkRecoRunPlaneEff->SetOn(DEFAULT_RECO_RUN_PLANE_EFF, kFALSE);
	fChkRecoWriteESDf->SetOn(DEFAULT_RECO_WRITE_ESDF, kFALSE);
	fChkRecoWriteAlignment->SetOn(DEFAULT_RECO_WRITE_ALIGN, kFALSE);
	fChkRecoCleanESD->SetOn(DEFAULT_RECO_CLEAN_ESD, kFALSE);
	
	// parse reco run detectors from string
	TString strRunDetectors( DEFAULT_RECO_DETECTORS );
	TObjArray* arrRecoDets = strRunDetectors.Tokenize(" ");
	
	for(int i=0; i<arrRecoDets->GetEntries(); ++i){
		TObjString* objStr = (TObjString*)arrRecoDets->At(i);
		if(objStr->GetString().BeginsWith("-")){ // detector is disabled
			for(int j=0; j<fDetectors->GetEntries();++j){
				TGCheckButton* bt = (TGCheckButton*)fDetectors->At(j);
				TString btTitle(Form("-%s",bt->GetTitle()));
				if(btTitle.CompareTo(objStr->GetString())==0) {// match
					bt->SetOn(kFALSE); // uncheck detector button
					break;
				}
			}
		}
		else if(objStr->GetString().CompareTo("ALL")==0){
			onCheckAllDetectors();
		}
		else{ // detector is enabled
			for(int j=0; j<fDetectors->GetEntries();++j){
				TGCheckButton* bt = (TGCheckButton*)fDetectors->At(j);
				TString btTitle(bt->GetTitle());
				if(btTitle.CompareTo(objStr->GetString())==0) {// match
					bt->SetOn(kTRUE); // check detector
					break;
				}
			}
		}
	}
	
	// logbook settings
	fEntryLogbookHost->SetText( DEFAULT_LOGBOOK_HOST, kFALSE);
	fEntryLogbookPort->SetText( Form("%d", DEFAULT_LOGBOOK_PORT), kFALSE);
	fEntryLogbookDB->SetText( DEFAULT_LOGBOOK_DB, kFALSE);
	fEntryLogbookUser->SetText( DEFAULT_LOGBOOK_USER, kFALSE);
	fEntryLogbookPass->SetText( DEFAULT_LOGBOOK_PASS, kFALSE);
}
