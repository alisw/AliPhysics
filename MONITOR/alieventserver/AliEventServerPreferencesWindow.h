// Author: Mihai Niculescu 2013

/**************************************************************************
 * Copyright(c) 1998-2013, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEventServerPreferencesWindow_H
#define AliEventServerPreferencesWindow_H


#include <TGFrame.h>

class TObjArray;
class TGCheckButton;
class TGTab;
class TGTextEntry;
class TGWindow;

class AliEventServerPreferencesWindow : public TGTransientFrame
{
public:
	AliEventServerPreferencesWindow(const TGWindow* parent, const char* window_title);
	virtual ~AliEventServerPreferencesWindow();
	
	// SLOTS
	void onApply();
	void onCancel();
	void onRestoreDefaults();
	void onCheckAllDetectors();
	void onUnCheckAllDetectors();
	
protected:
	Int_t ReadSettings();
	Int_t WriteSettings();
	void RestoreDefaults();
	
	TGTab* fTab;
	/* Server Settings Tab */
	TGTextEntry* fEntryServerHost;
	TGTextEntry* fEntryServerPort;
	
	/* Reconstruction Settings Tab */
	// CDB Manager widgets
	TGTextEntry* fEntryCDBDefaultStorage;
	TGTextEntry* fEntryCDBSpecStoragePath1;
	TGTextEntry* fEntryCDBSpecStorageValue1;
	TGTextEntry* fEntryCDBSpecStoragePath2;
	TGTextEntry* fEntryCDBSpecStorageValue2;
	TGTextEntry* fEntryCDBSpecStoragePath3;
	TGTextEntry* fEntryCDBSpecStorageValue3;
	// Reconstruction widgets
	TGTextEntry* fEntryRecoRunQA;
	TGTextEntry* fEntryRecoQARefDefStorage;
	TGCheckButton* fChkRecoRunGlobalQA;
	TGCheckButton* fChkRecoRunPlaneEff;
	TGCheckButton* fChkRecoWriteESDf;
	TGCheckButton* fChkRecoWriteAlignment;
	TGCheckButton* fChkRecoCleanESD;
	TObjArray* fDetectors; // array to pointers of TGCheckButtons for detectors
	
	/* Logbook Tab */
	TGTextEntry* fEntryLogbookHost;
	TGTextEntry* fEntryLogbookPort;
	TGTextEntry* fEntryLogbookDB;
	TGTextEntry* fEntryLogbookUser;
	TGTextEntry* fEntryLogbookPass;
	
private:
	void SetupServerTab(TGCompositeFrame* tab);
	void SetupRecoTab(TGCompositeFrame* tab);
	void SetupLogbookTab(TGCompositeFrame* tab);
				
	ClassDef(AliEventServerPreferencesWindow,0);
};
#endif
