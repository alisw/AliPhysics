#ifndef AliStorageAdministratorPanelListEvents_H
#define AliStorageAdministratorPanelListEvents_H

#include "AliStorageTypes.h"
#include "AliStorageEventManager.h"

#include <TGLabel.h>
#include <TGNumberEntry.h>
#include <TGListBox.h>
#include <TQObject.h>

class AliStorageAdministratorPanelListEvents : public TGMainFrame
{
public:
	static AliStorageAdministratorPanelListEvents* GetInstance();
    AliESDEvent* GetSelectedEvent(){return fCurrentEvent;}
	void onExit();
    void SelectedEvent(); //*SIGNAL*
private:
	AliStorageAdministratorPanelListEvents();
	virtual ~AliStorageAdministratorPanelListEvents();
	
	static AliStorageAdministratorPanelListEvents *fInstance;
	
    
    
	//gui components and methods
	TGLabel *fStatusLabel;
	TGNumberEntry *fRunNumberMinEntry;
	TGNumberEntry *fRunNumberMaxEntry;
	TGNumberEntry *fEventNumberMinEntry;
	TGNumberEntry *fEventNumberMaxEntry;
	TGNumberEntry *fMultiplicityMinEntry;
	TGNumberEntry *fMultiplicityMaxEntry;

	TGCheckButton *fProtonProtonCheckButton;
	TGCheckButton *fLeadLeadCheckButton;
	TGCheckButton *fTempCheckButton;
	TGCheckButton *fPermCheckButton;
	
	TGListBox *fEventsList;
	
	void InitWindow();
	void onCloseButton();
	void onGetListButton();
	void onMarkButton();
    void onLoadButton();

	std::vector<serverListStruct> fEventsListVector;

	virtual Bool_t ProcessMessage(Long_t msg, Long_t parm1, Long_t);
	void CloseWindow();

	storageSockets fServerSocket;
	AliStorageEventManager *fEventManager;
    
    AliESDEvent *fCurrentEvent;
	
	AliStorageAdministratorPanelListEvents(const AliStorageAdministratorPanelListEvents&);
	AliStorageAdministratorPanelListEvents& operator=(const AliStorageAdministratorPanelListEvents&);
	
	ClassDef(AliStorageAdministratorPanelListEvents,0);
};

#endif
