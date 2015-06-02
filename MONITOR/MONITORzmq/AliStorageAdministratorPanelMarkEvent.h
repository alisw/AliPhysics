#ifndef AliStorageAdministratorPanelMarkEvent_H
#define AliStorageAdministratorPanelMarkEvent_H

#include "AliZMQManager.h"

#include <TGLabel.h>
#include <TGNumberEntry.h>
#include <TG3DLine.h>

class AliStorageAdministratorPanelMarkEvent : public TGMainFrame
{
public:
	static AliStorageAdministratorPanelMarkEvent* GetInstance();
private:
	AliStorageAdministratorPanelMarkEvent();
	virtual ~AliStorageAdministratorPanelMarkEvent();

	static AliStorageAdministratorPanelMarkEvent *fInstance;
	
	//gui components and methods
	TGLabel *fStatusLabel;
	TGNumberEntry *fRunNumberEntry;
	TGNumberEntry *fEventNumberEntry;
	
	void InitWindow();
	void onCloseButton();
	void onMarkButton();
	void onExit();

	virtual Bool_t ProcessMessage(Long_t msg, Long_t parm1, Long_t);
	void CloseWindow();

	storageSockets fServerSocket;
	AliZMQManager *fEventManager;
	
	AliStorageAdministratorPanelMarkEvent(const AliStorageAdministratorPanelMarkEvent&);
	AliStorageAdministratorPanelMarkEvent& operator=(const AliStorageAdministratorPanelMarkEvent&);
	
	ClassDef(AliStorageAdministratorPanelMarkEvent,0);
};

#endif
