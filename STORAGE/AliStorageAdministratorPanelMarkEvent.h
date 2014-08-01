#ifndef AliStorageAdministratorPanelMarkEvent_H
#define AliStorageAdministratorPanelMarkEvent_H

#include "AliStorageEventManager.h"

#include <TGLabel.h>
#include <TGNumberEntry.h>
#include <TG3DLine.h>

namespace zmq
{
	class socket_t;
}

class AliStorageAdministratorPanelMarkEvent : public TGMainFrame
{
public:
	static AliStorageAdministratorPanelMarkEvent* GetInstance();
	void SetSocket(zmq::socket_t *socket);
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

	zmq::socket_t *fServerSocket;
	AliStorageEventManager *fEventManager;
	
	AliStorageAdministratorPanelMarkEvent(const AliStorageAdministratorPanelMarkEvent&);
	AliStorageAdministratorPanelMarkEvent& operator=(const AliStorageAdministratorPanelMarkEvent&);
	
	ClassDef(AliStorageAdministratorPanelMarkEvent,0);
};

#endif
