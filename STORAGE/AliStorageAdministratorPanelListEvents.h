#ifndef AliStorageAdministratorPanelListEvents_H
#define AliStorageAdministratorPanelListEvents_H

#include "AliStorageTypes.h"
#include "AliStorageEventManager.h"

#include <TGLabel.h>
#include <TGNumberEntry.h>
#include <TGListBox.h>

namespace zmq
{
	class socket_t;
}

class AliStorageAdministratorPanelListEvents : public TGMainFrame
{
public:
	static AliStorageAdministratorPanelListEvents* GetInstance();
	void SetSocket(zmq::socket_t *socket);
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
	void onExit();

	std::vector<serverListStruct> fEventsListVector;

	virtual Bool_t ProcessMessage(Long_t msg, Long_t parm1, Long_t);
	void CloseWindow();

	zmq::socket_t *fServerSocket;
	AliStorageEventManager *fEventManager;
	
	AliStorageAdministratorPanelListEvents(const AliStorageAdministratorPanelListEvents&);
	AliStorageAdministratorPanelListEvents& operator=(const AliStorageAdministratorPanelListEvents&);
	
	ClassDef(AliStorageAdministratorPanelListEvents,0);
};

#endif
