#ifndef AliStorageAdministratorPanelSetStorageParams_H
#define AliStorageAdministratorPanelSetStorageParams_H

#include "AliStorageEventManager.h"

#include <TGLabel.h>
#include <TGNumberEntry.h>
#include <TG3DLine.h>

namespace zmq
{
	class socket_t;
}

class AliStorageAdministratorPanelSetStorageParams : public TGMainFrame
{
public:
	static AliStorageAdministratorPanelSetStorageParams* GetInstance();
	void Setup(zmq::socket_t *socket,
		   int maxStorageSize,
		   int maxOccupation,
		   int removeEvents,
		   int eventsInChunk);
private:
	AliStorageAdministratorPanelSetStorageParams();
	virtual ~AliStorageAdministratorPanelSetStorageParams();

	static AliStorageAdministratorPanelSetStorageParams *fInstance;
	
	//gui components and methods
	TGLabel *fStatusLabel;
	TGNumberEntry *fMaxStorageSizeEntry;
	TGNumberEntry *fMaxOccupationEntry;
	TGNumberEntry *fRemoveEventsEntry;
	TGNumberEntry *fEventsInChunkEntry;
	
	void InitWindow();
	void onCloseButton();
	void onSetParamsButton();
	void onExit();

	virtual Bool_t ProcessMessage(Long_t msg, Long_t parm1, Long_t);
	void CloseWindow();
	
	zmq::socket_t *fClientSocket;
	AliStorageEventManager *fEventManager;

	AliStorageAdministratorPanelSetStorageParams(const AliStorageAdministratorPanelSetStorageParams&);
	AliStorageAdministratorPanelSetStorageParams& operator=(const AliStorageAdministratorPanelSetStorageParams&);
	
	ClassDef(AliStorageAdministratorPanelSetStorageParams,0);
};

#endif
