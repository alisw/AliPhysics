#ifndef AliStorageAdministratorPanelSetStorageParams_H
#define AliStorageAdministratorPanelSetStorageParams_H

#include "AliZMQManager.h"

#include <TGLabel.h>
#include <TGNumberEntry.h>
#include <TG3DLine.h>

class AliStorageAdministratorPanelSetStorageParams : public TGMainFrame
{
public:
	static AliStorageAdministratorPanelSetStorageParams* GetInstance();
	void Setup(storageSockets socket,
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
	
	storageSockets fClientSocket;
	AliZMQManager *fEventManager;

	AliStorageAdministratorPanelSetStorageParams(const AliStorageAdministratorPanelSetStorageParams&);
	AliStorageAdministratorPanelSetStorageParams& operator=(const AliStorageAdministratorPanelSetStorageParams&);
	
	ClassDef(AliStorageAdministratorPanelSetStorageParams,0);
};

#endif
