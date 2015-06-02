#ifndef AliStorageAdministratorPanel_H
#define AliStorageAdministratorPanel_H

#include "AliZMQManager.h"

#include <TGFrame.h>
#include <TGLabel.h>
#include <TThread.h>

class AliStorageAdministratorPanel : public TGMainFrame
{

public:
	AliStorageAdministratorPanel();
	virtual ~AliStorageAdministratorPanel();
private:
	bool fPanelQuited;
	
	//gui components and methods
	TGLabel *fConnectionLabel;
	TGLabel *fDataLabel;
	TGLabel *fSavingLabel;
	TGLabel *fCurrentSizeLabel;
	TGLabel *fMaxSizeLabel;
	TGLabel *fMaxOccupationLabel;
	TGLabel *fRemoveEventsLabel;
	TGLabel *fEventsInChunkLabel;
	
	void InitWindow();
	void SetupThreadsFrame();
	void SetupToolbar();
	void SetupFixedMenuBar();
	void SetupDockableMenuBar();

	//set labels method
	void SetLabel(TGLabel *label, int option);
	void SetLabelValue(TGLabel *label, long value, int option);
	
	//handle different actions
	void onExit();
	void onOption1();
	void onOption2();
	void onOption3();

	void onServerListEvents();
	void onServerMarkEvent();
    void onServerMarkAllEvents();
	void onServerGetEvent();
	void onClientSetParams();
	void onServerGetNextEvent();
	void onServerGetLastEvent();

	//client params
	int fMaxStorageSize;
	int fMaxOccupation;
	int fRemoveEvents;
	int fEventsInChunk;
	
	//handle different messages
	virtual Bool_t ProcessMessage(Long_t msg, Long_t parm1, Long_t);
	void CloseWindow();
	 
	//socket connection
	TThread *fCommunicationThread;
	storageSockets fCommunicationSocket;
	static void Dispatch(void *arg){static_cast<AliStorageAdministratorPanel*>(arg)->CheckStateHandle();}
	void CheckStateHandle();
	void CheckClientState(int option);
	
	storageSockets fServerSocket;//socket for two-way communication with AliStorageServerThread
	AliZMQManager *fEventManager;

	AliStorageAdministratorPanel(const AliStorageAdministratorPanel&);
	AliStorageAdministratorPanel& operator=(const AliStorageAdministratorPanel&);
	
	ClassDef(AliStorageAdministratorPanel,0);
};

#endif
