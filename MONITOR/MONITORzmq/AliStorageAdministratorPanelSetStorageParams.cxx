#include "AliStorageAdministratorPanelSetStorageParams.h"
#include "AliStorageTypes.h"

#include <iostream>

#include <TGFrame.h>
#include <TGButton.h>

using namespace std;

AliStorageAdministratorPanelSetStorageParams *AliStorageAdministratorPanelSetStorageParams::fInstance=0;

ClassImp(AliStorageAdministratorPanelSetStorageParams);

#define WINDOWS_WIDTH 200
#define WINDOWS_HEIGHT 200

enum BUTTON{
	BUTTON_CLOSE=1,
	BUTTON_SET 	
};

enum ENTRY{
	ENTRY_STORAGE_SIZE=1,
	ENTRY_OCCUPATION,
	ENTRY_REMOVE,
	ENTRY_EVENTS_CHUNK
};

AliStorageAdministratorPanelSetStorageParams::AliStorageAdministratorPanelSetStorageParams() :
	TGMainFrame(gClient->GetRoot(), 400, 400),
	fStatusLabel(0),
	fMaxStorageSizeEntry(0),
	fMaxOccupationEntry(0),
	fRemoveEventsEntry(0),
	fEventsInChunkEntry(0),
	fClientSocket(CLIENT_COMMUNICATION_REQ),
	fEventManager(0)
{
	fEventManager = AliZMQManager::GetInstance();
	InitWindow();
}

AliStorageAdministratorPanelSetStorageParams::~AliStorageAdministratorPanelSetStorageParams(){}

AliStorageAdministratorPanelSetStorageParams* AliStorageAdministratorPanelSetStorageParams::GetInstance()
{
	if(!fInstance){fInstance = new AliStorageAdministratorPanelSetStorageParams();}
	return fInstance;
}

void AliStorageAdministratorPanelSetStorageParams::Setup(storageSockets socket, int maxStorageSize, int maxOccupation, int removeEvents, int eventsInChunk)
{
	fClientSocket = socket;
	fMaxStorageSizeEntry->SetIntNumber(maxStorageSize);
	fMaxOccupationEntry->SetIntNumber(maxOccupation);
	fRemoveEventsEntry->SetIntNumber(removeEvents);
	fEventsInChunkEntry->SetIntNumber(eventsInChunk);
}

void AliStorageAdministratorPanelSetStorageParams::InitWindow()
{
	SetCleanup(kDeepCleanup);

	// max storage size
	AddFrame(new TGLabel(this,"Max storage size (MB):"),new TGLayoutHints(kLHintsLeft));

	fMaxStorageSizeEntry = new TGNumberEntry(this,
						     0,
						     6,
						     ENTRY_STORAGE_SIZE,
						     TGNumberFormat::kNESInteger,
						     TGNumberFormat::kNEAPositive,
						     TGNumberFormat::kNELNoLimits);
	AddFrame(fMaxStorageSizeEntry,new TGLayoutHints(kLHintsLeft));

	// max occupation
	AddFrame(new TGLabel(this,"Max occupation percent (%):"),new TGLayoutHints(kLHintsLeft));

	fMaxOccupationEntry = new TGNumberEntry(this,
						0,
						6,
						ENTRY_OCCUPATION,
						TGNumberFormat::kNESInteger,
						TGNumberFormat::kNEAPositive,
						TGNumberFormat::kNELNoLimits);
	AddFrame(fMaxOccupationEntry,new TGLayoutHints(kLHintsLeft));	

	// remove events percantage
	AddFrame(new TGLabel(this,"Remove events percentage (%):"),new TGLayoutHints(kLHintsLeft));

	fRemoveEventsEntry = new TGNumberEntry(this,
						0,
						6,
						ENTRY_OCCUPATION,
						TGNumberFormat::kNESInteger,
						TGNumberFormat::kNEAPositive,
						TGNumberFormat::kNELNoLimits);
	AddFrame(fRemoveEventsEntry,new TGLayoutHints(kLHintsLeft));

	// events in chunk
	AddFrame(new TGLabel(this,"Number of events in file:"),new TGLayoutHints(kLHintsLeft));

	fEventsInChunkEntry = new TGNumberEntry(this,
						0,
						6,
						ENTRY_OCCUPATION,
						TGNumberFormat::kNESInteger,
						TGNumberFormat::kNEAPositive,
						TGNumberFormat::kNELNoLimits);
	AddFrame(fEventsInChunkEntry,new TGLayoutHints(kLHintsLeft));	

	// status label
	fStatusLabel = new TGLabel(this,"");
	AddFrame(fStatusLabel,new TGLayoutHints(kLHintsExpandX | kLHintsLeft));

	AddFrame(new TGTextButton(this,"Close",BUTTON_CLOSE),
		 new TGLayoutHints(kLHintsLeft));

	AddFrame(new TGTextButton(this,"Set parameters",BUTTON_SET),
		 new TGLayoutHints(kLHintsRight));

	
	SetWindowName("Set Storage Parameters");
	MapSubwindows();
	Resize(WINDOWS_WIDTH,WINDOWS_HEIGHT);
	MapWindow();
}




void AliStorageAdministratorPanelSetStorageParams::onSetParamsButton()
{
	
	struct clientRequestStruct *requestMessage = new struct clientRequestStruct;	
        //get run and event number from TGNumberEntries
	requestMessage->messageType = REQUEST_SET_PARAMS;
	requestMessage->maxStorageSize = fMaxStorageSizeEntry->GetIntNumber()*1000000;
	requestMessage->maxOccupation = fMaxOccupationEntry->GetIntNumber();
	requestMessage->removeEvents = fRemoveEventsEntry->GetIntNumber();
	requestMessage->eventsInChunk = fEventsInChunkEntry->GetIntNumber();

	fEventManager->Send(requestMessage,fClientSocket);
    bool response;
    fEventManager->Get(&response,fClientSocket);
	
	if(response)
	{
		fStatusLabel->SetText("Params set");
		cout<<"ADMIN PANEL -- Params set succesfully"<<endl;
	}
	else
	{
		fStatusLabel->SetText("Couldn't set params");
		cout<<"ADMIN PANEL -- Couldn't set params"<<endl;
	}
}

void AliStorageAdministratorPanelSetStorageParams::onCloseButton(){onExit();}
void AliStorageAdministratorPanelSetStorageParams::CloseWindow(){onExit();}

void AliStorageAdministratorPanelSetStorageParams::onExit()
{
	cout<<"Quiting set storage params";
	if(fInstance){delete fInstance;fInstance=0;}
	cout<<" -- OK"<<endl;
}

Bool_t AliStorageAdministratorPanelSetStorageParams::ProcessMessage(Long_t msg, Long_t parm1, Long_t)
{
	switch (GET_MSG(msg))
	{
	case kC_COMMAND:
		switch (GET_SUBMSG(msg))
		{
		case kCM_BUTTON:
			switch(parm1)
			{
			case BUTTON_CLOSE:onCloseButton();break;
			case BUTTON_SET:onSetParamsButton();break;
			default:break;
			}
			break;
		default:break;
		}
		break;
	default:break;
	}

	return false;
}
