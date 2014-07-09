#include "AliStorageAdministratorPanelMarkEvent.h"
#include "AliStorageTypes.h"

#include <iostream>
#include "zmq.hpp"

#include <TGFrame.h>
#include <TGButton.h>

using namespace std;
using namespace zmq;

AliStorageAdministratorPanelMarkEvent *AliStorageAdministratorPanelMarkEvent::fInstance=0;

ClassImp(AliStorageAdministratorPanelMarkEvent);

#define WINDOWS_WIDTH 200
#define WINDOWS_HEIGHT 200

enum BUTTON{
	BUTTON_CLOSE=1,
	BUTTON_MARK 	
};

enum TEXTENTRY{
	TEXTENTRY_RUN=1,
	TEXTENTRY_EVENT
};

AliStorageAdministratorPanelMarkEvent::AliStorageAdministratorPanelMarkEvent() :
	TGMainFrame(gClient->GetRoot(), 400, 400),
	fStatusLabel(0),
	fRunNumberEntry(0),
	fEventNumberEntry(0),
	fServerSocket(0),
	fEventManager(0)
{
	fEventManager = new AliStorageEventManager();
	InitWindow();
}

AliStorageAdministratorPanelMarkEvent::~AliStorageAdministratorPanelMarkEvent()
{
	cout<<"ADMIN PANEL -- Mark Window descructor called";
//	if(fStatusLabel)delete fStatusLabel;
//	if(fRunNumberEntry)delete fRunNumberEntry;
//	if(fEventNumberEntry)delete fEventNumberEntry;
	//if(fServerSocket)delete fServerSocket;
//	DestroyWindow();
	cout<<" --- OK"<<endl;
}

AliStorageAdministratorPanelMarkEvent* AliStorageAdministratorPanelMarkEvent::GetInstance()
{
	if(!fInstance){fInstance = new AliStorageAdministratorPanelMarkEvent();}
	return fInstance;
}

void AliStorageAdministratorPanelMarkEvent::SetSocket(socket_t *socket)
{
	fServerSocket = socket;
}

void AliStorageAdministratorPanelMarkEvent::InitWindow()
{
	SetCleanup(kDeepCleanup);
	
	AddFrame(new TGLabel(this,"Run number:"),new TGLayoutHints(kLHintsLeft));

	fRunNumberEntry = new TGNumberEntry(this,
					    0,
					    6,
					    TEXTENTRY_RUN,
					    TGNumberFormat::kNESInteger,
					    TGNumberFormat::kNEAPositive,
					    TGNumberFormat::kNELNoLimits);
	AddFrame(fRunNumberEntry,new TGLayoutHints(kLHintsLeft));

	AddFrame(new TGLabel(this,"Event number:"),new TGLayoutHints(kLHintsLeft));

	fEventNumberEntry = new TGNumberEntry(this,
					      0,
					      6,
					      TEXTENTRY_EVENT,
					      TGNumberFormat::kNESInteger,
					      TGNumberFormat::kNEAPositive,
					      TGNumberFormat::kNELNoLimits);
	AddFrame(fEventNumberEntry,new TGLayoutHints(kLHintsLeft));	

	fStatusLabel = new TGLabel(this,"");
	AddFrame(fStatusLabel,new TGLayoutHints(kLHintsExpandX | kLHintsLeft));

	AddFrame(new TGTextButton(this,"Close",BUTTON_CLOSE),
		 new TGLayoutHints(kLHintsLeft));

	AddFrame(new TGTextButton(this,"Mark event",BUTTON_MARK),
		 new TGLayoutHints(kLHintsRight));

	
	SetWindowName("Mark Event");
	MapSubwindows();
	Resize(WINDOWS_WIDTH,WINDOWS_HEIGHT);
	MapWindow();
}




void AliStorageAdministratorPanelMarkEvent::onMarkButton()
{
	int runNumber;
	int eventNumber;

	//get run and event number from TGNumberEntries
	runNumber=fRunNumberEntry->GetIntNumber();
	eventNumber=fEventNumberEntry->GetIntNumber();
	
	struct serverRequestStruct *requestMessage = new struct serverRequestStruct;
	struct eventStruct mark;
	mark.runNumber = runNumber;
	mark.eventNumber = eventNumber;
	requestMessage->messageType = REQUEST_MARK_EVENT;
	requestMessage->event = mark;

	fEventManager->Send(requestMessage,fServerSocket);

	message_t *response = new message_t();
	fServerSocket->recv(response);
	char *result = (char*)response->data();
	
	if(!strcmp("true",result))
	{
		fStatusLabel->SetText("Event marked");
		cout<<"ADMIN PANEL -- Event marked succesfully"<<endl;
	}
	else
	{
		fStatusLabel->SetText("Couldn't mark this event");
		cout<<"ADMIN PANEL -- Could not matk event"<<endl;
	}
}

void AliStorageAdministratorPanelMarkEvent::onCloseButton(){onExit();}
void AliStorageAdministratorPanelMarkEvent::CloseWindow(){onExit();}

void AliStorageAdministratorPanelMarkEvent::onExit()
{
	cout<<"Quiting mark event";
	if(fInstance){delete fInstance;fInstance=0;}
	cout<<" -- OK"<<endl;
}

Bool_t AliStorageAdministratorPanelMarkEvent::ProcessMessage(Long_t msg, Long_t parm1, Long_t)
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
			case BUTTON_MARK:onMarkButton();break;
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
