#include "AliStorageAdministratorPanelListEvents.h"

#include <iostream>
#include <sstream>
#include <vector>

#include <TG3DLine.h>
#include <TGButton.h>
#include <TGFrame.h>

#include "zmq.hpp"

using namespace std;
using namespace zmq;

AliStorageAdministratorPanelListEvents *AliStorageAdministratorPanelListEvents::fInstance=0;

ClassImp(AliStorageAdministratorPanelListEvents);

#define WINDOWS_WIDTH 400
#define WINDOWS_HEIGHT 500

enum BUTTON{
	BUTTON_CLOSE=1,
	BUTTON_GET_LIST,
	BUTTON_CHECK_PP,
	BUTTON_CHECK_PBPB,
	BUTTON_CHECK_TEMP,
	BUTTON_CHECK_PERM,
	BUTTON_MARK_EVENT
};

enum TEXTENTRY{
	TEXTENTRY_RUN_MIN=1,
	TEXTENTRY_RUN_MAX,
	TEXTENTRY_EVENT_MIN,
	TEXTENTRY_EVENT_MAX,
	TEXTENTRY_MULTIPLICITY_MIN,
	TEXTENTRY_MULTIPLICITY_MAX
};

AliStorageAdministratorPanelListEvents::AliStorageAdministratorPanelListEvents() :
	TGMainFrame(gClient->GetRoot(),400,400),
	fStatusLabel(0),
	fRunNumberMinEntry(0),
	fRunNumberMaxEntry(0),
	fEventNumberMinEntry(0),
	fEventNumberMaxEntry(0),
	fMultiplicityMinEntry(0),
	fMultiplicityMaxEntry(0),
	fProtonProtonCheckButton(0),
	fLeadLeadCheckButton(0),
	fTempCheckButton(0),
	fPermCheckButton(0),
	fEventsList(0),
	fEventsListVector(0),
	fServerSocket(0),
	fEventManager(0)
{
	fEventManager = new AliStorageEventManager();
	InitWindow();
}

AliStorageAdministratorPanelListEvents::~AliStorageAdministratorPanelListEvents()
{
	cout<<"ADMIN PANEL -- List events descructor called";
//	DestroyWindow();
	cout<<" --- OK"<<endl;
}

AliStorageAdministratorPanelListEvents* AliStorageAdministratorPanelListEvents::GetInstance()
{
	if(!fInstance){fInstance = new AliStorageAdministratorPanelListEvents();}
	return fInstance;
}

void AliStorageAdministratorPanelListEvents::SetSocket(socket_t *socket)
{
	fServerSocket = socket;
}

void AliStorageAdministratorPanelListEvents::InitWindow()
{
	SetCleanup(kDeepCleanup);
	
	//min run number
	AddFrame(new TGLabel(this,"Minimum run number:"),new TGLayoutHints(kLHintsLeft));

	fRunNumberMinEntry = new TGNumberEntry(this,
					    0,
					    6,
					    TEXTENTRY_RUN_MIN,
					    TGNumberFormat::kNESInteger,
					    TGNumberFormat::kNEAPositive,
					    TGNumberFormat::kNELNoLimits);
	AddFrame(fRunNumberMinEntry,new TGLayoutHints(kLHintsLeft));

	//max run number
	AddFrame(new TGLabel(this,"Maximum run number:"),new TGLayoutHints(kLHintsLeft));

	fRunNumberMaxEntry = new TGNumberEntry(this,
					    999999,
					    6,
					    TEXTENTRY_RUN_MAX,
					    TGNumberFormat::kNESInteger,
					    TGNumberFormat::kNEAPositive,
					    TGNumberFormat::kNELNoLimits);
	AddFrame(fRunNumberMaxEntry,new TGLayoutHints(kLHintsLeft));

	//min event number
	AddFrame(new TGLabel(this,"Minimum event number:"),new TGLayoutHints(kLHintsLeft));

	fEventNumberMinEntry = new TGNumberEntry(this,
					      0,
					      6,
					      TEXTENTRY_EVENT_MIN,
					      TGNumberFormat::kNESInteger,
					      TGNumberFormat::kNEAPositive,
					      TGNumberFormat::kNELNoLimits);
	AddFrame(fEventNumberMinEntry,new TGLayoutHints(kLHintsLeft));	

	//max event number
	AddFrame(new TGLabel(this,"Maximum event number:"),new TGLayoutHints(kLHintsLeft));

	fEventNumberMaxEntry = new TGNumberEntry(this,
					      99999,
					      6,
					      TEXTENTRY_EVENT_MAX,
					      TGNumberFormat::kNESInteger,
					      TGNumberFormat::kNEAPositive,
					      TGNumberFormat::kNELNoLimits);
	AddFrame(fEventNumberMaxEntry,new TGLayoutHints(kLHintsLeft));	

	//min multiplicity
	AddFrame(new TGLabel(this,"Minimum multiplicity:"),new TGLayoutHints(kLHintsLeft));

	fMultiplicityMinEntry = new TGNumberEntry(this,
					      0,
					      6,
					      TEXTENTRY_MULTIPLICITY_MIN,
					      TGNumberFormat::kNESInteger,
					      TGNumberFormat::kNEAPositive,
					      TGNumberFormat::kNELNoLimits);
	AddFrame(fMultiplicityMinEntry,new TGLayoutHints(kLHintsLeft));	

	//max multiplicity
	AddFrame(new TGLabel(this,"Maximum multiplicity:"),new TGLayoutHints(kLHintsLeft));

	fMultiplicityMaxEntry = new TGNumberEntry(this,
					      9999,
					      6,
					      TEXTENTRY_MULTIPLICITY_MAX,
					      TGNumberFormat::kNESInteger,
					      TGNumberFormat::kNEAPositive,
					      TGNumberFormat::kNELNoLimits);
	AddFrame(fMultiplicityMaxEntry,new TGLayoutHints(kLHintsLeft));	

	//p-p check button
	AddFrame(new TGLabel(this,"System:"),new TGLayoutHints(kLHintsLeft));
	fProtonProtonCheckButton = new TGCheckButton(this,"p-p",BUTTON_CHECK_PP);
	fProtonProtonCheckButton->SetOn();
	AddFrame(fProtonProtonCheckButton,new TGLayoutHints(kLHintsLeft));

	//Pb-Pb check button
	fLeadLeadCheckButton = new TGCheckButton(this,"Pb-Pb",BUTTON_CHECK_PBPB);
	fLeadLeadCheckButton->SetOn();
	AddFrame(fLeadLeadCheckButton,new TGLayoutHints(kLHintsLeft));

	//temp check button
	AddFrame(new TGLabel(this,"Storage type:"),new TGLayoutHints(kLHintsLeft));
	fTempCheckButton = new TGCheckButton(this,"Temporary",BUTTON_CHECK_TEMP);
	fTempCheckButton->SetOn();
	AddFrame(fTempCheckButton,new TGLayoutHints(kLHintsLeft));

	//perm check button
	fPermCheckButton = new TGCheckButton(this,"Permanent",BUTTON_CHECK_PERM);
	fPermCheckButton->SetOn();
	AddFrame(fPermCheckButton,new TGLayoutHints(kLHintsLeft));
	
	// status label
	fStatusLabel = new TGLabel(this,"");
	AddFrame(fStatusLabel,new TGLayoutHints(kLHintsExpandX | kLHintsLeft));

	//buttons
	AddFrame(new TGTextButton(this,"Close",BUTTON_CLOSE),
		 new TGLayoutHints(kLHintsLeft));

	AddFrame(new TGTextButton(this,"Get event's list",BUTTON_GET_LIST),
		 new TGLayoutHints(kLHintsRight));

	AddFrame(new TGTextButton(this,"Mark selected event",BUTTON_MARK_EVENT),
		 new TGLayoutHints(kLHintsRight));


	//event's list
	fEventsList = new TGListBox(this,0);
	fEventsList->AddEntry(new TGString("Run   Event   System   Mult   Marked"),0);
	AddFrame(fEventsList,new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
	
	SetWindowName("List Events");
	MapSubwindows();
	Resize(WINDOWS_WIDTH,WINDOWS_HEIGHT);
	MapWindow();
}




void AliStorageAdministratorPanelListEvents::onGetListButton()
{
//prepare and send request message
	struct serverRequestStruct *requestMessage = new struct serverRequestStruct;
	struct listRequestStruct list; 

	//get listing parameters from somwhere
	list.runNumber[0]=fRunNumberMinEntry->GetIntNumber();
	list.runNumber[1]=fRunNumberMaxEntry->GetIntNumber();
	list.eventNumber[0]=fEventNumberMinEntry->GetIntNumber();
	list.eventNumber[1]=fEventNumberMaxEntry->GetIntNumber();
	if(fTempCheckButton->GetState()==1)
	{
		list.marked[0]=0;
	}
	else
	{
		list.marked[0]=-1;
	}
	if(fPermCheckButton->GetState()==1)
	{
		list.marked[1]=1;
	}
	else
	{
		list.marked[1]=-1;
	}
	list.multiplicity[0]=fMultiplicityMinEntry->GetIntNumber();
	list.multiplicity[1]=fMultiplicityMaxEntry->GetIntNumber();
	if(fProtonProtonCheckButton->GetState()==1)
	{
		strcpy(list.system[0],"p-p");
	}
	else
	{
		strcpy(list.system[0],"");
	}
	if(fLeadLeadCheckButton->GetState()==1)
	{
		strcpy(list.system[1],"Pb-Pb");
	}
	else
	{
		strcpy(list.system[1],"");
	}
	
	requestMessage->messageType = REQUEST_LIST_EVENTS;
	requestMessage->list = list;

	fEventManager->Send(requestMessage,fServerSocket);

	fEventsList->RemoveAll();
	fEventsList->AddEntry(new TGString("Run   Event   System   Mult   Marked"),0);
	
	vector<serverListStruct> receivedList = fEventManager->GetServerListVector(fServerSocket);
	
	cout<<"PANEL:"<<receivedList[0].runNumber<<endl;	
	cout<<"VECTOR SIZE:"<<receivedList.size()<<endl;
	
//do something with list of maching events
	cout<<"Received list of perm events"<<endl;

	for(unsigned int i=0;i<receivedList.size();i++)
	{
		fEventsList->InsertEntry(Form("%d   %d   %s   %d   %d   ",
					      receivedList[i].runNumber,
					      receivedList[i].eventNumber,
					      receivedList[i].system,
					      receivedList[i].multiplicity,
					      receivedList[i].marked),i+1,i);

		cout<<receivedList[i].runNumber<<receivedList[i].eventNumber<<endl;
	
	}

	fEventsListVector = receivedList;
	
	gClient->HandleInput();
	gClient->NeedRedraw(fEventsList, kTRUE);
	gClient->HandleInput();
	MapSubwindows();
	MapWindow();
	Layout();
}


void AliStorageAdministratorPanelListEvents::onMarkButton()
{
	int runNumber;
	int eventNumber;

	//get run and event number from selected row
	int selectedEventNumber = fEventsList->GetSelected()-1;

	cout<<"SELECTED:"<<selectedEventNumber<<endl;
	
	if(selectedEventNumber<0)return;
	
	runNumber=fEventsListVector[selectedEventNumber].runNumber;
	eventNumber=fEventsListVector[selectedEventNumber].eventNumber;
	
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
	//if(response)delete response;
	
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

void AliStorageAdministratorPanelListEvents::onCloseButton(){onExit();}
void AliStorageAdministratorPanelListEvents::CloseWindow(){onExit();}

void AliStorageAdministratorPanelListEvents::onExit()
{
	cout<<"Quiting list events";
	if(fInstance){delete fInstance;fInstance=0;}
	cout<<" -- OK"<<endl;
}

Bool_t AliStorageAdministratorPanelListEvents::ProcessMessage(Long_t msg, Long_t parm1, Long_t)
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
			case BUTTON_GET_LIST:onGetListButton();break;
			case BUTTON_MARK_EVENT:onMarkButton();break;
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
