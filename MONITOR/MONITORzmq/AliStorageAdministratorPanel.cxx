#include "AliStorageAdministratorPanel.h"
#include "AliStorageAdministratorPanelMarkEvent.h"
#include "AliStorageAdministratorPanelListEvents.h"
#include "AliStorageAdministratorPanelSetStorageParams.h"
#include "AliESDEvent.h"

#include <fstream>
#include <iostream>

#include <TGMenu.h>

using namespace std;

ClassImp(AliStorageAdministratorPanel);

#define WINDOWS_WIDTH 500
#define WINDOWS_HEIGHT 350

enum TOOLBUTTON{
	TOOLBUTTON_START=1,
	TOOLBUTTON_STOP,
	TOOLBUTTON_PREFERENCES,
	TOOLBUTTON_EXIT 	
};

enum MENUBAR{
	MENUBAR_CLIENT_SET_PARAMS=1,
	MENUBAR_SERVER_LIST_EVENTS,
    MENUBAR_SERVER_MARK_EVENT,
    MENUBAR_SERVER_MARK_ALL_EVENTS,
	MENUBAR_SERVER_GET_EVENT,
	MENUBAR_SERVER_GET_NEXT_EVENT,
	MENUBAR_SERVER_GET_LAST_EVENT
};

enum FRAME{
	CLIENT_FRAME=1,
	SERVER_FRAME
};

AliStorageAdministratorPanel::AliStorageAdministratorPanel() :
	TGMainFrame(gClient->GetRoot(), 400, 400),
	fPanelQuited(false),
	fConnectionLabel(0),
	fDataLabel(0),
	fSavingLabel(0),
	fCurrentSizeLabel(0),
	fMaxSizeLabel(0),
	fMaxOccupationLabel(0),
	fRemoveEventsLabel(0),
	fEventsInChunkLabel(0),
	fMaxStorageSize(0),
	fMaxOccupation(0),
	fRemoveEvents(0),
	fEventsInChunk(0),
	fCommunicationThread(0),
	fCommunicationSocket(CLIENT_COMMUNICATION_REQ),
	fServerSocket(SERVER_COMMUNICATION_REQ),
	fEventManager(0)
{
	InitWindow();
    
	//create event manager
	fEventManager = AliZMQManager::GetInstance();
    fEventManager->CreateSocket(CLIENT_COMMUNICATION_REQ);
    fEventManager->CreateSocket(SERVER_COMMUNICATION_REQ);

	// start communication with client thread
	fCommunicationThread = new TThread("fCommunicationThread",
                                       Dispatch,(void*)this);
	fCommunicationThread->Run();
}

AliStorageAdministratorPanel::~AliStorageAdministratorPanel()
{
	cout<<"ADMIN -- AliStorageAdministratorPanel descructor called";
	cout<<" --- OK"<<endl;
}

void AliStorageAdministratorPanel::CheckStateHandle()//ask client about its state
{
	while(!fPanelQuited)
	{
		CheckClientState(REQUEST_CONNECTION);
		sleep(1);
		CheckClientState(REQUEST_RECEIVING);
		sleep(1);
		CheckClientState(REQUEST_SAVING);
		sleep(1);
		CheckClientState(REQUEST_CURRENT_SIZE);
		sleep(1);
		CheckClientState(REQUEST_GET_PARAMS);
		sleep(1);
	}
}

void AliStorageAdministratorPanel::CheckClientState(int option)
{
	struct clientRequestStruct *request = new struct clientRequestStruct;
	request->messageType = option;
	if(!fEventManager->Send(request,fCommunicationSocket))
	{
		SetLabel(fConnectionLabel,STATUS_DOWN);
		SetLabel(fDataLabel,STATUS_DOWN);
		SetLabel(fSavingLabel,STATUS_DOWN);
		cout<<"ADMIN -- CLIENT IS DOWN"<<endl;
		return;
	}

	long response = -1;
	struct clientRequestStruct *responseParams = NULL;
	
	if(option == REQUEST_GET_PARAMS)
	{
        fEventManager->Get(responseParams,fCommunicationSocket);
	}
	else
	{
        fEventManager->Get(&response,fCommunicationSocket);
	}
	switch(option)
	{
	case REQUEST_CONNECTION:
		SetLabel(fConnectionLabel,response);
		break;
	case REQUEST_RECEIVING:
		SetLabel(fDataLabel,response);
		break;
	case REQUEST_SAVING:
		SetLabel(fSavingLabel,response);
		break;
	case REQUEST_CURRENT_SIZE:
		SetLabelValue(fCurrentSizeLabel,response,1);
		break;
	case REQUEST_GET_PARAMS:
		fMaxStorageSize = responseParams->maxStorageSize;
		fMaxOccupation = responseParams->maxOccupation;
		fRemoveEvents = responseParams->removeEvents;
		fEventsInChunk = responseParams->eventsInChunk;
		
		SetLabelValue(fMaxSizeLabel,fMaxStorageSize,1);
		SetLabelValue(fMaxOccupationLabel,fMaxOccupation,2);
		SetLabelValue(fRemoveEventsLabel,fRemoveEvents,2);
		SetLabelValue(fEventsInChunkLabel,fEventsInChunk,3);
		break;
	default:break;	
	}
	if(request){delete request;}
}

void AliStorageAdministratorPanel::InitWindow()
{
	SetCleanup(kDeepCleanup);
	// add main manubar on top of the window
	SetupFixedMenuBar();
	AddFrame(new TGHorizontal3DLine(this),
		 new TGLayoutHints(kLHintsNormal | kLHintsExpandX));
	// add toolbar with pictured buttons
	/*SetupToolbar();
	AddFrame(new TGHorizontal3DLine(this),
		 new TGLayoutHints(kLHintsNormal | kLHintsExpandX));
	// add dockable toolbar
	SetupDockableMenuBar();
	AddFrame(new TGHorizontal3DLine(this),
		 new TGLayoutHints(kLHintsNormal | kLHintsExpandX));
	*/
	
	// add frame with two segments for client and server
	SetupThreadsFrame();
	SetWindowName("ALICE Storage Manager");
	MapSubwindows();
	Resize(WINDOWS_WIDTH,WINDOWS_HEIGHT);
	MapWindow();
}

void AliStorageAdministratorPanel::SetupThreadsFrame()
{
	// create frame for both client and server
	TGHorizontalFrame *threadsFrame = new TGHorizontalFrame(this);

	//create client frame and add components to it
	TGCompositeFrame *clientThreadFrame = new TGCompositeFrame(threadsFrame,CLIENT_FRAME);
	fConnectionLabel = new TGLabel(clientThreadFrame,"======Waiting======");
	fDataLabel = new TGLabel(clientThreadFrame,"======Waiting======");
	fSavingLabel = new TGLabel(clientThreadFrame,"======Waiting======");
	fCurrentSizeLabel = new TGLabel(clientThreadFrame,"=========Waiting=========");
	fMaxSizeLabel = new TGLabel(clientThreadFrame,"==========Waiting=========");
	fMaxOccupationLabel = new TGLabel(clientThreadFrame,"==========Waiting=========");
	fRemoveEventsLabel = new TGLabel(clientThreadFrame,"==========Waiting=========");
	fEventsInChunkLabel = new TGLabel(clientThreadFrame,"==========Waiting=========");
	
	clientThreadFrame->AddFrame(new TGLabel(clientThreadFrame,
						"ALICE Storage Client"),
				    new TGLayoutHints(kLHintsCenterX));

	//connection label
	clientThreadFrame->AddFrame(new TGLabel(clientThreadFrame,
						"Connection:"),
				    new TGLayoutHints(kLHintsCenterX));
	clientThreadFrame->AddFrame(fConnectionLabel,
				    new TGLayoutHints(kLHintsCenterX));

	//data receiving label
	clientThreadFrame->AddFrame(new TGLabel(clientThreadFrame,
						"Receiving data:"),
				    new TGLayoutHints(kLHintsCenterX));
	clientThreadFrame->AddFrame(fDataLabel,
				    new TGLayoutHints(kLHintsCenterX));

	//saving label
	clientThreadFrame->AddFrame(new TGLabel(clientThreadFrame,
						"Saving:"),
				    new TGLayoutHints(kLHintsCenterX));
	clientThreadFrame->AddFrame(fSavingLabel,
				    new TGLayoutHints(kLHintsCenterX));

	//current size label
	clientThreadFrame->AddFrame(new TGLabel(clientThreadFrame,
						"Current storage size:"),
				    new TGLayoutHints(kLHintsCenterX));
	clientThreadFrame->AddFrame(fCurrentSizeLabel,
				    new TGLayoutHints(kLHintsCenterX));

	//max size label
	clientThreadFrame->AddFrame(new TGLabel(clientThreadFrame,
						"Max storage size:"),
				    new TGLayoutHints(kLHintsCenterX));
	clientThreadFrame->AddFrame(fMaxSizeLabel,
				    new TGLayoutHints(kLHintsCenterX));

	//max occupation label
	clientThreadFrame->AddFrame(new TGLabel(clientThreadFrame,
						"Max occupation percent:"),
				    new TGLayoutHints(kLHintsCenterX));
	clientThreadFrame->AddFrame(fMaxOccupationLabel,
				    new TGLayoutHints(kLHintsCenterX));

	//remove events label
	clientThreadFrame->AddFrame(new TGLabel(clientThreadFrame,
						"Remove events percentage:"),
				    new TGLayoutHints(kLHintsCenterX));
	clientThreadFrame->AddFrame(fRemoveEventsLabel,
				    new TGLayoutHints(kLHintsCenterX));

	//events in chunk label
	clientThreadFrame->AddFrame(new TGLabel(clientThreadFrame,
						"Number of events in chunk:"),
				    new TGLayoutHints(kLHintsCenterX));
	clientThreadFrame->AddFrame(fEventsInChunkLabel,
				    new TGLayoutHints(kLHintsCenterX));

	//create server frame and add components to it
	TGCompositeFrame *serverThreadFrame = new TGCompositeFrame(threadsFrame,SERVER_FRAME);
	serverThreadFrame->AddFrame(new TGLabel(serverThreadFrame,
						 "ALICE Storage Server"),
				     new TGLayoutHints(kLHintsCenterX));

	//add client and server frames to threads frame
	threadsFrame->AddFrame(clientThreadFrame,
				new TGLayoutHints(kLHintsLeft | kLHintsExpandX));
	threadsFrame->AddFrame(new TGVertical3DLine(threadsFrame),
		 new TGLayoutHints(kLHintsNormal | kLHintsExpandY));
	threadsFrame->AddFrame(serverThreadFrame,
				new TGLayoutHints(kLHintsRight | kLHintsExpandX));

	//add threads frame to main window
	AddFrame(threadsFrame,new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
}

void AliStorageAdministratorPanel::SetupFixedMenuBar()
{
	//create popup menu for client
	TGPopupMenu *clientPopup = new TGPopupMenu(fClient->GetRoot());
	clientPopup->AddEntry("&Set storage params",MENUBAR_CLIENT_SET_PARAMS);
	clientPopup->Associate(this);

	//create popup menu for server
	TGPopupMenu *serverPopup = new TGPopupMenu(fClient->GetRoot());
	serverPopup->AddEntry("&Get Events List",MENUBAR_SERVER_LIST_EVENTS);
	serverPopup->AddEntry("&Mark Event",MENUBAR_SERVER_MARK_EVENT);
    serverPopup->AddEntry("&Mark ALL Events",MENUBAR_SERVER_MARK_ALL_EVENTS);
	serverPopup->AddEntry("&Get Event (test)",MENUBAR_SERVER_GET_EVENT);
	serverPopup->AddEntry("&Get Next Event (test)",MENUBAR_SERVER_GET_NEXT_EVENT);
	serverPopup->AddEntry("&Get Last Event (test)",MENUBAR_SERVER_GET_LAST_EVENT);
	serverPopup->Associate(this);

	//create menubar
	TGMenuBar *menuBar = new TGMenuBar(this,1,1,kHorizontalFrame);
	menuBar->AddPopup("&Client",clientPopup,
			  new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0));
	menuBar->AddPopup("&Server", serverPopup,
			  new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0));

	AddFrame(menuBar,new TGLayoutHints(kLHintsTop | kLHintsExpandX));
}

//handle GUI actions
Bool_t AliStorageAdministratorPanel::ProcessMessage(Long_t msg, Long_t parm1, Long_t)
{
	switch (GET_MSG(msg))
	{

	case kC_COMMAND:
		switch (GET_SUBMSG(msg))
		{
		case kCM_BUTTON://when button is pressed
			break;
		case kCM_MENUSELECT://when mouse is over menu entry
			break;
		case kCM_MENU://when menu entry was clicked
			switch (parm1)
			{
			case MENUBAR_CLIENT_SET_PARAMS:onClientSetParams();break;
			case MENUBAR_SERVER_LIST_EVENTS:onServerListEvents();break;
			case MENUBAR_SERVER_MARK_EVENT:onServerMarkEvent();break;
            case MENUBAR_SERVER_MARK_ALL_EVENTS:onServerMarkAllEvents();break;
			case MENUBAR_SERVER_GET_EVENT:onServerGetEvent();break;
			case MENUBAR_SERVER_GET_NEXT_EVENT:onServerGetNextEvent();break;
			case MENUBAR_SERVER_GET_LAST_EVENT:onServerGetLastEvent();break;
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

void AliStorageAdministratorPanel::onServerListEvents()
{
	AliStorageAdministratorPanelListEvents *listEventsWindow =
		AliStorageAdministratorPanelListEvents::GetInstance();
	if(listEventsWindow){cout<<"List events window created"<<endl;}
}

void AliStorageAdministratorPanel::onServerMarkEvent()
{
	AliStorageAdministratorPanelMarkEvent *markEventWindow =
		AliStorageAdministratorPanelMarkEvent::GetInstance();
	if(markEventWindow){cout<<"Mark event window created"<<endl;}
}

void AliStorageAdministratorPanel::onServerMarkAllEvents()
{
    // get list of all events:
    struct serverRequestStruct *requestMessage = new struct serverRequestStruct;
    requestMessage->messageType = REQUEST_LIST_EVENTS;
//    requestMessage->list = list;

//    struct listRequestStruct list;
    
    requestMessage->runNumber[0]=0;
    requestMessage->runNumber[1]=999999;
    requestMessage->eventNumber[0]=0;
    requestMessage->eventNumber[1]=999999;
    requestMessage->marked[0]=1;
    requestMessage->marked[1]=0;
    requestMessage->multiplicity[0]=1;
    requestMessage->multiplicity[1]=999999;
    strcpy(requestMessage->system[0],"p-p");
    strcpy(requestMessage->system[1],"A-A");
    
    
    fEventManager->Send(requestMessage,SERVER_COMMUNICATION_REQ);
    vector<serverListStruct> *tmpVector;
    fEventManager->Get(tmpVector,SERVER_COMMUNICATION_REQ);
    vector<serverListStruct> &receivedList = *tmpVector;
    
    cout<<"ADMIN PANEL -- received list of marked events:"<<receivedList.size()<<endl;
    
    int failCounter=0;
    
    int i;
    for(i=0;i<receivedList.size();i++)
    {
        cout<<"marking:"<<i<<endl;
        struct serverRequestStruct *requestMessage = new struct serverRequestStruct;
        requestMessage->messageType = REQUEST_MARK_EVENT;
//        requestMessage->event = mark;

//        struct eventStruct mark;
        requestMessage->eventsRunNumber   = receivedList[i].runNumber;
        requestMessage->eventsEventNumber = receivedList[i].eventNumber;
     
        
        fEventManager->Send(requestMessage,SERVER_COMMUNICATION_REQ);
        bool response;
        fEventManager->Get(&response,SERVER_COMMUNICATION_REQ);
        delete requestMessage;
        if(!response){failCounter++;}
    }
    cout<<"Tried to mark: "<<i<<"events"<<endl;
    cout<<"Could not mark: "<<failCounter<<" events"<<endl;
}

void AliStorageAdministratorPanel::onClientSetParams()
{
	AliStorageAdministratorPanelSetStorageParams *setParamsWindow =
		AliStorageAdministratorPanelSetStorageParams::GetInstance();

	setParamsWindow->Setup(fCommunicationSocket,
			       fMaxStorageSize/1000000,
			       fMaxOccupation,
			       fRemoveEvents,
			       fEventsInChunk);
}

void AliStorageAdministratorPanel::onServerGetEvent()
{
	//this method is just for tests
	int runNumber=197669;
	int eventNumber=168;
	
	struct serverRequestStruct *requestMessage = new struct serverRequestStruct;
//	struct eventStruct mark;
	requestMessage->eventsRunNumber = runNumber;
	requestMessage->eventsEventNumber = eventNumber;
	requestMessage->messageType = REQUEST_GET_EVENT;
//	requestMessage->event = mark;

	fEventManager->Send(requestMessage,fServerSocket);
    AliESDEvent *resultEvent = NULL;
    fEventManager->Get(resultEvent,fServerSocket);
	
	if(resultEvent)
	{
		cout<<"ADMIN -- received event. Run no:"<<resultEvent->GetRunNumber()<<"\tEvent no in file:"<<resultEvent->GetEventNumberInFile()<<endl;
	}
	else
	{
		cout<<"ADMIN -- received no event"<<endl;
	}
}


void AliStorageAdministratorPanel::onServerGetNextEvent()
{
	//this method is just for tests
	int runNumber=197669;
	int eventNumber=33;
	
	struct serverRequestStruct *requestMessage = new struct serverRequestStruct;
//	struct eventStruct mark;
	requestMessage->eventsRunNumber = runNumber;
	requestMessage->eventsEventNumber = eventNumber;
	requestMessage->messageType = REQUEST_GET_NEXT_EVENT;
//	requestMessage->event = mark;

	fEventManager->Send(requestMessage,fServerSocket);
    AliESDEvent* resultEvent = NULL;
    fEventManager->Get(resultEvent,fServerSocket);
	if(resultEvent)
	{
		cout<<"ADMIN -- received event. Run no:"<<resultEvent->GetRunNumber()<<"\tEvent no in file:"<<resultEvent->GetEventNumberInFile()<<endl;
	}
	else
	{
		cout<<"ADMIN -- received no event"<<endl;
	}
}


void AliStorageAdministratorPanel::onServerGetLastEvent()
{	
	struct serverRequestStruct *requestMessage = new struct serverRequestStruct;
	requestMessage->messageType = REQUEST_GET_LAST_EVENT;

	fEventManager->Send(requestMessage,fServerSocket);
    AliESDEvent* resultEvent = NULL;
    fEventManager->Get(resultEvent,fServerSocket);

	if(resultEvent)
	{
		cout<<"ADMIN -- received event. Run no:"<<resultEvent->GetRunNumber()<<"\tEvent no in file:"<<resultEvent->GetEventNumberInFile()<<endl;
	}
	else
	{
		cout<<"ADMIN -- received no event"<<endl;
	}
}

void AliStorageAdministratorPanel::onExit()
{
	cout<<"ADMIN -- Quiting ALICE Storage Admin Panel"<<endl;
	fPanelQuited=true;
	gSystem->ExitLoop();
}

void AliStorageAdministratorPanel::CloseWindow(){onExit();}


//methods to change labels:
void AliStorageAdministratorPanel::SetLabel(TGLabel *label, int option)
{
	switch(option)
	{
	case STATUS_WAITING:
		label->SetText("Waiting");
		label->SetTextColor(kBlack);
		break;
	case STATUS_OK:
		label->SetText("OK");
		label->SetTextColor(kBlack);
		break;
	case STATUS_ERROR:
		label->SetText("ERROR");
		label->SetTextColor(kBlack);
		break;
	case STATUS_DOWN:
		label->SetText("CLIENT IS DOWN");
		label->SetTextColor(kRed);
		break;
	default:break;
	}
	gClient->HandleInput();
}

void AliStorageAdministratorPanel::SetLabelValue(TGLabel *label, long value,int option)
{
	// options:
	// 1 - MB
	// 2 - percentage
	// 3 - number

	switch(option)
	{
	case 1:
		label->SetText(Form("%lu B (%.2f MB)",value,(float)value/(1000.*1000.)));
		break;
	case 2:
		label->SetText(Form("%d %%",(int)value));
		break;
	case 3:
		label->SetText(Form("%d",(int)value));
		break;
	default:break;
	}
       	label->SetTextColor(kBlack);
	gClient->HandleInput();
}
//other methods that can be useful later

/*
void AliStorageAdministratorPanel::SetupToolbar()
{
	TGToolBar* toolBar = new TGToolBar(this);
	toolBar->AddButton(this,new TGPictureButton(toolBar,
						     Form("%s/MONITOR/icons/start.png",
						     gSystem->Getenv("ALICE_ROOT")),
						     TOOLBUTTON_START ) );
	toolBar->AddButton(this, new TGPictureButton(toolBar,
						      Form("%s/MONITOR/icons/stop.png",
						      gSystem->Getenv("ALICE_ROOT")),
						      TOOLBUTTON_STOP) ); 
	
	AddFrame(toolBar, new TGLayoutHints(kLHintsNormal | kLHintsExpandX));
}

void AliStorageAdministratorPanel::SetupDockableMenuBar()
{
	TGDockableFrame *menuDock = new TGDockableFrame(this);
	AddFrame(menuDock,new TGLayoutHints(kLHintsExpandX,0,0,1,0));
	menuDock->SetWindowName("ALICE Storage Dockable  MenuBar");
	menuDock->EnableUndock(kTRUE);
	menuDock->EnableHide(kTRUE);

	TGPopupMenu *popup1 = new TGPopupMenu(fClient->GetRoot());
	popup1->AddEntry("&Check client state",MENUBAR_1_OPTION_1);
	popup1->AddSeparator();
	popup1->AddEntry("&BBB",MENUBAR_1_OPTION_2);
	popup1->AddEntry("&CCC",MENUBAR_1_OPTION_3);

	popup1->DisableEntry(MENUBAR_1_OPTION_2);
	popup1->Associate(this);

	
	TGMenuBar *menuBar = new TGMenuBar(menuDock,1,1,kHorizontalFrame);
	menuBar->AddPopup("&File",
			  popup1,
			  new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0));

	menuDock->AddFrame(menuBar,
			   new TGLayoutHints(kLHintsTop | kLHintsExpandX));
			   }*/
