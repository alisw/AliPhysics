#include "AliStorageAdministratorPanelListEvents.h"

#include <iostream>
#include <sstream>
#include <vector>

#include <TG3DLine.h>
#include <TGButton.h>
#include <TGFrame.h>

using namespace std;

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
	BUTTON_MARK_EVENT,
    BUTTON_LOAD_EVENT
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
  TGMainFrame(gClient->GetRoot(),10,10,kMainFrame | kVerticalFrame),
	fStatusLabel(0),
	fRunMinEntry(0),
	fRunMaxEntry(0),
	fEventMinEntry(0),
	fEventMaxEntry(0),
	fMultiplicityMinEntry(0),
	fMultiplicityMaxEntry(0),
	fPPcheckbox(0),
	fPbPbcheckbox(0),
	fTemporaryCheckbox(0),
	fPermanentCheckbox(0),
	fListBox(0),
	fEventsListVector(0),
	fServerSocket(SERVER_COMMUNICATION_REQ),
	fEventManager(0)
{
	fEventManager = AliStorageEventManager::GetEventManagerInstance();
	fEventManager->CreateSocket(fServerSocket);

	SetName("List");
	SetLayoutBroken(kTRUE);

	InitWindow();
}

AliStorageAdministratorPanelListEvents::~AliStorageAdministratorPanelListEvents()
{
	cout<<"ADMIN PANEL -- List events descructor called";
	cout<<" --- OK"<<endl;
}

AliStorageAdministratorPanelListEvents* AliStorageAdministratorPanelListEvents::GetInstance()
{
	if(!fInstance){fInstance = new AliStorageAdministratorPanelListEvents();}
	return fInstance;
}


void AliStorageAdministratorPanelListEvents::SelectedEvent()
{
    Emit("SelectedEvent()");
}

void AliStorageAdministratorPanelListEvents::InitWindow()
{
   // "Run" group frame
   TGGroupFrame *fRunGroupFrame = new TGGroupFrame(this,"Run");
   fRunGroupFrame->SetLayoutBroken(kTRUE);
   
   fRunNumberSlider = new TGDoubleHSlider(fRunGroupFrame,312,kSlider1 | kScaleBoth,-1,kHorizontalFrame);
   fRunNumberSlider->SetRange(0,999999);
   fRunNumberSlider->SetPosition(0,999999);
   fRunGroupFrame->AddFrame(fRunNumberSlider, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fRunNumberSlider->MoveResize(10,20,312,24);
   fRunNumberSlider->Connect("PositionChanged()","AliStorageAdministratorPanelListEvents",this,"RunSliderPositionChanged()");

   fRunMinEntry = new TGNumberEntry(fRunGroupFrame, 
				    (Double_t) 0,
				    10,
				    TEXTENTRY_RUN_MIN,
				    (TGNumberFormat::EStyle) 0,
				    (TGNumberFormat::EAttribute) 1,
				    (TGNumberFormat::ELimit) 2,1,999999);
   fRunGroupFrame->AddFrame(fRunMinEntry, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fRunMinEntry->MoveResize(50,50,88,22);
   fRunMinEntry->Connect("ValueSet(Long_t)","AliStorageAdministratorPanelListEvents",this,"RunChanged()");
   
   fRunMaxEntry = new TGNumberEntry(fRunGroupFrame, 
				    (Double_t) 999999,
				    10,
				    TEXTENTRY_RUN_MAX,
				    (TGNumberFormat::EStyle) 0,
				    (TGNumberFormat::EAttribute) 1,
				    (TGNumberFormat::ELimit) 2,0,999999);
   fRunGroupFrame->AddFrame(fRunMaxEntry, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fRunMaxEntry->MoveResize(200,50,88,22);
   fRunMaxEntry->Connect("ValueSet(Long_t)","AliStorageAdministratorPanelListEvents",this,"RunChanged()");
   
   TGLabel *fRunNumberMinLabel = new TGLabel(fRunGroupFrame,"Min:");
   fRunNumberMinLabel->SetTextJustify(36);
   fRunNumberMinLabel->SetMargins(0,0,0,0);
   fRunNumberMinLabel->SetWrapLength(-1);
   fRunGroupFrame->AddFrame(fRunNumberMinLabel, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fRunNumberMinLabel->MoveResize(10,50,40,18);
   TGLabel *fRunNumberMaxLabel = new TGLabel(fRunGroupFrame,"Max:");
   fRunNumberMaxLabel->SetTextJustify(36);
   fRunNumberMaxLabel->SetMargins(0,0,0,0);
   fRunNumberMaxLabel->SetWrapLength(-1);
   fRunGroupFrame->AddFrame(fRunNumberMaxLabel, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fRunNumberMaxLabel->MoveResize(150,50,40,24);

   fRunGroupFrame->SetLayoutManager(new TGVerticalLayout(fRunGroupFrame));
   fRunGroupFrame->Resize(330,90);
   AddFrame(fRunGroupFrame, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fRunGroupFrame->MoveResize(8,8,330,90);

  // "Event" group frame
   TGGroupFrame *fEventGroupFrame = new TGGroupFrame(this,"Event");
   fEventGroupFrame->SetLayoutBroken(kTRUE);
   
   fEventSlider = new TGDoubleHSlider(fEventGroupFrame,312,kSlider1 | kScaleBoth,-1,kHorizontalFrame);
   fEventSlider->SetRange(0,999999);
   fEventSlider->SetPosition(0,999999);
   fEventGroupFrame->AddFrame(fEventSlider, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fEventSlider->MoveResize(10,20,312,24);
   fEventSlider->Connect("PositionChanged()","AliStorageAdministratorPanelListEvents",this,"EventSliderPositionChanged()");
   fEventMinEntry = new TGNumberEntry(fEventGroupFrame,0,10,TEXTENTRY_EVENT_MIN,(TGNumberFormat::EStyle) 0,(TGNumberFormat::EAttribute) 1,(TGNumberFormat::ELimit) 0,0,999999);
   fEventGroupFrame->AddFrame(fEventMinEntry, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fEventMinEntry->MoveResize(50,50,88,22);
   fEventMinEntry->Connect("ValueSet(Long_t)","AliStorageAdministratorPanelListEvents",this,"EventChanged()");
   
   fEventMaxEntry = new TGNumberEntry(fEventGroupFrame,999999,10,TEXTENTRY_EVENT_MAX,(TGNumberFormat::EStyle) 5,(TGNumberFormat::EAttribute) 1,(TGNumberFormat::ELimit) 0,0,999999);
   fEventGroupFrame->AddFrame(fEventMaxEntry, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fEventMaxEntry->MoveResize(200,50,88,22);
   fEventMaxEntry->Connect("ValueSet(Long_t)","AliStorageAdministratorPanelListEvents",this,"EventChanged()");
  
   TGLabel *fEventMinLabel = new TGLabel(fEventGroupFrame,"Min:");
   fEventMinLabel->SetTextJustify(36);
   fEventMinLabel->SetMargins(0,0,0,0);
   fEventMinLabel->SetWrapLength(-1);
   fEventGroupFrame->AddFrame(fEventMinLabel, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fEventMinLabel->MoveResize(10,50,40,18);
   
   TGLabel *fEventMaxLabel = new TGLabel(fEventGroupFrame,"Max:");
   fEventMaxLabel->SetTextJustify(36);
   fEventMaxLabel->SetMargins(0,0,0,0);
   fEventMaxLabel->SetWrapLength(-1);
   fEventGroupFrame->AddFrame(fEventMaxLabel, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fEventMaxLabel->MoveResize(150,50,40,24);

   fEventGroupFrame->SetLayoutManager(new TGVerticalLayout(fEventGroupFrame));
   fEventGroupFrame->Resize(330,100);
   AddFrame(fEventGroupFrame, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fEventGroupFrame->MoveResize(8,102,330,100);

  // "Multiplicity" group frame
   TGGroupFrame *fMultiplicityGroupFrame = new TGGroupFrame(this,"Multiplicity");
   fMultiplicityGroupFrame->SetLayoutBroken(kTRUE);
  
   fMultiplicitySlider = new TGDoubleHSlider(fMultiplicityGroupFrame,312,kSlider1 | kScaleBoth,-1,kHorizontalFrame);
   fMultiplicitySlider->SetRange(0,99999);
   fMultiplicitySlider->SetPosition(0,99999);
   fMultiplicityGroupFrame->AddFrame(fMultiplicitySlider, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fMultiplicitySlider->MoveResize(10,20,312,24);
   fMultiplicitySlider->Connect("PositionChanged()","AliStorageAdministratorPanelListEvents",this,"MultiplicitySliderPositionChanged()");
   
   fMultiplicityMinEntry = new TGNumberEntry(fMultiplicityGroupFrame,0,10,TEXTENTRY_MULTIPLICITY_MIN,(TGNumberFormat::EStyle) 0,(TGNumberFormat::EAttribute) 1,(TGNumberFormat::ELimit) 0,0,999999);
   fMultiplicityGroupFrame->AddFrame(fMultiplicityMinEntry, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fMultiplicityMinEntry->MoveResize(50,50,88,22);
   fMultiplicityMinEntry->Connect("ValueSet(Long_t)","AliStorageAdministratorPanelListEvents",this,"MultiplicityChanged()");

   fMultiplicityMaxEntry = new TGNumberEntry(fMultiplicityGroupFrame,99999,10,TEXTENTRY_MULTIPLICITY_MAX,(TGNumberFormat::EStyle) 0,(TGNumberFormat::EAttribute) 1,(TGNumberFormat::ELimit) 0,0,999999);
   fMultiplicityGroupFrame->AddFrame(fMultiplicityMaxEntry, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fMultiplicityMaxEntry->MoveResize(200,50,88,22);
   fMultiplicityMaxEntry->Connect("ValueSet(Long_t)","AliStorageAdministratorPanelListEvents",this,"MultiplicityChanged()");

   TGLabel *fMultiplicityMinLabel = new TGLabel(fMultiplicityGroupFrame,"Min:");
   fMultiplicityMinLabel->SetTextJustify(36);
   fMultiplicityMinLabel->SetMargins(0,0,0,0);
   fMultiplicityMinLabel->SetWrapLength(-1);
   fMultiplicityGroupFrame->AddFrame(fMultiplicityMinLabel, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fMultiplicityMinLabel->MoveResize(10,50,40,18);
   
   TGLabel *fMultiplicityMaxLabel = new TGLabel(fMultiplicityGroupFrame,"Max:");
   fMultiplicityMaxLabel->SetTextJustify(36);
   fMultiplicityMaxLabel->SetMargins(0,0,0,0);
   fMultiplicityMaxLabel->SetWrapLength(-1);
   fMultiplicityGroupFrame->AddFrame(fMultiplicityMaxLabel, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fMultiplicityMaxLabel->MoveResize(150,50,40,24);

   fMultiplicityGroupFrame->SetLayoutManager(new TGVerticalLayout(fMultiplicityGroupFrame));
   fMultiplicityGroupFrame->Resize(330,100);
   AddFrame(fMultiplicityGroupFrame, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fMultiplicityGroupFrame->MoveResize(8,196,330,100);


   // "Beam" group frame
   TGGroupFrame *fBeamGroupFrame = new TGGroupFrame(this,"Beam");
   fBeamGroupFrame->SetLayoutBroken(kTRUE);
   
   fPPcheckbox = new TGCheckButton(fBeamGroupFrame,"p-p",BUTTON_CHECK_PP);
   fPPcheckbox->SetTextJustify(36);
   fPPcheckbox->SetMargins(0,0,0,0);
   fPPcheckbox->SetWrapLength(-1);
   fBeamGroupFrame->AddFrame(fPPcheckbox, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fPPcheckbox->MoveResize(10,20,60,24);
   
   fPbPbcheckbox = new TGCheckButton(fBeamGroupFrame,"Pb-Pb",BUTTON_CHECK_PBPB);
   fPbPbcheckbox->SetTextJustify(36);
   fPbPbcheckbox->SetMargins(0,0,0,0);
   fPbPbcheckbox->SetWrapLength(-1);
   fBeamGroupFrame->AddFrame(fPbPbcheckbox, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fPbPbcheckbox->MoveResize(90,20,60,24);

   fBeamGroupFrame->SetLayoutManager(new TGVerticalLayout(fBeamGroupFrame));
   fBeamGroupFrame->Resize(330,60);
   AddFrame(fBeamGroupFrame, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fBeamGroupFrame->MoveResize(8,290,330,60);
  
   fGetListButton = new TGTextButton(this,"Get events' list",BUTTON_GET_LIST);
   fGetListButton->SetTextJustify(36);
   fGetListButton->SetMargins(0,0,0,0);
   fGetListButton->SetWrapLength(-1);
   fGetListButton->Resize(130,24);
   AddFrame(fGetListButton, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fGetListButton->MoveResize(8,408,130,24);
   
   fMarkButton = new TGTextButton(this,"Mark selected event",BUTTON_MARK_EVENT);
   fMarkButton->SetTextJustify(36);
   fMarkButton->SetMargins(0,0,0,0);
   fMarkButton->SetWrapLength(-1);
   fMarkButton->Resize(130,24);
   AddFrame(fMarkButton, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fMarkButton->MoveResize(8,432,130,24);
   
   fLoadButton = new TGTextButton(this,"Load selected event",BUTTON_LOAD_EVENT);
   fLoadButton->SetTextJustify(36);
   fLoadButton->SetMargins(0,0,0,0);
   fLoadButton->SetWrapLength(-1);
   fLoadButton->Resize(130,24);
   AddFrame(fLoadButton, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fLoadButton->MoveResize(8,456,130,24);

   // list box
   
   fListBox = new TGListBox(this);
   fListBox->SetName("fListBox");
   fListBox->AddEntry("Entry 1",0);
   fListBox->AddEntry("Entry 2",1);
   fListBox->AddEntry("Entry 3",2);
   fListBox->AddEntry("Entry 4",3);
   fListBox->AddEntry("Entry 5",4);
   fListBox->AddEntry("Entry 6",5);
   fListBox->AddEntry("Entry 7",6);
   fListBox->Resize(330,132);
   AddFrame(fListBox, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fListBox->MoveResize(8,488,330,132);

 
   // "Type" group frame
   TGGroupFrame *fTypeGroupFrame = new TGGroupFrame(this,"Type");
   fTypeGroupFrame->SetLayoutBroken(kTRUE);
   
   fTemporaryCheckbox = new TGCheckButton(fTypeGroupFrame,"Temporary",BUTTON_CHECK_TEMP);
   fTemporaryCheckbox->SetTextJustify(36);
   fTemporaryCheckbox->SetMargins(0,0,0,0);
   fTemporaryCheckbox->SetWrapLength(-1);
   fTypeGroupFrame->AddFrame(fTemporaryCheckbox, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fTemporaryCheckbox->MoveResize(10,20,90,24);
   
   fPermanentCheckbox = new TGCheckButton(fTypeGroupFrame,"Permanent",BUTTON_CHECK_PERM);
   fPermanentCheckbox->SetTextJustify(36);
   fPermanentCheckbox->SetMargins(0,0,0,0);
   fPermanentCheckbox->SetWrapLength(-1);
   fTypeGroupFrame->AddFrame(fPermanentCheckbox, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fPermanentCheckbox->MoveResize(110,20,90,24);

   fTypeGroupFrame->SetLayoutManager(new TGVerticalLayout(fTypeGroupFrame));
   fTypeGroupFrame->Resize(330,60);
   AddFrame(fTypeGroupFrame, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fTypeGroupFrame->MoveResize(8,344,330,60);

   SetMWMHints(kMWMDecorAll, kMWMFuncAll, kMWMInputModeless);
   MapSubwindows();
   Resize(GetDefaultSize());
   MapWindow();
   Resize(341,945);
}




void AliStorageAdministratorPanelListEvents::onGetListButton()
{
//prepare and send request message
	struct serverRequestStruct *requestMessage = new struct serverRequestStruct;
	struct listRequestStruct list; 

	//get listing parameters from somwhere
	list.runNumber[0]=fRunMinEntry->GetIntNumber();
	list.runNumber[1]=fRunMaxEntry->GetIntNumber();
	list.eventNumber[0]=fEventMinEntry->GetIntNumber();
	list.eventNumber[1]=fEventMaxEntry->GetIntNumber();
	if(fTemporaryCheckbox->GetState()==1)
	{
		list.marked[0]=0;
	}
	else
	{
		list.marked[0]=-1;
	}
	if(fPermanentCheckbox->GetState()==1)
	{
		list.marked[1]=1;
	}
	else
	{
		list.marked[1]=-1;
	}
	list.multiplicity[0]=fMultiplicityMinEntry->GetIntNumber();
	list.multiplicity[1]=fMultiplicityMaxEntry->GetIntNumber();
	if(fPPcheckbox->GetState()==1)
	{
		strcpy(list.system[0],"p-p");
	}
	else
	{
		strcpy(list.system[0],"");
	}
	if(fPbPbcheckbox->GetState()==1)
	{
		strcpy(list.system[1],"A-A");
	}
	else
	{
		strcpy(list.system[1],"");
	}
	
	requestMessage->messageType = REQUEST_LIST_EVENTS;
	requestMessage->list = list;

	fEventManager->Send(requestMessage,fServerSocket);

	fListBox->RemoveAll();
	fListBox->AddEntry(new TGString("Run   Event   System   Mult   Marked"),0);
	
	vector<serverListStruct> receivedList = fEventManager->GetServerListVector(fServerSocket);
    
	for(unsigned int i=0;i<receivedList.size();i++)
	{
		fListBox->InsertEntry(Form("%d   %d   %s   %d   %d   ",
					      receivedList[i].runNumber,
					      receivedList[i].eventNumber,
					      receivedList[i].system,
					      receivedList[i].multiplicity,
					      receivedList[i].marked),i+1,i);
	}

	fEventsListVector = receivedList;
	
	gClient->HandleInput();
	gClient->NeedRedraw(fListBox, kTRUE);
	gClient->HandleInput();
	MapSubwindows();
	MapWindow();
	Layout();
	delete requestMessage;
}


void AliStorageAdministratorPanelListEvents::onMarkButton()
{
	int runNumber;
	int eventNumber;

	//get run and event number from selected row
	int selectedEventNumber = fListBox->GetSelected()-1;
	
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

	bool response = fEventManager->GetBool(fServerSocket);
		
	if(response)
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

void AliStorageAdministratorPanelListEvents::onLoadButton()
{
    int selectedEventNumber = fListBox->GetSelected()-1;
    int runNumber=fEventsListVector[selectedEventNumber].runNumber;
    int eventNumber=fEventsListVector[selectedEventNumber].eventNumber;
    
    
    struct serverRequestStruct *requestMessage = new struct serverRequestStruct;
    struct eventStruct eventToLoad;
    eventToLoad.runNumber = runNumber;
    eventToLoad.eventNumber = eventNumber;
    requestMessage->messageType = REQUEST_GET_EVENT;
    requestMessage->event = eventToLoad;
    
    fEventManager->Send(requestMessage,fServerSocket);
    AliESDEvent *resultEvent = fEventManager->GetEvent(fServerSocket);
    
    if(resultEvent)
    {
        cout<<"ADMIN -- received event"<<endl;
        fCurrentEvent = resultEvent;
    }
    else
    {
        cout<<"ADMIN -- received no event"<<endl;
    }

    SelectedEvent();
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
            case BUTTON_LOAD_EVENT:onLoadButton();break;
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
	
void AliStorageAdministratorPanelListEvents::SetOfflineMode(Bool_t ison)
{

  if (ison) {
    fPPcheckbox->SetDisabledAndSelected(ison);
    fPbPbcheckbox->SetDisabledAndSelected(ison);
    fTemporaryCheckbox->SetDisabledAndSelected(ison);
    fPermanentCheckbox->SetDisabledAndSelected(ison);
  }
  else {
    fPPcheckbox->SetEnabled(!ison);
    fPbPbcheckbox->SetEnabled(!ison);
    fTemporaryCheckbox->SetEnabled(!ison);
    fPermanentCheckbox->SetEnabled(!ison);
    fPPcheckbox->SetOn();
    fPbPbcheckbox->SetOn();
    fTemporaryCheckbox->SetOn();
    fPermanentCheckbox->SetOn();
  }

  fRunMinEntry->SetState(!ison);
  fRunMaxEntry->SetState(!ison);
  fEventMinEntry->SetState(!ison);
  fEventMaxEntry->SetState(!ison);
  fMultiplicityMinEntry->SetState(!ison);
  fMultiplicityMaxEntry->SetState(!ison);

  // fCloseButton->SetEnabled(!ison);
  fGetListButton->SetEnabled(!ison);
  fMarkButton->SetEnabled(!ison);
  fLoadButton->SetEnabled(!ison);
}


void AliStorageAdministratorPanelListEvents::RunSliderPositionChanged()
{
  fRunMinEntry->SetNumber(fRunNumberSlider->GetMinPosition());
  fRunMaxEntry->SetNumber(fRunNumberSlider->GetMaxPosition());
}

void AliStorageAdministratorPanelListEvents::RunChanged()
{
  fRunNumberSlider->SetPosition(fRunMinEntry->GetIntNumber(),fRunMaxEntry->GetIntNumber());
}

void AliStorageAdministratorPanelListEvents::EventSliderPositionChanged()
{
  fEventMinEntry->SetNumber(fEventSlider->GetMinPosition());
  fEventMaxEntry->SetNumber(fEventSlider->GetMaxPosition());
}

void AliStorageAdministratorPanelListEvents::EventChanged()
{
  fEventSlider->SetPosition(fEventMinEntry->GetIntNumber(),fEventMaxEntry->GetIntNumber());
}

void AliStorageAdministratorPanelListEvents::MultiplicitySliderPositionChanged()
{
  fMultiplicityMinEntry->SetNumber(fMultiplicitySlider->GetMinPosition());
  fMultiplicityMaxEntry->SetNumber(fMultiplicitySlider->GetMaxPosition());
}

void AliStorageAdministratorPanelListEvents::MultiplicityChanged()
{
  fMultiplicitySlider->SetPosition(fMultiplicityMinEntry->GetIntNumber(),fMultiplicityMaxEntry->GetIntNumber());
}
