#include "AliStorageAdministratorPanelListEvents.h"
#include "AliOnlineReconstructionUtil.h"

#include <iostream>
#include <sstream>
#include <vector>
#include <set>

#include <TG3DLine.h>
#include <TGButton.h>
#include <TGFrame.h>

#include <TEnv.h>
#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>

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
    BUTTON_LOAD_EVENT,
    BUTTON_UPDATE_TC,
    TRIGGER_BOX
};

enum ENTRY{
    ENTRY_RUN_MIN=1,
    ENTRY_RUN_MAX,
    ENTRY_EVENT_MIN,
    ENTRY_EVENT_MAX,
    ENTRY_MULTIPLICITY_MIN,
    ENTRY_MULTIPLICITY_MAX
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
fGetListButton(0),
fMarkButton(0),
fLoadButton(0),
fListBox(0),
fEventsListVector(0),
fServerSocket(SERVER_COMMUNICATION_REQ),
fEventManager(0),
fCurrentEvent(0),
fOnlineMode(0)
{
    fEventManager = AliZMQManager::GetInstance();
    fEventManager->CreateSocket(fServerSocket);
    
    SetName("List");
    SetLayoutBroken(kTRUE);
    
    InitWindow();
    TriggerClassesFromCDB();
    fTriggerBox->Select(-1);
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
    
    fRunMinEntry = new TGNumberEntry(fRunGroupFrame,0,9,ENTRY_RUN_MIN,(TGNumberFormat::EStyle) 5);
    fRunGroupFrame->AddFrame(fRunMinEntry, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fRunMinEntry->MoveResize(38,22,80,22);
    
    TGLabel *fRunMinLabel = new TGLabel(fRunGroupFrame,"Min:");
    fRunMinLabel->SetTextJustify(36);
    fRunMinLabel->SetMargins(0,0,0,0);
    fRunMinLabel->SetWrapLength(-1);
    fRunGroupFrame->AddFrame(fRunMinLabel, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fRunMinLabel->MoveResize(8,22,30,18);
    
    TGLabel *fRunMaxLabel = new TGLabel(fRunGroupFrame,"Max:");
    fRunMaxLabel->SetTextJustify(36);
    fRunMaxLabel->SetMargins(0,0,0,0);
    fRunMaxLabel->SetWrapLength(-1);
    fRunGroupFrame->AddFrame(fRunMaxLabel, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fRunMaxLabel->MoveResize(120,22,30,24);
    
    fRunMaxEntry = new TGNumberEntry(fRunGroupFrame,999999,8,ENTRY_RUN_MAX,(TGNumberFormat::EStyle) 5);
    fRunGroupFrame->AddFrame(fRunMaxEntry, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fRunMaxEntry->MoveResize(150,22,72,22);
    
    fRunGroupFrame->SetLayoutManager(new TGVerticalLayout(fRunGroupFrame));
    AddFrame(fRunGroupFrame, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fRunGroupFrame->MoveResize(8,8,230,60);
    
    // "Event" group frame
    TGGroupFrame *fEventGroupFrame = new TGGroupFrame(this,"Event");
    fEventGroupFrame->SetLayoutBroken(kTRUE);
    
    fEventMinEntry = new TGNumberEntry(fEventGroupFrame,0,9,ENTRY_EVENT_MIN,(TGNumberFormat::EStyle) 5);
    fEventGroupFrame->AddFrame(fEventMinEntry, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fEventMinEntry->MoveResize(38,22,80,22);
    
    TGLabel *fEventMinLabel = new TGLabel(fEventGroupFrame,"Min:");
    fEventMinLabel->SetTextJustify(36);
    fEventMinLabel->SetMargins(0,0,0,0);
    fEventMinLabel->SetWrapLength(-1);
    fEventGroupFrame->AddFrame(fEventMinLabel, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fEventMinLabel->MoveResize(8,22,30,18);
    
    TGLabel *fEventMaxLabel = new TGLabel(fEventGroupFrame,"Max:");
    fEventMaxLabel->SetTextJustify(36);
    fEventMaxLabel->SetMargins(0,0,0,0);
    fEventMaxLabel->SetWrapLength(-1);
    fEventGroupFrame->AddFrame(fEventMaxLabel, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fEventMaxLabel->MoveResize(120,22,30,24);
    
    fEventMaxEntry = new TGNumberEntry(fEventGroupFrame,999999,8,ENTRY_EVENT_MAX,(TGNumberFormat::EStyle) 5);
    fEventGroupFrame->AddFrame(fEventMaxEntry, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fEventMaxEntry->MoveResize(150,22,72,22);
    
    fEventGroupFrame->SetLayoutManager(new TGVerticalLayout(fEventGroupFrame));
    AddFrame(fEventGroupFrame, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fEventGroupFrame->MoveResize(8,68,230,60);
    
    // "Multiplicity" group frame
    TGGroupFrame *fMultiplicityGroupFrame = new TGGroupFrame(this,"Multiplicity");
    fMultiplicityGroupFrame->SetLayoutBroken(kTRUE);
    
    fMultiplicityMinEntry = new TGNumberEntry(fMultiplicityGroupFrame,0,9,ENTRY_MULTIPLICITY_MIN,(TGNumberFormat::EStyle) 5);
    fMultiplicityGroupFrame->AddFrame(fMultiplicityMinEntry, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fMultiplicityMinEntry->MoveResize(38,22,80,22);
    
    TGLabel *fMultiplicityMinLabel = new TGLabel(fMultiplicityGroupFrame,"Min:");
    fMultiplicityMinLabel->SetTextJustify(36);
    fMultiplicityMinLabel->SetMargins(0,0,0,0);
    fMultiplicityMinLabel->SetWrapLength(-1);
    fMultiplicityGroupFrame->AddFrame(fMultiplicityMinLabel, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fMultiplicityMinLabel->MoveResize(8,22,30,18);
    
    TGLabel *fMultiplicityMaxLabel = new TGLabel(fMultiplicityGroupFrame,"Max:");
    fMultiplicityMaxLabel->SetTextJustify(36);
    fMultiplicityMaxLabel->SetMargins(0,0,0,0);
    fMultiplicityMaxLabel->SetWrapLength(-1);
    fMultiplicityGroupFrame->AddFrame(fMultiplicityMaxLabel, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fMultiplicityMaxLabel->MoveResize(120,22,30,24);
    
    fMultiplicityMaxEntry = new TGNumberEntry(fMultiplicityGroupFrame,999999,8,ENTRY_MULTIPLICITY_MAX,(TGNumberFormat::EStyle) 5);
    fMultiplicityGroupFrame->AddFrame(fMultiplicityMaxEntry, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fMultiplicityMaxEntry->MoveResize(150,22,72,22);
    
    fMultiplicityGroupFrame->SetLayoutManager(new TGVerticalLayout(fMultiplicityGroupFrame));
    AddFrame(fMultiplicityGroupFrame, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fMultiplicityGroupFrame->MoveResize(8,128,230,60);
    
    // "Beam" group frame
    TGGroupFrame *fBeamGroupFrame = new TGGroupFrame(this,"Beam");
    fBeamGroupFrame->SetLayoutBroken(kTRUE);
    
    fPPcheckbox = new TGCheckButton(fBeamGroupFrame,"p-p",BUTTON_CHECK_PP);
    fPPcheckbox->SetTextJustify(36);
    fPPcheckbox->SetMargins(0,0,0,0);
    fPPcheckbox->SetWrapLength(-1);
    fBeamGroupFrame->AddFrame(fPPcheckbox, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fPPcheckbox->MoveResize(8,22,40,19);
    
    fPbPbcheckbox = new TGCheckButton(fBeamGroupFrame,"Pb-Pb",BUTTON_CHECK_PBPB);
    fPbPbcheckbox->SetTextJustify(36);
    fPbPbcheckbox->SetMargins(0,0,0,0);
    fPbPbcheckbox->SetWrapLength(-1);
    fBeamGroupFrame->AddFrame(fPbPbcheckbox, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fPbPbcheckbox->MoveResize(60,22,55,19);
    
    fBeamGroupFrame->SetLayoutManager(new TGVerticalLayout(fBeamGroupFrame));
    AddFrame(fBeamGroupFrame, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fBeamGroupFrame->MoveResize(8,188,230,60);
    
    // "Storage" group frame
    TGGroupFrame *fStorageGroupFrame = new TGGroupFrame(this,"Storage");
    fStorageGroupFrame->SetLayoutBroken(kTRUE);
    
    fPermanentCheckbox = new TGCheckButton(fStorageGroupFrame,"Permanent",BUTTON_CHECK_PERM);
    fPermanentCheckbox->SetTextJustify(36);
    fPermanentCheckbox->SetMargins(0,0,0,0);
    fPermanentCheckbox->SetWrapLength(-1);
    fStorageGroupFrame->AddFrame(fPermanentCheckbox, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fPermanentCheckbox->MoveResize(8,22,80,19);
    
    fTemporaryCheckbox = new TGCheckButton(fStorageGroupFrame,"Temporary",BUTTON_CHECK_TEMP);
    fTemporaryCheckbox->SetTextJustify(36);
    fTemporaryCheckbox->SetMargins(0,0,0,0);
    fTemporaryCheckbox->SetWrapLength(-1);
    fStorageGroupFrame->AddFrame(fTemporaryCheckbox, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fTemporaryCheckbox->MoveResize(120,22,80,19);
    
    fStorageGroupFrame->SetLayoutManager(new TGVerticalLayout(fStorageGroupFrame));
    AddFrame(fStorageGroupFrame, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fStorageGroupFrame->MoveResize(8,248,230,60);
    
    // "Trigger" group frame
    TGGroupFrame *fTriggerGroupFrame = new TGGroupFrame(this,"Trigger");
    fTriggerGroupFrame->SetLayoutBroken(kTRUE);
    
    fTriggerBox = new TGComboBox(fTriggerGroupFrame,TRIGGER_BOX);
    fTriggerBox->AddEntry("No trigger selection",-1);
    fTriggerBox->Select(-1,kTRUE);
//    fTriggerBox->EnableTextInput(true);
    fTriggerGroupFrame->AddFrame(fTriggerBox, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fTriggerBox->MoveResize(8,22,200,19);
    
    fTriggerGroupFrame->SetLayoutManager(new TGVerticalLayout(fStorageGroupFrame));
    AddFrame(fTriggerGroupFrame, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fTriggerGroupFrame->MoveResize(8,310,230,60);
    
    
    // Status label group frame
    TGGroupFrame *fStatusGroupFrame = new TGGroupFrame(this,"Status");
    fStatusGroupFrame->SetLayoutBroken(kTRUE);
    
    fStatusLabel = new TGLabel(fStatusGroupFrame,"Status label");
    fStatusLabel->SetTextJustify(36);
    fStatusLabel->SetMargins(0,0,0,0);
    fStatusLabel->SetWrapLength(-1);
    fStatusGroupFrame->AddFrame(fStatusLabel, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fStatusLabel->MoveResize(8,22,80,19);
    
    fStatusGroupFrame->SetLayoutManager(new TGVerticalLayout(fStatusGroupFrame));
    AddFrame(fStatusGroupFrame, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fStatusGroupFrame->MoveResize(8,372,230,60);
    
    // buttons
    fMarkButton = new TGTextButton(this,"Mark event",BUTTON_MARK_EVENT);
    fMarkButton->SetTextJustify(36);
    fMarkButton->SetMargins(0,0,0,0);
    fMarkButton->SetWrapLength(-1);
    fMarkButton->Resize(100,24);
    AddFrame(fMarkButton, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fMarkButton->MoveResize(8,444,100,24);
    
    fLoadButton = new TGTextButton(this,"Load selected event",BUTTON_LOAD_EVENT);
    fLoadButton->SetTextJustify(36);
    fLoadButton->SetMargins(0,0,0,0);
    fLoadButton->SetWrapLength(-1);
    fLoadButton->Resize(130,24);
    AddFrame(fLoadButton, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLoadButton->MoveResize(108,444,130,24);
    
    fGetListButton = new TGTextButton(this,"List events",BUTTON_GET_LIST);
    fGetListButton->SetTextJustify(36);
    fGetListButton->SetMargins(0,0,0,0);
    fGetListButton->SetWrapLength(-1);
    fGetListButton->Resize(100,24);
    AddFrame(fGetListButton, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fGetListButton->MoveResize(8,468,100,24);
    
    fUpdateTriggersButton = new TGTextButton(this,"Update TCs",BUTTON_UPDATE_TC);
    fUpdateTriggersButton->SetTextJustify(36);
    fUpdateTriggersButton->SetMargins(0,0,0,0);
    fUpdateTriggersButton->SetWrapLength(-1);
    fUpdateTriggersButton->Resize(130,24);
    AddFrame(fUpdateTriggersButton, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fUpdateTriggersButton->MoveResize(108,468,130,24);
    
    // list box
    fListBox = new TGListBox(this);
    fListBox->SetName("fListBox");
    fListBox->AddEntry(new TGString("Run   Event   System   Mult   Marked"),0);
    AddFrame(fListBox, new TGLayoutHints(kLHintsLeft | kLHintsExpandY  | kLHintsTop,2,2,2,2));
    fListBox->MoveResize(8,510,230,436);
    
    MapSubwindows();
    MapWindow();
    Resize(252,809);
}

void AliStorageAdministratorPanelListEvents::TriggerClassesFromCDB()
{
    struct serverRequestStruct *requestMessage = new struct serverRequestStruct;
    requestMessage->messageType = REQUEST_GET_TRIGGER_LIST;
    
    if(!fEventManager->Send(requestMessage,fServerSocket))
    {
        cout<<"ADMIN PANEL -- couldn't send request for set of trigger classes"<<endl;
        return;
    }
    
    vector<string100> *tmpSet;
    if(!fEventManager->Get(tmpSet,fServerSocket))
    {
        cout<<"ADMIN PANEL -- problems getting server's response"<<endl;
        fEventManager->RecreateSocket(fServerSocket);
        return;
    }
    vector<string100> &receivedSet = *tmpSet;
    
    cout<<"ADMIN PANEL -- received set of size:"<<receivedSet.size()<<endl;
    
    vector<string100>::iterator it;
    int i=0;
    for (it = receivedSet.begin(); it != receivedSet.end(); ++it)
    {
        string100 cls = *it;
        fTriggerClasses.push_back(cls.data);
        fTriggerBox->AddEntry(cls.data,i++);
    }    
}


void AliStorageAdministratorPanelListEvents::onGetListButton()
{
    //prepare and send request message
    struct serverRequestStruct *requestMessage = new struct serverRequestStruct;
    
    //build query
    requestMessage->runNumber[0]=fRunMinEntry->GetIntNumber();
    requestMessage->runNumber[1]=fRunMaxEntry->GetIntNumber();
    requestMessage->eventNumber[0]=fEventMinEntry->GetIntNumber();
    requestMessage->eventNumber[1]=fEventMaxEntry->GetIntNumber();
    if(fTemporaryCheckbox->GetState()==1){requestMessage->marked[0]=0;}
    else{requestMessage->marked[0]=-1;}
    if(fPermanentCheckbox->GetState()==1){requestMessage->marked[1]=1;}
    else{requestMessage->marked[1]=-1;}
    requestMessage->multiplicity[0]=fMultiplicityMinEntry->GetIntNumber();
    requestMessage->multiplicity[1]=fMultiplicityMaxEntry->GetIntNumber();
    if(fPPcheckbox->GetState()==1){strcpy(requestMessage->system[0],"p-p");}
    else{strcpy(requestMessage->system[0],"");}
    if(fPbPbcheckbox->GetState()==1){strcpy(requestMessage->system[1],"A-A");}
    else{strcpy(requestMessage->system[1],"");}
    
    if(fTriggerBox->GetSelected()>=0)
    {
        const char* triggerClassName = fTriggerClasses[fTriggerBox->GetSelected()];
        cout<<"SELECTED TRIGGER CLASS:"<<triggerClassName<<endl;
        strcpy(requestMessage->triggerClass,triggerClassName);
    }
    else
    {
        strcpy(requestMessage->triggerClass,"No trigger selection");
    }
    requestMessage->messageType = REQUEST_LIST_EVENTS;
    
    if(!fEventManager->Send(requestMessage,fServerSocket))
    {
        cout<<"ADMIN PANEL -- couldn't send get list request"<<endl;
        return;
    }
    
    fListBox->RemoveAll();
    fListBox->AddEntry(new TGString("Run   Event   System   Mult   Marked"),0);
    
    vector<serverListStruct> *tmpVector;
    if(!fEventManager->Get(tmpVector,fServerSocket))
    {
        cout<<"ADMIN PANEL -- problems getting server's response"<<endl;
        fEventManager->RecreateSocket(fServerSocket);
        return;
    }
    vector<serverListStruct> &receivedList = *tmpVector;
    
    if(receivedList.size()==0){
        fStatusLabel->SetText("List is empty");
    }
    
    
    for(unsigned int i=0;i<receivedList.size();i++)
    {
        fListBox->AddEntry(new TGString(Form("%d   %d   %s   %d   %d   ",
                                             receivedList[i].runNumber,
                                             receivedList[i].eventNumber,
                                             receivedList[i].system,
                                             receivedList[i].multiplicity,
                                             receivedList[i].marked)),i+1);
    }
    fListBox->MoveResize(8,510,230,436);
    fEventsListVector = receivedList;
    
    gClient->HandleInput();
    gClient->NeedRedraw(fListBox, kTRUE);
    gClient->HandleInput();
    MapSubwindows();
    MapWindow();
    Layout();
    delete tmpVector;
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
    requestMessage->eventsRunNumber = runNumber;
    requestMessage->eventsEventNumber = eventNumber;
    requestMessage->messageType = REQUEST_MARK_EVENT;
    
    if(!fEventManager->Send(requestMessage,fServerSocket))
    {
        cout<<"ADMIN PANEL -- couldn't send mark request"<<endl;
        return;
    }
    
    bool response;
    if(!fEventManager->Get(&response,fServerSocket))
    {
        cout<<"ADMIN PANEL -- problems getting server's response"<<endl;
        fEventManager->RecreateSocket(fServerSocket);
        return;
    }
    
    if(response)
    {
        fStatusLabel->SetText("Event marked");
        cout<<"ADMIN PANEL -- Event marked succesfully"<<endl;
    }
    else
    {
        fStatusLabel->SetText("Error!");
        cout<<"ADMIN PANEL -- Could not matk event"<<endl;
    }
}

void AliStorageAdministratorPanelListEvents::onLoadButton()
{
    int selectedEventNumber = fListBox->GetSelected()-1;
    if(selectedEventNumber<0)return;
    int runNumber=fEventsListVector[selectedEventNumber].runNumber;
    int eventNumber=fEventsListVector[selectedEventNumber].eventNumber;
    
    
    struct serverRequestStruct *requestMessage = new struct serverRequestStruct;
    requestMessage->eventsRunNumber = runNumber;
    requestMessage->eventsEventNumber = eventNumber;
    requestMessage->messageType = REQUEST_GET_EVENT;
    
    if(!fEventManager->Send(requestMessage,fServerSocket))
    {
        cout<<"ADMIN PANEL -- couldn't send get event request"<<endl;
        return;
        
    }
    AliESDEvent *resultEvent;
    if(!fEventManager->Get(resultEvent,fServerSocket))
    {
        cout<<"ADMIN PANEL -- problems getting server's response"<<endl;
        fEventManager->RecreateSocket(fServerSocket);
        return;
    }
    
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
                case BUTTON_GET_LIST:onGetListButton();break;
                case BUTTON_MARK_EVENT:onMarkButton();break;
                case BUTTON_LOAD_EVENT:onLoadButton();break;
                case BUTTON_UPDATE_TC:TriggerClassesFromCDB();break;
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
    if(ison != fOnlineMode){fOnlineMode = ison;}
    else{return;}
    
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
        //    fPPcheckbox->SetOn();
        //    fPbPbcheckbox->SetOn();
        //    fTemporaryCheckbox->SetOn();
        //    fPermanentCheckbox->SetOn();
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
