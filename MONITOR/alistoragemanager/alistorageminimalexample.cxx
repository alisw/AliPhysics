#include "AliZMQManager.h"
#include <iostream>
#include <TFile.h>
#include <TThread.h>


using namespace std;

AliESDEvent *fCurrentEvent[2];
TTree *fCurrentTree[2];
TMutex fMutex;
int fEventInUse = 1;
int fWritingToEventIndex = 0;
bool fIsNewEventAvaliable = false;


int main()
{
    AliZMQManager *eventManager = AliZMQManager::GetInstance();
    
    fCurrentEvent[0]=0;
    fCurrentEvent[1]=0;
    fCurrentTree[0]=0;
    fCurrentTree[1]=0;
    AliESDEvent *tmpEvent = NULL;
    
    while(1)
    {
        eventManager->Get(tmpEvent,EVENTS_SERVER_SUB);
        
        if(tmpEvent)
        {
            if(tmpEvent->GetRunNumber()>=0)
            {
                fMutex.Lock();
                if(fEventInUse == 0){fWritingToEventIndex = 1;}
                else if(fEventInUse == 1){fWritingToEventIndex = 0;}
                cout<<"Received new event"<<endl;
                if(fCurrentEvent[fWritingToEventIndex])
                {
                    delete fCurrentEvent[fWritingToEventIndex];
                    fCurrentEvent[fWritingToEventIndex]=0;
                }
                fCurrentEvent[fWritingToEventIndex] = tmpEvent;
                fIsNewEventAvaliable = true;
                fMutex.UnLock();
            }
        }
    }
	return 0;
}
