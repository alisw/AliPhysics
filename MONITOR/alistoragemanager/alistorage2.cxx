#include "AliZMQManager.h"
#include "AliStorageTypes.h"
#include "AliESDEvent.h"

#include <TFile.h>

#include <iostream>

using namespace std;

int main(int argc, char **argv)
{
    AliZMQManager *eventManager = AliZMQManager::GetInstance();
    eventManager->CreateSocket(EVENTS_SERVER_SUB);
    
    AliESDEvent *event = NULL;
    TTree *tree = new TTree("esdTree","esdTree");
    
    int receiveStatus = false;
    
    TFile *outfile = new TFile();
    
    while (1)
    {
        cout<<"alistorage2 -- waiting for event..."<<endl;
        receiveStatus = eventManager->Get(event,EVENTS_SERVER_SUB);
        
        if (receiveStatus == 0){ // timeout
            continue;
        }
        else if(receiveStatus == -1){ // error, socket closed
            break;
        }
        else if(event && receiveStatus)
        {
            cout<<"alistorage2 -- received event"<<endl;
         
            event->WriteToTree(tree);
            tree->Fill();
            event->Reset();
            outfile->Open("/Users/Jerus/Desktop/AliESDs.root","UPDATE");
            tree->Write();
            outfile->Close();
//            delete event;event=0;
        }
        else
        {
            cout<<"alistorage2 -- ERROR -- NO DATA!"<<endl;
        }
    }

    
    return 0;
}
