#include "AliStorageClientThread.h"
    #include "AliZMQManager.h"

#include <signal.h>
#include <fstream>
#include <iostream>

using namespace std;

bool gClientQuit = false;
void GotSignalClient(int){gClientQuit = true;}

AliStorageClientThread::AliStorageClientThread() :
fDIMListenerThread(0),
fEventsCollectorThread(0),
fCommunicationThread(0),
fConnectionStatus(STATUS_WAITING),
fReceivingStatus(STATUS_WAITING),
fSavingStatus(STATUS_WAITING),
fCurrentStorageSize(0),
fMaximumStorageSize(0),
fStoragePath(""),
fNumberOfEventsInFile(0),
fStorageOccupationLevel(0),
fRemoveEventsPercentage(0)
{
    // make sure that when program is closed destructor will be called
    struct sigaction sa;
    memset(&sa,0,sizeof(sa));
    sa.sa_handler = GotSignalClient;
    sigfillset(&sa.sa_mask);
    sigaction(SIGINT,&sa,NULL);
   
    ifstream configFile (GetConfigFilePath());
    if (configFile.is_open())
    {
        string line;
        int from,to;
        while(configFile.good())
        {
            getline(configFile,line);
            from = line.find("\"")+1;
            to = line.find_last_of("\"");
            if(line.find("STORAGE_PATH=")==0){
                fStoragePath=line.substr(from,to-from);
            }
            else if(line.find("MAX_SIZE=")==0){
                fMaximumStorageSize=atoi(line.substr(from,to-from).c_str());
            }
            else if(line.find("MAX_OCCUPATION=")==0){
                fStorageOccupationLevel=atoi(line.substr(from,to-from).c_str());
            }
            else if(line.find("REMOVE_PERCENT=")==0){
                fRemoveEventsPercentage=atoi(line.substr(from,to-from).c_str());
            }
            else if(line.find("EVENTS_IN_FILE=")==0){
                fNumberOfEventsInFile=atoi(line.substr(from,to-from).c_str());
            }
        }
        if(configFile.eof()){configFile.clear();}
        configFile.close();
    }
    else{cout<<"CLIENT -- Unable to open config file"<<endl;}
    
    //create directory for storage if it doesn't exist
    gSystem->Exec(Form("mkdir -p %s",fStoragePath.c_str()));
    
    fDIMListenerThread = new AliDIMListenerThread();
    fEventsCollectorThread = new AliEventsCollectorThread(this);
    fCommunicationThread = new AliCommunicationThread(this);
}

AliStorageClientThread::~AliStorageClientThread()
{
    while(!gClientQuit){sleep(1);}
    cout<<"\n\nClosing\n\n"<<endl;
    AliZMQManager::GetInstance()->Close();
    if(fDIMListenerThread){delete fDIMListenerThread;}
    fEventsCollectorThread->Kill();
    fCommunicationThread->Kill();
}
