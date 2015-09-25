#include "AliDIMListenerThread.h"
#include "AliStorageTypes.h"

#include <iostream>
#include <sstream>
#include <fstream>

using namespace std;

AliDIMListenerThread::AliDIMListenerThread() : 
  fDimSORListener(0),
  fDimEORListener(0),
  fOnlineReconstructionHostname(""),
  fOnlineReconstructionUsername("")
{    
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
            if(line.find("EVENT_SERVER=")==0){fOnlineReconstructionHostname=line.substr(from,to-from);}
            else if(line.find("EVENT_SERVER_USER=")==0){fOnlineReconstructionUsername=line.substr(from,to-from);}
        }
        if(configFile.eof()){configFile.clear();}
        configFile.close();
    }
    else{cout<<"AliDIMListenerThread -- Unable to open config file"<<endl;}

    InitDIMListeners();

#ifdef ALI_DATE
    // check if there is a run in partition PHYSICS_1
    DimCurrentInfo SORrunNumber("/LOGBOOK/SUBSCRIBE/DAQ_SOR_PHYSICS_1",-1);
    DimCurrentInfo EORrunNumber("/LOGBOOK/SUBSCRIBE/DAQ_EOR_PHYSICS_1",-1);

    if(SORrunNumber.getData() && EORrunNumber.getData())
    {
        cout<<"DIM Listener -- current SOR signal:"<<SORrunNumber.getInt()<<endl;
        cout<<"DIM Listener -- current EOR signal:"<<EORrunNumber.getInt()<<endl;
        
        if(SORrunNumber.getInt() != EORrunNumber.getInt()){StartOfRun(SORrunNumber.getInt());}
    }
    else{cout<<"AliDIMListenerThread -- no data received from dim server"<<endl;}
#endif
}

AliDIMListenerThread::~AliDIMListenerThread()
{
  cout<<"AliDIMListenerThread -- destructor called...";
  if(fDimSORListener){delete fDimSORListener;fDimSORListener = 0;}
  if(fDimEORListener){delete fDimEORListener;fDimEORListener = 0;}

  // kill all running reconstructions (to be changed later)
//    gSystem->Exec(Form("ssh -n -f %s@%s \"killall alionlinereco\"",fOnlineReconstructionUsername.c_str(),fOnlineReconstructionHostname.c_str()));

  cout<<"OK"<<endl;
}

void AliDIMListenerThread::InitDIMListeners()
{

#ifdef ALI_DATE
  fDimSORListener = new AliDimIntNotifier("/LOGBOOK/SUBSCRIBE/DAQ_SOR_PHYSICS_1");
  fDimEORListener = new AliDimIntNotifier("/LOGBOOK/SUBSCRIBE/DAQ_EOR_PHYSICS_1");
  fDimSORListener->Connect("DimMessage(int)", "AliDIMListenerThread", this, "StartOfRun(int)");
  fDimEORListener->Connect("DimMessage(int)", "AliDIMListenerThread", this, "EndOfRun(int)");
#endif
}

void AliDIMListenerThread::StartOfRun(int run)
{
    cout<<"AliDIMListenerThread -- SOR signal received for run:"<<run<<endl;
   
    // Kill reconstruction and start new one
    gSystem->Exec(Form("ssh -n -f %s@%s \". ~/EventServerTesting/setEnv.sh;killall alionlinereco;alionlinereco %d\"",fOnlineReconstructionUsername.c_str(),fOnlineReconstructionHostname.c_str(),run));
}

void AliDIMListenerThread::EndOfRun(int run)
{
    cout<<"AliDIMListenerThread -- EOR signal received for run:"<<run<<endl;

    // kill all running reconstructions (to be changed later)
    gSystem->Exec(Form("ssh -n -f %s@%s \"killall alionlinereco\"",fOnlineReconstructionUsername.c_str(),fOnlineReconstructionHostname.c_str()));
}
