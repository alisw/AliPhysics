#include "AliDIMListenerThread.h"

#include <iostream>

using namespace std;

AliDIMListenerThread::AliDIMListenerThread()
{
    InitDIMListeners();
    
#ifdef ALI_DATE
    /* DimCurrentInfo SORrunNumber("/LOGBOOK/SUBSCRIBE/DAQ_SOR_PHYSICS_1",-1);
    DimCurrentInfo EORrunNumber("/LOGBOOK/SUBSCRIBE/DAQ_EOR_PHYSICS_1",-1);

    if(SORrunNumber.getData() && EORrunNumber.getData())
    {
        cout<<"DIM Listener -- current SOR signal:"<<SORrunNumber.getInt()<<endl;
        cout<<"DIM Listener -- current EOR signal:"<<EORrunNumber.getInt()<<endl;
        
        if(SORrunNumber.getInt() != EORrunNumber.getInt()){StartOfRun(SORrunNumber.getInt());}
    }
    else{cout<<"DIM Listener -- no data received from dim server"<<endl;}*/
#endif
}

AliDIMListenerThread::~AliDIMListenerThread()
{/*
    for (int i = 0; i < 5; ++i){
        if(fDimSORListener[i]) delete fDimSORListener[i];
        if(fDimEORListener[i]) delete fDimEORListener[i];
        
        fDimSORListener[i] = 0;
        fDimEORListener[i] = 0;
	}*/
}

void AliDIMListenerThread::InitDIMListeners()
{/*
    for (int i = 0; i < 5; ++i)
    {
#ifdef ALI_DATE
             if (i == 0)
        {
            fDimSORListener[i] = new AliDimIntNotifier("/LOGBOOK/SUBSCRIBE/DAQ_SOR_PHYSICS");
            fDimEORListener[i] = new AliDimIntNotifier("/LOGBOOK/SUBSCRIBE/DAQ_EOR_PHYSICS");
        }
        else
        {
            fDimSORListener[i] = new AliDimIntNotifier(Form("/LOGBOOK/SUBSCRIBE/DAQ_SOR_PHYSICS_%d", i));
            fDimEORListener[i] = new AliDimIntNotifier(Form("/LOGBOOK/SUBSCRIBE/DAQ_EOR_PHYSICS_%d", i));
        }
        
        fDimSORListener[i]->Connect("DimMessage(int)", "AliDIMListenerThread", this, "StartOfRun(int)");
        fDimEORListener[i]->Connect("DimMessage(int)", "AliDIMListenerThread", this, "EndOfRun(int)");
#else
        fDimSORListener[i]=0x0;
        fDimEORListener[i]=0x0;
#endif
    }
    */
}

void AliDIMListenerThread::StartOfRun(int run)
{
    cout<<"DIM Listener -- SOR signal received for run:"<<run<<endl;
}

void AliDIMListenerThread::EndOfRun(int run)
{
    cout<<"DIM Listener -- EOR signal received for run:"<<run<<endl;
}
