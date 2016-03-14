#include "AliEveMyoListener.h"

#include <AliEveMultiView.h>

#include <TEveViewer.h>
#include <TGLViewer.h>
#include <TEveManager.h>

#include <iostream>

using namespace std;

AliEveMyoListener::AliEveMyoListener()
{
    
    myo::Hub hub("ch.cern.AliEVE");
    myo::Myo* myo = hub.waitForMyo(10000);
    
    if (!myo) {
        cout<<"Unable to find a Myo!"<<endl;
        return;
    }
    cout << "Connected to a Myo armband!" << endl;
    
    DataCollector collector;
    hub.addListener(&collector);
    
    while (1) {
        hub.run(1000/20);
        collector.print();
    }

    
//    TEveViewer *view = AliEveMultiView::Instance()->Get3DView();
//    TGLViewer *gl = view->GetGLViewer();
//    
//    while(1)
//    {
//        gl->CurrentCamera().Rotate(0.0,1.0,false,false);
//        gEve->FullRedraw3D();
//        
//        sleep(1);
//    }
}

AliEveMyoListener::~AliEveMyoListener()
{
    
}
