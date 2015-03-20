#ifndef __AliDIMListenerThread__
#define __AliDIMListenerThread__

#include "AliDimIntNotifier.h"

class AliDimIntNotifier;

class AliDIMListenerThread
{
public:
    AliDIMListenerThread();
    ~AliDIMListenerThread();
    
    void StartOfRun(int run);
    void EndOfRun(int run);
    
private:
    void InitDIMListeners();
    
    AliDimIntNotifier *fDimSORListener; //now listening just in PHYSICS_1
    AliDimIntNotifier *fDimEORListener;

    std::string fOnlineReconstructionHostname;
    std::string fOnlineReconstructionUsername;
};

#endif
