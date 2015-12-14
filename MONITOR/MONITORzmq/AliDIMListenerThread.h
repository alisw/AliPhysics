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
    
    AliDIMListenerThread(const AliDIMListenerThread&);
    AliDIMListenerThread& operator=(const AliDIMListenerThread&);
    
    AliDimIntNotifier *fDimSORListener;//[5]; //now listening just in PHYSICS_1
    AliDimIntNotifier *fDimEORListener;//[5];

    std::string fOnlineReconstructionHostname;
    std::string fOnlineReconstructionUsername;
};

#endif
