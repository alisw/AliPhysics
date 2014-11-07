#ifndef __AliDIMListenerThread__
#define __AliDIMListenerThread__

#include "AliDimIntNotifier.h"
/*
#ifdef ALI_DATE
#include <dic.hxx>
#endif
*/
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
    
    AliDimIntNotifier *fDimSORListener[5];
    AliDimIntNotifier *fDimEORListener[5];
};

#endif /* defined(__AliDIMListenerThread__) */
