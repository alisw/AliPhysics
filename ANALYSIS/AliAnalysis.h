#ifndef ALIANALYSIS_H
#define ALIANALYSIS_H
//________________________________
///////////////////////////////////////////////////////////
//
// class AliAnalysis
//
// Base class for analysis
//
//
// Piotr.Skowronski@cern.ch
//
///////////////////////////////////////////////////////////

#include <TTask.h>

class AliESD;
class AliStack;
 
class AliAnalysis: public TTask
{
  public: 
    AliAnalysis();
    AliAnalysis(const char* name,const char* title);
    virtual ~AliAnalysis();
    
    virtual Int_t Init() = 0;
    virtual Int_t ProcessEvent(AliESD* esd, AliStack* stack = 0x0) = 0;
    virtual Int_t Finish() = 0;
    
    
    static Int_t GetDebug() {return fgkDebug;}
    static void  SetDebug(Int_t level) {fgkDebug = level;}
    
  protected:
    
  private:
    static Int_t fgkDebug;//! debug level
    ClassDef(AliAnalysis,1)
};

#endif
