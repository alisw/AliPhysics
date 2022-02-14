#ifndef __ALIANALYSISTASKEMCALTRIGGERBACKGROUND_H__
#define __ALIANALYSISTASKEMCALTRIGGERBACKGROUND_H__

#include <AliAnalysisTaskEmcal.h>
class THistManager;

namespace PWGJE {
    
namespace EMCALJetTasks { 

class AliAnalysisTaskEmcalTriggerBackground : public AliAnalysisTaskEmcal {
public:
    AliAnalysisTaskEmcalTriggerBackground();
    AliAnalysisTaskEmcalTriggerBackground(const char *name);
    virtual ~AliAnalysisTaskEmcalTriggerBackground();

    static AliAnalysisTaskEmcalTriggerBackground *AddTaskEmcalTriggerBackground(const char *nametag);

protected:
    virtual void UserCreateOutputObjects();
    virtual bool Run();

private:

    THistManager                *fHistos;               //!<! histograms

    ClassDef(AliAnalysisTaskEmcalTriggerBackground, 1);

};

}

}

#endif
