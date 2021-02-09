#ifndef ALIANALYSISTASKEMCALOUTLIERSGEN_H
#define ALIANALYSISTASKEMCALOUTLIERSGEN_H

#include "AliAnalysisTaskEmcalJet.h"
class TH1;

namespace PWGJE {

namespace EMCALJetTasks {

class AliAnalysisTaskEmcalOutliersGen : public AliAnalysisTaskEmcalJet {
public:
    AliAnalysisTaskEmcalOutliersGen();
    AliAnalysisTaskEmcalOutliersGen(const char *name);
    virtual ~AliAnalysisTaskEmcalOutliersGen() {}

    static AliAnalysisTaskEmcalOutliersGen *AddTaskEmcalOutliersGen(const char *name);


protected:
    virtual bool IsEventSelected() { return true; }
    virtual void UserCreateOutputObjects();
    virtual bool Run();

private:
    TH1         *fHistJetPt;            //

    ClassDef(AliAnalysisTaskEmcalOutliersGen, 1);
};
}
}
#endif
