/*
 * AliFemtoEventCuts.cxx
 *
 *  Created on: Nov 22, 2017
 *      Author: gu74req
 */

#include "AliLightNEventCuts.h"
ClassImp(AliLightNEventCuts)
AliLightNEventCuts::AliLightNEventCuts()
:fHist(0)
,fCutMinContrib(false)
,fMinContrib(2)
,fCutZVtx(false)
,fzVtxLow(0)
,fzVtxUp(0)
,fPileUpRejection(false)
,fUseMVPileUpRej(false)
,fCleanEvtMult(false)
,fUseSPDMult(false)
,fUseV0AMult(false)
,fUseV0CMult(false)
,fUseRef08Mult(false)
,fUseAliEvtCuts(false)
{
}

AliLightNEventCuts::~AliLightNEventCuts() {
    if (fHist) {
        delete fHist;
    }
}

bool AliLightNEventCuts::isSelected(AliLightNEvent *evt) {
    bool pass=true;
    fHist->FillEvtCounter(0);
    if (fUseAliEvtCuts) {
        if (!evt->PassAliEvtSelection()) {
            pass=false;
        } else {
            fHist->FillEvtCounter(1);
        }
    } else {
        //This needs to be always the case!
        if (!(evt->GetMagneticField()||evt->GetHasVertex())) {
            pass=false;
        } else {
            fHist->FillEvtCounter(2);
        }
        if (pass&&fCutMinContrib) {
            if (evt->GetNumberOfContributers()<fMinContrib) {
                pass=false;
            } else {
                fHist->FillEvtCounter(3);
            }
        }
        if (pass&&fCutZVtx) {
            if (evt->GetZVertex()<=fzVtxLow||evt->GetZVertex()>=fzVtxUp) {
                pass=false;
            } else {
                fHist->FillEvtCounter(4);
            }
        }
        if (pass&&fPileUpRejection) {
            if (evt->GetIsPileUp()) {
                pass=false;
            } else {
                fHist->FillEvtCounter(5);
            }
        }
        if (pass&&fCleanEvtMult) {
            if (pass&&fUseSPDMult) {
                if (!(evt->GetSPDMult()>0)) {
                    pass=false;
                } else {
                    fHist->FillEvtCounter(6);
                }
            }
            if (pass&&fUseV0AMult) {
                if (!(evt->GetV0AMult()>0)) {
                    pass=false;
                } else {
                    fHist->FillEvtCounter(7);
                    //Fill Hist
                }
            }
            if (pass&&fUseV0CMult) {
                if (!(evt->GetV0CMult()>0)) {
                    pass=false;
                } else {
                    fHist->FillEvtCounter(8);
                    //Fill Hist
                }
            }
            if (pass&&fUseRef08Mult) {
                if (!(evt->GetRefMult08()>0)) {
                    pass=false;
                } else {
                    fHist->FillEvtCounter(9);
                }
            }
        }
    }
    evt->SetSelectionStatus(pass);
    BookQA(evt);
    return pass;
}

void AliLightNEventCuts::InitQA() {
    fHist=new AliLightNEventHist();
    BookCuts();
}

AliLightNEventCuts* AliLightNEventCuts::StandardCutsRun1() {
    AliLightNEventCuts *evtCuts=new AliLightNEventCuts();
    evtCuts->SetCutMinContrib(1);
    evtCuts->SetZVtxPosition(-10.,10.);
    evtCuts->PileUpRejection(true);
    evtCuts->CleanUpMult(true,false,false,false);
    evtCuts->UseDontWorryEvtCuts(false);
    return evtCuts;
}

AliLightNEventCuts* AliLightNEventCuts::StandardCutsRun2() {
    AliLightNEventCuts *evtCuts=new AliLightNEventCuts();
    evtCuts->UseDontWorryEvtCuts(true);
    return evtCuts;
}
void AliLightNEventCuts::BookQA(AliLightNEvent *evt) {
    for (int i=0;i<2;++i) {
        if (i==0||(i==1&&evt->GetSelectionStatus())) {
            fHist->FillEvtNCont(i,evt->GetNumberOfContributers());
            fHist->FillEvtVtxX(i,evt->GetXVertex());
            fHist->FillEvtVtxY(i,evt->GetYVertex());
            fHist->FillEvtVtxZ(i,evt->GetZVertex());
            fHist->FillMultSPD(i,evt->GetSPDMult());
            fHist->FillMultV0A(i,evt->GetV0AMult());
            fHist->FillMultV0C(i,evt->GetV0CMult());
            fHist->FillMultRef08(i,evt->GetRefMult08());
        }
    }
}
void AliLightNEventCuts::BookCuts() {
    if (fCutMinContrib) {
        fHist->FillCuts(0,fMinContrib);
    } else {
        fHist->FillCuts(0,0);
    }
    if (fCutZVtx) {
        fHist->FillCuts(1,1);
        fHist->FillCuts(2,fzVtxLow);
        fHist->FillCuts(3,fzVtxUp);
    } else {
        fHist->FillCuts(1,0);
        fHist->FillCuts(2,0);
        fHist->FillCuts(3,0);
    }
    if (fPileUpRejection) {
        fHist->FillCuts(4,1);
    } else {
        fHist->FillCuts(4,0);
    }
    if (fUseMVPileUpRej) {
        fHist->FillCuts(5,1);
    } else {
        fHist->FillCuts(5,0);
    }
    if (fCleanEvtMult) {
        if (fUseSPDMult) {
            fHist->FillCuts(6,1);
        } else {
            fHist->FillCuts(6,0);
        }
        if (fUseV0AMult) {
            fHist->FillCuts(7,1);
        } else {
            fHist->FillCuts(7,0);
        }
        if (fUseV0CMult) {
            fHist->FillCuts(8,1);
        } else {
            fHist->FillCuts(8,0);
        }
        if (fUseRef08Mult) {
            fHist->FillCuts(9,1);
        } else {
            fHist->FillCuts(9,0);
        }
    } else {
        fHist->FillCuts(6,0);
        fHist->FillCuts(7,0);
        fHist->FillCuts(8,0);
        fHist->FillCuts(9,0);
    }
    if (fUseAliEvtCuts) {
        fHist->FillCuts(10,1);
    } else {
        fHist->FillCuts(10,0);
    }
}








