/*
 * AliFemtoEventCuts.cxx
 *
 *  Created on: Nov 22, 2017
 *      Author: gu74req
 */

#include "AliFemtoDreamEventCuts.h"
ClassImp(AliFemtoDreamEventCuts)
AliFemtoDreamEventCuts::AliFemtoDreamEventCuts()
    : fHist(0),
      fMinimalBooking(false),
      fCutMinContrib(false),
      fMinContrib(2),
      fCutZVtx(false),
      fzVtxLow(0),
      fzVtxUp(0),
      fPileUpRejection(false),
      fUseMVPileUpRej(false),
      fCleanEvtMult(false),
      fUseSPDMult(false),
      fUseV0AMult(false),
      fUseV0CMult(false),
      fUseRef08Mult(false),
      fUseMultPercentileCut(false),
      fMultPercentileMax(999.f),
      fUseAliEvtCuts(false),
      fCentVsMultPlots(false),
      fDoSpherCuts(false),
      fLowPtSpherCalc(0.5),
      fSpherCutsLow(0.f),
      fSpherCutsUp(1.f),
      fDoSpheroCuts(false),
      fSpheroCutsLow(0.f),
      fSpheroCutsUp(1.f){
}

AliFemtoDreamEventCuts::AliFemtoDreamEventCuts(
    const AliFemtoDreamEventCuts& cuts)
    : fHist(cuts.fHist),
      fMinimalBooking(cuts.fMinimalBooking),
      fCutMinContrib(cuts.fCutMinContrib),
      fMinContrib(cuts.fMinContrib),
      fCutZVtx(cuts.fCutZVtx),
      fzVtxLow(cuts.fzVtxLow),
      fzVtxUp(cuts.fzVtxUp),
      fPileUpRejection(cuts.fPileUpRejection),
      fUseMVPileUpRej(cuts.fUseMVPileUpRej),
      fCleanEvtMult(cuts.fCleanEvtMult),
      fUseSPDMult(cuts.fUseSPDMult),
      fUseV0AMult(cuts.fUseV0AMult),
      fUseV0CMult(cuts.fUseV0CMult),
      fUseRef08Mult(cuts.fUseRef08Mult),
      fUseMultPercentileCut(cuts.fUseMultPercentileCut),
      fMultPercentileMax(cuts.fMultPercentileMax),
      fUseAliEvtCuts(cuts.fUseAliEvtCuts),
      fCentVsMultPlots(cuts.fCentVsMultPlots),
      fDoSpherCuts(cuts.fDoSpherCuts),
      fLowPtSpherCalc(cuts.fLowPtSpherCalc),
      fSpherCutsLow(cuts.fSpherCutsLow),
      fSpherCutsUp(cuts.fSpherCutsUp),
      fDoSpheroCuts(cuts.fDoSpheroCuts),
      fSpheroCutsLow(cuts.fSpheroCutsLow),
      fSpheroCutsUp(cuts.fSpheroCutsUp){
}

AliFemtoDreamEventCuts& AliFemtoDreamEventCuts::operator=(
    const AliFemtoDreamEventCuts& cuts) {
  if (this == &cuts) {
    return *this;
  }
  this->fHist = cuts.fHist;
  this->fMinimalBooking = cuts.fMinimalBooking;
  this->fCutMinContrib = cuts.fCutMinContrib;
  this->fMinContrib = cuts.fMinContrib;
  this->fCutZVtx = cuts.fCutZVtx;
  this->fzVtxLow = cuts.fzVtxLow;
  this->fzVtxUp = cuts.fzVtxUp;
  this->fPileUpRejection = cuts.fPileUpRejection;
  this->fUseMVPileUpRej = cuts.fUseMVPileUpRej;
  this->fCleanEvtMult = cuts.fCleanEvtMult;
  this->fUseSPDMult = cuts.fUseSPDMult;
  this->fUseV0AMult = cuts.fUseV0AMult;
  this->fUseV0CMult = cuts.fUseV0CMult;
  this->fUseRef08Mult = cuts.fUseRef08Mult;
  this->fUseMultPercentileCut = cuts.fUseMultPercentileCut;
  this->fMultPercentileMax = cuts.fMultPercentileMax;
  this->fUseAliEvtCuts = cuts.fUseAliEvtCuts;
  this->fCentVsMultPlots = cuts.fCentVsMultPlots;
  this->fDoSpherCuts = cuts.fDoSpherCuts;
  this->fLowPtSpherCalc = cuts.fLowPtSpherCalc;
  this->fSpherCutsLow = cuts.fSpherCutsLow;
  this->fSpherCutsUp = cuts.fSpherCutsUp;
  this->fDoSpheroCuts = cuts.fDoSpheroCuts;
  this->fSpheroCutsLow = cuts.fSpheroCutsLow;
  this->fSpheroCutsUp = cuts.fSpheroCutsUp;
  return *this;
}

AliFemtoDreamEventCuts::~AliFemtoDreamEventCuts() {

}

bool AliFemtoDreamEventCuts::isSelected(AliFemtoDreamEvent *evt) {
  bool pass = true;
  if (!fMinimalBooking)
    fHist->FillEvtCounter(0);
  if (fUseAliEvtCuts) {
    if (!evt->PassAliEvtSelection()) {
      pass = false;
    } else {
      if (!fMinimalBooking)
        fHist->FillEvtCounter(1);
    }
  } else {
    //This needs to be always the case!
    if (!(evt->GetMagneticField() || evt->GetHasVertex())) {
      pass = false;
    } else {
      if (!fMinimalBooking)
        fHist->FillEvtCounter(2);
    }
    if (pass && fCutMinContrib) {
      if (evt->GetNumberOfContributers() < fMinContrib) {
        pass = false;
      } else {
        if (!fMinimalBooking)
          fHist->FillEvtCounter(3);
      }
    }
    if (pass && fCutZVtx) {
      if (evt->GetZVertex() <= fzVtxLow || evt->GetZVertex() >= fzVtxUp) {
        pass = false;
      } else {
        if (!fMinimalBooking)
          fHist->FillEvtCounter(4);
      }
    }
    if (pass && fPileUpRejection) {
      if (evt->GetIsPileUp()) {
        pass = false;
      } else {
        if (!fMinimalBooking)
          fHist->FillEvtCounter(5);
      }
    }
  }
  //we need to make sure that our evt mult estimator is >0 else we wont
  //be able to find a bin to put the event in.
  if (pass && fCleanEvtMult) {
    if (pass && fUseSPDMult) {
      if (!(evt->GetSPDMult() > 0)) {
        pass = false;
      } else {
        if (!fMinimalBooking)
          fHist->FillEvtCounter(6);
      }
    }
    if (pass && fUseV0AMult) {
      if (!(evt->GetV0AMult() > 0)) {
        pass = false;
      } else {
        if (!fMinimalBooking)
          fHist->FillEvtCounter(7);
        //Fill Hist
      }
    }
    if (pass && fUseV0CMult) {
      if (!(evt->GetV0CMult() > 0)) {
        pass = false;
      } else {
        if (!fMinimalBooking)
          fHist->FillEvtCounter(8);
        //Fill Hist
      }
    }
    if (pass && fUseRef08Mult) {
      if (!(evt->GetRefMult08() > 0)) {
        pass = false;
      } else {
        if (!fMinimalBooking)
          fHist->FillEvtCounter(9);
      }
    }
  }
  if (pass && fDoSpherCuts) {
    if (evt->GetSpher() <= fSpherCutsLow || evt->GetSpher() >= fSpherCutsUp) {
      pass = false;
    } else {
      if (!fMinimalBooking)
        fHist->FillEvtCounter(10);
    }
  }

  if (pass & fUseMultPercentileCut) {
    if (evt->GetV0MCentrality() > fMultPercentileMax) {
      pass = false;
    } else {
      if (!fMinimalBooking) {
        fHist->FillEvtCounter(11);
      }
    }
  }
  if (pass && fDoSpheroCuts) {
    if (evt->GetSphero() <= fSpheroCutsLow || evt->GetSphero() >= fSpheroCutsUp) {
      pass = false;
    } else {
      if (!fMinimalBooking)
        fHist->FillEvtCounter(12);
    }
  }


  evt->SetSelectionStatus(pass);
  BookQA(evt);
  return pass;
}

void AliFemtoDreamEventCuts::InitQA() {
  if (fMinimalBooking) {
    fHist = 0;
  } else {
    fHist = new AliFemtoDreamEventHist(fCentVsMultPlots);
    BookCuts();
  }
}

AliFemtoDreamEventCuts* AliFemtoDreamEventCuts::StandardCutsRun1() {
  AliFemtoDreamEventCuts *evtCuts = new AliFemtoDreamEventCuts();
  evtCuts->SetCutMinContrib(2);
  evtCuts->SetZVtxPosition(-10., 10.);
  evtCuts->PileUpRejection(true);
  evtCuts->CleanUpMult(true, false, false, false);
  evtCuts->UseDontWorryEvtCuts(false);
  return evtCuts;
}

AliFemtoDreamEventCuts* AliFemtoDreamEventCuts::StandardCutsRun2() {
  AliFemtoDreamEventCuts *evtCuts = new AliFemtoDreamEventCuts();
  evtCuts->UseDontWorryEvtCuts(true);
  return evtCuts;
}
void AliFemtoDreamEventCuts::BookQA(AliFemtoDreamEvent *evt) {
  if (!fMinimalBooking) {
    for (int i = 0; i < 2; ++i) {
      if (i == 0 || (i == 1 && evt->GetSelectionStatus())) {
        fHist->FillEvtNCont(i, evt->GetNumberOfContributers());
        fHist->FillEvtVtxX(i, evt->GetXVertex());
        fHist->FillEvtVtxY(i, evt->GetYVertex());
        fHist->FillEvtVtxZ(i, evt->GetZVertex());
        fHist->FillMultSPD(i, evt->GetSPDMult());
        fHist->FillMultV0A(i, evt->GetV0AMult());
        fHist->FillMultV0C(i, evt->GetV0CMult());
        fHist->FillMultV0M(i, evt->GetV0MMult());
        fHist->FillMultRef08(i, evt->GetRefMult08());
        fHist->FillMultPercentV0(i, evt->GetV0MCentrality());
        fHist->FillSPDTrackletsVsCluster(
            i, evt->GetSPDMult(),
            evt->GetSPDCluster(0) + evt->GetSPDCluster(1));
        fHist->FillSPDTrackletsLyVsCluster(i, 0, evt->GetSPDMult(),
                                           evt->GetSPDCluster(0));
        fHist->FillSPDTrackletsLyVsCluster(i, 1, evt->GetSPDMult(),
                                           evt->GetSPDCluster(1));
        fHist->FillEvtVtxZTrackvsSPD(i, evt->GetZVertexSPD(),
                                     evt->GetZVertexTracks());
        fHist->FillMagneticField(i, evt->GetBField());
        fHist->FillEvtSpher(i, evt->GetSpher());
        fHist->FillEvtSphero(i, evt->GetSphero());
        fHist->FillVZEROtiming(i, evt->GetV0ATime(), evt->GetV0CTime());
      }
    }

    if (fCentVsMultPlots) {
      fHist->FillCentVsMultV0A(evt->GetV0MCentrality(), evt->GetV0AMult());
      fHist->FillCentVsMultV0M(evt->GetV0MCentrality(), evt->GetV0MMult());
      fHist->FillCentVsMultV0C(evt->GetV0MCentrality(), evt->GetV0CMult());
      fHist->FillCentVsMultRef(evt->GetV0MCentrality(), evt->GetRefMult08());
    }
  }
}
void AliFemtoDreamEventCuts::BookCuts() {
  if (!fMinimalBooking) {
    if (fCutMinContrib) {
      fHist->FillCuts(0, fMinContrib);
    } else {
      fHist->FillCuts(0, 0);
    }
    if (fCutZVtx) {
      fHist->FillCuts(1, 1);
      fHist->FillCuts(2, fzVtxLow);
      fHist->FillCuts(3, fzVtxUp);
    } else {
      fHist->FillCuts(1, 0);
      fHist->FillCuts(2, 0);
      fHist->FillCuts(3, 0);
    }
    if (fPileUpRejection) {
      fHist->FillCuts(4, 1);
    } else {
      fHist->FillCuts(4, 0);
    }
    if (fUseMVPileUpRej) {
      fHist->FillCuts(5, 1);
    } else {
      fHist->FillCuts(5, 0);
    }
    if (fCleanEvtMult) {
      if (fUseSPDMult) {
        fHist->FillCuts(6, 1);
      } else {
        fHist->FillCuts(6, 0);
      }
      if (fUseV0AMult) {
        fHist->FillCuts(7, 1);
      } else {
        fHist->FillCuts(7, 0);
      }
      if (fUseV0CMult) {
        fHist->FillCuts(8, 1);
      } else {
        fHist->FillCuts(8, 0);
      }
      if (fUseRef08Mult) {
        fHist->FillCuts(9, 1);
      } else {
        fHist->FillCuts(9, 0);
      }
    } else {
      fHist->FillCuts(6, 0);
      fHist->FillCuts(7, 0);
      fHist->FillCuts(8, 0);
      fHist->FillCuts(9, 0);
    }
    if (fUseAliEvtCuts) {
      fHist->FillCuts(10, 1);
    } else {
      fHist->FillCuts(10, 0);
    }
    if (fDoSpherCuts) {
      fHist->FillCuts(11, fSpherCutsLow);
      fHist->FillCuts(12, fSpherCutsUp);
    } else {
      fHist->FillCuts(11, 0);
      fHist->FillCuts(12, 0);
    }
    if (fUseMultPercentileCut) {
      fHist->FillCuts(13, fMultPercentileMax);
    } else {
      fHist->FillCuts(13, 0);
    }
    if (fDoSpheroCuts) {
      fHist->FillCuts(14, fSpheroCutsLow);
      fHist->FillCuts(15, fSpheroCutsUp);
    } else {
      fHist->FillCuts(14, 0);
      fHist->FillCuts(15, 0);
    }
  }
}
