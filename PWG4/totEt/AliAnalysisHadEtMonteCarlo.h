//_________________________________________________________________________
//  Utility Class for transverse energy studies, charged hadrons
//  Base class for MC analysis
//  - MC output
// 
//Created by Christine Nattrass, Rebecca Scott, Irakli Martashvili
//University of Tennessee at Knoxville
//_________________________________________________________________________
#ifndef ALIANALYSISHADETMONTECARLO_H
#define ALIANALYSISHADETMONTECARLO_H

#include "AliAnalysisHadEt.h"
class AliVEvent;

class AliAnalysisHadEtMonteCarlo : public AliAnalysisHadEt
{

public:
   
  AliAnalysisHadEtMonteCarlo();
  virtual ~AliAnalysisHadEtMonteCarlo() {}
   
    virtual Int_t AnalyseEvent(AliVEvent* event);
    virtual Int_t AnalyseEvent(AliVEvent* event,AliVEvent* event2);

    //void FillHistograms();
    void CreateHistograms();
    virtual void Init();

    Float_t GetSimulatedHadronicEt(){return fSimHadEt;}
    Float_t GetSimulatedTotalEt(){return fSimTotEt;}

    void FillSimTotEtVsRecoTotEtFullAcceptanceTPC(Float_t et){FillHisto2D("SimTotEtVsRecoTotEtFullAcceptanceTPC",fSimTotEt,et,1.0);}
    void FillSimTotEtVsRecoTotEtFullAcceptanceITS(Float_t et){FillHisto2D("SimTotEtVsRecoTotEtFullAcceptanceITS",fSimTotEt,et,1.0);}
    void FillSimTotEtVsRecoTotEtEMCALAcceptanceTPC(Float_t et){FillHisto2D("SimTotEtVsRecoTotEtEMCALAcceptanceTPC",fSimTotEt,et,1.0);}
    void FillSimTotEtVsRecoTotEtEMCALAcceptanceITS(Float_t et){FillHisto2D("SimTotEtVsRecoTotEtEMCALAcceptanceITS",fSimTotEt,et,1.0);}
    void FillSimTotEtVsRecoTotEtPHOSAcceptanceTPC(Float_t et){FillHisto2D("SimTotEtVsRecoTotEtPHOSAcceptanceTPC",fSimTotEt,et,1.0);}
    void FillSimTotEtVsRecoTotEtPHOSAcceptanceITS(Float_t et){FillHisto2D("SimTotEtVsRecoTotEtPHOSAcceptanceITS",fSimTotEt,et,1.0);}
    void FillSimTotEtVsRecoTotEtFullAcceptanceTPCNoPID(Float_t et){FillHisto2D("SimTotEtVsRecoTotEtFullAcceptanceTPCNoPID",fSimTotEt,et,1.0);}
    void FillSimTotEtVsRecoTotEtFullAcceptanceITSNoPID(Float_t et){FillHisto2D("SimTotEtVsRecoTotEtFullAcceptanceITSNoPID",fSimTotEt,et,1.0);}
    void FillSimTotEtVsRecoTotEtEMCALAcceptanceTPCNoPID(Float_t et){FillHisto2D("SimTotEtVsRecoTotEtEMCALAcceptanceTPCNoPID",fSimTotEt,et,1.0);}
    void FillSimTotEtVsRecoTotEtEMCALAcceptanceITSNoPID(Float_t et){FillHisto2D("SimTotEtVsRecoTotEtEMCALAcceptanceITSNoPID",fSimTotEt,et,1.0);}
    void FillSimTotEtVsRecoTotEtPHOSAcceptanceTPCNoPID(Float_t et){FillHisto2D("SimTotEtVsRecoTotEtPHOSAcceptanceTPCNoPID",fSimTotEt,et,1.0);}
    void FillSimTotEtVsRecoTotEtPHOSAcceptanceITSNoPID(Float_t et){FillHisto2D("SimTotEtVsRecoTotEtPHOSAcceptanceITSNoPID",fSimTotEt,et,1.0);}
    void FillSimHadEtVsRecoHadEtFullAcceptanceTPC(Float_t et){FillHisto2D("SimHadEtVsRecoHadEtFullAcceptanceTPC",fSimHadEt,et,1.0);}
    void FillSimHadEtVsRecoHadEtFullAcceptanceITS(Float_t et){FillHisto2D("SimHadEtVsRecoHadEtFullAcceptanceITS",fSimHadEt,et,1.0);}
    void FillSimHadEtVsRecoHadEtEMCALAcceptanceTPC(Float_t et){FillHisto2D("SimHadEtVsRecoHadEtEMCALAcceptanceTPC",fSimHadEt,et,1.0);}
    void FillSimHadEtVsRecoHadEtEMCALAcceptanceITS(Float_t et){FillHisto2D("SimHadEtVsRecoHadEtEMCALAcceptanceITS",fSimHadEt,et,1.0);}
    void FillSimHadEtVsRecoHadEtPHOSAcceptanceTPC(Float_t et){FillHisto2D("SimHadEtVsRecoHadEtPHOSAcceptanceTPC",fSimHadEt,et,1.0);}
    void FillSimHadEtVsRecoHadEtPHOSAcceptanceITS(Float_t et){FillHisto2D("SimHadEtVsRecoHadEtPHOSAcceptanceITS",fSimHadEt,et,1.0);}
    void FillSimHadEtVsRecoHadEtFullAcceptanceTPCNoPID(Float_t et){FillHisto2D("SimHadEtVsRecoHadEtFullAcceptanceTPCNoPID",fSimHadEt,et,1.0);}
    void FillSimHadEtVsRecoHadEtFullAcceptanceITSNoPID(Float_t et){FillHisto2D("SimHadEtVsRecoHadEtFullAcceptanceITSNoPID",fSimHadEt,et,1.0);}
    void FillSimHadEtVsRecoHadEtEMCALAcceptanceTPCNoPID(Float_t et){FillHisto2D("SimHadEtVsRecoHadEtEMCALAcceptanceTPCNoPID",fSimHadEt,et,1.0);}
    void FillSimHadEtVsRecoHadEtEMCALAcceptanceITSNoPID(Float_t et){FillHisto2D("SimHadEtVsRecoHadEtEMCALAcceptanceITSNoPID",fSimHadEt,et,1.0);}
    void FillSimHadEtVsRecoHadEtPHOSAcceptanceTPCNoPID(Float_t et){FillHisto2D("SimHadEtVsRecoHadEtPHOSAcceptanceTPCNoPID",fSimHadEt,et,1.0);}
    void FillSimHadEtVsRecoHadEtPHOSAcceptanceITSNoPID(Float_t et){FillHisto2D("SimHadEtVsRecoHadEtPHOSAcceptanceITSNoPID",fSimHadEt,et,1.0);}
 private:
    //Declare it private to avoid compilation warning
    AliAnalysisHadEtMonteCarlo & operator = (const AliAnalysisHadEtMonteCarlo & g) ;//cpy assignment
    AliAnalysisHadEtMonteCarlo(const AliAnalysisHadEtMonteCarlo & g) ; // cpy ctor

    Float_t fSimHadEt;
    Float_t fSimTotEt;

    void ResetEventValues();
    ClassDef(AliAnalysisHadEtMonteCarlo, 1);
};

#endif // ALIANALYSISHADETMONTECARLO_H
