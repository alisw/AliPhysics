//Create by Christine Nattrass, Rebecca Scott, Irakli Martashvili
//University of Tennessee at Knoxville
//This is a container class for the correction factors for the hadronic component of transverse energy
//It is filled by the output of AliAnalysisTaskHadEt from spinning over Monte Carlo data (using AliAnalysisHadEtMonteCarlo)
//It is used by AliAnalysisTaskHadEt while spinning over reconstructed data (using AliAnalysisHadEtReconstructed)
//Please see https://twiki.cern.ch/twiki/bin/view/ALICE/ETCaloAnalysis
#ifndef ALIANALYSISHADETCORRECTIONS_H
#define ALIANALYSISHADETCORRECTIONS_H

#include "Rtypes.h"
#include "TString.h"
#include "TNamed.h"
#include "TH1D.h"

class AliAnalysisHadEtCorrections : public TNamed
{
public:
   
  AliAnalysisHadEtCorrections();
    virtual ~AliAnalysisHadEtCorrections();

    Float_t GetEtaCut() const {return fEtaCut;}
    Float_t GetAcceptanceCorrectionFull() const {return fAcceptanceCorrectionFull;}
    Float_t GetAcceptanceCorrectionEMCAL() const {return fAcceptanceCorrectionEMCAL;}
    Float_t GetAcceptanceCorrectionPHOS() const {return fAcceptanceCorrectionPHOS;}
    Float_t GetNeutralCorrection() const {return fNeutralCorrection;}
    Float_t GetNotHadronicCorrection() const {return fNotHadronicCorrection;}
    Float_t GetpTCutCorrectionTPC() const {return fpTcutCorrectionTPC;}
    Float_t GetpTCutCorrectionITS() const {return fpTcutCorrectionITS;}
    Float_t GetNeutralCorrectionLowBound() const {return fNeutralCorrectionLow;}
    Float_t GetNotHadronicCorrectionLowBound() const {return fNotHadronicCorrectionLow;}
    Float_t GetpTCutCorrectionTPCLowBound() const {return ffpTcutCorrectionTPCLow;}
    Float_t GetpTCutCorrectionITSLowBound() const {return ffpTcutCorrectionITSLow;}
    Float_t GetNeutralCorrectionHighBound() const {return fNeutralCorrectionHigh;}
    Float_t GetNotHadronicCorrectionHighBound() const {return fNotHadronicCorrectionHigh;}
    Float_t GetpTCutCorrectionTPCHighBound() const {return ffpTcutCorrectionTPCHigh;}
    Float_t GetpTCutCorrectionITSHighBound() const {return ffpTcutCorrectionITSHigh;}
    TH1D *GetNotIDCorrectionTPC() const {return fnotIDTPC;}
    TH1D *GetNotIDCorrectionITS() const {return fnotIDITS;}
    TH1D *GetNotIDCorrectionNoPID() const {return fnotIDNoID;}
    TH1D *GetEfficiencyPionTPC() const {return fEfficiencyPionTPC;}
    TH1D *GetEfficiencyKaonTPC() const {return fEfficiencyKaonTPC;}
    TH1D *GetEfficiencyProtonTPC() const {return fEfficiencyProtonTPC;}
    TH1D *GetEfficiencyHadronTPC() const {return fEfficiencyHadronTPC;}
    TH1D *GetEfficiencyPionITS() const {return fEfficiencyPionITS;}
    TH1D *GetEfficiencyKaonITS() const {return fEfficiencyKaonITS;}
    TH1D *GetEfficiencyProtonITS() const {return fEfficiencyProtonITS;}
    TH1D *GetEfficiencyHadronITS() const {return fEfficiencyHadronITS;}
    TH1D *GetBackgroundCorrectionTPC() const {return fBackgroundTPC;}
    TH1D *GetBackgroundCorrectionITS() const {return fBackgroundITS;}

    //This is stored as the inverse of the correction
    Float_t GetNotIDCorrectionTPC(float pT){return 1.0/(fnotIDTPC->GetBinContent(fnotIDTPC->FindBin(pT)));}
    Float_t GetNotIDCorrectionITS(float pT){return 1.0/(fnotIDITS->GetBinContent(fnotIDITS->FindBin(pT)));}
    Float_t GetNotIDCorrectionNoPID(float pT){return 1.0/(fnotIDNoID->GetBinContent(fnotIDNoID->FindBin(pT)));}
    //As is this...
    Float_t GetTPCEfficiencyCorrectionPion(float pT){return 1.0/(fEfficiencyPionTPC->GetBinContent(fEfficiencyPionTPC->FindBin(pT)));}
    Float_t GetTPCEfficiencyCorrectionKaon(float pT){return 1.0/(fEfficiencyKaonTPC->GetBinContent(fEfficiencyKaonTPC->FindBin(pT)));}
    Float_t GetTPCEfficiencyCorrectionProton(float pT){return 1.0/(fEfficiencyProtonTPC->GetBinContent(fEfficiencyProtonTPC->FindBin(pT)));}
    Float_t GetTPCEfficiencyCorrectionHadron(float pT){return 1.0/(fEfficiencyHadronTPC->GetBinContent(fEfficiencyHadronTPC->FindBin(pT)));}
    Float_t GetITSEfficiencyCorrectionPion(float pT){return 1.0/(fEfficiencyPionITS->GetBinContent(fEfficiencyPionITS->FindBin(pT)));}
    Float_t GetITSEfficiencyCorrectionKaon(float pT){return 1.0/(fEfficiencyKaonITS->GetBinContent(fEfficiencyKaonITS->FindBin(pT)));}
    Float_t GetITSEfficiencyCorrectionProton(float pT){return 1.0/(fEfficiencyProtonITS->GetBinContent(fEfficiencyProtonITS->FindBin(pT)));}
    Float_t GetITSEfficiencyCorrectionHadron(float pT){return 1.0/(fEfficiencyHadronITS->GetBinContent(fEfficiencyHadronITS->FindBin(pT)));}
    //...and these guys are too
    Float_t GetBackgroundCorrectionTPC(float pT){return 1.0/(fBackgroundTPC->GetBinContent(fBackgroundTPC->FindBin(pT)));}
    Float_t GetBackgroundCorrectionITS(float pT){return 1.0/(fBackgroundITS->GetBinContent(fBackgroundITS->FindBin(pT)));}


    void SetEtaCut(Float_t val){fEtaCut=val;}
    void SetAcceptanceCorrectionFull(Float_t val){fAcceptanceCorrectionFull=val;}
    void SetAcceptanceCorrectionEMCAL(Float_t val){fAcceptanceCorrectionEMCAL=val;}
    void SetAcceptanceCorrectionPHOS(Float_t val){fAcceptanceCorrectionPHOS=val;}
    void SetNeutralCorrection(Float_t val){fNeutralCorrection=val;}
    void SetNotHadronicCorrection(Float_t val){fNotHadronicCorrection=val;}
    void SetpTCutCorrectionTPC(Float_t val){fpTcutCorrectionTPC=val;}
    void SetpTCutCorrectionITS(Float_t val){fpTcutCorrectionITS=val;}
    void SetNeutralCorrectionLowBound(Float_t val){fNeutralCorrectionLow=val;}
    void SetNotHadronicCorrectionLowBound(Float_t val){fNotHadronicCorrectionLow=val;}
    void SetpTCutCorrectionTPCLowBound(Float_t val){ffpTcutCorrectionTPCLow=val;}
    void SetpTCutCorrectionITSLowBound(Float_t val){ffpTcutCorrectionITSLow=val;}
    void SetNeutralCorrectionHighBound(Float_t val){fNeutralCorrectionHigh=val;}
    void SetNotHadronicCorrectionHighBound(Float_t val){fNotHadronicCorrectionHigh=val;}
    void SetpTCutCorrectionTPCHighBound(Float_t val){ffpTcutCorrectionTPCHigh=val;}
    void SetpTCutCorrectionITSHighBound(Float_t val){ffpTcutCorrectionITSHigh=val;}
    void SetNotIDCorrectionTPC(TH1D *histo){fnotIDTPC=histo;}
    void SetNotIDCorrectionITS(TH1D *histo){fnotIDITS=histo;}
    void SetNotIDCorrectionNoPID(TH1D *histo){fnotIDNoID=histo;}
    void SetEfficiencyPionTPC(TH1D *histo){fEfficiencyPionTPC=histo;}
    void SetEfficiencyKaonTPC(TH1D *histo){fEfficiencyKaonTPC=histo;}
    void SetEfficiencyProtonTPC(TH1D *histo){fEfficiencyProtonTPC=histo;}
    void SetEfficiencyHadronTPC(TH1D *histo){fEfficiencyHadronTPC=histo;}
    void SetEfficiencyPionITS(TH1D *histo){fEfficiencyPionITS=histo;}
    void SetEfficiencyKaonITS(TH1D *histo){fEfficiencyKaonITS=histo;}
    void SetEfficiencyProtonITS(TH1D *histo){fEfficiencyProtonITS=histo;}
    void SetEfficiencyHadronITS(TH1D *histo){fEfficiencyHadronITS=histo;}
    void SetBackgroundCorrectionTPC(TH1D *histo){fBackgroundTPC=histo;}
    void SetBackgroundCorrectionITS(TH1D *histo){fBackgroundITS=histo;}


    AliAnalysisHadEtCorrections(const AliAnalysisHadEtCorrections *g) ; // cpy ctor
    //AliAnalysisHadEtCorrections & operator = (const AliAnalysisHadEtCorrections & g) ;//cpy assignment

protected:

    Float_t fEtaCut;//the eta cut used for this analysis
    Float_t fAcceptanceCorrectionFull;//the acceptance correction for full azimuthal acceptance
    Float_t fAcceptanceCorrectionEMCAL;//the acceptance correction for the EMCal azimuthal acceptance
    Float_t fAcceptanceCorrectionPHOS;//the acceptance correction for the PHOS azimuthal acceptance
    //Systematic errors - low value, mean value, high value
    Float_t fNeutralCorrection;//the correction for the fraction of energy from neutral particles (for using both the calorimeters and the tracking detectors)
    Float_t fNotHadronicCorrection;//the correction for the fraction of energy which is not measured by the tracking detectors 
    Float_t fpTcutCorrectionTPC;//the correction for the momentum cut for the tpc (150 MeV/c)
    Float_t fpTcutCorrectionITS;//the correction for the momentum cut for the ITS (100 MeV/c)
    Float_t fNeutralCorrectionLow;//the low bound on the neutral energy fraction correction
    Float_t fNotHadronicCorrectionLow;//the low bound on the hadronic energy fraction correction
    Float_t ffpTcutCorrectionTPCLow;//the low bound on the TPC momentum cut correction
    Float_t ffpTcutCorrectionITSLow;//the low bound on the ITS momentum cut correction
    Float_t fNeutralCorrectionHigh;//the high bound on the neutral energy correcton
    Float_t fNotHadronicCorrectionHigh;//the high bound on the hadronic energy correction
    Float_t ffpTcutCorrectionTPCHigh;//the high bound on the TPC momentum cut correction
    Float_t ffpTcutCorrectionITSHigh;//the high bound on the ITS momentum cut correction

    //Histograms with the pT dependent corrections
    TH1D *fnotIDTPC;//correction for unidentified tracks in the TPC
    TH1D *fnotIDITS;//correction for unidentified tracks in the ITS
    TH1D *fnotIDNoID;//correction for unidentified tracks assuming no ID
    TH1D *fEfficiencyPionTPC;//efficiency correction for pions in the TPC
    TH1D *fEfficiencyKaonTPC;//efficiency correction for kaons in the TPC
    TH1D *fEfficiencyProtonTPC;//efficiency correction for protons in the TPC
    TH1D *fEfficiencyHadronTPC;//efficiency correction for unidentified hadrons in the TPC
    TH1D *fEfficiencyPionITS;//efficiency correction for pions in the ITS
    TH1D *fEfficiencyKaonITS;//efficiency correction for kaons in the ITS
    TH1D *fEfficiencyProtonITS;//efficiency correction for protons in the ITS
    TH1D *fEfficiencyHadronITS;//efficiency correction for unidentified hadrons in the ITS
    TH1D *fBackgroundTPC;//background correction for the TPC
    TH1D *fBackgroundITS;//background correction for the ITS



private:
  //Declare it private to avoid compilation warning
    AliAnalysisHadEtCorrections & operator = (const AliAnalysisHadEtCorrections & g) ;//cpy assignment
    AliAnalysisHadEtCorrections(const AliAnalysisHadEtCorrections & g) ; // cpy ctor

    ClassDef(AliAnalysisHadEtCorrections, 1);
};

#endif // ALIANALYSISHADETCORRECTIONS_H
