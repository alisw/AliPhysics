//_________________________________________________________________________
//  Utility Class for transverse energy studies, charged hadrons
//  Base class for ESD analysis
//  - reconstruction output
//
//Created by Christine Nattrass, Rebecca Scott, Irakli Martashvili
//University of Tennessee at Knoxville
//_________________________________________________________________________
#ifndef ALIANALYSISHADETRECONSTRUCTED_H
#define ALIANALYSISHADETRECONSTRUCTED_H

#include "AliAnalysisHadEt.h"

class AliVParticle;
class AliAnalysisHadEtCorrections;
class TString;

class AliAnalysisHadEtReconstructed : public AliAnalysisHadEt
{

public:
   
    AliAnalysisHadEtReconstructed();
    virtual ~AliAnalysisHadEtReconstructed();
   
    virtual void SetConfigFile(const char *c) {fConfigFile = c;}
    virtual Int_t AnalyseEvent(AliVEvent* event);

    void CreateHistograms();
    virtual void Init();
    
protected:

    bool CheckGoodVertex(AliVParticle *track);
    AliAnalysisHadEtCorrections *corrections;

    TString       fConfigFile;        // the name of the ConfigFile
    //virtual bool TrackHitsCalorimeter(AliVParticle *track, Double_t magField) = 0;

    Float_t fCorrTotEtFullAcceptanceTPC;
    Float_t fCorrTotEtFullAcceptanceITS;
    Float_t fCorrHadEtFullAcceptanceTPC;
    Float_t fCorrHadEtFullAcceptanceITS;
    Float_t fCorrTotEtEMCALAcceptanceTPC;
    Float_t fCorrTotEtEMCALAcceptanceITS;
    Float_t fCorrHadEtEMCALAcceptanceTPC;
    Float_t fCorrHadEtEMCALAcceptanceITS;
    Float_t fCorrTotEtPHOSAcceptanceTPC;
    Float_t fCorrTotEtPHOSAcceptanceITS;
    Float_t fCorrHadEtPHOSAcceptanceTPC;
    Float_t fCorrHadEtPHOSAcceptanceITS;
    Float_t fRawEtFullAcceptanceTPC;
    Float_t fRawEtFullAcceptanceITS;
    Float_t fRawEtEMCALAcceptanceTPC;
    Float_t fRawEtEMCALAcceptanceITS;
    Float_t fRawEtPHOSAcceptanceTPC;
    Float_t fRawEtPHOSAcceptanceITS;


 private:
    //Declare it private to avoid compilation warning
    AliAnalysisHadEtReconstructed & operator = (const AliAnalysisHadEtReconstructed & g) ;//cpy assignment
    AliAnalysisHadEtReconstructed(const AliAnalysisHadEtReconstructed & g) ; // cpy ctor


    ClassDef(AliAnalysisHadEtReconstructed, 1);
};

#endif // ALIANALYSISHADETRECONSTRUCTED_H
