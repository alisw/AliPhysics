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
class AliESDtrack;

class AliAnalysisHadEtReconstructed : public AliAnalysisHadEt
{

public:
   
    AliAnalysisHadEtReconstructed();
    virtual ~AliAnalysisHadEtReconstructed();
   
    virtual void SetConfigFile(const char *c) {fConfigFile = c;}
    virtual Int_t AnalyseEvent(AliVEvent* event);

    //the "Corrected" variables are only corrected for the track-by-track corrections
    Float_t GetCorrectedHadEtFullAcceptanceTPC(){return fCorrHadEtFullAcceptanceTPC*fCorrectedHadEtFullAcceptanceTPC;}
    Float_t GetCorrectedHadEtFullAcceptanceITS(){return fCorrHadEtFullAcceptanceITS*(fCorrectedHadEtFullAcceptanceTPC+fCorrectedHadEtFullAcceptanceITS);}
    Float_t GetCorrectedHadEtEMCALAcceptanceTPC(){return fCorrHadEtEMCALAcceptanceTPC*fCorrectedHadEtEMCALAcceptanceTPC;}
    Float_t GetCorrectedHadEtEMCALAcceptanceITS(){return fCorrHadEtEMCALAcceptanceITS*(fCorrectedHadEtEMCALAcceptanceTPC+fCorrectedHadEtEMCALAcceptanceITS);}
    Float_t GetCorrectedHadEtPHOSAcceptanceTPC(){return fCorrHadEtPHOSAcceptanceTPC*fCorrectedHadEtPHOSAcceptanceTPC;}
    Float_t GetCorrectedHadEtPHOSAcceptanceITS(){return fCorrHadEtPHOSAcceptanceITS*(fCorrectedHadEtPHOSAcceptanceTPC+fCorrectedHadEtPHOSAcceptanceITS);}
    Float_t GetCorrectedTotEtFullAcceptanceTPC(){return fCorrTotEtFullAcceptanceTPC*fCorrectedHadEtFullAcceptanceTPC;}
    Float_t GetCorrectedTotEtFullAcceptanceITS(){return fCorrTotEtFullAcceptanceITS*(fCorrectedHadEtFullAcceptanceTPC+fCorrectedHadEtFullAcceptanceITS);}
    Float_t GetCorrectedTotEtEMCALAcceptanceTPC(){return fCorrTotEtEMCALAcceptanceTPC*fCorrectedHadEtEMCALAcceptanceTPC;}
    Float_t GetCorrectedTotEtEMCALAcceptanceITS(){return fCorrTotEtEMCALAcceptanceITS*(fCorrectedHadEtEMCALAcceptanceTPC+fCorrectedHadEtEMCALAcceptanceITS);}
    Float_t GetCorrectedTotEtPHOSAcceptanceTPC(){return fCorrTotEtPHOSAcceptanceTPC*fCorrectedHadEtPHOSAcceptanceTPC;}
    Float_t GetCorrectedTotEtPHOSAcceptanceITS(){return fCorrTotEtPHOSAcceptanceTPC*(fCorrectedHadEtPHOSAcceptanceTPC+fCorrectedHadEtPHOSAcceptanceITS);}
    Float_t GetRawEtFullAcceptanceTPC(){return fRawEtFullAcceptanceTPC;}
    Float_t GetRawEtFullAcceptanceITS(){return fRawEtFullAcceptanceITS+fRawEtFullAcceptanceTPC;}
    Float_t GetRawEtEMCALAcceptanceTPC(){return fRawEtEMCALAcceptanceTPC;}
    Float_t GetRawEtEMCALAcceptanceITS(){return fRawEtEMCALAcceptanceITS+fRawEtEMCALAcceptanceTPC;}
    Float_t GetRawEtPHOSAcceptanceTPC(){return fRawEtPHOSAcceptanceTPC;}
    Float_t GetRawEtPHOSAcceptanceITS(){return fRawEtPHOSAcceptanceITS+fRawEtPHOSAcceptanceTPC;}
    Float_t GetCorrectedHadEtFullAcceptanceTPCNoPID(){return fCorrHadEtFullAcceptanceTPC*fCorrectedHadEtFullAcceptanceTPCNoPID;}
    Float_t GetCorrectedHadEtFullAcceptanceITSNoPID(){return fCorrHadEtFullAcceptanceITS*(fCorrectedHadEtFullAcceptanceTPCNoPID+fCorrectedHadEtFullAcceptanceITSNoPID);}
    Float_t GetCorrectedHadEtEMCALAcceptanceTPCNoPID(){return fCorrHadEtEMCALAcceptanceTPC*fCorrectedHadEtEMCALAcceptanceTPCNoPID;}
    Float_t GetCorrectedHadEtEMCALAcceptanceITSNoPID(){return fCorrHadEtEMCALAcceptanceITS*(fCorrectedHadEtEMCALAcceptanceTPCNoPID+fCorrectedHadEtEMCALAcceptanceITSNoPID);}
    Float_t GetCorrectedHadEtPHOSAcceptanceTPCNoPID(){return fCorrHadEtPHOSAcceptanceTPC*fCorrectedHadEtPHOSAcceptanceTPCNoPID;}
    Float_t GetCorrectedHadEtPHOSAcceptanceITSNoPID(){return fCorrHadEtPHOSAcceptanceITS*(fCorrectedHadEtPHOSAcceptanceTPCNoPID+fCorrectedHadEtPHOSAcceptanceITSNoPID);}
    Float_t GetCorrectedTotEtFullAcceptanceTPCNoPID(){return fCorrTotEtFullAcceptanceTPC*fCorrectedHadEtFullAcceptanceTPCNoPID;}
    Float_t GetCorrectedTotEtFullAcceptanceITSNoPID(){return fCorrTotEtFullAcceptanceITS*(fCorrectedHadEtFullAcceptanceTPCNoPID+fCorrectedHadEtFullAcceptanceITSNoPID);}
    Float_t GetCorrectedTotEtEMCALAcceptanceTPCNoPID(){return fCorrTotEtEMCALAcceptanceTPC*fCorrectedHadEtEMCALAcceptanceTPCNoPID;}
    Float_t GetCorrectedTotEtEMCALAcceptanceITSNoPID(){return fCorrTotEtEMCALAcceptanceITS*(fCorrectedHadEtEMCALAcceptanceTPCNoPID+fCorrectedHadEtEMCALAcceptanceITSNoPID);}
    Float_t GetCorrectedTotEtPHOSAcceptanceTPCNoPID(){return fCorrTotEtPHOSAcceptanceTPC*fCorrectedHadEtPHOSAcceptanceTPCNoPID;}
    Float_t GetCorrectedTotEtPHOSAcceptanceITSNoPID(){return fCorrTotEtPHOSAcceptanceITS*(fCorrectedHadEtPHOSAcceptanceTPCNoPID+fCorrectedHadEtPHOSAcceptanceITSNoPID);}
    Float_t GetRawEtFullAcceptanceTPCNoPID(){return fRawEtFullAcceptanceTPCNoPID;}
    Float_t GetRawEtFullAcceptanceITSNoPID(){return fRawEtFullAcceptanceITSNoPID+fRawEtFullAcceptanceTPCNoPID;}
    Float_t GetRawEtEMCALAcceptanceTPCNoPID(){return fRawEtEMCALAcceptanceTPCNoPID;}
    Float_t GetRawEtEMCALAcceptanceITSNoPID(){return fRawEtEMCALAcceptanceITSNoPID+fRawEtEMCALAcceptanceTPCNoPID;}
    Float_t GetRawEtPHOSAcceptanceTPCNoPID(){return fRawEtPHOSAcceptanceTPCNoPID;}
    Float_t GetRawEtPHOSAcceptanceITSNoPID(){return fRawEtPHOSAcceptanceITSNoPID+fRawEtPHOSAcceptanceTPCNoPID;}

    void CreateHistograms();
    virtual void Init();
    
protected:

    Bool_t CheckGoodVertex(AliVParticle *track);
    AliAnalysisHadEtCorrections *corrections;

    TString       fConfigFile;        // the name of the ConfigFile
    //virtual bool TrackHitsCalorimeter(AliVParticle *track, Double_t magField) = 0;

    //correction factors
    Float_t fCorrTotEtFullAcceptanceTPC;//get the correction for total et for full acceptance, pt>0.15 GeV/c
    Float_t fCorrTotEtFullAcceptanceITS;//get the correction for total et for full acceptance, pt>0.10 GeV/c
    Float_t fCorrHadEtFullAcceptanceTPC;//get the correction for hadronic et for full acceptance, pt>0.15 GeV/c
    Float_t fCorrHadEtFullAcceptanceITS;//get the correction for hadronic et for full acceptance, pt>0.10 GeV/c
    Float_t fCorrTotEtEMCALAcceptanceTPC;//analogous to above for EMCal acceptance
    Float_t fCorrTotEtEMCALAcceptanceITS;//analogous to above for EMCal acceptance
    Float_t fCorrHadEtEMCALAcceptanceTPC;//analogous to above for EMCal acceptance
    Float_t fCorrHadEtEMCALAcceptanceITS;//analogous to above for EMCal acceptance
    Float_t fCorrTotEtPHOSAcceptanceTPC;//analogous to above for PHOS acceptance
    Float_t fCorrTotEtPHOSAcceptanceITS;//analogous to above for PHOS acceptance
    Float_t fCorrHadEtPHOSAcceptanceTPC;//analogous to above for PHOS acceptance
    Float_t fCorrHadEtPHOSAcceptanceITS;//analogous to above for PHOS acceptance
    //Et with various parameters...
    Float_t fCorrectedHadEtFullAcceptanceTPCNoPID;//get the corrected hadronic et for full acceptance, pt>0.15 GeV/c
    Float_t fCorrectedHadEtFullAcceptanceITSNoPID;//get the corrected hadronic et for full acceptance, pt>0.10 GeV/c
    Float_t fCorrectedHadEtEMCALAcceptanceTPCNoPID;//analogous to above for EMCal acceptance
    Float_t fCorrectedHadEtEMCALAcceptanceITSNoPID;//analogous to above for EMCal acceptance
    Float_t fCorrectedHadEtPHOSAcceptanceTPCNoPID;//analogous to above for PHOS acceptance
    Float_t fCorrectedHadEtPHOSAcceptanceITSNoPID;//analogous to above for PHOS acceptance
    Float_t fCorrectedHadEtFullAcceptanceTPC;//get the corrected hadronic et for full acceptance, pt>0.15 GeV/c
    Float_t fCorrectedHadEtFullAcceptanceITS;//get the corrected hadronic et for full acceptance, pt>0.10 GeV/c
    Float_t fCorrectedHadEtEMCALAcceptanceTPC;//analogous to above for EMCal acceptance
    Float_t fCorrectedHadEtEMCALAcceptanceITS;//analogous to above for EMCal acceptance
    Float_t fCorrectedHadEtPHOSAcceptanceTPC;//analogous to above for PHOS acceptance
    Float_t fCorrectedHadEtPHOSAcceptanceITS;//analogous to above for PHOS acceptance
    Float_t fRawEtFullAcceptanceTPC;//uncorrected Et for full acceptance, pT > 0.15 GeV/c
    Float_t fRawEtFullAcceptanceITS;//uncorrected Et for full acceptance, pT > 0.10 GeV/c
    Float_t fRawEtEMCALAcceptanceTPC;//uncorrected Et for EMCal acceptance, pT > 0.15 GeV/c
    Float_t fRawEtEMCALAcceptanceITS;//uncorrected Et for EMCal acceptance, pT > 0.10 GeV/c
    Float_t fRawEtPHOSAcceptanceTPC;//uncorrected Et for PHOS acceptance, pT > 0.15 GeV/c
    Float_t fRawEtPHOSAcceptanceITS;//uncorrected Et for PHOS acceptance, pT > 0.10 GeV/c
    Float_t fRawEtFullAcceptanceTPCNoPID;//uncorrected Et for full acceptance, pT > 0.15 GeV/c
    Float_t fRawEtFullAcceptanceITSNoPID;//uncorrected Et for full acceptance, pT > 0.10 GeV/c
    Float_t fRawEtEMCALAcceptanceTPCNoPID;//uncorrected Et for EMCal acceptance, pT > 0.15 GeV/c
    Float_t fRawEtEMCALAcceptanceITSNoPID;//uncorrected Et for EMCal acceptance, pT > 0.10 GeV/c
    Float_t fRawEtPHOSAcceptanceTPCNoPID;//uncorrected Et for PHOS acceptance, pT > 0.15 GeV/c
    Float_t fRawEtPHOSAcceptanceITSNoPID;//uncorrected Et for PHOS acceptance, pT > 0.10 GeV/c


 private:
    //Declare it private to avoid compilation warning
    AliAnalysisHadEtReconstructed & operator = (const AliAnalysisHadEtReconstructed & g) ;//cpy assignment
    AliAnalysisHadEtReconstructed(const AliAnalysisHadEtReconstructed & g) ; // cpy ctor

    void AddEt(Float_t rawEt, Float_t rawEtNoPID, Float_t corrEt, Float_t corrEtNoPID, Float_t pt, Bool_t IsTPC, Bool_t InPHOS, Bool_t InEMCAL);
    Bool_t IsInPHOS(AliESDtrack *track);
    Bool_t IsInEMCAL(AliESDtrack *track);

    void ResetEventValues();
    ClassDef(AliAnalysisHadEtReconstructed, 1);
};

#endif // ALIANALYSISHADETRECONSTRUCTED_H
