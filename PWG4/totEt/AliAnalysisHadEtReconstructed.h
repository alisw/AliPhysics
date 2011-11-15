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
//#include "PWG0/AliPWG0Helper.h"

class AliVParticle;
class AliAnalysisHadEtCorrections;
class TString;
class AliESDtrack;

class AliAnalysisHadEtReconstructed : public AliAnalysisHadEt
{

public:
   
    AliAnalysisHadEtReconstructed();
    virtual ~AliAnalysisHadEtReconstructed();
   
    virtual Int_t AnalyseEvent(AliVEvent* event, Int_t eventtype);
    virtual Int_t AnalyseEvent(AliVEvent* event){return AnalyseEvent(event,-1);}

    //the "Corrected" variables are only corrected for the track-by-track fCorrections
    Float_t GetCorrectedHadEtFullAcceptanceTPC() const {return fCorrHadEtFullAcceptanceTPC*fCorrectedHadEtFullAcceptanceTPC;}
    Float_t GetCorrectedHadEtFullAcceptanceITS() const {return fCorrHadEtFullAcceptanceITS*(fCorrectedHadEtFullAcceptanceTPC+fCorrectedHadEtFullAcceptanceITS);}
    Float_t GetCorrectedHadEtFullAcceptanceTPCAssumingPion() const {return fCorrHadEtFullAcceptanceTPC*fCorrectedHadEtFullAcceptanceTPCAssumingPion;}
    Float_t GetCorrectedHadEtFullAcceptanceITSAssumingPion() const {return fCorrHadEtFullAcceptanceITS*(fCorrectedHadEtFullAcceptanceTPCAssumingPion+fCorrectedHadEtFullAcceptanceITSAssumingPion);}
    Float_t GetCorrectedHadEtFullAcceptanceTPCAssumingProton() const {return fCorrHadEtFullAcceptanceTPC*fCorrectedHadEtFullAcceptanceTPCAssumingProton;}
    Float_t GetCorrectedHadEtFullAcceptanceITSAssumingProton() const {return fCorrHadEtFullAcceptanceITS*(fCorrectedHadEtFullAcceptanceTPCAssumingProton+fCorrectedHadEtFullAcceptanceITSAssumingProton);}
    Float_t GetCorrectedHadEtFullAcceptanceTPCAssumingKaon() const {return fCorrHadEtFullAcceptanceTPC*fCorrectedHadEtFullAcceptanceTPCAssumingKaon;}
    Float_t GetCorrectedHadEtFullAcceptanceITSAssumingKaon() const {return fCorrHadEtFullAcceptanceITS*(fCorrectedHadEtFullAcceptanceTPCAssumingKaon+fCorrectedHadEtFullAcceptanceITSAssumingKaon);}
    Float_t GetCorrectedHadEtEMCALAcceptanceTPC() const{return fCorrHadEtEMCALAcceptanceTPC*fCorrectedHadEtEMCALAcceptanceTPC;}
    Float_t GetCorrectedHadEtEMCALAcceptanceITS() const {return fCorrHadEtEMCALAcceptanceITS*(fCorrectedHadEtEMCALAcceptanceTPC+fCorrectedHadEtEMCALAcceptanceITS);}
    Float_t GetCorrectedHadEtPHOSAcceptanceTPC() const {return fCorrHadEtPHOSAcceptanceTPC*fCorrectedHadEtPHOSAcceptanceTPC;}
    Float_t GetCorrectedHadEtPHOSAcceptanceITS() const {return fCorrHadEtPHOSAcceptanceITS*(fCorrectedHadEtPHOSAcceptanceTPC+fCorrectedHadEtPHOSAcceptanceITS);}
    Float_t GetCorrectedTotEtFullAcceptanceTPC() const {return fCorrTotEtFullAcceptanceTPC*fCorrectedHadEtFullAcceptanceTPC;}
    Float_t GetCorrectedTotEtFullAcceptanceITS() const {return fCorrTotEtFullAcceptanceITS*(fCorrectedHadEtFullAcceptanceTPC+fCorrectedHadEtFullAcceptanceITS);}
    Float_t GetCorrectedTotEtEMCALAcceptanceTPC() const {return fCorrTotEtEMCALAcceptanceTPC*fCorrectedHadEtEMCALAcceptanceTPC;}
    Float_t GetCorrectedTotEtEMCALAcceptanceITS() const {return fCorrTotEtEMCALAcceptanceITS*(fCorrectedHadEtEMCALAcceptanceTPC+fCorrectedHadEtEMCALAcceptanceITS);}
    Float_t GetCorrectedTotEtPHOSAcceptanceTPC() const {return fCorrTotEtPHOSAcceptanceTPC*fCorrectedHadEtPHOSAcceptanceTPC;}
    Float_t GetCorrectedTotEtPHOSAcceptanceITS() const {return fCorrTotEtPHOSAcceptanceTPC*(fCorrectedHadEtPHOSAcceptanceTPC+fCorrectedHadEtPHOSAcceptanceITS);}
    Float_t GetRawEtFullAcceptanceTPC() const {return fRawEtFullAcceptanceTPC;}
    Float_t GetRawEtFullAcceptanceITS() const {return fRawEtFullAcceptanceITS+fRawEtFullAcceptanceTPC;}
    Float_t GetRawEtEMCALAcceptanceTPC() const {return fRawEtEMCALAcceptanceTPC;}
    Float_t GetRawEtEMCALAcceptanceITS() const {return fRawEtEMCALAcceptanceITS+fRawEtEMCALAcceptanceTPC;}
    Float_t GetRawEtPHOSAcceptanceTPC() const {return fRawEtPHOSAcceptanceTPC;}
    Float_t GetRawEtPHOSAcceptanceITS()const {return fRawEtPHOSAcceptanceITS+fRawEtPHOSAcceptanceTPC;}
    Float_t GetCorrectedHadEtFullAcceptanceTPCNoPID() const {return fCorrHadEtFullAcceptanceTPC*fCorrectedHadEtFullAcceptanceTPCNoPID;}
    Float_t GetCorrectedHadEtFullAcceptanceITSNoPID() const {return fCorrHadEtFullAcceptanceITS*(fCorrectedHadEtFullAcceptanceTPCNoPID+fCorrectedHadEtFullAcceptanceITSNoPID);}
    Float_t GetCorrectedHadEtEMCALAcceptanceTPCNoPID() const {return fCorrHadEtEMCALAcceptanceTPC*fCorrectedHadEtEMCALAcceptanceTPCNoPID;}
    Float_t GetCorrectedHadEtEMCALAcceptanceITSNoPID() const {return fCorrHadEtEMCALAcceptanceITS*(fCorrectedHadEtEMCALAcceptanceTPCNoPID+fCorrectedHadEtEMCALAcceptanceITSNoPID);}
    Float_t GetCorrectedHadEtPHOSAcceptanceTPCNoPID() const {return fCorrHadEtPHOSAcceptanceTPC*fCorrectedHadEtPHOSAcceptanceTPCNoPID;}
    Float_t GetCorrectedHadEtPHOSAcceptanceITSNoPID() const {return fCorrHadEtPHOSAcceptanceITS*(fCorrectedHadEtPHOSAcceptanceTPCNoPID+fCorrectedHadEtPHOSAcceptanceITSNoPID);}
    Float_t GetCorrectedTotEtFullAcceptanceTPCNoPID() const {return fCorrTotEtFullAcceptanceTPC*fCorrectedHadEtFullAcceptanceTPCNoPID;}
    Float_t GetCorrectedTotEtFullAcceptanceITSNoPID() const {return fCorrTotEtFullAcceptanceITS*(fCorrectedHadEtFullAcceptanceTPCNoPID+fCorrectedHadEtFullAcceptanceITSNoPID);}
    Float_t GetCorrectedTotEtEMCALAcceptanceTPCNoPID() const {return fCorrTotEtEMCALAcceptanceTPC*fCorrectedHadEtEMCALAcceptanceTPCNoPID;}
    Float_t GetCorrectedTotEtEMCALAcceptanceITSNoPID() const {return fCorrTotEtEMCALAcceptanceITS*(fCorrectedHadEtEMCALAcceptanceTPCNoPID+fCorrectedHadEtEMCALAcceptanceITSNoPID);}
    Float_t GetCorrectedTotEtPHOSAcceptanceTPCNoPID() const {return fCorrTotEtPHOSAcceptanceTPC*fCorrectedHadEtPHOSAcceptanceTPCNoPID;}
    Float_t GetCorrectedTotEtPHOSAcceptanceITSNoPID() const {return fCorrTotEtPHOSAcceptanceITS*(fCorrectedHadEtPHOSAcceptanceTPCNoPID+fCorrectedHadEtPHOSAcceptanceITSNoPID);}
    Float_t GetRawEtFullAcceptanceTPCNoPID() const {return fRawEtFullAcceptanceTPCNoPID;}
    Float_t GetRawEtFullAcceptanceITSNoPID() const {return fRawEtFullAcceptanceITSNoPID+fRawEtFullAcceptanceTPCNoPID;}
    Float_t GetRawEtEMCALAcceptanceTPCNoPID() const {return fRawEtEMCALAcceptanceTPCNoPID;}
    Float_t GetRawEtEMCALAcceptanceITSNoPID() const {return fRawEtEMCALAcceptanceITSNoPID+fRawEtEMCALAcceptanceTPCNoPID;}
    Float_t GetRawEtPHOSAcceptanceTPCNoPID() const {return fRawEtPHOSAcceptanceTPCNoPID;}
    Float_t GetRawEtPHOSAcceptanceITSNoPID() const {return fRawEtPHOSAcceptanceITSNoPID+fRawEtPHOSAcceptanceTPCNoPID;}
    Float_t GetCorrectedPiKPEtFullAcceptanceTPC() const {return fCorrectedHadEtFullAcceptanceTPC;}
    Float_t GetCorrectedPiKPEtFullAcceptanceITS() const {return fCorrectedHadEtFullAcceptanceITS+fCorrectedHadEtFullAcceptanceTPC;}
    Float_t GetCorrectedPiKPEtFullAcceptanceTPCNoPID() const {return fCorrectedHadEtFullAcceptanceTPCNoPID;}
    Float_t GetCorrectedPiKPEtFullAcceptanceITSNoPID() const {return fCorrectedHadEtFullAcceptanceITSNoPID+fCorrectedHadEtFullAcceptanceTPCNoPID;}
    void SetCorrections(AliAnalysisHadEtCorrections *corr){fCorrections = corr;}
    AliAnalysisHadEtCorrections *GetCorrections(){return fCorrections;}

    void CreateHistograms();
     virtual void Init();
    
protected:

    Bool_t CheckGoodVertex(AliVParticle *track);
    AliAnalysisHadEtCorrections *fCorrections;//corrections needed for hadronic et

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
    Float_t fCorrectedHadEtFullAcceptanceTPCAssumingPion;//get the corrected hadronic et for full acceptance, pt>0.15 GeV/c
    Float_t fCorrectedHadEtFullAcceptanceITSAssumingPion;//get the corrected hadronic et for full acceptance, pt>0.10 GeV/c
    Float_t fCorrectedHadEtFullAcceptanceTPCAssumingProton;//get the corrected hadronic et for full acceptance, pt>0.15 GeV/c
    Float_t fCorrectedHadEtFullAcceptanceITSAssumingProton;//get the corrected hadronic et for full acceptance, pt>0.10 GeV/c
    Float_t fCorrectedHadEtFullAcceptanceTPCAssumingKaon;//get the corrected hadronic et for full acceptance, pt>0.15 GeV/c
    Float_t fCorrectedHadEtFullAcceptanceITSAssumingKaon;//get the corrected hadronic et for full acceptance, pt>0.10 GeV/c
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

    void AddEt(Float_t rawEt, Float_t rawEtNoPID, Float_t corrEt, Float_t corrEtPion, Float_t corrEtProton, Float_t corrEtKaon, Float_t corrEtNoPID, Float_t pt, Bool_t IsTPC, Bool_t InPHOS, Bool_t InEMCAL);
    Bool_t IsInPHOS(AliESDtrack *track);
    Bool_t IsInEMCAL(AliESDtrack *track);

    void ResetEventValues();
    ClassDef(AliAnalysisHadEtReconstructed, 1);
};

#endif // ALIANALYSISHADETRECONSTRUCTED_H
