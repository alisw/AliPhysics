/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
//
//  Secondary vertexing construction Class
//  Construct secondary vertex from Beauty hadron with electron and
//  hadrons, then apply selection criteria
//

#ifndef ALIHFEELECBACKGROUND_H
#define ALIHFEELECBACKGROUND_H

#ifndef ROOT_TObject
//#include <TObject.h>
#endif

class AliESDEvent;
class AliESDpid;
class AliESDVertex;
class AliAODEvent;
class AliESDtrack;
class AliAODTrack;
class AliMCEvent;

class AliHFEpid;


//________________________________________________________________
class AliHFEelecbackground : public TObject {
  public: 
    AliHFEelecbackground();
    AliHFEelecbackground(const AliHFEelecbackground &p);
    AliHFEelecbackground &operator=(const AliHFEelecbackground &);
    virtual ~AliHFEelecbackground();
		virtual Bool_t Load(const Char_t *filename);
		virtual Bool_t Load(TList * const outputlist);

    void CreateHistograms(TList * const qaList);
		void Reset();

    Bool_t HasMCData() const { return TestBit(kHasMCData); };
    Bool_t IsAODanalysis() const { return TestBit(kAODanalysis); };
    Bool_t IsESDanalysis() const { return !TestBit(kAODanalysis); };

    void SetHasMCData(Bool_t hasMCdata = kTRUE) { SetBit(kHasMCData,hasMCdata); };
    void SetAODAnalysis() { SetBit(kAODanalysis, kTRUE); };
    void SetESDAnalysis() { SetBit(kAODanalysis, kFALSE); };
    void SetEvent(AliESDEvent* const ESD); 
    void SetEventAOD(AliAODEvent* const AOD){fAOD1=AOD;}; 
    void SetMCEvent(AliMCEvent* const mcEvent){fMCEvent=mcEvent;};  
    
    void SetOpeningAngleCut(Double_t openingAngleCut){fOpeningAngleCut = openingAngleCut;};
    void SetInvMassCut(Double_t invMassCut){fInvMassCut = invMassCut;};
    void SetChi2NdfCut(Double_t chi2NdfCut){fChi2NdfCut = chi2NdfCut;};
    void SetUseAliKFCode(Bool_t useAliKFCode){fUseAliKFCode = useAliKFCode;};    
    void SetSharedClusterCut(Bool_t sharedClusterCut){fSharedClusterCut = sharedClusterCut;};
    void SetRequireITSStandalone(Short_t requireITSStandalone){fRequireITSStandalone = requireITSStandalone;};
    void SetMinITSChi2(Double_t minITSChi2) {fMinITSChi2 = minITSChi2;};
    void SetMinNbCls(Int_t minNbCls){fMinNbCls = minNbCls;};
    void SetMinNbClsSDDSPD(Int_t minNbClsSDDSPD){fMinNbClsSDDSPD = minNbClsSDDSPD;};
    void SetPIDPartner();
    void SetPIDMethodPartner(AliHFEpid * const pid) {fPIDMethodPartner = pid;};
    void SetPIDMethodPartnerITS(AliESDpid * const pid) {fPIDMethodPartnerITS = pid;};
    void SetDebugLevel(Short_t debugLevel) { fDebugLevel = debugLevel;};
    
    Double_t GetOpeningAngleCut() const { return fOpeningAngleCut; };
    Double_t GetInvMassCut() const { return fInvMassCut; };
    Double_t GetChi2NdfCut() const { return fChi2NdfCut; };
    Bool_t GetUseAliKFCode() const { return fUseAliKFCode; };
    Bool_t GetSharedClusterCut() const { return fSharedClusterCut; };
    Short_t GetRequireITSStandalone() const { return fRequireITSStandalone; };
    Int_t GetMinNbCls() const { return fMinNbCls;};
    Double_t GetMinITSChi2() const { return fMinITSChi2; };
    Int_t GetMinNbClsSDDSPD() const { return fMinNbClsSDDSPD;};
    Bool_t GetPIDPartner() const { return fPIDPartner;};
    
    TList *GetList()  const           { return fList; };
    TList *GetListPostProcess() const { return fListPostProcess; };
    
    Bool_t SingleTrackCut(const AliESDtrack* const trackPart) const;
    Bool_t ShareCluster(AliESDtrack * const track1,AliESDtrack * const track2); 
    Bool_t PIDTrackCut(AliESDtrack* const trackPart);
    void PairAnalysis(AliESDtrack* const track, AliESDtrack* const trackpart); 
    void FillOutput(const Double_t *results,const Double_t *resultsr, Int_t sign); 
    void PostProcess();
    void Plot() const;
    
 private:
    enum{
      kHasMCData = BIT(15),             // bitset for mc data usage
	kAODanalysis = BIT(16)            // bitset for aod analysis
	};
    enum {kDatai=0, kDatar=1, kDatadca=2, kDatachi2Ndf=3, kMCo=4, kMCr=5, kMCdca=6, kMCchi2Ndf=7, kMCe=8, kMCcutPart0=9, kMCcutPart1=10, kMCcutPart2=11, kMCcutPart3=12};   // In the fList Data/MC
    enum {kOs=0, kPp=1, kNn=2, kR=3};            // In the last dimension Charge
    enum {kOos=0, kOss=1, kOr=2, kOdiff=3};      // outputs 
    enum {kNOutput=4,kNMCInfo=5};                // Nb of outputs
    enum{
      kElectronFromBackground = 0,
	kElectronFromGamma = 1,
	kElectronFromPi0 = 2,
	kElectronFromEta = 3,
	kElectronFromC = 4,
	kElectronFromB = 5		  
	};                                // MC: Origin 
    enum {kNotSplitted=0, kSplittedOs=1, kSplittedSs=2}; // MC: splitted 
    
    Bool_t CalculateMotherVariable(AliESDtrack* const track, AliESDtrack* const trackpart, Double_t *results);
    void CalculateMotherVariableR(AliESDtrack* const track, AliESDtrack* const trackpart, Double_t *results);
    Int_t IsMotherGamma(Int_t tr);
    Int_t IsMotherPi0(Int_t tr);
    Int_t IsMotherEta(Int_t tr);
    Int_t IsMotherC(Int_t tr);
    Int_t IsMotherB(Int_t tr);
    Int_t GetPdg(Int_t tr);
    Int_t GetLabMother(Int_t tr);
    
    static Double_t BetheBlochElectronITS(const Double_t *x, const Double_t * /*par*/);
    static Double_t BetheBlochMuonITS(const Double_t *x, const Double_t * /*par*/);
    static Double_t BetheBlochPionITS(const Double_t *x, const Double_t * /*par*/);
    static Double_t BetheBlochKaonITS(const Double_t *x, const Double_t * /*par*/);
    static Double_t BetheBlochProtonITS(const Double_t *x, const Double_t * /*par*/); 
    
    static Double_t BetheBlochElectronTPC(const Double_t *x, const Double_t * /*par*/);
    static Double_t BetheBlochMuonTPC(const Double_t *x, const Double_t * /*par*/);
    static Double_t BetheBlochPionTPC(const Double_t *x, const Double_t * /*par*/);
    static Double_t BetheBlochKaonTPC(const Double_t *x, const Double_t * /*par*/);
    static Double_t BetheBlochProtonTPC(const Double_t *x, const Double_t * /*par*/); 
    
    THnSparseF *fhtmp;     // Only to avoid coverity problem
    TH2F *fhtmpf;          // Only to avoid coverity problem
    TH1F *fhtmpp;          // Only to avoid coverity problem
    
    AliESDEvent* fESD1;              //! ESD pointer             
    AliAODEvent* fAOD1;              //! AOD pointer             
    AliMCEvent*  fMCEvent;           //! MC event             
    Double_t fBz;                    // Magnetic field 
    const AliESDVertex *fkVertex;    //! Primary vertex
    static const Double_t fgkMe;     //!  Mass of the electron
    Double_t fPtESD;                 //! pt of tagged electron
    Int_t fIndexTrack;               //! index track
    Int_t fPdg;                      //! pdg code track 
    Int_t fLabMother;                //! label first mother track 
    Int_t fIsFrom;                   //! is track from
    Int_t fMotherGamma;              //! Gamma, mother of track
    Int_t fMotherPi0;                //! Pi0, mother of track
    Int_t fMotherC;                  //! C, mother of track
    Int_t fMotherB;                  //! B, mother of track
    Int_t fMotherEta;                //! eta, mother of track
    Bool_t fIsPartner;               //! Are partners
    Bool_t fIsSplittedTrack;         //! Are splitted track
    
    Double_t fOpeningAngleCut;       //! Opening angle cut
    Double_t fInvMassCut;            //! Inv mass cut
    Double_t fChi2NdfCut;            //! chi2ndf cut for KF code
    
    Bool_t   fUseAliKFCode;          //! Use AliKF code to calculate the pair properties
    
    Bool_t  fSharedClusterCut;       //! Shared Cluster Cut
    Short_t fRequireITSStandalone;   //! ITS standalone: 1 and 2 (pureITSStandalone)
    Int_t  fMinNbCls;                //! Min Nb of clusters ITS or TPC
    Double_t fMinITSChi2;            //! ITS chi2 min
    Int_t  fMinNbClsSDDSPD;          //! Min Nb of clusters ITS SDD&SPD
    Bool_t fPIDPartner;              //! PID partner
    AliHFEpid *fPIDMethodPartner;    //! PID cuts
    AliESDpid *fPIDMethodPartnerITS; //! PID cuts ITS
    
    Short_t fDebugLevel;             //! Debug Level
    
    Bool_t   fCuts[10];               //! Cut passed already
    
    TList *fList;                    //! list for outputs
    TList *fListPostProcess;         //! list for postprocess
    
    static Bool_t  fgUseMCPID;       // flag to use MC PID for tagged electron
    
    ClassDef(AliHFEelecbackground,0);
};

#endif
