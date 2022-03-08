#ifndef AliAnalysisTaskR2P2multClass_H_Included
#define AliAnalysisTaskR2P2multClass_H_Included

#include "AliAnalysisTaskSE.h"
#include "TString.h"
#include "AliLog.h"

#include "AliPID.h"
#include "AliPIDResponse.h"

class AliAODEvent;
class AliESDEvent;
class AliInputEventHandler;
class TH1;
class TH2;
class TH2;
class TH3;
class TH1F;
class TH2F;
class TH2F;
class TH3F;
class TH1D;
class TH2D;
class TH2D;
class TH3D;
class TProfile;
class AliAnalysisUtils;
class AliEventCuts;
class AliHelperPID;
class AliPPVsMultUtils;
class AliVEvent;
class AliAODTrack;

class AliAnalysisTaskR2P2multClass : public AliAnalysisTaskSE
{
public:
    AliAnalysisTaskR2P2multClass();
    AliAnalysisTaskR2P2multClass(const TString & name);
 
  //PID functions
  //User should call ONLY the function GetParticleSpecies and set the PID strategy in the steering macro!
  Int_t TellParticleSpecies( AliVTrack * trk );//calculate the PID according to the slected method. // for pt cut analysis
  Int_t TellParticleSpecies_CircularCut( AliVTrack * trk );
  Int_t TellParticleSpecies_by_P( AliVTrack * trk );//calculate the PID according to the slected method. // for p cut analysis
  Int_t TellParticleSpecies_by_P_CircularCut( AliVTrack * trk );
  void CalculateNSigmas( AliVTrack * trk );   //Calcuate nsigma[ipart][idet], fill NSigma histos
  void CalculateTPCNSigmasElectron( AliVTrack * trk );
  void CheckTOF( AliVTrack * trk );   //check the TOF matching and set fHasTOFPID
  Double_t TOFBetaCalculation( AliVTrack * track ) const;
  Double_t massSquareCalculation( AliVTrack * track ) const;
  Float_t TPC_EventPlane(AliAODEvent *event);
  Bool_t Is2015PileUpEvent();
  Bool_t StoreEventMultiplicities(AliVEvent *event);
  Double_t CalculateSharedFraction(const TBits *triggerClusterMap,const TBits *assocClusterMap,const TBits *triggerShareMap,const TBits *assocShareMap);

  Bool_t ResonanceV0Cut(const Int_t pdgMother,
			  const Double_t px1, const Double_t py1,
			  const Double_t pz1, const Double_t px2,
			  const Double_t py2, const Double_t pz2);
  Bool_t photonConvCut(const Double_t px1, const Double_t py1,
		       const Double_t pz1, const Double_t px2,
		       const Double_t py2, const Double_t pz2);

 private:
    Double_t fnsigmas[4][2]; //nsigma values
    Bool_t fHasTOFPID;
    Double_t fNSigmaPID; // number of sigma for PID cut
    Double_t fNSigmaPID_veto;
    Double_t ptUpperLimit; //pt cut upper limit
    Double_t ptTOFlowerBoundary; // pt value which is the boundary between TPC & TOF.
    Double_t electronNSigmaVeto;
    Bool_t fRemoveTracksT0Fill;//if true remove tracks for which only StartTime from To-Fill is available (worst resolution)
    Double_t fSharedfraction_Pair_cut;
    
    AliAnalysisUtils *fUtils; //!
    AliPPVsMultUtils *fPPVsMultUtils; //!
    AliEventCuts *   fEventCut;  //!
    
    AliAnalysisTaskR2P2multClass(const  AliAnalysisTaskR2P2multClass&);
    const AliAnalysisTaskR2P2multClass& operator=(const  AliAnalysisTaskR2P2multClass&);
    
public:
    virtual ~AliAnalysisTaskR2P2multClass();
    
    // Implementation of interace methods
    //virtual void   ConnectInputData(Option_t *);
    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   FinishTaskOutput();
    virtual void   Terminate(Option_t* );
    virtual void   createHistograms();
    virtual void   finalizeHistograms();
    
    virtual void   addToList(TH1 *h);
    
    TH1D * createHisto1D(const TString &  name, const TString &  title,Int_t n, Double_t xmin, Double_t xmax,const TString &  xTitle, const TString &  yTitle);
    TH1D * createHisto1D(const TString &  name, const TString &  title,Int_t n, Double_t * bins,const TString &  xTitle, const TString &  yTitle);
    TH2D * createHisto2D(const TString &  name, const TString &  title,
                         Int_t nx, Double_t xmin, Double_t xmax, Int_t ny, Double_t ymin, Double_t ymax, const TString &  xTitle, const TString &  yTitle, const TString &  zTitle);
    TH2D * createHisto2D(const TString &  name, const TString &  title, Int_t nx, Double_t* xbins, Int_t ny, Double_t ymin, Double_t ymax,
                         const TString &  xTitle, const TString &  yTitle, const TString &  zTitle);
    
    TH1F * createHisto1F(const TString &  name, const TString &  title,Int_t n, Double_t xmin, Double_t xmax,const TString &  xTitle, const TString &  yTitle);
    TH1F * createHisto1F(const TString &  name, const TString &  title,Int_t n, Double_t * bins,const TString &  xTitle, const TString &  yTitle);
    TH2F * createHisto2F(const TString &  name, const TString &  title,
                         Int_t nx, Double_t xmin, Double_t xmax, Int_t ny, Double_t ymin, Double_t ymax, const TString &  xTitle, const TString &  yTitle, const TString &  zTitle);
    TH2F * createHisto2F(const TString &  name, const TString &  title, Int_t nx, Double_t* xbins, Int_t ny, Double_t ymin, Double_t ymax,
                         const TString &  xTitle, const TString &  yTitle, const TString &  zTitle);
    TH3F * createHisto3F(const TString &  name, const TString &  title,
                         Int_t nx, Double_t xmin, Double_t xmax, Int_t ny, Double_t ymin, Double_t ymax, Int_t nz, Double_t zmin, Double_t zmax,
                         const TString &  xTitle, const TString &  yTitle, const TString &  zTitle);
    TProfile * createProfile(const TString &  title,const TString &  description, Int_t n,  Double_t xMin,Double_t xMax,
                             const TString &  xTitle, const TString &  yTitle);
    TProfile * createProfile(const TString &  name,const TString &  description,
                             Int_t nx,  Double_t* bins,
                             const TString &  xTitle, const TString &  yTitle);
    
    //________________________________________________________________________
    
    Float_t  * getFloatArray(Int_t size, Float_t v);
    Double_t * getDoubleArray(Int_t size, Double_t v);
    void fillHistoWithArray(TH1 * h, Double_t * array, Int_t size);
    void fillHistoWithArray(TH2 * h, Double_t * array, Int_t size1, Int_t size2);
    void fillHistoWithArray(TH3 * h, Double_t * array, Int_t size1, Int_t size2, Int_t size3);
    void fillHistoWithArray(TH1 * h, Float_t * array, Int_t size);
    void fillHistoWithArray(TH2 * h, Float_t * array, Int_t size1, Int_t size2);
    void fillHistoWithArray(TH3 * h, Float_t * array, Int_t size1, Int_t size2, Int_t size3);
      
    virtual     void    SetDebugLevel( Int_t v )              { _debugLevel   = v; }
    virtual     void    SetSinglesOnly(Int_t v)               { _singlesOnly  = v; }
    virtual     void    SetPIDparticle( bool v )            { PIDparticle   = v; }
    virtual     void    SetUse_pT_cut( bool v )             { use_pT_cut   = v; }
    virtual     void    SetUse_AliHelperPID( bool v )       { useAliHelperPID   = v; }
    virtual     void    SetUse_CircularCutPID( bool v )     { useCircularCutPID = v; }
    virtual     void    SetIfContaminationInMC( bool v )    { NoContamination   = v; }
    virtual     void    SetIfContaminationWeakInMC( bool v )    { NoContaminationWeak   = v; }
    virtual     void    SetIfContaminationWeakMaterialInMC( bool v )    { NoContaminationWeakMaterial   = v; }
    virtual     void    SetIfMisIDWeakMaterialInMCClosure( bool v )     { Closure_NoMisIDWeakMaterial   = v; }
    virtual     void    SetUsePtEff(Int_t v)                { _usePtEff   = v; }
    virtual     void    SetUseWeights(Int_t v)                { _useWeights   = v; }
    virtual     void    SetUseRapidity(Int_t v)               { _useRapidity  = v; }
    virtual     void    SetEventPlane(bool v)               { _useEventPlane  = v; }
    virtual     void    SetEPmin( Double_t v)                 { EP_min          = v; }
    virtual     void    SetEPmax( Double_t v)                 { EP_max          = v; }
    virtual     void    SetSameFilter(Int_t v)                { _sameFilter   = v; }
    
    virtual     void    SetRejectPileup(Int_t v)              { _rejectPileup         = v; }
    virtual     void    SetRejectPairConversion(Int_t v)      { _rejectPairConversion = v; }
    virtual     void    SetVertexZMin(Double_t v)             { _vertexZMin           = v; }
    virtual     void    SetVertexZMax(Double_t v)             { _vertexZMax           = v; }
    virtual     void    SetVertexZWidth(Double_t v)           { _vertexZWidth         = v; }
    virtual     void    SetEtaWidth(Double_t v)               { _etaWidth             = v; }
    virtual     void    SetVertexXYMin(Double_t v)            { _vertexXYMin          = v; }
    virtual     void    SetVertexXYMax(Double_t v)            { _vertexXYMax          = v; }
    virtual     void    SetCentralityMethod(Int_t v)          { _centralityMethod     = v; }
    virtual     void    SetCentrality(Double_t centralityMin, Double_t centralityMax)
    {
        _centralityMin = centralityMin;
        _centralityMax = centralityMax;
    }


    virtual     void    SetRequestedCharge_1(Int_t v)     { _requestedCharge_1 = v; }
    virtual     void    SetRequestedCharge_2(Int_t v)     { _requestedCharge_2 = v; }
    virtual     void    SetPtMin1( Double_t v)            { _min_pt_1          = v; }

    virtual     void    SetPtMax1( Double_t v)            { _max_pt_1          = v; }
    virtual     void    SetPtBinWidth1( Double_t v)       { _width_pt_1        = v; }
    virtual     void    SetNPhiBins1( Int_t v)            { _nBins_phi_1       = v; }
    virtual     void    SetEtaMin1(Double_t v)            { _min_eta_1         = v; } // SetYMin1 acturally 
    virtual     void    SetEtaMax1(Double_t v)            { _max_eta_1         = v; } // SetYMax1 acturally
    virtual     void    SetPtMin2( Double_t v)            { _min_pt_2          = v; }
    virtual     void    SetPtMax2( Double_t v)            { _max_pt_2          = v; }
    virtual     void    SetPtBinWidth2( Double_t v)       { _width_pt_2        = v; }
    virtual     void    SetNPhiBins2( Int_t v)            { _nBins_phi_2       = v; }
    virtual     void    SetEtaMin2(Double_t v)            { _min_eta_2         = v; } // SetYMin2 acturally
    virtual     void    SetEtaMax2(Double_t v)            { _max_eta_2         = v; } // SetYMax2 acturally
    virtual     void    SetDcaZMin(Double_t v)            { _dcaZMin           = v; }
    virtual     void    SetDcaZMax(Double_t v)            { _dcaZMax           = v; }
    virtual     void    SetDcaXYMin(Double_t v)           { _dcaXYMin          = v; }
    virtual     void    SetDcaXYMax(Double_t v)           { _dcaXYMax          = v; }
    virtual     void    SetTPCNclus(Int_t v)              { _tpcnclus          = v; }
    virtual     void    SetChi2PerNDF(Double_t v)         { _chi2ndf           = v; }
    
    virtual     void    SetDedxMin(Double_t v)            { _dedxMin           = v; }
    virtual     void    SetDedxMax(Double_t v)            { _dedxMax           = v; }
    virtual     void    SetNClusterMin(Int_t v)           { _nClusterMin       = v; }
    virtual     void    SetTrackFilterBit(Int_t v)        { _trackFilterBit    = v; }

    virtual     void    SetPtEff_1(TH1F * v)           { _hPtEff_1          = v; }
    virtual     void    SetPtEff_2(TH1F * v)           { _hPtEff_2          = v; }
    
    virtual     void    SetWeigth_1(TH3F * v)           { _weight_1          = v; }
    virtual     void    SetWeigth_2(TH3F * v)           { _weight_2          = v; }

    AliHelperPID                   * GetHelperPID()          { return fHelperPID; }
    void SetHelperPID(AliHelperPID* pid)                     { fHelperPID = pid;  }

    void SetParticleSpecies( Int_t species )            { particleSpecies = species; }

    void SetAnalysisType( const char * analysisType ) { fAnalysisType = analysisType; }
    void SetSystemType( const char * systemType )     { fSystemType = systemType; }
    void SetResonancesCut( Bool_t NoResonances )      {  fExcludeResonancesInMC = NoResonances; }
    void SetElectronCut( Bool_t NoElectron )          { fExcludeElectronsInMC = NoElectron; }

    void SetNSigmaCut( Double_t nsigma )             { fNSigmaPID = nsigma; }
    void SetNSigmaCut_veto( Double_t nsigma )        { fNSigmaPID_veto = nsigma; }
    void SetPtCutUpperLimit( Double_t ptUpper )      { ptUpperLimit = ptUpper; }
    void SetPtTOFlowerBoundary( Double_t ptTPCTOFboundary )   { ptTOFlowerBoundary = ptTPCTOFboundary; }
    void SetElectronNSigmaVetoCut( Double_t electronVeto )   { electronNSigmaVeto = electronVeto; }
    void SetfRemoveTracksT0Fill( bool tof )     { fRemoveTracksT0Fill = tof; }    //fRemoveTracksT0Fill
    //void SetAliEventCuts(AliEventCuts * Event_Cut)     { fEventCut = Event_Cut; }
    void SetSharedFractionPairCut( Double_t v )   { fSharedfraction_Pair_cut = v; }
    
protected:
    
    // Handlers and events
    AliAODEvent*             fAODEvent;             //! AOD Event
    AliESDEvent*             fESDEvent;             //! ESD Event
    AliInputEventHandler*    fInputHandler;    //! Generic InputEventHandler
    AliPIDResponse*          fPIDResponse;  //!
    AliHelperPID* fHelperPID;       // points to class for PID
    
    // Histogram settings
    //TList*              _inputHistoList;
    TList*              _outputHistoList;   //!
    //Int_t _outputSlot;
    TClonesArray*                   fMCArray=nullptr; //!
    
    Double_t   _twoPi;
    long     _eventCount;
    
    //configuration variables and filters
    Int_t      _debugLevel;
    Int_t      _singlesOnly;
    bool      PIDparticle;
    bool      use_pT_cut;
    bool      useAliHelperPID;
    bool      useCircularCutPID;
    bool      NoContamination;
    bool      NoContaminationWeak;
    bool      NoContaminationWeakMaterial;
    bool      Closure_NoMisIDWeakMaterial;
    Int_t      _usePtEff;
    Int_t      _useWeights;
    Int_t      _useRapidity;
    bool     _useEventPlane;
    Double_t   EP_min;
    Double_t   EP_max;
    TH1F  *  _psi_EventPlane;
    Int_t      _sameFilter;
    Int_t      _rejectPileup;
    Int_t      _rejectPairConversion;
    Double_t   _vertexZMin;
    Double_t   _vertexZMax;
    Double_t   _vertexZWidth;
    Double_t   _etaWidth;
    Double_t   _vertexXYMin;
    Double_t   _vertexXYMax;
    Int_t      _centralityMethod;
    Double_t   _centralityMin;
    Double_t   _centralityMax;
    Int_t      _requestedCharge_1;
    Int_t      _requestedCharge_2;
    Double_t   _dcaZMin;
    Double_t   _dcaZMax;
    Double_t   _dcaXYMin;
    Double_t   _dcaXYMax;
    Double_t   _dedxMin;
    Double_t   _dedxMax;
    Int_t      _nClusterMin;
    Int_t      _trackFilterBit;
    Double_t particleSpecies;

    TString      fAnalysisType;
    TString      fSystemType;

    Bool_t fExcludeResonancesInMC;
    Bool_t fExcludeElectronsInMC;

    TFormula *f2015V0MtoTrkTPCout;
    TFormula *f2015V0MtoTrkTPCout_Upper;
    TFormula *fV0MtoTrkTPCout_lower_pp13TeV;
    TFormula *fV0MtoTrkTPCout_upper_pp13TeV;
    Int_t fV0Multiplicity;
    Int_t fV0AMultiplicity;
    Int_t fV0CMultiplicity;
    Int_t fV0Multiplicity_Victor;
    Int_t fNoOfTPCoutTracks;
    
    //----------------------------------------------------------------
    Double_t  fMultV0A;                                  //  mult. V0A
    Double_t  fMultV0C;                                  //  mult. V0C
    Double_t  fMultV0M;                                  //  mult. V0A+V0C

    //-------------------------------------------------------------


    
    Int_t _tpcnclus;
    Double_t _chi2ndf;
    
    // event and track wise variables
    
    Double_t _field;
    Int_t    _nTracks;
    Int_t    _nTracksTruth;
    Int_t nTracksMC;
    Int_t _nTpcCls;
    Double_t _mult0;
    Double_t _mult1;
    Double_t _mult2;
    Double_t _mult3;
    Double_t _mult4;
    Double_t _mult4a;
    Double_t _mult5;
    Double_t _mult6;
    Double_t _mult7;
    Double_t _mult8;
    
    //particle 1
    Int_t     arraySize;
    Int_t    *_id_1;               //!
    Int_t    *_charge_1;           //!
    //Int_t  *  _iPhi_1;            //!
    //Int_t  *  _iEta_1;            //!
    Int_t    *_iEtaPhi_1;            //!
    Int_t    *_iPt_1;             //!
    Float_t  *_pt_1;               //!
    Float_t  *_px_1;              //!
    Float_t  *_py_1;              //!
    Float_t  *_pz_1;              //!
    //Float_t * _phi_1;             //!
    //Float_t*  _eta_1;             //!
    Float_t  *_correction_1;           //!
    Float_t  *_dedx_1;           //!
    AliAODTrack ** _TrackArray;  //!  

    //particle 2
    Int_t    *_id_2;              //!
    Int_t    *_charge_2;            //!
    //Int_t    *_iPhi_2;            //!
    //Int_t    *_iEta_2;            //!
    Int_t    *_iEtaPhi_2;            //!
    Int_t    *_iPt_2;             //!
    Float_t  *_pt_2;              //!
    Float_t  *_px_2;              //!
    Float_t  *_py_2;              //!
    Float_t  *_pz_2;              //!
    //Float_t  *_phi_2;             //!
    //Float_t  *_eta_2;             //!
    Float_t  *_correction_2;           //!
    Float_t  *_dedx_2;           //!
    
    Float_t * _correctionPtEff_1;           //!
    Float_t * _correctionPtEff_2;           //!

    Float_t * _correctionWeight_1;           //!
    Float_t * _correctionWeight_2;           //!
    
    //histograming
    Int_t _nBins_M0;       Double_t _min_M0;       Double_t _max_M0;       Double_t _width_M0;
    Int_t _nBins_M1;       Double_t _min_M1;       Double_t _max_M1;       Double_t _width_M1;
    Int_t _nBins_M2;       Double_t _min_M2;       Double_t _max_M2;       Double_t _width_M2;
    Int_t _nBins_M3;       Double_t _min_M3;       Double_t _max_M3;       Double_t _width_M3;
    Int_t _nBins_M4;       Double_t _min_M4;       Double_t _max_M4;       Double_t _width_M4;
    Int_t _nBins_M5;       Double_t _min_M5;       Double_t _max_M5;       Double_t _width_M5;
    Int_t _nBins_M6;       Double_t _min_M6;       Double_t _max_M6;       Double_t _width_M6;
    Int_t _nBins_M7;       Double_t _min_M7;       Double_t _max_M7;       Double_t _width_M7;
    Int_t _nBins_M8;       Double_t _min_M8;       Double_t _max_M8;       Double_t _width_M8;
    
    Int_t _nBins_vertexZ;  Double_t _min_vertexZ;  Double_t _max_vertexZ;  Double_t _width_vertexZ;

    Int_t _nBins_pt_1;     Double_t _min_pt_1;     Double_t _max_pt_1;     Double_t _width_pt_1;
    Int_t _nBins_phi_1;    Double_t _min_phi_1;    Double_t _max_phi_1;    Double_t _width_phi_1;
    Int_t _nBins_eta_1;    Double_t _min_eta_1;    Double_t _max_eta_1;    Double_t _width_eta_1;
    Int_t _nBins_etaPhi_1;
    Int_t _nBins_etaPhiPt_1;
    Int_t _nBins_zEtaPhiPt_1;
    
    Int_t _nBins_pt_2;     Double_t _min_pt_2;     Double_t _max_pt_2;     Double_t _width_pt_2;
    Int_t _nBins_phi_2;    Double_t _min_phi_2;    Double_t _max_phi_2;    Double_t _width_phi_2;
    Int_t _nBins_eta_2;    Double_t _min_eta_2;    Double_t _max_eta_2;    Double_t _width_eta_2;
    Int_t _nBins_etaPhi_2;
    Int_t _nBins_etaPhiPt_2;
    Int_t _nBins_zEtaPhiPt_2;
    
    Int_t _nBins_etaPhi_12;
    
    Double_t __n1_1;
    Double_t __n1_2;
    Double_t __n2_12;
    Double_t __s1pt_1;
    Double_t __s1pt_2;
    Double_t __s2ptpt_12;
    Double_t __s2NPt_12;
    Double_t __s2PtN_12;
    
    Double_t __n1Nw_1;
    Double_t __n1Nw_2;
    Double_t __n2Nw_12;
    Double_t __s1ptNw_1;
    Double_t __s1ptNw_2;
    Double_t __s2ptptNw_12;
    Double_t __s2NPtNw_12;
    Double_t __s2PtNNw_12;
    

    Double_t * __n1_1_vsPt;   //!
    Double_t * __n1Nw_1_vsPt;   //!
    Double_t * __n1_1_vsPt_pdg;   //!
    Double_t * __n1_1_vsPt_pdg_Weak;   //!
    Double_t * __n1_1_vsPt_pdg_Weak_Material;   //!
    Double_t * __n1_1_vsPt_Weak;   //!
    Double_t * __n1_1_vsPt_Material;   //!
    Double_t * __n1_1_vsEtaPhi;     //!
    Double_t * __n1Nw_1_vsEtaPhi;     //!
    Double_t * __s1pt_1_vsEtaPhi;    //!
    Float_t  * __n1_1_vsZEtaPhiPt;    //!
    Float_t  * __wt_1_vsEtaPhi;    //! 
    Double_t * __n1_1_vsEta;   //!
    Double_t * __n1Nw_1_vsEta;   //!
    Double_t * __n1_1_vsPhi;   //!
    Double_t * __n1Nw_1_vsPhi;   //!

    Double_t * __n1_2_vsPt;   //!
    Double_t * __n1Nw_2_vsPt;   //!
    Double_t * __n1_2_vsPt_pdg;   //!
    Double_t * __n1_2_vsPt_pdg_Weak;   //!
    Double_t * __n1_2_vsPt_pdg_Weak_Material;   //!
    Double_t * __n1_2_vsPt_Weak;   //!
    Double_t * __n1_2_vsPt_Material;   //!
    Double_t * __n1_2_vsEtaPhi;     //!
    Double_t * __n1Nw_2_vsEtaPhi;     //!
    Double_t * __s1pt_2_vsEtaPhi;    //!
    Float_t  * __n1_2_vsZEtaPhiPt;    //!
    Float_t  * __wt_2_vsEtaPhi;    //! 
    Double_t * __n1_2_vsEta;   //!
    Double_t * __n1Nw_2_vsEta;   //!
    Double_t * __n1_2_vsPhi;   //!
    Double_t * __n1Nw_2_vsPhi;   //!

    //Double_t * __n2_12_vsPtPt;
    //Double_t * __n2_12_vsEtaPhi;
    //Double_t * __s2ptpt_12_vsEtaPhi;
    //Double_t * __s2PtN_12_vsEtaPhi;
    //Double_t * __s2NPt_12_vsEtaPhi;
    
    Double_t * __n2_12_vsPtPt;   //!
    Float_t  * __n2_12_vsEtaPhi;   //!
    Float_t  * __n2Nw_12_vsEtaPhi;   //!
    Float_t  * __s2ptpt_12_vsEtaPhi;   //!
    Float_t  * __s2PtN_12_vsEtaPhi;   //!
    Float_t  * __s2NPt_12_vsEtaPhi;   //!
    
    TH1F * _hPtEff_1;
    TH1F * _hPtEff_2;
    TH3F * _weight_1;
    TH3F * _weight_2;
    //    TProfile * _hProfPileupCut;
    TH1D * _eventDetails;    
    TH1D * _trackDetails;   	
    TH1D * _m0;
    TH1D * _m1;
    TH1D * _m2;
    TH1D * _m2DiffMultBeforeCut;
    TH1D * _m2DiffMult;
    TH1D * _m2RatioMult;
    TH1D * _m2Ratio2Mult;
    TH2F *multDiffVsTruth;
    TH2F *multDiffNegVsTruth;
    TH2F *multRecoVsTruth;

    TH2F *multDiff2VsTruth;
    TH2F *multDiff2NegVsTruth;
    TH2F *multReco2VsTruth;
    TH1D * _m3;
    TH1D * _m4;
    TH1D * _m5;
    TH1D * _m6;
    TH1D * _m7;
    TH1D * _m8;
    TH1D * _vertexZ;
    TH1D * _vertexZHisto_before;
    TH1D * _vertexZHisto;
    TH1D** fV0M;//!array of histograms for V0M flatness check

	
    
    /*TH1F * _Ncluster1;
    TH1F *_Ncluster1_NoMisID; 
    TH1F *_Ncluster1_NoMisID_NoWeakDecayed; 
    TH1F *_Ncluster1_NoMisID_NoWeakDecayed_NoMaterial;

    
    TH1D *h1_chiSqNDF_beforeCut;
    TH1D *h1_chiSqNDF;
    TH1D *h1_chiSqNDF_NoMisID; 
    TH1D *h1_chiSqNDF_NoMisID_NoWeakDecayed; 
    TH1D *h1_chiSqNDF_NoMisID_NoWeakDecayed_NoMaterial;
    */

    TH1F * _t0_1d;
    TH1F * _trackLength;
    TH1F * _trackLength_GetIntegratedLength;
    TH1F * _timeTOF_1d;
    TH1F * _realTOF_1d;
    TH1F * _t0_1d_POI;
    TH1F * _trackLength_POI;
    TH1F * _trackLength_GetIntegratedLength_POI;
    TH1F * _timeTOF_1d_POI;
    TH1F * _realTOF_1d_POI;
    
    TH1F * _etadis_POI_AliHelperPID;
    TH1F * _etadis_before_any_cuts;

    TH1F * _ydis_POI_AliHelperPID;
    
    //TH3F * _vZ_y_Pt_POI_AliHelperPID;
    //TH3F * _vZ_y_eta_POI_AliHelperPID;

    TH2F * _y_Pt_AllCh_MCAODTruth;
    TH2F * _y_Pt_POI_MCAODTruth;
    
    TH1F * _phidis_POI_AliHelperPID;
    TH1F * _phidis_before_any_cuts;

    TH1F * _ptMcTruth_all;
    TH1F * _ptMcTruth_phyPrimary;
    TH1F * _ptMcTruth_phyPrimary_noSecFromWeakDecay;
    TH1F * _ptMcTruth_phyPrimary_noSecFromWeakDecay_noSecFromMaterial;
    TH1F * _ptMcTruth_phyPrimary_noSecFromWeakDecay_noSecFromMaterial_noResonance;
  /*  TH1F * _dcazPos;
    TH1F * _dcazImpact;
    TH1F *_dcazImpact_NoMisID; 
    TH1F *_dcazImpact_NoMisID_NoWeakDecayed;   
    TH1F *_dcazImpact_NoMisID_NoWeakDecayed_NoMaterial; 
    
    TH1F * _dcaxyPos;
    TH1F * _dcaxyImpact;
    TH1F * _dcaxyPtDept;
    TH1F *_dcaxyImpact_NoMisID;  
    TH1F *_dcaxyImpact_NoMisID_NoWeakDecayed; 
    TH1F *_dcaxyImpact_NoMisID_NoWeakDecayed_NoMaterial;
*/
    
    TH2F *  _dedx_p;
    TH2F *  _dedx_p_POI_AliHelperPID;
    TH2F *  _dedx_p_AliHelperPID_no_Undefined;
    
    TH2F *  _beta_p;
    TH2F *  _beta_p_POI_AliHelperPID;
    TH2F *  _beta_p_AliHelperPID_no_Undefined;
    TH2F *  _nSigmaTOF_p_pion_before;
    TH2F *  _nSigmaTOF_p_pion_after;
    TH2F *  _nSigmaTOF_p_kaon_before;
    TH2F *  _nSigmaTOF_p_kaon_after;
    TH2F *  _nSigmaTOF_p_proton_before;
    TH2F *  _nSigmaTOF_p_proton_after;
    TH2F *  _nSigmaTPC_nSigmaTOF_Pion_before;
    TH2F *  _nSigmaTPC_nSigmaTOF_Pion_after;
    TH2F *  _nSigmaTPC_nSigmaTOF_Kaon_before;
    TH2F *  _nSigmaTPC_nSigmaTOF_Kaon_after;
    TH2F *  _nSigmaTPC_nSigmaTOF_Proton_before;
    TH2F *  _nSigmaTPC_nSigmaTOF_Proton_after;
    
    TH2F *  _inverse_beta_p;
    TH2F *  _inverse_beta_p_POI_AliHelperPID;
    TH2F *  _inverse_beta_p_AliHelperPID_no_Undefined;
    
    TH2F *  _msquare_p;
    TH2F *  _msquare_p_POI_AliHelperPID;
    TH2F *  _msquare_p_AliHelperPID_no_Undefined;
    TH2F *  _fhV0MvsTracksTPCout_before;
    TH2F *  _fhV0MvsTracksTPCout_after;
    TProfile *  _profV0MvsTPCout;
    TH2F *  _fV0MmultVsSpdTracklet_before;
    TH2F *  _fV0MmultVsSpdTracklet_after;

  
   TH1F *  _fhMultV0M;
   TH1F *  _fhMultV0A;
   TH1F *  _fhMultV0C;

   
   //-------------------------------------------------------------------------
    // PARTICLE 1 (satisfies filter 1)
    // Primary filled quantities

    TH1F      *  _n1_1_vsPt;
    TH1F      *  _n1Nw_1_vsPt;
    TH1F      *  _n1_1_vsPt_pdg;
    TH1F      *  _n1_1_vsPt_pdg_Weak;
    TH1F      *  _n1_1_vsPt_pdg_Weak_Material;

    TH1F      *  _n1_1_vsDCAzPos;
    TH1F      *  _n1_1_vsDCAzImpact_beforeCut;
    TH1F      *  _n1_1_vsDCAzImpact;
    TH1F      *  _n1_1_vsDCAzImpact_pdg;
    TH1F      *  _n1_1_vsDCAzImpact_pdg_Weak;
    TH1F      *  _n1_1_vsDCAzImpact_pdg_Weak_Material;

    TH1F      *  _n1_1_vsDCAxyPos;
    TH1F      *  _n1_1_vsDCAxyPtDept;
    TH1F      *  _n1_1_vsDCAxyImpact_beforeCut;
    TH1F      *  _n1_1_vsDCAxyImpact;
    TH1F      *  _n1_1_vsDCAxyImpact_pdg;
    TH1F      *  _n1_1_vsDCAxyImpact_pdg_Weak;
    TH1F      *  _n1_1_vsDCAxyImpact_pdg_Weak_Material;

    //   TH1F      *  _n1_1_vsNcluster2;
    TH1F * _Ncluster2;
    TH1F * _trackId; 
    TH1F * _trackIdTPCout;
    
    TH1F      *  _n1_1_vsNcluster1_beforeCut;
    TH1F      *  _n1_1_vsNcluster1;
    TH1F      *  _n1_1_vsNcluster1_pdg;
    TH1F      *  _n1_1_vsNcluster1_pdg_Weak;
    TH1F      *  _n1_1_vsNcluster1_pdg_Weak_Material;

    TH1F      *  _n1_1_vsChiSqPerNDF_beforeCut;
    TH1F      *  _n1_1_vsChiSqPerNDF;
    TH1F      *  _n1_1_vsChiSqPerNDF_pdg;
    TH1F      *  _n1_1_vsChiSqPerNDF_pdg_Weak;
    TH1F      *  _n1_1_vsChiSqPerNDF_pdg_Weak_Material;

    TH1F      *  _n1_1_vsPt_Weak;
    TH1F      *  _n1_1_vsPt_Material;
    TH1F      *  _n1_1_vsEta;
    TH1F      *  _n1Nw_1_vsEta;
    TH1F      *  _n1_1_vsPhi;
    TH1F      *  _n1Nw_1_vsPhi;
    //    TH1F      * h1f_wt1_vsEtaPhi;
    TH2F      *  _n1_1_vsEtaVsPhi;
    TH2F      *  _n1Nw_1_vsEtaVsPhi;
    TH2F      *  _s1pt_1_vsEtaVsPhi;
    TH3F      *  _n1_1_vsZVsEtaVsPhiVsPt;
    TH1F      *  _wt_1_vsEtaVsPhi;
    TProfile *  _n1_1_vsM;  // w/ weight
    TProfile *  _s1pt_1_vsM;
    TProfile *  _n1Nw_1_vsM; // w/o weight
    TProfile *  _s1ptNw_1_vsM;
    TH2D      *  _dedxVsP_1;
    TH2D      *  _corrDedxVsP_1;
    TH2F      *  _betaVsP_1;
    
    // PARTICLE 2 (satisfies filter 2)
    // Primary filled quantities
    TH1F      *  _n1_2_vsPt;
    TH1F      *  _n1Nw_2_vsPt;
    TH1F      *  _n1_2_vsPt_pdg;
    TH1F      *  _n1_2_vsPt_pdg_Weak;
    TH1F      *  _n1_2_vsPt_pdg_Weak_Material;

    TH1F      *  _n1_2_vsDCAzPos;
    TH1F      *  _n1_2_vsDCAzImpact_beforeCut;
    TH1F      *  _n1_2_vsDCAzImpact;
    TH1F      *  _n1_2_vsDCAzImpact_pdg;
    TH1F      *  _n1_2_vsDCAzImpact_pdg_Weak;
    TH1F      *  _n1_2_vsDCAzImpact_pdg_Weak_Material;

    TH1F      *  _n1_2_vsDCAxyPos;
    TH1F      *  _n1_2_vsDCAxyPtDept;
    TH1F      *  _n1_2_vsDCAxyImpact_beforeCut;
    TH1F      *  _n1_2_vsDCAxyImpact;
    TH1F      *  _n1_2_vsDCAxyImpact_pdg;
    TH1F      *  _n1_2_vsDCAxyImpact_pdg_Weak;
    TH1F      *  _n1_2_vsDCAxyImpact_pdg_Weak_Material;

    //   TH1F      *  _n1_2_vsNcluster2;

    TH1F      *  _n1_2_vsNcluster1_beforeCut;
    TH1F      *  _n1_2_vsNcluster1;
    TH1F      *  _n1_2_vsNcluster1_pdg;
    TH1F      *  _n1_2_vsNcluster1_pdg_Weak;
    TH1F      *  _n1_2_vsNcluster1_pdg_Weak_Material;

    TH1F      *  _n1_2_vsChiSqPerNDF_beforeCut;
    TH1F      *  _n1_2_vsChiSqPerNDF;
    TH1F      *  _n1_2_vsChiSqPerNDF_pdg;
    TH1F      *  _n1_2_vsChiSqPerNDF_pdg_Weak;
    TH1F      *  _n1_2_vsChiSqPerNDF_pdg_Weak_Material;

    TH1F      *  _n1_2_vsPt_Weak;
    TH1F      *  _n1_2_vsPt_Material;
    TH1F      *  _n1_2_vsEta;
    TH1F      *  _n1Nw_2_vsEta;
    TH1F      *  _n1_2_vsPhi;
    TH1F      *  _n1Nw_2_vsPhi;
    //    TH1F      * h1f_wt2_vsEtaPhi;
    TH2F      *  _n1_2_vsEtaVsPhi;
    TH2F      *  _n1Nw_2_vsEtaVsPhi;
    TH2F      *  _s1pt_2_vsEtaVsPhi;
    TH3F      *  _n1_2_vsZVsEtaVsPhiVsPt;
    TH1F      *  _wt_2_vsEtaVsPhi;
    TProfile *  _n1_2_vsM;
    TProfile *  _s1pt_2_vsM;
    TProfile *  _n1Nw_2_vsM; // w/o weight
    TProfile *  _s1ptNw_2_vsM;
    TH2D      *  _dedxVsP_2;
    TH2D      *  _corrDedxVsP_2;
    TH2F      *  _betaVsP_2;
    
    // Pairs 1 & 2
    //    TH1F      * h1f_wt12_vsEtaPhi;
    TH1F      * _n2_12_vsEtaPhi;
    TH1F      * _n2Nw_12_vsEtaPhi;
    TH2F      * _n2_12_vsPtVsPt;
    TH1F      * _s2PtPt_12_vsEtaPhi;
    TH1F      * _s2PtN_12_vsEtaPhi;
    TH1F      * _s2NPt_12_vsEtaPhi;
    
    TProfile * _n2_12_vsM;
    TProfile * _s2PtPt_12_vsM;
    TProfile * _s2PtN_12_vsM;
    TProfile * _s2NPt_12_vsM;
    TProfile * _n2Nw_12_vsM;
    TProfile * _s2PtPtNw_12_vsM;
    TProfile * _s2PtNNw_12_vsM;
    TProfile * _s2NPtNw_12_vsM;
    
    TH1F * _invMassRho0Before   ; 
    TH1F * _invMassK0sBefore    ; 
    TH1F * _invMassLambdaBefore ;  
    TH1F * _invMassRho0After    ;
    TH1F * _invMassK0sAfter     ;  
    TH1F * _invMassLambdaAfter  ;     
    TH1F     * _invMassKaon;
    TH1F     * _invMassKaonSq;
    TH1F     * _invMassElec_beforeCut;
    TH1F     * _invMassElec;
    TH1F     * _ClusterSharedFraction_beforeCut;
    TH1F     * _ClusterSharedFraction_afterCut;
    TH1F     * _ClusterSharedFraction_3by3Bins_beforeCut;
    TH1F     * _ClusterSharedFraction_3by3Bins_afterCut;

    ///thnsparse added by baidya on 25jul19 for correlation of detectors                             
    //    THnSparseD *fCorrDet_beforeCut;//
    // THnSparseD *fCorrDet_afterCut;   // 

    TString n1Name;
    TString n1NwName;
    TString n2Name;
    TString n2NwName;
    TString n3Name;
    TString n1n1Name;
    TString n1n1n1Name;
    TString n2n1Name;
    TString r1Name;
    TString r2Name;
    TString r3Name;
    TString r2r1Name;
    TString c2Name;
    TString c3Name;
    TString d3Name;
    TString p3Name;
    TString cName;
    
    TString intR2Name;
    TString binCorrName;
    TString intBinCorrName;
    
    TString countsName;
    TString part_1_Name;
    TString part_2_Name;
    TString part_3_Name;
    TString pair_12_Name;
    TString pair_13_Name;
    TString pair_23_Name;
    TString tripletName;
    
    TString avg;
    TString avgName;
    TString sumName;
    TString s1ptName;
    TString s1ptNwName;
    TString s1DptName;
    
    TString s2PtPtName;
    TString s2NPtName;
    TString s2PtNName;
    TString s2DptDptName;
    
    TString s2PtPtNwName;
    TString s2NPtNwName;
    TString s2PtNNwName;
    
    TString ptName;
    TString ptptName;
    TString pt1pt1Name;
    TString DptName;
    TString DptDptName;
    TString RDptDptName;
    TString nPtName;
    TString ptNName;
    TString seanName;
    
    TString _title_counts;
    
    TString _title_m0;
    TString _title_m1;
    TString _title_m2;
    TString _title_m3;
    TString _title_m4;
    TString _title_m5;
    TString _title_m6;
    TString _title_m7;
    TString _title_m8;
    
    TString _title_eta_1;
    TString _title_phi_1;
    TString _title_pt_1;
    TString _title_etaPhi_1;
    TString _title_n_1;
    TString _title_SumPt_1;
    TString _title_AvgPt_1;
    TString _title_AvgN_1;
    TString _title_AvgSumPt_1;
    
    TString _title_eta_2;
    TString _title_phi_2;
    TString _title_pt_2;
    TString _title_etaPhi_2;
    TString _title_n_2;
    TString _title_SumPt_2;
    TString _title_AvgPt_2;
    TString _title_AvgN_2;
    TString _title_AvgSumPt_2;
    
    TString _title_etaPhi_12;
    
    TString _title_wt;
    TString _title_AvgN2_12;
    TString _title_AvgSumPtPt_12;
    TString _title_AvgSumPtN_12;
    TString _title_AvgNSumPt_12;
    
    TString vsZ;
    TString vsM;
    TString vsPt;
    TString vsPhi; 
    TString vsEta; 
    TString vsEtaPhi; 
    TString vsPtVsPt;
    TString pdg;
    TString Weak;
    TString Material;
    
    ClassDef(AliAnalysisTaskR2P2multClass,1)
}; 


#endif

