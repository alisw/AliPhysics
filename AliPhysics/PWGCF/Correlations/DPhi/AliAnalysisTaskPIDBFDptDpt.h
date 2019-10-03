#ifndef AliAnalysisTaskPIDBFDptDpt_H_Included
#define AliAnalysisTaskPIDBFDptDpt_H_Included

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

class AliAnalysisTaskPIDBFDptDpt : public AliAnalysisTaskSE
{
public:
    AliAnalysisTaskPIDBFDptDpt();
    AliAnalysisTaskPIDBFDptDpt(const TString & name);
 
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
    AliEventCuts *   fEventCut;  //!
    
    AliAnalysisTaskPIDBFDptDpt(const  AliAnalysisTaskPIDBFDptDpt&);
    const AliAnalysisTaskPIDBFDptDpt& operator=(const  AliAnalysisTaskPIDBFDptDpt&);
    
public:
    virtual ~AliAnalysisTaskPIDBFDptDpt();
    
    // Implementation of interace methods
    //virtual void   ConnectInputData(Option_t *);
    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   FinishTaskOutput();
    virtual void   Terminate(Option_t* );
    virtual void   createHistograms();
    virtual void   finalizeHistograms();
    
    virtual void   addToList(TH1 *h);
    
    TH1D * createHisto1D(const TString &  name, const TString &  title,int n, double xmin, double xmax,const TString &  xTitle, const TString &  yTitle);
    TH1D * createHisto1D(const TString &  name, const TString &  title,int n, double * bins,const TString &  xTitle, const TString &  yTitle);
    TH2D * createHisto2D(const TString &  name, const TString &  title,
                         int nx, double xmin, double xmax, int ny, double ymin, double ymax, const TString &  xTitle, const TString &  yTitle, const TString &  zTitle);
    TH2D * createHisto2D(const TString &  name, const TString &  title, int nx, double* xbins, int ny, double ymin, double ymax,
                         const TString &  xTitle, const TString &  yTitle, const TString &  zTitle);
    
    TH1F * createHisto1F(const TString &  name, const TString &  title,int n, double xmin, double xmax,const TString &  xTitle, const TString &  yTitle);
    TH1F * createHisto1F(const TString &  name, const TString &  title,int n, double * bins,const TString &  xTitle, const TString &  yTitle);
    TH2F * createHisto2F(const TString &  name, const TString &  title,
                         int nx, double xmin, double xmax, int ny, double ymin, double ymax, const TString &  xTitle, const TString &  yTitle, const TString &  zTitle);
    TH2F * createHisto2F(const TString &  name, const TString &  title, int nx, double* xbins, int ny, double ymin, double ymax,
                         const TString &  xTitle, const TString &  yTitle, const TString &  zTitle);
    TH3F * createHisto3F(const TString &  name, const TString &  title,
                         int nx, double xmin, double xmax, int ny, double ymin, double ymax, int nz, double zmin, double zmax,
                         const TString &  xTitle, const TString &  yTitle, const TString &  zTitle);
    TProfile * createProfile(const TString &  title,const TString &  description, int n,  double xMin,double xMax,
                             const TString &  xTitle, const TString &  yTitle);
    TProfile * createProfile(const TString &  name,const TString &  description,
                             int nx,  double* bins,
                             const TString &  xTitle, const TString &  yTitle);
    
    //________________________________________________________________________
    
    float  * getFloatArray(int size, float v);
    double * getDoubleArray(int size, double v);
    void fillHistoWithArray(TH1 * h, double * array, int size);
    void fillHistoWithArray(TH2 * h, double * array, int size1, int size2);
    void fillHistoWithArray(TH3 * h, double * array, int size1, int size2, int size3);
    void fillHistoWithArray(TH1 * h, float * array, int size);
    void fillHistoWithArray(TH2 * h, float * array, int size1, int size2);
    void fillHistoWithArray(TH3 * h, float * array, int size1, int size2, int size3);
      
    virtual     void    SetDebugLevel( int v )              { _debugLevel   = v; }
    virtual     void    SetSinglesOnly(int v)               { _singlesOnly  = v; }
    virtual     void    SetPIDparticle( bool v )            { PIDparticle   = v; }
    virtual     void    SetUse_pT_cut( bool v )             { use_pT_cut   = v; }
    virtual     void    SetUse_AliHelperPID( bool v )       { useAliHelperPID   = v; }
    virtual     void    SetUse_CircularCutPID( bool v )     { useCircularCutPID = v; }
    virtual     void    SetIfContaminationInMC( bool v )    { NoContamination   = v; }
    virtual     void    SetIfContaminationWeakInMC( bool v )    { NoContaminationWeak   = v; }
    virtual     void    SetIfContaminationWeakMaterialInMC( bool v )    { NoContaminationWeakMaterial   = v; }
    virtual     void    SetIfMisIDWeakMaterialInMCClosure( bool v )     { Closure_NoMisIDWeakMaterial   = v; }
    virtual     void    SetUseWeights(int v)                { _useWeights   = v; }
    virtual     void    SetUseRapidity(int v)               { _useRapidity  = v; }
    virtual     void    SetEventPlane(bool v)               { _useEventPlane  = v; }
    virtual     void    SetEPmin( double v)                 { EP_min          = v; }
    virtual     void    SetEPmax( double v)                 { EP_max          = v; }
    virtual     void    SetSameFilter(int v)                { _sameFilter   = v; }
    
    virtual     void    SetRejectPileup(int v)              { _rejectPileup         = v; }
    virtual     void    SetRejectPairConversion(int v)      { _rejectPairConversion = v; }
    virtual     void    SetVertexZMin(double v)             { _vertexZMin           = v; }
    virtual     void    SetVertexZMax(double v)             { _vertexZMax           = v; }
    virtual     void    SetVertexZWidth(double v)           { _vertexZWidth         = v; }
    virtual     void    SetEtaWidth(double v)               { _etaWidth             = v; }
    virtual     void    SetVertexXYMin(double v)            { _vertexXYMin          = v; }
    virtual     void    SetVertexXYMax(double v)            { _vertexXYMax          = v; }
    virtual     void    SetCentralityMethod(int v)          { _centralityMethod     = v; }
    virtual     void    SetCentrality(double centralityMin, double centralityMax)
    {
        _centralityMin = centralityMin;
        _centralityMax = centralityMax;
    }
    
    virtual     void    SetRequestedCharge_1(int v)     { _requestedCharge_1 = v; }
    virtual     void    SetRequestedCharge_2(int v)     { _requestedCharge_2 = v; }
    virtual     void    SetPtMin1( double v)            { _min_pt_1          = v; }
    virtual     void    SetPtMax1( double v)            { _max_pt_1          = v; }
    virtual     void    SetPtBinWidth1( double v)       { _width_pt_1        = v; }
    virtual     void    SetNPhiBins1( int v)            { _nBins_phi_1       = v; }
    virtual     void    SetEtaMin1(double v)            { _min_eta_1         = v; } // SetYMin1 acturally 
    virtual     void    SetEtaMax1(double v)            { _max_eta_1         = v; } // SetYMax1 acturally
    virtual     void    SetPtMin2( double v)            { _min_pt_2          = v; }
    virtual     void    SetPtMax2( double v)            { _max_pt_2          = v; }
    virtual     void    SetPtBinWidth2( double v)       { _width_pt_2        = v; }
    virtual     void    SetNPhiBins2( int v)            { _nBins_phi_2       = v; }
    virtual     void    SetEtaMin2(double v)            { _min_eta_2         = v; } // SetYMin2 acturally
    virtual     void    SetEtaMax2(double v)            { _max_eta_2         = v; } // SetYMax2 acturally
    virtual     void    SetDcaZMin(double v)            { _dcaZMin           = v; }
    virtual     void    SetDcaZMax(double v)            { _dcaZMax           = v; }
    virtual     void    SetDcaXYMin(double v)           { _dcaXYMin          = v; }
    virtual     void    SetDcaXYMax(double v)           { _dcaXYMax          = v; }
    virtual     void    SetTPCNclus(int v)              { _tpcnclus          = v; }
    virtual     void    SetChi2PerNDF(double v)         { _chi2ndf           = v; }
    
    virtual     void    SetDedxMin(double v)            { _dedxMin           = v; }
    virtual     void    SetDedxMax(double v)            { _dedxMax           = v; }
    virtual     void    SetNClusterMin(int v)           { _nClusterMin       = v; }
    virtual     void    SetTrackFilterBit(int v)        { _trackFilterBit    = v; }
    virtual     void    SetWeigth_1(TH3F * v)           { _weight_1          = v; }
    virtual     void    SetWeigth_2(TH3F * v)           { _weight_2          = v; }

    AliHelperPID                   * GetHelperPID()          { return fHelperPID; }
    void SetHelperPID(AliHelperPID* pid)                     { fHelperPID = pid;  }

    void SetParticleSpecies( int species )            { particleSpecies = species; }

    void SetAnalysisType( const char * analysisType ) { fAnalysisType = analysisType; }
    void SetSystemType( const char * systemType )     { fSystemType = systemType; }
    void SetResonancesCut( Bool_t NoResonances )      { fExcludeResonancesInMC = NoResonances; }
    void SetElectronCut( Bool_t NoElectron )          { fExcludeElectronsInMC = NoElectron; }

    void SetNSigmaCut( double nsigma )             { fNSigmaPID = nsigma; }
    void SetNSigmaCut_veto( double nsigma )        { fNSigmaPID_veto = nsigma; }
    void SetPtCutUpperLimit( double ptUpper )      { ptUpperLimit = ptUpper; }
    void SetPtTOFlowerBoundary( double ptTPCTOFboundary )   { ptTOFlowerBoundary = ptTPCTOFboundary; }
    void SetElectronNSigmaVetoCut( double electronVeto )   { electronNSigmaVeto = electronVeto; }
    void SetfRemoveTracksT0Fill( bool tof )     { fRemoveTracksT0Fill = tof; }    //fRemoveTracksT0Fill
    //void SetAliEventCuts(AliEventCuts * Event_Cut)     { fEventCut = Event_Cut; }
    void SetSharedFractionPairCut( double v )   { fSharedfraction_Pair_cut = v; }
    
protected:
    
    // Handlers and events
    AliAODEvent*             fAODEvent;             //! AOD Event
    AliESDEvent*             fESDEvent;             //! ESD Event
    AliInputEventHandler*    fInputHandler;    //! Generic InputEventHandler
    
    AliPIDResponse*          fPIDResponse; //!
    AliHelperPID* fHelperPID;       // points to class for PID
    
    // Histogram settings
    //TList*              _inputHistoList;
    TList*              _outputHistoList;   //!
    //int _outputSlot;
    
    
    double   _twoPi;
    long     _eventCount;
    
    //configuration variables and filters
    int      _debugLevel;
    int      _singlesOnly;
    bool      PIDparticle;
    bool      use_pT_cut;
    bool      useAliHelperPID;
    bool      useCircularCutPID;
    bool      NoContamination;
    bool      NoContaminationWeak;
    bool      NoContaminationWeakMaterial;
    bool      Closure_NoMisIDWeakMaterial;
    int      _useWeights;
    int      _useRapidity;
    bool     _useEventPlane;
    double   EP_min;
    double   EP_max;
    TH1F  *  _psi_EventPlane;
    int      _sameFilter;
    int      _rejectPileup;
    int      _rejectPairConversion;
    double   _vertexZMin;
    double   _vertexZMax;
    double   _vertexZWidth;
    double   _etaWidth;
    double   _vertexXYMin;
    double   _vertexXYMax;
    int      _centralityMethod;
    double   _centralityMin;
    double   _centralityMax;
    int      _requestedCharge_1;
    int      _requestedCharge_2;
    double   _dcaZMin;
    double   _dcaZMax;
    double   _dcaXYMin;
    double   _dcaXYMax;
    double   _dedxMin;
    double   _dedxMax;
    int      _nClusterMin;
    int      _trackFilterBit;
    Double_t particleSpecies;

    TString      fAnalysisType;
    TString      fSystemType;

    Bool_t fExcludeResonancesInMC;
    Bool_t fExcludeElectronsInMC;

    TFormula *f2015V0MtoTrkTPCout;
    TFormula *f2015V0MtoTrkTPCout_Upper;
    Int_t fV0Multiplicity;
    Int_t fV0Multiplicity_Victor;
    Int_t fNoOfTPCoutTracks;
    
    int _tpcnclus;
    double _chi2ndf;
    
    // event and track wise variables
    
    double _field;
    int    _nTracks;
    double _mult0;
    double _mult1;
    double _mult2;
    double _mult3;
    double _mult4;
    double _mult4a;
    double _mult5;
    double _mult6;
    double _mult7;
    double _mult8;
    
    //particle 1
    int     arraySize;
    int    *_id_1;               //!
    int    *_charge_1;           //!
    //int  *  _iPhi_1;            //!
    //int  *  _iEta_1;            //!
    int    *_iEtaPhi_1;            //!
    int    *_iPt_1;             //!
    float  *_pt_1;               //!
    float  *_px_1;              //!
    float  *_py_1;              //!
    float  *_pz_1;              //!
    //float * _phi_1;             //!
    //float*  _eta_1;             //!
    float  *_correction_1;           //!
    float  *_dedx_1;           //!
    AliAODTrack ** _TrackArray;  //!
    
    //particle 2
    int    *_id_2;              //!
    int    *_charge_2;            //!
    //int    *_iPhi_2;            //!
    //int    *_iEta_2;            //!
    int    *_iEtaPhi_2;            //!
    int    *_iPt_2;             //!
    float  *_pt_2;              //!
    float  *_px_2;              //!
    float  *_py_2;              //!
    float  *_pz_2;              //!
    //float  *_phi_2;             //!
    //float  *_eta_2;             //!
    float  *_correction_2;           //!
    float  *_dedx_2;           //!
    
    float * _correctionWeight_1;           //!
    float * _correctionWeight_2;           //!
    
    //histograming
    int _nBins_M0;       double _min_M0;       double _max_M0;       double _width_M0;
    int _nBins_M1;       double _min_M1;       double _max_M1;       double _width_M1;
    int _nBins_M2;       double _min_M2;       double _max_M2;       double _width_M2;
    int _nBins_M3;       double _min_M3;       double _max_M3;       double _width_M3;
    int _nBins_M4;       double _min_M4;       double _max_M4;       double _width_M4;
    int _nBins_M5;       double _min_M5;       double _max_M5;       double _width_M5;
    int _nBins_M6;       double _min_M6;       double _max_M6;       double _width_M6;
    int _nBins_M7;       double _min_M7;       double _max_M7;       double _width_M7;
    int _nBins_M8;       double _min_M8;       double _max_M8;       double _width_M8;
    
    int _nBins_vertexZ;  double _min_vertexZ;  double _max_vertexZ;  double _width_vertexZ;
    
    int _nBins_pt_1;     double _min_pt_1;     double _max_pt_1;     double _width_pt_1;
    int _nBins_phi_1;    double _min_phi_1;    double _max_phi_1;    double _width_phi_1;
    int _nBins_eta_1;    double _min_eta_1;    double _max_eta_1;    double _width_eta_1;
    int _nBins_etaPhi_1;
    int _nBins_etaPhiPt_1;
    int _nBins_zEtaPhiPt_1;
    
    int _nBins_pt_2;     double _min_pt_2;     double _max_pt_2;     double _width_pt_2;
    int _nBins_phi_2;    double _min_phi_2;    double _max_phi_2;    double _width_phi_2;
    int _nBins_eta_2;    double _min_eta_2;    double _max_eta_2;    double _width_eta_2;
    int _nBins_etaPhi_2;
    int _nBins_etaPhiPt_2;
    int _nBins_zEtaPhiPt_2;
    
    int _nBins_etaPhi_12;
    
    double __n1_1;
    double __n1_2;
    double __n2_12;
    double __s1pt_1;
    double __s1pt_2;
    double __s2ptpt_12;
    double __s2NPt_12;
    double __s2PtN_12;
    
    double __n1Nw_1;
    double __n1Nw_2;
    double __n2Nw_12;
    double __s1ptNw_1;
    double __s1ptNw_2;
    double __s2ptptNw_12;
    double __s2NPtNw_12;
    double __s2PtNNw_12;
    
    double * __n1_1_vsPt;   //!
    double * __n1_1_vsPt_pdg;   //!
    double * __n1_1_vsPt_pdg_Weak;   //!
    double * __n1_1_vsPt_pdg_Weak_Material;   //!
    double * __n1_1_vsPt_Weak;   //!
    double * __n1_1_vsPt_Material;   //!
    double * __n1_1_vsEtaPhi;     //!
    double * __s1pt_1_vsEtaPhi;    //!
    float  * __n1_1_vsZEtaPhiPt;    //!
    
    double * __n1_2_vsPt;   //!
    double * __n1_2_vsPt_pdg;   //!
    double * __n1_2_vsPt_pdg_Weak;   //!
    double * __n1_2_vsPt_pdg_Weak_Material;   //!
    double * __n1_2_vsPt_Weak;   //!
    double * __n1_2_vsPt_Material;   //!
    double * __n1_2_vsEtaPhi;     //!
    double * __s1pt_2_vsEtaPhi;    //!
    float  * __n1_2_vsZEtaPhiPt;    //!
    
    //double * __n2_12_vsPtPt;
    //double * __n2_12_vsEtaPhi;
    //double * __s2ptpt_12_vsEtaPhi;
    //double * __s2PtN_12_vsEtaPhi;
    //double * __s2NPt_12_vsEtaPhi;
    
    double * __n2_12_vsPtPt;   //!
    float  * __n2_12_vsEtaPhi;   //!
    float  * __s2ptpt_12_vsEtaPhi;   //!
    float  * __s2PtN_12_vsEtaPhi;   //!
    float  * __s2NPt_12_vsEtaPhi;   //!
    
    TH3F * _weight_1;
    TH3F * _weight_2;
    TH1D * _eventAccounting;
    TH1D * _m0;
    TH1D * _m1;
    TH1D * _m2;
    TH1D * _m3;
    TH1D * _m4;
    TH1D * _m5;
    TH1D * _m6;
    TH1D * _m7;
    TH1D * _m8;
    TH1D * _vertexZ;
    
    TH1F * _Ncluster1;
    TH1F * _Ncluster2;

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
    
    TH1F * _dcaz;
    TH1F * _dcaxy;

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
    TH2F *  _fhV0MvsTracksTPCout_after;
    
    // PARTICLE 1 (satisfies filter 1)
    // Primary filled quantities
    TH1F      *  _n1_1_vsPt;
    TH1F      *  _n1_1_vsPt_pdg;
    TH1F      *  _n1_1_vsPt_pdg_Weak;
    TH1F      *  _n1_1_vsPt_pdg_Weak_Material;
    TH1F      *  _n1_1_vsPt_Weak;
    TH1F      *  _n1_1_vsPt_Material;
    TH2F      *  _n1_1_vsEtaVsPhi;
    TH2F      *  _s1pt_1_vsEtaVsPhi;
    TH3F      *  _n1_1_vsZVsEtaVsPhiVsPt;
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
    TH1F      *  _n1_2_vsPt_pdg;
    TH1F      *  _n1_2_vsPt_pdg_Weak;
    TH1F      *  _n1_2_vsPt_pdg_Weak_Material;
    TH1F      *  _n1_2_vsPt_Weak;
    TH1F      *  _n1_2_vsPt_Material;
    TH2F      *  _n1_2_vsEtaVsPhi;
    TH2F      *  _s1pt_2_vsEtaVsPhi;
    TH3F      *  _n1_2_vsZVsEtaVsPhiVsPt;
    TProfile *  _n1_2_vsM;
    TProfile *  _s1pt_2_vsM;
    TProfile *  _n1Nw_2_vsM; // w/o weight
    TProfile *  _s1ptNw_2_vsM;
    TH2D      *  _dedxVsP_2;
    TH2D      *  _corrDedxVsP_2;
    TH2F      *  _betaVsP_2;
    
    // Pairs 1 & 2
    TH1F      * _n2_12_vsEtaPhi;
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
    
    TH1F     * _invMassKaon;
    TH1F     * _invMassKaonSq;
    TH1F     * _invMassElec;
    TH1F     * _ClusterSharedFraction_beforeCut;
    TH1F     * _ClusterSharedFraction_afterCut;
    TH1F     * _ClusterSharedFraction_3by3Bins_beforeCut;
    TH1F     * _ClusterSharedFraction_3by3Bins_afterCut;
    
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
    
    ClassDef(AliAnalysisTaskPIDBFDptDpt,1)
}; 


#endif
