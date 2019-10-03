#ifndef AliAnalysisTask3PCorrelations_H_Included
#define AliAnalysisTask3PCorrelations_H_Included

#include "AliAnalysisTaskSE.h"
#include "TString.h"
#include "AliLog.h"

class AliAODEvent;
class AliESDEvent;
class AliInputEventHandler;
//class AliMCEvent;
//class AliMCEventHandler;
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

///This class is designed for 3 particle analysis.
///There are 3 sets of cuts 
/// #1 : 'trigger'
/// #2 : low pt associates
/// #3 : low pt associates
class AliAnalysisTask3PCorrelations : public AliAnalysisTaskSE
{
public:
  AliAnalysisTask3PCorrelations();  
  AliAnalysisTask3PCorrelations(const TString & name);
  
private:
  AliAnalysisTask3PCorrelations(const  AliAnalysisTask3PCorrelations&);
  const AliAnalysisTask3PCorrelations& operator=(const  AliAnalysisTask3PCorrelations&);

public:
  virtual ~AliAnalysisTask3PCorrelations();

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
  virtual     void    SetUseWeights(int v)                { _useWeights   = v; } 
  
  virtual     void    SetRejectPileup(int v)              { _rejectPileup         = v; } 
  virtual     void    SetRejectPairConversion(int v)      { _rejectPairConversion = v; } 
  virtual     void    SetVertexZMin(double v)             { _vertexZMin           = v; } 
  virtual     void    SetVertexZMax(double v)             { _vertexZMax           = v; } 
  virtual     void    SetVertexXYMin(double v)            { _vertexXYMin          = v; } 
  virtual     void    SetVertexXYMax(double v)            { _vertexXYMax          = v; } 
  virtual     void    SetCentralityMethod(int v)          { _centralityMethod     = v; }
  virtual     void    SetCentrality(double centralityMin, double centralityMax)
  {
    _centralityMin = centralityMin;
    _centralityMax = centralityMax;
  }
  virtual     void    SetPtMin1( double v)            { _min_pt_1          = v; } 
  virtual     void    SetPtMax1( double v)            { _max_pt_1          = v; } 
  virtual     void    SetEtaMin1(double v)            { _min_eta_1         = v; } 
  virtual     void    SetEtaMax1(double v)            { _max_eta_1         = v; } 
  virtual     void    SetPtMin2( double v)            { _min_pt_2          = v; } 
  virtual     void    SetPtMax2( double v)            { _max_pt_2          = v; } 
  virtual     void    SetEtaMin2(double v)            { _min_eta_2         = v; } 
  virtual     void    SetEtaMax2(double v)            { _max_eta_2         = v; } 
  virtual     void    SetPtMin3( double v)            { _min_pt_3          = v; } 
  virtual     void    SetPtMax3( double v)            { _max_pt_3          = v; } 
  virtual     void    SetEtaMin3(double v)            { _min_eta_3         = v; } 
  virtual     void    SetEtaMax3(double v)            { _max_eta_3         = v; } 
  virtual     void    SetDcaZMin(double v)            { _dcaZMin           = v; } 
  virtual     void    SetDcaZMax(double v)            { _dcaZMax           = v; } 
  virtual     void    SetDcaXYMin(double v)           { _dcaXYMin          = v; } 
  virtual     void    SetDcaXYMax(double v)           { _dcaXYMax          = v; } 
  virtual     void    SetDedxMin(double v)            { _dedxMin           = v; } 
  virtual     void    SetDedxMax(double v)            { _dedxMax           = v; } 
  virtual     void    SetNClusterMin(int v)           { _nClusterMin       = v; } 
  virtual     void    SetTrackFilterBit(int v)        { _trackFilterBit    = v; }
  virtual     void    SetWeigth_1P(TH3F * v)          { _weight_1P          = v; }
  virtual     void    SetWeigth_1M(TH3F * v)          { _weight_1M          = v; }
  virtual     void    SetWeigth_2P(TH3F * v)          { _weight_2P          = v; }
  virtual     void    SetWeigth_2M(TH3F * v)          { _weight_2M          = v; }
  virtual     void    SetWeigth_3P(TH3F * v)          { _weight_3P          = v; }
  virtual     void    SetWeigth_3M(TH3F * v)          { _weight_3M          = v; }
  
  
protected:      
  
  // Handlers and events
  AliAODEvent*             fAODEvent;             //! AOD Event 
  AliESDEvent*             fESDEvent;             //! ESD Event 
  AliInputEventHandler*    fInputHandler;    //! Generic InputEventHandler 
  
  // Histogram settings
  //TList*              _inputHistoList;
  TList*              _outputHistoList;
  //int _outputSlot;
  
  
  double   _twoPi;
  long     _eventCount;
  
  //configuration variables and filters
  int      _debugLevel;
  int      _singlesOnly; 
  int      _useWeights; 
  int      _rejectPileup; 
  int      _rejectPairConversion; 
  double   _vertexZMin; 
  double   _vertexZMax; 
  double   _vertexXYMin; 
  double   _vertexXYMax; 
  int      _centralityMethod;
  double   _centralityMin;
  double   _centralityMax;
  double   _dcaZMin; 
  double   _dcaZMax; 
  double   _dcaXYMin; 
  double   _dcaXYMax;
  double   _dedxMin;
  double   _dedxMax;
  int      _nClusterMin; 
  int      _trackFilterBit;
  
 
  // event and track wise variables
  
  double _field;
  int    _nTracks;
  double _mult0;
  double _mult1;
  double _mult2;
  double _mult3;
  double _mult4;
  double _mult5;
  double _mult6;
  
  //particle 1
  int     arraySize;
  
  int    *_id_1;              //!
  int    *_charge_1;          //! 
  int    *_iPhi_1;            //!
  float  *_pt_1;              //!
  float  *_px_1;              //!
  float  *_py_1;              //!
  float  *_pz_1;              //!
  float  *_correction_1;      //!
  float  *_dedx_1;            //!
  
  //particle 2
  int    *_id_2;              //!
  int    *_charge_2;          //!
  int    *_iPhi_2;            //!
  float  *_pt_2;              //!
  float  *_px_2;              //!
  float  *_py_2;              //!
  float  *_pz_2;              //!
  float  *_correction_2;      //!
  float  *_dedx_2;            //!
  
  //particle 3
  int    *_id_3;              //!
  int    *_charge_3;          //!
  int    *_iPhi_3;            //!
  float  *_pt_3;              //!
  float  *_px_3;              //!
  float  *_py_3;              //!
  float  *_pz_3;              //!
  float  *_correction_3;      //!
  float  *_dedx_3;            //!
  
  float * _correctionWeight_1P;           //!
  float * _correctionWeight_1M;           //!
  float * _correctionWeight_2P;           //!
  float * _correctionWeight_2M;           //!
  float * _correctionWeight_3P;           //!
  float * _correctionWeight_3M;           //!
  
  //histograming
  int _nBins_M0;       double _min_M0;       double _max_M0;       double _width_M0;
  int _nBins_M1;       double _min_M1;       double _max_M1;       double _width_M1;
  int _nBins_M2;       double _min_M2;       double _max_M2;       double _width_M2;
  int _nBins_M3;       double _min_M3;       double _max_M3;       double _width_M3;
  int _nBins_M4;       double _min_M4;       double _max_M4;       double _width_M4;
  int _nBins_M5;       double _min_M5;       double _max_M5;       double _width_M5;
  int _nBins_M6;       double _min_M6;       double _max_M6;       double _width_M6;
  int _nBins_vsM;      double _min_vsM;      double _max_vsM;      double _width_vsM;
  
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
  
  int _nBins_pt_3;     double _min_pt_3;     double _max_pt_3;     double _width_pt_3;
  int _nBins_phi_3;    double _min_phi_3;    double _max_phi_3;    double _width_phi_3;
  int _nBins_eta_3;    double _min_eta_3;    double _max_eta_3;    double _width_eta_3;
  int _nBins_etaPhi_3; 
  int _nBins_etaPhiPt_3;
  int _nBins_zEtaPhiPt_3;

  int _nBins_phi_12;
  int _nBins_phi_13;
  int _nBins_phi_23;
  int _nBins_phi_123;
  
  int _nBins_etaPhi_12;
  int _nBins_etaPhi_13;
  int _nBins_etaPhi_23;
  int _nBins_etaPhi_123;
    
  double __n1_1;
  double __n1_2;
  double __n1_3;
  double __n2_12;   
  double __n2_13;   
  double __n2_23;   
  double __n3_123;   
  
  double __n1Nw_1;
  double __n1Nw_2;
  double __n1Nw_3;
  double __n2Nw_12;   
  double __n2Nw_13;   
  double __n2Nw_23;   
  
  double __n3Nw_123;   
  
  double * __n1_1_vsPt;   //!
  double * __n1_1_vsPhi;     //! 
  float  * __n1_1P_vsZEtaPhiPt;    //!
  float  * __n1_1M_vsZEtaPhiPt;    //!
    
  double * __n1_2_vsPt;   //!
  double * __n1_2_vsPhi;     //! 
  float  * __n1_2P_vsZEtaPhiPt;    //!
  float  * __n1_2M_vsZEtaPhiPt;    //!
  
  double * __n1_3_vsPt;   //!
  double * __n1_3_vsPhi;     //! 
  float  * __n1_3P_vsZEtaPhiPt;    //!
  float  * __n1_3M_vsZEtaPhiPt;    //!
  
  float  * __n2_12_vsPhi;   //!
  float  * __n2_13_vsPhi;   //!
  float  * __n2_23_vsPhi;   //!
  float  * __n3_123_vsPhi;   //!
  
  TH3F * _weight_1P;
  TH3F * _weight_1M;
  TH3F * _weight_2P;
  TH3F * _weight_2M;
  TH3F * _weight_3P;
  TH3F * _weight_3M;
  TH1D * _eventAccounting;
  TH1D * _m0;
  TH1D * _m1;
  TH1D * _m2;
  TH1D * _m3;
  TH1D * _m4;
  TH1D * _m5;
  TH1D * _m6;
  TH1D * _vertexZ;
  
  // PARTICLE 1 (satisfies filter 1)
  // Primary filled quantities
  TH1F      *  _n1_1_vsPhi;
  TProfile  *  _n1_1_vsM;  // w/ weight
  TProfile  *  _n1Nw_1_vsM; // w/o weight

  TH1F      *  _n1_1_vsPt;         
  TH3F      *  _n1_1P_vsZVsEtaVsPhiVsPt;
  TH3F      *  _n1_1M_vsZVsEtaVsPhiVsPt;
  TH2F      *  _dedxVsP_1;
  TH2F      *  _corrDedxVsP_1;
  TH2F      *  _betaVsP_1;
  
  // PARTICLE 2 (satisfies filter 2)
  // Primary filled quantities
  TH1F      *  _n1_2_vsPhi;
  TProfile  *  _n1_2_vsM;  // w/ weight
  TProfile  *  _n1Nw_2_vsM; // w/o weight

  TH1F      *  _n1_2_vsPt;         
  TH3F      *  _n1_2P_vsZVsEtaVsPhiVsPt;
  TH3F      *  _n1_2M_vsZVsEtaVsPhiVsPt;
  TH2F      *  _dedxVsP_2;
  TH2F      *  _corrDedxVsP_2;
  TH2F      *  _betaVsP_2;

  // PARTICLE 3 (satisfies filter 3)
  // Primary filled quantities
  TH1F      *  _n1_3_vsPhi;
  TProfile  *  _n1_3_vsM;  // w/ weight
  TProfile  *  _n1Nw_3_vsM; // w/o weight

  TH1F      *  _n1_3_vsPt;         
  TH3F      *  _n1_3P_vsZVsEtaVsPhiVsPt;
  TH3F      *  _n1_3M_vsZVsEtaVsPhiVsPt;
  TH2F      *  _dedxVsP_3;
  TH2F      *  _corrDedxVsP_3;
  TH2F      *  _betaVsP_3;

  // Pairs 1 & 2  
  TH2F      * _n2_12_vsPhi;        
  TH2F      * _n2_13_vsPhi;        
  TH2F      * _n2_23_vsPhi;        

  // Triplets
  TH3F      * _n3_123_vsPhi;        
  
  TProfile * _n2_12_vsM;        
  TProfile * _n2_13_vsM;        
  TProfile * _n2_23_vsM;        
  TProfile * _n3_123_vsM;        
  TProfile * _n2Nw_12_vsM;        
  TProfile * _n2Nw_13_vsM;        
  TProfile * _n2Nw_23_vsM;        
  TProfile * _n3Nw_123_vsM;        
  
  TString n1Name;
  TString n2Name;
  TString n3Name;
  TString n1NwName;
  TString n2NwName;
  TString n3NwName;
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

  TString intR2Name;
  TString binCorrName;
  TString intBinCorrName;
  
  TString countsName;
  TString part_1_Name;
  TString part_1P_Name;
  TString part_1M_Name;
  TString part_2_Name;
  TString part_2P_Name;
  TString part_2M_Name;
  TString part_3_Name;
  TString part_3P_Name;
  TString part_3M_Name;
  TString pair_12_Name;
  TString pair_13_Name;
  TString pair_23_Name;
  TString triplet_Name;
  
  TString _title_counts;
  TString _title_m0;
  TString _title_m1;
  TString _title_m2;
  TString _title_m3;
  TString _title_m4;
  TString _title_m5;
  TString _title_m6;
  
  TString _title_phi_1;
  TString _title_phi_2;
  TString _title_phi_3;
  TString _title_phi_12;
  TString _title_phi_13;
  TString _title_phi_23;
  TString _title_phi_123;
   
  TString vsZ;
  TString vsM;
  TString vsPt;
  TString vsPhi; 
  TString vsEta; 
  TString vsEtaPhi; 
  
  ClassDef(AliAnalysisTask3PCorrelations,1)
}; 


#endif


