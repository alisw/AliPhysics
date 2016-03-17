#ifndef AliAnalysisTaskThermalGAFlow_h
#define AliAnalysisTaskThermalGAFlow_h

#include "AliAnalysisTaskSE.h"
//#include <vector>
class TNtupleD;
class TList;
class TH1F;
class TH2F;
class TH3F;
class TString;

class AliPHOSGeoUtils;
class AliPHOSGeometry;

class AliAODCaloCluster;
class AliESDCaloCluster;
#include <AliAODCaloCluster.h>
#include <AliESDCaloCluster.h>
#include <AliCaloPhoton.h>
#include <AliEventplane.h>

class AliAnalysisTaskThermalGAFlow : public AliAnalysisTaskSE {
  public:
    AliAnalysisTaskThermalGAFlow( const char *name = "AliAnalysisTaskThermalGAFlow");
    virtual ~AliAnalysisTaskThermalGAFlow();
    virtual void UserCreateOutputObjects();
    virtual void UserExec( Option_t *);
    virtual void Terminate( Option_t *);

    virtual void SetDebug(Int_t x) {fDebug = x;}
    virtual void SetMinCells(Int_t x) {fMinCells = x;}
    virtual void SetMinE(Double_t x) {fMinE = x;}
    virtual void SetMinTrackDr(Double_t x) {fMinTrackDr = x;}
    virtual void SetMaxVertexx(Double_t x) {fMaxVertexx = x;}
    virtual void SetMinCentrality(Double_t x) {fMinCentrality = x;}
    virtual void SetMaxCentrality(Double_t x) {fMaxCentrality = x;}
    virtual void SetCoreRadius(Double_t x) {fCoreRadius = x;}
    virtual void SetMinCoreEnergyRatio(Double_t x) {fMinCoreEnergyRatio = x;}
    virtual void SetMinLambdaDisp(Double_t x) {fMinLambdaDisp = x;}
    virtual void SetMinCPVStd(Double_t x) {fMinCPVStd = x;}

    virtual void SetMixVertxbins(Int_t x) {fMixVertxbins = x;}
    virtual void SetMixCentbins(Int_t x) {fMixCentbins = x;}
    virtual void SetMixEvbins(Int_t x) {fMixEvbins = x;}
    virtual void SetNptbins(Int_t x) {fNptbins = x;}

  protected:

    Int_t fRunNumber; //Run number of present run

    TList             *fAnalist; //The Analist - stores output data
    TList             *fAtlas; //The Event List for photon mixing.  Use an atlas to find where you are
    TObjArray         *fCompass; //Multidimensional Array of Event Lists for mixing.  Use a compass to find where you are going.
    TObjArray         *fSextant; //Array of Photons in an Event.  Use a Sextant to look at stars.
    AliCaloPhoton     *fphoton;  //Primary considered photon
    AliCaloPhoton     *fmixphoton; //A second considered photon used during mixing
    AliPHOSGeoUtils   *fPHOSgeomU; //Utilities to manipulate the PHOS geometry
    AliPHOSGeometry   *fPHOSgeom;  //Parameterization of the PHOS geometry
    AliVEvent         *fEvent; //Present Event
    AliESDEvent       *fesd; //fEvent, when esd
    AliAODEvent       *faod; //fEvent, when aod
    Double_t          fvertexx[3]; //The position of the event vertex in the global coordinate system
    Float_t           fcentrality; //Centrality information
    Double_t          fevPlane; //Event Plane Inforamtion -- A and C averaged
    Double_t          fevPlaneA; //Event Plane Inforamtion -- A only
    Double_t          fevPlaneC; //Event Plane Inforamtion -- C only
    Double_t          fevScaler; //The V0A and V0C difference squared.  Used for QA.
    TF1               *fNonLinearCorr; //Nonlinear energy correction for core energy correction

    Int_t fDebug;
    Int_t fMinCells;
    Double_t fMinE;
    Double_t fMinTrackDr;
    Double_t fMaxVertexx;
    Double_t fMinCentrality;
    Double_t fMaxCentrality;
    Double_t fCoreRadius;
    Double_t fMinCoreEnergyRatio;
    Double_t fMinLambdaDisp;
    Double_t fMinCPVStd;

    Int_t fMixVertxbins;
    Int_t fMixCentbins;
    Int_t fMixEvbins;
    Int_t fNptbins;

//    void FillHist(const char * key, Double_t x, Double_t y) const;
//    void FillHist(const char * key, Double_t x) const;
//    void FillHist(const char * key, Double_t x, Double_t y, Double_t z) const;
    void ConfigureEvent();
    void InitializeGeometry();
    void ScanBasicEventParameters();
    void CaptainsLog(Int_t s);
    void ScanClusters();
    void SavePhoton(AliESDCaloCluster* cluster);
    void SaveEvent();
    void MesonExclusion();
    Double_t PionMass(AliCaloPhoton* y1, AliCaloPhoton* y2);
    Double_t PionPt(AliCaloPhoton* y1, AliCaloPhoton* y2);
    TList *ReferenceAtlas(Double_t vertx, Double_t cent, Double_t evplane);
    Bool_t ApplyEventCuts();
    Bool_t ApplyClusterCuts(AliESDCaloCluster *cluster);
    Bool_t ApplyClusterCuts(AliAODCaloCluster *cluster);
    Bool_t ApplyCoreShapeCut(AliESDCaloCluster *cluster);
    Bool_t ApplyDispMatrixCut(AliESDCaloCluster *cluster);
    Bool_t ApplyCPCuts(AliESDCaloCluster *cluster);

    Int_t *GetPos(AliESDCaloCluster* cluster, Int_t x[4]);
    Int_t *GetPos(AliAODCaloCluster* cluster, Int_t x[4]);
    void Hypnotic();
    void IncFlow();
    Bool_t WeakCuts(AliESDCaloCluster *cluster);
    void FillHist(const char * key, Double_t x, Double_t y) const;
    void FillHist(const char * key, Double_t x) const;
    void FillHist(const char * key, Double_t x, Double_t y, Double_t z) const;

//*****************Add Task BS*****************//   
//We need to declare task to be private to avoid compilation warning
AliAnalysisTaskThermalGAFlow(const AliAnalysisTaskThermalGAFlow & g) ; //Copy Constructor
AliAnalysisTaskThermalGAFlow& operator=(const AliAnalysisTaskThermalGAFlow&); //Mystery Line

ClassDef(AliAnalysisTaskThermalGAFlow, 4);
};
#endif
