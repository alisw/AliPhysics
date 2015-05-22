#ifndef AliAnalysisTaskThermalGAFlow_h
#define AliAnalysisTaskThermalGAFlow_h

#include "AliAnalysisTaskSE.h"
#include <vector>
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

class AliAnalysisTaskThermalGAFlow : public AliAnalysisTaskSE {
  public:
    AliAnalysisTaskThermalGAFlow( const char *name = "AliAnalysisTaskThermalGAFlow");
    virtual ~AliAnalysisTaskThermalGAFlow();
    virtual void UserCreateOutputObjects();
    virtual void UserExec( Option_t *);
    virtual void Terminate( Option_t *);

//
    virtual void SetDebug(Int_t x) {fDebug = x;}
    virtual void SetMinCells(Int_t x) {fMinCells = x;}
    virtual void SetMinE(Double_t x) {fMinE = x;}
    virtual void SetMinTrackDr(Double_t x) {fMinTrackDr = x;}
    virtual void SetMaxVertexx(Double_t x) {fMaxVertexx = x;}
    virtual void SetMinCentrality(Double_t x) {fMinCentrality = x;}
    virtual void SetMaxCentrality(Double_t x) {fMaxCentrality = x;}
    virtual void SetCoreRadius(Double_t x) {fCoreRadius = x;}
    virtual void SetMinCoreEnergyRatio(Double_t x) {fMinCoreEnergyRatio = x;}

//

  protected:
    Int_t fRunNumber; //Run number of present run

    TList             *fAnalist;
    std::vector< std::vector<TObject*> >  fCompass;
    AliPHOSGeoUtils   *fPHOSgeomU;
    AliPHOSGeometry   *fPHOSgeom;
    AliVEvent         *fEvent; //Present Event
    AliESDEvent       *fesd; //fEvent, when esd
    AliAODEvent       *faod; //fEvent, when aod
    Double_t          fvertexx[3]; //The position of the event vertex in the global coordinate system
    Float_t           fcentrality;
    TF1               *fNonLinearCorr;

    Int_t fDebug;
    Int_t fMinCells;
    Double_t fMinE;
    Double_t fMinTrackDr;
    Double_t fMaxVertexx;
    Double_t fMinCentrality;
    Double_t fMaxCentrality;
    Double_t fCoreRadius;
    Double_t fMinCoreEnergyRatio;


    void FillHist(const char * key, Double_t x, Double_t y) const;
    void FillHist(const char * key, Double_t x) const;
    void FillHist(const char * key, Double_t x, Double_t y, Double_t z) const;
    void ConfigureEvent();
    void InitializeGeometry();
    void ScanBasicEventParameters();
    void CaptainsLog(Int_t s);
    void ScanClusters();
    void SavePhoton(AliESDCaloCluster* cluster);
    void ScanMod(Int_t* x, TLorentzVector p);
    Bool_t ApplyEventCuts();
    Bool_t ApplyClusterCuts(AliESDCaloCluster *cluster);
    Bool_t ApplyClusterCuts(AliAODCaloCluster *cluster);
    Bool_t ApplyShapeCuts(AliESDCaloCluster *cluster);
    Bool_t ApplyCPCuts(AliESDCaloCluster *cluster);
    Int_t *GetPos(AliESDCaloCluster* cluster, Int_t x[4]);
    Int_t *GetPos(AliAODCaloCluster* cluster, Int_t x[4]);

//    class histcache{
//      public:
//        void set_map (const char key, TObject*)
//        TObject* get_hist (const char) {return ;}
//      private:
//        const char hist_name;
//        TObject* hist_pointer;
//    } 
    

//*****************Add Task BS*****************//   
//We need to declare task to be private to avoid compilation warning
AliAnalysisTaskThermalGAFlow(const AliAnalysisTaskThermalGAFlow & g) ; //Copy Constructor
AliAnalysisTaskThermalGAFlow& operator=(const AliAnalysisTaskThermalGAFlow&); //Mystery Line

ClassDef(AliAnalysisTaskThermalGAFlow, 4);
};
#endif
