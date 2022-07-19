
//Done
#ifndef ALIANALYSISMEANPTDATA_H
#define ALIANALYSISMEANPTDATA_H

class TH1D;
class TH2F;
class TH3F;
class TProfile;
class TString;
class AliAODEvent;
class AliAODTrack;
class AliAODMCParticle;
class TList;
class AliESDtrackCuts;




#include "AliAnalysisTaskSE.h"
#include "THnSparse.h"
#include "TProfile.h"
#include "AliEventCuts.h"
#include "TProfile2D.h"

class AliAnalysisMeanPtdata : public AliAnalysisTaskSE {
 public:
    AliAnalysisMeanPtdata();
    AliAnalysisMeanPtdata(const char *name);
    Bool_t StoreEventMultiplicities(AliVEvent *event);
    virtual ~AliAnalysisMeanPtdata();
    
    virtual void     UserCreateOutputObjects();
    virtual void     UserExec(Option_t *option);
    virtual void     Terminate(Option_t *);

    //    virtual     void    SetUsePileupCut_PbPb5TeV( int v ){ _usePileupCut_PbPb5TeV = v;}
    virtual     void    SetPileupCut_PbPb5TeV( TProfile * v ){ _hProfPileupCut = v;}
    virtual     void    SetNClusterMin(int v){ _nClusterMin       = v; }
    virtual     void    SetFilterBit(int v)               { _fb  = v; }
    virtual     void    SetSinglesOnly(int v)               { _singlesOnly  = v; }
    virtual     void    SetDcaZMin(double v)            { _dcaZMin           = v; }
    virtual     void    SetDcaZMax(double v)            { _dcaZMax           = v; }
    virtual     void    SetDcaXYMin(double v)           { _dcaXYMin          = v; }
    virtual     void    SetDcaXYMax(double v)           { _dcaXYMax          = v; }
    virtual     void    SetChi2PerTPCCluster(double v)           { _chi2perTPC          = v; 
}
    virtual     void    SetChi2PerITSCluster(double v)           { _chi2perITS          = v; }

    virtual     void    SetVzMin(double v)            { _vZMin           = v; }
    virtual     void    SetVzMax(double v)            { _vZMax           = v; }
    virtual     void    SetNTPCCrossRows(double v)            { _nTPCCrossRows           = v; }
    virtual     void    SetPileUpEvent(double v)            { _pileUpEvent           = v; }




    AliEventCuts* GetAliEventCuts() const { return fAliEventCuts;}
    //    AliEventCuts fEventCuts;
    Int_t eventcount;
    void SetNContributors(Int_t nCont) {fNContributors = nCont;}
    void SetMaxVertexZDiff1(Float_t vZDiff1) {fMaxVertexZDiff1 = vZDiff1;}	
 private:
    TList           *fOutput;        // Output list
    AliPIDResponse  *fPIDResponse;	 // PID
	
    AliAODEvent     *fAOD;                   //! AOD event                        
    AliAODVertex    *fPrimaryVtx;            //! AOD vertex
    TH1D        *fEvents;                //
    TH1D        *fVtxZ;                // Vertex Z distribution 
    TH1F       *fVtxZCut;// Vertex Z dist. after vertex Z cut
    TH1F       *fVtxZCont;// Vertex Z dist. after vertex cut on nContributors
    TH1F       *fVtxZCutDiff;// Vertex Z dist. after vertex cut on vtx Z Difference
    TH1F*   fVtxZDiff1;// Difference 1 between vertex Z distributions
    TH1F*   fVtxZDiff2;// Difference 2 between vertex Z distributions
    TH1F*   fVtxZDiff3;// Difference 3 between vertex Z distributions
    TH1D        *fHistPt; //! Pt spectrum
    TH1D        *fHistEta;       // pseudorapidity spectrum
    TH1D        *fHistCentralityMultSelection; // centrality class selection by AliMultSelection

    TH1D        *fCount;
    TH1D        *fMultMeanQ;
    TH1D        *fMultTwopart;
    TH1D        *fMultThreepart;
    TH1D        *fMultFourpart;


    TH1D        *fMultQ1;
    TH1D        *fMultTwopartA;
    TH1D        *fMultThreepartA;
    TH1D        *fMultFourpartA;

    TH1D        *fMultA;
    TH1D        *fMultPairsA;
    TH1D        *fMultTripletsA;
    TH1D        *fMultQuadsA;
    TH1D        *fTpcNCrossedRows;
    TH1D        *fTpcNCluster;


    TProfile        *fMeanQ;
    TProfile        *fTwopart;
    TProfile        *fMeanQ10;
    TProfile        *fTwopart10;
    TProfile        *fMeanQ20;
    TProfile        *fTwopart20;
    TProfile        *fMeanQ50;
    TProfile        *fTwopart50;
    TProfile        *fMeanQ100;
    TProfile        *fTwopart100;
    TProfile        *fThreepart;
    TProfile        *fFourpart;

    TProfile        *fMeanQcent0p5;
    TProfile        *fTwopartcent0p5;
    TProfile        *fMeanQcent;
    TProfile        *fTwopartcent;
    TProfile        *fMeanQcent5;
    TProfile        *fTwopartcent5;
    //    TProfile        *_hProfPileupCut;
    TProfile        * _profV0MvsTPCout;
    TProfile2D            *fMeanQ_ss;
    TProfile2D        *fTwopart_ss;
    TProfile2D        *fThreepart_ss;
    TProfile2D        *fFourpart_ss;


    TProfile2D            *fMeanQcent_ss0p5;                                         
    TProfile2D        *fTwopartcent_ss0p5;                                                     TProfile2D            *fMeanQcent_ss;                                                     TProfile2D        *fTwopartcent_ss;                                                       TProfile2D            *fMeanQcent_ss5;                                                    TProfile2D        *fTwopartcent_ss5;




    TProfile        *fQ1A;
    TProfile        *fTwopartA;
    TProfile        *fThreepartA;
    TProfile        *fFourpartA;


    TProfile        *fA;
    TProfile        *fPairsA;
    TProfile        *fTripletsA;
    TProfile        *fQuadsA;

    TH1D    *fcentnpart;
    TH1D    *ftrack;
    TH2D    *fcentnpart2d;
    TH2D    *fcentnpart2d_1;
    
    TH2D    *h2d_tpc_ITSl1_before;
    TH2D    *h2d_tpc_ITSl2_before;
    TH2D    *h2d_ITSl1_ITSl2_before;

    TH2D    *h2d_tpc_ITSl1_after;
    TH2D    *h2d_tpc_ITSl2_after;
    TH2D    *h2d_ITSl1_ITSl2_after;

    TH2D    *_fhV0MvsTracksTPCout_before;
    TH2D    *_fhV0MvsTracksTPCout_after;

    TH2D    *_fV0MmultVsSpdTracklet_before;
    TH2D    *_fV0MmultVsSpdTracklet_after;


    TH1D        *fcentEvents;
    TH1D        *fcentEvents5;    
    TH1D        *fEventSee;    
    TH1F  *fHdcaxy;
    TH1F  *fHdcaz;


    Int_t        fNContributors;// Minimum contributors to the vertex   
    Float_t      fMaxVertexZDiff1;// Maximum value for Vertex Z difference TPC - global    
    TFormula *fV0MtoTrkTPCout_lower_PbPb5TeV;
    TFormula *fV0MtoTrkTPCout_upper_PbPb5TeV;
    //    int      _usePileupCut_PbPb5TeV;
    //    THnSparseD *fCorrDet_beforeCut;
    //    THnSparseD *fCorrDet_afterCut; 

    AliAnalysisMeanPtdata(const AliAnalysisMeanPtdata&); // not implemented
    AliAnalysisMeanPtdata& operator=(const AliAnalysisMeanPtdata&); // not implemented

 protected:


    //    AliAODEvent*  fAODEvent;
    AliEventCuts* fAliEventCuts; 
    TProfile        *_hProfPileupCut;
    Int_t fV0Multiplicity;
    Int_t fV0AMultiplicity;
    Int_t fV0CMultiplicity;
    Int_t fNoOfTPCoutTracks;    
    int      _fb;
    int      _nClusterMin;
    int      _singlesOnly;
    int      _nTPCCrossRows;
    int      _pileUpEvent;
    double   _dcaZMin;
    double   _dcaZMax;
    double   _dcaXYMin;
    double   _dcaXYMax;
    double   _chi2perTPC;
    double   _chi2perITS;
    double _vZMin;
    double _vZMax;


    ClassDef(AliAnalysisMeanPtdata, 1); // example of analysis
};

#endif

