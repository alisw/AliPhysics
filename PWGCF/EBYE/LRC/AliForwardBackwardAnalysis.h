#ifndef AliForwardBackwardAnalysis_cxx
#define AliForwardBackwardAnalysis_cxx

// analysis task which gets event tree with tracks

class TH1F;
class TH2F;
class TH1D;
class TH2D;
class TH3D;
class TList;
class TProfile;
class TTree;
class TF1;

class AliAODEvent;

#include "AliAnalysisTaskSE.h"


#define N_PT_BINS 3

#define DPT_N_ETAPHI_BINNINGS 1 //6 //5
#define DPT_N_CENTR_BINS 1 //10

#define FB_N_PHI_BINS 32
#define FB_N_ETAPHI_BINNINGS 2060 //1024*2+12 //12
#define FB_N_CENTR_BINNINGS 3




#define FB_max_tracks_in_win 2000

#include <iostream>
using namespace std;





class AliForwardBackwardAnalysis : public AliAnalysisTaskSE {
public:
    AliForwardBackwardAnalysis() : AliAnalysisTaskSE(){}//, fAOD(0), fOutputList(0)

    AliForwardBackwardAnalysis(const char *name);
    virtual ~AliForwardBackwardAnalysis() {}

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);

    void SetCentralityEstimator(const char* centralityEstimator) { fCentralityEstimator = centralityEstimator; }
    void SetAODtrackCutBit(Int_t bit){ fAODTrackCutBit = bit; }
private:
    //  AliESDEvent *fESD;    //! ESD object
    //  AliAODEvent *fAOD;    //! AOD object
    TList       *fOutputList; //! Output list

    TString fCentralityEstimator;//"V0M","TRK","TKL","ZDC","FMD"
    Int_t fAODTrackCutBit;//track cut bit from track selection (only used for AODs)

    Int_t fRunNumber, fNumberOfTracks; //
    // There are 3 centrlaity information for each event
    // "VOM","TRK","ZEMvsZDC"
    Float_t fCentPercentile[3]; //
    Float_t fVertexX,fVertexY,fVertexZ; //

    TH1F        *fHistVertexZ; //!
    TH2F        *fHistVertexXY; //!
    TH1F        *fHistCentrality; //!
    TH1F        *fHistTrackPt; //! Pt spectrum
    TH1F        *fHistTrackEta; //! Eta spectrum
    TH1F        *fHistTrackPhi; //! Phi spectrum
    TH1F        *fHistLog; //!


    TH2D * fHist2D_ESTIMATOR_VS_multTPC; //!
    TH2D * fHist2D_ESTIMATOR_VS_multTPC_afterCuts; //!
    TF1 *fBorderToCutOutliersLower; //!
    TF1 *fBorderToCutOutliersUpper; //!




    // ##### for FB 2D
    TH2D* fHistPtPt[FB_N_PHI_BINS][FB_N_PHI_BINS];	//! Pt-Pt 2D Profile
    TH2D* fHistNN[FB_N_PHI_BINS][FB_N_PHI_BINS];	//! N-N 2D Profile

    struct WinPairInfo
    {
        // ### win acceptance:
        //eta
        Double_t eta[4];
        //phi
        Double_t phi[4];
        //pt
        Double_t pt[4];

        // ### win pair filling:
        UShort_t nF;
        UShort_t nB;

        Double_t PtF;
        Double_t PtB;

        Double_t pipjF; //, sum_piF;
        Double_t pipjB; //, sum_piB;
        // to calc C:
        double arr_pF[FB_max_tracks_in_win];
        double arr_pB[FB_max_tracks_in_win];

        // for dptdpt:
        Double_t piFpjB;


        WinPairInfo()
        {
            reset();

            eta[0] = -5; eta[1] = 0;
            eta[2] = 0; eta[3] = 5;

            phi[0] = 0; phi[1] = TMath::TwoPi();
            phi[2] = 0; phi[3] = TMath::TwoPi();

            pt[0] = 0.2; pt[1] = 2.0;
            pt[2] = 0.2; pt[3] = 2.0;
        }

        void setEtaRanges( Double_t etaF_min, Double_t etaF_max, Double_t etaB_min, Double_t etaB_max )
        {
            eta[0] = etaF_min; eta[1] = etaF_max;
            eta[2] = etaB_min; eta[3] = etaB_max;
        }
        void setPhiRanges( Double_t phiF_min, Double_t phiF_max, Double_t phiB_min, Double_t phiB_max )
        {
            phi[0] = phiF_min; phi[1] = phiF_max;
            phi[2] = phiB_min; phi[3] = phiB_max;
        }
        void setPtRanges( Double_t ptF_min, Double_t ptF_max, Double_t ptB_min, Double_t ptB_max )
        {
            pt[0] = ptF_min; pt[1] = ptF_max;
            pt[2] = ptB_min; pt[3] = ptB_max;
        }

        void setFBacceptance( Double_t *_eta, Double_t *_phi, Double_t *_pt )
        {
            for( Int_t i=0; i < 4; i++ )
            {
                eta[i] = _eta[i];
                phi[i] = _phi[i];
                pt[i] = _pt[i];
            }
        }

        void reset()
        {
            if(0)cout << "reset(): eta: " << eta[0] << "..." << eta[1] << ", " << eta[2] << "..." << eta[3]
                      << ", phi: " << phi[0] << "..." << phi[1] << ", " << phi[2] << "..." << phi[3]
                      << ", pt: " << pt[0] << "..." << pt[1] << ", " << pt[2] << "..." << pt[3]
                      << ", nF=" << nF << ", nB=" << nB << ", PtF=" << PtF << ", PtB=" << PtB
                      << endl;

            nF = 0;
            nB = 0;

            PtF = 0;
            PtB = 0;

            pipjF = 0; //sum_piF = 0;
            pipjB = 0; //sum_piB = 0;

            piFpjB = 0;
        }
        void addTrack(double _pt, double _eta, double _phi, Short_t _charge, Int_t _pid)
        {
            // F:
            if( _pt > pt[0] && _pt < pt[1] )
                if( _eta > eta[0] && _eta < eta[1] )
                    if( _phi > phi[0] && _phi < phi[1] )
                    {
                        arr_pF[nF] = _pt;
                        nF++;
                        PtF += _pt;
                    }
            // B:
            if( _pt > pt[2] && _pt < pt[3] )
                if( _eta > eta[2] && _eta < eta[3] )
                    if( _phi > phi[2] && _phi < phi[3] )
                    {
                        arr_pB[nB] = _pt;
                        nB++;
                        PtB += _pt;
                    }
        }
        void finalActionsForEvent()
        {
            //F
            if ( nF > 0 )
                PtF /= nF;
            else
                PtF = -1;

            //B
            if ( nB > 0 )
                PtB /= nB;
            else
                PtB = -1;

            // to calc C:
            // F
            for( Int_t i = 0; i < nF; i++ )
                for( Int_t j = i+1; j < nF; j++ )
                    pipjF += arr_pF[i]*arr_pF[j];

            // B
            for( Int_t i = 0; i < nB; i++ )
                for( Int_t j = i+1; j < nB; j++ )
                    pipjB += arr_pB[i]*arr_pB[j];

            // for dpt-dpt
            for( Int_t i = 0; i < nF; i++ )
                for( Int_t j = 0; j < nB; j++ )
                    piFpjB += arr_pF[i]*arr_pB[j];

        }
    };

    // ##### for FB:
    WinPairInfo winPairs[N_PT_BINS][FB_N_ETAPHI_BINNINGS];
    TH3D* fHistFB_ingredients_etaPhi_centr[N_PT_BINS][FB_N_CENTR_BINNINGS]; //!

    // #####for DptDpt
    TH2D* fHistDptDpt_p1p2[N_PT_BINS][DPT_N_CENTR_BINS][DPT_N_ETAPHI_BINNINGS]; //!
    TH2D* fHistDptDpt_p1[N_PT_BINS][DPT_N_CENTR_BINS][DPT_N_ETAPHI_BINNINGS]; //!
    TH2D* fHistDptDpt_p2[N_PT_BINS][DPT_N_CENTR_BINS][DPT_N_ETAPHI_BINNINGS]; //!
    TH2D* fHistDptDpt_n1n2[N_PT_BINS][DPT_N_CENTR_BINS][DPT_N_ETAPHI_BINNINGS]; //!

    TH1D* fHistDptDpt_meanPt[N_PT_BINS][DPT_N_CENTR_BINS]; //!

    AliForwardBackwardAnalysis(const AliForwardBackwardAnalysis&); // not implemented
    AliForwardBackwardAnalysis& operator=(const AliForwardBackwardAnalysis&); // not implemented

    ClassDef(AliForwardBackwardAnalysis, 3); // example of analysis
};

#endif
