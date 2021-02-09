#ifndef AliAnalysisTaskEbyeChargeFlucpPb_H
#define AliAnalysisTaskEbyeChargeFlucpPb_H

class TH1F;
class TH1I;
class TH2;
class THn;
class TH1D;
class TList;
class AliAODtrack;
class AliAnalysisUtils;
class TTree;
class AliAODEvent;
class AliVEvent;
class TString;
class TObjArray;

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "THnSparse.h"
#include "THn.h"
#include "TTreeStream.h"


class AliAnalysisTaskEbyeChargeFlucpPb : public AliAnalysisTaskSE
{
public:
    AliAnalysisTaskEbyeChargeFlucpPb();
    AliAnalysisTaskEbyeChargeFlucpPb(const char *name);
    virtual                 ~AliAnalysisTaskEbyeChargeFlucpPb();
    
    virtual void            UserCreateOutputObjects();
    virtual void            UserExec(Option_t* option);
    virtual void            Terminate(Option_t* option);
    virtual void            doAODEvent();
    virtual void            doMCAODEvent();
    void                    SetAnalysisType(Bool_t IsMC = kTRUE){fAnalysisType = IsMC;};
    void                    FillMCEffMatrix();            // Prepare efficiency matrix
    AliEventCuts 	        fEventCuts;                   // Event cuts
//    void                    SetMaxTPCCluster(Int_t MaxTPCclus){fCutTPCMaxCls=MaxTPCclus;}
//    void                    SetNTPCCluster(Int_t TPCNclus){fCutTPCNCls=TPCNclus;}
//    void                    SetDCACut(Double_t DCAxyCut,Double_t DCAzCut){
//                            fCutDCAxy=DCAxyCut;
//                            fCutDCAz=DCAzCut;
//                             }
    void                    Setzvtxcut(Int_t zvtxcut){fzvtxcut=zvtxcut;};
    void                    SettrackBit(Int_t trackBit){ftrackBit=trackBit;};
    
private:
    Bool_t ProperVertex(AliVEvent *event) const;
    Bool_t AcceptTrack(AliAODTrack* aodtrack) const;
    Bool_t PassDCA(AliAODEvent *fAOD,AliAODTrack* aodtrack) const;
    
    
    AliAODEvent            *fAOD;                   // input event
    TClonesArray           *fArrayMC;               // array of MC particles
    TTree                  *fTree;                  // Tree for real Data
    TTree                  *fTreeMCrec;             // Tree for MC rec
    TTree                  *fTreeMCgen;             // Tree for MC gen
    TTreeSRedirector       *fTreeSRedirector;       // temp tree to dump output
    
    Bool_t                 fAnalysisType;          // "MC",  "AOD"
    TH1F                   *fhCent;                 // helper histogram for centrality
    TH1D                   *fHistCentralityMultSelection; // centrality class selection by AliMultSelection

    TList                  *fOutputList;            // output list
    Float_t                 fEta;                   // pseudo rapidity
    Float_t                 fpT;                    // transverse momentum
    Float_t                 fPhi;
    Short_t                 fCharge;                // charge information
    Float_t                 fCentrality;            // centrality information
    Int_t                   fPos;                   // positive charged particles
    Int_t                   fNeg;                   // negative charged particles
    Int_t                   fMomentsCross;          // moment of positive and negative charged particles
    Int_t                   fMomentsPos;            // moment of positive charged particles
    Int_t                   fMomentsNeg;            // moment of negative charged particles
    Int_t                   EventNumber;
    
    Float_t                 fEtaMCgen;              // MC gen pseudo rapidity
    Float_t                 fpTMCgen;               // MC gen transverse momentum
    Float_t                 fPhiMCgen;
    Short_t                 fChargegen;             // charge information
    Double_t                fVxMax;                 // vxmax
    Double_t                fVyMax;                 // vymax
    Double_t                fVzMax;                 // vzmax
    
    Int_t                   fRunNumber;
    Int_t                   genPos;                 // MC gen positive charged particles
    Int_t                   genNeg;                 // MC gen negative charged particles
    Int_t                   genNch;                 // MC gen negative charged particles
    Int_t                   genMomentsCross;        // MC gen moment of positive and negative charged particles
    Int_t                   genMomentsPos;          // MC gen moment of positive charged particles
    Int_t                   genMomentsNeg;          // MC gen moment of negative charged particles
    
    Int_t                   recPos;                 // MC rec positive charged particles
    Int_t                   recNeg;                 // MC rec negative charged particles
    Int_t                   recNch;                 // MC rec total charged particles
    Int_t                   Nch;
    Int_t                   recMomentsCross;        // MC rec moment of positive and negative charged particles
    Int_t                   recMomentsPos;          // MC rec moment of positive charged particles
    Int_t                   recMomentsNeg;          // MC rec moment of negative charged particles

    Int_t                   pdg;
    Int_t                   pdgMom;
    Int_t                   pdgMomPhysicalPrim;
    Int_t                   pdgMomPrim ;
    Int_t                   pdgPhysicalPrim ;
    Int_t                   pdgPrim ;
    
    Int_t                   pdggen;
    Int_t                   pdgMomgen;
    Int_t                   pdgMomPhysicalPrimgen;
    Int_t                   pdgMomPrimgen;
    Int_t                   pdgPhysicalPrimgen ;
    Int_t                   pdgPrimgen;
   
    Int_t                   fzvtxcut;
    Int_t                   ftrackBit;
//    Int_t                   fCutTPCMaxCls; // no. of TPC crossed rows
//    Int_t                   fCutTPCNCls;  // TPC n clusters
//    Double_t                fCutDCAxy ;   // DCA xy cut
//    Double_t                fCutDCAz ;    //DCA z
    
    
    Double_t                fEtaDown;
    Double_t                fEtaUp;
    
    // control and QA histograms
    //
    THnF                    *fHistPosEffMatrixRec;     // histogram efficiency matrix --> reconstructed traks
    THnF                    *fHistNegEffMatrixRec;     // histogram efficiency matrix --> generated traks
    THnF                    *fHistPosEffMatrixGen;     // histogram efficiency matrix --> reconstructed pions
    THnF                    *fHistNegEffMatrixGen;     // histogram efficiency matrix --> generated pions
    
    TH1I                    *fHistVertexNconributors;  // vertex contributors number
    TH1I                    *fHistVertexStats;         // Vertex reconstruction statistics
    TH1D                    *fHistVx;                  // Vx hist
    TH1D                    *fHistVy;                  // Vy hist
    TH1D                    *fHistVz;                  // Vz hist
    
    TH2F                    *fHistZVertexCent;
    TH1I                    *fEventStatistics;         // cut-by-cut counter of events
    TH1D                    *fHistClustersTPC;
    TH1D                    *fHistChi2perNDF;
    TH1F                    *fHistDCAz;
    TH1F                    *fHistDCAxy;
    TH1F                    *fHistMagneticField;
    TH1F                    *hGenPt;
    TH1F                    *hGenPhi;
    TH1F                    *hGenEta;
    TH1F                    *hTrackPt;
    TH1F                    *hTrackPhi;
    TH1F                    *hTrackEta;
    TH1F                    *hTrackPtallrec;
    TH1F                    *hTrackPhiallrec;
    TH1F                    *hTrackEtaallrec;
    
    AliAnalysisTaskEbyeChargeFlucpPb(const AliAnalysisTaskEbyeChargeFlucpPb&);
    AliAnalysisTaskEbyeChargeFlucpPb& operator=(const AliAnalysisTaskEbyeChargeFlucpPb&);
    ClassDef(AliAnalysisTaskEbyeChargeFlucpPb, 1);
};


#endif

