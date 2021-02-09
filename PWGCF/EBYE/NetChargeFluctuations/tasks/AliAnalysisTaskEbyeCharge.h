#ifndef AliAnalysisTaskEbyeCharge_H
#define AliAnalysisTaskEbyeCharge_H

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


class AliAnalysisTaskEbyeCharge : public AliAnalysisTaskSE
{
public:
    AliAnalysisTaskEbyeCharge();
    AliAnalysisTaskEbyeCharge(const char *name);
    virtual                 ~AliAnalysisTaskEbyeCharge();
    
    virtual void            UserCreateOutputObjects();
    virtual void            UserExec(Option_t* option);
    virtual void            Terminate(Option_t* option);
    virtual void            doAODEvent();
    virtual void            doMCAODEvent();
    void                    SetAnalysisType(Int_t IsMC){fAnalysisType = IsMC;};
    void                    FillMCEffMatrix();            // Prepare efficiency matrix
    AliEventCuts 	         fEventCuts;      /// Event cuts
    
    // Set the binning of centrality
    /*   void   SetCentralityBinning(const Int_t tmpCentbins, Float_t tmpfxCentBins[])
     {
     // Create the histograms to be used in the binning of cent
     fhCent =  new TH1F("fhCent","Centrality Bins",tmpCentbins-1 ,tmpfxCentBins );
     // ==========================================
     // prepare real data centrality bins
     fnCentbinsData = tmpCentbins;
     fnCentBinsMC   = tmpCentbins-1;
     
     fxCentBins = new Float_t[fnCentbinsData];
     for (Int_t i=0; i<fnCentbinsData; i++) fxCentBins[i] =  tmpfxCentBins[i];
     fcentDownArr = new Float_t[fnCentBinsMC];
     fcentUpArr   = new Float_t[fnCentBinsMC];
     
     for (Int_t i=0; i<fnCentbinsData-1; i++) {
     fcentDownArr[i] =  tmpfxCentBins[i];
     }
     for (Int_t i=1; i<fnCentbinsData; i++)  {
     fcentUpArr[i-1] =  tmpfxCentBins[i];
     }
     
     }*/
    
private:
    Bool_t ProperVertex(AliVEvent *event) const;
    Bool_t AcceptTrack(AliAODTrack* aodtrack) const;
    
    
    AliAODEvent            *fAOD;             // input event
    TTree                  *fTree;            // Tree for real Data
    TTree                  *fTreeMCrec;       // Tree for MC rec
    TTree                  *fTreeMCgen;       // Tree for MC gen
    TTree                  *fTreeTrackCuts;
    TTree                  *fTreeMCTrackCuts;
    TTreeSRedirector       *fTreeSRedirector; // temp tree to dump output
    
    Int_t                   fAnalysisType;    // "MC", "ESD", "AOD"
    TH1F                   *fhCent;          // helper histogram for TIdentity tree
    TH1D 	               *fHistCentralityMultSelection; // centrality class selection by AliMultSelection
    TList                  *fOutputList;      // output list
    Float_t                 fEta;             // pseudo rapidity
    Float_t                 fpT;              // transverse momentum
    Float_t                 fPhi;
    Short_t                 fCharge;          // charge information
    Float_t                 fCentrality;      // centrality information
    Int_t                   fPos;             // positive charged particles
    Int_t                   fNeg;             // negative charged particles
    Int_t                   fMomentsCross;    // moment of positive and negative charged particles
    Int_t                   fMomentsPos;      // moment of positive charged particles
    Int_t                   fMomentsNeg;      // moment of negative charged particles
    
    Float_t                 fEtaMCgen;        // MC gen pseudo rapidity
    Float_t                 fpTMCgen;         // MC gen transverse momentum
    Float_t                 fPhiMCgen;
    Float_t                 fEtaMC;           // MC rec pseudo rapidity
    Float_t                 fpTMC;            // MC rec transverse momentum
    
    Int_t                   genPos;           // MC gen positive charged particles
    Int_t                   genNeg;           // MC gen negative charged particles
    Int_t                   genMomentsCross;  // MC gen moment of positive and negative charged particles
    Int_t                   genMomentsPos;    // MC gen moment of positive charged particles
    Int_t                   genMomentsNeg;    // MC gen moment of negative charged particles
    Int_t                   recPos;           // MC gen positive charged particles
    Int_t                   recNeg;           // MC gen negative charged particles
    Int_t                   recMomentsCross;  // MC gen moment of positive and negative charged particles
    Int_t                   recMomentsPos;    // MC gen moment of positive charged particles
    Int_t                   recMomentsNeg;    // MC gen moment of negative charged particles
    
    Double_t                fVxMax;            //vxmax
    Double_t                fVyMax;             //vymax
    Double_t                fVzMax;             //vzmax
    
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
    Int_t                   pdgPrimgen ;
    
    Double_t                fEtaDown;
    Double_t                fEtaUp;
    /*   Int_t              fnCentBinsMC;
     Int_t                  fnCentbinsData;
     Float_t               *fcentDownArr;          //[fnCentBinsMC]
     Float_t               *fcentUpArr;            //[fnCentBinsMC]
     Float_t               *fxCentBins;  */        //[fnCentbinsData]
    
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
    TH1F                    *hGenPt;
    TH1F                    *hGenPhi;
    TH1F                    *hGenEta;
    TH1F                    *hTrackPt;
    TH1F                    *hTrackPhi;
    TH1F                    *hTrackEta;
    
    AliAnalysisTaskEbyeCharge(const AliAnalysisTaskEbyeCharge&);
    AliAnalysisTaskEbyeCharge& operator=(const AliAnalysisTaskEbyeCharge&);
    ClassDef(AliAnalysisTaskEbyeCharge, 1);
};


#endif

