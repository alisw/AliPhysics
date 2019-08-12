#ifndef AliAnalysisTaskEbyeNetChargeMCPbPbESD_H
#define AliAnalysisTaskEbyeNetChargeMCPbPbESD_H

class TH1F;
class TH1I;
class TH2;
class THn;
class TH1D;
class TList;
class AliESDtrack;
class AliAnalysisUtils;
class TTree;
class AliESDEvent;
class AliVEvent;
class TString;
class TObjArray;
class AliMCParticle;
//class AliStack;

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliVEvent.h" 
#include "THn.h"
#include "TTreeStream.h"


class AliAnalysisTaskEbyeNetChargeMCPbPbESD : public AliAnalysisTaskSE
{
public:
    AliAnalysisTaskEbyeNetChargeMCPbPbESD();
    AliAnalysisTaskEbyeNetChargeMCPbPbESD(const char *name);
    virtual                 ~AliAnalysisTaskEbyeNetChargeMCPbPbESD();
    
    virtual void            UserCreateOutputObjects();
    virtual void            UserExec(Option_t* option);
    virtual void            Terminate(Option_t* option);
    AliEventCuts 	        fEventCuts;                   // Event cuts
    
private:
  

    TClonesArray           *fArrayMC;               // array of MC particles
 
    TTree                  *fTreeMCgen;             // Tree for MC gen
    TTreeSRedirector       *fTreeSRedirector;       // temp tree to dump output
    AliStack               *fStack;                     //! MC stack
    TH1F                   *fhCent;                 // helper histogram for centrality
    TH1D                   *fHistCentralityMultSelection; // centrality class selection by AliMultSelection
    TH1F                   *fHistCentralityImpPar;         // control histogram for centrality

    
    TList                  *fOutputList;            // output list
 
    Int_t                   EventNumber;
    Int_t                   fpT;
    Int_t                   fPhi;
    Int_t                   fEta;
    
    Float_t                 fCentrality;             // centrality information


    Float_t                 fEtaMCgen;              // MC gen pseudo rapidity
    Float_t                 fpTMCgen;               // MC gen transverse momentum
    Float_t                 fPhiMCgen;
    Short_t                 fChargegen;             // charge information
 
    Double_t                fMCImpactParameter;
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
    
    Int_t                   pdggen;
    Int_t                   pdgMomgen;
    Int_t                   pdgMomPhysicalPrimgen;
    Int_t                   pdgMomPrimgen;
    Int_t                   pdgPhysicalPrimgen ;
    Int_t                   pdgPrimgen;
    
    Double_t                fEtaDown;
    Double_t                fEtaUp;
    
    TH1I                    *fHistVertexStats;         // Vertex reconstruction statistics
    TH1D                    *fHistVx;                  // Vx hist
    TH1D                    *fHistVy;                  // Vy hist
    TH1D                    *fHistVz;                  // Vz hist
    
    TH2F                    *fHistZVertexCent;
    TH1I                    *fEventStatistics;         // cut-by-cut counter of events
    TH1F                    *hGenPt;
    TH1F                    *hGenPhi;
    TH1F                    *hGenEta;

  
    AliAnalysisTaskEbyeNetChargeMCPbPbESD(const AliAnalysisTaskEbyeNetChargeMCPbPbESD&);
    AliAnalysisTaskEbyeNetChargeMCPbPbESD& operator=(const AliAnalysisTaskEbyeNetChargeMCPbPbESD&);
    ClassDef(AliAnalysisTaskEbyeNetChargeMCPbPbESD, 1);
};


#endif

