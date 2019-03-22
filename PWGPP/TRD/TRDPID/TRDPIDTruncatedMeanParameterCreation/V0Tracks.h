#ifndef V0TRACKS_h
#define V0TRACKS_h

#include "Rtypes.h"

#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TObject.h>

// class to read in the information from the PIDtree for the truncated mean method. More comments in VOtracks.initializeTree()

class V0Tracks {
private:
    TString inputFileString;
    TFile *inputFile;
    
public:
    TTree *tree;
    
    void setFile(TString);
    TString returnFile();
    
    void getTree();
    void initializeTree();
    
    Int_t getNumberOfEntries();
    Int_t getEntry(Int_t);
    
    // Return Functions
    Int_t       run();
    Int_t       trdNTracklets();
    Double_t    *trdSlices();
    Int_t       trdNCh();
    Int_t       trdNCls();
    Int_t       trdNClsGeneral();
    Double_t    trdSig();
    Double_t    *trdMom();
    Double_t    *trdY();
    Double_t    trdTheta();             // track
    Double_t    trdGlobalPhi();
    Double_t    *trdThetaLocal();       // layerwise    -> trdEtaLocal = trdThetaLocal (in LHC13bLHC13c_extrapolTheta)
    Double_t    *trdEtaLocal();
    Double_t    *trdPhiLocal();
    Double_t     centrality();
    Int_t       pId();
    Float_t     *nSigmaTPC();
    Float_t     *nSigmaTOF();
    Double_t    trdTPCtgl();
    
    
    // Declaration of leaf types
    Int_t           brun;
    Int_t           bTRDNtracklets;
    Double_t        bTRDslices[48];
    Double_t        bTRDMomentum[6];
    Double_t        bTRDY[6];
    Int_t           bTRDNcls;         // Anzahl der Cluster der gesamten Spur (Summe ueber alle Lagen in einem stack)
    Double_t        bTRDtheta;
    Double_t        bTRDglobalphi;
    Double_t        bTRDThetaLocal[6];
    Double_t        bTRDEtaLocal[6];
    Double_t        bTRDPhiLocal[6];
    Double_t        bCentrality;
    Double_t        bTRDsignal;         // truncated signal
    Int_t           bTRDnclsdEdx;       // Anzahl der cluster fuer truncated signal berechnung
    Int_t           bTRDnch;            // Anzahl der lagen truncated signal
    Int_t           bPDG;
    Float_t         bNSigmaTPC[3];      // tpc nsigma fuer elektronen, pionen, protonen
    Float_t         bNSigmaTOF[3];      // ensprechend tof
    Double_t        bTRDTPCtgl;
     
    // List of branches
    TBranch        *b_brun;
    TBranch        *b_bTRDNtracklets;
    TBranch        *b_bTRDslices;
    TBranch        *b_bTRDMomentum;
    TBranch        *b_bTRDY;
    TBranch        *b_bTRDNcls;
    TBranch        *b_bTRDtheta;
    TBranch        *b_bTRDglobalphi;
    TBranch        *b_bTRDthetalocal;
    TBranch        *b_bTRDetalocal;
    TBranch        *b_bTRDphilocal;
    TBranch        *b_bCentrality;
    TBranch        *b_bTRDsignal;
    TBranch        *b_bTRDnclsdEdx;
    TBranch        *b_bTRDnch;
    TBranch        *b_bPDG;
    TBranch        *b_bNSigmaTPC;
    TBranch        *b_bNSigmaTOF;
    TBranch        *b_bTRDTPCtgl;   //new TPC variable
};

#endif
