#ifndef ALIJETANALYSIS_H
#define ALIJETANALYSIS_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//---------------------------------------------------------------------
// JetAnalysis class 
// Perform Jet Analysis on already found jets
// Author: andreas.morsch@cern.ch, jgcn@mail.cern.ch
//         mercedes.lopez.noriega@cern.ch
//---------------------------------------------------------------------

#include <TObject.h> 
class AliLeading;
class AliJet;
class TH1;
class TH1F;
class TH2F;
class TProfile;
class TLorentzVector;

class AliJetAnalysis : public TObject
{
 public:
 
    AliJetAnalysis();
    virtual ~AliJetAnalysis();

    void Analyze();
    // define histograms
    void DefineHistograms();
    void DefineKineH();
    void DefineCorrH();
    void DefineCorr50H();
    void DefineShapH();
    void DefineFragH();
    void DefineTrigH();
    void DefineJtH();
    void DefinedNdxiH();
    // fill histograms
    void FillHistograms();
    void FillKineH();
    void FillCorrH();
    void FillCorr50H();
    void FillShapH(Float_t r);
    void FillFragH();
    void FillTrigH();
    void FillJtH();
    void FilldNdxiH();
    void FillBkgd(Int_t eventN, Int_t runN);
    // normalize histograms
    void NormHistograms();
    // plot histograms
    void PlotHistograms();
    void PlotKineH() const;
    void PlotCorrH() const;
    void PlotCorr50H() const;
    void PlotShapH() const;
    void PlotFragH() const;
    void PlotTrigH();
    // save histograms
    void SaveHistograms();
    void SaveKineH();
    void SaveCorrH();
    void SaveCorr50H();
    void SaveShapH();
    void SaveFragH();
    void SaveTrigH();
    void SaveJtH();
    void SavedNdxiH();
    // other functions
    void Shape(AliJet *j,TH1F* hs, TH1F* hr, TH1F* ha, TH2F* hd, TH2F* hp, TH1F* wd, Float_t r);
    void FragFun(AliJet *j,TH1F* hs, TH1F* hr, TH1F* ha);
    void Correlation(TLorentzVector *lv1,TLorentzVector *lv2,TH2F *h1, TH2F *h2, TH2F *h3, TH2F *h4);
    void Correlation50(AliJet *j,TLorentzVector *lv1,TLorentzVector *lv2,TH2F *h1, TH2F *h2, TH2F *h3, TH2F *h4);
    // setters
    void SetDirectory(char* directory) 
      {fDirectory = directory;}                      // directory where file with jets is
    void SetBkgdDirectory(char* directory) 
      {fBkgdDirectory = directory;}                  // directory where file with background is
    void SetOutputFile(const char* file) {fFile = file;}   // file where plots will be saved
    void SetPercentage(Float_t p) {fPercentage = p;} // minimum percentage of tracks coming from pythia (very aprox.)
    void SetEventRange(Int_t imin, Int_t imax) 
      {fEventMin = imin; fEventMax = imax;}           // first and last event
    void SetRunRange(Int_t imin, Int_t imax) 
      {fRunMin = imin; fRunMax = imax;}               // first and last run
    void SetMinimumMult(Int_t m){fminMult = m;}       // minimum multiplicity cut  
    void SetPythia(Bool_t f = kFALSE){fPythia = f;}   // If only pythia, to save everything...  
    void SetDoJt(Bool_t f = kTRUE){fDoJt = f;}        // To get j_T distribution
    void SetDodNdxi(Bool_t f = kTRUE){fDodNdxi = f;}  // To get #xi distribution
    void SetDoBkgd(Bool_t f = kTRUE) {fDoBkgd = f;}   // To get the bkgd for j_T, xi and dEdr in a hijing event
    void SetDoLeadPart(Bool_t f = kTRUE){fDoPart = f;}// To make plots for leading particle
    void SetDoGenJet(Bool_t f = kTRUE){fDoGenJ = f;}  // To make plots for generated jets
    void SetDoRecJet(Bool_t f = kTRUE){fDoRecJ = f;}  // To make plots for reconstructed jets
    void SetDoKinematics(Bool_t f = kTRUE){fDoKine = f;}       // To make the kine plots
    void SetDoCorrelations(Bool_t f = kTRUE){fDoCorr = f;}     // Correlation histograms 
    void SetDoCorr50(Bool_t f = kFALSE){fDoCorr50 = f;}        // Correlation histograms when one particle has more than 50% E
    void SetDoShape(Bool_t f = kTRUE){fDoShap = f;}            // Shape plots
    void SetDoFragmentations(Bool_t f = kTRUE){fDoFrag = f;}   // Fragmentation
    void SetDoTriggerBias(Bool_t f = kTRUE){fDoTrig = f;}      // Trigger bias plots
    void SetDivideEnergy(Float_t Efactor){fEfactor = Efactor;} // Divides E of rec.jet by Efactor
    void SetProperties(TH1* h,const char* x, const char* y) const;
    void SetReaderHeader(char *s="AliJetKineReaderHeader"){fReaderHeader = s;}
    void SetdEdrWeight();
    void SetPartPtCut(Float_t c){fPartPtCut = c;}
    void SetdrJt(Float_t r){fdrJt = r;}
    void SetdrdNdxi(Float_t r){fdrdNdxi = r;}
    void SetdrdEdr(Float_t r){fdrdEdr = r;}
    // getters
    Float_t GetdEdrWeight(Float_t eta, Float_t r);
    
 private:
    char*  fReaderHeader;    // Reader header
    char*  fDirectory;       // Directory
    char*  fBkgdDirectory;   // Directory for background
    const char*  fFile;      // Output file name
    Int_t  fEventMin;        // Minimum event number
    Int_t  fEventMax;        // Maximum event number
    Int_t  fRunMin;          // Minimum run number 
    Int_t  fRunMax;          // Maximum run number
    Int_t  fminMult;         // Minimum multiplicity for events
    Float_t fPercentage;     // percentage of pt from signal particles to accept a jet
    Float_t fPartPtCut;      // cut in the pt of particles in histos
    Float_t fdrJt;          // maximum dr for Jt plot
    Float_t fdrdNdxi;       // maximum dr for dN/dxi plot
    Float_t fdrdEdr;        // maximum dr for dE/dr plot
    Float_t fEfactor;        // factor by which energy the reconstructed jet will be divided

    Float_t fp0;    // percentage of tracks in reconstructed jet coming from pythia
                    // so far calculated in aprox. way, it needs to be improved!
    // for background from hijing events:
    Float_t fPtJ;     // P_T of the pythia jet
    Float_t fEJ;      // Energy of the pythia jet
    Float_t fEtaJ;    // Eta of the pythia jet
    Float_t fPhiJ;    // Phi of the pythia jet
    Float_t fjv3X, fjv3Y, fjv3Z;     // x,y,z of the pythia jet

    // user options    
    Bool_t fPythia;      // if pythia events
    Bool_t fDoPart;      // do analysis of leading particle
    Bool_t fDoGenJ;      // do analysis of leading generated jet
    Bool_t fDoRecJ;      // do analysis of leading rec jet
    Bool_t fDoKine;      // do kinematic plots
    Bool_t fDoCorr;      // do correlation plots
    Bool_t fDoCorr50;    // do correlation plots when one track more than 50% of jet energy
    Bool_t fDoShap;      // do shape plots
    Bool_t fDoFrag;      // do fragmentation plots 
    Bool_t fDoTrig;      // do trigger bias plots
    Bool_t fDoJt;        // do jt histo
    Bool_t fDodNdxi;     // do dN/dxi histo
    Bool_t fDoBkgd;      // get dN/dxi bkgd using hijing tracks only
    
    // weights
    Float_t fWeight;             // event weight
    Float_t fWShapR;             // weighted number of jets 
    Float_t fWFragR;             // weighted number of jets 
    Float_t fWeightdEdr[10][20]; // weight for acceptance of dE/dr histo
    Float_t fWdEdr;              // weighted number of events for dEdr histo
    Float_t fWJt;                // weight for Jt
    Float_t fWdNdxi;             // weight fro dNd#xi
    
    // leading hets and particles
    AliLeading* fPart;   // pointer to leading particle
    AliJet* fGenJ;       // pointer to leading generated jet
    AliJet* fRecJ;       // pointer to leading reconstructed jet
    AliJet* fRecB;       // pointer to leading reconstructed jet for background

    // kine histos
    TH1F *fRKineEneH;  // Reconstructed energy histo
    TH1F *fRKinePtH;   // Reconstructed Pt histo
    TH1F *fRKinePhiH;  // Reconstructed phi histo
    TH1F *fRKineEtaH;  // Reconstructed eta histo
    TH1F *fGKineEneH;  // Generated energy histo
    TH1F *fGKinePtH;   // Generated Pt histo
    TH1F *fGKinePhiH;  // Generated phi histo
    TH1F *fGKineEtaH;  // Generated eta histo
    TH1F *fPKineEneH;  // Pythia energy histo
    TH1F *fPKinePtH;   // Pythia Pt histo
    TH1F *fPKinePhiH;  // Pythia phi histo
    TH1F *fPKineEtaH;  // Pythia eta histo

    // correlation histograms
    TH2F *fPGCorrEneH;  // Energy correlation part-gen jet
    TH2F *fPGCorrPtH;   // Pt correlation part-gen jet
    TH2F *fPGCorrEtaH;  // Pseudorapidity correlation part-gen jet
    TH2F *fPGCorrPhiH;  // Azimuthal angle correlation part-gen jet
    TH2F *fPRCorrEneH;  // Energy correlation part-rec jet
    TH2F *fPRCorrPtH;   // Pt correlation part-rec jet
    TH2F *fPRCorrEtaH;  // Pseudorapidity correlation part-rec jet
    TH2F *fPRCorrPhiH;  // Azimuthal angle correlation part-rec jet
    TH2F *fRGCorrEneH;  // Energy correlation rec jet-gen jet
    TH2F *fRGCorrPtH;   // Pt correlation rec jet-gen jet
    TH2F *fRGCorrEtaH;  // Pseudorapidity correlation rec jet-gen jet
    TH2F *fRGCorrPhiH;  // Azimuthal angle correlation rec jet-gen jet
   
    // correlation histogramswhen one particle 
    // has more than 50% of the energy of the jet
    TH2F *fPRCorr50EneH;  // Energy correlation part-rec jet
    TH2F *fPRCorr50PtH;   // Pt correlation part-rec jet
    TH2F *fPRCorr50EtaH;  // Pseudorapidity correlation part-rec jet
    TH2F *fPRCorr50PhiH;  // Azimuthal angle correlation part-rec jet
    TH2F *fRGCorr50EneH;  // Energy correlation rec jet-gen jet
    TH2F *fRGCorr50PtH;   // Pt correlation rec jet-gen jet
    TH2F *fRGCorr50EtaH;  // Pseudorapidity correlation rec jet-gen jet
    TH2F *fRGCorr50PhiH;  // Azimuthal angle correlation rec jet-gen jet
   
    // fragmentation function and shape histos
    TH1F *fRFragSelH;  // Frag Fun of reconstructed jets (sel part)
    TH1F *fRFragRejH;  // Frag Fun of reconstructed jets (rej part)
    TH1F *fRFragAllH;  // Frag Fun of reconstructed jets (all part)
    TH1F *fRShapSelH;  // Shape of generated jets (sel part)
    TH1F *fRShapRejH;  // Shape of generated jets (rej part)
    TH1F *fRShapAllH;  // Shape of generated jets (all part)
    
    // trigger bias histos 
    TProfile *fGTriggerEneH;  // Generated energy (trigger bias)
    TProfile *fRTriggerEneH;  // Reconstructed energy (trigger bias)
    TProfile *fGPTriggerEneH; // Generated energy (trigger bias)
    TProfile *fPTriggerEneH;  // Leading particle energy (trigger bias)

    // dE/dr histo
    TH2F* fdEdrH;  // dE/dr histo
    TH2F* fdEdrB;  // dE/dr bkgdhisto
    TH2F* fPtEneH2;// fPtEneH2
    TH1F* fdEdrW;  // weights for dE/dr

    // Jt histo
    TH2F* fJtH;  // J_{T} histogram
    TH2F* fJtB;  // J_{T} bkgd histogram
    TH1F* fJtW;  // J_{T} weight

    // dN/dxi histo
    TH2F* fdNdxiH;  // dN/d#xi histo
    TH2F* fdNdxiB;  // dN/d#xi bkgd histo
    TH1F* fdNdxiW;  // dN/d#xi weight histo
    TH2F* fPtEneH;  // fPtEneH

protected:
    AliJetAnalysis(const AliJetAnalysis& rJetA);
    AliJetAnalysis& operator = (const AliJetAnalysis& rhsa);


    ClassDef(AliJetAnalysis,1)
};
 
#endif
