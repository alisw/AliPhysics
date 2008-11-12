#ifndef ALIJETDISTRIBUTIONS_H
#define ALIJETDISTRIBUTIONS_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
//---------------------------------------------------------------------
// JetDistributions class 
// Get different basic distributions
// Author: mercedes.lopez.noriega@cern.ch
//---------------------------------------------------------------------

#include <TObject.h> 
class AliLeading;
class AliJet;
class TH1;
class TH1F;
class TH2F;
class TProfile;
class TLorentzVector;

class AliJetDistributions : public TObject
{
 public:
 
    AliJetDistributions();
    ~AliJetDistributions(){;}

    void Analyze();
    void DefineHistograms();
    void FillHistograms();
    void FillDistributions(AliJet *j);
    void PlotHistograms();
    void SaveHistograms();
    
    // Setter
    void SetDirectory(char* directory) {fDirectory = directory;}
    void SetOutputFile(const char* file) {fFile = file;}
    void SetPercentage(Float_t p) { fPercentage = p;}
    void SetEventRange(Int_t imin, Int_t imax) {fEventMin = imin; fEventMax = imax;}
    void SetRunRange(Int_t imin, Int_t imax) {fRunMin = imin; fRunMax = imax;}
    void SetPythia(Bool_t f = kFALSE){fPythia = f;}    
    void SetProperties(TH1* h,const char* x, const char* y) const;
    void SetReaderHeader(const char *s="AliJetKineReaderHeader") {fReaderHeader = s;}
    void SetPartPtCut(Float_t c) { fPartPtCut = c; }

    void SetDoLeadPart(Bool_t f = kTRUE) {fDoPart = f;}
    void SetDoGenJet(Bool_t f = kTRUE) {fDoGenJ = f;}
    void SetDoRecJet(Bool_t f = kTRUE) {fDoRecJ = f;}
    
 private:
    char*  fReaderHeader;// Reader header
    char*  fDirectory;   // Directory
    const char*  fFile     ;   // Output file name
    Int_t  fEventMin;    // Minimum event number
    Int_t  fEventMax;    // Maximum event number
    Int_t  fRunMin;      // Minimum run number 
    Int_t  fRunMax;      // Maximum run number
    Float_t fPercentage; // percentage of pt from signal particles to accept a jet
    Float_t fPartPtCut;  // cut in the pt of particles in histos

    // user options
    Bool_t fPythia;      // if pythia events
    Bool_t fDoPart;      // do analysis of leading particle
    Bool_t fDoGenJ;      // do analysis of leading generated jet
    Bool_t fDoRecJ;      // do analysis of leading rec jet
   
    // leading hets and particles
    AliLeading* fPart;   // pointer to leading particle
    AliJet* fGenJ;       // pointer to leading generated jet
    AliJet* fRecJ;       // pointer to leading reconstructed jet

    // histos for reconstructed particles
    TH1F* fRetaH;   // Eta of reconstructed particle
    TH1F* fRphiH;   // Phi of reconstructed particle
    TH1F* fRptH;    // Pt of reconstructed particle
    TH2F* fRetaphiH;// Eta vs Phi of reconstructed particles

    TH1F* fMultH;   // Multiplicity

 protected:
    AliJetDistributions(const AliJetDistributions& rJetD);
    AliJetDistributions& operator = (const AliJetDistributions& rhsd);

    ClassDef(AliJetDistributions,1)
};
 
#endif
