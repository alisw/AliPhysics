#ifndef ALIJETBKG_H
#define ALIJETBKG_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 //---------------------------------------------------------------------
// Class to calculate the background per unit area
// manages the search for jets 
// Authors: Elena Bruna elena.bruna@yale.edu
//          Sevil Salur ssalur@lbl.gov
//          
//---------------------------------------------------------------------


class AliJetFinder;
class AliESDEvent;
class TTree;
class TChain;
class TString;
class AliAODEvent;
class AliJetHistos;
class AliFastJetInput;

class AliJetBkg : public TObject
{
 public:
    AliJetBkg();
    AliJetBkg(const AliJetBkg &input);
    AliJetBkg& operator=(const AliJetBkg& source);
    virtual ~AliJetBkg() {;}
    void SetHeader(AliJetHeader *header)  {fHeader=header;}
    void SetReader(AliJetReader *reader)  {fReader=reader;}
    void SetFastJetInput(AliFastJetInput *fjinput)  {fInputFJ=fjinput;}
    void BkgFastJetb(Double_t& x,Double_t& y, Double_t& z);
    void BkgFastJetWoHardest(Double_t& x,Double_t& y, Double_t& z);
    Float_t BkgFastJet();
    Float_t BkgChargedFastJet();
    Float_t BkgStat();
    Float_t BkgFastJetCone(TClonesArray* fAODJets);
    Float_t BkgRemoveJetLeading(TClonesArray* fAODJets);
    Float_t BkgRemoveJetLeadingFromUArray(TClonesArray* fAODJets);
    Float_t EtaToTheta(Float_t arg);
    Bool_t EmcalAcceptance(const Float_t eta, const Float_t phi, const Float_t radius);
    static Double_t BkgFunction(Double_t *x,Double_t *par);
    
 private:
    Double_t CalcRho(vector<fastjet::PseudoJet> input_particles,Double_t RparamBkg,TString method);
    void CalcRhob(Double_t& median, Double_t& sigma, Double_t& 
meanarea,vector<fastjet::PseudoJet> input_particles,Double_t RparamBkg,TString 
method);
    void CalcRhoWoHardest(Double_t& median, Double_t& sigma, Double_t& 
meanarea,vector<fastjet::PseudoJet> input_particles,Double_t RparamBkg,TString 
method);
    AliJetReader *fReader;   //! reader
    AliJetHeader *fHeader;   //! header
    AliFastJetInput *fInputFJ; //! input particles

  ClassDef(AliJetBkg, 1); // Analysis task for standard jet analysis
};
 
#endif
