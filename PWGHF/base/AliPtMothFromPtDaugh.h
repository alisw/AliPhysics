#ifndef ALIPTMOTHFROMPTDAUGH_H
#define ALIPTMOTHFROMPTDAUGH_H
/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

////////////////////////////////////////////////////////////////////////////
//  Class to perform pt-spectra (and ptMin-spectra) extraction of mothers //
//  particles starting from measured pt-spectra of daughters particles    //
//  that come from inclusive decays.                                      //
//                                                                        //  
//  Contact: Giuseppe.Bruno@ba.infn.it  &  Fiorella.Fionda@ba.infn.it     //
////////////////////////////////////////////////////////////////////////////

class TH1F;
class TNtuple;
class AliStack;
class TParticle;
class TArrayI;
#include <TNamed.h>
#include <TString.h>

//________________________________________________________________
class AliPtMothFromPtDaugh : public TNamed {
  
public:

  typedef enum {
          kUserAnalysis,
          kBtoJPSI,
          kBtoEle,
          kBtoMuon,
          kBtoD0
          } Analysis_mode;
  
  AliPtMothFromPtDaugh();
  AliPtMothFromPtDaugh(const char* name, const char* title);
  virtual ~AliPtMothFromPtDaugh();
  
  Bool_t CreateWeights(); // method which create the containers for the corrections for the deconvolution
  void DeleteWeights();   // reset the containers
  Bool_t ReadHistoPtDaught(const TH1F *hist); // this is input, i.e. the spectrum of the daughter.
  //
  // setters  
  //
  void SetDefaultAnalysis(Analysis_mode mode){fAnalysisMode = mode; InitDefaultAnalysis();}
  void SetPdgDaugh(Int_t pdgD); // set the pdg of the daughter particle
  void SetBeautyMothers(); // all the Beauty particles are considered as "mothers" 
  void SetPdgMothers(Int_t n_mothers,Int_t *pdgM);  // define which particles are considered as possible mothers 
  
  // (here you can set the pt binning of the final spectrum, i.e. that of the mothers )
  void SetBinsPtMoth(Double_t ptmin, Double_t ptmax,Int_t nbins,Double_t alpha=1.0); // binning of histo pt(Mothers)
  void SetBinsPtMoth(Int_t nbins, Double_t *edgeBins);                          //   "      "    "      " 
  void SetBinsPtMinMoth(Double_t ptmin, Double_t ptmax,Int_t nbins,Double_t alpha=1.0); // histo of pt_min(Mothers)
  void SetBinsPtMinMoth(Int_t nbins, Double_t *edgeBins);                          // histo of pt_min(Mothers)
  //
  void SetYmothers(Double_t yMin, Double_t yMax){fyMothMin = yMin;                     // define rapidity range 
                                                  fyMothMax = yMax; SetUseEta(kFALSE);} // of the mothers
  void SetYdaughter(Double_t yMin, Double_t yMax){fyDaughMin = yMin;                     // define rapidity range
                                                   fyDaughMax = yMax; SetUseEta(kFALSE);} // of the daughter

  void SetEtaMothers(Double_t etaMin, Double_t etaMax){fyMothMin = etaMin;   // define pseudo-rapidity range of mothers
                                                        fyMothMax = etaMax; SetUseEta(kTRUE);}
  void SetEtaDaughter(Double_t etaMin, Double_t etaMax){fyDaughMin = etaMin; // define pseudo-rapidity range of daughter
                                                        fyDaughMax = etaMax; SetUseEta(kTRUE);}
  void SetDecayNtupla(TNtuple *DecKine){fDecayKine = DecKine;}               // define TNtupla with kinematic informations 
  //getters
  Double_t* GetBinsSize(const TH1F *hist, Int_t &n) const;
  Bool_t GetEtaMothers(Double_t &etaMin, Double_t &etaMax) const;  // get pseudorapidity edges for mothers if pseudorapidity
                                                                    // is used    (return kFALSE if rapidity has been used) 
  Bool_t GetEtaDaughter(Double_t &etaMin, Double_t &etaMax) const; // get pseudorapidity edges for daughters if pseudorapidity
                                                                    // is used    (return kFALSE if rapidity has ben used)
  Bool_t GetYMothers(Double_t &yMin, Double_t &yMax) const;    // get rapidity edges for mothers if rapidity is used
  Bool_t GetYDaughter(Double_t &yMin, Double_t &yMax) const;   // get rapidity edges for daughters if rapidity is used

  Int_t GetPdgDaugh() const {return fDaughter;} // return the pdg of the daughter particle
  Int_t* GetPdgMothers(Int_t &n_mothers) const; // return the pdg codes of the mothers particles 
 
  // main method to evaluate mothers spectra (for pt and ptMin)
   Bool_t EvaluatePtMoth();
   void WritePtMothHistoToFile(TString fileOutName="Mothers.root");

  // return values of correction factors for pt-mothers ditribution
  Double_t GetW(Int_t i,Int_t j) const;
  Double_t GetF(Int_t i) const;

  // return values of correction factors for ptMin-mothers ditribution
  Double_t GetWmin(Int_t i,Int_t j) const;
  Double_t GetFmin(Int_t i) const;
   
  // return errors on correction factors for pt-mothers distribution
  Double_t GetStatErrW(Int_t i,Int_t j) const; 
  Double_t GetStatErrF(Int_t i) const;
 
  // return errors on correction factors for ptMin-mothers distribution
  Double_t GetStatErrWmin(Int_t i,Int_t j) const;
  Double_t GetStatErrFmin(Int_t i) const;

   // return pointers to pt and ptMin mothers histograms 
   TH1F* GetHistoPtMother() const {return fHistoPtMothers;}
   TH1F* GetHistoPtMinMother() const {return fHistoPtMinMothers;}

   //method to read kinematic
   Int_t GiveBinIndex(Double_t Ptpart,const TH1F *ptHist) const;
   Bool_t Rapidity(const TParticle *particle, Double_t &y);
   Bool_t IsMothers(Int_t pdgCode);
   Bool_t IsSelectedDaugh(const TParticle *part, Int_t &labelDaugh, AliStack * const stack);
   Bool_t CutDaugh(Double_t yD,Double_t ptD);

private:
  void InitDefaultAnalysis(); // set the default analysis for one of the specific Analysis_mode, apart kUserAnalysis
  // methods to evaluate correction factors
  Bool_t EvaluateWij(); // 
  Bool_t EvaluateFi();  // 
 
  void SetUseEta(Bool_t useEta){fUseEta=useEta;}   // decide whether to use rapidity or pseudo-rapidity cut
  Double_t* SetBinsSize(Double_t ptmin, Double_t ptmax,Int_t nbins,Double_t alpha);  
  // method to evaluate raw-mothers spectra (for pt and ptMin)
  Bool_t EvaluatePtMothRaw(TH1F *histoPt, TH1F *histoPtMin); 
  //method to evaluate statistical errors on extracted distribution
  Bool_t EvaluateErrPt(Double_t *erStat);
  Bool_t EvaluateErrPtMin(Double_t *erStat);
  void SetPdgMothersPrivate(Int_t n_mothers,Int_t *pdgM);  // define which particles are considered as possible mothers 

  TNtuple*             fDecayKine; //Ntupla to store kinematic information of Decay 
  Double_t**           fWij; //pointer to store correction factors
  Double_t*            fFi;  //"            "                "
  Double_t**           fWijMin; //"        "                "
  Double_t*            fFiMin; //"         "                "
  TH1F*                fHistoPtDaughter; //pointers to pt-histogram of Daughther 
  TH1F*                fHistoPtMothers; //pointers to ptMin-histogram of mothers
  TH1F*                fHistoPtMinMothers; //pointers for pt_min-histogram of Mothers
  TArrayI*             fMothers; //Array with pdg codes of mothers
  Int_t                fDaughter; //pdg code of daughter 
  Double_t             fyMothMax; //max rapidity (or pseudorapidity) of mothers
  Double_t             fyMothMin; //min  "             "                 " 
  Double_t             fyDaughMax; //max "             "             of daughters  
  Double_t             fyDaughMin; //min "             "                 "
  Bool_t               fUseEta; //kTRUE if pseudorapidity range is used 
  Analysis_mode        fAnalysisMode; //analysis mode
  
  AliPtMothFromPtDaugh(const AliPtMothFromPtDaugh &c);
  AliPtMothFromPtDaugh &operator=(const AliPtMothFromPtDaugh &c);
  
  ClassDef(AliPtMothFromPtDaugh,1); // mother pt spectrum from daughter pt spectrum
};


#endif
