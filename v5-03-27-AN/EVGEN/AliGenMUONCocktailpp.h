#ifndef ALIGENMUONCOCKTAILPP_H
#define ALIGENMUONCOCKTAILPP_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//
// Class to create the coktail for physics with muons for pp collisions
// using the The followoing sources: 
// jpsi, psiP, upsilon, upsilonP, upsilonPP, open charm and open beauty
//

#include "AliGenCocktail.h"
#include "AliDecayer.h"

class AliGenCocktailEntry;
class AliGenParam;

class AliGenMUONCocktailpp : public AliGenCocktail
{
 public:

    AliGenMUONCocktailpp();
    enum CMSEnergyCode { kCMS07TeV, kCMS10TeV, kCMS14TeV, kCMS09TeVpPb, kCMS09TeVPbp, kCMS03TeVPbPb, kNCMSEs };    

    virtual ~AliGenMUONCocktailpp();    
    virtual void Init();
    virtual void CreateCocktail();
    virtual void Generate();    
    Int_t   GetNSucceded()         const {return fNSucceded;}    
    Int_t   GetNGenerated()        const {return fNGenerated;}
    Int_t   GetCentralityBin()     const {return fCentralityBin;}
    Int_t   GetMuonMultiplicity()  const {return fMuonMultiplicity;}
    Float_t GetMuonPtCut()         const {return fMuonPtCut;}
    Float_t GetMuonPCut()          const {return fMuonPCut;}    
    Float_t GetMuonThetaMin()      const {return fMuonThetaMinCut;}
    Float_t GetMuonThetaMax()      const {return fMuonThetaMaxCut;}	    
    Float_t GetMuonOriginCut()     const {return fMuonOriginCut;}	    
    Float_t GetDecayModeResonance()const {return fDecayModeResonance;}
    Float_t GetDecayModePythia()   const {return fDecayModePythia;}
    
    void    SetCentralityBin(Int_t bin) { fCentralityBin = bin;}
    void    SetMuonMultiplicity(Int_t MuonMultiplicity) { fMuonMultiplicity = MuonMultiplicity;}
    void    SetMuonPtCut(Float_t PtCut) { fMuonPtCut = PtCut;}
    void    SetMuonPCut(Float_t PCut) { fMuonPCut = PCut;}    
    void    SetMuonOriginCut(Float_t originCut) { fMuonOriginCut = originCut;}
    void    SetMuonThetaRange(Float_t ThetaMin, Float_t ThetaMax){
	fMuonThetaMinCut=ThetaMin;
	fMuonThetaMaxCut=ThetaMax; }    
    void    SetDecayer(AliDecayer* const decayer){fDecayer = decayer;}
    void    SetDecayModeResonance(Decay_t decay){ fDecayModeResonance = decay;}
    void    SetDecayModePythia(Decay_t decay){ fDecayModePythia = decay;}
    void    SetResPolarization(Double_t JpsiPol, Double_t PsiPPol, Double_t UpsPol, 
    				Double_t UpsPPol, Double_t UpsPPPol, char *PolFrame);

    void    SetCMSEnergy(CMSEnergyCode cmsEnergy);
    void    SetSigmaSilent() { fSigmaSilent = kTRUE; }
    
 protected:

    //
 private:
    AliGenMUONCocktailpp(const AliGenMUONCocktailpp &cocktail); 
    AliGenMUONCocktailpp & operator=(const AliGenMUONCocktailpp &cocktail); 

    void AddReso2Generator(Char_t *nameReso, AliGenParam* const genReso, Double_t sigmaReso, Double_t polReso);
    
    AliDecayer* fDecayer; // External decayer
    Decay_t fDecayModeResonance; //decay mode in which resonances are forced to decay, default: kAll
    Decay_t fDecayModePythia; //decay mode in which particles in Pythia are forced to decay, default: kAll
    Int_t   fMuonMultiplicity; // Muon multiplicity for the primordial trigger
    Float_t fMuonPtCut;// Transverse momentum cut for muons
    Float_t fMuonPCut;// Momentum cut for muons    
    Float_t fMuonThetaMinCut;// Minimum theta cut for muons
    Float_t fMuonThetaMaxCut; // Maximum theta cut for muons
    Float_t fMuonOriginCut; //use only muons whose "part->Vz()" value is larger than fMuonOrigin
    Int_t   fNSucceded;// Number of Succes in the (di)-muon generation in the acceptance
    Int_t   fNGenerated;// Number of generated cocktails
    Int_t   fCentralityBin;// Collision centrality bin number
    Double_t fJpsiPol, fChic1Pol, fChic2Pol, fPsiPPol, fUpsPol, fUpsPPol, fUpsPPPol;//Resonances polarization parameters
    Int_t    fPolFrame;//Resonances polarization frame (Collins-Soper / Helicity)
//    Int_t fCMSEnergy; // CMS beam energy
    
    Double_t fCMSEnergyTeV;                 // energy
    Double_t fCMSEnergyTeVArray[kNCMSEs];   //!
    Double_t fSigmaReaction;                // xsec pp
    Double_t fSigmaReactionArray[kNCMSEs];  //!
    Double_t fSigmaJPsi;                    // xsec JPsi
    Double_t fSigmaJPsiArray[kNCMSEs];      //!
    Double_t fSigmaChic1;                   // xsec Chic1 
    Double_t fSigmaChic1Array[kNCMSEs];     //!
    Double_t fSigmaChic2;                   // xsec Chic2 
    Double_t fSigmaChic2Array[kNCMSEs];     //!
    Double_t fSigmaPsiP;                    // xsec Psi-prime
    Double_t fSigmaPsiPArray[kNCMSEs];      //!
    Double_t fSigmaUpsilon;                 // xsec Upsilon
    Double_t fSigmaUpsilonArray[kNCMSEs];   //!
    Double_t fSigmaUpsilonP;                // xsec Upsilon-prime 
    Double_t fSigmaUpsilonPArray[kNCMSEs];  //!
    Double_t fSigmaUpsilonPP;               // xsec Upsilon-double-prime
    Double_t fSigmaUpsilonPPArray[kNCMSEs]; //!
    Double_t fSigmaCCbar;                   // xsec correlated charm
    Double_t fSigmaCCbarArray[kNCMSEs];     //!
    Double_t fSigmaBBbar;                   // xsec correlated beauty
    Double_t fSigmaBBbarArray[kNCMSEs];     //!
    
    Bool_t   fSigmaSilent;    // hide values of cross-sections in output

    ClassDef(AliGenMUONCocktailpp,5)  //  cocktail for physics in the Alice
};

#endif



