
#ifndef ALISPECTRAAODHISTOMANAGER_H
#define ALISPECTRAAODHISTOMANAGER_H

/*  See cxx source for full Copyright notice */

//-------------------------------------------------------------------------
//                      AliSpectraAODHistoManager
//
//
//
//
// Authors: Michele Floris, CERN, Philip Versteeg, UU, Redmer Bertens, UU
//-------------------------------------------------------------------------

class AliAODEvent;
class TH1F;
class TH2F;
class TList;
#include "TNamed.h"

namespace AliSpectraNameSpace
{
   enum AODPtHist_t
   {
     // MF 22/02/2012
     // Add histograms 2D DCA_xy (-3,3) vs pt (0, 3)
     // For Rec data/MC, Rec MC primaries, Rec MC secondaries weak decay, Rec MC secondaries material (x6, each particle hypothesis)

      // 6 Pt Generated True Primary
      kHistPtGenTruePrimaryProtonPlus=0,          // Pt histo for protons +, generated tracks, true ID, primary Event
      kHistPtGenTruePrimaryKaonPlus,            // Pt histo for kaons +, generated tracks, true ID, primary Event
      kHistPtGenTruePrimaryPionPlus,            // Pt histo for pions +, generated tracks, true ID, primary Event
      kHistPtGenTruePrimaryProtonMinus,         // Pt histo for protons -, generated tracks, true ID, primary Event
      kHistPtGenTruePrimaryKaonMinus,           // Pt histo for kaons -, generated tracks, true ID, primary Event
      kHistPtGenTruePrimaryPionMinus,           // Pt histo for pions -, generated tracks, true ID, primary Event
      kNPtGenHist = kHistPtGenTruePrimaryPionMinus,                    // Number of ptGen-likehistos histos
      
      // 6 Pt Reconstructed Sigma
      kHistPtRecSigmaProtonPlus,            // Pt histo for protons +, reconstructed tracks, sigma ID
      kHistPtRecSigmaKaonPlus,                  // Pt histo for kaons +, reconsructed tracks, sigma ID
      kHistPtRecSigmaPionPlus,                  // Pt histo for pions +, reconstructed tracks, sigma ID
      kHistPtRecSigmaProtonMinus,               // Pt histo for protons -, reconstructed tracks, sigma ID
      kHistPtRecSigmaKaonMinus,                 // Pt histo for kaons -, reconstructed tracks, sigma ID
      kHistPtRecSigmaPionMinus,                 // Pt histo for pions -, reconstructed tracks, sigma ID
      
      // 6 Pt Reconstructed True
      kHistPtRecTrueProtonPlus,                 // Pt histo for protons +, reconstructed tracks, true ID
      kHistPtRecTrueKaonPlus,                   // Pt histo for kaons +, reconsructed tracks, true ID
      kHistPtRecTruePionPlus,                   // Pt histo for pions +, reconstructed tracks, true ID
      kHistPtRecTrueProtonMinus,                // Pt histo for protons -, reconstructed tracks, true ID
      kHistPtRecTrueKaonMinus,                  // Pt histo for kaons -, reconstructed tracks, true ID
      kHistPtRecTruePionMinus,                  // Pt histo for pions -, reconstructed tracks, true ID
            
      // 6 Pt Reconstructed Sigma Primary
      kHistPtRecSigmaPrimaryProtonPlus,         // Pt histo for protons +, reconstructed tracks, sigma ID, primary Event
      kHistPtRecSigmaPrimaryKaonPlus,           // Pt histo for kaons +, reconsructed tracks, sigma ID, primary Event
      kHistPtRecSigmaPrimaryPionPlus,           // Pt histo for pions +, reconstructed tracks, sigma ID, primary Event
      kHistPtRecSigmaPrimaryProtonMinus,        // Pt histo for protons -, reconstructed tracks, sigma ID, primary Event
      kHistPtRecSigmaPrimaryKaonMinus,          // Pt histo for kaons -, reconstructed tracks, sigma ID, primary Event
      kHistPtRecSigmaPrimaryPionMinus,          // Pt histo for pions -, reconstructed tracks, sigma ID, primary Event
            
      // 6 Pt Reconstructed Sigma Secondary Material
      kHistPtRecSigmaSecondaryMaterialProtonPlus,       // Pt histo for protons +, reconstructed tracks, sigma ID, secondary Event
      kHistPtRecSigmaSecondaryMaterialKaonPlus,         // Pt histo for kaons +, reconsructed tracks, sigma ID, secondary Event
      kHistPtRecSigmaSecondaryMaterialPionPlus,         // Pt histo for pions +, reconstructed tracks, sigma ID, secondary Event
      kHistPtRecSigmaSecondaryMaterialProtonMinus,      // Pt histo for protons -, reconstructed tracks, sigma ID, secondary Event
      kHistPtRecSigmaSecondaryMaterialKaonMinus,        // Pt histo for kaons -, reconstructed tracks, sigma ID, secondary Event
      kHistPtRecSigmaSecondaryMaterialPionMinus,        // Pt histo for pions -, reconstructed tracks, sigma ID, secondary Event

      // 6 Pt Reconstructed Sigma Secondary WeakDecay
      kHistPtRecSigmaSecondaryWeakDecayProtonPlus,       // Pt histo for protons +, reconstructed tracks, sigma ID, secondary Event
      kHistPtRecSigmaSecondaryWeakDecayKaonPlus,         // Pt histo for kaons +, reconsructed tracks, sigma ID, secondary Event
      kHistPtRecSigmaSecondaryWeakDecayPionPlus,         // Pt histo for pions +, reconstructed tracks, sigma ID, secondary Event
      kHistPtRecSigmaSecondaryWeakDecayProtonMinus,      // Pt histo for protons -, reconstructed tracks, sigma ID, secondary Event
      kHistPtRecSigmaSecondaryWeakDecayKaonMinus,        // Pt histo for kaons -, reconstructed tracks, sigma ID, secondary Event
      kHistPtRecSigmaSecondaryWeakDecayPionMinus,        // Pt histo for pions -, reconstructed tracks, sigma ID, secondary Event

      // 6 Pt Reconstructed True Primary
      kHistPtRecTruePrimaryProtonPlus,          // Pt histo for protons +, reconstructed tracks, true ID, primary event
      kHistPtRecTruePrimaryKaonPlus,            // Pt histo for kaons +, reconsructed tracks, true ID, primary event
      kHistPtRecTruePrimaryPionPlus,            // Pt histo for pions +, reconstructed tracks, true ID, primary event
      kHistPtRecTruePrimaryProtonMinus,         // Pt histo for protons -, reconstructed tracks, true ID, primary event
      kHistPtRecTruePrimaryKaonMinus,           // Pt histo for kaons -, reconstructed tracks, true ID, primary event
      kHistPtRecTruePrimaryPionMinus,           // Pt histo for pions -, reconstructed tracks, true ID, primary event
      
      // Rest
      kHistPtRec,                               // Pt histo for all particles, reconstructed tracks
      kHistPtGen,                               // Pt histo for all particles, generated tracks
      kNPtRecHist = kHistPtGen,                    // Number of ptRec-likehistos histos
      
      kHistPIDTPC,                              // Particle Identification histo
      kHistPIDTOF,                              
      kNHistPID =kHistPIDTOF,                           
      
      kHistNSigProtonTPC,                       // NSigma separation plot    
      kHistNSigKaonTPC,                              
      kHistNSigPionTPC,                              
      kHistNSigProtonPtTPC,                              
      kHistNSigKaonPtTPC,                              
      kHistNSigPionPtTPC,                              
      
      kHistNSigProtonTOF,                              
      kHistNSigKaonTOF,                              
      kHistNSigPionTOF,                              
      kHistNSigProtonPtTOF,                              
      kHistNSigKaonPtTOF,                              
      kHistNSigPionPtTOF,                              
     
      kHistNSigProtonTPCTOF,                             
      kHistNSigKaonTPCTOF,                              
      kHistNSigPionTPCTOF,                              
      kHistNSigProtonPtTPCTOF,                              
      kHistNSigKaonPtTPCTOF,                              
      kHistNSigPionPtTPCTOF,
      kNHistNSig=kHistNSigPionPtTPCTOF,                              
      
      kHistqVecPos,
      kHistqVecNeg,
      kNHist,                                   // Total number of histos
   };  // Type of events plotted in Pt Histogram

   const char * kHistName[] =
   {
      // 6 Pt Reconstructed Sigma Primary
      "histPtGenTruePrimaryProtonPlus",         // Pt histo for protons +, generated tracks, sigma ID, primary Event
      "histPtGenTruePrimaryKaonPlus",           // Pt histo for kaons +, generated tracks, sigma ID, primary Event
      "histPtGenTruePrimaryPionPlus",           // Pt histo for pions +, generated tracks, sigma ID, primary Event
      "histPtGenTruePrimaryProtonMinus",          // Pt histo for protons -, generated tracks, sigma ID, primary Event
      "histPtGenTruePrimaryKaonMinus",            // Pt histo for kaons -, generated tracks, sigma ID, primary Event
      "histPtGenTruePrimaryPionMinus",            // Pt histo for pions -, generated tracks, sigma ID, primary Event
      
      // 6 Pt Reconstructed Sigma
      "histPtRecSigmaProtonPlus",               // Pt histo for protons +, reconstructed tracks, sigma ID
      "histPtRecSigmaKaonPlus",                 // Pt histo for kaons +, reconsructed tracks, sigma ID
      "histPtRecSigmaPionPlus",                 // Pt histo for pions +, reconstructed tracks, sigma ID
      "histPtRecSigmaProtonMinus",              // Pt histo for protons -, reconstructed tracks, sigma ID
      "histPtRecSigmaKaonMinus",                // Pt histo for kaons -, reconstructed tracks, sigma ID
      "histPtRecSigmaPionMinus",                // Pt histo for pions -, reconstructed tracks, sigma ID
      
      // 6 Pt Reconstructed True
      "histPtRecTrueProtonPlus",                // Pt histo for protons +, reconstructed tracks, true ID
      "histPtRecTrueKaonPlus",                  // Pt histo for kaons +, reconsructed tracks, true ID
      "histPtRecTruePionPlus",                  // Pt histo for pions +, reconstructed tracks, true ID
      "histPtRecTrueProtonMinus",               // Pt histo for protons -, reconstructed tracks, true ID
      "histPtRecTrueKaonMinus",                 // Pt histo for kaons -, reconstructed tracks, true ID
      "histPtRecTruePionMinus",                 // Pt histo for pions -, reconstructed tracks, true ID

      // 6 Pt Reconstructed Sigma Primary
      "histPtRecSigmaPrimaryProtonPlus",        // Pt histo for protons +, reconstructed tracks, sigma ID, primary Event
      "histPtRecSigmaPrimaryKaonPlus",          // Pt histo for kaons +, reconsructed tracks, sigma ID, primary Event
      "histPtRecSigmaPrimaryPionPlus",          // Pt histo for pions +, reconstructed tracks, sigma ID, primary Event
      "histPtRecSigmaPrimaryProtonMinus",       // Pt histo for protons -, reconstructed tracks, sigma ID, primary Event
      "histPtRecSigmaPrimaryKaonMinus",         // Pt histo for kaons -, reconstructed tracks, sigma ID, primary Event
      "histPtRecSigmaPrimaryPionMinus",         // Pt histo for pions -, reconstructed tracks, sigma ID, primary Event
      
      // 6 Pt Reconstructed Sigma Seconday
      "histPtRecSigmaSecondaryMaterialProtonPlus",      // Pt histo for protons +, reconstructed tracks, sigma ID, secondary Event
      "histPtRecSigmaSecondaryMaterialKaonPlus",        // Pt histo for kaons +, reconsructed tracks, sigma ID, secondary Event
      "histPtRecSigmaSecondaryMaterialPionPlus",        // Pt histo for pions +, reconstructed tracks, sigma ID, secondary Event
      "histPtRecSigmaSecondaryMaterialProtonMinus",     // Pt histo for protons -, reconstructed tracks, sigma ID, secondary Event
      "histPtRecSigmaSecondaryMaterialKaonMinus",       // Pt histo for kaons -, reconstructed tracks, sigma ID, secondary Event
      "histPtRecSigmaSecondaryMaterialPionMinus",       // Pt histo for pions -, reconstructed tracks, sigma ID, secondary Event
      
      // 6 Pt Reconstructed Sigma Seconday
      "histPtRecSigmaSecondaryWeakDecayProtonPlus",      // Pt histo for protons +, reconstructed tracks, sigma ID, secondary Event
      "histPtRecSigmaSecondaryWeakDecayKaonPlus",        // Pt histo for kaons +, reconsructed tracks, sigma ID, secondary Event
      "histPtRecSigmaSecondaryWeakDecayPionPlus",        // Pt histo for pions +, reconstructed tracks, sigma ID, secondary Event
      "histPtRecSigmaSecondaryWeakDecayProtonMinus",     // Pt histo for protons -, reconstructed tracks, sigma ID, secondary Event
      "histPtRecSigmaSecondaryWeakDecayKaonMinus",       // Pt histo for kaons -, reconstructed tracks, sigma ID, secondary Event
      "histPtRecSigmaSecondaryWeakDecayPionMinus",       // Pt histo for pions -, reconstructed tracks, sigma ID, secondary Event
        
      // 6 Pt Reconstructed True
      "histPtRecTruePrimaryProtonPlus",         // Pt histo for protons +, reconstructed tracks, true ID, primary event
      "histPtRecTruePrimaryKaonPlus",           // Pt histo for kaons +, reconsructed tracks, true ID, primary event
      "histPtRecTruePrimaryPionPlus",           // Pt histo for pions +, reconstructed tracks, true ID, primary event
      "histPtRecTruePrimaryProtonMinus",        // Pt histo for protons -, reconstructed tracks, true ID, primary event
      "histPtRecTruePrimaryKaonMinus",          // Pt histo for kaons -, reconstructed tracks, true ID, primary event
      "histPtRecTruePrimaryPionMinus",          // Pt histo for pions -, reconstructed tracks, true ID, primary event
      
      // Rest
      "histPtRec",                              // Pt histo for all particles, reconstructed tracks
      "histPtGen",                              // Pt histo for all particles, generated tracks
     
      "histPIDTPC",                             // Particle Identification histo
      "histPIDTOF",                             
     
      "histNSigProtonTPC",                      // NSigma Separation plot
      "histNSigKaonTPC",
      "histNSigPionTPC",
      "histNSigProtonPtTPC",
      "histNSigKaonPtTPC",
      "histNSigPionPtTPC",
      
      "histNSigProtonTOF",
      "histNSigKaonTOF",
      "histNSigPionTOF",
      "histNSigProtonPtTOF",
      "histNSigKaonPtTOF",
      "histNSigPionPtTOF",
      
      "histNSigProtonTPCTOF",
      "histNSigKaonTPCTOF",
      "histNSigPionTPCTOF",
      "histNSigProtonPtTPCTOF",
      "histNSigKaonPtTPCTOF",
      "histNSigPionPtTPCTOF",
      
      "histqPos",                             // qVecVsCentrality
      "histqNeg"
   };
   
   enum AODParticleSpecies_t
   {
     kProton = 0,
     kKaon,
     kPion,
     kNSpecies = kPion,
   }; // Particle species used in plotting

   const char * kParticleSpecies[] =
   {
       "ProtonPlus",
       "ProtonMinus",
       "KaonPlus",
       "KaonMinus",
       "PionPlus",
       "PionMinus",
   };
   enum AODHistoType_t
   {
       kHistSpectraRec = 0,
       kHistSpectraGen,
       kHNHistoTypes = kHistSpectraGen,
   }; // Types of histos

}

using namespace AliSpectraNameSpace;

class AliSpectraAODHistoManager : public TNamed
{
public:
   AliSpectraAODHistoManager() :  TNamed(), fOutputList(0) {}
   AliSpectraAODHistoManager(const char *name);
   virtual  ~AliSpectraAODHistoManager() {}


   TH2F*   BookPtGenHistogram(const char * name);
   TH2F*   BookPtRecHistogram(const char * name);
   TH2F*   BookPIDHistogram(const char * name);
   TH2F*   BookNSigHistogram(const char * name);
   TH2F*   BookqVecHistogram(const char * name);
   TH1F*   GetPtHistogram1D(const char * name,Double_t minDCA,Double_t maxDCA);
   TH1F*   GetDCAHistogram1D(const char * name,Double_t minPt,Double_t maxPt);
   TH2*     GetHistogram(UInt_t id)      {      return (TH2*) fOutputList->At(id);   }
   TH2*     GetPtHistogram(UInt_t id)    {      return (TH2*) fOutputList->At(id);   }
   TH2*     GetPtHistogram(const char * name)   {      return (TH2*) fOutputList->FindObject(name);   }
   TH2*     GetPtHistogramByName(UInt_t id)     {      return (TH2*) fOutputList->FindObject(kHistName[id]); }  // Use this if you want to read a file saved with a different histo list   
   TH2*     GetPIDHistogram(UInt_t id)   {      return (TH2*) fOutputList->At(id);   }
   TH2*     GetPIDHistogram(const char * name)  {      return (TH2*) fOutputList->FindObject(name);   }
   TH2*     GetPIDHistogramByName(UInt_t id)    {      return (TH2*) fOutputList->FindObject(kHistName[id]);  }// Use this if you want to read a file saved with a different histo list
   TH2*     GetNSigHistogram(UInt_t id)   {      return (TH2*) fOutputList->At(id);   }
   TH2*     GetNSigHistogram(const char * name)  {      return (TH2*) fOutputList->FindObject(name);   }
   TH2*     GetNSigHistogramByName(UInt_t id)    {      return (TH2*) fOutputList->FindObject(kHistName[id]);  }// Use this if you want to read a file saved with a different histo list
   TH2*     GetqVecHistogram(UInt_t id)   {      return (TH2*) fOutputList->At(id);   }
   TH2*     GetqVecHistogram(const char * name)  {      return (TH2*) fOutputList->FindObject(name);   }
   TH2*     GetqVecHistogramByName(UInt_t id)    {      return (TH2*) fOutputList->FindObject(kHistName[id]);  }// Use this if you want to read a file saved with a different histo list
  
   //TH1F*   GetTH1F(UInt_t id)            {      return (TH1F*) GetPtHistogram(id);   }
   //TH2F*   GetTH2F(UInt_t id)            {      return (TH2F*) GetPIDHistogram(id);   }

  TList * GetOutputList() {return fOutputList;}

  Long64_t Merge(TCollection* list);


private:
   TList     *fOutputList;  // List of Pt Histo's

   AliSpectraAODHistoManager(const AliSpectraAODHistoManager&);
   AliSpectraAODHistoManager& operator=(const AliSpectraAODHistoManager&);

   ClassDef(AliSpectraAODHistoManager, 1);

};
#endif

