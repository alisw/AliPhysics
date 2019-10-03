/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

// AliFlowOnTheFlyEventGenerator
// origin: Redmer Alexander Bertens (rbertens@cern.ch, rbertens@nikhef.nl, rbertens@uu.nl)

#ifndef ALIFLOWONTHEFLYEVENTGENERATOR_H
#define ALIFLOWONTHEFLYEVENTGENERATOR_H

#include "TObject.h"
#include "TClonesArray.h"
#include "TRandom.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TParticle.h"

class AliFlowEventSimple;
class TVirtualMCDecayer;

class AliFlowOnTheFlyEventGenerator : public TObject {
    public:
        // general
        AliFlowOnTheFlyEventGenerator();                
        AliFlowOnTheFlyEventGenerator(Bool_t qa, Int_t ff, Int_t mult, TVirtualMCDecayer* decayer, Bool_t a, Bool_t b, Bool_t c, Bool_t d);
        virtual ~AliFlowOnTheFlyEventGenerator();
        class NaiveFlowAndSpectrumGenerator : public TObject { // small nested helper class. use methods of AliFlowOnTheFlyEventGenerator to get access to members
            public: 
                NaiveFlowAndSpectrumGenerator(Short_t pdg, Bool_t qa, Int_t ff) : fPdg(pdg), fQA(qa), fFF(ff), fpt(0), fv2(0), fv3(0), fQApt(0), fQAv2(0), fQAv3(0) {
                    TParticle* t = new TParticle(); t->SetPdgCode(fPdg);
                    fpt = new TF1(Form("pt_%i", pdg), Form("x/TMath::Power(1+2*TMath::Sqrt(%f*%f+x*x),4)", t->GetMass(), t->GetMass()),0.,20);
                    fv2 = new TF1(Form("v2_%i", pdg), Form("TMath::Log(x+1)*.2/TMath::Power(%f*%f+x*x,.2)", t->GetMass(), t->GetMass()), 0., 20);
                    fv3 = new TF1(Form("v3_%i", pdg), Form("TMath::Log(x+1)*.1/TMath::Power(%f*%f+x*x,.2)", t->GetMass(), t->GetMass()), 0., 20);
                    if(fQA) { fQApt = new TH1F(Form("pt_%i", fPdg),Form("pt_%i", fPdg),400,0,20);
                              fQAv2 = new TH2F(Form("v2_%i", fPdg),Form("v2_%i", fPdg),400,0,20, 400, -.5, .5);
                              fQAv3 = new TH2F(Form("v3_%i", fPdg),Form("v3_%i", fPdg),400,0,20, 400, -.5, .5);      }
                    delete t;
                }
                virtual ~NaiveFlowAndSpectrumGenerator() {if(fpt) delete fpt; if(fv2) delete fv2; if(fv3) delete fv3; if(fQA) {delete fQApt; delete fQAv2; delete fQAv3;};}
                Short_t         GetPDGCode()            const {return fPdg;}
                Double_t        GetPt()                 const {double _pt(fpt->GetRandom()); if(fQA) fQApt->Fill(_pt); return _pt;}
                Double_t        GetV2(Double_t pt)      const {double _v2(fv2->Eval(pt)); if(fQA&&fFF==0) fQAv2->Fill(pt, _v2); return _v2;} 
                Double_t        GetV3(Double_t pt)      const {double _v3(fv3->Eval(pt)); if(fQA&&fFF==0) fQAv3->Fill(pt, _v3); return _v3;}
                TF1*            GetPtSpectrum()         const {return fpt;}                     // return pt fSpectrum
                TF1*            GetDifferentialV2()     const {return fv2;}                     // return fDifferential v2
                TF1*            GetDifferentialV3()     const {return fv3;}                     // return fDifferential v3
                TH1*            GetQAType(Int_t t)      const { if(t==0) return (TH1*)fQApt; 
                                                                if(t==1) return (TH1*)fQAv2; 
                                                                if(t==2) return (TH1*)fQAv3; 
                                                                return 0x0; }                   // return base pointer to QA histo class
                void            FillV2(Double_t p, Double_t v)  {fQAv2->Fill(p,v);}             // fill QA histo in case of fluctuations
                void            SetPtSpectrum(TF1* s)           {fpt = s;}                      // set custom fSpectrum
                void            SetDifferentialV2(TF1* v2)      {fv2 = v2; }                    // set custom fDifferential v2
                void            SetDifferentialV3(TF1* v3)      {fv3 = v3; }                    // set custom fDifferential v3
            private:
                Short_t                 fPdg;                                                   // pdg value of track
                Bool_t                  fQA;                                                    // make fQA plots for all generated species
                Int_t                   fFF;                                                    // introduce e-by-e flow fluctuations
                TF1*                    fpt;                                                    // !pt fSpectrum
                TF1*                    fv2;                                                    // !fDifferential v2
                TF1*                    fv3;                                                    // !fDifferential v3
                TH1F*                   fQApt;                                                  // !pt fSpectrum for fQA
                TH2F*                   fQAv2;                                                  // !v2 for fQA
                TH2F*                   fQAv3;                                                  // !v3 for fQA
                NaiveFlowAndSpectrumGenerator(const NaiveFlowAndSpectrumGenerator& dummy);              // not implemented
                NaiveFlowAndSpectrumGenerator& operator =(const NaiveFlowAndSpectrumGenerator& dummy);  // not implemented
        };      // end of NaiveFlowAndSpectrumGenerator
        // access to some members of the hidden nested class
        NaiveFlowAndSpectrumGenerator*  Find(Short_t pdg, Bool_t make);
        TObjArray*                      GetGenerators()  {return fGenerators; }
        void                            SetPtSpectrum(const char* func, Short_t pdg);
        void                            SetPtDependentV2(const char* func, Short_t pdg);  
        void                            SetPtDependentV3(const char* func, Short_t pdg);
        TF1*                            GetPtSpectrum(Short_t pdg)              {return Find(pdg, kTRUE)->GetPtSpectrum();}
        TF1*                            GetDifferentialV2(Short_t pdg)          {return Find(pdg, kTRUE)->GetDifferentialV2();}
        TF1*                            GetDifferentialV3(Short_t pdg)          {return Find(pdg, kTRUE)->GetDifferentialV3();}
        TH1*                            GetQAType(Short_t pdg, Int_t type)      {return Find(pdg, kTRUE)->GetQAType(type);}
        // event generator 
        void                    AddV2(TParticle* particle, Double_t v2, Double_t fluc);
        void                    AddV2(TClonesArray* event);
        void                    SetAfterBurnerPrecision(Double_t a, Int_t b)    { fPrecisionPhi = a; fMaxNumberOfIterations = b; }
        void                    GenerateOnTheFlyTracks(Int_t mult, Int_t pid, TClonesArray* event, Double_t fluc);
        void                    DecayOnTheFlyTracks(TClonesArray* event); 
        void                    ForceGammaDecay(TClonesArray* arr, TParticle* part); 
        AliFlowEventSimple*     GenerateOnTheFlyEvent(TClonesArray* event, Int_t nSpecies, Int_t species[], Int_t mult[], Int_t bg, Bool_t fluc);
        void                    EmbedEvent(TClonesArray* embedMe)               {fEmbedMe = embedMe;}
        AliFlowEventSimple*     ConvertTClonesToFlowEvent(TClonesArray* event, Int_t totalMultiplicity);
        void                    AddV2Mothers(Bool_t b)                          { fAddV2Mothers = b; }
        void                    AddV3Mothers(Bool_t b)                          { fAddV3Mothers = b; }
        void                    AddV2Daughters(Bool_t b)                        { fAddV2Daughters = b; }
        void                    AddV3Daughters(Bool_t b)                        { fAddV3Daughters = b; }
        // 'bookkeeping'
        void                    InitGenerators();
        void                    PrintGenerators();
        void                    DoGeneratorQA(Bool_t v2, Bool_t v3);

    private:
        // (transient) class members
        TObjArray*              fGenerators;                    // array with generators
        TClonesArray*           fEmbedMe;                       // tclones array for embedding
        AliFlowEventSimple*     fFlowEvent;                     //! flow event simple for output
        TVirtualMCDecayer*      fDecayer;                       // virtual decayer, needs to be set in macro
        Bool_t                  fAddV2Mothers;                  // add v2 to mother tracks
        Bool_t                  fAddV3Mothers;                  // add v3 to mother tracks
        Bool_t                  fAddV2Daughters;                // add v2 to daughter tracks (not implemented)
        Bool_t                  fAddV3Daughters;                // add v3 to daughter tracks (not implemented)
        Double_t                fPsi2;                          // 2nd order symmetry plane
        Double_t                fPsi3;                          // 3rd order symmetry plane  (not implemented)
        Double_t                fPrecisionPhi;                  // afterburner convergence precision
        Int_t                   fMaxNumberOfIterations;         // afterburner convergence precision
        Bool_t                  fQA;                            // save qa histograms for all generated species
        Int_t                   fFF;                            // introduce e-by-e flow fluctuations
        // assignment operator and copy constructor
        AliFlowOnTheFlyEventGenerator(const AliFlowOnTheFlyEventGenerator& dummy);              // not implemented
        AliFlowOnTheFlyEventGenerator& operator =(const AliFlowOnTheFlyEventGenerator& dummy);  // not implemented

        ClassDef(AliFlowOnTheFlyEventGenerator, 0)
};

#endif
