/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/**************************************************************************
 *                                                                        *
 * QA class of Heavy Flavor quark and fragmeted/decayed particles         *
 *                                                                        *
 **************************************************************************/

#ifndef _ALIHFEMCQA_H_
#define _ALIHFEMCQA_H_

#ifndef ROOT_TObject
#include <TObject.h>
#endif

class TH1F;
class TH2F;
class TParticle;
class TString;
class AliStack;

//________________________________________________________________
class AliHFEmcQA: public TObject {

        public: 
                enum heavyType {fkCharm=4, fkBeauty=5};
                enum qType {fkQuark, fkantiQuark, fkElectron, fkElectron2nd, fkeHadron, fkDeHadron};
                AliHFEmcQA();
                AliHFEmcQA(const AliHFEmcQA &p); // copy constructor
                AliHFEmcQA &operator=(const AliHFEmcQA &); // assignment operator

                virtual ~AliHFEmcQA();

                void PostAnalyze();
                void SetVerbosity(Bool_t verb=kTRUE) {fVerbos=verb;}

        protected:
                void SetStack(AliStack* stack){fStack=stack;} // set stack pointer
                void CreateHistograms(const Int_t kquark, TString hnopt=""); // create histograms for mc qa analysis
                void Init();
                void GetQuarkKine(Int_t iTrack, const Int_t kquark); // get heavy quark kinematics distribution
                void GetDecayedKine(Int_t iTrack, const Int_t kquark, const Int_t kdecayed, Bool_t isbarrel=kFALSE); // get decay electron kinematics distribution
                void EndOfEventAna(const Int_t kquark); // run analysis which should be done at the end of the event loop
                void IdentifyMother(Int_t mother_label, Int_t &mother_pdg, Int_t &grandmother_label); // 
                void HardScattering(const Int_t kquark, Int_t &motherID, Int_t &mothertype, Int_t &motherlabel); // check if the quark is produced from hard scattering
                void reportStrangeness(Int_t &motherID, Int_t &mothertype, Int_t &motherlabel); // report if the quark production process is unknown
                Bool_t IsFromInitialShower(Int_t inputmotherlabel, Int_t &motherID, Int_t &mothertype, Int_t &motherlabel); // check if the quark is produced from initial parton shower 
                Bool_t IsFromFinalParton(Int_t inputmotherlabel, Int_t &motherID, Int_t &mothertype, Int_t &motherlabel); // check if the quark is produced from final parton shower
                Float_t GetRapidity(TParticle *part); // return rapidity

                Bool_t fVerbos; //

                AliStack* fStack; //            

                static const Int_t fgkGluon; //
                static const Int_t fgkPDGInterested; //
                static const Int_t fgkMaxGener; //
                static const Int_t fgkMaxIter; //
                static const Int_t fgkqType; //


                enum ProcessType_t
                        {
                        fkPairCreationFromq,  fkPairCreationFromg,  fkFlavourExitation,  fkGluonSplitting, fkInitialPartonShower, fkLightQuarkShower
                        };

                TString fkqTypeLabel[6]; //

                struct hists{
                        TH1F *fPdgCode; // histogram to store particle pdg code
                        TH1F *fPt; // histogram to store pt
                        TH1F *fY; // histogram to store rapidity
                        TH1F *fEta; // histogram to store eta
                };
                struct histsComm {
                        TH1F *fNq; // histogram to store number of quark
                        TH1F *fProcessID; // histogram to store process id 
                        TH2F *fePtRatio; //
                        TH2F *fDePtRatio; //
                        TH2F *feDistance; //
                        TH2F *fDeDistance; //
                };

                hists fHist[2][6]; //
                histsComm fHistComm[2]; //

                TParticle *fHeavyQuark[50]; //
                Int_t fIsHeavy[2]; //
                Int_t fNparents; //
                Int_t fParentSelect[2][7]; //


		ClassDef(AliHFEmcQA,0); // HFE Monte Carlo quality assurance
};

#endif
