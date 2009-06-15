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
 * Secondary vertexing construction Class                                 *
 *                                                                        *
 **************************************************************************/


#ifndef _ALIHFESECVTX_H_
#define _ALIHFESECVTX_H_

#ifndef ROOT_TObject
#include <TObject.h>
#endif

class TH1F;
class TH2F;
class TString;
class AliESDEvent;
class AliESDtrack;
class AliStack;

//________________________________________________________________
class AliHFEsecVtx : public TObject {

        public: 
                AliHFEsecVtx();
                AliHFEsecVtx(const AliHFEsecVtx &p); // copy constructor
                AliHFEsecVtx &operator=(const AliHFEsecVtx &); // assignment operator
                virtual ~AliHFEsecVtx();


        protected:
                void CreateHistograms(TString hnopt="");

                void Init();
                void InitAnaPair(); // initialize default parameters 
                void SetEvent(AliESDEvent* ESD){fESD1=ESD;}; // set ESD pointer
                void SetStack(AliStack* stack){fStack=stack;} // set stack pointer
                void AnaPair(AliESDtrack* ESDtrack1, AliESDtrack* ESDtrack2, Int_t sourcePart, Int_t index2); // do e-h analysis
                void RunSECVTX(AliESDtrack *track); // run secondary vertexing algorithm
                void FindSECVTXCandid4Tracks(AliESDtrack* track1); // find secondary vertex for 4 tracks
                void FindSECVTXCandid3Tracks(AliESDtrack* track1); // find secondary vertex for 3 tracks        
                void FindSECVTXCandid2Tracks(AliESDtrack* track1); // find secondary vertex for 2 tracks
                void CalcSECVTXProperty(AliESDtrack* track1, AliESDtrack* track2, AliESDtrack* track3); // calculated distinctive variables
                void ApplyPairTagCut(); // apply strict tagging condition for 1 pair secondary vertex
                void ApplyVtxTagCut(Double_t chi2); // apply tagging condition for 3 track secondary vertex
                //void ApplyTagCut(Double_t chi2, AliESDtrack* track1, AliESDtrack* track2, AliESDtrack* track3); // apply tagging condition for 3 track secondary vertex
                void ResetTagVar(); // reset tagging variables

                Bool_t ApplyPairTagCut(Int_t id); // apply strict tagging condition for 1 pair secondary vertex
                Bool_t ApplyTagCut(Double_t chi2); // apply tagging condition
                Bool_t SingleTrackCut(AliESDtrack* track1); // single track cut

                Int_t GetMCPID(AliESDtrack *track); // return MC pid
                Int_t PairOrigin(AliESDtrack* track1, AliESDtrack* track2); // return pair origin as a pdg code
                Int_t PairCode(AliESDtrack* track1,AliESDtrack* track2); // return corresponding pair code to pdg code
                Int_t GetElectronSource(Int_t mclabel); // return origin of the electron

                Double_t GetSecVtxChi2(AliESDtrack* track1, AliESDtrack* track2, AliESDtrack* track3, AliESDtrack* track4); // get secondary vertex chi2 for 4 tracks
                Double_t GetSecVtxChi2(AliESDtrack* track1, AliESDtrack* track2, AliESDtrack* track3); // get secondary vertex chi2 for 3 tracks
                Double_t GetSecVtxChi2(AliESDtrack* track1, AliESDtrack* track2); // get secondary vertex chi2 for 2 tracks


        private:

                AliESDEvent* fESD1; // ESD pointer             
                AliStack* fStack; // stack pointer              

                Int_t fParentSelect[2][7]; //
                Int_t fNparents; // 

                TString fkSourceLabel[10]; //

                enum {fkAll, fkDirectCharm, fkDirectBeauty, fkBeautyCharm, fkGamma, fkPi0, fkElse, fkBeautyGamma, fkBeautyPi0, fkBeautyElse};
                enum {fkCharm=4, fkBeauty=5};

                struct histspair{
                        TH2F *fInvMass; // histogram to store pair invariant mass
                        TH2F *fInvMassCut1; // histogram to store pair invariant mass after cut1
                        TH2F *fInvMassCut2; // histogram to store pair invariant mass after cut2
                        TH1F *fKFChi2; // histogram to store pair vertex chi2
                        TH1F *fCosOpenAngle; // histogram to store pair opening angle
                        TH2F *fSignedLxy; // histogram to store signed Lxy
                        TH1F *fKFIP; // histogram to store pair impact parameter
                };

                struct histstag{
                        TH1F *fPtBeautyElec; // histogram for electrons of single track cut passed 
                        TH1F *fPtTaggedElec; // histogram for total b tagged electrons
                        TH1F *fPtRightTaggedElec; // histogram for right b tagged electrons
                        TH1F *fPtWrongTaggedElec; // histogram for wrong b tagged electrons
                };

                histspair fHistPair[10]; // 
                histstag fHistTagged; //

                Int_t fPairTagged; //
                Int_t fpairedTrackID[20]; //
                Double_t fpairedChi2[20]; //
                Double_t fpairedInvMass[20]; //
                Double_t fpairedSignedLxy[20]; //

                Int_t fid[4][3]; //
                Int_t fia[4][3][2]; //

                Double_t fdistTwoSecVtx; // 
                Double_t fcosPhi; // 
                Double_t fsignedLxy; // signed Lxy of secondary vertex
                Double_t finvmass; // invariant mass of secondary vertex
                Double_t finvmassSigma; // invariant mass sigma of secondary vertex

                Bool_t fBTagged; // if b tagged, set true

		ClassDef(AliHFEsecVtx,0); // secondary vertex for electrons with AliKFParticle package
};

#endif
