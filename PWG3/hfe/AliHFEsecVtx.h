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
//
//  Secondary vertexing construction Class
//  Construct secondary vertex from Beauty hadron with electron and
//  hadrons, then apply selection criteria
//

#ifndef ALIHFESECVTX_H
#define ALIHFESECVTX_H

#ifndef ROOT_TObject
//#include <TObject.h>
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

                void CreateHistograms(TString hnopt="");

                void SetEvent(AliESDEvent* const ESD){fESD1=ESD;}; // set ESD pointer
                void SetStack(AliStack* const stack){fStack=stack;} // set stack pointer
                void InitAnaPair(); // initialize default parameters 
                void AnaPair(AliESDtrack* ESDtrack1, AliESDtrack* ESDtrack2, Int_t index2); // do e-h analysis
                void RunSECVTX(AliESDtrack *track); // run secondary vertexing algorithm

                Int_t GetMCPID(AliESDtrack *track); // return MC pid
                Int_t PairOrigin(AliESDtrack* track1, AliESDtrack* track2); // return pair origin as a pdg code
                Int_t PairCode(AliESDtrack* track1,AliESDtrack* track2); // return corresponding pair code to pdg code
                Int_t GetElectronSource(Int_t mclabel); // return origin of the electron

                Bool_t SingleTrackCut(AliESDtrack* track1) const; // single track cut

        protected:

                void Init();
                void FindSECVTXCandid4Tracks(AliESDtrack* track1); // find secondary vertex for 4 tracks
                void FindSECVTXCandid3Tracks(AliESDtrack* track1); // find secondary vertex for 3 tracks        
                void FindSECVTXCandid2Tracks(AliESDtrack* track1); // find secondary vertex for 2 tracks
                void CalcSECVTXProperty(AliESDtrack* track1, AliESDtrack* track2, AliESDtrack* track3); // calculated distinctive variables
                void ApplyPairTagCut(); // apply strict tagging condition for 1 pair secondary vertex
                void ApplyVtxTagCut(Double_t chi2, Int_t id1, Int_t id2); // apply tagging condition for 3 track secondary vertex
                void ResetTagVar(); // reset tagging variables

                Bool_t ApplyPairTagCut(Int_t id); // apply strict tagging condition for 1 pair secondary vertex
                Bool_t ApplyTagCut(Double_t chi2); // apply tagging condition


                Double_t GetSecVtxChi2(AliESDtrack* track1, AliESDtrack* track2, AliESDtrack* track3, AliESDtrack* track4); // get secondary vertex chi2 for 4 tracks
                Double_t GetSecVtxChi2(AliESDtrack* track1, AliESDtrack* track2, AliESDtrack* track3); // get secondary vertex chi2 for 3 tracks
                Double_t GetSecVtxChi2(AliESDtrack* track1, AliESDtrack* track2); // get secondary vertex chi2 for 2 tracks


        private:

                AliESDEvent* fESD1; // ESD pointer             
                AliStack* fStack; // stack pointer              

                Int_t fParentSelect[2][7]; // heavy hadron species
                Int_t fNparents; // number of heavy hadrons to be considered

                TString fkSourceLabel[10]; // electron source label

                enum {kAll, kDirectCharm, kDirectBeauty, kBeautyCharm, kGamma, kPi0, kElse, kBeautyGamma, kBeautyPi0, kBeautyElse};
                enum {kCharm=4, kBeauty=5};

                struct AliHistspair{
                        TH2F *fInvMass; // histogram to store pair invariant mass
                        TH2F *fInvMassCut1; // histogram to store pair invariant mass after cut1
                        TH2F *fInvMassCut2; // histogram to store pair invariant mass after cut2
                        TH1F *fKFChi2; // histogram to store pair vertex chi2
			TH1F *fOpenAngle; // histogram to store pair opening angle
                        TH1F *fCosOpenAngle; // histogram to store cosine of the pair opening angle
                        TH2F *fSignedLxy; // histogram to store signed Lxy
                        TH1F *fKFIP; // histogram to store pair impact parameter
                        TH1F *fIPMax; // histogram to store larger impact parameter among pair tracks

			AliHistspair()
			: fInvMass()
			, fInvMassCut1()
			, fInvMassCut2()
			, fKFChi2()
			, fOpenAngle()
			, fCosOpenAngle()
			, fSignedLxy()
			, fKFIP()
			, fIPMax()
			{
			  // define constructor
			}
			AliHistspair(const AliHistspair & p)
			: fInvMass(p.fInvMass)
			, fInvMassCut1(p.fInvMassCut1)
			, fInvMassCut2(p.fInvMassCut2)
			, fKFChi2(p.fKFChi2)
			, fOpenAngle(p.fOpenAngle)
			, fCosOpenAngle(p.fCosOpenAngle)
			, fSignedLxy(p.fSignedLxy)
			, fKFIP(p.fKFIP)
			, fIPMax(p.fIPMax)
			{
			  // copy constructor
			}
			AliHistspair &operator=(const AliHistspair &)
			{
			  // assignment operator, not yet implemented
			  return *this;
			}
                };

                struct AliHiststag{
                        TH1F *fPtBeautyElec; // histogram for electrons of single track cut passed 
                        TH1F *fPtTaggedElec; // histogram for total b tagged electrons
                        TH1F *fPtRightTaggedElec; // histogram for right b tagged electrons
                        TH1F *fPtWrongTaggedElec; // histogram for wrong b tagged electrons
                        TH1F *fInvmassBeautyElecSecVtx;  // histogram for right-tagged b invariant mass
                        TH1F *fInvmassNotBeautyElecSecVtx; // histogram for mis-tagged b invariant mass
                        TH1F *fSignedLxyBeautyElecSecVtx; // histogram for right-tagged b signed Lxy
                        TH1F *fSignedLxyNotBeautyElecSecVtx; // histogram for mis-tagged b signed Lxy
                        TH1F *fDistTwoVtxBeautyElecSecVtx; // histogram for right-tagged b two vertex distance
                        TH1F *fDistTwoVtxNotBeautyElecSecVtx; // histogram for mis-tagged b two vertex distance
                        TH1F *fCosPhiBeautyElecSecVtx; // histogram for right-tagged b cos of opening angle
                        TH1F *fCosPhiNotBeautyElecSecVtx; // histogram for mis-tagged b cos of opening angle
                        TH1F *fChi2BeautyElecSecVtx; // histogram for right-tagged b chi2 of secondary vertex 
                        TH1F *fChi2NotBeautyElecSecVtx; // histogram for mis-tagged b chi2 of secondary vertex
                        TH1F *fInvmassBeautyElec2trkVtx; // histogram for right-tagged b invariant mass of two track vertex
                        TH1F *fInvmassNotBeautyElec2trkVtx; // histogram for mis-tagged b invariant mass of two track vertex
                        TH1F *fSignedLxyBeautyElec2trkVtx; // histogram for right-tagged b signed Lxy of two track vertex
                        TH1F *fSignedLxyNotBeautyElec2trkVtx; // histogram for mis-tagged b signed Lxy of two track vertex
                        TH1F *fChi2BeautyElec2trkVtx; // histogram for right-tagged b chi2 of two track vertex
                        TH1F *fChi2NotBeautyElec2trkVtx; // histogram for mis-tagged b chi2 of two track vertex

			AliHiststag()
			: fPtBeautyElec()
			, fPtTaggedElec()
			, fPtRightTaggedElec()
			, fPtWrongTaggedElec()
			, fInvmassBeautyElecSecVtx()
			, fInvmassNotBeautyElecSecVtx()
			, fSignedLxyBeautyElecSecVtx()
			, fSignedLxyNotBeautyElecSecVtx()
			, fDistTwoVtxBeautyElecSecVtx()
			, fDistTwoVtxNotBeautyElecSecVtx()
			, fCosPhiBeautyElecSecVtx()
			, fCosPhiNotBeautyElecSecVtx()
			, fChi2BeautyElecSecVtx()
			, fChi2NotBeautyElecSecVtx()
			, fInvmassBeautyElec2trkVtx()
			, fInvmassNotBeautyElec2trkVtx()
			, fSignedLxyBeautyElec2trkVtx()
			, fSignedLxyNotBeautyElec2trkVtx()
			, fChi2BeautyElec2trkVtx()
			, fChi2NotBeautyElec2trkVtx()
			{
			  // define constructor
			}
			AliHiststag(const AliHiststag & p)
			: fPtBeautyElec(p.fPtBeautyElec)
			, fPtTaggedElec(p.fPtTaggedElec)
			, fPtRightTaggedElec(p.fPtRightTaggedElec)
			, fPtWrongTaggedElec(p.fPtWrongTaggedElec)
			, fInvmassBeautyElecSecVtx(p.fInvmassBeautyElecSecVtx)
			, fInvmassNotBeautyElecSecVtx(p.fInvmassNotBeautyElecSecVtx)
			, fSignedLxyBeautyElecSecVtx(p.fSignedLxyBeautyElecSecVtx)
			, fSignedLxyNotBeautyElecSecVtx(p.fSignedLxyNotBeautyElecSecVtx)
			, fDistTwoVtxBeautyElecSecVtx(p.fDistTwoVtxBeautyElecSecVtx)
			, fDistTwoVtxNotBeautyElecSecVtx(p.fDistTwoVtxNotBeautyElecSecVtx)
			, fCosPhiBeautyElecSecVtx(p.fCosPhiBeautyElecSecVtx)
			, fCosPhiNotBeautyElecSecVtx(p.fCosPhiNotBeautyElecSecVtx)
			, fChi2BeautyElecSecVtx(p.fChi2BeautyElecSecVtx)
			, fChi2NotBeautyElecSecVtx(p.fChi2NotBeautyElecSecVtx)
			, fInvmassBeautyElec2trkVtx(p.fInvmassBeautyElec2trkVtx)
			, fInvmassNotBeautyElec2trkVtx(p.fInvmassNotBeautyElec2trkVtx)
			, fSignedLxyBeautyElec2trkVtx(p.fSignedLxyBeautyElec2trkVtx)
			, fSignedLxyNotBeautyElec2trkVtx(p.fSignedLxyNotBeautyElec2trkVtx)
			, fChi2BeautyElec2trkVtx(p.fChi2BeautyElec2trkVtx)
			, fChi2NotBeautyElec2trkVtx(p.fChi2NotBeautyElec2trkVtx)
			{
			  // copy constructor
			}
			AliHiststag &operator=(const AliHiststag &)
			{
			  // assignment operator, not yet implemented
			  return *this;
			}
                };

                AliHistspair fHistPair[10]; // struct of above given histogram for different electron sources
                AliHiststag fHistTagged; // struct of above given histogram

                Int_t fPairTagged; // number of tagged e-h pairs
                Int_t fpairedTrackID[20]; // paird hadron track track 
                Double_t fpairedChi2[20]; // pair chi2
                Double_t fpairedInvMass[20]; // pair invariant mass
                Double_t fpairedSignedLxy[20]; // pair signed Lxy

                Int_t fid[4][3]; // index to store sorting result
                Int_t fia[4][3][2]; // index to store sorting result

                Double_t fdistTwoSecVtx; // distance between two pair vertex
                Double_t fcosPhi; // cos of opening angle of two pair vertex
                Double_t fsignedLxy; // signed Lxy of secondary vertex
                Double_t finvmass; // invariant mass of secondary vertex
                Double_t finvmassSigma; // invariant mass sigma of secondary vertex

                Bool_t fBTagged; // if b tagged, set true
                Bool_t fBElec; // if set true for b electron, set true

        ClassDef(AliHFEsecVtx,0);
};

#endif
