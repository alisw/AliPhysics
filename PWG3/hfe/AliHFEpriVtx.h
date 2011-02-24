#ifndef ALIHFEPRIVTX_H
#define ALIHFEPRIVTX_H

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

/* $Id$ */ 

//
// QA class of primary vertex study for Heavy Flavor electrons
// this has functionality to reject electrons from primary vertex
// and check primary vertex characteristics
//

#ifndef ROOT_TObject
//#include <TObject.h>
#endif

class TH1F;
class TH1I;
class TH2F;
class TString;
class AliESDEvent;
class AliESDtrack;
class AliMCEvent;
class AliKFParticle;

//________________________________________________________________
class AliHFEpriVtx : public TObject {

        public: 
                AliHFEpriVtx();
                AliHFEpriVtx(const AliHFEpriVtx &p); // copy constructor
                AliHFEpriVtx &operator=(const AliHFEpriVtx &); // assignment operator
                virtual ~AliHFEpriVtx();

                void CreateHistograms(TString hnopt=""); // create histograms
                void Init();
                void SetEvent(AliESDEvent * const ESD){fESD1=ESD;}; // set ESD pointer
                void SetMCEvent(AliMCEvent * const mcEvent){fMCEvent=mcEvent;} // set stack pointer
                void CountNtracks(Int_t sourcePart, Int_t recpid, Double_t recprob); // count number of tracks passed certain cut
                void FillNtracks(); // fill counted number of tracks
                void CountPriVxtElecContributor(AliESDtrack *ESDelectron, Int_t sourcePart, Int_t recpid, Double_t recprob); 
                void FillNprimVtxContributor() const;
                void RecalcPrimvtx(Int_t nkftrk, const Int_t * const, const AliKFParticle * const); //recalculate primary vertex after removing given tracks
                void RecalcPrimvtx(AliESDtrack * const ESDelectron); //recalculate primary vertex after removing given track
                void GetRecalcPrimvtx(Double_t privtx[3]) const {
                    privtx[0]=fPVxRe; privtx[1]=fPVyRe; privtx[2]=fPVzRe;
                }
                void GetNPriVxtContributor();
                Double_t GetDistanceFromRecalVertexXY(AliESDtrack * const ESDelectron);
		            Int_t GetNsectrk2prim() const {return fNsectrk2prim;}; 
                Int_t GetMCPID(AliESDtrack const *track); // return mc pid


        private:

                AliESDEvent* fESD1; // ESD event 
                AliMCEvent* fMCEvent; // MC Event

                TString fkSourceLabel[10]; // storing source label

                enum kSources {kAll, kDirectCharm, kDirectBeauty, kBeautyCharm, kGamma, kPi0, kElse, kBeautyGamma, kBeautyPi0, kBeautyElse};

                struct AliHists{
                        TH1F *fNtracks; // histogram to fill number of counted tracks for different sources
                        TH1F *fNprimVtxContributor; // histogram to fill number of tracks contributing primary vertex 
                        TH1F *fPtElec; // histogram to fill pt of electron tracks
                        TH1F *fPtElecContributor; // histogram to fill pt of electron tracks contributing primary vertex
                        Int_t fNtrackCount; // number of counted track
                        Int_t fNprimVtxContributorCount; // number of tracks contributing primary vertex

			AliHists()
			: fNtracks()
			, fNprimVtxContributor()
			, fPtElec()
			, fPtElecContributor()
			, fNtrackCount(0)
			, fNprimVtxContributorCount(0)
			{
			  // default constructor
			}

			AliHists(const AliHists & p)
			: fNtracks(p.fNtracks)
			, fNprimVtxContributor(p.fNprimVtxContributor)
			, fPtElec(p.fPtElec)
			, fPtElecContributor(p.fPtElecContributor)
			, fNtrackCount(p.fNtrackCount)
			, fNprimVtxContributorCount(p.fNprimVtxContributorCount)
			{
			  // copy constructor
			}
			AliHists &operator=(const AliHists &)
			{
			  // assignment operator, not yet implemented
			  return *this;
			}
                };
                AliHists fPrimVtx[10]; // define structure of histograms

                Int_t fNtrackswoPid; //  number of track counted
                TH1F *fHNtrackswoPid; // histogram to fill number of track counted
                TH1I *fNESDprimVtxContributor; // histogram to fill number of primary vertex contributor for given event 
                TH1I *fNESDprimVtxIndices; // histogram to fill number of primary vertex indices for given event
                TH2F *fDiffDCAvsPt; // histogram to fill DCA difference as a function of pT
                TH2F *fDiffDCAvsNt; // histogram to fill DCA difference as a function of pT

                Int_t fNsectrk2prim; // # of secvtx tracks contributing to primvtx calculation
                Double_t fPVxRe;     // recalculated primary vertex x 
                Double_t fPVyRe;     // recalculated primary vertex y 
                Double_t fPVzRe;     // recalculated primary vertex z  

        ClassDef(AliHFEpriVtx,0);
};

#endif
