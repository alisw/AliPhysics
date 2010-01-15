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

#ifndef ALIHFEELECBACKGROUND_H
#define ALIHFEELECBACKGROUND_H

#ifndef ROOT_TObject
//#include <TObject.h>
#endif

class AliESDEvent;
class AliESDVertex;
class AliAODEvent;
class AliESDtrack;
class AliAODTrack;
class AliMCEvent;

//________________________________________________________________
class AliHFEelecbackground : public TObject {

        public: 
                AliHFEelecbackground();
                AliHFEelecbackground(const AliHFEelecbackground &p);
                AliHFEelecbackground &operator=(const AliHFEelecbackground &);
                virtual ~AliHFEelecbackground();
		virtual Bool_t Load(const Char_t *filename);

                void CreateHistograms(TList* const qaList);
		void Reset();

                Bool_t HasMCData() const { return TestBit(kHasMCData); };
                Bool_t IsAODanalysis() const { return TestBit(kAODanalysis); };
                Bool_t IsESDanalysis() const { return !TestBit(kAODanalysis); };

		void SetHasMCData(Bool_t hasMCdata = kTRUE) { SetBit(kHasMCData,hasMCdata); };
                void SetAODAnalysis() { SetBit(kAODanalysis, kTRUE); };
                void SetESDAnalysis() { SetBit(kAODanalysis, kFALSE); };
                void SetEvent(AliESDEvent* const ESD); 
                void SetEventAOD(AliAODEvent* const AOD){fAOD1=AOD;}; 
                void SetMCEvent(AliMCEvent* const mcEvent){fMCEvent=mcEvent;};  

		TList *GetList()  const           { return fList; };
		TList *GetListPostProcess() const { return fListPostProcess; };

		Bool_t SingleTrackCut(AliESDtrack* const track1) const; 
                void PairAnalysis(AliESDtrack* const track, AliESDtrack* const trackpart); 
		void FillOutput(Double_t *results, Double_t *resultsr, Int_t sign); 
                void PostProcess();
		void Plot() const;


 private:

                enum{
                    kHasMCData = BIT(15),     // bitset for mc data usage
		    kAODanalysis = BIT(16)    // bitset for aod analysis
		    };
		enum {kPp=0, kNn=1, kSs=2, kR=3, kOs=4};
		enum {kOos=0, kOss=1, kOr=2, kOdiff=3};  // outputs 
		enum{
		  kElectronFromGamma = 0,
		    kElectronFromPi0 = 1,
		    kElectronFromC = 2,
		    kElectronFromB = 3,
		    kElectronFromEta = 4,
		    kSplittedTrackss =5,
		    kSplittedTrackos =6
		    };
		enum {kNOutput=4, kNSignComb=5, kNMCInfo=7}; 
		
		void CalculateMotherVariable(AliESDtrack* const track, AliESDtrack* const trackpart, Double_t *results);
		void CalculateMotherVariableR(AliESDtrack* const track, AliESDtrack* const trackpart, Double_t *results);
		Int_t IsMotherGamma(Int_t tr);
		Int_t IsMotherPi0(Int_t tr);
		Int_t IsMotherEta(Int_t tr);
		Int_t IsMotherC(Int_t tr);
		Int_t IsMotherB(Int_t tr);
		Int_t GetPdg(Int_t tr);
		Int_t GetLabMother(Int_t tr);


                AliESDEvent* fESD1;              // ESD pointer             
                AliAODEvent* fAOD1;              // AOD pointer             
                AliMCEvent*  fMCEvent;           // MC event             
		Double_t fBz;                    // Magnetic field 
		const AliESDVertex *fkVertex;     // Primary vertex
		static const Double_t    fgkMe;  //  Mass of the electron
		
		Double_t fPtESD;                 // pt of tagged electron
		Int_t fIndexTrack;               // index track
		Int_t fPdg;                      // pdg code track 
		Int_t fLabMother;                // label first mother track 
		Int_t fIsFrom;                   // is track from
		Int_t fMotherGamma;              // Gamma, mother of track
		Int_t fMotherPi0;                // Pi0, mother of track
		Int_t fMotherC;                  // C, mother of track
		Int_t fMotherB;                  // B, mother of track
		Int_t fMotherEta;                // eta, mother of track
		Bool_t fIsPartner;               // Are partners
		Bool_t fIsSplittedTrack;         // Are splitted track
   
		
                TList *fList;                    // list for outputs
		TList *fListPostProcess;         // list for postprocess

		static Bool_t  fgUseMCPID;    // flag to use MC PID for tagged electron
    
    ClassDef(AliHFEelecbackground,0);
};

#endif
