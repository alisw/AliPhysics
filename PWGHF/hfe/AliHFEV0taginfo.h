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
#ifndef ALIHFEV0TAGINFO_H
#define ALIHFEV0TAGINFO_H

#ifndef ROOT_TNamed
#include <TNamed.h>
#endif


#include "AliPID.h"
#include "AliESDv0KineCuts.h"
#include "AliAODv0KineCuts.h"

class AliVEvent;

class AliHFEV0taginfo: public TNamed{

    public:

        AliHFEV0taginfo();
        AliHFEV0taginfo(const char* name);
        AliHFEV0taginfo(const AliHFEV0taginfo &ref);
        AliHFEV0taginfo& operator=(const AliHFEV0taginfo &ref);
        virtual ~AliHFEV0taginfo();

        void Reset(); //resets the fTaggedTracks TList
        void TagV0Tracks(AliVEvent *event); // loops over V0s in event and fills fTaggedTracks with V0 tracks
        AliPID::EParticleType GetV0Info(Int_t trackID); //check for V0 information from track ID 
        Float_t GetV0ProdR(Int_t trackID); //check for V0 information from track ID 
        void SetIsAODana() { fIsAODana = kTRUE; } // Setter AOD analysis
        void SetIsESDana() { fIsAODana = kFALSE; } // Setter ESD analysis
        void AddTrack(Int_t TrackID, Int_t pdgCode, Double_t prodR); //Translates the pdg code to AliPID enum and adds track to tagged list (currently only e, pion and proton tracks)

        class AliHFEV0tag: public TObject{
            public:
                AliHFEV0tag();
                AliHFEV0tag(Int_t TrackID, AliPID::EParticleType Pinfo, Double_t ProdR);
                AliHFEV0tag(const AliHFEV0tag &ref);
                AliHFEV0tag &operator=(const AliHFEV0tag &ref);
                virtual ~AliHFEV0tag();

                //TObject virtual functions
                virtual Bool_t IsSortable() const { return kTRUE; };
                virtual Bool_t IsEqual(const TObject *ref) const;
                virtual Int_t Compare(const TObject *ref) const;

                //Setter
                void SetTrack(Int_t trackID, AliPID::EParticleType Pinfo); //Set track ID
                void SetProdR(Int_t trackID, Double_t prodR); //Set V0 daughter production vertex radius
                //Getter
                Int_t GetTrackID() const { return fTrackID; };
                AliPID::EParticleType GetPinfo() const { return fPinfo; };
                Double_t GetProdR() const { return fProdR; };
            private:
                Int_t fTrackID;
                AliPID::EParticleType fPinfo;  
                Double_t fProdR; // V0 daughter production vertex radius in x-y plane

                ClassDef(AliHFEV0taginfo::AliHFEV0tag,1);
        };
    private:
        Bool_t fIsAODana;
        TList *fTaggedTracks;
        AliESDv0KineCuts *fV0finder;
        AliAODv0KineCuts *fAODV0finder;

        ClassDef(AliHFEV0taginfo,1)
};
#endif
