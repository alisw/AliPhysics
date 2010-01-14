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
//  Container class to store pair characteristics
//  for secondary vertex analysis
//  from example, qusi-invariant mass, signed Lxy are stored
//

#ifndef ALIHFEPAIRS_H
#define ALIHFEPAIRS_H

#ifndef ROOT_TObject
#include <TObject.h>
#endif

//________________________________________________________________
class AliHFEpairs : public TObject {

        public: 
                AliHFEpairs();
                AliHFEpairs(const AliHFEpairs &p); // copy constructor
                AliHFEpairs &operator=(const AliHFEpairs &); // assignment operator
                virtual ~AliHFEpairs();

                Int_t GetTrkLabel() const {return fTrkLabel;}
                Int_t GetPairCode() const {return fPairCode;}
                Double_t GetInvmass() const {return fInvmass;}
                Double_t GetKFChi2() const {return fKFChi2;}
                Double_t GetOpenangle() const {return fOpenangle;}
                Double_t GetCosOpenangle() const {return fCosOpenangle;}
                Double_t GetSignedLxy() const {return fSignedLxy;}
                Double_t GetKFIP() const {return fKFIP;}

                void SetTrkLabel(Int_t label) {fTrkLabel = label;}
                void SetInvmass(Double_t invmass) {fInvmass = invmass;}
                void SetKFChi2(Double_t kfchi2) {fKFChi2 = kfchi2;}
                void SetOpenangle(Double_t openangle) {fOpenangle = openangle;}
                void SetCosOpenangle(Double_t cosopenangle) {fCosOpenangle = cosopenangle;}
                void SetSignedLxy(Double_t signedlxy) {fSignedLxy = signedlxy;}
                void SetKFIP(Double_t kfip) {fKFIP = kfip;}
                void SetPairCode(Int_t paircode) {fPairCode = paircode;}

        protected:

                Int_t fTrkLabel;        // paired track label
                Int_t fPairCode;        // paired track mc code
                Double_t fInvmass;      // pair invariant mass 
                Double_t fKFChi2;       // pair kf vertex chi2 
                Double_t fOpenangle;    // pair opening angle 
                Double_t fCosOpenangle; // pair cos(opening angle)
                Double_t fSignedLxy;    // pair signed Lxy
                Double_t fKFIP;         // impact parameter of the pair

        private:

        ClassDef(AliHFEpairs,0);
};

#endif
