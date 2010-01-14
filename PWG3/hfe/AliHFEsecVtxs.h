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
//  Secondary vertexing container to store secondary vertex characteristics of 
//  2 or 3 particle sec vertex
//  from example, qusi-invariant mass, signed Lxy are stored
//

#ifndef ALIHFESECVTXS_H
#define ALIHFESECVTXS_H

#ifndef ROOT_TObject
#include <TObject.h>
#endif

//________________________________________________________________
class AliHFEsecVtxs : public TObject {

        public: 
                AliHFEsecVtxs();
                AliHFEsecVtxs(const AliHFEsecVtxs &p); // copy constructor
                AliHFEsecVtxs &operator=(const AliHFEsecVtxs &); // assignment operator
                virtual ~AliHFEsecVtxs();

                Int_t GetTrkLabel1() const {return fTrkLabel1;}
                Int_t GetTrkLabel2() const {return fTrkLabel2;}
                Double_t GetInvmass() const {return fInvmass;}
                Double_t GetKFChi2() const {return fKFChi2;}
                Double_t GetSignedLxy() const {return fSignedLxy;}
                Double_t GetKFIP() const {return fKFIP;}

                void SetTrkLabel1(Int_t label) {fTrkLabel1 = label;}
                void SetTrkLabel2(Int_t label) {fTrkLabel2 = label;}
                void SetInvmass(Double_t invmass) {fInvmass = invmass;}
                void SetKFChi2(Double_t kfchi2) {fKFChi2 = kfchi2;}
                void SetSignedLxy(Double_t signedlxy) {fSignedLxy = signedlxy;}
                void SetKFIP(Double_t kfip) {fKFIP = kfip;}

        protected:
                Int_t fTrkLabel1;    // track 1 label associated to secvtx 
                Int_t fTrkLabel2;    // track 2 label associated to secvtx
                Double_t fInvmass;   // secvtx invariant mass
                Double_t fKFChi2;    // secvtx chi2 
                Double_t fSignedLxy; // secvtx signed Lxy
                Double_t fKFIP;      // secvtx impact parameter 

        private:

        ClassDef(AliHFEsecVtxs,0);
};

#endif
