#ifndef ALIHFESECVTXS_H
#define ALIHFESECVTXS_H

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
//  Secondary vertexing container to store secondary vertex characteristics of 
//  2 or 3 particle sec vertex
//  from example, qusi-invariant mass, signed Lxy are stored
//

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
                Int_t GetMCCode() const {return fMCCode;}
                Double_t GetInvmass() const {return fInvmass;}
                Double_t GetKFChi2() const {return fKFChi2;}
                Double_t GetSignedLxy() const {return fSignedLxy;}
                Double_t GetSignedLxy2() const {return fSignedLxy2;}
                Double_t GetKFIP() const {return fKFIP;}
                Double_t GetKFIP2() const {return fKFIP2;}

                void SetTrkLabel1(Int_t label) {fTrkLabel1 = label;}
                void SetTrkLabel2(Int_t label) {fTrkLabel2 = label;}
                void SetInvmass(Double_t invmass) {fInvmass = invmass;}
                void SetKFChi2(Double_t kfchi2) {fKFChi2 = kfchi2;}
                void SetSignedLxy(Double_t signedlxy) {fSignedLxy = signedlxy;}
                void SetSignedLxy2(Double_t signedlxy2) {fSignedLxy2 = signedlxy2;}
                void SetKFIP(Double_t kfip) {fKFIP = kfip;}
                void SetKFIP2(Double_t kfip2) {fKFIP2 = kfip2;}
                void SetMCCode(Int_t mccode) {fMCCode = mccode;}

        protected:
                Int_t fTrkLabel1;    // track 1 label associated to secvtx 
                Int_t fTrkLabel2;    // track 2 label associated to secvtx
                Int_t fMCCode;       // track mc code
                Double_t fInvmass;   // secvtx invariant mass
                Double_t fKFChi2;    // secvtx chi2 
                Double_t fSignedLxy; // secvtx signed Lxy
                Double_t fSignedLxy2; // recalculated secvtx signed Lxy
                Double_t fKFIP;       // secvtx impact parameter 
                Double_t fKFIP2;      // recalculated secvtx impact parameter 

        private:

        ClassDef(AliHFEsecVtxs,1);
};

#endif
