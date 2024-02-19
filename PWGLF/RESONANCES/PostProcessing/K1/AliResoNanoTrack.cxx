/*************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

/* AliResoNanoTrack.h
 * Header file for the AliResoNanoTrack class, which is a simple class
 * to store track information in a reduced format
 *
 * Author: Bong-Hwi Lim
 *
 */
#include "AliResoNanoTrack.h"
class AliResoNanoTrack;

AliResoNanoTrack::AliResoNanoTrack() : TLorentzVector(),
                                       fID(-1),
                                       fParticleType(-999),
                                       fCharge(-999),
                                       fLabel(-1),
                                       fStatus(-999),
                                       fFlags(-999),
                                       fDCAxy(-999),
                                       fDCAz(-999),
                                       fTOFNSigmaPi(-999),
                                       fTPCNSigmaPi(-999),
                                       fTOFNSigmaKa(-999),
                                       fTPCNSigmaKa(-999),
                                       fTOFNSigmaPr(-999),
                                       fTPCNSigmaPr(-999),
                                       fMCPDGCode(-999),
                                       fMCMotherPDGCode(-999),
                                       fMCMotherID(-1) { ; } // default constructor

AliResoNanoTrack::AliResoNanoTrack(float px, float py, float pz, float e, Int_t id, Short_t ptype, Char_t charge) : TLorentzVector(px, py, pz, e),
                                                                                                                    fID(id),
                                                                                                                    fParticleType(ptype),
                                                                                                                    fCharge(charge),
                                                                                                                    fLabel(-1),
                                                                                                                    fStatus(-999),
                                                                                                                    fFlags(-999),
                                                                                                                    fDCAxy(-999),
                                                                                                                    fDCAz(-999),
                                                                                                                    fTOFNSigmaPi(-999),
                                                                                                                    fTPCNSigmaPi(-999),
                                                                                                                    fTOFNSigmaKa(-999),
                                                                                                                    fTPCNSigmaKa(-999),
                                                                                                                    fTOFNSigmaPr(-999),
                                                                                                                    fTPCNSigmaPr(-999),
                                                                                                                    fMCPDGCode(-999),
                                                                                                                    fMCMotherPDGCode(-999),
                                                                                                                    fMCMotherID(-1) { ; } // constructor

AliResoNanoTrack::AliResoNanoTrack(const AliResoNanoTrack &a) : TLorentzVector(a.Px(), a.Py(), a.Pz(), a.E()),
                                                                fID(a.GetID()),
                                                                fParticleType(a.GetParticleType()),
                                                                fCharge(a.GetCharge()),
                                                                fLabel(a.GetLabel()),
                                                                fStatus(a.GetStatus()),
                                                                fFlags(a.GetFlags()),
                                                                fDCAxy(a.GetDCAxy()),
                                                                fDCAz(a.GetDCAz()),
                                                                fTOFNSigmaPi(a.GetTOFNSigmaPi()),
                                                                fTPCNSigmaPi(a.GetTPCNSigmaPi()),
                                                                fTOFNSigmaKa(a.GetTOFNSigmaKa()),
                                                                fTPCNSigmaKa(a.GetTPCNSigmaKa()),
                                                                fTOFNSigmaPr(a.GetTOFNSigmaPr()),
                                                                fTPCNSigmaPr(a.GetTPCNSigmaPr()),
                                                                fMCPDGCode(a.GetMCPDGCode()),
                                                                fMCMotherPDGCode(a.GetMCMotherPDGCode()),
                                                                fMCMotherID(a.GetMCMotherID()) { ; } // copy constructor

AliResoNanoTrack::AliResoNanoTrack(const TLorentzVector &a) : TLorentzVector(a.Px(), a.Py(), a.Pz(), a.E()),
                                                              fID(-1),
                                                              fParticleType(-999),
                                                              fCharge(-999),
                                                              fLabel(-1),
                                                              fStatus(-999),
                                                              fFlags(-999),
                                                              fDCAxy(-999),
                                                              fDCAz(-999),
                                                              fTOFNSigmaPi(-999),
                                                              fTPCNSigmaPi(-999),
                                                              fTOFNSigmaKa(-999),
                                                              fTPCNSigmaKa(-999),
                                                              fTOFNSigmaPr(-999),
                                                              fTPCNSigmaPr(-999),
                                                              fMCPDGCode(-999),
                                                              fMCMotherPDGCode(-999),
                                                              fMCMotherID(-1) { ; } // copy constructor

AliResoNanoTrack::~AliResoNanoTrack() { ; } // destructor
