// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef ALIEVE_JetPlane_H
#define ALIEVE_JetPlane_H

#include <TEveUtil.h>
#include <TEveElement.h>
#include <TEveTrans.h>

#include <TAtt3D.h>
#include <TAttBBox.h>

#include <AliAODJet.h>
#include <AliAODTrack.h>

#include <vector>


class AliEveJetPlane : public TEveElementList,
		       public TAtt3D,
		       public TAttBBox
{
  friend class AliEveJetPlaneGL;

private:
  AliEveJetPlane(const AliEveJetPlane&);            // Not implemented
  AliEveJetPlane& operator=(const AliEveJetPlane&); // Not implemented

protected:
  Float_t fMinEta;    // Min eta for display.
  Float_t fMaxEta;    // Max eta for display.
  Float_t fMinPhi;    // Min phi for display.
  Float_t fMaxPhi;    // Max phi for display.
  Int_t   fNEtaDiv;   // Number of eta divisions for display.
  Int_t   fNPhiDiv;   // Number of phi divisions for display.

  Float_t fEtaScale;          // Multiplier for eta.
  Float_t fPhiScale;          // Multiplier for phi.
  Float_t fEnergyScale;       // Multiplier for energy.
  Float_t fEnergyColorScale;  // Multiplier for energy color.

  Color_t fGridColor; // Color of coordinate grid.

  TEveTrans  fHMTrans;// Transformation matrix.

  std::vector<AliAODJet>   fJets;   // Jets to display.
  std::vector<AliAODTrack> fTracks; // Tracks to display.

  Bool_t                 fRnrJets;       // Show jets.
  Bool_t                 fRnrTracks;     // Show tracks.
  Bool_t                 fOneSelection;  // One object selected.
  Bool_t                 fTwoSelection;  // Two objects selected.

  AliAODJet             *fJet1,   *fJet2;    // Selection jets.
  AliAODTrack           *fTrack1, *fTrack2;  // Selection tracks.

  Int_t                  fSelectionFlag; // Selection state, handled by GL renderer.

  // Common settings:
  static Bool_t fgOneMomentumXYZ;       // Display momentum as coordinates.
  static Bool_t fgOneMomentumPhiTheta;  // Display momentum as phi/theta.
  static Bool_t fgOneEta;               // Display eta.
  static Bool_t fgOneE;                 // Display energy.
  static Bool_t fgOneChgMass;           // Display charge and mass.

public:
  AliEveJetPlane(Int_t iev);
  virtual ~AliEveJetPlane() {}

  void AddJet(AliAODJet jet);
  void AddTrack(AliAODTrack track);

  Int_t GetNEtaDiv() const  { return fNEtaDiv; }
  void  SetNEtaDiv(Int_t r) { fNEtaDiv = r; }

  Int_t GetNPhiDiv() const  { return fNPhiDiv; }
  void  SetNPhiDiv(Int_t r) { fNPhiDiv = r; }

  Bool_t GetRnrJets() const   { return fRnrJets; }
  void   SetRnrJets(Bool_t r) { fRnrJets = r; }

  Bool_t GetRnrTracks() const   { return fRnrTracks; }
  void   SetRnrTracks(Bool_t r) { fRnrTracks = r; }

  Bool_t GetOneSelection() const   { return fOneSelection; }
  void   SetOneSelection(Bool_t r) { fOneSelection = r; }

  Bool_t GetTwoSelection() const   { return fTwoSelection; }
  void   SetTwoSelection(Bool_t r) { fTwoSelection = r; }

  Float_t GetEnergyScale() const { return fEnergyScale; }
  void    SetEnergyScale(Float_t s) { fEnergyScale = s; }

  Float_t GetEnergyColorScale() const { return fEnergyColorScale; }
  void    SetEnergyColorScale(Float_t s) { fEnergyColorScale = s; }

  const AliAODJet& GetJet1() const { return *fJet1; }
  const AliAODJet& GetJet2() const { return *fJet2; }
  const AliAODTrack& GetTrack1() const { return *fTrack1; }
  const AliAODTrack& GetTrack2() const { return *fTrack2; }

  void    SetJet1(AliAODJet* s) { fJet1 = s; }
  void    SetJet2(AliAODJet* s) { fJet2 = s; }
  void    SetTrack1(AliAODTrack* s) { fTrack1 = s; }
  void    SetTrack2(AliAODTrack* s) { fTrack2 = s; }

  void    SetSelectionFlag(Int_t s) { fSelectionFlag = s;}

  virtual Bool_t  CanEditMainColor()   { return kTRUE; }

  virtual Bool_t  CanEditMainHMTrans() { return kTRUE; }
  virtual TEveTrans* PtrMainHMTrans()  { return &fHMTrans; }

  virtual void ComputeBBox();
  virtual void Paint(Option_t* option = "");

  TEveTrans& RefHMTrans()                    { return fHMTrans; }
  void SetTransMatrix(Double_t* carr)        { fHMTrans.SetFrom(carr); }
  void SetTransMatrix(const TGeoMatrix& mat) { fHMTrans.SetFrom(mat);  }

  ClassDef(AliEveJetPlane, 1); // Show jets and tracks in eta-phi plane.
}; // endclass AliEveJetPlane

#endif
