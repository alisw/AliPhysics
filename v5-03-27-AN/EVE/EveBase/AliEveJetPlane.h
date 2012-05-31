// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveJetPlane_H
#define AliEveJetPlane_H

#include <TEveElement.h>

#include <TAtt3D.h>
#include <TAttBBox.h>

#include <AliAODJet.h>
#include <AliAODTrack.h>

#include <vector>

//==============================================================================
//
// AliEveJetPlane
//
// Class for display of jets and tracks in eta-phi plane.

class AliEveJetPlane : public TEveElementList,
		       public TAtt3D,
		       public TAttBBox
{
  friend class AliEveJetPlaneGL;

public:
  AliEveJetPlane(Int_t iev);
  virtual ~AliEveJetPlane();

  void AddJet(AliAODJet* jet);
  void AddTrack(AliAODTrack* track);

  void CreateArrows();

  Int_t GetNEtaDiv() const  { return fNEtaDiv; }
  void  SetNEtaDiv(Int_t r) { fNEtaDiv = r; }

  Int_t GetNPhiDiv() const  { return fNPhiDiv; }
  void  SetNPhiDiv(Int_t r) { fNPhiDiv = r; }

  Bool_t GetRnrJets() const   { return fRnrJets; }
  void   SetRnrJets(Bool_t r) { fRnrJets = r; CreateArrows(); }

  Bool_t GetRnrTracks() const   { return fRnrTracks; }
  void   SetRnrTracks(Bool_t r) { fRnrTracks = r; CreateArrows(); }

  Bool_t GetOneSelection() const   { return fOneSelection; }
  void   SetOneSelection(Bool_t r) { fOneSelection = r; }

  Bool_t GetTwoSelection() const   { return fTwoSelection; }
  void   SetTwoSelection(Bool_t r) { fTwoSelection = r; }

  Float_t GetEnergyScale() const { return fEnergyScale; }
  void    SetEnergyScale(Float_t s) { fEnergyScale = s; CreateArrows(); }

  Float_t GetArrowJetScale() const { return fArrowJetScale; }
  void    SetArrowJetScale(Float_t s) { fArrowJetScale = s; CreateArrows(); }
	
  Float_t GetArrowTrackScale() const { return fArrowTrackScale; }
  void    SetArrowTrackScale(Float_t s) { fArrowTrackScale = s; CreateArrows(); }

  const AliAODJet& GetJet1() const { return *fJet1; }
  const AliAODJet& GetJet2() const { return *fJet2; }
  const AliAODTrack& GetTrack1() const { return *fTrack1; }
  const AliAODTrack& GetTrack2() const { return *fTrack2; }

  void    SetJet1(AliAODJet* s) { fJet1 = s; }
  void    SetJet2(AliAODJet* s) { fJet2 = s; }
  void    SetTrack1(AliAODTrack* s) { fTrack1 = s; }
  void    SetTrack2(AliAODTrack* s) { fTrack2 = s; }

  void    SetSelectionFlag(Int_t s) { fSelectionFlag = s;}
  void    SelectionAdded(TEveElement* el);

  virtual Bool_t  CanEditMainColor()const { return kTRUE; }

  virtual void ComputeBBox();
  virtual void Paint(Option_t* option = "");

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
  Float_t fArrowJetScale;     // Multiplier for jet arrow dim.
  Float_t fArrowTrackScale;   // Multiplier for track arrow dim.
	
  Color_t fGridColor; // Color of coordinate grid.

  std::vector<AliAODJet>   fJets;   // Jets to display.
  std::vector<AliAODTrack> fTracks; // Tracks to display.

  Bool_t                 fRnrJets;       // Show jets.
  Bool_t                 fRnrTracks;     // Show tracks.
  Bool_t                 fOneSelection;  // One object selected.
  Bool_t                 fTwoSelection;  // Two objects selected.
  Bool_t                 fSelConnected;  // Connected to EVE selection.

  AliAODJet             *fJet1,   *fJet2;    // Selection jets.
  AliAODTrack           *fTrack1, *fTrack2;  // Selection tracks.

  Int_t                  fSelectionFlag; // Selection state, handled by GL renderer.

  // Common settings:
  static Bool_t fgOneMomentumXYZ;       // Display momentum as coordinates.
  static Bool_t fgOneMomentumPhiTheta;  // Display momentum as phi/theta.
  static Bool_t fgOneEta;               // Display eta.
  static Bool_t fgOneE;                 // Display energy.
  static Bool_t fgOneChgMass;           // Display charge and mass.

private:
  AliEveJetPlane(const AliEveJetPlane&);            // Not implemented
  AliEveJetPlane& operator=(const AliEveJetPlane&); // Not implemented

  ClassDef(AliEveJetPlane, 0); // Show jets and tracks in eta-phi plane.
};

#endif
