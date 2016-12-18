#ifndef ALIANALYSISTASKEMCALJETCDF_H
#define ALIANALYSISTASKEMCALJETCDF_H
/// \file AliAnalysisTaskEmcalJetCDF.h
/// \brief Declaration of class AliAnalysisTaskEmcalJetCDF
///
/// \author Adrian SEVCENCO <Adrian.Sevcenco@cern.ch>, Institute of Space Science, Romania
/// \date Mar 23, 2015

/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TH1D.h>
#include <TH2D.h>
#include <TArrayD.h>
#include <TVector2.h>
#include "AliAnalysisTaskEmcalJet.h"
#include <THistManager.h>

/// \class AliAnalysisTaskEmcalJetCDF
/// \brief Analysis of jet shapes and FF of all jets and leading jets
///
/// Analysis of leading jets distribution of pt and multiplicity, R distribution
/// N80 and Pt80 and Toward, Transverse, Away UE histos
class AliAnalysisTaskEmcalJetCDF : public AliAnalysisTaskEmcalJet {
public:

    AliAnalysisTaskEmcalJetCDF();
    AliAnalysisTaskEmcalJetCDF ( const char *name );
    virtual ~AliAnalysisTaskEmcalJetCDF();

    void                        UserCreateOutputObjects();
    void                        Terminate ( Option_t *option );

    THistManager                fHistManager   ;///< Histogram manager

protected:
    void                        ExecOnce();
    Bool_t                      Run() ;

    /// Filling of histograms
    /// \return kTRUE if filling is succesful
    Bool_t                      FillHistograms()   ;

    /// Get pointer to a histogram from manager
    /// \param histName Name of the histogram
    /// \return TObject*
    TObject* GetHistogram ( const char* histName );

private:
    AliAnalysisTaskEmcalJetCDF ( const AliAnalysisTaskEmcalJetCDF& );           // not implemented
    AliAnalysisTaskEmcalJetCDF &operator= ( const AliAnalysisTaskEmcalJetCDF& ); // not implemented

    /// \cond CLASSIMP
    ClassDef ( AliAnalysisTaskEmcalJetCDF, 8 );
    /// \endcond

};

namespace NS_AliAnalysisTaskEmcalJetCDF {
  /// (pt,index) pair
  typedef std::pair<Double_t, Int_t> ptidx_pair;

  /// functional for sorting pair by first element - descending
  struct sort_descend
    {
    bool operator () ( const ptidx_pair &p1, const ptidx_pair &p2 )  { return p1.first > p2.first ; }
    };

  /// Computing the Mag2 of an AliVParticle derived object
  /// \param trk AliVParticle
  /// \return Mag2
  inline Double_t Mag2 (const AliVParticle& trk) 
    { return trk.Px()*trk.Px() + trk.Py()*trk.Py() + trk.Pz()*trk.Pz(); }

  /// Computing the Mag of an AliVParticle derived object
  /// \param trk AliVParticle
  /// \return Mag2
  inline Double_t Mag(const AliVParticle& trk) { return TMath::Sqrt(Mag2(trk)); }

  /// Computing the Dot product of two AliVParticle derived objects
  /// \param trk1 AliVParticle 1
  /// \param trk2 AliVParticle 2
  /// \return Mag2
  inline Double_t Dot  (const AliVParticle& trk1, const AliVParticle& trk2 )
    { return trk1.Px()*trk2.Px() + trk1.Py()*trk2.Py() + trk1.Pz()*trk2.Pz(); }

  /// Computing the transverse component squared of trk1 w.r.t trk2
  /// \param trk1 AliVParticle 1
  /// \param trk2 AliVParticle 2
  /// \return Perp2 of trk1 w.r.t. trk2
  inline Double_t Perp2( const AliVParticle& trk1, const AliVParticle& trk2) {
    Double_t mag1  = Mag2(trk1);
    Double_t mag2  = Mag2(trk2);
    Double_t dotp  = Dot(trk1,trk2);
    if (mag2 > 0.0) { mag1 -= dotp*dotp/mag2; }
    if (mag1 <= 0)  { mag1 = 0; }
    return mag1;
    }

  /// Computing the transverse component of trk1 w.r.t trk2
  /// \param trk1 AliVParticle 1
  /// \param trk2 AliVParticle 2
  /// \return Perp of trk1 w.r.t. trk2
  inline Double_t Perp (const AliVParticle& trk1, const AliVParticle& trk2 ) { return TMath::Sqrt(Perp2(trk1,trk2)); }

  /// Sorting of tracks in the event by pt (descending)
  /// \param event AliVEvent*
  /// \return vector of indexes of constituents
  std::vector<Int_t>          SortTracksPt ( AliVEvent* event );

  /// Sorting of tracks in the event by pt (descending) - using a particle container
  /// \param track_container AliParticleContainer*
  /// \return vector of indexes of constituents
  std::vector<Int_t>          SortTracksPt ( AliParticleContainer* track_container );

  /// Get P() fraction of constituent to jet
  /// \param jet AliEmcalJet*
  /// \param trk AliVParticle*
  /// \return Z = trk->P()/jet->P()
  inline Double_t             Z_ptot ( const AliEmcalJet* jet, const AliVParticle* trk ) // Get Z of constituent trk ; p total
    {
    if (trk->P() < 1e-6) return 0.;
    return (trk != 0) ? trk->P()/ jet->P() : 0.;
    }

  /// Get Pt() fraction of constituent to jet
  /// \param jet AliEmcalJet*
  /// \param trk AliVParticle*
  /// \return Z = trk->Pt()/jet->Pt()
  inline Double_t             Z_pt   ( const AliEmcalJet* jet, const AliVParticle* trk ) // Get Z of constituent trk ; pt
    {
    if (trk->P() < 1e-6) return 0.;
    return (trk != 0) ? trk->Pt() / jet->Pt() : 0.;
    }

  /// Get Xi for a double value z
  /// \return Xi
  inline Double_t              Xi ( Double_t z ) { return TMath::Log ( 1/z ); } // Get Xi of value z

  /// Return dR dinstance in eta,phi plane between 2 AliVParticle derived objects
  /// \param part1 AliVParticle*
  /// \param part2 AliVParticle*
  /// \return #eta-#phi distance from jet axis
  inline Double_t              DeltaR ( const AliVParticle* part1, const AliVParticle* part2 )
    {
    Double_t dPhi = part1->Phi() - part2->Phi();
    Double_t dEta = part1->Eta() - part2->Eta();
    dPhi = TVector2::Phi_mpi_pi ( dPhi );
    return TMath::Sqrt ( dPhi * dPhi + dEta * dEta );
    }

  /// Search for index(int) in array of ints
  /// \param index the index to be searched
  /// \param array array of ints
  /// \return kTRUE if found
  inline Int_t                IdxInArray ( Int_t index, TArrayI &array )
    {
    for ( Int_t i = 0; i < array.GetSize(); i++ ) { if ( index == array[i] ) { return i; } }
    return -1;
    }

  /// Add a AliAnalysisTaskEmcalJetCDF task - detailed signature
  /// \param ntracks name of tracks collection
  /// \param nclusters name of clusters collection
  /// \param ncells name of EMCAL cell collection
  /// \param tag tag name of analysis task
  /// \return AliAnalysisTaskEmcalJetCDF* task
  AliAnalysisTaskEmcalJetCDF* AddTaskEmcalJetCDF (
                                const char* ntracks    = "usedefault",
                                const char* nclusters  = "usedefault",
                                const char* ncells     = "usedefault",
                                const char* tag        = "CDF"
                              );

  /// Set parameters of a jet container
  /// \param jetCont AliJetContainer*
  /// \param jetptmin : min pt of jets in this container (default = 1.)
  /// \param jetptmax : max pt of jets in this container (default = 500.)
  /// \param jetareacutperc : cut jets under percentage of area given by algo radius (default = 0.)
  /// \param leadhadtype : 0 = charged, 1 = neutral, 2 = both (default = 2)
  /// \param nLeadJets : how many jets are to be considered the leading jet(s) (default = 1)
  /// \param mintrackpt : min track constituent pt to accept the jet (default = 0.15)
  /// \param maxtrackpt : max track constituent pt to accept the jet (default = 1000.)
  /// \return
  void jetContSetParams (
                          AliJetContainer* jetCont,
                          Float_t jetptmin = 1.,
                          Float_t jetptmax = 500.,
                          Float_t jetareacutperc = 0.,
                          Int_t leadhadtype = 2,
                          Int_t nLeadJets = 1,
                          Float_t mintrackpt = 0.15,
                          Float_t maxtrackpt = 1000.
                          );

} // end of NS_AliAnalysisTaskEmcalJetCDF

#endif // end of #ifndef ALIANALYSISTASKEMCALJETCDF_H

// kate: indent-mode none; indent-width 2; replace-tabs on;
