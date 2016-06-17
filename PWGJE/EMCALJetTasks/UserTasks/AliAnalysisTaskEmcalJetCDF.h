#ifndef ALIANALYSISTASKEMCALJETCDF_H
#define ALIANALYSISTASKEMCALJETCDF_H
/// \file AliAnalysisTaskEmcalJetCDF.h
/// \brief Declaration of class AliAnalysisTaskEmcalJetCDF
///
/// \author Adrian SEVCENCO <Adrian.Sevcenco@cern.ch>, Institute of Space Science, Romania
/// \date Mar 23, 2015

/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskEmcalJet.h"
#include <THistManager.h>

/// \class AliAnalysisTaskEmcalJetCDF
/// \brief Analysis of jet shapes and FF of all jets and leading jets
///
/// Analysis of leading jets distribution of pt and multiplicity, R distribution
/// N80 and Pt80 and Toward, Transverse, Away UE histos
class AliAnalysisTaskEmcalJetCDF : public AliAnalysisTaskEmcalJet
  {
  public:

    AliAnalysisTaskEmcalJetCDF();
    AliAnalysisTaskEmcalJetCDF ( const char *name );
    virtual ~AliAnalysisTaskEmcalJetCDF();

    void                        UserCreateOutputObjects();
    void                        Terminate ( Option_t *option );

    THistManager                fHistManager                                      ;///< Histogram manager

    /// Get P() fraction of constituent to jet
    /// \param AliEmcalJet* jet
    /// \param AliVParticle* trk
    /// \return Z = trk->P()/jet->P()
    Double_t                    Z_ptot ( const AliEmcalJet* jet, const AliVParticle* trk )  const; // Get Z of constituent trk ; p total

    /// Get Pt() fraction of constituent to jet
    /// \param AliEmcalJet* jet
    /// \param AliVParticle* trk
    /// \return Z = trk->Pt()/jet->Pt()
    Double_t                    Z_pt   ( const AliEmcalJet* jet, const AliVParticle* trk )  const; // Get Z of constituent trk ; pt

    /// Get Xi for a double value z
    /// \return Xi
    Double_t                    Xi ( Double_t z )  const { return TMath::Log ( 1/z ); } // Get Xi of value z

  protected:
    void                        ExecOnce();
    Bool_t                      Run() ;

    /// Filling of histograms
    /// \return kTRUE if filling is succesful
    Bool_t                      FillHistograms()   ;

    /// Sorting of tracks in the event by pt (descending)
    /// \param AliVEvent*
    /// \return vector of indexes of constituents
    std::vector<Int_t>          SortTracksPt ( AliVEvent *event ) const;

    /// Sorting of tracks in the event by pt (descending) - using a particle container
    /// \param AliParticleContainer*
    /// \return vector of indexes of constituents
    std::vector<Int_t>          SortTracksPt ( AliParticleContainer *track_container ) const;

//     /// Return dR dinstance in eta,phi plane between 2 AliVParticle derived objects
//     /// \param AliVParticle* particle1
//     /// \param AliVParticle* particle2
//     /// \return distance
//     Double_t                    DeltaR ( const AliVParticle *part1, const AliVParticle *part2 );

//     /// Search for index(int) in array of ints
//     /// \param index - the int to be searched
//     /// \param array of ints
//     /// \return kTRUE if found
//     Bool_t                      IdxInArray ( Int_t index, TArrayI &array );

    /// Get pointer to a histogram from manager
    /// \param Name of histogram
    /// \return TObject*
    TObject* GetHistogram ( const char* histName );

  private:
    /// (pt,index) pair
    typedef std::pair<Double_t, Int_t> ptidx_pair;

    /// functional for sorting pair by first element - descending
    struct sort_descend
      {
      bool operator () ( const ptidx_pair &p1, const ptidx_pair &p2 )  { return p1.first > p2.first ; }
      };

    AliAnalysisTaskEmcalJetCDF ( const AliAnalysisTaskEmcalJetCDF& );           // not implemented
    AliAnalysisTaskEmcalJetCDF &operator= ( const AliAnalysisTaskEmcalJetCDF& ); // not implemented

    /// \cond CLASSIMP
    ClassDef ( AliAnalysisTaskEmcalJetCDF, 6 );
    /// \endcond

  };
#endif

// kate: indent-mode none; indent-width 2; replace-tabs on;
