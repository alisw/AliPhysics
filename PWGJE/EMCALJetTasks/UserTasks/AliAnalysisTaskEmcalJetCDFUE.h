/// \class AliAnalysisTaskEmcalJetCDFUE
/// \brief Analysis of leading jets distribution of pt and multiplicity
///
/// Analysis of UE - Toward, Transverse, Away histos
///
/// \author Adrian SEVCENCO <Adrian.Sevcenco@cern.ch>, Institute of Space Science, Romania
/// \date Mar 23, 2015


#ifndef ALIANALYSISTASKEMCALJETCDFUE_H
#define ALIANALYSISTASKEMCALJETCDFUE_H

#include <cstdio>
#include <Rtypes.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TProfile.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TProfile.h>
#include <TMath.h>

#include "AliLog.h"
#include "AliEmcalJet.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"

#include "AliAnalysisTaskEmcalJet.h"

class TH1;
class TH2;
class TProfile;
class AliAnalysisUtils;
class AliJetContainer;
class AliParticleContainer;
class AliClusterContainer;

class AliAnalysisTaskEmcalJetCDFUE : public AliAnalysisTaskEmcalJet
  {
  public:

    AliAnalysisTaskEmcalJetCDFUE();
    AliAnalysisTaskEmcalJetCDFUE ( const char *name );
    virtual ~AliAnalysisTaskEmcalJetCDFUE();

    void                        UserCreateOutputObjects();
    void                        Terminate ( Option_t *option );

  protected:

    /// Filling of histograms
    /// \return kTRUE if filling is succesful
    Bool_t                      FillHistograms()   ;

    void                        ExecOnce();
    Bool_t                      Run() ;

    /// Sorting of tracks in the event by pt (descending)
    /// \param AliVEvent*
    /// \return vector of indexes of constituents
    std::vector<Int_t>          SortTracksPt ( AliVEvent *event ) const;

    /// Sorting of tracks in the event by pt (descending) - using a particle container
    /// \param AliParticleContainer*
    /// \return vector of indexes of constituents
    std::vector<Int_t>          SortTracksPt ( AliParticleContainer *track_container ) const;

    /// Return dR dinstance in eta,phi plane between 2 AliVParticle derived objects
    /// \param AliVParticle* particle1
    /// \param AliVParticle* particle2
    /// \return distance
    Double_t                    DeltaR ( const AliVParticle *part1, const AliVParticle *part2 );

    /// Search for index(int) in array of ints
    /// \param index - the int to be searched
    /// \param array of ints
    /// \return kTRUE if found
    Bool_t                      IdxInArray ( Int_t index, TArrayI &array );

    TProfile  *fH21;               //!<!  N(in the event - including jet1) vs P_{T}(jet1) - full phi
    TProfile  *fH21Toward;         //!<!  N(in the event - including jet1) vs P_{T}(jet1) - toward
    TProfile  *fH21Transverse_min; //!<!  N(in the event - including jet1) vs P_{T}(jet1) - trans min
    TProfile  *fH21Transverse_max; //!<!  N(in the event - including jet1) vs P_{T}(jet1) - trans max
    TProfile  *fH21Away;           //!<!  N(in the event - including jet1) vs P_{T}(jet1) - away

    TH1D      *fH21_bin;               //!<!  N(in the event - including jet1) vs P_{T}(jet1) - per bin sum - full phi
    TH1D      *fH21Toward_bin;         //!<!  N(in the event - including jet1) vs P_{T}(jet1) - per bin sum - toward
    TH1D      *fH21Transverse_min_bin; //!<!  N(in the event - including jet1) vs P_{T}(jet1) - per bin sum - trans min
    TH1D      *fH21Transverse_max_bin; //!<!  N(in the event - including jet1) vs P_{T}(jet1) - per bin sum - trans max
    TH1D      *fH21Away_bin;           //!<!  N(in the event - including jet1) vs P_{T}(jet1) - per bin sum - away

    TH1D      *fH21_bin_wojet1;               //!<!  N(in the event - without jet1) vs P_{T}(jet1) - per bin sum - full phi
    TH1D      *fH21Toward_bin_wojet1;         //!<!  N(in the event - without jet1) vs P_{T}(jet1) - per bin sum - toward
    TH1D      *fH21Transverse_min_bin_wojet1; //!<!  N(in the event - without jet1) vs P_{T}(jet1) - per bin sum - trans min
    TH1D      *fH21Transverse_max_bin_wojet1; //!<!  N(in the event - without jet1) vs P_{T}(jet1) - per bin sum - trans max
    TH1D      *fH21Away_bin_wojet1;           //!<!  N(in the event - without jet1) vs P_{T}(jet1) - per bin sum - away

    TProfile  *fH22;               //!<!  PT_{sum}(in the event - including jet1) vs P_{T}(jet1) - full phi
    TProfile  *fH22Toward;         //!<!  PT_{sum}(in the event - including jet1) vs P_{T}(jet1) - toward
    TProfile  *fH22Transverse_min; //!<!  PT_{sum}(in the event - including jet1) vs P_{T}(jet1) - trans min
    TProfile  *fH22Transverse_max; //!<!  PT_{sum}(in the event - including jet1) vs P_{T}(jet1) - trans max
    TProfile  *fH22Away;           //!<!  PT_{sum}(in the event - including jet1) vs P_{T}(jet1) - away

    TH1D  *fH22_bin;               //!<!  PT_{sum}(in the event - including jet1) vs P_{T}(jet1) - per bin sum - full phi
    TH1D  *fH22Toward_bin;         //!<!  PT_{sum}(in the event - including jet1) vs P_{T}(jet1) - per bin sum - toward
    TH1D  *fH22Transverse_min_bin; //!<!  PT_{sum}(in the event - including jet1) vs P_{T}(jet1) - per bin sum - trans min
    TH1D  *fH22Transverse_max_bin; //!<!  PT_{sum}(in the event - including jet1) vs P_{T}(jet1) - per bin sum - trans max
    TH1D  *fH22Away_bin;           //!<!  PT_{sum}(in the event - including jet1) vs P_{T}(jet1) - away

    TH1D  *fH22_bin_wojet1;               //!<!  PT_{sum}(in the event - without jet1) vs P_{T}(jet1) - per bin sum - full phi
    TH1D  *fH22Toward_bin_wojet1;         //!<!  PT_{sum}(in the event - without jet1) vs P_{T}(jet1) - per bin sum - toward
    TH1D  *fH22Transverse_min_bin_wojet1; //!<!  PT_{sum}(in the event - without jet1) vs P_{T}(jet1) - per bin sum - trans min
    TH1D  *fH22Transverse_max_bin_wojet1; //!<!  PT_{sum}(in the event - without jet1) vs P_{T}(jet1) - per bin sum - trans max
    TH1D  *fH22Away_bin_wojet1;           //!<!  PT_{sum}(in the event - without jet1) vs P_{T}(jet1) - per bin sum - away

    TH1D      *fH23;               //!<!  Event Pt Distribution of particles
    TH1D      *fH23Toward;         //!<!  'Toward' Pt Distribution of particles
    TH1D      *fH23Transverse_min; //!<!  'Transverse' MIN Pt Distribution of particles
    TH1D      *fH23Transverse_max; //!<!  'Transverse' MAX Pt Distribution of particles
    TH1D      *fH23Away;           //!<!  'Away' Pt Distribution of particles

    TProfile  *fH40;           //!<!  total particles fNPart w.r.t PTmax (pt of leading particle from leading jet) - full phi
    TProfile  *fH40toward;     //!<!  total particles fNPart w.r.t PTmax (pt of leading particle from leading jet) - toward
    TProfile  *fH40away;       //!<!  total particles fNPart w.r.t PTmax (pt of leading particle from leading jet) - away
    TProfile  *fH40transmin;   //!<!  total particles fNPart w.r.t PTmax (pt of leading particle from leading jet) - trans min
    TProfile  *fH40transmax;   //!<!  total particles fNPart w.r.t PTmax (pt of leading particle from leading jet) - trans max

    TH1D  *fH40_bin;           //!<!  total particles fNPart w.r.t PTmax (pt of leading particle from leading jet) - per bin sum - full phi
    TH1D  *fH40toward_bin;     //!<!  total particles fNPart w.r.t PTmax (pt of leading particle from leading jet) - per bin sum - toward
    TH1D  *fH40away_bin;       //!<!  total particles fNPart w.r.t PTmax (pt of leading particle from leading jet) - per bin sum - away
    TH1D  *fH40transmin_bin;   //!<!  total particles fNPart w.r.t PTmax (pt of leading particle from leading jet) - per bin sum - trans min
    TH1D  *fH40transmax_bin;   //!<!  total particles fNPart w.r.t PTmax (pt of leading particle from leading jet) - per bin sum - trans max

    TProfile  *fH41;           //!<!  PTsum w.r.t PTmax - full phi
    TProfile  *fH41toward;     //!<!  PTsum w.r.t PTmax - toward
    TProfile  *fH41away;       //!<!  PTsum w.r.t PTmax - away
    TProfile  *fH41transmin;   //!<!  PTsum w.r.t PTmax - trans min
    TProfile  *fH41transmax;   //!<!  PTsum w.r.t PTmax - trans max

    TH1D  *fH41_bin;           //!<!  PTsum w.r.t PTmax - per bin sum - full phi
    TH1D  *fH41toward_bin;     //!<!  PTsum w.r.t PTmax - per bin sum - toward
    TH1D  *fH41away_bin;       //!<!  PTsum w.r.t PTmax - per bin sum - away
    TH1D  *fH41transmin_bin;   //!<!  PTsum w.r.t PTmax - per bin sum - trans min
    TH1D  *fH41transmax_bin;   //!<!  PTsum w.r.t PTmax - per bin sum - trans max

    AliJetContainer           *fJetsCont;                   //!<! Jets Container
    AliParticleContainer      *fTracksCont;                 //!<! Tracks Container
    AliClusterContainer       *fCaloClustersCont;           //!<! Clusters Container
    TClonesArray              *fTracksContArray;            //!<! the array of tracks from the tracks container
    TClonesArray              *fCaloClustContArray;         //!<! the array of clusters from the tracks container

    AliEmcalJet               *fJet1;                       //!<! Leading Jet
    UInt_t                     fNJets_accepted;             ///<  Number of Jets found in event - accepted cuts applied by JetContainer
    UInt_t                     fNaccPart;                   ///<  Multiplicity in event - accepted tracks in tracks container
    Int_t                      fNaccClus;                   //!<! Multiplicity in event - accepted clusters in cluster container
    Double_t                   fEvPt;                       ///<  Scalar sum of pt off all accepted tracks in events

    std::vector<Int_t>         fJet1_sorted_idxvec;         ///< vector of sorted indexes of particles in leading jet
    std::vector<Int_t>         fEvent_sorted_idxvec;        ///< vector of sorted indexes of accepted tracks in the event

  private:
    /// (pt,index) pair
    typedef std::pair<Double_t, Int_t> ptidx_pair;

    /// functional for sorting pair by first element - descending
    struct sort_descend
      {
      bool operator () ( const ptidx_pair &p1, const ptidx_pair &p2 )  { return p1.first > p2.first ; }
      };

    AliAnalysisTaskEmcalJetCDFUE ( const AliAnalysisTaskEmcalJetCDFUE & );           // not implemented
    AliAnalysisTaskEmcalJetCDFUE &operator= ( const AliAnalysisTaskEmcalJetCDFUE & ); // not implemented

    /// \cond CLASSIMP
    ClassDef ( AliAnalysisTaskEmcalJetCDFUE, 0 );
    /// \endcond

  };
#endif

// kate: indent-mode none; indent-width 2; replace-tabs on;
