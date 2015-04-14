/// \class AliAnalysisTaskEmcalJetCDF
/// \brief Analysis of leading jets distribution of pt and multiplicity
///
/// Analysis of leading jets distribution of pt and multiplicity, R distribution
/// N80 and Pt80 and Toward, Transverse, Away UE histos
///
/// \author Adrian SEVCENCO <Adrian.Sevcenco@cern.ch>, Institute of Space Science, Romania
/// \date Mar 23, 2015


#ifndef ALIANALYSISTASKEMCALJETCDF_H
#define ALIANALYSISTASKEMCALJETCDF_H

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

class AliAnalysisTaskEmcalJetCDF : public AliAnalysisTaskEmcalJet
  {
  public:

    AliAnalysisTaskEmcalJetCDF();
    AliAnalysisTaskEmcalJetCDF ( const char *name );
    virtual ~AliAnalysisTaskEmcalJetCDF();

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

    TH1D       *fH1;               //!<!  Pt distribution of all jets
    TH1D       *fH2;               //!<!  Eta distribution of all jets
    TH1D       *fH3;               //!<!  Phi distribution of all jets
    TH1D       *fH4;               //!<!  Multiplicity of all jets // 1 unit of multiplicity /bin
    TH1D       *fH5;               //!<!  Distribution of jets in events
    TH1D       *fH6;               //!<!  Jet1 Multiplicity Distribution
    TProfile   *fH7;               //!<!  N(jet1) vs P_{T}(jet1)

    TH1D       *fH8;               //!<!  Momentum distribution for leading jet (fragmentation function)
    TH1D       *fH8xi;             //!<!  Xi distribution for leading jet (fragmentation function)

    TH1D       *fH8_all;           //!<!  Momentum distribution for jets (fragmentation function)
    TH1D       *fH8xi_all;         //!<!  Xi distribution for jets (fragmentation function)

    TProfile   *fH9;               //!<!  N vs the Azimuthal Angle from Jet1
    TProfile   *fH10;              //!<!  P_{T} sum vs the Azimuthal Angle from Jet1

    TH1D       *fH9_bin;           //!<!  N vs the Azimuthal Angle from Jet1 - per bin sum
    TH1D       *fH10_bin;          //!<!  P_{T} sum vs the Azimuthal Angle from Jet1 - per bin sum

    TH1D       *fH9_bin_wojet1;    //!<!  N vs the Azimuthal Angle from Jet1 - per bin sum ; WITHOUT JET1 constituents
    TH1D       *fH10_bin_wojet1;   //!<!  P_{T} sum vs the Azimuthal Angle from Jet1 - per bin sum ; WITHOUT JET1 constituents

    TProfile  *fH15;               //!<!  <p_{T}> track vs the Distance R from Jet1
    TProfile  *fH15_n80;           //!<!  <p_{T}> track vs the Distance R from Jet1 - 80% of particles
    TProfile  *fH15_pt80;          //!<!  <p_{T}> track vs the Distance R from Jet1 - 80% of Pt

    TH1D      *fH15_bin;           //!<!  p_{T} track vs the Distance R from Jet1
    TH1D      *fH15_bin_n80;       //!<!  p_{T} track vs the Distance R from Jet1 - 80% of particles
    TH1D      *fH15_bin_pt80;      //!<!  p_{T} track vs the Distance R from Jet1 - 80% of Pt

    TH1D      *fH15_bin_all;           //!<!  p_{T} track vs the Distance R from owner jet
    TH1D      *fH15_bin_n80_all;       //!<!  p_{T} track vs the Distance R from owner jet - 80% of particles
    TH1D      *fH15_bin_pt80_all;      //!<!  p_{T} track vs the Distance R from owner jet - 80% of Pt

    TH1D      *fH20;               //!<!  Distribution of R in leading jet
    TH1D      *fH20_n80;           //!<!  Distribution of R in leading jet - 80% of particles
    TH1D      *fH20_pt80;          //!<!  Distribution of R in leading jet - 80% of Pt

    TH1D      *fH20_all;               //!<!  Distribution of R in jets
    TH1D      *fH20_n80_all;           //!<!  Distribution of R in jets - 80% of particles
    TH1D      *fH20_pt80_all;          //!<!  Distribution of R in jets - 80% of Pt

    TH1D      *fH23jet1;           //!<!  Jet1 Pt Distribution of particles

    TProfile  *fH24;               //!<!  Jet1 Size vs P_{T}(jet1) - 80% of particles
    TProfile  *fH25;               //!<!  Jet1 Size vs P_{T}(jet1) - 80% of Pt

    TProfile  *fH26;               //!<!  N_{chg} vs the Distance R from Jet1
    TProfile  *fH26_n80;           //!<!  N_{chg} vs the Distance R from Jet1 - 80% of particles
    TProfile  *fH26_pt80;          //!<!  N_{chg} vs the Distance R from Jet1 - 80% of Pt
    TProfile  *fH26jet1;           //!<!  N_{chg}(jet1) vs the Distance R from Jet1
    TProfile  *fH26jet1_n80;       //!<!  N_{chg}(jet1) vs the Distance R from Jet1 - 80% of particles
    TProfile  *fH26jet1_pt80;      //!<!  N_{chg}(jet1) vs the Distance R from Jet1 - 80% of Pt

    TH1D  *fH26_bin;               //!<!  N_{chg} vs the Distance R from Jet1
    TH1D  *fH26_n80_bin;           //!<!  N_{chg} vs the Distance R from Jet1 - 80% of particles
    TH1D  *fH26_pt80_bin;          //!<!  N_{chg} vs the Distance R from Jet1 - 80% of Pt
    TH1D  *fH26jet1_bin;           //!<!  N_{chg}(jet1) vs the Distance R from Jet1
    TH1D  *fH26jet1_n80_bin;       //!<!  N_{chg}(jet1) vs the Distance R from Jet1 - 80% of particles
    TH1D  *fH26jet1_pt80_bin;      //!<!  N_{chg}(jet1) vs the Distance R from Jet1 - 80% of Pt

    TProfile  *fH27;           //!<!  PT_{sum} vs the Distance R from Jet1
    TProfile  *fH27_n80;       //!<!  PT_{sum} vs the Distance R from Jet1 - 80% of particles
    TProfile  *fH27_pt80;      //!<!  PT_{sum} vs the Distance R from Jet1 - 80% of Pt
    TProfile  *fH27jet1;       //!<!  PT_{sum}(jet1) vs the Distance R from Jet1
    TProfile  *fH27jet1_n80;   //!<!  PT_{sum}(jet1) vs the Distance R from Jet1 - 80% of particles
    TProfile  *fH27jet1_pt80;  //!<!  PT_{sum}(jet1) vs the Distance R from Jet1 - 80% of Pt

    TH1D  *fH27_bin;           //!<!  PT_{sum} vs the Distance R from Jet1
    TH1D  *fH27_n80_bin;       //!<!  PT_{sum} vs the Distance R from Jet1 - 80% of particles
    TH1D  *fH27_pt80_bin;      //!<!  PT_{sum} vs the Distance R from Jet1 - 80% of Pt
    TH1D  *fH27jet1_bin;       //!<!  PT_{sum}(jet1) vs the Distance R from Jet1
    TH1D  *fH27jet1_n80_bin;   //!<!  PT_{sum}(jet1) vs the Distance R from Jet1 - 80% of particles
    TH1D  *fH27jet1_pt80_bin;  //!<!  PT_{sum}(jet1) vs the Distance R from Jet1 - 80% of Pt


    //_________________________________________________________________________
    AliJetContainer           *fJetsCont;                   //!<! Jets Container
    AliParticleContainer      *fTracksCont;                 //!<! Tracks Container
    AliClusterContainer       *fCaloClustersCont;           //!<! Clusters Container
    TClonesArray              *fTracksContArray;            //!<! the array of tracks from the tracks container

    AliEmcalJet               *fJet1;                       //!<! Leading Jet
    UInt_t                     fNJets_accepted;             ///<  Number of Jets found in event - accepted cuts applied by JetContainer
    UInt_t                     fNaccPart;                   ///<  Multiplicity in event - accepted tracks in tracks container
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

    AliAnalysisTaskEmcalJetCDF ( const AliAnalysisTaskEmcalJetCDF & );           // not implemented
    AliAnalysisTaskEmcalJetCDF &operator= ( const AliAnalysisTaskEmcalJetCDF & ); // not implemented

    /// \cond CLASSIMP
    ClassDef ( AliAnalysisTaskEmcalJetCDF, 2 );
    /// \endcond

  };
#endif

// kate: indent-mode none; indent-width 2; replace-tabs on;
