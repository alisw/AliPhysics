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
#include <iostream>
#include <vector>
#include <utility>
#include <algorithm>

#include <Rtypes.h>
#include <TH1.h>
#include <TH1D.h>
#include <TProfile.h>
#include <TMath.h>

#include <AliEmcalJet.h>
#include <AliJetContainer.h>
#include <AliParticleContainer.h>
#include <AliClusterContainer.h>
#include <THistManager.h>

#include "AliAnalysisTaskEmcalJet.h"

using std::cout;
using std::endl;
using std::vector;
using std::pair;

class TH1;
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

    /// Setting of index of the jet container to be processed
    /// \param index of jet container
    void                        SetProcessJetCont (Int_t i) {idx_jetcont = i;}

    /// Get what jet container is set to be processed
    /// \return index of jer container
    Int_t                       GetProcessJetCont () {return idx_jetcont;}

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
    /// \param index of jet container
    /// \return true/false for succesful processing
    Bool_t                      ProcessJetContainer(Int_t idx_jet_container = 0);

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

    // Histogram list
    TH1D      *fH1;                //!<!  Pt distribution of all jets
    TH1D      *fH2;                //!<!  Eta distribution of all jets
    TH1D      *fH3;                //!<!  Phi distribution of all jets

//________________________________________________________________________
    TH1D      *fH4;                //!<!  Multiplicity of all jets - tracks // 1 unit of multiplicity /bin
    TH1D      *fH4c;               //!<!  Multiplicity of all jets - all constituents // 1 unit of multiplicity /bin

//________________________________________________________________________
    TH1D      *fH5;                //!<!  Distribution of jets in events

//________________________________________________________________________
    TH1D      *fH6;                //!<!  Jet1 Multiplicity Distribution - charged tracks
    TH1D      *fH6c;               //!<!  Jet1 Multiplicity Distribution - all constituents

//________________________________________________________________________
    TProfile  *fH7;                //!<!  N(jet1) vs P_{T}(jet1)
    TProfile  *fH7all;             //!<!  N(jet1) vs P_{T} - all jets

//________________________________________________________________________
    TH1D      *fH8;                //!<!  Momentum distribution for leading jet (fragmentation function)
    TH1D      *fH8xi;              //!<!  Xi distribution for leading jet (fragmentation function)

//________________________________________________________________________
    TH1D      *fH8_p;                //!<!  Momentum distribution for leading jet (fragmentation function) - Z_ptot
    TH1D      *fH8xi_p;              //!<!  Xi distribution for leading jet (fragmentation function) - Z_ptot

//________________________________________________________________________
    TH1D      *fH8_pt;                //!<!  Momentum distribution for leading jet (fragmentation function) - Z_pt
    TH1D      *fH8xi_pt;              //!<!  Xi distribution for leading jet (fragmentation function) - Z_pt

//________________________________________________________________________
    TH1D      *fH8_all;            //!<!  Momentum distribution for jets (fragmentation function)
    TH1D      *fH8xi_all;          //!<!  Xi distribution for jets (fragmentation function)

//________________________________________________________________________
    TH1D      *fH8_all_p;            //!<!  Momentum distribution for jets (fragmentation function) - Z_ptot
    TH1D      *fH8xi_all_p;          //!<!  Xi distribution for jets (fragmentation function) - Z_ptot

//________________________________________________________________________
    TH1D      *fH8_all_pt;            //!<!  Momentum distribution for jets (fragmentation function) - Z_pt
    TH1D      *fH8xi_all_pt;          //!<!  Xi distribution for jets (fragmentation function) - Z_pt

//________________________________________________________________________
    TProfile  *fH15;               //!<!  <p_{T}> track vs the Distance R from Jet1

    TProfile  *fH15_n90;           //!<!  <p_{T}> track vs the Distance R from Jet1 - 90% of particles
    TProfile  *fH15_n85;           //!<!  <p_{T}> track vs the Distance R from Jet1 - 85% of particles
    TProfile  *fH15_n80;           //!<!  <p_{T}> track vs the Distance R from Jet1 - 80% of particles
    TProfile  *fH15_n75;           //!<!  <p_{T}> track vs the Distance R from Jet1 - 75% of particles
    TProfile  *fH15_n70;           //!<!  <p_{T}> track vs the Distance R from Jet1 - 70% of particles

    TProfile  *fH15_pt90;          //!<!  <p_{T}> track vs the Distance R from Jet1 - 80% of Pt
    TProfile  *fH15_pt85;          //!<!  <p_{T}> track vs the Distance R from Jet1 - 85% of Pt
    TProfile  *fH15_pt80;          //!<!  <p_{T}> track vs the Distance R from Jet1 - 80% of Pt
    TProfile  *fH15_pt75;          //!<!  <p_{T}> track vs the Distance R from Jet1 - 75% of Pt
    TProfile  *fH15_pt70;          //!<!  <p_{T}> track vs the Distance R from Jet1 - 70% of Pt

//________________________________________________________________________
    TProfile  *fH15all;            //!<!  <p_{T}> track vs the Distance R from Jet1

    TProfile  *fH15all_n90;        //!<!  <p_{T}> track vs the Distance R from Jet1 - 90% of particles
    TProfile  *fH15all_n85;        //!<!  <p_{T}> track vs the Distance R from Jet1 - 85% of particles
    TProfile  *fH15all_n80;        //!<!  <p_{T}> track vs the Distance R from Jet1 - 80% of particles
    TProfile  *fH15all_n75;        //!<!  <p_{T}> track vs the Distance R from Jet1 - 75% of particles
    TProfile  *fH15all_n70;        //!<!  <p_{T}> track vs the Distance R from Jet1 - 70% of particles

    TProfile  *fH15all_pt90;       //!<!  <p_{T}> track vs the Distance R from Jet1 - 90% of Pt
    TProfile  *fH15all_pt85;       //!<!  <p_{T}> track vs the Distance R from Jet1 - 85% of Pt
    TProfile  *fH15all_pt80;       //!<!  <p_{T}> track vs the Distance R from Jet1 - 80% of Pt
    TProfile  *fH15all_pt75;       //!<!  <p_{T}> track vs the Distance R from Jet1 - 75% of Pt
    TProfile  *fH15all_pt70;       //!<!  <p_{T}> track vs the Distance R from Jet1 - 70% of Pt

//________________________________________________________________________
    TH1D      *fH15_bin;           //!<!  p_{T} SUM track vs the Distance R from Jet1

    TH1D      *fH15_bin_n90;       //!<!  p_{T} SUM track vs the Distance R from Jet1 - 90% of particles
    TH1D      *fH15_bin_n85;       //!<!  p_{T} SUM track vs the Distance R from Jet1 - 85% of particles
    TH1D      *fH15_bin_n80;       //!<!  p_{T} SUM track vs the Distance R from Jet1 - 80% of particles
    TH1D      *fH15_bin_n75;       //!<!  p_{T} SUM track vs the Distance R from Jet1 - 75% of particles
    TH1D      *fH15_bin_n70;       //!<!  p_{T} SUM track vs the Distance R from Jet1 - 70% of particles

    TH1D      *fH15_bin_pt90;      //!<!  p_{T} SUM track vs the Distance R from Jet1 - 80% of Pt
    TH1D      *fH15_bin_pt85;      //!<!  p_{T} SUM track vs the Distance R from Jet1 - 80% of Pt
    TH1D      *fH15_bin_pt80;      //!<!  p_{T} SUM track vs the Distance R from Jet1 - 80% of Pt
    TH1D      *fH15_bin_pt75;      //!<!  p_{T} SUM track vs the Distance R from Jet1 - 80% of Pt
    TH1D      *fH15_bin_pt70;      //!<!  p_{T} SUM track vs the Distance R from Jet1 - 80% of Pt

//________________________________________________________________________
    TH1D      *fH15all_bin;       //!<!  p_{T} SUM track vs the Distance R from owner jet

    TH1D      *fH15all_bin_n90;   //!<!  p_{T} SUM track vs the Distance R from owner jet - 80% of particles
    TH1D      *fH15all_bin_n85;   //!<!  p_{T} SUM track vs the Distance R from owner jet - 80% of particles
    TH1D      *fH15all_bin_n80;   //!<!  p_{T} SUM track vs the Distance R from owner jet - 80% of particles
    TH1D      *fH15all_bin_n75;   //!<!  p_{T} SUM track vs the Distance R from owner jet - 80% of particles
    TH1D      *fH15all_bin_n70;   //!<!  p_{T} SUM track vs the Distance R from owner jet - 80% of particles

    TH1D      *fH15all_bin_pt90;  //!<!  p_{T} SUM track vs the Distance R from owner jet - 80% of Pt
    TH1D      *fH15all_bin_pt85;  //!<!  p_{T} SUM track vs the Distance R from owner jet - 80% of Pt
    TH1D      *fH15all_bin_pt80;  //!<!  p_{T} SUM track vs the Distance R from owner jet - 80% of Pt
    TH1D      *fH15all_bin_pt75;  //!<!  p_{T} SUM track vs the Distance R from owner jet - 80% of Pt
    TH1D      *fH15all_bin_pt70;  //!<!  p_{T} SUM track vs the Distance R from owner jet - 80% of Pt

//________________________________________________________________________
    TH1D      *fH20;               //!<!  Distribution of R in leading jet

    TH1D      *fH20_n90;           //!<!  Distribution of R in leading jet - 90% of particles
    TH1D      *fH20_n85;           //!<!  Distribution of R in leading jet - 85% of particles
    TH1D      *fH20_n80;           //!<!  Distribution of R in leading jet - 80% of particles
    TH1D      *fH20_n75;           //!<!  Distribution of R in leading jet - 75% of particles
    TH1D      *fH20_n70;           //!<!  Distribution of R in leading jet - 70% of particles

    TH1D      *fH20_pt90;          //!<!  Distribution of R in leading jet - 90% of Pt
    TH1D      *fH20_pt85;          //!<!  Distribution of R in leading jet - 85% of Pt
    TH1D      *fH20_pt80;          //!<!  Distribution of R in leading jet - 80% of Pt
    TH1D      *fH20_pt75;          //!<!  Distribution of R in leading jet - 75% of Pt
    TH1D      *fH20_pt70;          //!<!  Distribution of R in leading jet - 70% of Pt


//________________________________________________________________________
    TH1D      *fH20all;           //!<!  Distribution of R in jets

    TH1D      *fH20all_n90;       //!<!  Distribution of R in jets - 80% of particles
    TH1D      *fH20all_n85;       //!<!  Distribution of R in jets - 80% of particles
    TH1D      *fH20all_n80;       //!<!  Distribution of R in jets - 80% of particles
    TH1D      *fH20all_n75;       //!<!  Distribution of R in jets - 80% of particles
    TH1D      *fH20all_n70;       //!<!  Distribution of R in jets - 80% of particles

    TH1D      *fH20all_pt90;      //!<!  Distribution of R in jets - 80% of Pt
    TH1D      *fH20all_pt85;      //!<!  Distribution of R in jets - 80% of Pt
    TH1D      *fH20all_pt80;      //!<!  Distribution of R in jets - 80% of Pt
    TH1D      *fH20all_pt75;      //!<!  Distribution of R in jets - 80% of Pt
    TH1D      *fH20all_pt70;      //!<!  Distribution of R in jets - 80% of Pt

//________________________________________________________________________
    TH1D      *fHg;                //!<!  Distribution of girth (radial girth) g = sum_jet_parts ( r_i * ( pt_i/pt_jet ) )

    TH1D      *fHg_n90;            //!<!  Distribution of girth (radial girth) g
    TH1D      *fHg_n85;            //!<!  Distribution of girth (radial girth) g
    TH1D      *fHg_n80;            //!<!  Distribution of girth (radial girth) g
    TH1D      *fHg_n75;            //!<!  Distribution of girth (radial girth) g
    TH1D      *fHg_n70;            //!<!  Distribution of girth (radial girth) g

    TH1D      *fHg_pt90;           //!<!  Distribution of girth (radial girth) g
    TH1D      *fHg_pt85;           //!<!  Distribution of girth (radial girth) g
    TH1D      *fHg_pt80;           //!<!  Distribution of girth (radial girth) g
    TH1D      *fHg_pt75;           //!<!  Distribution of girth (radial girth) g
    TH1D      *fHg_pt70;           //!<!  Distribution of girth (radial girth) g

//________________________________________________________________________
    TH1D      *fHptd;                //!<!  Distribution of dispersion d pt_D = sqrt ( sum (pt_i^2) )/sum (pt_i)

    TH1D      *fHptd_n90;            //!<!  Distribution of dispersion pt_D
    TH1D      *fHptd_n85;            //!<!  Distribution of dispersion pt_D
    TH1D      *fHptd_n80;            //!<!  Distribution of dispersion pt_D
    TH1D      *fHptd_n75;            //!<!  Distribution of dispersion pt_D
    TH1D      *fHptd_n70;            //!<!  Distribution of dispersion pt_D

    TH1D      *fHptd_pt90;           //!<!  Distribution of dispersion pt_D
    TH1D      *fHptd_pt85;           //!<!  Distribution of dispersion pt_D
    TH1D      *fHptd_pt80;           //!<!  Distribution of dispersion pt_D
    TH1D      *fHptd_pt75;           //!<!  Distribution of dispersion pt_D
    TH1D      *fHptd_pt70;           //!<!  Distribution of dispersion pt_D

//________________________________________________________________________
    TProfile  *fH23;               //!<!  Jet1 Size vs P_{T}(jet1)
    TProfile  *fH24;               //!<!  Jet1 Size vs P_{T}(jet1) - 80% of particles
    TProfile  *fH25;               //!<!  Jet1 Size vs P_{T}(jet1) - 80% of Pt

//________________________________________________________________________
    TProfile  *fH23all;            //!<!  Jet1 Size vs P_{T}(jet)
    TProfile  *fH24all;            //!<!  Jet1 Size vs P_{T}(jet) - 80% of particles
    TProfile  *fH25all;            //!<!  Jet1 Size vs P_{T}(jet) - 80% of Pt


    //_________________________________________________________________________
    AliJetContainer           *fJetsCont;                   //!<! Jets Container
    AliParticleContainer      *fTracksCont;                 //!<! Tracks Container
    AliClusterContainer       *fCaloClustersCont;           //!<! Clusters Container
    TClonesArray              *fTracksContArray;            //!<! the array of tracks from the tracks container
    TClonesArray              *fCaloClustContArray;         //!<! the array of clusters from the tracks container

    Int_t                      idx_jetcont;                 ///< index of jet container to be processed
    Int_t                      fNJets_accepted;             ///< Number of Jets found in event - accepted cuts applied by JetContainer
    Int_t                      fNaccPart;                   ///< Multiplicity in event - accepted tracks in tracks container
    Int_t                      fNaccClus;                   ///< Multiplicity in event - accepted clusters in cluster container
    THistManager               fHistManager;                ///< Histogram manager

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
    ClassDef ( AliAnalysisTaskEmcalJetCDF, 4 );
    /// \endcond

  };
#endif

// kate: indent-mode none; indent-width 2; replace-tabs on;
