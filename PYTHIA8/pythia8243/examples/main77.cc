// main77.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Boris Blok,       blok@physics.technion.ac.il
//          Paolo Gunnellini, Paolo.Gunnellini@cern.ch

// The program calculates the 4 jet DPS cross section as a function of
// standard observables according to the model described in
// "Dynamical approach to MPI four-jet production in Pythia",
// B. Blok , P. Gunnellini
// Eur.Phys.J. C75 (2015) no.6, 282  [arXiv:1503.08246 [hep-ph]]

// The model is simulated by a dynamical reweighting of Pythia events.
// The notations and observables calculated are as in Fig. 4 of the article.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

//==========================================================================

// Method to calculate DPS weight according to the model above, 
// relative to the default one in Pythia.

double Reweighting(double x1, double x2, double x3, double x4,
  double scale1, double scale2){

  double newsigmaeff=2.*3.1415*(4.1+0.28*log(0.0012/x1)+4.1
    +0.28*log(0.0012/x2)+4.1+0.28*log(0.0012/x3)+4.1+0.28*log(0.0012/x4));

  newsigmaeff=newsigmaeff*0.389; //in mb
  //rescaling function

  double y=1.0;//0.5, 1 or 2

  double c1=0.18527-0.2058534*log(2*y)+0.0736003*log(2*y)*log(2*y);
  double c2=-0.016428+0.0254615*log(2*y)-0.0097517*log(2*y)*log(2*y);
  double c3=0.001444-0.00073161*log(2*y)+0.000314591*log(2*y)*log(2*y);

  double R=c1*log(scale2*scale2/y)+c2*log(scale2*scale2/y)*log(scale2*scale2/y)
    +c3*log(scale2*scale2/y)*log(scale2*scale2/y)*log(scale2*scale2/y);

  double ReweightFactor=R+(0.03-0.014427*log(y))*log((scale1*scale1)
    /(scale2*scale2));

  newsigmaeff=newsigmaeff/(1+ReweightFactor);

  // The final weight is the sigma effective implemented in the tune
  // reweighted with the value of the new sigma effective that you
  // obtained after the x-dependence and the Q^2 dependence.

  // Sigma effective value from our paper https://arxiv.org/pdf/1503.08246.pdf
  double FinalReweight=29.719/newsigmaeff;//only valid for 7 TeV

  return FinalReweight;

}

//==========================================================================

int main(){

  // Number of events, generated and listed ones (for jets).
  int nEvent    = 1000;

  // Select common parameters for SlowJet and FastJet analyses.
  int    power   = -1;     // -1 = anti-kT; 0 = C/A; 1 = kT.
  double R       = 0.5;    // Jet size.
  double pTMin   = 15.0;    // Min jet pT.
  double etaMax  = 5.0;    // Pseudorapidity range of detector.
  int    select  = 2;      // Which particles are included?
  int    massSet = 2;      // Which mass are they assumed to have?

  // Generator.
  Pythia pythia;

  // Process selection.
  pythia.readString("HardQCD:all = on");
  pythia.readString("PhaseSpace:pTHatMin = 15.");

  // pThat reweighting.
  pythia.readString("PhaseSpace:bias2Selection = on");
  pythia.readString("PhaseSpace:bias2SelectionPow = 4.");
  pythia.readString("PhaseSpace:bias2SelectionRef = 15.");

  // No event record printout.
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");

  // Tune from our paper https://arxiv.org/pdf/1503.08246.pdf
  pythia.readString("Tune:pp=5");
  pythia.readString("Tune:ee=3");
  pythia.readString("MultipartonInteractions:bProfile=1");
  pythia.readString("MultipartonInteractions:pT0Ref=2.659");
  pythia.readString("ColourReconnection:range=3.540");

  // LHC initialization.
  pythia.readString("Beams:eCM = 7000.");
  pythia.init();

  //Histograms
  Hist _hist_DeltaS_norm("DeltaS between the two jet systems", 14, 0., 3.1415);
  Hist _hist_DeltaPtRelSoft_norm("Normalized pT balance between the two jet"
    " systems", 10, 0., 1.);

  // Set up SlowJet jet finder in native mode.
  SlowJet slowJet( power, R, pTMin, etaMax, select, massSet, 0, false);
  //setup generation

  double sumWeights = 0.;

  // Generate events.
  for(int iEvent=0;iEvent<nEvent;++iEvent){
    if (!(pythia.next())){
        continue;//skip event in case of problem.
      }

    // Initial values for current event.
    unsigned int num_jets=0;
    unsigned int num_Softjets=0;
    double x1=pythia.info.x1();
    double x2=pythia.info.x2();
    double x3=-1.;
    double x4=-1.;
    bool hasMPI=false;

    // Identify incoming state of second MPI, if present.
    for (int i=7;i<pythia.event.size();++i){
      if(fabs(pythia.event[i].status())==31){
        x3=pythia.event[i].e()/pythia.event[1].e();
        x4=pythia.event[i+1].e()/pythia.event[2].e();
        hasMPI=true;
        break;
      }
    }

    // Find reweighting factor; check if the MPI scale is higher than 15 GeV.
    double scale1=pythia.info.pTMPI(0);
    double scale2=(hasMPI)? pythia.info.pTMPI(1):1.;
    double ReweightFactor= (scale2>15 && x3>0 && x4>0)
      ? Reweighting(x1,x2,x3,x4,scale1,scale2) : 1.;
    double weightFirst = pythia.info.weight();
    double weight=weightFirst*ReweightFactor;

    // Analyze jet content of event.
    slowJet.analyze( pythia.event );
    int nSlow = slowJet.sizeJet();
    vector <Vec4> myJets;
    for (int i = 0; i < nSlow; ++i) {
      if(slowJet.pT(i)>=50 && abs(slowJet.y(i))<4.7) num_jets+=1;
      if(slowJet.pT(i)>=20 && abs(slowJet.y(i))<4.7){
        num_Softjets+=1;
        myJets.push_back(slowJet.p(i));
      }
    }

    // Events with at least two hard and exactly four soft jets.
    if (num_jets >= 2 && num_Softjets == 4) {
      sumWeights+=weight;
      double SptSoft = (myJets.at(2)+myJets.at(3)).pT()
        / (myJets.at(2).pT()+myJets.at(3).pT());
      _hist_DeltaPtRelSoft_norm.fill(SptSoft, weight);
      double DeltaS = abs( (myJets.at(0)+myJets.at(1)).phi()
        - (myJets.at(2)+myJets.at(3)).phi() );
      if (DeltaS > M_PI) DeltaS = 2. * M_PI - DeltaS;
      _hist_DeltaS_norm.fill(DeltaS,weight);
    }
  }

  // Print run statistics and normalized histograms.
  pythia.stat();
  double binWidthDeltaPtRel  = 0.1;
  double binWidthDeltaS      = M_PI/14;
  _hist_DeltaPtRelSoft_norm /= (sumWeights*binWidthDeltaPtRel);
  _hist_DeltaS_norm         /= (sumWeights*binWidthDeltaS);
  cout << _hist_DeltaPtRelSoft_norm << _hist_DeltaS_norm;

  // Histogram in format to allow Python plotting.
  HistPlot hpl("main77plot");
  hpl.plotFrame("out77plot1", _hist_DeltaS_norm, "$\\Delta$ S (rad)",
                "1/$\\sigma$ d$\\sigma$/d$\\Delta$S",
                "1/$\\sigma$ d$\\sigma$/d$\\Delta$S",
                "h", "proton-proton collisions, $\\sqrt{s}$ = 7 TeV", true);
  hpl.plotFrame("out77plot2", _hist_DeltaPtRelSoft_norm, "$\\Delta$ p_T",
                "1/$\\sigma$ d$\\sigma$/d$\\Delta$p_T",
                "1/$\\sigma$ d$\\sigma$/d$\\Delta$p_T",
                "h", "proton-proton collisions, $\\sqrt{s}$ = 7 TeV", true);
  hpl.plot();

  // Done.
  return 0;

}
