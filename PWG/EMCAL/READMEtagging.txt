/*! \page READMEtagging Jet Tagging

AddTaskEmcalJetTagger.C contains two methods:
- one takes as input already existing jet arrays
- the other takes care of running the jet algorithm according to the settings.

The actual tagging is performed in the method AliAnalysisTaskEmcalJetTagger::MatchJetsGeo. The distance between pair of jets in the two inputs array is compared and those who are reciprocally closest are matched.

There are two tagging types:

~~~{.cxx}
enum JetTaggingType {
  kTag      = 0,
  kClosest  = 1
};
~~~

- `kTag` is used for "hard core" tagging and switches the TagStatus of the AliEmcalJet to 1 and set the TaggedJet to the companion jet.
- `kClosest` is used for the matching of the reconstructed and particle level jets in PYTHIA. It sets the ClosestJet (AliEmcalJet::fClosestJets[0]) and the distance.

## Hard core tagging

Used to reject combinatorial jets without applying an explicit high-pT cut on the constituents. Two jet finders are needed:

- normal jet finder used in the analysis (e.g. AKT, p<sub>T,const</sub> > 0.15 GeV/c, E<sub>min</sub> > 0.3 GeV) --> **jet**
-  jet finder with harder constituent cut (e.g. AKT,  p<sub>T,const</sub> > 4 GeV /c, E<sub>min</sub> > 0.3 GeV) --> **tag**

The **jets** are matched geometrically to the **tags** and flagged as "tagged".

JetFinderAKTCharged_R04_Escheme
~~~{.cxx}
gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C");
AliEmcalJetTask *jet = AddTaskEmcalJet(kUsedTracks, "", 1, 0.4, 1, kTrackPtCut, kClusPtCut, 0.005, 0, "Jet", 0., kFALSE, kFALSE, kFALSE);

//optional: background subtraction utils
AliEmcalJetUtilityGenSubtractor* genUtil = jet->AddUtility(new AliEmcalJetUtilityGenSubtractor("GenSubtractor"));
genUtil->SetUseExternalBkg(kTRUE);
genUtil->SetRhoName("RhoSparseR040");
genUtil->SetRhomName("RhoMassSparseR040");
genUtil->SetGenericSubtractionJetMass(kTRUE);
~~~

JetFinderAKTCharged_R04_Escheme_4GeV
~~~{.cxx}
AliEmcalJetTask *jettag = AddTaskEmcalJet(kUsedTracks, "", 1, 0.4, 1, 4., 4., 0.005, 0, "Jet", 0., kFALSE, kFALSE, kFALSE);
~~~

JetTaggerChargedR040_v2
~~~{.cxx}
gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJetTagger.C");
AliAnalysisTaskEmcalJetTagger *jetTagger = AddTaskEmcalJetTagger("Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme", "Jet_AKTChargedR040_PicoTracks_pT4000_E_scheme", 0.4, "", "", kUsedTracks, "", "TPC", "", kPhysSel);
jetTagger->SetNCentBins(1);
jetTagger->SetIsPythia(kIsPythia);
jetTagger->SetJetTaggingType(AliAnalysisTaskEmcalJetTagger::kTag);
jetTagger->SetJetTaggingMethod(AliAnalysisTaskEmcalJetTagger::kGeo);
AliJetContainer *contBase = jetTagger->GetJetContainer(0);
contBase->SetJetPtCut(0.15);
AliJetContainer *contTag = jetTagger->GetJetContainer(1);
contTag->SetJetPtCut(4.);
~~~

## Matching PYTHIA detector and PYTHIA particle level jets
It is done geometrically, as explained before. It is possible to request that the detector level jet has a given fraction of p<sub>T</sub> of the particle level jet. This doesn't affect the tagging but the result is shown in one of the output histograms (fh2PtJet2VsFraction) and the cut is applied to the histograms related to the "tagged" jets.
~~~{.cxx}
AliAnalysisTaskEmcalJetTagger *taskJetTagger = AddTaskEmcalJetTagger(nJetsThrm.Data(),nJetsPyt.Data(),0.4,rhoName.Data(),rhoMassName,"ThrmTracksEmb","","TPC",kCentEst,pSel,"","");
     taskJetTagger->SelectCollisionCandidates(AliVEvent::kAny);
     taskJetTagger->SetNCentBins(1);
     taskJetTagger->SetJetTaggingType(AliAnalysisTaskEmcalJetTagger::kClosest);
     taskJetTagger->SetMinFractionShared(fracShared);
~~~

### Minimum p<sub>T</sub> fraction shared

Note that the method:
~~~{.cxx}
AliJetContainer::GetFractionSharedPt(const AliEmcalJet *jet1, AliParticleContainer *cont2)
~~~
can be used without specifying `cont2` if the track containers of the two jets is the same.
If `cont2` is specified, a geometrical matching of the constituents of `jet1` and `jet2 = jet1->ClosestJet()` is performed. The tracks of `jet2` contained in `cont2` that have very similar (equal)  #eta, #phi, and p<sub>T</sub> as the tracks of `jet1` in `cont1` are accounted for the shared p<sub>T</sub> fraction.

*/
