/*! \page READMEtracks Track selection
# Default track selection

To select the correct track filtering call
~~~{.cxx}
AliTrackContainer::SetDefTrackCutsPeriod(cRunPeriod);
~~~

where `cRunPeriod` is a string (e.g. `"lhc11h"`). This function has to be called only once. All the AliTrackContainer objects created after this, will have the correct settings to properly filter tracks in that particular dataset. In LEGO trains this line can be added in the global variables field (managed by the LEGO train operators).

In order to take advantage of all the features go the EMCal Jet Framework user tasks should derive from AliAnalysisTaskEmcal or AliAnalysisTaskEmcalJet. Tracks should be accessed using an AliTrackContainer object. This line should be used to set up a AliTrackContainer object when configuring your task (e.g. in the AddTask macro):

~~~{.cxx}
userTask->AddTrackContainer("tracks");
~~~

where `userTask` is the pointer to the user task object. This will work in the same way for both ESD and AOD (remember that `"Tracks"` should be used for ESD).

By default, hybrid track selection will be performed by AliTrackContainer. Other track selections are possibile. TPC-only tracks can be selected doing:

~~~{.cxx}
userTask->GetTrackContainer(0)->SetTrackFilterType(AliEmcalTrackSelection::kTPCOnlyTracks);
~~~

# Custom track cuts
Custom track cuts can also be used with the AliTrackContainer:

~~~{.cxx}
userTask->GetTrackContainer(0)->SetTrackFilterType(AliEmcalTrackSelection::kCustomTrackFilter);
~~~

To set custom track cuts we have to distinguish between ESD and AOD.

## ESD

For ESD analysis one can add AliESDtrackCuts objects to AliTrackContainer. This is done via:
~~~{.cxx}
userTask->GetTrackContainer(0)->AddTrackCuts(trackCuts);
~~~

where trackCuts is a AliESDtrackCuts object, e.g.:
~~~{.cxx}
AliESDtrackCuts* myCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kFALSE,1);
~~~

_Tip_: it is actually possibile to use any class derived from AliVCuts.

_Tip 2_: use AliTrackContainer::SetSelectionModeAny() and AliTrackContainer::SetSelectionModeAll() to switch between a track filtering that requires **all** cuts to accept the track (i.e. AND) or a track filtering that requires that **any** of the cuts accept the track (i.e. OR).

## AOD

For AOD analysis one can add a filter bit mask to AliTrackContainer. This is done via:
~~~{.cxx}
userTask->GetTrackContainer(0)->SetAODFilterBits(myFilterBits);
~~~

where `myFilterBits` is an integer containing a filter bit mask, e.g.:

~~~{.cxx}
UInt_t myFilterBits = 1<<8 | 1<<9;
~~~

The filter bit mask is always applied as an "OR" of the filter bits, namely track is accepted if **any** of the bits that are set in the test mask match with track filter bit mask.

_Tip_: you can also add AliVCuts (and AliESDtrackCuts) to the AliTrackContainer object even when analyzing AOD, the same way as it is done for ESD. Be aware that applying AliESDtrackCuts to AOD may have a different behavior compared to when they are applied to ESD.

# MC Particles

The AliMCParticleContainer class is able to directly filter primary particles in AOD events. Therefore in the case of ESD the AliEmcalMCTrackSelector task should be used to create a collection of AliAODMCParticle objects that can be dealt with the AliMCParticleContainer class (se below details about the AliEmcalMCTrackSelector task). This may change in the future (proper notification will be done in the relevant mailing list).

To access and filter primary particles in a MC event add an AliMCParticleContainer object to your task (assuming it derives from AliAnalysisTaskEmcal or AliAnalysisTaskEmcalJet):

~~~{.cxx}
AliMCParticleContainer *partCont = userTask->AddMCParticleContainer("mcparticles");
~~~

Sometimes one also needs to select only particles coming from a particle generator (this is the case of "cocktail" MC productions). This can also be done by selecting:

~~~{.cxx}
partCont->SetGeneratorIndex(index);
~~~

where `index` is a short integer corresponding to the requested generator. A frequent case is HIJING productions with injected rare probes. In this case usually `index=0` to select only HIJING particles.

# AliEmcalMCTrackSelector task (ESD only)
<strong>Class name</strong>: AliEmcalMCTrackSelector (PWG/EMCAL).

<strong>Add task macro</strong>: PWG/EMCAL/macros/AddTaskMCTrackSelector.C

~~~{.cxx}
AliEmcalMCTrackSelector* AddTaskMCTrackSelector(
  const char *outname    = "mcparticles",
  Bool_t      nk         = kFALSE,
  Bool_t      ch         = kFALSE,
  Double_t    etamax     = 1,
  Bool_t      physPrim   = kTRUE
)
~~~

<strong>Suggested parameter list</strong>:
~~~{.cxx}
"mcparticles", kFALSE, kFALSE, -1, kFALSE
~~~

When this task is executed at the beginning of an analysis train running on ESD, the MC particles can be filtered using the AliMCParticleContainer class as done for AOD.

## Using virtual track selection     {#EMCALVirtualTrackSelection}

In case a more strict track selection needs to be applied, users 
should use the virtual track selection. The virtual track selection provides
an interface handling the track selection in a seemless way both for ESDs and
AODs. The selection criteria have to be implemented by the user.

~~~{.cxx}
AliEmcalTrackSelectionESD *sel = new AliEmcalTrackSelectionESD;
sel->AddTrackCuts(AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE, 1));
~~~

Inside the task the virtual track selection can be used in two ways: Either one checks for
every track separate, via the function ``IsTrackAccepted()``, or one requests a TObjArray of all 
accepted tracks via ``GetAcceptedTracks()``. Everything needed at this step is provided by
AliEmcalTrackSelection, the virtual base class. 

*/
