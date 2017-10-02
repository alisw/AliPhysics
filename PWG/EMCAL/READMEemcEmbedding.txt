/*! \page READMEemcEmbedding EMCal Embedding Framework

# EMCal Embedding Framework

The EMCal Embedding Framework provides the means to embed an already reconstructed signal event into an
existing event. For example, the user may embed a Pythia signal event into a Pb-Pb event in an effort to
study the jet energy scale. To do so, a new class, AliAnalysisTaskEmcalEmbeddingHelper, manages access to
the embedded event, while users usually access the new input objects via AliEmcalContainer derived classes
(of course, direct access to the input objects is possible for those who are really interested - see
[below](\ref emcEmbeddingAdvancedTopics)). Note that the use of AliEmcalContainer objects does __not__
require use of AliAnalysisTaskEmcal, although using it does simplify the process.

Note that this framework was presented at an %Analysis Tutorial in Nov 2016. The presentation is a supplement
to the documentation presented here (it is highly recommended to read this documentation __first__) and
provides a practical example of embedding charged and full jets. For more information, see
[below](\ref emcEmbeddingAnalysisTutorial).

# Basics of embedding via the Embedding Helper

Access to the embedded event is managed via the embedding helper, AliAnalysisTaskEmcalEmbeddingHelper. The
embedding helper should be the first task after AliTaskCDBconnect! It is very straightforward to configure
the task:

~~~{.cxx}
// Create the Embedding Helper
AliAnalysisTaskEmcalEmbeddingHelper * embeddingHelper = AliAnalysisTaskEmcalEmbeddingHelper::AddTaskEmcalEmbeddingHelper();
// Set the file pattern. This example uses ptHardBin 4 of LHC12a15e_fix.
// The pT hard bin and anchor run can also be set by adding a printf() wild card to the string (ie %d)
// See the documentation of AliAnalysisTaskEmcalEmbeddingHelper::GetFilenames()
embeddingHelper->SetFilePattern("/alice/sim/2012/LHC12a15e_fix/169838/%d/AOD149");
// If the embedded file is an ESD, then set:
embeddingHelper->SetESD();
// Add additional configure as desired.
// For some examples...
// ... randomly select which file to start from:
embeddingHelper->SetRandomFileAccess(kTRUE);
// ... Start from a random event within each file
embeddingHelper->SetRandomEventNumberAccess(kTRUE);
// ... Set pt hard bin properties
embeddingHelper->SetPtHardBin(4);
embeddingHelper->SetNPtHardBins(11);
// etc..
// As your last step, always initialize the helper!
embeddingHelper->Initialize();
~~~

There are a large number of configuration options, including event selection and MC outlier rejection, which
are available as options for the embedding helper. See AliAnalysisTaskEmcalEmbeddingHelper. Regardless of the
options set, always be certain to call AliAnalysisTaskEmcalEmbeddingHelper::Initialize() when you are done
configuring the task!

Once the embedding helper has been created, it will manage access to the file. Then, the basic ideas is that
users access the embedded input objects via AliEmcalContainer derived classes (AliClusterContainer,
AliParticleContainer, AliTrackContainer, etc). For example, accessing tracks in an embedded AOD file via a
track container would look like:

~~~{.cxx}
// Create the new track container
AliTrackContainer * trackCont = new AliTrackContainer("tracks");
// Tell the track container that it should get the tracks from the embedded event
trackCont->SetIsEmbedding(kTRUE);
~~~

Although there are some details to understand, that's roughly all there is to it from the user perspective - it
should just work! (Recall again that AliEmcalContainers are __not__ required, although they substantially simplify
the user experience - see [below](\ref emcEmbeddingAdvancedTopics) for more information if you don't want to use them.)

With these basics in mind, there are a few main topics to discuss:
- [Embedding charged input objects (charged tracks)](\ref emcEmbeddingChargedInputObjects)
- [Embedding full (charged + neutral) input objects (cells, clusters, and tracks)](\ref emcEmbeddingFullInputObjects)
- [Recording Embedded Event properties](\ref emcEmbeddingQA)
- [Supporting multiple EMCal containers](\ref emcEmbeddingMultipleContainerSupport)

# Charged input objects (charged tracks)            {#emcEmbeddingChargedInputObjects}

To access charged tracks, it is as simple as described above: Simply create a track container and set that it
is embedding. To repeat from above:

~~~{.cxx}
// Create the new track container
AliTrackContainer * trackCont = new AliTrackContainer("tracks");
// Tell the track container that it should get the tracks from the embedded event
trackCont->SetIsEmbedding(kTRUE);
~~~

That is all there is to it - the embedded tracks are now available!

Although this process is quite straightforward, for this to be useful, we need to check on some additional
details. In particular, to actually embed the tracks in another event, that means that we will need to add at
least two track containers - one corresponding the input event, and one to the embedded event. Thus, it is
imperative that whatever task you are adding the containers to actually supports multiple containers explicitly.
Many base classes already do. For more information, see [below](\ref emcEmbeddingMultipleContainerSupport).

# Full (charged + neutral) input objects (cells, clusters, and tracks)      {#emcEmbeddingFullInputObjects}

In the case of embedding full objects, everything described above still applies! The Embedding Helper still
needs to be configured in the same way, and tracks can be accessed in exactly the same way! We are just
building on previous steps described here.

To use full objects, we need to also configure the cells and clusters, which involves configuring the
EMCal Corrections Framework. Fortunately, the correction framework will handle all of the steering, so
so we will just use advanced configuration techniques to setup the correction task (please review the
[correction task](\ref READMEemcCorrections) documentation if you are not familiar with these techniques).

Consider an example where we run the cell corrections separately for the input and embedded events, then
combine those cells to run the clusterizer, and then run the cluster level corrections. A sample
configuration (with an explanation below) is as follows:

~~~
inputObjects:
    cells:
        cells:
            branchName: "usedefault"
        cells_embed:
            branchName: "usedefault"
            embedding: true
        cells_combined:
            branchName: "emcalCellsCombined"
    clusterContainers:
        # Used for the clusterizer
        baseClusterContainer:
            branchName: "caloClustersCombined"
        # Used after clusterizer
        baseClusterContainer_1:
            minE: 0.0
            minPt: 0.0
        # Used for cluster-track matcher and after
        baseClusterContainer_2:
            branchName: "caloClustersCombined"
            minE: 0.0
            minPt: 0.0
            clusNonLinCorrEnergyCut: 0.15
    trackContainers:
        trackContainerName:
            # Sets the branch name
            branchName: "usedefault"
            minPt: 0.15
        trackContainerName_embed:
            embedding: true
CellEnergy_data:
    enabled: true
    cellsNames:
        - cells
CellBadChannel_data:
    enabled: true
    cellsNames:
        - cells
CellBadChannel_embed:
    enabled: true
    cellsNames:
        - cells_embed
CellTimeCalib_data:
    enabled: true
    cellsNames:
        - cells
CellCombineCollections_combined:
    enabled: true
    # Name of the cells branch in the embedded (external) event.
    externalCellsBranchName: "usedefault"
    # Name of the cells branch to be created for the combined cells in the input event.
    combinedCellsBranchName: "emcalCellsCombined"
    # Name of the cells to connect from the input event
    cellsNames:
        - cells
Clusterizer_combined:
    enabled: true
    cellsNames:
        - cells_combined
    # By selecting the cluster container here, we set where it will be output
    clusterContainersNames:
        - baseClusterContainer
ClusterExotics_combined:
    enabled: true
    cellsNames:
        - cells_combined
    clusterContainersNames:
        - baseClusterContainer_1
ClusterTrackMatcher_combined:
    enabled: true
    cellsNames:
        - cells_combined
    clusterContainersNames:
        - baseClusterContainer_2
    trackContainersNames:
        - trackContainerName
        - trackContainerName_embed

# Add additional correction tasks from here...
# ...
~~~

And then the add task configuration will look something like:

~~~{.cxx}
// To store the 3 corrections tasks
TObjArray correctionTasks;

// Create the 3 needed correction tasks
// These handle the cell corrections to the input data
// Selecting all corrections with "_data" in the name
correctionTasks.Add(AliEmcalCorrectionTask::AddTaskEmcalCorrectionTask("data"));
// These handle the cell corrections to the embedded data
// Selecting all corrections with "_embed" in the name
correctionTasks.Add(AliEmcalCorrectionTask::AddTaskEmcalCorrectionTask("embed"));
// These handle the cluster corrections to the combined (input + embedded) data
// Selecting all corrections with "_combined" in the name
// It is important that combined is last!
correctionTasks.Add(AliEmcalCorrectionTask::AddTaskEmcalCorrectionTask("combined"));

// Loop over all of the correction tasks to configure them the same
AliEmcalCorrectionTask * tempCorrectionTask = 0;
TIter next(&correctionTasks);
while (( tempCorrectionTask = static_cast<AliEmcalCorrectionTask *>(next())))
{
    tempCorrectionTask->SelectCollisionCandidates(kPhysSel);
    // Local configuration
    tempCorrectionTask->SetUserConfigurationFilename("userConfiguration.yaml");
    // Initialize the configuration
    tempCorrectionTask->Initialize();
}
~~~

## Cells

The cell corrections are split up into multiple sets of corrections because you often need to treat MC and
data cells differently. In particular, you'll often need to disable the energy and time calibration tasks
for the MC production. Once these individual cell corrections are completed, the cell collections will
need to be combined to be an input to the clusterizer.

To combine the cells, we utilize the task AliEmcalCorrectionCellCombineCollections which will handle the
combining the two cell collections. To ensure that each cell branch is corrected before they all of the
cells are combined, it is imperative that the correction tasks running the cells corrections on the input
and embedded events are run before the combined correction task! It is done properly in the example above.

## Clusters

After the cells have been combined, only one correction task needs to continue. This correction task is
labeled as "combined". Since all corrections are disabled by default, the other correction tasks will
automatically stop after running the cell corrections. All cluster corrections are then run as normal,
using the combined cluster (which are in cluster branch "caloClustersCombined" and accessible in the input
objects "baseClusterContainer", "baseClusterContainer_1", and "baseClusterContainer_2").

In the case that the user _does not_ want to run the clusterizer again (this should be fairly uncommon),
then multiple cluster containers can simply be added to each correction task, just as one would handle
tracks.

## Tracks (in cluster corrections)

Simply add each relevant track container to the relevant task. Make sure to mark the proper track container
as embedded for cluster corrections which require tracks. Everything will just work, as in the charged track
example above.

__Note on the Hadronic Correction__: Due to technical limitations when working with AOD files, the track
selection of the first particle container is applied to __all__ particles from __all__ particle containers.
This should not matter for most users, but it is important to be aware of it. A warning will be thrown when
running with such a configuration.

# Multiple Container Support                       {#emcEmbeddingMultipleContainerSupport}

To actually embed into another event, the user must add at least two EMCal containers to the task:
one corresponding the input event, and one to the embedded event. Thus, that task must support multiple containers.
There are two parts to this support:

 1. Adding multiple containers to the task.
 2. Properly utilizing multiple containers.

Fortunately, adding multiple containers often comes for free: any classes that inherit from
AliAnalysisTaskEmcal (and by extension AliAnalysisTaskEmcalJet) supports adding multiple containers
automatically.

Properly utilizing multiple containers can be a bit more complicated. There are two relevant groups here:

1. Core framework and embedding relevant classes.
2. The user's own analysis task.

In the case of the core framework and embedding relevant classes, most of them already support
multiple containers. Examples include
- AliEmcalJetTask
- All EMCal cluster and track corrections (except for the clusterizer) - see the
  [section on full input objects](\ref emcEmbeddingFullInputObjects) for more.
- AliAnalysisTaskEmcalQGTagging
- AliJetResponseMaker

In the case of the user's own task needing access to multiple containers, it is often as simple as adding
just a few lines of code to loop over container collections. For AliAnalysisTaskEmcal derived tasks, it should
look something like this:

~~~{.cxx}
// Assuming that we want to loop over clusters
AliClusterContainer * clusCont = 0;
TIter nextClusCont(&fClusterCollArray);
while ((clusCont = static_cast<AliClusterContainer*>(nextClusCont()))) {
    // A cluster container called cont is available.
    // It can be used as normal. That might look like:
    for (auto clusterIter = clusCont->all_momentum()) {
        // Get the cluster
        AliVCluster * cluster = clusterIter.second;
    }

    // Alternatively, if the index of the cluster is needed, use:
    auto clusItCont = clusCont->accepted_momentum();
    for (AliClusterIterableMomentumContainer::iterator clusIterator = clusItCont.begin(); clusIterator != clusItCont.end(); ++clusIterator) {
        // Get the cluster
        AliVCluster * cluster = clusterIter.second;
        // Get the cluster index
        int clusterIndex = clusIterator.current_index();
    }
}
~~~

Looping over multiple track containers is extremely similar. Once you've setup a loop for both, your class
is ready to go! It is up to you on how to handle it from here.

# Recording of Embedded Events Properties                                   {#emcEmbeddingQA}

The embedding helper will automatically store the properties of the embedded events such as cross section
and number of trials. These properties are recorded before event selection is applied to either the embedded
events or the data. However, it is often helpful to record those same properties after both the embedded and
data event selections. For such a case, there is an additional class, AliEmcalEmbeddingQA.

To use this case requires a few simple steps. First, an AliEmcalEmbeddingQA object should be added to your task.
Next, it must be initialized via AliEmcalEmbeddingQA::Initialize(). This should usually be done when other
histograms are being created and configured, such as in ``UserCreateOutputObjects()``. This histograms that are
created in this class should be added to your tasks output list via AliEmcalEmbeddingQA::AddQAPlotsToList().
Lastly, the properties should be recorded by calling AliEmcalEmbeddingQA::RecordEmbeddedEventProperties() after
the event selection has been performed in your task. The properties of the embedded events will then be available
in your tasks output.

Note that for compatibility reasons, AliAnalysisTaskEmcal::IsPythia() should be set to __false__ if you task
inherits from that class. Otherwise, histogram names will conflict.

# Note on jets and jet finding

When handling jet finding, a bit more care needs to be applied, especially if apply an artificial tracking
efficiency. This is due to the fact that the jet finder keeps track of constituents by index, which can cause
some ambiguity in which object applies to which container (if you are interested, see the
[advanced topics](\ref emcEmbeddingAdvancedTopics) for further details). Thus, the recommend approach is as
follows:

To apply an artificial efficiency to the embedded input objects, it is best to do so via the jet finder:

~~~{.cxx}
// Create the finder jet task as usual (called "jetTask")
// Set the track efficiency as desired
jetTask->SetTrackEfficiency(0.94);
// Tell it to apply to embedding only
jetTask->SetTrackEfficiencyOnlyForEmbedding(kTRUE);
~~~

Whether or not an artificial efficiency is applied, it is best to access jet constituents through new, dedicated
functions in AliEmcalJet which properly handle returning the appropriate object.

~~~{.cxx}
// Assuming that you have an AliEmcalJet called "jet"
// And you want to access the constituent tracks
for (int i = 0; i < jet->GetNumberOfTracks(); ++i) {
    // Access the track via the jet
    AliVTrack * track = jet->Track(i);
}
~~~

Previous functions such as AliEmcalJet::TrackAt(Int_t index, TCloensArray * arr) still work, but with multiple
containers, it is impossible to disambiguate which object comes from which container without external help
(this is handled by AliEmcalContainerIndexMap). AliEmcalJet::Track(Int_t index) handles such situations properly
and is provided as a convenience. And of course the use is still welcome to get the index directly via
AliEmcalJet::TrackAt(Int_t index) and then retrieve the object themselves (using AliEmcalContainerIndexMap will
be required if there are multiple TClonesArrays or containers available).

The above procedure also applies equally to clusters. The function AliEmcalContainer::Cluster(Int_t index)
replaces AliEmcalJet::ClusterAt(Int_t index, TClonesArray * arr) in the same way that
AliEmcalJet::Track(Int_t index) replaces AliEmcalJet::TrackAt(Int_t index, TClonesArray * arr).

Note that AliEmcalJet::Track() and AliEmcalJet::Cluster() will work regardless of which track or cluster
containers were added to the jet container that was used to access the jet!

# Analysis Tutorial Presentation               {#emcEmbeddingAnalysisTutorial}

The framework was presented at the %Analysis Tutorial during the
[November 2016 ALICE Week](https://indico.cern.ch/event/586577/). Note that a work in progress version
of the framework was presented. However, the slides were updated afterwards, and the interface to the
framework was already nearly complete, so the information should be still be useful. Also note the
presentation on the EMCal corrections framework. Since it is an integral part of steering the embedding
framework in the EMCal, it also worth watching that tutorial if you are unfamiliar with the framework
and are interested in embedding cells or clusters.

# Advanced topics                       {#emcEmbeddingAdvancedTopics}

## Direct access to the embedded event

If users want to access the events directly instead of using AliEmcalContainer derived classes - this is
__not recommended__ - then the embedded event can be retrieved via:

~~~{.cxx}
AliAnalysisTaskEmcalEmbeddingHelper::GetInstance()->GetExternalEvent()
~~~

From here, the user is then responsible for retrieving the information they are interested in.

## Framework details

These details are intended for experts - users can safely ignore them!

The framework was designed with the goal of making it as separate (ie encapsulated) as possible from all
other classes. By leveraging existing code as much as possible and keeping implementation details hidden
from the user, it ensures that the framework will maintain compatibility as the code base evolves, as well
as keeping each task focused on its particular task.

The implementation details for a variety of classes are as follows:

### AliAnalysisTaskEmcalEmbeddingHelper

This class is documented fairly extensively in the class, AliAnalysisTaskEmcalEmbeddingHelper.
Please see there for additional details!

### AliEmcalContainer

AliEmcalContainer::SetArray() checks whether the container is embedded and automatically selects the proper
event. Thus, this works seamlessly for the user. Plus, since it is in the base class, it will work for all
derived classes - nothing additional code is required.

### AliEmcalContainerIndexMap

AliEmcalContainerIndexMap is an important class that manages the mapping between all of the input object
arrays and the indices that are often stored in classes to identify other input objects. Effectively, it
can map from a local index in each array to a global index. For an instance in which this matters,
the cluster-track matcher stores the track index number in the cluster and the cluster index number in
the track under certain configurations. Without this mapping, it would be impossible to know which input
object array the index belongs. Note that although this map is used internally, it should be transparent
to most uses. If they only add one TClonesArray of a particular input object type, it will never be
apparent that it is used (ie indices will just look as normal).

As a concrete example, consider two input arrays, A and B. A has 200 entries and B has 150 entires. Since
the size of the arrays can change between events and we need a consistent mapping for A and B, a large
offset is required selected between the two arrays. For example, 100000. Thus, the third element in array
A is index 2, while the third element in array B is index 100002. Then, when 100002 is given to
AliEmcalContainerIndexMap, it will automatically return the third element in array B.

Note that this mapping can be between TCloensArrays and indices or AliEmcalContainer derived classes and
indices. AliClusterContainer and AliParticleContainer automatically store the mapping between accessed
TClonesArrays and indices. These mappings can be copied and automatically adapted to AliEmcalContainer
derived classes for use in another class by calling AliEmcalContainerMap::CopyMappingFrom().

Additional information is available in the class documentation (AliEmcalContainerIndexMap).

## Appendix

Development of this framework was documented in a [JIRA ticket](https://alice.its.cern.ch/jira/browse/ALPHY-53).

In case of problems in the code or documentation, please contact:
[Raymond Ehlers](mailto:raymond.ehlers@cern.ch)

*/
