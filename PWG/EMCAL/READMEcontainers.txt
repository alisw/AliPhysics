/*! \page READMEcontainers EMCal Containers

EMCAL containers are used to handle arrays of objects (particles, clusters, jets) shared
among different tasks within the input ESD or AOD event via an easy interface hiding the
direct access from the user. Each task can handle several containers as needed. Once they
are defined, the AliAnalysisTaskEmcal will automatically connect them to the content they
are supposed to handle, so the content can be used from that moment on. Normally EMCAL 
containers are defined in the add macro of the task.

To better understand containers, consider the specific example of a track container:

# Create a container   {#createEMCalContainer}

Create a container in your ``AddTask`` macro using the name of the object of interest:
 
~~~{.cxx}
// The name of the tracks provided by the user (in this case, "tracks") must match the name in the input event!
AliTrackContainer * exampleTracks = new AliTrackContainer("tracks", "myTrackContainer");
// Make the tracks available in your task which inherits from AliAnalysisTaskEmcal
myTask->AdoptParticleContainer(exampleTracks);
~~~
  
As an alternative, the user can create the container and add it to the task in one step by using AliAnalysisTaskEmcal::AddTrackContainer(). Either approach yields equivalent results. 

Note the name given as the first argument to the ``AliTrackContainer`` is the name of the tracks branch in the event! The names are documented at the end of [this page](#emcalContainerBranchNames), but it is not recommended to use these names directly! Instead, it is very highly recommended to implement the "usedefault" pattern in your Add Task! For an example implementation of this pattern, see ``AddTaskEmcalJet`` (this is still a good example, even if you are not using jets!). Note that many of the shared Add Task macros use this pattern. In the rare case that the use needs to use the names directly, they are documented in the [section below](#emcalContainerBranchNames).

The user may create as many containers as desired, and each one will be **treated independently**. In general, the container interface is extremely flexible, allowing the user to use the contained objects however is desired. 

# Apply cuts

A user may also apply any desired kinematic cuts:
  
~~~{.cxx}
exampleTracks->SetMinPt(2);
exampleTracks->SetEtaLimits(-0.9, 0.9);
// ...
~~~

A wide variety of cuts are available, reducing code duplication and allowing the analyzer to focus on their analysis. One may also apply track selection. See \ref READMEtracks for more information.

Containers must be configured **independently** for each task. That is, for each created container (say, via the ``AddContainer()`` function), you must specify the cuts you would like to apply — the cuts on the previous containers are not preserved. 

# Iterate through the collection        {#emcalContainerIterateTechniques}

The containers are automatically loaded and ready to use in the ``Run()`` method in the user task. All you need to do is retrieve the tracks and iterate through the available tracks subject to the specified cuts:

~~~{.cxx}
// Retrieve the tracks
// Can also retrieve tracks by index
AliTrackContainer * tracks = dynamic_cast<AliTrackContainer *>(GetParticleContainer("myTrackContainer"));

// Iterable approach (using C++11)
AliTLorentzVector track;
for (auto trackIterator : tracks->accepted_momentum() )
{
  // trackIterator is a std::map of AliTLorentzVector and AliVTrack

  track.Clear();
  track = trackIterator.first;                   // Get the four-momentum
  AliVTrack * fullTrack = trackIterator.second;  // Get the full track

}
~~~

There are two recommended ways to iterate, depending on the need:

- Iterates through tracks that pass the specified kinematic cuts
  ~~~{.cxx}
  tracks->accepted_momentum() 
  ~~~

- Iterates through all tracks
  ~~~{.cxx}
  tracks->all_momentum() 
  ~~~

Both of these make available the correct four-momentum, and the full track object.  

Alternately (not recommended), one can get just the ``AliVTrack`` object with ``AliEmcalIterableContainer::accepted()`` and ``AliEmcalIterableContainer::all()``. If doing so with clusters, see the warnings [below](\ref emcalContainerClusterEnergyCorrections) about retrieving the proper cluster energy.

For more details on this issue, as well as more generally on iteration techniques, see ``AliEmcalIterableContainer``.

For more information on the containers, see the base class, ``AliEmcalContainer``, as well as the particular containers, ``AliClusterContainer``, ``AliParticleContainer``, ``AliTrackContainer``, and ``AliJetContainer``.

# Accessing corrected cluster energy            {#emcalContainerClusterEnergyCorrections}

The new framework uses a single cluster container for all levels of energy correction — the different levels of correction are stored as fields in each cluster. As there are now multiple possibilities for the cluster energy in one cluster object and in the container, you explicitly need to specify which energy you want to use in your analysis. 

From the cluster, energy can be retrieved as:

- Raw energy
  ~~~{.cxx}
  cluster->E()
  ~~~
  This energy usually:
  - implements cell-level energy/time calibrations and bad channel removal
  - does not include cluster-level non-linearity correction or any analysis-specific corrections (such as the "hadronic correction")

- Energy after non-linearity correction
  ~~~{.cxx}
  cluster->GetNonLinCorrEnergy()
  ~~~

- Energy after hadronic correction
  ~~~{.cxx}
  cluster->GetHadCorrEnergy()
  ~~~

One can set the default energy in the container, to be used for a given task by setting:
~~~{.cxx}
container->SetDefaultClusterEnergy(AliVCluster::VCluUserDefEnergy i)
~~~

where the available options for ``VCluUserDefEnergy`` are:
~~~{.cxx}
AliVCluster::kNonLinCorr
AliVCluster::kHadCorr
AliVCluster::kUserDefEnergy1
AliVCluster::kUserDefEnergy2
AliVCluster::kLastUserDefEnergy
~~~

This will then return the default energy in the ``AliTLorentzVector`` if the clusters are iterated over using ``clusters->accepted_momentum()`` or ``clusters->all_momentum()``, as above. Setting the default energy type will also return the default energy when using the ``GetMomentum(…)`` function of **cluster container**. In particular,
~~~{.cxx}
container->GetMomentum(AliTLorentzVector& mom, const AliVCluster* vcl)
~~~

Note that setting the default does not change the meaning of the cluster energy functions ``E()``, ``GetNonLinCorrEnergy()``, ``GetHadCorrEnergy()``.

The energy is also available via the **cluster** ``GetMomentum(…)`` functions, but it requires substantially more care to use it correctly. In general, it is recommended to use the options described above. However, if access directly from the cluster is required, the options are:
- Uses the cluster energy type specified when ``GetMomentum(…)`` is called:
  ~~~{.cxx}
  cluster->GetMomentum(AliTLorentzVector& mom, const Double_t* vertex, AliVCluster::VCluUserDefEnergy_t i)
  ~~~

- Will **always** give you the raw energy:
  ~~~{.cxx}
  cluster->GetMomentum(AliTLorentzVector& mom, Double_t* vertex)
  ~~~

# Branch names of cells, clusters, and tracks   {#emcalContainerBranchNames}

It is highly recommended not to use these names directly! Instead, use the "usedefault" pattern describe in the [section above](#createEMCalContainer).

For nearly all uses, they should use the following names:

| Object   | AOD            | ESD            |
| -------- | -------------- | -------------- |
| Cells    | "EMCALCells"   | "emcalCells"   |
| Clusters | "caloClusters" | "CaloClusters" |
| Tracks   | "tracks"       | "Tracks"       |

Default containers as listed above can be accessed directly passing the argument *usedefault* as name of the 
input containers. The example below demonstrates this for the cluster container:

~~~{.cxx}
AliClusterContainer *clustercont = new AliClusterContainer("usedefault", "myclustercont");
~~~

or

~~~{.cxx}
AliAnalysisTaskEmcal *emcaltask = new AliAnalysisTaskEmcal("emcaltask", kTRUE);
AliClusterContainer *clustercont = emcaltask->AddClusterContainer("usedefault");
~~~

*/ 
