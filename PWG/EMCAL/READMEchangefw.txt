/*! \page READMEchangefw Transitioning from the old to new EMCal framework

# Introduction

This section lists the main differences between the old and the new framework.
While the old framework produced a new cluster container for each correction
that was applied to the clusters, the new framework works with only one cluster
container. In contrast to the old framework in the new framework specific energy fields of the cluster container are set to keep track of modifications to the cluster energy.

# How to switch from old to new framework?

The following things need to be changed in order to work with the new EMCal framework:

- [Your run macro](\ref changeFrameworkRunMacro)
- [Your add task macro](\ref changeFrameworkAddTask)
- [Potentially your analysis task](\ref changeFrameworkAnalysisTask)

Generally, it is recommended to read the documentation on containers, available [here](\ref READMEcontainers).
 
# Run Macro            {#changeFrameworkRunMacro}

The main changes in the run macro account for the input and output cluster names which are used in the different tasks.

## Old framework

Each task needs the specific name of the container that contains the applied energy corrections, as described in the table below

| After the Task              | Function to get the correct energy <br> (new framework) | Output ClusterName <br> (old framework)  |
| --------------------------- | ------------------------------------------------------- | ---------------------------------------- |
| AddTaskClusterizerFast      | E(), raw energy                                         | CaloClusters                             |
| AddTaskEmcalClusterMaker    | GetNonLinCorrEnergy()                                   | EmcCaloClusters                          |
| AddTaskHadCorr              | GetHadCorrEnergy()                                      | CaloClustersCorr                         |

## New framework

Use containers to manage access to cells, clusters, and tracks. For more information on the containers, see the [containers](\ref READMEcontainers) page. Generally, use the "usedefault" pattern (described on the containers page) to manage the names of the various objects. Note that for shared tasks which have been updated to the new framework, the run macro should use "usedefault" as input name for tasks and leave output name (if present) empty "". 

In addition to generally using containers, the track container needs to be properly intialized. To properly intialiaze the track container, the macro must set:
~~~{.cxx}
AliTrackContainer::SetDefTrackCutsPeriod(sRunPeriod)
~~~
This function only has to be called once. All track containers (created via AliParticleContainer) created after this will have the correct settings to properly filter tracks in that particular dataset.

For an example of an updated run macro, see ``runEMCalJetSampleTask.C``.

In addition, Don't forget to change the cluster names in your wagons in the train as well!

# AddTask macro                    {#changeFrameworkAddTask}

Since the new framework corrects clusters in place (adding the correction into a new field in the cluster), you will not need to keep track of changing cluster names. Instead, there is a name for each object (clusters, tracks) that will be consistent through a run macro. Consequently, the best way to set the names is by using the "usedefault" pattern as described in the [containers page](\ref READMEcontainers).

Your AddTask macro should also set the default energy you want to use in the container (it could also be set in the run macro). This needs to be specified, otherwise the raw energy is used by default in the ClusterContainer. A possible configuration example is below

~~~{.cxx}
task->GetClusterContainer(clusName)->SetClusECut(0);  //by default set to 0.15
task->GetClusterContainer(clusName)->SetClusPtCut(0); //by default set to 0
task->GetClusterContainer(clusName)->SetClusUserDefEnergyCut(AliVCluster::kHadCorr, 0.15);
task->GetClusterContainer(clusName)->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
~~~

Possible default energy options are available [here](\ref emcalContainerClusterEnergyCorrections).

After making changes to the Add Task, remember to make the appropriate changes to your train wagons!

# Changes in the Analysis Task                    {#changeFrameworkAnalysisTask}

As there are now multiple possibilities for the cluster energy in one cluster object and
in the container, you explicitly need to specify which energy you want to use in your analysis.
Specifically, 
~~~{.cxx}
cluster->E() // Will always give you the raw energy!
cluster->GetMomentum(Vector,Vertex) // Will give you the raw energy!
~~~
**No matter what default energy you have set for your cluster container!**

For the various approaches to getting the desired cluster energy, see [here](\ref emcalContainerClusterEnergyCorrections).

For information on accessing the containers in your analysis task (including best practices), see [container access techniques](\ref emcalContainerIterateTechniques).

## Important things to keep in mind about the new framework

Now you carefully have to reevaluate which energy you use in your tasks. This
might not be the one you think it is. Old functions to retrieve the energy might not
give you back what you expect. See information above.

Some other notes:

- The `GetAccCluster()` function (which calls `AliClusterContainer::ApplyClusterCuts()`)
will loop over all possible energy cuts, no matter which default you have specified in your task.
So if you want to cut only on, for example, the hadcorr energy, make sure the other ones are set to zero.

- Containers are **configured independently for each task**. You always need to specify cuts and default energy each time!

- It is strongly encouraged to the container by name. If you get the container by index, be careful as ``AliAnalysisTaskEmcal::AddTrackContainer()`` and ``AliAnalysisTaskEmcal::AddParticleContainer()`` will both increase the same index.

## Additional Material:

Mini tutorials and slides:

- https://indico.cern.ch/event/525454/contributions/2152155/attachments/1267654/1877240/2016-05-03-EMCalJetFrameworkTutorial.pdf
- https://indico.cern.ch/event/463551/contributions/1138526/attachments/1193228/1732592/2015-11-24-EMCalJetFrameworkUpdate.pdf

E-mails:

- Salvatore 9th of December 2015<br>
- Salvatore 1st February 2016<br>
- Salvatore 17th February 2016<br>
- Joel 26 April 2016<br>
- Salvatore 5th May 2016<br>

*/

