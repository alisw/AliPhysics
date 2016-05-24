/*! \page READMEchangefw The “new EMCal framework”

# Introduction

This section lists the main differences between the old and the new framework.
While the old framework produced a new cluster container for each correction
that was applied to the clusters, the new framework works with only one cluster
container. In contrast to the old framework in the new framework specific energy fields of the cluster container are set to keep track of modifications to the cluster energy.


| After the Task        | Function to get the correct energy <br> (new framework) | Output ClusterName <br> (old framework)  |
| -------------               |-------------|                       -----|
| AddTaskClusterizerFast      | E(), raw energy       | CaloClusters     |
| AddTaskEmcalClusterMaker    | GetNonLinCorrEnergy() |  EmcCaloClusters |
| AddTaskHadCorr              | GetHadCorrEnergy()    | CaloClustersCorr |


# How to switch from old to new framework?

The following things need to be changed in order to work with the “new EMCal framework”.<br>
Adapt your run macro, the wagons in your train, your add task macro, and maybe your analysis task.
 
### Run Macro

The main changes in the run macro account for the input and output cluster names which are used in the different tasks.<br>
Old:<br>
Each task needs the specific name of the container that contains the applied energy corrections. See table above.<br>
New:<br>
Use "usedefault" as input name for tasks and leave output name (if present) empty "".<br>
For example see:<br>
NewRunMacro: `runEMCalJetSampleTask.C`

In your run macro you have to call now:<br>
`AliTrackContainer::SetDefTrackCutsPeriod(sRunPeriod);` <br>
This function has to be called only once. All the AliParticleContainer objects created after this, will have the correct settings to properly filter tracks in that particular dataset.

Don't forget to change the cluster names in your wagons in the train as well!

### AddTask macro

Since the new framework works with only one cluster container, you don't need to set a specific name in the input arguments of your Addtask macro. Your Addtask macro should now use the "usedefault" containers as input arguments. In your macro they then can be se to "tracks" and "caloClusters" (in most of the cases). See as an example `AddTaskEmcalJetSample.C` around line 44.

Second thing you should set in your AddTask macro is the default energy you want to use in the container (could also be set in the run macro).
This needs to be specified, otherwise the raw energy is used by default in the ClusterContainer.

Eg.<br>
~~~{.cxx}
Task->GetClusterContainer(clusName)->SetClusECut(0);  //by default set to 0.15 always
Task->GetClusterContainer(clusName)->SetClusPtCut(0); //by default set to 0
Task->GetClusterContainer(clusName)->SetClusUserDefEnergyCut(AliVCluster::kHadCorr,clusptcut);
Task->GetClusterContainer(clusName)->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
~~~

Available options are (`VCluUserDefEnergy`)
~~~{.cxx}
AliVCluster:: kNonLinCorr
AliVCluster:: kHadCorr
AliVCluster:: kUserDefEnergy1
AliVCluster:: kUserDefEnergy2
AliVCluster:: kLastUserDefEnergy
~~~

### Changes in the Analysis Task

As there are now multiple possibilities for the cluster energy in one cluster object and
in the container, you explicitly need to specify which energy you want to use in your analysis.
Eg.<br>
`cluster->E()` will always give you the raw energy<br>
`cluster->GetMomentum(Vector,Vertex)` will give you the raw energy.<br>
No matter what default energy you have set for your cluster container!

Better to use the two following functions:<br>
`cluster->GetMomentum(Vector,Vertex, VCluUserDefEnergy_t t)`  (see available options above for VCluUserDefEnergy).<br>
`clusterContainer->GetMomentum(Vector,Cluster)` (This option uses the default cluster energy)

New useful methods of the cluster container:<br>
~~~{.cxx}
SetDefaultClusterEnergy(AliVCluster::ToBeSpecified);
SetClusUserDefEnergyCut(AliVCluster::ToBeSpecified);
~~~

New useful methods of the clusters:
~~~{.cxx}
GetUserDefEnergy(AliVCluster::ToBeSpecified)
GetMomentum (AliTLorentzVector&, const Double_t*, AliVCluster::VCluUserDefEnergy_t);
GetHadCorrEnergy()
GetNonLinCorrEnergy() 
~~~

# Important: Things to keep in mind about the new framework

Now you carefully have to reevaluate which energy you use in your tasks. This
might not be the one you think it is. Old functions to retrieve the energy might not
give you back what you expect. See information above.

The `GetAccCluster()` function (calls `AliClusterContainer::ApplyClusterCuts`)
will loop over all energy cuts, no matter which default you have specified in your task.
So if you want to cut only on the hadcorr energy eg. make sure the other ones are set to zero.

Iterator: Strongly discouraged are: `GetNext{Particle,Cluster,Jet}`.
A better method to use is `clusters->GetCluster(index)` or `clusters->GetAcceptCluster(index)`.
Even better is the new C++11 syntax:<br>
<b>Needs to be updated here</b>

If you get the container by index, be careful as `AddTrackContainer()` and `AddParPcleContainer()` will both increase the same index.

Containers are configured independently for each task. You always need to specify cuts and default energy each time!

It is encouraged to use the `GetContainer(name/index)` function with name and not with index.
 
If you add a container with<br>
`task->AddParticleContainer(…)`
you can not retrieve it with<br>
`mytask->GetTrackContainer(0)`
If you add it with<br>
`task->AddTrackContainer(…)`
you can nevertheless retrieve it with<br>
`mytask->GetParticleContainer(0)`
 

### Additional Material:

Mini tutorials and slides<br>
https://indico.cern.ch/event/525454/contributions/2152155/attachments/1267654/1877240/2016-05-03-EMCalJetFrameworkTutorial.pdf<br>
https://indico.cern.ch/event/463551/contributions/1138526/attachments/1193228/1732592/2015-11-24-EMCalJetFrameworkUpdate.pdf<br>
E-mails:<br>
Salvatore 9th of December 2015<br>
Salvatore 1st February 2016<br>
Salvatore 17th February 2016<br>
Joel 26 April 2016<br>
Salvatore 5th May 2016<br>
*/

