/** 

@defgroup pwg_muondep_compact Quick Acc x Eff computation

Acc x Eff computation based on a single compacted simulation

@ingroup pwg_muondep

Using the AliMuonCompact* classes to compute a quick J/psi acceptance times efficiency

The starting point is a (pseudo-) ideal simulation of J/psi, producing (muon) ESDs (and Kinematics).

Ref "how to make pseudo-ideal simulations"

Then those ESDs are "compacted" into much smaller objects by the `AliMuonCompactTreeMaker` class.

~~~{.cxx}
.x runCompactTreeMaker("list.esd.txt","compacttreemaker.root");
~~~

The size reduction is very significant. For instance, for a 800k J/psi simulation (pp 13 TeV input distributions), the
ESDs alone are 6.9 Gbytes (the Kinematics and TrackRefs are 4.0 GB and 3.9 GB respectively). The corresponding output
`compacttreemaker.root` file is 69 Mbytes in comparison.  

A binary file must be created to hold the status of all the MCH manus for all the runs you intend to get an AccxEff for.
This implies a "scan" of the OCDB, and as such is better done using `cvmfs` OCDB or even a local copy of the (relevant
parts of) OCDB if you have one at hand.

```{.cxx}
AliMuonCompactManuStatus ms;
ms.WriteToBinaryFile("runlist.lhc15pp.txt","manustatus.lhc15pp.bis.dat","local:///alice/data/2015/OCDB");
```

Finally the compacted objects are used, together with the manu status file, to compute a quick AccxEff. The principle is
to loop over all tracks. For each track, we loop on its clusters and reject the ones which are "sitting on" a dead manu
(where dead can mean different things, e.g. dead because of HV, dead just because of occupancy, etc...). Once the
clusters have been rejected, we then check if the track still survive the tracking cuts. With the surviving tracks we
then do a simple invariant mass counting to get the number of J/psi corresponding to the given detector status.

```{.cxx}
int numberOfEventsToUse=100000; // use 0 for all
AliMuonCompactQuickAccEff q(numberOfEventsToUse);
q.ComputeEvolutionFromManuStatus("compacttreemaker.root","runlist.lhc15pp.txt","lhc15pp.allowing.monocathodes.root","manustatus.lhc15pp.dat","local:///alice/data/2015/OCDB",0);
```

The `AliMuonCompactQuickAccEffChecker` has been used to validate the method using a full simulation (aka regular one)
made by Hugo and Astrid for 2015 pp periods.

```
AliMuonCompactQuickAccEffChecker::CompareWithFullAccEff("EfficiencyJPsiRun_from_astrid.root","lhc15pp.allowing.monocathodes.root","lhc15pp.rejecting.monocathodes.root")
```
*/

