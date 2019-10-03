# MC benchmark studies #
These are the analysis files for the MC benchmark studies done by Christian Bourjau in 2015 (Supervisor Michele Floris) for the HMTF. Often, observables are binned in narrow bins of "multiplicity". However, the choice of the multiplicity estimator introduces a great bias on the observable in question and thus needs to be understood. This is the goal of this service task. It produces a number of "standard plots" from MC truth data for various "multiplicity" estimators and observables. Also, this Task is agnostic of the event generator allowing for comparison between different models.

As a matter of fact, the concept of a multiplicity estimator was generalized in this code to an "event classifier". The code is organized around the concept of classifiers and observables and is modular in this regard. A classifier is an algorithm which calculates one "classifier value" per event. The task has many such classifiers. One example of such a classifier is the number of charged tracks N_ch within |eta| < 0.5 or the number of primary interactions nMPI.

The other logical unit of this service task are "observables". An example of an observable is the distribution of charged particles with respect to eta, binned in a given classifier. Thus, every observable has (at least) one classifier attached to it. The histogram associated with an observable is included within each observable class.

Furthermore, it is possible to define a "global trigger" for the task. The trigger is just another classifier, so it produces a classifier value for a given event. That trigger is currently implemented in a very simple way in `AliAnalysisTaskHMTFMCMultEst.cxx`. When adding the task, the user can call `SetGlobalTrigger` in order to activate a global trigger for this task. Possible choices of parameters are: "" (default), "InelGt0" and "SphericityGt09". Then, these triggers are initialized accordingly in  `UserCreateOutputObjects` without the need for altering the task's code. However, if one wants to add a new global trigger, one will have to define it in `UserCreateOutputObjects` (see the existing ones as a reference). 

## The execution logic of this task

In UserCreateOutputObjects():

 - Create the desired classifier objects
 - Create the desired observables.
 Since each observable needs at least one classifier, this is best done in a loop over the classifiers

In UserExec():

 - Reset the classifier values before doing anything!
 - `Fill()` each observable with the current event

## How to create a new classifier
Every classifier needs to inherit from `AliEventClassifierBase`. All the logic of the new classifier should be confined to its constructor and to the function `CalculateClassifierValue`. The result of that function should than be written to `fClassifierValue`.

## How to create a new observable
Similar to the classifiers, every observable needs to inherit from `AliObservableBase`. All the logic of the observable should be confined to its constructor (where the histogram is defined) and to the `Fill()` function. When requiring the classifier value, one should call the public function of the classifier `GetClassifierValue()`, which takes care of caching the value for successive calls.

## How to add these file to aliphysics

 - Add the cxx files to the `SRCS` variable in the file `CMakeList.txt` in this folder.
 - Add the files to the `PWGHMTFLinkDef.h`

