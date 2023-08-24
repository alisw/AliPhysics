/*! \page READMEtrgfw The EMCAL trigger framework

# Reconstruction of EMCAL trigger patches from FastOR / FEE data

The EMCAL Level1 trigger decision is based on reconstructed trigger patches. Events are selected once a 
trigger patch of a given type is found above threshold. Trigger patches are reconstructed using a sliding 
window.

Currently two Level1 triggers are implemented: A gamma trigger and a jet trigger. The setup of the gamma 
trigger is unique for all data sets while the setup of the jet trigger might differ for various data sets.

| Gamma trigger            | Jet trigger                    |
|--------------------------|--------------------------------|
| Subregion size: 1        | Subregion size: 4              |
| Number of subregions: 2  | Number of subregions: 2 or 4   |
| Total patch size: 2x2    | Total patch size: 8x8 or 16x16 |

Trigger patches are reconstructed by the trigger maker wagon consisting of the task AliEmcalTriggerMakerTask.
The trigger maker task configures and runs the trigger maker kernel for every event. The trigger maker 
kernel contains the trigger algorithms for gamma and jet trigger. Trigger algorithms are run for EMCAL
and DCAL separately.

## Configuration of the trigger maker 

The trigger maker can be added to the train by using its add macro AddEmcalTriggerMakerNew.C

~~~{.cxx}
gROOT->LoadMacro(Form("%s/PWG/EMCAL/macros/AddTaskTriggerMakerNew.C", gSystem->ExpandPathName(gSystem->Getenv("ALICE_PHYSICS"))));
AliEmcalTriggerMakerTask *triggermaker = AddTaskTriggerMakerNew();
triggermaker->SelectCollisionCandidates(AliVEvent::kINT7 | AliVEvent::kEMC7 | AliVEvent::kEMCEJE | AliVEvent::kEMCEGA);
~~~

For existing datasets, identified by the run numbers in the dataset, the trigger maker can configure the trigger algorithms 
itself and no interaction has to be done by the user. Thus it is possible to provide a configuration by the user. This is 
particularly of relevance if a desired configuration for a future dataset needs to be tested before on existing data. The
following configurations are implemented so far:

| Configuration       | Data sets:                | Detectors    | Remarks                                                           |
|---------------------|---------------------------|--------------|-------------------------------------------------------------------|
|ConfigureForPP2011   | pp, 2011                  | EMCAL        | Only level0 algorithm                                             |       
|ConfigureForPbPb2011 | Pb-Pb, 2011               | EMCAL        | L1 for gamma and jet trigger (16x16), only high-threshold online bits     |
|ConfigureForPP2011   | pp, 2012                  | EMCAL        | L1 for gamma and jet trigger (16x16), only high-threshold online bits     |
|ConfigureForPPb2013  | p-Pb, 2013                | EMCAL        | L1 for gamma and jet trigger (16x16), high- and low-threshold online bits |
|ConfigureForPP2015   | pp, 2015-2016, p-Pb 2016  | EMCAL + DCAL | L1 for gamma and jet trigger (16x16), high- and low-threshold online bits |
|ConfigureForPbPb2015 | Pb-Pb 2015                | EMCAL + DCAL | L1 for gamma and jet trigger (8x8), only high-threshold online bits, background patches |

More configurations need to implemented before running by the experts.

The following example configures the trigger maker for Pb-Pb data from 2015:

~~~{.cxx}
AliEmcalTriggerMakerTask *triggermaker = AddTaskTriggerMakerNew();
triggermaker->GetTriggerMaker()->ConfigureForPbPb2015();
~~~

## Masking of bad cells or fastors

The trigger maker reads in cell and fastor data from the input event. Handling of bad cells is typically done by the 
EMCAL tender or the [EMCAL correction task](\ref READMEemcCorrections). Thus the trigger maker provies the possibility
to mask additional bad cells:

- AddOfflineBadChannel (Short_t absId)
- ReadOfflineBadChannelFromStream (std::istream &stream)
- ReadOfflineBadChannelFromFile (const char *fname)

These are functions of the [trigger maker kernel](\ref AliEmcalTriggerMakerKernel) Cells with the given IDs are ignored when 
reading the cell data. The following example configures the trigger maker masking cel 12303:

~~~{.cxx}
AliEmcalTriggerMakerTask *triggermaker = AddTaskTriggerMakerNew();
triggermaker->GetTriggerMaker()->AddOfflineBadChannel(12303);
~~~

FastORs can only be handled by the trigger maker. Also for fastors it is possible to define certain fastors directly

- AddFastORBadChannel (Short_t absId)
- ReadFastORBadChannelFromStream (std::istream &stream)
- ReadFastORBadChannelFromFile (const char *fname)

For masked fastor channels the L1 ADC (L1 time sums) are ignored when reading the FastOR data, trigger bits however are still
handled.

Apart from applying the masking manually, it is also possible to read the masking from OCDB or an OADB container. 
The masking in the OCDB corresponds to the masking which was applied at hardware level. It is of relevance in simulation 
for fee / offlinepatches in order to have the same acceptance map as it was applied online. The trigger maker will
also take care about mapping issues for run1 and run2. In order to apply the masking from the OCDB, users have to do the 
following:

~~~{.cxx}
AliEmcalTriggerMakerTask *triggermaker = AddTaskTriggerMakerNew();
triggermaker->SetLoadFastORMaskingFromOCDB();
~~~

__Note__: The CDBconnect task is required when applying the FastOR masking from the OCDB.
 
A second option is to provide an OADB container with a masking. This option is focused on maskin fastors which
turned out to be bad during offline QA and were not handled online. In this case the triggered sample is cleaned
up from noisy triggers. The OADB container has to follow the definitions

- The name of the OADB container has to be *AliEmcalMaskedFastors*
- The content is a TObjArray with FastOR Absolute IDs, stored as TParameter<int>

The following code example masks FastORs 1211, 1560, 2300 for the run range 195000 to 200000:

~~~{.cxx}
AliOADBContainer fastorContainer("AliEmcalMaskedFastors");
TObjArray *fastorData = new TObjArray;
std::array<int, 3> masked = {1211, 1560, 2300};
for(auto m : masked) fastorData->Add(new TParameter<int>(Form("%d", m), m));
fastorContainer.AppendObject(fastorData, 195000, 200000);
fastorContainer.WriteToFile("EmcalMaskedFastors.root");
~~~

The new container needs to be specified when adding the trigger maker:

~~~{.cxx}
AliEmcalTriggerMakerTask *triggermaker = AddTaskTriggerMakerNew();
triggermaker->SetMaskedFastorOADBContainer("EmcalMaskedFastors");
~~~

In case users want to apply the FastOR masking also to cell data, this can be 
specified via SetApplyTRUMaskingToFEE();

__Important notice for simulations__: In order to have a realistic simulation of the acceptance, the only bad 
channel map on fastor level should be used, not the one on cell level. As tender or correction task will handle bad
cells, the offline bad cell map will be applied. Therefor the trigger maker should run before the tender / correction
task in case of simulations.

# Simulation of noise

The trigger maker supports generating noise when calculating FastOR signals from 
smeared FEE signals. Noise is simulated for each FastOR (run2 - except PHOS hole)
and every event, either as constant value or via a gaussian model. In order to switch
on noise simulation one needs to call

~~~.{cxx}
// for constant noise
triggermaker->SetConstNoiseFEESmear(0.005);   // Add 5 MeV noise to each FastOR
// for gaussian noise
triggermaker->SetGaussianNoiseFEESmear(0.005, 0.001); // Add noise simulated as gaussian noise with mean 5 MeV and sigma 1 MeV
~~~

# Selection of trigger patches by the user

Trigger patches are attached to the input event as new list object of type TClonesArray. The list object
has the name *EmcalTriggers*. The AliAnalysisTaskEmcal can import the object automatically once the list 
object is specified in the constructor:

~~~{.cxx}
this->SetCaloTriggerPatchInfoName("EmcalTriggers");
~~~

The following example code selects events with EMCAL gamma trigger patches with a minimum FEE energy of
12 GeV:

~~~{.cxx}
TClonesArray *triggerpatches = fInputEvent->FindListObject("EmcalTriggers");
bool triggered = false;
for(auto p : *triggerpatches){
  AliEMCALTriggerPatchInfo *recpatch = static_cast<AliEMCALTriggerPatchInfo *>(p);
  if(!recpatch->IsEMCAL()) continue;            // patch on EMCAL site                
  if(!recpatch->IsGammaLowSimple()) continue;   // patch is a gamma (2x2) patch, reconstructed from cells
  if(recpatch->GetPatchE() > 12.){
    triggered = true;
    break;
  }
}
~~~

## Selection of triggered events using trigger patches

In some case, i. e. in simulation or in order to select noise, EMCAL-triggered events need to be selected
using trigger patches. The wagon *PWG::EMCAL::AliAnalysisTaskEmcalTriggerSelection* performs the selection
of EMCAL-triggered events based on period-dependent configurations and adds a container with all the trigger
decisions for the various supported Level1-triggers to the event. Users can access the trigger decision
container in their task and use it for their event selection.

### 1. Add task to the analysis

In a simple user analysis do

~~~{.cxx}
gROOT->LoadMacro(Form("%s/PWG/EMCAL/AddEmcalTriggerSelectionTask.C", gSystem->Getenv("ALICE_PHYSICS)));
PWG::EMCAL::AliAnalysisTaskEmcalTriggerSelection *trgseltask = AddEmcalTriggerSelectionTask();
trgseltask->SetGlobalDecisionContainerName("EmcalTriggerDecision");
trgseltask->AutoConfigure(kDataset);
~~~

Within the train wagon, call the macro *PWG/EMCAL/AddEmcalTriggerSelectionTask.C* with the arguments (no arguments).
In the wagon configuration add

~~~{.cxx}
__R_ADDTASK__->SetGlobalDecisionContainerName("EmcalTriggerDecision");
__R_ADDTASK__->AutoConfigure(kDataset);
~~~

kDataset should be the name of a supported dataset. Please refer to the description of class PWG::EMCAL::AliAnalysisTaskEmcalTriggerSelection
for supported datasets.

### 2. Use the trigger selection in your analysis

In case the trigger selection task is added to the analysis, users can access it and use it for event selection. The
trigger seletion results are store in a trigger decision container (PWG::EMCAL::AliEmcalTriggerDecisionContainer). The
trigger decision container hold the trigger decision objects, consisting of the final trigger decision (fired yes/no),
the maximum trigger patch firing the trigger and all trigger patches firing the trigger. If one is just interested in
the event selection result it is sufficient to query the IsEventSelected method of the trigger decision container with
the valid name of the Level1 trigger to be selected.

The following code demonstrates how to select EG1 events using the EMCAL trigger decision container

~~~{.cxx}
auto triggercont = static_cast<PWG::EMCAL::AliEmcalTriggerDecisionContainer *>(fInputEvent->FindListObject(fNameTriggerDecisionContainer));
if(!triggercont){
  AliErrorStream() <<  "Trigger decision container not found in event - not possible to select EMCAL triggers" << std::endl;    }
  return false;
if(!triggercont->IsEventSelected(fTriggerSelectionString)) return false;
~~~


## Trigger bit setup

Trigger bits are used to store online trigger selection by the STU as well as offline trigger selection. 
Trigger patches can be offline or recalc patches.

| Online               | Offline Simple       | Recalc                 | 
|----------------------|----------------------|------------------------|
| Patches firing the <br/> trigger in the STU | Recalculated patch <br/> based on FEE energy | Recalculated patch <br/> based on FastOR energy  |


Accoring to the new definition, where there exists only one trigger algorithm per L1 trigger type,  trigger 
patches are always offline simple and recalc patches at the same time and store both online ADC values and
offline calibrated FEE energy. Recalc patches with an ADC value above online threshold applied for the dataset
must fullfil per definiton the trigger requirement in the STU and should therefor be online patches as well.

A certain bit number is assigned to each patch type in a bitmap in the AliEMCALTriggerPatchInfo object. Ranges
are reserved for certain patch types:

- Bit 0 - Bit 15 for online patches
- Bit 16 - Bit 23 for recalc patches 
- Bit 24 - Bit 32 for offline patches

Note that for online bits the bits are shifted in case of data by the amount of online trigger types. The order
of the bits is defined by the configuration used in the reconstruction.

| 3-bit configuration | 5-bit configuration |
|---------------------|---------------------|
| 0) Level0           | 0) Level0           |
| 1) Gamma            | 1) Gamma-high       |
| 2) Jet              | 2) Gamma-low        |
|                     | 3) Jet-high         |
|                     | 4) Jet-low          |

The access to the information is generalized by the functions IsLevel0(), IsJetHigh(), IsJetLow(), IsGammaHigh(),
IsGammaLow(). Users should not access the trigger bits directly but use the helper functions. In case of a
reconstruction with the 3-bit (referred to as *old*) configuration jet patches will be both jet-high and jet-low
patches at the same time, similarly for gamma patches. In case of the 5-bit (referred to as *new*) configuration
it is mandatory to distinguish between high and low threshold.

Offline and recalc bits are set by the trigger maker and independent of the software version used in the reconstruction.
The order of the bits follows the 5-bit configuration. Also for offline and recalc patches general getter functions
are available making the check for a patch type easy for the user (IsJetLowSimple(), IsJetHighSimple(), IsGammaLowSimple(),
IsGammaHighSimple() in case of offline simple patches and IsJetLowRecalc(), IsJetHighRecalc(), IsGammaLowRecalc() and
IsGammaHighRecalc() in case of recalc patches). For both offline simple and recalc patches the set of high-threshold
patches is automatically a subset of the set of low-threshold patches as the trigger decision to a patch is applied
in case the patch ADC is above threshold. Thresholds can be set via

~~~{.cxx}
AliEmcalTriggerMakerTask *trgtask = new AliEmcalTriggerMakerTask("NewTriggerMaker", kTRUE);
trgtask->SetTriggerThresholdJetLow(0, 0, 120);      // V0-centrality independent threshold at 120 ADC counts
~~~

In case the thresholds are set to 0, no trigger decision to patches is applied, and the output contains all possible
gamma (2x2) and jet (8x8 or 16x16) patches. In this case the user is responsible for the trigger selection.
*/