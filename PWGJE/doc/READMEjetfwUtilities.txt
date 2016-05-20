/*! \page READMEjetfwUtilities Jet Framework Utilities

# Jet Framework Utilities (e.g. FJ contrib)

Additional utilities can be attached to the AliEmcalJetTask object, in a similar fashion as it is done for the AliTender class. The utility classes have to derive from the abstract class AliEmcalJetUtility. An EMCal jet utility class can implement any of the following four methods declared as virtual in AliEmcalJetUtility:

~~~{.cxx}
virtual void Init() = 0;                                                        // Executed only once at the end of AliEmcalJetTask::DoInit()
virtual void Prepare(AliFJWrapper& fjw) = 0;                                    // Executed for each event at the beginning of AliEmcalJetTask::FillJetBranch()
virtual void ProcessJet(AliEmcalJet* jet, Int_t ij, AliFJWrapper& fjw) = 0;     // Executed for each jet in the loop in AliEmcalJetTask::FillJetBranch()
virtual void Terminate(AliFJWrapper& fjw) = 0;                                  // Executed for each event at the end of AliEmcalJetTask::FillJetBranch()
~~~

At the moment two utilities are available, which make use of the FastJet contribs: the generic subtractor implemented in AliEmcalJetUtilityGenSubtractor and the constituent subtractor implemented in AliEmcalJetUtilityConstSubtractor. For example, to add the generic subtractor to a previously defined AliEmcalJetTask object named jetTask one can use the following code.

~~~{.cxx}
AliEmcalJetUtilityGenSubtractor* genSub = jetTask->AddUtility(new AliEmcalJetUtilityGenSubtractor("GenSubtractor"));
genSub->SetUseExternalBkg(kTRUE);
genSub->SetRhoName(rhoname);
genSub->SetRhomName(rhomname);
genSub->SetGenericSubtractionJetMass(kTRUE);
~~~

or for the constituent subtraction:

~~~{.cxx}
AliEmcalJetUtilityConstSubtractor* constUtil = jetTask->AddUtility(new AliEmcalJetUtilityConstSubtractor("ConstSubtractor"));
constUtil->SetUseExternalBkg(kTRUE);
constUtil->SetRhoName(rhoname);
constUtil->SetRhomName(rhomname);
constUtil->SetJetsSubName(Form("%sConstSub",jetTask->GetName()));
constUtil->SetParticlesSubName(Form("%s_%sConstSub",kTracksName.Data(),jetTask->GetName()));
~~~

_Note:_ A description of the AliFWWrapper should be added

*/
