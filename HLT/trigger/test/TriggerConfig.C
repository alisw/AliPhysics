// $Id:  $
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Artur Szostak <artursz@iafrica.com>                   *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/// @file   TriggerConfig.C
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   28 Oct 2009
/// @brief  Implementation of the global trigger menu configuration routines.
///
/// The TriggerConfig routine generates a global trigger menu configuration used
/// to test the AliHLTGlobalTriggerComponent class.

/**
 * Initialise a trigger menu to test a configuration where all trigger menu entries
 * have different priority values. Thus only the first trigger menu item that
 * matches is activated.
 */
void PriorityGroupTestConfig()
{
	AliHLTGlobalTriggerConfig config("Priority group test config");
	config.AddItem(3, "triggerTPC && triggerSSD", "triggerSSD", "Fast SSD trigger");
	config.AddItem(2, "triggerMUON || triggerSSD", "triggerMUON | triggerSSD", "MUON trigger");
	config.AddItem(1, "triggerTPC", "triggerTPC | triggerMUON | triggerSSD", "TPC trigger");
	config.SetDefaultTriggerDescription("No trigger");
}

/**
 * Initialise a trigger menu to test a configuration where all trigger menu entries
 * have the same priority values. Thus all entries form a single priority group and
 * multiple items can be matched with the trigger input.
 */
void SingleGroupTestConfig()
{
	AliHLTGlobalTriggerConfig config("Single group test config");
	config.AddItem("triggerTPC", "triggerTPC", "TPC trigger");
	config.AddItem("triggerMUON", "triggerMUON", "MUON trigger");
	config.AddItem("triggerSSD", "triggerSSD", "SSD trigger");
	config.SetDefaultConditionOperator("||");
	config.SetDefaultDomainOperator("|");
}

/**
 * Create a trigger menu configuration to test the prescalar mechanism.
 */
void PrescalarTestConfig()
{
	AliHLTGlobalTriggerConfig config("Prescalar test config");
	config.AddItem(2, "triggerTPC", "triggerTPC", 3, "TPC trigger");
	config.AddItem(1, "triggerMUON", "triggerMUON", "MUON trigger");
	config.SetDefaultConditionOperator("||");
	config.SetDefaultTriggerDescription("No trigger");
	config.SetDefaultTriggerDomain(AliHLTTriggerDomain("*******:HLT "));
}

/**
 * Create a trigger menu configuration to test the trigger menu symbol mechanism.
 */
void SymbolTestConfig()
{
	AliHLTGlobalTriggerConfig config("Symbol test config");
	config.AddSymbol("domain-All", "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"*******:***,-DAQRDOUT:TST\")");
	config.AddSymbol("PHOSclusters", "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"CLUSTERS:PHOS\")");
	config.AddSymbol("PHOStracks", "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"TRACKS:PHOS\")");
	config.AddSymbol("domainPHOS", "AliHLTTriggerDomain", "", "PHOSclusters | PHOStracks");
	config.AddSymbol("triggerTPC", "bool", "this->Result()", "false", "AliHLTTriggerDecision");
	config.AddSymbol("trigClasses", "Double_t", "this->GetScalar(\"TrigClass\")", "0", "AliHLTScalars");
	config.AddItem(2, "true", "domain-All", 4, "Pass through");
	config.AddItem(1, "trigClasses == 0x2", "triggerTPC | domainPHOS", "Trigger class 2");
}

/**
 * This will make a complex configuration that uses multiple features of the
 * trigger menu configuration system.
 */
void ComplexTestConfig()
{
	AliHLTGlobalTriggerConfig config("Complex test config");
	config.AddSymbol("domain-All", "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"*******:***,-DAQRDOUT:TST\")");
	config.AddSymbol("PHOSclusters", "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"CLUSTERS:PHOS\")");
	config.AddSymbol("PHOStracks", "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"TRACKS:PHOS\")");
	config.AddSymbol("domainPHOS", "AliHLTTriggerDomain", "", "PHOSclusters | PHOStracks");
	config.AddItem(4, "true", "domain-All", 7, "Pass through");
	config.AddItem(3, "triggerSSD", "triggerSSD |", "SSD trigger 1");
	config.AddItem(3, "triggerMUON", "triggerMUON", "MUON trigger 1");
	config.AddItem(2, "triggerSSD ||", "triggerSSD | domainPHOS", "SSD trigger 2");
	config.AddItem(2, "triggerMUON ||", "triggerMUON", "MUON trigger 2");
	config.AddItem(1, "triggerTPC", "triggerTPC | triggerSSD | triggerMUON", "Slow trigger");
	config.SetDefaultTriggerDomain(AliHLTTriggerDomain("*******:MUON"));
	config.SetDefaultConditionOperator("&&");
}

/**
 * This will make a configuration that will test software triggers.
 */
void SoftwareTriggersTestConfig()
{
	AliHLTGlobalTriggerConfig config("Software triggers test config");
	config.AddSymbol("domainPHOS", "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"*******:PHOS\")");
	config.AddSymbol("domainSPD", "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"*******:SPD\")");
	config.AddItem(5, "START_OF_DATA", "START_OF_DATA", "Start of data");
	config.AddItem(4, "END_OF_DATA", "END_OF_DATA", "End of data");
	config.AddItem(3, "SOFTWARE", "domainSPD | SOFTWARE", "Software trigger");
	config.AddItem(2, "CALIBRATION", "domainPHOS", "Calibration trigger");
	config.AddItem(1, "triggerMUON", "triggerMUON", "MUON trigger");
}

/**
 * Top level configuration routine for the global trigger component tests.
 */
void TriggerConfig(int version = 0)
{
	switch (version)
	{
	case 0: PriorityGroupTestConfig(); break;
	case 1: SingleGroupTestConfig(); break;
	case 2: PrescalarTestConfig(); break;
	case 3: SymbolTestConfig(); break;
	case 4: ComplexTestConfig(); break;
	case 5: SoftwareTriggersTestConfig(); break;
	default: AliHLTGlobalTriggerConfig config("Empty test config");
	}
}
