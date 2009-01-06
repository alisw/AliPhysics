
void TriggerConfig()
{
	AliHLTGlobalTriggerConfig config("test Config");
	config.AddSymbol("domainAll", "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"*******:***,-DAQRDOUT:TST\")");
	config.AddItem("true", "domainAll", 5, "Trigger Type: pass through");
	config.AddItem("Trigger1", "Trigger1 | Trigger2", 3, "Trigger Type: 1");
	config.AddItem("Trigger2", "Trigger2", 0, "Trigger Type: 2");
	config.Print();
}

