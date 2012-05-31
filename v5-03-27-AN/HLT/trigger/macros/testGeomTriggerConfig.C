void testGeomTriggerConfig()
{
  AliHLTGlobalTriggerConfig config("Global Trigger Test Config");
  config.AddItem("BarrelGeomMultiplicityTrigger", "BarrelGeomMultiplicityTrigger", "Barrel Geom Multiplicity Trigger");
  config.Print();
}
