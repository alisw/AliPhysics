void
ShowFMDITS(const char* name="ALIC")
{
  AliLog::SetModuleDebugLevel("FMD", 1);
  gAlice->InitMC("$ALICE_ROOT/FMD/Config.C");

  const char* inv[] = { "RB24", "RB26Pipe", 
			"ICYL", 
			"ICU0", "ICU1", "ICU2", "ICU3", "ICU4", 
			"ICU5", "ICU6", "ICU7", "ICU8", "ICU9",
			"ICC0", "ICC1", "ICC2", "ICC3", "ICC4", 
			"ICC5", "ICC6", "ICC7", "ICC8", "ICC9",
			"IHK1", "IHK2", 
			"ISR1", "ISR2", "ISR3", "ISR6",
			"ITSSPD", 
			"ITSsddLayer3", "ITSsddLayer4", 
			"ITSssdLayer5", "ITSssdLayer6", 
			"ITSsddForward3Pos", "ITSsddForward3Neg", 
			"ITSsddForward4Pos", "ITSsddForward4Neg", 
			"Lay5LadderSupportRing", "Lay6LadderSupportRing",
			"EndCapSupportSystemLayer5Sx",
			"EndCapSupportSystemLayer5Dx",
			"EndCapSupportSystemLayer6Sx",
			"EndCapSupportSystemLayer6Dx",
			"SDDCarbonFiberCylinder",
			"SDDCarbonFiberCone", 
			"SSDexternalcylinder", 
			"SSDfoamcylinder",
			"SSDCableLay5RightSidePCon1", 
			"SSDCableLay5RightSidePCon2", 
			"SSDCableLay5RightSidePCon3", 
			"SSDCableLay5RightSidePConWater1", 
			"SSDCableLay5RightSidePConWater2", 
			"SSDCableLay5RightSidePConWater3", 
			"SSDCableLay5RightSideTube", 
			"SSDCableLay5RightSideWaterTube", 
			"SSDCableLay6MaterialBudgetPCon", 
			"SSDCableLay6MaterialBudgetTube", 
			"SSDCableLay6RightSidePCon", 
			"SSDCableLay6RightSideTube", 
			"SSDCableLay6RightSideWaterPCon", 
			"SSDCableLay6RightSideWaterTube", 
			"ITSssdCone", 
			"ITScablesSDDpcon1Plast",
			"ITScablesSDDpcon3Plast",
			"vSddCableInterCyl",
			"vpcon1container",
			"vpcon2container",
			"vpcon3container",
			0 };
  const char** ptr = inv;
  TGeoVolume*  vol = 0;
  while (*ptr) {
    const char* n = *ptr++;
    vol = gGeoManager->GetVolume(n);
    if (!vol) { 
      std::cerr << "Volume " << n << " not found" << std::endl;
      continue;
    }
    // std::cout << "Processing " << n << std::endl;
    vol->SetVisDaughters(kFALSE);
    vol->SetVisContainers(kTRUE);
    vol->SetVisibility(kFALSE);
    vol->InvisibleAll(kTRUE);
  }

  vol = gGeoManager->GetVolume("SPDcentralomega");
  if (vol) vol->SetTransparency(64);
  vol = gGeoManager->GetVolume("SPDendcabomaga");
  if (vol) vol->SetTransparency(64);
  vol = gGeoManager->GetVolume("SPDconeshieldH1");
  if (vol) vol->SetTransparency(64);
  vol = gGeoManager->GetVolume("SPDconeshieldH2");
  if (vol) vol->SetTransparency(64);

  vol = gGeoManager->GetVolume(name);
  if (!vol) return;
  new TBrowser("geoBrowser", gGeoManager);
  vol->Draw("ogl");
}

