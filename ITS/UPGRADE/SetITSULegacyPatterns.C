void SetITSULegacyPatterns()
{
  // set the ITSU patterns to "legacy mode" to be able to read old geometry
  gSystem->Load("libITSUpgradeBase");
  AliITSUGeomTGeo::SetITSVolPattern("ITSV");
  AliITSUGeomTGeo::SetITSLayerPattern("ITSULayer");
  AliITSUGeomTGeo::SetITSWrapVolPattern("ITSUWrapVol");
  AliITSUGeomTGeo::SetITSStavePattern("ITSULadder");
  AliITSUGeomTGeo::SetITSHalfStavePattern("");
  AliITSUGeomTGeo::SetITSModulePattern("");
  AliITSUGeomTGeo::SetITSChipPattern("ITSUModule");
  AliITSUGeomTGeo::SetITSSensorPattern("ITSUSensor");
  //
}
