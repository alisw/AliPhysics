#!/bin/bash

rm -f default_geo.root
aliroot <<"FNORD"
  TFile::Open("$ALICE_ROOT/OCDB/GRP/Geometry/Data/Run0_999999999_v0_s0.root")
  gFile->Get("AliCDBEntry");
  gGeoManager->DefaultColors()
  gGeoManager->Export("default_geo.root", "", "")
  .q
FNORD
