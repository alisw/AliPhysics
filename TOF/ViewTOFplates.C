void ViewTOFplates()
{
//=====> Level 1
  // Level 1 for TOF volumes
  gMC->Gsatt("B077","seen",0);


//==========> Level 2
  // Level 2
  gMC->Gsatt("B076","seen",-1); // all B076 sub-levels skipped -
  gMC->Gsatt("B071","seen",0);
  gMC->Gsatt("B074","seen",0);
  gMC->Gsatt("B075","seen",0);
  gMC->Gsatt("B080","seen",0); // B080 does not has sub-level

  // Level 2 of B071
  gMC->Gsatt("B063","seen",-1); // all B063 sub-levels skipped   -
  gMC->Gsatt("B065","seen",-1); // all B065 sub-levels skipped   -
  gMC->Gsatt("B067","seen",-1); // all B067 sub-levels skipped   -
  gMC->Gsatt("B069","seen",-1); // all B069 sub-levels skipped   -
  gMC->Gsatt("B056","seen",0);  // B056 does not has sub-levels  -
  gMC->Gsatt("B059","seen",-1); // all B059 sub-levels skipped   -
  gMC->Gsatt("B072","seen",-1); // all B072 sub-levels skipped   -
  gMC->Gsatt("BTR1","seen",0);  // BTR1 do not have sub-levels   -
  gMC->Gsatt("BTO1","seen",0);  


  // Level 2 of B074
  gMC->Gsatt("BTR2","seen",0); // BTR2 does not has sub-levels -
  gMC->Gsatt("BTO2","seen",0);

  // Level 2 of B075
  gMC->Gsatt("BTR3","seen",0); // BTR3 do not have sub-levels -
  gMC->Gsatt("BTO3","seen",0);

// ==================> Level 3
  // Level 3 of B071 / Level 2 of BTO1
  gMC->Gsatt("FTOC","seen",-2);
  gMC->Gsatt("FTOB","seen",-2);
  gMC->Gsatt("FTOA","seen",-2);

  // Level 3 of B074 / Level 2 of BTO2
  // -> cfr previous settings

  // Level 3 of B075 / Level 2 of BTO3
  // -> cfr previous settings

  
}
