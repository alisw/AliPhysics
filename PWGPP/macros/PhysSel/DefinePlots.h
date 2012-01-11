enum {
  kGraphACCOverAllAC,
  kGraphACCOverAllACS2,
  kGraphACCOverAllSMH,

  kGraphACCOverAllACRej,
  kGraphACCOverAllACS2Rej,
  kGraphACCOverAllSMHRej,

  kGraphV0BGOverAllAC,
  kGraphV0BGOverAllAC2,
  kGraphV0BGOverAllSMH,

  kGraphV0BGOverAllACRej,
  kGraphV0BGOverAllAC2Rej,
  kGraphV0BGOverAllSMHRej,


  kGraphNevAC,
  kGraphNevAC2, 
  kGraphNevSMH,

  kGraphNevACBefore, 
  kGraphNevAC2Before,
  kGraphNevSMHBefore,

  kGraphNevACAfter, 
  kGraphNevAC2After,
  kGraphNevSMHAfter,

  kGraphNevACratioAfter, 
  kGraphNevAC2ratioAfter,
  kGraphNevSMHratioAfter,

  kGraphNevACratioAfterRej, 
  kGraphNevAC2ratioAfterRej,
  kGraphNevSMHratioAfterRej,


  kGraphBGOverAllAC,
  kGraphBGOverAllAC2,
  kGraphBGOverAllSMH,

  kGraphBGOverAllACRej,
  kGraphBGOverAllAC2Rej,
  kGraphBGOverAllSMHRej,

  kGraphV0AOverV0CAC,
  kGraphV0AOverV0CAC2,
  kGraphV0AOverV0CSMH,

  kGraphV0AOverV0CACRej,
  kGraphV0AOverV0CAC2Rej,
  kGraphV0AOverV0CSMHRej,


  kGraphV0ABGOverAllAC,
  kGraphV0ABGOverAllAC2,
  kGraphV0ABGOverAllSMH,

  kGraphV0ABGOverAllACRej,
  kGraphV0ABGOverAllAC2Rej,
  kGraphV0ABGOverAllSMHRej,

  kGraphV0CBGOverAllAC,
  kGraphV0CBGOverAllAC2,
  kGraphV0CBGOverAllSMH,

  kGraphV0CBGOverAllACRej,
  kGraphV0CBGOverAllAC2Rej,
  kGraphV0CBGOverAllSMHRej,

  kGraphFO1OverAllAC,
  kGraphFO1OverAllAC2,
  kGraphFO1OverAllSMH,

  kGraphFO1OverAllACRej,
  kGraphFO1OverAllAC2Rej,
  kGraphFO1OverAllSMHRej,

  kGraphFO2OverAllAC,
  kGraphFO2OverAllAC2,
  kGraphFO2OverAllSMH,

  kGraphFO2OverAllACRej,
  kGraphFO2OverAllAC2Rej,
  kGraphFO2OverAllSMHRej,





  kGraphNev90, 
  kGraphTemp,
  kNGraphs
};
const char * ylabels[] = {
  "Accepted / All [CEMC7-B]",
  "Accepted / All [CMUSH7-B]",
  "Accepted / All [CINT7-B]",

  "Accepted / All Rej [CEMC7-B]",
  "Accepted / All Rej [CMUSH7-B]",
  "Accepted / All Rej [CINT7-B]",

  "V0BG / All [CEMC7-B]", 
  "V0BG / All [CMUSH7-B]", 
  "V0BG / All [CINT7-B]", 

  "V0BG / All Rej [CEMC7-B]", 
  "V0BG / All Rej [CMUSH7-B]", 
  "V0BG / All Rej [CINT7-B]", 


  "Nev [Accepted, CEMC7-B]" ,
  "Nev [Accepted, CMUSH7-B]" ,
  "Nev [Accepted, CINT7-B]" ,

  "Nev [All Before, CEMC7-B]" ,
  "Nev [All Before, CMUSH7-B]" ,
  "Nev [All Before, CINT7-B]" ,
  
  "Nev [All After, CEMC7-B]" ,
  "Nev [All After, CMUSH7-B]" ,
  "Nev [All After, CINT7-B]" ,

  "Nev_AFS / Nev_BPS [CEMC7-B]" ,
  "Nev_AFS / Nev_BPS [CMUSH7-B]" ,
  "Nev_AFS / Nev_BPS [CINT7-B]" ,

  "Nev_AFS / Nev_BPS  [CEMC7-B] Rej" ,
  "Nev_AFS / Nev_BPS  [CMUSH7-B] Rej" ,
  "Nev_AFS / Nev_BPS  [CINT7-B] Rej" ,

  "BG / All [CEMC7-B]", 
  "BG / All [CMUSH7-B]", 
  "BG / All [CINT7-B]", 

  "BG / All Rej [CEMC7-B]", 
  "BG / All Rej [CMUSH7-B]", 
  "BG / All Rej [CINT7-B]", 

  "Nev_V0A / Nev_V0C [CEMC7-B]",
  "Nev_V0A / Nev_V0C [CMUSH7-B]",
  "Nev_V0A / Nev_V0C [CINT7-B]",

  "Nev_V0A / Nev_V0C  [CEMC7-B] Rej",
  "Nev_V0A / Nev_V0C  [CMUSH7-B] Rej",
  "Nev_V0A / Nev_V0C  [CINT7-B] Rej",

  "BGV0A / All [CEMC7-B]",
  "BGV0A / All [CMUSH7-B]",
  "BGV0A / All [CINT7-B]",

  "BGV0A / All  [CEMC7-B] Rej",
  "BGV0A / All  [CMUSH7-B] Rej",
  "BGV0A / All  [CINT7-B] Rej",

  "BGV0C / All [CEMC7-B]",
  "BGV0C / All [CMUSH7-B]",
  "BGV0C / All [CINT7-B]",

  "BGV0C / All  [CEMC7-B] Rej",
  "BGV0C / All  [CMUSH7-B] Rej",
  "BGV0C / All  [CINT7-B] Rej",

  "FO1 / All [CEMC7-B]",
  "FO1 / All [CMUSH7-B]",
  "FO1 / All [CINT7-B]",

  "FO1 / All  [CEMC7-B] Rej",
  "FO1 / All  [CMUSH7-B] Rej",
  "FO1 / All  [CINT7-B] Rej",


  "FO2 / All [CEMC7-B]",
  "FO2 / All [CMUSH7-B]",
  "FO2 / All [CINT7-B]",

  "FO2 / All  [CEMC7-B] Rej",
  "FO2 / All  [CMUSH7-B] Rej",
  "FO2 / All  [CINT7-B] Rej",


  "Events above 90% (centrality fit V0)",
  "Temp"
};
const char * gnames[] = {
  "grGraphACCOverAllAC",
  "grGraphACCOverAllACS2",
  "grGraphACCOverAllSMH", 

  "grGraphACCOverAllACRej", 
  "grGraphACCOverAllACS2Rej", 
  "grGraphACCOverAllSMHRej",

  "grGraphV0BGOverAllAC",
  "grGraphV0BGOverAllAC2",
  "grGraphV0BGOverAllSMH",

  "grGraphV0BGOverAllACRej",
  "grGraphV0BGOverAllAC2Rej",
  "grGraphV0BGOverAllSMHRej",

  "grGraphNevAC",
  "grGraphNevAC2",
  "grGraphNevSMH",


  "grGraphNevACAll",
  "grGraphNevAC2All",
  "grGraphNevSMHAll",

  "grGraphNevACAfter",
  "grGraphNevAC2After",
  "grGraphNevSMHAfter",

  "grGraphNevACratioAfter",
  "grGraphNevAC2ratioAfter",
  "grGraphNevSMHratioAfter",

  "grGraphNevACratioAfterRej",
  "grGraphNevAC2ratioAfterRej",
  "grGraphNevSMHratioAfterRej",

  "grGraphBGOverAllAC",
  "grGraphBGOverAllAC2",
  "grGraphBGOverAllSMH",

  "grGraphBGOverAllACRej",
  "grGraphBGOverAllAC2Rej",
  "grGraphBGOverAllSMHRej",

  "grGraphV0AOverV0CAC",
  "grGraphV0AOverV0CAC2",
  "grGraphV0AOverV0CSMH",

  "grGraphV0AOverV0CACRej",
  "grGraphV0AOverV0CAC2Rej",
  "grGraphV0AOverV0CSMHRej",


  "grGraphV0ABGOverAllAC",
  "grGraphV0ABGOverAllAC2",
  "grGraphV0ABGOverAllSMH",

  "grGraphV0ABGOverAllACRej",
  "grGraphV0ABGOverAllAC2Rej",
  "grGraphV0ABGOverAllSMHRej",


  "grGraphV0CBGOverAllAC",
  "grGraphV0CBGOverAllAC2",
  "grGraphV0CBGOverAllSMH",

  "grGraphV0CBGOverAllACRej",
  "grGraphV0CBGOverAllAC2Rej",
  "grGraphV0CBGOverAllSMHRej",

  "grGraphFO1OverAllAC",
  "grGraphFO1OverAllAC2",
  "grGraphFO1OverAllSMH",

  "grGraphFO1OverAllACRej",
  "grGraphFO1OverAllAC2Rej",
  "grGraphFO1OverAllSMHRej",


  "grGraphFO2OverAllAC",
  "grGraphFO2OverAllAC2",
  "grGraphFO2OverAllSMH",

  "grGraphFO2OverAllACRej",
  "grGraphFO2OverAllAC2Rej",
  "grGraphFO2OverAllSMHRej",


  "grNev90",
  "grTemp"
};
