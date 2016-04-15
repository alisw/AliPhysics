/**
 * @file   GRP.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Wed Aug 20 10:00:39 2014
 * 
 * @brief  Class that hold summary of GRP data.
 * 
 * 
 */
// #include <fstream>

/**
 * Class that hold summary of GRP data, and has methods to obtain
 * these, either from a previously dumped file, or directly from OCDB.
 * 
 * This is used by the simulation setup to ensure we have the proper
 * beam settings for an anchor run. 
 */
struct GRPData
{
  /** 
   * A beam 
   */
  struct Beam
  {
    UInt_t z;
    UInt_t a;
    /** 
     * Get the per nucleon beam energy given a per charge beam energy
     * 
     * @param e Per charge beam energy 
     * 
     * @return Per nucleon beam energy 
     */
    Float_t GetEnergy(Float_t e) const 
    {
      return e * float(z)/a;
    }
    /** 
     * Get the per nucleon beam momentum given a per charge beam energy
     * 
     * @param e Per charge beam energy 
     * 
     * @return Per nucleon beam momentum
     */
    Float_t GetMomentum(Float_t e) const
    {
      Double_t be = GetEnergy(e);
      Double_t mp =TDatabasePDG::Instance()->GetParticle(2212)->Mass(); 
      return TMath::Sqrt(be*be-mp*mp);
    }
    /** 
     * Check if this is a proton
     * 
     * @return True if z and a are 1 
     */
    Bool_t IsP() const { return z == 1 && a == 1; }
    /** 
     * Check if this is a nucleus 
     * 
     * @return true if not P
     */
    Bool_t IsA() const { return !IsP(); }
    /** 
     * Get the beam type name 
     * 
     * @return P or A
     */
    const char* Name() { return (IsP() ? "P" : IsA() ? "A" : ""); }
    /** 
     * Set the beam type from single beam spec from LHC 
     * 
     * @param b Beam type from LHC 
     */
    void Set(UInt_t b)
    {
      z = b % 1000;
      a = b / 1000;
    }
    /** 
     * Set from either a full LHC spec or from ALICE spec 
     * 
     * @param b Beam
     */
    void Set(const TString& b) 
    {
      b.ToUpper();
      a = 0; 
      z = 0;
      if (b.EqualTo("P") || b.EqualTo("PROTON")) { 
	z = a = 1;
      }
      else if (b.EqualTo("A") || b.BeginsWidth("PB")) {
	z = 82;
	a = 208;
      }
    }
  };
  UInt_t  beamEnergy; // Total energy in center of mass
  UInt_t  energy; // Center of mass energy per nucleon
  TString period; // The period 
  UInt_t  run;   // The run number 
  Beam    beam1; // Target beam 
  Beam    beam2; // Projectile beam
  /** 
   * Constructor. 
   * 
   * @param r Run number
   */
  GRPData(UInt_t r) 
    : beamEnergy(0), energy(0), period(""), beam1(), beam2(), run(r)
  {
    Info("GRP", "Will try from file %s", FileName());
    if (!ReadFromFile()) { 
      Warning("GRP", "Failed to read from file, trying from OCDB");
      if (!ReadFromOCDB(run)) 
	Fatal("GRP", "Failed to get GRP values");
    }
    Print();
  }
  /** 
   * Print information 
   * 
   */
  void Print()
  {
    Printf("%s run %d\n"
	   "  Beam energy: %d GeV\n"
	   "  sqrt(s_NN):  %d GeV\n"
	   "  Beam 1:      %s (%d %d)\n"
	   "  Beam 2:      %s (%d %d)",
	   period.Data(), run, beamEnergy, energy, 
	   beam1.Name(), beam1.a, beam1.z,
	   beam2.Name(), beam2.a, beam2.z);
  }
  /** 
   * Check if this run was a p-p run
   * 
   * @return Return true if both beams are p
   */
  Bool_t IsPP() const { return beam1.IsP() && beam2.IsP(); }
  /** 
   * Check if this run was a p-A run
   * 
   * @return Return true beam 1 is p and 2 is A
   */
  Bool_t IsPA() const { return beam1.IsP() && beam2.IsA(); }
  /** 
   * Check if this run was a A-p run
   * 
   * @return Return true beam 1 is A and 2 is p
   */
  Bool_t IsAP() const { return beam1.IsA() && beam2.IsP(); }
  /** 
   * Check if this run was a A-A run
   * 
   * @return Return true if both beams are A
   */
  Bool_t IsAA() const { return beam1.IsA() && beam2.IsA(); }

  UInt_t CMEnergy(UShort_t how=0) const
  {
    switch (how) {
    case 1:
      return UInt_t(2*TMath::Sqrt(Float_t(beam1.z*beam2.z)/beam1.a/beam2.a)
		    * beamEnergy);
    case 2:
      return
	UInt_t(TMath::Sqrt(TMath::Power(BeamEnergy(1)+BeamEnergy(2),2)-
			   TMath::Power(BeamMomentum(1)-BeamMomentum(2),2)));

    }
    return energy;
  }
  Double_t BeamEnergy(UShort_t which)
  {
    switch (which) {
    case 1: return beam1.GetEnergy(beamEnergy);
    case 2: return beam2.GetEnergy(beamEnergy);
    }
    return beamEnergy;
  }
  Double_t BeamMomentum(UShort_t which)
  {
    switch (which) {
    case 1: return beam1.GetMomentum(beamEnergy);
    case 2: return beam2.GetMomentum(beamEnergy);
    }
    return beamEnergy;
  }
  /** 
   * Get the year 
   *
   * @return Year 
   */
  UInt_t Year() const 
  {
    TString tmp = period(3,2);
    return tmp.Atoi() + 2000;
  }
  const char* FileName() const { return "grp.dat"; }
  /** 
   * Get unsigned int value from string value of TObjString 
   * 
   * @param o Object
   * 
   * @return value 
   */
  static UInt_t Obj2UInt(const TObject* o) 
  {
    return (static_cast<TObjString*>(o))->String().Atoi();
  }
  /** 
   * Read values from a file 
   * 
   * @return true on success
   */
  Bool_t ReadFromFile() 
  {
    if (gSystem->AccessPathName(FileName())) {
      Info("GRP", "Cannot open file %s", FileName());
      return false;
    }
    
    std::ifstream* pin = new std::ifstream("grp.dat");
    std::ifstream& in  = *pin;
    TString line;
    TString env;
    do {
      line.ReadLine(in);
      if (line.IsNull()) continue;
      if (line[0] == '#') continue;
      env = line;
      break;
    } while (!in.eof());
    if (env.IsNull()) {
      Info("GRP", "No line found");
      return false;
    }
    Info("GRP", "Got the line %s", env.Data());
    TObjArray* tokens = env.Tokenize(":");
    if (tokens->GetEntries() < 8) {
      Warning("GRP", "Failed to get enough field from GRP_SUMMARY=%s",
	      env.Data());
      return false;
    }
    period     = tokens->At(0)->GetName();
    run        = Obj2UInt(tokens->At(1));
    beamEnergy = Obj2UInt(tokens->At(2));
    energy     = Obj2UInt(tokens->At(3));
    beam1.a    = Obj2UInt(tokens->At(4));
    beam1.z    = Obj2UInt(tokens->At(5));
    beam2.a    = Obj2UInt(tokens->At(6));
    beam2.z    = Obj2UInt(tokens->At(7));
    pin->close();
    delete pin;
    return true;
  }
  /** 
   * Read values from OCDB.  Also dumps values to file. 
   * 
   * @param r run number 
   * 
   * @return true on success
   */
  Bool_t ReadFromOCDB(UInt_t r)
  {
     AliCDBManager* cdb = AliCDBManager::Instance();
     cdb->SetDefaultStorageFromRun(r);
     cdb->SetRun(r);
     AliCDBEntry*   ent = cdb->Get("GRP/GRP/Data");
     if (!ent) {
       Warning("GRP", "No GRP data found");
       cdb->SetRun(-1);
       return false;
     }
     // Info("GRP", "Got GRP");
     // ent->PrintMetaData();

     AliGRPObject*  obj        = static_cast<AliGRPObject*>(ent->GetObject());
     // obj->Dump();
     run                       = r;
     period                    = obj->GetLHCPeriod();
     beamEnergy                = TMath::Ceil(obj->GetBeamEnergy());
     TString        beam1T     = obj->GetSingleBeamType(0);
     TString        beam2T     = obj->GetSingleBeamType(1);

     if (!beam1T.IsNull() && !beam2T.IsNull()) {
       beam1.Set(beam1T.Atoi());
       beam2.Set(beam2T.Atoi());
     }
     else {
       TString beamType = obj->GetBeamType();
       if (beamType.IsNull()) {
	 Warning("GRP", "No beam type");
	 cdb->SetRun(-1);
	 return false;
       }
       
       TObjArray* ab  = beamType.Tokenize("-");
       ab->ls();
       beam1T = ab->At(0)->GetName();
       beam2T = ab->At(1)->GetName();
       beam1.Set(beam1T);
       beam2.Set(beam2T);
     }
     Info("", "Energy=%d",  beamEnergy);
     if (beam1.IsA() && beam2.IsA()) {
       // Revert calculation done in AliGRPObject - sigh
       beamEnergy *= 208./82;
       Info("", "Beam energy now %d", beamEnergy);
     }
     // Massage the beam energy in case we had sqrt{s_NN} instead of
     // beam energy.
     if (TMath::Abs(beamEnergy - 1379.8) < 10 && beam1.IsA() && beam2.IsA()) 
       beamEnergy = 3500;
     Float_t  sNN = 2*beamEnergy *
       TMath::Sqrt(Float_t(beam1.z*beam2.z)/(beam1.a*beam2.a));
     energy = (TMath::Abs(sNN -  2760) < 10 ?  2760 :
	       TMath::Abs(sNN -  5023) < 10 ?  5023 :
	       TMath::Abs(sNN -  2360) < 10 ?  2360 :
	       TMath::Abs(sNN -   900) < 10 ?   900 :
	       TMath::Abs(sNN -  7000) < 10 ?  7000 :
	       TMath::Abs(sNN -  8000) < 10 ?  8000 :
	       TMath::Abs(sNN - 13000) < 10 ? 13000 :
	       0);
     
     const char* fn = FileName();
     std::ofstream* pout = new std::ofstream(fn);
     std::ofstream& out  = *pout;
     out << "# PERIOD:RUN:BEAMENERGY:ENERGY:BEAM1A:BEAM1Z:BEAM2A:BEAM2Z\n" 
	 << period     << ":" 
	 << run        << ":"
	 << beamEnergy << ":" 
	 << energy     << ":" 
	 << beam1.a    << ":"
	 << beam1.z    << ":"
	 << beam2.a    << ":"
	 << beam2.z    << "\n"
	 << "# EOF "   << std::endl;
     pout->close();
     delete pout;
     cdb->SetRun(-1);
     return true;
  }
  /** 
   * Get the value of an environment variable as a unsigned int 
   * 
   * @param envName Enviroment variable name 
   * @param def     Default value 
   *
   * @return As unsigned int, or default
   */
  static UInt_t Env2UInt(const char* envName, UInt_t def)
  {
    TString val(gSystem->Getenv(envName));
    if (val.IsNull()) return def;
    return UInt_t(val.Atoll());
  }
};
GRPData* grp;
void GRP(UInt_t run=0)
{
  grp = new GRPData(run);
}
// 
// EOF
// 



  
