///////////////////////////////////////////////////////////////////////////////
// TTherminator: global interface class to the Therminator model.            //
// Initialized global variables and runs the event generation                //
///////////////////////////////////////////////////////////////////////////////
/******************************************************************************
 *                      T H E R M I N A T O R                                 *
 *                   THERMal heavy-IoN generATOR                              *
 *                           version 1.0                                      *
 *                                                                            *
 * Authors of the model: Wojciech Broniowski, Wojciech.Broniowski@ifj.edu.pl, *
 *                       Wojciech Florkowski, Wojciech.Florkowski@ifj.edu.pl  *
 * Authors of the code:  Adam Kisiel, kisiel@if.pw.edu.pl                     *
 *                       Tomasz Taluc, ttaluc@if.pw.edu.pl                    *
 * Code designers: Adam Kisiel, Tomasz Taluc, Wojciech Broniowski,            *
 *                 Wojciech Florkowski                                        *
 *                                                                            *
 * For the detailed description of the program and furhter references         * 
 * to the description of the model plesase refer to: nucl-th/0504047,         *
 * accessibile at: http://www.arxiv.org/nucl-th/0504047                       *
 *                                                                            *
 * Homepage: http://hirg.if.pw.edu.pl/en/therminator/                         *
 *                                                                            *
 * This code can be freely used and redistributed. However if you decide to   *
 * make modifications to the code, please contact the authors, especially     *
 * if you plan to publish the results obtained with such modified code.       *
 * Any publication of results obtained using this code must include the       *
 * reference to nucl-th/0504047 and the published version of it, when         *
 * available.                                                                 *
 *                                                                            *
 *****************************************************************************/
#ifndef TTHERMINATOR_H
#define TTHERMINATOR_H

#include <TGenerator.h>
#include <map>
#include <THGlobal.h>
#include "ReadPar.h"
#include "Parser.h"
#include "ParticleDB.h"
#include "ParticleType.h"
#include "DecayTable.h"
#include "Therminator/Event.h"

class TTherminator: public TGenerator {
 public:
  TTherminator();
  TTherminator(const TTherminator & therm);
  TTherminator& operator=(const TTherminator & therm);
  virtual ~TTherminator();
  
  virtual void        ReadParameters();
  
  virtual void        Initialize();
  
  virtual void        GenerateEvent();
  
  virtual Int_t       ImportParticles(TClonesArray *particles, Option_t *option="");
  virtual TObjArray*  ImportParticles(Option_t *option="");

 private:
  int Run();
  Integrator *fCalka; // Integrator class
  Event      *fEvent; // The therminator event
  ParticleDB *fPartDB;// Particle properties database
  TString     fInputDir; // Name of directory with SHARE input files

  ClassDef(TTherminator,1) // Hijing parametrisation generator
};

#endif
