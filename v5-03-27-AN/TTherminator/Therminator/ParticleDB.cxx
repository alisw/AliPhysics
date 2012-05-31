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
 * to the description of the model plesase refer to: nucl-th/0011222,         *
 * accessibile at: http://www.arxiv.org/nucl-th/0011222                       *
 *                                                                            *
 * Homepage: http://hirg.if.pw.edu.pl/en/therminator/                         *
 *                                                                            *
 * This code can be freely used and redistributed. However if you decide to   *
 * make modifications to the code, please contact the authors, especially     *
 * if you plan to publish the results obtained with such modified code.       *
 * Any publication of results obtained using this code must include the       *
 * reference to nucl-ex/0011222 and the published version of it, when         *
 * available.                                                                 *
 *                                                                            *
 *****************************************************************************/
#include "ParticleDB.h"

ParticleDB::ParticleDB()
{
  mParticleTable.clear();
  mParticleNames.clear();
}

ParticleDB::~ParticleDB()
{
}

int           
ParticleDB::AddParticleType(ParticleType *aPartType)
{
  mParticleTable.push_back(*aPartType);
  mParticleNames[aPartType->GetName()] = mParticleTable.size()-1;
  return  mParticleTable.size()-1;
}

ParticleType* 
ParticleDB::GetParticleType(int aIndex)
{
  return &(mParticleTable[aIndex]);
}

ParticleType* 
ParticleDB::GetParticleType(std::string aName)
{
  return &(mParticleTable[mParticleNames[aName]]);
}

int           
ParticleDB::GetParticleTypeIndex(std::string aName)
{
  return mParticleNames[aName];
}

int           
ParticleDB::GetParticleTypeCount()
{
  return mParticleTable.size();
}

int           
ParticleDB::ExistsParticleType(std::string aName)
{
  return mParticleNames.count(aName);
}
