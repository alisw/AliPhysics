// $Id$
// Category: run
//
// Author: I. Hrivnacova
//
// Enum AliPrimaryGenerator
// ------------------------
// Enumaration of available primary generators.

#ifndef ALI_PRIMARY_GENERATOR_H
#define ALI_PRIMARY_GENERATOR_H

// available primary generators

enum AliPrimaryGenerator {
    kGun,               // gun (can be set interactively) 
    kGeantino,          // geantino with random direction
    kChargedGeantino,   // chargedgeantino with random direction
    kStack              // AliGenerator from MC stack
};  

#endif //ALI_PRIMARY_GENERATOR_H
