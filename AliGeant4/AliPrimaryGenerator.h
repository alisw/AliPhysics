// $Id$

#ifndef ALI_PRIMARY_GENERATOR_H
#define ALI_PRIMARY_GENERATOR_H

// available primary generators

enum AliPrimaryGenerator {
    kGun,               // gun (can be set interactively) 
    kGeantino,          // geantino with random direction
    kChargedGeantino,   // chargedgeantino with random direction
    kAliGenerator       // AliGenerator from AliRoot
};  

#endif //ALI_PRIMARY_GENERATOR_H
