#ifndef ALIEMCALTRIGGEREMULATION_H
#define ALIEMCALTRIGGEREMULATION_H
/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TF1.h>
#include <TObject.h>

class AliEmcalContainer;

namespace PWGJE {
	
namespace EMCALJetTasks {

/**
 * @class AliEmcalTriggerEmulation
 * @brief Emulate trigger decision offline
 * @ingroup PWGJETASKS
 * @author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * @since May 4, 2016
 *
 * This class emulates the EMCAL trigger decision by selecting events with at least on
 * object containing energy (cluster, trigger patch) above energy threshold. The energy
 * threshold is evaluated event-by-event, smeared with the resolution.
 *
 * This class is intended for the usage on MC events in order to evaluate the EMCAL trigger
 * efficiency there.
 */
class AliEmcalTriggerEmulation : public TObject {
public:
  /**
   * Dummy constructor, for ROOT I/O, not for normal usage
   * Setting default values for energy threshold and resolution
   */
	AliEmcalTriggerEmulation();
	/**
	 * Constructor, setting energy threshold and resolution used for the smearing
	 * @param[in] energy Energy threshold (units should match the target object [GeV/ ADC counts])
	 * @param[in] resolution Energy resolution (same units as energy threshold)
	 */
	AliEmcalTriggerEmulation(double energy, double resolution);
	/**
	 * Destructor
	 */
	virtual ~AliEmcalTriggerEmulation() {}

	/**
	 * Select event according to objects above energy threshold. The threshold
	 * is evaluated for each event separately smearing the nominal energy
	 * threshold with a resolution. Objects are taken from the input container.
	 * The cut will be implemented for
	 * - Cluster container
	 * - Trigger patch container
	 * @param[in] cont Container with input objects (clusters / trigger patches)
	 * @return True if at least one object above threshold is found, false otherwise
	 */
	bool SelectEvent(const AliEmcalContainer * const cont) const;

	/**
	 * Set the energy threshold used in the event selection
	 * @param[in] energy Energy threshold (units should match the target object [GeV/ ADC counts])
	 */
	void SetEnergy(double energy) { fThresholdEngine.SetParameter(1, energy); }
	/**
	 * Set the energy resolution used for the smearing of the threshold
	 * @param[in] resolution Energy resolution (same units as energy threshold)
	 */
	void SetResolution(double resolution) { fThresholdEngine.SetParameter(2, resolution); }
	/**
	 * Set the energy threshold and resolution used for the smearing
	 * @param[in] energy Energy threshold (units should match the target object [GeV/ ADC counts])
	 * @param[in] resolution Energy resolution (same units as energy threshold)
	 */
	void Set(double energy, double resolution) {
	  fThresholdEngine.SetParameter(1, energy);
	  fThresholdEngine.SetParameter(2, resolution);
	}

protected:
	TF1                       fThresholdEngine;   ///< Engine for energy threshold calculation

	ClassDef(AliEmcalTriggerEmulation, 1);
};

} /* namespace EMCALJetTasks */

} /* namespace PWGJE */

#endif /* ALIEMCALTRIGGEREMULATION_H */
