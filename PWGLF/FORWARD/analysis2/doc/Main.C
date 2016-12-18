/** 
 * @if pwglf_fwd_main
 * @mainpage ALICE PWGLF Forward Multiplcity Analysis 
 * @else 
 * @page PWGLF_FORWARD PWGLF Forward Multiplcity Analysis 
 * @endif
 * 
 * This is the analysis code for analysis of the Forward data. 
 * 
 * @par Sections 
 *
 * - @subpage pwglf_fwd_mult_doc
 * - @subpage pwglf_fwd_density_doc
 * - @subpage pwglf_fwd_dist_doc
 * - @subpage pwglf_fwd_flow_doc
 * - @subpage pwglf_fwd_mid_doc
 * - @subpage pwglf_forward_tracklets
 *   - @subpage pwglf_fwd_spd_tracklet_1
 *   - @subpage pwglf_fwd_spd_tracklet_2
 * - @subpage pwglf_fwd_sim
 * - @subpage pwglf_fwd_fastsim
 * 
 * @par External Information 
 *
 * - <a href="https://aliceinfo.cern.ch/Notes/node/107">Pb-Pb 2.76TeV analysis note</a> (satellite-main collisions)
 * - <a href="https://aliceinfo.cern.ch/Notes/node/384">Pb-Pb 2.76TeV analysis note</a> (all centralities)
 * - <a href="https://aliceinfo.cern.ch/Notes/node/455">Pb-Pb 5.02TeV analysis note</a>
 *
 * @par Observations on implementing tasks 
 *
 * - Any object propagated to the output must be allocated on the heap
 *   and not deleted by the task.  The reason is, that PROOF cleans-up
 *   tasks _before_ the output is flushed to disk
 */

/** 
 * @defgroup pwglf_forward PWGLF Forward analysis
 *
 * Code to do the multiplicity analysis in the forward psuedo-rapidity
 * regions
 *
 */
/** 
 * @defgroup pwglf_forward_tasks Tasks
 *
 * Code to do the multiplicity analysis in the forward psuedo-rapidity
 * regions
 *
 * @ingroup pwglf_forward 
 */
/** 
 * @defgroup pwglf_forward_topical Topical
 *
 * The code divided according to topic
 *
 * @ingroup pwglf_forward 
 */
