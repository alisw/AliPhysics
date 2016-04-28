/** 
 * @mainpage ALICE PWGLF Forward Multiplcity Analysis 
 * 
 * This is the analysis code for analysis of the Forward data. 
 * 
 * @par Sections 
 *
 * - @subpage train_setup_doc
 * - @subpage mult_doc
 * - @subpage density_doc
 * - @subpage dist_doc
 * - @subpage flow_doc
 * - @subpage mid_doc
 * - <a href="modules.html"><b>Modules</b></a>
 * 
 * @par External Information 
 *
 * - <a href="https://aliceinfo.cern.ch/Notes/node/25">Analysis Note</a>
 *
 * @par Observations on implementing tasks 
 *
 * - Any object propagated to the output must be allocated on the heap
 *   and not deleted by the task.  The reason is, that PROOF cleans-up
 *   tasks _before_ the output is flushed to disk
 */
