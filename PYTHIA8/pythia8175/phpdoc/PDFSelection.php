<html>
<head>
<title>PDF Selection</title>
<link rel="stylesheet" type="text/css" href="pythia.css"/>
<link rel="shortcut icon" href="pythia32.gif"/>
</head>
<body>

<script language=javascript type=text/javascript>
function stopRKey(evt) {
var evt = (evt) ? evt : ((event) ? event : null);
var node = (evt.target) ? evt.target :((evt.srcElement) ? evt.srcElement : null);
if ((evt.keyCode == 13) && (node.type=="text"))
{return false;}
}

document.onkeypress = stopRKey;
</script>
<?php
if($_POST['saved'] == 1) {
if($_POST['filepath'] != "files/") {
echo "<font color='red'>SETTINGS SAVED TO FILE</font><br/><br/>"; }
else {
echo "<font color='red'>NO FILE SELECTED YET.. PLEASE DO SO </font><a href='SaveSettings.php'>HERE</a><br/><br/>"; }
}
?>

<form method='post' action='PDFSelection.php'>

<h2>PDF Selection</h2>

This page contains five subsections. The first deals with how to 
pick  the parton distribution set for protons, including from LHAPDF, 
to be used for all proton and antiproton beams. The second is a special
option that allows a separate PDF set to be used for the hard process
only, while the first choice would still apply to everything else.
The third and fourth give access to pion and Pomeron PDF's, respectively,
the latter being used to describe diffractive systems.
The fifth gives the possibility to switch off the lepton 
"parton density". More information on PDF classes is found 
<?php $filepath = $_GET["filepath"];
echo "<a href='PartonDistributions.php?filepath=".$filepath."' target='page'>";?>here</a>.

<h3>Parton densities for protons</h3>

The selection of parton densities is made once and then is propagated 
through the program. It is essential to make an informed choice, 
for several reasons [<a href="Bibliography.php" target="page">Kas10</a>]: 
<br/><b>Warning 1:</b> the choice of PDF set affects a number of
properties of events. A change of PDF therefore requires a complete 
retuning e.g.  of the multiparton-interactions model for minimum-bias and 
underlying events.
<br/><b>Warning 2:</b> People often underestimate the differences 
between different sets on the market. The sets for the same order are 
constructed to behave more or less similarly at large <i>x</i> and 
<i>Q^2</i>, while the multiparton interactions are dominated by the 
behaviour in the region of small <i>x</i> and <i>Q^2</i>. A good 
PDF parametrization ought to be sensible down to <i>x = 10^-6</i> 
(<i>x = 10^-7</i>) and <i>Q^2 = 1</i> GeV^2 for Tevatron (LHC) 
applications. Unfortunately there are distributions on the market that 
completely derail in that region. The <code>main51.cc</code> and 
<code>main52.cc</code> programs in the <code>examples</code> 
subdirectory provide some examples of absolutely minimal sanity checks 
before a new PDF set is put in production.
<br/><b>Warning 3:</b> NLO and LO sets tend to have quite different
behaviours, e.g. NLO ones have less gluons at small x, which then is 
compensated by positive corrections in the NLO matrix elements.
Therefore do not blindly assume that an NLO tune has to be better than 
an LO one when combined with the LO matrix elements in PYTHIA. There are 
explicit examples where such thinking can lead you down the wrong alley,
especially if you study low-<i>pT</i> physics. In the list below you 
should therefore be extra cautious when using set 6 or set 9.

<p/>
The simplest option is to pick one 
of the distributions available internally:

<br/><br/><table><tr><td><strong>PDF:pSet  </td><td>  &nbsp;&nbsp;(<code>default = <strong>2</strong></code>; <code>minimum = 1</code>; <code>maximum = 12</code>)</td></tr></table>
Parton densities to be used for proton beams (and, by implication,
antiproton ones):
<br/>
<input type="radio" name="1" value="1"><strong>1 </strong>: GRV 94L, LO <ei>alpha_s(M_Z) = 0.128</ei> (this set is out of date, but retained for historical comparisons).<br/>
<input type="radio" name="1" value="2" checked="checked"><strong>2 </strong>: CTEQ 5L, LO <ei>alpha_s(M_Z) = 0.127</ei> (this set is also out of date, but not badly so, and many tunes  are based on it).<br/>
<input type="radio" name="1" value="3"><strong>3 </strong>: MRST LO* (2007),  NLO <ei>alpha_s(M_Z) = 0.12032</ei>.<br/>
<input type="radio" name="1" value="4"><strong>4 </strong>: MRST LO** (2008),  NLO <ei>alpha_s(M_Z) = 0.11517</ei>.<br/>
<input type="radio" name="1" value="5"><strong>5 </strong>: MSTW 2008 LO (central member),  LO <ei>alpha_s(M_Z) = 0.13939</ei>.<br/>
<input type="radio" name="1" value="6"><strong>6 </strong>: MSTW 2008 NLO (central member),  NLO <ei>alpha_s(M_Z) = 0.12018</ei> (NLO, see Warning 3 above).<br/>
<input type="radio" name="1" value="7"><strong>7 </strong>: CTEQ6L, NLO <ei>alpha_s(M_Z) = 0.1180</ei>.<br/>
<input type="radio" name="1" value="8"><strong>8 </strong>: CTEQ6L1, LO <ei>alpha_s(M_Z) = 0.1298</ei>.<br/>
<input type="radio" name="1" value="9"><strong>9 </strong>: CTEQ66.00 (NLO, central member),  NLO <ei>alpha_s(M_Z) = 0.1180</ei> (NLO, see Warning 3 above).<br/>
<input type="radio" name="1" value="10"><strong>10 </strong>: CT09MC1, LO <ei>alpha_s(M_Z) = 0.1300</ei>.<br/>
<input type="radio" name="1" value="11"><strong>11 </strong>: CT09MC2, NLO <ei>alpha_s(M_Z) = 0.1180</ei>.<br/>
<input type="radio" name="1" value="12"><strong>12 </strong>: CT09MCS, NLO <ei>alpha_s(M_Z) = 0.1180</ei>.<br/>
<br/><b>Note:</b> the <i>alpha_s(M_Z)</i> values and the order of the
running in the description above is purely informative, and does not 
affect any other parts of the program. Instead you have the freedom to
set <i>alpha_s(M_Z)</i> value and running separately for 
<?php $filepath = $_GET["filepath"];
echo "<a href='CouplingsAndScales.php?filepath=".$filepath."' target='page'>";?>hard processes</a> 
(including resonance decays),
<?php $filepath = $_GET["filepath"];
echo "<a href='MultipartonInteractions.php?filepath=".$filepath."' target='page'>";?>multiparton interactions</a>,
<?php $filepath = $_GET["filepath"];
echo "<a href='SpacelikeShowers.php?filepath=".$filepath."' target='page'>";?>initial-state radiation</a>, and
<?php $filepath = $_GET["filepath"];
echo "<a href='TimelikeShowers.php?filepath=".$filepath."' target='page'>";?>final-state radiation</a>.

<p/>
This is a reasonably complete list of recent LO fits, both
ones within the normal LO context and ones with modifications for better
matching to event generators. In addition two older sets are 
included for backwards reference (most studies to date are based on 
CTEQ 5L). If you link to the 
<a href="http://projects.hepforge.org/lhapdf/" target="page">LHAPDF 
library</a> [<a href="Bibliography.php" target="page">Wha05</a>] you get access to a much wider selection.
<br/><b>Warning 1:</b> owing to previous problems with the behaviour 
of PDF's beyond the <i>x</i> and <i>Q^2</i> boundaries of a set, 
you should only use LHAPDF <b>version 5.3.0 or later</b>.
<br/><b>Warning 2:</b> the behaviour of the LHAPDF sets need not be
identical with the implementation found in PYTHIA. Specifically we
are aware of the following points that may influence a comparison.
<br/>(a) CTEQ 5L in PYTHIA is the parametrization, in LHAPDF the grid
interpolation. 
<br/>(b) MRST LO* and LO** in PYTHIA is based on an updated edition,
where one makes use of the expanded MSTW grid format, while LHAPDF
is based on the original smaller grid. 
<br/>(c) The CTEQ 6 and CT09MC sets in PYTHIA are frozen at the 
boundaries of the grid, by recommendation of the authors, while 
LHAPDF also offers an option with a smooth extrapolation outside 
the grid boundaries. 

<br/><br/><strong>PDF:useLHAPDF</strong>  <input type="radio" name="2" value="on"><strong>On</strong>
<input type="radio" name="2" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
If off then the choice of proton PDF is based on <code>PDF:pSet</code>
above. If on then it is instead based on the choice of 
<code>PDF:LHAPDFset</code> and <code>PDF:LHAPDFmember</code> below.
<br/><b>Note:</b> in order for this option to work you must have 
compiled PYTHIA appropriately and have set the <code>LHAPATH</code> 
environment variable to provide the data-files directory of your local 
LHAPDF installation. See the README file in the <code>examples</code> 
directory for further instructions. 
  

<br/><br/><table><tr><td><strong>PDF:LHAPDFset  </td><td></td><td> <input type="text" name="3" value="MRST2004FF4lo.LHgrid" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>MRST2004FF4lo.LHgrid</strong></code>)</td></tr></table>
Name of proton PDF set from LHAPDF to be used. You have to choose 
from the 
<a href="http://projects.hepforge.org/lhapdf/pdfsets" target="page">
list of available sets</a>. Examples of some fairly recent ones 
(but still less recent than found above) would be 
cteq61.LHpdf, cteq61.LHgrid, cteq6l.LHpdf, cteq6ll.LHpdf, 
MRST2004nlo.LHpdf, MRST2004nlo.LHgrid, MRST2004nnlo.LHgrid and 
MRST2004FF3lo.LHgrid. If you pick a LHpdf set it will require some 
calculation the first time it is called. 
<br/><b>Technical note:</b> if you provide a name beginning with a 
slash (/) it is assumed you want to provide the full file path and then
<code>initPDFsetM(name)</code> is called, else the correct path is assumed 
already set and <code>initPDFsetByNameM(name)</code> is called.
   

<br/><br/><table><tr><td><strong>PDF:LHAPDFmember  </td><td></td><td> <input type="text" name="4" value="0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>)</td></tr></table>
Further choice of a specific member from the set picked above. Member 0
should normally correspond to the central value, with higher values
corresponding to different error PDF's somewhat off in different 
directions. You have to check from set to set which options are open.
<br/><b>Note:</b> you can only use one member in a run, so if you
want to sweep over many members you either have to do many separate
runs or, as a simplification, save the 
<?php $filepath = $_GET["filepath"];
echo "<a href='EventInformation.php?filepath=".$filepath."' target='page'>";?>pdf weights</a> at the hard scattering
and do an offline reweighting of events.
     

<br/><br/><strong>PDF:extrapolateLHAPDF</strong>  <input type="radio" name="5" value="on"><strong>On</strong>
<input type="radio" name="5" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Parton densities have a guaranteed range of validity in <i>x</i>
and <i>Q^2</i>, and what should be done beyond that range usually is 
not explained by the authors of PDF sets. Nevertheless these boundaries
very often are exceeded, e.g. minimum-bias studies at LHC may sample
<i>x</i> values down to <i>10^-8</i>, while many PDF sets stop
already at <i>10^-5</i>. The default behaviour is then that the 
PDF's are frozen at the boundary, i.e. <i>xf(x,Q^2)</i> is fixed at
its value at <i>x_min</i> for all values <i>x &lt; x_min</i>,
and so on. This is a conservative approach. Alternatively, if you
switch on extrapolation, then parametrizations will be extended beyond
the boundaries, by some prescription. In some cases this will provide a
more realistic answer, in others complete rubbish. Another problem is 
that some of the PDF-set codes will write a warning message anytime the
limits are exceeded, thus swamping your output file. Therefore you should 
study a set seriously before you run it with this switch on.
  

<p/> 
If you want to use PDF's not found in LHAPDF, or you want to interface
LHAPDF another way, you have full freedom to use the more generic 
<?php $filepath = $_GET["filepath"];
echo "<a href='PartonDistributions.php?filepath=".$filepath."' target='page'>";?>interface options</a>.

<h3>Parton densities for protons in the hard process</h3>

The above options provides a PDF set that will be used everywhere:
for the hard process, the parton showers and the multiparton interactions
alike. As already mentioned, therefore a change of PDF should be
accompanied by a <b>complete</b> retuning of the whole MPI framework,
and maybe more. There are cases where one may want to explore 
different PDF options for the hard process, but would not want to touch 
the rest. If several different sets are to be compared, a simple
reweighting based on the <?php $filepath = $_GET["filepath"];
echo "<a href='EventInformation.php?filepath=".$filepath."' target='page'>";?>originally 
used</a> flavour, <i>x</i>, <i>Q^2</i> and PDF values may offer the 
best route. The options in this section allow a choice of the PDF set
for the hard process alone, while the choice made in the previous section
would still be used for everything else. The hardest interaction
of the minimum-bias process is part of the multiparton-interactions
framework and so does not count as a hard process here. 

<p/>
Of course it is inconsistent to use different PDF's in different parts 
of an event, but if the <i>x</i> and <i>Q^2</i> ranges mainly accessed 
by the components are rather different then the contradiction would not be
too glaring. Furthermore, since standard PDF's are one-particle-inclusive
we anyway have to 'invent' our own PDF modifications to handle configurations
where more than one parton is kicked out of the proton [<a href="Bibliography.php" target="page">Sjo04</a>]. 

<p/>
The PDF choices that can be made are the same as above, so we do not 
repeat the detailed discussion.

<br/><br/><strong>PDF:useHard</strong>  <input type="radio" name="6" value="on"><strong>On</strong>
<input type="radio" name="6" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
If on then select a separate PDF set for the hard process, using the 
variables below. If off then use the same PDF set for everything,
as already chosen above.   
  

<br/><br/><table><tr><td><strong>PDF:pHardSet  </td><td>  &nbsp;&nbsp;(<code>default = <strong>2</strong></code>; <code>minimum = 1</code>; <code>maximum = 12</code>)</td></tr></table>
Parton densities to be used for proton beams (and, by implication,
antiproton ones):
<br/>
<input type="radio" name="7" value="1"><strong>1 </strong>: GRV 94L, LO <ei>alpha_s(M_Z) = 0.128</ei> (out of date).<br/>
<input type="radio" name="7" value="2" checked="checked"><strong>2 </strong>: CTEQ 5L, LO <ei>alpha_s(M_Z) = 0.127</ei> (slightly out of date; many tunes are based on it).<br/>
<input type="radio" name="7" value="3"><strong>3 </strong>: MRST LO* (2007),  NLO <ei>alpha_s(M_Z) = 0.12032</ei>.<br/>
<input type="radio" name="7" value="4"><strong>4 </strong>: MRST LO** (2008),  NLO <ei>alpha_s(M_Z) = 0.11517</ei>.<br/>
<input type="radio" name="7" value="5"><strong>5 </strong>: MSTW 2008 LO (central member),  LO <ei>alpha_s(M_Z) = 0.13939</ei>.<br/>
<input type="radio" name="7" value="6"><strong>6 </strong>: MSTW 2008 NLO (central member),  LO <ei>alpha_s(M_Z) = 0.12018</ei>.<br/>
<input type="radio" name="7" value="7"><strong>7 </strong>: CTEQ6L, NLO <ei>alpha_s(M_Z) = 0.1180</ei>.<br/>
<input type="radio" name="7" value="8"><strong>8 </strong>: CTEQ6L1, LO <ei>alpha_s(M_Z) = 0.1298</ei>.<br/>
<input type="radio" name="7" value="9"><strong>9 </strong>: CTEQ66.00 (NLO, central member),  NLO <ei>alpha_s(M_Z) = 0.1180</ei>.<br/>
<input type="radio" name="7" value="10"><strong>10 </strong>: CT09MC1, LO <ei>alpha_s(M_Z) = 0.1300</ei>.<br/>
<input type="radio" name="7" value="11"><strong>11 </strong>: CT09MC2, NLO <ei>alpha_s(M_Z) = 0.1180</ei>.<br/>
<input type="radio" name="7" value="12"><strong>12 </strong>: CT09MCS, NLO <ei>alpha_s(M_Z) = 0.1180</ei>.<br/>

<br/><br/><strong>PDF:useHardLHAPDF</strong>  <input type="radio" name="8" value="on"><strong>On</strong>
<input type="radio" name="8" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
If off then the choice of proton PDF is based on <code>hardpPDFset</code>
above. If on then it is instead based on the choice of 
<code>hardLHAPDFset</code> and <code>hardLHAPDFmember</code> below.
Note that if you want to use LHAPDF here, and you also use LHAPDF
for the "normal" PDF set, then LHAPDF must have been compiled so as to 
handle (at least) two concurrent sets, with the configure statement
<code>--with-max-num-pdfsets=2</code>. 
  

<br/><br/><table><tr><td><strong>PDF:hardLHAPDFset  </td><td></td><td> <input type="text" name="9" value="MRST2004FF4lo.LHgrid" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>MRST2004FF4lo.LHgrid</strong></code>)</td></tr></table>
Name of proton PDF set from LHAPDF to be used. 
   

<br/><br/><table><tr><td><strong>PDF:hardLHAPDFmember  </td><td></td><td> <input type="text" name="10" value="0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>)</td></tr></table>
Further choice of a specific member from the set picked above. 
     

<p/>
Note that there is no separate equivalent of the 
<code>PDF:extrapolateLHAPDF</code> flag specifically for the hard
PDF. Since LHAPDF only has one global flag for extrapolation or not,
the choice for the normal PDF's also applies to the hard ones.

<h3>Parton densities for pions</h3>

The parton densities of the pion are considerably less well known than
those of the proton. There are only rather few sets on the market,
and none particularly recent. Only one comes built-in, but others can 
be accessed from LHAPDF. Input parametrizations are for the <i>pi+</i>.
>From this the <i>pi-</i> is obtained by charge conjugation and the 
<i>pi0</i> from averaging (half the pions have <i>d dbar</i> 
valence quark content, half <i>u ubar</i>.

<p/>
Much of the switches are taken over from the proton case, with obvious
modifications; therefore the description is briefer. Currently we have 
not seen the need to allow separate parton densities for hard processes. 
When using LHAPDF the <code>PDF:extrapolateLHAPDF</code> switch of the 
proton also applies to pions. 
 
<br/><br/><table><tr><td><strong>PDF:piSet  </td><td>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = 1</code>; <code>maximum = 1</code>)</td></tr></table>
Internal parton densities that can be used for pion beams, currently with 
only one choice.
<br/>
<input type="radio" name="11" value="1" checked="checked"><strong>1 </strong>: GRV 92 L.<br/>

<br/><br/><strong>PDF:piUseLHAPDF</strong>  <input type="radio" name="12" value="on"><strong>On</strong>
<input type="radio" name="12" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
If off then the choice of proton PDF is based on <code>PDF:piSet</code>
above. If on then it is instead based on the choice of 
<code>PDF:piLHAPDFset</code> and <code>PDF:piLHAPDFmember</code> below.
  

<br/><br/><table><tr><td><strong>PDF:piLHAPDFset  </td><td></td><td> <input type="text" name="13" value="OWPI.LHgrid" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>OWPI.LHgrid</strong></code>)</td></tr></table>
Name of pion PDF set from LHAPDF to be used. You have to choose from the 
<a href="http://projects.hepforge.org/lhapdf/pdfsets" target="page">
list of available sets</a>. 
   

<br/><br/><table><tr><td><strong>PDF:piLHAPDFmember  </td><td></td><td> <input type="text" name="14" value="0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>)</td></tr></table>
Further choice of a specific member from the set picked above.
     

<h3>Parton densities for Pomerons</h3>

The Pomeron is introduced in the description of diffractive events, 
i.e. a diffractive system is viewed as a Pomeron-proton collision at a 
reduced CM energy. Here the PDF's are even less well known. 
Most experimental parametrizations are NLO, which makes them less
well suited for Monte Carlo applications. Furthermore note that 
the momentum sum is arbitrarily normalized to a non-unity value.

<br/><br/><table><tr><td><strong>PDF:PomSet  </td><td>  &nbsp;&nbsp;(<code>default = <strong>6</strong></code>; <code>minimum = 1</code>; <code>maximum = 6</code>)</td></tr></table>
Parton densities that can be used for Pomeron beams. 
<br/>
<input type="radio" name="15" value="1"><strong>1 </strong>: <ei>Q^2</ei>-independent parametrizations <ei>xf(x) = N_ab x^a (1 - x)^b</ei>, where <ei>N_ab</ei> ensures unit momentum sum. The <ei>a</ei> and <ei>b</ei> parameters can be  set separately for the gluon and the quark distributions. The momentum fraction of gluons and quarks can be freely mixed, and  production of <ei>s</ei> quarks can be suppressed relative to  that of <ei>d</ei> and <ei>u</ei> ones, with antiquarks as likely  as quarks. See further below how to set the six parameters of this  approach. <br/>
<input type="radio" name="15" value="2"><strong>2 </strong>: <ei>pi0</ei> distributions, as specified in the  section above. <br/>
<input type="radio" name="15" value="3"><strong>3 </strong>: the H1 2006 Fit A NLO <ei>Q^2</ei>-dependent  parametrization, based on a tune to their data <ref>H1P06</ref>, rescaled by the factor <code>PomRescale</code> below. <br/>
<input type="radio" name="15" value="4"><strong>4 </strong>: the H1 2006 Fit B NLO <ei>Q^2</ei>-dependent  parametrization, based on a tune to their data <ref>H1P06</ref>, rescaled by the factor <code>PomRescale</code> below. <br/>
<input type="radio" name="15" value="5"><strong>5 </strong>: the H1 2007 Jets NLO <ei>Q^2</ei>-dependent  parametrization, based on a tune to their data <ref>H1P07</ref>, rescaled by the factor <code>PomRescale</code> below. <br/>
<input type="radio" name="15" value="6" checked="checked"><strong>6 </strong>: the H1 2006 Fit B LO <ei>Q^2</ei>-dependent  parametrization, based on a tune to their data <ref>H1P06</ref>, rescaled by the factor <code>PomRescale</code> below. <br/>

<br/><br/><table><tr><td><strong>PDF:PomGluonA </td><td></td><td> <input type="text" name="16" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>; <code>minimum = -0.5</code>; <code>maximum = 2.</code>)</td></tr></table>
the parameter <i>a</i> in the ansatz <i>xg(x) = N_ab x^a (1 - x)^b</i>
for option 1 above.
  

<br/><br/><table><tr><td><strong>PDF:PomGluonB </td><td></td><td> <input type="text" name="17" value="3." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>3.</strong></code>; <code>minimum = 0.</code>; <code>maximum = 10.</code>)</td></tr></table>
the parameter <i>b</i> in the ansatz <i>xg(x) = N_ab x^a (1 - x)^b</i>
for option 1 above.
  

<br/><br/><table><tr><td><strong>PDF:PomQuarkA </td><td></td><td> <input type="text" name="18" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>; <code>minimum = -0.5</code>; <code>maximum = 2.</code>)</td></tr></table>
the parameter <i>a</i> in the ansatz <i>xq(x) = N_ab x^a (1 - x)^b</i>
for option 1 above.
  

<br/><br/><table><tr><td><strong>PDF:PomQuarkB </td><td></td><td> <input type="text" name="19" value="3." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>3.</strong></code>; <code>minimum = 0.</code>; <code>maximum = 10.</code>)</td></tr></table>
the parameter <i>b</i> in the ansatz <i>xq(x) = N_ab x^a (1 - x)^b</i>
for option 1 above.
  

<br/><br/><table><tr><td><strong>PDF:PomQuarkFrac </td><td></td><td> <input type="text" name="20" value="0.2" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.2</strong></code>; <code>minimum = 0.</code>; <code>maximum = 1.</code>)</td></tr></table>
the fraction of the Pomeron momentum carried by quarks 
for option 1 above, with the rest carried by gluons.
  

<br/><br/><table><tr><td><strong>PDF:PomStrangeSupp </td><td></td><td> <input type="text" name="21" value="0.5" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.5</strong></code>; <code>minimum = 0.</code>; <code>maximum = 1.</code>)</td></tr></table>
the suppression of the <i>s</i> quark density relative to that of the 
<i>d</i> and <i>u</i> ones for option 1 above.
  

<br/><br/><table><tr><td><strong>PDF:PomRescale </td><td></td><td> <input type="text" name="22" value="1.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.0</strong></code>; <code>minimum = 0.5</code>; <code>maximum = 5.0</code>)</td></tr></table>
Rescale the four H1 fits above by this uniform factor, e.g. to bring 
up their momentum sum to around unity. By default all three have
a momentum sum of order 0.5, suggesting that a factor around 2.0
should be used. You can use <code>examples/main51.cc</code> to get
a more precise value. Note that also other parameters in the 
<?php $filepath = $_GET["filepath"];
echo "<a href='Diffraction.php?filepath=".$filepath."' target='page'>";?>diffraction</a> framework may need to
be retuned when this parameter is changed.
  

<h3>Parton densities for leptons</h3>

For electrons/muons/taus there is no need to choose between different 
parametrizations, since only one implementation is available, and 
should be rather uncontroversial (apart from some technical details).
However, insofar as e.g. <i>e^+ e^-</i> data often are corrected 
back to a world without any initial-state photon radiation, it is 
useful to have a corresponding option available here.

<br/><br/><strong>PDF:lepton</strong>  <input type="radio" name="23" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="23" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Use parton densities for lepton beams or not. If off the colliding
leptons carry the full beam energy, if on part of the energy is 
radiated away by initial-state photons. In the latter case the
initial-state showers will generate the angles and energies of the
set of photons that go with the collision. In addition one collinear
photon per beam carries any leftover amount of energy not described
by shower emissions. If the initial-state showers are switched off 
these collinear photons will carry the full radiated energy.  
   

<p/>
Neutrinos are always taken pointlike. Do note that the phase space 
selection machinery currently does not allow one resolved and one 
unresolved beam. For lepton-neutrino collisions to work you must 
therefore set <code>PDF:lepton = off</code>.

<h3>Incoming parton selection</h3>

There is one useful degree of freedom to restrict the set of incoming 
quark flavours for hard processes. It does not change the PDF's as such, 
only which quarks are allowed to contribute to the hard-process cross 
sections. Note that separate but similarly named modes are available 
for multiparton interactions and spacelike showers.

<br/><br/><table><tr><td><strong>PDFinProcess:nQuarkIn  </td><td></td><td> <input type="text" name="24" value="5" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>5</strong></code>; <code>minimum = 0</code>; <code>maximum = 5</code>)</td></tr></table>
Number of allowed incoming quark flavours in the beams; a change 
to 4 would thus exclude <i>b</i> and <i>bbar</i> as incoming 
partons, etc.
  

<input type="hidden" name="saved" value="1"/>

<?php
echo "<input type='hidden' name='filepath' value='".$_GET["filepath"]."'/>"?>

<table width="100%"><tr><td align="right"><input type="submit" value="Save Settings" /></td></tr></table>
</form>

<?php

if($_POST["saved"] == 1)
{
$filepath = $_POST["filepath"];
$handle = fopen($filepath, 'a');

if($_POST["1"] != "2")
{
$data = "PDF:pSet = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "off")
{
$data = "PDF:useLHAPDF = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "MRST2004FF4lo.LHgrid")
{
$data = "PDF:LHAPDFset = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "0")
{
$data = "PDF:LHAPDFmember = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "off")
{
$data = "PDF:extrapolateLHAPDF = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "off")
{
$data = "PDF:useHard = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
if($_POST["7"] != "2")
{
$data = "PDF:pHardSet = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
if($_POST["8"] != "off")
{
$data = "PDF:useHardLHAPDF = ".$_POST["8"]."\n";
fwrite($handle,$data);
}
if($_POST["9"] != "MRST2004FF4lo.LHgrid")
{
$data = "PDF:hardLHAPDFset = ".$_POST["9"]."\n";
fwrite($handle,$data);
}
if($_POST["10"] != "0")
{
$data = "PDF:hardLHAPDFmember = ".$_POST["10"]."\n";
fwrite($handle,$data);
}
if($_POST["11"] != "1")
{
$data = "PDF:piSet = ".$_POST["11"]."\n";
fwrite($handle,$data);
}
if($_POST["12"] != "off")
{
$data = "PDF:piUseLHAPDF = ".$_POST["12"]."\n";
fwrite($handle,$data);
}
if($_POST["13"] != "OWPI.LHgrid")
{
$data = "PDF:piLHAPDFset = ".$_POST["13"]."\n";
fwrite($handle,$data);
}
if($_POST["14"] != "0")
{
$data = "PDF:piLHAPDFmember = ".$_POST["14"]."\n";
fwrite($handle,$data);
}
if($_POST["15"] != "6")
{
$data = "PDF:PomSet = ".$_POST["15"]."\n";
fwrite($handle,$data);
}
if($_POST["16"] != "0.")
{
$data = "PDF:PomGluonA = ".$_POST["16"]."\n";
fwrite($handle,$data);
}
if($_POST["17"] != "3.")
{
$data = "PDF:PomGluonB = ".$_POST["17"]."\n";
fwrite($handle,$data);
}
if($_POST["18"] != "0.")
{
$data = "PDF:PomQuarkA = ".$_POST["18"]."\n";
fwrite($handle,$data);
}
if($_POST["19"] != "3.")
{
$data = "PDF:PomQuarkB = ".$_POST["19"]."\n";
fwrite($handle,$data);
}
if($_POST["20"] != "0.2")
{
$data = "PDF:PomQuarkFrac = ".$_POST["20"]."\n";
fwrite($handle,$data);
}
if($_POST["21"] != "0.5")
{
$data = "PDF:PomStrangeSupp = ".$_POST["21"]."\n";
fwrite($handle,$data);
}
if($_POST["22"] != "1.0")
{
$data = "PDF:PomRescale = ".$_POST["22"]."\n";
fwrite($handle,$data);
}
if($_POST["23"] != "on")
{
$data = "PDF:lepton = ".$_POST["23"]."\n";
fwrite($handle,$data);
}
if($_POST["24"] != "5")
{
$data = "PDFinProcess:nQuarkIn = ".$_POST["24"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>

<!-- Copyright (C) 2013 Torbjorn Sjostrand -->
