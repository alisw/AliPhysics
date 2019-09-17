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
<ol id="toc">
  <li><a href="#section0">Parton densities for protons</a></li>
  <li><a href="#section1">Parton densities for protons in the hard process</a></li>
  <li><a href="#section2">Nuclear modifications of parton densities</a></li>
  <li><a href="#section3">Parton densities for pions</a></li>
  <li><a href="#section4">Parton densities for Pomerons</a></li>
  <li><a href="#section5">Parton densities for photons</a></li>
  <li><a href="#section6">Parton densities for leptons</a></li>
  <li><a href="#section7">Incoming parton selection</a></li>
</ol>

 
This page contains several subsections. The first deals with how to 
pick the parton distribution set for protons, including from LHAPDF, 
to be used for all proton and antiproton beams. The second is a special 
option that allows a separate PDF set to be used for the hard process 
only, while the first choice would still apply to everything else. 
The third introduces the possibility of nuclear modifications. 
Further sections give access to pion, Pomeron and photon PDF's, 
respectively, the second being used to describe diffractive systems. 
Towards the end comes the possibility to switch off the lepton 
"parton density", and photons from lepton beams. More information 
on PDF classes is found <?php $filepath = $_GET["filepath"];
echo "<a href='PartonDistributions.php?filepath=".$filepath."' target='page'>";?>here</a>. 
 
<br/><br/><strong>PDF:extrapolate</strong>  <input type="radio" name="1" value="on"><strong>On</strong>
<input type="radio" name="1" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Allow PDF sets to be extrapolated to small <i>x</i> values, instead of 
being frozen at <i>x_min</i>. This is a global flag that affects all 
PDF sets used, whenever extrapolation has been implemented. Among 
internal PDFs, all Pomeron sets are affected by this flag, as are 
the CTEQ6/CT09 proton ones, the NNPDF 3.1 ones and others accessed 
by the <code>LHAGrid1</code> approach. For the rest some by default 
extrapolate to small <i>x</i> (GRV 94 L, MRST/MSTW) while others are 
frozen at the border (CTEQ 5 L, NNPDF 2.3). When in doubt, check whether 
and how the behaviour depends on the choice made for your region of 
interest. When LHAPDF (5 or 6) is used, the extrapolation switch is set 
according to the choice here, and the behaviour is according to the 
rules of the respective program. 
<br/>To put the issue in context, parton densities have a guaranteed 
range of validity in <i>x</i> and <i>Q^2</i>, and what should be done 
beyond that range usually is not explained by the authors of PDF sets. 
Nevertheless these boundaries very often are exceeded, e.g. minimum-bias 
studies at LHC may sample <i>x</i> values down to <i>10^-8</i>, while 
many PDF sets stop already at <i>10^-5</i>. The default behaviour is 
then that the PDF's are frozen at the boundary, i.e. <i>xf(x,Q^2)</i> is 
fixed at its value at <i>x_min</i> for all values <i>x &lt; 
x_min</i>, and so on. This is a conservative approach. Alternatively, 
if you switch on extrapolation, then parametrizations will be extended 
beyond the boundaries, by some prescription. In some cases this will 
provide a more realistic answer, in others complete rubbish. Another 
problem is that some of the PDF-set codes will write a warning message 
anytime the limits are exceeded, thus swamping your output 
file. Therefore you should study a set seriously before you run it 
with this switch on. 
<br/><b>Warning:</b>It has been found out that the LHAPDF program by 
default uses a damping of PDFs at low <i>Q</i> scales, below 
<i>Q_min</i>, based on an anomalous dimension ansatz. This overlaps 
with the damping imposed in the MPI framework by the <i>p_T0</i> 
parameter, and to have both would probably imply doublecounting of 
effects. Therefore, as of version 8.227, PDFs are frozen below 
<i>Q_min</i>. This change affects the LHAPDF 5 interface. The native 
LHAPDF 6 interface already contained this restriction, as does the PDFs 
that come with PYTHIA. Also limits at <i>Q_max</i> and <i>x_max</i> 
are checked and PDFs frozen outside them, so the extrapolate option now 
is strictly a choice of low-<i>x</i> behaviour. 
   
 
<a name="section0"></a> 
<h3>Parton densities for protons</h3> 
 
PYTHIA comes with a reasonably complete list of recent LO fits built-in, 
both ones within the normal LO context and ones with modifications for 
better matching to event generators. In addition two older sets are 
included for backwards reference (most studies to date are based on 
CTEQ 5L). Therefore there is no real need to link any external PDF sets. 
 
<p/> 
If the internal PDF sets are not sufficient, the 
<a href="http://projects.hepforge.org/lhapdf/" target="page">LHAPDF 
library</a> [<a href="Bibliography.php#refWha05" target="page">Wha05</a>,<a href="Bibliography.php#refBuc15" target="page">Buc15</a>] gives you access to a much wider 
selection. 
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
 
<p/> 
If you do not want to install LHAPDF, it is possible to use LHAPDF6 
data grids natively in PYTHIA. This is based on a simplified 
implementation of interpolation in a <code>.dat</code> "lhagrid1" 
file, and so does not give fully identical results, and also is not 
foolproof. 
 
<p/> 
The selection of parton densities is made once and then is propagated 
through the program. It is essential to make an informed choice, 
for several reasons [<a href="Bibliography.php#refKas10" target="page">Kas10</a>]: 
<br/><b>Warning 1:</b> the choice of PDF set affects a number of 
properties of events. A change of PDF therefore requires a complete 
retuning e.g.  of the multiparton-interactions model for minimum-bias and 
underlying events. Conversely, the 
<?php $filepath = $_GET["filepath"];
echo "<a href='Tunes.php?filepath=".$filepath."' target='page'>";?>pp physics tunes</a> are all made for a specific 
PDF tune, and the chosen (or default) tune will therefore overwrite 
the <code>PDF:pSet</code> default value described below. If you want 
to set <code>PDF:pSet</code> differently it should be done <i>after</i> 
the <code>Tune:pp</code> value, if any, has been set. 
<br/><b>Warning 2:</b> People often underestimate the differences 
between different sets on the market. The sets for the same order are 
constructed to behave  or less similarly at large <i>x</i> and 
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
especially if you study low-<i>pT</i> physics. A longer discussion on 
this point can be found in <a href="../pdfdoc/pdfwarning.pdf">this note</a>. 
In the list below you should therefore be extra cautious when using 
set 6 or set 9. 
 
<br/><br/><table><tr><td><strong>PDF:pSet  </td><td></td><td> <input type="text" name="2" value="13" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>13</strong></code>)</td></tr></table>
Parton densities to be used for proton beams (and, by implication, 
antiproton ones). Note that the choice of a string input (rather than 
e.g. an integer) allows to pick either an internal, LHAPDF5 or LHAPDF6 
set in one single setting, by some behind-the-scenes machinations. 
<br/><code>option </code><strong> 1</strong> : GRV 94L, LO <i>alpha_s(M_Z) = 0.128</i> 
(this set is out of date, but retained for historical comparisons).   
<br/><code>option </code><strong> 2</strong> : CTEQ 5L, LO <i>alpha_s(M_Z) = 0.127</i> 
(this set is also out of date, but not badly so, and many tunes 
are based on it).   
<br/><code>option </code><strong> 3</strong> : MRST LO* (2007), 
NLO <i>alpha_s(M_Z) = 0.12032</i>.   
<br/><code>option </code><strong> 4</strong> : MRST LO** (2008), 
NLO <i>alpha_s(M_Z) = 0.11517</i>.   
<br/><code>option </code><strong> 5</strong> : MSTW 2008 LO (central member), 
LO <i>alpha_s(M_Z) = 0.13939</i>.   
<br/><code>option </code><strong> 6</strong> : MSTW 2008 NLO (central member), 
NLO <i>alpha_s(M_Z) = 0.12018</i> (NLO, see Warning 3 above).   
<br/><code>option </code><strong> 7</strong> : CTEQ6L, NLO <i>alpha_s(M_Z) = 0.1180</i>.   
<br/><code>option </code><strong> 8</strong> : CTEQ6L1, LO <i>alpha_s(M_Z) = 0.1298</i>.   
<br/><code>option </code><strong> 9</strong> : CTEQ66.00 (NLO, central member), 
NLO <i>alpha_s(M_Z) = 0.1180</i> (NLO, see Warning 3 above).   
<br/><code>option </code><strong> 10</strong> : CT09MC1, LO <i>alpha_s(M_Z) = 0.1300</i>.   
<br/><code>option </code><strong> 11</strong> : CT09MC2, NLO <i>alpha_s(M_Z) = 0.1180</i>.   
<br/><code>option </code><strong> 12</strong> : CT09MCS, NLO <i>alpha_s(M_Z) = 0.1180</i>.   
<br/><code>option </code><strong> 13</strong> : NNPDF2.3 QCD+QED LO <i>alpha_s(M_Z) = 0.130</i>. 
   
<br/><code>option </code><strong> 14</strong> : NNPDF2.3 QCD+QED LO <i>alpha_s(M_Z) = 0.119</i>. 
   
<br/><code>option </code><strong> 15</strong> : NNPDF2.3 QCD+QED NLO <i>alpha_s(M_Z) = 0.119</i>. 
   
<br/><code>option </code><strong> 16</strong> : NNPDF2.3 QCD+QED NNLO <i>alpha_s(M_Z) = 0.119</i>. 
   
<br/><b>Warning :</b>the following four NNPDF 3.1 sets are quite 
different from the NNPDF 2.3 ones, and cannot be used interchangeably, 
but need retuning of the MPI framework. Some also do not contain QED 
evolution. 
<br/><code>option </code><strong> 17</strong> : NNPDF3.1 QCD LO <i>alpha_s(M_Z) = 0.130</i>. 
   
<br/><code>option </code><strong> 18</strong> : NNPDF3.1 QCD LO <i>alpha_s(M_Z) = 0.118</i>. 
   
<br/><code>option </code><strong> 19</strong> : NNPDF3.1 QCD+LUXQED NLO <i>alpha_s(M_Z) = 0.118</i>. 
   
<br/><code>option </code><strong> 20</strong> : NNPDF3.1 QCD+LUXQED NNLO <i>alpha_s(M_Z) = 0.118</i>. 
   
<br/><code>option </code><strong> 21</strong> : NNPDF3.1sx+LHCb NLO+NLLx LUXQED 
<i>alpha_s(M_Z) = 0.118</i> [<a href="Bibliography.php#refBer18" target="page">Ber18</a>]. While at NLO, the 
additional small-<i>x</i> resummation, anchored by LHC-b data, offers 
a  reasonable small-<i>x</i> behaviour than most NLO PDFs, as 
required for the successful usage e.g. with traditional "improved LL" 
parton showers. The photon part is unchanged from the earlier NNPDF 3.1 
QED analysis [<a href="Bibliography.php#refBer17" target="page">Ber17</a>]. 
<br/><b>Warning :</b>in version 8.235 the 21 identifier was used to 
denote and earlier attempt to obtain a  reasonable small-<i>x</i> 
behaviour. This PDF set is superseded by the new 21 and 22 sets, and has 
been removed, as was forewarned. 
   
<br/><code>option </code><strong> 22</strong> : NNPDF3.1sx+LHCb NNLO+NLLx LUXQED 
<i>alpha_s(M_Z) = 0.118</i>. Comments as for 21, but this set is at 
NNLO rather than NLO. 
   
<br/><code>option </code><strong> LHAPDF5:set/member</strong> : Use an external LHAPDF set 
where <code>set</code> is the name of the set to use 
and <code>member</code> is the member of the set to use. The value 
for <code>set</code> is the name of the PDF set to use while the value 
for <code>member</code> must be an integer and is the member of the 
set to use. If member is not supplied, then <code>0</code> is assumed. 
   
<br/><code>option </code><strong> LHAPDF6:set/member</strong> : Same as 
for <code>LHAPDF5:set/member</code> but now the LHAPDF6 library is 
used instead. 
   
<br/><code>option </code><strong> LHAGrid1:filename</strong> : Use the internal implementation 
of interpolation in <code>.dat</code> files in the default "lhagrid1" 
LHAPDF6 format. This is a simplified implementation, with cubic 
interpolation in <i>ln(x)</i> and in <i>ln(Q2)</i>. If 
there are several <i>Q^2</i> subgrids they have to have the same 
<i>x</i> grid. (Linear interpolation in <i>ln(Q2)</i> is used, 
should a subgrid contain fewer than four <i>Q2</i> values.) 
Other restrictions may also apply, so use with caution. 
If the <code>filename</code> begins with a / it is supposed to contain 
the absolute path to the file, and if not the file is supposed to be 
located in the standard <code>share/Pythia8/xmldoc</code> directory. 
Note that, unlike LHAPDF, there is no explicit hierarchy of a set 
containing separate members; each <code>.dat</code> file can be used 
without any reference to the set it is a member of. 
   
<br/><b>Warning 1:</b> the <i>alpha_s(M_Z)</i> values and the order of the 
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
<br/><b>Warning 2:</b> in order for <code>LHAPDF</code> PDF sets to work 
you must have compiled the appropriate LHAPDF plugins for PYTHIA and 
have set the <code>LHAPATH</code> environment variable 
(or <code>LHAPDF_DATA_PATH</code>) to provide the data-files directory 
of your local LHAPDF installation. See the README file in 
the <code>examples</code> directory for further instructions. 
<br/><b>Warning 3:</b> it is technically possible to simultaneously 
use <code>LHAPDF5</code> and <code>LHAPDF6</code> PDF sets at the same 
time for the two beams, but such a configuration is not officially 
supported and strongly discouraged. 
   
 
<br/><br/><table><tr><td><strong>PDF:pSetB  </td><td></td><td> <input type="text" name="3" value="void" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>void</strong></code>)</td></tr></table>
Parton densities to be used by <i>proton beam B</i>, with the same 
options available as for <code>PDF:pSet</code>. If this option is set 
to <code>void</code> then the same PDF set as <code>PDF:pSet</code> is 
used. 
   
 
<p/> 
If you want to use PDF's not found in LHAPDF, or you want to interface 
LHAPDF another way, you have full freedom to use the more generic 
<?php $filepath = $_GET["filepath"];
echo "<a href='PartonDistributions.php?filepath=".$filepath."' target='page'>";?>interface options</a>. 
 
<a name="section1"></a> 
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
where  than one parton is kicked out of the proton [<a href="Bibliography.php#refSjo04" target="page">Sjo04</a>]. 
 
<p/> 
The PDF choices that can be made are the same as above, so we do not 
repeat the detailed discussion. 
 
<br/><br/><strong>PDF:useHard</strong>  <input type="radio" name="4" value="on"><strong>On</strong>
<input type="radio" name="4" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
If on then select a separate PDF set for the hard process, using the 
variables below. If off then use the same PDF set for everything, 
as already chosen above. 
   
 
<br/><br/><table><tr><td><strong>PDF:pHardSet  </td><td></td><td> <input type="text" name="5" value="void" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>void</strong></code>)</td></tr></table>
Parton densities to be used by the proton beams of the hard process, 
with the same options available as for <code>PDF:pSet</code>. If this 
option is set to <code>void</code> then the same PDF set 
as <code>PDF:pSet</code> is used. 
   
 
<br/><br/><table><tr><td><strong>PDF:pHardSetB  </td><td></td><td> <input type="text" name="6" value="void" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>void</strong></code>)</td></tr></table>
Parton densities to be used by <i>proton beam B</i> of the hard 
process, with the same options available as 
for <code>PDF:pSet</code>. If this option is set to <code>void</code> 
then the same PDF set as <code>PDF:pHardSet</code> is used. 
   
 
<a name="section2"></a> 
<h3>Nuclear modifications of parton densities</h3> 
 
<p/> 
Nuclear modifications of the PDFs are implemented for the hard-process 
generation only. The final PDF value is calculated for an average nucleon 
within given nucleus, i.e. 
<i>f_i^A(x,Q^2) = (Z/A)*f_i^(p/A) + ((A-Z)/A)*f_i^(n/A)</i>, where 
<i>A</i> is the nuclear mass number and <i>Z</i> the number of protons, 
set using the PDG code for nucleus. The neutron PDFs are obtained by 
applying isospin symmetry, e.g. <i>f_u^(n/A)(x,Q^2) = f_d^(p/A)(x,Q^2)</i>. 
The nuclear PDFs implemented provide only the nuclear modification so the 
full PDF is calculated by multiplying the selected free proton PDF with the 
modification. 
<p/> 
 
<br/><br/><strong>PDF:useHardNPDFA</strong>  <input type="radio" name="7" value="on"><strong>On</strong>
<input type="radio" name="7" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
If on, the hard processes are generated with nuclear modifications for 
beam A. 
   
 
<br/><br/><strong>PDF:useHardNPDFB</strong>  <input type="radio" name="8" value="on"><strong>On</strong>
<input type="radio" name="8" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
If on, the hard processes are generated with nuclear modifications for 
beam B. 
   
 
<br/><br/><table><tr><td><strong>PDF:nPDFSetA  </td><td>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = 0</code>; <code>maximum = 3</code>)</td></tr></table>
The nuclear modication to be used for beam A if enabled with the switch above. 
<br/>
<option value = "0"> Only Isospin effect.</option> 
<option value = "1"> EPS09, LO <ref>Esk09</ref>.</option> 
<option value = "2"> EPS09, NLO <ref>Esk09</ref>. The grid files can be 
found from 
<a href="https://www.jyu.fi/science/en/physics/research/highenergy/urhic/npdfs/eps09"> 
here</a> and are to be stored in the same folder as other PDF grid files 
(usually share/Pythia8/xmldoc/). 
</option> 
<option value = "3"> EPPS16, NLO <ref>Esk16</ref>. The grid files can be 
found from 
<a href="https://www.jyu.fi/science/en/physics/research/highenergy/urhic/npdfs/epps16-nuclear-pdfs"> 
here</a>. 
</option> 
 
<br/><br/><table><tr><td><strong>PDF:nPDFSetB  </td><td>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = 0</code>; <code>maximum = 3</code>)</td></tr></table>
The nuclear modication to be used for beam B. Same options as above. 
<br/>
<option value = "0"> Only Isospin effect.</option> 
<option value = "1"> EPS09, LO. </option> 
<option value = "2"> EPS09, NLO.</option> 
<option value = "3"> EPPS16, NLO.</option> 
 
<br/><br/><table><tr><td><strong>PDF:nPDFBeamA  </td><td></td><td> <input type="text" name="11" value="100822080" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>100822080</strong></code>)</td></tr></table>
The PDG code for nuclear beam A, provides the number of protons and 
neutrons. Default code for Pb. 
   
 
<br/><br/><table><tr><td><strong>PDF:nPDFBeamB  </td><td></td><td> <input type="text" name="12" value="100822080" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>100822080</strong></code>)</td></tr></table>
The PDG code for nucleus B. 
   
 
<a name="section3"></a> 
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
 
<br/><br/><table><tr><td><strong>PDF:piSet  </td><td></td><td> <input type="text" name="13" value="1" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>)</td></tr></table>
Parton densities that can be used for pion beams, currently with 
only one internal choice. 
<br/><code>option </code><strong> 1</strong> : GRV 92 L.   
<br/><code>option </code><strong> LHAPDF5:set/member</strong> : Use an external LHAPDF set 
where <code>set</code> is the name of the set to use 
and <code>member</code> is the member of the set to use. The value 
for <code>set</code> can either be a relative path to the LHAPDF path, 
or an absolute path. The value for <code>member</code> must be an 
integer. 
   
<br/><code>option </code><strong> LHAPDF6:set/member</strong> : Same as 
for <code>LHAPDF5:set/member</code> but now the LHAPDF6 library is 
used instead. 
   
<br/><code>option </code><strong> LHAGrid1:filename</strong> : Use the internal implementation 
of interpolation in <code>.dat</code> files in the default "lhagrid1" 
LHAPDF6 format, cf. the corresponding proton option. 
   
   
 
<br/><br/><table><tr><td><strong>PDF:piSetB  </td><td></td><td> <input type="text" name="14" value="void" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>void</strong></code>)</td></tr></table>
Parton density for <i>pion beam B</i>. If this option is set 
to <code>void</code> then the same PDF set as <code>PDF:piSet</code> 
is used. 
   
 
<a name="section4"></a> 
<h3>Parton densities for Pomerons</h3> 
 
The Pomeron is introduced in the description of diffractive events, 
i.e. a diffractive system is viewed as a Pomeron-proton collision at a 
reduced CM energy. Here the PDF's are even less well known. 
Most experimental parametrizations are NLO, which makes them less 
well suited for Monte Carlo applications. Furthermore note that 
the momentum sum is arbitrarily normalized to a non-unity value. 
 
<br/><br/><table><tr><td><strong>PDF:PomSet  </td><td></td><td> <input type="text" name="15" value="6" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>6</strong></code>)</td></tr></table>
Parton densities that can be used for Pomeron beams. 
<br/><code>option </code><strong> 1</strong> : <i>Q^2</i>-independent parametrizations 
<i>xf(x) = N_ab x^a (1 - x)^b</i>, where <i>N_ab</i> ensures 
unit momentum sum. The <i>a</i> and <i>b</i> parameters can be 
set separately for the gluon and the quark distributions. The 
momentum fraction of gluons and quarks can be freely mixed, and 
production of <i>s</i> quarks can be suppressed relative to 
that of <i>d</i> and <i>u</i> ones, with antiquarks as likely 
as quarks. See further below how to set the six parameters of this 
approach. 
   
<br/><code>option </code><strong> 2</strong> : <i>pi0</i> distributions, as specified in the 
section above. 
   
<br/><code>option </code><strong> 3</strong> : the H1 2006 Fit A NLO <i>Q^2</i>-dependent 
parametrization, based on a tune to their data [<a href="Bibliography.php#refH1P06" target="page">H1P06</a>], 
rescaled by the factor <code>PomRescale</code> below. 
   
<br/><code>option </code><strong> 4</strong> : the H1 2006 Fit B NLO <i>Q^2</i>-dependent 
parametrization, based on a tune to their data [<a href="Bibliography.php#refH1P06" target="page">H1P06</a>], 
rescaled by the factor <code>PomRescale</code> below. 
   
<br/><code>option </code><strong> 5</strong> : the H1 2007 Jets NLO <i>Q^2</i>-dependent 
parametrization, based on a tune to their data [<a href="Bibliography.php#refH1P07" target="page">H1P07</a>], 
rescaled by the factor <code>PomRescale</code> below. 
   
<br/><code>option </code><strong> 6</strong> : the H1 2006 Fit B LO <i>Q^2</i>-dependent 
parametrization, based on a tune to their data [<a href="Bibliography.php#refH1P06" target="page">H1P06</a>], 
rescaled by the factor <code>PomRescale</code> below. 
   
<br/><code>option </code><strong> 7</strong> : the ACTW B NLO <i>Q^2</i>-dependent 
parametrization with <i>epsilon=0.14</i>, 
based on a tune to H1 and ZEUS data [<a href="Bibliography.php#refAlv99" target="page">Alv99</a>], 
rescaled by the factor <code>PomRescale</code> below. 
   
<br/><code>option </code><strong> 8</strong> : the ACTW D NLO <i>Q^2</i>-dependent 
parametrization with <i>epsilon=0.14</i>, 
based on a tune to H1 and ZEUS data [<a href="Bibliography.php#refAlv99" target="page">Alv99</a>], 
rescaled by the factor <code>PomRescale</code> below. 
   
<br/><code>option </code><strong> 9</strong> : the ACTW SG NLO <i>Q^2</i>-dependent 
parametrization with <i>epsilon=0.14</i>, 
based on a tune to H1 and ZEUS data [<a href="Bibliography.php#refAlv99" target="page">Alv99</a>], 
rescaled by the factor <code>PomRescale</code> below. 
   
<br/><code>option </code><strong> 10</strong> : the ACTW D NLO <i>Q^2</i>-dependent 
parametrization with <i>epsilon=0.19</i>, 
based on a tune to H1 and ZEUS data [<a href="Bibliography.php#refAlv99" target="page">Alv99</a>], 
rescaled by the factor <code>PomRescale</code> below. 
   
<br/><code>option </code><strong> 11</strong> : a rescaling of the proton PDF, 
<i>xf^pom(x)=xf^p(x x_pom)</i>, used in the 
<code>Angantyr</code> model for Heavy Ion collisions. For high <i>x</i> 
there is an additional suppression by <i>(1-x)^p</i>, where the power is 
given by <code>PDF:PomHixSupp</code> below. 
   
<br/><code>option </code><strong> 12</strong> : The GKG18-DPDF LO Fit A central member 
<i>Q^2</i>-dependent parametrization based on a tune to 
H1 and ZEUS data [<a href="Bibliography.php#refGoh18" target="page">Goh18</a>]. 
   
<br/><code>option </code><strong> 13</strong> : The GKG18-DPDF LO Fit B central member 
<i>Q^2</i>-dependent parametrization based on a tune to 
H1 and ZEUS data [<a href="Bibliography.php#refGoh18" target="page">Goh18</a>]. 
   
<br/><code>option </code><strong> 14</strong> : The GKG18-DPDF NLO Fit A central member 
<i>Q^2</i>-dependent parametrization based on a tune to 
H1 and ZEUS data [<a href="Bibliography.php#refGoh18" target="page">Goh18</a>]. 
   
<br/><code>option </code><strong> 15</strong> : The GKG18-DPDF NLO Fit B central member 
<i>Q^2</i>-dependent parametrization based on a tune to 
H1 and ZEUS data [<a href="Bibliography.php#refGoh18" target="page">Goh18</a>]. 
   
<br/><code>option </code><strong> LHAPDF5:set/member</strong> : Use an external LHAPDF5 set, 
cf. the corresponding proton option. 
   
<br/><code>option </code><strong> LHAPDF6:set/member</strong> : Use an external LHAPDF6 set, 
cf. the corresponding proton option. 
   
<br/><code>option </code><strong> LHAGrid1:filename</strong> : Use the internal implementation 
for a LHAPDF6 set, cf. the corresponding proton option. 
   
   
 
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
Rescale several of the fits above by this uniform factor, e.g. to bring 
up their momentum sum to around unity. By default many of the sets have 
a momentum sum of order 0.5, suggesting that a factor around 2.0 
should be used. You can use <code>examples/main51.cc</code> to get 
a more precise value. Note that also other parameters in the 
<?php $filepath = $_GET["filepath"];
echo "<a href='Diffraction.php?filepath=".$filepath."' target='page'>";?>diffraction</a> framework may need to 
be retuned when this parameter is changed. Specifically 
<code>Diffraction:PomFluxRescale</code> should be set to the inverse 
of <code>PDF:PomRescale</code> to preserve the cross section for hard 
diffractive processes. 
   
 
<br/><br/><table><tr><td><strong>PDF:PomHixSupp </td><td></td><td> <input type="text" name="23" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>; <code>minimum = 0.</code>; <code>maximum = 10.</code>)</td></tr></table>
the power in the suppression of the high-x PDF for option 11 above. 
   
 
<a name="section5"></a> 
<h3>Parton densities for photons</h3> 
 
Photon PDFs describe the partonic content of the resolved photons and 
can be used to generate any process initiated by quarks and gluons. 
 
<p/> 
There are several PDF sets available for photons, although there have not 
been much activity recently. Currently one internal set is included, but 
more sets are available from LHAPDF5. The sets from LHAPDF5 can only be 
used as PDFs in the hard process (see <code>PDF:GammaHardSet</code> below). 
In case of photons the parton shower and beam remnant generation 
require additional methods that are provided only for internal sets. 
Currently no photon PDFs have been included in LHAPDF6. 
 
<br/><br/><table><tr><td><strong>PDF:GammaSet  </td><td>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = 1</code>; <code>maximum = 1</code>)</td></tr></table>
Parton densities that can be used for resolved photon beams. 
<br/>
<input type="radio" name="24" value="1" checked="checked"><strong>1 </strong>:  CJKL, based on <ref>Cor03</ref> but the rescaling  for heavy quarks due to kinematic constraints in DIS is undone to obtain  correct behaviour for photon-photon/hadron collisions.<br/>
 
<br/><br/><table><tr><td><strong>PDF:GammaHardSet  </td><td></td><td> <input type="text" name="25" value="void" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>void</strong></code>)</td></tr></table>
Parton densities to be used by the beams of the hard process. For photons 
the other options are the ones provided by LHAPDF5. If this option is set 
to <code>void</code> then the same PDF set as <code>PDF:GammaSet</code> is 
used. 
   
 
<a name="section6"></a> 
<h3>Parton densities for leptons</h3> 
 
For electrons/muons/taus there is no need to choose between different 
parametrizations, since only one implementation is available, and 
should be rather uncontroversial (apart from some technical details). 
However, insofar as e.g. <i>e^+ e^-</i> data often are corrected 
back to a world without any initial-state photon radiation, it is 
useful to have a corresponding option available here. 
 
<br/><br/><strong>PDF:lepton</strong>  <input type="radio" name="26" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="26" value="off"><strong>Off</strong>
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
unresolved lepton beam. For lepton-neutrino collisions to work you must 
therefore set <code>PDF:lepton = off</code>. 
 
<h4>Photons from lepton beams</h4> 
 
Lepton beams can emit photons and therefore may have partonic 
content. The PDFs describing these can be obtained by convoluting 
the photon flux with the selected photon PDFs. The photon flux 
is modelled according to equivalent photon approximation (EPA) which 
gives the flux of bremsstrahlung photons. 
 
<br/><br/><strong>PDF:lepton2gamma</strong>  <input type="radio" name="27" value="on"><strong>On</strong>
<input type="radio" name="27" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Gives photon beams from leptons. Both, unresolved (direct) and resolved 
contributions are included, see <?php $filepath = $_GET["filepath"];
echo "<a href='Photoproduction.php?filepath=".$filepath."' target='page'>";?> 
Photoproduction</a> for details. Can be used only with charged leptons. 
The applied photon PDF set is selected with the <code>PDF:GammaSet</code> 
and <code>PDF:GammaHardSet</code> options above. Events with two unresolved 
photon initiators can be generated also with the <code>PDF:lepton = on</code> 
but then additional phase-space cuts (e.g. cut on the invariant mass of the 
photon-photon pair) are not applied. 
   
 
<br/><br/><table><tr><td><strong>PDF:lepton2gammaSet  </td><td>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = 1</code>; <code>maximum = 2</code>)</td></tr></table>
The type of photon flux. 
<br/>
<input type="radio" name="28" value="1" checked="checked"><strong>1 </strong>:  Convolute the photon flux from EPA with the selected photon  PDF set. Convolution integral is performed "on the fly", meaning that the  actual integral is not computed but the <ei>x_gamma</ei> is sampled  event-by-event. Since the final PDF value depends on the sampled value for  <ei>x_gamma</ei> the phase-space sampling is set up using an overestimate for  the PDFs. This makes the process selection somewhat less efficient compared  to the case where the PDFs are fixed (e.g. for protons).<br/>
<input type="radio" name="28" value="2"><strong>2 </strong>:  Uses an approximation of the photon flux to sample  processes and corrects this later with an externally provided flux. For  leptons a bit less efficient than option 1 but allows straightforward  implementation of photon fluxes from different particles. To use this option  user has to provide the external photon flux using method  <code>Pythia::setPhotonFluxPtr(PDF*, PDF*)</code> as demostrated in the  sample program <code>main70.cc</code>.<br/>
 
<br/><br/><table><tr><td><strong>PDF:lepton2gammaApprox  </td><td>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = 1</code>; <code>maximum = 2</code>)</td></tr></table>
Controls which type of overestimate is used for photon flux sampling. 
<br/>
<input type="radio" name="29" value="1" checked="checked"><strong>1 </strong>:  Estimate optimized for photons from leptons. Works  reasonable well also for photoproduction in p+p.<br/>
<input type="radio" name="29" value="2"><strong>2 </strong>:  Estimate optimized for ultraperipheral heavy-ion collisions  as presented in <code>main70.cc</code>. Here the estimate is divided into two  regions, a power-law in <ei>x_gamma</ei> below <ei>x_cut</ei> and an  exponential fall-off above, see the related parameters below. Default values  are optimized for p+Pb collisions where Pb-nucleus provide the photons but  they should work reasonably well also for other similar  configurations.<br/>
<br/><b>Note:</b> Parameters do not affect the flux itself, only the sampling 
efficiency. 
 
<br/><br/><table><tr><td><strong>PDF:gammaFluxApprox2bMin </td><td></td><td> <input type="text" name="30" value="7.336" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>7.336</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
Minimal allowed impact parameter for which the flux is considered. Units in 
<code>fm</code>. Should match the flux provided by user. 
   
 
<br/><br/><table><tr><td><strong>PDF:gammaFluxApprox2mBeam </td><td></td><td> <input type="text" name="31" value="0.9314" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.9314</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
Per-nucleon mass used for the overestimate. Units in 
<code>GeV</code> and should again match to the user-provided flux. 
   
 
<br/><br/><table><tr><td><strong>PDF:gammaFluxApprox2xPow </td><td></td><td> <input type="text" name="32" value="1.15" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.15</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
Value of the exponent of the power law. The default value should work well 
for the foreseen cases, so vary with caution. 
   
 
<br/><br/><table><tr><td><strong>PDF:gammaFluxApprox2xCut </td><td></td><td> <input type="text" name="33" value="0.01" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.01</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
Value that defines at which <i>x_gamma</i> different approximations are 
used. As above, vary with caution. 
   
 
<a name="section7"></a> 
<h3>Incoming parton selection</h3> 
 
There is one useful degree of freedom to restrict the set of incoming 
quark flavours for hard processes. It does not change the PDF's as such, 
only which quarks are allowed to contribute to the hard-process cross 
sections. Note that separate but similarly named modes are available 
for multiparton interactions and spacelike showers. 
 
<br/><br/><table><tr><td><strong>PDFinProcess:nQuarkIn  </td><td></td><td> <input type="text" name="34" value="5" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>5</strong></code>; <code>minimum = 0</code>; <code>maximum = 5</code>)</td></tr></table>
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

if($_POST["1"] != "off")
{
$data = "PDF:extrapolate = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "13")
{
$data = "PDF:pSet = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "void")
{
$data = "PDF:pSetB = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "off")
{
$data = "PDF:useHard = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "void")
{
$data = "PDF:pHardSet = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "void")
{
$data = "PDF:pHardSetB = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
if($_POST["7"] != "off")
{
$data = "PDF:useHardNPDFA = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
if($_POST["8"] != "off")
{
$data = "PDF:useHardNPDFB = ".$_POST["8"]."\n";
fwrite($handle,$data);
}
if($_POST["9"] != "1")
{
$data = "PDF:nPDFSetA = ".$_POST["9"]."\n";
fwrite($handle,$data);
}
if($_POST["10"] != "1")
{
$data = "PDF:nPDFSetB = ".$_POST["10"]."\n";
fwrite($handle,$data);
}
if($_POST["11"] != "100822080")
{
$data = "PDF:nPDFBeamA = ".$_POST["11"]."\n";
fwrite($handle,$data);
}
if($_POST["12"] != "100822080")
{
$data = "PDF:nPDFBeamB = ".$_POST["12"]."\n";
fwrite($handle,$data);
}
if($_POST["13"] != "1")
{
$data = "PDF:piSet = ".$_POST["13"]."\n";
fwrite($handle,$data);
}
if($_POST["14"] != "void")
{
$data = "PDF:piSetB = ".$_POST["14"]."\n";
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
if($_POST["23"] != "0.")
{
$data = "PDF:PomHixSupp = ".$_POST["23"]."\n";
fwrite($handle,$data);
}
if($_POST["24"] != "1")
{
$data = "PDF:GammaSet = ".$_POST["24"]."\n";
fwrite($handle,$data);
}
if($_POST["25"] != "void")
{
$data = "PDF:GammaHardSet = ".$_POST["25"]."\n";
fwrite($handle,$data);
}
if($_POST["26"] != "on")
{
$data = "PDF:lepton = ".$_POST["26"]."\n";
fwrite($handle,$data);
}
if($_POST["27"] != "off")
{
$data = "PDF:lepton2gamma = ".$_POST["27"]."\n";
fwrite($handle,$data);
}
if($_POST["28"] != "1")
{
$data = "PDF:lepton2gammaSet = ".$_POST["28"]."\n";
fwrite($handle,$data);
}
if($_POST["29"] != "1")
{
$data = "PDF:lepton2gammaApprox = ".$_POST["29"]."\n";
fwrite($handle,$data);
}
if($_POST["30"] != "7.336")
{
$data = "PDF:gammaFluxApprox2bMin = ".$_POST["30"]."\n";
fwrite($handle,$data);
}
if($_POST["31"] != "0.9314")
{
$data = "PDF:gammaFluxApprox2mBeam = ".$_POST["31"]."\n";
fwrite($handle,$data);
}
if($_POST["32"] != "1.15")
{
$data = "PDF:gammaFluxApprox2xPow = ".$_POST["32"]."\n";
fwrite($handle,$data);
}
if($_POST["33"] != "0.01")
{
$data = "PDF:gammaFluxApprox2xCut = ".$_POST["33"]."\n";
fwrite($handle,$data);
}
if($_POST["34"] != "5")
{
$data = "PDFinProcess:nQuarkIn = ".$_POST["34"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>
 
<!-- Copyright (C) 2019 Torbjorn Sjostrand --> 
