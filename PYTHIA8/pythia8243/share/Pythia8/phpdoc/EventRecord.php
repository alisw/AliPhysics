<html>
<head>
<title>The Event Record</title>
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

<form method='post' action='EventRecord.php'>
 
<h2>The Event Record</h2> 
<ol id="toc">
  <li><a href="#section0">Basic output methods</a></li>
  <li><a href="#section1">Further output methods</a></li>
  <li><a href="#section2">Constructors and modifications of the event record</a></li>
  <li><a href="#section3">The Junction Class</a></li>
  <li><a href="#section4">Subsystems</a></li>
</ol>

 
A <code>Pythia</code> instance contains two members of the 
<code>Event</code> class. The one called <code>process</code> provides 
a brief summary of the main steps of the hard process, while the 
one called <code>event</code> contains the full history. The 
user would normally interact mainly with the second one, so 
we will exemplify primarily with that one. 
 
<p/> 
The <code>Event</code> class to first approximation is a vector of 
<code>Particle</code>s, so that it can expand to fit the current 
event size. The index operator is overloaded, so that e.g. 
<code>event[i]</code> corresponds to the <i>i</i>'th particle 
of the object <code>event</code>. Thus <code>event[i].id()</code> 
returns the identity of the <i>i</i>'th particle, and so on. 
Therefore the methods of the 
<code><?php $filepath = $_GET["filepath"];
echo "<a href='ParticleProperties.php?filepath=".$filepath."' target='page'>";?>Particle</a></code> class 
are at least as essential as those of the <code>Event</code> class 
itself. 
 
<p/> 
As used inside PYTHIA, some conventions are imposed on the structure 
of the event record. Entry 0 of the <code>vector&lt;Particle&gt;</code> 
is used to represent the event as a whole, with its total four-momentum 
and invariant mass, but does not form part of the event history. 
Lines 1 and 2 contains the two incoming beams, and only from here on 
history tracing works as could be expected. That way unassigned mother 
and daughter indices can be put 0 without ambiguity. Depending on the 
task at hand, a loop may therefore start at index 1 rather than 0 
without any loss. Specifically, for translation to other event record 
formats such as HepMC [<a href="Bibliography.php#refDob01" target="page">Dob01</a>], where the first index is 1, the 
Pythia entry 0 definitely ought to be skipped in order to minimize the 
danger of indexing errors. 
 
<p/> 
In the following we will list the methods available. 
Only a few of them have a function to fill in normal user code. 
 
<a name="section0"></a> 
<h3>Basic output methods</h3> 
 
Some methods are available to read out information on the 
current event record: 
 
<a name="anchor1"></a>
<p/><strong> Particle& Event::operator[](int i) &nbsp;</strong> <br/>
   
<a name="anchor2"></a>
<strong> const Particle& Event::operator[](int i) &nbsp;</strong> <br/>
returns a (<code>const</code>) reference to the <i>i</i>'th particle 
in the event record, which can be used to get (or set) all the 
<?php $filepath = $_GET["filepath"];
echo "<a href='ParticleProperties.php?filepath=".$filepath."' target='page'>";?>properties</a> of this particle. 
   
 
<a name="anchor3"></a>
<p/><strong> Particle& Event::front() &nbsp;</strong> <br/>
   
<a name="anchor4"></a>
<strong> Particle& Event::at(int i) &nbsp;</strong> <br/>
   
<a name="anchor5"></a>
<strong> Particle& Event::back() &nbsp;</strong> <br/>
   
returns a reference to the zeroth, <i>i</i>'th or last particle 
in the event record, as an alternative to the methods above. 
 
<a name="anchor6"></a>
<p/><strong> int Event::size() &nbsp;</strong> <br/>
The event size, i.e. the size of the <code>vector&lt;Particle&gt;</code>. 
Thus valid particles, to be accessed by the above indexing operator, 
are stored in the range <i>0 &lt;= i &lt; size()</i>. See comment 
above about the (ir)relevance of entry 0. 
   
 
<a name="anchor7"></a>
<p/><strong> void Event::list(bool showScaleAndVertex = false, bool showMothersAndDaughters = false, int precision = 3) &nbsp;</strong> <br/>
Provide a listing of the whole event, i.e. of the 
<code>vector&lt;Particle&gt;</code>. The basic identity 
code, status, mother, daughter, colour, four-momentum and mass data 
are always given, but the methods can also be called with a few 
optional arguments for further information: 
<br/><code>argument</code><strong> showScaleAndVertex </strong> (<code>default = <strong>off</strong></code>) :  optionally give a 
second line for each particle, with the production scale (in GeV), 
the particle polarization (dimensionless), the production vertex 
(in mm or mm/c) and the invariant lifetime (also in mm/c). 
   
<br/><code>argument</code><strong> showMothersAndDaughters </strong> (<code>default = <strong>off</strong></code>) :  
gives a list of all daughters and mothers of a particle, as defined by 
the <code>motherList(i)</code> and <code>daughterList(i)</code> methods 
described below. It is mainly intended for debug purposes. 
   
<br/><code>argument</code><strong> precision </strong> (<code>default = <strong>3</strong></code>) :  the number of digits to the right 
of the decimal point shown for momenta, energies andf masses. Can be set 
above 3, but reducing it below 3 will have no effect. This option is 
intended for expert users, e.g. for debugging purposes, and so no effort 
has been made to stretch header and footer to match. 
   
   
 
<a name="section1"></a> 
<h3>Further output methods</h3> 
 
Many event properties are accessible via the <code>Info</code> class, 
<?php $filepath = $_GET["filepath"];
echo "<a href='EventInformation.php?filepath=".$filepath."' target='page'>";?>see here</a>. Since they are used 
directly in the event generation, a few are stored directly in the 
<code>Event</code> class, however. 
 
<a name="anchor8"></a>
<p/><strong> void Event::scale( double scaleIn) &nbsp;</strong> <br/>
   
<a name="anchor9"></a>
<strong> double Event::scale() &nbsp;</strong> <br/>
set or get the scale (in GeV) of the hardest process in the event. 
Matches the function of the <code>scale</code> variable in the 
<?php $filepath = $_GET["filepath"];
echo "<a href='LesHouchesAccord.php?filepath=".$filepath."' target='page'>";?>Les Houches Accord</a>. 
   
 
<a name="anchor10"></a>
<p/><strong> void Event::scaleSecond( double scaleSecondIn) &nbsp;</strong> <br/>
   
<a name="anchor11"></a>
<strong> double Event::scaleSecond() &nbsp;</strong> <br/>
set or get the scale (in GeV) of a second hard process in the event, 
in those cases where such a one 
<?php $filepath = $_GET["filepath"];
echo "<a href='SecondHardProcess.php?filepath=".$filepath."' target='page'>";?>has been requested</a>. 
   
 
<p/> 
One data member in an <code>Event</code> object is used to keep track 
of the largest <code>col()</code> or <code>acol()</code> colour tag set 
so far, so that new ones do not clash. 
 
<br/><br/><table><tr><td><strong>Event:startColTag  </td><td></td><td> <input type="text" name="1" value="100" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>100</strong></code>; <code>minimum = 0</code>; <code>maximum = 1000</code>)</td></tr></table>
This sets the initial colour tag value used, so that the first one 
assigned is <code>startColTag + 1</code>, etc. The Les Houches accord 
[<a href="Bibliography.php#refBoo01" target="page">Boo01</a>] suggests this number to be 500, but 100 works equally 
well. 
   
 
<a name="anchor12"></a>
<p/><strong> void Event::initColTag(int colTag = 0) &nbsp;</strong> <br/>
forces the current colour tag value to be the larger of the input 
<code>colTag</code> and the above <code>Event:startColTag</code> 
values. 
   
 
<a name="anchor13"></a>
<p/><strong> int Event::lastColTag() &nbsp;</strong> <br/>
returns the current maximum colour tag. 
   
 
<a name="anchor14"></a>
<p/><strong> int Event::nextColTag() &nbsp;</strong> <br/>
increases the current maximum colour tag by one and returns this 
new value. This method is used whenever a new colour tag is needed. 
   
 
<a name="section2"></a> 
<h3>Constructors and modifications of the event record</h3> 
 
Although you would not normally need to create your own 
<code>Event</code> instance, there may be times where that 
could be convenient. The typical example would be if you want to 
create a new event record that is the sum of a few different ones, 
e.g. if you want to simulate pileup events. There may also be cases 
where you want to add one or a few particles to an existing event 
record. 
 
<a name="anchor15"></a>
<p/><strong> Event::Event(int capacity = 100) &nbsp;</strong> <br/>
creates an empty event record, but with a reserved size 
<i>capacity</i> for the <code>Particle</code> vector. 
   
 
<a name="anchor16"></a>
<p/><strong> Event& Event::operator=(const Event& oldEvent) &nbsp;</strong> <br/>
copies the input event record. 
   
 
<a name="anchor17"></a>
<p/><strong> Event& Event::operator+=(const Event& addEvent) &nbsp;</strong> <br/>
appends an event to an existing one. For the appended particles 
mother, daughter and colour tags are shifted to make a consistent 
record. The zeroth particle of the appended event is not copied, 
but the zeroth particle of the combined event is updated to the 
full energy-momentum content. 
   
 
<a name="anchor18"></a>
<p/><strong> void Event::init(string headerIn = &quot;&quot;, ParticleData* particleDataPtrIn = 0, int startColTagIn = 100) &nbsp;</strong> <br/>
initializes colour, the pointer to the particle database, and the 
header specification used for the event listing. We remind that a 
<code>Pythia</code> object contains two event records 
<code>process</code> and <code>event</code>. Thus one may e.g. 
call either  <code>pythia.process.list()</code> or 
<code>pythia.event.list()</code>. To distinguish those two rapidly 
at visual inspection, the <code>"Pythia Event Listing"</code> header 
is printed out differently, in one case adding 
<code>"(hard process)"</code> and in the other 
<code>"(complete event)"</code>. When <code>+=</code> is used to 
append an event, the modified event is printed with 
<code>"(combination of several events)"</code> as a reminder. 
   
 
<a name="anchor19"></a>
<p/><strong> void Event::clear() &nbsp;</strong> <br/>
empties event record. Specifically the <code>Particle</code> vector 
size is reset to zero. 
   
 
<a name="anchor20"></a>
<p/><strong> void Event::free() &nbsp;</strong> <br/>
empties event record, like <code>clear()</code> above, but also frees 
up the memory of the <code>Particle</code> vector. 
   
 
<a name="anchor21"></a>
<p/><strong> void Event::reset() &nbsp;</strong> <br/>
empties the event record, as <code>clear()</code> above, but then 
fills the zero entry of the <code>Particle</code> vector with the 
pseudoparticle used to represent the event as a whole. At this point 
the pseudoparticle is not assigned any momentum or mass. 
   
 
<a name="anchor22"></a>
<p/><strong> void Event::popBack(int n = 1) &nbsp;</strong> <br/>
removes the last <i>n</i> particle entries; must be a positive 
number. History (and other) information of remaning entries is 
untouched, and so may be internally inconsistent. 
   
 
<a name="anchor23"></a>
<p/><strong> void Event::remove(int iFirst, int iLast, bool shiftHistory = true) &nbsp;</strong> <br/>
removes particles in the range between indices <code>iFirst</code> 
and <code>iLast</code>, including the endpoints. By default all mother 
and daughter indices above the removed range are shifted down by the 
number of removed entries, while indices in the removed range are put 
zero. Optionally these shifts can be omitted. Other information remains 
unchanged, which may lead to inconsistencies. If the decay products of 
a particle are removed, e.g., the mother particle status should be set 
positive, cf. <code>Particle::undoDecay()</code>. 
   
 
<a name="anchor24"></a>
<p/><strong> int Event::append(Particle entryIn) &nbsp;</strong> <br/>
appends a particle to the bottom of the event record and 
returns the index of this position. 
   
 
<a name="anchor25"></a>
<p/><strong> int Event::append(int id, int status, int mother1, int mother2, int daughter1, int daughter2, int col, int acol, double px, double py, double pz,  double e, double m = 0., double scale = 0., double pol = 9.) &nbsp;</strong> <br/>
appends a particle to the bottom of the event record and 
returns the index of this position; 
<?php $filepath = $_GET["filepath"];
echo "<a href='ParticleProperties.php?filepath=".$filepath."' target='page'>";?>see here</a> for the meaning 
of the various particle properties. 
   
 
<a name="anchor26"></a>
<p/><strong> int Event::append(int id, int status, int mother1, int mother2, int daughter1, int daughter2, int col, int acol, Vec4 p, double m = 0., double scale = 0., double pol = 9.) &nbsp;</strong> <br/>
appends a particle to the bottom of the event record and 
returns the index of this position, as above but with four-momentum 
as a <code>Vec4</code>. 
   
 
<a name="anchor27"></a>
<p/><strong> int Event::append(int id, int status, int col, int acol, double px, double py, double pz, double e, double m = 0., double scale = 0., double pol = 9.) &nbsp;</strong> <br/>
   
<a name="anchor28"></a>
<strong> int Event::append(int id, int status, int col, int acol, Vec4 p, double m = 0., double scale = 0., double pol = 9.) &nbsp;</strong> <br/>
appends a particle to the bottom of the event record and 
returns the index of this position, as above but with vanishing 
(i.e. zero) mother and daughter indices. 
   
 
<a name="anchor29"></a>
<p/><strong> int Event::setEvtPtr(int iSet = -1) &nbsp;</strong> <br/>
send in the <code>this</code> pointer of the current <code>Event</code> 
itself to the particle <code>iSet</code>, by default the most recently 
appended particle. Also generates a pointer to the 
<code>ParticleDataEntry</code> object of the identity code of the particle. 
   
 
<a name="anchor30"></a>
<p/><strong> int Event::copy(int iCopy, int newStatus = 0) &nbsp;</strong> <br/>
copies the existing particle in entry <code>iCopy</code> to the 
bottom of the event record and returns the index of this position. 
By default, i.e. with <code>newStatus = 0</code>, everything is 
copied precisely as it is, which means that history information 
has to be modified further by hand to make sense. With a positive 
<code>newStatus</code>, the new copy is set up to be the daughter of 
the old, with status code <code>newStatus</code>, while the status 
code of <code>iCopy</code> is negated. With a negative 
<code>newStatus</code>, the new copy is instead set up to be the 
mother of <code>iCopy</code>. An attempt to copy an out-of-range 
entry will return -1. 
   
 
<a name="anchor31"></a>
<p/><strong> void Event::restorePtrs() &nbsp;</strong> <br/>
each particle in the event record has a pointer to the event itself 
and another to the particle species it belongs to. The latter pointer 
is automatically set/changed whenever the particle identity is 
set/changed by one of the normal methods. Of course the pointer values 
are specific to the memory locations of the current run, and so it has 
no sense to save them if events are written to file. Should you use some 
persistency scheme that bypasses the normal methods when the event is 
read back in, you can use <code>restorePtrs()</code> afterwards to set 
these pointers appropriately. 
   
 
<a name="anchor32"></a>
<p/><strong> int Event::nFinal(bool chargedOnly = false) &nbsp;</strong> <br/>
return the number of final-state particles in the current event record, 
optionally counting charged particles only. 
   
 
<a name="anchor33"></a>
<p/><strong> double Event::dyAbs(int i1, int i2) &nbsp;</strong> <br/>
   
<a name="anchor34"></a>
<strong> double Event::detaAbs(int i1, int i2) &nbsp;</strong> <br/>
   
<a name="anchor35"></a>
<strong> double Event::dphiAbs(int i1, int i2) &nbsp;</strong> <br/>
   
<a name="anchor36"></a>
<strong> double Event::RRapPhi(int i1, int i2) &nbsp;</strong> <br/>
   
<a name="anchor37"></a>
<strong> double Event::REtaPhi(int i1, int i2) &nbsp;</strong> <br/>
return the separation between two particles in the event record, 
in true rapidity, in pseudorapidity, in phi angle, and in the <i>R</i> 
combination of true or pseudorapidity with phi, when required with 
absolute sign so as to avoid negative numbers. 
   
 
<p/> 
A few methods exist to rotate and boost events. These derive from the 
<?php $filepath = $_GET["filepath"];
echo "<a href='FourVectors.php?filepath=".$filepath."' target='page'>";?>Vec4</a> methods, and affect both the 
momentum and the vertex (position) components of all particles. 
 
<a name="anchor38"></a>
<p/><strong> void Event::rot(double theta, double phi) &nbsp;</strong> <br/>
rotate all particles in the event by this polar and azimuthal angle 
(expressed in radians). 
   
 
<a name="anchor39"></a>
<p/><strong> void Event::bst(double betaX, double betaY, double betaZ) &nbsp;</strong> <br/>
   
<a name="anchor40"></a>
<strong> void Event::bst(double betaX, double betaY, double betaZ, double gamma) &nbsp;</strong> <br/>
   
<a name="anchor41"></a>
<strong> void Event::bst(const Vec4& vec) &nbsp;</strong> <br/>
boost all particles in the event by this three-vector. 
Optionally you may provide the <i>gamma</i> value as a fourth argument, 
which may help avoid roundoff errors for big boosts. You may alternatively 
supply a <code>Vec4</code> four-vector, in which case the boost vector 
becomes <i>beta = p/E</i>. 
   
 
<a name="anchor42"></a>
<p/><strong> void Event::rotbst(const RotBstMatrix& M, bool boostVertices = true) &nbsp;</strong> <br/>
rotate and boost by the combined action encoded in the 
<code><?php $filepath = $_GET["filepath"];
echo "<a href='FourVectors.php?filepath=".$filepath."' target='page'>";?>RotBstMatrix</a> M</code>. 
If the optional second argument is false only the four-momenta are 
boosted, and not the production vertices. 
   
 
<a name="section3"></a> 
<h3>The Junction Class</h3> 
 
The event record also contains a vector of junctions, which often 
is empty or else contains only a very few per event. Methods are 
available to add further junctions or query the current junction list. 
This is only for the expert user, however, and is not discussed 
further here, but only the main points. 
 
<p/> 
A junction stores the properties associated with a baryon number that 
is fully resolved, i.e. where three different colour indices are 
involved. There are two main applications, 
<ol> 
<li>baryon beams, where at least two valence quarks are kicked out, 
and so the motion of the baryon number is nontrivial;</li> 
<li>baryon-number violating processes, e.g. in SUSY with broken 
<i>R</i>-parity.</li> 
</ol> 
Information on junctions is set, partly in the process generation, 
partly in the beam remnants machinery, and used by the fragmentation 
routines, but the normal user does not have to know the details. 
 
<p/> 
For each junction, information is stored on the kind of junction, and 
on the three (anti)colour indices that are involved in the junction. 
The possibilities foreseen are: 
<ul> 
<li><code>kind = 1</code> : incoming colourless particle to three 
outgoing colours (e.g. baryon beam remnant or 
<i>neutralino &rarr; q q q</i>);</li> 
<li><code>kind = 2</code> : incoming colourless particle to three 
outgoing anticolours;</li> 
<li><code>kind = 3</code> : one incoming anticolour (stored first) 
and two outgoing  colours (e.g. antisquark decaying to two quarks, or 
  gluino decay to three quarks);</li> 
<li><code>kind = 4</code> : one incoming colour (stored first) and two 
outgoing anticolours (e.g. squark decaying to two antiquarks, or 
  gluino decaying to three antiquarks);</li> 
<li><code>kind = 5</code> : two incoming anticolours (stored first) 
and one outgoing colour (e.g. resonant squark production through RPV);</li> 
<li><code>kind = 6</code> : two incoming colours (stored first) 
and one outgoing anticolour (e.g. resonant antisquark production 
  through RPV); 
</li> 
</ul> 
The odd (even) <code>kind</code> codes corresponds to a +1 (-1) change in 
baryon number across the junction. 
 
<p/> 
The kind and colour information in the list of junctions can be set 
or read with methods of the <code>Event</code> class, but are not of 
common interest and so not described here. 
 
<p/> 
A listing of current junctions can be obtained with the 
<code>listJunctions()</code> method. 
 
<a name="section4"></a> 
<h3>Subsystems</h3> 
 
Separate from the event record as such, but closely tied to it is the 
<code><?php $filepath = $_GET["filepath"];
echo "<a href='AdvancedUsage.php?filepath=".$filepath."' target='page'>";?>PartonSystems</a></code> class, 
which mainly stores the parton indices of incoming and outgoing partons, 
classified by collision subsystem. Such information is needed to 
interleave multiparton interactions, initial-state showers and final-state 
showers, and append beam remnants. It could also be used in other places. 
It is intended to be accessed only by experts, such as implementors of 
<?php $filepath = $_GET["filepath"];
echo "<a href='ImplementNewShowers.php?filepath=".$filepath."' target='page'>";?>new showering models</a>. 
 
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

if($_POST["1"] != "100")
{
$data = "Event:startColTag = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>
 
<!-- Copyright (C) 2019 Torbjorn Sjostrand --> 
