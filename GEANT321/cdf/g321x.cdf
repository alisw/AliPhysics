*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:40  cernlib
* Geant
*
*
*CMZ :  3.21/02 29/03/94  15.41.35  by  S.Giani
*-- Author :
>Menu GEANT
>Guidance
GEANT specific commands.
 
>Name GKDRAW
 
>Menu /GEANT/CVOL
>Guidance
Clipping commands.
The hidden line removal technique is necessary to visualize properly
very complex detectors. At the same time, it can be useful to visualize
the inner elements of a detector in detail. For this purpose, the
commands menu CVOL has been developed: these commands allow
subtractions (via boolean operation) of given shapes from any part of
the detector, therefore showing its inner contents. It is possible
to clip each different volume by means of a different shape (BOX ,
TUBE, CONE, SPHE are available). If '*' is given as the name of the
volume to be clipped, all volumes are clipped by the given shape.
A volume can be clipped at most twice (even by
different shapes); if a volume is explicitely clipped
twice, the '*' will not act on it anymore. Giving '.' as the name
of the volume to be clipped will reset the clipping.
 
>Command BOX
>Parameters
CNNV  ' Name of volume to be clipped          ' C  D='*   '
+
XMIN  ' Lower limit of the Shape X coordinate ' R  D=-10000.
XMAX  ' Upper limit of the Shape X coordinate ' R  D=-9999.
YMIN  ' Lower limit of the Shape Y coordinate ' R  D=-10000.
YMAX  ' Upper limit of the Shape Y coordinate ' R  D=-9999.
ZMIN  ' Lower limit of the Shape Z coordinate ' R  D=-10000.
ZMAX  ' Upper limit of the Shape Z coordinate ' R  D=-9999.
>Guidance
This command performs a boolean subtraction between the volume
CNVV and a box placed in the MARS according the values of the given
coordinates. See also CVOL.
The following commands will clip by a box,
with a vertex at the origin, the volume specified by NAME (a valid
string for the NAME of the volume can be found using the DTREE command).
 EXAMPLE -
 dopt hide on
 satt * seen -2
 draw NAME 40 40 0 10 10 .01 .01
 next
 box NAME 0 1000 0 1000 0 1000
 draw NAME 40 40 0 10 10 .01 .01
 box .
 
>Action GXDRAW
 
>Command TUBE
>Parameters
CNVV  ' Name of volume to be clipped          ' C  D='*   '
+
RMAX  ' External radius of tube               ' R  D=0.1
ZDEM  ' Half length of tube axis              ' R  D=0.1
XMED  ' Center X coordinate                   ' R  D=-10000.
YMED  ' Center Y coordinate                   ' R  D=-10000.
ZMED  ' Center Z coordinate                   ' R  D=-10000.
>Guidance
This command performs a boolean subtraction between the volume
CNVV and a tube; the tube has the given parameters and is placed in
the MARS according the given coordinates of its center.
See also CVOL.
The following commands will clip, by a tube,
positioned according to the given parameters, the volume specified
by NAME (a valid string for the NAME of the volume
can be found using the DTREE command).
 EXAMPLE -
 dopt hide on
 satt * seen -2
 draw NAME 40 40 0 10 10 .01 .01
 next
 tube * 500 1000 500 0 0
 draw NAME 40 40 0 10 10 .01 .01
 box .
 
>Action GXDRAW
 
>Command CONE
>Parameters
CNVV  ' Name of volume to be clipped          ' C  D='*   '
+
RMAX1 ' Min external radius                   ' R  D=0.1
RMAX2 ' Max external radius                   ' R  D=0.1
ZDEM  ' Half length of cone axis              ' R  D=0.1
XMED  ' Center X coordinate                   ' R  D=-10000.
YMED  ' Center Y coordinate                   ' R  D=-10000.
ZMED  ' Center Z coordinate                   ' R  D=-10000.
>Guidance
This command performs a boolean subtraction between the volume
CNVV and a cone; the cone has the given parameters and is placed in
the MARS according to the given coordinates of its center.
See also CVOL.
The following commands will clip by a cone,
positioned according the given parameters, the volume specified
by NAME (a valid string for the NAME of the volume
can be found using the DTREE command).
 EXAMPLE -
 dopt hide on
 satt * seen -2
 draw NAME 40 40 0 10 10 .01 .01
 next
 cone * 1 750 1000 0 0 1000
 draw NAME 40 40 0 10 10 .01 .01
 box .
 
>Action GXDRAW
 
>Command SPHE
>Parameters
CNVV  ' Name of volume to be clipped          ' C  D='*   '
+
RMAX  ' External radius of sphere             ' R  D=0.1
XMED  ' Center X coordinate                   ' R  D=-10000.
YMED  ' Center Y coordinate                   ' R  D=-10000.
ZMED  ' Center Z coordinate                   ' R  D=-10000.
>Guidance
This command performs a boolean subtraction between the volume
CNVV and a sphere; the sphere has the given parameters and is placed in
the MARS according to the given coordinates of its center.
See also CVOL. The following commands clip by a sphere,
positioned according to the given parameters, the volume specified
by NAME (a valid string for the NAME of the volume
can be found using the DTREE command).
EXAMPLE -
 dopt hide on
 satt * seen -2
 draw NAME 40 40 0 10 10 .01 .01
 next
 sphe * 500 0 0 500
 draw NAME 40 40 0 10 10 .01 .01
 box .
 
>Action GXDRAW
 
>Command VALCUT
>Parameters
XCUT 'x coordinate of cutted value' R D=0.
YCUT 'y coordinate of cutted value' R D=0.
ZCUT 'z coordinate of cutted value' R D=0.
>Guidance
It allows the cutting in the ray-tracing. All the volumes are cutted
from XCUT to +BIG along the x axis, from YCUT to +BIG along the y axis
and from ZCUT to +BIG along the z axis.
 
>Action GXDRAW
 
>Menu /GEANT/DRAWING
>Guidance
Drawing commands. These commands allow the visualization in several ways
of the volumes defined in the geometrical data structure. It is possible
to draw the logical tree of volumes belonging to the detector (DTREE),
to show their geometrical specification (DSPEC,DFSPC), to draw them
and their cut views (DRAW, DCUT). Moreover, it is possible to execute
these commands when the hidden line removal option is activated; in
this case, the volumes can be also either translated in the space
(SHIFT), or clipped by boolean operation (CVOL). In addition, it is
possible to fill the surfaces of the volumes
with solid colours when the shading option (SHAD) is activated.
Several tools (ZOOM, LENS) have been developed to zoom detailed parts
of the detectors or to scan physical events as well.
Finally, the command MOVE will allow the rotation, translation and zooming
on real time parts of the detectors or tracks and hits of a simulated event.
Ray-tracing commands. In case the command (DOPT RAYT ON) is executed,
the drawing is performed by the Geant ray-tracing;
automatically, the color is assigned according to the tracking medium of each
volume and the volumes with a density lower/equal than the air are considered
transparent; if the option (USER) is set (ON) (again via the command (DOPT)),
the user can set color and visibility for the desired volumes via the command
(SATT), as usual, relatively to the attributes (COLO) and (SEEN).
The resolution can be set via the command (SATT * FILL VALUE), where (VALUE)
is the ratio between the number of pixels drawn and 20 (user coordinates).
Parallel view and perspective view are possible (DOPT PROJ PARA/PERS); in the
first case, we assume that the first mother volume of the tree is a box with
dimensions 10000 X 10000 X 10000 cm and the view point (infinetely far) is
5000 cm far from the origin along the Z axis of the user coordinates; in the
second case, the distance between the observer and the origin of the world
reference system is set in cm by the command (PERSP NAME VALUE); grand-angle
or telescopic effects can be achieved changing the scale factors in the command
(DRAW). When the final picture does not occupy the full window,
mapping the space before tracing can speed up the drawing, but can also
produce less precise results; values from 1 to 4 are allowed in the command
(DOPT MAPP VALUE), the mapping being more precise for increasing (VALUE); for
(VALUE = 0) no mapping is performed (therefore max precision and lowest speed).
The command (VALCUT) allows the cutting of the detector by three planes
ortogonal to the x,y,z axis. The attribute (LSTY) can be set by the command
SATT for any desired volume and can assume values from 0 to 7; it determines
the different light processing to be performed for different materials:
0 = dark-matt, 1 = bright-matt, 2 = plastic, 3 = ceramic, 4 = rough-metals,
5 = shiny-metals, 6 = glass, 7 = mirror. The detector is assumed to be in the
dark, the ambient light luminosity is 0.2 for each basic hue (the saturation
is 0.9) and the observer is assumed to have a light source (therefore he will
produce parallel light in the case of parallel view and point-like-source
light in the case of perspective view).
 
>Command DRAW
>Parameters
NAME   'Volume name' C
+
THETA  'Viewing angle theta (for 3D projection)' R R=0.:180.
PHI    'Viewing angle phi (for 3D projection)' R R=0.:360.
PSI    'Viewing angle psi (for 2D rotation)' R R=0.:360.
U0     'U-coord. (horizontal) of volume origin' R
V0     'V-coord. (vertical) of volume origin' R
SU     'Scale factor for U-coord.' R
SV     'Scale factor for V-coord.' R
>Guidance
 CALL GDRAW(name,theta,phi,psi,u0,v0,su,sv)
If optional parameters are missing, the corresponding values are
taken from the common /GCDRAW/. This command will draw the volumes,
selected with their graphical attributes, set by the SATT
facility. The drawing may be performed with hidden line removal
and with shading effects according to the value of the options HIDE
and SHAD; if the option SHAD is ON, the contour's edges can be
drawn or not. If the option HIDE is ON, the detector can be
exploded (BOMB), clipped with different shapes (CVOL), and some
of its parts can be shifted from their original
position (SHIFT). When HIDE is ON, if
the drawing requires more than the available memory, the program
will evaluate and display the number of missing words
(so that the user can increase the
size of its ZEBRA store). Finally, at the end of each drawing (with HIDE on),
the program will print messages about the memory used and
statistics on the volumes' visibility.
The following commands will produce the drawing of a green
volume, specified by NAME, without using the hidden line removal
technique, using the hidden line removal technique,
with different linewidth and colour (red), with
solid colour, with shading of surfaces, and without edges.
Finally, some examples are given for the ray-tracing. (A possible
string for the NAME of the volume can be found using the command DTREE).
 EXAMPLE -
 satt * seen -2
 satt NAME colo 3
 draw NAME 40 40 0 10 10 .01 .01
 next
 dopt hide on
 draw NAME 40 40 0 10 10 .01 .01
 next
 satt NAME colo 2
 satt NAME lwid 4
 draw NAME 40 40 0 10 10 .01 .01
 next
 dopt shad on
 satt * lwid 1
 satt NAME fill 1
 draw NAME 40 40 0 10 10 .01 .01
 next
 satt NAME fill 3
 draw NAME 40 40 0 10 10 .01 .01
 next
 dopt edge off
 draw NAME 40 40 0 10 10 .01 .01
 dopt rayt on
 satt * fill 20
 dopt mapp 1
 draw NAME 40 40 0 10 10 .01 .01
 dopt proj pers
 persp NAME 500
 draw NAME 40 40 0 10 10 1 1
 valcut 100 100 100
 dopt mapp 0
 dopt user on
 satt NAM1 seen 0
 satt NAM2 colo 2
 draw NAME 40 40 0 10 10 5 5
 
>Action GXDRAW
 
>Command SPOT
>Parameters
XLPOS 'x coordinate of light source' R
YLPOS 'y coordinate of light source' R
ZLPOS 'z coordinate of light source' R
INTEN 'intensity of light source' I
>Guidance
This point-like light source can be moved in the space and its intensity
can be changed (INTEN going from 0 to 10) relatively to the ambience light.
>Action GXDRAW
 
>Command VAR5D
>Parameters
TSEQTO 'total sequential time' R
NPROC  'number of processors' I
NMPTOT 'number of message passing' I
TOTMBY 'total megabytes transfert' R
TSEQ   'not parallelized code' R
TLAT   'latency time' R
TNET   'network speed in Mbytes/sec' R
>Guidance
It sets the values of the parameters expressed in the formula and
specify which variables must be assumed as x,y,z (setting their value
to 1001,1002,1003, respectively).
>Action GXDRAW
 
>Command RANG5D
>Parameters
X1MIN 'x coordinate min' R
X1MAX 'x coordinate max' R
Y1MIN 'y coordinate min' R
Y1MAX 'y coordinate max' R
Z1MIN 'z coordinate min' R
Z1MAX 'z coordinate max' R
>Guidance
It sets the range for the x,y,z variables.
>Action GXDRAW
 
>Command DVOLUME
>Parameters
N      'Number of elements in arrays LNAMES and LNUMBS' I D=1
NAMNUM 'Volume names and numbers (ex. "NAME1,NR1,NAME2,NR2")' C
CHNRS    'Reference system used' C D='MARS' R='MARS,DRS'
+
THETA  'Viewing angle theta (for 3D projection)' R R=0.:360.
PHI    'Viewing angle phi (for 3D projection)' R R=0.:360.
PSI    'Viewing angle psi (for 2D rotation)' R R=0.:180.
U0     'U-coord. (horizontal) of volume origin' R
V0     'V-coord. (vertical) of volume origin' R
SU     'Scale factor for U-coord.' R
SV     'Scale factor for V-coord.' R
>Guidance
 CALL GDRVOL(n,lnames,lnumbs,nrs,theta,phi,psi,u0,v0,su,sv)
N is the number of levels from the top of the geometry structure
to the volume lnames(n),lnumbs(n) to be drawn.
NAMNUM contain the arrays lnames and lnumbs,
identifying the path, in pairs and separated by commas; for
example (with n=2) :
'lname(1),lnumbs(1),lname(2),lnumbs(2) '
CHNRS is the name of the reference system used: MARS for MAster Reference
System or DRS for Daughter Reference System.
NRS=0 for MARS or NRS<>0 for DRS
If optional parameters are missing, the current values in /GCDRAW/
are taken.
>Action GXDRAW
 
>Command DCUT
>Parameters
NAME   'Volume name' C
CAXIS  'Axis value' C R='X,Y,Z'
CUTVAL 'Cut plane distance from the origin along the axis' R
+
U0     'U-coord. (horizontal) of volume origin' R
V0     'V-coord. (vertical) of volume origin' R
SU     'Scale factor for U-coord.' R
SV     'Scale factor for V-coord.' R
>Guidance
 CALL GDRAWC(name,iaxis,cutval,u0,v0,su,sv)
The cut plane is normal to caxis (X,Y,Z), corresponding to iaxis (1,2,3),
and placed at the distance cutval from the origin.
The resulting picture is seen from the the same axis.
If optional parameters are missing, the current values in /GCDRAW/
are taken.
When HIDE Mode is ON, it is possible to get the same effect with
the CVOL/BOX command.
>Action GXDRAW
 
>Command DXCUT
>Parameters
NAME   'Volume name' C
CUTTHE 'Theta angle of the line normal to cut plane' R R=0.:360.
CUTPHI 'Phi angle of the line normal to cut plane' R R=0.:360.
CUTVAL 'Cut plane distance from the origin along the axis' R
+
THETA  'Viewing angle theta (for 3D projection)' R R=0.:360.
PHI    'Viewing angle phi (for 3D projection)' R R=0.:360.
U0     'U-coord. (horizontal) of volume origin' R
V0     'V-coord. (vertical) of volume origin' R
SU     'Scale factor for U-coord.' R
SV     'Scale factor for V-coord.' R
>Guidance
 CALL GDRAWX(name,cutthe,cutphi,cutval,theta,phi,u0,v0,su,sv)
The cut plane is normal to the line given by the cut angles
cutthe and cutphi and placed at the distance cutval from the origin.
The resulting picture is seen from the viewing angles theta,phi.
If optional parameters are missing, the current values in /GCDRAW/
are taken.
>Action GXDRAW
 
>Command SHIFT
>Parameters
CNVN  ' Name of volume to be shifted        ' C  D='*'
XXXX  ' Shift along X axis                  ' R  D=0.
YYYY  ' Shift along Y axis                  ' R  D=0.
ZZZZ  ' Shift along Z axis                  ' R  D=0.
>Guidance
To draw a volume shifted from its initial position when hidden
line removal is ON. It can be useful if you want to extract a
volume or some volumes from the detector to show them more clearly.
The last requested SHIFT for each volume
NAME is performed. Moreover, the SHIFT of
each volume will be performed starting from where its mother has
been shifted, so that it's easier to SHIFT nicely sets
of volumes using the mother-daughter relationships.
If '.' is given as the name of the volume
to be shifted, the shifts for all volumes will be reset.
The following commands will produce the translation along
the Z-axis of the previously drawn volume:
 EXAMPLE -
 dopt hide on
 satt * seen -2
 draw NAME 40 40 0 10 10 .01 .01
 shift NAME 0 0 10
 
>Action GXDRAW
 
>Command BOMB
>Parameters
BOOM  ' Exploding factor for volumes position ' R  D=0. R=-10.:10.
>Guidance
To 'explode' the detector. If BOOM is positive (values smaller
than 1. are suggested, but any value is possible)
all the volumes are shifted by a distance
proportional to BOOM along the direction between their centre
and the origin of the MARS; the volumes which are symmetric
with respect to this origin are simply not shown.
BOOM equal to 0 resets the normal mode.
A negative (greater than -1.) value of
BOOM will cause an 'implosion'; for even lower values of BOOM
the volumes' positions will be reflected respect to the origin.
This command can be useful to improve the 3D effect for very
complex detectors. The following commands will make explode the
detector:
 EXAMPLE -
 dopt hide on
 satt * seen 1
 draw NAME 40 40 0 10 10 .01 .01
 bomb 1
 next
 draw NAME 40 40 0 10 10 .01 .01
 
>Action GXDRAW
 
>Command DTREE
>Parameters
+
NAME   'Volume name' C D=' '
LEVMAX 'Depth level' I D=3 R=-15:15
ISELT  'Options    ' I D=111
>Guidance
This command allows the drawing of the logical tree,
displaying the name, the multiplicity and other information about the volumes,
via a call to GDTREE(name,levmax,isel):
if the third parameter is not given (default), the command will
produce the drawing of the tree displaying, for each volume, the
number of the following levels (red arrows) and of the preceeding
levels (green arrows); then the control is automatically given to the
mouse: clicking on the left button when the cursor is inside a volume's
pave will perform a DSPEC for that volume; doing the same when the cursor
is on a red arrow, will perform a DTREE for the relative volume (the
number of levels displayed depending on the clicked arrow); doing the
same for the 'i-th' green arrow of a given volume, will perform a DTREE
for its mother-volume staying 'i' levels before.
If running with X-windows, the drawing of the specification (DSPEC)
is performed
in a different window to speed up the scanning of the tree.
Iterating this procedure it is possible to analyse very easily and quickly
any kind of tree. Clicking the right button of the mouse will return
the control to the command mode.
If the ISELT parameter is given,
then the TREE will work as in the
previous version, with ISELT up to 10001.
The following command will perform a drawing of the tree and give the
control to the user via the mouse:
 EXAMPLE -
 dtree NAME 3
 
>Action GXDRAW
 
>Command DSPEC
>Parameters
NAME   'Volume name' C
>Guidance
Trough a call to GDSPEC(name), this command allows one to show three
views of the volume (two cut-views and a 3D view), together with
its geometrical specifications. The 3D drawing will
be performed according the current values of the options HIDE and
SHAD and according the current CVOL clipping parameters for that
volume.
>Action GXDRAW
 
>Command D3DSPEC
>Parameters
NAME   'Volume name' C
+
TETA3  'Theta angle' R D=40. R=0.:180.
PHI3   'Phi angle'   R D=40. R=0.:360.
PSI3   'Psi angle'   R D=0.  R=0.:360.
U03    'U-coord. (horizontal) of volume origin' R D=10. R=-40.:40.
V03    'V-coord. (vertical) of volume origin' R D=10. R=-40.:40.
ZM3    'Zoom factor for current size factors' R D=1. R=0.00001:10.
>Guidance
Trough a call to GSPE3D, this command allows one to show
the volume (3D views in real time), together with
its geometrical specifications (if using MOTIF). The 3D drawing will
be performed according the current values of the options HIDE and
SHAD and according the current CVOL clipping parameters for that
volume.
>Action GXDRAW
 
>Command DFSPC
>Parameters
NAME   'Volume name' C
+
CSORT  'Alphabetic sorting flag' C D='N' R='Y,N,0,1'
CINTER 'Interactive/Batch version' C D='I' R='I,B,0,1'
>Guidance
 CALL GDFSPC(name,isort,inter)
Same as DSPEC, but it will draw the specifications for all the volumes.
If the alphabetic sorting flag is YES, all pictures will be drawn in ascending
alphabetic order; isort is set to 1.
If INTERACTIVE, (inter=1), the routine will prompt the user at each plot
before doing a clear screen, otherwise it will clear automatically
the screen before starting a new frame.
>Action GXDRAW
 
>Command DTEXT
>Parameters
X0     'X-coord. (horizontal) of text string' R D=10. R=0.:20.
Y0     'Y-coord. (vertical) of text string' R D=10. R=0.:20.
TEXT   'Text string' C D='GEANT'
SIZE   'Character size (cm)' R D=.5
ANGLE  'Rotation angle (deg)' R D=0. R=0.:360.
LWID   'Line width' I D=4
CENT   'Centering option' C D='CENT' R='CENT,LEFT,RIGH'
>Guidance
 CALL GDRAWT(x0,y0,text,size,angle,lwid,opt)
It allows one to draw some text in the current picture.
Now more than 160 colours are available. The text colour
must be set via the command IGSET. The size of the
text will follow the zooming factors in the view banks.
>Action GXDRAW
 
>Command DVECTOR
>Parameters
XVECT  'Vector containing X-coord. (horizontal)' C
YVECT  'Vector containing Y-coord. (vertical)' C
NPOINT 'Number of coord.' I
>Guidance
Draw a polyline of 'npoint' point via
a call to GDRAWV(xvect,yvect,npoint)
where xvect and yvect are two KUIP vectors
>Action GXDRAW
 
>Command DSCALE
>Parameters
U      'U-coord. (horizontal) of the centre of scale' R
V      'V-coord. (vertical) of the centre of scale' R
>Guidance
 CALL GDSCAL(u,v)
It draws a scale centered in U,V.
>Action GXDRAW
 
>Command DAXIS
>Parameters
X0     'X-coord. of axis origin' R
Y0     'Y-coord. of axis origin' R
Z0     'Z-coord. of axis origin' R
DX     'Axis size' R
>Guidance
 CALL GDAXIS(x0,y0,z0,dx)
This commmand superimposes the axis of the MARS on the
current picture. It is useful for finding immediately the
orientation of the current drawing of the detector in the space.
>Action GXDRAW
 
>Command DMAN
>Parameters
U      'U-coord. (horizontal) of the centre of man' R
V      'V-coord. (vertical) of the centre of man' R
TYPE   'Man, Wm1, Wm2, Wm3' C D='MAN' R='MAN,WM1,WM2,WM3'
>Guidance
 CALL GDMAN(u,v),CALL GDWMN1(u,v),CALL GDWMN2(u,v),CALL GDWMN2(u,v)
It superimposes the picure of a man or of a woman, chosen among
three different ones, with the same scale factors as the detector
in the current drawing.
>Action GXDRAW
 
>Command DHEAD
>Parameters
+
ISEL   'Option flag' I D=111110
NAME   'Title' C D=' '
CHRSIZ 'Character size (cm) of title NAME' R D=0.6
>Guidance
 CALL GDHEAD(isel,name,chrsiz)
ISEL =
 0      to have only the header lines
 xxxxx1 to add the text name centered on top of header
 xxxx1x to add global detector name (first volume) on left
 xxx1xx to add date on right
 xx1xxx to select thick characters for text on top of header
 x1xxxx to add the text 'EVENT NR x' on top of header
 1xxxxx to add the text 'RUN NR x' on top of header
NOTE that ISEL=x1xxx1 or ISEL=1xxxx1 are illegal choices,
i.e. they generate overwritten text.
NAME is the title
and CHRSIZ the character size in cm of text name.
>Action GXDRAW
 
>Command MEASURE
>Guidance
Position the cursor on the first point (u1,v1) and hit the space bar(GKS).
Position the cursor on the second point (u2,v2) and hit the space bar(GKS).
Clicking the left button of the mouse (X11) will have the same effect as
hiting the space bar (GKS).
The command will compute and print the distance in space separating
the two points on the projection view. It can be useful to measure
distances either between volumes or between tracks or hits.
>Action GXDRAW
 
>Command PICK
>Parameters
>Guidance
Activates graphic input to identify detector elements
in a cut view. Clicking on the left button of the mouse when
the cursor is in a given point of the drawing and clicking again
(outside the detector) will produce the following effect:
a line joininig the two points will be drawn together with
the name and the medium number of the volume picked
with the first clicking close to the second point.
>Action GXPICK
 
>Command MOVE
>Parameters
NAME   'Volume name' C D='    '
+
NOPT   'S=sample mode,T=tracks,H=hits' C D='    '
>Guidance
Positioning some daughter volumes inside a 'mother', it can be
important to check if overlaps between such volumes have occurred.
Instead of putting the drawing in a view bank, zooming, and iterating
the process for different viewing angles of the same detector, the
MOVE facility has been developed (for machines running with X11):
it is sufficient to draw a view of the volumes to be analysed (after
setting the proper SEEN, COLO, etc. attributes) and then to enter
'MOVE' followed by the same 'NAME' used for the last command DRAW.
The detector will appear in a panel with five buttons at the
bottom: THETA, PHI, TRASL, ZOOM, OFF. Clicking on the left button
of the mouse, when the cursor is inside the THETA area, will rotate the
detector along the polar angle theta according to the
backward-to-forward movement of the mouse
(clicking up and down the left button if
not in sample mode); clicking on the right button of
the mouse will stop the rotation; clicking now on the
left button of the mouse when inside the PHI area will activate a
rotation along the polar angle phi. In the same way, activating the
TRASL button, the detector can be translated in the u,v plane
of the screen according to the 2D-movement of the mouse. Finally,
activating the ZOOM button, the detector will be zoomed (or unzoomed)
according to the backward-to-forward movement of the mouse. Clicking on the
OFF button will return the control to the 'command mode'. The MOVE
command will work also with hidden line removal and shading options
(when SHAD is on the background will be black);
moreover, if the volumes are clipped, exploded, shifted, etc., they
will be 'MOVED' with these features as well.
Tracks and hits of a previously stored physical event can be moved
together with the detector, allowing a dynamical 3-D analysis of the
simulated events. Clicking the central button of the mouse when a good
view of the event is found, will stop any movement and the mouse will
allow the normal picking capabilities first for the tracks and then for
the hits. After clicking of the right button, the normal
movement will restart to find another interesting view of the event
and to iterate the process.
The MOVE is also available in sample mode.
The following commands will produce a drawing of a volume
and then will give the control to the MOVE panel; try the following
possibilities:
 EXAMPLE 1 -
 dopt hide off
 satt * seen -2
 draw NAME 40 40 0 10 10 .01 .01
 move NAME
 EXAMPLE 2 -
 dopt hide on
 satt * seen -2
 draw NAME 40 40 0 10 10 .01 .01
 move NAME
 EXAMPLE 3 -
 dopt shad on
 satt * colo 3
 satt * fill 2
 dopt edge off
 draw NAME 40 40 0 10 10 .01 .01
 move NAME
 
>Action GXDRAW
 
>Command MOVE3D
>Parameters
NAME   'Volume name' C D='    '
+
THETA  'Viewing angle theta (for 3D projection)' R D=40. R=0.:180.
PHI    'Viewing angle phi (for 3D projection)' R D=40. R=0.:360.
PSI    'Viewing angle psi (for 2D rotation)' R D=0. R=0.:180.
U0     'U-coord. (horizontal) of volume origin' R D=10. R=0.:20.
V0     'V-coord. (vertical) of volume origin' R D=10. R=0.:20.
SU     'Scale factor for U-coord.' R D=0.01
SV     'Scale factor for V-coord.' R D=0.01
SZ     'Scale zoom factor' R D=1. R=0.1:10.
NOPT   'T=tracks,H=hits' C D='    ' R='T,H'
>Guidance
Same functionality of the command MOVE interfaced with MOTIF.
>Action GXDRAW
 
>Command PERSP
>Parameters
NAME   'Volume name' C D='    '
DISTT  'Volume distance from observer' R D=1000.
+
SAMP   'Control to the mouse' C D='OFF '
>Guidance
To control the perspective according to the variation of the distance
between the observer and the object (if PROJ has the value PERS).
If SAMP is ON the control of the distance is given via the mouse.
>Action GXDRAW
 
>Command LENS
>Parameters
KNUM   'View bank identifier' I D=1
+
KSAM   'Sample mode         ' C D='OFF '
>Guidance
Interactive zooming for detectors and events when running
with X-windows. Using this command, when showing the contents of a
view bank, it is possible to click (left button) in two points of the
drawing (which will represent the left upper corner and the right
bottom corner of the part to be zoomed). After the second click
a new 'window' will appear to fit the frame defined
by the two clicks and it will show a zoomed view as seen from a
lens with those dimensions. Clicking now the central button will
translate the lens over the drawing, while clicking the right button
will stop it. Moreover, clicking the left button of the
mouse, the lens will increase (or decrease) its magnification
power according to the backward-to-forward movement of the mouse.
A click on the right button will stop this action and it is possible
to restart the translation of the lens or, clicking
on the right button again, to make the lens disappear. It is then possible
to open another 'window-lens' with different dimensions. Thus,
this command can be useful to scan detailed parts of a detector or
to scan hits and showers for events. Clicking the right
button when no lens is displayed will return the control to the
'command mode'. The LENS is also available in sample mode when KSAM is
'ON'.
The following commands will fill a view bank and will
allow to scan the detector and an event previously stored
via the use of LENS (when running
with X-windows):
 EXAMPLE -
 satt * seen 1
 dopen 1
 draw NAME 40 40 0 10 10 .01 .01
 dxyz 0
 dhits * * 0 0 .2
 dclose
 dsh 1
 lens 1 on
 
>Action GXDRAW
 
>Command ZOOM
>Parameters
+
ZFU    'Zoom factor for U-coord. (horizontal)' R D=2.
ZFV    'Zoom factor for V-coord. (vertical)' R D=2.
ISEL   'Options' I D=1
UZ0    'U-coord. of the centre of zoom rectangle' R R=0.:20. D=10.
VZ0    'V-coord. of the centre of zoom rectangle' R R=0.:20. D=10.
U0     'U-coord. of the centre of resulting zoomed rectangle' R R=0.:20. D=10.
V0     'V-coord. of the centre of resulting zoomed rectangle' R R=0.:20. D=10.
>Guidance
 CALL GDZOOM(zfu,zfv,uz0,vz0,u0,v0)
This command sets the zoom parameters that will be used by
subsequent calls to the drawing routines. Each zoom operation is always
relative to the status of the current zoom parameters.
The scale factors in u,v are respectively  zfu,zfv.
zfu=0 (or zfv=0) will act as a reset (i.e. unzoomed viewing).
The zoom is computed around uz0,vz0 (user coordinates),
and the resulting picture will be centered at u0,v0.
The use of the space bar is replaced by the left button of the mouse
running with X11:
 
If isel=0 :
 1. position the cursor at (uz0,vz0)
 2. type the space bar (GKS)
(u0,v0 are chosen at centre of screen)
 
If isel=1 :
 1. position the cursor at first corner of zoom rectangle
 2. type the space bar (GKS)
 3. position the cursor at second corner of zoom rectangle
 4. type the space bar (GKS)
(zfu,zfv are chosen according to the zoom rectangle;
uz0,vz0 are chosen at the centre of the zoom rectangle;
u0,v0 are chosen at centre of screen)
 
If isel=2 :
 1. position the cursor at (uz0,vz0)
 2. type the space bar (GKS)
 3. position the cursor at (u0,v0)
 4. type the space bar (GKS)
 
If isel=1000+n and running with X-windows:
 1. n must be the identifier of an active view bank
 2. clicking on the left button of the mouse will display
    a zoomed view (computed around the cursor position) of
    the previous drawing in a new window
 3. it is now possible to iterate the zooming from the new window
 4. clicking on the right button will return the control to the
    main window
 5. clicking on the left button it is possible to open new windows
    zooming in other points of the detector
 6. clicking on the right button when the main window is active
    will return the control to the 'command mode'.
>Action GXDRAW
 
>Command DXYZ
>Parameters
+
ITRA   'Track number' I D=0
>Guidance
 CALL GDXYZ(itra)
Draw tracks previously stored via GSXYZ.
>Action GXDRAW
 
>Command KXYZ
>Parameters
+
EPSILO 'Delta angle' R D=0.25
>Guidance
 CALL GKXYZ(epsilo)
The picking of track points requires the JXYZ data structure
and is  repeated until the character typed is 'Q' or 'q' (GKS)
or the right button of the mouse is clicked (X11).
EPSILO is the delta angle used for picking; if EPSILO=0
there is no optimization performed and
over all the track points the one nearest to the pick
point is taken.
>Action GXDRAW
 
>Command DPART
>Parameters
+
ITRA   'Track number' I D=0
ISEL   'Option flag' I D=11
SIZE   'Character size (cm) for particle names' R D=0.25
>Guidance
 CALL GDPART(itra,isel,size)
 isel=x1 to draw the track number
 isel=1x to draw the particle name
>Action GXDRAW
 
>Command DHITS
>Parameters
+
CHUSET  'User set identifier' C D='*'
CHUDET  'User detector identifier' C D='*'
ITRA   'Number of the selected track' I D=0
ISYMB  'Character selection number' I D=0
SSYMB  'Size of characters (cm)' R D=0.1
>Guidance
CALL GDHITS(chuset,chudet,itra,isymb,ssymb).
The character plotted at each hit point may be chosen by isymb :
      -1   (small) hardware points             (fast)
       0   software crosses                    (default)
   840,850   empty/full circles                  (slow)
   841,851   empty/full squares                  (slow)
   842,852   empty/full triangles (up)           (slow)
   843,853   empty diamond/full triangle (down)  (slow)
   844,854   empty/full stars                    (slow)
Except for isymb=-1, the size of the character on the screen can be
chosen by SSYMB cm. The hit colour will follow the value of TXCI (text
colour) for isymb>0, the value of PMCI (polymarkers colour) for isymb<0,
the value of PLCI (polyline colour) for isymb=0.
>Action GXDRAW
 
>Command KHITS
>Parameters
+
CHUSET  'User set identifier' C D='*'
CHUDET  'User detector identifier' C D='*'
EPSILO 'Pick aperture' R D=0.1
>Guidance
 CALL GKHITS(chuset,chudet,epsilo)
The picking of hit points requires the appropriate JSET data structure
have been filled
and is  repeated until the character typed is 'Q' or 'q' (GKS) or the
right button of the mouse is clicked (X11).
If the character typed to pick is 'K' or 'k' then the
kinematics of the corresponding track is also printed.
The search is made of all the hits of all tracks in
detector CHUDET of set CHUSET.
EPSILO is the pick aperture; if EPSILO<0 its absolute value is taken
and in addition the pick aperture is drawn; if EPSILO=0
there is an infinite pick aperture and
over all the hits the one nearest to the pick point is taken.
>Action GXDRAW
 
>Command DCHIT
>Parameters
+
CHUSET  'User set identifier' C D='*'
CHUDET  'User detector identifier' C D='*'
ITRA   'Number of the selected track' I D=0
ISYMB  'Character selection number' I D=0
SIZMAX 'Maximum character size (cm)' R D=1
IHIT   'Index of array HITS' I D=4
HITMIN 'Lower boundary of HITS(IHIT)' R D=0
HITMAX 'Upper boundary of HITS(IHIT)' R D=0
>Guidance
 CALL GDCHIT(chuset,chudet,itra,isymb,sizmax,ihit,hitmin,hitmax)
The character plotted at each hit point may be chosen via
CSYMB; isymb is composed as:
      -1   (small) hardware points             (fast)
       0   software crosses                    (default)
 840,850   empty/full circles                  (slow)
 841,851   empty/full squares                  (slow)
 842,852   empty/full triangles (up)           (slow)
 843,853   empty diamond/full triangle (down)  (slow)
 844,854   empty/full stars                    (slow)
Except for isymb=-1 the SIZE of the character on the screen
is a function of HITS(IHIT), the array containing the calorimeter
quantity, with HITMIN and HITMAX defining its range.
The maximum character size (used in overflow) is SIZMAX.
 SIZE = SIZMAX * ( HITS(IHIT) - HITMIN ) / HITMAX
>Action GXDRAW
 
>Command DUVIEW
>Parameters
NAME   'Detector name' C
TYPE   'View name' C
CPXTYP 'Complexity name' C
+
IVIEW  'View number where picture is stored' I D=0
>Guidance
 CALL GUVIEW(name,type,cpxtyp,iview)
>Action GXDRAW
 
>Name GKGCON
 
>Menu /GEANT/GRAPHICS_CONTROL
>Guidance
Graphics control commands.
 
>Command DOPEN
>Parameters
IVIEW  'View number' I
>Guidance
 CALL GDOPEN(iview)
When a drawing is very complex and requires a long time to be
executed, it can be useful to store it in a view bank: after a
call to DOPEN and the execution of the drawing (nothing will
appear on the screen), and after a necessary call to DCLOSE,
the contents of the bank can be displayed in a very fast way
through a call to DSHOW; therefore, the detector can be easily
zoomed many times in different ways. Please note that the pictures
with solid colours can now be stored in a view bank or in 'PICTURE FILES'.
>Action GXGCON
 
>Command DSHOW
>Parameters
+
IVIEW  'View number' I
>Guidance
 CALL GDSHOW(iview)
It shows on the screen the contents of a view bank. It
can be called after a view bank has been closed.
>Action GXGCON
 
>Command DELETE
>Parameters
IVIEW  'View number' I
>Guidance
 CALL GDELET(iview)
It deletes a view bank from memory.
>Action GXGCON
 
>Command DCLOSE
>Guidance
 CALL GDCLOS
It closes the currently open view bank; it must be called after the
end of the drawing to be stored.
>Action GXGCON
 
>Command CHANGEWK
>Guidance
CALL GCHNWK
It open a new workstation (if not already opened) and activate it
(deactivating the default one).
>Action GXGCON
 
>Command RESETWK
>Guidance
CALL GRESWK
It deactivate the previosly activated workstation and reactivate the
default one.
>Action GXGCON
 
>Command SSETAT
>Parameters
IOPT   'Attribute name' C
>Guidance
Set current attribute.
>Action GXGCON
 
>Command SSETVA
>Parameters
+
RVAL   'Attribute value' R D=1. R=-10.:10.
>Guidance
Set current attribute value.
>Action GXGCON
 
>Command SATT
>Parameters
+
NAME   'Volume name' C D='*   '
IOPT   'Name of the attribute to be set' C D='DEFA'
IVAL   'Value to which the attribute is to be set' I D=10000
>Guidance
 CALL GSATT(name,iopt,ival)
name='*' stands for all the volumes.
iopt can be chosen among the following :
 
 'WORK'   0=volume name is inactive for the tracking
          1=volume name is active for the tracking (default)
 
 'SEEN'   0=volume name is invisible
          1=volume name is visible (default)
         -1=volume invisible with all its descendants in the tree
         -2=volume visible but not its descendants in the tree
 
 'LSTY'   line style 1,2,3,... (default=1)
          LSTY=7 will produce a very precise approximation for
          revolution bodies.
 
 'LWID'   line width -7,...,1,2,3,..7 (default=1)
          LWID<0 will act as abs(LWID) was set for the volume
          and for all the levels below it. When SHAD is 'ON', LWID
          represent the linewidth of the scan lines filling the surfaces
          (whereas the FILL value represent their number). Therefore
          tuning this parameter will help to obtain the desired
          quality/performance ratio.
 
 'COLO'   colour code -166,...,1,2,..166 (default=1)
          n=1=black
          n=2=red;    n=17+m, m=0,25, increasing luminosity according to 'm';
          n=3=green;  n=67+m, m=0,25, increasing luminosity according to 'm';
          n=4=blue;   n=117+m, m=0,25, increasing luminosity according to 'm';
          n=5=yellow; n=42+m, m=0,25, increasing luminosity according to 'm';
          n=6=violet; n=142+m, m=0,25, increasing luminosity according to 'm';
          n=7=lightblue; n=92+m, m=0,25, increasing luminosity according to 'm';
          colour=n*10+m, m=1,2,...9, will produce the same colour
          as 'n', but with increasing luminosity according to 'm';
          COLO<0 will act as if abs(COLO) was set for the volume
          and for all the levels below it.
          When for a volume the attribute FILL is > 1 (and the
          option SHAD is on), the ABS of its colour code must be < 8
          because an automatic shading of its faces will be
          performed.
 
 'FILL'   (1992) fill area  -7,...,0,1,...7 (default=0)
          when option SHAD is 'on' the FILL attribute of any
          volume can be set different from 0 (normal drawing);
          if it is set to 1, the faces of such volume will be filled
          with solid colours; if ABS(FILL) is > 1, then a light
          source is placed along the observer line, and the faces of
          such volumes will be painted by colours whose luminosity
          will depend on the amount of light reflected;
          if ABS(FILL) = 1, then it is possible to use all the 166
          colours of the colour table, becouse the automatic shading
          is not performed;
          for increasing values of FILL the drawing will be performed
          with higher and higher resolution improving the quality (the
          number of scan lines used to fill the faces increases with FILL);
          it is possible to set different values of FILL
          for different volumes, in order to optimize at the same time
          the performance and the quality of the picture;
          FILL<0 will act as if abs(FILL) was set for the volume
          and for all the levels below it.
          This kind of drawing can be saved in 'picture files'
          or in view banks.
          0=drawing without fill area
          1=faces filled with solid colours and resolution = 6
          2=lowest resolution (very fast)
          3=default resolution
          4=.................
          5=.................
          6=.................
          7=max resolution
          Finally, if a coloured background is desired, the FILL
          attribute for the first volume of the tree must be set
          equal to -abs(colo), colo being >0 and <166.
 
 'SET '   set number associated to volume name
 'DET '   detector number associated to volume name
 'DTYP'   detector type (1,2)
>Action GXGCON
 
>Command SCALE
>Parameters
GSCU   'Scale factor for U-coord.' R
GSCV   'Scale factor for V-coord.' R
>Guidance
Change the scale factors GSCU and GSCV in /GCDRAW/.
>Action GXGCON
 
>Command COLOR
>Parameters
ICOL   'Colour code' I D=1
>Guidance
 CALL GDCOL(-abs(icol))
>Action GXGCON
 
>Command LWID
>Parameters
LWIDTH 'Line width code' I D=1
>Guidance
 CALL GDLW(-abs(lwidth))
>Action GXGCON
 
>Command NEXT
>Guidance
Clear screen (start a new picture on graphics file, if opened).
>Action GXGCON
 
>Command DOPT
>Parameters
+
IOPT   'Option name' C D='*'
IVAL   'Option value' C D='*'
>Guidance
 CALL GDOPT(iopt,ival)
To set/modify the drawing options.
   IOPT   IVAL      Action
 
   THRZ    ON       Draw tracks in R vs Z
           OFF (D)  Draw tracks in X,Y,Z
           180
           360
   PROJ    PARA (D) Parallel projection
           PERS     Perspective
   TRAK    LINE (D) Trajectory drawn with lines
           POIN       " " with markers
   HIDE    ON       Hidden line removal using the CG package
           OFF (D)  No hidden line removal
   SHAD    ON       Fill area and shading of surfaces.
           OFF (D)  Normal hidden line removal.
   RAYT    ON       Ray-tracing on.
           OFF (D)  Ray-tracing off.
   EDGE    OFF      Does not draw contours when shad is on.
           ON  (D)  Normal shading.
   MAPP    1,2,3,4  Mapping before ray-tracing.
           0   (D)  No mapping.
   USER    ON       User graphics options in the raytracing.
           OFF (D)  Automatic graphics options.
>Action GXGCON
 
 
>Command SIZE
>Parameters
+
XSIZE 'Size along X' R D=20.
YSIZE 'Size along Y' R D=20.
>Guidance
Set the size of the picture.
On the terminal, the pictures will have the ratio YSIZE/XSIZE, and,
if a metafile is produced, pictures will be YSIZE by XSIZE cm.
This command sets the parameters for the normalisation transformation
number 1 to [0-XSIZE], [0-YSIZE].
>Action GXGCON
 
>Command SPERS
>Parameters
DPERS  'Distance from the origin' R
>Guidance
Set the variable dpers in /GCDRAW/, representing
the distance from the origin when using option PERSpective.
>Action GXGCON
 
>Command MAP_COLOR
>Parameters
+
ICADD  'Colour table index' I D=0
ICVAL  'Colour table value' I D=0
>Guidance
Sets the color table LOOKTB(ICADD)=ICVAL.
If ICADD=0 then LOOKTB(1:16) is taken.
If ICVAL is omitted the current value of LOOKTB(ICADD) is shown.
>Action GXGCON
 
>Name GKLIST
>Menu /GEANT/LISTS
>Guidance
 
 
>Command HSTA
>Parameters
+
LHSTA_1  'user word' C
LHSTA_2  'user word' C
LHSTA_3  'user word' C
LHSTA_4  'user word' C
LHSTA_5  'user word' C
LHSTA_6  'user word' C
LHSTA_7  'user word' C
LHSTA_8  'user word' C
LHSTA_9  'user word' C
LHSTA_10  'user word' C
LHSTA_11  'user word' C
LHSTA_12  'user word' C
LHSTA_13  'user word' C
LHSTA_14  'user word' C
LHSTA_15  'user word' C
LHSTA_16  'user word' C
LHSTA_17  'user word' C
LHSTA_18  'user word' C
LHSTA_19  'user word' C
LHSTA_20  'user word' C
>Guidance
The command HSTA is similar to the HSTA data records. It can accept
up to 20 4-character words. If the first argument is '.', the number
of words is reset to 0 and all the words to four blanks.
>Action GXLIST
 
>Command GET
>Parameters
+
LGET_1  'user word' C
LGET_2  'user word' C
LGET_3  'user word' C
LGET_4  'user word' C
LGET_5  'user word' C
LGET_6  'user word' C
LGET_7  'user word' C
LGET_8  'user word' C
LGET_9  'user word' C
LGET_10  'user word' C
LGET_11  'user word' C
LGET_12  'user word' C
LGET_13  'user word' C
LGET_14  'user word' C
LGET_15  'user word' C
LGET_16  'user word' C
LGET_17  'user word' C
LGET_18  'user word' C
LGET_19  'user word' C
LGET_20  'user word' C
>Guidance
The command GET is similar to the GET data records. It can accept
up to 20 4-character words. If the first argument is '.', the number
of words is reset to 0 and all the words to four blanks.
>Action GXLIST
 
>Command SAVE
>Parameters
+
LSAVE_1  'user word' C
LSAVE_2  'user word' C
LSAVE_3  'user word' C
LSAVE_4  'user word' C
LSAVE_5  'user word' C
LSAVE_6  'user word' C
LSAVE_7  'user word' C
LSAVE_8  'user word' C
LSAVE_9  'user word' C
LSAVE_10  'user word' C
LSAVE_11  'user word' C
LSAVE_12  'user word' C
LSAVE_13  'user word' C
LSAVE_14  'user word' C
LSAVE_15  'user word' C
LSAVE_16  'user word' C
LSAVE_17  'user word' C
LSAVE_18  'user word' C
LSAVE_19  'user word' C
LSAVE_20  'user word' C
>Guidance
The command SAVE is similar to the SAVE data records. It can accept
up to 20 4-character words. If the first argument is '.', the number
of words is reset to 0 and all the words to four blanks.
>Action GXLIST
 
>Command SETS
>Parameters
+
LSETS_1  'user word' C
LSETS_2  'user word' C
LSETS_3  'user word' C
LSETS_4  'user word' C
LSETS_5  'user word' C
LSETS_6  'user word' C
LSETS_7  'user word' C
LSETS_8  'user word' C
LSETS_9  'user word' C
LSETS_10  'user word' C
LSETS_11  'user word' C
LSETS_12  'user word' C
LSETS_13  'user word' C
LSETS_14  'user word' C
LSETS_15  'user word' C
LSETS_16  'user word' C
LSETS_17  'user word' C
LSETS_18  'user word' C
LSETS_19  'user word' C
LSETS_20  'user word' C
>Guidance
The command SETS is similar to the SETS data records. It can accept
up to 20 4-character words. If the first argument is '.', the number
of words is reset to 0 and all the words to four blanks.
>Action GXLIST
 
>Command LPRIN
>Parameters
+
LPRIN_1  'user word' C
LPRIN_2  'user word' C
LPRIN_3  'user word' C
LPRIN_4  'user word' C
LPRIN_5  'user word' C
LPRIN_6  'user word' C
LPRIN_7  'user word' C
LPRIN_8  'user word' C
LPRIN_9  'user word' C
LPRIN_10  'user word' C
LPRIN_11  'user word' C
LPRIN_12  'user word' C
LPRIN_13  'user word' C
LPRIN_14  'user word' C
LPRIN_15  'user word' C
LPRIN_16  'user word' C
LPRIN_17  'user word' C
LPRIN_18  'user word' C
LPRIN_19  'user word' C
LPRIN_20  'user word' C
>Guidance
The command PRIN is similar to the PRIN data records. It can accept
up to 20 4-character words. If the first argument is '.', the number
of words is reset to 0 and all the words to four blanks.
>Action GXLIST
 
>Command GEOM
>Parameters
+
LGEOM_1  'user word' C
LGEOM_2  'user word' C
LGEOM_3  'user word' C
LGEOM_4  'user word' C
LGEOM_5  'user word' C
LGEOM_6  'user word' C
LGEOM_7  'user word' C
LGEOM_8  'user word' C
LGEOM_9  'user word' C
LGEOM_10  'user word' C
LGEOM_11  'user word' C
LGEOM_12  'user word' C
LGEOM_13  'user word' C
LGEOM_14  'user word' C
LGEOM_15  'user word' C
LGEOM_16  'user word' C
LGEOM_17  'user word' C
LGEOM_18  'user word' C
LGEOM_19  'user word' C
LGEOM_20  'user word' C
>Guidance
The command GEOM is similar to the GEOM data records. It can accept
up to 20 4-character words. If the first argument is '.', the number
of words is reset to 0 and all the words to four blanks.
>Action GXLIST
 
>Command VIEW
>Parameters
+
LVIEW_1  'user word' C
LVIEW_2  'user word' C
LVIEW_3  'user word' C
LVIEW_4  'user word' C
LVIEW_5  'user word' C
LVIEW_6  'user word' C
LVIEW_7  'user word' C
LVIEW_8  'user word' C
LVIEW_9  'user word' C
LVIEW_10  'user word' C
LVIEW_11  'user word' C
LVIEW_12  'user word' C
LVIEW_13  'user word' C
LVIEW_14  'user word' C
LVIEW_15  'user word' C
LVIEW_16  'user word' C
LVIEW_17  'user word' C
LVIEW_18  'user word' C
LVIEW_19  'user word' C
LVIEW_20  'user word' C
>Guidance
The command VIEW is similar to the VIEW data records. It can accept
up to 20 4-character words. If the first argument is '.', the number
of words is reset to 0 and all the words to four blanks.
>Action GXLIST
 
>Command PLOT
>Parameters
+
LPLOT_1  'user word' C
LPLOT_2  'user word' C
LPLOT_3  'user word' C
LPLOT_4  'user word' C
LPLOT_5  'user word' C
LPLOT_6  'user word' C
LPLOT_7  'user word' C
LPLOT_8  'user word' C
LPLOT_9  'user word' C
LPLOT_10  'user word' C
LPLOT_11  'user word' C
LPLOT_12  'user word' C
LPLOT_13  'user word' C
LPLOT_14  'user word' C
LPLOT_15  'user word' C
LPLOT_16  'user word' C
LPLOT_17  'user word' C
LPLOT_18  'user word' C
LPLOT_19  'user word' C
LPLOT_20  'user word' C
>Guidance
The command PLOT is similar to the PLOT data records. It can accept
up to 20 4-character words. If the first argument is '.', the number
of words is reset to 0 and all the words to four blanks.
>Action GXLIST
 
>Command STAT
>Parameters
+
LSTAT_1  'user word' C
LSTAT_2  'user word' C
LSTAT_3  'user word' C
LSTAT_4  'user word' C
LSTAT_5  'user word' C
LSTAT_6  'user word' C
LSTAT_7  'user word' C
LSTAT_8  'user word' C
LSTAT_9  'user word' C
LSTAT_10  'user word' C
LSTAT_11  'user word' C
LSTAT_12  'user word' C
LSTAT_13  'user word' C
LSTAT_14  'user word' C
LSTAT_15  'user word' C
LSTAT_16  'user word' C
LSTAT_17  'user word' C
LSTAT_18  'user word' C
LSTAT_19  'user word' C
LSTAT_20  'user word' C
>Guidance
The command STAT is similar to the STAT data records. It can accept
up to 20 4-character words. If the first argument is '.', the number
of words is reset to 0 and all the words to four blanks.
>Action GXLIST
 
>Command RGET
>Parameters
+
LRGET_1  'user word' C
LRGET_2  'user word' C
LRGET_3  'user word' C
LRGET_4  'user word' C
LRGET_5  'user word' C
LRGET_6  'user word' C
LRGET_7  'user word' C
LRGET_8  'user word' C
LRGET_9  'user word' C
LRGET_10  'user word' C
LRGET_11  'user word' C
LRGET_12  'user word' C
LRGET_13  'user word' C
LRGET_14  'user word' C
LRGET_15  'user word' C
LRGET_16  'user word' C
LRGET_17  'user word' C
LRGET_18  'user word' C
LRGET_19  'user word' C
LRGET_20  'user word' C
>Guidance
The command RGET is similar to the RGET data records. It can accept
up to 20 4-character words. If the first argument is '.', the number
of words is reset to 0 and all the words to four blanks.
>Action GXLIST
 
>Command RSAV
>Parameters
+
LRSAVE_1  'user word' C
LRSAVE_2  'user word' C
LRSAVE_3  'user word' C
LRSAVE_4  'user word' C
LRSAVE_5  'user word' C
LRSAVE_6  'user word' C
LRSAVE_7  'user word' C
LRSAVE_8  'user word' C
LRSAVE_9  'user word' C
LRSAVE_10  'user word' C
LRSAVE_11  'user word' C
LRSAVE_12  'user word' C
LRSAVE_13  'user word' C
LRSAVE_14  'user word' C
LRSAVE_15  'user word' C
LRSAVE_16  'user word' C
LRSAVE_17  'user word' C
LRSAVE_18  'user word' C
LRSAVE_19  'user word' C
LRSAVE_20  'user word' C
>Guidance
The command RSAV is similar to the RSAV data records. It can accept
up to 20 4-character words. If the first argument is '.', the number
of words is reset to 0 and all the words to four blanks.
>Action GXLIST
 
>Name GKGEOM
>Menu /GEANT/GEOMETRY
>Guidance
Geometry commands.
 
>Command OPTI
>Parameters
IOPTI  'GSORD optimisation level' I D=0 R=-1,2
>Guidance
This flag controls the tracking optimisation performed via the
GSORD routine:
    1 no optimisation at all; GSORD calls disabled;
    0 no optimisation; only user calls to GSORD kept;
    1 all non-GSORDered volumes are ordered along the best axis;
    2 all volumes are ordered along the best axis.
>Action GXGEOM
 
>Command SVOLU
>Parameters
NAME   'Volume name' C
SHAPE  'Volume type' C
NUMED  'Tracking medium number' I
NPAR   'Number of shape parameters' I
PAR    'Vector containing shape parameters' C
>Guidance
 CALL GSVOLU(name,shape,numed,par,npar,ivolu)
where par is a KUIP vector.
It creates a new volume in the JVOLUM data structure.
>Action GXGEOM
 
>Command SPOS
>Parameters
NAME   'Volume name' C
NUMBER 'Copy number of the volume' I
MOTHER 'Mother volume name' C
X0     'X coord. of the volume in mother ref. sys.' R
Y0     'Y coord. of the volume in mother ref. sys.' R
Z0     'Z coord. of the volume in mother ref. sys.' R
IROT   'Rotation matrix number w.r.t. mother ref. sys.' I
ONLY   'ONLY/MANY flag' C
>Guidance
 CALL GSPOS(name,number,mother,x0,y0,z0,irot,only)
It positions a previously defined volume in the mother.
>Action GXGEOM
 
>Command SDVN
>Parameters
NAME   'Volume name' C
MOTHER 'Mother volume name' C
NDIV   'Number of divisions' I
CAXIS  'Axis value' C R='X,Y,Z,1,2,3'
>Guidance
 CALL GSDVN(name,mother,ndiv,iaxis)
X,Y,Z of CAXIS will be translated to 1,2,3 for IAXIS.
It divides a previously defined volume.
>Action GXGEOM
 
>Command PVOLU
>Parameters
NUMB   'Volume ID' I
>Guidance
 CALL GPVOLU(numb)
Prints volumes' specifications.
>Action GXGEOM
 
>Command SROTM
>Parameters
IROT   'Rotation matrix number' I
THETA1 'Polar angle for axis I'        R D=0. R=0.:180.
PHI1   'Azimuthal angle for axis I'    R D=0. R=0.:360.
THETA2 'Polar angle for axis II'       R D=0. R=0.:180.
PHI2   'Azimuthal angle for axis II'   R D=0. R=0.:360.
THETA3 'Polar angle for axis III'      R D=0. R=0.:180.
PHI3   'Azimuthal angle for axis III'  R D=0. R=0.:360.
>Guidance
 CALL GSROTM(irot,theta1,phi1,theta2,phi2,theta3,phi3)
It defines the rotation matrix number IROT.
>Action GXGEOM
 
>Command PROTM
>Parameters
NUMB   'Matrix ID' I
>Guidance
 CALL GPROTM(numb)
Print matrixes' specifications.
>Action GXGEOM
 
 
>Command STMED
>Parameters
NTMED  'Tracking medium number' I D=1
NAME   'Tracking medium name' C
NMAT   'Material number' I D=1
ISVOL  'Sensitive volume flag' I D=0
IFIELD 'Magnetic field' I D=0
FIELDM 'Max. field value (Kilogauss)' R D=0
TMAXFD 'Max. angle due to field (deg/step)' R D=0.01
STEMAX 'Max. step allowed'  R D=1.E+10
DEEMAX 'Max. fraction of energy lost in a step' R D=0.01
EPSIL  'Tracking precision (cm)' R D=0.01
STMIN  'Min. step due to continuos processes (cm)' R D=0.1
>Guidance
      CALL GSTMED(ntmed,name,nmat,isvol,ifield,fieldm,tmaxfd,
     +            stemax,deemax,epsil,stmin,0,0)
IFIELD = 0 if no magnetic field; IFIELD = -1 if user decision in GUSWIM;
IFIELD = 1 if tracking performed with GRKUTA; IFIELD = 2 if tracking
performed with GHELIX; IFIELD = 3 if tracking performed with GHELX3.
>Action GXGEOM
 
>Command PTMED
>Parameters
NUMB   'Medium ID' I
>Guidance
 CALL GPTMED(numb)
Print tracking media's specifications.
>Action GXGEOM
 
>Command EDITV
>Parameters
+
ISEL   'Options' I D=0
NAME   'Volume name' C D='   '
>Guidance
 CALL GEDITV(isel,name)
When the routine prompts for input parameters that do not need
to be changed, type return.
ISEL is used to select the editing operation to be performed:
 ISEL=0, CALL GGCLOS
 ISEL=1, to modify shape parameters PAR given by GSVOLU
 ISEL=2, to modify NAME given by GSVOLU
 ISEL=3, to delete NAME given by GSVOLU
 ISEL=4, to unlink NAME,NR given by GSPOS/GSDVN/GSDV..
 ISEL=5, to modify X0,Y0,Z0 of NAME,NR given by GSPOS
 ISEL=6, to modify IROT of NAME,NR given by GSPOS
 ISEL=7, to modify NDIV given by GSDVN
 ISEL=8, to modify IAXIS given by GSDVN
>Action GXGEOM
 
>Command CADINT
>Parameters
FNAME   'Name of the SET file'                        C D='example.set'
ANAME   'Name of the volume'                          C
NBINS   'Number of the instances'                     I D=1
LUNIT   'Logical unit number for SET file'            I D=66
LUNIT2  'Logical unit number for material file'       I D=67
INST    'Name of your institute'                      C D='CERN'
SITE    'Name of site'                                C D='MEYRIN'
DEPT    'Name of departement'                         C D='CN'
RESP    'Name of sender'                              C D='god_knows_who'
>Guidance
 CALL GTXSET(fname,aname,nbins,lunit,lunit2,inst,site,dept,resp)
This command produces a SET file describing the given volume with
the contents currently set visible. (Use the visibility attribute,
see SATT SEEN.) The description is given as a flat assembly
related to the global coordinate system.
The ouput can be read into CAD systems (EUCLID-IS) trough a SET interface.
A list of materials of the volumes in the SET file and the GEANT tree
is written into a file with the same filename as the SET file,
but with extension .mat.
>Action GXGEOM
 
>Command REUCLID
>Parameters
LUN    'Logical unit of the file to be read'        I R=1:100
FNAME  'Name of the EUCLID file to be read'         C
>Guidance
          CALL GREUCL(LUN,FNAME)
Calls the routine to read into GEANT a geometry from an ASCII file
written by the EUCLID-GEANT interface.
>Action GXGEOM
 
>Command WEUCLID
>Parameters
LUN    'Logical unit of the file to be written'       I R=1:100
FNAME  'Name of the EUCLID file to be written'        C
TOPVOL 'Volume name of the starting node'             C
+
NUMBER 'Copy number of TOPVOL (relevant for GSPOSP)'  I D=1
NLEVEL 'Number of levels in the tree structure'       I D=15
 
>Guidance
          CALL GWEUCL(LUN,FNAME)
Calls the routine to write the current GEANT geometry into an ASCII file
in EUCLID compatible format.
>Action GXGEOM
 
>Menu /GEANT/CREATE
>Guidance
It creates volumes of the given shape interactively.
CALL GSVOLU(name,shape,numed,par,npar,ivolu)
where par is a KUIP vector
 
>Command SBOX
>Parameters
NAME   'Volume name' C
NUMED  'Tracking medium number' I
HALFX  'Half X length' R
HALFY  'Half Y length' R
HALFZ  'Half Z length' R
+
YESNO  'GSPOSP option' C D='NO' R='YES,NO'
>Guidance
>Action GXGEOM
 
>Command STRD1
>Parameters
NAME   'Volume name' C
NUMED  'Tracking medium number' I
HLFDWX 'Half X length in Lower Z Surface' R
HLFUPX 'Half X length in Upper Z Surface' R
HALFY  'Half Y length' R
HALFZ  'Half Z length' R
+
YESNO  'GSPOSP option' C D='NO' R='YES,NO'
>Guidance
>Action GXGEOM
 
 
>Command STRD2
>Parameters
NAME    'Volume name' C
NUMED   'Tracking medium number' I
HLFDWX  'Half X length in Lower Z Surface' R
HLFUPX  'Half X length in Upper Z Surface' R
HLFDWY  'Half Y length in Lower Z Surface' R
HLFUPY  'Half Y length in Upper Z Surface' R
HALFZ   'Half Z length' R
+
YESNO   'GSPOSP option' C D='NO' R='YES,NO'
>Guidance
>Action GXGEOM
 
 
>Command STUBE
>Parameters
NAME   'Volume name' C
NUMED  'Tracking medium number' I
INRAD  'Inside Radius' R
OUTRAD 'Outside Radius' R
HALFZ  'Half Z length' R
+
YESNO  'GSPOSP option' C D='NO' R='YES,NO'
>Guidance
>Action GXGEOM
 
 
>Command STUBS
>Parameters
NAME   'Volume name' C
NUMED  'Tracking medium number' I
INRAD  'Inside Radius' R
OUTRAD 'Outside Radius' R
HALFZ  'Half Z length' R
SPHI   'Start of section PHI' R R=0.:360.
EPHI   'End of section PHI' R R=0.:360.
+
YESNO  'GSPOSP option' C D='NO' R='YES,NO'
>Guidance
>Action GXGEOM
 
 
>Command SCONE
>Parameters
NAME   'Volume name' C
NUMED  'Tracking medium number' I
INRDW  'Inside Radius in Lower Z Surface' R
OUTRDW 'Outside Radius in Lower Z Surface' R
INRUP  'Inside Radius in Upper Z Surface' R
OUTRUP 'Outside Radius in Upper Z Surface' R
HALFZ  'Half Z length' R
+
YESNO  'GSPOSP option' C D='NO' R='YES,NO'
>Guidance
>Action GXGEOM
 
 
>Command SCONS
>Parameters
NAME   'Volume name' C
NUMED  'Tracking medium number' I
INRDW  'Inside Radius in Lower Z Surface' R
OUTRDW 'Outside Radius in Lower Z Surface' R
INRUP  'Inside Radius in Upper Z Surface' R
OUTRUP 'Outside Radius in Upper Z Surface' R
HALFZ  'Half Z length' R
SPHI   'Start of section PHI' R R=0.:360.
EPHI   'End of section PHI' R R=0.:360.
+
YESNO  'GSPOSP option' C D='NO' R='YES,NO'
>Guidance
>Action GXGEOM
 
 
>Command SSPHE
>Parameters
NAME   'Volume name' C
NUMED  'Tracking medium number' I
INRAD  'Inside Radius' R
OUTRAD 'Outside Radius' R
SPHI   'Start of section PHI' R R=0.:360.
EPHI   'End of section PHI' R R=0.:360.
STHETA 'Start of section THETA' R
ETHETA 'End of section THETA' R
+
YESNO  'GSPOSP option' C D='NO' R='YES,NO'
>Guidance
>Action GXGEOM
 
 
>Command SPARA
>Parameters
NAME   'Volume name' C
NUMED  'Tracking medium number' I
HALFX  'Half X length' R
HALFY  'Half Y length' R
HALFZ  'Half Z length' R
AXIS   'Angle of Y mid-faces segment to Y axis' R R=0.:360.
PHI    'PHI angle of Low Z mid-face to High Z mid-face segment' R R=0.:360.
THETA  'THETA angle of mid-low-Z-face to mid-high-Z-face segment' R R=0.:360.
+
YESNO  'GSPOSP option' C D='NO' R='YES,NO'
>Guidance
>Action GXGEOM
 
 
>Name GKCONT
 
>Menu /GEANT/CONTROL
>Guidance
Control commands.
 
>Command KINE
>Parameters
IKINE     'IKINE'     I D=1
+
PKINE1  'PKINE(1)'  R
PKINE2  'PKINE(2)'  R
PKINE3  'PKINE(3)'  R
PKINE4  'PKINE(4)'  R
PKINE5  'PKINE(5)'  R
PKINE6  'PKINE(6)'  R
PKINE7  'PKINE(7)'  R
PKINE8  'PKINE(8)'  R
PKINE9  'PKINE(9)'  R
PKINE10 'PKINE(10)' R
>Guidance
Set the variables in /GCFLAG/ IKINE, PKINE(10)
>Action GXCONT
 
>Command RUNG
>Parameters
IDRUN  'User run number'             I
IDEVT  'User starting event number'  I
>Guidance
Set the run number and the starting value for the user event number.
>Action GXCONT
 
>Command SORD
>Parameters
ISTORD  'Flag to control user ordering of the stack' I D=1 R=1,0
>Guidance
If ISTORD is set to 1, the particle with the highest value of the
user weight UPWGHT will be selected to be tracked next.
>Action GXCONT
 
>Command GTIME
>Parameters
TIMINT  'Total time after initialisation'          R
TIMEND  'Time reserved for the termination phase'  R
ITIME   'Frequency of control printing'            I
>Guidance
These commands have limited use in the interactive version. In
particular the value of TIMINT is disregarded by GEANT.
>Action GXCONT
 
>Command TRACK
>Guidance
Restart tracking, clearing the track and hit
banks, but keeping the kinematics.
>Action GXCONT
 
>Command TRIGGER
>Parameters
+
N      'Number of events' I D=1
>Guidance
Start one or more new events.
>Action GXCONT
 
>Command RNDM
>Parameters
+
ISEED1  'First seed for the random number generator' I
ISEED2  'Second seed for the random number generator' I
>Guidance
Set the seeds for the random number generator. If no numbers are
given, the currents seeds are printed.
>Action GXCONT
 
>Command SWITCH
>Parameters
ISWI   'Switch number' I
IVAL   'New switch value' I
>Guidance
Change one element of array ISWIT(10) in /GCFLAG/
>Action GXCONT
 
 
>Command MZLOGL
>Parameters
LEVEL  'MZ log level' I D=0
>Guidance
Set the log level for the MZ package of ZEBRA: CALL MZLOGL(0,level)
 LEVEL = -3   no messages at all
         -2   error messages only
         -1   terse logging
          0   normal
         +1   log rare events
         +2   log calls to MZ routines
>Action GXCONT
 
 
>Command PRINT
>Parameters
NAME   'Name' C
NUMBER 'Number' I D=0
>Guidance
 CALL GPRINT(name,number)
>Action GXCONT
 
>Command OUTPUT_LP
>Parameters
LOUT   'New output unit' I
>Guidance
To change lout in /GCUNIT/
Note: unit numbers 5,11,12,13,14,15 are reserved and cannot be used.
>Action GXCONT
 
>Command PHITS
>Parameters
+
CHUSET  'User set' C D='*'
CHUDET  'User detector' C D='*'
NUMHI  'Hit number' I D=0
>Guidance
 CALL GPHITS(chuset,chudet)
>Action GXCONT
 
>Command PDIGI
>Parameters
+
CHUSET  'User set' C D='*'
CHUDET  'User detector' C D='*'
>Guidance
 CALL GPDIGI(chuset,chudet)
>Action GXCONT
 
>Command SMATE
>Parameters
IMAT   'Material number' I
NAMATE 'Material name' C
A      'Atomic weight' R
Z      'Atomic number' R
DENS   'Density' R
RADL   'Radiation lenght' R
ABSL   'Absorption lenght' R
UBUF   ' ' R
NWBUF  ' ' I
>Guidance
 CALL GSMATE(imat,namate,a,z,dens,radl,absl,ubuf,nwbuf)
>Action GXCONT
 
>Command SMIXT
>Parameters
IMAT   'Material number' I
NAMATE 'Material name' C
A      'Atomic weight' R
Z      'Atomic number' R
DENS   'Density' R
NLMAT  'Flag for WMAT' I
WMAT   'Relative weights or n. of atoms in molecule' R
>Guidance
 CALL GSMIXT(imat,namate,a,z,dens,nlmat,wmat)
>Action GXCONT
 
>Command PMATE
>Parameters
NUMB   'Material number' I
>Guidance
 CALL GPMATE(numb)
>Action GXCONT
 
>Command PRMAT
>Parameters
IMATE  'Material number' I
IPART  'Particle number' I
MECAN  'Mechanism' C
>Guidance
 CALL GPRMAT(imate,ipart,mecan,nekbin,elow)
>Action GXCONT
 
>Command PLMAT
>Parameters
IMATE  'Material number' I
IPART  'Particle number' I
MECAN  'Mechanism' C
+
IDM    'ID mode option' I D=0
>Guidance
CALL GPLMAT(imate,ipart,mecan,nekbin,elow,idm)
 IDM convention for histogramming mode :
 IDM.gt.0  fill, print,   keep   histogram(s)
 IDM.eq.0  fill, print,   delete histogram(s)
 IDM.lt.0  fill, noprint, keep   histogram(s)
If MECAN = 'ALL' all the mechanisms are histogrammed. If the material number
is negative, the cross sections relative to material ABS(IMATE) will
be histogrammed in barns rather than in 1/cm.
>Action GXCONT
 
>Command DRMAT
>Parameters
IMATE  'Material number' I
IPART  'Particle number' I
+
MECAN  'List of mechanism' C D='ALL'
>Guidance
CALL GDRMAT(imate,ipart,mecan,nmec)
If MECAN = 'ALL' all the mechanisms are plotted. If the material number
is negative, the cross sections relative to material ABS(IMATE) will
be plotted in barns rather than in 1/cm.
Note that it is not possible to plot anything if GSTMED has not been called
for the material number IMATE.
>Action GXCONT
 
>Command STPAR
>Parameters
ITMED  'Medium number' I
CHPAR  'Cut or mechanism' C
PARVAL 'Value' R
>Guidance
CALL GSTPAR(itmed,chpar,parval)
>Action GXCONT
 
>Command SPART
>Parameters
IPART  'Particle number' I
NAPART 'Particle name' C
ITRTYP ' ' I
AMASS  'Mass' R
CHARGE 'Charge' R
TLIFE  'Lifetime' R
UBUF   ' ' R
NWBUF  ' ' I
BRATIO 'Branching ratios' R
MODE   'Decay mode' I
>Guidance
CALL GSPART(ipart,napart,itrtyp,amass,charge,tlife,ubuf,nwbuf);
CALL GSDK(ipart,bratio,mode)
>Action GXCONT
 
>Command PPART
>Parameters
NUMB   'Particle number' I
>Guidance
CALL GPPART(numb)
>Action GXCONT
 
>Command PRKINE
>Parameters
NUMB   'Track number' I
>Guidance
CALL GPKINE(numb)
>Action GXCONT
 
>Command DEBUG
>Parameters
+
IDEB   'Debug option' C D='ON' R='ON,OFF'
>Guidance
If ideb='ON  ' then :
 idebug=1, idemin=1, idemax=1000000, itime=1
else :
 idebug=0, idemin=0, idemax=0
>Action GXCONT
 
>Name GKDZ
 
>Menu /GEANT/DZ
>Command SURV
>Parameters
NAME 'Bank name' C
+
NUMBER 'Bank number' I D=1
>Guidance
Print a survey of the structure identified by NAME, NUMBER.
>Action GXDZ
 
>Command SHOW
>Parameters
NAME 'Bank name' C
+
NUMBER 'Bank number' I D=1
CHOPT 'Options' C D='BSV'
>Guidance
Display the contents of a bank or a data structure
identified by its NAME and NUMBER.
The output format of the data part is controlled by the internal
or external I/O characteristic.
 CHOPT='B' Print the bank.
 CHOPT='S' Print the bank contents from left to right Sideways
           with up to ten elements per line.
 CHOPT='V' Print the vertical (down) structure.
 CHOPT='D' Print the bank contents from top to bottom Downwards
           with five elements per line.
 CHOPT='L' Print the linear structure.
 CHOPT='Z' Print the data part of each bank in hexadecimal format
>Action GXDZ
 
>Command SNAP
>Parameters
+
IDIV 'Division number ' I D=2 R=0:24
CHOPT 'Options' C D='M'
>Guidance
Snap of one or more divisions.
Provides a snapshot of one or more divisions in a ZEBRA store.
The kind of information provided is controlled by CHOPT.
 CHOPT='M' Print Map entry for each bank
 CHOPT='E' Extend map entry to dump all links of each bank
           (otherwise only as many links as will fit on a line)
 CHOPT='F' Full. Dump all active banks, links and data
 CHOPT='K' Kill. Dropped banks to be treated as active
           (dropped banks are not normally dumped under D or F option)
 CHOPT='L' Dump all Link areas associated with the store
 CHOPT='W' Dump the Working space, links and data
 CHOPT='Z' Dump the information in hexadecimal.
>Action GXDZ
 
>Command VERIFY
>Parameters
+
IDIV 'Division number ' I D=0 R=0:24
CHOPT 'Options' C D='CLSU'
>Guidance
Check the structure of one or more ZEBRA divisions.
The verification detail depends on the settings in CHOPT.
 CHOPT='C' Check chaining of banks only
 CHOPT='L' Check validity of the structural links (implies 'C')
 CHOPT='S' Check the store parameters
 CHOPT='U' Check the validity of the up and origin (implies 'C')
 CHOPT='F' Errors are considered fatal and generate a call to ZFATAL
>Action GXDZ
 
>Command STORE
>Parameters
+
IXSTOR 'Store number' I D=0  R=0:24
>Guidance
Display the structure of the ZEBRA store IXSTOR.
Output the parameters characterizing the store, followed by a
list of all divisions and all link areas associated with the store in
question.
>Action GXDZ
 
>Command DDIV
>Parameters
+
IDIV   'Division number'      I D=2
PATH   'Name of the doc file' C D=' '
>Guidance
Facility to display the layout of stores and divisions.
 
 CALL DZDDIV(idiv,LDUMMY,path,'IN',1,0,1,IWTYPE)
 
>Action GXDZ
 
>Command DISP
>Parameters
BANK   'Name of the bank' C
+
PATH   'Name of the doc file' C D=' '
NUMBER 'Number of the bank' I D=1
>Guidance
Interactive bank display tool.
 
 CALL DZDISP(IXSTOR,LBANK,path,'N',1,0,1,IWTYPE)
 
>Action GXDZ
 
>Command DIRZ
>Parameters
+
PATH   'Name of the RZ directory to analyse' C
>Guidance
Facility to display RZ directory trees.
 
 CALL DZDIRZ(0,LDUMMY,0,path,'N',1,0,1)
 
>Action GXDZ
 
>Name GKFZ
>Menu /GEANT/FZ
>Guidance
ZEBRA/FZ commands
 
>Command FZIN
>Parameters
LUN         'Fortran unit of the FZ file' I
KEYSU       'Name of the data structure to be retrieved' C
+
IDENT       'Version of the data structure to be retrieved' I D=0
>Guidance
Equivalent to a call to:
 
       CALL GFIN(LUN,KEYSU,1,IDENT,' ',IER)
 
>Action GXFZ
 
>Command FZOPEN
>Parameters
LUN         'Fortran unit with which to open the file' I
FILE        'Name of the file to be opened' C
LUNTYP      'Type of FZ file to be opened by GOPEN' C D='XI'
LEN         'Recordlenght of the file' I D=0
+
CHOPT       'Optional parameter to specify the action' C D=' '
>Guidance
Equivalent to a call to:
 
       CALL GOPEN(LUN,FILE,LUNTYP,LEN,IER)
 
If CHOPT = I then a call to GFIN or GFOUT will be performed in addition
according to the value of LUNTYP, with the key INIT to save or retrieve
the whole initialization data structure.
>Action GXFZ
 
>Command FZOUT
>Parameters
LUN         'Fortran unit of the FZ file' I
KEYSU       'Name of the data structure to be saved' C
+
IDENT       'Version of the data structure to be saved' I D=1
>Guidance
Equivalent to a call to:
 
       CALL GFOUT(LUN,KEYSU,1,IDENT,' ',IER)
 
>Action GXFZ
 
>Command FZCLOSE
>Parameters
LUN         'Fortran unit of the FZ to close' I
>Guidance
Equivalent to a call to:
 
       CALL GCLOSE(LUN,IER)
 
>Action GXFZ
 
>Name GKRZ
>Menu /GEANT/RZ
>Guidance
ZEBRA/RZ commands.
 
>Command PQUEST
>Parameters
+
IQ1    'Lower limit for IQ index' I D=1
IQ2    'Upper limit for IQ index' I D=20
>Guidance
Print the array IQUEST in /QUEST/.
>Action GXRZ
 
>Command FILE
>Parameters
LUN    'Logical unit number' I
FNAME 'File name' C
+
CHOPT 'Options' C D=' ' R=' ,A,N,U'
>Guidance
Open a GRZ file.
 CHOPT=' ' readonly mode
 CHOPT='U' update mode
 CHOPT='N' create new file
 CHOPT='I' Read all structures from existing file
 CHOPT='O' Write all structures on file
>Action GXRZ
 
>Command REND
>Parameters
LUNRZ  'Logical unit number' I
>Guidance
Close an RZ file opened by GRFILE on logical unit LUNRZ.
 CALL GREND(LUNRZ)
>Action GXRZ
 
 
>Command MDIR
>Parameters
CHDIR  'Directory name' C
+
CHOPT  'Options' C D=' '
>Guidance
To create a new RZ directory below the current directory.
with
 RZTAGS(1)='Object'
 RZTAGS(2)='Idvers-NR '
>Action GXRZ
 
>Command CDIR
>Parameters
+
CHPATH 'Path name' C D=' '
CHOPT  'CHOPT' C D=' '
>Guidance
Change or print the current directory.
 Ex.  CD dir1         ; make DIR1 the new CWD
      CD //file1/dir2 ; make //FILE1/DIR2 the new CWD
      CD              ; print the name of the CWD
>Action GXRZ
 
>Command IN
>Parameters
OBJECT 'Structure name' C
+
IDVERS 'Version number' I D=1
CHOPT  'Option'         C D=' '
>Guidance
Read data structure identified by OBJECT,IDVERS into memory.
  MATE read JMATE structure
  TMED read JTMED structure
  VOLU read JVOLUM structure
  ROTM read JROTM structure
  SETS read JSET  structure
  PART read JPART structure
  SCAN read LSCAN structure
  INIT read all above data structures
>Action GXRZ
 
>Command OUT
>Parameters
OBJECT 'Structure name' C
+
IDVERS 'Version number' I D=1
CHOPT  'Option'         C D=' '
>Guidance
Write data structure identified by OBJECT,IDVERS to RZ file.
  MATE write JMATE structure
  TMED write JTMED structure
  VOLU write JVOLUM structure
  ROTM write JROTM structure
  SETS write JSET  structure
  PART write JPART structure
  SCAN write LSCAN structure
  INIT write all above data structures
>Action GXRZ
 
>Command LDIR
>Parameters
+
CHPATH 'Path name' C D=' '
CHOPT  'CHOPT' C D=' '
>Guidance
List the contents of a directory (memory or disk).
To list all RZ files currently open, type 'LD //'.
>Action GXRZ
 
>Command PURGE
>Parameters
+
NKEEP  'Number of cycles to keep' I D=1
>Guidance
Purge an RZ directory.
>Action GXRZ
 
>Command SCR
>Parameters
OBJECT 'Structure name' C
+
IDVERS 'Version number' I D=1
>Guidance
Delete entry identified by OBJECT,IDVERS on RZ file.
OBJECT may be : MATE,TMED,VOLU,ROTM,SETS,PART,SCAN, *
If OBJECT= *    delete all entries with IDVERS.
>Action GXRZ
 
>Command LOCK
>Parameters
CHDIR  'Lock identifier' C D='RZFILE'
>Guidance
Lock an RZ directory.
>Action GXRZ
 
>Command FREE
>Parameters
CHDIR  'Lock identifier' C D='RZFILE'
>Guidance
Free an RZ directory.
>Action GXRZ
 
>Name GKSCAN
>Menu /GEANT/SCAN
>Guidance
To define parameters for the SCAN geometry. If the routine GUSTEP
and GUKINE are properly instrumented (see examples in GEANX),
when the TRI command is entered NTETA Geantinos will be
tracked through the real detector starting at the vertex position
defined by the command vertex. A simplified version of the geometry
is automatically generated in (ETA,PHI) or (THETA,PHI) following
the option given in the command TETA. The data structure LSCAN
generated may be saved on an RZ file for subsequent processing.
This data structure may be used for fast parametrization techniques.
 
>Command PHI
>Parameters
NPHI 'Number of PHI divisions' I D=90
+
PHIMIN 'Minimum PHI in degrees' R D=0. R=0.:360.
PHIMAX 'Maximum PHI in degrees' R D=360. R=0.:360.
>Guidance
To specify number of divisions along PHI. If no parameter is
given, the current values of the parameters are displayed.
>Action GXSCAN
 
>Command TETA
>Parameters
NTETA 'Number of TETA divisions' I D=90
+
TETMIN 'Minimum value of TETA' R
TETMAX 'Maximum value of TETA' R
DIVTYP 'Type of TETA division' I R=1:3
>Guidance
To specify number of divisions along TETA.
If DIVTYP=1 divisions in pseudo-rapidity ETA.
If DIVTYP=2 divisions in degrees following the THETA angle.
If DIVTYP=3 divisions in cos(TETA).
If no parameter is given, the current values of the parameters
are displayed.
>Action GXSCAN
 
>Command SLIST
>Parameters
LIST 'List of master volumes' C
>Guidance
Only boundary crossings of volumes given in LIST will be seen
in the SCAN geometry. If no parameters are given, the current
SCAN volumes will be listed. If a full stop (.) is given, the list
of scan volumes will be erased.
>Action GXSCAN
 
>Command VERTEX
>Parameters
VX 'Scan X-origin' R D=0.
VY 'Scan Y-origin' R D=0.
VZ 'Scan Z-origin' R D=0.
>Guidance
All Geantinos tracked will start from position VX,VY,VZ.
>Action GXSCAN
 
>Command SFACTORS
>Parameters
FACTX0 'Scale factor for SX0' R D=100.
FACTL  'Scale factor for SL' R D=1000.
FACTR  'Scale factor for R'  R D=100.
>Guidance
Set scale factors for SX0,SL and R. The given scale factors must be
such that:
  SX0*FACTX0 < 2**15-1 (32767)
  SL*FACTL   < 2**10-1 (1023)
  SR*FACTR   < 2**17-1 (131071)
>Action GXSCAN
 
>Command STURN
>Parameters
CHOPT 'SCAN mode setting' C R='ON,OFF,INIT'
>Guidance
Switch on/off SCAN mode. If SCAN mode is on, SCAN geantinos
are generated and tracked to fill (or complete) the current
scan data structure. If SCAN mode is off, normal kinematics
generation and tracking will take place. If INIT is given,
the current SCAN data structure (if any) will be dropped
and SCAN mode will be turned on.
>Action GXSCAN
 
>Command PCUTS
>Parameters
+
IPARAM   'Parametrization Flag'                        I R=0:1
PCUTGA   'Parametrization Cut for gammas'              R
PCUTEL   'Parametrization Cut for electrons'           R
PCUTHA   'Parametrization Cut for charged hadrons'     R
PCUTNE   'Parametrization Cut for neutral hadrons'     R
PCUTMU   'Parametrization Cut for muons'               R
>Guidance
Control parametrization at tracking time.
 
     IPARAM=0       No parametrization is performed
     IPARAM=1       Parametrization is performed
 
If parametrization is active and a particle falls below its
parametrization cut, then the particle will be replaced by
a parametrized shower which will be tracked in the SCAN
geometry.
>Action GXSCAN
 
>Command LSCAN
>Parameters
ID 'Lego plot identifier' I D=2000
+
VOLUME 'Volume name' C D='XXXX'
CHOPT 'List of options' C D='OPX' R=' ,O,P,I,X,L'
>Guidance
Generates and plot a table of physics quantities such as
the total number of radiation lengths or interaction lengths
in function of the SCAN parameters TETA,PHI.
  CHOPT='O' table is generated at Exit  of VOLUME.
  CHOPT='I' table is generated at Entry of VOLUME.
  CHOPT='X' radiation lengths
  CHOPT='L' Interaction lengths
  CHOPT='P' Plot the table
If VOLUME='XXXX' Mother volume is used.
>Action GXSCAN
 
>Command HSCAN
>Parameters
IDPHI 'Histogram/phi identifier' I D=1000
+
VOLUME 'Volume name' C D='XXXX'
CHOPT 'List of options' C D='OPX' R=' ,O,P,I,X,L'
>Guidance
Generates and plots an histogram of physics quantities such as
the total number of radiation lengths or interaction lengths
as a function of the SCAN parameter TETA for a given value of PHI.
  CHOPT='O' histogram is generated at Exit  of VOLUME.
  CHOPT='I' histogram is generated at Entry of VOLUME.
  CHOPT='X' radiation lengths
  CHOPT='L' Interaction lengths
  CHOPT='P' Plot the histogram
If VOLUME='XXXX' Mother volume is used.
The histogram identifier IDPHI is used to also identify which
PHI division to plot: IPHI=MOD(IDPHI,1000).
If IPHI=0, then all PHI divisions are generated (not plotted)
with histogram identifiers IDPHI+PHI division number.
>Action GXSCAN
 
>Name GKPHYS
>Menu /GEANT/PHYSICS
>Guidance
Commands to set physics parameters.
 
>Command ANNI
>Parameters
+
IANNI 'Flag IANNI' I D=1 R=0,1,2
>Guidance
To control positron annihilation.
 IANNI=0 no annihilation
      =1 annihilation. Decays processed.
      =2 annihilation. No decay products stored.
>Action GXPHYS
 
>Command AUTO
>Parameters
+
IAUTO 'Flag IAUTO' I D=1 R=0,1
>Guidance
To control automatic calculation of tracking medium parameters:
 IAUTO=0 no automatic calculation;
      =1 automati calculation.
>Action GXPHYS
 
>Command BREM
>Parameters
+
IBREM 'Flag IBREM' I D=1 R=0,1,2
>Guidance
To control bremstrahlung.
 IBREM=0 no bremstrahlung
      =1 bremstrahlung. Photon processed.
      =2 bremstrahlung. No photon stored.
>Action GXPHYS
 
>Command CKOV
>Parameters
+
ICKOV 'Flag ICKOV' I D=0 R=0,1,2
>Guidance
To control Cerenkov production
 ICOMP=0 no Cerenkov;
      =1 Cerenkov;
      =2 Cerenkov with primary stopped at each step.
>Action GXPHYS
 
>Command COMP
>Parameters
+
ICOMP 'Flag ICOMP' I D=1 R=0,1,2
>Guidance
To control Compton scattering
 ICOMP=0 no Compton
      =1 Compton. Electron processed.
      =2 Compton. No electron stored.
>Action GXPHYS
 
>Command DCAY
>Parameters
+
IDCAY 'Flag IDCAY' I D=1 R=0,1,2
>Guidance
To control Decay mechanism.
 IDCAY=0 no decays.
      =1 Decays. secondaries processed.
      =2 Decays. No secondaries stored.
>Action GXPHYS
 
>Command DRAY
>Parameters
+
IDRAY 'Flag IDRAY' I D=1 R=0,1,2
>Guidance
To control delta rays mechanism.
 IDRAY=0 no delta rays.
      =1 Delta rays. secondaries processed.
      =2 Delta rays. No secondaries stored.
>Action GXPHYS
 
>Command ERAN
>Parameters
+
EKMIN  'Minimum energy of the tables' R D=1E-5
EKMAX  'Maximum energy of the tables' R D=1E+4
NEKBIN 'Number of bins in the tables' I D=90 R=1:200
>Guidance
To define the range and binning of internal tables.
>Action GXPHYS
 
>Command HADR
>Parameters
+
IHADR 'Flag IHADR' I D=1
>Guidance
To control hadronic interactions.
 IHADR=0 no hadronic interactions.
      =1 Hadronic interactions. secondaries processed.
      =2 Hadronic interactions. No secondaries stored.
>Action GXPHYS
 
>Command LABS
>Parameters
+
LABS   'Flag LABS' I D=0
>Guidance
To control absorbtion of Cerenkov photons:
    LABS=0 no absorbtion of photons;
    LABS=1 absorbtion of photons;
>Action GXPHYS
 
>Command LOSS
>Parameters
+
ILOSS 'Flag ILOSS' I D=2 R=0,1,2,3,4
>Guidance
To control energy loss.
 ILOSS=0 no energy loss;
      =1 restricted energy loss fluctuations;
      =2 complete energy loss fluctuations;
      =3 same as 1;
      =4 no energy loss fluctuations.
If the value ILOSS is changed, then cross-sections and energy loss
tables must be recomputed via the command 'PHYSI'.
>Action GXPHYS
 
>Command MULS
>Parameters
+
IMULS 'Flag IMULS' I D=1 R=0,1,2,3
>Guidance
To control multiple scattering.
 IMULS=0 no multiple scattering.
      =1 Moliere or Coulomb scattering.
      =2 Moliere or Coulomb scattering.
      =3 Gaussian scattering.
>Action GXPHYS
 
>Command MUNU
>Parameters
+
IMUNU 'Flag IMUNU' I D=1 R=0,1,2
>Guidance
To control muon nuclear interactions.
 IMUNU=0 no muon-nuclear interactions.
      =1 Nuclear interactions. Secondaries processed.
      =2 Nuclear interactions. Secondaries not processed.
>Action GXPHYS
 
>Command PAIR
>Parameters
+
IPAIR 'Flag IPAIR' I D=1 R=0,1,2
>Guidance
To control pair production mechanism.
 IPAIR=0 no pair production.
      =1 Pair production. secondaries processed.
      =2 Pair production. No secondaries stored.
>Action GXPHYS
 
>Command PFIS
>Parameters
+
IPFIS 'Flag IPFIS' I D=1 R=0,1,2
>Guidance
To control photo fission mechanism.
 IPFIS=0 no photo fission.
      =1 Photo fission. secondaries processed.
      =2 Photo fission. No secondaries stored.
>Action GXPHYS
 
>Command PHOT
>Parameters
+
IPHOT 'Flag IPHOT' I D=1 R=0,1,2
>Guidance
To control Photo effect.
 IPHOT=0 no photo electric effect.
      =1 Photo effect. Electron processed.
      =2 Photo effect. No electron stored.
>Action GXPHYS
 
>Command RAYL
>Parameters
+
IRAYL 'Flag IRAYL' I D=1 R=0,1
>Guidance
To control Rayleigh scattering.
 IRAYL=0 no Rayleigh scattering.
      =1 Rayleigh.
>Action GXPHYS
 
>Command STRA
>Parameters
+
ISTRA 'Flag ISTRA' I D=0 R=0,1,2
>Guidance
To control energy loss fluctuation model:
 ISTRA=0 Urban model;
      =1 PAI model;
      =2 PAI+ASHO model (not active at the moment).
>Action GXPHYS
 
>Command SYNC
>Parameters
+
ISYNC 'Flag ISYNC' I D=1 R=0,1
>Guidance
To control synchrotron radiation:
 ISYNC=0 no synchrotron radiation;
      =1 synchrotron radiation.
>Action GXPHYS
 
>Command CUTS
>Parameters
+
CUTGAM   'Cut for gammas'              R D=0.001
CUTELE   'Cut for electrons'           R D=0.001
CUTHAD   'Cut for charged hadrons'     R D=0.01
CUTNEU   'Cut for neutral hadrons'     R D=0.01
CUTMUO   'Cut for muons'               R D=0.01
BCUTE    'Cut for electron brems.'     R D=-1.
BCUTM    'Cut for muon brems.'         R D=-1.
DCUTE    'Cut for electron delta-rays' R D=-1.
DCUTM    'Cut for muon delta-rays'     R D=-1.
PPCUTM   'Cut for e+e- pairs by muons' R D=0.01
TOFMAX   'Time of flight cut'          R D=1.E+10
GCUTS    '5 user words'                R D=0.
>Guidance
To change physics cuts. If no parameter is given, the list
of the current cuts is printed.
 If the default values (-1.) for       BCUTE ,BCUTM ,DCUTE ,DCUTM
 are not modified, they will be set to CUTGAM,CUTGAM,CUTELE,CUTELE
 respectively.
If one of the parameters from CUTGAM to PPCUTM included
is modified, cross-sections and energy loss tables must be
recomputed via the command 'PHYSI'.
>Action GXPHYS
 
>Command DRPRT
>Parameters
IPART    'GEANT particle number'         I
IMATE    'GEANT material number'         I
STEP     'step length in centimeters'    R
+
NPOINT 'number of logarithmically spaced energy points' I D=10 R=2:100
>Guidance
This routine prints the relevant parameters linked with the energy loss
fluctuation.
>Action GXPHYS
 
>Command PHYSI
>Guidance
Call the GEANT initialisation routine GPHYSI to recompute
the tables of cross-sections and energy loss. This command
must be invoked after CUTS, LOSS or ERAN commands.
>Action GXPHYS
 
>Name GKFORT
>Menu FORTRAN
 
>Command FORTRAN
>Parameters
FNAME 'File name' C
>Guidance
The routines in the file FNAME will be compiled by COMIS.
If routines with names: UGEOM,GUKINE,GUOUT,UGLAST are found,
then they will be automatically called by GXINT instead of
the routines with the same names compiled with the standard
Fortran compiler and linked with the application.
The user callable routines from the GEANT library as well as
routines from PACKLIB (HBOOK,HPLOT,HIGZ,ZEBRA) may be called
from these user routines. All GEANT common blocks may be
referenced.
In case where the routine UGEOM is called several times,
it is important to DROP all the initialisation data structures
JVOLUM,JMATE,JTMED,etc already in memory by using the routine GIDROP.
 Example of an interactive session where the routine UGEOM is modified:
.
   GEANT > Edit ugeom.for
   GEANT > Fortran ugeom.for
   GEANT > Call GIDROP
   GEANT > Call UGEOM
   GEANT > Dtree
   GEANT > Edit ugeom.for
   GEANT > Fortran ugeom.for
   GEANT > Call GIDROP
   GEANT > Call UGEOM
   GEANT > Dtree
 
If FNAME='-', calls to user routines is reset and standard
routines called instead.
>Action GXFORT
 
 
