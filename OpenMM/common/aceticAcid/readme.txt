The files in this directory are the result of an automated
process to parameterise acetic acid using the AMBER tool
ANTECHAMBER. The aim is to generate an example AMBER topology/
parameter assignment.

The process starts with a hand drawn molecule in the ./raw
directory and results in an AMBER prmtop and inpcrd. The
assignment of charges and atomtypes is done using Antechamber
and the system is built using tleap. Antechambe assigns 
GAFF parameters which are compatible with parm94.
Finally, the energy of the system is calculated using SANDER.

The entire process is controlled by a Make file.




This is the final topology of the molecule:


GAFF atom types		Atom ids
--------------          --------

hc         o                H2        O
  \       /                  \       /
   \     /                    \     /
hc--c3--c                   H--C--C1
   /     \                    /     \
  /       \                  /       \
hc         oh--ho          H1         O1-H3



The antechamber generated mol2 file contains GAFF atomtype
assignments and charges:


@<TRIPOS>MOLECULE
ACE
    8     7     1     0     0
SMALL
bcc


@<TRIPOS>ATOM
      1 C          -3.9150   -0.4930    0.0000 c3        1 ACE     -0.200100
      2 C1         -2.7230    0.4150    0.0490 c         1 ACE      0.623100
      3 H          -4.8240    0.1060   -0.0980 hc        1 ACE      0.086700
      4 H1         -3.9740   -1.0750    0.9240 hc        1 ACE      0.068700
      5 H2         -3.8440   -1.1540   -0.8680 hc        1 ACE      0.068700
      6 O          -2.7880    1.6330    0.0830 o         1 ACE     -0.485000
      7 O1         -1.5520   -0.2590    0.0650 oh        1 ACE     -0.584100
      8 H3         -1.6580   -1.2280    0.0390 ho        1 ACE      0.422000
@<TRIPOS>BOND
     1    1    2 1
     2    1    3 1
     3    1    4 1
     4    1    5 1
     5    6    2 2
     6    2    7 1
     7    7    8 1
@<TRIPOS>SUBSTRUCTURE
     1 ACE         1 TEMP              0 ****  ****    0 ROOT


Dissecting the prmtop and inpcrd files, one can see that the following
parameters were used:


Atom 
-----

    MASS          ATPOL
c  12.01         0.616               Sp2 C carbonyl group
c3 12.01         0.878               Sp3 C
hc 1.008         0.135               H bonded to aliphatic carbon without electrwd. group
o  16.00         0.434               Oxygen with one connected atom
oh 16.00         0.465               Oxygen in hydroxyl group
ho 1.008         0.135               Hydroxyl group


         AMASS      Atomic mass of the center having the symbol "KNDSYM".

         ATPOL      The atomic polarizability for each atom (in A**3)
                    This is the type of polarizability used in sander
                    and gibbs. No parameters are supplied for this since
                    the feature is still in development (Amber 4.1).


Bond
-----
         RK       REQ
c3-hc  337.3    1.0920       SOURCE3    2815    0.0059
ho-oh  369.6    0.9740       SOURCE3     367    0.0105
c -c3  328.3    1.5080       SOURCE1    2949    0.0060
c -o   648.0    1.2140       SOURCE1    3682    0.0165
c -oh  466.4    1.3060       SOURCE1     271    0.0041

Where

         RK         The harmonic force constant for the bond "IBT"-"JBT".
                    The unit is kcal/mol/(A**2).

         REQ        The equilibrium bond length for the above bond in Angstroms

Angle
-----

             TK         TEQ
c -c3-hc   47.200     109.680   SOURCE3          614    0.6426
c -oh-ho   51.190     107.370   SOURCE3           34    1.6830
hc-c3-hc   39.430     108.350   SOURCE3         2380    0.9006
c3-c -o    68.030     123.110   SOURCE3          267    3.0977
c3-c -oh   69.840     112.200   SOURCE3           14    1.8324
o -c -oh   77.380     122.880   SOURCE3           33    2.1896

Where:

         TK         The harmonic force constants for the angle "ITT"-"JTT"-
                    "KTT" in units of kcal/mol/(rad**2) (radians are the
                    traditional unit for angle parameters in force fields).

         TEQ        The equilibrium bond angle for the above angle in degrees.


Dihedrals
---------
[c3-c -oh-ho]

             IDIVF  PK           PHASE             PN 
X -c -oh-X    2    4.600       180.000           2.000      Junmei et al, 1999

[hc-c3-c -o ]
hc-c3-c -o    1    0.80          0.0            -1.         Junmei et al, 1999
hc-c3-c -o    1    0.08        180.0             3.         Junmei et al, 1999


[o-c-oh-ho]
ho-oh-c -o    1    2.30        180.0            -2.         Junmei et al, 1999
ho-oh-c -o    1    1.90          0.0             1.         Junmei et al, 1999

[c3-o-c-oh]
c3-o -c -oh         1.1          180.          2.

where:

         IDIVF      The factor by which the torsional barrier is divided.
                    Consult Weiner, et al., JACS 106:765 (1984) p. 769 for
                    details. Basically, the actual torsional potential is

                           (PK/IDIVF) * (1 + cos(PN*phi - PHASE))

         PK         The barrier height divided by a factor of 2.

         PHASE      The phase shift angle in the torsional function.

                    The unit is degrees.

         PN         The periodicity of the torsional barrier.
                    NOTE: If PN .lt. 0.0 then the torsional potential
                          is assumed to have more than one term, and the
                          values of the rest of the terms are read from the
                          next cards until a positive PN is encountered.  The
                          negative value of pn is used only for identifying
                          the existence of the next term and only the
                          absolute value of PN is kept.


Note for more information about dihedrals, please see page 254
of http://ambermd.org/doc8/amber8.pdf


VDW
===
		R      EDEP
  hc          1.4870  0.0157             OPLS
  c3          1.9080  0.1094             OPLS
  c           1.9080  0.0860             OPLS
  o           1.6612  0.2100             OPLS
  oh          1.7210  0.2104             OPLS
  ho          0.0000  0.0000             OPLS Jorgensen, JACS,110,(1988),1657

where: 
	R 	= The van der Waals radius of the atoms having the symbol (angstroms)
	EDEP	= The 6-12 potential well depth. (kcal/mol)


Refs
====
http://ambermd.org/formats.html#parm.dat
http://ambermd.org/formats.html#topology
