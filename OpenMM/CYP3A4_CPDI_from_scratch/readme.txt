This is an example whereby CYP3A4 is built from scratch, 
using Modeller and parameters for the CPDI state and
the modified cysteine residues that have been converted
using JAMBER.

1) CYP must be added to /usr/local/lib/python2.6/dist-packages/simtk/openmm/app/data/residues.xml

 <Residue name="CYP">
  <Bond from="-C" to="N"/>
  <Bond from="C" to="CA"/>
  <Bond from="C" to="O"/>
  <Bond from="C" to="OXT"/>
  <Bond from="CA" to="CB"/>
  <Bond from="CA" to="HA"/>
  <Bond from="CA" to="N"/>
  <Bond from="CB" to="HB2"/>
  <Bond from="CB" to="HB3"/>
  <Bond from="CB" to="SG"/>
  <Bond from="H" to="N"/>
  <Bond from="H2" to="N"/>
  <Bond from="H3" to="N"/>
 </Residue>

<Residue name="HEM">
	<Bond from="NC" to="C1C" />
	<Bond from="NC" to="C4C" />
	<Bond from="NC" to="FE" />
	<Bond from="C1C" to="C2C" />
	<Bond from="C1C" to="CHC" />
	<Bond from="C4C" to="C3C" />
	<Bond from="C4C" to="CHD" />
	<Bond from="C2C" to="C3C" />
	<Bond from="C2C" to="CMC" />
	<Bond from="C3C" to="CAC" />
	<Bond from="CHD" to="HHD" />
	<Bond from="CHD" to="C1D" />
	<Bond from="C1D" to="ND" />
	<Bond from="C1D" to="C2D" />
	<Bond from="ND" to="C4D" />
	<Bond from="ND" to="FE" />
	<Bond from="C4D" to="C3D" />
	<Bond from="C4D" to="CHA" />
	<Bond from="C3D" to="C2D" />
	<Bond from="C3D" to="CAD" />
	<Bond from="C2D" to="CMD" />
	<Bond from="CHA" to="HHA" />
	<Bond from="CHA" to="C1A" />
	<Bond from="C1A" to="C2A" />
	<Bond from="C1A" to="NA" />
	<Bond from="C2A" to="C3A" />
	<Bond from="C2A" to="CAA" />
	<Bond from="C4A" to="NA" />
	<Bond from="C4A" to="C3A" />
	<Bond from="C4A" to="CHB" />
	<Bond from="NA" to="FE" />
	<Bond from="C3A" to="CMA" />
	<Bond from="CHB" to="HHB" />
	<Bond from="CHB" to="C1B" />
	<Bond from="C1B" to="C2B" />
	<Bond from="C1B" to="NB" />
	<Bond from="C2B" to="C3B" />
	<Bond from="C2B" to="CMB" />
	<Bond from="NB" to="C4B" />
	<Bond from="NB" to="FE" />
	<Bond from="C4B" to="CHC" />
	<Bond from="C4B" to="C3B" />
	<Bond from="CHC" to="HHC" />
	<Bond from="FE" to="O1" />
	<Bond from="C3B" to="CAB" />
	<Bond from="CAB" to="HAB" />
	<Bond from="CAB" to="CBB" />
	<Bond from="CBB" to="HBB1" />
	<Bond from="CBB" to="HBB2" />
	<Bond from="CAC" to="HAC" />
	<Bond from="CAC" to="CBC" />
	<Bond from="CBC" to="HBC1" />
	<Bond from="CBC" to="HBC2" />
	<Bond from="CMB" to="HMB1" />
	<Bond from="CMB" to="HMB2" />
	<Bond from="CMB" to="HMB3" />
	<Bond from="CMC" to="HMC1" />
	<Bond from="CMC" to="HMC2" />
	<Bond from="CMC" to="HMC3" />
	<Bond from="CMD" to="HMD1" />
	<Bond from="CMD" to="HMD2" />
	<Bond from="CMD" to="HMD3" />
	<Bond from="CMA" to="HMA1" />
	<Bond from="CMA" to="HMA2" />
	<Bond from="CMA" to="HMA3" />
	<Bond from="CAA" to="HAA1" />
	<Bond from="CAA" to="HAA2" />
	<Bond from="CAA" to="CBA" />
	<Bond from="CAD" to="HAD1" />
	<Bond from="CAD" to="HAD2" />
	<Bond from="CAD" to="CBD" />
	<Bond from="CBA" to="HBA1" />
	<Bond from="CBA" to="HBA2" />
	<Bond from="CBA" to="CGA" />
	<Bond from="CBD" to="HBD1" />
	<Bond from="CBD" to="HBD2" />
	<Bond from="CBD" to="CGD" />
	<Bond from="CGD" to="O1D" />
	<Bond from="CGD" to="O2D" />
	<Bond from="CGA" to="O1A" />
	<Bond from="CGA" to="O2A" />
</Residue>


2) CYP and HEM hydrogen atoms must be added to /usr/local/lib/python2.6/dist-packages/simtk/openmm/app/data/hydrogens.xml

 <Residue name="CYP">
  <Variant name="CYP"/>
  <H name="H" parent="N"/>
  <H name="H2" parent="N" terminal="N"/>
  <H name="H3" parent="N" terminal="N"/>
  <H name="HA" parent="CA"/>
  <H name="HB2" parent="CB"/>
  <H name="HB3" parent="CB"/>
 </Residue>

<Residue name="HEM">
	<H name="HHD" parent="CHD" />
	<H name="HHA" parent="CHA" />
	<H name="HHB" parent="CHB" />
	<H name="HHC" parent="CHC" />
	<H name="HAB" parent="CAB" />
	<H name="HBB1" parent="CBB" />
	<H name="HBB2" parent="CBB" />
	<H name="HAC" parent="CAC" />
	<H name="HBC1" parent="CBC" />
	<H name="HBC2" parent="CBC" />
	<H name="HMB1" parent="CMB" />
	<H name="HMB2" parent="CMB" />
	<H name="HMB3" parent="CMB" />
	<H name="HMC1" parent="CMC" />
	<H name="HMC2" parent="CMC" />
	<H name="HMC3" parent="CMC" />
	<H name="HMD1" parent="CMD" />
	<H name="HMD2" parent="CMD" />
	<H name="HMD3" parent="CMD" />
	<H name="HMA1" parent="CMA" />
	<H name="HMA2" parent="CMA" />
	<H name="HMA3" parent="CMA" />
	<H name="HAA1" parent="CAA" />
	<H name="HAA2" parent="CAA" />
	<H name="HAD1" parent="CAD" />
	<H name="HAD2" parent="CAD" />
	<H name="HBA1" parent="CBA" />
	<H name="HBA2" parent="CBA" />
	<H name="HBD1" parent="CBD" />
	<H name="HBD2" parent="CBD" />
</Residue>





2) OXT must be added to THR
	ATOM   3767  OXT THR A 499     -27.630 -61.169 -23.602  1.00 88.35           O

3) Bug must be fixed in internal/pdbstructure.py and pdbfile.py so that the atom's
element gets passed up to PDBFile's topology.

(REPORTED: https://simtk.org/tracker/index.php?func=detail&aid=1792&group_id=161&atid=435 )

For example, within the PDB 1TQN:

....
HETATM 3772  CHD HEM A 508     -15.552 -23.643  -7.998  1.00 36.27           C
HETATM 3773  NA  HEM A 508     -15.862 -21.317 -12.447  1.00 33.91           N
HETATM 3774  C1A HEM A 508     -15.597 -20.050 -12.014  1.00 31.80           C
HETATM 3775  C2A HEM A 508     -15.544 -19.183 -13.121  1.00 31.63           C
HETATM 3776  C3A HEM A 508     -15.780 -19.956 -14.223  1.00 32.28           C
....

atom 3373 will be classed as a Sodium (Na) and not Nitrogen.


Within the pdbfile.py:__init__() function, a heuristic is used to guess an atom type's
element based upon its PDB name, whilst trying to build the PDBFile topology:

                    if element is None:
                        # Try to guess the element.

                        upper = atomName.upper()
                        if upper.startswith('CL'):
                            element = elem.chlorine
                        elif upper.startswith('NA'):



Surely, an attempt should be made to use information in internal/pdbstructure.py 
before falling back to this heuristic, i.e. another method should be added to
internal/pdbstructure.py:

    def get_name_with_spaces(self):
        return self._name_with_spaces
    name_with_spaces = property(get_name_with_spaces, set_name_with_spaces, doc='four-character residue name including spaces')

    def get_element(self):
        return self.element

    def get_name(self):
        return self._name
    name = property(get_name, doc='residue name')


Then in pdbfile.py:__init__(), the following should be done:

                for atom in residue.atoms:
                    atomName = atom.get_name()
                    if atomName in atomReplacements:
                        atomName = atomReplacements[atomName]
                    atomName = atomName.strip()

                    # Wrap this in some form of try

                    element = atom.get_element()

                    #Then, perhaps fallback to heuristic
                    # ...etc...

                    newAtom = top.addAtom(atomName, element, r)
                    atomByNumber[atom.serial_number] = newAtom
                    pos = atom.get_position().value_in_unit(nanometers)
                    coords.append(Vec3(pos[0], pos[1], pos[2]))




This has been fixed in OpenMM 5.0


4.1) What about:
	/usr/local/lib/python2.6/dist-packages/simtk/openmm/app/data/residues.xml

5) We are still missing two improper dihedrals from the simulation:

	  i- j -k -l

	CBA-O1A-CGA-O2A
	 61- 69- 68- 71
	 c3-  o-  c-  o

	CBD-O1D-CGD-O2D
	 64- 70- 67- 72
	 c3-  o-  c-  o

By AMBER convention, the central file is in the k position of the improper.
These are the propionate impropers:

        O2A
          \
          CGA-CBA
          /
        O1A

For some reason, the gaff term, is not being matched.
	X -o -c -o          1.1          180.          2.           JCC,7,(1986),230

This is present in the XML:
	<Improper class1="" class2="o" class3="c" class4="o" periodicity1="2" k1="4.6024" phase1="3.141592653589793"/>

but is not being matched. Forcefield.createSystem() generates a list of impropers from a topology. However, it seems to storer the dihedrals in a different order:

Forcefield.data.impropers: 68 71 61 69
Forcefield.data.impropers: 67 64 72 70

  It seems that with all impropers in converting from parm to OpenMMXML:

	1) Atomtypes in positions 1 and 3 are swapped 

	2) Then Atomtypes in positions 2 and 3 are swapped

i.e.
	C -CT-N -O 1.1 180. 2. Junmei et al.1999
	 
	C -CT-N -O 
		=> N -CT-C -O 
			=> N -C -CT-O


         Finally to :
	 <Improper class1="N" class2="C" class3="CT" class4="O" periodicity1="2" phase1="3.14159265359" k1="4.6024"/>


6) Tests with the CPDI only residue (OpenMM/CPDI) have verified that the modeller.addHydrogens is working ok. However,
 with the complete model (OpenMM/CYP3A4_CPDI_from_scratch) is not connecting Fe to SH and also, some hydrogen atoms
are being added incorrectly to the HEM. This is manifest after a minimisation:
	HBD2
	HBB2
	HMC2
	HAA2

There are also some other issues seen with other hydrogens as well.
In fact, the problem is seen (in modeller.pdb) after modeller has operated on the structure:
	PDBFile.writeFile(modeller.getTopology(), modeller.getPositions(),open('modeller.pdb', 'w'))
Problems with the HB2 CYP residue are also seen.

Could it be an issue with:
	HEM's ExternalBond
	CYP's ExternalBonds
	Connectivity information in simtk/openmm/app/data/residues.xml

		

Fixed; two issues. The type name, i.e. foo_0 here

	<Type name="foo_0" class="nc" element="N" mass="14.01"/>

has to unique across ALL xml loaded in.

Secondly, here: 

	forceField = ForceField('CPDI_CYP.xml', 'amber99sb.xml', 'tip3p.xml')

# if amber99sb is
# first, CYP will be matched against CYX in amber99sb.xml
# and not CYP in CPDI_CYP.xml


Notes
=====

OpenMM will not add any zero k valued angles or torsions, whereas leap will.

Why is there a duplicate term in CPDI_CYP.xml:
	$ grep SH CPDI_CYP.xml  | uniq -d
 <Bond class1="fe" class2="SH" length="0.2565" k="32635.200000000004"/>

