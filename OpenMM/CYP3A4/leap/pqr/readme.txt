The following residue assignments have been made by PDB2PQR:

# grep ^ATOM ../leap/1TQN.pdb | awk '{print $4, $6}'  | uniq | grep -v WAT  > seqres.pdb
# grep ^ATOM         1TQN.pqr | awk '{print $4, $5}'  | uniq | grep -v WAT  > seqres.pqr
# diff seqres.p{db,qr} | grep ">" | sed s/"> "//

HIE 30
HID 54
HIE 65
HID 267
HID 287
HIE 324
HID 402


However, using PROPKA (2), the following is predicted:

HIP 30
HIP 54
HIE 65
HID 267
HID 287
HIE 324
HID 402
CYM 468

Using PROPKA (2) on a PDB generated from the end of stage 1 (extra residues (282->285) filtered out here)

HIE 28
HIE 30
HIE 54
HIE 65
HIE 267
HIE 287
HIE 324
HIE 402
CYM 468


Using PROPKA (2) on 1TQN_with_missing_residues.pdb  (extra residues (282->285) filtered out in results here)

HIP 30
HIP 54
HIE 65
HID 267
HID 287
HIE 324
HID 402
CYM 468



