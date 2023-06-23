from amberti.utils import run_command
from amberti.amber import antechamber

"""
antechamber -i M7G.sdf -fi sdf -o M7G.gaff.mol2 -fo mol2 -rn M7G -at gaff -an yes -dr no -pf yes -c bcc -nc -2
parmchk2 -i M7G.gaff.mol2 -f mol2 -o M7G.gaff.frcmod -s gaff -a yes

tleap -f - <<_EOF
source leaprc.gaff
loadamberparams M7G.gaff.frcmod
M7G = loadmol2 M7G.gaff.mol2 
saveoff M7G M7G.lib
savepd b M7G M7G.pdb
quit
_EOF

"""



