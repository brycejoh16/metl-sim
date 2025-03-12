#!/bin/bash

pos=$1
mut=$2

dock_xml="mutate_dock_relax.xml"

temp_xml="temp.xml"

content=$(<"mutate_dock_relax.xml")

echo "$content" > "$temp_xml"

sed -i "s/<MutateResidue name=\"mutant1\" target=\"A\" new_res=\"MUT\"\/>/<MutateResidue name=\"mutant1\" target=\"${pos}\" new_res=\"${mut}\"\/>/" "$temp_xml"


#transfer files from /software

cp /software/alira3/rosetta3.13/main/source/bin/rosetta_scripts.static.linuxgccrelease .
 
cp -R /software/alira3/rosetta3.13/main/database .

#set up file


#Run Rosetta

tar -xzf database.tar.gz
chmod a+x rosetta_scripts.static.linuxgccrelease


./rosetta_scripts.static.linuxgccrelease -database ./database -parser:protocol temp.xml -s SadA_NSLeu_Corrected_3701_best_structure_0044_correct_seq.pdb -restore_pre_talaris_2013_behavior true -in:auto_setup_metals -extra_res_fa AKG.params -extra_res_fa NEU.params -ex1 -ex2 -no_optH false -flip_HNQ true -ignore_ligand_chi true -nstruct 1 -out:pdb -overwrite

#./relax.static.linuxgccrelease $1







