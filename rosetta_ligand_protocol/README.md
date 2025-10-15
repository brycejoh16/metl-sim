# Rosetta Ligand Protocol Exploration

## Docker Image
1. Dowload docker image (9 GB total)
```angular2html
docker pull arnvsharma/metl-sim:latest
```
2. Give docker containers permission to access metl-sim folder on docker app.
3. Start a mounted docker container 
```angular2html
docker run -it -v /Users/brycejohnson/Desktop/proteins/metl-sim/:/rosetta arnvsharma/metl-sim:latest /bin/bash  
```

4. Cd into `/rosetta` in the docker container. 


## Protocols 

I would like to explore a few of the parameters. The parameters I would like to explore include: 
- nstructs 
- how the best term is chosen from nstructs 
- fixing the bb and chi coordinates of metal coordinating atoms, and then also removing the 
auto metal movers setup. 
- Will do a second round where I look at the effects of the autometal movers set up. 

# Runnign the Original protocol
I will run a brief protocol making the following mutations : 
```
SadA_rosetta_2024_3_6_p.pdb R83A
SadA_rosetta_2024_3_6_p.pdb E95A
```
Then I will run the ligand protocol. 
```angular2html
python code/sadA_rosetta_ligand.py --save_wd --log_dir_base rosetta_ligand_protocol/output --num_structs 1 --variants_fn variant_lists/SadA_rosetta_2024_3_6_p_example.txt
```

Great! It is making both of the mutations. Now I can start making changes! 

## Experiments to run:

- protocols.qsar.scoring_grid.GridManager: [ WARNING ] option -qsar:grid_dir is not set.  Use this flag to specify a directory to store scoring grids.  This will save you a huge amount of time

- Lowest interface delta X 
- Clustering 
- Increasing the ligand radius for the HighResolutionDocker 
- FastRelax with constraints

- Transform or no transform 
- Conformers or no conformers for the NSL 
- Number of Nstructs  




Where I left of Oct 13 , 2025 

I have to still get the container working CHTC again. I can't seem to figure it out. But the 
aptainer wouldn't let me run it in an interactive environment. Will need to open up an issue. 



## Increasing 

# Increasing nstructs

I will start with a single mutation. 