import os
os.system("mkdir Prokka_annotations")# Making a folder to store all the annotations from prokka
os.chdir("Assemblies")# Changing directory to Assemblies to make Prokka work
os.system('ls *.fna | parallel --verbose "prokka {} --prefix {.}_out"')# Run prokka and start annotate all the .fna input 

os.chdir("..")# Back to the previous menu to organize output files
command1 = "mv ./Assemblies/*v1_genomic_out Prokka_annotations"# move out outputfiles
os.system(command1)

command2 = "find ./Prokka_annotations -name \*.faa -exec cp {} Proteinseq \;"#find out .faa protein files and store them into Proteinseq file
os.system(command2)
