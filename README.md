VMD Extension Functions
==============

This is a collection of TCL-VMD functions that support extraction of
structural data from large-scale simulations. These functions are
currently meant for TCL-VMD programmers.  Features easy semantics to

 * Iterate a block of code over frames
 * Iterate a block of code over trajectory files
 * Compute the number and fraction of native contacts
 * Compute distance matrices
 * ...and more

See the full documentation at http://tonigi.github.io/vmd_extensions/


Citation
--------

If you use this software, please cite:

> Giorgino T. Analysis Libraries for Molecular Trajectories: A Cross-Language Synopsis. In: Bonomi M, Camilloni C, editors. Biomolecular Simulations: Methods and Protocols [Internet]. New York, NY: Springer; 2019 [cited 2020 Jun 29]. p. 503–27. (Methods in Molecular Biology). Available from: https://doi.org/10.1007/978-1-4939-9608-7_20 [Preprint](https://github.com/giorginolab/analysis_libraries_chapter)


```
@incollection{giorgino_analysis_2019,
	address = {New York, NY},
	series = {Methods in {Molecular} {Biology}},
	title = {Analysis {Libraries} for {Molecular} {Trajectories}: {A} {Cross}-{Language} {Synopsis}},
	isbn = {978-1-4939-9608-7},
	shorttitle = {Analysis {Libraries} for {Molecular} {Trajectories}},
	url = {https://doi.org/10.1007/978-1-4939-9608-7_20},
	abstract = {Analyzing the results of molecular dynamics (MD)-based simulations usually entails extensive manipulations of file formats encoding both the topology (e.g., the chemical connectivity) and configurations (the trajectory) of the simulated system. This chapter reviews a number of software libraries developed to facilitate interactive and batch analysis of MD results with scripts written in high-level, interpreted languages. It provides a beginners’ introduction to MD analysis presenting a side-by-side comparison of major scripting languages used in MD and shows how to perform common analysis tasks within the Visual Molecular Dynamics (VMD), Bio3D, MDTraj, MDAnalysis, and High-Throughput Molecular Dynamics (HTMD) environments.},
	language = {en},
	urldate = {2020-06-29},
	booktitle = {Biomolecular {Simulations}: {Methods} and {Protocols}},
	publisher = {Springer},
	author = {Giorgino, Toni},
	editor = {Bonomi, Massimiliano and Camilloni, Carlo},
	year = {2019},
	doi = {10.1007/978-1-4939-9608-7_20},
	keywords = {Bio3D, HTMD, MDAnalysis, MDTraj, Molecular dynamics, Scripting languages, Trajectory analysis, VMD},
	pages = {503--527}
}

```
