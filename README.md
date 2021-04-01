# 2MMN40---Molecular-Modelling
A simple but fairly powerful Molecular Dynamics simulator in Python, developed for the course 2MMN40 at Eindhoven University of Technology (hence the 'week' structure; everything for the simulator can be found in 'Final Project' and 'FINAL_CODE'). The result consists of a topology file generator for water and ethanol, simulation, and analysis of results. The project was of limited length; there is plenty of room for improvement, as described in the report. 


To run the simulation, do the following: 
- Run 'GenerateTopology.py' for the desired substance and box size
- Enter the substance and box size (either using the function 'setSimulation' or manually) and run simulation.py for the desired simulation time
- To analyse energy levels, use 'energyAnalysis.py'; for plots of radial distribution functions, use 'RDF.py'

Output files are written to .xyz, which can be loaded into VMD (http://www.ks.uiuc.edu/Research/vmd/) to visualise the result. 
