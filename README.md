# Greenopt

Why use greenopt?
-	Greenopt was designed to help decision-makers incorporate green infrastructure (GI) into their water resources management plans while considering conflicting management objectives such as: 
  o	minimizing nutrient, sediment, and heavy metal loads (referred to here as pollutants)
  o	minimize negative impacts of stormwater runoff such as erosion, and
  o	minimize capital costs and maintenance costs. 
-	More precisely, greenopt facilitates a comprehensive understanding of the tradeoffs between conflicting management objectives by depicting the extent to which different candidate management plans are able to satisfy such objectives.  
  o	In other words, the inflection points and tradeoffs in a solution set (of management plan options) represent valuable information for decision-makers. For instance, greenopt can help answer questions such as, “At what degree of GI implementation is reducing nutrient load most cost-effective?” Or, “On which hydrologic response units (HRUs) should GI be implemented to achieve highest loading reductions per unit GI?” 

Simulation / optimization for decision-support:
-	This decision-support tool simulates changes in runoff and pollutant load with implementation of GI to evaluate the impact of different candidate plans on the management objectives. 
-	It uses multi-objective optimization to generate a Pareto optimal set of candidate management plans. These management plans are considered Pareto optimal because a plan’s performance in one objective cannot be improved without degrading its performance in another objective (Coello Coello, 1999).

Green infrastructure defined for the purposes of this tool:
-	By green infrastructure (GI), we are referring to structural stormwater controls, such as rain gardens, bioretention cells, vegetative swales, infiltration trenches, green roofs, etc., which have been shown in many cases to reduce stormwater runoff and improve water quality.  

What other programs do you need to use greenopt?
-	Greenopt was designed to compliment EPA’s Watershed Management Decision Support Tool (WMOST). Both WMOST and greenopt promote integrated water resources management, the concept that achieving sustainable, efficient water use requires consideration of the interconnectivity of our water resources.
-	Greenopt relies on WMOST for its user-interface and data import capabilities. On the Excel-based WMOST interface, users can indicate their interest in different types of GI, locations for GI implementation (in terms of lumped HRUs), bounds on amounts of implementation, as well as other parameters. In terms of data import capabilities, WMOST hosts a database of watershed runoff and associated loading (TN, TP, TSS, Zn) timeseries for many watersheds in the northeast. 
  o	Such runoff and loading timeseries are from calibrated HSPF, SWAT, or in some cases SWMM models.
-	For simulation, greenopt uses the EPA’s Stormwater Management Model (SWMM) to simulate changes in runoff and pollutant loading in the watershed due to implementation of GI.
-	For optimization, greenopt uses the Borg Multi-Objective Evolutionary Algorithm (MOEA). The Borg MOEA is freely available for noncommercial use, but it must be requested here: http://borgmoea.org/

What are the differences between WMOST and greenopt? 
-	Single versus multi-objective optimization (respectively)
-	All infrastructure versus green infrastructure (respectively)

What is provided in my Github repo?
Python script files (.py) 
-	Routines in runprep.py, runopt.py
-	Functions in read.py, write.py, calcs.py, checks.py
-	Global parameters in par.py 
-	One option in greenopt uses the script file swmmtoolbox.py, developed by Tim Cera, to read the binary output file produced by SWMM.

