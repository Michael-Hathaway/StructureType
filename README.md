<h1> StructureType Module </h1>
<p> Project for Hendrix Lab | Oregon State University </p>

<h3>About the StructureType Module</h3>
<p> The StructureType module is a python module that provides a user friendly mechanism for working with
RNA structure type files in the python programming language. </p>

<h3>Turner Parameters</h3>
<p>Source: https://rna.urmc.rochester.edu/NNDB/turner04/index.html</p>
<p>These are the current set of nearest neighbor parameters for RNA folding compiled by the Turner group. Both free energy changes at 37 ºC and enthalpy changes have been estimated, allowing for structure prediction at arbitrary temperature. These parameters are used by the StructureType module to calculate free energy values for the RNA molecules</p>
<p>The TurnerParameters directory contains three subdirectories: parameterTextFiles, scripts, and parameters. The parameterTextFiles are the text file found at the url above and provide the parameter value for RNA folding. The scripts directory contains several python scripts used to parse these parameters into python dictionaries and write them to python files so they can be imported by the StructureType Module. The parameters directory contains the .py files produced by the scripts.</p>
