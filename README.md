<h1> bpRNAStructure Module </h1>
<p> Project for Hendrix Lab | Oregon State University </p>

<h3>Usage</h3>
<p>To use this module clone this directory onto your machine and add the path to the directory to your PYTHONPATH environment varibale. This can be done by adding the following lines to your .bashrc or .zshrc file:
<ul>
    <li>export PYTHONPATH="${PYTHONPATH}:/path/to/bpRNAStructure"</li>
    <li>export PYTHONPATH="${PYTHONPATH}:/path/to/bpRNAStructure/TurnerParameters/parameters"</li>
</ul>
The Structure and StructureComponent modules can them be imported using:
<ul>
    <li>'import Structure as ST'</li>
    <li>'import StructureComponents as STC'</li>
</ul>
</p>

<h3>About the bpRNAStructure Module</h3>
<p> The bpRNAStructure Module is a python module that provides a user friendly mechanism for working with
RNA structure type files in the python programming language. </p>

<h4>Structure.py</h4>
<p>This file defines the Structure object and includes functionality for parsing the Structure Type file, as well as for accessing all the information stored in it</p>

<h4>Structure Components.py</h4>
<p>This file defines classes for all the secondary structures that are characterized in the Structure Type file. These secondary structures include: Stems, Bulges, Hairpins, InnerLoops, MultiLoops, ExternalLoops, PseudoKnots, Ends, and NCBPs. Each class provides specific functionality for accessing the information about each structure, as well as functionality for calculating the energy associated with each structure.</p>

<h3>Turner Parameters</h3>
<p>Source: https://rna.urmc.rochester.edu/NNDB/turner04/index.html</p>
<p>These are the current set of nearest neighbor parameters for RNA folding compiled by the Turner group. Both free energy changes at 37 ºC and enthalpy changes have been estimated, allowing for structure prediction at arbitrary temperature. These parameters are used by the StructureType module to calculate free energy values for the RNA molecules</p>
<p>The TurnerParameters directory contains three subdirectories: parameterTextFiles, scripts, and parameters. The parameterTextFiles are the text file found at the url above and provide the parameter value for RNA folding. The scripts directory contains several python scripts used to parse these parameters into python dictionaries and write them to python files so they can be imported by the StructureType Module. The parameters directory contains the .py files produced by the scripts.</p>
