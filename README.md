#GHSR_archive

This archive provides python scripts and partially processed data to produce Figs. 1 and 3-7 in the paper:

"Analysis of the dynamics of the human growth hormone secretagogue receptor reveals insights into the energy landscape of the molecule"

by

Albert A. Smith, Emelyne M. Pacull, Sabrina Stecher, Peter W. Hildebrand, Alexander Vogel, Daniel Huster 

The following scripts can be run at the command line to produce components of figures in the paper and SI. These will also produce the source data files.
Fig1_overall.py,
Fig3_CACB.py,
Fig4_ghrelinCC.py,
Fig5_aromatics.py,
Fig6_PCA.py,
Fig7_PCA_path.py,
SI_Fig2_CCcount.py,
SI_Fig4_PCA_path.py

The following scripts are provided for transparency in the data processing, but cannot be run without access to the full MD trajectories. 
These scripts create the five project folders found in Projects/, and also the folder PCA_results
setup_direct.py,
setup_directCACB.py,
setup_iRED_aromatic.py,
setup_iRED_ghrelin.py,
setup_PCA.py

There is NO INSTALLATION required for the provided code. Just place everything in a folder, navigate there, and run with python3.
However, python3 and the following modules must be installed from other sources (these are the tested versions, although other versions may work).

Python v. 3.7.3,
numpy v. 1.19.2,
scipy v. 1.5.2,
MDAnalysis v. 1.0.0,
matplotlib v. 3.4.2

Additionally, ChimeraX must be installed and its location provided to pyDR. How to do this is provided at a comment at the top of the provided scripts.

We recommend installing Anaconda: https://docs.continuum.io/anaconda/install/
The Anaconda installation includes Python, numpy, scipy, pandas, and matplotlib. 
(I also highly recommend using Spyder, which comes with Anaconda, for running the provided scripts interactively, such that one may stop to understand each step in the overall analysis)

MDAnalysis is installed by running:
conda config --add channels conda-forge
conda install mdanalysis
(https://www.mdanalysis.org/pages/installation_quick_start/)

All files are copyrighted under the GNU General Public License. A copy of the license has been provided in the file LICENSE

Copyright 2023 Albert Smith-Penzel

This research was supported by Deutsche Forschungsgemeinschaft (DFG) through CRC 1423, project number 421152132, subproject A02 and through DFG grant 450148812 (awarded to Albert A. Smith)
