<img src="docs/README_pics/logo_fingernat.png" width="500" class="center" />


Welcome to fingeRNAt's README
===================

fingeRNAt is a software to calculate Structural Interaction Fingerprint (SIFt) in nucleic acids - ligands complexes.

[![CI (conda)](https://github.com/n-szulc/fingeRNAt/workflows/CI%20(conda)/badge.svg?branch=master)](https://github.com/n-szulc/fingeRNAt/actions?query=workflow%3A%22CI+%28conda%29%22)
[![Python 3.8](https://img.shields.io/badge/python-3.8-blue.svg)](https://www.python.org/downloads/release/python-360/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Project Status: Active - The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
![Last modified](https://img.shields.io/github/last-commit/n-szulc/fingeRNAt?color=green)

<!-- TOC START min:1 max:6 link:true asterisk:false update:true -->
- [Overview](#overview)
	- [What is the Structural Interaction Fingerprint (SIFt)?](#what-is-the-structural-interaction-fingerprint-sift)
	- [What are the Structural Interaction Fingerprint (SIFt) applications?](#what-are-the-structural-interaction-fingerprint-sift-applications)
- [Installation](#installation)
	- [Recommended method](#recommended-method)
	- [Manual installation](#manual-installation)
- [Usage](#usage)
	- [Quick start :zap:](#quick-start-zap)
	- [Parameters description](#parameters-description)
	- [Inputs](#inputs)
	- [Structural Interaction Fingerprint (SIFt) types](#structural-interaction-fingerprints-sifs-types)
	- [Interactions with inorganic ions](#interactions-with-inorganic-ions)
	- [User defined thresholds](#user-defined-thresholds)
	- [Outputs](#outputs)
		- [`FULL`](#full)
		- [`SIMPLE`](#simple)
		- [`PBS`](#pbs)
	- [Wrappers](#wrappers)
	- [Heatmap plotting](#heatmap-plotting)
	- [Usage examples](#usage-examples)
	- [Graphical User Interface](#graphical-user-interface)
	- [Debugging mode](#debugging-mode)
	- [Warnings/Errors](#warningserrors)
- [Frequently Asked Questions (FAQ)](#frequently-asked-questions-faq)
- [fingerDISt :straight_ruler:](#fingerdist-straight_ruler)
	- [Installation](#installation-1)
	- [Usage](#usage-1)
		- [Quick start :zap:](#quick-start-zap-1)
		- [Parameters description](#parameters-description-1)
		- [Inputs](#inputs-1)
		- [Distance Metrics](#distance-metrics)
		- [Outputs](#outputs-1)
		- [Usage examples](#usage-examples-1)
- [Documentation](#documentation)
- [Unit test](#unit-test)
- [Detection of interactions](#detection-of-interactions)
	- [Nucleic acid properties](#nucleic-acid-properties)
	- [Ligand properties](#ligand-properties)
	- [Molecular Interactions' Geometric Rules](#molecular-interactions-geometric-rules)
		- [1. Hydrogen bonds](#1-hydrogen-bonds)
		- [2. Halogen bonds](#2-halogen-bonds)
		- [3. Cation - anion interactions](#3-cation---anion-interactions)
		- [4. Pi - cation & 5. Pi - anion interactions](#4-pi---cation-&-5-pi---anion-interactions)
		- [6. Pi - stacking interactions](#6-pi---stacking-interactions)
		- [7. Ion-mediated interactions](#7-ion-mediated-interactions)
		- [8. Water-mediated interactions](#8-water-mediated-interactions)
		- [9. Lipophilic interactions](#9-lipophilic-interactions)
- [Contributors](#contributors)
- [Feedback](#feedback)
- [Acknowledgments](#acknowledgments)
- [How to cite](#how-to-cite)
- [License](#license)
<!-- TOC END -->

# Overview

fingeRNAt is a Python 3.8 script calculating Structural Interaction Fingerprint (**SIFt**) in nucleic acid - ligand complexes.

It detects different non-covalent interactions in the input complex and returns long binary string - **SIFt**, describing if particular interaction between each nucleic acid residue and ligand occurred or not. **SIFt** can be calculated for the following complexes:

<p align="center">
<img src="docs/README_pics/detected_interactions.png" width="900" />
</p>

fingeRNAt runs under Python 3.5 - 3.8 on Linux, Mac OS and Windows.

## What is the Structural Interaction Fingerprint (SIFt)?

Structural Interaction Fingerprint (**SIFt**) is a binary string, describing existence (1/0) of specified molecular interactions between all target's residues and the ligand (*Deng et al.*, 2004).

<p align="center">
<img src="docs/README_pics/SIFs.png" width="800" />
</p>

<br/>
<br/>

<p align="center">
<img src="docs/README_pics/SIFs_merging.png" width="500" />
</p>

<br/>

**SIFt** translate information about 3D interactions in the target-ligand complex into a string, where the respective bit in the fingerprint is set to 1 in case of detecting particular interaction, and to 0 otherwise.

Therefore, **the interactions are represented in a unified fashion, thus allowing for easy high throughput computational analysis**, as they provide a full picture of all interactions within the complex.

## What are the Structural Interaction Fingerprint (SIFt) applications?

* Machine Learning
* Clustering
* Comparison of interaction modes

# Installation

Recommended fingeRNAt usage is in conda environment.

## Recommended method

1. Install conda

  Please refer to [conda manual](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) and install conda version with Python 3.x according to your operating system.


2. Download fingeRNAt repository

      Manually - click on the green field `Clone or download`, then `Download ZIP`

      **or**

      Clone it into the desired location [requires [git](https://git-scm.com/downloads) installation]

      `git clone --depth=1 https://github.com/n-szulc/fingernat.git`

3. Restore conda environment

      `conda env create -f fingeRNAt/env/fingeRNAt_env.yml`

## Manual installation

Required dependencies are:

- Python 3.8
- openbabel 3.1.1
- numpy  
- pandas
- matplotlib
- rdkit
- tk
- sphinx

# Usage

## Quick start :zap:

To call fingeRNAt with example inputs:

```bash
conda activate fingernat

cd fingeRNAt

python code/fingeRNAt.py -r example_inputs/1aju_model1.pdb -l example_inputs/ligands.sdf
```

## Parameters description

fingeRNAt accepts the following parameters:

`-r` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; path to RNA/DNA structure; see -> [Inputs](#Inputs)

`[-l]` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;path to ligands file; see -> [Inputs](#Inputs)

`[-f]` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; optional Structural Interactions Fingerprint (SIFt) type;

 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;   available types are: `FULL` [default], &nbsp;&nbsp;`SIMPLE`, &nbsp;&nbsp;`PBS`

  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; see -> [SIFt types](#structural-interaction-fingerprints-sift-types)

`[-o]` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; optional path to save output

`[-h2o]`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; optional detection of water-mediated interactions; applies only to SIFt type `FULL`. If not passed, all columns containing information about water-mediated interactions are empty (`None`) (applies also for wrappers)

`[-dha]`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; optional Donor-Hydrogen-Acceptor angle calculation when detecting hydrogen bonds;

  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
see -> [1. Hydrogen Bonds](#1-hydrogen-bonds)

`[-print]`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; print detected interactions for each nucleic acid - ligand complex on screen

`[-detail]`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; generate an additional file with detailed data on detected interactions see -> [PyMOL visualization](#pymol-visualization)


`[-vis]`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; optional SIFt results heatmap visualization; see -> [Heatmap plotting](#heatmap-plotting)

`[-wrapper]` &nbsp;&nbsp;&nbsp;&nbsp; optional SIFt results wrapper, see -> [Wrappers](#Wrappers)


 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;   available types are: `ACUG`, &nbsp;&nbsp;`PuPy`, &nbsp;&nbsp;`Counter`

 `[-verbose]` &nbsp;&nbsp;&nbsp;&nbsp; provides additional information about performed calculations at the given moment

 `[-debug]` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; enters debug mode, see -> [Debbuging mode](#debugging-mode)

 `[-h]` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; show help message

## Inputs

1. `-r `: path to receptor - RNA/DNA structure
    - supported file types: pdb, mol2
    - **only 1 model** of RNA/DNA structure
         - if there are more models, you have to choose only one (e.g. manually delete remaining models)
    - **no ligands**
    - may contain water & ions

   &#x1F534; **Hydrogens need to be added**

2. `-l`: optional path to ligand - small molecule, RNA/DNA/LNA or protein
	- **small molecule ligands**
         - supported file types: sdf
         - possible multiple poses of ligands in one file
    - **RNA/DNA structure**
    	- supported file types: sdf
    	- possible multiple models of RNA/DNA structure
    	- **only** RNA/DNA chains
        	- no water, ions, ligands

   &#x1F535; if `-l` is not specified, fingeRNAt **will find all inorganic ions in the receptor file and treat them as ligands**, see -> [Interactions with inorganic ions](#interactions-with-inorganic-ions).

**Additional notes:**
- Input ligand molecules should have assigned desired protonation state and formal charges.
- In the receptor molecule, charges on the phosphate groups do not need to be assigned (fingeRNAt always treats OP1 and OP2 atoms as negatively charged anions).
- Please pay attention to sdf ligand files converted from pdbqt/mol2 files, if the *formal charges* are preserved in the sdf files.
-	All the missing ligands' hydrogens will be automatically added.

## Structural Interaction Fingerprint (SIFt) types

fingeRNAt allows to calculate the following SIFt types:

- `FULL`

    Calculates **nine non-covalent interactions** for each RNA/DNA residue - ligand  pair:

	* hydrogen bondings (HB)
	* halogen bondings (HAL)
	* cation - anion (CA)
	* Pi - cation (Pi\_Cation)
	* Pi - anion (Pi\_anion)
	* Pi - stacking (Pi\_Stacking) interactions
	* ion-mediated; distinguishes between:
		* Magnessium-mediated (Mg\_mediated)
		* Potassium-mediated (K\_mediated),
		* Sodium-mediated (Na\_mediated),
		* Other ion-mediated (Other\_mediated)
	* water-mediated (Water\_mediated); **only if `-hoh` parameter was passed**; otherwise this interaction is assigned as `None`
	* lipophilic (lipophilic\_mediated)

	&#x1F536;Returns nine 0/1 values for each residue.

- `SIMPLE`

    **Calculates distances** between each RNA/DNA residue and the ligand; returns 1 if the distance does not exceed the declared threshold (default = 4.0 &#8491;), 0 otherwise. Does not take into account distances between hydrogens or hydrogen - heavy atom.

	&#x1F536;Returns one 0/1 value for each residue.

- `PBS`

    Divides each RNA/DNA residue into three groups: **Phosphate, Base, Sugar**. Then, **for each group separately**, calculates distance to the ligand; returns 1 if the distance does not exceed the declared threshold (default = 4.0 &#8491;), 0 otherwise. Does not take into account distances between hydrogens or hydrogen - heavy atom.

	&#x1F536;Returns three 0/1 values for each residue.
    > **_NOTE:_** Only for RNA/DNA with canonical residues.

## Interactions with inorganic ions

It is possible to calculate contacts for nucleic acid - ions from the same pdb/mol2 file.

In such case, parameter `-l` should be ommited, which allows fingeRNAt to treat all inorganic ions as ligands and calculate SIFt for each residue - ion pair.

As the aforementioned interactions are detected based on  contact, **only SIFt type SIMPLE or PBS may be calculated**.

  **Usage example**
	`python code/fingeRNAt.py -r example_inputs/3d2v.mol2 -f PBS`


## User defined thresholds

All the default thresholds can be changed in `code/config.py`

## Outputs

Calculated SIFt are saved to tsv files - a simple text format similar to csv, except for the data being tab-separated instead of comma-separated.

If fingeRNAT was run without optional parameter `-o`, script will create `outputs/` directory in the working directory and save there the output in tsv format. Otherwise fingeRNAt will save outputs in the user-specified location.

Example outputs for SIFt types and their wrappers are available from `fingeRNAt/example_outputs`.

### `FULL`
<br/>

<p align="center">
<img src="docs/README_pics/full.png" width="900" />
</p>

Sample extract of output of running `python code/fingeRNAt.py -r example_inputs/1aju_model1.pdb -l example_inputs/ligands.sdf -h2o`

> **_NOTE:_** If fingeRNAt was called without `-h2o` parameter, all columns containing information about water-mediated interactions are empty (`None`) (applies also for wrappers; (see -> ['Parameters description'](#parameters-description))).

### `SIMPLE`

<p align="center">
<img src="docs/README_pics/SIMPLE.png" width="900" />
</p>

Sample extract of output of running `python code/fingeRNAt.py -r example_inputs/1aju_model1.pdb -l example_inputs/ligands.sdf -f SIMPLE`


<br/>

### `PBS`

<p align="center">
<img src="docs/README_pics/PBS.png" width="900" />
</p>

Sample extract of output of running `python code/fingeRNAt.py -r example_inputs/1aju_model1.pdb -l example_inputs/ligands.sdf -f PBS`

<br/>


## Wrappers

Calculated SIFt (of any type) can be wrapped, which allows for representing it in lower resolution.

The results for the SIFt calculations and all passed wrappers are saved to separate tsv files.
Multiple wrappers may be passed at once (comma-separated; see -> ['Usage examples'](#usage-examples)).

The following three types of wrappers are available:

- `ACUG`

	Wraps calculated results according to nucleotide. Provides information if particular kind of interaction between e.g. any adenine from RNA/DNA and ligand occurred.

	<p align="center">
	<img src="docs/README_pics/PBS_ACUG.png" width="900" />
	</p>

	Sample extract of output of running `python code/fingeRNAt.py -r example_inputs/1aju_model1.pdb -l example_inputs/ligands.sdf -f PBS -wrapper ACUG`


- `PuPy`

	Wraps calculated results according to nucleobase type (purine or pyrimidyne). Provides information if particular kind of interaction between e.g. any purine from RNA/DNA and ligand occurred.

	<p align="center">
	<img src="docs/README_pics/FULL_PuPy.png" width="900" />
	</p>

	Sample extract of output of running `python code/fingeRNAt.py -r example_inputs/1aju_model1.pdb -l example_inputs/ligands.sdf -wrapper PuPy`

> **_NOTE:_**  As `-h2o` parameter was not passed, the columns containing information about water-mediated interactions are empty (`None`) (see -> ['Parameters description'](#parameters-description)).

- `Counter`

	Counts total number of given interaction type.

	<p align="center">
	<img src="docs/README_pics/FULL_Counter.png" width="900" />
	</p>

	Sample extract of output of running `python code/fingeRNAt.py -r example_inputs/1aju_model1.pdb -l example_inputs/ligands.sdf -wrapper Counter`

> **_NOTE:_**  As -h2o parameter was not passed, the columns containing information about water-mediated interactions are empty (`None`) (see -> ['Parameters description'](#parameters-description))

## Heatmap plotting

All SIFt outputs can be visualized as heatmap and saved as png files with the same name as tsv output. For large SIFt, it may be necessary to manually change the figure size in fingeRNAt.py to generate the plot.


<p align="center">
<img src="docs/README_pics/heatmap_PuPy.png" width="1000" />
</p>

Heatmap for SIFt type `FULL` with wrapper `PuPy` obtained from running `python code/fingeRNAt.py -r example_inputs/1aju_model1.pdb -l example_inputs/ligands.sdf -vis -wrapper PuPy`


<p align="center">
<img src="docs/README_pics/heatmap_Counter.png" width="600" />
</p>

Heatmap for SIFt type `FULL` with wrapper `Counter` obtained from running `python code/fingeRNAt.py -r example_inputs/1aju_model1.pdb -l example_inputs/ligands.sdf -wrapper Counter -h2o -vis`

## PyMOL visualization

Dedicated PyMOL plugin was created to visualize detected interactions based on `-detail` outputs (see -> [Detail mode](#detail-mode)).

## Usage examples

- Calculate SIFt `FULL`, print detected interactions on screen and save the output in the default location with the default filename.

`python code/fingeRNAt.py -r example_inputs/1aju_model1.pdb -l example_inputs/ligands.sdf -print`

- Calculate SIFt `SIMPLE` and save the output in the user-declared location with the default filename.

`python code/fingeRNAt.py -r example_inputs/1aju_model1.pdb -l example_inputs/ligands.sdf -f SIMPLE -o /path/to/my_output/`

- Calculate SIFt `PBS`, see what is being calculated using verbose mode, and save the output with the default filename in the `outputs` directory together with heatmap.

`python code/fingeRNAt.py -r example_inputs/1aju_model1.pdb -l example_inputs/ligands.sdf -f PBS -verbose -vis`

- Calculate default SIFt `FULL` and save it's output along with table containing details on each detected interaction with the default filenames in the `outputs` directory.

`python code/fingeRNAt.py -r example_inputs/1aju_model1.pdb -l example_inputs/ligands.sdf -detail`

- Calculate default SIFt `FULL`, consider water-mediated interactions, and save it's output and three wrapped outputs with the default filenames in the `outputs` directory.

`python code/fingeRNAt.py -r example_inputs/1aju_model1.pdb -l example_inputs/ligands.sdf -h2o -wrapper ACUG,PuPy,Counter`

- Calculate SIFt for nucleic acid - inorganic ions from the same pdb/mol2 file (see -> [Interactions with inorganic ions](#interactions-with-inorganic-ions))

`python code/fingeRNAt.py -r example_inputs/3d2v.mol2 -f PBS`

## Graphical User Interface

To use Graphical User Interface (GUI), simply run

`python code/gui.py`

GUI is user-friendly and has all aforementioned functionalities.

<p align="center">
<img src="docs/README_pics/gui.png" width="600" />
</p>

## Detail mode

fingeRNAt has additional mode `-detail` which saves all detected interactions, for any SIFt type, to separate tsv file (with prefix `DETAIL_`).

The interactions saved in `-detail` mode follow the same convention - **each row contains information about interaction between given ligand - receptor**, thus allowing for their easy, high-throughput processing.

The abovementioned data serve also as input to dedicated Pymol plugin, created by us, to visualize all interactions detected by fingeRNAt (see -> [PyMOL visualization](#pymol-visualization)).

> **_NOTE:_** In case of ion- and water-mediated interactions, two rows will correspond to one such interaction: (i) ligand - ion/water molecule and (ii) ion/water molecule - RNA/DNA. Therefore in case (i) ion/water molecule is treated as the receptor and in case (ii) as the ligand.

> **_NOTE 2:_** We refer to ligands as structures from ligands' file `-l`, as water molecules and inorganic ions are supposed to be provided in the RNA/DNA pdb/mol2 file.

The following data are saved in `-detail` mode:

* **Ligand_name**
  * for ligands from sdf file: ligand's name
 * for inorganic ions: ion's name
 * for water molecules: HOH
* **Ligand_pose**
 * for ligands from sdf file: ligand's pose number (indexed from 1)
 * for inorganic ions/water molecules: 0
* **Ligand_occurrence_in_sdf**
 * if `-l` ligands' file was passed
   * for ligands from sdf file: ligand's occurrence from the beginning of sdf file
	* for inorganic ions/water molecules: ligand's occurrence (from the beginning of sdf file), for which they mediate the interaction with nucleic acid
 * if `-l` ligands' file was not passed (see -> [Interactions with inorganic ions](#interactions-with-inorganic-ions)): 0
* **Interaction**: interaction type
* **Ligand_Atom**
  * for ligands from sdf file: ligand's atom index
 * for inorganic ions/water molecules: residue number : chain
* **Ligand_X/Ligand_Y/Ligand_Z**: coordinates of ligand's/ion's/water molecule's atom
* **Receptor_Residue_Name/Receptor_Number/Receptor_Chain**:
 * for interactions ligand/inorganic ions/water molecule - nucleic acid: receptor's nucleotide name/receptor's residue number/receptor's chain
 * for interactions ligand - ion: ion's name/ion's residue number/ion's chain
 * for interactions ligand - water molecule: HOH/water molecule's residue number/water molecule's chain
* **Receptor_Atom**
 * for interactions ligand/inorganic ions/water molecule - nucleic acid: receptor's atom ID
 * for interactions ligand - ion: ion's name
 * for interactions ligand - water molecule: O
* **Receptor_X/Receptor_Y/Receptor_Z**: coordinates of ligand's/ion's/water molecule's atom
* **Distance**: distance between ligand's and residue's atoms [&#8491;]

**Sample output**
see -> [DETAIL_3d2v.mol2_redocked.sdf_FULL.tsv](example_outputs/DETAIL_3d2v.mol2_redocked.sdf_FULL.tsv)

## Debugging mode

fingeRNAt has debugging mode in which it prints on screen exhaustive information about detected interactions.

Debugging mode may be used with each SIFt type and provides the following information:

- `FULL`

  Prints the following **properties of each ligand**:
  1. atom indices of hydrogen bonds acceptors & donors
  2. atom indices of halogen bonds donors
  3. atom indices of cations & anions
  4. number of ligand's aromatic rings
  5. IDs of inorganic ions in electrostatic contact with ligand_name_detail
  6. IDs of water molecules in contact with ligand
  7. atom indices of lipophilic atoms

  Prints the following **properties for each residue** of nucleic acid:
  1. atom IDs of hydrogen bonds acceptors & donors
  2. atom IDs of anions
  3. atom IDs of aromatic rings

  Prints **detected interactions** for pair of each nucleic acid residue - ligand:
  1. atoms creating hydrogen bond with their distance and angle (if parameter `-dha` was passed), separately for cases when nucleic acid is hydrogen bond acceptor and hydrogen bond donor
  2. atoms creating halogen bond with their distance and angles
  3. atoms creating cation-anion interaction with their distance (note that nucleic acid's atoms are anions)
  4. atoms creating Pi-cation interaction with their distance and angle
  5. atoms creating Pi-anion interaction with their distance and angle (note that ligand's atoms are anions)
  6. atoms creating Pi-stacking interaction type Sandwich/Displaced with their distance, offset and angle
  7. atoms creating Pi-stacking interaction type T-shaped with their distance, offset and angle
	8. atoms creating ion-mediated interaction with their distances (nucleic acid - ion & ion - ligand)
	9. atoms creating water-mediated interaction with their distances (nucleic acid - water & water - ligand)
	10. atoms creating lipophilic interaction with their distance

- `SIMPLE`

    For each detected contact of nucleic acid residue - ligand prints information about atoms IDs and the distance between them.

- `PBS`

    For each detected contact of nucleic acid residue's group - ligand prints information about atoms IDs and the distance between them.

> **_NOTE:_** In all above cases, only the first detected interaction of given type is printed, as fingeRNAt stops searching for more once it detected one particular interaction.

## Warnings/Errors

Please pay attention to the following types of errors: **Could not sanitize molecule ending on line ...**. This means that RDKit library used by the fingeRNAt cannot properly read the molecule.

e.g.

```
[08:46:45] non-ring atom 39 marked aromatic
[08:46:45] ERROR: Could not sanitize molecule ending on line 168
[08:46:45] ERROR: non-ring atom 39 marked aromatic
```

- Error: **non-ring atom ... marked aromatic**
  - Solution: please make sure that the mentioned molecule has a proper aromatic ring representation.
- Error: **Explicit valence for atom # ..., is greater than permitted** (eg., *Explicit valence for atom # 18 O, 3, is greater than permitted*)
  - Solution: please make sure that the indicated atom(s) have a proper valence number, i.e., they form a correct number of bonds.

# Frequently Asked Questions (FAQ)

**What happens when I have a non-canonical nucleotide in my nucleic acid?**

* If you have a residue with **only** non-canonical name (all atom names are canonical), e.g. X

|   | `FULL`  | `SIMPLE`  | `PBS`  |
|:-:|:-:|:-:|:-:|
| No wrapper  |  OK |  OK  | OK   |
| `ACUG`  | Omits interaction for residue with non-canonical name  | Omits interaction for residue with non-canonical name  |  Omits interaction for residue with non-canonical name |
|  `PuPy` |  Omits interaction for residue with non-canonical name | Omits interaction for residue with non-canonical name  |  Omits interaction for residue with non-canonical name |
| `Counter`  |  OK |  OK | OK  |

* If you have a residue with a canonical name but with non-canonical atom name e.g. P9

|   | `FULL`  | `SIMPLE`  | `PBS`  |
|:-:|:-:|:-:|:-:|
| No wrapper  |  OK |  OK  | Does not work   |
| `ACUG`  |  OK   | OK  |  Does not work|
|  `PuPy` |   OK  | OK |  Does not work|
| `Counter`  |  OK |  OK | Does not work  |

Please note, we consider both oxygens from phosphate group (OP1 and OP2) of nucleic acid as negatively charged, therefore fingeRNAt will not consider differently named atoms as anions.

* If you have a residue with non-canonical name and non-canonical atom name e.g. P9

|   | `FULL`  | `SIMPLE`  | `PBS`  |
|:-:|:-:|:-:|:-:|
| No wrapper  |  OK |  OK  | Does not work   |
| `ACUG`  |  Omits interaction for residue with non-canonical name | Omits interaction for residue with non-canonical name  |  Does not work |
|  `PuPy` | Omits interaction for residue with non-canonical name  | Omits interaction for residue with non-canonical name  |  Does not work |
| `Counter`  |  OK |  OK | Does not work |

Please note,  we consider both oxygens from phosphate group (OP1 and OP2) of nucleic acid as negatively charged, therefore fingeRNAt will not consider differently named atoms as anions.

---

**What happens when I have nucleic acid with two residues with the same number e.g. due to errors in structure?**

In case of SIFt types `SIMPLE` and `PBS`, the only difference is that their outputs will have two columns with the same name in output. SIFt and the wrapped results are correct.

However in case of SIFt type `FULL`, there will be two columns with the same name but their Pi-interactions may be swapped and their SIFt as well as the wrapped results may be unreliable.

# fingerDISt :straight_ruler:

<img src="docs/README_pics/logo_fingerdist.png" width="300" class="center" />

fingerDISt is an additional, standalone tool which calculates different Distance Metrics based on Structural Interaction Fingerprint (SIFt) - outputs of fingeRNAt.

It calculates the selected Distance Metric for all SIFt vs. all SIFt from input file - creates a matrix of scores and saves it to tsv file.

## Installation

fingerDISt, similarly like fingeRNAt, requires Python 3.5 - 3.8, and may be run from within fingeRNAt's environment, but it is not obligatory. No external modules are needed.

## Usage

### Quick start :zap:

```bash
cd fingeRNAt

python code/fingerDISt.py -i example_outputs/FULL/1aju_model1.pdb_ligands.sdf_FULL.tsv -m tanimoto
```

### Parameters description

fingerDISt accepts the following parameters:

`-i` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; path to tsv/csv file with calculated SIFt; see -> [Inputs](#Inputs-1)

`-m` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; types of desired Distance Metrics; see -> [Distance Metrics](#distance-metrics)

`[-o]` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; optional path to save output

`[-verbose]` &nbsp;&nbsp;&nbsp;&nbsp; prints calculated Distance Metrics on the screen

`[-h]` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; show help message

### Inputs

1. `-i `: path to tsv/csv file with calculated SIFs

  **fingeRNAt outputs are fingerDISt inputs**

  Example inputs may be found in `fingeRNAt/tests/expected_outputs`.

### Distance Metrics

fingerDISt calculates the following Distance Metrics:

* Tanimoto coefficient
* Cosine similarity
* Manhattan
* Euclidean
* Square Euclidean
* Half Square Euclidean
* Soergel Distance

The abovementioned Distance Metrics calculations were implemented based on code [github.com/varun-parthasarathy/crux-fr-sprint/blob/master/DistanceMetrics.py](github.com/varun-parthasarathy/crux-fr-sprint/blob/master/DistanceMetrics.py) under MIT license.

  > **_NOTE:_** **Tanimoto coefficient works only for SIFt with binary values**, therefore it may not work on input SIFt wrapped with `Counter` wrapper.


> **_NOTE 2:_** It automatically replaces `None` with 0, meaning that Distance Metrics can be calculated for SIFt type `FULL`, which was called without `-h2o` parameter.

### Outputs

fingerDISt saves scores for each selected Distance Metric to separate tsv files - a simple text format similar to csv, except for the data being tab-separated instead of comma-separated.

If fingerDISt was run without optional parameter `-o`, script will create `outputs/` directory in the working directory and save there the output in tsv format. Otherwise fingerDISt will save outputs in the user-specified location.

<p align="center">
<img src="docs/README_pics/Tanimoto.png" width="900" />
</p>

Sample output of running `python code/fingerDISt.py -i tests/expected_outputs/1aju_model1.pdb_ligands.sdf_FULL.tsv -m tanimoto`

### Usage examples

- Calculate all available Distance Metrics on SIFs inputs type `FULL` and save the output with the default filename in the `outputs` directory.

`python code/fingerDISt.py -i tests/expected_outputs/1aju_model1.pdb_ligands.sdf_FULL.tsv -m manhattan,square_euclidean,euclidean,half_square_euclidean,cosine_similarity,tanimoto`

- Calculate two Distance Metrics on SIFs inputs type `PBS` wrapped with `ACUG` wrapper and save the output to user-specified location with the default filename.

`python code/fingerDISt.py -i tests/expected_outputs/1aju_model1.pdb_ligands.sdf_PBS_ACUG.tsv -m manhattan,square_euclidean -o my_dir/`

- Calculate one Distance Metric on SIFs inputs type `SIMPLE`, print it on the screen and save the output to the user-specified location with the desired filename.

`python code/fingerDISt.py -i tests/expected_outputs/1aju_model1.pdb_ligands.sdf_SIMPLE.tsv -m tanimoto -verbose -o my_dir/my_filename`


# Documentation

To generate fingeRNAt documentation file using sphinx:

```bash
cd docs
make html
```

The documentation will be available from `_build/html`.

# Unit test

To run a unit test:

```bash
cd tests
python fingeRNAt_test.py
```

# Detection of interactions

Inspired by [PLIP](https://github.com/pharmai/plip) implementation. Applies to SIFt type `FULL`.

## Nucleic acid properties

The nucleic acid's detected properties are as follows:

* Hydrogen Bonds Acceptors & Donors - detected with OpenBabel
* Halogen Bonds Acceptors - detected with OpenBabel (same as Hydrogen Bonds Acceptors)
* Negative charges - assigned to *both* oxygens (OP1 & OP2) of each residue's phosphate group
* Aromatic rings -  detected with OpenBabel

## Ligand properties

The ligand's detected properties are as follows:

* Hydrogen Bonds Acceptors & Donors - detected with OpenBabel
* Halogen Bonds Donors - detected with OpenBabel
* Positive & Negative charges - detected with OpenBabel
* Aromatic rings - detected with RDKit
* Lipophilic atoms - detected with RDKit

## Molecular Interactions' Geometric Rules

fingeRNAt detects the following nine non-covalent interactions:

<img src="docs/README_pics/Interactions.png" width="900"/>


## 1. Hydrogen bonds

  Geometric rule:
  * |D - A| < 3.9 &#8491;
	(*Torshin, Weber, & Harrison*, 2002)


  > **_NOTE:_** If hydrogens are present in RNA/DNA structure, fingeRNAt can be run with parameter `-dha`, that additionaly calculates Donor-Hydrogen-Acceptor angle  used as supplementary criteria in hydrogen bonds detection:&nbsp;&nbsp;&nbsp;&nbsp;
  100&deg; <   D-H-A angle < 260&deg; <br/>
  Parameter `-dha` applies only to `FULL` SIFt type, as `SIMPLE` & `PBS` do not calculate hydrogen bonds.

## 2. Halogen bonds

  Geometric rules:
	- |X - A| < 4.0 &#8491;
	(*Auffinger et al.*, 2004)
	- C-X-A angle ~ 165&deg; &#177; 30&deg;
	(*Auffinger et al.*, 2004)
	- X-A-A' angle ~ 120&deg; &#177; 30&deg;
	(*Auffinger et al.*, 2004)

## 3. Cation - anion interactions

  Geometric rule:
  - 0.5 &#8491; < |cation - anion| < 5.5 &#8491;
  (*Barlow and Thornton*, 1983)

  > **_NOTE:_** fingeRNAt considers both oxygens from phosphate group (OP1 and OP2) of RNA/DNA as negatively charged.

## 4. Pi - cation** & 5. Pi - anion interactions

  Geometric rules:

- |cation/anion - aromatic ring center| < 6.0 &#8491;
(*Gallivan and Dougherty*, 1999)
- angle between the ring plane and the line between cation/anion - ring center ~ 90&deg; &#177; 30&deg;

## 6. Pi - stacking interactions

  Geometric rules:

  * Common rules for all Pi - stacking interactions' types:

    * |rings' centroids| < 5.5 &#8491;
		(*McGaughey*, 1998)
    * rings' outset < 2.0 &#8491;

  * For sandwich & parallel - displaced:
    * angle between the ring planes < 30&deg;

  * For T - shaped:
    * angle between the ring planes ~ 90&deg; &#177; 30&deg;

> **_NOTE:_** fingeRNAt considers all three abovementioned Pi - stacking interactions' types.

##  7. Ion-mediated interactions

  Geometric rules:
	* |ligand's nitrogen/oxygen/sulphur - ion| <= X
	* |nucleic acid's nitrogen/oxygen - ion| <= X
	where:
	  * X = 3.2 &#8491; for magnesium ion
	  (*Zheng et al.*, 2015)
	  * X = 3.9 &#8491; for potassium ion
		(*Zheng et al.*, 2008)
		* X = 3.6 &#8491; for sodium ion
		(*Zheng et al.*, 2008)
	  * X = 3.5 &#8491; for other ions
		(*Zheng et al.*, 2008)

## 8. Water-mediated interactions

> **_NOTE:_**  Only if fingeRNAt was called with `-h2o` parameter.

  Geometric rules:
  * |ligand's hydrogen bond donor/acceptor - water (oxygen)| <= 3.5 &#8491;
	(*Poornima & Dean*, 1995)
  * |nucleic acid's hydrogen bond donor/acceptor - water (oxygen)| <= 3.5 &#8491;
	(*Poornima & Dean*, 1995)

## 9. Lipophilic interactions

  Geometric rule:
  * |nucleic acid's carbon - ligand's lipophilic atom| <= 4.0 &#8491;
	(*Padroni et al.*, 2020)

# Contributors

Natalia Szulc, @n-szulc ![](https://img.shields.io/badge/nszulc-%40iimcb.gov.pl-brightgreen)

Filip Stefaniak, @filipsPL ![](https://img.shields.io/badge/fstefaniak-%40genesilico.pl-brightgreen)

# Feedback

We welcome any feedback, please send an email to Natalia Szulc  ![](https://img.shields.io/badge/nszulc-%40iimcb.gov.pl-brightgreen)


# Acknowledgments

Special thanks of gratitude to [Masoud Farsani](https://github.com/mafarsani), [Pritha Ghosh](https://github.com/prithaghosh) and Tomasz Wirecki for their invaluable feedback as well as to Prof. Janusz M. Bujnicki and the entire [Bujnicki Lab](http://genesilico.pl) for all the support and project guidelines.

Extensive script testing provided by Zuzanna Mackiewicz has been a great help in developing this tool.

Assistance provided by [Open Babel Community](http://openbabel.org) was greatly appreciated.


# How to cite

If you use this software, please cite:

*fingeRNAt - a novel tool for high-throughput analysis of nucleic acid-ligand interactions*

Natalia A. Szulc, Zuzanna Mackiewicz, Janusz M. Bujnicki, Filip Stefaniak
[in preparation]

# License

fingeRNAt is licensed under the GNU General Public License v3.0.
