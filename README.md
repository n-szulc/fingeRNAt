
<img src="docs/README_pics/logo_fingernat.png" width="500" class="center" />


# Welcome to fingeRNAt's README

fingeRNAt is a software to calculate Structural Interaction Fingerprints in nucleic acids - ligands complexes.

[![CI (conda)](https://github.com/n-szulc/fingeRNAt/workflows/CI%20(conda)/badge.svg?branch=master)](https://github.com/n-szulc/fingeRNAt/actions?query=workflow%3A%22CI+%28conda%29%22)

<!-- TOC depthFrom:1 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 -->

- [Overview](#overview)
- [Installation](#installation)
	- [Recommended installation instructions](#recommended-installation-instructions)
	- [Manual installation](#manual-installation)
- [Usage](#usage)
	- [Quick start](#quick-start)
	- [Parametres description](#parametres-description)
	- [Inputs](#inputs)
	- [Structural Interactions Fingerprints' (SIFs) types](#structural-interactions-fingerprints-sifs-types)
	- [Molecular Interactions' Geometric Rules](#molecular-interactions-geometric-rules)
		- [1. Hydrogen Bonds](#1-hydrogen-bonds)
		- [2. Halogen Bonds](#2-halogen-bonds)
		- [3. Cation - Anion](#3-cation---anion)
		- [4. Pi - Cation & 5. Pi - Anion](#4-pi---cation-5-pi---anion)
		- [6. Pi - Stacking](#6-pi---stacking)
	- [Defining own thresholds](#defining-own-thresholds)
	- [Outputs](#outputs)
		- [`FULL`](#full)
		- [`SIMPLE`](#simple)
		- [`PBS`](#pbs)
		- [`XP`](#xp)
	- [Wrappers](#wrappers)
	- [Visualization](#visualization)
	- [Usage examples](#usage-examples)
	- [Documentation](#documentation)
	- [Unit test](#unit-test)
	- [GUI](#gui)
- [Feedback](#feedback)
- [Acknowledgments](#acknowledgments)
- [Citing](#citing)
- [License](#license)

<!-- /TOC -->

# Overview

fingeRNAt is a Python 3.8 script calculating Structural Interactions Fingerprints (SIFs) in complexes of:

| Nucleic acid |Ligand|
|:---:|:---:|
| RNA | small molecule ligand |
| RNA | RNA  |
| RNA | protein |
| DNA | small molecule ligand |
| DNA | DNA  |
| DNA | protein |


fingeRNAt calculates different interaction between input RNA/DNA structure and ligand, returns long binary string describing if particular interaction type occurred between given nucleic acid residue and ligand or not.

fingeRNAt runs under Python 3.5 - 3.8 on Linux, Mac OS and Windows.

# Installation

Recommended fingeRNAt usage is in conda environment.

## Recommended installation instructions

1. Install conda

  Please refer to [conda manual](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) and install conda version according to your operating system. Please use Python3 version.


2. Download fingeRNAt repository

      Manually - click on the green field `Clone or download`, then `Download ZIP`
      
      **or** 
      
      Clone it into the desired location [requires [git](https://git-scm.com/downloads) installation] `git clone https://github.com/n-szulc/fingernat.git`

3. Restore conda environment

      `conda env create -f fingeRNAt/env/fingeRNAt_env.yml`


## Manual installation

Required dependencies are:

- Python 3.8
- openbabel 3.1.1
- numpy  
- pandas
- matplotlib
- tk
- sphinx


# Usage

## Quick start

To call fingeRNAt with example inputs:

```bash

conda activate fingernat

cd fingeRNAt

python code/fingeRNAt.py -r example_inputs/1aju_model1.pdb -l example_inputs/ligands.sdf -f SIMPLE
```

## Parametres description

<br/>
where:

`-r` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; path to RNA/DNA structure; see -> [Inputs](#Inputs)

`-l` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; path to ligands file; see -> [Inputs](#Inputs)

`[-f]` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; optional Structural Interaction Fingerprint (SIFt) type;

 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;   available types are: `FULL` [default], &nbsp;&nbsp;`SIMPLE`, &nbsp;&nbsp;`PBS`, &nbsp;&nbsp;`XP`

  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; see -> [SIFs types](#structural-interactions-fingerprints-sifs-types)

`[-o]` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; optional output file name

`[-dha]`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; optional donor-hydrogen-acceptor angle calculation when detecting hydrogen bonds; see -> [1. Hydrogen Bonds](#1-hydrogen-bonds)

`[-vis]`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; optional calculated Structural Interaction Fingerprints results heatmap visualization; see -> [Visualization](#Heatmap-visualization)

`[-wrapper]` optional calculated Structural Interaction Fingerprints results wrapper, see -> [Wrappers](#Wrappers)

 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;   available types are: `ACUG`, &nbsp;&nbsp;`PuPy`, &nbsp;&nbsp;`Counter`

 `[-h]` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; show help message

## Inputs

1. `-r `: path to RNA/DNA structure
    - supported file types: pdb, mol2
    - **only 1 model** of RNA/DNA structure
         - if there are more models, you have to choose only one (e.g. manually delete remaining models)
    - **only** RNA/DNA chains
        - no water, ions, ligands
    - optionally with hydrogens
        - if RNA/DNA structure was obtained from NMR, hydrogens are already there
        - if RNA/DNA structure was obtained from XR or cryo-EM, hydrogens may be added (e.g. in PyMol, VMD, Chimera)
2. `-l`: path to small molecule ligands **OR** RNA/DNA structure
	- **small molecule ligands**
         - supported file types: sdf, mol2, pdb for small molecule ligands
         - possible multiple ligands in multiple poses in one file
    - **RNA/DNA structure**
    	- supported file types: pdb, mol2
    	- possible multiple models of RNA/DNA structure
    	- **only** RNA/DNA chains
        	- no water, ions, ligands
    	- optionally with hydrogens
        	- if RNA/DNA structure was obtained from NMR, hydrogens are already there
        	- if RNA/DNA structure was obtained from XR or cryo-EM, hydrogens may be added (e.g. in PyMol, VMD, Chimera)


## Structural Interactions Fingerprints' (SIFs) types

*Structural interactions fingerprint (SIFt)* is a binary string, describing existence (1/0) of specified molecular interactions between all RNA/DNA residues and ligand (small molecule compound or another RNA/DNA structure).

<img src="docs/README_pics/SIFs.png" width="1000" />

<br/>
<br/>

 Each 1AJU residue has calculated six different molecular interactions with it's small molecule ligand.

**Possible SIFs types `[-f]`**

- `FULL`

    Calculates **six different molecular interactions** for each RNA/DNA residue - ligand (small molecule compound or another RNA/DNA structure) pair: hydrogen bondings (HB), halogen bondings (HAL), cation - anion (CA), Pi - cation (Pi\_Cation), Pi - anion (Pi\_anion) & Pi - stacking (Pi\_Stacking) interactions; returns six 0/1 values for each RNA residue.

- `SIMPLE`

    Calculates **distances** between each RNA/DNA residue and ligand (small molecule compound or another RNA/DNA structure); returns 1 if the distance is less than declared CUTOFF (default = 4.0 A), 0 otherwise.

- `PBS`

    Divides each RNA/DNA residue in **3 groups (Phosphate, Base, Sugar)** and for each group **calculates distance** to the ligand (small molecule compound or another RNA/DNA structure); returns three 0/1 values for each group within RNA/DNA residue - 1 if the distance is less than declared CUTOFF (default = 4.0 A), 0 otherwise.

    > **_NOTE:_** Only for RNA/DNA with canonical residues.

- `XP`

    Calculates the same six different molecular interactions for each RNA/DNA residue as `FULL`, however it is of no binary type - it **calculates total number of each potential interactions occurrence** (except Pi - interactions) for each RNA/DNA residue - ligand (small molecule compound or another RNA/DNA structure) pair, therefore being an **extra precision hologram**.

		> **_NOTE:_** It returns total number of potential interactions between given residue and ligand pair, e.g. if the residue has one hydrogen bond donor and the ligand has two hydrogen bond acceptors, both fulfilling hydrogen bonding geometrical rules (see -> [Hydrogen Bonds](###1.Hydrogen Bonds)), XP will return 2 for the given residue-ligand pair, despite the fact that one hydrogen bond donor may interact with only one hydrogen bond acceptor.

    - In case of hydrogen bonds, it not only calculates total number of its potential occurrence in each RNA/DNA - ligand pair, but **also assigns each hydrogen bond as strong/moderate/weak type** and calculates total number of each type potential occurrence in each RNA/DNA - ligand pair.

    - **It does not calculate total number of potential Pi - interactions** due to considering both purine's rings separately - if it calculated total numbers of Pi - interactions, it would calculate RNA/DNA purine - ligand's aromatic ring interaction as two independent interactions, which would not be true.

## Molecular Interactions' Geometric Rules

### 1. Hydrogen Bonds

<img src="docs/README_pics/hb.png" width="250" alt="Torshin, Weber, & Harrison, 2002" />

<br/>
<br/>

**Geometric rules:**

- D - A distance < 3.9 &#8491;

> **_NOTE:_**
If hydrogens are present in RNA/DNA structure and in ligand, fingeRNAt can be run with additional flag `-dha`, that additionaly calculates angle between donor, hydrogen and acceptor and also takes it's value into account in hydrogen bonds calculation
- 100&deg; < D - H - A angle < 260&deg;
Applies only to FULL/XP SIFt type, as SIMPLE & PBS do not calculate hydrogen bonds.


(*Torshin, Weber, & Harrison*, 2002)

In case of `XP` hologram, there is additional assignment of each hydrogen bond type. Depending on D - A distance, each hydrogen bond can be assigned as strong/moderate/weak (*Jeffrey*, 1997).

- 2.2 &#8491; < D - A distance < 2.5 &#8491;: strong
- 2.5 &#8491; < D - A distance < 3.5 &#8491;: moderate, mostly electrostatic
- 3.5 &#8491; < D - A distance < 4.0 &#8491;: weak, electrostatic

### 2. Halogen Bonds

<img src="docs/README_pics/hal.png" width="300" alt="Auffinger et al., 2004" />

<br/>
<br/>

**Geometric rules:**


- X - O distance < 4.0 &#8491;
- C - X - O angle ~ 165&deg; &#177; 30&deg;
- X - O - Y angle ~ 120&deg; &#177; 30&deg;

(*Auffinger et al.*, 2004)

### 3. Cation - Anion

<img src="docs/README_pics/ca.png" width="200" alt="Barlow and Thornton, 1983" />

<br/>
<br/>

**Geometric rule:**

- 0.5 &#8491; < cation - anion distance  < 5.5 &#8491;

(*Barlow and Thornton*, 1983)

### 4. Pi - Cation & 5. Pi - Anion

<img src="docs/README_pics/pi-ion.png" width="200" alt="Wikimedia Commons, modified" />

<br/>
<br/>

**Geometric rules:**

- cation/anion - aromatic ring center distance < 6.0 &#8491; (*Gallivan and Dougherty*, 1999)
- angle between the ring plane and the line between cation/anion - ring center ~ 90&deg; &#177; 30&deg;

### 6. Pi - Stacking

<img src="docs/README_pics/pi-stacking.png" width="500" alt="Wikimedia Commons, modified" />

<br/>
<br/>

> **_NOTE:_** All above interactions' types are considered by fingeRNAt.

**Geometric rules:**

Common rules for all Pi - stacking interactions' types:

- ring centroids distance < 5.5 &#8491; (*McGaughey*, 1998)
- rings' outset < 2.0 &#8491;

For Sandwich & Parallel - displaced:
- angle between the ring planes < 30&deg;

For T - shaped:
- angle between the ring planes ~ 90&deg; &#177; 30&deg;

## Defining own thresholds

All the default thresholds can be changed in `code/config.py`

## Outputs

If fingeRNAT was run without optional parameter `-o`, script will create `outputs/` directory in the working directory and save there the output in tsv format.

### `FULL`

<img src="docs/README_pics/full-explanation.png" width="900" />

Sample extract of output of running `python code/fingeRNAt.py -r example_inputs/1aju_model1.pdb -l example_inputs/ligands.sdf`

<br/>


### `SIMPLE`

<img src="docs/README_pics/simple-explanation.png" width="900" />

Sample extract od output of running `python code/fingeRNAt.py -r example_inputs/1aju_model1.pdb -l example_inputs/ligands.sdf -f SIMPLE`


<br/>

### `PBS`

<img src="docs/README_pics/pbs-explanation.png" width="900" />

Sample extract of output of running `python code/fingeRNAt.py -r example_inputs/1aju_model1.pdb -l example_inputs/ligands.sdf -f PBS`

<br/>

### `XP`

<img src="docs/README_pics/xp-explanation.png" width="1000" />

Sample extract of output of running `python code/fingeRNAt.py -r example_inputs/1aju_model1.pdb -l example_inputs/ligands.sdf -f XP`


## Wrappers

Calculated SIFs of all 4 types can be wrapped, representing SIFs in decreasing resolutions. Multiple wrappers, comma-separated, may be passed at once. The results for the fingerprint calculations and all passed wrappers are saved as separate tsv files.

There are 3 types of wrappers:

- `ACUG`

	Wraps calculated results according to nucleotide, gives information if particular kind of interaction between e.g. any adenine from RNA/DNA and ligand occurred (SIFt types: `SIMPLE`, `PBS`, `FULL`) or returns number of these interactions with all adenines (SIFt type `XP`; see -> [`XP`](#xp)).

	<img src="docs/README_pics/acug_full.png" width="900" />

	Sample extract of output of running `python code/fingeRNAt.py -r example_inputs/1aju_model1.pdb -l example_inputs/ligands.sdf -wrapper ACUG`

	<br/>

	<img src="docs/README_pics/acug_xp.png" width="1000" />

	Sample extract of output of running `python code/fingeRNAt.py -r example_inputs/1aju_model1.pdb -l example_inputs/ligands.sdf -f XP -wrapper ACUG`

- `PuPy`

	Wraps calculated results according to nucleobase type (purine or pyrimidyne), gives information if particular kind of interaction between e.g. any purine from RNA/DNA and ligand occurred (`SIMPLE`, `PBS`, `FULL`) or returns number of these interactions with all purines (`XP`; see -> [`XP`](#xp)).

	<img src="docs/README_pics/pupy_full.png" width="900" />

	Sample extract of output of running `python code/fingeRNAt.py -r example_inputs/1aju_model1.pdb -l example_inputs/ligands.sdf -wrapper PuPy`

	<br/>

	<img src="docs/README_pics/pupy_xp.png" width="1000" />

	Sample extract of output of running `python code/fingeRNAt.py -r example_inputs/1aju_model1.pdb -l example_inputs/ligands.sdf -f XP -wrapper PuPy`

- `Counter`

	Counts total number of given interaction type for any SIFt type. In `SIMPLE`, `PBS` & `FULL` SIFs types sums all binary interactions values, but if ran together with `XP`, calculates total number of non-Pi interaction & total number of binary interactions values for Pi interactions (see -> [`XP`](#xp)).

	<img src="docs/README_pics/counter_full.png" width="900" />

	Sample extract of output of running `python code/fingeRNAt.py -r example_inputs/1aju_model1.pdb -l example_inputs/ligands.sdf -wrapper Counter`

	<br/>

	<img src="docs/README_pics/counter_xp.png" width="1000" />

	Sample extract of output of running `python code/fingeRNAt.py -r example_inputs/1aju_model1.pdb -l example_inputs/ligands.sdf -f XP -wrapper Counter`

## Visualization

All SIFs outputs can be visualized as heatmap and saved as png files with the same name as tsv output.

<img src="docs/README_pics/simple_heatmap.png" width="1200" />

Heatmap for SIFt type `SIMPLE` obtained from running `python code/fingeRNAt.py -r example_inputs/1aju_model1.pdb -l example_inputs/ligands.sdf -f SIMPLE -vis`

<br/>

<img src="docs/README_pics/full_acug_heatmap.png" width="1000" />

Heatmap for SIFt type `FULL` with wrapper `ACUG` obtained from running `python code/fingeRNAt.py -r example_inputs/1aju_model1.pdb -l example_inputs/ligands.sdf -vis -wrapper ACUG`

<br/>

<img src="docs/README_pics/xp_pupy_heatmap.png" width="800" />

Heatmap for SIFt type `XP` with wrapper `PuPy` obtained from running `python code/fingeRNAt.py -r example_inputs/1aju_model1.pdb -l example_inputs/ligands.sdf -f XP -vis -wrapper PuPy`

<br/>

<img src="docs/README_pics/xp_counter_heatmap.png" width="400" />

Heatmap for SIFt type `XP` with wrapper `Counter` obtained from running `python code/fingeRNAt.py -r example_inputs/1aju_model1.pdb -l example_inputs/ligands.sdf -f XP -vis -wrapper Counter`

## Usage examples

`python code/fingeRNAt.py -r example_inputs/1aju_model1.pdb -l example_inputs/ligands.sdf -f SIMPLE -o /path/to/my_output`

Calculates fingerprint SIMPLE and saves the SIFs output in the declared location.

`python code/fingeRNAt.py -r example_inputs/1aju_model1.pdb -l example_inputs/ligands.sdf -f PBS -vis`

Calculates fingerprint PBS and saves the SIFs output with the deafult name in the current directory together with results heatmap.

`python code/fingeRNAt.py -r example_inputs/1aju_model1.pdb -l example_inputs/ligands.sdf -wrapper ACUG,PuPy,Counter`

Calculates default fingerprint FULL and saves the SIFs output and 3 SIFs wrapped outputs with the deafult names in the current directory.

`python code/fingeRNAt.py -r example_inputs/1aju_model1.pdb -l example_inputs/ligands.sdf -f XP -dha -o /path/to/my_output -vis -wrapper ACUG`

Calculates fingerprint XP considering donor-hydrogen-acceptor angle calculation when detecting hydrogen bonds and saves the SIFs output, 1 SIFs wrapped output and 2 results heatmaps in the declared location.

## Documentation

To generate documentation file using sphinx:

`cd docs`

`make html`

The documentation will be available from `_build/html`.

## Unit test

To run a unit test:

`cd tests`

`python fingeRNAt_test.py`

## GUI

To use GUI, simply run

`python code/gui.py`

GUI is user-friendly and has all aforementioned functionalities.

<img src="docs/README_pics/gui.png" width="600" />


# Feedback

We welcome any feedback, please send an email to Natalia Szulc <img src="docs/README_pics/nszulc_mail.png" width="130" />


# Acknowledgments

Special thanks of gratitude to [Masoud Farsani](https://github.com/mafarsani), [Pritha Ghosh](https://github.com/prithaghosh) and [Tomasz Wirecki](https://github.com/fryzjergda) for their invaluable feedback as well as to Prof. Janusz M. Bujnicki and the entire [Bujnicki Lab](http://genesilico.pl) for all the support and project guidelines.

Extensive script testing provided by Zuzanna Mackiewicz has been a great help in developing this tool.

Assistance provided by [Open Babel Community](http://openbabel.org) was greatly appreciated.


# Citing

Authors:

Natalia A. Szulc,
<img src="docs/README_pics/nszulc_mail.png" width="130" />

Filip Stefaniak, &nbsp;&nbsp;
<img src="docs/README_pics/fstefaniak_mail.png" width="135" />

<br />

If you use this software, please cite:

**fingeRNAt - a software for analysis of nucleic acids-ligand complexes. Design and applications.**

Natalia A. Szulc, Zuzanna Mackiewicz, Janusz M. Bujnicki, Filip Stefaniak
[in preparation]

# License

fingeRNAt is licensed under the GNU General Public License v3.0
