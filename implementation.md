fingeRNAt implementation details
=========


<!-- TOC START min:1 max:6 link:true asterisk:false update:true -->
- [Detection of interactions](#detection-of-interactions)
  - [Nucleic acid properties](#nucleic-acid-properties)
  - [Ligand properties](#ligand-properties)
  - [Molecular Interactions' Geometric Rules](#molecular-interactions-geometric-rules)
  - [1. Hydrogen bonds](#1-hydrogen-bonds)
  - [2. Halogen bonds](#2-halogen-bonds)
  - [3. Cation - anion interactions](#3-cation---anion-interactions)
  - [4. Pi - cation & 5. Pi - anion interactions](#4-pi---cation--5-pi---anion-interactions)
  - [6. Pi - stacking interactions](#6-pi---stacking-interactions)
  - [7. Ion-mediated interactions](#7-ion-mediated-interactions)
  - [8. Water-mediated interactions](#8-water-mediated-interactions)
  - [9. Lipophilic interactions](#9-lipophilic-interactions)
<!-- TOC END -->


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

<img src="docs/README_pics/interactions-explanation.png" width="900"/>


## 1. Hydrogen bonds

  Geometric rule:
  - |D - A| < 3.9 &#8491;
	(*Torshin, Weber, & Harrison*, 2002)


  > **_NOTE:_** If hydrogens are present in RNA/DNA structure, fingeRNAt can be run with parameter `-dha`, that additionaly calculates Donor-Hydrogen-Acceptor angle  used as supplementary criteria in hydrogen bonds detection:<br/>100&deg; < D-H-A angle < 260&deg; (*Adasme et al.*, 2021)<br/>
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

## 4. Pi - cation & 5. Pi - anion interactions

  Geometric rules:

  - |cation/anion - aromatic ring center| < 6.0 &#8491;
  (*Gallivan and Dougherty*, 1999)
  - angle between the ring plane and the line between cation/anion - ring center ~ 90&deg; &#177; 30&deg;

 > **_NOTE:_** In case of Pi - cation interaction, aromatic ring is only from the nucleic acid side (as nucleic acids are negatively charged). However, Pi - anion interaction is considered both ways: (i) nucleic acid's aromatic ring - ligands's anion and (ii) nucleic acid's anion (from phosphate group; see above) - ligand's aromatic ring.

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
  * |nucleic acid's nitrogen/oxygen - ion| <= X, where:
    * X = 3.2 &#8491; for magnesium ion (*Zheng et al.*, 2015)
    * X = 3.9 &#8491; for potassium ion (*Zheng et al.*, 2008)
    * X = 3.6 &#8491; for sodium ion (*Zheng et al.*, 2008)
    * X = 3.5 &#8491; for other ions (*Zheng et al.*, 2008)

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
