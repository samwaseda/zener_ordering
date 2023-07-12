# Order-disorder phase diagram of Zener ordering

This repository contains a workflow that allows the user to study the order-disorder transition phase diagram using metadynamics for varying temperature and carbon concentration.

## Introduction

Zener ordering is a phenomenon in which carbon atoms in ferrite (bcc Fe) occupy one of the three sub-lattices of octahedral interstitial sites. It has been first theorized in the middle of the last century, but it has so far not been observed on the microscopic level.

## Theory

### Order parameter $z$

Nomenclature:

- $i, j, k, l$: Coordinate indices
- $\alpha$: Enumeration of atoms
- $p$: Index for 0 and 1

We define the order parameter $z$ via:

```math
z = \sqrt{\frac{3}{2}\frac{\sum_in_i^2}{n_{\mathrm C}^2}-\frac{1}{2}}
```

where $n_{\mathrm C}$ is the total number of carbon atoms and $n_i$ is the number of carbon atoms in the sub-lattice $i$, which is defined by:

```math
n_i = \sum_\alpha\frac{\sum_p \cos^2\left(\sum_kE_{pik}x_{\alpha k}\right)}{\sum_{pj} \cos^2\left(\sum_kE_{pjk}x_{\alpha k}\right)}
```

where we defined the matrices $E$ via:

```math
E_{pij} = \begin{pmatrix} 0 & 1 & (-1)^p \\ (-1)^p & 0 & 1 \\ 1 & (-1)^p & 0 \\\end{pmatrix}
```

The main point is $E$, where it used to have only the first matrix, but this one was not enough to have a symmetric distribution in non-octahedral sites. For the energy calculation, the gradient is given by:

```math
\frac{\partial z}{\partial x_{\alpha l}}= \frac{3}{2zn_{\mathrm C}^2}\sum_in_i\left(-\frac{\sum_p E_{pil}\sin\left(2\sum_kE_{pik}x_{\alpha k}\right)}{\sum_{pj} \cos^2\left(\sum_kE_{pjk}x_{\alpha k}\right)}+\frac{\sum_{pj} E_{pjl}\sin\left(2\sum_kE_{pik}x_{\alpha k}\right)\sum_p\cos^2\left(\sum_kE_{pik}x_{\alpha k}\right)}{\left(\sum_{pj} \cos^2\left(\sum_kE_{pjk}x_{\alpha k}\right)\right)^2}\right)
```

This definition of order parameter therefore allows for the calculation of derivatives for each atom.

### Metadynamics

We define the histogram of metadynamics via:

```math
B(z) = w\sum_i e^{-\frac{-(z-z_i)^2}{2\sigma^2}}
```
