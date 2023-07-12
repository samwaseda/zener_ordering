# Order-disorder phase diagram of Zener ordering

This repository contains a workflow that allows the user to study the order-disorder transition phase diagram using metadynamics for varying temperature and carbon concentration.

## Theory


We define the order parameter $z$ via:

$$z = \sqrt s = \sqrt{\frac{3}{2}\frac{\sum_in_i^2}{n_{\mathrm C}^2}-\frac{1}{2}}$$ 
$$n_i = \sum_a\frac{\sum_p \cos^2\left(\sum_kE_{pik}x_{ak}\right)}{\sum_{pj} \cos^2\left(\sum_kE_{pjk}x_{ak}\right)}$$
where we defined the matrices $E$ via:

```math
E_{kij} = \begin{pmatrix} 0 & 1 & (-1)^k \\ (-1)^k & 0 & 1 \\ 1 & (-1)^k & 0 \\\end{pmatrix}
```

The main point is $E$, where it used to have only the first matrix, but this one was not enough to have a symmetric distribution in non-octahedral sites. For the energy calculation, the gradient is given by:

$$\frac{\partial z}{\partial x_{al}}= \frac{3}{2zn_{\mathrm C}^2}\sum_in_i\left(-\frac{\sum_p E_{pil}\sin\left(2\sum_kE_{pik}x_{ak}\right)}{\sum_{pj} \cos^2\left(\sum_kE_{pjk}x_{ak}\right)}+\frac{\sum_{pj} E_{pjl}\sin\left(2\sum_kE_{pik}x_{ak}\right)\sum_p\cos^2\left(\sum_kE_{pik}x_{ak}\right)}{\left(\sum_{pj} \cos^2\left(\sum_kE_{pjk}x_{ak}\right)\right)^2}\right)$$

\end{document}

