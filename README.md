# KaneMele

This is a simple MATLAB implementation of the Kane-Mele model[^1].

## How to Use

First, run `graphene_kane_mele.m` in MATLAB. This will generate a file named `wan_basis.mat`,  
which contains all the required information for the hopping model.  
Then, run `wannier_solver.m` to compute and plot the full band structure.

If needed, you can modify the parameters at the beginning of `graphene_kane_mele.m`:

```Matlab
hopping_on_site_A = 0;
hopping_on_site_B = 0;
hopping_nn_t = -3;

hopping_nnn_kane_mele = -0.05;
hopping_nnn_t = 0;
```

## File List

| File                 | Description                                                                 |
|----------------------|-----------------------------------------------------------------------------|
| `graphene_kane_mele.m` | Generates all hopping constants for the Kane-Mele model                    |
| `wannier_solver.m`     | Solves a generalized Wannier hopping model and plots the band structure    |
| `KPOINTS`              | Specifies the band structure path, using the [VASP format](https://www.vasp.at/wiki/index.php/KPOINTS) |

## Example Output

<img src="https://github.com/angushphys/KaneMele/blob/main/result/kane_mele_result.png" width="560" height="600">

---

[^1]: C. L. Kane and E. J. Mele, [Phys. Rev. Lett. 95, 226801 (2005)](https://doi.org/10.1103/PhysRevLett.95.226801)
