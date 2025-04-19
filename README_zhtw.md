KaneMele
====================
這是一個簡單實作 Kane-Mele model [^1] 的 Matlab 程式。

如何使用
--------------------
先用 Matlab 執行 `graphene_kane_mele.m`，這會生成一個 `wan_basis.mat`的檔案，
裡面會包含 hopping model 所有需要的訊息，再使用 `wannier_solver.m`，
就可以畫出完整的 band。

需要的話則可以修改 `graphene_kane_mele.m` 開頭的參數

```Matlab
hopping_on_site_A = 0;
hopping_on_site_B = 0;
hopping_nn_t = -3;

hopping_nnn_kane_mele = -0.05;
hopping_nnn_t = 0;
```


檔案列表
--------------------

| File                 | Descript  |
|----------------------|-----------|
| graphene_kane_mele.m | 生成 Kane-Mele model 的所有 hopping costatant                                             |
| wannier_solver.m     | 用來解一般化的 Wannier hopping model，並畫出band structure                                  |
| KPOINTS              | 描述 band structure 路徑，使用[VASP的格式](https://www.vasp.at/wiki/index.php/KPOINTS)      |

參考結果
--------------------
<img src="https://github.com/angushphys/KaneMele/blob/main/result/kane_mele_result.png" width="560" height="600">


[^1]: C. L. Kane and E. J. Mele, [Phys. Rev. Lett. 95, 226801 (2005)](https://doi.org/10.1103/PhysRevLett.95.226801)
