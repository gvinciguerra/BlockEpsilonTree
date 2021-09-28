<p align="center"><img alt="Block-ε tree" src="http://pages.di.unipi.it/vinciguerra/img/block_epsilon_tree.svg" height="150"></p>

The [block-ε tree](http://pages.di.unipi.it/vinciguerra/publication/repetition-and-linearity-aware-rank-select-dictionaries/) is a compressed rank/select dictionary that achieves new space-time trade-offs by exploiting the approximate linearity and the repetitiveness of the data.
It is based on a combination of the _LA-vector_ ([paper](https://doi.org/10.1137/1.9781611976472.4), [code](https://github.com/gvinciguerra/la_vector)) and the _block tree_ ([paper](https://doi.org/10.1016/j.jcss.2020.11.002), [code](https://github.com/elarielcl/BlockTrees)).

## Usage

This is a header-only library. To compile the [example](example.cpp), use the following commands:

```sh
git clone https://github.com/gvinciguerra/BlockEpsilonTree.git
cd BlockEpsilonTree
cmake . -DCMAKE_BUILD_TYPE=Release
make -j8
```

## License

This project is released for academic purposes under the terms of the GNU General Public License v3.0. Some methods implemented in this project are **patent pending**.

If you use this code for your research, please cite:

> Paolo Ferragina, Giovanni Manzini, and Giorgio Vinciguerra. Repetition- and linearity-aware rank/select dictionaries.  In: Proceedings of the 32nd International Symposium on Algorithms and Computation (ISAAC), 2021.

```bibtex
@inproceedings{Ferragina:2021isaac,
  author = {Ferragina, Paolo and Manzini, Giovanni and Vinciguerra, Giorgio},
  booktitle = {Proceedings of the 32nd International Symposium on Algorithms and Computation (ISAAC)},
  title = {Repetition- and linearity-aware rank/select dictionaries},
  year = {2021}}
```