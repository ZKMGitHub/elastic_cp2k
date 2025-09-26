# elastic_cp2k
Automating the Calculation of Crystal Elastic Constants Using CP2K
## Install from source code
```bash
git clone https://github.com/ZKMGitHub/elastic_cp2k.git
cd elastic_cp2k
pip install .
```
## Preparation
* The executable, compiled, parallel version of CP2K. (pspm or popt)
* Input file for performing cell optimization tasks in CP2K. (Should be named `cellopt.inp`)
* Input file for performing geometry optimization tasks in CP2K. (Should be named `geoopt.inp`)
* Initial structure file of the calculation system. (Format: .xyz)
* Basis set and pseudopotential files. (e.g. `BASIS_MOLOPT` and `POTENTIAL`)
* `config.yaml`. (Translatable parameters that can be modified for program invocation)
```bash
# configs/config.yaml.example

num_cores: 52
cp2k_exec_path: "/home/zkm/CP2K/cp2k-2024.3/exe/local/cp2k.psmp"
cellopt_inp: "cellopt.inp"
cellopt_out: "cellopt.out"
cellopt_jobname: "cellopt"
geoopt_inp: "geoopt.inp"
geoopt_out: "geoopt.out"
geoopt_jobname: "geoopt"
potential_file: "POTENTIAL"
basis_set_file: "BASIS_MOLOPT"
dftd3_dat: "dftd3.dat"
up: 0.01
```
## Example
For the supercell calculation of elastic constants for Jennite with a $1\times2\times1$ supercell.
- Activate CP2k Environment. (Need to set up the environment yourself)
  ```bash
  source ~/.bashrc
  conda activate CP2K
  ```
- Configure environment variables. (Replace with your own executable path)
  ```bash
  source /home/zkm/CP2K/cp2k-2024.3/tools/toolchain/install/setup
  ```
- Define the parameters in `config.yaml` file.
  ```bash
  cd ./example
  ```
- Run the calculation process.
  ```bash
   elastic_cp2k or elastic_cp2k -c path/to/your/config.yaml
  ```
After the run is completed, a ```Cij.dat``` file will be generated in the current folder to record the values of the elastic constants.
- e.g. Jennite $1\times2\times1$
```bash
113.931 41.606 37.737 4.471 1.955 5.386
41.606 124.829 42.071 6.368 -9.861 -8.352
37.737 42.071 65.579 -2.692 -0.739 -2.444
4.471 6.368 -2.692 25.023 1.085 -7.662
1.955 -9.861 -0.739 1.085 16.353 3.684
5.386 -8.352 -2.444 -7.662 3.684 37.282
```
## Description
- The method for constructing the crystal cell is:
  - cell parameter *a* is parallel to *X*-axis and *b* is in the *XY* plane
- $\epsilon_4$, $\epsilon_5$, and $\epsilon_6$ are the Engineering shear strains between *YZ*, *XZ*, and *XY* planes.

$$
  \begin{bmatrix}
   \sigma_1 \rule{0pt}{0.75em}\\
   \sigma_2 \rule{0pt}{0.75em}\\
   \sigma_3 \rule{0pt}{0.75em}\\
   \sigma_4 \rule{0pt}{0.75em}\\
   \sigma_5 \rule{0pt}{0.75em}\\
   \sigma_6 \rule{0pt}{0.75em}
  \end{bmatrix}=
  \begin{bmatrix}
   C_{11} & C_{12} & C_{13} & C_{14} & C_{15} & C_{16} \rule{0pt}{0.75em}\\
   & C_{22} & C_{23} & C_{24} & C_{25} & C_{26} \rule{0pt}{0.75em}\\
   & & C_{33} & C_{34} & C_{35} & C_{36} \rule{0pt}{0.75em}\\
   & & & C_{44} & C_{45} & C_{46} \rule{0pt}{0.75em}\\
   & & & & C_{55} & C_{56} \rule{0pt}{0.75em}\\
   & & & & & C_{66} \rule{0pt}{0.75em}
  \end{bmatrix}
  \begin{bmatrix}
   \epsilon_1 \rule{0pt}{0.75em}\\
   \epsilon_2 \rule{0pt}{0.75em}\\
   \epsilon_3 \rule{0pt}{0.75em}\\
   \epsilon_4 \rule{0pt}{0.75em}\\
   \epsilon_5 \rule{0pt}{0.75em}\\
   \epsilon_6 \rule{0pt}{0.75em}
  \end{bmatrix}
$$
