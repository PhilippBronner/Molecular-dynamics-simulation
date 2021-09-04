# z0-precisionmeasurement

Precision measurements and tests of the Standard Model using OPAL data at LEP. \\
The project consists of several parts. :

- Cut selection in the Monte Carlo (MC) data, to separate the different particle types
- Test of the selected cuts on the cummulative MC data.
- t- and s- channel seperation for the electrons
- Calculation of the total cross section of for the Z-boson, $ e^+e^- -> f \bar{f} $, for resulting in electrons, muons, tauons or hadronic particles.
- Evaluations following from the crossestion, depending on the energy and fitted with a Breit-Wegner function.
- Tests on the forward-backward asymmetry and $\sin^2(\theta_\text{W})$ in muon final states
- Tests on the lepton universality.

## Installing Python

### with conda (on windows)

1. Install [anaconda](https://docs.anaconda.com/anaconda/install/) (miniconda probably works too?).
2. `conda env create` reads the `environment.yml` file in this
    repository, creates a new env and installs all necessary packages
    into it.
3. Activate the new env: `conda activate z0-env`

### on linux

1. Installing python

```
sudo apt-get install python3
```

2. Insall pip

```
sudo apt-get install python3-pip
```

3. Using pip to install jupyter

```
python3 -m pip install jupyter
```

## Installing python libraries

1. Install the necessary python libraries using pip

```
python3 -m pip install uproot
python3 -m pip install awkward
python3 -m pip install mplhep
python3 -m pip install numpy
python3 -m pip install matplotlib
python3 -m pip install scipy
```

## Clone the Repository

1. Clone the repository from github

```
git clone https://github.com/lukasg96/z0-precisionmeasurement.git
```

# Usage

1. Start `jupyter-lab` (or `jupyter-notebook` if you prefer) in the
    environment. This will usually open your browesr automatically.

2. navigate to the file `z0_experiment_Philipp_Bronner_Lukas_Grunwald.ipynb` in the repository and open in jupyter.

3. Run the code cells from top to bottom (in order)

# Contact

FPII group nr: 8

Philipp Bronner: [philippbronner-at-t-online.de](http://philippbronner-at-t-online.de)

Lukas Grunwald: [lukas-at-grunwald-elzach.de](http://lukas-at-grunwald-elzach.de)

Project Link: https://github.com/lukasg96/z0-precisionmeasurement.git

