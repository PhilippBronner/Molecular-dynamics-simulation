# Water droplett

Molecular dynamics simulation of a water droplet. \\
The evolution of the water molecules from an initial grid constellation to an equilibrium constellation can be observed. 
Furthermore the motion of the molecules in the stable droplet can be observed and material parameters and functions can be extracted from this data.

The project consists of several parts. :

- Initialization of the melecular setup on a 3D grid
- Equilibarion of the initial particles until the equilibrium constellation of a water droplet is reached.
- The stable droplet envolves over a period of time and the trajectory is stored
- Analysis: The collected data can be used to determine material parameters and functions as the radial-distribution function.
- Visualization: The created "trajectroy.dump" files can be visulaized with programms as "OVITO".

## Installing Python

### with conda (on windows)

1. Install [anaconda](https://docs.anaconda.com/anaconda/install/) (miniconda probably works too?).
2. `conda env create` reads the `environment.yml` file in this
    repository, creates a new env and installs all necessary packages
    into it.
3. Activate the new env: `conda activate MD-sim-env`

### on linux

1. Installing python

```
sudo apt-get install python3
```

2. Insall pip

```
sudo apt-get install python3-pip
```


3. Using pip to install spyder

```
python3 -m pip install spyder
```

## Installing python libraries

1. Install the necessary python libraries using pip

```
python3 -m pip install tqdm
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

1. Start `spyder` and execute the file "Main_Water_Droplet_Philipp_Bronner.py". This can also be done in every other python compiler.

2. The created "trejectory.dump" files can be opend with OVITO and the trajectory can be visulaized.


# Contact

Philipp Bronner: [philippbronner-at-t-online.de] (http://philippbronner-at-t-online.de)


Project Link: https://github.com/PhilippBronner/Molecular-dynamics-simulation.git
