# Application for calculate partial atomic charges with empirical method Molecular Graph Charge (MGC) and Orbital Graph Charge (OGC).
by Šimon Zouvala  "445475 (at) mail (dot) muni (dot) cz"
github link: https://github.com/SimonZouvala/MGC-OGC

*Requirements*
- Matplotlib, NumPy libraries

*Compilation*
- Simply run on terminal write python3.6 main.py + argument for run part of application what you want 

*Arguments for run*
---------------------------------------------------------------
usage: main.py calculation [-h] [--eem] [--mgc] [--ogc]
                           [--parameters PARAMETERS] [--output OUTPUT]
                           [calculate]

positional arguments:
  calculate             Give one file with molecules (.sdf) for mgcm and one
                        file with parameters for EEM (.xml)

optional arguments:
  -h, --help            show this help message and exit
  --eem                 Give this argument, if you want calculate with EEM
  --mgc                 Give this argument, if you want calculate with MGC
  --ogc                 Give this argument, if you want calculate with OGC
  --parameters PARAMETERS
                        Give this argument for parameters, if you want
                        calculate with EEM
  --output OUTPUT       Give a name file, for output calculate

---------------------------------------------------------------
usage: main.py structure [-h] [--parameters PARAMETERS]
                         [--molecules MOLECULES] [--nobond]

optional arguments:
  -h, --help            show this help message and exit
  --parameters PARAMETERS
                        Give a file with parameters (EEM) (.xml)
  --molecules MOLECULES
                        Give a file with molecules (.sdf)
  --nobond              Give this argument, if you want not type bond.

---------------------------------------------------------------
usage: main.py graph [-h] draw_graph draw_graph


positional arguments:
  draw_graph  Give two charges results file, if you have not it, you can used
              CALCULATION

optional arguments:
  -h, --help  show this help message and exit

