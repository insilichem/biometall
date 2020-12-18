.. BioMetAll: Identifying metal-binding sites in proteins from backbone preorganization

   https://github.com/insilichem/biometall

   Copyright 2020 José-Emilio Sánchez-Aparicio, Laura Tiessler-Sala,
   Lorea Velasco-Carneros, Lorena Roldán-Martín, Giuseppe Sciortino,
   Jean-Didier Maréchal


BioMetAll
=========

.. image:: https://readthedocs.org/projects/biometall/badge/?version=latest
   :target: https://biometall.readthedocs.io/en/latest/

.. image:: https://img.shields.io/pypi/v/biometall
   :target: https://pypi.org/project/biometall/

.. image:: https://img.shields.io/badge/python-3.7-green
   :target: https://www.python.org/downloads/release/python-377/

.. image:: https://img.shields.io/pypi/l/biometall
   :target: https://opensource.org/licenses/BSD-3-Clause

.. image:: https://img.shields.io/badge/doi-https%3A%2F%2Fdoi.org%2F10.1021%2Facs.jcim.0c00827-blue
   :target: https://doi.org/10.1021/acs.jcim.0c00827

BioMetAll is a command line application to allow the identification of metal-binding
sites in proteins from backbone preorganization.

.. image:: docs/images/logo_biometall.png
    :alt: BioMetAll logo

Features
--------

**Different options to customize the search:**

- sites with a minimum of coordinating amino acids
- sites with a custom combination of coordinating amino acids
- taking into account coordinations with backbone oxygen atoms (or only with the side chains)
- scanning the whole protein or only a region of interest

**Search for a structural motif and propose mutations:**

- scan the structure to search a specific motif (e.g. HIS,HIS,ASP/GLU)
- propose mutations to complete the motif (e.g. having a HIS,HIS motif, search for mutations to complete with either a GLU or an ASP)

**Possible applications**

- screening of a pool of .pdb structures
- identification of conformational changes that alter the formation of metal-binding sites
- metalloenzyme design

A more complete overview of the capabilities of the program, illustrated with several case studies, is
available in our `ChemRxiv preprint <https://doi.org/10.26434/chemrxiv.12668651.v1>`_.

Documentation and support
-------------------------

Documentation is available at `this link <https://biometall.readthedocs.io/en/latest/>`_.

Installation instructions are available at `this webpage <https://biometall.readthedocs.io/en/latest/installation.html>`_.

If you need help with BioMetAll, please use the `issues page <https://github.com/insilichem/biometall/issues>`_ of this GitHub repository. You can drop me a message at `joseemilio.sanchez@uab.cat <mailto:joseemilio.sanchez@uab.cat>`_ too.

License
-------

BioMetAll is an open-source software licensed under the BSD-3 Clause License. Check the details in the `LICENSE <https://github.com/insilichem/biometall/blob/master/LICENSE>`_ file.

History of versions
-------------------

- **v1.0:** Release version used in the preparation of the JCIM article.

- **v0.2:** New *backbone_clashes* and *sidechain_clashes* parameters, which allow to customize the filtering of probes with clashes.

- **v0.1:** Release version used in the preparation of the ChemRxiv preprint.

OS Compatibility
----------------

BioMetAll is compatible with Linux, macOS and Windows.

If you find some dificulties when installing it in a concrete distribution, please use the `issues page <https://github.com/insilichem/biometall/issues>`_ to report them.

How to cite this software
-------------------------

Sánchez-Aparicio, J.-E.; Tiessler-Sala, L.; Velasco-Carneros, L.; Roldán-Martín, L.; Sciortino, G.; Maréchal, J.-D.. Biometall: Identifying Metal-binding Sites in Proteins from Backbone Preorganization, *J. Chem. Inf. Model.*, **2020**, https://dx.doi.org/10.1021/acs.jcim.0c00827.

Acknowledgements
----------------

Project template based on the
`Computational Molecular Science Python Cookiecutter <https://github.com/molssi/cookiecutter-cms>`_ version 1.2.

Standalone executables have been created with `PyInstaller <https://www.pyinstaller.org/>`_














