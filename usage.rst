Usage
=====

.. _installation:

Installation
------------

To use MIRISim Coronagraphy, first install MIRISim: https://wiki.miricle.org/Public/MirisimInstallation.

You should now be able to ``conda activate mirisim``.

Next:

1. git clone the following MIRISim Coronagraphy repository into wherever you want it to live on your computer: https://github.com/spacegal-spiff/mirisim_coronagraphy

2. cd to the location of the original MIRISim installation on your computer

3. Rename this new mirisim folder (e.g. ``mv mirisim mirisim_old``) so you can create a new symbolic link to the new Coronagraphy repository

4. Create a symbolic link to the Coronagraphy repository (e.g. ``ln -s /location/of/mirisim_coronagraphy mirisim``)

5. Within the mirisim environment:

    a. Install webbpsf and poppy using these instructions: https://webbpsf.readthedocs.io/en/latest/installation.html#installing-the-required-data-files
    b. Install webbpsf_ext using these instructions: https://pypi.org/project/webbpsf-ext/
    c. Set the required paths outlined in the instructions above (you can also add these paths to your .bashrc/.zshrc/equivalent files so you don't have to set the paths before every session)


.. _example_files:

Example files
------------

To run MIRISim Coronagraphy after installing it, you'll need need to cd into a folder containing the following three files:

1. A configuration .yaml file
2. A setup .ini file
2. The actual python script to run the simulation

An example configuration file:
.. literalinclude:: example_setup_config.yaml
  :language: YAML 

An example setup file:
.. literalinclude:: scene.ini

A script to run the simulation:
.. literalinclude:: run_mirisim_coron.py