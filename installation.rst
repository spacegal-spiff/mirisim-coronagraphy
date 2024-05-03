Installation
=====


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

