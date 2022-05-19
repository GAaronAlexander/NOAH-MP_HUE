# hrldas

This is the offical HRLDAS/Noah-MP unified Github repository for code downloading and contribution. Note that the HRLDAS/Noah-MP model is a community model contributed by the whole Noah-MP scientific community. For maintenance and release of this GitHub, please contact: Cenlin He (cenlinhe@ucar.edu) and Fei Chen (feichen@ucar.edu).

HRLDAS (High Resolution Land Data Assimilation System) contains the Noah-MP Land Surface Model and can
be executed in offline uncoupled mode for a 2-D gridded domain or a single point. 

As of May 6, 2020, this repository does not explicitly contain the Noah-MP code and parameter tables, but
connects to the Noah-MP repository (https://github.com/NCAR/noahmp) through git submodules. To clone
this repository and populate with the Noah-MP code/parameters, use the following command:

git clone --recurse-submodules https://github.com/NCAR/hrldas


Some changes have been made to the structure of archiving the stand-alone version of HRLDAS/Noah-MP code in the Github repository. Now, it separately archives the core Noah-MP source code in another core Noah-MP Github repository (https://github.com/NCAR/noahmp) and the rest of the HRLDAS-related files in this HRLDAS Github repository. The HRLDAS Github repo is already linked to the core Noah-MP code repository through git submodules. This new archiving structure will allow different host model platforms/systems (e.g., HRLDAS, WRF, UFS, NWM, LIS, etc.) to connect to the core Noah-MP source code and develop their own host model drivers. 


Model developers can make code development based on the develop branch and create pull request to the develop branch. The pull request will be reviewed by Noah-MP model release team and if everything looks good, the new code development will be merged to the official 'develop' branch. Eventually, the updates in develop branch will be merged to the master branch for official Noah-MP model release.

Branch structures of this Noah-MP repository:

1. "master" branch: most stable & recent version (we suggest users to use this branch), updated whenever there are bug fixes or major model update/release, automatically connected to the master branch of the Noah-MP core repository through git submodule;

2. "develop" branch: used for continuous HRLDAS/NoahMP development, keep being updated by including bug fixes and code updates (e.g., new physics options, processes, etc.), automatically connected to the develop branch of the Noah-MP core repository through git submodule; 

3. "release-v4.3-WRF" branch: stable version consistent with the WRF version 4.3 release;

4. "wrfv4.2" branch: stable version consistent with the WRF version 4.2 release;

*Note: The HRLDAS branches and codes will always be consistent with WRF and currently include the gecros crop scheme.

Some suggestions for model developers to contribute to Noah-MP code through the Github repository (typical procedures):

1. Step (1) Create a fork of this official HRLDAS repository to your own Github account; 

2. Step (2) Make code updates based on the "develop" branch of the forked repository under your own account; 

3. Step (3) Finalize and test the code updates you make; 

4. Step (4) Submit a pull request for your model updates from your own forked Github repository to the "develop" branch of the official repository;

5. Step (5) The Noah-MP release team reviews and tests the model updates in the submitted pull request and discusses with the developer if there is any problem; 

6. Step (6) The Noah-MP release team confirms the pull request and merges the updated code to the "develop" branches in the official repository;

7. Step (7) The Noah-MP release team will merge the updated "develop" branch to the master branch and version-release branch during the annual model release.

Note: If updates are made to both the core NoahMP source code and other HRLDAS-related files, two separate pull requests need to be submitted to the core NoahMP repository (https://github.com/NCAR/noahmp) and this HRLDAS repository, respectively, regarding the changes in the Noah-MP code file and other HRLDAS files. This could be done by using the same title for the two pull requests to ensure that the submitted code changes are handled together by the release team in the two repositories.
