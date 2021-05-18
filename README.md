<b>Description</b>
<p>Python-based system builder for constant-pH simulations in Gromacs.</p>

<b>Install instructions</b>
1. If you have a GPU and want to use GPU-acceleration, make sure you first install <a href="https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html#pre-installation-actions">CUDA</a>.
2. Clone <a href="https://bitbucket.org/berkhess/gromacs-constantph/branch/clean-cpHMD-branch">clean-cpHMD-branch</a>.
3. Install using the instructions <a href="https://manual.gromacs.org/documentation/current/install-guide/index.html">here</a>. I personally use the following CMake flags:
`cmake .. -DGMX_BUILD_OWN_FFTW=ON -DGMX_GPU=CUDA -DGMX_CUDA_TARGET_COMPUTE=60 -DGMX_USE_RDTSCP=ON -DGMX_SIMD=AVX2_256 -DCMAKE_INSTALL_PREFIX=/usr/local/gromacs_constantph`
4. Clone phbuilder (this) repository.
5. Add phbuilder directory to PYTHONPATH by adding `export PYTHONPATH=$PYTHONPATH:/path/to/phbuilder` to your `~/.bashrc` file (and reload terminal(s)).

<b>Design principles</b>
<p>I hardly use classes for the developing this builder as I don't particularly like them in Python. Instead, I decided to use the shelve module, which allows for storing and loading python data structures in a binary file. Here, this binary file is called 'universe'. The idea is to just store any data structure that might be of value to other builder/analysis functions in this 'universe', whether it be a float or a long list of Residue objects.</p>
<p>To reduce coupling, one thing I think is very useful is the ability to run any command at any stage of the building process. If a function requires a variable that was not added to the universe yet (because some other builder function was skipped), the function will simply prompt the user to input the missing data instead of exiting (any data structure except for d_residues can be manually added this way). </p>
<p>All (encapsulated) gromacs commands executed by the various builder functions will be logged in the "builder.log" file (this is helpful for debugging purposes). </p>

| Reference list of builder functions |
| :------------- |
| `universe.add(varName, value)` <br /> Add variable to universe. |
| `universe.has(varName)` <br /> Check whether variable with varName is in universe. |
| `universe.get(varName)` <br /> Retrieve variable with varName from universe. |
| `universe.inspect()` <br /> Print all the data members stored in universe. This is very useful for debugging or checking the parameters used for building the simulation in question. |
| `protein.process(fname, d_model=1, d_ALI='A', d_chain=[], resetResId=False):` <br /> Processes a .pdb file. This should be the first step in building the simulation. `fname` is the name of the .pdb file, `d_model` is the MODEL number, `d_ALI` is the alternative location indicator, `d_chain` is a list of chains (default = all), and `resetResId` resets the residue numbering. Note that all atoms in the input protein MUST belong to a chain (must have a chain letter), otherwise the building process won't work. |
| `protein.countres(resName)` <br /> Counts and returns the number of residues of a certain `resName`. |
| `protein.add_box(d_boxMargin, d_boxType='cubic')` <br /> Add periodic box. This function encapsulates `gmx editconf`. |
| `protein.add_buffer(ph_bufpdbName, ph_bufitpName, ph_bufMargin=2.0, ph_bufnmol=-1, attempts=10000)` <br /> Add buffer particle(s). This function encapsulates `gmx insert-molecules`. Will be skipped if either `ph_constantpH` or `ph_restrainpH` is False. Requires the structure (.pdb) file of the buffer particle as well as the topology (.itp) file. One can optionally set the minimum distance between the buffer particle(s) (themselves) and the protein, as well as the number of attempts. The number of buffer particles can also be set. By default, `ph_bufnmol` = #titratable groups. |
| `protein.add_water()` <br /> Solvate the system. This function encapsulates `gmx solvate`. |
| `add_ions(neutral=True, conc=0, pname='NA', nname='CL')` <br /> Neutralize the system and/or add a specified salt concentration (`conc` in mmol/ml) using `pname` and `nname` ion types (default NA and CL). This function encapsulates `gmx genion`. |
| `topol.generate(d_modelFF, d_modelWater, d_terministring="")` <br /> Generate the protein topology. This function encapsulates `gmx pdb2gmx`, and makes sure all protonatable residues are put in their protonated state when `ph_constantpH` is True. `d_modelFF` is the force field, while `d_modelWater` is the water model. `d_terministring` is a string consisting of two letters, specifying how `pdb2gmx` should treat the termini of the (various chains of) the protein. Commonly used: 00 = charged termini (`pdb2gmx` default), 11 = neutral termini, 34 = use this with capped tripeptide. |
| `topol.restrain_dihedrals(resName, atomNameList, Type, phi, dphi, fc)` <br /> Restrain dihedrals of a certain residue. This is an optional step that should be run after `topol.generate()`. `resName` is the residue name, `atomNameList` is a list containing four strings, indicating with four atoms form the dihedral (e.g. `[' OE1', ' CD ', ' OE2', ' HE2']` to restrain the carboxyl dihedral of GLU), and `Type`, `phi`, `dphi`, and `fc` are restraining parameters (see Gromacs manual). |
| `topol.restrain_dihedrals_by_idx(indices, Type, phi, dphi, fc)` <br /> Restrain dihedral between atom indices. This is an optional step that should be run after `topol.generate()`. `Type`, `phi`, `dphi`, and `fc` are restraining parameters (see Gromacs manual). |
| `md.gen_mdp(Type, nsteps=25000, nstxout=0, posres=False)` <br /> Generate .mdp file. Possible Types are: 'EM', 'NVT', 'NPT', 'MD'. If `posres=True`, the entire protein will be position restrained (not just COM). If you use energy_minimize(), tcouple(), and pcoule(), this function only needs to be called explicitly to generate MD.mdp.
| `md.gen_constantpH(ph_pH, ph_lambdaM, ph_nstout, ph_barrierE, cal=False, lambdaInit=0.5)` <br /> Generate constant-pH input parameters. Will be skipped if `d_constantpH` is False. This function currently works in four steps: 1) write general constant-pH parameters to MD.mdp, 2) write different lambda group types to MD.mdp, 3) write different lambda residues to MD.mdp, 4) write lambda groups and corresponding atom numbers to index.ndx.  `ph_pH` is the pH of the simulation/buffer, `ph_lambdaM` is the mass of the lambda particle, `ph_nstout` is the output frequency for the lambda_xxx.dat files, `ph_barrierE` is the barrier energy in kJ/mol. `cal` and `lambdaInit` should be left alone unless you're doing a calibration. |
| `md.energy_minimize()` <br /> Generate EM.mdp and perform energy minimization. Encapsulates `gmx grompp` and `gmx mdrun`. |
| `md.energy_tcouple()` <br /> Generate NVT.mdp and perform temperature coupling. Encapsulates `gmx grompp` and `gmx mdrun`. Velocities are generated in this step, and it should not be skipped. |
| `md.energy_pcouple()` <br /> Generate NPT.mdp and perform pressure coupling. Encapsulates `gmx grompp` and `gmx mdrun`. Can sometimes be skipped for small systems to save time. |
| `write.run(gmxPath="/usr/local/gromacs", options="")` <br /> Write executable bash file called run.sh which calls `gmx grompp` and `gmx mdrun` for running the simulation. `options` can be used to specify a string for additional parameters for mdrun such as `-pme cpu`, which is currently required when using GPU acceleration. |
| `write.reset()` <br /> Write executable bash file called reset.sh which can be used to clean up / delete all generate files. Proceed with caution! |

<br />

| Reference list of analysis functions  |
| :------------- |
| `analyze.fitCalibration(order=5, compare=[])` <br /> Fit dV/dl coefficients and plot the fit. `compare` can be used to specify different dV/dl coefficients for a comparison. |
| `analyze.glicphstates()` <br /> Analyze pH states in GLIC. |
| `analyze.plotlambda(plotBUF=False)` <br /> Plot lambda trajectories of all lambda coordinates. Also plot buffer trajectory if `plotBUF=True`.|
| `analyze.plotpotentials(pKa)` <br /> Plot the correction and pH potentials for given pKa. Uses lambda_dwp.dat for obtaining the bias potential. |
| `analyze.plotforces(pKa)` <br /> Plot the correction and pH forces for given pKa. Uses lambda_dwp.dat for obtaining the bias potential. |
| `analyze.plothistogram(fname, bins=200)` <br /> Plot histogram for a given lambda_xxx.dat file. |

<br />

<b>Reference list of data members used/stored in universe:</b>
* All data members have either a `d_` or a `ph_` prefix to distinguish them.
* `d_` means it's a general data member while `ph_` means it's specific to some constant-pH functionality.

<table class="tg">
<tbody>
  <tr>
    <td class="tg-0pky">d_ALI</td>
    <td class="tg-0pky">Alternative Location Indicator in .pdb files. Default = 'A'.</td>
  </tr>
  <tr>
    <td class="tg-0pky">d_box</td>
    <td class="tg-0pky">The (CRYST1) line a .pdb file defining the periodic box.</td>
  </tr>
  <tr>
    <td class="tg-0pky">d_boxMargin</td>
    <td class="tg-0pky">Box size margin / option -d for gmx editconf.</td>
  </tr>
  <tr>
    <td class="tg-0pky">d_chain</td>
    <td class="tg-0pky">List of different chains in protein in .pdb file.</td>
  </tr>
  <tr>
    <td class="tg-0pky">d_dt</td>
    <td class="tg-0pky">MD time step size (ps). Currently hardcoded at 0.002 ps.</td>
  </tr>
  <tr>
    <td class="tg-0pky">d_model</td>
    <td class="tg-0pky">The MODEL number that was specified when processing .pdb file. Default = 1.</td>
  </tr>
  <tr>
    <td class="tg-0pky">d_modelFF</td>
    <td class="tg-0pky">Force field.</td>
  </tr>
  <tr>
    <td class="tg-0pky">d_modelWater</td>
    <td class="tg-0pky">Water model.</td>
  </tr>
  <tr>
    <td class="tg-0pky">d_nameList</td>
    <td class="tg-0pky">List containing succession of .pdb names when performing various steps. Usually a tool will take the last name contained in this list and us that as input.</td>
  </tr>
  <tr>
    <td class="tg-0pky">d_nsteps</td>
    <td class="tg-0pky">Number of MD steps to perform during production (MD) run.</td>
  </tr>
  <tr>
    <td class="tg-0pky">d_nstxout</td>
    <td class="tg-0pky">Interval in which to output a trajectory frame. E.g. d_nstxout=10000 means output a frame every d_dt * d_nstxout = 20ps.</td>
  </tr>
  <tr>
    <td class="tg-0pky">d_pdbName</td>
    <td class="tg-0pky">Name of the input protein minus .pdb extention.</td>
  </tr>
  <tr>
    <td class="tg-0pky">d_residues</td>
    <td class="tg-0pky">List of (all) Residue objects. Basically it's the entire .pdb file loaded into my own data structure.</td>
  </tr>
  <tr>
    <td class="tg-0pky">d_terministring</td>
    <td class="tg-0pky">String consisting of two letters, specifying how gmx pdb2gmx should treat the termini of the (various chains of) the protein. Commonly used: 00 = charged termini (pdb2gmx default), 11 = neutral termini, 34 = use this with capped tripeptide.</td>
  </tr>
  <tr>
    <td class="tg-0pky">d_title</td>
    <td class="tg-0pky">TITLE of .pdb file. If none was set, d_pdbName will be used.</td>
  </tr>
  <tr>
    <td class="tg-0pky">ph_ASP_dvdl</td>
    <td class="tg-0pky">Contains dV/dl calibration coefficients for ASP. Has to be added manually to universe before you start to build the system.</td>
  </tr>
  <tr>
    <td class="tg-0pky">ph_BUF_dvdl</td>
    <td class="tg-0pky">Containd dV/dl calibration coefficients for BUF. Has to be added manually to universe before you start to build the system.</td>
  </tr>
  <tr>
    <td class="tg-0pky">ph_GLU_dvdl</td>
    <td class="tg-0pky">Containd dV/dl calibration coefficients for GLU. Has to be added manually to universe before you start to build the system.</td>
  </tr>
  <tr>
    <td class="tg-0pky">ph_barrier_E</td>
    <td class="tg-0pky">Barrier energy (kJ/mol). Good value is 5.0 or 7.5 kJ/mol.</td>
  </tr>
  <tr>
    <td class="tg-0pky">ph_bufMargin</td>
    <td class="tg-0pky">Minimum distance (nm) between the buffer molecule(s) (themselves) and the protein. This should be at least 1.5nm to avoid unphysical interactions.</td>
  </tr>
  <tr>
    <td class="tg-0pky">ph_bufitpname</td>
    <td class="tg-0pky">Name of the topology file of the buffer molecule.</td>
  </tr>
  <tr>
    <td class="tg-0pky">ph_bufnmol</td>
    <td class="tg-0pky">Number of buffer molecules.</td>
  </tr>
  <tr>
    <td class="tg-0pky">ph_pdbname</td>
    <td class="tg-0pky">Name of the structure file of the buffer molecule.</td>
  </tr>
  <tr>
    <td class="tg-0pky">ph_constantpH</td>
    <td class="tg-0pky">Turn constant-pH on or off. Has to be added manually to universe before you start to build the system.</td>
  </tr>
  <tr>
    <td class="tg-0pky">ph_dvdl_initList</td>
    <td class="tg-0pky">List of lambda coordinates for which to compute mean dV/dl during calibration.</td>
  </tr>
  <tr>
    <td class="tg-0pky">ph_dvdl_meanList</td>
    <td class="tg-0pky">List of mean dV/dl values at said points.</td>
  </tr>
  <tr>
    <td class="tg-0pky">ph_dvdl_stdList</td>
    <td class="tg-0pky">List of standard deviations of dV/dl values at said points.</td>
  </tr>  
  <tr>
    <td class="tg-0pky">ph_lambdaM</td>
    <td class="tg-0pky">Mass (Da) of lambda particle. Good value is 5.0 Da.</td>
  </tr>
  <tr>
    <td class="tg-0pky">ph_nstout</td>
    <td class="tg-0pky">Similar to d_nstxout. Interval in which to output a line of lambda data to the various lambda_x.dat files.</td>
  </tr>
  <tr>
    <td class="tg-0lax">ph_pH</td>
    <td class="tg-0lax">pH of the simulation/buffer.</td>
  </tr>
  <tr>
    <td class="tg-0lax">ph_restrainpH</td>
    <td class="tg-0lax">Turn the charge-leveling scheme "charge-restraining" on or off. Has to be added manually to universe before you start to build the system.</td>
  </tr>
</tbody>
</table>

<b>Known issues</b>
* You can currently use `protein.add_buf()` to add less buffer molecules than you have titratable groups, however for some reason (either phbuilder or gromacs build), this will give nonsensical results.

<b>To-do</b>
* Enable the charge leveling scheme termed "charge coupling", i.e. put charge on an "ion".
* Enable multistate.
