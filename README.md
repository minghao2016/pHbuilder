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
    <td class="tg-0pky">MD time step size. Currently hardcoded at 0.002 ps.</td>
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
    <td class="tg-0pky">Barrier energy in kJ/mol. Good value is 5.0 or 7.5 kJ/mol.</td>
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
    <td class="tg-0pky">Mass of lambda particle. Good value is 5.0 Da.</td>
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
* You can currently use add_buf() to add less buffer molecules than you have titratable groups, however for some reason (either phbuilder or gromacs build), this does not work yet.

<b>To-do</b>
* Enable the charge leveling scheme termed "charge coupling", i.e. put charge on an "ion".
* Put the lambda charges on virtual sites instead of the actual atoms themselves.
* Enable multistate.
