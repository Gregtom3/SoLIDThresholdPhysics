# SoLIDThresholdPhysics
Generators for studying vector meson production at SoLID

The project at this stage consists of the following generators...

1. (TO BE IMPLEMENTED) Bethe-Heitler e-e+ production off a proton target
2. Coherent Bethe-Heitler e-e+ production off a deuteron target
3. Coherent exclusive vector meson production of J/Psi off a deuteron target
4. (TO BE IMPLEMENTED) Coherent exclusive vector meson production of J/Psi off a proton target
5. (TO BE IMPLEMENTED) Incoherent exclusive vector meson production of J/Psi off a deuteron target (p/n)
6. (TO BE IMPLEMENTED) Exclusive vector meson productions of Phi and Psi(2S)

# Instructions
The generator can be run remotely with the proper python3 packages installed. However, for users wishing to run parallel computing jobs on Jefferson Lab's *ifarm* , one must have a computing account already set up. The program also assumes you have PYROOT enabled in your python3 installation.

1. Clone the repository
```
git clone https://github.com/Gregtom3/SoLIDThresholdPhysics
```
2. Install the following python3 packages
```
pip3 install tqdm numpy matplotlib scipy uproot
```
3. To run either the Bethe-Heitler or DVMP event generators, a runcard must be created. Sample runcards are given in **./runcards/** for Bethe-Heitler, Electroproduction, and Photoproduction at SoLID. To use a runcard on the event generator, one must be in the home directory of the repository. As an example, to run the BH simulation...
```
python3 ./evtgen/run_bh.py --runcard ./runcards/ed_bh_8.8GeV.card
```
Similarly, for running the exclusive vector meson production generators...
```
python3 ./evtgen/run_dvmp.py --runcard ./runcards/ed_dvmp_8.8GeV_electroproduction.card
```
The **run_bh.py** and **run_dvmp.py** scripts will event-by-event generate the corresponding physics processes, caclulating event weights, and smearing final state particle kinematics. The results of the simulation are stored in TTrees within the `output_file_location` directory. Two TTrees are saved for each simulation. The first, titled `tree`, assumes perfect 4pi acceptance and saves each generated event. The second, titled `treeAcc` contains fewer events, and only those which passed the 3-fold coincidence demanded by the production process.

4. If you wish to execute multiple simulations on *ifarm*, please run the following script from the _SoLIDThresholdPhysics_ directory and follow the prompts
```
bash ./scripts/create_project.sh
```
This will create multiple runcards within a new subdirectory of _./runcards/_ which can be run using their corresponding slurm scripts in a new subdirectory of _./scripts/_ with the **sbatch** command (again, from the _SoLIDThresholdPhysics_ directory)

# Note

1. jpsi-d coherent photoproduction is checked, jpsi-d electroproduction doesn't run somehow, BH-d photoproduction is not checked (by Zhiwen, 2024/04/10)
2. degree of freedom and phasespace which the generator chooses. this choice means there is no physics forbidden region in the simple phasespace
```
   jpsi photoproduction
        solid angle 4pi from gamma+d-->J/Psi+d' in CM frame
        solid angle 4pi from J/Psi-->ee in CM frame
        Egmax - Egmin from real photon of EPA and Bremsstrahlung
   jpsi electroproduction
        solid angle 4pi from gamma+d-->J/Psi+d' in CM frame
        solid angle 4pi from J/Psi-->ee in CM frame
        Egmax - Egmin from virtual photon of scattered e-
        solid angle 4pi from scattered e-
```
