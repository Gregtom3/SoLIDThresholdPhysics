#!/bin/bash
USER=$USER
PWD=$PWD
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Check if the user is inside the "SoLIDThresholdPhysics" directory
current_dir=$(basename "$PWD")

if [ "$current_dir" != "SoLIDThresholdPhysics" ]; then
    echo -e "${RED}Please run this script inside the 'SoLIDThresholdPhysics' directory. Exiting.${NC}"
    exit 1
fi

# Prompt for simulation name
read -p "Enter the name of the simulation: " simulation_name

# Check if the project already exists
base_dir="/volatile/clas12/users/$USER/solid.data/projects/"
sim_dir="$base_dir/$simulation_name"

if [ ! -d "$base_dir" ]; then
    echo "The directory '$base_dir' does not exist. Please create it and try again."
    exit 1
fi

if [ -d "$sim_dir" ]; then
    echo -e "${RED}The project '$simulation_name' already exists.${NC}"
    read -p "Do you want to delete the existing project? (y/n): " confirm_delete
    case $confirm_delete in
        [Yy]* ) 
            rm -rf "$sim_dir"
            echo -e "${GREEN}Deleted the existing project '$simulation_name'.${NC}"
            ;;
        [Nn]* ) 
            echo "Exiting."
            exit 1
            ;;
        * ) 
            echo "Invalid choice. Exiting."
            exit 1
            ;;
    esac
fi

# Create BH directories
data_dir_bh="$sim_dir/bh/data"
out_dir_bh="$sim_dir/bh/out"
err_dir_bh="$sim_dir/bh/err"

mkdir -p "$out_dir_bh" "$err_dir_bh" "$data_dir_bh"

echo -e "${GREEN}Created directories:${NC}"
echo -e "${GREEN}$data_dir_bh${NC}"
echo -e "${GREEN}$out_dir_bh${NC}"
echo -e "${GREEN}$err_dir_bh${NC}"

# Create electroproduction directories
data_dir_electro="$sim_dir/electroproduction/data"
out_dir_electro="$sim_dir/electroproduction/out"
err_dir_electro="$sim_dir/electroproduction/err"

mkdir -p "$out_dir_electro" "$err_dir_electro" "$data_dir_electro"

echo -e "${GREEN}Created directories for Electroproduction:${NC}"
echo -e "${GREEN}$data_dir_electro${NC}"
echo -e "${GREEN}$out_dir_electro${NC}"
echo -e "${GREEN}$err_dir_electro${NC}"

# Create photoproduction directories
data_dir_photo="$sim_dir/photoproduction/data"
out_dir_photo="$sim_dir/photoproduction/out"
err_dir_photo="$sim_dir/photoproduction/err"

mkdir -p "$out_dir_photo" "$err_dir_photo" "$data_dir_photo"

echo -e "${GREEN}Created directories for Photoproduction:${NC}"
echo -e "${GREEN}$data_dir_photo${NC}"
echo -e "${GREEN}$out_dir_photo${NC}"
echo -e "${GREEN}$err_dir_photo${NC}"

# Prompt the user for runcard details
echo "************** GENERAL RUNCARD INFO **************"

#read -p "num_events: " num_events
read -p "beam_energy (GeV): " beam_energy
read -p "beam_current (uA): " beam_current
read -p "target_length (cm): " target_length
read -p "target_type ('p' or 'd'): " target_type
read -p "days: " days
read -p "detector ('SoLID' or 'CLAS12'): " detector

echo "************** BETHE-HEITLER INFO **************"
read -p "num_events (For BH simulation): " num_events_bh
read -p "num_batches (For BH simulation): " num_batches_bh
read -p "photon_energy_min (GeV): " photon_energy_min_bh
read -p "photon_energy_max (GeV): " photon_energy_max_bh
read -p "Mll_min (GeV): " Mll_min
read -p "Mll_max (GeV): " Mll_max

# Calculate the squared quantities
Mll2_min_bh=$(echo "$Mll_min^2" | bc)
Mll2_max_bh=$(echo "$Mll_max^2" | bc)

while true; do
    read -p "t_min (GeV^2): " t_min_bh
    read -p "t_max (GeV^2): " t_max_bh

    if (( t_min_bh <= 0 )) && (( t_max_bh <= 0 )) && (( t_min_bh < t_max_bh )); then
        break
    else
        echo "Please make sure both t_min and t_max are negative (or zero) and that t_min < t_max."
    fi
done

read -p "Q2_max (GeV^2): " Q2_max_bh

echo "************** DVMP INFO **************"
read -p "num_events (For DVMP Photoproduction simulation): " num_events_photo
read -p "num_batches (For DVMP Photoproduction simulation): " num_batches_photo
read -p "num_events (For DVMP Electroproduction simulation): " num_events_electro
read -p "num_batches (For DVMP Electroproduction simulation): " num_batches_electro
read -p "model_type ('PomeronLQCD' or '23g'): " model_type_dvmp
read -p "photon_energy_min (GeV): " photon_energy_min_dvmp
read -p "photon_energy_max (GeV): " photon_energy_max_dvmp
read -p "Q2_max (GeV^2): " Q2_max_dvmp
read -p "t_min (GeV^2):" t_min_dvmp


# Create the runcard directory if it doesn't exist
runcard_dir="./runcards/$simulation_name"
mkdir -p "$runcard_dir"

##################################################################
# Generate the runcard for BETHE-HEITLER
runcard_path_bh="$runcard_dir/bh.card"

cat > "$runcard_path_bh" <<EOF
# Runcard for DVMP generator
#############################
num_events: $num_events_bh
output_file_location: $data_dir_bh
output_file_prefix: bh

# Experiment configuration
#############################
beam_energy: $beam_energy
beam_current: $beam_current
target_length: $target_length
target_type: $target_type
days: $days
detector: $detector

# Kinematic Limits
#############################
photon_energy_min: $photon_energy_min_bh
photon_energy_max: $photon_energy_max_bh

Mll2_min: $Mll2_min_bh
Mll2_max: $Mll2_max_bh

t_min: $t_min_bh
t_max: $t_max_bh
Q2_max: $Q2_max_bh
EOF

##################################################################
# Generate the runcard for PHOTOPRODUCTION
runcard_path_photoproduction="$runcard_dir/photoproduction.card"

cat > "$runcard_path_photoproduction" <<EOF
# Runcard for DVMP generator
#############################
num_events: $num_events_photo
output_file_location: $data_dir_photo
output_file_prefix: photo
process: photoproduction
model_type: $model_type_dvmp

# Experiment configuration
#############################
beam_energy: $beam_energy
beam_current: $beam_current
target_length: $target_length
target_type: $target_type
days: $days
detector: $detector

# Kinematic Limits
#############################
photon_energy_min: $photon_energy_min_dvmp
photon_energy_max: $photon_energy_max_dvmp
Q2_max: $Q2_max_dvmp
t_min: $t_min_dvmp

EOF

##################################################################
# Generate the runcard for ELECTROPRODUCTION
runcard_path_electroproduction="$runcard_dir/electroproduction.card"

cat > "$runcard_path_electroproduction" <<EOF
# Runcard for DVMP generator
#############################
num_events: $num_events_electro
output_file_location: $data_dir_electro
output_file_prefix: electro
process: electroproduction
model_type: $model_type_dvmp

# Experiment configuration
#############################
beam_energy: $beam_energy
beam_current: $beam_current
target_length: $target_length
target_type: $target_type
days: $days
detector: $detector

# Kinematic Limits
#############################
photon_energy_min: $photon_energy_min_dvmp
photon_energy_max: $photon_energy_max_dvmp
Q2_max: $Q2_max_dvmp
t_min: $t_min_dvmp
EOF

# Create slurm folder
mkdir scripts/$simulation_name/


# Generate the Slurm file for BETHE-HEITLER
cat > ./scripts/$simulation_name/bh.slurm <<EOF
#!/bin/bash
#SBATCH --account=clas12
#SBATCH --partition=production
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=bh_evtgen_${simulation_name}
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --output=$out_dir_bh/%x-%j-%N.out
#SBATCH --error=$err_dir_bh/%x-%j-%N.err
#SBATCH --array=0-$(($num_batches_bh-1))
/apps/python3/3.9.7/bin/python3 evtgen/run_bh.py --runcard $runcard_path_bh --batch \$SLURM_ARRAY_TASK_ID
EOF

echo -e "${GREEN}Slurm file ./scripts/$simulation_name/bh.slurm has been generated.${NC}"




# Generate the Slurm file for PHOTOPRODUCTION
cat > ./scripts/$simulation_name/photoproduction.slurm <<EOF
#!/bin/bash
#SBATCH --account=clas12
#SBATCH --partition=production
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=photoproduction_evtgen_${simulation_name}
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --output=$out_dir_photo/%x-%j-%N.out
#SBATCH --error=$err_dir_photo/%x-%j-%N.err
#SBATCH --array=0-$(($num_batches_photo-1))

/apps/python3/3.9.7/bin/python3 evtgen/run_dvmp.py --runcard $runcard_path_photoproduction --batch \$SLURM_ARRAY_TASK_ID
EOF

echo -e "${GREEN}Slurm file ./scripts/$simulation_name/photoproduction.slurm has been generated.${NC}"


# Generate the Slurm file for ELECTROPRODUCTION
cat > ./scripts/$simulation_name/electroproduction.slurm <<EOF
#!/bin/bash
#SBATCH --account=clas12
#SBATCH --partition=production
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=electroproduction_evtgen_${simulation_name}
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --output=$out_dir_electro/%x-%j-%N.out
#SBATCH --error=$err_dir_electro/%x-%j-%N.err
#SBATCH --array=0-$(($num_batches_electro-1))

/apps/python3/3.9.7/bin/python3 evtgen/run_dvmp.py --runcard $runcard_path_electroproduction --batch \$SLURM_ARRAY_TASK_ID
EOF

echo -e "${GREEN}Slurm file ./scripts/$simulation_name/electroproduction.slurm has been generated.${NC}"


