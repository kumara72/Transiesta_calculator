#!/usr/bin/env python3

import os
import sys
import subprocess
import shutil
import re
from textwrap import dedent
from string import Template

# libraries
try:
    import matplotlib
    matplotlib.use('Agg') 
    import matplotlib.pyplot as plt
    import numpy as np
except ImportError:
    print("Error: Matplotlib and Numpy are required.")
    sys.exit(1)

# PLOT SCRIPT
PLOT_SCRIPT_CONTENT = """
import sisl
import matplotlib.pyplot as plt
import sys
import os

# --- Handle Arguments from the Automation Script ---
# The automation script passes the filename as the first argument
if len(sys.argv) < 2:
    tbt_file = "siesta.TBT.nc" # Default if run manually
    label_text = tbt_file
else:
    tbt_file = sys.argv[1]
    label_text = sys.argv[2] if len(sys.argv) > 2 else tbt_file

# --- Setup the plot ---
plt.figure(figsize=(8, 6))

# --- Load the file ---
try:
    if not os.path.exists(tbt_file):
        print(f"Error: The file '{tbt_file}' was not found.")
        sys.exit(0)

    tbt = sisl.get_sile(tbt_file)

    # --- Plot the transmission for that file ---
    # Using 'b-' for a solid blue line as requested
    plt.plot(tbt.E, tbt.transmission(), label=label_text, color='b', linestyle='-')

    # --- Labels and formatting ---
    plt.legend()
    plt.ylabel('Transmission')
    plt.xlabel('Energy [eV]')
    plt.ylim([0, None])  # Start y-axis at 0
    plt.title("Transmission Spectrum")
    plt.grid(True)

    # --- Save the figure as a PNG file ---
    # We create a unique name based on the input label so they don't overwrite each other
    clean_label = label_text.replace(' ', '_').replace('=', '').replace('.', 'p')
    output_filename = f"transmission_{clean_label}.png"
    
    plt.savefig(output_filename)

    print(f"Success! Graph saved as '{output_filename}'")

except FileNotFoundError:
    print(f"Error: The file '{tbt_file}' was not found.")
except Exception as e:
    print(f"An error occurred: {e}")
"""

def create_plot_script(directory):
    path = os.path.join(directory, "plot_transmission.py")
    with open(path, "w") as f: f.write(PLOT_SCRIPT_CONTENT)

# HELPER
def clean_file_remove_kgrid(input_path, output_path):
    try:
        with open(input_path, 'r') as f_in: lines = f_in.readlines()
        with open(output_path, 'w') as f_out:
            skipping = False
            for line in lines:
                check = line.lower().strip()
                if check.startswith("%block kgrid"): skipping = True; continue
                if skipping and check.startswith("%endblock kgrid"): skipping = False; continue
                if not skipping: f_out.write(line)
            f_out.write("\n")
        return True
    except Exception as e:
        print(f"Error cleaning file: {e}"); sys.exit(1)

def extract_kpoints_block_from_file(fdf_path):
    try:
        with open(fdf_path, 'r') as f: lines = f.readlines()
        block = []; in_b = False; found = False
        for line in lines:
            check = line.lower().strip()
            if check.startswith("%block kgrid"): in_b = True; found = True; block.append("%block kgrid_Monkhorst_Pack"); continue
            if in_b and check.startswith("%endblock kgrid"): in_b = False; block.append("%endblock kgrid_Monkhorst_Pack"); break
            if in_b: block.append(line.strip())
        if found: return "\n".join(block)
        return None
    except: return None

def format_manual_k_block(k_str):
    try:
        k = [int(x) for x in k_str.split()]
        return dedent(f"""
            %block kgrid_Monkhorst_Pack
            {k[0]} 0 0 0.0
            0 {k[1]} 0 0.0
            0 0 {k[2]} 0 0.0
            %endblock kgrid_Monkhorst_Pack
            """).strip()
    except:
        return dedent("""
            %block kgrid_Monkhorst_Pack
            1 0 0 0.0
            0 1 0 0.0
            0 0 1 0.0
            %endblock kgrid_Monkhorst_Pack""").strip()

def create_transiesta_run_file(device_geom_path, output_path):
    print(f"--- Using user provided device file: {device_geom_path} ---")
    try:
        shutil.copy(device_geom_path, output_path)
    except Exception as e:
        print(f"Error copying device file: {e}"); sys.exit(1)

# TEMPLATE

# 1. MESH
MESH_SCRIPT_TEMPLATE = Template("""
#!/bin/bash
echo "--- Starting: Mesh Cutoff Convergence Test ---"
module load siesta || echo "Warning: 'module load siesta' failed."

CLEAN_FDF="$clean_fdf_name"
MESH_SEQUENCE="$mesh_sequence"
SIESTA_COMMAND="$siesta_command"
TBTRANS_COMMAND="$tbtrans_command"

rm -rf cont && mkdir cont
for mesh in $$MESH_SEQUENCE
do
    DIR_NAME="mesh_$${mesh}Ry"
    echo "--- Preparing: $$DIR_NAME ---"
    mkdir -p $$DIR_NAME
    [ -f cont/siesta.DM ] && cp cont/siesta.DM $$DIR_NAME/
    cd $$DIR_NAME
    
    cp ../*psf . 2>/dev/null
    cp ../electrodel.fdf . 2>/dev/null
    cp ../electroder.fdf . 2>/dev/null
    cp ../*.TSHS . 2>/dev/null
    cp ../*.TSDE . 2>/dev/null
    cp ../plot_transmission.py .
    
    cp ../$$CLEAN_FDF siesta.fdf
    sed -i "s/^[[:space:]]*SystemLabel.*/SystemLabel siesta/" siesta.fdf
    sed -i "s/^[[:space:]]*MeshCutoff.*/MeshCutoff $$mesh Ry/" siesta.fdf
    
    # REMOVE OLD SETTINGS TO AVOID CONFLICTS
    sed -i "/SCF.Mixer.Weight/d" siesta.fdf
    sed -i "/SCF.Mixer.History/d" siesta.fdf
    sed -i "/MaxSCFIterations/d" siesta.fdf
    
cat <<EOF >> siesta.fdf

MaxSCFIterations 2000
SCF.Mixer.Weight    0.1
SCF.Mixer.History   15
$fixed_k_block
DM.UseSaveDM .true.
EOF
    
    echo "--- Running Siesta (Mesh: $$mesh) ---"
    $$SIESTA_COMMAND < siesta.fdf > siesta.out
    
    # CHECK AND RUN TBTRANS
    if [ -f siesta.TSHS ]; then
        echo "--- Running TBTrans ---"
        $$TBTRANS_COMMAND < siesta.fdf > siesta.tbt.out
        if [ -f siesta.TBT.nc ]; then
            python3 plot_transmission.py siesta.TBT.nc "Mesh=$${mesh}Ry"
        fi
    else
        echo "!!! ERROR: siesta.TSHS not found. Siesta failed. !!!"
        tail -n 15 siesta.out
    fi
    
    cd ..
    rm -rf cont && mkdir cont
    if [ -f ./$$DIR_NAME/siesta.DM ]; then
        cp ./$$DIR_NAME/siesta.DM cont/
    else
        echo "--- WARNING: siesta.DM missing. ---"
    fi
done
""")

# 2. KPOINT
KPOINT_SCRIPT_TEMPLATE = Template("""
#!/bin/bash
echo "--- Starting: K-Point Convergence Test ---"
module load siesta || echo "Warning: 'module load siesta' failed."

CLEAN_FDF="$clean_fdf_name"
KPOINT_SEQUENCE="$kpoint_sequence"
FIXED_MESH_CUTOFF="$fixed_mesh_cutoff"
SIESTA_COMMAND="$siesta_command"
TBTRANS_COMMAND="$tbtrans_command"

rm -rf cont_kpoint && mkdir cont_kpoint
for k in $$KPOINT_SEQUENCE
do
    DIR_NAME="$kpoint_test_dir_name"
    echo "--- Preparing: $$DIR_NAME ---"
    mkdir -p $$DIR_NAME
    [ -f cont_kpoint/siesta.DM ] && cp cont_kpoint/siesta.DM $$DIR_NAME/
    cd $$DIR_NAME
    cp ../*psf . 2>/dev/null
    cp ../electrodel.fdf . 2>/dev/null
    cp ../electroder.fdf . 2>/dev/null
    cp ../*.TSHS . 2>/dev/null
    cp ../*.TSDE . 2>/dev/null
    cp ../plot_transmission.py .
    
    cp ../$$CLEAN_FDF siesta.fdf
    sed -i "s/^[[:space:]]*SystemLabel.*/SystemLabel siesta/" siesta.fdf
    sed -i "s/^[[:space:]]*MeshCutoff.*/MeshCutoff $$FIXED_MESH_CUTOFF Ry/" siesta.fdf

    sed -i "/SCF.Mixer.Weight/d" siesta.fdf
    sed -i "/SCF.Mixer.History/d" siesta.fdf
    sed -i "/MaxSCFIterations/d" siesta.fdf
    
    echo "MaxSCFIterations 2000" >> siesta.fdf
    echo "SCF.Mixer.Weight 0.1" >> siesta.fdf
    echo "SCF.Mixer.History 15" >> siesta.fdf
    
    echo "" >> siesta.fdf
    $kpoint_echo_command >> siesta.fdf
    echo "" >> siesta.fdf
    echo "DM.UseSaveDM .true." >> siesta.fdf
    
    echo "--- Running Siesta (K: Varying) ---"
    $$SIESTA_COMMAND < siesta.fdf > siesta.out
    
    if [ -f siesta.TSHS ]; then
        echo "--- Running TBTrans ---"
        $$TBTRANS_COMMAND < siesta.fdf > siesta.tbt.out
        if [ -f siesta.TBT.nc ]; then
            python3 plot_transmission.py siesta.TBT.nc "Kpoint=$$k"
        fi
    else
        echo "!!! ERROR: siesta.TSHS not found. !!!"
        tail -n 15 siesta.out
    fi

    cd ..
    rm -rf cont_kpoint && mkdir cont_kpoint
    if [ -f ./$$DIR_NAME/siesta.DM ]; then
        cp ./$$DIR_NAME/siesta.DM cont_kpoint/
    else
        echo "--- WARNING: siesta.DM missing. ---"
    fi
done
""")

# 3. LATTICE
LATTICE_SCRIPT_TEMPLATE = Template("""
#!/bin/bash
echo "--- Starting: Z-Lattice Convergence Test ---"
module load siesta || echo "Warning: 'module load siesta' failed."

CLEAN_FDF="$clean_fdf_name"
LATTICE_SEQUENCE="$lattice_sequence"
FIXED_MESH_CUTOFF="$fixed_mesh_cutoff"
SIESTA_COMMAND="$siesta_command"
TBTRANS_COMMAND="$tbtrans_command"
FIXED_K_BLOCK='$fixed_k_block'

rm -rf cont_lattice && mkdir cont_lattice
for lat_z in $$LATTICE_SEQUENCE
do
    DIR_NAME="lattice_z_$${lat_z}A"
    echo "--- Preparing: $$DIR_NAME ---"
    mkdir -p $$DIR_NAME
    [ -f cont_lattice/siesta.DM ] && cp cont_lattice/siesta.DM $$DIR_NAME/
    cd $$DIR_NAME
    cp ../*psf . 2>/dev/null
    cp ../electrodel.fdf . 2>/dev/null
    cp ../electroder.fdf . 2>/dev/null
    cp ../*.TSHS . 2>/dev/null
    cp ../*.TSDE . 2>/dev/null
    cp ../plot_transmission.py .
    
    cp ../$$CLEAN_FDF siesta.fdf
    
    sed -i "s/^[[:space:]]*SystemLabel.*/SystemLabel siesta/" siesta.fdf
    sed -i "s/^[[:space:]]*MeshCutoff.*/MeshCutoff $$FIXED_MESH_CUTOFF Ry/" siesta.fdf
    
cat <<EOF >> siesta.fdf

$fixed_k_block

EOF
    
    sed -i "s/^[[:space:]]*0.00000000[[:space:]]*0.00000000[[:space:]].*# Z-vector/  0.00000000  0.00000000  $$lat_z # Z-vector/" siesta.fdf

    echo "DM.UseSaveDM .true." >> siesta.fdf
    echo "--- Running Siesta (Lat: $$lat_z) ---"
    $$SIESTA_COMMAND < siesta.fdf > siesta.out
    
    if [ -f siesta.TSHS ]; then
        echo "--- Running TBTrans ---"
        $$TBTRANS_COMMAND < siesta.fdf > siesta.tbt.out
        if [ -f siesta.TBT.nc ]; then
            python3 plot_transmission.py siesta.TBT.nc "Lat=$$lat_z"
        fi
    fi

    cd ..
    rm -rf cont_lattice && mkdir cont_lattice
    if [ -f ./$$DIR_NAME/siesta.DM ]; then
        cp ./$$DIR_NAME/siesta.DM cont_lattice/
    else
        echo "--- WARNING: siesta.DM missing. ---"
    fi
done
""")

# 4. FINAL CALC
FINAL_CALC_SCRIPT_TEMPLATE = Template("""
#!/bin/bash
echo "--- Starting: Final Converged Calculation ---"
module load siesta || echo "Warning: 'module load siesta' failed."

CLEAN_FDF="$clean_fdf_name"
CONV_MESH="$conv_mesh"
CONV_LAT_Z="$conv_lat_z"
SIESTA_COMMAND="$siesta_command"
TBTRANS_COMMAND="$tbtrans_command"
SYSTEM_LABEL="$system_label"
DO_LATTICE_CONV="$do_lattice_conv"
FIXED_K_BLOCK='$fixed_k_block'

DIR_NAME="optimize"
mkdir -p $$DIR_NAME
cd $$DIR_NAME

cp ../*psf . 2>/dev/null
cp ../electrodel.fdf . 2>/dev/null
cp ../electroder.fdf . 2>/dev/null
cp ../*.TSHS . 2>/dev/null
cp ../*.TSDE . 2>/dev/null
cp ../plot_transmission.py .

cp ../$$CLEAN_FDF siesta.fdf
sed -i "s/^[[:space:]]*MeshCutoff.*/MeshCutoff $$CONV_MESH Ry/" siesta.fdf
sed -i "s/^[[:space:]]*SystemLabel.*/SystemLabel siesta/" siesta.fdf

if [ "$$DO_LATTICE_CONV" = "true" ]; then
    sed -i "s/^[[:space:]]*0.00000000[[:space:]]*0.00000000[[:space:]].*# Z-vector/  0.00000000  0.00000000  $$CONV_LAT_Z # Z-vector/" siesta.fdf
fi

sed -i "/SCF.Mixer.Weight/d" siesta.fdf
sed -i "/SCF.Mixer.History/d" siesta.fdf
sed -i "/MaxSCFIterations/d" siesta.fdf

cat <<EOF >> siesta.fdf

MaxSCFIterations 2000
SCF.Mixer.Weight    0.1
SCF.Mixer.History   15
$fixed_k_block

EOF

echo "--- Running Siesta ---"
$$SIESTA_COMMAND < siesta.fdf > siesta.out

if [ -f siesta.TSHS ]; then
    echo "--- Running TBTrans ---"
    $$TBTRANS_COMMAND < siesta.fdf > siesta.tbt.out
    if [ -f siesta.TBT.nc ]; then
        python3 plot_transmission.py siesta.TBT.nc "Device_0.0V"
    fi
else
    echo "!!! CRITICAL: Siesta failed to produce TSHS !!!"
    tail -n 20 siesta.out
    exit 1
fi

if [ -f siesta.TSHS ] || [ -f siesta.DM ]; then
    echo "--- Renaming files ---"
    for f in siesta.*; do mv "$$f" "$${SYSTEM_LABEL}.$${f#siesta.}"; done
fi
""")

# 5. VOLTAGE SWEEP (FIXED: Copies TSHS/DM/TSDE for device & electrodes from PREVIOUS folder)
VOLTAGE_SWEEP_SCRIPT_TEMPLATE = Template("""
#!/bin/bash
echo "--- Starting: Voltage Sweep Workflow ---"
module load siesta || echo "Warning: 'module load siesta' failed."

START_DIR=$$(pwd)
SUCCESS_MSG="Job completed"
SIESTA_COMMAND="mpirun -np $mpi_procs siesta"
TBTRANS_COMMAND="$tbtrans_command"

# Variables passed from Python to guarantee correctness
OPTIMIZED_MESH="$opt_mesh"
OPTIMIZED_K_BLOCK='$opt_k_block'

# Uses the optimized FDF provided by the main script
TEMPLATE_FDF="$template_fdf_name"
PSF_FILES="*.psf"

print_step() {
    echo ""
    echo "=================================================================="
    echo "  $$1"
    echo "=================================================================="
}

print_step "1. Running 0V Device Calculation"
VOLT_DIR="volt_0.0V"
LOG_FILE="device_0.0V.out"
mkdir -p $$VOLT_DIR
cd $$VOLT_DIR

cp ../$$TEMPLATE_FDF .
cp ../$$PSF_FILES . 2>/dev/null
cp ../electrodel.TSHS .
cp ../electrodel.TSDE .
cp ../electroder.TSHS .
cp ../electroder.TSDE .
cp ../electrodel.fdf . 2>/dev/null
cp ../electroder.fdf . 2>/dev/null
cp ../plot_transmission.py .

# --- COPY INPUTS STRICTLY ---
# Copy the optimized TSHS and DM from parent folder (setup by python script)
cp ../siesta.TSHS .
if [ -f ../siesta.DM ]; then cp ../siesta.DM .; fi

# --- FORCE OPTIMIZED PARAMETERS ---
sed -i "s/^[[:space:]]*MeshCutoff.*/MeshCutoff $$OPTIMIZED_MESH Ry/" $$TEMPLATE_FDF

# Append the correct K-block
cat <<EOF >> $$TEMPLATE_FDF

$opt_k_block

EOF

# Ensure SystemLabel is siesta
sed -i "/SystemLabel/d" $$TEMPLATE_FDF
sed -i "/SystemName/d" $$TEMPLATE_FDF
echo "SystemLabel siesta" >> $$TEMPLATE_FDF

# Ensure DM saving is ON
sed -i "/DM.UseSaveDM/d" $$TEMPLATE_FDF
echo "DM.UseSaveDM .true." >> $$TEMPLATE_FDF

# Ensure Voltage is 0.0 eV
sed -i "/TS.Voltage/d" $$TEMPLATE_FDF
echo "" >> $$TEMPLATE_FDF
echo "TS.Voltage 0.0 eV" >> $$TEMPLATE_FDF

echo "Running 0V siesta..."
$$SIESTA_COMMAND < $$TEMPLATE_FDF > $$LOG_FILE

if grep -q "$$SUCCESS_MSG" "$$LOG_FILE"; then
    echo "SUCCESS: 0V Device calculation completed."
    
    if [ -f siesta.TSHS ]; then 
        echo "Running 0V TBTrans..."
        $$TBTRANS_COMMAND < $$TEMPLATE_FDF > tbtrans_0.0V.out
        
        if [ -f siesta.TBT.nc ]; then
            mv siesta.TBT.nc device_0.0V.TBT.nc
            echo "Plotting 0.0V transmission..."
            python3 plot_transmission.py device_0.0V.TBT.nc "Voltage=0.0V"
        fi
    else
        echo "!!! ERROR: siesta.TSHS missing after 0V run !!!"
        exit 1
    fi
else
    print_step "!!! ERROR: 0V DEVICE FAILED !!!"
    tail -n 10 $$LOG_FILE
    exit 1
fi
cd $$START_DIR

VOLTAGES=($voltage_list)
PREV_VOLTAGES=("0.0" $${VOLTAGES[@]::$$(( $${#VOLTAGES[@]} - 1 ))})

for i in $${!VOLTAGES[@]}; do
    V=$${VOLTAGES[$$i]}
    PREV_V=$${PREV_VOLTAGES[$$i]}
    VOLT_DIR="volt_$${V}V"
    PREV_VOLT_DIR="volt_$${PREV_V}V"
    LOG_FILE="device_$${V}V.out"
    
    print_step "2. Running $$V V Device Calculation"
    
    mkdir -p $$VOLT_DIR
    cd $$VOLT_DIR

    # 1. Basic files from ROOT (FDF, PSF, Plotting)
    cp ../$$TEMPLATE_FDF .
    cp ../$$PSF_FILES . 2>/dev/null
    cp ../plot_transmission.py .
    
    # 2. VITAL: Copy Restart Files from PREVIOUS VOLTAGE Folder
    SOURCE_DIR="../volt_$${PREV_V}V"
    echo "--- Copying TSHS, DM, and TSDE from $$SOURCE_DIR ---"
    
    # Copy Core files
    cp $$SOURCE_DIR/siesta.DM .
    cp $$SOURCE_DIR/siesta.TSHS .
    
    # Copy Electrodes (TSHS and TSDE) from previous run
    cp $$SOURCE_DIR/electrodel.TSHS .
    cp $$SOURCE_DIR/electroder.TSHS .
    
    if [ -f $$SOURCE_DIR/electrodel.TSDE ]; then cp $$SOURCE_DIR/electrodel.TSDE .; fi
    if [ -f $$SOURCE_DIR/electroder.TSDE ]; then cp $$SOURCE_DIR/electroder.TSDE .; fi

    # Copy Siesta TSDE if available (Device Delta Energy)
    if [ -f $$SOURCE_DIR/siesta.TSDE ]; then cp $$SOURCE_DIR/siesta.TSDE .; fi

    # 3. Force Optimized Parameters (Mesh/K-Points) onto the template again
    sed -i "s/^[[:space:]]*MeshCutoff.*/MeshCutoff $$OPTIMIZED_MESH Ry/" $$TEMPLATE_FDF
    cat <<EOF >> $$TEMPLATE_FDF

$opt_k_block

EOF

    # 4. Update FDF for new Voltage (using eV)
    sed -i "/SystemLabel/d" $$TEMPLATE_FDF; echo "SystemLabel siesta" >> $$TEMPLATE_FDF
    sed -i "/DM.UseSaveDM/d" $$TEMPLATE_FDF; echo "DM.UseSaveDM .true." >> $$TEMPLATE_FDF
    sed -i "/TS.Voltage/d" $$TEMPLATE_FDF
    
    echo "" >> $$TEMPLATE_FDF
    echo "TS.Voltage $$V eV" >> $$TEMPLATE_FDF

    echo "Running $$V V siesta..."
    $$SIESTA_COMMAND < $$TEMPLATE_FDF > $$LOG_FILE

    if grep -q "$$SUCCESS_MSG" "$$LOG_FILE"; then
        echo "SUCCESS: $$V V Device calculation completed."
        
        if [ -f siesta.TSHS ]; then
            echo "Running $$V V TBTrans..."
            $$TBTRANS_COMMAND < $$TEMPLATE_FDF > tbtrans_$${V}V.out
            
            if [ -f siesta.TBT.nc ]; then
                mv siesta.TBT.nc device_$${V}V.TBT.nc
                echo "Plotting $${V}V transmission..."
                python3 plot_transmission.py device_$${V}V.TBT.nc "Voltage=$${V}V"
            fi
        fi
    else
        print_step "!!! ERROR: $$V V DEVICE FAILED !!!"
        tail -n 10 $$LOG_FILE
        exit 1
    fi
    cd $$START_DIR
done
print_step "SCRIPT COMPLETE"
echo "All scattering calculations finished successfully."
""")

# HELPER

def ask_question(prompt, default=None):
    val = input(f"{prompt} [{default}]: ").strip()
    return val if val else default

def run_command(command, cwd='.'):
    print(f"\n[RUNNING IN {cwd}]: $ {command}\n")
    try:
        subprocess.run(command, cwd=cwd, shell=True, check=True, executable='/bin/bash')
    except subprocess.CalledProcessError as e:
        print(f"ERROR: Command failed with exit code {e.returncode}")
        sys.exit(1)

def check_and_tag_z_vector(fdf_path):
    print(f"Checking Z-vector tag in {fdf_path}...")
    try:
        with open(fdf_path, 'r') as f: lines = f.readlines()
        tagged = any("# Z-vector" in line for line in lines)
        if tagged: 
            print("  - OK: Z-vector tag already present.")
            return True
        new_lines = []
        in_vectors = False
        count = 0
        for line in lines:
            if "%block LatticeVectors" in line: in_vectors = True
            if "%endblock LatticeVectors" in line: in_vectors = False
            if in_vectors and line.strip() and not line.strip().startswith("%"):
                count += 1
                if count == 3: line = line.strip() + " # Z-vector\n"
            new_lines.append(line)
        with open(fdf_path, 'w') as f: f.writelines(new_lines)
        return True
    except Exception as e:
        print(f"Error checking Z-vector: {e}")
        return False

def parse_energy(filepath):
    if not os.path.exists(filepath): return None
    try:
        with open(filepath, 'r') as f:
            for line in f:
                if "siesta:" in line and "Total =" in line: return float(line.split('=')[1].strip())
    except: return None
    return None

def plot_convergence_graph(x, y, xl, title, fname, cwd):
    try:
        plt.figure(figsize=(10, 6)); plt.plot(x, y, 'o-'); plt.title(title); plt.xlabel(xl)
        plt.ylabel("Total Energy (eV)"); plt.grid(True, ls='--', alpha=0.7)
        plt.savefig(os.path.join(cwd, fname)); plt.close()
    except: pass

def find_converged_value_by_delta(params, engs, th):
    valid = [(p, e) for p, e in zip(params, engs) if e is not None]
    if not valid: return params[-1]
    p, e = zip(*valid)
    for i in range(len(e) - 1):
        if abs(e[i+1] - e[i]) < th: return p[i]
    return p[-1]

def find_converged_value_by_min(params, engs):
    valid = [(p, e) for p, e in zip(params, engs) if e is not None]
    if not valid: return params[-1]
    p, e = zip(*valid)
    return p[e.index(min(e))]

def copy_psf_files(src, dest):
    try:
        if os.path.abspath(src) != os.path.abspath(dest):
            os.makedirs(dest, exist_ok=True)
        for i in os.listdir(src):
            if i.endswith(".psf"): shutil.copy(os.path.join(src, i), dest)
    except: pass

# CORE
def run_convergence_stage(stage_dir, base_fdf_path, siesta_command, plot_title, k_z_val, system_label, conv_threshold, base_dir,
                          mesh_sequence, kpoint_sequence, lattice_sequence,
                          k_dim, k_fixed_x, k_fixed_y, tbtrans_command,
                          initial_k_block=None, do_lattice_convergence=True):
    
    os.makedirs(stage_dir, exist_ok=True)
    copy_psf_files(base_dir, stage_dir)
    create_plot_script(stage_dir)
    
    current_k_block = None
    if initial_k_block: current_k_block = initial_k_block
    else:
        extracted = extract_kpoints_block_from_file(base_fdf_path)
        print(f"\n--- K-Point Selection for {plot_title} ---")
        if extracted:
            print(f"Found:\n{extracted.strip()}")
            if ask_question("Use these? (y/n)", "y").lower() == 'y': current_k_block = extracted
            else: current_k_block = format_manual_k_block(ask_question("Enter manual K-points", "1 1 50"))
        else: current_k_block = format_manual_k_block(ask_question("Enter manual K-points", "1 1 50"))

    clean_name = "clean_template.fdf"
    if not clean_file_remove_kgrid(base_fdf_path, os.path.join(stage_dir, clean_name)): sys.exit(1)

    # 1. Mesh
    print(f"\n--- Stage: Mesh Convergence ({plot_title}) ---")
    ms = os.path.join(stage_dir, "run_mesh_conv.sh")
    with open(ms, 'w') as f: f.write(dedent(MESH_SCRIPT_TEMPLATE.substitute(
        clean_fdf_name=clean_name, mesh_sequence=mesh_sequence, siesta_command=siesta_command,
        tbtrans_command=tbtrans_command, fixed_k_block=current_k_block)))
    os.chmod(ms, 0o755); run_command(f"./{os.path.basename(ms)}", cwd=stage_dir)
    
    m_vals = [int(m) for m in mesh_sequence.split()]
    engs = [parse_energy(os.path.join(stage_dir, f"mesh_{m}Ry", "siesta.out")) for m in m_vals]
    plot_convergence_graph(m_vals, engs, "Mesh (Ry)", f"Mesh {plot_title}", "mesh_plot.png", cwd=stage_dir)
    conv_mesh = find_converged_value_by_delta(m_vals, engs, conv_threshold)
    print(f"--- WINNER: Mesh = {conv_mesh} Ry ---")
    
    # 2. K-Point
    print(f"\n--- Stage: K-Point Convergence ({plot_title}) ---")
    ks = os.path.join(stage_dir, "run_kpoint_conv.sh")
    if k_dim == '1D': k_echo = f'echo -e "%block kgrid_Monkhorst_Pack\\n {k_fixed_x} 0 0 0.0\\n 0 {k_fixed_y} 0 0.0\\n 0 0 $k 0.0\\n%endblock kgrid_Monkhorst_Pack"'; k_dir = "kpoint_z_$k"
    elif k_dim == '3D': k_echo = f'echo -e "%block kgrid_Monkhorst_Pack\\n $k 0 0 0.0\\n 0 $k 0 0.0\\n 0 0 $k 0.0\\n%endblock kgrid_Monkhorst_Pack"'; k_dir = "kpoint_${k}x${k}x${k}"
    else: k_echo = f'echo -e "%block kgrid_Monkhorst_Pack\\n $k 0 0 0.0\\n 0 $k 0 0.0\\n 0 0 {k_z_val} 0.0\\n%endblock kgrid_Monkhorst_Pack"'; k_dir = "kpoint_${k}x${k}"

    with open(ks, 'w') as f: f.write(dedent(KPOINT_SCRIPT_TEMPLATE.substitute(
        clean_fdf_name=clean_name, kpoint_sequence=kpoint_sequence, fixed_mesh_cutoff=conv_mesh,
        siesta_command=siesta_command, tbtrans_command=tbtrans_command,
        kpoint_echo_command=k_echo, kpoint_test_dir_name=k_dir)))
    os.chmod(ks, 0o755); run_command(f"./{os.path.basename(ks)}", cwd=stage_dir)
    
    k_vals = [int(k) for k in kpoint_sequence.split()]
    engs = []
    for k in k_vals:
        if k_dim == '1D': fname = f"kpoint_z_{k}"
        elif k_dim == '3D': fname = f"kpoint_{k}x{k}x{k}"
        else: fname = f"kpoint_{k}x{k}"
        engs.append(parse_energy(os.path.join(stage_dir, fname, "siesta.out")))
    plot_convergence_graph(k_vals, engs, "K-Point", f"K-Point {plot_title}", "kpoint_plot.png", cwd=stage_dir)
    conv_k = find_converged_value_by_delta(k_vals, engs, conv_threshold)
    print(f"--- WINNER: K-Point = {conv_k} ---")

    if k_dim == '1D': k_blk = dedent(f"""%block kgrid_Monkhorst_Pack\n{k_fixed_x} 0 0 0.0\n0 {k_fixed_y} 0 0.0\n0 0 {conv_k} 0.0\n%endblock kgrid_Monkhorst_Pack""").strip()
    elif k_dim == '3D': k_blk = dedent(f"""%block kgrid_Monkhorst_Pack\n{conv_k} 0 0 0.0\n0 {conv_k} 0 0.0\n0 0 {conv_k} 0.0\n%endblock kgrid_Monkhorst_Pack""").strip()
    else: k_blk = dedent(f"""%block kgrid_Monkhorst_Pack\n{conv_k} 0 0 0.0\n0 {conv_k} 0 0.0\n0 0 {k_z_val} 0.0\n%endblock kgrid_Monkhorst_Pack""").strip()

    # 3. Lattice
    if do_lattice_convergence:
        print(f"\n--- Stage: Z-Lattice Convergence ({plot_title}) ---")
        ls = os.path.join(stage_dir, "run_lattice_conv.sh")
        with open(ls, 'w') as f: f.write(dedent(LATTICE_SCRIPT_TEMPLATE.substitute(
            clean_fdf_name=clean_name, lattice_sequence=lattice_sequence, fixed_mesh_cutoff=conv_mesh,
            siesta_command=siesta_command, tbtrans_command=tbtrans_command, fixed_k_block=k_blk)))
        os.chmod(ls, 0o755); run_command(f"./{os.path.basename(ls)}", cwd=stage_dir)
        
        l_vals = [float(l) for l in lattice_sequence.split()]
        engs = [parse_energy(os.path.join(stage_dir, f"lattice_z_{l}A", "siesta.out")) for l in l_vals]
        plot_convergence_graph(l_vals, engs, "Lattice Z", f"Lattice {plot_title}", "lattice_plot.png", cwd=stage_dir)
        conv_lat_z = find_converged_value_by_min(l_vals, engs)
        print(f"--- WINNER: Lattice = {conv_lat_z} A ---")
    else:
        print("--- SKIPPING lattice convergence (Device Mode) ---")
        conv_lat_z = "0.0"
        try:
            with open(base_fdf_path, 'r') as f:
                for line in f:
                    if "# Z-vector" in line: conv_lat_z = line.split()[2]; break
        except: pass

    # 4. Final
    print(f"\n--- Stage: Final Calculation ({plot_title}) ---")
    fs = os.path.join(stage_dir, "run_final_converged.sh")
    with open(fs, 'w') as f: f.write(dedent(FINAL_CALC_SCRIPT_TEMPLATE.substitute(
        clean_fdf_name=clean_name, conv_mesh=conv_mesh, conv_lat_z=conv_lat_z,
        siesta_command=siesta_command, tbtrans_command=tbtrans_command, system_label=system_label,
        fixed_k_block=k_blk, do_lattice_conv="true" if do_lattice_convergence else "false")))
    os.chmod(fs, 0o755); run_command(f"./{os.path.basename(fs)}", cwd=stage_dir)
    
    return {"mesh": conv_mesh, "kpoint": conv_k, "lattice_z": conv_lat_z, "final_k_block": k_blk, "optimize_dir": os.path.join(stage_dir, "optimize")}

def main():
    print("\n" + "="*60 + "\n TranSIESTA Workflow (Sequential) \n" + "="*60)
    base_dir = ask_question("Enter project directory", os.path.abspath("."))
    is_sym = ask_question("Symmetrical electrodes? (y/n)", "y").lower() == 'y'
    
    if is_sym:
        el_fdf = ask_question("Electrode FDF", os.path.join(base_dir, "electrode.fdf"))
        el_L = el_R = el_fdf
    else:
        el_L = ask_question("Left Electrode FDF", os.path.join(base_dir, "electrodel.fdf"))
        el_R = ask_question("Right Electrode FDF", os.path.join(base_dir, "electroder.fdf"))
        
    dev_fdf = ask_question("Device FDF", os.path.join(base_dir, "device.fdf"))
    mpi_procs = ask_question("MPI Processors", "32")
    tbtrans_cmd = ask_question("TBTrans Command", "tbtrans")
    conv_th = float(ask_question("Convergence Threshold (eV)", "0.01"))
    
    mesh_seq = ask_question("Mesh sequence", "150 200")
    k_seq = ask_question("K-point sequence", "1 2 3")
    lat_seq = ask_question("Lattice sequence", "2.80 2.85")
    
    k_dim = ask_question("K-point dimension (1D/2D/3D)", "2D").upper()
    kx, ky, kz_elec, kz_dev = "1", "1", "100", "1"
    if k_dim == '1D': kx = ask_question("Fixed k_x", "1"); ky = ask_question("Fixed k_y", "1")
    elif k_dim == '2D': kz_elec = ask_question("Fixed k_z (Electrode)", "100"); kz_dev = ask_question("Fixed k_z (Device)", "1")
    volt_list = ask_question("Voltage list", "0.01 0.05")

    check_and_tag_z_vector(el_L); 
    if not is_sym: check_and_tag_z_vector(el_R)
    check_and_tag_z_vector(dev_fdf)

    cmd = f"mpirun -np {mpi_procs} siesta --electrode"
    
    if is_sym:
        res = run_convergence_stage(os.path.join(base_dir, "2_electrode_conv"), el_L, cmd, "Electrode", kz_elec, "electrode", conv_th, base_dir, mesh_seq, k_seq, lat_seq, k_dim, kx, ky, tbtrans_cmd, do_lattice_convergence=True)
        res_L = res_R = res
    else:
        res_L = run_convergence_stage(os.path.join(base_dir, "2_elec_L_conv"), el_L, cmd, "Left", kz_elec, "electrodel", conv_th, base_dir, mesh_seq, k_seq, lat_seq, k_dim, kx, ky, tbtrans_cmd, do_lattice_convergence=True)
        res_R = run_convergence_stage(os.path.join(base_dir, "2_elec_R_conv"), el_R, cmd, "Right", kz_elec, "electroder", conv_th, base_dir, mesh_seq, k_seq, lat_seq, k_dim, kx, ky, tbtrans_cmd, do_lattice_convergence=True)

    # DEVICE
    print("\n" + "-"*60 + "\n Stage: Device Convergence \n" + "-"*60)
    dev_dir = os.path.join(base_dir, "3_device_conv")
    os.makedirs(dev_dir, exist_ok=True)
    copy_psf_files(base_dir, dev_dir)
    
    if is_sym:
        for f in os.listdir(res['optimize_dir']):
            if f.startswith("electrode."):
                ext = f.split("electrode")[1]
                shutil.copy(os.path.join(res['optimize_dir'], f), os.path.join(dev_dir, f"electrodel{ext}"))
                shutil.copy(os.path.join(res['optimize_dir'], f), os.path.join(dev_dir, f"electroder{ext}"))
        shutil.copy(el_L, os.path.join(dev_dir, "electrodel.fdf"))
        shutil.copy(el_L, os.path.join(dev_dir, "electroder.fdf"))
    else:
        shutil.copy(el_L, os.path.join(dev_dir, "electrodel.fdf"))
        shutil.copy(el_R, os.path.join(dev_dir, "electroder.fdf"))
        for f in os.listdir(res_L['optimize_dir']):
             if f.startswith("electrodel."): shutil.copy(os.path.join(res_L['optimize_dir'], f), dev_dir)
        for f in os.listdir(res_R['optimize_dir']):
             if f.startswith("electroder."): shutil.copy(os.path.join(res_R['optimize_dir'], f), dev_dir)

    dev_run_fdf = os.path.join(dev_dir, "device_run.fdf")
    create_transiesta_run_file(dev_fdf, dev_run_fdf)

    # Build Device K-block
    k_dev_val = res_L['kpoint']
    if k_dim == '1D': k_blk = dedent(f"""%block kgrid_Monkhorst_Pack\n{kx} 0 0 0.0\n0 {ky} 0 0.0\n0 0 {k_dev_val} 0.0\n%endblock kgrid_Monkhorst_Pack""").strip()
    else: k_blk = dedent(f"""%block kgrid_Monkhorst_Pack\n{k_dev_val} 0 0 0.0\n0 {k_dev_val} 0 0.0\n0 0 {kz_dev} 0.0\n%endblock kgrid_Monkhorst_Pack""").strip()

    dev_cmd = f"mpirun -np {mpi_procs} siesta"
    dev_res = run_convergence_stage(dev_dir, dev_run_fdf, dev_cmd, "Device", kz_dev, "device", conv_th, base_dir, mesh_seq, k_seq, lat_seq, k_dim, kx, ky, tbtrans_cmd, initial_k_block=k_blk, do_lattice_convergence=False)

    # SWEEP
    print("\n" + "-"*60 + "\n Stage: Voltage Sweep \n" + "-"*60)
    sweep_dir = os.path.join(base_dir, "4_voltage_sweep")
    os.makedirs(sweep_dir, exist_ok=True)
    
    # Copy Electrodes and Device files
    shutil.copy(os.path.join(dev_dir, "electrodel.TSHS"), sweep_dir)
    shutil.copy(os.path.join(dev_dir, "electrodel.TSDE"), sweep_dir)
    shutil.copy(os.path.join(dev_dir, "electroder.TSHS"), sweep_dir)
    shutil.copy(os.path.join(dev_dir, "electroder.TSDE"), sweep_dir)
    shutil.copy(os.path.join(dev_dir, "electrodel.fdf"), sweep_dir)
    shutil.copy(os.path.join(dev_dir, "electroder.fdf"), sweep_dir)
    copy_psf_files(base_dir, sweep_dir)
    create_plot_script(sweep_dir)

    # --- FIND OPTIMIZED TSHS/DM AND FORCE NAME 'siesta.*' ---
    opt_dir = dev_res['optimize_dir']
    
    if os.path.exists(os.path.join(opt_dir, "siesta.TSHS")):
        shutil.copy(os.path.join(opt_dir, "siesta.TSHS"), os.path.join(sweep_dir, "siesta.TSHS"))
    elif os.path.exists(os.path.join(opt_dir, "device.TSHS")):
        shutil.copy(os.path.join(opt_dir, "device.TSHS"), os.path.join(sweep_dir, "siesta.TSHS"))
    else:
        print("!!! WARNING: Optimized TSHS not found in Optimize directory !!!")

    if os.path.exists(os.path.join(opt_dir, "siesta.DM")):
        shutil.copy(os.path.join(opt_dir, "siesta.DM"), os.path.join(sweep_dir, "siesta.DM"))
    elif os.path.exists(os.path.join(opt_dir, "device.DM")):
        shutil.copy(os.path.join(opt_dir, "device.DM"), os.path.join(sweep_dir, "siesta.DM"))

    # USE OPTIMIZED FDF AND TSHS
    opt_fdf_src = os.path.join(dev_res['optimize_dir'], "siesta.fdf")
    opt_fdf_dst = os.path.join(sweep_dir, "device_optimized.fdf")
    
    if os.path.exists(opt_fdf_src):
        print(f"--- Found Optimized FDF: {opt_fdf_src} ---")
        clean_file_remove_kgrid(opt_fdf_src, opt_fdf_dst)
    else:
        print("!!! WARNING: Optimized FDF not found, reverting to raw input. !!!")
        clean_file_remove_kgrid(dev_run_fdf, opt_fdf_dst)

    # Pass Converged Values to Bash Script
    v_script = os.path.join(sweep_dir, "run_sweep.sh")
    with open(v_script, 'w') as f: f.write(dedent(VOLTAGE_SWEEP_SCRIPT_TEMPLATE.substitute(
        mpi_procs=mpi_procs, 
        voltage_list=volt_list, 
        template_fdf_name="device_optimized.fdf", 
        tbtrans_command=tbtrans_cmd,
        opt_mesh=dev_res['mesh'],            
        opt_k_block=dev_res['final_k_block'] 
    )))
    os.chmod(v_script, 0o755)
    run_command(f"./{os.path.basename(v_script)}", cwd=sweep_dir)

if __name__ == "__main__":
    try: main()
    except KeyboardInterrupt: sys.exit(1)