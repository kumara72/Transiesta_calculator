#!/usr/bin/env python3

import os
import sys
import subprocess
import shutil
import re
from textwrap import dedent
from string import Template

# Libraries
try:
    import matplotlib
    matplotlib.use('Agg') 
    import matplotlib.pyplot as plt
    import numpy as np
except ImportError:
    print("Error: Matplotlib and Numpy are required.")
    sys.exit(1)

# plotting
# This script will be injected into every folder to plot the .TBT.nc file
PLOT_SCRIPT_CONTENT = """
import sisl
import matplotlib.pyplot as plt
import sys
import os

if len(sys.argv) < 2:
    print("Usage: python plot_transmission.py <TBT_FILE> [LABEL]")
    sys.exit(1)

tbt_file = sys.argv[1]
label_text = sys.argv[2] if len(sys.argv) > 2 else tbt_file

if not os.path.exists(tbt_file):
    print(f"Skipping plot: {tbt_file} not found (Calculation might have failed).")
    sys.exit(0)

try:
    tbt = sisl.get_sile(tbt_file)
    E = tbt.E
    T = tbt.transmission()

    plt.figure(figsize=(8, 6))
    plt.plot(E, T, label=label_text, color='b', linewidth=1.5)
    plt.axvline(x=0, color='k', linestyle='--', linewidth=0.8) # Fermi level
    
    plt.title(f"Transmission: {label_text}")
    plt.xlabel("Energy (eV-Ef)")
    plt.ylabel("Transmission")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.ylim(bottom=0)

    out_name = f"transmission_{label_text.replace(' ', '_').replace('=', '').replace('.', 'p')}.png"
    plt.savefig(out_name, dpi=100)
    print(f"   -> Plot saved: {out_name}")

except Exception as e:
    print(f"   -> Plotting failed for {tbt_file}: {e}")
"""

def create_plot_script(directory):
    """Writes the python plotting script to the target directory."""
    path = os.path.join(directory, "plot_transmission.py")
    with open(path, "w") as f:
        f.write(PLOT_SCRIPT_CONTENT)
    return "plot_transmission.py"

# file cleaner
def clean_file_remove_kgrid(input_path, output_path):
    try:
        with open(input_path, 'r') as f_in:
            lines = f_in.readlines()
        
        with open(output_path, 'w') as f_out:
            skipping = False
            for line in lines:
                check_line = line.lower().strip()
                if check_line.startswith("%block kgrid"): 
                    skipping = True
                    continue 
                if skipping and check_line.startswith("%endblock kgrid"):
                    skipping = False
                    continue
                if not skipping:
                    f_out.write(line)
            f_out.write("\n")
        return True
    except Exception as e:
        print(f"ERROR cleaning file: {e}")
        sys.exit(1)

def extract_kpoints_block_from_file(fdf_path):
    try:
        with open(fdf_path, 'r') as f:
            lines = f.readlines()
        block_content = []
        in_block = False
        found = False
        for line in lines:
            check_line = line.lower().strip()
            if check_line.startswith("%block kgrid"):
                in_block = True
                found = True
                block_content.append("%block kgrid_Monkhorst_Pack")
                continue
            if in_block and check_line.startswith("%endblock kgrid"):
                in_block = False
                block_content.append("%endblock kgrid_Monkhorst_Pack")
                break
            if in_block:
                block_content.append(line.strip())
        if found: return "\n".join(block_content)
        else: return None
    except Exception: return None

def format_manual_k_block(k_str):
    try:
        k = [int(x) for x in k_str.split()]
        if len(k) != 3: raise ValueError
        return dedent(f"""
            %block kgrid_Monkhorst_Pack
            {k[0]}   0   0   0.0
             0  {k[1]}   0   0.0
             0   0 {k[2]}  0.0
            %endblock kgrid_Monkhorst_Pack
            """).strip()
    except:
        return dedent("""
            %block kgrid_Monkhorst_Pack
            1   0   0   0.0
             0  1   0   0.0
             0   0  1  0.0
            %endblock kgrid_Monkhorst_Pack
            """).strip()

# Transiesta block
def create_transiesta_run_file(device_geom_path, output_path):
    print(f"--- Creating device_run.fdf using geometry from {os.path.basename(device_geom_path)} ---")
    try:
        with open(device_geom_path, 'r') as f:
            geom_content = f.read()
        
        ts_block = """
# ==========================================
# --- Advanced TranSIESTA Settings ---
# ==========================================

SolutionMethod          transiesta
TS.Voltage              0.000 eV
TS.SaveHS               .true.

# --- Mixing & Convergence ---
SCF.Mixer.Weight        0.05
SCF.Mixer.History       12
SCF.DM.Tolerance        0.0001
MaxSCFIterations        5000

# --- Complex Contour Integration ---
%block TS.ChemPots
  Left
  Right
%endblock TS.ChemPots

%block TS.ChemPot.Left
  mu V/2
  contour.eq
    begin
      C-Left
      T-Left
    end
%endblock TS.ChemPot.Left

%block TS.ChemPot.Right
  mu -V/2
  contour.eq
    begin
      C-Right
      T-Right
    end
%endblock TS.ChemPot.Right

TS.Contours.Eq.Pole 2.5 eV

%block TS.Contour.C-Left
  part circle
   from -40. eV + V/2 to -10 kT + V/2
     points 50
      method g-legendre
%endblock TS.Contour.C-Left

%block TS.Contour.T-Left
  part tail
   from prev to inf
     points 10
      method g-fermi
%endblock TS.Contour.T-Left

%block TS.Contour.C-Right
  part circle
   from -40. eV -V/2 to -10 kT -V/2
     points 50
      method g-legendre
%endblock TS.Contour.C-Right

%block TS.Contour.T-Right
  part tail
   from prev to inf
     points 10
      method g-fermi
%endblock TS.Contour.T-Right

%block TS.Contours.nEq
  neq-1
%endblock TS.Contours.nEq

%block TS.Contour.nEq.neq-1
  part line
   from -|V|/2 - 5 kT to |V|/2 + 5 kT
     delta 0.005 eV
      method mid-rule
%endblock TS.Contour.nEq.neq-1

# --- Electrode Definitions ---
%block TS.Elecs
  Left
  Right
%endblock TS.Elecs

%block TS.Elec.Left
  HS electrodel.TSHS
  semi-inf-direction -a3
  electrode-position 1
%endblock TS.Elec.Left

%block TS.Elec.Right
  HS electroder.TSHS
  semi-inf-direction +a3
  electrode-position end -1
%endblock TS.Elec.Right

# --- TBTrans Options ---
%block TBT.k
  diag 1 1 1
%endblock TBT.k
TBT.DOS.A T

%block TBT.Contours
    line
%endblock TBT.Contours

%block TBT.Contour.line
  part line
     from -1. eV to 1. eV
      delta 0.01 eV
        method mid-rule
%endblock TBT.Contour.line
"""
        with open(output_path, 'w') as f:
            f.write(geom_content)
            f.write(ts_block)
            
    except Exception as e:
        print(f"Error creating TranSIESTA file: {e}")
        sys.exit(1)

# Script Templates

MESH_SCRIPT_TEMPLATE = Template("""
#!/bin/bash
echo "--- Starting: Mesh Cutoff Convergence Test ---"
module load siesta || echo "Warning: 'module load siesta' failed."

CLEAN_FDF="$clean_fdf_name"
MESH_SEQUENCE="$mesh_sequence"
SIESTA_COMMAND="$siesta_command"

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
    
    # Copy plotting script
    cp ../plot_transmission.py .
    
    cp ../$$CLEAN_FDF siesta.fdf
    sed -i "s/^[[:space:]]*SystemLabel.*/SystemLabel siesta/" siesta.fdf
    sed -i "s/^[[:space:]]*MeshCutoff.*/MeshCutoff $$mesh Ry/" siesta.fdf
    
    # Force MaxSCFIterations
    if grep -q "MaxSCFIterations" siesta.fdf; then
        sed -i "s/^[[:space:]]*MaxSCFIterations.*/MaxSCFIterations 5000/" siesta.fdf
    else
        echo "MaxSCFIterations 5000" >> siesta.fdf
    fi
    
    # Force Safe Mixing
    sed -i "/SCF.Mixer.Weight/d" siesta.fdf
    sed -i "/SCF.Mixer.History/d" siesta.fdf
    
cat <<EOF >> siesta.fdf

SCF.Mixer.Weight    0.01
SCF.Mixer.History   20
$fixed_k_block
DM.UseSaveDM .true.
EOF
    
    echo "--- Running calculation in $$DIR_NAME (Mesh: $$mesh) ---"
    $$SIESTA_COMMAND < siesta.fdf > siesta.out
    
    # PLOT TRANSMISSION
    if [ -f siesta.TBT.nc ]; then
        python3 plot_transmission.py siesta.TBT.nc "Mesh=$${mesh}Ry"
    fi
    
    cd ..
    rm -rf cont && mkdir cont
    if [ -f ./$$DIR_NAME/siesta.DM ]; then
        cp ./$$DIR_NAME/siesta.DM cont/
    else
        echo "--- ERROR: siesta.DM was not created in $$DIR_NAME. ---"
        exit 1
    fi
done
""")

KPOINT_SCRIPT_TEMPLATE = Template("""
#!/bin/bash
echo "--- Starting: K-Point Convergence Test ---"
module load siesta || echo "Warning: 'module load siesta' failed."

CLEAN_FDF="$clean_fdf_name"
KPOINT_SEQUENCE="$kpoint_sequence"
FIXED_MESH_CUTOFF="$fixed_mesh_cutoff"
SIESTA_COMMAND="$siesta_command"

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
    
    # Copy plotting script
    cp ../plot_transmission.py .
    
    cp ../$$CLEAN_FDF siesta.fdf
    sed -i "s/^[[:space:]]*SystemLabel.*/SystemLabel siesta/" siesta.fdf
    sed -i "s/^[[:space:]]*MeshCutoff.*/MeshCutoff $$FIXED_MESH_CUTOFF Ry/" siesta.fdf

    if grep -q "MaxSCFIterations" siesta.fdf; then
        sed -i "s/^[[:space:]]*MaxSCFIterations.*/MaxSCFIterations 5000/" siesta.fdf
    else
        echo "MaxSCFIterations 5000" >> siesta.fdf
    fi
    
    sed -i "/SCF.Mixer.Weight/d" siesta.fdf
    sed -i "/SCF.Mixer.History/d" siesta.fdf
    echo "SCF.Mixer.Weight 0.01" >> siesta.fdf
    echo "SCF.Mixer.History 20" >> siesta.fdf
    
    echo "" >> siesta.fdf
    $kpoint_echo_command >> siesta.fdf
    echo "" >> siesta.fdf
    echo "DM.UseSaveDM .true." >> siesta.fdf
    
    echo "--- Running calculation in $$DIR_NAME (K: Varying) ---"
    $$SIESTA_COMMAND < siesta.fdf > siesta.out
    
    # PLOT TRANSMISSION
    if [ -f siesta.TBT.nc ]; then
        python3 plot_transmission.py siesta.TBT.nc "Kpoint=$$k"
    fi

    cd ..
    rm -rf cont_kpoint && mkdir cont_kpoint
    if [ -f ./$$DIR_NAME/siesta.DM ]; then
        cp ./$$DIR_NAME/siesta.DM cont_kpoint/
    else
        echo "--- ERROR: siesta.DM was not created in $$DIR_NAME. ---"
        exit 1
    fi
done
""")

LATTICE_SCRIPT_TEMPLATE = Template("""
#!/bin/bash
echo "--- Starting: Z-Lattice Convergence Test ---"
module load siesta || echo "Warning: 'module load siesta' failed."

CLEAN_FDF="$clean_fdf_name"
LATTICE_SEQUENCE="$lattice_sequence"
FIXED_MESH_CUTOFF="$fixed_mesh_cutoff"
SIESTA_COMMAND="$siesta_command"
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
    
    cp ../$$CLEAN_FDF siesta.fdf
    
    sed -i "s/^[[:space:]]*SystemLabel.*/SystemLabel siesta/" siesta.fdf
    sed -i "s/^[[:space:]]*MeshCutoff.*/MeshCutoff $$FIXED_MESH_CUTOFF Ry/" siesta.fdf
    
cat <<EOF >> siesta.fdf

$fixed_k_block

EOF
    
    sed -i "s/^[[:space:]]*0.00000000[[:space:]]*0.00000000[[:space:]].*# Z-vector/  0.00000000  0.00000000  $$lat_z # Z-vector/" siesta.fdf

    echo "DM.UseSaveDM .true." >> siesta.fdf
    echo "--- Running calculation in $$DIR_NAME (Lat: $$lat_z) ---"
    $$SIESTA_COMMAND < siesta.fdf > siesta.out
    cd ..
    rm -rf cont_lattice && mkdir cont_lattice
    if [ -f ./$$DIR_NAME/siesta.DM ]; then
        cp ./$$DIR_NAME/siesta.DM cont_lattice/
    else
        echo "--- ERROR: siesta.DM was not created in $$DIR_NAME. ---"
        exit 1
    fi
done
""")

FINAL_CALC_SCRIPT_TEMPLATE = Template("""
#!/bin/bash
echo "--- Starting: Final Converged Calculation ---"
module load siesta || echo "Warning: 'module load siesta' failed."

CLEAN_FDF="$clean_fdf_name"
CONV_MESH="$conv_mesh"
CONV_LAT_Z="$conv_lat_z"
SIESTA_COMMAND="$siesta_command"
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

# Copy plotting script
cp ../plot_transmission.py .

cp ../$$CLEAN_FDF siesta.fdf
sed -i "s/^[[:space:]]*MeshCutoff.*/MeshCutoff $$CONV_MESH Ry/" siesta.fdf
sed -i "s/^[[:space:]]*SystemLabel.*/SystemLabel siesta/" siesta.fdf

if [ "$$DO_LATTICE_CONV" = "true" ]; then
    echo "--- Applying converged lattice constant: $$CONV_LAT_Z ---"
    sed -i "s/^[[:space:]]*0.00000000[[:space:]]*0.00000000[[:space:]].*# Z-vector/  0.00000000  0.00000000  $$CONV_LAT_Z # Z-vector/" siesta.fdf
fi

cat <<EOF >> siesta.fdf

$fixed_k_block

EOF

echo "--- Running final calculation in $$DIR_NAME ---"
$$SIESTA_COMMAND < siesta.fdf > siesta.out

# PLOT TRANSMISSION
if [ -f siesta.TBT.nc ]; then
    python3 plot_transmission.py siesta.TBT.nc "Device_0.0V"
fi

if [ -f siesta.TSHS ] || [ -f siesta.DM ]; then
    echo "--- Renaming files from siesta.* to $$SYSTEM_LABEL.* ---"
    for f in siesta.*; do
        mv "$$f" "$${SYSTEM_LABEL}.$${f#siesta.}"
    done
else
    echo "--- ERROR: Final calculation failed. Check $$DIR_NAME/siesta.out ---"
    exit 1
fi
""")

VOLTAGE_SWEEP_SCRIPT_TEMPLATE = Template("""
#!/bin/bash
echo "--- Starting: Voltage Sweep Workflow ---"
module load siesta || echo "Warning: 'module load siesta' failed. Assuming siesta is in PATH."

START_DIR=$$(pwd)
SUCCESS_MSG="Job completed"
SIESTA_COMMAND="mpirun -np $mpi_procs siesta"

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

# Copy template and files
cp ../$$TEMPLATE_FDF .
cp ../$$PSF_FILES . 2>/dev/null
cp ../electrodel.TSHS .
cp ../electrodel.TSDE .
cp ../electroder.TSHS .
cp ../electroder.TSDE .
cp ../electrodel.fdf . 2>/dev/null
cp ../electroder.fdf . 2>/dev/null
cp ../plot_transmission.py .

echo "Running 0V siesta command..."
$$SIESTA_COMMAND < $$TEMPLATE_FDF > $$LOG_FILE

if grep -q "$$SUCCESS_MSG" "$$LOG_FILE"; then
    echo "SUCCESS: 0V Device calculation completed."
    
    # RENAME TBT FILE
    if [ -f siesta.TBT.nc ]; then
        mv siesta.TBT.nc device_0.0V.TBT.nc
        echo "Plotting 0.0V transmission..."
        python3 plot_transmission.py device_0.0V.TBT.nc "Voltage=0.0V"
    fi
else
    print_step "!!! ERROR: 0V DEVICE FAILED !!!"
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
    
    print_step "2. Running $$V V Device Calculation (Chained from $$PREV_V V)"
    
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
    
    cp $$START_DIR/$$PREV_VOLT_DIR/device.DM .
    cp $$START_DIR/volt_0.0V/device.TSHS .

    sed -i.bak "s/^[[:space:]]*TS.Voltage.*/TS.Voltage $$V V/" $$TEMPLATE_FDF
    if ! grep -q "DM.UseSaveDM" $$TEMPLATE_FDF; then
        echo "DM.UseSaveDM .true." >> $$TEMPLATE_FDF
    else
        sed -i.bak "s/^[[:space:]]*DM.UseSaveDM.*/DM.UseSaveDM .true./" $$TEMPLATE_FDF
    fi

    $$SIESTA_COMMAND < $$TEMPLATE_FDF > $$LOG_FILE

    if grep -q "$$SUCCESS_MSG" "$$LOG_FILE"; then
        echo "SUCCESS: $$V V Device calculation completed."
        
        # RENAME TBT FILE
        if [ -f siesta.TBT.nc ]; then
            mv siesta.TBT.nc device_${V}V.TBT.nc
            echo "Plotting ${V}V transmission..."
            python3 plot_transmission.py device_${V}V.TBT.nc "Voltage=${V}V"
        fi
    else
        print_step "!!! ERROR: $$V V DEVICE FAILED !!!"
        exit 1
    fi
    cd $$START_DIR
done
print_step "SCRIPT COMPLETE"
echo "All scattering calculations finished successfully."
""")

# Helper Functions

def print_header(message):
    print("\n" + "=" * 70)
    print(f" {message} ".center(70))
    print("=" * 70)

def print_subheader(message):
    print("\n" + "-" * 70)
    print(f" {message} ")
    print("-" * 70)

def ask_question(prompt, default=None):
    val = input(f"{prompt} [{default}]: ").strip()
    return val if val else default

def run_command(command, cwd='.'):
    print(f"\n[RUNNING IN {cwd}]: $ {command}\n")
    try:
        subprocess.run(command, cwd=cwd, shell=True, check=True, executable='/bin/bash')
    except subprocess.CalledProcessError as e:
        print_header(f"ERROR: Command failed with exit code {e.returncode}")
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
                if count == 3:
                    line = line.strip() + " # Z-vector\n"
            new_lines.append(line)
            
        with open(fdf_path, 'w') as f: f.writelines(new_lines)
        print("  - Tagged Z-vector.")
        return True
    except Exception as e:
        print(f"Error checking Z-vector: {e}")
        return False

def parse_energy(filepath):
    if not os.path.exists(filepath): return None
    try:
        with open(filepath, 'r') as f:
            for line in f:
                if "siesta:" in line and "Total =" in line:
                    return float(line.split('=')[1].strip())
    except: return None
    return None

def plot_convergence_graph(x_values, y_values, x_label, title, plot_filename, cwd='.'):
    plot_data = [(x, y) for x, y in zip(x_values, y_values) if y is not None]
    if not plot_data:
        return
    plot_x, plot_y = zip(*plot_data)
    try:
        plt.figure(figsize=(10, 6))
        plt.plot(x_values, y_values, marker='o', linestyle='-')
        plt.title(title)
        plt.xlabel(x_label)
        plt.ylabel("Total Energy (eV)")
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.savefig(os.path.join(cwd, plot_filename))
        plt.close()
    except: pass

def find_converged_value_by_delta(parameters, energies, threshold, param_name="Parameter"):
    valid = [(p, e) for p, e in zip(parameters, energies) if e is not None]
    if not valid: return parameters[-1]
    params, engs = zip(*valid)
    for i in range(len(engs) - 1):
        if abs(engs[i+1] - engs[i]) < threshold: return params[i]
    return params[-1]

def find_converged_value_by_min(parameters, energies, param_name="Parameter"):
    valid = [(p, e) for p, e in zip(parameters, energies) if e is not None]
    if not valid: return parameters[-1]
    params, engs = zip(*valid)
    return params[engs.index(min(engs))]

def copy_psf_files(src_dir, dest_dir):
    try:
        if os.path.abspath(src_dir) != os.path.abspath(dest_dir):
            os.makedirs(dest_dir, exist_ok=True)
        for item in os.listdir(src_dir):
            if item.endswith(".psf"): shutil.copy(os.path.join(src_dir, item), dest_dir)
    except: pass

# Just Copy User File 
def create_transiesta_run_file(device_geom_path, output_path):
    print(f"--- Using user provided device file: {device_geom_path} ---")
    try:
        shutil.copy(device_geom_path, output_path)
    except Exception as e:
        print(f"Error copying device file: {e}")
        sys.exit(1)

# core function 
def run_convergence_stage(stage_dir, base_fdf_path, siesta_command, plot_title, k_z_val, system_label, conv_threshold, base_dir,
                          mesh_sequence, kpoint_sequence, lattice_sequence,
                          k_dim, k_fixed_x, k_fixed_y,
                          initial_k_block=None, 
                          do_lattice_convergence=True):
    
    os.makedirs(stage_dir, exist_ok=True)
    copy_psf_files(base_dir, stage_dir)
    
    # WRITE THE PLOTTING SCRIPT
    create_plot_script(stage_dir)
    
    current_k_block = None
    if initial_k_block:
        current_k_block = initial_k_block
        print(f"--- Using inherited K-points from previous stage ---")
    else:
        extracted_block = extract_kpoints_block_from_file(base_fdf_path)
        print(f"\n--- K-Point Selection for {plot_title} Mesh Test ---")
        if extracted_block:
            print(f"Found in file:\n{extracted_block.strip()}")
            choice = ask_question("Use these K-points for the Mesh Test? (y/n)", "y") 
            if choice.lower() == 'y':
                current_k_block = extracted_block
            else:
                default_k = "1 1 50"
                user_k = ask_question(f"Enter manual K-points (e.g., '{default_k}')", default_k)
                current_k_block = format_manual_k_block(user_k)
        else:
            default_k = "1 1 50"
            user_k = ask_question(f"Enter manual K-points (e.g., '{default_k}')", default_k)
            current_k_block = format_manual_k_block(user_k)

    # Pre-clean
    clean_template_name = "clean_template.fdf"
    clean_template_path = os.path.join(stage_dir, clean_template_name)
    if not clean_file_remove_kgrid(base_fdf_path, clean_template_path):
        sys.exit(1)

    # 1. Mesh
    print_subheader(f"Stage: Mesh Convergence ({plot_title})")
    mesh_script_path = os.path.join(stage_dir, "run_mesh_conv.sh")
    with open(mesh_script_path, 'w') as f: f.write(dedent(MESH_SCRIPT_TEMPLATE.substitute(
        clean_fdf_name=clean_template_name, mesh_sequence=mesh_sequence,
        siesta_command=siesta_command, fixed_k_block=current_k_block)))
    os.chmod(mesh_script_path, 0o755)
    run_command(f"./{os.path.basename(mesh_script_path)}", cwd=stage_dir)
    
    mesh_vals = [int(m) for m in mesh_sequence.split()]
    energies = [parse_energy(os.path.join(stage_dir, f"mesh_{m}Ry", "siesta.out")) for m in mesh_vals]
    plot_convergence_graph(mesh_vals, energies, "Mesh (Ry)", f"Mesh {plot_title}", "mesh_plot.png", cwd=stage_dir)
    conv_mesh = find_converged_value_by_delta(mesh_vals, energies, conv_threshold)
    print(f"--- WINNER: Mesh = {conv_mesh} Ry ---")
    
    # 2. K-Point
    print_subheader(f"Stage: K-Point Convergence ({plot_title})")
    kpoint_script_path = os.path.join(stage_dir, "run_kpoint_conv.sh")
    
    if k_dim == '1D':
        k_echo_cmd = f'echo -e "%block kgrid_Monkhorst_Pack\\n {k_fixed_x} 0 0 0.0\\n 0 {k_fixed_y} 0 0.0\\n 0 0 $k 0.0\\n%endblock kgrid_Monkhorst_Pack"'
        kpoint_dir_name = "kpoint_z_$k"
    elif k_dim == '3D':
        k_echo_cmd = f'echo -e "%block kgrid_Monkhorst_Pack\\n $k 0 0 0.0\\n 0 $k 0 0.0\\n 0 0 $k 0.0\\n%endblock kgrid_Monkhorst_Pack"'
        kpoint_dir_name = "kpoint_${k}x${k}x${k}"
    else:
        k_echo_cmd = f'echo -e "%block kgrid_Monkhorst_Pack\\n $k 0 0 0.0\\n 0 $k 0 0.0\\n 0 0 {k_z_val} 0.0\\n%endblock kgrid_Monkhorst_Pack"'
        kpoint_dir_name = "kpoint_${k}x${k}"

    with open(kpoint_script_path, 'w') as f: f.write(dedent(KPOINT_SCRIPT_TEMPLATE.substitute(
        clean_fdf_name=clean_template_name, kpoint_sequence=kpoint_sequence,
        fixed_mesh_cutoff=conv_mesh, siesta_command=siesta_command,
        kpoint_echo_command=k_echo_cmd, kpoint_test_dir_name=kpoint_dir_name)))
    os.chmod(kpoint_script_path, 0o755)
    run_command(f"./{os.path.basename(kpoint_script_path)}", cwd=stage_dir)
    
    k_vals = [int(k) for k in kpoint_sequence.split()]
    # Logic to construct folder name for parsing
    energies = []
    for k in k_vals:
        if k_dim == '1D': fname = f"kpoint_z_{k}"
        elif k_dim == '3D': fname = f"kpoint_{k}x{k}x{k}"
        else: fname = f"kpoint_{k}x{k}"
        energies.append(parse_energy(os.path.join(stage_dir, fname, "siesta.out")))
        
    plot_convergence_graph(k_vals, energies, "K-Point", f"K-Point {plot_title}", "kpoint_plot.png", cwd=stage_dir)
    conv_k = find_converged_value_by_delta(k_vals, energies, conv_threshold)
    print(f"--- WINNER: K-Point = {conv_k} ---")

    # block for final selection
    if k_dim == '1D':
        conv_k_block = dedent(f"""
            %block kgrid_Monkhorst_Pack
            {k_fixed_x} 0 0 0.0
            0 {k_fixed_y} 0 0.0
            0 0 {conv_k} 0.0
            %endblock kgrid_Monkhorst_Pack""").strip()
    elif k_dim == '3D':
        conv_k_block = dedent(f"""
            %block kgrid_Monkhorst_Pack
            {conv_k} 0 0 0.0
            0 {conv_k} 0 0.0
            0 0 {conv_k} 0.0
            %endblock kgrid_Monkhorst_Pack""").strip()
    else:
        conv_k_block = dedent(f"""
            %block kgrid_Monkhorst_Pack
            {conv_k} 0 0 0.0
            0 {conv_k} 0 0.0
            0 0 {k_z_val} 0.0
            %endblock kgrid_Monkhorst_Pack""").strip()

    # 3. Lattice
    if do_lattice_convergence:
        print_subheader(f"Stage: Z-Lattice Convergence ({plot_title})")
        lat_script_path = os.path.join(stage_dir, "run_lattice_conv.sh")
        with open(lat_script_path, 'w') as f: f.write(dedent(LATTICE_SCRIPT_TEMPLATE.substitute(
            clean_fdf_name=clean_template_name, lattice_sequence=lattice_sequence,
            fixed_mesh_cutoff=conv_mesh, siesta_command=siesta_command,
            fixed_k_block=conv_k_block)))
        os.chmod(lat_script_path, 0o755)
        run_command(f"./{os.path.basename(lat_script_path)}", cwd=stage_dir)
        
        lat_vals = [float(l) for l in lattice_sequence.split()]
        energies = [parse_energy(os.path.join(stage_dir, f"lattice_z_{l}A", "siesta.out")) for l in lat_vals]
        plot_convergence_graph(lat_vals, energies, "Lattice Z", f"Lattice {plot_title}", "lattice_plot.png", cwd=stage_dir)
        conv_lat_z = find_converged_value_by_min(lat_vals, energies)
        print(f"--- WINNER: Lattice = {conv_lat_z} A ---")
    else:
        print("--- SKIPPING lattice convergence (Device Mode) ---")
        conv_lat_z = "0.0"
        # Try to read existing lattice from file just for logging
        try:
            with open(base_fdf_path, 'r') as f:
                for line in f:
                    if "# Z-vector" in line:
                        conv_lat_z = line.split()[2]
                        break
        except: pass

    # 4. Final
    print_subheader(f"Stage: Final Calculation ({plot_title})")
    final_script_path = os.path.join(stage_dir, "run_final_converged.sh")
    with open(final_script_path, 'w') as f: f.write(dedent(FINAL_CALC_SCRIPT_TEMPLATE.substitute(
        clean_fdf_name=clean_template_name, conv_mesh=conv_mesh, conv_lat_z=conv_lat_z,
        siesta_command=siesta_command, system_label=system_label,
        fixed_k_block=conv_k_block, do_lattice_conv="true" if do_lattice_convergence else "false")))
    os.chmod(final_script_path, 0o755)
    run_command(f"./{os.path.basename(final_script_path)}", cwd=stage_dir)
    
    return {
        "mesh": conv_mesh, "kpoint": conv_k, "lattice_z": conv_lat_z,
        "final_k_block": conv_k_block, "optimize_dir": os.path.join(stage_dir, "optimize")
    }

def main():
    print_header("TranSIESTA Workflow (Sequential)")
    
    base_dir = ask_question("Enter project directory", os.path.abspath("."))
    is_symmetrical = ask_question("Symmetrical electrodes? (y/n)", "y").lower() == 'y'
    
    if is_symmetrical:
        elec_fdf = ask_question("Electrode FDF", os.path.join(base_dir, "electrode.fdf"))
        el_L, el_R = elec_fdf, elec_fdf
    else:
        el_L = ask_question("Left Electrode FDF", os.path.join(base_dir, "electrodel.fdf"))
        el_R = ask_question("Right Electrode FDF", os.path.join(base_dir, "electroder.fdf"))
        
    dev_fdf = ask_question("Device FDF", os.path.join(base_dir, "device.fdf"))
    mpi_procs = ask_question("MPI Processors", "32")
    conv_th = float(ask_question("Convergence Threshold (eV)", "0.01"))
    
    mesh_seq = ask_question("Mesh sequence", "150 200")
    k_seq = ask_question("K-point sequence", "1 2 3")
    lat_seq = ask_question("Lattice sequence", "2.80 2.85")
    
    k_dim = ask_question("K-point dimension (1D/2D/3D)", "2D").upper()
    kx, ky, kz_elec, kz_dev = "1", "1", "100", "1"
    if k_dim == '1D':
        kx = ask_question("Fixed k_x", "1")
        ky = ask_question("Fixed k_y", "1")
    elif k_dim == '2D':
        kz_elec = ask_question("Fixed k_z (Electrode)", "100")
        kz_dev = ask_question("Fixed k_z (Device)", "1")
        
    volt_list = ask_question("Voltage list", "0.01 0.05")

    # Prep files
    check_and_tag_z_vector(el_L)
    if not is_symmetrical: check_and_tag_z_vector(el_R)
    check_and_tag_z_vector(dev_fdf)

    cmd = f"mpirun -np {mpi_procs} siesta --electrode"
    
    # ELECTRODE
    if is_symmetrical:
        res = run_convergence_stage(os.path.join(base_dir, "2_electrode_conv"), el_L, cmd, "Electrode", kz_elec, "electrode", conv_th, base_dir, mesh_seq, k_seq, lat_seq, k_dim, kx, ky, do_lattice_convergence=True)
        res_L, res_R = res, res
    else:
        res_L = run_convergence_stage(os.path.join(base_dir, "2_elec_L_conv"), el_L, cmd, "Left", kz_elec, "electrodel", conv_th, base_dir, mesh_seq, k_seq, lat_seq, k_dim, kx, ky, do_lattice_convergence=True)
        res_R = run_convergence_stage(os.path.join(base_dir, "2_elec_R_conv"), el_R, cmd, "Right", kz_elec, "electroder", conv_th, base_dir, mesh_seq, k_seq, lat_seq, k_dim, kx, ky, do_lattice_convergence=True)

    # DEVICE
    dev_dir = os.path.join(base_dir, "3_device_conv")
    os.makedirs(dev_dir, exist_ok=True)
    copy_psf_files(base_dir, dev_dir)
    
    # Copy electrode results
    if is_symmetrical:
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
        # Copy results logic...
        for f in os.listdir(res_L['optimize_dir']):
             if f.startswith("electrodel."): shutil.copy(os.path.join(res_L['optimize_dir'], f), dev_dir)
        for f in os.listdir(res_R['optimize_dir']):
             if f.startswith("electroder."): shutil.copy(os.path.join(res_R['optimize_dir'], f), dev_dir)

    dev_run_fdf = os.path.join(dev_dir, "device_run.fdf")
    create_transiesta_run_file(dev_fdf, dev_run_fdf)

    # Update lattice Z in device file
    try:
        run_command(f"sed -i 's/^[[:space:]]*0.00000000[[:space:]]*0.00000000[[:space:]].*# Z-vector/  0.00000000  0.00000000  {res_L['lattice_z']} # Z-vector/' {dev_run_fdf}", cwd=dev_dir)
    except: pass

    # Build Device K-block
    k_dev_val = res_L['kpoint']
    if k_dim == '1D':
        dev_k_block = dedent(f"""
            %block kgrid_Monkhorst_Pack
            {kx} 0 0 0.0
            0 {ky} 0 0.0
            0 0 {k_dev_val} 0.0
            %endblock kgrid_Monkhorst_Pack""").strip()
    else:
        dev_k_block = dedent(f"""
            %block kgrid_Monkhorst_Pack
            {k_dev_val} 0 0 0.0
            0 {k_dev_val} 0 0.0
            0 0 {kz_dev} 0.0
            %endblock kgrid_Monkhorst_Pack""").strip()

    dev_cmd = f"mpirun -np {mpi_procs} siesta"
    dev_res = run_convergence_stage(dev_dir, dev_run_fdf, dev_cmd, "Device", kz_dev, "device", conv_th, base_dir, mesh_seq, k_seq, lat_seq, k_dim, kx, ky, initial_k_block=dev_k_block, do_lattice_convergence=False)

    # SWEEP
    sweep_dir = os.path.join(base_dir, "4_voltage_sweep")
    os.makedirs(sweep_dir, exist_ok=True)
    for f in os.listdir(dev_res['optimize_dir']):
        if f.startswith("device."): shutil.copy(os.path.join(dev_res['optimize_dir'], f), sweep_dir)
    
    shutil.copy(os.path.join(dev_dir, "electrodel.TSHS"), sweep_dir)
    shutil.copy(os.path.join(dev_dir, "electrodel.TSDE"), sweep_dir)
    shutil.copy(os.path.join(dev_dir, "electroder.TSHS"), sweep_dir)
    shutil.copy(os.path.join(dev_dir, "electroder.TSDE"), sweep_dir)
    shutil.copy(os.path.join(dev_dir, "electrodel.fdf"), sweep_dir)
    shutil.copy(os.path.join(dev_dir, "electroder.fdf"), sweep_dir)
    shutil.copy(dev_run_fdf, sweep_dir)
    copy_psf_files(base_dir, sweep_dir)
    create_plot_script(sweep_dir)
    
    v_script = os.path.join(sweep_dir, "run_sweep.sh")
    with open(v_script, 'w') as f: f.write(dedent(VOLTAGE_SWEEP_SCRIPT_TEMPLATE.substitute(
        mpi_procs=mpi_procs, voltage_list=volt_list, template_fdf_name="device_run.fdf")))
    os.chmod(v_script, 0o755)
    run_command(f"./{os.path.basename(v_script)}", cwd=sweep_dir)

if __name__ == "__main__":
    try: main()
    except KeyboardInterrupt: sys.exit(1)
