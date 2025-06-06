import os
import shutil
import h5py
import numpy as np

# Part O. Locate current script
script_dir = os.path.dirname(os.path.abspath(__file__))


# Part I. Read input card coupling.dat for the mesh information
coupling_path = os.path.join(script_dir, 'coupling.dat')
with open(coupling_path, 'r') as coupling_file:
    lines = coupling_file.readlines()
    for i, line in enumerate(lines):
        if line.strip() == '#DIM':
            dim = lines[i+1].split()
        elif line.strip() == '#OUTDIMF':
            out_F_dim = lines[i+1].split()
        elif line.strip() == '#OUTDIMC':
            out_C_dim = lines[i+1].split()
        elif line.strip() == '#OUTDIMM':
            out_M_dim = lines[i+1].split()
        elif line.strip() == '#OUTDIMR':
            out_R_dim = lines[i+1].split()
        elif line.strip() == '#P_bndry in [cm] {xmin xmax ymin ymax zmin zmax}':
            P_bndry = lines[i+1].split()
        elif line.strip() == '#F_bndry in [m] {xmin xmax ymin ymax zmin zmax}':
            F_bndry = lines[i+1].split()
        elif line.strip() == '#C_bndry in [m] {xmin xmax ymin ymax zmin zmax}':
            C_bndry = lines[i+1].split()
        elif line.strip() == '#M_bndry in [m] {xmin xmax ymin ymax zmin zmax}':
            M_bndry = lines[i+1].split()
        elif line.strip() == '#R_bndry in [m] {xmin xmax ymin ymax zmin zmax}':
            R_bndry = lines[i+1].split()
        elif line.startswith('#Multilevel Flag'):
            flag = int(lines[i].split()[-1])
        elif line.startswith('#Coupling mode'):
            mode = int(lines[i+1])


# Part II. Update macro definitions in udf.c, and modify paths
udf_path = os.path.join(script_dir, 'libudf', 'src', 'udf.c')  # udf.c path
temp_file_path = udf_path + '.tmp'  # Temporal file for modifying

with open(udf_path, 'r', encoding='utf-8') as udf_file, open(temp_file_path, 'w', encoding='utf-8') as temp_file:
    for line in udf_file:

        # Mesh number definition

        if line.startswith('#define DIM0'):
            temp_file.write(f'#define DIM0 {dim[0]}\n')
        elif line.startswith('#define DIM1'):
            temp_file.write(f'#define DIM1 {dim[1]}\n')
        elif line.startswith('#define DIM2'):
            temp_file.write(f'#define DIM2 {dim[2]}\n')
        elif line.startswith('#define OUTDIMF0'):
            temp_file.write(f'#define OUTDIMF0 {out_F_dim[0]}\n')
        elif line.startswith('#define OUTDIMF1'):
            temp_file.write(f'#define OUTDIMF1 {out_F_dim[1]}\n')
        elif line.startswith('#define OUTDIMF2'):
            temp_file.write(f'#define OUTDIMF2 {out_F_dim[2]}\n')
        elif line.startswith('#define OUTDIMC0'):
            temp_file.write(f'#define OUTDIMC0 {out_C_dim[0]}\n')
        elif line.startswith('#define OUTDIMC1'):
            temp_file.write(f'#define OUTDIMC1 {out_C_dim[1]}\n')
        elif line.startswith('#define OUTDIMC2'):
            temp_file.write(f'#define OUTDIMC2 {out_C_dim[2]}\n')
        elif line.startswith('#define OUTDIMM0'):
            temp_file.write(f'#define OUTDIMM0 {out_M_dim[0]}\n')
        elif line.startswith('#define OUTDIMM1'):
            temp_file.write(f'#define OUTDIMM1 {out_M_dim[1]}\n')
        elif line.startswith('#define OUTDIMM2'):
            temp_file.write(f'#define OUTDIMM2 {out_M_dim[2]}\n')
        elif line.startswith('#define OUTDIMR0'):
            temp_file.write(f'#define OUTDIMR0 {out_R_dim[0]}\n')
        elif line.startswith('#define OUTDIMR1'):
            temp_file.write(f'#define OUTDIMR1 {out_R_dim[1]}\n')
        elif line.startswith('#define OUTDIMR2'):
            temp_file.write(f'#define OUTDIMR2 {out_R_dim[2]}\n')

        # Mesh boundary definition

        elif line.startswith('#define F_xmin'):
            temp_file.write(f'#define F_xmin {float(F_bndry[0])}\n')
        elif line.startswith('#define F_xmax'):
            temp_file.write(f'#define F_xmax {float(F_bndry[1])}\n')
        elif line.startswith('#define F_ymin'):
            temp_file.write(f'#define F_ymin {float(F_bndry[2])}\n')
        elif line.startswith('#define F_ymax'):
            temp_file.write(f'#define F_ymax {float(F_bndry[3])}\n')
        elif line.startswith('#define F_zmin'):
            temp_file.write(f'#define F_zmin {float(F_bndry[4])}\n')
        elif line.startswith('#define F_zmax'):
            temp_file.write(f'#define F_zmax {float(F_bndry[5])}\n')

        elif line.startswith('#define C_xmin'):
            temp_file.write(f'#define C_xmin {float(C_bndry[0])}\n')
        elif line.startswith('#define C_xmax'):
            temp_file.write(f'#define C_xmax {float(C_bndry[1])}\n')
        elif line.startswith('#define C_ymin'):
            temp_file.write(f'#define C_ymin {float(C_bndry[2])}\n')
        elif line.startswith('#define C_ymax'):
            temp_file.write(f'#define C_ymax {float(C_bndry[3])}\n')
        elif line.startswith('#define C_zmin'):
            temp_file.write(f'#define C_zmin {float(C_bndry[4])}\n')
        elif line.startswith('#define C_zmax'):
            temp_file.write(f'#define C_zmax {float(C_bndry[5])}\n')

        elif line.startswith('#define M_xmin'):
            temp_file.write(f'#define M_xmin {float(M_bndry[0])}\n')
        elif line.startswith('#define M_xmax'):
            temp_file.write(f'#define M_xmax {float(M_bndry[1])}\n')
        elif line.startswith('#define M_ymin'):
            temp_file.write(f'#define M_ymin {float(M_bndry[2])}\n')
        elif line.startswith('#define M_ymax'):
            temp_file.write(f'#define M_ymax {float(M_bndry[3])}\n')
        elif line.startswith('#define M_zmin'):
            temp_file.write(f'#define M_zmin {float(M_bndry[4])}\n')
        elif line.startswith('#define M_zmax'):
            temp_file.write(f'#define M_zmax {float(M_bndry[5])}\n')

        elif line.startswith('#define R_xmin'):
            temp_file.write(f'#define R_xmin {float(R_bndry[0])}\n')
        elif line.startswith('#define R_xmax'):
            temp_file.write(f'#define R_xmax {float(R_bndry[1])}\n')
        elif line.startswith('#define R_ymin'):
            temp_file.write(f'#define R_ymin {float(R_bndry[2])}\n')
        elif line.startswith('#define R_ymax'):
            temp_file.write(f'#define R_ymax {float(R_bndry[3])}\n')
        elif line.startswith('#define R_zmin'):
            temp_file.write(f'#define R_zmin {float(R_bndry[4])}\n')
        elif line.startswith('#define R_zmax'):
            temp_file.write(f'#define R_zmax {float(R_bndry[5])}\n')

        elif line.startswith('int Multilevel_flag'):
            temp_file.write(f'int Multilevel_flag = {flag};\n')

        else:
            temp_file.write(line)

shutil.copyfile(temp_file_path, udf_path)  # Forcefully replace original file

os.remove(temp_file_path)  # Delete temporal file

makefile_path = os.path.join(script_dir, 'libudf', 'src', 'makefile')  # makefile path
temp_file_path = makefile_path + '.tmp'  # Temporal file for modifying

with open(makefile_path, 'r', encoding='utf-8') as make_file, open(temp_file_path, 'w', encoding='utf-8') as temp_file:
    for line in make_file:
        if line.startswith('EX_LIB'):
            temp_file.write(f'EX_LIB={script_dir}/libh5rw.so\n')
        else:
            temp_file.write(line)

shutil.copyfile(temp_file_path, makefile_path)  # Forcefully replace original file

os.remove(temp_file_path)  # Delete temporal file


# Part III. Update macro definitions in h5rw.cpp
h5rw_path = os.path.join(script_dir, 'lib_h5rw', 'src', 'h5rw.cpp')  # h5rw.cpp path
temp_file_path = h5rw_path + '.tmp'  # Temporal file for modifying

with open(h5rw_path, 'r', encoding='utf-8') as h5rw_file, open(temp_file_path, 'w', encoding='utf-8') as temp_file:
    for line in h5rw_file:

        # Mesh number definition

        if line.startswith('#define DIM0'):
            temp_file.write(f'#define DIM0 {dim[0]}\n')
        elif line.startswith('#define DIM1'):
            temp_file.write(f'#define DIM1 {dim[1]}\n')
        elif line.startswith('#define DIM2'):
            temp_file.write(f'#define DIM2 {dim[2]}\n')
        elif line.startswith('#define OUTDIMF0'):
            temp_file.write(f'#define OUTDIMF0 {out_F_dim[0]}\n')
        elif line.startswith('#define OUTDIMF1'):
            temp_file.write(f'#define OUTDIMF1 {out_F_dim[1]}\n')
        elif line.startswith('#define OUTDIMF2'):
            temp_file.write(f'#define OUTDIMF2 {out_F_dim[2]}\n')
        elif line.startswith('#define OUTDIMC0'):
            temp_file.write(f'#define OUTDIMC0 {out_C_dim[0]}\n')
        elif line.startswith('#define OUTDIMC1'):
            temp_file.write(f'#define OUTDIMC1 {out_C_dim[1]}\n')
        elif line.startswith('#define OUTDIMC2'):
            temp_file.write(f'#define OUTDIMC2 {out_C_dim[2]}\n')
        elif line.startswith('#define OUTDIMM0'):
            temp_file.write(f'#define OUTDIMM0 {out_M_dim[0]}\n')
        elif line.startswith('#define OUTDIMM1'):
            temp_file.write(f'#define OUTDIMM1 {out_M_dim[1]}\n')
        elif line.startswith('#define OUTDIMM2'):
            temp_file.write(f'#define OUTDIMM2 {out_M_dim[2]}\n')
        elif line.startswith('#define OUTDIMR0'):
            temp_file.write(f'#define OUTDIMR0 {out_R_dim[0]}\n')
        elif line.startswith('#define OUTDIMR1'):
            temp_file.write(f'#define OUTDIMR1 {out_R_dim[1]}\n')
        elif line.startswith('#define OUTDIMR2'):
            temp_file.write(f'#define OUTDIMR2 {out_R_dim[2]}\n')

        else:
            temp_file.write(line)

shutil.copyfile(temp_file_path, h5rw_path)  # Forcefully replace original file

os.remove(temp_file_path)  # Delete temporal file


# Part IV. Update mesh tally in .rmc input card
rmc_files = [filename for filename in os.listdir(script_dir) if filename.endswith('.rmc')]

if len(rmc_files) == 1:
    rmc_path = os.path.join(script_dir, rmc_files[0])  # .rmc path
    print("The only existing .rmc file found:", rmc_path)
else:
    print("No or multiple .rmc files found in the directory.")

temp_file_path = rmc_path + '.tmp'  # emporal file for modifying

with open(rmc_path, 'r', encoding='utf-8') as rmc_file, open(temp_file_path, 'w', encoding='utf-8') as temp_file:
    for line in rmc_file:

        # Mesh number definition

        if line.startswith('                    Bound ='):
            temp_file.write(f'                    Bound = {P_bndry[0]} {P_bndry[1]} {P_bndry[2]} {P_bndry[3]} {P_bndry[4]} {P_bndry[5]}\n')
        elif line.startswith('                    Scope ='):
            temp_file.write(f'                    Scope = {dim[0]} {dim[1]} {dim[2]}\n')
        else:
            temp_file.write(line)

shutil.copyfile(temp_file_path, rmc_path)  # Forcefully replace original file

os.remove(temp_file_path)  # Delete temporal file

if mode == 0:
    # Part V. Create initial .h5 files

    # Create initial .h5 files for fuel
    filename_fuel = os.path.join(script_dir, 'info_fuel.h5')
    file_fuel = h5py.File(filename_fuel, 'w')

    group_fuel = file_fuel.create_group('Geometry')
    group_fuel.attrs['MeshType'] = np.array([1], dtype=np.int32)

    bin_number_fuel = np.array([int(out_F_dim[0]), int(out_F_dim[1]), int(out_F_dim[2])])
    group_fuel.create_dataset('BinNumber', data=bin_number_fuel, dtype=np.float64)

    boundary_fuel = np.array([[float(F_bndry[0]), float(F_bndry[1])], [float(F_bndry[2]), float(F_bndry[3])], [float(F_bndry[4]), float(F_bndry[5])]])
    boundary_fuel *= 100
    group_fuel.create_dataset('Boundary', data=boundary_fuel, dtype=np.float64)

    data_fuel = np.full((int(out_F_dim[0]), int(out_F_dim[1]), int(out_F_dim[2])), float(lines[6]), dtype=np.float64)
    file_fuel.create_dataset('/temp_fuel', data=data_fuel, dtype=np.float64)

    file_fuel.create_dataset('/temp_fuel_max', data=float(lines[6]), dtype=np.float64)
    file_fuel.create_dataset('/temp_fuel_average', data=float(lines[6]), dtype=np.float64)

    if flag == 1:
        for i in range(1, 6):
            dataset_name = f'/temp_fuel{i}'
            data_fuel_multilevel = np.full((int(out_F_dim[0]), int(out_F_dim[1]), int(out_F_dim[2])), float(lines[6]), dtype=np.float64)
            file_fuel.create_dataset(dataset_name, data=data_fuel_multilevel, dtype=np.float64)

    file_fuel.close()

    # Create initial .h5 files for coolant
    filename_coolant = os.path.join(script_dir, 'info_coolant.h5')
    file_coolant = h5py.File(filename_coolant, 'w')

    group_coolant = file_coolant.create_group('Geometry')
    group_coolant.attrs['MeshType'] = np.array([1], dtype=np.int32)

    bin_number_coolant = np.array([int(out_C_dim[0]), int(out_C_dim[1]), int(out_C_dim[2])])
    group_coolant.create_dataset('BinNumber', data=bin_number_coolant, dtype=np.float64)

    boundary_coolant = np.array([[float(C_bndry[0]), float(C_bndry[1])], [float(C_bndry[2]), float(C_bndry[3])], [float(C_bndry[4]), float(C_bndry[5])]])
    boundary_coolant *= 100
    group_coolant.create_dataset('Boundary', data=boundary_coolant, dtype=np.float64)

    data_coolant_temp = np.full((int(out_C_dim[0]), int(out_C_dim[1]), int(out_C_dim[2])), float(lines[9]), dtype=np.float64)
    file_coolant.create_dataset('/temp_coolant', data=data_coolant_temp, dtype=np.float64)

    data_coolant_density = np.full((int(out_C_dim[0]), int(out_C_dim[1]), int(out_C_dim[2])), float(lines[10]) / 1000, dtype=np.float64)
    file_coolant.create_dataset('/r_coolant', data=data_coolant_density, dtype=np.float64)

    file_coolant.close()

    # Create initial .h5 files for moderator
    filename_moderator = os.path.join(script_dir, 'info_moderator.h5')
    file_moderator = h5py.File(filename_moderator, 'w')

    group_moderator = file_moderator.create_group('Geometry')
    group_moderator.attrs['MeshType'] = np.array([1], dtype=np.int32)

    bin_number_moderator = np.array([int(out_M_dim[0]), int(out_M_dim[1]), int(out_M_dim[2])])
    group_moderator.create_dataset('BinNumber', data=bin_number_moderator, dtype=np.float64)

    boundary_moderator = np.array([[float(M_bndry[0]), float(M_bndry[1])], [float(M_bndry[2]), float(M_bndry[3])], [float(M_bndry[4]), float(M_bndry[5])]])
    boundary_moderator *= 100
    group_moderator.create_dataset('Boundary', data=boundary_moderator, dtype=np.float64)

    data_moderator = np.full((int(out_M_dim[0]), int(out_M_dim[1]), int(out_M_dim[2])), float(lines[7]), dtype=np.float64)
    file_moderator.create_dataset('/temp_moderator', data=data_moderator, dtype=np.float64)

    file_moderator.close()

    # Create initial .h5 files for reflector
    filename_reflector = os.path.join(script_dir, 'info_reflector.h5')
    file_reflector = h5py.File(filename_reflector, 'w')

    group_reflector = file_reflector.create_group('Geometry')
    group_reflector.attrs['MeshType'] = np.array([1], dtype=np.int32)

    bin_number_reflector = np.array([int(out_R_dim[0]), int(out_R_dim[1]), int(out_R_dim[2])])
    group_reflector.create_dataset('BinNumber', data=bin_number_reflector, dtype=np.float64)

    boundary_reflector = np.array([[float(R_bndry[0]), float(R_bndry[1])], [float(R_bndry[2]), float(R_bndry[3])], [float(R_bndry[4]), float(R_bndry[5])]])
    boundary_reflector *= 100
    group_reflector.create_dataset('Boundary', data=boundary_moderator, dtype=np.float64)

    data_reflector = np.full((int(out_R_dim[0]), int(out_R_dim[1]), int(out_R_dim[2])), float(lines[8]), dtype=np.float64)
    file_reflector.create_dataset('/temp_reflector', data=data_reflector, dtype=np.float64)

    file_reflector.close()

    # Part VI. Save .h5 files
    source_files = [
        "info_fuel.h5",
        "info_coolant.h5",
        "info_moderator.h5",
        "info_reflector.h5",
    ]
    destination_folder = os.path.join(script_dir, "histories")
    os.makedirs(destination_folder, exist_ok=True)

    for file in source_files:
        first_dot_index = file.find('.')
        file_name, file_ext = (file[:first_dot_index], file[first_dot_index:]) if first_dot_index != -1 else (file, '')
        new_file_name = f"{file_name}_iter0{file_ext}"
        source_path = os.path.join(script_dir, file)
        destination_path = os.path.join(destination_folder, new_file_name)
        shutil.copyfile(source_path, destination_path)
