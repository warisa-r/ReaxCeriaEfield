
#Function to extract lattice constant from dump file
def lattice_constant_dump(file_path):
    lx_values = []
    ly_values = []
    lz_values = []
    timesteps = []

    with open(file_path, 'r') as file:
        lx = ly = lz = None
        timestep = None
        
        for line in file:
            if "ITEM: TIMESTEP" in line:
                timestep = int(file.readline().strip())
                timesteps.append(timestep)

            elif "ITEM: BOX BOUNDS" in line:
                lx_line = file.readline().strip().split()
                ly_line = file.readline().strip().split()
                lz_line = file.readline().strip().split()

                lx = float(lx_line[1]) - float(lx_line[0])
                ly = float(ly_line[1]) - float(ly_line[0])
                lz = float(lz_line[1]) - float(lz_line[0])

                lx_values.append(lx)
                ly_values.append(ly)
                lz_values.append(lz)

    return timesteps, lx_values, ly_values, lz_values

#Function to extract lattice constant from log file
def lattice_constant_log(file_path):
    lx_values = []
    ly_values = []
    lz_values = []
    timesteps = []

    with open(file_path, 'r') as file:
        data_section = False
        for line in file:
            if 'Step' in line:
                data_section = True
                continue  
            
            if not data_section or line.strip() == '' or line.startswith('Loop time') or line.startswith('Performance'):
                data_section = False
                continue
            
            parts = line.split()
            
            if len(parts) >= 6:
                try:
                    timestep = int(parts[0])
                    lx = float(parts[3])
                    ly = float(parts[4])
                    lz = float(parts[5])
                    
                    timesteps.append(timestep)
                    lx_values.append(lx)
                    ly_values.append(ly)
                    lz_values.append(lz)
                except ValueError as e:
                    print(f"ValueError: {e} for line: {line.strip()}")
                    continue
    
    return timesteps, lx_values, ly_values, lz_values
