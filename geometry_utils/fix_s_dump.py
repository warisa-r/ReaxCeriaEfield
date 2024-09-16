def fix_s_lammpstrj(input_file, output_file):
    # modify the lammpstrj with ss condition, make it visualizable in vmd, by adding vacuum to everysteps.

    with open(input_file,'r') as file:
        lines = file.readlines()

    _temp_line = lines[7]
    _item = [float(num) for num in _temp_line.split()]
    vacuum = _item[0]

    index = 0
    for i, line in enumerate(lines):
        if 'BOX' in line:
            temp_line = lines[i+3]
            item = [float(num) for num in temp_line.split()]
            item[0] -= vacuum
            lines[i+3] = f'{item[0]} {item[1]} {item[2]} \n'
    with open(output_file,'w') as file:
        file.writelines(lines)                   
