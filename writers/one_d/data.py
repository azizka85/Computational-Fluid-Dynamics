from utils.fs import ensure_path


def write_data(u, t, dy, m, out_dir):
    file_path = f'{out_dir}/data/data.{m:03}.vtk'

    ensure_path(file_path)

    ny = len(u)

    with open(file_path, 'w') as file:
        file.write('# vtk DataFile Version 3.0\n')
        file.write(f'TIME {t:.3f}\n')
        file.write('ASCII\n')
        file.write('DATASET STRUCTURED_GRID\n')
        file.write(f'DIMENSIONS 1 {ny} 1\n')
        file.write(f'POINTS {ny} float\n')

        file.writelines(f'0.0 {(i*dy):.3f} 0.0\n' for i in range(ny))

        file.write('FIELD FieldData 1\n')
        file.write('Time 1 1 float\n')
        file.write(f'{t:.3f}\n')
        file.write(f'POINT_DATA {ny}\n')
        file.write('SCALARS u float\n')
        file.write('LOOKUP_TABLE default\n')

        file.writelines(f'{u[i]:.3f}\n' for i in range(ny))

def write_statistics(stats, out_dir):
    file_path = f'{out_dir}/convergence.csv'

    ensure_path(file_path)

    with open(file_path, 'w') as file:
        file.write('n,  t,   max_diff\n')

        file.writelines(f'{stat[0]}, {stat[1]:.3f}, {stat[2]:.5f}\n' for stat in stats)
