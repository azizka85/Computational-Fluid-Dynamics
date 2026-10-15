from utils.fs import ensure_path

def write_data(w, t, dx, dy, m, nx, ny, out_dir):
    file_path = f'{out_dir}/data/data.{m:03}.vtk'

    ensure_path(file_path)

    with open(file_path, 'w') as file:
        file.write('# vtk DataFile Version 3.0\n')
        file.write(f'TIME {t:.3f}\n')
        file.write('ASCII\n')
        file.write('DATASET STRUCTURED_GRID\n')
        file.write(f'DIMENSIONS {nx} {ny} 1\n')
        file.write(f'POINTS {nx*ny} float\n')

        for j in range(ny):
            for i in range(nx):
                x = dx * i
                y = dy * j

                file.write(f'{x:.3f} {y:.3f} 0.0\n')

        file.write('FIELD FieldData 1\n')
        file.write('Time 1 1 float\n')
        file.write(f'{t:.3f}\n')
        file.write(f'POINT_DATA {nx*ny}\n')
        file.write('SCALARS w float\n')
        file.write('LOOKUP_TABLE default\n')

        for j in range(ny):
            for i in range(nx):
                p = j + i*ny

                file.write(f'{w[p]:.3f}\n')
