from argparse import ArgumentParser

import pandas as pd

def get_manifest(manifest_path: str, info) -> str:
    all_files = set(sum([list(info[col]) for col in info.columns.values if col != 'total_size'], []))

    manifest_lines = []
    with open(manifest_path, 'r') as fp:
        for i, line in enumerate(fp):
            if i == 0 or line.split('\t')[1] in all_files:
                manifest_lines.append(line)

    return ''.join(manifest_lines)

if __name__ == '__main__':
    parser = ArgumentParser('get_manifest', description='Get GDC manifest file for info index')
    parser.add_argument('manifest', type=str)
    parser.add_argument('info', type=str)
    parser.add_argument('--case_id', type=str)
    parser.add_argument('--out', type=str, help='Filtered manifest output')

    args = parser.parse_args()

    info = pd.read_csv(args.info, index_col='case_id')
    if args.case_id:
        info = info.loc[[args.case_id]]
    
    manifest = get_manifest(args.manifest, info)

    if args.out:
        with open(args.out, 'w') as fp:
            fp.write(manifest)
    else:
        print(manifest, end='')
