import os
import glob
import tqdm
import argparse
import subprocess

def solve(json_path, solve_command):
    try:
        subprocess.run(
            solve_command.format(input=json_path),
            cwd='../src', shell=True, check=True)
    except subprocess.CalledProcessError:
        pass

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('solve_command')
    parser.add_argument('json_paths', nargs='*')
    args = parser.parse_args()

    if len(args.json_paths) == 0:
        args.json_paths = list(glob.glob('../data/problems/problem-*.json'))

    for json_path in tqdm.tqdm(args.json_paths):
        solve(json_path, args.solve_command)

if __name__ == '__main__':
    main()