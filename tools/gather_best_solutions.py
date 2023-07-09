#!/usr/bin/env python3

import os
import re
import glob
import json
import shutil
import argparse
import subprocess

BASE_DIR = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
DATA_DIR = os.path.join(BASE_DIR, 'data')
SOLUTIONS_DIR = os.path.join(DATA_DIR, 'solutions')


def eval_solution(solution_path):
    if os.name == 'posix':
        cwd = os.path.join(BASE_DIR, 'src')
        solver_path = './solver'
    elif os.name == 'nt':
        # FileNotFoundError...
        cwd = os.path.join(BASE_DIR, 'vs')
        solver_path = os.path.join('x64', 'Release', 'solver.exe')
    cmds = [solver_path, 'eval', solution_path]
    res = subprocess.run(cmds, cwd=cwd, check=True, capture_output=True, text=True)
    return int(res.stdout.split('=')[1].strip().split(' ')[0])

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('target_dir', nargs='?', type=str, default=os.path.join(SOLUTIONS_DIR, 'bests'))
    parser.add_argument('--solutions-dir', type=str, default=SOLUTIONS_DIR)
    parser.add_argument('--set-constant-volumes', type=float, default=-1.0)
    args = parser.parse_args()

    assert os.path.isdir(args.solutions_dir)

    solution_pattern = os.path.join(args.solutions_dir, '**', 'solution-*.json')

    id_to_paths = {}

    for filepath in glob.glob(solution_pattern, recursive=True):
        filename = os.path.basename(filepath)
        m = re.search(r'solution-([0-9]+)(_|-)?.*.json', filename)
        if not m:
            print(f"Can't parse filename: {filename}. {m}")
            continue
        id = int(m[1])
        if id < 1:
            print(f'Invalid problem ID ({id}) for {filename}')
            continue
        if id not in id_to_paths:
            id_to_paths[id] = []
        id_to_paths[id].append(filepath)
    
    id_to_best_score = {}
    id_to_best_path = {}

    # compute_score 自体が並列化されているのでそのまま
    for id, paths in id_to_paths.items():
        if id not in id_to_best_score:
            id_to_best_score[id] = 0
        for path in paths:
            score = eval_solution(path)
            if id_to_best_score[id] < score:
                id_to_best_score[id] = score
                id_to_best_path[id] = path
        print(id, id_to_best_score[id], id_to_best_path[id])
    
    if not os.path.exists(args.target_dir):
        os.makedirs(args.target_dir)

    for id, path in id_to_best_path.items():
        basename = f'solution-{id}.json'
        target_path = os.path.join(args.target_dir, basename)
        if args.set_constant_volumes < 0:
            shutil.copy2(path, target_path)
        else:
            with open(path, 'r', encoding='utf-8') as f:
                solution = json.load(f)
            num_musicians = len(solution['placements'])
            solution['volumes'] = [args.set_constant_volumes] * num_musicians
            with open(target_path, 'w', encoding='utf-8') as f:
                json.dump(solution, f, indent=4)


if __name__ == '__main__':
    main()
