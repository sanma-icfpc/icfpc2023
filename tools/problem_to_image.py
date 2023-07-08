#!/usr/bin/env python3

import os
import glob
import argparse
import subprocess


def problem_to_image(json_path, png_path):
    subprocess.run(
        ['./solver', 'problem-to-png', json_path, png_path],
        cwd='../src', check=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('json_paths', nargs='*')
    parser.add_argument('--output-dir', default='../data/problems-image')
    args = parser.parse_args()

    if len(args.json_paths) == 0:
        args.json_paths = list(glob.glob('../data/problems/problem-*.json'))

    os.makedirs(args.output_dir, exist_ok=True)

    for json_path in args.json_paths:
        png_path = os.path.join(args.output_dir, os.path.splitext(
            os.path.basename(json_path))[0] + '.png')
        problem_to_image(json_path, png_path)


if __name__ == '__main__':
    main()
