#!/usr/bin/env python3

import os
import sys
import requests

URL_DOMAIN = "api.icfpcontest.com"
BASE_DIR = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
DATA_DIR = f"{BASE_DIR}/data/problems"


def command_problem(ids):
    # problem?problem_id=[problem_id:u32]
    if not isinstance(ids, list):
        ids = [ids]
    ids = [int(x) for x in ids]
    print(BASE_DIR)

    for id in ids:
        filepath = f"{DATA_DIR}/problem-{id}.json"
        if os.path.isfile(filepath):
            continue

        url = f"https://{URL_DOMAIN}/problem?problem_id={id}"
        print(f'Downloading {url}', file=sys.stderr)
        response = requests.get(url)
        data = response.text
        data = data.replace("\\n", '\n')
        data = data.replace('\\"', '"')
        with open(filepath, 'w') as f:
            f.write(data)

    return True


def main():
    subcommand = "" if len(sys.argv) < 2 else sys.argv[1]
    if subcommand == "problem":
        if command_problem(sys.argv[2:]):
            return

    print("Usage: " + sys.argv[0] + " subcommand", file=sys.stderr)


if __name__ == '__main__':
    main()
