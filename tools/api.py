#!/usr/bin/env python3

import re
import os
import sys
import json
import requests
import dateutil.parser

URL_DOMAIN = "api.icfpcontest.com"
BASE_DIR = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
DATA_DIR = f"{BASE_DIR}/data"
PROBLEMS_DIR = f"{DATA_DIR}/problems"
SCOREBOARD_DIR = f"{DATA_DIR}/scoreboard"

# Do not share the token with others during the contest.
TOKEN = 'eyJ0eXAiOiJKV1QiLCJhbGciOiJIUzI1NiJ9.eyJ1aWQiOiI2NGEzNjZmMDhjNjg1MzEzZDFjNjBkOGYiLCJpYXQiOjE2ODg4MDA0MjYsImV4cCI6MTY4OTgwMDQyNn0.ppnYGzlvv6GSsDNQUuvmVHdjoaCd46jXt860E0OqdUs'

USAGE = """
$ api.py <command> [other_arguments]

  where
    <command> := "problem" | "problems" | "scoreboard"

  Commands:
    problem:     $ api.py problem <problem_id>
      Downloads JSON file for a problem <problem_id>. <problem_id> is
      a number. <problem_id> can be a list of IDs.

    problems:    $ api.py problems
      Displays the number of available problems, and downloads all problems
      you don't have in your local problems directory.

    scoreboard:  $ api.py scoreboard
      Downloads and saves the current scoreboard.

    submit:      $ api.py submit <json_file | directory>
      Submits solution(s) to the official server. A name of JSON file
      should be either of "problem-XX.json" or "problem-XX_yyyy.json"
      where XX represents its problem ID.
"""


def command_problem(ids):
    # problem?problem_id=[problem_id:u32]
    if not isinstance(ids, list):
        ids = [ids]
    ids = [int(x) for x in ids]

    for id in ids:
        filepath = f"{PROBLEMS_DIR}/problem-{id}.json"
        if os.path.isfile(filepath):
            continue

        url = f"https://{URL_DOMAIN}/problem?problem_id={id}"
        print(f'Downloading {url}', file=sys.stderr)
        response = requests.get(url)
        data = response.text
        data = json.loads(data)
        data = data['Success']
        data = data.replace("\\n", '\n')
        data = data.replace('\\"', '"')
        # print(data)
        with open(filepath, 'w') as f:
            f.write(data)

    return True


def command_problems(options):
    # problems
    # api.py problems --get # get new problems
    url = f"https://{URL_DOMAIN}/problems"
    response = requests.get(url)
    data = response.text
    data = data.replace("\\n", '\n')
    data = data.replace('\\"', '"')
    data = json.loads(data)
    number = data['number_of_problems']
    print(f"Number of problems: {number}")

    current = sum(os.path.isfile(os.path.join(PROBLEMS_DIR, name))
                  for name in os.listdir(PROBLEMS_DIR))
    if current < number:
        print(
            f"Currently we have only {current} problems. Download unavailable problems.")
        command_problem(list(range(1, number + 1)))

    return True


def command_scoreboard(_options):
    # scoreboard
    url = f"https://{URL_DOMAIN}/scoreboard"
    response = requests.get(url)
    data = response.text
    data = data.replace("\\n", '\n')
    data = data.replace('\\"', '"')
    data = json.loads(data)

    updated_at = dateutil.parser.parse(data['updated_at'])
    start_at = dateutil.parser.parse('2023-07-07T12:00:00Z')
    duration = (updated_at - start_at).seconds // 60
    filepath = f"{SCOREBOARD_DIR}/{duration:05d}.json"
    with open(filepath, 'w') as f:
        json.dump(data, f)
    print(f"Save scoreboard {duration:05d}.json at {updated_at}.")

    # Dump some data
    if data['frozen']:
        print("** Frozen **")

    print("Top 10 scores are:")
    for (rank, team) in enumerate(data['scoreboard']):
        name = team['username']
        if rank >= 10 and name != 'sanma':
            continue

        if len(name) > 20:
            name = name[:17] + '...'
        score = team['score']
        print("{:2d} : {:20s} : {:13.1f}".format(rank+1, name, score))

    return True


def command_submit(path):
    filepaths = path_to_filepaths(path)

    for filepath in filepaths:
        filename = os.path.basename(filepath)
        m = re.search('[^\\d]-(\\d+)([-_\\d]*)\\.json', filename)
        if not m:
            print(f"Can't parse filename: {filename}. {m}")
            continue
        id = int(m[1])
        if id < 1 or id > 999:
            print(f'Invalid problem ID ({id}) for {filename}')
            continue
        data = {"problem_id": id}
        with open(filepath, 'r') as f:
            content = json.load(f)
            data["contents"] = json.dumps(content)
        data = json.dumps(data)

        # Submit with a POST method
        url = f"https://{URL_DOMAIN}/submission"
        headers = {'Content-Type': 'application/json',
                   'Authorization': f'Bearer {TOKEN}'}
        response = requests.post(url=url, headers=headers, data=data)
        print(f"Submit {filename} for {id}: {response.text}", file=sys.stderr)
    return True


def path_to_filepaths(input_path):
    paths = [input_path] if isinstance(input_path, str) else input_path
    for path in paths:
        if os.path.isfile(path):
            yield path
        else:
            names = os.listdir(path)
            for name in names:
                yield os.path.join(path, name)


def main():
    subcommand = "" if len(sys.argv) < 2 else sys.argv[1]
    if subcommand == "problem":
        if command_problem(sys.argv[2:]):
            return
    elif subcommand == "problems":
        if command_problems(sys.argv[2:]):
            return
    elif subcommand == "scoreboard":
        if command_scoreboard(sys.argv[2:]):
            return
    elif subcommand == "submit":
        if command_submit(sys.argv[2:]):
            return

    print("Usage: " + sys.argv[0] + " subcommand", file=sys.stderr)


if __name__ == '__main__':
    main()
