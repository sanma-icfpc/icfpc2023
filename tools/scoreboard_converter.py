#!/usr/bin/env python3

import json
import os
import sys
import time

SCOREBOARD_DIR = 'data/scoreboard/'
DOCS_DIR = 'docs/'
START_UNIXTIME = time.mktime(time.strptime('Fri Jul 7 12:00:00 2023'))


def time2min(timestamp):
    timestamp = timestamp[:timestamp.rfind('.')]
    unixtime = time.mktime(time.strptime(timestamp, "%Y-%m-%dT%H:%M:%S"))
    return int(unixtime - START_UNIXTIME) // 60


def main():
    scoreboard_dir = SCOREBOARD_DIR
    if len(sys.argv) >= 2:
        scoreboard_dir = sys.argv[1]
    scoreboard_dir = os.path.abspath(scoreboard_dir)
    scoreboards = dict()
    for filename in os.listdir(scoreboard_dir):
        filepath = os.path.join(scoreboard_dir, filename)
        try:
            with open(filepath) as f:
                print(filepath)
                j = json.load(f)
                elapsed = time2min(j['updated_at'])
                scoreboards[elapsed] = j
        except:
            print('Fail to parse ' + filename)

    teams = set()
    for snapshot in scoreboards.values():
        names = map(lambda x: x['username'], snapshot['scoreboard'])
        teams |= set(names)
    timestamps = sorted(scoreboards.keys())

    scoreboard2 = dict()
    for team in teams:
        scoreboard2[team] = list()
    for minutes in timestamps:
        scoreboard = scoreboards[minutes]
        for team in scoreboard['scoreboard']:
            name = team['username']
            score = team['score']
            scoreboard2[name].append({'x': minutes, 'y': score})

    scoreboard_filepath = os.path.join(DOCS_DIR, 'scoreboard.json')
    with open(scoreboard_filepath, 'w') as f:
        json.dump(scoreboard2, f)

    # Dump top 10 teams on frozen
    last_score = scoreboards[timestamps[-1]]['scoreboard']
    last_score = sorted(last_score, key=lambda x: x['score'], reverse=True)
    print(list(map(lambda x: x['username'], last_score[:10])))


if __name__ == '__main__':
    main()
