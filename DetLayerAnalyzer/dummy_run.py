import os
import csv
import subprocess
import sys

file = sys.argv[1]
cmd = []
with open(file) as f:
    reader = csv.reader(f,delimiter=",")
    for row in reader:
        cmd.append("./run.sh %s %s %s %s %s %s %s"%(row[0],row[1],row[2],row[3],row[4],row[5], row[6]))
        print(row)
#subprocess.run(["ls"],shell=True)
for c in cmd:
    subprocess.run(c,shell=True)

