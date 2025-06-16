#!/usr/bin/env python
import subprocess
import time
import sys

if len(sys.argv) != 2:
    print("Usage: python compare_speed.py <input_file>")
    sys.exit(1)

input_file = sys.argv[1]

start = time.time()
subprocess.run(["./build/msa", input_file], check=True, stdout=subprocess.DEVNULL)
end = time.time()
my_time_UPGMA = end - start
print(f"Moj algoritam, UPGMA: {my_time_UPGMA:.2f} sekundi")

start = time.time()
subprocess.run(["./build/msa", "-n", input_file], check=True, stdout=subprocess.DEVNULL)
end = time.time()
my_time_NJ = end - start
print(f"Moj algoritam, NJ: {my_time_NJ:.2f} sekundi")

start = time.time()
subprocess.run(["clustalw", input_file], check=True, stdout=subprocess.DEVNULL)
end = time.time()
clustal_time = end - start
print(f"ClustalW: {clustal_time:.2f} sekundi")
