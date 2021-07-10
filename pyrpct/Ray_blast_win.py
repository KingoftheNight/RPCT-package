import ray
import subprocess
import time
import sys
import os
from Visual import visual_longcommand_raywindows
now_path = os.getcwd()
fuction_path = os.path.dirname(__file__)
# psiblast_path = os.path.join(os.path.join(fuction_path, 'bin'), 'psiblast.exe')


commands = sys.argv[1]
db_path = os.path.join(os.path.join(os.path.join(now_path, 'pyrpct'), 'blastDB'), 'nr')
with open(commands, 'r', encoding='UTF-8') as f:
    data = f.readlines()
file_box = []
out_box = []
for line in data:
    line = line.strip('\n')
    file_box.append(line.split('@@')[0])
    out_box.append(line.split('@@')[1])

ray.init()
start = time.time()


@ray.remote
def blast(file, database_path, out):
    if out not in os.listdir(os.path.split(out)[0]):
        command = visual_longcommand_raywindows('psiblast', file, database_path, out)
        subprocess.Popen(command, shell=True).communicate()


results_id = []
for i in range(len(file_box)):
    results_id.append(blast.remote(file_box[i], db_path, out_box[i]))

ray.get(results_id)
ray.shutdown()
print("等待时间: {}s".format(time.time()-start))
