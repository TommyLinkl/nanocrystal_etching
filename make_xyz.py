import numpy as np
import sys

# Takes config file output from nanocrystal class and converts it to an xyz file for VMD

direc = "TestFiles/"
trial_number = 0
temperature = 0.1
mu = float(sys.argv[1])
radius = 20.0
config_file = direc + "t%dctest%.2f_%.2f_%.2f.txt" % (trial_number, temperature, mu, radius)

f = open(config_file, 'r')
g = open(config_file.replace(".txt", ".xyz"), "w")

count_line = 0
max_particles = 75000
count_particles = 0
step_number = 0
g.write(str(max_particles) + "\n")
g.write("\n")

for line in f:
    if count_line >= 3:
    
        line = line.strip()
        line_list = line.split()

        if len(line_list) == 4:
            for i in range(max_particles - count_particles):
                g.write("1 0 0 0\n")
            time = line[2]
            g.write(str(max_particles) + "\n")
            g.write("\n")
            
            count_particles = 0
            step_number += 1
            print(step_number)
                
        else:
            g.write(line + "\n")
            count_particles += 1
    count_line += 1

f.close()
g.close()
