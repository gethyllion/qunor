
import os, argparse
argus = argparse.ArgumentParser()
argus.add_argument('files','-f',)
file = [i for i in os.listdir() if i.endswith('.txt') ]
try:

    data = open(file[0]).read()

except Exception as e:
    print(e)

data = data.split('\n')


for i in data:
	if i =='':
		data.remove(i)
		
print(data)

with open('task_q.txt','w') as f:
    f.write(repr(data))
