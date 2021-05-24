from sys import argv

keys = argv[1]
key_file = open(keys, 'r')
key_list = []
for line in key_file:
    key_list.append(line.strip()[:-4])

key_file.close()

key_file = open(keys, 'w')
for key in key_list:
    key_file.write(key)
    key_file.write('\n')

key_file.close()
